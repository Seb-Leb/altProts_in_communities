
import pickle
import csv
import numpy as np
import networkx as nx
from multiprocessing import Pool
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles

from altProts_in_communities.utils import is_alt

def get_bioplex_network(baits=None):
    G_b = nx.Graph()
    with open('BioPlex_interactionList_v4a.tsv', 'r') as f:
        for n,l in enumerate(f):
            ls = l.strip().split('\t')
            if n==0:
                keys = ls
                continue
            line = dict(zip(keys, ls))
            if baits is not None:
                if line['SymbolA'] in baits or line['SymbolB'] in baits:
                    G_b.add_edge(line['SymbolA'], line['SymbolB'])
            else:
                G_b.add_edge(line['SymbolA'], line['SymbolB'])

    G_b_edges = set()
    for e in G_b.edges:
        G_b_edges.add(frozenset(e))
    return G_b, G_b_edges

G_b, G_b_edges = get_bioplex_network()

def get_BP3_network(file_path):
    bp3_net = set()
    with open(file_path, 'r') as f:
        reader = csv.reader(f, delimiter='\t', quotechar='"')
        for n,l in enumerate(reader):
            if n==0:
                keys = l
                continue
            line = dict(zip(keys, l))
            bp3_net.add(frozenset((line['SymbolA'], line['SymbolB'])))
    return bp3_net

bp3_net = get_BP3_network('BioPlex_293T_Network_10K_Dec_2019.tsv')

class ThresholdSelect:
    def __init__(self, features_pkl_path, BP_edges, thresholds):
        self.predictions = pickle.load(open(features_pkl_path, 'rb'))
        self.BP_edges = BP_edges # set of frozensets
        self.thresholds = thresholds

    def draw(self):
        scores = [x for x in self.scores if type(x)==tuple]
        jacc, recall, precision, fscore = zip(*self.scores)
        thresholds = self.thresholds[:len(self.scores)]

        fig, ax = plt.subplots(figsize=(5,5))
        ax.plot(thresholds, fscore, label='Fscore', color='#7A3350')
        ax.plot(thresholds, precision,  label='precision', color='k')
        ax.plot(thresholds, recall, label='recall', color='#093E63')
        ax.plot(thresholds, jacc, label='Jaccard', color='#E4A737')
        ax.set_xlabel('NB score threshold')
        ax.set_ylabel('scores')
        plt.axvline(0.045, c='gray', linestyle='dashed')
        plt.legend()
        plt.show()

    def explore_thresholds(self):
        scores = []
        with Pool(8) as p:
            scores = p.map(self.compute_metrics_at_thres, self.thresholds)
        self.scores = scores

    def compute_recall(self, HCIP, ref_edges_set):
        TP = HCIP.intersection(ref_edges_set)
        FN = ref_edges_set.difference(HCIP)
        return len(TP)/(len(TP) + len(FN))

    def compute_precision(self, HCIP, ref_edges_set):
        TP = HCIP.intersection(ref_edges_set)
        FP = HCIP.difference(ref_edges_set)
        if len(TP) + len(FP) == 0:
            return 0.
        return len(TP)/(len(TP) + len(FP))

    def compute_f1(self, HCIP, ref_edges_set):
        recall    = self.compute_recall(HCIP, ref_edges_set)
        precision = self.compute_precision(HCIP, ref_edges_set)
        if recall + precision == 0 :
            return 0.
        return (2 * recall * precision) / (recall + precision)

    def compute_metrics_at_thres(self, thres):
        HCIP = []
        for pred in self.predictions:
            batch, bait, prey, target, pnoint_pint = pred
            if bait==prey:continue
            if pnoint_pint[1] < thres: continue # select only interactions over pint threshold of input thres
            HCIP.append(pred)
        
        HCIP_net = []
        for batch, bait, prey, target, pnoint_pint in HCIP:
            HCIP_net.append((bait, prey, pnoint_pint[1]))
        
        G_o_noalt = nx.Graph()
        G_o_noalt.add_edges_from([(x[0], x[1]) for x in HCIP_net if not is_alt(x[1])])
        G_o_noalt.remove_edges_from(nx.selfloop_edges(G_o_noalt))
        G_o_edges = set(frozenset(e) for e in G_o_noalt.edges)
        intersection = G_o_edges.intersection(self.BP_edges)
        union = G_o_edges.union(self.BP_edges)
        
        jacc      = len(intersection)/len(union)
        recall    = self.compute_recall(G_o_edges, self.BP_edges)
        precision = self.compute_precision(G_o_edges, self.BP_edges)
        fscore    = self.compute_f1(G_o_edges, self.BP_edges)
        
        return jacc, recall, precision, fscore

def tt_type_node_attrs(G, HCIP_baits):    
    node_attrs = {}
    for n in G.nodes():
        node_attrs[n] = {}
        if is_alt(n):
            node_attrs[n]['tt_type'] = 'alt'
            node_attrs[n]['node_type'] = 'prey'
        else:
            node_attrs[n]['tt_type'] = 'ref'
            if n in HCIP_baits:
                node_attrs[n]['node_type'] = 'bait'
            else:
                node_attrs[n]['node_type'] = 'prey'

    return node_attrs