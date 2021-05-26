
import pickle
import numpy as np
from multiprocessing import Pool
from matplotlib import pyplot as plt

class ThresholdSelect:
    def __init__(self, features_pkl_path):
        self.predictions = pickle.load(open(features_pkl_path, 'rb'))

    def draw(self):
        scores = [x for x in self.scores if type(x)==tuple]
        jacc, recall, precision, fscore = zip(*self.scores)
        thresholds = thresholds[:len(self.scores)]

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

    def explore_thresholds(self, thresholds):
        scores = []
        with Pool(8) as p:
            scores = p.map(self.compute_metrics_at_thres, thresholds)
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
        recall    = compute_recall(HCIP, ref_edges_set)
        precision = compute_precision(HCIP, ref_edges_set)
        if recall + precision == 0 :
            return 0.
        return (2 * recall * precision) / (recall + precision)

    def compute_metrics_at_thres(self, thres):
        HCIP = []
        for pred in bp_predictions_9d:
            batch, bait, prey, target, pnoint_pint = pred
            if bait==prey:continue
            if pnoint_pint[1] < thres: continue # select only intercations over pint threshold of input thres
            HCIP.append(pred)
        
        HCIP_net = []
        for batch, bait, prey, target, pnoint_pint in HCIP:
            HCIP_net.append((bait, prey, pnoint_pint[1]))
        
        G_o_noalt = nx.Graph()
        G_o_noalt.add_edges_from([(x[0], x[1]) for x in HCIP_net if not is_alt(x[1])])
        G_o_noalt.remove_edges_from(nx.selfloop_edges(G_o_noalt))
        G_o_edges = set(frozenset(e) for e in G_o_noalt.edges)
        intersection = G_o_edges.intersection(G_b_edges)
        union = G_o_edges.union(G_b_edges)
        
        jacc      = len(intersection)/len(union)
        recall    = compute_recall(G_o_edges, G_b_edges)
        precision = compute_precision(G_o_edges, G_b_edges)
        fscore    = compute_f1(G_o_edges, G_b_edges)
        
        return jacc, recall, precision, fscore