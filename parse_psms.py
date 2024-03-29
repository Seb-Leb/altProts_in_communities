import os
import re
import pickle
import csv
import sys
maxInt = sys.maxsize

while True:
    # decrease the maxInt value by factor 10 
    # as long as the OverflowError occurs.
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)

from collections import Counter, namedtuple
import numpy as np
from tqdm import tqdm
import pandas as pd
import itertools as itt

from multiprocessing import Pool, Process, Queue
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles

from Bio import SeqIO
from altProts_in_communities.hierarchical_report_parser import *
from altProts_in_communities.utils import *

# Make a list of ref seq
refseqs = []
for line in SeqIO.parse("human-openprot-r1_6-refprots-+uniprot2019_03_01.fasta", "fasta"):
    refseqs.append(str(line.seq))

# paths for all MS reports
if os.path.exists('bioplex_PSreport_paths.pkl'):
    exp_report_paths = pickle.load(open('bioplex_PSreport_paths.pkl', 'rb'))
else:
    exp_report_paths = {}
    search_paths = ['/home/xroucou_group/analysis/mass_spec/reanalysis_r1.6/BioPlex_1/', '/home/xroucou_group/analysis/mass_spec/reanalysis_r1.6/BioPlex_2/']
    for search_path in search_paths:
        for dirpath, subdirs, files in os.walk(search_path):
            for f in files:
                if 'BioPlex' in dirpath and 'Default_Hierarchical_Report' in f:
                    exp_id = f.split('_')[0]
                    exp_report_paths[exp_id] = os.path.join(dirpath, f)
    pickle.dump(exp_report_paths, open('bioplex_PSreport_paths.pkl', 'wb'))

def get_protgrp_psms(fpath):
    report = parse(fpath)
    prot_grps = {}
    for protein_grp in report:
        
        prot_grp_acc = protein_grp.accession.split('|')[0]
        for peptide in protein_grp.peptide_rows:
            prot_grp_id = protein_grp.row[0][0] # re-define for every pep
            pep_seq = peptide.row[2]
            prot_accs_genes = []
            for p in peptide.peptide_names:
                fasta_header = parse_fasta_header(p)
                if 'PA' in fasta_header:
                    prot_accs_genes.append((fasta_header['PA'], prot_gene_dict[fasta_header['PA']]))
                else:
                    prot_accs_genes.append((prot_grp_acc, prot_gene_dict[prot_grp_acc]))
            
            tt_type_isalt = [is_alt(acc) for acc, gene in prot_accs_genes]

            if any(tt_type_isalt) and not all(tt_type_isalt):                                         # if any but not all members of pep_grp are alt,
                prot_accs_genes = list(itt.compress(prot_accs_genes, [not x for x in tt_type_isalt])) # remove them 
            
            protgrp_isalt = False
            if all(tt_type_isalt):
                protgrp_isalt = True
                if 'alt' not in prot_grp_id:
                    prot_grp_id = prot_grp_id+'_alt'

            for psm in peptide.psm_rows:
                if 'Confident' not in psm.validation : continue
                if protgrp_isalt:
                    if prot_grp_id not in prot_grps:
                        prot_grps[prot_grp_id] = {
                                                  'prot_accs_genes':set(prot_accs_genes), 
                                                  'psms':[(pep_seq, psm.spectrum_title)],
                                                 }
                    else:
                        prot_grps[prot_grp_id]['prot_accs_genes'] = prot_grps[prot_grp_id]['prot_accs_genes'].intersection(set(prot_accs_genes))
                        prot_grps[prot_grp_id]['psms'].append((pep_seq, psm.spectrum_title))
                else:
                    if prot_grp_id not in prot_grps:
                        prot_grps[prot_grp_id] = {
                                                  'prot_accs_genes':set(prot_accs_genes), 
                                                  'psms':[(pep_seq, psm.spectrum_title)],
                                                 }
                    else:
                        prot_grps[prot_grp_id]['prot_accs_genes'] = prot_grps[prot_grp_id]['prot_accs_genes'].intersection(set(prot_accs_genes))
                        prot_grps[prot_grp_id]['psms'].append((pep_seq, psm.spectrum_title))
    return prot_grps

def remove_invalid_grp(protgrp_psms, valid_prot_accs):
    item_ls = protgrp_psms.copy().items()
    for grp_id, grp in item_ls:
        if len(grp['prot_accs_genes'].intersection(valid_prot_accs))==0:
            del protgrp_psms[grp_id]
    return protgrp_psms

def get_single_gene_psms(protgrp_psms):
    single_gene_psms = {}
    for grp_id, grp in protgrp_psms.items():
        prot_accs, genes = zip(*grp['prot_accs_genes'])
        if len(set(genes)) == 1:
            gene = genes[0]
            if all(is_alt(x) for x in prot_accs):
                alt_acc = sorted(prot_accs)[0]
                gene = '|'.join([gene, alt_acc])
            single_gene_psms[grp_id] = {'gene':gene, 'psms':grp['psms']}
            
    # group single gene prot groups under gene
    psms = {}
    single_genes = set()
    for gene, psms_grp in itt.groupby(sorted(single_gene_psms.values(), key=lambda x: x['gene']), key=lambda x: x['gene']):
        psms_ls = [x for y in [x['psms'] for x in psms_grp] for x in y]
        n_upep = len(set(x[0] for x in psms_ls))
        if n_upep>1 or is_alt(gene):   # check for minimum of 2 unique peptide if not alt
            psms[gene] = {'psms':psms_ls, 'n_upep':n_upep}
            single_genes.add(gene)
    return psms, single_genes, set(single_gene_psms.keys())

def get_multigene_psms(protgrp_psms, bait):
    rep_psms, singles_genes, singles_grpids = get_single_gene_psms(protgrp_psms)
    for grp_id, grp in protgrp_psms.items():
        if grp_id in singles_grpids: continue
        prot_accs, genes = zip(*grp['prot_accs_genes']) 
        already_detected = set(genes).intersection(singles_genes)
        if len(already_detected)>=1:
            for gene in list(already_detected):
                rep_psms[gene]['psms'].extend(grp['psms'])

    return rep_psms

def clean_all_grps(protgrp_psms, valid_prot_accs):
    for grp_id, prot_grp in protgrp_psms.items():
        prot_accs_genes = prot_grp['prot_accs_genes'].intersection(valid_prot_accs)
        if len(prot_accs_genes)>0:
            protgrp_psms[grp_id]['prot_accs_genes'] = prot_accs_genes
        else: # prot group is not valid
            del protgrp_psms[grp_id]
    return protgrp_psms

def get_all_psms(bait):
    psms = {
        bait:{
            'exps':{exp_id:{'protgrp_psms':get_protgrp_psms(exp_report_paths[exp_id])} for exp_id in bait_rep_dict[bait] if exp_id in exp_report_paths},
            'valid_prot_accs':set(),
        }
    }
    for exp_id, exp in psms[bait]['exps'].items(): # Obtain a set of all proteins detected for each rep, take the intersection
        exp['all_prots'] = set(x for y in [list(x['prot_accs_genes']) for x in exp['protgrp_psms'].values()] for x in y)
        
    valid_prot_accs = set.intersection(*[x['all_prots'] for x in psms[bait]['exps'].values()])
    psms[bait]['valid_prot_accs'] = valid_prot_accs
    
    for exp_id, exp in psms[bait]['exps'].items(): # Remove protein groups containing solely proteins seen in only 1 rep
        exp['protgrp_psms_valid'] = remove_invalid_grp(exp['protgrp_psms'], valid_prot_accs)

    for exp_id, exp in psms[bait]['exps'].items():
        exp['protgrp_psms_valid_clean'] = clean_all_grps(exp['protgrp_psms_valid'], valid_prot_accs)
            
    for exp_id, exp in psms[bait]['exps'].items(): # assign PSMs to prot groups, single gene and multigene
        exp['psms'] = get_multigene_psms(exp['protgrp_psms_valid_clean'], bait)
        
    #pickle.dump(psms, open('bioplex_psms/{}_psms.pkl'.format(bait), 'wb'))
    return psms

def get_HPP_peps(alt_acc, peps):
    prot_seq = altprotseq_dict[alt_acc]
    peps = [p for p in list(peps) if len(p)>8] # min 9aas
    
    pep_locs = []
    for pep in peps:
        pep_start = prot_seq.index(pep)
        pep_end   = pep_start + len(pep)
        pep_locs.append((pep_start, pep_end))
    
    pep_ranges = []
    pep_start, pep_end = 0, 0
    for n, pep_loc in enumerate(sorted(pep_locs)):
        pep_start, pep_end = pep_loc
        if n==0 or pep_start>=prev_pep_end:
            pep_range = [pep_loc]
        elif pep_start<prev_pep_end:
            pep_range.append(pep_loc)
            prev_pep_start, prev_pep_end = pep_loc
            continue
        pep_ranges.append(pep_range)
        prev_pep_start, prev_pep_end = pep_loc
        
    HPP_peps = []
    for pep_range in pep_ranges:
        pep_start, pep_end = max(pep_range, key=lambda x: x[0]-x[1])
        HPP_peps.append(prot_seq[pep_start:pep_end])
    return HPP_peps

def get_nonnested_peps(alt_acc, peps, margin=2):
    prot_seq = altprotseq_dict[alt_acc]
    peps = [p for p in list(peps) if len(p)>8] # min 9aas
    
    pep_len_locs = []
    for pep in peps:
        pep_start = prot_seq.index(pep)
        pep_end   = pep_start + len(pep)
        pep_len_locs.append({
            'pep':pep, 
            'len':len(pep), 
            'start':pep_start, 
            'end':pep_end
        })
    
    pep_len_locs = sorted(pep_len_locs, key=lambda x: x['len'])
    non_nested = []
    for n, pep_len_loc in enumerate(pep_len_locs):
        nested = False
        for longer_pep in pep_len_locs[n+1:]:
            if np.abs(np.diff([pep_len_loc['start'], longer_pep['start']]))<margin or np.abs(np.diff([pep_len_loc['end'], longer_pep['end']]))<margin:
                nested = True
                break
        if not nested:
            non_nested.append(pep_len_loc['pep'])
    return non_nested

def get_coverage(alt_acc, peps):
    prot_seq = altprotseq_dict[alt_acc]
    seq_aa = np.zeros(len(prot_seq))
    for pep in peps:
        pep_start = prot_seq.index(pep)
        pep_end   = pep_start + len(pep)
        seq_aa[pep_start:pep_end] += np.array([1]*len(pep))
    return sum(aa>0 for aa in seq_aa)/len(prot_seq)

def pepset_is_HPP(alt_acc, peps):
    prot_seq = altprotseq_dict[alt_acc]
    peps = [p for p in list(peps) if len(p)>8] # min 9aas
    
    pep_locs = {}
    for pep in peps:
        pep_start = prot_seq.index(pep)
        pep_locs[pep] = (pep_start, pep_start + len(pep))
    
    pep_pairs = [tuple(x) for x in list(set(frozenset(x) for x in itt.permutations(peps, 2))) if len(x)>1]
    pep_pairs_locs = [sorted([pep_locs[pep] for pep in pep_pair]) for pep_pair in pep_pairs]
    
    no_overlap_pairs, overlap_pairs = [], []
    for pep_pair in pep_pairs_locs:
        if pep_pair[0][1]<pep_pair[1][0]:# non-overlapping pairs
            no_overlap_pairs.append(pep_pair)
        else:# overlapping pairs
            overlap_pairs.append(pep_pair)
    valid_overlap_pairs = [pep_pair for pep_pair in overlap_pairs if np.abs(np.diff([pep_pair[0][0], pep_pair[1][1]]))>17]
    
    if len(no_overlap_pairs)>0 or len(valid_overlap_pairs)>0:
        return True
    
    return False

def get_pep_bait_rep_cnts(psms):
    pep_bait_rep_cnts = {}
    pep_reps = [(pep, spec.split('.')[0]) for pep, spec in psms]
    pep_bait_reps = sorted([(pep, rep_bait_dict[rep], rep) for pep, rep in pep_reps])
    for pep, pep_bait_rep_grp in itt.groupby(pep_bait_reps, lambda x: x[0]):
        pep_bait_rep_ls = list(pep_bait_rep_grp)
        pep_bait_rep_cnts[pep] = Counter([bait for pep, bait, rep in pep_bait_rep_ls])
    return pep_bait_rep_cnts

def get_bait_peps(psms):
    bait_peps = {}
    pep_reps = [(pep, spec.split('.')[0]) for pep, spec in psms]
    bait_pep_reps = sorted([(rep_bait_dict[rep], pep, rep) for pep, rep in pep_reps])
    for bait, bait_pep_rep_grp in itt.groupby(bait_pep_reps, key=lambda x: x[0]):
        bait_peps[bait] = set(pep for bait, pep, rep in bait_pep_rep_grp)
    return bait_peps