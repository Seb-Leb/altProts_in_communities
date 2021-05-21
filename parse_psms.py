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
from scipy.stats import binom_test
import statsmodels.api as sm
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
                    #if fasta_header['PA'] not in prot_gene_dict: continue
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
        
    #return psms
    #pickle.dump(psms, open('bioplex_psms/{}_psms.pkl'.format(bait), 'wb'))
    return psms