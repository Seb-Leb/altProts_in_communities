import os
import pickle
import itertools as itt
from Bio import SeqIO

aa_weights = { 
    'X': 110, 'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121, 'E': 147, 'Q': 146,
    'G': 75, 'H': 155, 'I': 131, 'L':131,'K': 146, 'M': 149, 'F': 165, 'P': 115,
    'S': 105, 'T': 119, 'W': 204, 'Y': 181,'V': 117
}

def is_alt(acc):
    if 'IP_' in acc:
        return True
    return False

fasta_fields = ['OS', 'GN', 'TA', 'PA']
def parse_fasta_header(h):
    h = h.split()
    acc = h[0].split('|')[0]
    res = {}
    for f in h[1:]:
        for field in fasta_fields:
            if f[:2] == field:
                res[field] = f[3:]
    if 'GN' not in res:
        res['GN'] = 'unknown'
    if 'PA' in res  and ',' in res['PA']:
        res['PA'] = res['PA'].split(',')[0]
    return res 

def sanitize(s):
    return s.replace('"', '')

def enzyme_digest(protseq, enzyme, miss_cleavage=2, **kwargs): 
    if enzyme == 'trypsin':
        protseq = re.sub("K(?!P)", "K,", protseq)
        protseq = re.sub("R(?!P)", "R,", protseq)
    if enzyme == 'lysc':
        raw_peps_0mc = re.sub("K", "K,", protseq)
    if enzyme == 'gluc':
        raw_peps_0mc = re.sub("E", "E,", protseq)
    if enzyme == 'chymotrypsin':
        if 'chymo_aa' in kwargs.keys():
            chymo_aa = kwargs['chymo_aa']
        else:
            chymo_aa = 'YWFML'
        for a in chymo_aa:
            protseq = protseq.replace(a, a+',')
    raw_peps_0mc = protseq.split(",")
    raw_peps = [] 
    raw_peps.extend(raw_peps_0mc)
    if miss_cleavage == 0:
        return raw_peps
    raw_peps_1mc = [raw_peps_0mc[n] + raw_peps_0mc[n+1] for n in range(len(raw_peps_0mc) - 1)]
    raw_peps.extend(raw_peps_1mc)
    if miss_cleavage == 1:
        return raw_peps
    raw_peps_2mc = [raw_peps_1mc[n] + raw_peps_0mc[n+2] for n in range(len(raw_peps_1mc) - 1)]
    raw_peps.extend(raw_peps_2mc)
    if miss_cleavage == 2:
        return raw_peps
    raise 'miss_cleavage must be 0, 1, or 2.'

def pep_size_filter(raw_peps):
    pred_peps = [raw_pep for raw_pep in raw_peps if len(raw_pep) >= 7 and sum(aa_weights[aa] for aa in raw_pep) < 4600]
    return pred_peps

def check_pep_unicity(pep, refseqs):
    for refseq in refseqs:
        if pep in refseq: return False
    return True

if os.path.exists('prot_gene_dict.pkl'):
    prot_gene_dict = pickle.load(open('prot_gene_dict.pkl', 'rb'))
else:
    prot_gene_dict = {}
    with open('human-openprot-r1_6-refprots+altprots+isoforms-+uniprot2019_03_01.tsv', 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for n,row in enumerate(reader):
            if n==0:
                continue
            elif n==1:
                cols=row
                continue
            line=dict(zip(cols,row))
            prot_gene_dict[line['protein accession numbers'].split('.')[0]] = line['gene symbol']
    pickle.dump(prot_gene_dict, open('prot_gene_dict.pkl', 'wb'))

if os.path.exists('protgene_len_dict.pkl'):
    protgene_len_dict = pickle.load(open('protgene_len_dict.pkl', 'rb'))
else:
    protgene_len_dict = {}
    for record in SeqIO.parse('uman-openprot-r1_6-refprots+altprots+isoforms-+uniprot2019_03_01.fasta', 'fasta'):
        header = parse_fasta_header(record.description)
        gene = header['GN']
        prot_len = len(str(record.seq))
        prot_acc = record.name.split('|')[0]
        if is_alt(prot_acc):
            protgene_len_dict[prot_acc] = prot_len
            continue
        if gene not in protgene_len_dict:
            protgene_len_dict[gene] = prot_len
        elif protgene_len_dict[gene] < prot_len:
            protgene_len_dict[gene] = prot_len
    pickle.dump(protgene_len_dict, open('protgene_len_dict.pkl', 'wb'))

if os.path.exists('altprotseq_dict.pkl'):
    altprotseq_dict = pickle.load(open('altprotseq_dict.pkl', 'rb'))
else:
    altprotseq_dict = {}
    for record in SeqIO.parse('human-openprot-r1_6-refprots+altprots+isoforms-+uniprot2019_03_01.fasta', 'fasta'):
        prot_acc = record.name.split('|')[0]
        if is_alt(prot_acc):
            altprotseq_dict[prot_acc] = str(record.seq)
    pickle.dump(altprotseq_dict, open('altprotseq_dict.pkl', 'wb'))

exp_id_name = {}
with open('BioPlex_exp_id_name.csv', 'r') as f:
    for l in f:
        exp_id, exp_name = l.strip().split(',')
        exp_id_name[exp_name] = exp_id

baits_exps = []
with open('archiveSummary_BP2.tsv','r') as f:
    for n,l in enumerate(f):
        if n==0: continue
        ls = l.split('\t')
        file_details = ls[0].split('_')
        exp_name, plate, well, rep  = file_details[:4]
        bait   = ls[2].split(' ')[0][6:].strip()
        baits_exps.append((exp_name, plate, well, rep, bait, exp_id_name[exp_name]))
        
bp1_baits = set()
with open('archiveSummary_BP1.tsv','r') as f:
    for n,l in enumerate(f):
        if n==0: continue
        ls = l.split('\t')
        file_details = ls[0].split('_')
        exp_name, plate, well, rep = file_details[:4]
        bait = ls[2].split(' ')[0][6:]
        bp1_baits.add(bait)
        if exp_name not in set(x[0] for x in baits_exps):
            baits_exps.append((exp_name, plate, well, rep, bait, exp_id_name[exp_name]))

bait_rep_dict = {bait:[x[0] for x in grp] for bait,grp in itt.groupby(baits_exps, lambda x: x[4])}
rep_bait_dict = {exp_id:bait for exp_id, _, _, _, bait, _ in baits_exps}
bait_num_dict = {b:n for b,n in zip(bait_rep_dict.keys(), range(len(bait_rep_dict.keys())))}
batch_n_wells = {batch:len(set(x[2] for x in grp)) for batch,grp in itt.groupby(sorted(baits_exps, key=lambda x: x[1]), key=lambda x: x[1])}
bait_batch_dict = {bait:batch for _,batch,_,_,bait,_ in baits_exps}
all_baits = list(bait_rep_dict.keys())

if all(os.path.exists(x) for x in ['BP_noFilters.pkl', 'BP_preys_per_bait.pkl', 'BP2_unflt_psms.pkl']):
    BP_noFilters      = pickle.load(open('BP_noFilters.pkl', 'rb'))
    BP_preys_per_bait = pickle.load(open('BP_preys_per_bait.pkl', 'rb'))
    BP2_unflt_psms    = pickle.load(open('BP2_unflt_psms.pkl', 'rb'))
else:
    with open('BaitPreyPairs_noFilters_BP2a.tsv', 'r') as f:
        for n,l in enumerate(f):
            ls = l.strip().split('\t')
            if n==0:
                keys = ls
                continue
            line = dict(zip(keys, ls))
            
            for i in line:
                if not line[i]:continue
                line[i] = dtypes[i](line[i])
                
            BP_noFilters.append(line)
            if line['bait_symbol'] not in BP_preys_per_bait:
                BP_preys_per_bait[line['bait_symbol']] = [line['symbol'],]
            else:
                BP_preys_per_bait[line['bait_symbol']].append(line['symbol'])
    
    BP2_unflt_psms = dict()
    grps = itt.groupby(sorted(BP_noFilters, key=lambda x: x['bait_symbol']), key=lambda x: x['bait_symbol'])
    for bait, bait_grp in grps:
        BP2_unflt_psms[bait] = {bp['symbol']:bp['total_psms'] for bp in bait_grp}
    
    pickle.dump(BP_noFilters, open('BP_noFilters.pkl', 'wb'))
    pickle.dump(BP_preys_per_bait, open('BP_preys_per_bait.pkl', 'wb'))
    pickle.dump(BP2_unflt_psms, open('BP2_unflt_psms.pkl', 'wb'))