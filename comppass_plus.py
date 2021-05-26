import os
import csv
import sys
import pandas as pd
import pickle
import itertools as itt
import numpy as np
from sklearn.naive_bayes import GaussianNB
from orderedset import OrderedSet
from multiprocessing import Pool

# column data types for feature table
columns_dtypes = {
    'Batch': str,
    'Prey': str,
    'AvePSM': float,
    'Bait': str,
    'Entropy': float,
    'Experiment_ID': str,
    'Prey_Gene': str,
    'WD': float,
    'Z': float,
    'cnt_peps': int,
    'cnt_upeps': int,
    'excluded_by_binom_test': lambda x: 0 if x == 'False' else 1,
    'is_alt': lambda x: 0 if x == 'False' else 1,
    'Ratio': float,
    'total_PSMs': lambda x: int(float(x)),
    'ratio_total_PSMs': float,
    'pep_bin': int,
    'pep_ratio': float,
    'batch_Z': float,
    'Z_binned': int,
    'WD_binned': int,
    'batch_Z_binned': int,
    'Entropy_binned': int,
    'Ratio_binned': int,
    'total_PSMs_binned': int,
    'ratio_total_PSMs_binned': int,
    'pep_ratio_binned': int,
    'HCIP_label': int
}

def get_batch_data(batch, training=False):
    batch_data = pickle.load(open('bioplex_data_batches/{}.pkl'.format(batch), 'rb'))
    baits_preys, v, t = zip(*batch_data)
    if training:
        mask_neg = np.array([frozenset(bp) in BP_neg for bp in baits_preys])
        mask_pos = t
        mask = mask_pos + mask_neg
        baits_preys, v, t = [list(itt.compress(x, mask)) for x in [baits_preys, v, t]]
    return baits_preys, v, t

def get_training_data(batch):
    X,y = [], []
    for b in all_batches:
        if batch == b: continue # exclude a batch's data from its training set
        bp, v, t = get_batch_data(b, training=False)
        mask = [not is_alt(x[1]) for x in bp]
        v = list(itt.compress(v, mask))
        t = list(itt.compress(t, mask))
        X.extend(v)
        y.extend(t)
    return X, y

def train(batch):
    X, y = get_training_data(batch)
    classifier   = GaussianNB()
    model        = classifier.fit(X, y)
    return {batch : model}

def compute_models(batches, serial=True):
    batch_models_dict = {}
    if not serial:
        with Pool(10) as p:
            models_by_batch = p.map(train, batches)
    else:
        models_by_batch = []
        for b in batches:
            models_by_batch.append(train(b))
            
    for m in models_by_batch:
        batch_models_dict.update(m)
    return batch_models_dict

def predict(batch):
    bp, v, t = get_batch_data(batch)
    clf = batch_models_dict[batch]
    res = clf.predict_proba(v)
    return [(batch, *bp, t, res) for bp, t, res in zip(bp, t, res)]

def predict_all(all_batches):
    with Pool(12) as p:
        results = p.map(predict, all_batches)
    return [x for y in results for x in y]