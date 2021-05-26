import re
import sys
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy.stats import binom_test
import statsmodels.api as sm
from multiprocessing import Pool

feats = ['Z', 'WD', 'batch_Z', 'Entropy', 'Ratio', 'total_PSMs', 'ratio_total_PSMs', 'pep_ratio']

def batch_Z(batch, prey, psm, dist_dict):
    u, se = dist_dict[batch, prey]
    return (psm - u)/se

def bin_feature(feat_vals):
    bins = np.linspace(min(feat_vals)-1e-5, max(feat_vals), 1000)
    return bins