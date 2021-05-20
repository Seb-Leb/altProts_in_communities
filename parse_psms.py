import os
import re
import pickle
import csv
import sys
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
from hierarchical_report_parser import *

aa_weights = { 
    'X': 110, 'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121, 'E': 147, 'Q': 146,
    'G': 75, 'H': 155, 'I': 131, 'L':131,'K': 146, 'M': 149, 'F': 165, 'P': 115,
    'S': 105, 'T': 119, 'W': 204, 'Y': 181,'V': 117
}

