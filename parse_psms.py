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
from utils import *

# Make a list of ref seq
refseqs = []
for line in SeqIO.parse("human-openprot-r1_6-refprots-+uniprot2019_03_01.fasta", "fasta"):
    refseqs.append(str(line.seq))

