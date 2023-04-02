#Global imports
import pandas as pd
import numpy as np

from pybiomart import Dataset, Server

import glob
import os

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.path as mpath
import matplotlib.patches as mpatches

from scipy.stats import ttest_ind, ttest_rel, f_oneway, t

import seaborn as sns

sns.set()

from plot_mutation import plot_mutation
from plot_data import plot_exons

#Local imports
from mafWrapper import mafWrapper
from help_functions.median_of_ratios_normalization import median_of_ratios_normalization
from help_functions.find_expression import find_expressions
from plotting.plot_binned_mutations import plot_binned_mutations



import argparse

parser = argparse.ArgumentParser(
    prog='SpliceMEAN',
    description='Expression analyses '
)

parser.add_argument('--MAF', type=str, nargs=1, help="path folder containing MAF file for all samples and a map.txt file containing sample name and file name for sample. Sample name must correlate with sample name in count file.")
parser.add_argument('--features', type=str, nargs=1, help="path to csv file containing all feature counts as produced på SGSeq.")
parser.add_argument('--sep', default="\t", type=str, nargs=1, help="If a different seperation is used in features csv. Deafult \\t")
parser.add_argument("--out", default="./", type=str, nargs=1, help="Path to folder where outputs will be placed. Will create a new folder output at path given.")
parser.add_argument("--save", default=True, type=bool, help="Arguemnt to write tables and figures to file. Deafult true. Will overwrite old fiures in outfolder.")
parser.add_argument("--plot_filetype", default="pdf", type=str, help="Sets the filetype of exported plots. All formats accepted by matplotlib as accepted.")

args = parser.parse_args()

if not os.path.isdir(args.MAF):
    raise FileNotFoundError("Path to MAF folder not found")

if not os.path.isfile(args.features):
    raise FileNotFoundError("Path to features not found")

if not os.path.isdir(args.out):
    raise FileNotFoundError("Outpath not found for argument given in --out")

if not os.path.isdir(os.path.join(args.out, "output")):
    os.makedirs(os.path.join(args.out, "output"))


#Read in features from csv

df = pd.read_csv(args.features, sep=args.sep)

#Drop ending of tumor type, not needed as we dont have multiple samples from same patients
df.columns = df.columns.str.replace('-01A','')

#Drop duplicate rows and remove features with unknown gene location.
df.drop_duplicates(inplace=True)
df = df[df["geneName"].str.contains("nan") == False]


#Save name of all sample, maybe this should be read from a file
s = df.columns[12:]


#Read maf_path and initatie the wrapper for searching mutations in MAFs
dfs = mafWrapper(args.MAF)


median_ratio_norm = median_of_ratios_normalization(df, s)
median_ratio_norm = median_ratio_norm[median_ratio_norm["geneName"].notna()]
median_sorted = find_expressions(median_ratio_norm, df, s, show=False)

plot_binned_mutations(median_sorted, "is_related_mutated", "Distribution of samples with mutations, where mutations is close to the exon in question. Left is features with highest expression.", path=args.out, filetype=args.plot_type)