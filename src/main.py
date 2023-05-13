# Global imports
import argparse
from plotting.plotters import vulcano_plot, plot_binned_mutations, ma_plot, plot_gene_heatmap
from help_functions.differently_expressed_exons import differently_expressed_linked_exons, differently_expressed_exons
from help_functions.find_expression import find_expressions
from help_functions.median_of_ratios_normalization import median_of_ratios_normalization
from mafWrapper import mafWrapper
import pandas as pd
import numpy as np


import glob
import os

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.path as mpath
import matplotlib.patches as mpatches

from scipy.stats import ttest_ind, ttest_rel, f_oneway, t

import seaborn as sns

sns.set()

# from plotting import plot_mutation
# from plot_data import plot_exons

# Local imports


parser = argparse.ArgumentParser(
    prog='SpliceExp',
    description='Expression analyses '
)

parser.add_argument('--MAF', type=str, nargs=1,
                    help="path folder containing MAF file for all samples and a map.txt file containing sample name and file name for sample. Sample name must correlate with sample name in count file.")
parser.add_argument('--features', type=str, nargs=1,
                    help="path to csv file containing all feature counts as produced på SGSeq.")
parser.add_argument('--sep', default="\t", type=str, nargs=1,
                    help="If a different seperation is used in features csv. Deafult \\t")
parser.add_argument("--out", default="." + os.path.sep, type=str, nargs=1,
                    help="Path to folder where outputs will be placed. Will create a new folder output at path given.")
parser.add_argument("--save", default=True, action='store_true',
                    help="Arguemnt to write tables and figures to file. Deafult true. Will overwrite old fiures in outfolder.")
parser.add_argument("--plot_filetype", default="pdf", type=str,
                    help="Sets the filetype of exported plots. All formats accepted by matplotlib as accepted.")

args = parser.parse_args()


if not os.path.isdir(args.MAF[0]):
    raise FileNotFoundError("Path to MAF folder not found")

if not os.path.isfile(args.features[0]):
    raise FileNotFoundError("Path to features not found")

if not os.path.isdir(args.out[0]):
    raise FileNotFoundError("Outpath not found for argument given in --out")

if not os.path.isdir(os.path.join(args.out, "output")):
    os.makedirs(os.path.join(args.out, "output"))

outpath = os.path.abspath(os.path.join(args.out[0], "output"))


if args.save:
    save = True
else:
    save = False

show = False

# Read in features from csv

df = pd.read_csv(args.features[0], sep=args.sep[0])

# Drop ending of tumor type, not needed as we dont have multiple samples from same patients
df.columns = df.columns.str.replace('-01A', '')
s = df.columns[12:]

# Drop duplicate rows and remove features with unknown gene location.
df.drop_duplicates(inplace=True)
df = df[df["geneName"].str.contains("nan") == False]
df = df[df[df.columns[12:]].sum(axis=1) > 5]

# Save name of all sample, maybe this should be read from a file


# Read maf_path and initatie the wrapper for searching mutations in MAFs
dfs = mafWrapper(os.path.abspath(args.MAF[0]))



median_ratio_norm = median_of_ratios_normalization(df, s)
median_ratio_norm = median_ratio_norm[median_ratio_norm["geneName"].notna()]
median_ratio_norm_centerd = median_ratio_norm.copy()
median_ratio_norm_centerd[s] = median_ratio_norm_centerd[s].apply(lambda x: x-x.mean(), axis=1)
median_sorted = find_expressions(
    df=median_ratio_norm_centerd, dfs=dfs, raw_df=df, col=s, show=False)

plot_binned_mutations(median_sorted, "is_related_mutated", "Distribution of samples with mutations",
                      path=outpath, filetype=args.plot_filetype, write=save)


differently_exons = differently_expressed_exons(median_ratio_norm, dfs, s)

linked_expression = differently_expressed_linked_exons(median_ratio_norm, dfs)


ma_plot(differently_exons, "mean_expression", "fold_change", dfs, filename=os.path.join(outpath, f"ma_single_exons.{args.plot_filetype}"), show=show, save=True)


vulcano_plot(linked_expression, x='fold_change', y='p_value', dfs=dfs, save=save,
             filename=os.path.join(outpath, f"linked_expression_vulcano.{args.plot_filetype}", ))


if save:
    median_ratio_norm.to_csv(os.path.join(outpath, "normalized_counts.csv"), sep=";")
    median_sorted.to_csv(os.path.join(outpath, "stacked_normalized_counts.csv"), sep=";")
    median_ratio_norm_centerd.head(50).to_csv(os.path.join(outpath, "normalized_counts_centerd_top_50.csv"), sep=";")
    median_ratio_norm_centerd.tail(50).to_csv(os.path.join(outpath, "normalized_counts_centerd_bottom_50.csv"), sep=";")

