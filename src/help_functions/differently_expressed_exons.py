from help_functions.fold_change import calc_fold_change

import pandas as pd
from scipy.stats import ttest_ind
import numpy as np


def differently_expressed_linked_exons(use_df, dfs):
    """
    Not completed function
    """
    differential_df = pd.DataFrame(
        columns=["geneName", "start", "end", "sample", "fold_change", "p_value"])

    p_values = {}
    fold_change = {}

    for geneIDs in use_df["geneName"].dropna().unique():
        for geneID in geneIDs.split(","):
            exons = use_df[(use_df["geneName"] == geneID)]
            junctions = exons[exons["type"] == "J"]
            exons = exons[exons["type"] == "E"]
            muta_samp = dfs.get_mutated_samples(geneID)

            norm_samples = [i for i in exons.columns[12:]
                            if i not in muta_samp]

            for row in exons.iterrows():
                for i in muta_samp:
                    if dfs.is_mutation_close_to_feature(geneID, i, row[1]):
                        junction_starting_in_exon = junctions[(
                            row[1]["start"] <= junctions["start"]) & (junctions["start"] <= row[1]["end"])]
                        for single_junction in junction_starting_in_exon.iterrows():

                            mut_exons = pd.concat(
                                [exons.loc[[row[0]]], exons[(exons["start"] == single_junction[1]["end"])]])

                            t_test = ttest_ind(mut_exons[[i]].mean(axis=1).to_numpy(
                            ), mut_exons[norm_samples].mean(axis=1).to_numpy(), axis=None, equal_var=True)
                            if geneID in p_values:
                                p_values[geneID].append(t_test)
                            else:
                                p_values[geneID] = [t_test]

                            fchg = calc_fold_change(
                                mut_exons[[i]], mut_exons[norm_samples])

                            if geneID in fold_change:
                                fold_change[geneID].append(fchg)
                            else:
                                fold_change[geneID] = [fchg]

                            differential_df = pd.concat([differential_df, pd.DataFrame(
                                {"geneName": row[1]["geneName"], "start": row[1]["start"], "end": row[1]["end"], "sample": i, "fold_change": fchg, "p_value": t_test.pvalue}, index=[0])])

    return differential_df


def differently_expressed_exons(use_df, dfs, samples):
    """
    Not completed function
    """
    differential_df = pd.DataFrame(
        columns=["geneName", "start", "end", "sample", "fold_change", "mean_expression"])

    p_values = {}
    fold_change = {}

    for geneIDs in use_df["geneName"].dropna().unique():
        for geneID in geneIDs.split(","):
            exons = use_df[(use_df["geneName"] == geneID)]
            junctions = exons[exons["type"] == "J"]
            exons = exons[exons["type"] == "E"]
            muta_samp = dfs.get_mutated_samples(geneID)

            norm_samples = [i for i in exons[samples]
                            if i not in muta_samp]

            for row in exons.iterrows():
                for i in muta_samp:
                    if dfs.is_mutation_close_to_feature(geneID, i, row[1]):
                        junction_starting_in_exon = junctions[(
                            row[1]["start"] <= junctions["start"]) & (junctions["start"] <= row[1]["end"])]
                        for single_junction in junction_starting_in_exon.iterrows():

                            mut_exons = pd.concat(
                                [exons.loc[[row[0]]], exons[(exons["start"] == single_junction[1]["end"])]])
                            

                            mean_expression = row[1][samples].average()

                            fchg = calc_fold_change(
                                mut_exons[[i]], mut_exons[norm_samples])

                            if geneID in fold_change:
                                fold_change[geneID].append(fchg)
                            else:
                                fold_change[geneID] = [fchg]

                            differential_df = pd.concat([differential_df, pd.DataFrame(
                                {"geneName": row[1]["geneName"], "start": row[1]["start"], "end": row[1]["end"], "sample": i, "fold_change": fchg, "mean_expression": mean_expression}, index=[0])])

    return differential_df