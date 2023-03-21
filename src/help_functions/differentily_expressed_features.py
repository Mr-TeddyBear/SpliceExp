""""
This script should be converted in to a proper function and is not complete

"""
geneID = "7273"





differential_df = pd.DataFrame(columns=["geneID", "start", "end", "sample", "fold_change", "p_value"])

p_values = {}
fold_change = {}

use_df = median_ratio_norm

for geneIDs in use_df["geneName"].dropna().unique():
    for geneID in geneIDs.split(","):
        exons = use_df[(use_df["geneName"] == geneID)]
        junctions = exons[exons["type"] == "J"]
        exons =  exons[exons["type"] == "E"]
        muta_samp = dfs.get_mutated_samples(geneID)
        

        norm_samples = [i for i in exons.columns[12:] if i not in muta_samp]

        for row in exons.iterrows():
            for i in muta_samp:
                if dfs.is_mutation_close_to_feature(geneID, i, row[1]):
                    junction_starting_in_exon = junctions[(row[1]["start"] <= junctions["start"]) & (junctions["start"] <= row[1]["end"])]
                    for single_junction in junction_starting_in_exon.iterrows():
                        
                        mut_exons = pd.concat([exons.loc[[row[0]]],exons[(exons["start"] == single_junction[1]["end"])]])

                        #Perform ttest on subsequent exons in gene.
                        t_test = ttest_ind(mut_exons[[i]].mean(axis=1).to_numpy(), mut_exons[norm_samples].mean(axis=1).to_numpy(), axis=None, equal_var=True)
                        print(t_test, mut_exons[norm_samples])
                        print(mut_exons[[i]].mean(axis=1).to_numpy(), mut_exons[norm_samples].mean(axis=1).to_numpy())
                        if geneID in p_values:
                            p_values[geneID].append(t_test)
                        else:
                            p_values[geneID] = [t_test]

                        # Calculate the fold change

                        #calc_fold_change = lambda group1, group2: np.average( np.log2(group1) ) - np.average( np.log2(group2) )
                        
                        fchg = calc_fold_change(mut_exons[[i]], mut_exons[norm_samples])


                        if geneID in fold_change:
                            fold_change[geneID].append(fchg)
                        else:
                            fold_change[geneID] = [fchg]

                        print(mut_exons[[i]], mut_exons[norm_samples]) 
                        print(fchg, "\n")

                        differential_df = pd.concat([differential_df, pd.DataFrame({"geneID": row[1]["geneName"], "start": row[1]["start"], "end": row[1]["end"], "sample": i, "fold_change": fchg, "p_value": t_test.pvalue}, index=[0])])
                break