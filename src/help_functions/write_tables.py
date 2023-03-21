def write_tables(dataframe, filename_prepend):
    dataframe = dataframe.copy()

    maf_df = dfs.get_df()

    tmp = [(maf_df[maf_df["Entrez_Gene_Id"] == int(i)])["Hugo_Symbol"] if len(i.split(",")) == 1 else "region" for i in dataframe["geneName"] ]
    tmp = [i[0] if len(i) != 0 else None for i in tmp]
    dataframe["hugo_symbol"] = tmp


    tmp = [(dfs.get_information_about_mutations(row[1]["geneName"], row[1]["sample"], row[1])) for row in dataframe.iterrows()]
    tmp = [i[['Hugo_Symbol','Start_Position', 'End_Position', 'Variant_Type', 'Variant_Classification']] if isinstance(i,pd.DataFrame) else pd.DataFrame(data={"Hugo_Symbol": None , "Start_Position": None, "End_Position" : None, "Variant_Type":None,"Variant_Classification":None}, index=[0]) for i in tmp]

    tmp = pd.concat(tmp)

    for i in tmp.columns:
        dataframe[i] = tmp[i].to_numpy()

    print(dataframe.columns)
    for i in ["txName", "geneID", "featureID"]:
        if i in dataframe.columns:
            dataframe.drop(columns=[i])

    (dataframe[dataframe["type"] == "E"].head(50)).to_csv(f"{filename_prepend}_50_highest_score_exon.csv", sep=";")
    (dataframe[dataframe["type"] == "E"].tail(50)).to_csv(f"{filename_prepend}_50_lowest_score_exon.csv", sep=";")

    (dataframe[dataframe["type"] == "J"].head(50)).to_csv(f"{filename_prepend}_50_highest_score_junction.csv", sep=";")
    (dataframe[dataframe["type"] == "J"].tail(50)).to_csv(f"{filename_prepend}_50_lowest_junction.csv", sep=";")



    (dataframe[(dataframe["type"] == "E") & (dataframe["is_related_mutated"])].head(50)).to_csv(f"{filename_prepend}_50_highest_score_exon_m_mut.csv", sep=";")
    (dataframe[(dataframe["type"] == "E") & (dataframe["is_related_mutated"])].tail(50)).to_csv(f"{filename_prepend}_50_lowest_score_exon_m_mut.csv", sep=";")

    (dataframe[(dataframe["type"] == "J") & (dataframe["is_related_mutated"])].head(50)).to_csv(f"{filename_prepend}_50_highest_score_junction_m_mut.csv", sep=";")
    (dataframe[(dataframe["type"] == "J") & (dataframe["is_related_mutated"])].tail(50)).to_csv(f"{filename_prepend}_50_lowest_junction_m_mut.csv", sep=";")


    return dataframe