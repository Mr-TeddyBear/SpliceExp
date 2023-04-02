import pandas as pd


def write_tables(dataframe, dfs, filename_prepend, path, ends=None, n=None):
    """
    This function writes dataframes to csv files. 

    Args:

    dataframe (pandas.DataFrame) : A stacked features dataframe
    dfs (mafWrapper) ; A mafwarpper instanace
    filename_prepend (str) : string that is put before .csv in the filename.
    path (str) : 
    ends (str) : None, head or tail. Choose of only part of the dataframe should be written to file.
    n (int) : Only used of ends is set. Choose n tail or n head entries.
    
    """
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



    match ends:
        case "tail":
            (dataframe[dataframe["type"] == "E"].tail(n)).to_csv(f"{filename_prepend}_exon_tail_{n}.csv", sep=";")
            (dataframe[dataframe["type"] == "J"].tail(n)).to_csv(f"{filename_prepend}_junction_tail_{n}.csv", sep=";")

        case "head":
            (dataframe[dataframe["type"] == "E"].head(50)).to_csv(f"{filename_prepend}_exon_head_{n}.csv", sep=";")
            (dataframe[dataframe["type"] == "J"].head(50)).to_csv(f"{filename_prepend}_junction_head_{n}.csv", sep=";")

        case None:
            (dataframe[dataframe["type"] == "E"]).to_csv(f"{filename_prepend}_exon.csv", sep=";")
            (dataframe[dataframe["type"] == "J"]).to_csv(f"{filename_prepend}_junction.csv", sep=";")



    return dataframe