import numpy as np


def median_of_ratios_normalization(exon_junction_df, samples):
    """
    Function that performs median of ratios normalization.

    Args:
        exon_junction_df (pd.DataFrame): Dataframe containgen non-normalized data
        samples (list): list of samples to be normalized in dataframe

    Returns:
        pd.DataFrame: a normalized dataframe. Contains all the same fields as in input
                      dataframe, but the counts are normalized.

    """
    ratio_df = exon_junction_df.copy()
    # Calcualte the psedua referance sample value
    ratio_df["pseudo_reference_sample"] = exon_junction_df[samples].apply(lambda x: np.power(x[x != 0], 1/len(x[x != 0])).prod(), axis=1)
    #display(ratio_df)

    #Calcualte normalisation facotr for rows
    ratio_df[samples] = ratio_df.apply(lambda x: x[samples]/x["pseudo_reference_sample"], axis=1)

    norm_df = exon_junction_df.copy()

    #display(ratio_df.head())

    for i in (exon_junction_df.geneName.dropna()).unique():
        # Calcualte the sample normalisation factor
        tmp_exon_norm_factor = np.median(ratio_df[(ratio_df["geneName"] == i) & (ratio_df["type"] == "E")][samples], axis=0)
        tmp_junction_norm_factor = np.median(ratio_df[(ratio_df["geneName"] == i) & (ratio_df["type"] == "J")][samples], axis=0)

        # Do the actual normalisation of samples
        norm_df.loc[(norm_df["geneName"] == i) & (norm_df["type"] == "E"), samples] = norm_df.loc[(norm_df["geneName"] == i) & (norm_df["type"] == "E")][samples]/tmp_exon_norm_factor
        norm_df.loc[(norm_df["geneName"] == i) & (norm_df["type"] == "J"), samples] = norm_df.loc[(norm_df["geneName"] == i) & (norm_df["type"] == "J")][samples]/tmp_junction_norm_factor


    return norm_df