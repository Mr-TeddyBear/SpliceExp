import matplotlib.pyplot as plt
import os

def plot_gene_heatmap(geneName, dataframe, mutations, center=0,show = False, save = False, filename=None):
    """
    Plot the heatmap of expression levels for given gene.

    Args:
        geneName (str): Hugo Symbol for given to plot
        dataframe (pd.DataFrame): A dataframe containing information about given gene
        mutations (mafWrapper): An instace of the mafWrapper class.
        center (int, str): Middle value for heatmap color scale.
        show (bool): If figures should be shown
        save (bool): If figures should be saved
        filename (str): Figure name, not used if save if False
    """
    geneName = mutations.get_entrez_id(geneName)

    data = dataframe.loc[dataframe["geneName"] == str(geneName)]

    exons = data[data["type"] == "E"]
    junctions = data[data["type"] == "J"]

    fig, axs = plt.subplots(nrows=4,ncols=1, figsize=(10,15), gridspec_kw={"height_ratios": [1,1,3,3]})


    #Plot postion of mutations in gene
    sample_muts = mutations.get_geneID(geneName)
    for key in sample_muts:
        value = sample_muts[key]
        if not value.empty:
            #print(value["Start_Position"])
            for ax in axs[:2]:
                for i in range(value["Start_Position"].values[0],value["End_Position"].values[0]+1):
                    ax.axvline(i, label=key, alpha=0.5)

    mut_exon = np.zeros_like(exons["type"], dtype=bool)

    #Plot exons
    count = 0
    for index, row in exons.iterrows():
        axs[0].plot([row["start"],row["end"]], [count,count], c="b")
        for key in sample_muts:
            value = sample_muts[key]
            if not value.empty:
                for s,e in zip(value["Start_Position"], value["End_Position"]):
                    if row["start"] - 10 < s and e < row["end"] + 10:
                        mut_exon[count] = True
        count += 1


    mut_junction = np.zeros_like(junctions["type"], dtype=bool)
    #Plot junctions
    count = 0
    for index, row in junctions.iterrows():
        axs[1].plot([row["start"], row["end"]], [count,count], c="r")
        for key in sample_muts:
            value = sample_muts[key]
            if not value.empty:
                for s,e in zip(value["Start_Position"], value["End_Position"]):
                    if row["start"] - 10 < s and e < row["end"] + 10:
                        mut_junction[count] = True
        count += 1

       
    
    mut_pos = {"E": mut_exon, "J":mut_junction}

    sns.set(rc = {'figure.figsize':(20,5)})

    for i,j in enumerate(["E", "J"]):
        tmp_data = data[data["type"] == j]
        if center == "mean" : 
            g = sns.heatmap((tmp_data.iloc[:,12:]), ax=axs[i+2], yticklabels=tmp_data.iloc[:,5], cmap="Spectral_r", linecolor="k", linewidths=0.009)
        else:
            g = sns.heatmap((tmp_data.iloc[:,12:]), ax=axs[i+2], yticklabels=tmp_data.iloc[:,5], center=0, cmap="Spectral_r", linecolor="k", linewidths=0.009)

        #print(g.get_yticklabels())
        g.set_yticklabels([f"{j.get_text()}{i}" for i,j in enumerate(g.get_yticklabels())], rotation=0)
            

        for tick_label in g.get_xticklabels():

            tick_text = tick_label.get_text()
            if not sample_muts[tick_text].empty:
                tick_label.set_color("r")

        for p, tick_label in enumerate(g.get_yticklabels()):
            if mut_pos[j][p]:
                tick_label.set_color("r")


    fig.suptitle(f"{mutations.get_df()[mutations.get_df()['Entrez_Gene_Id'] == geneName]['Hugo_Symbol'][0]}")
    
    if save:
        if filename:
            plt.savefig(filename)
        else:
            plt.savefig(os.path.join(filename, f"{geneName}_heatmap.pdf"))

    if show:
        plt.show()