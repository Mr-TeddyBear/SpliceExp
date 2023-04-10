import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os


def vulcano_plot(df, x, y, dfs, y_function=lambda x: x, show=False, save=False, filename=None):

    plt.figure()

    plt.scatter(x=df[x], y=df[y].apply(lambda x: -np.log10(x)),
                s=1, label="Not significant")
    muts = dfs.get_df()
    # highlight down- or up- regulated genes
    down = df[(df[x] <= -2) & (df[y] <= 0.01)]
    up = df[(df[x] >= 2) & (df[y] <= 0.01)]
    print("\n\n\n", up.head(), "\n\n\n")
    print("\n\n\n", down.head(), "\n\n\n")
    plt.scatter(x=down[x], y=down[y].apply(
        lambda x: -np.log10(x)), s=3, label="Down-regulated", color="blue")
    plt.scatter(x=up[x], y=up[y].apply(lambda x: -np.log10(x)),
                s=3, label="Up-regulated", color="red")

    texts = []
    for i, r in up.iterrows():
        texts.append(plt.text(x=r[x], y=-np.log10(r[y]), s=muts[muts["Entrez_Gene_Id"]
                     == int(r["geneName"])]["Hugo_Symbol"].unique()[0]))
    for i, r in down.iterrows():
        texts.append(plt.text(x=r[x], y=-np.log10(r[y]), s=muts[muts["Entrez_Gene_Id"]
                     == int(r["geneName"])]["Hugo_Symbol"].unique()[0]))

    # adjust_text(texts,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    plt.xlabel("log2FC")
    plt.ylabel("-logFDR")
    plt.axvline(-2, color="grey", linestyle="--")
    plt.axvline(2, color="grey", linestyle="--")
    plt.axhline(2, color="grey", linestyle="--")
    plt.legend()
    plt.autoscale()

    if save:
        if filename:
            plt.savefig(filename)
        else:
            plt.savefig(f"vulcano_plot_{filename}.pdf")

    if show:
        plt.show()

def ma_plot(df, x, y, dfs, show=False, save=False, filename=None):

    plt.figure()

    plt.scatter(x=np.log2(df[x]), y=df[y],
                s=1, label="Not significant")
    #muts = dfs.get_df()
    # highlight down- or up- regulated genes
    #down = df[(df[x] <= -2) & (df[y] <= 0.01)]
    #up = df[(df[x] >= 2) & (df[y] <= 0.01)]
    #print("\n\n\n", up.head(), "\n\n\n")
    #print("\n\n\n", down.head(), "\n\n\n")
    #plt.scatter(x=down[x], y=down[y].apply(
    #    lambda x: -np.log10(x)), s=3, label="Down-regulated", color="blue")
    #plt.scatter(x=up[x], y=up[y].apply(lambda x: -np.log10(x)),
    #            s=3, label="Up-regulated", color="red")

    #texts = []
    #for i, r in up.iterrows():
    #    texts.append(plt.text(x=r[x], y=-np.log10(r[y]), s=muts[muts["Entrez_Gene_Id"]
    #                 == int(r["geneName"])]["Hugo_Symbol"].unique()[0]))
    #for i, r in down.iterrows():
    #    texts.append(plt.text(x=r[x], y=-np.log10(r[y]), s=muts[muts["Entrez_Gene_Id"]
    #                 == int(r["geneName"])]["Hugo_Symbol"].unique()[0]))

    # adjust_text(texts,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    plt.xlabel("Mean expression")
    plt.ylabel("-logFDR")
    #plt.axvline(-2, color="grey", linestyle="--")
    #plt.axvline(2, color="grey", linestyle="--")
    #plt.axhline(2, color="grey", linestyle="--")
    #plt.legend()
    plt.autoscale()

    if save:
        if filename:
            plt.savefig(filename)
        else:
            plt.savefig(f"vulcano_plot_{filename}.pdf")

    if show:
        plt.show()
        




def plot_gene_heatmap(geneName, dataframe, mutations, center=0,show = False, save = False, filename=None):
    """
    Plot the heatmap of expression levels for given gene.

    Args:
        geneName (str): Entrez gene Id
        dataframe (pd.DataFrame): A dataframe containing information about given gene
        mutations (pd.DataFrame): A dataframe containging information about 
                                  mutations for samples.
        center (int, str): Middle value for heatmap color scale.
        show (bool): If figures should be shown
        save (bool): If figures should be saved
        filename (str): Figure name, not used if save if False
    """
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
            plt.savefig(f"{looking_for}_{gene.symbol}.pdf")

    if show:
        plt.show()
        


def plot_binned_mutations(df, target_column, title, n_bins = 10, write=False, path=None, filetype="pdf"):
    exons = df[df["type"] == "E"]
    junction = df[df["type"] == "J"]

    fig, ax = plt.subplots(2,1)


    skip = len(exons[target_column])//n_bins
    bins = [sum(exons[target_column][i:i+skip]) for i in range(0, len(exons[target_column]), skip)]
    ax[0].bar(x = range(len(bins)), height = bins)
    ax[0].set_title("Exon")
    ax[0].set_xlim([-1,len(bins)])


    skip = len(junction[target_column])//n_bins
    bins = [sum(junction[target_column][i:i+skip]) for i in range(0, len(junction[target_column]), skip)]
    ax[1].bar(range(len(bins)), bins)
    ax[1].set_title("Junction")
    ax[1].set_xlim([-1,len(bins)])

    print(len(bins), skip)

    fig.suptitle(title)

    if write:
        fig.savefig(os.path.join(path, f"distribution_mutated_features.{filetype}"))

    return fig, ax