import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os
from adjustText import adjust_text

def vulcano_plot(df, x, y, dfs, y_function=lambda x: x, show=False, save=False, filename=None):
    """
    Creates a vulcano plot from a pandas dataframe.
    Args:
        df (pandas.DataFrame)   : A pandas dataframe containg all values that should be plottet.
        x (str)                 : Name of the column in the dataframe used for x values.
        y (str)                 : Name of column in the dataframe used for y values.
        dfs (mafWrapper)        : An instace of a mafWrapper with MAF file loaded.
        y_function (function)   : An optional function used to transform the y values.
        show (bool)             : If matplotlib should show the plot.
        save (bool)             : If the plot should be written to file.
        filename (str)          : Full path to where the file should be stored and the file name.
    """

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

    adjust_text(texts)#,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    plt.xlabel("log2FC")
    plt.ylabel("-log P-value")
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
    """
    Function creates a MA-plot for data given.
    Args:
        df (pandas.DataFrame)   : A pandas dataframe containg all values that should be plottet.
        x (str)                 : Name of the column in the dataframe used for x values.
        y (str)                 : Name of column in the dataframe used for y values.
        dfs (mafWrapper)        : An instace of a mafWrapper with MAF file loaded.
        show (bool)             : If matplotlib should show the plot.
        save (bool)             : If the plot should be written to file.
        filename (str)          : Full path to where the file should be stored and the file name.
    """


    plt.figure()

    plt.scatter(x=np.log2(df[x]), y=df[y],
                s=1, label="Not significant")
    muts = dfs.get_df()
    # highlight down- or up- regulated genes
    down = df[(df[y] <= -2)]
    up = df[(df[y] >= 2)]
    plt.scatter(x=down[x].apply(lambda x: np.log2(x)), y=down[y], s=3, label="Down-regulated", color="blue")
    plt.scatter(x=up[x].apply(lambda x: np.log2(x)), y=up[y],
                s=3, label="Up-regulated", color="red")

    texts = []
    for i, r in up.iterrows():
        texts.append(plt.text(x=np.log2(r[x]), y=r[y], s=muts[muts["Entrez_Gene_Id"]
                     == int(r["geneName"])]["Hugo_Symbol"].unique()[0]))
    for i, r in down.iterrows():
        texts.append(plt.text(x=np.log2(r[x]), y=r[y], s=muts[muts["Entrez_Gene_Id"]
                     == int(r["geneName"])]["Hugo_Symbol"].unique()[0]))

    adjust_text(texts)#,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    plt.xlabel("Mean normalized expression")
    plt.ylabel("log2FC")
    mean_fdr = df[(df[y] != np.inf) & (df[y]  != -np.inf) & (~df[y].isna())][y].mean()
    plt.axhline(mean_fdr, color="red", linestyle="--",label=f"Mean {mean_fdr:.4f}")
    plt.axhline(2, color="grey", linestyle="--")
    plt.axhline(-2, color="grey", linestyle="--")
    plt.legend()
    plt.autoscale()

    if save:
        if filename:
            plt.savefig(filename)
        else:
            plt.savefig(f"vulcano_plot_{filename}.pdf")
        
        up.to_csv(os.path.join(os.path.split(filename)[0], "upregulated_genes.csv"), sep=";")
        down.to_csv(os.path.join(os.path.split(filename)[0], "downregulated_genes.csv"), sep=";")

    if show:
        plt.show()
        




def plot_gene_heatmap(geneName_hugo, dataframe, mutations, center=0,show = False, save = False, filename=None, sub_selection = None):
    """
    Plot the heatmap of expression levels for given gene.

    Args:
        geneName_hugo (str): Entrez gene Id
        dataframe (pd.DataFrame): A dataframe containing information about given gene
        mutations (pd.DataFrame): A dataframe containging information about 
                                  mutations for samples.
        center (int, str): Middle value for heatmap color scale.
        show (bool): If figures should be shown
        save (bool): If figures should be saved
        filename (str): Figure name, not used if save if False
    """
    geneName = mutations.get_entrez_id(geneName_hugo)[0]
    data = dataframe.loc[dataframe["geneName"] == str(geneName)]
    print(geneName)

    exons = data[data["type"] == "E"]
    junctions = data[data["type"] == "J"]

    fig, axs = plt.subplots(nrows=2,ncols=1, figsize=(10,10), gridspec_kw={"height_ratios": [3,3]})


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
#        axs[0].plot([row["start"],row["end"]], [count,count], c="b")
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
    #    axs[1].plot([row["start"], row["end"]], [count,count], c="r")
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
        y_ticks = [f"{j}{k+1}" for k in range(tmp_data.shape[0])]
        tmp_data["type"] = y_ticks

        if sub_selection:                        
            tmp_data = pd.DataFrame(tmp_data[(tmp_data["start"] >  sub_selection[0]) & (tmp_data["end"] < sub_selection[-1])])
        
        print(tmp_data)
        
        if center == "mean" : 
            g = sns.heatmap((tmp_data.iloc[:,12:]), ax=axs[i], yticklabels=tmp_data.iloc[:,5], cmap="Spectral_r", linecolor="k", linewidths=0.009)
        else:
            tmp_fig, axs[i] = plt.subplots(tmp_data.shape[0],1)
            for k in range(tmp_data.shape[0]):
                sns.heatmap((tmp_data.iloc[:,12:]), ax=axs[i+2][k], yticklabels=tmp_data.iloc[k,12:], center=0, cmap="Spectral_r", linecolor="k", linewidths=0.009)#, vmin=min(tmp_data.iloc[4,12:]), vmax=max(tmp_data.iloc[4,12:]))

        print(g)
        g.set_yticklabels(tmp_data["type"].to_numpy(), rotation=0)#[f"{j.get_text()}{i+1+sub_selection[0]}" for i,j in enumerate(g.get_yticklabels())], rotation=0)
            

        for tick_label in g.get_xticklabels():

            tick_text = tick_label.get_text()
            if not sample_muts[tick_text].empty:
                tick_label.set_color("darkred")

        for p, tick_label in enumerate(g.get_yticklabels()):
            if mut_pos[j][p]:
                tick_label.set_color("darkred")

    axs[0].set_xticklabels(["" for i in axs[0].get_xticklabels()])

    #for i in [0,1]:
    #    axs[i].ticklabel_format(useOffset=False, style='plain')

    if sub_selection:
        fig.suptitle(f"$\it{{ {mutations.get_df()[mutations.get_df()['Entrez_Gene_Id'] == geneName]['Hugo_Symbol'][0]} }}$, Genomic region: {tmp_data['seqnames'].to_numpy()[0]} {sub_selection[0]:,} to {sub_selection[-1]:,}")
    else:
        fig.suptitle(f"{mutations.get_df()[mutations.get_df()['Entrez_Gene_Id'] == geneName]['Hugo_Symbol'][0]}")
    
    if save:
        if filename:
            plt.savefig(filename)
        else:
            plt.savefig(os.path.join(filename, f"{geneName_hugo}_heatmap.pdf"))

    if show:
        plt.show()
        


def plot_binned_mutations(df, target_column, title, n_bins = 10, write=False, path=None, filetype="pdf"):
    """
    This function create and plots binned count of mutations. 
    Dynamicly creates bins depending on number of binns wanted
    Args:
        df (pandas.DataFrame): A data frame with sorted stacked count values.
                               This data frame need to be from the find_expression function.
        target_column (str) :  Name of column with count values
        title (str)         :  Title of the plot
        n_bins (int)        :  Number of bins to use, defaults to 10
        write (bool)        :  If the function should be written to file
        path (str)          :  Path to folder if the figure should be written to file
        filetype (str)      :  Filtype of plot. All matplotlib formats supported. Default is pdf. 
    """
    exons = df[df["type"] == "E"]
    junction = df[df["type"] == "J"]

    fig, ax = plt.subplots(2,1)


    skip = len(exons[target_column])//n_bins
    bins = [sum(exons[target_column][i:i+skip]) for i in range(0, len(exons[target_column]), skip)]
    ax[0].bar(x = range(len(bins)), height = bins)
    ax[0].set_title("Exon", y = -0.2)
    ax[0].set_xlim([-1,len(bins)-1])


    skip = len(junction[target_column])//n_bins
    bins = [sum(junction[target_column][i:i+skip]) for i in range(0, len(junction[target_column]), skip)]
    ax[1].bar(range(len(bins)), bins)
    ax[1].set_title("Junction", y=-0.2)
    ax[1].set_xlim([-1,len(bins)-1])

    labels = ["" for i in range(len(ax[1].get_xticklabels()))]
    ax[0].set_xticklabels(labels)
    labels[1] = "High expression"
    labels[-1] = "Low expression"


    ax[1].set_xticklabels(labels)
    fig.suptitle(title)

    if write:
        fig.savefig(os.path.join(path, f"distribution_mutated_features.{filetype}"))

    return fig, ax