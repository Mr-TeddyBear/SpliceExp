import matplotlib.pyplot as plt
import numpy as np

def ma_plot(df, x, y, dfs, show=False, save=False, filename=None):

    plt.figure()

    plt.scatter(x=df[x], y=df[y],
                s=1, label="Not significant")
    muts = dfs.get_df()
    # highlight down- or up- regulated genes
    #down = df[(df[x] <= -2) & (df[y] <= 0.01)]
    #up = df[(df[x] >= 2) & (df[y] <= 0.01)]
    #print("\n\n\n", up.head(), "\n\n\n")
    #print("\n\n\n", down.head(), "\n\n\n")
    #plt.scatter(x=down[x], y=down[y].apply(
    #    lambda x: -np.log10(x)), s=3, label="Down-regulated", color="blue")
    #plt.scatter(x=up[x], y=up[y].apply(lambda x: -np.log10(x)),
                s=3, label="Up-regulated", color="red")

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
    