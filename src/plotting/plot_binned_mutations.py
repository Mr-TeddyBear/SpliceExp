import matplotlib.pyplot as plt
import os


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