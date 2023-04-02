import numpy as np


def calc_fold_change(gp1, gp2):
    gp1 = [np.log2(i) if i != 0 else 1 for i in gp1.to_numpy().ravel()]
    gp2 = [np.log2(i) if i != 0 else 1 for i in gp2.to_numpy().ravel()]
    return np.average(gp1) - np.average(gp2)
