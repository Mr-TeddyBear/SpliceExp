import pandas as pd

up = pd.read_csv("downregulated_genes.csv", sep=";")
down = pd.read_csv("upregulated_genes.csv", sep=";")




with pd.option_context('mode.use_inf_as_na',True):
    up = (up[~up["fold_change"].isna()].iloc[2:])
    down = down[~down["fold_change"].isna()]


    (down.style.format(precision=2)).to_latex("upregulated_table.tex", column_format="l|r|r|r|r|l|r|r", label="tab:singel_upregulated", hrules=True)
    (up.style.format(precision=2)).to_latex("downregulated_table.tex", column_format="l|r|r|r|r|l|r|r", label="tab:single_downregulated", hrules=True)
