from plotting.plot_binned_mutations import plot_binned_mutations
import pandas as pd
import matplotlib.pyplot as plt

def ravel_dataframe(df, columns_to_flatten, expression_column_name = "expression"):
        """
        This function stacks all the reads in to a single columns.
        """
        sub_df = []
        for i in columns_to_flatten:
                tmp = df.iloc[:,:12]
                tmp[expression_column_name] = df[i]
                tmp["sample"] = i

                sub_df.append(tmp)
        return(pd.concat(sub_df))


def is_mutated(df, dfs):
    """
    Create a columns in a ravel dataframe that indicates if the gene the feature is in have an mutation.

    Args:
        df (pd.DataFrame): A raveld dataframe
    """
    df["is_mutated"] = [dfs.is_mutated(row[1]["geneName"], row[1]["sample"]) for row in df.iterrows()]
    df.astype({"is_mutated": bool})

def is_related_mutation(df, dfs):
    """
    Create a column in a ravel dataframe that indicates if the feature has a mutation in it.

    Args:
        df (pd.DataFrame): A raveld dataframe
    """

    df["is_related_mutated"] = [dfs.is_mutation_close_to_feature(row[1]["geneName"], row[1]["sample"], row[1]) for row in df.iterrows()]
    df.astype({"is_related_mutated": bool})




def find_expressions(df, dfs, raw_df, col, show, title1 = "Mutated samples", title2 = "Mutation is close to feature"):
    df = df.copy()
    df = ravel_dataframe(df, col)

    stacked_raw_expression = ravel_dataframe(raw_df, col, "raw_expression")

    #display(df.head())
    #display(stacked_raw_expression.head())

    df = pd.concat([df, stacked_raw_expression["raw_expression"]], ignore_index=True)
    
    df.sort_values("expression", inplace=True, ascending=False)
    df = df[df["expression"].notna()]
    
    is_mutated(df, dfs)
    is_related_mutation(df, dfs)
    if show:
        plot_binned_mutations(df, "is_mutated", title1)
        plt.show()
        plot_binned_mutations(df, "is_related_mutated", title2)
        plt.show()
    


    return df