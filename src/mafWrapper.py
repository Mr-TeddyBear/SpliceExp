import pandas as pd
import glob
import os

class mafWrapper:
    def __init__(self, path: str):
        """
        This class is a tool for finding mutations in genes. Reads in MAF files

        Args:
            path (str) -- Path to folder containing MAF files and a map file
                          Map file contains a mapping between sample name and
                          filename. Map file should contain two columns (filename, samplename) and be tab sperated. 

        """
        self.mapMAF = pd.read_csv(os.path.join(path,"map.txt"), sep="    ", header=None, names=["file", "barcode"])
        dfs = {}
        files = glob.glob(os.path.join(f"{path}","*","*.maf.gz"))
        
        for i in files:
            bcode = (self.mapMAF.loc[self.mapMAF["file"] == os.path.split(i)[-1]]["barcode"]).values[0]
            dfs[bcode] = pd.read_csv(i, delimiter="\t", comment="#")
        
        self.dfs = dfs
        
    def __getitem__(self, params):
        tmp_rtrn = {}
        for i in self.dfs:
            tmp_rtrn[i] = self.dfs[i][params]

        return tmp_rtrn
    
    def get_df(self) -> pd.DataFrame:
        """
        Returns a concatatantion of all sub dataframe for all samples.
        """
        return pd.concat(self.dfs)
    
    def get_geneID(self, geneID: int) -> dict:
        """
        Get all features from all samples for a given gene

        Args:
            geneID (str, int): Entrez Gene Id
        """
        tmp_rtrn = {}
        for i in self.dfs:
            tmp_rtrn[i] = self.dfs[i][self.dfs[i]["Entrez_Gene_Id"] == int(geneID)]

        return tmp_rtrn


    def get_mutated_samples(self, geneID : int) -> list:
        """
        Returns a list of all samples that have a mutation in geneID

        Args:
            geneID (int) : Entrez gene Id

        Returns:
            list : list of mutated samples
        
        """
        tmp_rtrn = []
        for i in self.dfs:
            tmp = self.dfs[i][self.dfs[i]["Entrez_Gene_Id"] == int(geneID)]
            if not tmp.empty:
                tmp_rtrn.append(i)
        

        return tmp_rtrn

    def is_mutated(self, geneID: int, sample: str) -> bool:
        """
        Boolean method to test of a sample has a mutation in geneID

        Args:
            geneID (int) : Entrez gene Id
            sample (str) : Name of sample to test

        Returns:
            bool: True or false

        """
        for i in (str(geneID)).split(","):
            if not self.dfs[sample][self.dfs[sample]["Entrez_Gene_Id"] == int(i)].empty:
                return True
            else:
                return False
        return False

    def is_mutation_close_to_feature(self, geneID: int, sample: str, feature: pd.DataFrame| pd.Series, out_of_feature_boundry: int = 10) -> bool:
        """
        Boolean function to test if a mutation is close a given feature in sample


        Args:
            geneID (int) : Entrez gene Id
            sample (str) : Name of sample to test
            feature (pd.DataFrame, pd.Series) : row from pandas dataframe
                                                out_of_feature_boundry (int): how many positions a mutation can be from the edge of a feature to be counted as close.

        Returns:
            bool: True or False

        """
        for i in (str(geneID)).split(","):
            if not (tmp := self.dfs[sample][self.dfs[sample]["Entrez_Gene_Id"] == int(i)]).empty:
                for s,e in zip(tmp["Start_Position"], tmp["End_Position"]):
                    if feature["start"] - out_of_feature_boundry < s  and  e < out_of_feature_boundry + feature["end"]:
                        return True
        return False

    def get_information_about_mutations(self, geneID: int, sample: str, feature: pd.DataFrame | pd.Series, out_of_feature_boundry: int = 10):
        """
        Args:
            geneID (int) : Entrez gene Id
            sample (str) : Name of sample to test
            feature (pd.DataFrame, pd.Series) : row from pandas dataframe
                                                out_of_feature_boundry (int): how many positions a
                                                mutation can be from the edge of a feature to be counted
                                                as close.

        Returns:
            (None, pd.DataFrame): Returns all mutatuions close to feature in sample, None if there are
                                  none mutations.

       
        """
        for i in (str(geneID)).split(","):
            if not (tmp := self.dfs[sample][self.dfs[sample]["Entrez_Gene_Id"] == int(i)]).empty:
                for i,(s,e) in enumerate(zip(tmp["Start_Position"], tmp["End_Position"])):
                    if feature["start"] - out_of_feature_boundry < s  and  e < out_of_feature_boundry + feature["end"]:
                        return tmp.iloc[[i]]
        return None
    

    def get_Hugo_symbol(self, entrez_id):
        """
        Converts Entrez id to Hugo symbol.
        Args:
            entrez_id (str): The entrez id for a gene, will only return Hugo symbol for
                             the gene if it is in the mutations files.
        """
        tmp = self.get_df()
        return tmp[tmp["Entrez_Gene_Id"] == int(entrez_id)]["Hugo_Symbol"].unique()
    
    def get_entrez_id(self, hugo_symbol):
        """
        Converts Hugo symbol to Entrez Id.
        Args:
            hugo_symbol (str): Hugo symbol for a gene found in the mutation files.
        """
        tmp = self.get_df()
        return tmp[tmp["Hugo_Symbol"] == hugo_symbol]["Entrez_Gene_Id"].unique()
