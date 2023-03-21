from pybiomart import Dataset, Server


class mapGeneID:
    """
    Wrapper class for pybiomart.
    Makes finding switching between IDs easier. Currently got methods for
    ensembleID, HugoSymbol and entrezID
    """
    def __init__(self, filters, wanted_fileds=['ensembl_gene_id', 'external_gene_name', "entrezgene_id"]):
        self.wanted_feilds = wanted_fileds
        self.filters = filters
        self.geneMap = self._callBiomart(wanted_fileds, filters)        

    def _callBiomart(self, wanted_fileds, filters):
        dataset = Dataset(name='hsapiens_gene_ensembl',
                        host='http://www.ensembl.org/')

        geneMapID = dataset.query(attributes=wanted_fileds,
              filters=filters)
        return geneMapID

    def get_ensembleId(self, gID, from_type):
        match from_type:
            case "name":
                tmp = self.geneMap[self.geneMap["NCBI gene (formerly Entrezgene) ID"] == str(gID)]["Gene stable ID"]
            case "entrez":
                tmp = self.geneMap[self.geneMap["NCBI gene (formerly Entrezgene) ID"] == str(gID)]["Gene stable ID"]
        if not tmp.empty:
            return tmp[0]
        else:
            return None

    def get_geneName(self, gID, from_type):
        match from_type:
            case "entrez":
                tmp = self.geneMap[self.geneMap["NCBI gene (formerly Entrezgene) ID"] == str(gID)]["Gene name"]
            case "ensembl":
                tmp = self.geneMap[self.geneMap["Gene stable ID"] == str(gID)]["Gene name"]
        if not tmp.empty:
            return tmp[0]
        else:
            return None

    def get_entrezId(self, gID, from_type):
        match from_type:
            case "ensembl":
                tmp = self.geneMap[self.geneMap["Gene stable ID"] == str(gID)]["NCBI gene (formerly Entrezgene) ID"]
            case "name":
                tmp = self.geneMap[self.geneMap["Gene name"] == str(gID)]["NCBI gene (formerly Entrezgene) ID"]
        if not tmp.empty:
            return tmp[0]
        else:
            return None

    def __call__(self):
        return self.geneMap
    
    def __str__(self):
        string = ""

        for i in self.geneMap.iterrows():
            string += str(i)
        
        return string


class mapGeneID:
    """
    Wrapper class for pybiomart.
    Makes finding switching between IDs easier. Currently got methods for
    ensembleID, HugoSymbol and entrezID
    """
    def __init__(self, filters, wanted_fileds=['ensembl_gene_id', 'external_gene_name', "entrezgene_id"]):
        self.wanted_feilds = wanted_fileds
        self.filters = filters
        self.geneMap = self._callBiomart(wanted_fileds, filters)        

    def _callBiomart(self, wanted_fileds, filters):
        dataset = Dataset(name='hsapiens_gene_ensembl',
                        host='http://www.ensembl.org/')

        geneMapID = dataset.query(attributes=wanted_fileds,
              filters=filters)
        return geneMapID

    def get_ensembleId(self, gID, from_type):
        match from_type:
            case "name":
                tmp = self.geneMap[self.geneMap["NCBI gene (formerly Entrezgene) ID"] == str(gID)]["Gene stable ID"]
            case "entrez":
                tmp = self.geneMap[self.geneMap["NCBI gene (formerly Entrezgene) ID"] == str(gID)]["Gene stable ID"]
        if not tmp.empty:
            return tmp[0]
        else:
            return None

    def get_geneName(self, gID, from_type):
        match from_type:
            case "entrez":
                tmp = self.geneMap[self.geneMap["NCBI gene (formerly Entrezgene) ID"] == str(gID)]["Gene name"]
            case "ensembl":
                tmp = self.geneMap[self.geneMap["Gene stable ID"] == str(gID)]["Gene name"]
        if not tmp.empty:
            return tmp[0]
        else:
            return None

    def get_entrezId(self, gID, from_type):
        match from_type:
            case "ensembl":
                tmp = self.geneMap[self.geneMap["Gene stable ID"] == str(gID)]["NCBI gene (formerly Entrezgene) ID"]
            case "name":
                tmp = self.geneMap[self.geneMap["Gene name"] == str(gID)]["NCBI gene (formerly Entrezgene) ID"]
        if not tmp.empty:
            return tmp[0]
        else:
            return None

    def __call__(self):
        return self.geneMap
    
    def __str__(self):
        string = ""

        for i in self.geneMap.iterrows():
            string += str(i)
        
        return string


class mapGeneID:
    """
    Wrapper class for pybiomart.
    Makes finding switching between IDs easier. Currently got methods for
    ensembleID, HugoSymbol and entrezID
    """
    def __init__(self, filters, wanted_fileds=['ensembl_gene_id', 'external_gene_name', "entrezgene_id"]):
        self.wanted_feilds = wanted_fileds
        self.filters = filters
        self.geneMap = self._callBiomart(wanted_fileds, filters)        

    def _callBiomart(self, wanted_fileds, filters):
        dataset = Dataset(name='hsapiens_gene_ensembl',
                        host='http://www.ensembl.org/')

        geneMapID = dataset.query(attributes=wanted_fileds,
              filters=filters)
        return geneMapID

    def get_ensembleId(self, gID, from_type):
        match from_type:
            case "name":
                tmp = self.geneMap[self.geneMap["NCBI gene (formerly Entrezgene) ID"] == str(gID)]["Gene stable ID"]
            case "entrez":
                tmp = self.geneMap[self.geneMap["NCBI gene (formerly Entrezgene) ID"] == str(gID)]["Gene stable ID"]
        if not tmp.empty:
            return tmp[0]
        else:
            return None

    def get_geneName(self, gID, from_type):
        match from_type:
            case "entrez":
                tmp = self.geneMap[self.geneMap["NCBI gene (formerly Entrezgene) ID"] == str(gID)]["Gene name"]
            case "ensembl":
                tmp = self.geneMap[self.geneMap["Gene stable ID"] == str(gID)]["Gene name"]
        if not tmp.empty:
            return tmp[0]
        else:
            return None

    def get_entrezId(self, gID, from_type):
        match from_type:
            case "ensembl":
                tmp = self.geneMap[self.geneMap["Gene stable ID"] == str(gID)]["NCBI gene (formerly Entrezgene) ID"]
            case "name":
                tmp = self.geneMap[self.geneMap["Gene name"] == str(gID)]["NCBI gene (formerly Entrezgene) ID"]
        if not tmp.empty:
            return tmp[0]
        else:
            return None

    def __call__(self):
        return self.geneMap
    
    def __str__(self):
        string = ""

        for i in self.geneMap.iterrows():
            string += str(i)
        
        return string