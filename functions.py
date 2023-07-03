import numpy as np
import scanpy as sc

def quality_metrics(adata, gene_subtypes=["mt"]):
    
    #qc metrics
    sc.pp.calculate_qc_metrics(adata,percent_top=(3,),inplace=True)

    #mt fraction
    for gene_subtype in gene_subtypes:
        adata.obs[gene_subtype+"_fraction"] = np.array(adata.X[:,[i.startswith(gene_subtype) for i in adata.var["Gene"].index.values]]).mean(dim=1)