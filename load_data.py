import scirpy as ir
import tarfile
import warnings
from glob import glob

import anndata
import matplotlib.pyplot as plt
import muon as mu
import pandas as pd
import scanpy as sc
import scirpy as ir


# mdata = ir.datasets.wu2020_3k()
# print(mdata)

# airr_data = mdata["airr"].obsm["airr"]
#
# print(type(airr_data))
# print(airr_data)
#
# cell0 = mdata["airr"].obsm["airr"][0]
# print(cell0[0].fields)

# gex = mdata["gex"]

# print(gex.obs.head())
# print(gex.var.head())

# cell_id = gex.obs_names[0]
#
# print(cell_id)
# print(mdata["gex"].obs.loc[cell_id])
# print ("---------")
# print(mdata["airr"].obs.loc[cell_id])
#
# cell_id = "LN1_GTAGGCCAGCGTAGTG-1"
#
# for gene in ["CD4", "CD8A", "FOXP3", "GZMB"]:
#     value = mdata["gex"][cell_id, gene].X
#     print(gene, value)

# for field in [
#     "locus",
#     "v_call",
#     "d_call",
#     "j_call",
#     "c_call",

#     "cdr3_aa",
#     "productive"
# ]:
#     print(field, ":", cell0[0][field])

# gex = mdata["gex"]
#
# print(gex.obs.head())
# print(gex.var.head())

import tarfile
import warnings
from glob import glob


sc.set_figure_params(figsize=(4, 4))
sc.settings.verbosity = 2  # verbosity: errors (0), warnings (1), info (2), hints (3)

adata_tcr = ir.io.read_10x_vdj("example_data/liao-2019-covid19/GSM4385993_C144_filtered_contig_annotations.csv.gz")

# Load the associated transcriptomics data
adata = sc.read_10x_h5("example_data/liao-2019-covid19/GSM4339772_C144_filtered_feature_bc_matrix.h5")
adata.var_names_make_unique()

print (adata)
