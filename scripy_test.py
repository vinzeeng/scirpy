import muon as mu
import numpy as np
import scanpy as sc
import scirpy as ir
import seaborn as sns
import pandas as pd#
from cycler import cycler
from matplotlib import cm as mpl_cm
from matplotlib import pyplot as plt



# temporary fix for deprecated matplotlib functionality
import IPython.display
from matplotlib_inline.backend_inline import set_matplotlib_formats

IPython.display.set_matplotlib_formats = set_matplotlib_formats

sc.set_figure_params(figsize=(4, 4))
sc.settings.verbosity = 2  # verbosity: errors (0), warnings (1), info (2), hints (3)

print("Hello Analysis of 3k T cells from cancer#")
sc.logging.print_header()

mdata = ir.datasets.wu2020_3k()
#print(mdata)

airr = mdata["airr"]
#print(airr)

tcr = airr.obsm["airr"]


print("tcr")

for i in range(5):
     print("Cell", i)
     for i in range(5):
         print("Cell ID:", airr.obs_names[i])
         print(tcr[i].to_list())

df = pd.DataFrame([
    {
        "cell_id": "AAAGCAAGTACGGTAC-1",
        "locus": "TRA",
        "junction_aa": "CAYGTGNQFYF",
    },
    {
        "cell_id": "AAAGCAAGTACGGTAC-1",
        "locus": "TRB",
        "junction_aa": "CSAVGTGLPRDIQYF",
    },
])

adata_tcr = ir.io.read_airr(df)

#print(adata_tcr)
#print(adata_tcr.obsm["airr"][0].to_list())
