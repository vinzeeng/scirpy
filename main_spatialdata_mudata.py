from helper.slide_tags_reader import read_slidetags
from helper.spatialdata_scirpy_converter import spatialdata_to_anndata
from helper.spatialdata_scirpy_converter import spatialdata_to_mudata

import scirpy as ir
import matplotlib.pyplot as plt

print("hello scirpy spatial mutdata")

sdata = read_slidetags("SCP2176")
mdata     = spatialdata_to_mudata(sdata)

ir.pp.index_chains(
    mdata,
    filter=["require_junction_aa"]
)

ir.tl.chain_qc(mdata)

ir.pp.ir_dist(
    mdata,
    sequence="aa",
    metric="tcrdist"
)

ir.tl.define_clonotype_clusters(
    mdata,
    sequence="aa",
    metric="tcrdist",
    receptor_arms="any",
    dual_ir="primary_only"
)

ir.pl.plot_spatial(mdata)
plt.show()
