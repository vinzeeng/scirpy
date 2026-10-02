from helper.slide_tags_reader import read_slidetags
from helper.spatialdata_scirpy_converter import spatialdata_to_anndata

import scirpy as ir
import matplotlib.pyplot as plt

print("hello scirpy spatial data")

sdata = read_slidetags("SCP2176")
adata     = spatialdata_to_anndata(sdata)

ir.pp.index_chains(
    adata,
    filter=["require_junction_aa"]
)

ir.tl.chain_qc(adata)

ir.pp.ir_dist(
    adata,
    sequence="aa",
    metric="tcrdist"
)

ir.tl.define_clonotype_clusters(
    adata,
    sequence="aa",
    metric="tcrdist",
    receptor_arms="any",
    dual_ir="primary_only"
)

# ir.pl.plot_spatial(
#     adata,
#     color="cc_aa_tcrdist",
#     filter_by="cell_type",
#     filter_value="tumour_1",
# )

ir.pl.plot_spatial_all(
    adata,
    color="cc_aa_tcrdist",
    filter_by="cell_type",
    filter_value="tumour_1",
)



plt.show()
