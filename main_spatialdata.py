import matplotlib.pyplot as plt
import squidpy as sq



from helper.slide_tags_reader import read_slidetags

import scirpy as ir


print("Scirpy analysis directly on SpatialData")


sdata = read_slidetags("SCP2176")

ir.pp.index_chains(
    sdata,
    filter=["require_junction_aa"],
)

ir.tl.chain_qc(
    sdata
)

ir.pp.ir_dist(
    sdata,
    sequence="aa",
    metric="tcrdist",
)

ir.tl.define_clonotype_clusters(
    sdata,
    sequence="aa",
    metric="tcrdist",
    receptor_arms="any",
    dual_ir="primary_only",
)

plt.figure()

ir.pl.plot_spatial(
    sdata,
    color="cc_aa_tcrdist",
    filter_by="cell_type",
    filter_value="tumour_1",
)

plt.show()

plt.figure()

ir.pl.plot_spatial_all(
    sdata,
    color="cc_aa_tcrdist",
    filter_by="cell_type",
    filter_value="tumour_1",
)

plt.show()


adata = sdata.tables["adata"]

print(adata)
