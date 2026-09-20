import matplotlib.pyplot as plt
from helper import melanoma_data_loader as mda

import scirpy as ir

print("hello scirpy")
adata = mda.import_melanoma_data()

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

ir.pl.plot_spatial(adata)
plt.show()



