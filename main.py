import scirpy as ir

from spatial_data_reader import read_spatial_data

print("Hello Spatial Data in Scirpy")

sdata = read_spatial_data("SCP2176")

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

ir.tl.spatial_clonotype_distance(sdata)

ir.pl.plot_spatial()
