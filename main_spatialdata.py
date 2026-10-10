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

ir.tl.spatial_clonotype_distance(sdata)

adata = sdata.tables["adata"]
result = adata.uns["spatial_clonotype_distance"]

plot_data = (
    result[
        (result["n_cells_spatial"] >= 5)
        & result["mean_pairwise_distance"].notna()
    ]
    .sort_values("mean_pairwise_distance")
)

fig, ax = plt.subplots(figsize=(12, 5))

bars = ax.bar(
    plot_data.index.astype(str),
    plot_data["mean_pairwise_distance"],
)
ax.bar_label(
    bars,
    labels=[f"n={int(n)}" for n in plot_data["n_cells_spatial"]],
    padding=3,
    rotation=90,
    fontsize=8,
)

ax.set_xlabel("Clonotype cluster")
# The reader does not specify a physical unit for the supplied X/Y coordinates.
ax.set_ylabel("Mean pairwise distance (coordinate units)")
ax.set_title("Mean spatial distance within TCR clonotypes\nn = spatially mapped cells")
ax.tick_params(axis="x", labelrotation=90)
ax.margins(y=0.2)
fig.tight_layout()

plt.show()

# adata = sdata.tables["adata"]
#
# result = adata.uns["spatial_clonotype_distance"]
# result = adata.uns["spatial_clonotype_distance"]
#
# print(
#     result[
#         [
#             "n_cells",
#             "n_cells_spatial",
#             "mean_nearest_neighbor_distance",
#         ]
#     ].head(20)
# )

# ir.pl.plot_spatial(
#     sdata,
#     color="cc_aa_tcrdist",
#     filter_by="cell_type",
#     filter_value="tumour_1",
# )
#
#
# ir.pl.plot_spatial_all(
#     sdata,
#     color="cc_aa_tcrdist",
#     filter_by="cell_type",
#     filter_value="tumour_1",
# )
#
# plt.show()


# adata = sdata.tables["adata"]
# print (adata)
#
# print(adata.obs["cc_aa_tcrdist"].head())
# print(adata.obsm["spatial"][:5])
