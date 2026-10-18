import matplotlib.pyplot as plt
from spatialdata import SpatialData


def plot_spatial_clonotype_distance(
    sdata: SpatialData,
    top_n: int = 20
):
    adata = sdata.tables["adata"]

    distances = adata.uns["spatial_clonotype_distance"]

    # Nur die 20 Cluster mit den meisten räumlichen Zellen
    distances = distances.sort_values(
        "n_cells_spatial",
        ascending=False
    ).head(top_n)

    fig, ax = plt.subplots(figsize=(12, 6))

    ax.bar(
        distances.index.astype(str),
        distances["mean_pairwise_distance"]
    )

    ax.set_xlabel("Clonotype Cluster")
    ax.set_ylabel("Mean Pairwise Distance (µm)")
    ax.set_title("Spatial Distribution of TCR Clonotype Clusters")

    ax.tick_params(axis="x", rotation=90)

    plt.tight_layout()
    plt.show()

import matplotlib.pyplot as plt
from spatialdata import SpatialData



import matplotlib.pyplot as plt
from spatialdata import SpatialData


def plot_spatial_data(
    sdata: SpatialData,
    cell_type: str | None = None
):
    points = sdata.points["spatial"].compute()

    # Optional nach einem Zelltyp filtern
    if cell_type is not None:
        points = points[points["cell_type"] == cell_type]

    fig, ax = plt.subplots(figsize=(12, 9))

    # Zellen nach Zelltyp gruppieren
    for type_name, cells in points.groupby("cell_type"):

        ax.scatter(
            cells["x"],
            cells["y"],
            s=8,
            alpha=0.7,
            label=str(type_name)
        )

    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
    ax.set_title("Spatial distribution by cell type")

    ax.set_aspect("equal")
    ax.invert_xaxis()
    ax.invert_yaxis()

    ax.grid(True, alpha=0.2)
    ax.legend(
        title="Cell types",
        bbox_to_anchor=(1.05, 1),
        loc="upper left"
    )

    plt.tight_layout()
    plt.show()
