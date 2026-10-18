import math

from spatialdata import SpatialData
import pandas as pd

def spatial_clonotype_distance (sdata: SpatialData, min_cells: int = 2):

    points = sdata.points["spatial"].compute()

    adata = sdata.tables["adata"]
    clonotypes = adata.obs["cc_aa_tcrdist"]

    groups = {}

    for cell_id, clonotype in clonotypes.items():
        if pd.isna(clonotype):
            continue

        if clonotype not in groups:
            groups[clonotype] = []

        groups[clonotype].append(cell_id)

    results = []

    for clonotype_id, cell_ids in groups.items():

        if (len(cell_ids) < min_cells):
            continue

        spatial_cell_ids = []

        for cell_id in cell_ids:
            if cell_id not in points.index:
                continue

            cell = points.loc[cell_id]

            if pd.isna(cell["x"]) or pd.isna(cell["y"]):
                continue

            spatial_cell_ids.append(cell_id)

        if len(spatial_cell_ids) < max(2, min_cells):
            continue

        sum_of_all_distances = 0
        number_of_distances = 0

        for i, first_cell_id in enumerate(spatial_cell_ids):

            for second_cell_id in spatial_cell_ids[i + 1:]:
                distance = calc_distance(
                    points,
                    first_cell_id,
                    second_cell_id
                )

                sum_of_all_distances += distance
                number_of_distances += 1

        mean_distance = (
            sum_of_all_distances / number_of_distances
        )

        results.append({
            "clonotype_id": clonotype_id,
            "n_cells": len(cell_ids),
            "n_cells_spatial": len(spatial_cell_ids),
            "n_distances": number_of_distances,
            "mean_pairwise_distance": mean_distance
        })

    result_df = pd.DataFrame(
        results,
        columns=[
            "clonotype_id",
            "n_cells",
            "n_cells_spatial",
            "n_distances",
            "mean_pairwise_distance"
        ]
    )
    result_df = result_df.set_index("clonotype_id")

    adata.uns["spatial_clonotype_distance"] = result_df

    return result_df

def calc_distance(
    points: pd.DataFrame,
    cell_id_1: str,
    cell_id_2: str
) -> float:

    cell_1 = points.loc[cell_id_1]
    cell_2 = points.loc[cell_id_2]

    dx = cell_1["x"] - cell_2["x"]
    dy = cell_1["y"] - cell_2["y"]

    return math.sqrt(dx ** 2 + dy ** 2)
