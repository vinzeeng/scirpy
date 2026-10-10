import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist

from scirpy.util import DataHandler


def spatial_clonotype_distance(
    data: DataHandler.TYPE,
    *,
    clonotype_key: str = "cc_aa_tcrdist",
    spatial_key: str = "spatial",
    key_added: str = "spatial_clonotype_distance",
) -> None:
    """Describe the average spatial distance within each TCR clonotype cluster.

    The primary metric, ``mean_pairwise_distance``, is the average Euclidean
    distance between all unordered cell pairs of the same clonotype with valid
    spatial coordinates. Smaller values describe a more compact spatial
    distribution; larger values describe a more dispersed distribution.
    ``mean_nearest_neighbor_distance`` provides an additional measure of local
    density within the clonotype.

    Parameters
    ----------
    data
        AnnData, MuData, SpatialData, or DataHandler.
        The AnnData representation selected by DataHandler is used.
    clonotype_key
        Column in ``adata.obs`` containing clonotype assignments.
        Cells with missing clonotype assignments are ignored.
    spatial_key
        Key in ``adata.obsm`` containing numeric x/y coordinates of shape
        ``(n_cells, 2)``.

        Cells with missing or non-finite spatial coordinates are excluded
        from the spatial distance calculation.
    key_added
        Key in ``adata.uns`` where the resulting per-clonotype table is stored.

    Returns
    -------
    None
        Stores the per-clonotype DataFrame in ``adata.uns[key_added]`` on the
        AnnData representation selected by DataHandler (AIRR modality for MuData).
        SpatialData results remain in its original table.

    Notes
    -----
    The result contains one row per observed clonotype.

    ``n_cells``
        Total number of cells assigned to the clonotype.

    ``n_cells_spatial``
        Number of clonotype cells with valid spatial coordinates.

    ``mean_pairwise_distance``
        Primary metric: mean Euclidean distance across all unordered pairs
        of spatially mapped cells belonging to the clonotype. Each pair is
        counted once; self-pairs are excluded.

    ``median_pairwise_distance``
        Median Euclidean distance across all unordered pairs of spatially
        mapped cells belonging to the clonotype.

    ``mean_nearest_neighbor_distance``
        Additional local density metric: mean distance from each spatially
        mapped cell to its nearest spatially mapped cell of the same clonotype.
        Smaller values indicate closer local neighbors.

    Distance values are ``NaN`` when fewer than two cells with valid spatial
    coordinates are available.

    Distances are expressed in the units of the supplied spatial coordinates.
    Pairwise summaries require quadratic time and memory in the number of mapped
    cells per clonotype. Nearest neighbors use a spatial tree without allocating
    a square distance matrix. No statistical tests are performed.
    """
    adata = DataHandler(data, "airr").adata

    if clonotype_key not in adata.obs.columns:
        raise KeyError(f"{clonotype_key!r} not found in adata.obs.")

    if spatial_key not in adata.obsm:
        raise KeyError(f"{spatial_key!r} not found in adata.obsm.")

    clonotypes = adata.obs[clonotype_key]

    coordinates = np.asarray(adata.obsm[spatial_key])

    if coordinates.ndim != 2 or coordinates.shape != (adata.n_obs, 2):
        raise ValueError(f"adata.obsm[{spatial_key!r}] must have shape (n_cells, 2).")

    try:
        coordinates = pd.DataFrame(coordinates).apply(pd.to_numeric).to_numpy(dtype=float, na_value=np.nan)
    except (TypeError, ValueError) as e:
        raise ValueError("Spatial coordinates must be numeric.") from e

    assigned = clonotypes.notna().to_numpy()

    groups = pd.Series(
        clonotypes.to_numpy()[assigned],
        index=np.flatnonzero(assigned),
    )

    rows = []
    labels = []

    for label, cell_indices in groups.groupby(
        groups,
        sort=False,
        observed=True,
    ).groups.items():
        cell_indices = np.asarray(
            list(cell_indices),
            dtype=int,
        )

        n_cells = len(cell_indices)

        clonotype_coordinates = coordinates[cell_indices]

        valid_coordinates = np.isfinite(clonotype_coordinates).all(axis=1)

        clonotype_coordinates = clonotype_coordinates[valid_coordinates]

        n_cells_spatial = len(clonotype_coordinates)

        mean_pairwise_distance = np.nan
        median_pairwise_distance = np.nan
        mean_nearest_neighbor_distance = np.nan

        if n_cells_spatial >= 2:
            distances = pdist(
                clonotype_coordinates,
                metric="euclidean",
            )

            mean_pairwise_distance = distances.mean()

            median_pairwise_distance = np.median(distances)

            # k=2 excludes the query cell itself; coincident cells still have distance zero.
            nearest_neighbor_distances = cKDTree(clonotype_coordinates).query(clonotype_coordinates, k=2)[0][:, 1]
            mean_nearest_neighbor_distance = nearest_neighbor_distances.mean()

        labels.append(label)

        rows.append(
            (
                n_cells,
                n_cells_spatial,
                mean_pairwise_distance,
                median_pairwise_distance,
                mean_nearest_neighbor_distance,
            )
        )

    result = pd.DataFrame(
        rows,
        index=pd.Index(
            labels,
            name=clonotype_key,
        ),
        columns=[
            "n_cells",
            "n_cells_spatial",
            "mean_pairwise_distance",
            "median_pairwise_distance",
            "mean_nearest_neighbor_distance",
        ],
    )

    result = result.astype(
        {
            "n_cells": int,
            "n_cells_spatial": int,
            "mean_pairwise_distance": float,
            "median_pairwise_distance": float,
            "mean_nearest_neighbor_distance": float,
        }
    )

    adata.uns[key_added] = result
