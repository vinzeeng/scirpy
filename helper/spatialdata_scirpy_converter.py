import pandas as pd
import scirpy as ir

from anndata import AnnData
from mudata import MuData
from spatialdata import SpatialData


def spatialdata_to_anndata(
    sdata: SpatialData,
    table_key: str | None = None,
) -> AnnData:
    """
    Convert SpatialData into a Scirpy-compatible AnnData object.

    Spatial coordinates are stored in obsm["spatial"].

    If alpha/beta TCR sequences are available in obs,
    they are converted into Scirpy AIRR format and stored
    in obsm["airr"].
    """

    adata = _get_table(
        sdata,
        table_key,
    ).copy()

    spatial = _get_spatial_coordinates(
        sdata,
        adata,
    )

    adata.obsm["spatial"] = (
        spatial.to_numpy()
    )

    if _contains_tcr_data(adata):

        airr_df = _create_airr_dataframe(
            adata.obs
        )

        adata_airr = ir.io.read_airr(
            airr_df
        )

        barcode_to_index = {
            barcode: index
            for index, barcode
            in enumerate(
                adata_airr.obs_names
            )
        }

        airr_index = [
            barcode_to_index.get(
                barcode,
                None,
            )
            for barcode
            in adata.obs_names
        ]

        adata.obsm["airr"] = (
            adata_airr.obsm["airr"][
                airr_index
            ]
        )

    return adata


def spatialdata_to_mudata(
    sdata: SpatialData,
    table_key: str | None = None,
) -> MuData:
    """
    Convert SpatialData into a Scirpy-compatible MuData object.

    Result
    ------
    MuData
    ├── mod["gex"]
    │   └── obsm["spatial"]
    │
    └── mod["airr"]
        └── obsm["airr"]

    If no TCR data are available, only the gex modality
    is created.
    """

    gex = _get_table(
        sdata,
        table_key,
    ).copy()

    spatial = _get_spatial_coordinates(
        sdata,
        gex,
    )

    gex.obsm["spatial"] = (
        spatial.to_numpy()
    )

    if not _contains_tcr_data(gex):

        return MuData(
            {
                "gex": gex,
            }
        )

    airr_df = _create_airr_dataframe(
        gex.obs
    )

    airr = ir.io.read_airr(
        airr_df
    )

    return MuData(
        {
            "gex": gex,
            "airr": airr,
        }
    )


def _get_table(
    sdata: SpatialData,
    table_key: str | None,
) -> AnnData:
    """
    Return a table from SpatialData.

    If exactly one table exists, it is selected
    automatically.
    """

    if table_key is not None:

        if table_key not in sdata.tables:
            raise KeyError(
                f"Table '{table_key}' not found."
            )

        return sdata.tables[
            table_key
        ]

    if len(sdata.tables) == 1:

        return next(
            iter(
                sdata.tables.values()
            )
        )

    raise ValueError(
        "SpatialData contains multiple tables. "
        "Please specify table_key."
    )


def _get_spatial_coordinates(
    sdata: SpatialData,
    adata: AnnData,
) -> pd.DataFrame:
    """
    Extract spatial x/y coordinates.

    Supported SpatialData representations:
        - Shapes linked through spatialdata_attrs
        - Points linked through observation IDs
    """

    # --------------------------------------------------
    # Shapes
    #
    # Example: Visium
    #
    # tables["adata"]
    #       |
    #       | region / instance_key
    #       v
    # shapes["spots"]
    # --------------------------------------------------

    if "spatialdata_attrs" in adata.uns:

        attrs = adata.uns[
            "spatialdata_attrs"
        ]

        region = attrs.get(
            "region"
        )

        instance_key = attrs.get(
            "instance_key"
        )

        if region is not None:

            if isinstance(
                region,
                (list, tuple),
            ):

                if len(region) != 1:
                    raise NotImplementedError(
                        "Multiple spatial regions "
                        "are not supported yet."
                    )

                region = region[0]

            if region in sdata.shapes:

                return _get_coordinates_from_shapes(
                    sdata,
                    adata,
                    region,
                    instance_key,
                )

    # --------------------------------------------------
    # Points
    #
    # Example: Slide-tags
    #
    # points["spatial"]
    #     NAME
    #     X
    #     Y
    #     cell_type
    # --------------------------------------------------

    if len(sdata.points) == 1:

        return _get_coordinates_from_points(
            sdata,
            adata,
        )

    raise ValueError(
        "Could not determine spatial coordinates "
        "from SpatialData."
    )


def _get_coordinates_from_shapes(
    sdata: SpatialData,
    adata: AnnData,
    region: str,
    instance_key: str | None,
) -> pd.DataFrame:
    """
    Extract x/y coordinates from a SpatialData Shapes element.

    For polygons and circles, their centroid is used
    as the spatial coordinate.
    """

    if instance_key is None:
        raise ValueError(
            "SpatialData table does not define "
            "an instance_key."
        )

    if instance_key not in adata.obs.columns:
        raise ValueError(
            f"Instance key '{instance_key}' "
            f"not found in adata.obs."
        )

    shapes = sdata.shapes[
        region
    ]

    centroids = (
        shapes
        .geometry
        .centroid
    )

    coordinates = pd.DataFrame(
        {
            "X": centroids.x,
            "Y": centroids.y,
        },
        index=shapes.index,
    )

    instance_ids = adata.obs[
        instance_key
    ]

    spatial = coordinates.reindex(
        instance_ids
    )

    # The reindex operation currently uses the spatial
    # instance IDs. Restore the AnnData observation index.
    spatial.index = (
        adata.obs_names
    )

    return spatial


def _get_coordinates_from_points(
    sdata: SpatialData,
    adata: AnnData,
) -> pd.DataFrame:
    """
    Extract x/y coordinates from a SpatialData Points element.

    Supports the original Slide-tags column names:

        NAME
        X
        Y

    and an already normalized representation:

        cell_id
        x
        y

    Additional point annotations such as cell_type are
    copied into adata.obs.
    """

    points = next(
        iter(
            sdata.points.values()
        )
    ).compute()

    # --------------------------------------------------
    # Determine observation ID column
    # --------------------------------------------------

    if "NAME" in points.columns:

        id_column = "NAME"

    elif "cell_id" in points.columns:

        id_column = "cell_id"

    else:

        raise ValueError(
            "Could not determine the observation ID "
            "column of the SpatialData Points element."
        )

    # --------------------------------------------------
    # Determine coordinate columns
    # --------------------------------------------------

    if (
        "X" in points.columns
        and "Y" in points.columns
    ):

        x_column = "X"
        y_column = "Y"

    elif (
        "x" in points.columns
        and "y" in points.columns
    ):

        x_column = "x"
        y_column = "y"

    else:

        raise ValueError(
            "Could not determine x/y coordinate columns "
            "of the SpatialData Points element."
        )

    # --------------------------------------------------
    # Validate IDs
    # --------------------------------------------------

    if points[id_column].duplicated().any():

        raise ValueError(
            "SpatialData Points contain duplicate "
            "observation IDs."
        )

    points = points.set_index(
        id_column
    )

    # --------------------------------------------------
    # Spatial coordinates
    # --------------------------------------------------

    spatial = points[
        [
            x_column,
            y_column,
        ]
    ].reindex(
        adata.obs_names
    )

    spatial.columns = [
        "X",
        "Y",
    ]

    # --------------------------------------------------
    # Additional point annotations -> adata.obs
    #
    # Example:
    # cell_type
    # --------------------------------------------------

    annotation_columns = [
        column
        for column in points.columns
        if column not in {
            x_column,
            y_column,
        }
    ]

    for column in annotation_columns:

        adata.obs[column] = (
            points[column]
            .reindex(
                adata.obs_names
            )
        )

    return spatial


def _contains_tcr_data(
    adata: AnnData,
) -> bool:
    """
    Check whether raw alpha/beta TCR sequences
    are available.
    """

    return (
        "alpha" in adata.obs.columns
        or "beta" in adata.obs.columns
    )


def _create_airr_dataframe(
    obs: pd.DataFrame,
) -> pd.DataFrame:
    """
    Convert alpha/beta TCR sequences into a minimal
    AIRR-compatible DataFrame.
    """

    records = []

    for cell_id, row in obs.iterrows():

        alpha = row.get(
            "alpha"
        )

        beta = row.get(
            "beta"
        )

        if pd.notna(alpha):

            records.append(
                {
                    "cell_id": cell_id,
                    "locus": "TRA",
                    "junction_aa": alpha,
                }
            )

        if pd.notna(beta):

            records.append(
                {
                    "cell_id": cell_id,
                    "locus": "TRB",
                    "junction_aa": beta,
                }
            )

    return pd.DataFrame(
        records
    )
