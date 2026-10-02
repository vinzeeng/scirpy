from pathlib import Path

import pandas as pd
import scanpy as sc

from anndata import AnnData
from spatialdata import SpatialData
from spatialdata.models import PointsModel, TableModel


def read_slidetags(
    path: str | Path,
) -> SpatialData:
    """
    Read the SCP2176 Slide-tags dataset into SpatialData.

    Result
    ------
    SpatialData
    ├── points["spatial"]
    │   ├── NAME
    │   ├── X
    │   ├── Y
    │   └── cell_type
    │
    └── tables["adata"]
        └── AnnData containing
            - gene expression
            - metadata
            - raw TCR information
    """

    path = Path(path)

    spatial_path = (
        path
        / "cluster"
        / "HumanMelanomaMultiome_spatial.csv"
    )

    metadata_path = (
        path
        / "metadata"
        / "HumanMelanomaMultiome_metadata.csv"
    )

    tcr_path = (
        path
        / "other"
        / "slidetags_multiome_tcr.csv"
    )

    expression_directory = _find_10x_expression_directory(
        path / "expression"
    )

    _validate_file(spatial_path)
    _validate_file(metadata_path)
    _validate_file(tcr_path)

    spatial = _read_spatial_points(
        spatial_path
    )

    adata = _read_table(
        expression_directory,
        metadata_path,
        tcr_path,
    )

    return SpatialData(
        points={
            "spatial": spatial,
        },
        tables={
            "adata": adata,
        },
    )


def _read_spatial_points(
    spatial_path: Path,
):
    """
    Read the spatial coordinates from the original
    Slide-tags spatial CSV.

    Original column names are preserved.
    """

    df = pd.read_csv(
        spatial_path
    )

    df = _drop_unnamed_columns(
        df
    )

    required_columns = {
        "NAME",
        "X",
        "Y",
    }

    missing_columns = (
        required_columns
        - set(df.columns)
    )

    if missing_columns:
        raise ValueError(
            f"Spatial file is missing columns: "
            f"{missing_columns}"
        )

    # The TYPE row contains column metadata and does not
    # represent a biological observation.
    df = df[
        df["NAME"] != "TYPE"
    ].copy()

    df["X"] = pd.to_numeric(
        df["X"],
        errors="raise",
    )

    df["Y"] = pd.to_numeric(
        df["Y"],
        errors="raise",
    )

    if df["NAME"].duplicated().any():
        raise ValueError(
            "Spatial file contains duplicate NAME values."
        )

    return PointsModel.parse(
        df,
        coordinates={
            "x": "X",
            "y": "Y",
        },
    )


def _read_table(
    expression_directory: Path,
    metadata_path: Path,
    tcr_path: Path,
) -> AnnData:
    """
    Read gene expression and attach metadata and
    raw TCR information using the original cell barcodes.
    """

    # --------------------------------------------------
    # Gene expression
    # --------------------------------------------------

    adata = sc.read_10x_mtx(
        expression_directory
    )

    # --------------------------------------------------
    # Metadata
    # --------------------------------------------------

    metadata = pd.read_csv(
        metadata_path
    )

    metadata = _drop_unnamed_columns(
        metadata
    )

    metadata = _prepare_metadata(
        metadata
    )

    if metadata is not None:
        adata.obs = adata.obs.join(
            metadata,
            how="left",
        )

    # --------------------------------------------------
    # TCR
    # --------------------------------------------------

    tcr = pd.read_csv(
        tcr_path
    )

    tcr = _drop_unnamed_columns(
        tcr
    )

    if "CB" not in tcr.columns:
        raise ValueError(
            "TCR file does not contain "
            "the expected 'CB' column."
        )

    if tcr["CB"].duplicated().any():
        raise ValueError(
            "TCR file contains duplicate CB values."
        )

    # CB contains the same cell barcodes used by
    # the 10x gene-expression AnnData.
    tcr = tcr.set_index(
        "CB"
    )

    adata.obs = adata.obs.join(
        tcr,
        how="left",
    )

    return TableModel.parse(
        adata
    )


def _prepare_metadata(
    metadata: pd.DataFrame,
) -> pd.DataFrame | None:
    """
    Prepare metadata for joining with the GEX AnnData.
    """

    if "NAME" not in metadata.columns:
        return None

    metadata = metadata[
        metadata["NAME"] != "TYPE"
    ].copy()

    if metadata["NAME"].duplicated().any():
        raise ValueError(
            "Metadata contains duplicate NAME values."
        )

    # NAME contains the same cell barcodes used
    # by AnnData.obs_names.
    metadata = metadata.set_index(
        "NAME"
    )

    return metadata


def _find_10x_expression_directory(
    expression_path: Path,
) -> Path:
    """
    Find the directory containing the standard
    10x gene-expression matrix files.
    """

    if not expression_path.exists():
        raise FileNotFoundError(
            f"Expression directory not found: "
            f"{expression_path}"
        )

    required_files = {
        "matrix.mtx.gz",
        "barcodes.tsv.gz",
        "features.tsv.gz",
    }

    for directory in expression_path.iterdir():

        if not directory.is_dir():
            continue

        files = {
            file.name
            for file in directory.iterdir()
        }

        if required_files.issubset(
            files
        ):
            return directory

    raise FileNotFoundError(
        "Could not find a directory containing "
        "matrix.mtx.gz, barcodes.tsv.gz and "
        "features.tsv.gz."
    )


def _drop_unnamed_columns(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Remove automatically generated CSV index columns
    such as 'Unnamed: 0'.
    """

    unnamed_columns = [
        column
        for column in df.columns
        if str(column).startswith(
            "Unnamed:"
        )
    ]

    if unnamed_columns:
        df = df.drop(
            columns=unnamed_columns
        )

    return df


def _validate_file(
    path: Path,
) -> None:
    """
    Check whether a required dataset file exists.
    """

    if not path.is_file():
        raise FileNotFoundError(
            f"Required file not found: {path}"
        )
