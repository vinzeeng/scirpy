from pathlib import Path

import pandas as pd
import scanpy as sc

from anndata import AnnData
from spatialdata import SpatialData
from spatialdata.models import PointsModel, TableModel


def read_slidetags(path: str | Path) -> SpatialData:
    """
    Read the SCP2176 Slide-tags dataset into a SpatialData object.

    Imported data:
        - spatially mapped nuclei
        - gene expression
        - metadata
        - TCR data

    ATAC data is currently not imported.

    Parameters
    ----------
    path
        Path to the SCP2176 root directory.

    Returns
    -------
    SpatialData
        SpatialData object containing:
        - points["nuclei"]
        - tables["gex"]
        - tables["tcr"]
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

    nuclei = _read_nuclei(
        spatial_path
    )

    gex = _read_gex(
        expression_directory,
        metadata_path,
    )

    tcr = _read_tcr(
        tcr_path
    )

    return SpatialData(
        points={
            "nuclei": nuclei,
        },
        tables={
            "gex": gex,
            "tcr": tcr,
        },
    )


def _read_nuclei(
    spatial_path: Path,
):
    """
    Read spatial coordinates and cell annotations.
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

    # TYPE is a technical metadata row,
    # not a biological observation.
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

    df = df.rename(
        columns={
            "NAME": "cell_id",
            "X": "x",
            "Y": "y",
        }
    )

    if df["cell_id"].duplicated().any():
        raise ValueError(
            "Spatial file contains duplicate cell IDs."
        )

    return PointsModel.parse(
        df
    )


def _read_gex(
    expression_directory: Path,
    metadata_path: Path,
) -> AnnData:
    """
    Read 10x gene-expression data and attach metadata.
    """

    gex = sc.read_10x_mtx(
        expression_directory
    )

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
        gex.obs = gex.obs.join(
            metadata,
            how="left",
        )

    return TableModel.parse(
        gex
    )


def _read_tcr(
    tcr_path: Path,
) -> AnnData:
    """
    Read Slide-tags TCR data.
    """

    df = pd.read_csv(
        tcr_path
    )

    df = _drop_unnamed_columns(
        df
    )

    if "CB" not in df.columns:
        raise ValueError(
            "TCR file does not contain the expected 'CB' column."
        )

    # Use a consistent name for cell identifiers.
    df = df.rename(
        columns={
            "CB": "cell_id"
        }
    )

    if df["cell_id"].duplicated().any():
        raise ValueError(
            "TCR file contains duplicate cell IDs."
        )

    df = df.set_index(
        "cell_id"
    )

    tcr = AnnData(
        obs=df
    )

    return TableModel.parse(
        tcr
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

    metadata = metadata.rename(
        columns={
            "NAME": "cell_id"
        }
    )

    if metadata["cell_id"].duplicated().any():
        raise ValueError(
            "Metadata contains duplicate cell IDs."
        )

    metadata = metadata.set_index(
        "cell_id"
    )

    return metadata


def _find_10x_expression_directory(
    expression_path: Path,
) -> Path:
    """
    Find the directory containing the 10x matrix files.
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
    Check whether a required file exists.
    """

    if not path.is_file():
        raise FileNotFoundError(
            f"Required file not found: {path}"
        )
