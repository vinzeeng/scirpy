from pathlib import Path

import scanpy as sc
import awkward as ak

from anndata import AnnData
from scirpy.io import from_airr_cells

import pandas as pd
from numba.core.cgutils import false_bit
from spatialdata import SpatialData
from spatialdata.models import PointsModel, TableModel

from scirpy.io import AirrCell


def _read_tcr (path: str| Path) -> AnnData:

    cells = []
    path = Path(path)

    df = pd.read_csv(path)

    df.set_index("CB")
    df = df.drop(columns=["Unnamed: 0"])
    #df = df.ilo

    for index, row in df.iterrows():

        cb = row["CB"]
        alpha = row["alpha"]
        beta = row["beta"]

        alpha_chain = False
        if pd.notna(alpha) and alpha != "":
            alpha_chain = True

        beta_chain = False
        if pd.notna(beta) and beta != "":
            beta_chain = True

        if alpha_chain or beta_chain:
            cell = AirrCell(cb)

            if alpha_chain:
                alpha_chain = AirrCell.empty_chain_dict()

                alpha_chain.update({
                    "locus": "TRA",
                    "junction_aa": alpha
                })

                cell.add_chain(alpha_chain)

            if beta_chain:
                beta_chain = AirrCell.empty_chain_dict()

                beta_chain.update({
                    "locus": "TRB",
                    "junction_aa": beta
                })

                cell.add_chain(beta_chain)

            cells.append(cell)

    adata = from_airr_cells(cells)
    return adata


def read_spatial_data (path: str| Path) -> SpatialData:
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

    points =  read_points(spatial_path)
    adata = sc.read_10x_mtx(expression_directory)

    metadata = _read_metadata(metadata_path)
    adata.obs = adata.obs.join(metadata, how="left")

    tcr_ann_data = _read_tcr(tcr_path)

    print(len(tcr_ann_data.obsm["airr"]))
    print(tcr_ann_data.obsm)
    print(tcr_ann_data.obs_names[:5])
    print(adata.obs_names[:5])

    # Position jeder GEX-Zelle im AIRR-AnnData ermitteln
    indices = tcr_ann_data.obs_names.get_indexer(adata.obs_names)

    # AIRR-Array um eine leere Kettenliste erweitern
    airr = tcr_ann_data.obsm["airr"]
    extended = ak.concatenate([airr, airr[:1, :0]], axis=0)

    # Fehlende Zellen auf die leere Liste verweisen lassen
    indices[indices == -1] = len(airr)

    # In das GEX-AnnData übernehmen
    adata.obsm["airr"] = extended[indices]

    print(tcr_ann_data.shape)

    # tcr = tcr.drop(columns=["Unnamed: 0"])
    #
    # print(tcr[["alpha", "beta"]].head(10).to_string())
    # print(tcr[["alpha", "beta"]].dtypes)
    #
    # adata.obs = adata.obs.join(tcr, how="left")

    table = TableModel.parse(
        adata
    )

    return SpatialData(
        points={
            "spatial": points,
        },
        tables={
            "adata": table,
        },
    )


def read_points (file_path: str
                 ) -> PointsModel:
    df = pd.read_csv(file_path)

    df = df.iloc[1:].copy()
    df = df.set_index("NAME")
    df["X"] = pd.to_numeric(df["X"])
    df["Y"] = pd.to_numeric(df["Y"])

    return PointsModel.parse(df,
                             coordinates={"x": "X", "y": "Y"}
                             )

def _find_10x_expression_directory (file_path: str
    )-> str:
    directory = Path(file_path)

    for subdir in directory.iterdir():
        if not subdir.is_dir():
            continue

        required_files = {
            "matrix.mtx.gz",
            "barcodes.tsv.gz",
            "features.tsv.gz",
        }


        files = set()

        for file in subdir.iterdir():
            files.add(file.name)

        if required_files.issubset(files):
            return str(subdir)

    raise FileNotFoundError(
        "Could not find a directory containing "
        "matrix.mtx.gz, barcodes.tsv.gz and "
        "features.tsv.gz."
    )

def _read_metadata(file_path: Path) -> pd.DataFrame:
    df =  pd.read_csv(file_path)
    df = df.set_index("NAME")
    df = df.iloc[1:].copy()

    return df

