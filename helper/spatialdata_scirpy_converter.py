import pandas as pd
import scirpy as ir

from anndata import AnnData
from spatialdata import SpatialData

import pandas as pd
import scirpy as ir

from mudata import MuData
from spatialdata import SpatialData


def spatialdata_to_mudata(
    sdata: SpatialData,
) -> MuData:
    """
    Convert Slide-tags SpatialData into a
    Scirpy-compatible MuData object.

    Result
    ------
    MuData
    ├── mod["gex"]
    │   └── obsm["spatial"]
    └── mod["airr"]
        └── obsm["airr"]
    """

    # --------------------------------------------------
    # 1. GEX
    # --------------------------------------------------

    gex = sdata.tables["gex"].copy()

    # --------------------------------------------------
    # 2. Spatial coordinates
    # --------------------------------------------------

    nuclei = (
        sdata.points["nuclei"]
        .compute()
        .set_index("cell_id")
    )

    spatial = nuclei[
        ["x", "y"]
    ].reindex(
        gex.obs_names
    )

    gex.obsm["spatial"] = (
        spatial.to_numpy()
    )

    # Spatial cell-type annotation
    if "cell_type" in nuclei.columns:
        gex.obs["cell_type"] = (
            nuclei["cell_type"]
            .reindex(gex.obs_names)
        )

    # --------------------------------------------------
    # 3. TCR -> AIRR
    # --------------------------------------------------

    tcr = (
        sdata.tables["tcr"]
        .obs
        .copy()
    )

    airr_df = _create_airr_dataframe(
        tcr
    )

    airr = ir.io.read_airr(
        airr_df
    )

    # --------------------------------------------------
    # 4. Create MuData
    # --------------------------------------------------

    mdata = MuData(
        {
            "gex": gex,
            "airr": airr,
        }
    )

    return mdata


def _create_airr_dataframe(
    tcr: pd.DataFrame,
) -> pd.DataFrame:
    """
    Convert the Slide-tags alpha/beta columns
    into a minimal AIRR-compatible table.
    """

    records = []

    for cell_id, row in tcr.iterrows():

        alpha = row.get("alpha")
        beta = row.get("beta")

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

def spatialdata_to_anndata(
    sdata: SpatialData,
) -> AnnData:
    """
    Convert the Slide-tags SpatialData object into
    a Scirpy-compatible AnnData object.

    The resulting AnnData contains:

        X
            Gene expression.

        obs
            Cell metadata.

        var
            Gene metadata.

        obsm["spatial"]
            Spatial x/y coordinates.

        obsm["airr"]
            AIRR-formatted TCR data.
    """

    # --------------------------------------------------
    # 1. Gene expression as base AnnData
    # --------------------------------------------------

    adata = sdata.tables["gex"].copy()

    # --------------------------------------------------
    # 2. Spatial coordinates
    # --------------------------------------------------

    nuclei_df = (
        sdata.points["nuclei"]
        .compute()
        .set_index("cell_id")
    )

    spatial_df = nuclei_df[
        ["x", "y"]
    ].reindex(
        adata.obs_names
    )

    adata.obsm["spatial"] = (
        spatial_df.to_numpy()
    )

    # Cell type from the spatial file
    if "cell_type" in nuclei_df.columns:
        adata.obs["cell_type"] = (
            nuclei_df["cell_type"]
            .reindex(adata.obs_names)
        )

    # --------------------------------------------------
    # 3. Raw TCR data
    # --------------------------------------------------

    tcr_df = (
        sdata.tables["tcr"]
        .obs
        .copy()
    )

    # --------------------------------------------------
    # 4. Convert alpha/beta sequences to AIRR format
    # --------------------------------------------------

    airr_df = _create_airr_dataframe(
        tcr_df
    )

    # Scirpy creates an AnnData containing obsm["airr"]
    adata_tcr = ir.io.read_airr(
        airr_df
    )

    # --------------------------------------------------
    # 5. Align AIRR data to the GEX cells
    # --------------------------------------------------

    barcode_to_index = {
        barcode: index
        for index, barcode
        in enumerate(adata_tcr.obs_names)
    }

    airr_index = [
        barcode_to_index.get(
            barcode,
            None,
        )
        for barcode in adata.obs_names
    ]

    adata.obsm["airr"] = (
        adata_tcr.obsm["airr"][
            airr_index
        ]
    )

    return adata


def _create_airr_dataframe(
    tcr_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Convert Slide-tags alpha/beta columns into
    a minimal AIRR-compatible DataFrame.
    """

    airr_records = []

    for cell_id, row in tcr_df.iterrows():

        alpha = row.get("alpha")
        beta = row.get("beta")

        if pd.notna(alpha):
            airr_records.append(
                {
                    "cell_id": cell_id,
                    "locus": "TRA",
                    "junction_aa": alpha,
                }
            )

        if pd.notna(beta):
            airr_records.append(
                {
                    "cell_id": cell_id,
                    "locus": "TRB",
                    "junction_aa": beta,
                }
            )

    return pd.DataFrame(
        airr_records
    )
