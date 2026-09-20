from anndata import AnnData
import pandas as pd
import scanpy as sc

import scirpy as ir

def create_airr_df(filename: str) -> pd.DataFrame:
    raw_df = pd.read_csv("../SCP2176/other/slidetags_multiome_tcr.csv")

    tcr_df = pd.DataFrame(columns=["cell_id", "locus", "junction_aa"])

    for index in raw_df.index:
        line = raw_df.loc[index]

        cell_id = line["CB"]
        alpha = line["alpha"]
        beta = line["beta"]

        if pd.notna(alpha):
            tcr_df.loc[len(tcr_df)] = {
                "cell_id": cell_id,
                "locus": "TRA",
                "junction_aa": alpha
            }

        if pd.notna(beta):
            tcr_df.loc[len(tcr_df)] = {
                "cell_id": cell_id,
                "locus": "TRB",
                "junction_aa": beta
            }

    return tcr_df


def get_obs_index_mapping(source_obs_names, target_obs_names):
    barcode2idx = {
        barcode: i
        for i, barcode in enumerate(source_obs_names)
    }

    return [
        barcode2idx.get(barcode, None)
        for barcode in target_obs_names
    ]

def import_melanoma_data() -> AnnData:
    # GEX laden
    adata = sc.read_10x_mtx(
        "../SCP2176/expression/6426588be2c9c436276e802d/"
    )

    # TCR laden
    tcr_df = create_airr_df(
        "../SCP2176/other/slidetags_multiome_tcr.csv"
    )
    adata_tcr = ir.io.read_airr(tcr_df)
    tcr_idx = get_obs_index_mapping(
        adata_tcr.obs_names,
        adata.obs_names
    )
    adata.obsm["airr"] = adata_tcr.obsm["airr"][tcr_idx]

    # Spatial laden
    adata_spatial = ir.io.read_spatial_data(
        "../SCP2176/cluster/HumanMelanomaMultiome_spatial.csv"
    )
    spatial_df = pd.DataFrame(
        adata_spatial.obsm["spatial"],
        index=adata_spatial.obs_names,
        columns=["X", "Y"]
    )
    spatial_df = spatial_df.reindex(adata.obs_names)
    adata.obsm["spatial"] = spatial_df.to_numpy()
    adata.obs["cell_type"] = adata_spatial.obs["cell_type"].reindex(adata.obs_names)

    return adata
