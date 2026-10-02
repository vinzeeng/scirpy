import pandas as pd
import matplotlib.pyplot as plt

from anndata import AnnData
from mudata import MuData
from matplotlib.ticker import MultipleLocator

def plot_spatial_all(
    data: AnnData | MuData,
    color: str,
    filter_by: str | None = None,
    filter_value: str | None = None,
):
    spatial, obs = _get_spatial_data(data)

    df = pd.DataFrame(
        spatial,
        index=obs.index,
        columns=["X", "Y"]
    )

    mask = df.notna().all(axis=1)
    df = df[mask].copy()

    color_values = _get_obs_column(
        data,
        color
    )

    df[color] = color_values.reindex(
        df.index
    )

    if filter_by is not None and filter_value is not None:
        filter_values = _get_obs_column(
            data,
            filter_by
        )

        df[filter_by] = filter_values.reindex(
            df.index
        )

        df = df[
            df[filter_by] == filter_value
        ]

    ax = plt.gca()

    ax.xaxis.set_major_locator(
        MultipleLocator(1000)
    )

    ax.yaxis.set_major_locator(
        MultipleLocator(1000)
    )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    ax.set_title(
        f"Spatial distribution by {color}"
    )

    ax.grid(True)

    # ---------------------------------
    # Hintergrund: alle gefilterten Zellen
    # ---------------------------------

    ax.scatter(
        df["X"],
        df["Y"],
        alpha=0.25,
        label="all cells"
    )

    # ---------------------------------
    # Vordergrund: Zellen mit Farbwert
    # ---------------------------------

    colored_df = df[
        df[color].notna()
    ]

    for value in colored_df[color].unique():

        subset = colored_df[
            colored_df[color] == value
        ]

        ax.scatter(
            subset["X"],
            subset["Y"],
            label=value
        )

    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.set_aspect("equal")

    ax.legend()

    plt.show()

def plot_spatial(
    data: AnnData | MuData,
    color: str,
    filter_by: str | None = None,
    filter_value: str | None = None,
):
    spatial, obs = _get_spatial_data(
        data
    )

    df = pd.DataFrame(
        spatial,
        index=obs.index,
        columns=["X", "Y"]
    )

    # Nur Zellen mit gültigen Spatial-Koordinaten
    mask = df.notna().all(axis=1)
    df = df[mask].copy()

    # Werte für die Einfärbung holen
    color_values = _get_obs_column(
        data,
        color
    )

    df[color] = color_values.reindex(
        df.index
    )

    # Optionaler Filter, z. B. nach Gewebetyp / Zelltyp
    if filter_by is not None and filter_value is not None:

        filter_values = _get_obs_column(
            data,
            filter_by
        )

        df[filter_by] = filter_values.reindex(
            df.index
        )

        df = df[
            df[filter_by] == filter_value
        ]

    ax = plt.gca()

    ax.xaxis.set_major_locator(
        MultipleLocator(1000)
    )

    ax.yaxis.set_major_locator(
        MultipleLocator(1000)
    )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    ax.invert_yaxis()
    ax.invert_xaxis()

    ax.set_title(
        f"Spatial distribution by {color}"
    )

    ax.grid(True)

    for value in df[color].dropna().unique():

        subset = df[
            df[color] == value
        ]

        ax.scatter(
            subset["X"],
            subset["Y"],
            label=value
        )

    ax.legend()

    plt.show()


def _get_spatial_data(
    data: AnnData | MuData,
):
    if isinstance(data, AnnData):
        return (
            data.obsm["spatial"],
            data.obs,
        )

    if isinstance(data, MuData):
        return (
            data["gex"].obsm["spatial"],
            data["gex"].obs,
        )

    raise TypeError(
        "Expected AnnData or MuData."
    )


def _get_obs_column(
    data: AnnData | MuData,
    column: str,
) -> pd.Series:

    if isinstance(data, AnnData):

        if column not in data.obs.columns:
            raise KeyError(
                f"Column '{column}' not found."
            )

        return data.obs[column]

    if isinstance(data, MuData):

        # First try global MuData obs
        if column in data.obs.columns:
            return data.obs[column]

        # Then try gex obs
        if column in data["gex"].obs.columns:
            return data["gex"].obs[column]

        # Then try airr obs
        if column in data["airr"].obs.columns:
            return data["airr"].obs[column]

        # Scirpy often prefixes AIRR columns
        airr_column = f"airr:{column}"

        if airr_column in data.obs.columns:
            return data.obs[airr_column]

        raise KeyError(
            f"Column '{column}' not found."
        )

    raise TypeError(
        "Expected AnnData or MuData."
    )
