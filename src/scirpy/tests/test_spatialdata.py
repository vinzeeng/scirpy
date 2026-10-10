import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from anndata import AnnData
from mudata import MuData

import scirpy as ir
from scirpy.util import DataHandler
from scirpy.util._spatialdata import prepare_spatialdata_for_scirpy, spatialdata_to_anndata, spatialdata_to_mudata


@pytest.mark.parametrize("container", [AnnData, MuData])
def test_data_handler_preserves_existing_input(adata_tra, container):
    adata_tra = DataHandler.default(adata_tra).adata
    data = adata_tra if container is AnnData else MuData({"airr": adata_tra})
    handler = DataHandler.default(data)
    assert handler.data is data
    assert handler.adata is adata_tra
    assert handler.airr is adata_tra.obsm["airr"]
    assert handler.chain_indices is adata_tra.obsm["chain_indices"]
    with pytest.raises(AttributeError, match="without SpatialData"):
        _ = handler.spatialdata


@pytest.fixture
def spatial_input():
    spatialdata = pytest.importorskip("spatialdata")
    from spatialdata.models import PointsModel

    table = AnnData(
        X=np.ones((3, 2)),
        obs=pd.DataFrame(
            {"alpha": ["CAVRDSNYQLIW", None, None], "beta": ["CASSLGQETQYF", "CASSIRSSYEQYF", None]},
            index=["c1", "c2", "c3"],
        ),
    )
    points = PointsModel.parse(
        pd.DataFrame(
            {"NAME": ["c2", "c3", "c1"], "X": [2.0, 5.0, 1.0], "Y": [4.0, 6.0, 3.0], "cell_type": ["B", "B", "T"]}
        ),
        coordinates={"x": "X", "y": "Y"},
    )
    return spatialdata.SpatialData(tables={"table": table}, points={"spatial": points})


def test_data_handler_spatialdata_identity_and_repeated_calls(spatial_input, monkeypatch):
    from scirpy.util import _spatialdata

    table = spatial_input.tables["table"]
    original_points = spatial_input.points["spatial"].compute()
    original_x = table.X
    calls = []
    read_airr = ir.io.read_airr
    coordinates = _spatialdata._get_spatial_coordinates

    def counted_read_airr(*args, **kwargs):
        calls.append("airr")
        return read_airr(*args, **kwargs)

    def counted_coordinates(*args, **kwargs):
        calls.append("spatial")
        return coordinates(*args, **kwargs)

    monkeypatch.setattr(ir.io, "read_airr", counted_read_airr)
    monkeypatch.setattr(_spatialdata, "_get_spatial_coordinates", counted_coordinates)
    handler = DataHandler.default(spatial_input)
    assert handler.data is handler.adata is table
    assert handler.spatialdata is spatial_input
    assert table.X is original_x
    assert table.n_vars == 2
    npt.assert_equal(table.obsm["spatial"], [[1, 3], [2, 4], [5, 6]])
    assert list(table.obs["cell_type"]) == ["T", "B", "B"]
    wrapped = DataHandler(handler)
    assert wrapped.data is table
    assert wrapped.spatialdata is spatial_input
    airr = table.obsm["airr"]
    spatial = table.obsm["spatial"]
    ir.tl.chain_qc(spatial_input)
    ir.pp.ir_dist(spatial_input, sequence="aa", metric="identity")
    ir.tl.define_clonotype_clusters(spatial_input, sequence="aa", metric="identity", receptor_arms="any")
    assert "receptor_type" in table.obs
    assert "cc_aa_identity" in table.obs
    assert "ir_dist_aa_identity" in table.uns
    assert DataHandler.default(spatial_input).data is table
    assert table.obsm["airr"] is airr
    assert table.obsm["spatial"] is spatial
    assert calls == ["spatial", "airr"]
    pdt.assert_frame_equal(spatial_input.points["spatial"].compute(), original_points)
    with pytest.raises(AttributeError):
        _ = handler.mdata


def test_data_handler_spatialdata_airr(spatial_input):
    handler = DataHandler.default(spatial_input)
    assert ak.to_list(handler.airr.locus) == [["TRA", "TRB"], ["TRB"], []]
    assert list(handler.adata.obs_names) == ["c1", "c2", "c3"]


def test_data_handler_spatialdata_chain_indices(spatial_input):
    table = spatial_input.tables["table"]
    ir.pp.index_chains(spatial_input, filter=["require_junction_aa"])
    assert "chain_indices" in table.obsm
    assert "chain_indices" in table.uns
    handler = DataHandler.default(spatial_input)
    assert ak.to_list(handler.chain_indices.VJ) == [[0, None], [None, None], [None, None]]
    assert ak.to_list(handler.chain_indices.VDJ) == [[1, None], [0, None], [None, None]]
    assert list(ir.get.airr(spatial_input, "junction_aa", "VDJ_1")[:2]) == ["CASSLGQETQYF", "CASSIRSSYEQYF"]


def test_data_handler_spatialdata_existing_arrays(spatial_input, adata_tra):
    table = spatial_input.tables["table"]
    adata_tra = DataHandler.default(adata_tra).adata
    table.obsm["airr"] = adata_tra.obsm["airr"][:3]
    table.obsm["chain_indices"] = adata_tra.obsm["chain_indices"][:3]
    table.obsm["spatial"] = np.zeros((3, 2))
    airr, indices, spatial = (table.obsm[key] for key in ("airr", "chain_indices", "spatial"))
    handler = DataHandler.default(spatial_input)
    assert handler.airr is airr
    assert handler.chain_indices is indices
    assert handler.adata.obsm["spatial"] is spatial


@pytest.mark.parametrize("geometry_type", ["circles", "polygons"])
def test_spatialdata_shapes(geometry_type):
    spatialdata = pytest.importorskip("spatialdata")
    import geopandas as gpd
    from shapely.geometry import Point, box
    from spatialdata.models import ShapesModel, TableModel

    geometry = [Point(1, 3), Point(2, 4)] if geometry_type == "circles" else [box(0, 2, 2, 4), box(1, 3, 3, 5)]
    shapes = gpd.GeoDataFrame({"geometry": geometry}, index=[10, 20])
    if geometry_type == "circles":
        shapes["radius"] = 1.0
    shapes = ShapesModel.parse(shapes)
    table = AnnData(obs=pd.DataFrame({"region": ["spots", "spots"], "instance": [20, 10]}, index=["c2", "c1"]))
    table = TableModel.parse(table, region="spots", region_key="region", instance_key="instance")
    sdata = spatialdata.SpatialData(tables={"table": table}, shapes={"spots": shapes})
    assert prepare_spatialdata_for_scirpy(sdata) is sdata.tables["table"]
    npt.assert_equal(table.obsm["spatial"], [[2, 4], [1, 3]])


def test_data_handler_spatialdata_ambiguous_tables(spatial_input):
    spatial_input.tables["second"] = spatial_input.tables["table"].copy()
    with pytest.raises(ValueError, match="multiple tables.*table_key"):
        DataHandler.default(spatial_input)
    assert prepare_spatialdata_for_scirpy(spatial_input, "second") is spatial_input.tables["second"]
    with pytest.raises(KeyError, match="not found"):
        prepare_spatialdata_for_scirpy(spatial_input, "missing")


def test_spatialdata_no_sequences(spatial_input):
    table = spatial_input.tables["table"]
    table.obs[["alpha", "beta"]] = None
    ir.pp.index_chains(spatial_input)
    assert ak.to_list(table.obsm["airr"]) == [[], [], []]
    ir.tl.chain_qc(spatial_input)
    assert list(table.obs["receptor_type"]) == ["no IR"] * 3


@pytest.mark.parametrize("plot", [ir.pl.plot_spatial, ir.pl.plot_spatial_all])
@pytest.mark.parametrize("container", ["spatialdata", "anndata", "mudata"])
def test_spatial_plot(spatial_input, monkeypatch, plot, container):
    ir.tl.chain_qc(spatial_input)
    table = spatial_input.tables["table"]
    data = spatial_input if container == "spatialdata" else table
    if container == "mudata":
        data = MuData({"gex": table, "airr": table.copy()})
    monkeypatch.setattr(plt, "show", lambda: None)
    plt.figure()
    try:
        plot(data, color="receptor_type", filter_by="cell_type", filter_value="T")
        assert len(plt.gca().collections) > 0
    finally:
        plt.close()


@pytest.mark.parametrize("convert", [spatialdata_to_anndata, spatialdata_to_mudata])
def test_explicit_spatialdata_converters_are_independent(spatial_input, convert):
    table = spatial_input.tables["table"]
    result = convert(spatial_input)
    assert result is not table
    if isinstance(result, MuData):
        assert result.mod["gex"] is not table
        assert "airr" in result.mod["airr"].obsm
    else:
        assert "airr" in result.obsm
    assert "airr" not in table.obsm
    assert "spatial" not in table.obsm
