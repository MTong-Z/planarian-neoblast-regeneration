"""Check compact-data integrity and agreement with saved course-project results."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def validate():
    data, tables = ROOT / "data/processed", ROOT / "results/tables"
    manifest = json.loads((data / "provenance.json").read_text())
    for item in manifest["compact_files"]:
        assert hashlib.sha256((data / item["file"]).read_bytes()).hexdigest() == item["sha256"], item["file"]
    overall = pd.read_csv(tables / "overall_neoblast_proportion.csv", index_col=0)
    assert overall.total_cells.sum() == 55014
    assert overall.neoblast_cells.sum() == 14601
    # Benchmarks transcribed from the original executed notebook, before refactoring.
    expected = [23.776371, 30.412371, 28.602267, 30.134887, 24.297823, 27.201721, 24.570490, 21.446154]
    np.testing.assert_allclose(overall.neoblast_percent, expected, atol=1e-6)
    subtype = pd.read_csv(tables / "subtype_percent.csv", index_col=0)
    assert subtype.shape == (8, 9)
    np.testing.assert_allclose(subtype.sum(axis=1), 100)
    assert round(subtype.loc["cut2d", "neuron neoblast"], 1) == 15.6
    assert round(subtype.loc["cut7d", "tspan-1+ neoblast"], 1) == 23.7
    domains = pd.read_csv(tables / "spatial_domain_percent.csv", index_col=0)
    assert round(domains.loc["cut0d", "Anterior body"], 2) == 18.86
    assert round(domains.loc["cut7d", "Anterior body"], 2) == 42.33
    assert round(domains.loc["cut7d", "Wounds"], 2) == 21.82
    assert np.isnan(domains.loc["cut0d", "Wounds"]), "Missing domains must not become zeros"
    assert (domains.stack() >= 0).all() and (domains.stack() <= 100).all()
    markers = pd.read_csv(tables / "marker_top_genes.csv")
    assert len(markers) == 25 and markers.subtype.nunique() == 5
    assert markers.loc[markers.subtype == "neuron neoblast", "names"].iloc[0] == "SMED30035424"
    assert markers.pvals_adj.between(0, 1).all()
    panel = pd.read_csv(data / "marker_expression_summary.csv")
    assert np.isfinite(panel[["mean_processed_expression", "fraction_above_zero"]]).all().all()
    assert panel.fraction_above_zero.between(0, 1).all()
    figures = list((ROOT / "results/figures").glob("*.png"))
    assert len(figures) == 5 and all(p.stat().st_size > 1000 for p in figures)
    print("PASS: input checksums, 55,014 cells, 14,601 neoblasts, temporal benchmarks,")
    print("      subtype composition, spatial benchmarks/missingness, markers, and five figures.")


if __name__ == "__main__":
    validate()
