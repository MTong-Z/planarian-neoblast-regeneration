"""Recompute original subtype-vs-rest Wilcoxon marker rankings from expression data."""
import argparse
import hashlib
import json
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from neoblast_analysis.io import read_neoblast_matrix

KEY_SUBTYPES = ["tspan-1+ neoblast", "neuron neoblast", "protonephridia neoblast",
                "epidermal neoblast", "muscle/parenchymal neoblast"]


def compute(expression_path, annotation_path, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = read_neoblast_matrix(expression_path)
    matrix_values = adata.X.data if hasattr(adata.X, "toarray") else adata.X
    if not np.isfinite(matrix_values).all():
        raise ValueError("The input expression matrix has non-finite values. Retrieve an intact processed file before marker analysis.")
    # Suppress only floating-point overflow/invalid warnings from Scanpy's
    # log-fold-change calculation. On this processed matrix those LFCs are
    # not interpretable; retain Wilcoxon scores and BH-adjusted p-values.
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        sc.tl.rank_genes_groups(adata, groupby="Annotation", method="wilcoxon",
                               reference="rest", corr_method="benjamini-hochberg", use_raw=False)
    all_markers = pd.concat([
        sc.get.rank_genes_groups_df(adata, group=group).assign(subtype=group)
        for group in sorted(adata.obs.Annotation.unique())
    ], ignore_index=True)[["subtype", "names", "scores", "pvals", "pvals_adj"]]
    all_markers.to_csv(output_dir / "marker_rankings.csv.gz", index=False,
                       compression={"method": "gzip", "mtime": 0})
    top = pd.concat([all_markers.loc[all_markers.subtype == s].head(5) for s in KEY_SUBTYPES], ignore_index=True)
    # The course supplied this identifier-to-human-symbol lookup. Retain only
    # rows used here, rather than claiming to have rerun BLAST or validated orthology.
    annotations = pd.read_csv(annotation_path)
    top = top.merge(annotations[["smed_id", "gene_symbol"]], left_on="names", right_on="smed_id", how="left")
    top["gene_symbol"] = top.gene_symbol.fillna("no BLAST symbol")
    top.to_csv(output_dir / "marker_top_genes.csv", index=False)
    top[["names", "gene_symbol"]].drop_duplicates().rename(columns={"names": "smed_id"}).to_csv(
        output_dir / "marker_gene_annotations.csv", index=False)
    genes = top.names.drop_duplicates().tolist()
    values = adata[:, genes].X
    if hasattr(values, "toarray"):
        values = values.toarray()
    panels = []
    for group in sorted(adata.obs.Annotation.unique()):
        block = np.asarray(values[adata.obs.Annotation.to_numpy() == group])
        panels.append(pd.DataFrame({"Annotation": group, "gene": genes,
                                    "mean_processed_expression": block.mean(axis=0),
                                    "fraction_above_zero": (block > 0).mean(axis=0)}))
    pd.concat(panels, ignore_index=True).to_csv(output_dir / "marker_expression_summary.csv", index=False)
    note = {
        "matrix_file": Path(expression_path).name, "n_cells": adata.n_obs, "n_genes": adata.n_vars,
        "method": "Wilcoxon rank-sum, subtype versus all other neoblasts",
        "multiple_testing": "Benjamini-Hochberg within each subtype",
        "matrix_min": float(np.min(matrix_values)), "matrix_max": float(np.max(matrix_values)),
        "has_negative_values": bool(np.any(matrix_values < 0)),
        "expression_file_sha256": hashlib.sha256(Path(expression_path).read_bytes()).hexdigest(),
        "annotation_file_sha256": hashlib.sha256(Path(annotation_path).read_bytes()).hexdigest(),
        "scanpy_version": version("scanpy"),
        "note": "Uses all 3,000 genes retained in the provided processed neoblast matrix. Above-zero fractions are not transcript-detection rates; LFCs are excluded. This analysis pools time points and does not model biological replicates.",
    }
    (output_dir / "marker_provenance.json").write_text(json.dumps(note, indent=2, allow_nan=False) + "\n")
    print(f"Ranked {adata.n_vars:,} genes in {adata.n_obs:,} neoblasts; saved results to {output_dir}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression-file", type=Path, default=Path("data/raw/singlecell_h5ad/adata_Neoblast_recovered.h5ad"))
    parser.add_argument("--annotation-file", type=Path, default=Path("data/processed/marker_gene_annotations.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed"))
    args = parser.parse_args()
    compute(args.expression_file, args.annotation_file, args.output_dir)
