"""Recover invalid rows of the supplied scaled matrix from published count data.

Apply the study's normalize_total(1e4), log1p, and scale(max_value=10)
recipe using the mean and std retained in the neoblast object. Before replacing
any rows, require agreement with at least 95% of supplied profiles to an absolute
tolerance of 2e-5. Reconstruct disagreeing profiles; source files are never edited.
"""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse


def recover(counts_path, neoblast_path, output_path, report_path):
    counts_path, neoblast_path, output_path = map(Path, [counts_path, neoblast_path, output_path])
    if output_path.resolve() in [counts_path.resolve(), neoblast_path.resolve()]:
        raise ValueError("Recovery output must be a different file from both inputs")
    with h5py.File(neoblast_path, "r") as neo, h5py.File(counts_path, "r") as full:
        supplied = neo["X"][:]
        bad_rows = ((~np.isfinite(supplied)) | (np.abs(supplied) > 1000)).any(axis=1)
        obs_lookup = {value: i for i, value in enumerate(full["obs/_index"].asstr()[:])}
        gene_lookup = {value: i for i, value in enumerate(full["var/_index"].asstr()[:])}
        selected_rows = np.array([obs_lookup[value] for value in neo["obs/_index"].asstr()[:]])
        selected_genes = np.array([gene_lookup[value] for value in neo["var/_index"].asstr()[:]])
        mean, std = neo["var/mean"][:], neo["var/std"][:]
        if not (np.isfinite(mean).all() and np.isfinite(std).all() and (std > 0).all()):
            raise ValueError("Invalid stored scaling parameters")
        matrix = full["X"]
        pointers = matrix["indptr"][:]
        n_genes = int(matrix.attrs["shape"][1])
        recovered = supplied.copy()
        max_intact_difference = 0.0
        repaired_rows = np.zeros(len(supplied), dtype=bool)
        for start in range(0, len(selected_rows), 256):
            stop = min(start + 256, len(selected_rows))
            rows = selected_rows[start:stop]
            first, last = int(rows.min()), int(rows.max()) + 1
            left, right = int(pointers[first]), int(pointers[last])
            raw_values = matrix["data"][left:right]
            if not (np.isfinite(raw_values).all() and (raw_values >= 0).all()):
                raise ValueError("Counts contain invalid or negative values")
            block = sparse.csr_matrix((raw_values, matrix["indices"][left:right],
                                       pointers[first:last+1] - left), shape=(last-first, n_genes))
            block = block[rows-first]
            totals = np.asarray(block.sum(axis=1)).ravel()
            if (totals <= 0).any():
                raise ValueError("Count matrix contains an empty selected cell")
            values = block[:, selected_genes].toarray().astype(np.float64)
            values *= 1e4 / totals[:, None]
            values = np.minimum((np.log1p(values) - mean) / std, 10).astype(np.float32)
            with np.errstate(invalid="ignore", over="ignore"):
                difference_by_row = np.max(np.abs(values - supplied[start:stop]), axis=1)
            intact = (~bad_rows[start:stop]) & (difference_by_row <= 2e-5)
            if intact.any():
                difference = float(np.max(np.abs(values[intact] - supplied[start:stop][intact])))
                max_intact_difference = max(max_intact_difference, difference)
            repair = ~intact
            repaired_rows[start:stop] = repair
            recovered[start:stop][repair] = values[repair]
        if repaired_rows.mean() > .05:
            raise ValueError(f"Reconstruction disagrees with more than 5% of profiles: {repaired_rows.sum()} rows. Investigate the inputs.")
        if not np.isfinite(recovered).all():
            raise ValueError("Recovered matrix is not finite")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(neoblast_path, output_path)
    with h5py.File(output_path, "r+") as handle:
        handle["X"][:] = recovered
    report = {
        "counts_source": counts_path.name, "scaled_source": neoblast_path.name,
        "counts_source_sha256": hashlib.file_digest(counts_path.open("rb"), "sha256").hexdigest(),
        "scaled_source_sha256": hashlib.file_digest(neoblast_path.open("rb"), "sha256").hexdigest(),
        "recovered_source_sha256": hashlib.file_digest(output_path.open("rb"), "sha256").hexdigest(),
        "n_cells": int(len(supplied)), "n_genes": int(supplied.shape[1]),
        "rows_with_nonfinite_or_extreme_values": int(bad_rows.sum()),
        "rows_reconstructed": int(repaired_rows.sum()), "matching_rows_validated": int((~repaired_rows).sum()),
        "max_absolute_difference_on_intact_rows": max_intact_difference,
        "agreement_tolerance": 2e-5,
        "method": "Library-size normalization to 10000; log1p; supplied per-gene mean/std; upper clipping at 10. Replace only profiles disagreeing with reconstruction by more than 2e-5 or containing non-finite/extreme values. Require agreement on at least 95% of profiles.",
        "study_recipe": "https://github.com/yangteam/Spatiotemporal-transcriptomic-atlas-of-planarian-regeneration/blob/main/Analysis/5.%20scRNA-seq%20Neoblast%20analysis.py",
    }
    report_path = Path(report_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts-file", type=Path, default=Path("data/raw/singlecell_h5ad/adata_scRNA_Annotated.h5ad"))
    parser.add_argument("--neoblast-file", type=Path, default=Path("data/raw/singlecell_h5ad/adata_Neoblast.h5ad"))
    parser.add_argument("--output-file", type=Path, default=Path("data/raw/singlecell_h5ad/adata_Neoblast_recovered.h5ad"))
    parser.add_argument("--report-file", type=Path, default=Path("data/processed/matrix_recovery.json"))
    args = parser.parse_args()
    recover(args.counts_file, args.neoblast_file, args.output_file, args.report_file)
