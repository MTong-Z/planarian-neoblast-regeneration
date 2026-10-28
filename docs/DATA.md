# Data and reproducibility

## Study and upstream provenance

This project reanalyzes processed data from Cui et al. (2023),
[doi:10.1038/s41467-023-39016-0](https://doi.org/10.1038/s41467-023-39016-0).
The official [OMIX003867 record](https://ngdc.cncb.ac.cn/omix/release/OMIX003867)
provides both the single-cell H5AD and Scanpy Visium archives. The original study's
raw reads are under [CRA007941](https://ngdc.cncb.ac.cn/gsa/browse/CRA007941).
Downloading raw reads is not necessary for this project.

| Archive | Official download | Contents used |
| --- | --- | --- |
| OMIX003867-01 | [singlecell_h5ad ZIP](https://download.cncb.ac.cn/OMIX/OMIX003867/OMIX003867-01.zip) | `adata_scRNA_Annotated.h5ad`, `adata_Neoblast.h5ad` |
| OMIX003867-03 | [Scanpy_Visium ZIP](https://download.cncb.ac.cn/OMIX/OMIX003867/OMIX003867-03.zip) | Six `Visium_*.h5ad` files |

The archive sizes listed by OMIX are approximately 444 MB and 761 MB; extracted
files require several GB. Network speed and archive structure are controlled by
the upstream repository. Downloads are optional because compact inputs are included.

The course also supplied `smed_blast.csv`, a planarian-ID-to-human-symbol lookup.
Only the lookup rows used for the displayed markers are retained in
`marker_gene_annotations.csv`. The project did not generate this alignment or
independently validate orthology. This lookup is not a claim that those genes have
identical functions in the two species.

## Compact inputs included in GitHub

| File under `data/processed/` | Unit | Columns / purpose |
| --- | --- | --- |
| `cell_metadata.csv.gz` | One of 55,014 sampled cells | `TimePoint`, `BigCellType`, `batch` |
| `neoblast_metadata.csv.gz` | One of 14,601 neoblasts | `TimePoint`, `Annotation`, `batch` |
| `visium_spots.csv.gz` | One spatial spot | `TimePoint`, `library_id`, `Section_id`, `STAGATE`, supplied `cell2location_*` estimates, `spatial_x`, `spatial_y` |
| `marker_rankings.csv.gz` | One gene/subtype comparison | Wilcoxon score, unrounded p-value and BH-adjusted p-value; 3,000 genes per subtype |
| `marker_top_genes.csv` | One selected marker/subtype pair | Five markers each for five selected subtypes, with supplied symbol lookup |
| `marker_expression_summary.csv` | One selected gene/subtype pair | Mean value and fraction above zero in the supplied processed matrix |
| `marker_gene_annotations.csv` | One selected gene | Minimal subset of the course-provided symbol lookup |
| `provenance.json` | Source/derived-file inventory | Source filenames, sizes, observation counts, retained fields, compact-input SHA-256 checksums |
| `marker_provenance.json` | Marker-analysis provenance | Input checksums, matrix dimensions/range, method, Scanpy version |
| `matrix_recovery.json` | Expression-data audit | Source checksums, reconstructed profiles, tolerance, and agreement against the remaining profiles |

These compact tables contain no full gene-expression matrix or tissue image.
Rows in cell/spot tables preserve source order; observation identifiers are omitted
because this analysis does not need them. Existing annotations and model estimates
are retained, not recomputed. `batch` and `Section_id` are included for context,
but the original descriptive analysis does not model them as biological replicates.

## Rebuild from full processed files

From the repository root, after `uv sync --frozen`:

```bash
# Download both official processed-data archives and extract the required files.
uv run python scripts/download_data.py --include-spatial

# Re-extract only the metadata, coordinates, and supplied spatial estimates.
uv run python scripts/prepare_data.py --source-root data/raw

# Recompute gene-level marker discovery from the neoblast expression matrix.
uv run python scripts/recover_neoblast_matrix.py
uv run python scripts/compute_markers.py

# Rebuild figures and validate them against the original course-project results.
uv run python scripts/run_analysis.py
uv run python scripts/validate_results.py
```

For marker recomputation alone, omit `--include-spatial` when downloading, then run
`scripts/recover_neoblast_matrix.py` and `scripts/compute_markers.py`; the included compact spot table is sufficient for the
rest of the project. `--annotation-file` accepts an alternative complete lookup.
The included lookup covers the original displayed top markers, not all 3,000 genes.

If the archives have already been downloaded, arrange the required files as:

```text
data/raw/
  singlecell_h5ad/
    adata_scRNA_Annotated.h5ad
    adata_Neoblast.h5ad
  Scanpy_Visium/
    Visium_0hpa.h5ad
    Visium_6hpa.h5ad
    Visium_12hpa.h5ad
    Visium_1dpa.h5ad
    Visium_3dpa.h5ad
    Visium_7dpa.h5ad
```

Large raw/processed expression files and ZIP archives are ignored by Git.
`scripts/prepare_data.py` also accepts the original course-directory structure via
`--source-root`, so the source files do not need to be moved or modified.

## Interpretation of the exported quantities

1. **Overall single-cell proportion:** cells labeled `Neoblast` divided by all
   sampled cells at that time point.
2. **Subtype proportion:** cells with a subtype annotation divided by all neoblasts
   at that time point. Fold change is relative to `cut0d` with this same denominator.
3. **Marker ranking:** each annotated subtype versus all other neoblasts, pooling
   time points. No time-specific differential-expression model is fit. The matrix
   already contains processed expression for 3,000 retained genes; log-fold-change
   estimates are excluded because the matrix includes negative values. The
   repository edition uses the recovered matrix described below. Numerically zero
   p-values can result from floating-point underflow; they are not exactly zero probabilities.
4. **Marker-panel size:** fraction of processed values above zero, not RNA detection
   rate. Panel color is the subtype's mean processed expression, rescaled per gene
   across subtypes to 0–1 for display. Rescaling is only a visualization step.
5. **Spatial fraction:** supplied `cell2location_Neoblast` divided by the sum of all
   available `cell2location_*` columns at that spot. This denominator covers the
   supplied seven major cell types, not a direct count of every cell in the tissue.
6. **Domain average:** unweighted mean of spot-level fractions, pooling all spots
   assigned to a broad domain at that time point. Gut, body-wall, and epidermis
   label variants are merged; the source spelling `Secertory` is corrected.
7. **Missing domains:** a missing domain/time combination remains NaN. In particular,
   `Wounds` is labeled only at day 7 in these inputs; that does not mean wounds were
   biologically absent at earlier time points.

## Validation scope

The analysis was rerun locally from the original processed files during repository
preparation. The compact workflow and the executed notebook were then checked against
the saved May 2026 notebook's cell counts, proportions, subtype peaks, and selected
spatial-domain values. Marker rankings were recomputed after the matrix audit and
correction detailed below. The download links were checked against the official
OMIX page; fresh multi-GB downloads were not required for the metadata/spatial rerun.
Checks do not establish biological
replicability or validate the upstream annotations and abundance models.

## Expression-matrix recovery

During repository preparation, the supplied scaled neoblast object had 164 profiles
with non-finite or extreme values, plus 11 additional profiles inconsistent with
the original count-based preprocessing. An initial portion streamed from the
official OMIX ZIP matched the corresponding local bytes, including the invalid
region; the issue was not treated as a valid biological signal.

`scripts/recover_neoblast_matrix.py` matches cells and genes to
`adata_scRNA_Annotated.h5ad`, normalizes counts to 10,000 per cell, applies log1p,
then uses the stored neoblast gene means and standard deviations and clips the
upper values at 10. This follows the study's
[neoblast preprocessing recipe](https://github.com/yangteam/Spatiotemporal-transcriptomic-atlas-of-planarian-regeneration/blob/main/Analysis/5.%20scRNA-seq%20Neoblast%20analysis.py).

The reconstruction agrees with **14,426 profiles** within a maximum absolute
difference of **1.90735e-6**. Only the **175 disagreeing profiles** are replaced in a
separate recovered file; the originals are never edited. The script stops if more
than 5% of profiles disagree, preventing automatic application to an unrelated or
substantially different dataset. Checksums and the audit are recorded in
`matrix_recovery.json`. Marker computation refuses non-finite matrices. The
repository's recomputed marker results should be interpreted as a documented
data-quality correction to the course analysis, rather than byte-identical reuse
of its original marker output.
