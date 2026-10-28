"""Reproduce figures and tables from compact, study-derived input tables.

Marker discovery is a separate expression-matrix step in compute_markers.py.
The marker panel here uses its exported summary, without rerunning the test.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TIME_ORDER = ["cut0d", "cut6h", "cut12h", "cut24h", "cut2d", "cut3d", "cut5d", "cut7d"]
SPATIAL_ORDER = ["cut0d", "cut6h", "cut12h", "cut24h", "cut3d", "cut7d"]


def load_data(data_dir):
    """Load the compact inputs included in the repository."""
    data_dir = Path(data_dir)
    names = {
        "cells": "cell_metadata.csv.gz",
        "neoblasts": "neoblast_metadata.csv.gz",
        "spots": "visium_spots.csv.gz",
        "markers": "marker_top_genes.csv",
        "marker_panel": "marker_expression_summary.csv",
    }
    return {key: pd.read_csv(data_dir / name) for key, name in names.items()}


def simplify_domain(value):
    value = str(value)
    if value.startswith("Gut"):
        return "Gut"
    if value.startswith("Body wall"):
        return "Body wall"
    if "epidermis" in value.lower():
        return "Epidermis"
    return "Secretory" if value == "Secertory" else value


def summarize_cells(cells, neoblasts):
    counts = pd.crosstab(cells.TimePoint, cells.BigCellType).reindex(TIME_ORDER)
    overall = pd.DataFrame({"total_cells": counts.sum(axis=1), "neoblast_cells": counts["Neoblast"]})
    overall["neoblast_percent"] = 100 * overall.neoblast_cells / overall.total_cells
    subtype_counts = pd.crosstab(neoblasts.TimePoint, neoblasts.Annotation).reindex(TIME_ORDER)
    subtype_percent = 100 * subtype_counts.div(subtype_counts.sum(axis=1), axis=0)
    fold_change = subtype_percent.div(subtype_percent.loc["cut0d"], axis=1)
    change = pd.DataFrame({
        "cut0d_percent": subtype_percent.loc["cut0d"],
        "min_percent": subtype_percent.min(),
        "max_percent": subtype_percent.max(),
        "time_of_max": subtype_percent.idxmax(),
        "max_abs_log2_fold_change": np.log2(fold_change).abs().max(),
    }).sort_values("max_abs_log2_fold_change", ascending=False)
    return overall, subtype_counts, subtype_percent, fold_change, change


def summarize_spots(spots):
    spots = spots.copy()
    abundance_cols = [c for c in spots if c.startswith("cell2location_")]
    denominator = spots[abundance_cols].sum(axis=1)
    spots["neoblast_proportion"] = spots.cell2location_Neoblast.div(denominator.where(denominator > 0))
    spots["domain"] = spots.STAGATE.map(simplify_domain)
    summary = spots.groupby(["TimePoint", "domain"], observed=True).neoblast_proportion.agg(
        mean_neoblast_proportion="mean", n_spots="count"
    ).reset_index()
    summary["mean_neoblast_percent"] = 100 * summary.mean_neoblast_proportion
    domain_percent = summary.pivot(index="TimePoint", columns="domain", values="mean_neoblast_percent")
    domain_percent = domain_percent.reindex(SPATIAL_ORDER).sort_index(axis=1)
    return spots, summary, domain_percent


def _finish(fig, path):
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_overall(overall, path):
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(overall.index, overall.neoblast_percent, marker="o", lw=2, color="#1f5c4d")
    ax.set(xlabel="Time after amputation", ylabel="Neoblasts (% of sampled cells)",
           title="Early increase in sampled neoblast proportion")
    ax.tick_params(axis="x", rotation=45)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    _finish(fig, path)


def plot_subtypes(percent, change, path):
    fig, axes = plt.subplots(1, 2, figsize=(17, 5))
    percent.plot.bar(stacked=True, colormap="tab20", width=.8, ax=axes[0])
    axes[0].set(title="Neoblast subtype composition", xlabel="Time after amputation", ylabel="Percent of neoblasts")
    axes[0].legend(title="Subtype", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    for subtype in change.head(5).index:
        axes[1].plot(percent.index, percent[subtype], marker="o", label=subtype)
    axes[1].set(title="Subtypes with largest relative changes", xlabel="Time after amputation", ylabel="Percent of neoblasts")
    axes[1].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    for ax in axes:
        ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    _finish(fig, path)


def plot_markers(panel, markers, path):
    """Plot means and fractions above zero in the supplied processed matrix.

    The original neoblast matrix contains negative values. A positive value
    is not a raw transcript-detection event, so the legend says above zero.
    """
    genes = markers.names.drop_duplicates().tolist()
    groups = sorted(panel.Annotation.unique())
    mean = panel.pivot(index="Annotation", columns="gene", values="mean_processed_expression").reindex(index=groups, columns=genes)
    frac = panel.pivot(index="Annotation", columns="gene", values="fraction_above_zero").reindex(index=groups, columns=genes)
    relative = mean.subtract(mean.min(), axis=1).div((mean.max() - mean.min()).replace(0, 1), axis=1)
    xx, yy = np.meshgrid(np.arange(len(genes)), np.arange(len(groups)))
    fig, ax = plt.subplots(figsize=(14, 6))
    points = ax.scatter(xx.ravel(), yy.ravel(), s=frac.to_numpy().ravel() * 160,
                        c=relative.to_numpy().ravel(), cmap="Reds", vmin=0, vmax=1,
                        edgecolors="#777777", linewidths=.25)
    ax.set_xticks(range(len(genes)), genes, rotation=90, fontsize=8)
    ax.set_yticks(range(len(groups)), groups, fontsize=9)
    ax.invert_yaxis()
    ax.set_title("Subtype-enriched markers in the supplied processed matrix")
    fig.colorbar(points, ax=ax, label="Relative mean processed expression", shrink=.65)
    for fraction in [.2, .6, 1.0]:
        ax.scatter([], [], s=fraction*160, color="gray", label=f"{fraction:.0%}")
    ax.legend(title="Cells above zero\nin supplied matrix", bbox_to_anchor=(1.22, 1), loc="upper left")
    fig.tight_layout()
    _finish(fig, path)


def plot_domains(domain_percent, path):
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))
    cmap = plt.cm.YlGnBu.copy()
    cmap.set_bad("#eeeeee")
    im = axes[0].imshow(np.ma.masked_invalid(domain_percent.to_numpy(float)), aspect="auto", cmap=cmap)
    axes[0].set_xticks(range(len(domain_percent.columns)), domain_percent.columns, rotation=45, ha="right")
    axes[0].set_yticks(range(len(domain_percent)), domain_percent.index)
    axes[0].set(title="Mean estimated neoblast fraction by tissue domain", xlabel="Provided STAGATE domain", ylabel="Time after amputation")
    fig.colorbar(im, ax=axes[0], label="Mean estimated fraction (%)")
    for domain in ["Anterior body", "Body wall", "Epidermis", "Secretory", "Wounds"]:
        if domain in domain_percent:
            axes[1].plot(domain_percent.index, domain_percent[domain], marker="o", label=domain)
    axes[1].set(title="Selected tissue domains", xlabel="Time after amputation", ylabel="Mean estimated neoblast fraction (%)")
    axes[1].tick_params(axis="x", rotation=45)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    _finish(fig, path)


def plot_spatial(spots, path):
    chosen = []
    for time in ["cut0d", "cut24h", "cut7d"]:
        at_time = spots.loc[spots.TimePoint == time]
        library = at_time.library_id.value_counts().idxmax()
        chosen.append((time, library, at_time.loc[at_time.library_id == library]))
    vmax = np.nanpercentile(np.concatenate([df.neoblast_proportion for _, _, df in chosen]), 98)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for ax, (time, library, df) in zip(axes, chosen):
        im = ax.scatter(df.spatial_x, df.spatial_y, c=df.neoblast_proportion,
                        s=5, cmap="viridis", vmin=0, vmax=vmax)
        ax.set(title=f"{time}\nLibrary: {library}", aspect="equal", xticks=[], yticks=[])
        ax.invert_yaxis()
    fig.colorbar(im, ax=axes, label="Estimated neoblast fraction", shrink=.8)
    _finish(fig, path)
    return pd.DataFrame([{"TimePoint": t, "library_id": lib, "n_spots": len(df)} for t, lib, df in chosen])


def run_analysis(data_dir, output_dir):
    """Regenerate all panels and summary tables from the compact inputs."""
    data = load_data(data_dir)
    output_dir = Path(output_dir)
    figures, tables = output_dir / "figures", output_dir / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    overall, counts, percent, fc, change = summarize_cells(data["cells"], data["neoblasts"])
    spots, domains, domain_percent = summarize_spots(data["spots"])
    plot_overall(overall, figures / "01_overall_neoblast_proportion.png")
    plot_subtypes(percent, change, figures / "02_neoblast_subtype_dynamics.png")
    plot_markers(data["marker_panel"], data["markers"], figures / "03_marker_expression.png")
    plot_domains(domain_percent, figures / "04_spatial_domain_summary.png")
    libraries = plot_spatial(spots, figures / "05_spatial_maps.png")
    outputs = {
        "overall_neoblast_proportion": overall, "subtype_counts": counts,
        "subtype_percent": percent, "subtype_fold_change": fc, "subtype_change_summary": change,
        "spatial_domain_percent": domain_percent,
    }
    for name, frame in outputs.items():
        frame.to_csv(tables / f"{name}.csv")
    domains.to_csv(tables / "spatial_domain_summary.csv", index=False)
    libraries.to_csv(tables / "representative_libraries.csv", index=False)
    data["markers"].to_csv(tables / "marker_top_genes.csv", index=False)
    return {**outputs, "markers": data["markers"], "spatial_domain_summary": domains, "libraries": libraries}
