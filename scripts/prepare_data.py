"""Extract compact, non-expression inputs from the study's processed files."""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import pandas as pd

from neoblast_analysis.io import raw_paths, read_obs

VISIUM = {
    "cut0d": "Visium_0hpa.h5ad", "cut6h": "Visium_6hpa.h5ad",
    "cut12h": "Visium_12hpa.h5ad", "cut24h": "Visium_1dpa.h5ad",
    "cut3d": "Visium_3dpa.h5ad", "cut7d": "Visium_7dpa.h5ad",
}


def prepare(source_root, output_dir):
    source_root, output_dir = Path(source_root), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    singlecell, spatial, _ = raw_paths(source_root)
    records = []
    frames = {}
    for filename, columns, dest in [
        ("adata_scRNA_Annotated.h5ad", ["TimePoint", "BigCellType", "batch"], "cell_metadata.csv.gz"),
        ("adata_Neoblast.h5ad", ["TimePoint", "Annotation", "batch"], "neoblast_metadata.csv.gz"),
    ]:
        path = singlecell / filename
        with h5py.File(path, "r") as handle:
            frames[dest] = read_obs(handle, columns)
        records.append({"source_file": filename, "source_bytes": path.stat().st_size,
                        "n_observations": len(frames[dest]), "retained_columns": columns})
    spatial_frames = []
    for time, filename in VISIUM.items():
        path = spatial / filename
        with h5py.File(path, "r") as handle:
            columns = ["library_id", "Section_id", "STAGATE"] + sorted(
                c for c in handle["obs"] if c.startswith("cell2location_")
            )
            frame = read_obs(handle, columns)
            coords = handle["obsm/spatial"][:]
            frame["spatial_x"], frame["spatial_y"] = coords[:, 0], coords[:, 1]
            frame.insert(0, "TimePoint", time)
        spatial_frames.append(frame)
        records.append({"source_file": filename, "source_bytes": path.stat().st_size,
                        "n_observations": len(frame), "retained_columns": list(frame.columns)})
        print(f"Extracted {filename}: {len(frame):,} spots", flush=True)
    frames["visium_spots.csv.gz"] = pd.concat(spatial_frames, ignore_index=True)
    compact = []
    for filename, frame in frames.items():
        path = output_dir / filename
        frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})
        compact.append({"file": filename, "rows": len(frame), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    manifest = {
        "study_doi": "10.1038/s41467-023-39016-0", "accession": "OMIX003867",
        "source_files": records, "compact_files": compact,
        "note": "These tables retain supplied cell labels, tissue labels, spatial coordinates, and cell2location estimates. They contain no gene-expression matrix.",
    }
    (output_dir / "provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=Path("data/raw"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed"))
    args = parser.parse_args()
    prepare(args.source_root, args.output_dir)
