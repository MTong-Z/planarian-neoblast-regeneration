"""Regenerate the committed result tables and figures from compact inputs."""
import argparse
from pathlib import Path

from neoblast_analysis import run_analysis

if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=root / "data/processed")
    parser.add_argument("--output-dir", type=Path, default=root / "results")
    args = parser.parse_args()
    result = run_analysis(args.data_dir, args.output_dir)
    print(result["overall_neoblast_proportion"].round(3).to_string())
    print(f"Saved figures and tables to {args.output_dir}")
