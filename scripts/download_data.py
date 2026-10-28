"""Download processed expression files from the study's official OMIX record.

By default download the single-cell archive. --include-spatial additionally
downloads the six Visium objects needed to rebuild the compact spot tables.
"""
import argparse
import shutil
import urllib.request
from pathlib import Path
from zipfile import ZipFile

BASE_URL = "https://download.cncb.ac.cn/OMIX/OMIX003867"
SINGLECELL = {"adata_scRNA_Annotated.h5ad", "adata_Neoblast.h5ad"}
VISIUM = {f"Visium_{time}.h5ad" for time in ["0hpa", "6hpa", "12hpa", "1dpa", "3dpa", "7dpa"]}


def fetch(archive_name, destination, wanted):
    destination.mkdir(parents=True, exist_ok=True)
    archive = destination.parent / archive_name
    if not archive.exists():
        partial = archive.with_suffix(".zip.part")
        url = f"{BASE_URL}/{archive_name}"
        print(f"Downloading {url}", flush=True)
        with urllib.request.urlopen(url, timeout=120) as response, partial.open("wb") as target:
            expected_bytes = response.headers.get("Content-Length")
            shutil.copyfileobj(response, target, length=1024 * 1024)
        if expected_bytes and partial.stat().st_size != int(expected_bytes):
            raise IOError(f"Incomplete download: {partial}. Retrieve the complete archive from the official OMIX page and retry.")
        # Check the central directory before treating an archive as complete.
        with ZipFile(partial) as bundle:
            bundle.infolist()
        partial.replace(archive)
    with ZipFile(archive) as bundle:
        members = {Path(item.filename).name: item for item in bundle.infolist()
                   if not item.is_dir() and not item.filename.startswith("__MACOSX/")}
        missing = wanted - members.keys()
        if missing:
            raise ValueError(f"Missing expected files in {archive.name}: {sorted(missing)}")
        for name in sorted(wanted):
            target = destination / name
            if target.exists():
                if target.stat().st_size != members[name].file_size:
                    raise FileExistsError(f"Existing file has a different size: {target}")
                print(f"Already present: {target}", flush=True)
                continue
            partial = target.with_suffix(".h5ad.part")
            with bundle.open(members[name]) as source, partial.open("wb") as output:
                shutil.copyfileobj(source, output, length=1024 * 1024)
            partial.replace(target)
            print(f"Extracted {target}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw"))
    parser.add_argument("--include-spatial", action="store_true")
    args = parser.parse_args()
    fetch("OMIX003867-01.zip", args.raw_dir / "singlecell_h5ad", SINGLECELL)
    if args.include_spatial:
        fetch("OMIX003867-03.zip", args.raw_dir / "Scanpy_Visium", VISIUM)
