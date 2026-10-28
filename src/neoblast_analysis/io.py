"""Read selected fields from current or legacy AnnData HDF5 files."""
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def read_array(dataset):
    if h5py.check_string_dtype(dataset.dtype) is not None:
        return dataset.asstr()[:]
    return dataset[:]


def read_obs(handle, columns):
    obs = handle["obs"]
    values = {}
    for column in columns:
        field = obs[column]
        if isinstance(field, h5py.Group) and "codes" in field:
            codes, categories = field["codes"][:], read_array(field["categories"])
            values[column] = pd.Categorical.from_codes(codes, categories)
        elif "__categories" in obs and column in obs["__categories"]:
            values[column] = pd.Categorical.from_codes(field[:], read_array(obs["__categories"][column]))
        else:
            values[column] = read_array(field)
    return pd.DataFrame(values)


def raw_paths(source_root):
    """Accept the documented raw-data layout or the original course layout."""
    source_root = Path(source_root)
    if (source_root / "homeworks").exists():
        return (source_root / "homeworks/hw5/data/singlecell_h5ad",
                source_root / "homeworks/hw4/data/OMIX003867/Scanpy_Visium",
                source_root / "homeworks/hw5/data/smed_blast.csv")
    return (source_root / "singlecell_h5ad", source_root / "Scanpy_Visium", source_root / "smed_blast.csv")


def read_neoblast_matrix(path):
    """Load only the provided neoblast matrix, labels, and gene identifiers."""
    import anndata as ad
    from scipy import sparse

    with h5py.File(path, "r") as handle:
        node = handle["X"]
        if isinstance(node, h5py.Dataset):
            matrix = node[:]
        else:
            shape = tuple(node.attrs["shape"])
            kind = node.attrs.get("encoding-type", "csr_matrix")
            constructor = sparse.csc_matrix if kind == "csc_matrix" else sparse.csr_matrix
            matrix = constructor((node["data"][:], node["indices"][:], node["indptr"][:]), shape=shape)
        obs = read_obs(handle, ["Annotation", "TimePoint"])
        var = handle["var"]
        index_key = var.attrs.get("_index", "_index")
        genes = read_array(var[index_key])
    obs.index = pd.Index([f"cell_{i}" for i in range(len(obs))])
    return ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=pd.Index(genes.astype(str))))
