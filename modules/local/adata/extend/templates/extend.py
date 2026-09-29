#!/usr/bin/env python3
"""
Extend an AnnData object with pickled obs, uns, etc. content produced by other
modules; one directory per slot, with the file name as the key.

Every slot index must match the base object exactly, and existing columns or
keys are not replaced. `allow_missing` re-indexes the slot index to the base
adata object instead, and `overwrite` allows existing columns and keys to be
replaced.
"""

# Disable OpenMP CPU topology detection for macOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import logging
import pickle
import platform
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from threadpoolctl import threadpool_limits

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Limit BLAS/OpenMP threads to the allocated CPUs
threadpool_limits(int("${task.cpus}"))


def load_pickle(path):
    """Load a pickle file."""
    if path.suffix == ".pkl":
        with open(path, "rb") as f:
            loaded = pickle.load(f)
        return loaded
    else:
        raise ValueError(f"Unsupported file extension: `{path}`")


def extend_adata(adata, allow_missing, overwrite):
    """Extend an adata object with pickled obs, uns, etc."""

    obs_paths = sorted(Path("obs/").glob("*"))
    uns_paths = sorted(Path("uns/").glob("*"))

    for path in obs_paths:
        name = path.stem
        df = load_pickle(path)

        # Check that the indices are identical
        if not df.index.equals(adata.obs_names):
            if allow_missing:
                df = df.reindex(adata.obs_names)
                logger.info(f"Re-indexed {name} with `adata.obs_names` index")
            else:
                raise ValueError(f"Index for obs.{name} differs from adata.obs")

        # Check if columns already exists in the adata object
        for col in df.columns:
            if col in adata.obs.columns:
                if overwrite:
                    del adata.obs[col]
                    logger.warning(f"Overwrite existing `adata.obs[{col}]`")
                else:
                    raise ValueError(f"Column `{col}` already exists")

        adata.obs = pd.concat([adata.obs, df], axis=1)
        logger.info(f"Extended `adata.obs` with {name} data")

    for path in uns_paths:
        name = path.stem
        loaded = load_pickle(path)

        # Check if the named data is already present in the adata object
        if name in adata.uns:
            if overwrite:
                logger.warning(f"Overwrite existing `adata.uns[{name}]`")
            else:
                raise ValueError(f"data for `uns.{name}` already exists")

        adata.uns[name] = loaded
        logger.info(f"Extended `adata.uns` with {name} data")

    return adata


def write_versions(process_name):
    """Write software versions to a YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "anndata": importlib.metadata.version("anndata"),
            "pandas": pd.__version__,
            "numpy": np.__version__,
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Extend an adata object."""

    # Template variables
    base = "${base}"
    allow_missing = "${allow_missing}" == "true"
    overwrite = "${overwrite}" == "true"
    prefix = "${prefix}"
    process_name = "${task.process}"

    adata_base = ad.read_h5ad(base)
    adata = extend_adata(adata_base, allow_missing, overwrite)

    adata.write_h5ad(f"{prefix}.h5ad")

    write_versions(process_name)

if __name__ == "__main__":
    main()
