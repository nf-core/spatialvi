#!/usr/bin/env python3
"""
Extend an AnnData object with pickled `obs`, `var`, `obsm`, `varm`, `obsp`,
`uns` and H5AD layer content produced by other modules; one directory per slot,
with the file name as the key, except for `obs` and `var`, which use the column
names.

Every slot index must match the input adata object exactly, and existing columns
or keys are not replaced; `overwrite` allows existing columns and keys to be
replaced. With `align`, slot data is matched to the input adata object by name
instead: entries are selected and reordered, and names without slot data become
missing values where the slot can hold them (`obs`, `var`, `obsm`, `varm`).
Graphs (`obsp`) are never aligned. Layers must contain every name of the input
adata object.
"""

# Disable OpenMP CPU topology detection for macOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

# Keep caches in the task's work directory, which is always writable and
# private to the task
os.environ["NUMBA_CACHE_DIR"] = os.path.join(os.getcwd(), ".cache", "numba")
os.environ["MPLCONFIGDIR"] = os.path.join(os.getcwd(), ".cache", "matplotlib")
os.environ["XDG_CACHE_HOME"] = os.path.join(os.getcwd(), ".cache")

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

logging.basicConfig(
    level=logging.INFO,
    format="%(name)s - %(levelname)s: %(message)s"
)
logger = logging.getLogger(__name__)

# Limit BLAS/OpenMP threads to the allocated CPUs
threadpool_limits(int("${task.cpus}"))


def check_suffix(slot, path):
    """
    Check that slot data has the file type of its slot: H5AD for `layers`,
    pickle for all other slots.
    """
    expected_suffix = ".h5ad" if slot == "layers" else ".pkl"
    if path.suffix != expected_suffix:
        raise ValueError(
            f"Slot data `{path}` must be a `{expected_suffix}` file"
        )


def load_pickle(path):
    """Load slot data from a pickle file."""
    with open(path, "rb") as f:
        loaded = pickle.load(f)
    return loaded


def load_layer(path):
    """
    Load layer data from an H5AD file.

    Only `X`, `obs_names` and `var_names` may be present; the minimal amount of
    data required for extending with a layer.
    """
    layer_adata = ad.read_h5ad(path)
    non_empty = []
    for slot in ["obs", "var", "obsm", "varm", "obsp", "varp", "uns", "layers"]:
        content = getattr(layer_adata, slot)
        # `obs` and `var` always have rows (the names), so check columns
        if slot in ("obs", "var"):
            n_entries = len(content.columns)
        else:
            n_entries = len(content)
        if n_entries > 0:
            non_empty.append(f"`{slot}`")
    if non_empty:
        raise ValueError(
            f"Layer data `{path}` has content in {', '.join(non_empty)}; only "
            "`X`, `obs_names` and `var_names` are allowed"
        )
    return layer_adata


def handle_index_mismatches(adata, data, slot, name, align):
    """
    Check that the index of slot data matches the corresponding adata axis.

    `obs` and `obsm` are aligned to `obs_names`, `var` and `varm` to
    `var_names`. A mismatch (different values or order) raises an error, unless
    `align` is set, in which case the data is re-indexed to the axis and
    unmatched entries become missing values.
    """
    idx_name = slot[0:3]
    slot_idx = getattr(adata, idx_name).index
    if not data.index.equals(slot_idx):
        if align:
            # Re-indexing should fail with duplicated names
            if not (slot_idx.is_unique and data.index.is_unique):
                raise ValueError(
                    f"Can't re-index slot data `{slot}/{name}.pkl`: names in "
                    f"`adata.{idx_name}_names` or the slot data are duplicated"
                )
            # Log number of missing entries
            n_missing = (~slot_idx.isin(data.index)).sum()
            data = data.reindex(slot_idx)
            logger.warning(
                f"Re-indexed slot data `{slot}/{name}.pkl` to "
                f"`adata.{idx_name}_names`; {n_missing} of {len(slot_idx)} "
                "entries have no slot data"
            )
        else:
            raise ValueError(
                f"Index of slot data `{slot}/{name}.pkl` differs from "
                f"`adata.{idx_name}_names`; set `align` to match it by name"
            )
    return data


def handle_name_collisions(adata, slot, name, names, overwrite):
    """
    Handle overlapping names between data and `adata.{slot}`.

    Overwrites an already existing column/key in `adata.{slot}` if `overwrite`
    is set, otherwise raises an error.
    """
    if name in names:
        if overwrite:
            del getattr(adata, slot)[name]
            logger.warning(f'Overwrote existing `adata.{slot}["{name}"]`')
        else:
            raise ValueError(f'`adata.{slot}["{name}"]` already exists')
    return adata


def extend_frame(adata, data, slot, name, align, overwrite):
    """
    Add the columns of a DataFrame to `adata.obs` or `adata.var`.
    """
    data = handle_index_mismatches(adata, data, slot, name, align)

    adata_cols = getattr(adata, slot).columns
    for data_col in data.columns:
        adata = handle_name_collisions(
            adata,
            slot,
            data_col,
            adata_cols,
            overwrite
        )

    setattr(adata, slot, pd.concat([getattr(adata, slot), data], axis=1))
    logger.info(
        f"Added columns {list(data.columns)} to `adata.{slot}` from slot data "
        f"`{slot}/{name}.pkl`"
    )
    return adata


def extend_matrix(adata, data, slot, name, align, overwrite):
    """Add a matrix to `adata.obsm` or `adata.varm` under the key `name`."""
    data = handle_index_mismatches(adata, data, slot, name, align)
    names = getattr(adata, slot).keys()
    adata = handle_name_collisions(adata, slot, name, names, overwrite)
    getattr(adata, slot)[name] = data.to_numpy()
    logger.info(f'Added `adata.{slot}["{name}"]`')
    return adata


def extend_pairwise(adata, data, name, overwrite):
    """
    Add a pairwise observation annotation to `adata.obsp` with `name` key.

    The `data` is required to be a dictionary with the "matrix" and "index"
    keys, storing the sparse matrix and the corresponding `obs_names`.
    Differences between the index and `adata.obs_names` are not allowed.
    """
    if not adata.obs_names.equals(data["index"]):
        raise ValueError(
            f"Index of slot data `obsp/{name}.pkl` differs from "
            "`adata.obs_names`; graphs can't be re-indexed"
        )
    names = adata.obsp.keys()
    adata = handle_name_collisions(adata, "obsp", name, names, overwrite)
    adata.obsp[name] = data["matrix"]
    logger.info(f'Added `adata.obsp["{name}"]`')
    return adata


def extend_uns(adata, data, name, overwrite):
    """Add unstructured data to `adata.uns` under the key `name`."""
    adata = handle_name_collisions(adata, "uns", name, adata.uns, overwrite)
    adata.uns[name] = data
    logger.info(f'Added `adata.uns["{name}"]`')
    return adata


def extend_layer(base_adata, layer_adata, name, align, overwrite):
    """
    Add `layer_adata.X` to `base_adata.layers` under the key `name`.

    Duplicate obs/var names in either object is not allowed, though `base_adata`
    is allowed to be a subset of `layer_adata` when `align` is set; layers are
    never re-indexed.
    """
    path_label = f"layers/{name}.h5ad"
    layer_label = f'`adata.layers["{name}"]`'

    # Check both axes for duplication/equivalency before taking any actions
    needs_subset = False
    for axis in ("obs", "var"):
        base_names = getattr(base_adata, f"{axis}_names")
        layer_names = getattr(layer_adata, f"{axis}_names")

        if not base_names.is_unique:
            raise ValueError(
                f"Can't add {layer_label}: names in `adata.{axis}_names` are "
                "duplicated"
            )
        if not layer_names.is_unique:
            raise ValueError(
                f"Can't add {layer_label}: names in `{axis}_names` of "
                f"`{path_label}` are duplicated"
            )

        if base_names.equals(layer_names):
            continue
        if not align:
            raise ValueError(
                f"Index of layer data `{path_label}` differs from "
                f"`adata.{axis}_names`; set `align` to select the matching "
                "names"
            )
        n_missing = (~base_names.isin(layer_names)).sum()
        if n_missing > 0:
            raise ValueError(
                f"Can't add {layer_label}: {n_missing} of {len(base_names)} "
                f"names in `adata.{axis}_names` are missing from `{path_label}`"
            )
        needs_subset = True

    base_adata = handle_name_collisions(
        base_adata,
        "layers",
        name,
        base_adata.layers,
        overwrite
    )
    if needs_subset:
        subset = layer_adata[base_adata.obs_names, base_adata.var_names].copy()
        base_adata.layers[name] = subset.X
    else:
        base_adata.layers[name] = layer_adata.X
    logger.info(f"Added {layer_label}")

    return base_adata


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
    """Extend an adata object with slot data."""

    # Template variables
    h5ad = "${h5ad}"
    align = "${align}" == "true"
    overwrite = "${overwrite}" == "true"
    prefix = "${prefix}"
    process_name = "${task.process}"

    slot_dirs = {
        "obs": Path("obs/"),
        "var": Path("var/"),
        "obsm": Path("obsm/"),
        "varm": Path("varm/"),
        "obsp": Path("obsp/"),
        "uns": Path("uns/"),
        "layers": Path("layers/")
    }
    slot_paths = {
        slot: sorted(directory.glob("*"))
        for slot, directory in slot_dirs.items()
    }

    # Abort if no slot data is given
    if not any(slot_paths.values()):
        raise ValueError("No slot data given; there is nothing to extend")

    # Abort if any slot data has the wrong file type
    for slot, paths in slot_paths.items():
        for path in paths:
            check_suffix(slot, path)

    # `X` isn't used, so it stays on disk until the output is written
    adata = ad.read_h5ad(h5ad, backed="r")
    logger.info(f"Read base AnnData with shape {adata.shape}: {h5ad}")

    for slot, paths in slot_paths.items():
        for path in paths:
            if slot == "layers":
                data = load_layer(path)
            else:
                data = load_pickle(path)
            name = path.stem
            if slot in ("obs", "var"):
                adata = extend_frame(
                    adata,
                    data,
                    slot,
                    name,
                    align,
                    overwrite
                )
            elif slot in ("obsm", "varm"):
                adata = extend_matrix(
                    adata,
                    data,
                    slot,
                    name,
                    align,
                    overwrite
                )
            elif slot == "obsp":
                adata = extend_pairwise(
                    adata,
                    data,
                    name,
                    overwrite
                )
            elif slot == "uns":
                adata = extend_uns(
                    adata,
                    data,
                    name,
                    overwrite
                )
            elif slot == "layers":
                adata = extend_layer(
                    adata,
                    data,
                    name,
                    align,
                    overwrite
                )

    adata.write_h5ad(f"{prefix}.h5ad")
    logger.info(f"Written extended AnnData to: {prefix}.h5ad")

    write_versions(process_name)

if __name__ == "__main__":
    main()
