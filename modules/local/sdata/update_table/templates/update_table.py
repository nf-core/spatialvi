#!/usr/bin/env python3
"""
Update SpatialData tables with processed AnnData objects. Also updates
associated spatial elements to match filtered observations.

Supports:
- Single-sample: replace table entirely.
- Multi-sample: update tables from concatenated AnnData using library_key.
"""

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

# Fix numba caching issue in read-only containers
os.environ['NUMBA_CACHE_DIR'] = '/tmp/numba_cache'
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib'
os.environ['XDG_CACHE_HOME'] = '/tmp/cache'

import importlib.metadata
import platform
import shutil
import sys

import anndata as ad
import numpy as np
import spatialdata
import yaml
from spatialdata.models import TableModel

# -----------------------------------------------------------------------------
# Single-sample operations
# -----------------------------------------------------------------------------


def find_table_name(sdata, adata, sample_id):
    """Determine which table to update for single-sample mode."""

    # Keep only alphanumeric characters, underscores and hyphens in sample ID
    sample_id_clean = "".join(
        c for c in sample_id if c.isalnum() or c in ["_", "-"]
    )

    if "table_name" in adata.uns:
        return adata.uns["table_name"]
    elif f"{sample_id_clean}_table" in sdata.tables:
        return f"{sample_id_clean}_table"
    elif len(sdata.tables) > 0:
        return list(sdata.tables.keys())[0]
    else:
        print("ERROR: No tables found in SpatialData object")
        sys.exit(1)


def replace_table(sdata, adata, table_name):
    """Replace a table entirely, preserving SpatialData metadata."""

    print(f"Replacing table '{table_name}'")
    # print(f"Updating table: {table_name}")
    print(f"Original table shape: {sdata.tables[table_name].shape}")
    print(f"New AnnData shape: {adata.shape}")

    original_table = sdata.tables[table_name]
    spatialdata_attrs = original_table.uns.get("spatialdata_attrs", {})
    region = spatialdata_attrs.get("region")
    region_key = spatialdata_attrs.get("region_key")
    instance_key = spatialdata_attrs.get("instance_key")
    # region = spatialdata_attrs.get("region", None)
    # region_key = spatialdata_attrs.get("region_key", None)
    # instance_key = spatialdata_attrs.get("instance_key", None)
    print(f"Original region: {region}")
    print(f"Original region_key: {region_key}")
    print(f"Original instance_key: {instance_key}")

    # Handle region being a numpy array (convert to list/string)
    if isinstance(region, np.ndarray):
        region = region.tolist()
    if isinstance(region, list) and len(region) == 1:
        region = region[0]

    # Ensure instance_key column exists in the new AnnData
    if instance_key and instance_key not in adata.obs.columns:
        if instance_key in original_table.obs.columns:
            common_idx = adata.obs.index.intersection(original_table.obs.index)
            if len(common_idx) == len(adata.obs.index):
                adata.obs[instance_key] = original_table.obs.loc[adata.obs.index, instance_key]
        else:
            print(f"WARNING: Could not match all indices for {instance_key}")
    else:
        print(f"WARNING: instance_key '{instance_key}' not found in original table")

    # Restore region_key column if missing
    if region_key and region_key not in adata.obs.columns:
        if region_key in original_table.obs.columns:
            common_idx = adata.obs.index.intersection(original_table.obs.index)
            if len(common_idx) == len(adata.obs.index):
                adata.obs[region_key] = original_table.obs.loc[adata.obs.index, region_key]
        elif region:
            region_value = region if isinstance(region, str) else region[0]
            adata.obs[region_key] = region_value

    # Create new table with proper metadata
    try:
        if region and instance_key:
            new_table = TableModel.parse(
                adata,
                region=region,
                region_key=region_key,
                instance_key=instance_key
            )
            print("Created new table using TableModel.parse()")
        else:
            # Fallback: copy uns from original
            print("WARNING: Missing region/instance_key metadata, copying from original")
            adata.uns["spatialdata_attrs"] = spatialdata_attrs.copy()
            new_table = adata
    except Exception as e:
        print(f"WARNING: TableModel.parse failed ({e}), copying from original")
        adata.uns["spatialdata_attrs"] = spatialdata_attrs.copy()
        new_table = adata

    # Replace table
    del sdata.tables[table_name]
    try:
        sdata.tables[table_name] = new_table
        print(f"Successfully updated table '{table_name}'")
    except Exception as e:
        print(f"WARNING: Failed to set table via sdata.tables: {e}")
        print("Attempting to set table using internal method...")

    # try:
    #     # Try setting directly on the internal dict
    #     if hasattr(sdata, '_tables'):
    #         sdata._tables[table_name] = new_table
    #         print("Set table using _tables dict")
    #     else:
    #         raise AttributeError("No _tables attribute found")
    # except Exception as e2:
    #     print(f"ERROR: All methods failed: {e2}")
    #     sys.exit(1)

    # Update spatial elements if observations were filtered
    if region and adata.shape[0] < original_table.shape[0]:
        print(f"Observations were filtered: {original_table.shape[0]} -> {adata.shape[0]}")
        region_name = region if isinstance(region, str) else region[0]
        # region_name = region if isinstance(region, str) else (region[0] if isinstance(region, list) else None)
        if region_name in sdata.shapes:
        # if region_name and region_name in sdata.shapes:
            try:
                matched, _ = spatialdata.match_element_to_table(
                    sdata,
                    element_name=region_name,
                    table_name=table_name
                )
                sdata.shapes[region_name] = matched[region_name]
                print(f"Updated shapes '{region_name}' to match filtered table")
            except Exception as e:
                print(f"  WARNING: Could not update shapes: {e}")
    else:
        print("No filtering detected or no region to update")


# -----------------------------------------------------------------------------
# Multi-sample operations
# -----------------------------------------------------------------------------


def update_tables_from_concat(sdata, adata_concat, library_key):
    """Update multiple tables from concatenated AnnData."""
    if library_key not in adata_concat.obs.columns:
        print(f"ERROR: Column '{library_key}' not found in AnnData")
        sys.exit(1)

    library_ids = adata_concat.obs[library_key].unique()
    print(f"Found {len(library_ids)} libraries: {library_ids.tolist()}")

    for table_name in sdata.tables.keys():
        mask = adata_concat.obs[library_key] == table_name
        if not mask.any():
            print(f"WARNING: No observations for '{table_name}'")
            continue

        adata_subset = adata_concat[mask].copy()
        table = sdata.tables[table_name]
        print(f"Updating '{table_name}' ({adata_subset.shape[0]} obs)")

        # Build index mapping (remove concatenation suffix)
        index_map = {idx: idx.rsplit("-", 1)[0] for idx in adata_subset.obs_names}

        # Add obsm fields
        for key in adata_subset.obsm.keys():
            n_features = adata_subset.obsm[key].shape[1]
            if key not in table.obsm:
                table.obsm[key] = np.full((table.n_obs, n_features), np.nan, dtype=np.float32)

            for i, new_idx in enumerate(adata_subset.obs_names):
                orig_idx = index_map[new_idx]
                if orig_idx in table.obs_names:
                    pos = table.obs_names.get_loc(orig_idx)
                    table.obsm[key][pos] = adata_subset.obsm[key][i]

            print(f"  Added obsm['{key}']")

        # Add obs columns
        skip = {library_key}
        new_cols = set(adata_subset.obs.columns) - set(table.obs.columns) - skip

        for col in new_cols:
            table.obs[col] = np.nan
            for i, new_idx in enumerate(adata_subset.obs_names):
                orig_idx = index_map[new_idx]
                if orig_idx in table.obs_names:
                    table.obs.loc[orig_idx, col] = adata_subset.obs.iloc[i][col]

            print(f"  Added obs['{col}']")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def write_versions():
    """Write software versions to YAML."""
    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "spatialdata": importlib.metadata.version("spatialdata"),
            "anndata": importlib.metadata.version("anndata"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Synchronize AnnData back to SpatialData."""

    # Template variables
    input_sdata = "${sdata}"
    input_adata = "${adata}"
    output_sdata = "${prefix}" + "_updated.zarr"
    sample_id = "${meta.id}"
    library_key = "${library_key}"

    print(f"Input SpatialData: {input_sdata}")
    print(f"Input AnnData: {input_adata}")
    print(f"Output: {output_sdata} SpatialData")
    print(f"Mode: {'multi-sample' if library_key else 'single-sample'}")

    # Validate
    if os.path.abspath(input_sdata) == os.path.abspath(output_sdata):
        print("ERROR: Input and output paths are the same!")
        sys.exit(1)

    if os.path.exists(output_sdata):
        print(f"Removing existing output directory: {output_sdata}")
        shutil.rmtree(output_sdata)

    # Read data
    sdata = spatialdata.read_zarr(input_sdata)
    adata = ad.read_h5ad(input_adata)
    print(f"Available tables: {list(sdata.tables.keys())}")
    print(f"Available shapes: {list(sdata.shapes.keys())}")

    # Process based on mode
    if library_key:
        update_tables_from_concat(sdata, adata, library_key)
    else:
        table_name = find_table_name(sdata, adata, sample_id)
        replace_table(sdata, adata, table_name)

    # Write output
    sdata.write(output_sdata)
    print(f"Written: {output_sdata}")

    write_versions()


if __name__ == "__main__":
    main()
