#!/usr/bin/env python3
"""
Merge multiple SpatialData objects into one.
"""

import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import platform

import spatialdata
import yaml


def read_sdatas(file_paths):
    """Read multiple SpatialData objects from Zarr directories."""
    sdata_list = []
    for file in sorted(file_paths):
        print(f"Reading: {file}")
        sdata = spatialdata.read_zarr(file)
        sdata_list.append(sdata)
    return sdata_list


def write_versions(process_name):
    """Write software versions to a YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "spatialdata": importlib.metadata.version("spatialdata"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Merge multiple SpatialData objects into one."""
    # Template variables
    zarrs = "${sdata}".split()
    output_zarr = "${prefix}.zarr"
    process_name = "${task.process}"

    sdata_list = read_sdatas(zarrs)
    print(f"Merging {len(zarrs)} SpatialData objects")

    # TODO: Should some of these be template variables that the user can change?
    output_sdata = spatialdata.concatenate(
        sdata_list,
        region_key=None,
        instance_key=None,
        concatenate_tables=False,
        obs_names_make_unique=True,
        modify_tables_inplace=False,
    )

    output_sdata.write(output_zarr, overwrite=True)
    print(f"Written merged SpatialData to: {output_zarr}")

    write_versions(process_name)


if __name__ == "__main__":
    main()
