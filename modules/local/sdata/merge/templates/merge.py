#!/usr/bin/env python3
"""
Merge multiple SpatialData objects into one.
"""

import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import yaml
import spatialdata


def read_sdatas(file_paths):
    """Read multiple SpatialData objects from Zarr directories."""
    sdata_list = []
    for file in sorted(file_paths):
        print(f"Reading: {file}")
        sdata = spatialdata.read_zarr(file)
        sdata_list.append(sdata)
    return sdata_list


def write_versions(process_name):
    """Write software versions to YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "spatialdata": spatialdata.__version__,
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Merge multiple SpatialData objects into one."""
    input_files = "${sdata}".split()
    output_path = "${prefix}.zarr"

    # Read all SpatialData objects
    print(f"Merging {len(input_files)} SpatialData objects")
    sdata_list = read_sdatas(input_files)

    # Merge the objects into one
    output_sdata = spatialdata.concatenate(
        sdata_list,
        region_key=None,
        instance_key=None,
        concatenate_tables=False,
        obs_names_make_unique=True,
        modify_tables_inplace=False,
    )

    # Write output
    output_sdata.write(output_path, overwrite=True)
    print(f"Written merged SpatialData to: {output_path}")

    write_versions("${task.process}")


if __name__ == "__main__":
    main()
