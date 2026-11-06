#!/usr/bin/env python3
"""
Group Visium HD 2 µm bins into cells with bin2cell, using a segmentation
(label mask or GeoJSON) in full-resolution image coordinates. The cell table
becomes the main table of the SpatialData object.
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
import json
import logging
import platform
import re

import bin2cell as b2c
import numpy as np
import scanpy as sc
import shapely
import spatialdata
import tifffile
import yaml
from spatialdata.models import ShapesModel, TableModel
from spatialdata.transformations import get_transformation
from threadpoolctl import threadpool_limits

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Limit BLAS/OpenMP threads to the allocated CPUs
threadpool_limits(int("${task.cpus}"))

NO_CELLS_MESSAGE = (
    "No bins were assigned to a cell. The segmentation must cover the tissue "
    "and be in full-resolution image pixel coordinates."
)


def is_geojson(segmentation):
    """Return True if the segmentation file is GeoJSON rather than a mask."""
    return segmentation.lower().endswith((".geojson", ".json"))


def read_geojson_features(geojson_file):
    """Return the features of a GeoJSON file."""
    with open(geojson_file) as f:
        return json.load(f)["features"]


def labels_from_geojson(features, xy):
    """Assign each bin centroid to the polygon containing it."""
    polygons = [shapely.geometry.shape(feat["geometry"]) for feat in features]
    tree = shapely.STRtree(polygons)
    bin_idx, poly_idx = tree.query(shapely.points(xy), predicate="within")
    labels = np.zeros(len(xy), dtype=int)
    labels[bin_idx] = poly_idx + 1
    return labels


def geojson_object_class(features):
    """Return the classification name of the first feature, if any."""
    if not features:
        return None
    properties = features[0].get("properties") or {}
    classification = properties.get("classification") or {}
    return classification.get("name")


def labels_from_mask(mask_file, xy):
    """Look up the label under each bin centroid."""
    if mask_file.endswith(".npy"):
        mask = np.load(mask_file)
    else:
        mask = tifffile.imread(mask_file).squeeze()
    rows = np.clip(xy[:, 1].astype(int), 0, mask.shape[0] - 1)
    cols = np.clip(xy[:, 0].astype(int), 0, mask.shape[1] - 1)
    return mask[rows, cols].astype(int)


def group_bins(spaceranger_dir, segmentation):
    """Run bin2cell on the 2 µm bins and return a cell-level AnnData."""
    adata = b2c.read_visium(os.path.join(spaceranger_dir, "binned_outputs", "square_002um"))
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_counts=1)

    xy = np.asarray(adata.obsm["spatial"])
    object_class = None
    if is_geojson(segmentation):
        features = read_geojson_features(segmentation)
        object_class = geojson_object_class(features)
        labels = labels_from_geojson(features, xy)
    else:
        labels = labels_from_mask(segmentation, xy)
    n_labelled = int((labels > 0).sum())
    logger.info(f"Bins with a label: {n_labelled} / {adata.n_obs}")
    if n_labelled == 0:
        raise ValueError(NO_CELLS_MESSAGE)
    adata.obs["labels"] = labels

    # Nuclei are expanded into cells; whole-cell segmentations are not. As in
    # Space Ranger, only a GeoJSON classified as "cell" is taken as whole
    # cells; masks and unclassified GeoJSON are treated as nuclei.
    whole_cells = object_class == "cell"
    if whole_cells:
        logger.info("Segmentation is whole cells: not expanding labels")
        labels_key = "labels"
    else:
        logger.info("Segmentation is nuclei: expanding labels into cells")
        b2c.expand_labels(adata, labels_key="labels", expanded_labels_key="labels_expanded")
        labels_key = "labels_expanded"
    cells = b2c.bin_to_cell(adata, labels_key=labels_key, spatial_keys=["spatial"])
    logger.info(f"Cells: {cells.n_obs}")
    if cells.n_obs == 0:
        raise ValueError(NO_CELLS_MESSAGE)

    library_id = list(adata.uns["spatial"].keys())[0]
    mpp = adata.uns["spatial"][library_id]["scalefactors"]["microns_per_pixel"]
    return cells, mpp


def add_cells_to_sdata(sdata, cells, mpp, sample_id, hd_bin_size):
    """Add cells as circles and the cell table as the main table."""
    bins_name = f"{sample_id}_square_{hd_bin_size:03d}um"
    cells_name = f"{sample_id}_cells"
    table_name = f"{sample_id}_table"

    # Circle radius from the number of 2 µm bins per cell
    radius = np.sqrt(cells.obs["bin_count"].to_numpy() * 4 / np.pi) / mpp
    sdata.shapes[cells_name] = ShapesModel.parse(
        np.asarray(cells.obsm["spatial"]),
        geometry=0,
        radius=radius,
        index=cells.obs_names.astype(int),
        transformations=get_transformation(sdata.shapes[bins_name], get_all=True),
    )

    cells.obs["cell_id"] = cells.obs_names.astype(int)
    cells.obs["region"] = cells_name
    cells.obs["region"] = cells.obs["region"].astype("category")
    # bin2cell copies the images into uns; SpatialData already holds them
    cells.uns.pop("spatial", None)

    sdata.tables[f"{sample_id}_bin_table"] = sdata.tables[table_name]
    del sdata.tables[table_name]
    sdata.tables[table_name] = TableModel.parse(
        cells, region=cells_name, region_key="region", instance_key="cell_id"
    )


def write_versions(process_name):
    """Write software versions to a YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "bin2cell": importlib.metadata.version("bin2cell"),
            "spatialdata": importlib.metadata.version("spatialdata"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Group bins into cells and update the SpatialData object."""

    # Template variables
    input_sdata = "${sdata}"
    spaceranger_dir = "${spaceranger_dir}"
    segmentation = "${segmentation}"
    hd_bin_size = int("${hd_bin_size}")
    sample_id = "${meta.id}"
    output_sdata = "${prefix}.zarr"
    process_name = "${task.process}"

    # Sample ID must only contain alphanumerics, underscores and dashes
    sample_id = re.sub(r"[^a-zA-Z0-9_-]", "", sample_id)

    cells, mpp = group_bins(spaceranger_dir, segmentation)

    sdata = spatialdata.read_zarr(input_sdata)
    add_cells_to_sdata(sdata, cells, mpp, sample_id, hd_bin_size)

    logger.info(f"Writing SpatialData to: {output_sdata}")
    sdata.write(output_sdata, overwrite=True)

    write_versions(process_name)

if __name__ == "__main__":
    main()
