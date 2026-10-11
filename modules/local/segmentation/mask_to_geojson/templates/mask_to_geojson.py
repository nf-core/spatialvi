#!/usr/bin/env python3
"""
Convert a label mask to GeoJSON polygons in full-resolution pixel coordinates.
The output follows the Space Ranger segmentation GeoJSON format.
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

import cv2
import numpy as np
import tifffile
import yaml
from skimage.measure import regionprops

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def label_to_polygon(mask, prop, scale):
    """Return the outer contour of one label as a closed coordinate list."""
    min_row, min_col, max_row, max_col = prop.bbox
    binary = (mask[min_row:max_row, min_col:max_col] == prop.label).astype(np.uint8)
    contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    contour = max(contours, key=cv2.contourArea).reshape(-1, 2)
    if len(contour) < 3:
        return None
    # Contour points are pixel indices; map pixel centres back to full
    # resolution (identity when scale is 1)
    coords = (contour + [min_col, min_row] + 0.5) * scale - 0.5
    coords = coords.tolist()
    coords.append(coords[0])
    return coords


def write_versions(process_name):
    """Write software versions to a YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "opencv": cv2.__version__,
            "scikit-image": importlib.metadata.version("scikit-image"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Convert mask to GeoJSON."""

    # Template variables
    mask_file = "${mask}"
    # Factor the image was shrunk by before segmentation; polygons are scaled back up
    scale = float("${downscale}")
    object_class = "${object_class}"
    output_geojson = "${prefix}_segmentation.geojson"
    process_name = "${task.process}"

    logger.info(f"Reading mask: {mask_file}")
    mask = tifffile.imread(mask_file)
    if mask.ndim > 2:
        mask = mask.squeeze()

    features = []
    for prop in regionprops(mask):
        coords = label_to_polygon(mask, prop, scale)
        if coords is None:
            continue
        features.append({
            "type": "Feature",
            "geometry": {"type": "Polygon", "coordinates": [coords]},
            "properties": {
                "cell_id": int(prop.label),
                "classification": {"name": object_class},
            },
        })
    logger.info(f"Polygons: {len(features)}")

    with open(output_geojson, "w") as f:
        json.dump({"type": "FeatureCollection", "features": features}, f)
    logger.info(f"Written: {output_geojson}")

    write_versions(process_name)

if __name__ == "__main__":
    main()
