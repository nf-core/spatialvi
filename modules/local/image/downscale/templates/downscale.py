#!/usr/bin/env python3
"""
Downscale a full-resolution tissue image for segmentation.
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
import platform

import cv2
import numpy as np
import tifffile
import yaml

logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def to_rgb(img):
    """Return an (H, W, 3) uint8-compatible array."""
    if img.ndim == 3 and img.shape[0] <= 4 and img.shape[0] < img.shape[2]:
        img = np.moveaxis(img, 0, -1)
    if img.ndim == 3 and img.shape[2] == 4:
        img = img[..., :3]
    if img.ndim == 2:
        img = np.stack([img] * 3, axis=-1)
    return img


def write_versions(process_name):
    """Write software versions to a YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "opencv": cv2.__version__,
            "tifffile": importlib.metadata.version("tifffile"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Downscale the image by an integer factor."""

    # Template variables
    image = "${image}"
    downscale = int("${downscale}")
    output_image = "${prefix}_downscaled.tif"
    process_name = "${task.process}"

    logger.info(f"Reading image: {image}")
    if image.lower().endswith((".tif", ".tiff", ".btf")):
        img = tifffile.imread(image)
    else:
        # JPEG or PNG: OpenCV reads colour images as BGR
        img = cv2.imread(image, cv2.IMREAD_UNCHANGED)
        if img is None:
            raise ValueError(f"Could not read image: {image}")
        if img.ndim == 3:
            img = cv2.cvtColor(img, cv2.COLOR_BGRA2RGBA if img.shape[2] == 4 else cv2.COLOR_BGR2RGB)
    img = to_rgb(img)
    height, width = img.shape[:2]
    logger.info(f"Image size: {width} x {height}, downscale factor: {downscale}")

    out = cv2.resize(
        img, (width // downscale, height // downscale), interpolation=cv2.INTER_AREA
    )
    tifffile.imwrite(output_image, out, compression="zlib")
    logger.info(f"Written: {output_image} ({out.shape[1]} x {out.shape[0]})")

    write_versions(process_name)

if __name__ == "__main__":
    main()
