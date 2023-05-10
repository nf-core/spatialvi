#!/usr/bin/env python

# Load packages
import argparse
import json
import pandas as pd
from anndata import AnnData
from matplotlib.image import imread
from pathlib import Path
from scanpy import read_10x_mtx
from typing import Union, Optional


# Function to read MTX
def read_visium_mtx(
    path: Union[str, Path],
    genome: Optional[str] = None,
    *,
    count_file: str = "filtered_feature_bc_matrix.h5",
    library_id: str = None,
    load_images: Optional[bool] = True,
    source_image_path: Optional[Union[str, Path]] = None,
) -> AnnData:
    """\
    Read 10x-Genomics-formatted visum dataset.
    In addition to reading regular 10x output,
    this looks for the `spatial` folder and loads images,
    coordinates and scale factors.
    Based on the `Space Ranger output docs`_.
    See :func:`~scanpy.pl.spatial` for a compatible plotting function.
    .. _Space Ranger output docs: https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/overview
    Parameters
    ----------
    path
        Path to directory for visium datafiles.
    genome
        Filter expression to genes within this genome.
    count_file
        Which file in the passed directory to use as the count file. Typically would be one of:
        'filtered_feature_bc_matrix.h5' or 'raw_feature_bc_matrix.h5'.
    library_id
        Identifier for the visium library. Can be modified when concatenating multiple adata objects.
    source_image_path
        Path to the high-resolution tissue image. Path will be included in
        `.uns["spatial"][library_id]["metadata"]["source_image_path"]`.
    Returns
    -------
    Annotated data matrix, where observations/cells are named by their
    barcode and variables/genes by gene name. Stores the following information:
    :attr:`~anndata.AnnData.X`
        The data matrix is stored
    :attr:`~anndata.AnnData.obs_names`
        Cell names
    :attr:`~anndata.AnnData.var_names`
        Gene names
    :attr:`~anndata.AnnData.var`\\ `['gene_ids']`
        Gene IDs
    :attr:`~anndata.AnnData.var`\\ `['feature_types']`
        Feature types
    :attr:`~anndata.AnnData.uns`\\ `['spatial']`
        Dict of spaceranger output files with 'library_id' as key
    :attr:`~anndata.AnnData.uns`\\ `['spatial'][library_id]['images']`
        Dict of images (`'hires'` and `'lowres'`)
    :attr:`~anndata.AnnData.uns`\\ `['spatial'][library_id]['scalefactors']`
        Scale factors for the spots
    :attr:`~anndata.AnnData.uns`\\ `['spatial'][library_id]['metadata']`
        Files metadata: 'chemistry_description', 'software_version', 'source_image_path'
    :attr:`~anndata.AnnData.obsm`\\ `['spatial']`
        Spatial spot coordinates, usable as `basis` by :func:`~scanpy.pl.embedding`.
    """

    path = Path(path)
    adata = read_10x_mtx(path / "raw_feature_bc_matrix")

    adata.uns["spatial"] = dict()

    if library_id is None:
        library_id = "library_id"

    adata.uns["spatial"][library_id] = dict()

    if load_images:
        files = dict(
            tissue_positions_file=path / "spatial/tissue_positions_list.csv",
            scalefactors_json_file=path / "spatial/scalefactors_json.json",
            hires_image=path / "spatial/tissue_hires_image.png",
            lowres_image=path / "spatial/tissue_lowres_image.png",
        )

        # Check if files exist; continue if images are missing
        for f in files.values():
            if not f.exists():
                if any(x in str(f) for x in ["hires_image", "lowres_image"]):
                    print("You seem to be missing an image file.")
                    print("Could not find '{f}'.")
                else:
                    raise OSError(f"Could not find '{f}'")

        # Check for existance of images
        adata.uns["spatial"][library_id]["images"] = dict()
        for res in ["hires", "lowres"]:
            try:
                adata.uns["spatial"][library_id]["images"][res] = imread(str(files[f"{res}_image"]))
            except Exception:
                raise OSError(f"Could not find '{res}_image'")

        # Read JSON scale factors
        adata.uns["spatial"][library_id]["scalefactors"] = json.loads(files["scalefactors_json_file"].read_bytes())
        adata.uns["spatial"][library_id]["metadata"] = {k: "NA" for k in ("chemistry_description", "software_version")}

        # Read coordinates
        positions = pd.read_csv(files["tissue_positions_file"], header=None)
        positions.columns = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_col_in_fullres",
            "pxl_row_in_fullres",
        ]
        positions.index = positions["barcode"]
        adata.obs = adata.obs.join(positions, how="left")
        adata.obsm["spatial"] = adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy()
        adata.obs.drop(
            columns=["barcode", "pxl_row_in_fullres", "pxl_col_in_fullres"],
            inplace=True,
        )

        # Put absolute image path in uns
        if source_image_path is not None:
            source_image_path = str(Path(source_image_path).resolve())
            adata.uns["spatial"][library_id]["metadata"]["source_image_path"] = str(source_image_path)

    return adata


# Parse command-line arguments
parser = argparse.ArgumentParser(description="Load spatial transcriptomics data from MTX matrices and aligned images.")
parser.add_argument(
    "--SRCountDir", metavar="SRCountDir", type=str, default=None, help="Input directory with Spaceranger data."
)
parser.add_argument("--outAnnData", metavar="outAnnData", type=str, default=None, help="Output h5ad file path.")
args = parser.parse_args()

# Read Visium data
st_adata = read_visium_mtx(args.SRCountDir, library_id=None, load_images=True, source_image_path=None)

# Write raw anndata to file
st_adata.write(args.outAnnData)
