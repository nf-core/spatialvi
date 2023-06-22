#!/usr/bin/env python

# Load packages
import argparse
import json
import pandas as pd
from anndata import AnnData
from matplotlib.image import imread
from pathlib import Path
from scanpy import read_10x_h5
from typing import Union, Optional


# Function to read MTX
def read_visium_mtx(
    path: Union[str, Path],
    *,
    load_images: bool = True,
    library_id: Optional[str] = None,
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
        Path to a spaceranger output directory
    load_images:
        Whether or not to load images
    library_id
        Identifier for the visium library. Can be modified when concatenating multiple adata objects.

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
    adata = read_10x_h5(path / "raw_feature_bc_matrix.h5")

    adata.uns["spatial"] = dict()

    if library_id is None:
        library_id = "library_id"

    adata.uns["spatial"][library_id] = dict()

    if load_images:
        files = dict(
            tissue_positions_file=path / "tissue_positions.csv",
            scalefactors_json_file=path / "scalefactors_json.json",
            hires_image=path / "tissue_hires_image.png",
            lowres_image=path / "tissue_lowres_image.png",
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
        positions = pd.read_csv(files["tissue_positions_file"], index_col="barcode", dtype={"in_tissue": bool})
        adata.obs = adata.obs.join(positions, how="left")
        adata.obsm["spatial"] = adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy()
        adata.obs.drop(
            columns=["pxl_row_in_fullres", "pxl_col_in_fullres"],
            inplace=True,
        )

    return adata


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Load spatial transcriptomics data from MTX matrices and aligned images."
    )
    parser.add_argument(
        "--SRCountDir", metavar="SRCountDir", type=str, default=None, help="Input directory with Spaceranger data."
    )
    parser.add_argument("--outAnnData", metavar="outAnnData", type=str, default=None, help="Output h5ad file path.")
    args = parser.parse_args()

    # Read Visium data
    st_adata = read_visium_mtx(args.SRCountDir, library_id=None, load_images=True)

    # Write raw anndata to file
    st_adata.write(args.outAnnData)
