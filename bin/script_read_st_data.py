#!/opt/conda/bin/python

# Load packages 
import os
import sys
import argparse
import scanpy as sc
import numpy as np

from scanpy import read_10x_mtx
from pathlib import Path
from typing import Union, Dict, Optional
import json
import numpy as np
import pandas as pd
from matplotlib.image import imread
import anndata
from anndata import AnnData, read_csv


# Parse command-line arguments
parser = argparse.ArgumentParser(description='Load spatial traqnscriptomics data from MTX or HDF5 count matrices and aligned images.')

parser.add_argument('--outsPath', metavar='outspath', type=str, default=None, help='Path to Space Range outs directory, etc.')
parser.add_argument('--saveFile', metavar='savefile', type=str, default=None, help='Path to a file to save h5ad data into.')
parser.add_argument('--countsFile', metavar='countsfile', type=str, default='raw_feature_bc_matrix.h5', help='Name of the HDF5 file.')
parser.add_argument('--npCountsOutputName', metavar='npzoutput', type=str, default=None, help='Name of the NPZ file.')

parser.add_argument('--minCounts', metavar='cutoff', type=int, default=1, help='Min counts per spot.')
parser.add_argument('--minCells', metavar='cutoff', type=int, default=1, help='Min cells per gene.')

args = parser.parse_args()


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
    #adata = read_10x_h5(path / count_file, genome=genome)
    adata = read_10x_mtx(path / 'raw_feature_bc_matrix')
    
    adata.uns["spatial"] = dict()

    #from h5py import File

    #with File(path / 'raw_feature_bc_matrix' / count_file, mode="r") as f:
    #    attrs = dict(f.attrs)
    #if library_id is None:
    #    library_id = str(attrs.pop("library_ids")[0], "utf-8")
    if library_id is None:
        library_id = 'library_id'

    adata.uns["spatial"][library_id] = dict()

    if load_images:
        files = dict(
            tissue_positions_file=path / 'spatial/tissue_positions_list.csv',
            scalefactors_json_file=path / 'spatial/scalefactors_json.json',
            hires_image=path / 'spatial/tissue_hires_image.png',
            lowres_image=path / 'spatial/tissue_lowres_image.png',
        )

        # check if files exists, continue if images are missing
        for f in files.values():
            if not f.exists():
                if any(x in str(f) for x in ["hires_image", "lowres_image"]):
                    print("You seem to be missing an image file.")
                    print("Could not find '{f}'.")
                    #logg.warning(
                    #    f"You seem to be missing an image file.\n"
                    #    f"Could not find '{f}'."
                    #)
                else:
                    raise OSError(f"Could not find '{f}'")

        adata.uns["spatial"][library_id]['images'] = dict()
        for res in ['hires', 'lowres']:
            try:
                adata.uns["spatial"][library_id]['images'][res] = imread(
                    str(files[f'{res}_image'])
                )
            except Exception:
                raise OSError(f"Could not find '{res}_image'")

        # read json scalefactors
        adata.uns["spatial"][library_id]['scalefactors'] = json.loads(
            files['scalefactors_json_file'].read_bytes()
        )

        #adata.uns["spatial"][library_id]["metadata"] = {
        #    k: (str(attrs[k], "utf-8") if isinstance(attrs[k], bytes) else attrs[k])
        #    for k in ("chemistry_description", "software_version")
        #    if k in attrs
        #}
        
        adata.uns["spatial"][library_id]["metadata"] = {k: "NA" for k in ("chemistry_description", "software_version")}

        # read coordinates
        positions = pd.read_csv(files['tissue_positions_file'], header=None)
        positions.columns = [
            'barcode',
            'in_tissue',
            'array_row',
            'array_col',
            'pxl_col_in_fullres',
            'pxl_row_in_fullres',
        ]
        positions.index = positions['barcode']

        adata.obs = adata.obs.join(positions, how="left")

        adata.obsm['spatial'] = adata.obs[
            ['pxl_row_in_fullres', 'pxl_col_in_fullres']
        ].to_numpy()
        adata.obs.drop(
            columns=['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'],
            inplace=True,
        )

        # put image path in uns
        if source_image_path is not None:
            # get an absolute path
            source_image_path = str(Path(source_image_path).resolve())
            adata.uns["spatial"][library_id]["metadata"]["source_image_path"] = str(
                source_image_path
            )

    return adata

for fname in os.listdir(outsPath):
    if args.countsFile in fname:
        args.countsFile = fname
        break


# Main script
if args.countsFile in os.listdir(args.outsPath):
    st_adata = sc.read_visium(args.outsPath, count_file=args.countsFile, library_id=None, load_images=True, source_image_path=None)
else:
    st_adata = read_visium_mtx(args.outsPath)

st_adata.var_names_make_unique()
sc.pp.filter_cells(st_adata, min_counts=args.minCounts)
sc.pp.filter_genes(st_adata, min_cells=args.minCells)

if not os.path.exists(os.path.dirname(args.saveFile)):
    os.makedirs(os.path.dirname(args.saveFile))

st_adata.write(args.saveFile)

X = np.array(st_adata[st_adata.obs['in_tissue']==1].X.todense()).T
np.savez_compressed(os.path.dirname(args.saveFile) + '/' + args.npCountsOutputName, X)

exit(0)
