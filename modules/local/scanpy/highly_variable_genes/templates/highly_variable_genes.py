#!/usr/bin/env python3
"""
Identify highly variable genes (HVGs) in the dataset.
"""

import os
import sys

# Fix numba caching issue in read-only containers
os.environ['NUMBA_CACHE_DIR'] = '/tmp/numba_cache'
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib'
os.environ['XDG_CACHE_HOME'] = '/tmp/cache'

import platform
import importlib.metadata
import yaml
import scanpy as sc
import anndata as ad
import numpy as np

# Parameters from Nextflow
input_adata = "${adata}"
output_adata = "${prefix}" + ".h5ad"
n_top_genes = int("${n_top_genes}")
flavor = "${flavor}"

# Read AnnData
adata = ad.read_h5ad(input_adata)

print(f"Finding highly variable genes in: {input_adata}")
print(f"Shape: {adata.shape}")
print(f"Number of top genes requested: {n_top_genes}")
print(f"Flavor: {flavor}")

n_cells, n_genes = adata.shape

# Validate input
if n_cells == 0:
    print("ERROR: AnnData has 0 cells/spots.")
    sys.exit(1)

if n_genes == 0:
    print("ERROR: AnnData has 0 genes.")
    sys.exit(1)

# Handle edge case: too few genes
if n_genes < 10:
    print(f"WARNING: AnnData has only {n_genes} genes. This is too few for meaningful HVG selection.")
    print("Marking all genes as highly variable.")
    
    # Mark all genes as highly variable
    adata.var["highly_variable"] = True
    adata.var["highly_variable_rank"] = np.arange(n_genes)
    adata.var["means"] = np.array(adata.X.mean(axis=0)).flatten()
    adata.var["dispersions"] = np.zeros(n_genes)
    adata.var["dispersions_norm"] = np.zeros(n_genes)
    
    n_hvgs = n_genes
    
    # Store parameters in uns
    adata.uns["hvg"] = {
        "flavor": flavor,
        "n_top_genes": n_top_genes,
        "n_hvgs_found": int(n_hvgs),
        "warning": f"Only {n_genes} genes available, all marked as HVG"
    }

else:
    # Adjust n_top_genes if necessary
    actual_n_top_genes = min(n_top_genes, n_genes)
    
    if actual_n_top_genes < n_top_genes:
        print(f"WARNING: Requested {n_top_genes} HVGs but only {n_genes} genes available.")
        print(f"Adjusting to select top {actual_n_top_genes} genes.")
    
    # Find highly variable genes
    try:
        sc.pp.highly_variable_genes(
            adata,
            flavor=flavor,
            n_top_genes=actual_n_top_genes,
            inplace=True
        )
        
        n_hvgs = adata.var["highly_variable"].sum()
        print(f"Identified {n_hvgs} highly variable genes")
        
    except Exception as e:
        print(f"WARNING: HVG selection with flavor '{flavor}' failed: {e}")
        print("Attempting with flavor 'cell_ranger'...")
        
        try:
            sc.pp.highly_variable_genes(
                adata,
                flavor="cell_ranger",
                n_top_genes=actual_n_top_genes,
                inplace=True
            )
            
            n_hvgs = adata.var["highly_variable"].sum()
            print(f"Identified {n_hvgs} highly variable genes with cell_ranger flavor")
            flavor = "cell_ranger"
            
        except Exception as e2:
            print(f"WARNING: HVG selection with 'cell_ranger' also failed: {e2}")
            print("Falling back to marking top genes by mean expression as HVG...")
            
            # Fallback: use mean expression to select top genes
            mean_expr = np.array(adata.X.mean(axis=0)).flatten()
            top_indices = np.argsort(mean_expr)[::-1][:actual_n_top_genes]
            
            adata.var["highly_variable"] = False
            adata.var.iloc[top_indices, adata.var.columns.get_loc("highly_variable")] = True
            adata.var["means"] = mean_expr
            adata.var["highly_variable_rank"] = np.nan
            adata.var.loc[adata.var["highly_variable"], "highly_variable_rank"] = np.arange(actual_n_top_genes)
            
            n_hvgs = actual_n_top_genes
            print(f"Selected top {n_hvgs} genes by mean expression")
    
    # Store parameters in uns
    adata.uns["hvg"] = {
        "flavor": flavor,
        "n_top_genes": n_top_genes,
        "n_top_genes_actual": actual_n_top_genes,
        "n_hvgs_found": int(n_hvgs),
    }

# Ensure highly_variable column exists and is boolean
if "highly_variable" not in adata.var.columns:
    adata.var["highly_variable"] = True
else:
    adata.var["highly_variable"] = adata.var["highly_variable"].astype(bool)

# Print summary
n_hvgs = adata.var["highly_variable"].sum()
print(f"Final HVG count: {n_hvgs}")
print(f"Percentage of genes: {n_hvgs / n_genes * 100:.1f}%")

# Write output
adata.write_h5ad(output_adata)
print(f"Written AnnData with HVG annotations to: {output_adata}")

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": importlib.metadata.version("scanpy"),
        "anndata": importlib.metadata.version("anndata"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)