#!/usr/bin/env python3
"""
Calculate QC metrics for AnnData using scanpy.

Adds the following annotations:
- var: 'mt', 'ribo', 'hb' (boolean flags for gene types)
- obs: 'n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb'
"""

# Fix numba caching issue in read-only containers
import os
os.environ['NUMBA_CACHE_DIR'] = '/tmp/numba_cache'
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib'
os.environ['XDG_CACHE_HOME'] = '/tmp/cache'

import importlib.metadata
import platform
import sys
import yaml

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse


def validate_adata(adata, input_path):
    """Validate that AnnData has sufficient data for QC calculation."""
    n_cells, n_genes = adata.shape

    if n_cells == 0:
        print("ERROR: AnnData has 0 cells/spots. Cannot calculate QC metrics.")
        sys.exit(1)

    if n_genes == 0:
        print("ERROR: AnnData has 0 genes. Cannot calculate QC metrics.")
        sys.exit(1)

    if n_genes < 10:
        print(f"WARNING: AnnData has only {n_genes} genes. This may indicate a problem with the input data.")

    print(f"Calculating QC metrics for: {input_path}")
    print(f"Shape: {adata.shape}")


def annotate_gene_types(adata):
    """Annotate gene types (mitochondrial, ribosomal, haemoglobin) in var."""
    # Mitochondrial genes (MT- prefix, case-insensitive)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

    # Ribosomal genes (RPS/RPL prefix, case-insensitive)
    adata.var["ribo"] = adata.var_names.str.upper().str.match(r"^RP[SL]")

    # Haemoglobin genes (HBA/HBB prefix, case-insensitive)
    adata.var["hb"] = adata.var_names.str.upper().str.match(r"^HB[^P]")

    # Print gene type counts
    n_mt = adata.var["mt"].sum()
    n_ribo = adata.var["ribo"].sum()
    n_hb = adata.var["hb"].sum()

    print("Gene type counts:")
    print(f"  MT genes: {n_mt}")
    print(f"  Ribo genes: {n_ribo}")
    print(f"  Hb genes: {n_hb}")

    return {"mt": n_mt, "ribo": n_ribo, "hb": n_hb}


def determine_qc_vars(gene_counts):
    """Determine which qc_vars to use based on gene counts."""
    qc_vars = []
    for var_name, count in gene_counts.items():
        if count > 0:
            qc_vars.append(var_name)
        else:
            print(f"WARNING: No {var_name} genes found in the dataset.")
    return qc_vars


def determine_percent_top(n_genes):
    """Determine percent_top parameter based on number of genes."""
    if n_genes >= 500:
        return [50, 100, 200, 500]
    elif n_genes >= 200:
        return [50, 100, 200]
    elif n_genes >= 100:
        return [50, 100]
    elif n_genes >= 50:
        return [50]
    elif n_genes >= 10:
        return [n_genes]
    else:
        return []


def calculate_qc_metrics(adata):
    """Calculate QC metrics for AnnData object."""
    # Store raw counts layer if not already present
    if "raw" not in adata.layers:
        adata.layers["raw"] = adata.X.copy()

    # Ensure X is in sparse format for compatibility
    if not scipy.sparse.issparse(adata.X):
        adata.X = scipy.sparse.csr_matrix(adata.X)

    # Make var names unique before processing
    adata.var_names_make_unique()

    # Annotate gene types
    gene_counts = annotate_gene_types(adata)

    # Determine QC parameters
    qc_vars = determine_qc_vars(gene_counts)
    percent_top = determine_percent_top(adata.shape[1])

    print(f"Using percent_top: {percent_top if percent_top else 'disabled'}")
    print(f"Using qc_vars: {qc_vars if qc_vars else 'none'}")

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=qc_vars if qc_vars else [],
        percent_top=percent_top if percent_top else None,
        inplace=True,
        log1p=False
    )

    # Ensure all expected columns exist (add zeros if missing)
    for var_name in ["mt", "ribo", "hb"]:
        if f"pct_counts_{var_name}" not in adata.obs.columns:
            adata.obs[f"pct_counts_{var_name}"] = 0.0
        if f"total_counts_{var_name}" not in adata.obs.columns:
            adata.obs[f"total_counts_{var_name}"] = 0.0

    # Ensure basic columns exist
    if "total_counts" not in adata.obs.columns:
        adata.obs["total_counts"] = np.array(adata.X.sum(axis=1)).flatten()
    if "n_genes_by_counts" not in adata.obs.columns:
        adata.obs["n_genes_by_counts"] = np.array((adata.X > 0).sum(axis=1)).flatten()

    return adata


def print_qc_summary(adata):
    """Print summary of QC metrics."""
    print("QC metrics summary:")
    print(f"  Median total counts: {adata.obs['total_counts'].median():.1f}")
    print(f"  Median genes per spot: {adata.obs['n_genes_by_counts'].median():.1f}")
    print(f"  Median MT %: {adata.obs['pct_counts_mt'].median():.2f}")
    print(f"  Median Ribo %: {adata.obs['pct_counts_ribo'].median():.2f}")
    print(f"  Median Hb %: {adata.obs['pct_counts_hb'].median():.2f}")
    print(f"Total counts range: [{adata.obs['total_counts'].min():.0f}, {adata.obs['total_counts'].max():.0f}]")
    print(f"Genes per spot range: [{adata.obs['n_genes_by_counts'].min():.0f}, {adata.obs['n_genes_by_counts'].max():.0f}]")


def write_versions(process_name):
    """Write software versions to YAML file."""
    versions = {
        process_name: {
            "python": platform.python_version(),
            "scanpy": importlib.metadata.version("scanpy"),
            "anndata": importlib.metadata.version("anndata"),
        }
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


def main():
    """Calculate QC metrics for an AnnData object."""

    # Template variables
    input_adata = "${adata}"
    output_adata = "${prefix}.h5ad"
    process_name = "${task.process}"

    # Read AnnData
    adata = ad.read_h5ad(input_adata)

    # Validate input
    validate_adata(adata, input_adata)

    # Calculate QC metrics
    adata = calculate_qc_metrics(adata)

    # Print summary
    print_qc_summary(adata)

    # Write output
    adata.write_h5ad(output_adata)
    print(f"Written AnnData with QC metrics to: {output_adata}")

    # Write versions
    write_versions(process_name)


if __name__ == "__main__":
    main()
