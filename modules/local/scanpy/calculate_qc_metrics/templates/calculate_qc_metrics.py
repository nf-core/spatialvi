#!/usr/bin/env python3
"""
Calculate QC metrics for AnnData using scanpy.

Adds the following annotations:
- var: 'mt', 'ribo', 'hb' (boolean flags for gene types)
- obs: 'n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb'
"""

# Required before importing Numba-dependent packages in read-only containers
import os
os.environ["NUMBA_CACHE_DIR"] = "/tmp/numba_cache"
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib"
os.environ["XDG_CACHE_HOME"] = "/tmp/cache"

import importlib.metadata
import platform

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse
import yaml


def validate_adata(adata):
    """Validate that AnnData has sufficient data for QC calculation."""
    n_obs, n_genes = adata.shape
    if n_obs == 0:
        raise ValueError("AnnData has 0 obs; cannot calculate QC metrics.")
    if n_genes == 0:
        raise ValueError("AnnData has 0 genes; cannot calculate QC metrics.")
    if n_genes < 10:
        print(f"WARNING: AnnData has only {n_genes} genes. This may indicate a problem with the input data.")


def annotate_gene_types(adata):
    """Annotate gene types (mitochondrial, ribosomal, haemoglobin) in var."""
    # Case-insensitive annotation
    var_names_upper = adata.var_names.str.upper()
    adata.var["mt"] = var_names_upper.str.startswith("MT-")
    adata.var["ribo"] = var_names_upper.str.match(r"^RP[SL]")
    adata.var["hb"] = var_names_upper.str.match(r"^HB[^P]")

    gene_counts = {
        "mt": adata.var["mt"].sum(),
        "ribo": adata.var["ribo"].sum(),
        "hb": adata.var["hb"].sum(),
    }

    print("Gene type counts:")
    print(f"  MT genes: {gene_counts['mt']}")
    print(f"  Ribo genes: {gene_counts['ribo']}")
    print(f"  Hb genes: {gene_counts['hb']}")

    return gene_counts


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
    percent_top = [t for t in [500, 200, 100, 50] if n_genes >= t]
    if not percent_top and n_genes >= 10:
        percent_top = [n_genes]
    return percent_top


def ensure_qc_columns_exist(adata):
    """Ensure all expected QC columns exist, adding zeros if missing."""
    for var_name in ["mt", "ribo", "hb"]:
        if f"pct_counts_{var_name}" not in adata.obs.columns:
            adata.obs[f"pct_counts_{var_name}"] = 0.0
        if f"total_counts_{var_name}" not in adata.obs.columns:
            adata.obs[f"total_counts_{var_name}"] = 0.0

    if "total_counts" not in adata.obs.columns:
        adata.obs["total_counts"] = np.array(adata.X.sum(axis=1)).flatten()
    if "n_genes_by_counts" not in adata.obs.columns:
        adata.obs["n_genes_by_counts"] = np.array((adata.X > 0).sum(axis=1)).flatten()

    return adata


def calculate_qc_metrics(adata):
    """Calculate QC metrics for AnnData object."""
    # Store raw counts layer if not already present
    if "raw" not in adata.layers:
        adata.layers["raw"] = adata.X.copy()

    # Ensure X is in sparse format for compatibility
    if not scipy.sparse.issparse(adata.X):
        adata.X = scipy.sparse.csr_matrix(adata.X)

    adata.var_names_make_unique()

    gene_counts = annotate_gene_types(adata)
    qc_vars = determine_qc_vars(gene_counts)
    percent_top = determine_percent_top(adata.shape[1])

    print(f"Using percent_top: {percent_top if percent_top else 'disabled'}")
    print(f"Using qc_vars: {qc_vars if qc_vars else 'none'}")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=qc_vars,
        percent_top=percent_top if percent_top else None,
        inplace=True,
        log1p=False,
    )

    adata = ensure_qc_columns_exist(adata)

    return adata


def print_qc_summary(adata):
    """Print summary of QC metrics."""
    obs = adata.obs

    print("QC metrics summary:")
    print(f"  Median total counts: {obs['total_counts'].median():.1f}")
    print(f"  Median genes per obs: {obs['n_genes_by_counts'].median():.1f}")
    print(f"  Median MT %: {obs['pct_counts_mt'].median():.2f}")
    print(f"  Median Ribo %: {obs['pct_counts_ribo'].median():.2f}")
    print(f"  Median Hb %: {obs['pct_counts_hb'].median():.2f}")
    print(
        f"Total counts range: "
        f"[{obs['total_counts'].min():.0f}, {obs['total_counts'].max():.0f}]"
    )
    print(
        f"Genes per obs range: "
        f"[{obs['n_genes_by_counts'].min():.0f}, "
        f"{obs['n_genes_by_counts'].max():.0f}]"
    )


def write_versions(process_name):
    """Write software versions to a YAML file."""
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

    print(f"Calculating QC metrics for sample: {input_adata}")
    adata = ad.read_h5ad(input_adata)
    print(f"AnnData shape: {adata.shape}")

    validate_adata(adata)
    adata = calculate_qc_metrics(adata)
    print_qc_summary(adata)

    adata.write_h5ad(output_adata)
    print(f"Written AnnData with QC metrics to: {output_adata}")

    write_versions(process_name)


if __name__ == "__main__":
    main()
