#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# v3 was the version from claude


# In[ ]:


# v4 is the version from scGPT


# In[1]:


from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import mmread
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt

# ----------------------------
# CONFIG
# ----------------------------
SHOW = True  # True in Jupyter; set False for headless/batch runs on cluster

base_dir = Path("/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/JK4488_01_combined/pooled")
velocyto_dir = base_dir / "JK4488_01_starsolo_Velocyto_filteredMatrix"
output_dir = base_dir / "velocity_analysis_starsolo_output_v4"
output_dir.mkdir(parents=True, exist_ok=True)

# Make scanpy/scvelo save figures into output_dir/figures by default
sc.settings.figdir = str(output_dir)
scv.settings.figdir = str(output_dir)

scv.settings.verbosity = 3
scv.settings.presenter_view = True

# Notebook-friendly defaults (optional)
# sc.settings.set_figure_params(dpi=120, facecolor="white")

import matplotlib as mpl
mpl.rcParams["figure.dpi"] = 120
mpl.rcParams["savefig.dpi"] = 300
mpl.rcParams["figure.facecolor"] = "white"
mpl.rcParams["axes.facecolor"] = "white"

# If running headless, use a non-interactive matplotlib backend
if not SHOW:
    import matplotlib
    matplotlib.use("Agg")

# ----------------------------
# LOAD MATRICES (cells x genes)
# ----------------------------
spliced   = mmread(velocyto_dir / "spliced.mtx").T.tocsr()
unspliced = mmread(velocyto_dir / "unspliced.mtx").T.tocsr()
ambig     = mmread(velocyto_dir / "ambiguous.mtx").T.tocsr()

barcodes = pd.read_csv(velocyto_dir / "barcodes.tsv", header=None)[0].astype(str).values

# features.tsv can be 1-col or 3-col; handle both
feat = pd.read_csv(velocyto_dir / "features.tsv", header=None, sep="\t")
if feat.shape[1] >= 2:
    genes = feat[1].astype(str).values  # gene_name
else:
    genes = feat[0].astype(str).values  # fallback

print(f"Spliced:   {spliced.shape[0]:,} cells x {spliced.shape[1]:,} genes")
print(f"Unspliced: {unspliced.shape[0]:,} cells x {unspliced.shape[1]:,} genes")
print(f"Ambig:     {ambig.shape[0]:,} cells x {ambig.shape[1]:,} genes")

# ----------------------------
# BUILD ANNDATA
# ----------------------------
adata = sc.AnnData(X=spliced)
adata.obs_names = barcodes
adata.var_names = genes
adata.var_names_make_unique()  # important if gene names repeat

adata.layers["spliced"] = spliced
adata.layers["unspliced"] = unspliced
adata.layers["ambiguous"] = ambig

# ----------------------------
# BASIC QC (customize thresholds)
# ----------------------------
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

print("Before filter:", adata.n_obs, "cells;", adata.n_vars, "genes")

# explicit thresholds
adata = adata[adata.obs["n_genes_by_counts"] >= 200, :].copy()
adata = adata[adata.obs["n_genes_by_counts"] <= 12000, :].copy()
adata = adata[adata.obs["pct_counts_mt"] < 20, :].copy()

print("After filter: ", adata.n_obs, "cells;", adata.n_vars, "genes")

# ----------------------------
# scVelo velocity pipeline (recommended minimal)
# ----------------------------
# Note: filter_and_normalize already does some HVG selection; it will modify adata.X.
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# stochastic is fine; try dynamical later if needed
scv.tl.velocity(adata, mode="stochastic")
scv.tl.velocity_graph(adata)

# embedding + clustering for visualization
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# ----------------------------
# PLOTS (display in notebook + save)
# ----------------------------
# IMPORTANT: scanpy/scvelo "save=" expects a suffix string.
# It will save into: sc.settings.figdir / "figures" / "<plotname><suffix>"
# Example: umap => "figures/umap_umap_clusters.png"
#
# If you *really* want exact filenames, we can switch to sc.pl.umap(..., show=False) then savefig manually.

sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="on data",
    show=SHOW,
    save="_umap_clusters.png",
)

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="leiden",
    legend_loc="right margin",
    show=SHOW,
    save="_velocity_stream.png",
)

scv.pl.velocity_embedding(
    adata,
    basis="umap",
    color="leiden",
    arrow_length=3,
    arrow_size=2,
    show=SHOW,
    save="_velocity_arrows.png",
)

scv.tl.velocity_confidence(adata)
scv.pl.scatter(
    adata,
    basis="umap",
    c="velocity_confidence",
    show=SHOW,
    save="_velocity_confidence.png",
)

# If you add ANY matplotlib plots (ax=...), do:
# fig.savefig(...)
# if SHOW: plt.show()
# plt.close(fig)

# ----------------------------
# SAVE RESULTS
# ----------------------------
adata.write(output_dir / "velocity_analysis.h5ad")
print("Saved:", output_dir / "velocity_analysis.h5ad")

# Where the figures actually go (scanpy/scvelo create a "figures/" subfolder)
fig_dir = Path(sc.settings.figdir) / "figures"
print("Figures in:", fig_dir)


# In[ ]:





# In[ ]:




