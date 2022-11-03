# -*- coding: utf-8 -*-
# @Author: Jarning
# @Date:   2022-01-08 16:18:28
# @Last Modified by:   Jarning
# @Last Modified time: 2022-04-15 10:52:02

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from matplotlib.pyplot import rc_context


def scvi_main(species):
    inData = "../results/02.preprocess_{species}/TCA_{species}.processed.h5ad".format(**locals())
    outData = "../results/04.scvi/TCA_{species}.scvi.h5ad".format(**locals())
    print("Loading and processing dataset...")
    adata = sc.read(inData)
    sc.pp.filter_genes(adata, min_counts=3)
    adata.layers["counts"] = adata.X.tocsr()
    adata.raw = adata  # keep full dimension safe

    print("Finding highly variable genes...")
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=2000,
        layer="counts",
        batch_key="GSE_ID",
        subset=True
    )
    print("SCVI modeling...")
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="GSM_ID")
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()
    # vae.save("../results/03.scvi/TCA_{species}.vae.model")
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    print("Building graph...")
    sc.pp.neighbors(adata, use_rep="X_scVI")
    print("Clustering...")
    sc.tl.leiden(adata)
    print("Calculating UMAP...")
    sc.tl.umap(adata)

    # save results
    print("Saving results...")
    adata.write_h5ad(outData)
    adata = sc.read(outData)
    # plots
    with rc_context({"figure.figsize": (8, 8), "figure.dpi": (300)}):
        sc.pl.umap(adata, color=["GSE_ID"], frameon=False, save="_{species}_gseid.png".format(**locals()))
        sc.pl.umap(adata, color=["Cell_type_published"], frameon=False, save="_{species}_published_celltype.png".format(**locals()))

    emb_umap = pd.DataFrame(adata.obsm['X_umap'])
    emb_umap.columns = ["UMAP_1", "UMAP_2"]
    emb_umap.index = adata.obs.index

    emb_scvi = pd.DataFrame(adata.obsm["X_scVI"])
    emb_scvi.columns = ["scVI_"+str(i+1) for i in range(emb_scvi.shape[1])]
    emb_scvi.index = adata.obs.index

    emb_umap.to_csv("../results/04.scvi/TCA_{species}.emb_umap.csv".format(**locals()))
    emb_scvi.to_csv("../results/04.scvi/TCA_{species}.emb_scvi.csv".format(**locals()))
    adata.obs.to_csv("../results/04.scvi/TCA_{species}.cellmeta.csv".format(**locals()))


if __name__ == '__main__':
    outDIR = "../results/04.scvi"
    os.path.exists(outDIR) or os.makedirs(outDIR)
    scvi_main("human")
    scvi_main("mouse")

