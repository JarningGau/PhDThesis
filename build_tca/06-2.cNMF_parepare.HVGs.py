# -*- coding: utf-8 -*-
# @Author: Jarning
# @Date:   2022-04-10 15:15:07
# @Last Modified by:   Jarning
# @Last Modified time: 2022-04-15 22:15:51
import os
import anndata
import numpy as np
import pandas as pd
import scanpy as sc

def mouse_main():
    OUTDIR = "../results/06.cNMF_mouse"
    os.path.exists(OUTDIR) or os.makedirs(OUTDIR)
    print("Loading data...")
    adata = sc.read("{OUTDIR}/mTCA.counts.h5ad".format(**locals()))
    ## HVGs from all genes
    print("Calculating HVGs from all genes...")
    n_top_genes = [2000, 4000, 6000]
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=max(n_top_genes),
        batch_key="GSE_ID",
        subset=True
    )
    hvg_df = adata.var
    hvg_df.sort_values(by=['variances_norm'], inplace=True, ascending=False)
    hvg_df.to_csv("{OUTDIR}/HVGs.all.info.csv".format(**locals()))
    for topN in n_top_genes:
        f = open("{OUTDIR}/HVGs.{topN}.all.txt".format(**locals()), 'w')
        for i in hvg_df.name.tolist()[0:topN]:
            f.write(i+'\n')
    ## HVGs from coding genes
    print("Calculating HVGs from coding genes...")
    adata = sc.read("{OUTDIR}/mTCA.counts.protein_coding.h5ad".format(**locals()))
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=max(n_top_genes),
        batch_key="GSE_ID",
        subset=True
    )
    hvg_df = adata.var
    hvg_df.sort_values(by=['variances_norm'], inplace=True, ascending=False)
    hvg_df.to_csv("{OUTDIR}/HVGs.protein_coding.info.csv".format(**locals()))
    for topN in n_top_genes:
        f = open("{OUTDIR}/HVGs.{topN}.protein_coding.txt".format(**locals()), 'w')
        for i in hvg_df.name.tolist()[0:topN]:
            f.write(i+'\n')


def human_main():
    OUTDIR = "../results/06.cNMF_human"
    os.path.exists(OUTDIR) or os.makedirs(OUTDIR)
    print("Loading data...")
    adata = sc.read("{OUTDIR}/hTCA.counts.h5ad".format(**locals()))
    ## HVGs from all genes
    print("Calculating HVGs from all genes...")
    n_top_genes = [2000, 4000, 6000]
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=max(n_top_genes),
        batch_key="GSE_ID",
        subset=True
    )
    hvg_df = adata.var
    hvg_df.sort_values(by=['variances_norm'], inplace=True, ascending=False)
    hvg_df.to_csv("{OUTDIR}/HVGs.all.info.csv".format(**locals()))
    for topN in n_top_genes:
        f = open("{OUTDIR}/HVGs.{topN}.all.txt".format(**locals()), 'w')
        for i in hvg_df.name.tolist()[0:topN]:
            f.write(i+'\n')
    ## HVGs from coding genes
    print("Calculating HVGs from coding genes...")
    adata = sc.read("{OUTDIR}/hTCA.counts.protein_coding.h5ad".format(**locals()))
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=max(n_top_genes),
        batch_key="GSE_ID",
        subset=True
    )
    hvg_df = adata.var
    hvg_df.sort_values(by=['variances_norm'], inplace=True, ascending=False)
    hvg_df.to_csv("{OUTDIR}/HVGs.protein_coding.info.csv".format(**locals()))
    for topN in n_top_genes:
        f = open("{OUTDIR}/HVGs.{topN}.protein_coding.txt".format(**locals()), 'w')
        for i in hvg_df.name.tolist()[0:topN]:
            f.write(i+'\n')


if __name__ == '__main__':
    mouse_main()
    # human_main()
