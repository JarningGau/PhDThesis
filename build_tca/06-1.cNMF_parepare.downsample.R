library(here)
library(Seurat, lib.loc = "/home/software/R/seuratV3_env/")
library(tidyverse)
options(stringsAsFactors = FALSE)
source(here("scripts/R/utils.R"))


## mouse
if (T) {
  outpath <- here("results/06.cNMF_mouse")
  safe_mkdir(outpath)
  seu <- sceasy::convertFormat(here("results/02.preprocess_mouse/TCA_mouse.processed.h5ad"), from = "anndata", to = "seurat")
  cellmeta <- read.csv(here("results/04.scvi/TCA_mouse.cellmeta.csv"), row.names = 1)
  seu$Leiden <- factor(cellmeta[rownames(seu@meta.data), ]$leiden)

  ## from ~210,000 down to ~35,000 cells
  Idents(seu) <- seu$Leiden
  seu.ds <- subset(seu, downsample = 1000)

  ## gene sets
  expr.in.cells <- Matrix::rowSums(seu.ds[["RNA"]]@counts > 0)
  features <- readRDS(here("data/gene_set/features.mm10.ens98.rds"))
  rownames(features) <- sub("|", "-", features$uniqName, fixed = T)
  features$expr_in_cells <- expr.in.cells[rownames(features)]
  ribo.genes <- readLines(here("data/gene_set/MGI_GO_ribosomal_subunit.txt"))
  features$is_ribo <- features$Symbol %in% ribo.genes
  features$is_mito <- startsWith(features$Symbol, "mt-")
  features.use <- features %>%
    subset(!is_mito & !is_ribo) %>%
    subset(expr_in_cells >= 10)
  coding.genes <- subset(features.use, biotype == "protein_coding")

  ## save *.h5ad
  seu <- seu[rownames(features.use), ]
  sceasy::convertFormat(seu, from = "seurat", to = "anndata", main_layer = "counts",
                        outFile = file.path(outpath, "mTCA.counts.h5ad"))
  seu.coding <- seu[rownames(coding.genes), ]
  sceasy::convertFormat(seu.coding, from = "seurat", to = "anndata", main_layer = "counts",
                        outFile = file.path(outpath, "mTCA.counts.protein_coding.h5ad"))

  seu.ds <- seu.ds[rownames(features.use), ]
  sceasy::convertFormat(seu.ds, from = "seurat", to = "anndata", main_layer = "counts",
                        outFile = file.path(outpath, "mTCA_ds1k.counts.h5ad"))

  seu.ds.coding <- seu.ds[rownames(coding.genes), ]
  sceasy::convertFormat(seu.ds.coding, from = "seurat", to = "anndata", main_layer = "counts",
                        outFile = file.path(outpath, "mTCA_ds1k.counts.protein_coding.h5ad"))
}

## human
if (F) {
  outpath <- here("results/06.cNMF_human/")
  safe_mkdir(outpath)
  seu <- sceasy::convertFormat(here("results/02.preprocess_human/TCA_human.processed.h5ad"), from = "anndata", to = "seurat")
  cellmeta <- read.csv(here("results/04.scvi/TCA_human.cellmeta.csv"), row.names = 1)
  seu$Leiden <- factor(cellmeta[rownames(seu@meta.data), ]$leiden)

  ## from ~167,000 down to ~29,000 cells
  Idents(seu) <- seu$Leiden
  seu.ds <- subset(seu, downsample = 1000)

  ## gene sets
  expr.in.cells <- Matrix::rowSums(seu.ds[["RNA"]]@counts > 0)
  features <- readRDS(here("data/gene_set/features.GRCh38.ens98.rds"))
  rownames(features) <- sub("|", "-", features$uniqName, fixed = T)
  features$expr_in_cells <- expr.in.cells[rownames(features)]
  ribo.genes <- readLines(here("data/gene_set/MSigDB_GO_ribosomal_subunit.txt"))
  features$is_ribo <- features$Symbol %in% ribo.genes
  features$is_mito <- startsWith(features$Symbol, "MT-")
  features.use <- features %>%
    subset(!is_mito & !is_ribo) %>%
    subset(expr_in_cells >= 10)
  coding.genes <- subset(features.use, biotype == "protein_coding")

  ## save *.h5ad
  seu <- seu[rownames(features.use), ]
  sceasy::convertFormat(seu, from = "seurat", to = "anndata", main_layer = "counts",
                        outFile = file.path(outpath, "hTCA.counts.h5ad"))
  seu.coding <- seu[rownames(coding.genes), ]
  sceasy::convertFormat(seu.coding, from = "seurat", to = "anndata", main_layer = "counts",
                        outFile = file.path(outpath, "hTCA.counts.protein_coding.h5ad"))

  seu.ds <- seu.ds[rownames(features.use), ]
  sceasy::convertFormat(seu.ds, from = "seurat", to = "anndata", main_layer = "counts",
                        outFile = file.path(outpath, "hTCA_ds1k.counts.h5ad"))
  seu.ds.coding <- seu.ds[rownames(coding.genes), ]
  sceasy::convertFormat(seu.ds.coding, from = "seurat", to = "anndata", main_layer = "counts",
                        outFile = file.path(outpath, "hTCA_ds1k.counts.protein_coding.h5ad"))
}
