library(tidyverse)
library(Seurat, lib.loc = "/home/software/R/seuratV3_env/")
library(here)
source(here("scripts/R/utils.R"))

########################### Mouse ###################################
if (T){
  inpath <- here("results/05.metacell_mouse")
  outpath <- here("results/06.cNMF_mouse")
  datasets <- list.dirs("../data/matrix", recursive = F, full.names = F)
  datasets <- datasets[grep("*mouse*", datasets)]
  safe_mkdir(outpath)

  ## merge all metacells
  metacells_mat <- lapply(datasets, function(gse_dataset) {
    obj.names <- paste0(gse_dataset, ".metacell.rds")
    metacells <- readRDS(file.path(inpath, obj.names))
    metacells$mat
  }) %>% do.call(cbind, .)
  ## gene sets
  expr.in.cells <- rowSums(metacells_mat > 0)
  features <- readRDS(here("data/gene_set/features.mm10.ens98.rds"))
  rownames(features) <- sub("|", "-", features$uniqName, fixed = T)
  features$expr_in_cells <- expr.in.cells[rownames(features)]
  ribo.genes <- readLines(here("data/gene_set/MGI_GO_ribosomal_subunit.txt"))
  features$is_ribo <- features$Symbol %in% ribo.genes
  features$is_mito <- startsWith(features$Symbol, "mt-")
  features.use <- features %>%
    subset(!is_mito & !is_ribo) %>%
    subset(expr_in_cells >= 3)
  ## save in h5ad format
  metacells_mat <- metacells_mat[rownames(features.use), ]
  seu <- CreateSeuratObject(counts = metacells_mat)
  sceasy::convertFormat(seu, from = "seurat", to = "anndata", main_layer = "counts",
                        outFile = file.path(outpath, "mTCA_metacell.counts.h5ad"))
  coding.genes <- subset(features.use, biotype == "protein_coding")
  seu.coding <- seu[rownames(coding.genes), ]
  sceasy::convertFormat(seu.coding, from = "seurat", to = "anndata", main_layer = "counts",
                        outFile = file.path(outpath, "mTCA_metacell.counts.protein_coding.h5ad"))
}

########################### Human ###################################
if (T){
  inpath <- here("results/05.metacell_human")
  outpath <- here("results/06.cNMF_human")
  datasets <- list.dirs("../data/matrix", recursive = F, full.names = F)
  datasets <- datasets[grep("*human*", datasets)]
  safe_mkdir(outpath)

  ## merge all metacells
  metacells_mat <- lapply(datasets, function(gse_dataset) {
    obj.names <- paste0(gse_dataset, ".metacell.rds")
    metacells <- readRDS(file.path(inpath, obj.names))
    metacells$mat
  }) %>% do.call(cbind, .)
  ## gene sets
  expr.in.cells <- rowSums(metacells_mat > 0)
  features <- readRDS(here("data/gene_set/features.GRCh38.ens98.rds"))
  rownames(features) <- sub("|", "-", features$uniqName, fixed = T)
  features$expr_in_cells <- expr.in.cells[rownames(features)]
  ribo.genes <- readLines(here("data/gene_set/MSigDB_GO_ribosomal_subunit.txt"))
  features$is_ribo <- features$Symbol %in% ribo.genes
  features$is_mito <- startsWith(features$Symbol, "MT-")
  features.use <- features %>%
    subset(!is_mito & !is_ribo) %>%
    subset(expr_in_cells >= 3)
  ## save in h5ad format
  metacells_mat <- metacells_mat[rownames(features.use), ]
  seu <- CreateSeuratObject(counts = metacells_mat)
  sceasy::convertFormat(seu, from = "seurat", to = "anndata", main_layer = "counts",
                        outFile = file.path(outpath, "hTCA_metacell.counts.h5ad"))
  coding.genes <- subset(features.use, biotype == "protein_coding")
  seu.coding <- seu[rownames(coding.genes), ]
  sceasy::convertFormat(seu.coding, from = "seurat", to = "anndata", main_layer = "counts",
                        outFile = file.path(outpath, "hTCA_metacell.counts.protein_coding.h5ad"))
}
