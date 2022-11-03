library(Seurat, lib.loc = "/home/software/R/seuratV3_env/")
library(tidyverse)
library(magrittr)
library(here)
source(here("scripts/R/utils.R"))

makeMetaCells <- function(seu, min.cells=10, reduction="umap", dims=1:2, k.param=10, cores=20) {
  seu <- seu %>%
    FindNeighbors(reduction = reduction, dims=dims, k.param=k.param) %>%
    FindClusters(res=50)
  metadata <- seu@meta.data
  metadata$METACELL_ID <- factor(metadata$seurat_clusters)
  dge_mat <- seu[["RNA"]]@counts

  dge_mat_mc <- parallel::mclapply(levels(metadata$METACELL_ID), function(xx) {
    cells <- rownames(subset(metadata, METACELL_ID == xx))
    Matrix::rowSums(dge_mat[, cells])
  }, mc.cores = cores)
  dge_mat_mc <- do.call(cbind, dge_mat_mc)

  metacell_metadata <-
    metadata[["METACELL_ID"]] %>%
    table() %>%
    as.data.frame() %>%
    set_colnames(c("METACELL_ID", "CELL_COUNT"))
  rownames(metacell_metadata) <- metacell_metadata[["METACELL_ID"]]
  kept.cells <- subset(metacell_metadata, CELL_COUNT>=min.cells)[["METACELL_ID"]]
  metacells <- list(
    mat = dge_mat_mc[, kept.cells],
    metadata = metacell_metadata[kept.cells, ]
  )
  colnames(metacells$mat) <- paste0(seu@project.name, ".METACELL_", kept.cells)
  rownames(metacells$metadata) <- colnames(metacells$mat)
  metacells
}

sumN <- function(seu, field="Cell_type_published", N=20, min.cells=10, seed=1024, cores=20) {
  dge_mat <- seu[["RNA"]]@counts
  metadata <- seu@meta.data
  cell.types <- unique(metadata[[field]])

  ## calculate pool.IDs
  set.seed(seed)
  metadata <- lapply(cell.types, function(cc) {
    new.metadata <- subset(metadata, get(field) == cc)
    new.metadata <- new.metadata[permute::shuffle(rownames(new.metadata)), ]
    new.metadata$PoolID <- floor(0:(nrow(new.metadata)-1) / N)
    new.metadata$METACELL_ID <- paste(new.metadata[[field]], new.metadata$PoolID, sep = ".")
    new.metadata
  }) %>% do.call(rbind, .)
  metadata$METACELL_ID %<>% as.factor()

  ## Pool counts
  dge_mat_mc <- parallel::mclapply(levels(metadata$METACELL_ID), function(xx) {
    cells <- rownames(subset(metadata, METACELL_ID == xx))
    Matrix::rowSums(dge_mat[, cells])
  }, mc.cores = cores)
  dge_mat_mc <- do.call(cbind, dge_mat_mc)

  ## Metacell metadata
  metacell_metadata <-
    metadata[["METACELL_ID"]] %>%
    table() %>%
    as.data.frame() %>%
    set_colnames(c("METACELL_ID", "CELL_COUNT"))
  rownames(metacell_metadata) <- metacell_metadata[["METACELL_ID"]]
  kept.cells <- subset(metacell_metadata, CELL_COUNT>=min.cells)[["METACELL_ID"]]
  metacells <- list(
    mat = dge_mat_mc[, kept.cells],
    metadata = metacell_metadata[kept.cells, ]
  )
  colnames(metacells$mat) <- paste0(seu@project.name, ".METACELL_", kept.cells)
  rownames(metacells$metadata) <- colnames(metacells$mat)
  metacells
}

MC_main <- function(gse_dataset, ...) {
  obj.name <- paste0(gse_dataset, ".seu.hc.rds")
  message(paste("Loading", obj.name, "..."))
  seu <- readRDS(file.path(inpath, obj.name))

  message("Compute metacell assignments ...")
  metacells <- makeMetaCells(seu, ...)

  message("Saving metacells ...")
  out.obj <- paste0(gse_dataset, ".metacell.rds")
  saveRDS(metacells, file.path(outpath, out.obj))
  invisible()
}

sumN_main <- function(gse_dataset, ...) {
  obj.name <- paste0(gse_dataset, ".seu.hc.rds")
  message(paste("Loading", obj.name, "..."))
  seu <- readRDS(file.path(inpath, obj.name))

  message("Compute metacell assignments ...")
  metacells <- sumN(seu, ...)

  message("Saving metacells ...")
  out.obj <- paste0(gse_dataset, ".metacell.rds")
  saveRDS(metacells, file.path(outpath, out.obj))
  invisible()
}

################################################################################
## Mouse
if (TRUE) {
  inpath <- "../results/02.preprocess_mouse"
  outpath <- "../results/05.metacell_mouse"
  datasets <- list.dirs("../data/matrix", recursive = F, full.names = F)
  datasets <- datasets[grep("*mouse*", datasets)]
  safe_mkdir(outpath)

  labeled.datasets <- c("Cell.Res_2018_mouse_GSE107644", "Nat.Com_2021_mouse_GSE148032") ## smart-seq2
  for (gse_dataset in datasets) {
    if (gse_dataset %in% labeled.datasets) {
      sumN_main(gse_dataset)
    } else {
      MC_main(gse_dataset)
    }
  }
}


################################################################################
## Human
if (FALSE) {
  inpath <- here("results/02.preprocess_human")
  outpath <- here("results/05.metacell_human")
  datasets <- list.dirs(here("data/matrix"), recursive = F, full.names = F)
  datasets <- datasets[grep("*human*", datasets)]
  safe_mkdir(outpath)
  ## make metacells
  labeled.datasets <- c("Cell.Stem.Cell_2017_human_GSE86146", "Cell.Stem.Cell_2017_human_GSE86146") ## smart-seq2
  for (gse_dataset in datasets) {
    if (gse_dataset %in% labeled.datasets) {
      sumN_main(gse_dataset)
    } else {
      MC_main(gse_dataset)
    }
  }
}


