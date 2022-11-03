library(celda)
library(SingleCellExperiment)
library(Seurat, lib.loc = "/home/software/R/seuratV3_env/")
library(pbapply)
library(here)

# https://bioconductor.org/packages/release/bioc/manuals/celda/man/celda.pdf

## functions
#'@title read_solo_sparse
read_solo_sparse <- function(cells, regions, mtx) {
  mtx <- Matrix::readMM(mtx)
  regions <- read.table(regions, sep = "\t", header = F)
  cells <- read.table(cells, sep = "\t", header = F)
  rownames(mtx)<- regions$V1
  colnames(mtx)<- cells$V1
  return(mtx)
}


saveMM <- function(mtx, outpath) {
  if (!dir.exists(outpath)) {
    dir.create(outpath, recursive = T)
  }
  Matrix::writeMM(mtx, paste(outpath, "matrix.mtx", sep = "/"))
  writeLines(rownames(mtx), paste(outpath, "features.tsv", sep = "/"))
  writeLines(colnames(mtx), paste(outpath, "barcodes.tsv", sep = "/"))
  R.utils::gzip(paste(outpath, "matrix.mtx", sep = "/"))
}


#' @title remove ambient RNAs
#' @description only for droplet-based methods, like 10x, dropseq, and inDrops.
remove_ambient_RNAs <- function(gse_dataset) {
  path = here("data/matrix", gse_dataset)
  gsms <- list.dirs(path, recursive = F, full.names = F)

  pblapply(gsms, function(gsmID) {
    barcodes = paste(path, gsmID, "raw", "barcodes.tsv", sep = "/")
    features = paste(path, gsmID, "raw", "features.tsv", sep = "/")
    mtx = paste(path, gsmID, "raw", "matrix.mtx.gz", sep = "/")
    message(paste("\nloading", gsmID, "..."))
    raw.mtx <- read_solo_sparse(barcodes, features, mtx)
    # divide cells and droplets
    barcodes <- colnames(raw.mtx)
    cells = paste(path, gsmID, "filtered", "barcodes.tsv", sep = "/")
    cells <- read.table(cells, sep = "\t", header = F)$V1
    droplets <- setdiff(barcodes, cells)
    cells.mtx <- raw.mtx[, cells]
    bg.mtx <- raw.mtx[, droplets]
    colnames(cells.mtx) <- paste(gsmID, colnames(cells.mtx), sep="_")
    colnames(bg.mtx) <- paste(gsmID, colnames(bg.mtx), sep="_")
    sce <- SingleCellExperiment(list(counts = cells.mtx))
    sce.bg <- SingleCellExperiment(list(counts = bg.mtx))
    # celda::decontX for ambient RNA removal
    sce <- decontX(sce, background = sce.bg)
    # save corrected counts
    outpath <- file.path(path, gsmID, "corrected")
    message(paste("saving corrected counts to", outpath))
    saveMM(mtx = round(decontXcounts(sce), 0), outpath)
    write.table(x = colData(sce),
                file = file.path(outpath, "decontX.info.tsv"),
                sep = "\t", quote = F)
  })
}

## main
inpath = here("data/matrix/")
datasets <- list.dirs(inpath, recursive = F, full.names = F)

## human
remove_ambient_RNAs("Cell.Rep_2018_human_GSE109033")
remove_ambient_RNAs("Cell.Rep_2019_human_GSE124263")
remove_ambient_RNAs("Cell.Res_2018_human_GSE120508")
remove_ambient_RNAs("Cell.Stem.Cell_2017_human_GSE86146")
remove_ambient_RNAs("Cell.Stem.Cell_2018_human_GSE106487")
remove_ambient_RNAs("Cell.Stem.Cell_2020_human_GSE134144")
remove_ambient_RNAs("Nat.Com_2020_human_GSE149512")

## mouse
remove_ambient_RNAs("Cell.Rep_2018_mouse_GSE109037")
remove_ambient_RNAs("Cell.Res_2018_mouse_GSE107644")
remove_ambient_RNAs("Dev_2020_mouse_GSE130593")
remove_ambient_RNAs("Dev.Cell_2018_mouse_GSE112393")
remove_ambient_RNAs("Elife_2019_mouse_GSE113293")
remove_ambient_RNAs("Nat.Com_2019_mouse_E-MTAB-6946")
remove_ambient_RNAs("Nat.Com_2019_mouse_GSE124904")
remove_ambient_RNAs("Plos.Gen_2019_mouse_GSE121904")


