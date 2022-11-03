library(here)

#### load functions ####

# load the kallisto count matrix
read_kb_sparse <- function(sample){
  cells.file <- here("data/benchmark/kb_matrix", sample, "cells_x_genes.barcodes.txt.gz")
  regions.file <- here("data/benchmark/kb_matrix", sample, "cells_x_genes.genes.txt.gz")
  mtx.file <- here("data/benchmark/kb_matrix", sample, "cells_x_genes.mtx.gz")
  mtx <- Matrix::readMM(mtx.file)
  # the sparse matrix with rows are cells and columns are peaks/features
  mtx<- Matrix::t(mtx)
  regions <- readLines(regions.file)
  cells <- readLines(cells.file)
  # remove version of ensembl gene id: ENSG00000243485.5 -> ENSG00000243485
  rownames(mtx)<- sub("\\.[0-9]+", "", regions)
  # cellranger add -1 to the cell barcode, I add it for later compare with cellranger output
  colnames(mtx)<- paste0(cells, "-1")
  return(mtx)
}


# load the cellranger count matrix
read_cr_sparse <- function(sample){
  cells.file <- here("data/benchmark/10X_matrix", sample, "barcodes.tsv.gz")
  regions.file <- here("data/benchmark/10X_matrix", sample, "features.tsv.gz")
  mtx.file <- here("data/benchmark/10X_matrix", sample, "matrix.mtx.gz")

  mtx <- Matrix::readMM(mtx.file)
  regions <- read.table(regions.file, sep = "\t", header = FALSE)
  cells <- readLines(cells.file)
  rownames(mtx)<- regions$V1
  colnames(mtx)<- cells
  return(mtx)
}

# load the STARsolo count matrix
read_solo_sparse <- function(sample) {
  cells.file <- here("data/benchmark/solo_matrix", sample, "barcodes.tsv.gz")
  regions.file <- here("data/benchmark/solo_matrix", sample, "features.tsv.gz")
  mtx.file <- here("data/benchmark/solo_matrix", sample, "matrix.mtx.gz")

  mtx <- Matrix::readMM(mtx.file)
  regions <- read.table(regions.file, sep = "\t", header = FALSE)
  cells <- readLines(cells.file)
  rownames(mtx)<- regions$V1
  # cellranger add -1 to the cell barcode, I add it for later compare with cellranger output
  colnames(mtx)<- paste0(cells, "-1")
  return(mtx)
}


#### compare cell barcodes ####

kb_cr_barcodes <- function(sample) {
  kb_path <- here("data/benchmark/kb_matrix", sample, "cells_x_genes.barcodes.txt.gz")
  cr_path <- here("data/benchmark/10X_matrix", sample, "barcodes.tsv.gz")
  kb_bc <- readLines(kb_path)
  kb_bc <- paste0(kb_bc, "-1")
  cr_bc <- readLines(cr_path)
  ab_bc <- intersect(kb_bc, cr_bc) |> length()
  a_bc <- setdiff(kb_bc, cr_bc) |> length()
  b_bc <- setdiff(cr_bc, kb_bc) |> length()
  data.frame(
    sample = sample,
    set_name = c("kb_uniq", "common", "cr_uniq"),
    set_size = c(a_bc, ab_bc, b_bc),
    total = sum(c(a_bc, ab_bc, b_bc)),
    set_perc = c(a_bc, ab_bc, b_bc) / sum(c(ab_bc, b_bc)) * 100
  )
}

solo_cr_barcodes <- function(sample) {
  solo_path <- here("data/benchmark/solo_matrix", sample, "barcodes.tsv.gz")
  cr_path <- here("data/benchmark/10X_matrix", sample, "barcodes.tsv.gz")
  solo_path <- readLines(solo_path)
  solo_path <- paste0(solo_path, "-1")
  cr_bc <- readLines(cr_path)
  ab_bc <- intersect(solo_path, cr_bc) %>% length()
  a_bc <- setdiff(solo_path, cr_bc) %>% length()
  b_bc <- setdiff(cr_bc, solo_path) %>% length()
  data.frame(
    sample = sample,
    set_name = c("solo_uniq", "common", "cr_uniq"),
    set_size = c(a_bc, ab_bc, b_bc),
    total = sum(c(a_bc, ab_bc, b_bc)),
    set_perc = c(a_bc, ab_bc, b_bc) / sum(c(ab_bc, b_bc)) * 100
  )
}

#### compare UMI ####
## total UMI
kb_cr_umi <- function(sample, kb_mm, cr_mm) {
  bc.used <- intersect(colnames(kb_mm), colnames(cr_mm))
  kb_mm <- kb_mm[, bc.used]
  cr_mm <- cr_mm[, bc.used]
  a_umi <- Matrix::colSums(kb_mm)
  b_umi <- Matrix::colSums(cr_mm)
  data.frame(
    sample = sample,
    barcode = bc.used,
    kb = a_umi,
    cr = b_umi
  )
}

solo_cr_umi <- function(sample, solo_mm, cr_mm) {
  bc.used <- intersect(colnames(solo_mm), colnames(cr_mm))
  solo_mm <- solo_mm[, bc.used]
  cr_mm <- cr_mm[, bc.used]
  a_umi <- Matrix::colSums(solo_mm)
  b_umi <- Matrix::colSums(cr_mm)
  data.frame(
    sample = sample,
    barcode = bc.used,
    solo = a_umi,
    cr = b_umi
  )
}

## gene level UMI
# calculating PCCs of the profile of the same cells from different pipe
# ref to: https://divingintogeneticsandgenomics.rbind.io/post/compare-kallisto-bustools-and-cellranger-for-single-nuclei-seqencing-data/
pair_column_cor <- function(x, y) {
  bc.used <- intersect(colnames(x), colnames(y))
  gene.used <- intersect(rownames(x), rownames(y))
  x <- x[gene.used, bc.used] %>% as.matrix()
  y <- y[gene.used, bc.used] %>% as.matrix()

  sqr = function(x) x*x
  if(!is.matrix(x)||!is.matrix(y)||any(dim(x)!=dim(y)))
    stop("Please supply two matrices of equal size.")
  x   = sweep(x, 2, colMeans(x))
  y   = sweep(y, 2, colMeans(y))
  cor = colSums(x*y) /  sqrt(colSums(sqr(x))*colSums(sqr(y)))
  return(cor)
}

