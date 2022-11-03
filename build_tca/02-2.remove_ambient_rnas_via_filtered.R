library(celda)
library(SingleCellExperiment)
library(Seurat, lib.loc = "/home/software/R/seuratV3_env/")
library(pbapply)
library(here)

# https://bioconductor.org/packages/release/bioc/manuals/celda/man/celda.pdf

## functions
remove_ambient_RNAs <- function(gse_dataset, do.batch=TRUE) {
  message(paste("Loading", gse_dataset, "..."))
  objname <- paste(inpath, paste0(gse_dataset, ".seu.rds"), sep="/")
  seu <- readRDS(objname)
  sce <- sceasy::convertFormat(seu, from="seurat", to="sce")
  if (do.batch) {
    sce <- decontX(sce, z=sce$RNA_snn_res.1, batch=sce$GSM_ID)
  } else {
    sce <- decontX(sce, z=sce$Cell_type_published)
  }
  seu.backup <- seu
  seu[["RNA"]]@counts <- decontXcounts(sce)
  seu <- NormalizeData(seu)
  seu[["RNA"]]@counts <- round(seu[["RNA"]]@counts, 0)
  seu$decontX_contamination <- colData(sce)$decontX_contamination

  objname <- paste(outpath, paste0(gse_dataset, ".h5ad"), sep="/")
  sceasy::convertFormat(seu, from = "seurat", to = "anndata", main_layer = "counts", outFile = objname)
}

## human
inpath <- here("results/02.preprocess_human/")
outpath <- here("results/03.decontX_human/")
safe_mkdir(outpath)

remove_ambient_RNAs("Cell.Rep_2018_human_GSE109037")
remove_ambient_RNAs("Cell.Rep_2019_human_GSE124263")
remove_ambient_RNAs("Cell.Res_2018_human_GSE120508")
remove_ambient_RNAs("Cell.Stem.Cell_2017_human_GSE86146", do.batch = FALSE)
remove_ambient_RNAs("Cell.Stem.Cell_2018_human_GSE106487", do.batch = FALSE)
remove_ambient_RNAs("Cell.Stem.Cell_2020_human_GSE134144")
remove_ambient_RNAs("Nat.Com_2020_human_GSE149512")

## mouse
inpath <- here("results/02.preprocess_mouse/")
outpath <- here("results/03.decontX_mouse/")
safe_mkdir(outpath)

remove_ambient_RNAs("Cell.Rep_2018_mouse_GSE109033")
remove_ambient_RNAs("Dev_2020_mouse_GSE130593")
remove_ambient_RNAs("Dev.Cell_2018_mouse_GSE112393")
remove_ambient_RNAs("Elife_2019_mouse_GSE113293")
remove_ambient_RNAs("Nat.Com_2019_mouse_E-MTAB-6946")
remove_ambient_RNAs("Nat.Com_2019_mouse_GSE124904")
remove_ambient_RNAs("Plos.Gen_2019_mouse_GSE121904")
remove_ambient_RNAs("Cell.Res_2018_mouse_GSE107644", do.batch = FALSE)
remove_ambient_RNAs("Nat.Com_2021_mouse_GSE148032", do.batch = FALSE)
