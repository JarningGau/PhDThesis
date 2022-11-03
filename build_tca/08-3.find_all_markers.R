library(here)
library(Seurat, lib.loc = "/home/software/R/seuratV3_env/")
source(here("scripts/R/mcFindAllMarkers.R"))

OUTDIR = here("results/08.iterative_clustering/")

## human
if (F) {
  seu <- readRDS(here("results/08.iterative_clustering/hTCA.core.seurat.rds"))
  Idents(seu) <- factor(seu$Cell_type_final)
  seu.ds <- subset(seu, downsample = 200)
  seu.ds <- NormalizeData(seu.ds)
  all.markers <- mcFindAllMarkers(seu.ds)
  readr::write_tsv(all.markers, file = file.path(OUTDIR, "celltype_markers.hTCA.tsv"))

  # FeaturePlot(seu.ds, features = c("TEX19"))
  saveRDS(seu.ds, file.path(OUTDIR, "hTCA.seurat.sub200.rds"))
}

## mouse
if (T) {
  seu <- readRDS(here("results/08.iterative_clustering/mTCA.core.seurat.rds"))
  Idents(seu) <- factor(seu$Cell_type_final)
  seu.ds <- subset(seu, downsample = 200)
  seu.ds <- NormalizeData(seu.ds)
  all.markers <- mcFindAllMarkers(seu.ds)
  readr::write_tsv(all.markers, file = file.path(OUTDIR, "celltype_markers.mTCA.tsv"))

  FeaturePlot(seu.ds, features = c("Prdm14"))
  saveRDS(seu.ds, file.path(OUTDIR, "mTCA.seurat.sub200.rds"))
}

