#### Load data ####
read_solo_sparse <- function(cells, regions, mtx) {
  mtx <- Matrix::readMM(mtx)
  regions <- read.table(regions, sep = "\t", header = F)
  cells <- read.table(cells, sep = "\t", header = F)
  rownames(mtx)<- regions$V1
  colnames(mtx)<- cells$V1
  return(mtx)
}

read_solo_dataset <- function(gse_dataset, assay="raw") {
  path = here("data/matrix", gse_dataset)
  gsms <- list.dirs(path, recursive = F, full.names = F)
  pbapply::pblapply(gsms, function(gsmID) {
    cells = paste(path, gsmID, assay, "barcodes.tsv", sep = "/")
    features = paste(path, gsmID, assay, "features.tsv", sep = "/")
    mtx = paste(path, gsmID, assay, "matrix.mtx.gz", sep = "/")
    mm <- read_solo_sparse(cells = cells, regions = features, mtx = mtx)
    if (assay != "corrected") {
      colnames(mm) <- paste(gsmID, colnames(mm), sep="_")
    }
    mm
  }) %>% do.call(cbind, .)
}

#' read ambient rna ratio calculated by decontX
read_metadata_dataset <- function(gse_dataset, assay="corrected") {
  path = paste("../data/matrix", gse_dataset, sep = "/")
  gsms <- list.dirs(path, recursive = F, full.names = F)
  pbapply::pblapply(gsms, function(gsmID) {
    meta.file = paste(path, gsmID, assay, "decontX.info.tsv", sep = "/")
    read.table(meta.file, header = T, row.names = 1)
  }) %>% do.call(rbind, .)
}

load_dataset <- function(gse_dataset, species=NULL) {
  message(paste("Loading", gse_dataset, "..."))
  path = here("data/matrix", gse_dataset)
  ## load gene set
  if (species == "human") {
    ribo.genes <- readLines(here("data/gene_set/MSigDB_GO_ribosomal_subunit.txt"))
    features <- readRDS(here("data/gene_set/features.GRCh38.ens98.rds"))
  } else if (species == "mouse") {
    ribo.genes <- readLines(here("data/gene_set/MGI_GO_ribosomal_subunit.txt"))
    features <- readRDS(here("data/gene_set/features.mm10.ens98.rds"))
  } else {
    stop("Invalid 'species'. Please choose one of 'human' and 'mouse'.")
  }
  rownames(features) <- features$Ensembl.ID
  chrX.genes <- sub("|", "-", subset(features, chromosome == "X")$uniqName, fixed = T)
  chrY.genes <- sub("|", "-", subset(features, chromosome == "Y")$uniqName, fixed = T)
  ## load counts
  mm <- read_solo_dataset(gse_dataset, assay = "filtered")
  rownames(mm) <- features[rownames(mm), ]$uniqName
  ## new seurat
  seu <- CreateSeuratObject(counts = mm, project = gse_dataset)
  mito.pattern <- switch(species, "human" = "MT-", "mouse" = "mt-")
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = mito.pattern)
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, features = ribo.genes)
  seu[["percent.chrX"]] <- PercentageFeatureSet(seu, features = chrX.genes)
  seu[["percent.chrY"]] <- PercentageFeatureSet(seu, features = chrY.genes)
  celda.meta <- read_metadata_dataset(gse_dataset, assay = "corrected")
  seu[["decontX_contamination"]] <- celda.meta[rownames(seu@meta.data), ]$decontX_contamination
  # seu[["decontX_contamination"]] <- 0 ## for Nat.Com_2021_mouse_GSE148032, no decontX results.
  return(seu)
}

load_processed_dataset <- function(gse_dataset, species="human") {
  if (species == "human") {
    obj.name <- paste0(gse_dataset, ".seu.rds")
  }
  if (species == "mouse") {
    obj.name <- paste0(gse_dataset, ".seu.hc.rds")
  }
  message(paste("Loading", obj.name, "..."))
  seu <- readRDS(here(paste0("results/02.preprocess_", species), obj.name))
  seu
}

#### preprocess ####

seurat_one_step <- function(seu, outpath, do.plot=TRUE) {
  if (!dir.exists(outpath)) {
    stop(paste0(outpath, "does not exist."))
  }
  seu.project.name <- seu@project.name
  message("Norm & HGVs ...")
  seu.list <- SplitObject(seu, split.by = "orig.ident")
  for (i in 1:length(seu.list)) {
    seu.list[[i]] <- NormalizeData(seu.list[[i]], verbose = FALSE)
    seu.list[[i]] <- FindVariableFeatures(seu.list[[i]], selection.method = "vst",
                                          nfeatures = 2000, verbose = FALSE)
  }
  message("MNN correct ...")
  seu <- RunFastMNN(object.list = seu.list)
  message("UMAP ...")
  seu <- RunUMAP(seu, reduction = "mnn", dims = 1:30)
  K.auto <- signif(sqrt(ncol(seu)), 0)
  message("Cluster ...")
  seu <- seu %>%
    FindNeighbors(dims = 1:30, reduction = "mnn", k.param = K.auto, verbose = FALSE) %>%
    FindClusters(resolution = 1, verbose = FALSE)
  if (do.plot) {
    DimPlot(seu, group.by = "RNA_snn_res.1", label = T, label.size = 7) + NoLegend() + ggtitle("res=1")
    plot.file <- paste0(seu.project.name, ".cluster.png")
    ggsave(file.path(outpath, plot.file), width = 7, height = 5)
  }
  seu@project.name <- seu.project.name
  return(seu)
}

add_metadata <- function(seu) {
  gse_dataset <- seu@project.name
  path = here("data/matrix", gse_dataset)
  seu$seurat_clusters <- NULL
  gsm.meta <- read.table(paste(path, "meta.tsv", sep = "/"), sep = "\t", header = T)
  metadata.old <- seu@meta.data
  metadata.new <- metadata.old %>%
    mutate(
      GSM_ID = orig.ident,
      GSE_ID = mapvalues.2(x = GSM_ID, from = gsm.meta$GSM_ID, gsm.meta$GSE_ID),
      Sample_name = mapvalues.2(x = GSM_ID, from = gsm.meta$GSM_ID, gsm.meta$Sample_name),
      Species = mapvalues.2(x = GSM_ID, from = gsm.meta$GSM_ID, gsm.meta$Species),
      Platform = mapvalues.2(x = GSM_ID, from = gsm.meta$GSM_ID, gsm.meta$Platform),
      Cell_enrichment = mapvalues.2(x = GSM_ID, from = gsm.meta$GSM_ID, gsm.meta$Cell_enrichment),
      Condition = mapvalues.2(x = GSM_ID, from = gsm.meta$GSM_ID, gsm.meta$Condition),
      Genotype = mapvalues.2(x = GSM_ID, from = gsm.meta$GSM_ID, gsm.meta$Genotype),
      Age = mapvalues.2(x = GSM_ID, from = gsm.meta$GSM_ID, gsm.meta$Age),
      Sex = mapvalues.2(x = GSM_ID, from = gsm.meta$GSM_ID, gsm.meta$Sex),
      Cell_type_published = NA,
      Cell_type_manual = NA,
      Cell_type_predicted = NA
    )
  rownames(metadata.new) <- rownames(metadata.old)
  seu@meta.data <- metadata.new
  return(seu)
}


# library(celda)
# calculate_ambient_rna_ratio <- function(seu, do.batch=TRUE) {
#   sce <- sceasy::convertFormat(seu, from="seurat", to="sce")
#   if (do.batch) {
#     sce <- decontX(sce, z=sce$RNA_snn_res.1, batch=sce$GSM_ID)
#   } else {
#     sce <- decontX(sce, z=sce$Cell_type_published)
#   }
#   seu$decontX_contamination <- SingleCellExperiment::colData(sce)$decontX_contamination
#   seu[["correct"]] <- CreateAssayObject(counts = round(decontXcounts(sce), 0))
#   seu$nCount_correct <- NULL
#   seu$nFeature_correct <- NULL
#   return(seu)
# }


#### Doublets ####
library(DoubletFinder)

PreprocessSeurat <- function(seu, PCs=1:10) {
  message("  Normalize ...")
  seu <- NormalizeData(seu, verbose = F)
  message("  Find variable genes ...")
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = F)
  message("  Scale data ...")
  seu <- ScaleData(seu, verbose = F)
  message("  Run PCA ...")
  seu <- RunPCA(seu, verbose = F)
  message("  Find neighbors ...")
  k <- round(ncol(seu)*0.02)
  k <- ifelse(k < 20, 20, k)
  seu <- FindNeighbors(seu, dims = 1:10, reduction = "pca", k.param = k, verbose = F)
  message("  Find clusters ...")
  seu <- FindClusters(seu, resolution = 0.4, verbose = F)
  seu
}

FindOptimalpK <- function(seu, PCs=1:10, num.cores = 8) {
  sweep.res.list <- paramSweep_v3(seu, PCs = PCs, sct = F, num.cores = num.cores)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn[which.max(bcmvn$BCmetric), ]$pK
  as.numeric(as.character(pK))
}

DF <- function(seu, PCs=1:10, auto.pK=TRUE, auto.cluster=TRUE) {
  message("1. Preprocessing ...")
  # this step is for getting clusters
  if (auto.cluster) {
    seu <- PreprocessSeurat(seu, PCs = PCs)
  }

  message("2. Find optimal pK ...")
  if (auto.pK) {
    optimal.pK <- FindOptimalpK(seu, PCs = PCs, num.cores = 8)
    message(paste0("   Optimal pK = ", optimal.pK))
  } else {
    optimal.pK <- 0.09 # default
  }

  message("3. Mark doublets ...")
  homotypic.prop <- modelHomotypic(seu$seurat_clusters)
  nExp_poi <- round(0.075*ncol(seu))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  pN = 0.25 # default 0.25
  pANN_name <- paste("pANN", pN, optimal.pK, nExp_poi, sep = "_")
  DF_name <- paste("DF.classifications", pN, optimal.pK, nExp_poi, sep = "_")
  DF_name.adj <- paste("DF.classifications", pN, optimal.pK, nExp_poi.adj, sep = "_")

  seu <- doubletFinder_v3(seu, PCs = PCs, pN = pN, pK = optimal.pK, nExp = nExp_poi, reuse.pANN = F, sct = F)
  seu <- doubletFinder_v3(seu, PCs = PCs, pN = pN, pK = optimal.pK, nExp = nExp_poi.adj, reuse.pANN = pANN_name, sct = F)
  # summarize results
  seu$DF.classifications <- ifelse(seu@meta.data[[DF_name.adj]] == "Doublet", "Doublet.hc",
                                   ifelse(seu@meta.data[[DF_name]] == "Doublet", "Doublet.lc", "Singlet"))
  seu@meta.data[[pANN_name]] <- NULL
  seu@meta.data[[DF_name]] <- NULL
  seu@meta.data[[DF_name.adj]] <- NULL
  seu
}


DF_main <- function(seu, split.by=NULL) {
  if (is.null(split.by)) {
    seu.list <- list(seu)
  } else {
    seu.list <- SplitObject(seu, split.by = split.by)
  }
  names(seu.list) <- NULL
  metadata.new <- pbapply::pblapply(seu.list, function(xx) {
    seu.tmp <- DF(xx, PCs=1:10, auto.pK=TRUE, auto.cluster=TRUE)
    seu.tmp@meta.data
  }) %>% do.call(rbind, .)
  # rownames(metadata.new) <- sub("^GSM[0-9]+\\.", "", rownames(metadata.new))
  seu$DF.classifications <- metadata.new[rownames(seu@meta.data), ]$DF.classifications
  seu
}
