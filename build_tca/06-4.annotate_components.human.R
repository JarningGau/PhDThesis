library(here)
library(readr)
library(magrittr)
library(ggplot2)
options(stringsAsFactors = FALSE)
source(here("scripts/R/plotRegulonRank.R"))
source(here("scripts/R/utils.R"))

K = 80
run.name = paste0("K", K)
INDIR = here("results/06.cNMF_human/hTCA_mc.all.cnmf.hvg2000")


#### geneset analysis ####
top.genes <- readr::read_tsv(file.path(INDIR, run.name, paste0("top100_genes.k_",K,".dt_0_2.txt")), show_col_types = FALSE)
term2gene <- readRDS(here("data/gene_set/HCL_TableS2_full.human_symbol.rds"))
head(term2gene)

e.res <- parallel::mclapply(seq_along(top.genes), function(ii) {
  y <- clusterProfiler::enricher(top.genes[[ii]], TERM2GENE=term2gene, minGSSize=20)
  y %<>% as.data.frame()
  if (nrow(y)) {
    y$program <- names(top.genes)[ii]
  }
  return(y)
}, mc.cores = 20) %>% dplyr::bind_rows(.)

readr::write_tsv(e.res, file.path(INDIR, run.name, paste0("top100_genes.k_",K,".gsea.txt")))

#### regulon-specific-score (RSS) ####
## load dataset
cellmeta <- read.csv(here("results/04.scvi/TCA_human.cellmeta.csv"), row.names = 1)
aucMat <- readRDS(file.path(INDIR, run.name, "hTCA.AUCell.rds"))
## cell type matrix
cell.types <- names(table(cellmeta$leiden))
ctMat <- lapply(cell.types, function(i) {
  as.numeric(cellmeta$leiden == i)
}) %>% do.call(cbind, .)
colnames(ctMat) <- cell.types
rownames(ctMat) <- rownames(cellmeta)
## rss matrix
rssMat <- parallel::mclapply(colnames(aucMat), function(i) {
  sapply(colnames(ctMat), function(j) {
    1 - philentropy::JSD(rbind(aucMat[, i], ctMat[, j]), unit = 'log2', est.prob = "empirical")
  })
}, mc.cores = 20)
rssMat <- do.call(rbind, rssMat)
rownames(rssMat) <- colnames(aucMat)
colnames(rssMat) <- colnames(ctMat)

## plots
plot.list <- lapply(colnames(rssMat), function(xx) PlotRegulonRank(rssMat, xx))
plotDIR <- file.path(INDIR, run.name, "RSS_figure")
safe_mkdir(plotDIR)

for (i in seq_along(plot.list)) {
  print(plot.list[[i]])
  plot.file <- paste0("Leiden_", colnames(rssMat)[i], ".png")
  ggsave(file.path(plotDIR, plot.file), width = 3, height = 5, units = "in", dpi = 300)
}

