library(here)
library(tidyverse)
library(scProject)
source(here("scripts/R/00.theme.R"))
source(here("scripts/R/utils.R"))

INFILE = here("results/06.cNMF_human/hTCA_metacell.counts.h5ad")
HVGFILE = here("results/06.cNMF_human/HVGs.2000.all.txt")
OUTPATH = here("results/06.cNMF_human/hTCA_mc.all.cnmf.hvg2000")

## step1
if (F) {
  FindOptimalK(counts.fn = INFILE,
               components = seq(20,100,5),
               out.path = OUTPATH,
               run.name = "K_20_100_by5",
               n.iter = 20,
               genes.fn = HVGFILE,
               cores = 40)
}

## step2
if (T) {
  K = 80
  run.name = paste0("K",K)
  RunCNMF(counts.fn = INFILE,
          K = K,
          out.path = OUTPATH,
          run.name = run.name,
          n.iter = 100,
          genes.fn = HVGFILE,
          cores = 40,
          n.top.genes = 100,
          local.density.cutoff = 0.2,
          show.clustering = T)


  ## compute gene module score
  top.genes <- read.table(file.path(OUTPATH, run.name, paste0("top100_genes.k_",K,".dt_0_2.txt")), header = T)
  seu <- sceasy::convertFormat(here("results/02.preprocess_human/TCA_human.processed.h5ad"),
                               from = "anndata", to = "seurat")
  ras.mat <- ComputeModuleScore(seu[["RNA"]]@counts, gene.sets = top.genes, cores = 20)
  newdata <- Matrix::t(ras.mat) %>% as.data.frame()
  saveRDS(newdata, file.path(OUTPATH, run.name, "hTCA.AUCell.rds"))

  ## plots
  newdata <- readRDS(file.path(OUTPATH, run.name, "hTCA.AUCell.rds"))
  umap.emb <- read.csv(here("results/04.scvi/TCA_human.emb_umap.csv"), row.names = 1)
  sel.cells <- sample(rownames(newdata), size = 1e4)
  data.plot <- cbind(newdata[sel.cells, ], umap.emb[sel.cells, ])
  data.plot <- data.plot %>%
    pivot_longer(cols = 1:ncol(newdata), names_to = "component", values_to = "score")

  plot.list <- data.plot %>%
    group_split(component) %>%
    map(
      ~ggplot(., aes(UMAP_1, UMAP_2, color = score)) +
        geom_point(size = .5) +
        scale_color_viridis_c() +
        facet_grid(~ component, labeller = function(x) label_value(x, multi_line = FALSE))
    )

  safe_mkdir(file.path(OUTPATH, run.name, "figure"))
  prog.names <- sort(unique(data.plot$component))
  for (i in seq_along(plot.list)) {
    print(plot.list[[i]])
    plot.file <- paste0(prog.names[i], ".png")
    ggsave(file.path(OUTPATH, run.name, "figure", plot.file), width = 6, height = 5, units = "in", dpi = 300)
  }
}

