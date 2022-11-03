calculate.consensus.metrics <- function(cellmeta, cn) {
  ## cell type entropy
  clusters <- sort(unique(cellmeta[[cn]]))
  tmp.data <- table(cellmeta$Cell_type_corrected, cellmeta[[cn]]) %>% as.matrix()
  ct.entropy <- apply(tmp.data, 2, function(xx) entropy::entropy.empirical(xx, unit = "log2"))
  cellmeta$Cell_type_entropy <- ct.entropy[as.character(cellmeta[[cn]])]
  ## sample diversity
  clusters <- sort(unique(cellmeta[[cn]]))
  tmp.data.2 <- table(cellmeta$Sample_name, cellmeta[[cn]]) %>% as.matrix()
  sn.entropy <- apply(tmp.data.2, 2, function(xx) entropy::entropy.empirical(xx, unit = "log2"))
  cellmeta$Subject_diversity <- sn.entropy[as.character(cellmeta[[cn]])]
  return(cellmeta)
}

plot.consensus.metrics <- function(cellmeta, cn) {
  p1 <- ggplot(cellmeta, aes(UMAP_1, UMAP_2, color=factor(get(cn)))) +
    geom_point(size = .1, alpha=.1, show.legend = F) +
    geom_text(inherit.aes = F, data = get_label_pos(cellmeta, emb = "UMAP", group.by=cn),
              aes(x,y,label=label), size=3) +
    theme_void()

  p2 <- ggplot(cellmeta, aes(UMAP_1, UMAP_2, color=Cell_type_entropy)) +
    geom_point(size = .1, alpha=.1) +
    scale_color_viridis_c() +
    geom_text(inherit.aes = F, data = get_label_pos(cellmeta, emb = "UMAP", group.by=cn),
              aes(x,y,label=label), size=3) +
    theme_void()

  p3 <- ggplot(cellmeta, aes(UMAP_1, UMAP_2, color=Subject_diversity)) +
    geom_point(size = .1, alpha=.1) +
    scale_color_viridis_c() +
    geom_text(inherit.aes = F, data = get_label_pos(cellmeta, emb = "UMAP", group.by=cn),
              aes(x,y,label=label), size=3) +
    theme_void()

  data.plot <- cellmeta %>%
    select(which(cn == colnames(cellmeta)), Cell_type_entropy, Subject_diversity) %>%
    distinct()

  p4 <- ggplot(data.plot, aes(Cell_type_entropy, Subject_diversity)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label=get(cn))) +
    theme_classic(base_size = 15)

  p2 + p3 + p1 + p4
}
