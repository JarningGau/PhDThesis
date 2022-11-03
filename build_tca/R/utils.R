#' @description Make directory if it not exists.
safe_mkdir <- function(dn) {
  if (!dir.exists(dn)) {
    dir.create(dn, recursive = T)
  }
}


#' @description Negate operation of '%in%'.
`%notin%` <- Negate(`%in%`)


#' @description Item in 'x' that match items 'from' will be replaced by items in 'to', matched by position.
#' If items in 'x' were not in 'from', NA would replace at the matched position.
#' @param x <vector> The vector to modifiy
#' @param from <vector> keys
#' @param to <vector> values
mapvalues.2 <- function(x, from ,to) {
  z <- plyr::mapvalues(x, from, to, warn_missing = F)
  ifelse(x %in% from, z, NA)
}

#' .replace function like pandas.DataFrame.replace()
.replace <- function(col.dat, h) {
  h2 <- unique(col.dat)
  h2 <- setdiff(h2, names(h))
  names(h2) <- h2
  h <- c(h, h2)
  h[col.dat]
}

#' calculating the position of cluster labels
get_label_pos <- function(data, emb = "tSNE", group.by="ClusterID") {
  new.data <- data[, c(paste(emb, 1:2, sep = "_"), group.by)]
  colnames(new.data) <- c("x","y","cluster")
  clusters <- names(table(new.data$cluster))
  new.pos <- lapply(clusters, function(i) {
    tmp.data = subset(new.data, cluster == i)
    data.frame(
      x = median(tmp.data$x),
      y = median(tmp.data$y),
      label = i)
  })
  do.call(rbind, new.pos)
}


#' @description downsample the data.frame by group
downsample <- function(data, batch, min.size=2000, seed=1024) {
  set.seed(seed)
  cates <- unique(data[[batch]])
  cates <- cates[!is.na(cates)]
  rns.ds <- lapply(cates, function(xx) {
    rns <- rownames(data[data[[batch]] == xx & !is.na(data[[batch]]), ])
    nr <- length(rns)
    if(nr > min.size) {
      rns <- sample(rns, size = min.size, replace = F)
    }
    rns
  }) %>% do.call(c, .)
  data[rns.ds, ]
}
