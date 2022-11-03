library(ggplot2)

theme_set(
  theme_bw(base_size = 15) +
    theme(axis.line = element_line(color="black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = .5, face = "bold"),
          panel.grid = element_blank(),
          legend.background = element_rect(fill=alpha('white', 0))
    )
)

