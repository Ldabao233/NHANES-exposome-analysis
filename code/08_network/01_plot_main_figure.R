# 04_plot_main_Figure7.R
# Style follows the user's previous enrichment.R:
# clean theme_bw/classic, blue-to-red enrichment scale, horizontal bars,
# ggalluvial/networkD3-style Sankey, and compact publication-ready PDF export.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggalluvial)
library(patchwork)
library(igraph)

repo_root <- normalizePath(Sys.getenv("NHANES_REPO", "."), mustWork = TRUE)
setwd(file.path(repo_root, "results", "network"))
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

save_pub_r <- function(plot, filename, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = width_in, height = height_in)
  print(plot)
  dev.off()
  grDevices::pdf(
    paste0(filename, ".pdf"), width = width_in, height = height_in,
    family = "Helvetica", useDingbats = FALSE
  )
  print(plot)
  dev.off()
  ragg::agg_tiff(
    paste0(filename, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, compression = "lzw"
  )
  print(plot)
  dev.off()
  ragg::agg_png(
    paste0(filename, ".png"), width = width_in, height = height_in,
    units = "in", res = 180
  )
  print(plot)
  dev.off()
}

# ==========================================================
# A. Chemical-associated targets ∩ CVD-associated genes
# ==========================================================
chemical_targets <- read_csv(file.path(repo_root, "data", "processed", "network", "chemical_targets.csv"),
                             show_col_types = FALSE)$Gene
cvd_genes <- read_csv(file.path(repo_root, "data", "processed", "network", "cvd_genes_unique.csv"),
                      show_col_types = FALSE)$Gene

venn_list <- list(
  "LE8-associated chemical targets" = unique(chemical_targets),
  "CVD-associated genes" = unique(cvd_genes)
)

# Native ggplot2 two-set Venn, avoiding an optional package dependency.
venn_circle <- function(cx, cy, r, label, fill) {
  theta <- seq(0, 2 * pi, length.out = 361)
  tibble(
    x = cx + r * cos(theta),
    y = cy + r * sin(theta),
    Set = label,
    Fill = fill
  )
}

only_a <- length(setdiff(venn_list[[1]], venn_list[[2]]))
only_b <- length(setdiff(venn_list[[2]], venn_list[[1]]))
overlap <- length(intersect(venn_list[[1]], venn_list[[2]]))

venn_shapes <- bind_rows(
  venn_circle(-0.38, 0, 0.82, "Targets", "#8ba1cb"),
  venn_circle(0.38, 0, 0.82, "CVD", "#c7d2e6")
)

pA <- ggplot() +
  geom_polygon(
    data = venn_shapes,
    aes(x = x, y = y, group = Set, fill = Fill),
    alpha = 0.62,
    colour = "black",
    linewidth = 0.45
  ) +
  scale_fill_identity() +
  annotate("text", x = -0.78, y = 1.02,
           label = "Chemical-associated\ntargets", size = 3.1, fontface = "bold") +
  annotate("text", x = 0.78, y = 1.02,
           label = "CVD-associated\ngenes", size = 3.1, fontface = "bold") +
  annotate("text", x = -0.59, y = 0, label = only_a, size = 4.5, fontface = "bold") +
  annotate("text", x = 0, y = 0, label = overlap, size = 4.8, fontface = "bold") +
  annotate("text", x = 0.59, y = 0, label = only_b, size = 4.5, fontface = "bold") +
  coord_fixed(xlim = c(-1.35, 1.35), ylim = c(-1.25, 1.28), expand = FALSE) +
  theme_void()

ggsave("figures/Figure7A_Targets_CVD_Venn.pdf",
       pA, width = 6.5, height = 5.5)

# ==========================================================
# B. GO-BP enrichment
# ==========================================================
go <- read_csv("tables/Figure7B_GO_BP_representative10.csv",
               show_col_types = FALSE) %>%
  arrange(LogQ) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

pB <- ggplot(go, aes(x = MinusLog10FDR, y = Description)) +
  geom_col(fill = "#8ba1cb", width = 0.7) +
  labs(
    x = "-log10(FDR)",
    y = "GO Biological Processes"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

ggsave("figures/Figure7B_GO_BP.pdf",
       pB, width = 8, height = 6.5)

# ==========================================================
# C. KEGG enrichment bubble plot
# ==========================================================
kegg <- read_csv("tables/Figure7C_KEGG_top15.csv", show_col_types = FALSE) %>%
  arrange(LogQ, desc(GeneRatio), Description) %>%
  mutate(
    PlotLabel = stringr::str_wrap(Description, width = 38),
    PlotLabel = factor(PlotLabel, levels = rev(PlotLabel))
  )

pC <- ggplot(kegg, aes(x = GeneRatio, y = PlotLabel)) +
  geom_point(
    aes(size = Count, fill = MinusLog10FDR),
    shape = 21, colour = "black", stroke = 0.35, alpha = 0.95
  ) +
  scale_fill_gradientn(
    colors = c("#2166ac", "#fddbc7", "#b2182b"),
    name = "-log10(FDR)"
  ) +
  scale_size_continuous(range = c(3.5, 7.5), name = "Gene count") +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0.08, 0.12))
  ) +
  labs(x = "Gene ratio", y = "KEGG pathways") +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = "black", size = 8.5),
    axis.text.y = element_text(color = "black", size = 8, lineheight = 0.86),
    axis.title = element_text(face = "bold", size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7.5)
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      barwidth = grid::unit(32, "mm"),
      barheight = grid::unit(3.2, "mm")
    ),
    size = guide_legend(
      title.position = "top",
      override.aes = list(fill = "grey75")
    )
  )

ggsave("figures/Figure7C_KEGG.pdf",
       pC, width = 8.5, height = 7.2)

# ==========================================================
# E. LE8-associated chemicals → candidate CVD-related targets
#    Keep top 20 genes by number of linked chemicals to avoid a "hairball".
# ==========================================================
links <- read_csv("tables/TableS_chemical_candidate_gene_links.csv",
                  show_col_types = FALSE)

top20 <- links %>%
  distinct(Chemical, Gene) %>%
  count(Gene, name = "N_chemicals") %>%
  arrange(desc(N_chemicals), Gene) %>%
  slice_head(n = 20) %>%
  pull(Gene)

sankey <- links %>%
  filter(Gene %in% top20) %>%
  distinct(Chemical, Gene) %>%
  mutate(value = 1)

# Reproduce the manuscript's original two-column target-link style: fixed-size
# chemical and target nodes connected by thin, semi-transparent S-curves.
chemical_labels <- c(
  "1-Aminonaphthalene" = "1-AN",
  "2-Aminonaphthalene" = "2-AN",
  "2,5-Dimethylfuran" = "2,5-DMF",
  "2,6-Dimethylaniline" = "2,6-DMA",
  "4-Aminobiphenyl" = "4-ABP",
  "A-alpha-C" = "A-alpha-C",
  "Acrolein" = "Acrolein",
  "Acrylamide" = "Acrylamide",
  "Acrylonitrile" = "Acrylonitrile",
  "Benzene" = "Benzene",
  "Cadmium" = "Cd",
  "Ethylbenzene" = "Ethylbenzene",
  "Fluorene" = "Fluorene",
  "HEMA mercapturic acid" = "HEMA",
  "Isoprene" = "Isoprene",
  "Manganese" = "Mn",
  "MeA-alpha-C" = "MeA-alpha-C",
  "Mono-benzyl phthalate" = "MBzP",
  "N,N-Dimethylformamide" = "N,N-DMF",
  "Naphthalene" = "Naphthalene",
  "Phenanthrene" = "Phenanthrene",
  "Styrene" = "Styrene",
  "Toluene" = "Toluene",
  "Xylene" = "Xylene"
)

# High-contrast qualitative colours. Opposing hues are interleaved so that
# vertically adjacent chemicals do not receive visually similar colours.
chemical_colour_wheel <- grDevices::hcl(
  h = seq(10, 370, length.out = 25)[-25], c = 78, l = 60
)
chemical_contrast_order <- as.vector(rbind(1:12, 13:24))
chemical_colours <- chemical_colour_wheel[chemical_contrast_order]

target_colour_wheel <- grDevices::hcl(
  h = seq(20, 380, length.out = 21)[-21], c = 58, l = 70
)
target_contrast_order <- as.vector(rbind(1:10, 11:20))
target_colours <- target_colour_wheel[target_contrast_order]

chemical_nodes <- sankey %>%
  count(Chemical, name = "Degree") %>%
  arrange(Chemical) %>%
  mutate(
    y = seq(24, 1, length.out = n()),
    Label = unname(chemical_labels[Chemical]),
    Label = if_else(is.na(Label), Chemical, Label),
    Fill = chemical_colours[seq_len(n())]
  )

target_nodes <- sankey %>%
  count(Gene, name = "Degree") %>%
  arrange(desc(Degree), Gene) %>%
  mutate(
    y = seq(24, 1, length.out = n()),
    Fill = target_colours[seq_len(n())]
  )

link_positions <- sankey %>%
  left_join(chemical_nodes %>% select(Chemical, y_left = y, Colour = Fill),
            by = "Chemical") %>%
  left_join(target_nodes %>% select(Gene, y_right = y), by = "Gene") %>%
  mutate(LinkID = row_number())

curve_data <- bind_rows(lapply(seq_len(nrow(link_positions)), function(i) {
  t <- seq(0, 1, length.out = 80)
  u <- 1 - t
  row <- link_positions[i, ]
  tibble(
    LinkID = row$LinkID,
    x = 3 * u^2 * t * 0.34 + 3 * u * t^2 * 0.66 + t^3,
    y = u^3 * row$y_left + 3 * u^2 * t * row$y_left +
      3 * u * t^2 * row$y_right + t^3 * row$y_right,
    Colour = row$Colour
  )
}))

node_half_width <- 0.028
node_half_height <- 0.28

pE <- ggplot() +
  geom_path(
    data = curve_data,
    aes(x = x, y = y, group = LinkID, colour = Colour),
    linewidth = 0.42, alpha = 0.50, lineend = "round"
  ) +
  geom_rect(
    data = chemical_nodes,
    aes(xmin = -node_half_width, xmax = node_half_width,
        ymin = y - node_half_height, ymax = y + node_half_height,
        fill = Fill),
    colour = "#555555", linewidth = 0.35
  ) +
  geom_rect(
    data = target_nodes,
    aes(xmin = 1 - node_half_width, xmax = 1 + node_half_width,
        ymin = y - node_half_height, ymax = y + node_half_height,
        fill = Fill),
    colour = "#555555", linewidth = 0.35
  ) +
  geom_text(
    data = chemical_nodes,
    aes(x = -0.055, y = y, label = Label),
    hjust = 1, size = 2.7, colour = "black"
  ) +
  geom_text(
    data = target_nodes,
    aes(x = 1.055, y = y, label = Gene),
    hjust = 0, size = 2.7, colour = "black"
  ) +
  annotate("text", x = -0.28, y = 24.65, label = "E",
           fontface = "bold", size = 4) +
  scale_colour_identity() +
  scale_fill_identity() +
  coord_cartesian(xlim = c(-0.30, 1.23), ylim = c(0.45, 24.75), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(4, 5, 4, 5))

ggsave("figures/Figure7E_Chemical_Target_Sankey.pdf",
       pE, width = 7.0, height = 8.5)

# Save exact links used in Panel E.
write_csv(sankey, "tables/Figure7E_Sankey_links_top20_targets.csv")

# ==========================================================
# F. Functional MCODE modules in the PPI network
#    Six small-multiple networks avoid implying that distances between
#    disconnected modules have biological meaning. Node size represents
#    degree in the full 102-protein PPI network.
# ==========================================================
edges <- read_csv("tables/TableS_Metascape_PPI_edges.csv", show_col_types = FALSE)
mcode <- read_csv("tables/TableS_Metascape_MCODE_nodes.csv", show_col_types = FALSE)

g <- graph_from_data_frame(edges %>% select(Source, Target), directed = FALSE)

node_stat <- tibble(
  Gene = V(g)$name,
  Degree = degree(g),
  Betweenness = betweenness(g, normalized = TRUE)
) %>%
  left_join(mcode, by = "Gene") %>%
  arrange(desc(Degree), desc(Betweenness))

# Candidate hubs were defined conservatively as proteins in the top decile
# for both degree (local connectedness) and betweenness (bridging role) in the
# complete PPI network. This rule is independent of MCODE membership.
hub_degree_cutoff <- unname(quantile(node_stat$Degree, 0.90, type = 7))
hub_betweenness_cutoff <- unname(quantile(node_stat$Betweenness, 0.90, type = 7))

node_stat <- node_stat %>%
  mutate(
    Candidate_Hub = Degree >= hub_degree_cutoff &
      Betweenness >= hub_betweenness_cutoff,
    Hub_Criterion = if_else(
      Candidate_Hub,
      sprintf(
        "Degree >= %.0f and Betweenness >= %.5f (top decile for both)",
        hub_degree_cutoff, hub_betweenness_cutoff
      ),
      "Not retained"
    )
  )

write_csv(node_stat, "tables/TableS_PPI_topology_and_MCODE.csv")
write_csv(
  node_stat %>%
    filter(Candidate_Hub) %>%
    select(Gene, Degree, Betweenness, MCODE_Cluster, MCODE_Score, Hub_Criterion),
  "tables/TableS_candidate_hub_genes.csv"
)

module_info <- tribble(
  ~MCODE_Cluster, ~Module_title,
  1, "Lipoprotein response and\nmetabolic regulation",
  2, "Cardiovascular morphogenesis and\nTGF-beta/NOTCH signaling",
  3, "Receptor tyrosine kinase-PI3K\nsignaling",
  4, "MAPK/ERK signaling and\nfeedback regulation",
  5, "Cytokine-VEGF/PI3K-Akt\nsignaling",
  6, "Complement and coagulation\ncascades"
)

module_colours <- c(
  `1` = "#D55E00",
  `2` = "#009E73",
  `3` = "#3478B4",
  `4` = "#8E63B0",
  `5` = "#CC6677",
  `6` = "#B79F00"
)

degree_limits <- range(node_stat$Degree, na.rm = TRUE)

plot_mcode_module <- function(cluster_id) {
  module_nodes <- mcode %>%
    filter(MCODE_Cluster == cluster_id) %>%
    distinct(Gene, .keep_all = TRUE)

  module_edges <- edges %>%
    filter(Source %in% module_nodes$Gene,
           Target %in% module_nodes$Gene,
           Source != Target) %>%
    distinct(Source, Target)

  module_graph <- graph_from_data_frame(
    module_edges %>% select(Source, Target),
    directed = FALSE,
    vertices = module_nodes %>% transmute(name = Gene)
  )

  set.seed(7200 + cluster_id)
  module_layout <- if (vcount(module_graph) <= 3) {
    layout_in_circle(module_graph)
  } else {
    layout_with_kk(module_graph)
  }
  module_layout <- norm_coords(module_layout, xmin = -1, xmax = 1,
                               ymin = -1, ymax = 1)

  node_xy <- tibble(
    Gene = V(module_graph)$name,
    x = module_layout[, 1],
    y = module_layout[, 2]
  ) %>%
    left_join(
      node_stat %>% select(Gene, Degree, Betweenness, Candidate_Hub),
      by = "Gene"
    )

  edge_xy <- as_data_frame(module_graph, what = "edges") %>%
    rename(Source = from, Target = to) %>%
    left_join(node_xy %>% select(Gene, x_source = x, y_source = y),
              by = c("Source" = "Gene")) %>%
    left_join(node_xy %>% select(Gene, x_target = x, y_target = y),
              by = c("Target" = "Gene"))

  info <- module_info %>% filter(MCODE_Cluster == cluster_id)
  module_score <- unique(module_nodes$MCODE_Score)[1]
  panel_title <- sprintf(
    "M%d  %s\nn = %d; MCODE score = %.2f",
    cluster_id, info$Module_title, nrow(module_nodes), module_score
  )

  ggplot() +
    geom_segment(
      data = edge_xy,
      aes(x = x_source, y = y_source, xend = x_target, yend = y_target),
      colour = "#A99BC3", alpha = 0.55, linewidth = 0.55,
      lineend = "round"
    ) +
    geom_point(
      data = node_xy %>% filter(!Candidate_Hub),
      aes(x = x, y = y, size = Degree, shape = Candidate_Hub),
      fill = unname(module_colours[as.character(cluster_id)]),
      colour = "#263238", stroke = 0.5, alpha = 0.96,
      show.legend = cluster_id == 1
    ) +
    geom_point(
      data = node_xy %>% filter(Candidate_Hub),
      aes(x = x, y = y, size = Degree, shape = Candidate_Hub),
      fill = "#D73027", colour = "#5A1010", stroke = 0.6, alpha = 0.98,
      show.legend = cluster_id == 1
    ) +
    ggrepel::geom_text_repel(
      data = node_xy,
      aes(x = x, y = y, label = Gene),
      size = 2.35, family = "Helvetica", colour = "#1F1F1F",
      box.padding = 0.22, point.padding = 0.26,
      min.segment.length = 0, segment.colour = "#A0A0A0",
      segment.size = 0.25, max.overlaps = Inf,
      force = 1.8, force_pull = 0.25, max.time = 2,
      seed = 7200 + cluster_id
    ) +
    scale_size_continuous(
      limits = degree_limits, range = c(2.7, 7.2),
      breaks = c(5, 15, 25, 35), name = "Degree"
    ) +
    scale_shape_manual(
      values = c(`FALSE` = 21, `TRUE` = 21),
      breaks = c(FALSE, TRUE),
      labels = c("Other protein", "Candidate hub"),
      name = "Node type"
    ) +
    guides(
      size = guide_legend(
        title.position = "left", title.hjust = 0.5,
        override.aes = list(
          fill = "#AEB7BF", colour = "#263238", alpha = 1, shape = 21
        )
      ),
      shape = guide_legend(
        title.position = "left", title.hjust = 0.5,
        override.aes = list(
          fill = c("#AEB7BF", "#D73027"),
          colour = c("#263238", "#5A1010"),
          alpha = 1, size = 4, shape = c(21, 21)
        )
      )
    ) +
    scale_x_continuous(expand = expansion(mult = 0.28)) +
    scale_y_continuous(expand = expansion(mult = 0.28)) +
    coord_equal(clip = "off") +
    labs(title = panel_title) +
    theme_void(base_family = "Helvetica") +
    theme(
      plot.title = element_text(
        colour = unname(module_colours[as.character(cluster_id)]),
        face = "bold", size = 7.2, lineheight = 0.96,
        hjust = 0, margin = margin(b = 2.5)
      ),
      panel.border = element_rect(colour = "#D7D7D7", fill = NA,
                                  linewidth = 0.4),
      plot.margin = margin(5, 7, 4, 7),
      legend.position = "bottom",
      legend.title = element_text(size = 7, face = "bold"),
      legend.text = element_text(size = 6.5),
      legend.key.width = grid::unit(6, "mm")
    )
}

module_plots <- lapply(1:6, plot_mcode_module)

pF <- wrap_plots(module_plots, ncol = 3, guides = "collect") +
  plot_annotation(
    tag_levels = list("F"),
    theme = theme(
      plot.tag = element_text(family = "Helvetica", face = "bold", size = 10),
      plot.tag.position = c(0.006, 0.994),
      plot.margin = margin(5, 4, 4, 4)
    )
  ) &
  theme(legend.position = "bottom")

# ==========================================================
# Optional combined layout of A/B/C only.
# E and F are often clearer as full-width panels.
# ==========================================================
pB_combined <- pB + theme(
  axis.text = element_text(color = "black", size = 5.5),
  axis.title = element_text(face = "bold", size = 6.5)
)
pC_combined <- pC + theme(
  axis.text = element_text(color = "black", size = 5.5),
  axis.title = element_text(face = "bold", size = 6.5),
  legend.text = element_text(size = 5),
  legend.title = element_text(size = 5.5)
)
pA_combined <- pA
pA_combined$layers[[2]]$aes_params$label <- "Targets"
pA_combined$layers[[3]]$aes_params$label <- "CVD genes"
pA_combined$layers[[2]]$aes_params$size <- 1.7
pA_combined$layers[[3]]$aes_params$size <- 1.7
pA_combined$layers[[4]]$aes_params$size <- 2.6
pA_combined$layers[[5]]$aes_params$size <- 2.6
pA_combined$layers[[6]]$aes_params$size <- 2.6
pABC <- (pA_combined | pB_combined) / pC_combined
ggsave("figures/Figure7_ABC_preview.pdf",
       pABC, width = 14, height = 11)

# Publication exports and a lower-resolution visual-QA preview, all rendered
# with R graphics devices.
save_pub_r(pA, "figures/Figure7A_Targets_CVD_Venn", 90, 75)
save_pub_r(pB, "figures/Figure7B_GO_BP", 120, 95)
save_pub_r(pC, "figures/Figure7C_KEGG", 120, 120)
save_pub_r(pE, "figures/Figure7E_Chemical_Target_Sankey", 135, 165)
save_pub_r(pF, "figures/Figure7F_PPI_candidate_hubs", 183, 140)
save_pub_r(pABC, "figures/Figure7_ABC_preview", 183, 155)
