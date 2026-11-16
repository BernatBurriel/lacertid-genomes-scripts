#!/usr/bin/env Rscript
# Micro-derived chromosome composition analysis.
#
# Input TSV columns: species, chr, label, psa_color, size_mb, gc, gene_density
# TREE_FILE: time-calibrated squamate species tree in Newick format (ultrametric)
#
# Each panel is saved separately.
# Panels:
#   panel_A_gc_scatter        — per-species GC vs log10(size) grid
#   panel_B_gd_scatter        — per-species gene density vs log10(size) grid
#   panel_C_lacertid_strip    — within-lacertid strip: macro vs micro-derived vs true-micro
#   panel_D_time_decay        — GC and gene density vs divergence time from Lacertidae crown

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ape)        # for read.tree / node.depth
})

# ── paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE))
TSV        <- file.path(SCRIPT_DIR, "outgroups_composition.tsv")
TREE_FILE  <- file.path(SCRIPT_DIR, "squamate_species_tree.tree")
OUTDIR     <- SCRIPT_DIR

# ── species metadata ──────────────────────────────────────────────────────────
SP_INFO <- data.frame(
  species = c("Gallus_gallus","Dibamus_smithi","Iguana_iguana","Varanus_acanthurus",
              "Vipera_latastei","Aspidocelis_tigris","Rhineura_floridana",
              "Psammodromus_algirus","Acanthodactylus_erythrurus","Eremias_argus",
              "Timon_lepidus","Lacerta_agilis","Iberolacerta_aurelioi",
              "Zootoca_vivipara","Scelarcis_perspicillata","Podarcis_muralis"),
  label   = c("G. gallus","D. smithi","I. iguana","V. acanthurus",
              "V. latastei","A. tigris","R. floridana",
              "P. algirus","A. erythrurus","E. argus",
              "T. lepidus","L. agilis","I. aurelioi",
              "Z. vivipara","S. perspicillata","P. muralis"),
  kind    = c(rep("outgroup",7), rep("lacertid",9)),
  stringsAsFactors = FALSE
)
OUTGROUP_SPS <- SP_INFO$species[SP_INFO$kind == "outgroup"]
LACERTID_SPS <- SP_INFO$species[SP_INFO$kind == "lacertid"]

# ── colour scheme ─────────────────────────────────────────────────────────────
# psa_color → category
COLOR_CAT <- c(
  "#C0C0C0" = "macro",    "#B2D78A" = "macro",
  "#E96456" = "micro_mac","#FDBC6B" = "micro_mac","#42956B" = "micro_mac",
  "#BB8F6E" = "micro_mac","#F6979A" = "micro_mac","#8ED0C4" = "micro_mac",
  "#E36EAA" = "micro_mac","#157AB7" = "true_micro"
)

CAT_ORDER  <- c("macro","micro_mac","true_micro")
CAT_LABELS <- c("Macro-derived","Micro-derived\n(in macro)","True micro\n(Chr18)")
CAT_COLORS <- c(macro = "#C0C0C0", micro_mac = "#888888", true_micro = "#157AB7")

# ── load data ─────────────────────────────────────────────────────────────────
df <- read.table(TSV, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df <- df[df$species %in% SP_INFO$species, ]
df$category <- COLOR_CAT[df$psa_color]
df$kind     <- SP_INFO$kind[match(df$species, SP_INFO$species)]
df$sp_label <- SP_INFO$label[match(df$species, SP_INFO$species)]

# ── divergence times from Lacertidae crown (requires ultrametric tree) ────────
tree <- read.tree(TREE_FILE)
tree <- keep.tip(tree, intersect(tree$tip.label, SP_INFO$species))

# Node age = distance from node to any leaf (ultrametric)
node_age <- function(tr, node) {
  tips <- tr$tip.label
  leaf <- tips[1]
  dist_mat <- dist.nodes(tr)
  dist_mat[node, which(tr$tip.label == leaf)]
}

lac_tips  <- intersect(LACERTID_SPS, tree$tip.label)
lac_crown <- getMRCA(tree, lac_tips)
lac_crown_age <- node_age(tree, lac_crown)

div_times <- sapply(OUTGROUP_SPS, function(sp) {
  if (!sp %in% tree$tip.label) return(NA_real_)
  sp_node <- which(tree$tip.label == sp)
  mrca    <- getMRCA(tree, c(sp_node, lac_tips[1]))
  node_age(tree, mrca)
})

df$div_time <- div_times[df$species]
message(sprintf("Lacertidae crown age: %.1f Mya", lac_crown_age))

# ── helpers ───────────────────────────────────────────────────────────────────
sig_label <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
}

save_plot <- function(p, name, w = 10, h = 7) {
  ggsave(file.path(OUTDIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUTDIR, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
}

# ── PANEL A: per-species GC scatter grid ─────────────────────────────────────
make_sp_grid <- function(y_key, y_label, panel_name) {
  plots <- lapply(seq_len(nrow(SP_INFO)), function(i) {
    sp  <- SP_INFO$species[i]
    lbl <- SP_INFO$label[i]
    sub <- df[df$species == sp, ]
    bg  <- if (SP_INFO$kind[i] == "lacertid") "#eef4ee" else "#f8f8f8"
    ggplot(sub, aes(x = log10(size_mb), y = .data[[y_key]])) +
      geom_point(aes(colour = psa_color), size = 1.5, alpha = 0.9, shape = 16) +
      scale_colour_identity() +
      labs(title = bquote(italic(.(lbl))), x = "log10 size (Mb)", y = y_label) +
      theme_bw(base_size = 7) +
      theme(plot.background = element_rect(fill = bg, colour = NA),
            plot.title = element_text(size = 6.5))
  })
  wrap_plots(plots, ncol = 4)
}

pA <- make_sp_grid("gc",           "GC content",             "panel_A_gc_scatter")
pB <- make_sp_grid("gene_density", "Gene density (genes/Mb)","panel_B_gd_scatter")
save_plot(pA, "panel_A_gc_scatter",  w = 12, h = 10)
save_plot(pB, "panel_B_gd_scatter",  w = 12, h = 10)

# ── PANEL C: within-lacertid strip plot ───────────────────────────────────────
lac_df <- df[df$kind == "lacertid" & !is.na(df$category), ]
lac_df$cat_f <- factor(lac_df$category, levels = CAT_ORDER, labels = CAT_LABELS)

make_strip <- function(y_key, y_label) {
  v_mac <- lac_df[[y_key]][lac_df$category == "macro"]
  v_mic <- lac_df[[y_key]][lac_df$category == "micro_mac"]
  v_tru <- lac_df[[y_key]][lac_df$category == "true_micro"]
  p1    <- wilcox.test(v_mac, v_mic, exact = FALSE)$p.value
  p2    <- wilcox.test(v_mic, v_tru, exact = FALSE)$p.value

  ymax <- max(lac_df[[y_key]], na.rm = TRUE)
  ann  <- data.frame(
    x    = c(1.5, 2.5),
    xend = c(1.5, 2.5),
    y    = c(ymax * 1.04, ymax * 1.08),
    label = c(sig_label(p1), sig_label(p2))
  )

  set.seed(0)
  ggplot(lac_df, aes(x = cat_f, y = .data[[y_key]])) +
    geom_jitter(aes(colour = psa_color), width = 0.18, size = 2, alpha = 0.85,
                shape = 16, show.legend = FALSE) +
    stat_summary(fun = median, geom = "crossbar", width = 0.44,
                 colour = "#222", linewidth = 0.8) +
    geom_text(data = ann, aes(x = x, y = y, label = label),
              hjust = 0.5, size = 3.5, colour = "#555") +
    scale_colour_identity() +
    labs(x = NULL, y = y_label) +
    theme_bw(base_size = 10)
}

pC1 <- make_strip("gc",           "GC content")
pC2 <- make_strip("gene_density", "Gene density (genes/Mb)")
save_plot(pC1 + pC2, "panel_C_lacertid_strip", w = 10, h = 4.5)

# ── PANEL D: time-decay scatter (GC and gene density vs divergence time) ──────
out_micro <- df[df$kind == "outgroup" & df$category %in% c("micro_mac","true_micro") &
                !is.na(df$div_time), ]

lac_medians <- lac_df %>%
  group_by(category) %>%
  summarise(gc = median(gc, na.rm = TRUE),
            gene_density = median(gene_density, na.rm = TRUE),
            .groups = "drop")

make_decay <- function(y_key, y_label) {
  sub_reg <- out_micro[out_micro$category == "micro_mac", ]
  fit     <- lm(sub_reg[[y_key]] ~ sub_reg$div_time)
  r       <- cor(sub_reg$div_time, sub_reg[[y_key]])
  p_val   <- summary(fit)$coefficients[2, 4]

  x_range <- range(out_micro$div_time, na.rm = TRUE)
  fit_line <- data.frame(
    div_time = seq(x_range[1], x_range[2], length.out = 100)
  )
  fit_line[[y_key]] <- predict(fit, newdata = data.frame(
    `sub_reg$div_time` = fit_line$div_time))

  lac_med_mic <- lac_medians[[y_key]][lac_medians$category == "micro_mac"]
  lac_med_tru <- lac_medians[[y_key]][lac_medians$category == "true_micro"]

  ggplot(out_micro, aes(x = div_time, y = .data[[y_key]])) +
    annotate("rect", xmin = 0, xmax = lac_crown_age,
             ymin = -Inf, ymax = Inf, alpha = 0.06, fill = "darkgreen") +
    geom_hline(yintercept = lac_med_mic, linetype = "dashed",
               colour = "#888888", linewidth = 0.8) +
    geom_hline(yintercept = lac_med_tru, linetype = "dashed",
               colour = "#157AB7", linewidth = 0.8) +
    geom_line(data = fit_line, aes(x = div_time, y = .data[[y_key]]),
              colour = "#333", linewidth = 1.2, alpha = 0.7) +
    geom_point(aes(colour = psa_color), size = 3, alpha = 0.85, shape = 16) +
    scale_colour_identity() +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("r = %.2f, %s", r, sig_label(p_val)),
             hjust = 1.1, vjust = 1.5, size = 3.5, colour = "#333") +
    labs(x = "Divergence from Lacertidae crown (Mya)", y = y_label) +
    theme_bw(base_size = 10)
}

pD1 <- make_decay("gc",           "GC content")
pD2 <- make_decay("gene_density", "Gene density (genes/Mb)")
save_plot(pD1 + pD2, "panel_D_time_decay", w = 10, h = 4.5)

message("Done. Panels saved to ", OUTDIR)
