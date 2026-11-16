#!/usr/bin/env Rscript
# Sensitivity of macro-fraction to minspan threshold.
# 6 squamate outgroups only (no Gallus).
# Macro-derived autosomes → grey; Chr2, Chr8, Chr13 retain paper colours.
#
# Input:  synteny_blocks_detailed_ms<N>.tsv files in BLOCKS_DIR
#         (produced by parse_synteny_blocks.py for each minspan value)
# Output: synteny_sensitivity_minspan.pdf / .png

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

BLOCKS_DIR <- "synteny_blocks"   # directory with synteny_blocks_detailed_ms*.tsv files
OUTDIR     <- "."

# 6 squamate outgroups — Gallus excluded
OUTGROUP_SP <- c("Vipera latastei", "Aspidocelis tigris", "Varanus acanthurus",
                 "Iguana iguana", "Rhineura floridana", "Dibamus smithi")

CHR_ORDER  <- c(paste0("SUPER_", 1:18), "SUPER_Z")
CHR_LABELS <- c(as.character(1:17), "18 (true micro)", "Z")

# Panel assignment: fixed by P. algirus biology + squamate reclassification
# SUPER_1-12 = P. algirus macros → Macro-derived
# SUPER_13   = Macro-derived without Gallus (p_macro >> 0.5)
# SUPER_14-17, Z = micro-derived; SUPER_18 = true micro
CHR_CLASS <- c(
  "1"="Macro-derived",  "2"="Macro-derived",  "3"="Macro-derived",
  "4"="Macro-derived",  "5"="Macro-derived",  "6"="Macro-derived",
  "7"="Macro-derived",  "8"="Micro-derived",  "9"="Macro-derived",
  "10"="Macro-derived", "11"="Macro-derived", "12"="Macro-derived",
  "13"="Macro-derived", "17"="Macro-derived",
  "14"="Micro-derived", "15"="Micro-derived", "16"="Micro-derived",
  "Z"="Micro-derived",
  "18 (true micro)"="True micro"
)

GREY_MACRO <- "#949495"
CHR_COLORS <- c(
  "1"  = GREY_MACRO, "2"  = "#E96456", "3"  = GREY_MACRO,
  "4"  = GREY_MACRO, "5"  = GREY_MACRO, "6"  = GREY_MACRO,
  "7"  = GREY_MACRO, "8"  = "#FDBC6B", "9"  = GREY_MACRO,
  "10" = GREY_MACRO, "11" = GREY_MACRO, "12" = GREY_MACRO,
  "13" = "#42956B", "17" = "#B2D78A",
  "14" = "#BB8F6E", "15" = "#F6979A", "16" = "#8ED0C4",
  "Z"  = "#E36EAA",
  "18 (true micro)" = "#157AB7"
)

DRAW_ORDER <- c(
  as.character(c(1,3:7,9:12)), "17",
  "2", "13",
  "8", "14","15","16","Z",
  "18 (true micro)"
)

# ── helpers ───────────────────────────────────────────────────────────────────
load_tsv <- function(suffix) {
  f  <- file.path(BLOCKS_DIR, paste0("synteny_blocks_detailed_", suffix, ".tsv"))
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df %>%
    filter(species %in% OUTGROUP_SP) %>%
    group_by(species, psam_chr, size_class) %>%
    summarize(n_blocks = sum(n_blocks), n_genes = sum(n_genes), .groups = "drop") %>%
    as.data.frame()
}

pmacro_mean <- function(df) {
  macro <- df[df$size_class == "macro", c("species","psam_chr","n_blocks")]
  micro <- df[df$size_class == "micro", c("species","psam_chr","n_blocks")]
  names(macro)[3] <- "n_macro"; names(micro)[3] <- "n_micro"
  m <- merge(macro, micro, by = c("species","psam_chr"), all = TRUE)
  m[is.na(m)] <- 0
  m$p_macro <- with(m, ifelse(n_macro + n_micro > 0,
                               n_macro / (n_macro + n_micro), NA_real_))
  m$psam_chr <- factor(m$psam_chr, levels = CHR_ORDER, labels = CHR_LABELS)
  m %>%
    group_by(psam_chr) %>%
    summarize(
      mean_pm   = mean(p_macro, na.rm = TRUE),
      n_species = sum(!is.na(p_macro)),
      .groups   = "drop"
    ) %>%
    as.data.frame()
}

# ── discover available files ──────────────────────────────────────────────────
ALL_MS   <- seq(5, 70, by = 5)
SUFFIXES <- Filter(function(x) {
  file.exists(file.path(BLOCKS_DIR, paste0("synteny_blocks_detailed_ms", x, ".tsv")))
}, ALL_MS)

if (length(SUFFIXES) == 0)
  stop("No synteny_blocks_detailed_ms*.tsv files found in: ", BLOCKS_DIR)

# ── build sensitivity table ───────────────────────────────────────────────────
sens_list <- lapply(SUFFIXES, function(ms) {
  pmacro_mean(load_tsv(paste0("ms", ms))) %>% mutate(minspan = ms)
})
sens <- do.call(rbind, sens_list)

sens$chr_label <- as.character(sens$psam_chr)
sens$chr_class <- CHR_CLASS[sens$chr_label]

sens$mean_pm <- ifelse(
  sens$chr_class == "True micro" & sens$n_species < length(OUTGROUP_SP),
  NA_real_,
  sens$mean_pm
)

sens$chr_label <- factor(sens$chr_label, levels = DRAW_ORDER)
sens$chr_class <- factor(sens$chr_class,
                         levels = c("Macro-derived","Micro-derived","True micro"))
sens$lwd <- ifelse(sens$chr_label %in% c("2","13"), 1.1, 0.65)

sens_pts <- sens[!is.na(sens$mean_pm), ]

# ── plot ──────────────────────────────────────────────────────────────────────
x_breaks <- seq(10, 70, by = 10)
x_minor  <- seq(5,  70, by = 5)

p <- ggplot(sens, aes(x = minspan, y = mean_pm,
                      group = chr_label, color = chr_label)) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_line(aes(linewidth = lwd), alpha = 0.9, na.rm = FALSE) +
  geom_point(data = sens_pts, size = 2.2, shape = 16) +
  scale_linewidth_identity() +
  scale_color_manual(
    values = CHR_COLORS,
    breaks = CHR_LABELS,
    labels = paste("Chr", CHR_LABELS),
    name   = "P. algirus chr"
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1),
    name   = "Mean fraction macro (across outgroup species)"
  ) +
  scale_x_continuous(
    name         = "Minspan threshold",
    breaks       = x_breaks,
    minor_breaks = x_minor
  ) +
  facet_wrap(~ chr_class, ncol = 3) +
  guides(color = guide_legend(ncol = 1, keyheight = unit(0.38, "cm"),
                              override.aes = list(linewidth = 1.5, size = 0))) +
  labs(
    title    = "Sensitivity of macro-fraction to minspan threshold",
    subtitle = paste0(
      "Lines show mean p_macro across ", length(OUTGROUP_SP),
      " outgroup species;  flat = robust classification.\n",
      "True micro: points omitted when any species lacks data."
    )
  ) +
  theme_bw(base_size = 9) +
  theme(
    strip.text         = element_text(size = 8.5, face = "bold"),
    legend.position    = "right",
    legend.text        = element_text(size = 7),
    legend.key.size    = unit(0.35, "cm"),
    legend.title       = element_text(size = 7.5, face = "bold"),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.3),
    panel.grid.minor.x = element_line(color = "grey93", linewidth = 0.2),
    panel.grid.minor.y = element_blank(),
    plot.title         = element_text(size = 9, face = "bold"),
    plot.subtitle      = element_text(size = 7, color = "grey40")
  )

out_pdf <- file.path(OUTDIR, "synteny_sensitivity_minspan.pdf")
out_png <- file.path(OUTDIR, "synteny_sensitivity_minspan.png")

ggsave(out_pdf, p, width = 9.5, height = 4, units = "in", device = cairo_pdf)
ggsave(out_png, p, width = 9.5, height = 4, units = "in", dpi = 220)
message("Saved: ", out_pdf)
