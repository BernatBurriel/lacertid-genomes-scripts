#!/usr/bin/env Rscript
# GC content and gene density enrichment analysis.
#
# Input TSV columns: species, chr, label, psa_color, size_mb, gc, gene_density
# CLF_TXT: plain-text chromosome classification file (see _parse_classification below)
#
# Each panel is saved as a separate PDF/PNG.
# Panels:
#   panel_AB_allspecies_scatter  — GC and gene density vs log10(size), all species
#   panel_CD_outgroup_enrichment — enrichment relative to remaining macros, outgroups
#   panel_EF_lacertid_scatter    — GC and gene density vs log10(size), lacertids only
#   panel_GH_lacertid_enrichment — enrichment relative to remaining macros, lacertids

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# ── paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE))
TSV        <- file.path(SCRIPT_DIR, "outgroups_composition.tsv")
OUTDIR     <- SCRIPT_DIR

# ── species sets ──────────────────────────────────────────────────────────────
LACERTID_SP <- c(
  "Psammodromus_algirus", "Acanthodactylus_erythrurus", "Eremias_argus",
  "Timon_lepidus", "Lacerta_agilis", "Zootoca_vivipara",
  "Scelarcis_perspicillata", "Podarcis_muralis"
)
ALL_OUTGROUP_SP <- c(
  "Gallus_gallus", "Dibamus_smithi",
  "Vipera_latastei", "Iguana_iguana", "Varanus_acanthurus",
  "Rhineura_floridana", "Aspidocelis_tigris"
)
ENRICH_OUT_SP <- setdiff(ALL_OUTGROUP_SP, "Gallus_gallus")

# ── colours ───────────────────────────────────────────────────────────────────
MACRO_PSA      <- "#C0C0C0"
CHR2_PSA       <- "#E96456"
CHR13_PSA      <- "#42956B"
CHR17_PSA      <- "#B2D78A"
CHR18_PSA      <- "#157AB7"
MICRO_DER_COLS <- c("#FDBC6B","#BB8F6E","#F6979A","#8ED0C4","#E36EAA")
MACRO_REM_COLS <- c(MACRO_PSA, CHR17_PSA)

MACRO_FILL      <- "#D0D0D0"
MICRO_DER_FILL  <- "#F68C1F"
TRUE_MICRO_FILL <- "#1477A8"

LAC_COL        <- "#44B2A8"
TOXICOFERA_COL <- "#AA66A5"
RHINEURA_COL   <- "#7492CB"
ASPIDO_COL     <- "#EDE756"
GREY_DARK      <- "#555555"
GREY_MID       <- "#ABABAB"

SP_COLOR <- c(
  Gallus_gallus      = GREY_DARK,   Dibamus_smithi     = GREY_MID,
  Vipera_latastei    = TOXICOFERA_COL, Iguana_iguana  = TOXICOFERA_COL,
  Varanus_acanthurus = TOXICOFERA_COL, Rhineura_floridana = RHINEURA_COL,
  Aspidocelis_tigris = ASPIDO_COL,
  Psammodromus_algirus = LAC_COL, Acanthodactylus_erythrurus = LAC_COL,
  Eremias_argus = LAC_COL, Timon_lepidus = LAC_COL,
  Lacerta_agilis = LAC_COL, Zootoca_vivipara = LAC_COL,
  Scelarcis_perspicillata = LAC_COL, Podarcis_muralis = LAC_COL
)

LAC_SP_COL <- c(
  Psammodromus_algirus       = "#1B7837",
  Acanthodactylus_erythrurus = "#762A83",
  Eremias_argus              = "#4393C3",
  Timon_lepidus              = "#D6604D",
  Lacerta_agilis             = "#F4A582",
  Zootoca_vivipara           = "#92C5DE",
  Scelarcis_perspicillata    = "#878787",
  Podarcis_muralis           = "#E08214"
)

# ── 1. Load data ──────────────────────────────────────────────────────────────
df_raw <- read.table(TSV, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df     <- df_raw[df_raw$species %in% names(SP_COLOR), ]
df$log_size <- log10(df$size_mb)
df$group    <- ifelse(df$species %in% LACERTID_SP, "lacertid", "outgroup")

df_lac <- df[df$species %in% LACERTID_SP, ]
df_out <- df[df$species %in% ENRICH_OUT_SP, ]

# ── 2. Assign chromosome categories ──────────────────────────────────────────
classify_outgroup <- function(psa, is_macro) {
  ifelse(!is_macro, "true_micro",
    ifelse(psa == CHR13_PSA, "chr13",
      ifelse(psa == CHR2_PSA, "chr2", "macro_rem")))
}

classify_lacertid <- function(psa) {
  ifelse(psa == CHR18_PSA, "true_micro",
    ifelse(psa %in% MICRO_DER_COLS, "micro_est",
      ifelse(psa == CHR13_PSA, "chr13",
        ifelse(psa == CHR2_PSA, "chr2",
          ifelse(psa %in% MACRO_REM_COLS, "macro_rem", NA_character_)))))
}

# For outgroups, is_macro must be supplied from the classification file or TSV
# Here we derive it from psa_color: micro = Chr18 or any outgroup micro color
# (adjust if your TSV already has an is_macro column)
df_out$is_macro  <- !(df_out$psa_color == CHR18_PSA)
df_out$enr_cat   <- classify_outgroup(df_out$psa_color, df_out$is_macro)
df_lac$enr_cat   <- classify_lacertid(df_lac$psa_color)

# ── 3. Per-species enrichment relative to remaining macros ────────────────────
add_enrichment <- function(df_sp) {
  ref <- df_sp %>%
    filter(enr_cat == "macro_rem") %>%
    group_by(species) %>%
    summarise(ref_gc = mean(gc), ref_gd = mean(gene_density), .groups = "drop")
  df_sp %>%
    left_join(ref, by = "species") %>%
    mutate(
      gc_enr = 100 * (gc - ref_gc),
      gd_enr = gene_density - ref_gd
    )
}

df_out <- add_enrichment(df_out)
df_lac <- add_enrichment(df_lac)

# ── helpers ───────────────────────────────────────────────────────────────────
save_plot <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(OUTDIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUTDIR, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
}

XLABEL <- expression(log[10] ~ "chromosome size (Mb)")

enrich_panel <- function(df_sp, enr_col, y_label, cats, cat_labels, sp_col_map) {
  df_sp <- df_sp %>% filter(!is.na(.data[[enr_col]]) & !is.na(enr_cat) &
                              enr_cat %in% cats)
  df_sp$enr_cat <- factor(df_sp$enr_cat, levels = cats)

  sp_means <- df_sp %>%
    group_by(species, enr_cat) %>%
    summarise(val = mean(.data[[enr_col]], na.rm = TRUE), .groups = "drop") %>%
    mutate(sp_col = sp_col_map[species])

  set.seed(42)
  ggplot(df_sp, aes(x = enr_cat, y = .data[[enr_col]])) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
    geom_jitter(aes(colour = psa_color), width = 0.22, size = 1.5, alpha = 0.35,
                show.legend = FALSE) +
    geom_jitter(data = sp_means, aes(x = enr_cat, y = val, colour = sp_col),
                width = 0.12, size = 3.5, shape = 16, alpha = 0.92,
                show.legend = FALSE) +
    scale_colour_identity() +
    scale_x_discrete(labels = setNames(cat_labels, cats)) +
    labs(x = NULL, y = y_label) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank())
}

# ── PANEL A–B: all-species scatter ───────────────────────────────────────────
make_AB <- function(y_key, y_label) {
  df$sp_col <- SP_COLOR[df$species]
  df$is_lac <- df$species %in% LACERTID_SP
  ggplot(df, aes(x = log_size, y = .data[[y_key]])) +
    geom_point(data = df[df$is_lac, ], colour = LAC_COL,
               size = 2.2, alpha = 1, shape = 16) +
    geom_point(data = df[!df$is_lac, ], aes(colour = sp_col),
               size = 2.2, alpha = 0.9, shape = 16) +
    scale_colour_identity() +
    labs(x = XLABEL, y = y_label) +
    theme_bw(base_size = 10)
}

pA <- make_AB("gc",           "GC content")
pB <- make_AB("gene_density", "Gene density (genes/Mb)")
save_plot(pA + pB, "panel_AB_allspecies_scatter", w = 11, h = 4.5)

# ── PANEL C–D: outgroup enrichment ───────────────────────────────────────────
OUT_CATS   <- c("true_micro","chr13","chr2","macro_rem")
OUT_LABELS <- c("Retained\nmicro","Chr 13\nhomologues","Chr 2\nhomologues","Remaining\nmacros")

pC <- enrich_panel(df_out, "gc_enr",  "GC enrichment (pp)",         OUT_CATS, OUT_LABELS, SP_COLOR)
pD <- enrich_panel(df_out, "gd_enr",  "Gene-density enrichment\n(genes/Mb)", OUT_CATS, OUT_LABELS, SP_COLOR)
save_plot(pC + pD, "panel_CD_outgroup_enrichment", w = 11, h = 4.5)

# ── PANEL E–F: lacertid scatter ───────────────────────────────────────────────
make_EF <- function(y_key, y_label) {
  df_lac_plot <- df_lac
  df_lac_plot$is_macro_rem <- df_lac_plot$psa_color %in% MACRO_REM_COLS
  ggplot(df_lac_plot, aes(x = log_size, y = .data[[y_key]])) +
    geom_point(data = df_lac_plot[df_lac_plot$is_macro_rem, ],
               aes(colour = psa_color), size = 3, alpha = 1, shape = 16) +
    geom_point(data = df_lac_plot[!df_lac_plot$is_macro_rem, ],
               aes(colour = psa_color), size = 3, alpha = 1, shape = 16) +
    scale_colour_identity() +
    labs(x = XLABEL, y = y_label) +
    theme_bw(base_size = 10)
}

pE <- make_EF("gc",           "GC content")
pF <- make_EF("gene_density", "Gene density (genes/Mb)")
save_plot(pE + pF, "panel_EF_lacertid_scatter", w = 11, h = 4.5)

# ── PANEL G–H: lacertid enrichment ───────────────────────────────────────────
LAC_CATS   <- c("true_micro","micro_est","chr13","chr2","macro_rem")
LAC_LABELS <- c("Retained\nmicro","Micro-\nderived","Chr 13\nhomologues","Chr 2\nhomologues","Remaining\nmacros")

pG <- enrich_panel(df_lac, "gc_enr",  "GC enrichment (pp)",         LAC_CATS, LAC_LABELS, LAC_SP_COL)
pH <- enrich_panel(df_lac, "gd_enr",  "Gene-density enrichment\n(genes/Mb)", LAC_CATS, LAC_LABELS, LAC_SP_COL)
save_plot(pG + pH, "panel_GH_lacertid_enrichment", w = 11, h = 4.5)

message("Done. Panels saved to ", OUTDIR)
