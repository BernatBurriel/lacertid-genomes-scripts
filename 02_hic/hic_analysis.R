## ============================================================================
## hic_analysis.R
## Hi-C nuclear-positioning pipeline in R, reading ICE-balanced .cool files
## directly via rhdf5.
##
## Dependencies:
##   BiocManager::install("rhdf5")
##   install.packages(c("data.table","Matrix","ggplot2","patchwork",
##                       "minpack.lm","lme4","lmerTest","viridis","scales"))
##
## Input:  ICE-balanced .cool files at 50 kb resolution (one per species)
## Output: TSV tables + PNG/PDF figures
## ============================================================================

library(rhdf5)
library(data.table)
library(Matrix)
library(ggplot2)
library(patchwork)
library(minpack.lm)
library(lme4)
library(lmerTest)
library(viridis)
library(scales)

## ── 0. Configuration ─────────────────────────────────────────────────────────

COOL_DIR <- "cool_files"   # directory with .cool files (one per species)
OUT_DIR  <- "HiC_figures"
DATA_DIR <- file.path(OUT_DIR, "data")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIES <- list(
  Psam_alg        = file.path(COOL_DIR, "contacts.Psam_alg.50000.cool"),
  acanthodactylus = file.path(COOL_DIR, "contacts.acanthodactylus.50000.cool"),
  timon           = file.path(COOL_DIR, "contacts.timon.50000.cool"),
  scelarcis       = file.path(COOL_DIR, "contacts.scelarcis.50000.cool")
)

OUTGROUPS <- list(
  aspidocelis = file.path(COOL_DIR, "contacts.aspidocelis.50000.cool"),
  rhineura    = file.path(COOL_DIR, "contacts.rhineura.50000.cool")
)

# Chromosome classification per species (homology-based, not size-based)
# true_micro = genuine microchromosome; micro_derived = fused; macro_derived = ancestral macro
CHROM_CLASS <- list(
  Psam_alg = data.table(
    chrom = paste0("SUPER_", 1:20),
    class = c("macro_derived","macro_derived","macro_derived","macro_derived",
              "macro_derived","macro_derived","macro_derived","micro_derived",  # 8
              "macro_derived","macro_derived","macro_derived","macro_derived",
              "macro_derived","micro_derived","micro_derived","micro_derived",  # 14,15,16
              "macro_derived","true_micro","Z","W")
  ),
  acanthodactylus = data.table(
    chrom = paste0("SUPER_", 1:19),
    class = c("macro_derived","macro_derived","macro_derived","macro_derived",
              "macro_derived","macro_derived","macro_derived","micro_derived",  # 8
              "macro_derived","macro_derived","macro_derived","macro_derived",
              "micro_derived","macro_derived","micro_derived","macro_derived",  # 13,15
              "micro_derived","true_micro","Z")                                  # 17
  ),
  timon = data.table(
    chrom = paste0("SUPER_", 1:18),
    class = c("macro_derived","macro_derived","macro_derived","macro_derived",
              "macro_derived","macro_derived","micro_derived","macro_derived",  # 7
              "macro_derived","macro_derived","macro_derived","micro_derived",  # 12
              "macro_derived","micro_derived","micro_derived","macro_derived",  # 14,15
              "true_micro","Z")
  ),
  scelarcis = data.table(
    chrom = paste0("SUPER_", 1:19),
    class = c("macro_derived","macro_derived","macro_derived","macro_derived",
              "macro_derived","macro_derived","macro_derived","micro_derived",  # 8
              "macro_derived","macro_derived","macro_derived","macro_derived",
              "micro_derived","macro_derived","micro_derived","micro_derived",  # 13,15,16
              "macro_derived","true_micro","Z")
  )
)

# Fusion class for micro-derived chromosomes (triple = chr8; double = others)
FUSION_CLASS <- list(
  Psam_alg        = c(SUPER_8  = "triple", SUPER_14 = "double",
                      SUPER_15 = "double", SUPER_16 = "double"),
  acanthodactylus = c(SUPER_8  = "triple", SUPER_13 = "double",
                      SUPER_15 = "double", SUPER_17 = "double"),
  timon           = c(SUPER_7  = "triple", SUPER_12 = "double",
                      SUPER_14 = "double", SUPER_15 = "double"),
  scelarcis       = c(SUPER_8  = "triple", SUPER_13 = "double",
                      SUPER_15 = "double", SUPER_16 = "double")
)

# Outgroup chromosome classes (homology-based)
OUTGROUP_CLASS <- list(
  aspidocelis = data.table(
    chrom = c(paste0("CM10788", 0:9), "CM107890",
              paste0("CM10789", 1:9), "CM107900","CM107901","CM107902"),
    class = c(rep("macro", 11), rep("micro", 12))
  ),
  rhineura = data.table(
    chrom = paste0("Chromosome_", 1:22),
    class = c(rep("macro", 11), rep("micro", 11))
  )
)

CHR_COLORS <- c(macro_derived = "#4A90D9", micro_derived = "#E8762B",
                true_micro = "#7B3294", Z = "#999999")

## ── 1. Helper: read ICE-balanced .cool file ───────────────────────────────────

read_cool <- function(cool_path, resolution = 50000) {
  chroms <- data.table(
    name   = as.character(h5read(cool_path, "chroms/name")),
    length = as.numeric(h5read(cool_path, "chroms/length"))
  )
  bins <- data.table(
    chrom  = as.character(h5read(cool_path, "bins/chrom")),
    start  = as.numeric(h5read(cool_path, "bins/start")),
    end    = as.numeric(h5read(cool_path, "bins/end")),
    weight = as.numeric(h5read(cool_path, "bins/weight"))
  )
  bins[, bin_id := .I - 1L]

  pixels <- data.table(
    bin1_id = as.integer(h5read(cool_path, "pixels/bin1_id")),
    bin2_id = as.integer(h5read(cool_path, "pixels/bin2_id")),
    count   = as.numeric(h5read(cool_path, "pixels/count"))
  )
  w <- bins$weight
  pixels[, balanced := count * w[bin1_id + 1L] * w[bin2_id + 1L]]
  pixels <- pixels[is.finite(balanced)]

  list(chroms = chroms, bins = bins, pixels = pixels)
}

## ── 2. Helper: chromosome × chromosome contact matrix ────────────────────────

chrom_contact_matrix <- function(cool_data, chroms_keep = NULL) {
  bins   <- cool_data$bins
  pixels <- cool_data$pixels

  pix <- pixels[, .(bin1_id, bin2_id, balanced)]
  pix[, chrom1 := bins$chrom[bin1_id + 1L]]
  pix[, chrom2 := bins$chrom[bin2_id + 1L]]

  if (!is.null(chroms_keep))
    pix <- pix[chrom1 %in% chroms_keep & chrom2 %in% chroms_keep]

  agg <- pix[, .(contacts = sum(balanced)), by = .(chrom1, chrom2)]

  chr_names <- unique(c(agg$chrom1, agg$chrom2))
  if (!is.null(chroms_keep)) chr_names <- chroms_keep[chroms_keep %in% chr_names]
  n   <- length(chr_names)
  mat <- matrix(0, n, n, dimnames = list(chr_names, chr_names))
  for (i in seq_len(nrow(agg))) {
    c1 <- agg$chrom1[i]; c2 <- agg$chrom2[i]; v <- agg$contacts[i]
    mat[c1, c2] <- mat[c1, c2] + v
    if (c1 != c2) mat[c2, c1] <- mat[c2, c1] + v
  }
  mat
}

## ── 3. Helper: log₂(O/E) ─────────────────────────────────────────────────────

log2oe <- function(mat) {
  grand  <- sum(mat)
  marg   <- rowSums(mat)
  expect <- outer(marg, marg) / grand
  lo     <- log2(mat / expect)
  lo[!is.finite(lo)] <- NA
  lo
}

## ── 4. Helper: trans-contact fraction ────────────────────────────────────────

trans_fraction <- function(mat) {
  total <- rowSums(mat)
  cis   <- diag(mat)
  (total - cis) / total
}

## ── 5. Helper: PCA on log₂(O/E) matrix ───────────────────────────────────────

run_pca <- function(lo_mat) {
  lo_notrans <- lo_mat
  diag(lo_notrans) <- 0
  lo_notrans[is.na(lo_notrans)] <- 0
  lo_scaled <- t(scale(t(lo_notrans)))
  pca    <- prcomp(lo_scaled, center = TRUE, scale. = FALSE)
  scores <- as.data.frame(pca$x[, 1:2])
  scores$chrom <- rownames(lo_mat)
  colnames(scores)[1:2] <- c("PC1", "PC2")
  list(scores = scores, pca = pca,
       var_explained = summary(pca)$importance[2, 1:2])
}

## ── 6. Main loop: per-species analysis ───────────────────────────────────────

results_list <- list()
log2oe_pairs <- list()
trans_list   <- list()

for (sp in names(SPECIES)) {
  message("Processing: ", sp)
  cool_data  <- read_cool(SPECIES[[sp]])
  class_tbl  <- CHROM_CLASS[[sp]]

  keep_chroms   <- class_tbl$chrom
  keep_analysis <- class_tbl[class != "Z" & class != "W", chrom]

  mat_all  <- chrom_contact_matrix(cool_data, chroms_keep = keep_chroms)
  mat_anal <- mat_all[keep_analysis, keep_analysis]
  lo       <- log2oe(mat_anal)

  idx      <- which(upper.tri(lo, diag = FALSE), arr.ind = TRUE)
  pairs_dt <- data.table(
    species = sp,
    chrom1  = rownames(lo)[idx[, 1]],
    chrom2  = colnames(lo)[idx[, 2]],
    log2_oe = lo[idx]
  )
  log2oe_pairs[[sp]] <- pairs_dt

  mat_full_anal <- mat_all[keep_chroms[keep_chroms %in% rownames(mat_all)],
                            keep_chroms[keep_chroms %in% colnames(mat_all)]]
  tf       <- trans_fraction(mat_full_anal)
  sizes_mb <- cool_data$bins[chrom %in% names(tf),
                               .(size_mb = sum(end - start) / 1e6), by = chrom]
  setkey(sizes_mb, chrom)
  tf_dt <- data.table(species = sp, chrom = names(tf), trans_frac = as.numeric(tf))
  tf_dt <- merge(tf_dt, class_tbl, by = "chrom")
  tf_dt <- merge(tf_dt, sizes_mb,  by = "chrom")
  trans_list[[sp]] <- tf_dt

  pca_res <- run_pca(lo)
  scores  <- as.data.table(pca_res$scores)

  micro_mean <- scores[chrom %in% class_tbl[class == "micro_derived", chrom], mean(PC1)]
  macro_mean <- scores[chrom %in% class_tbl[class == "macro_derived", chrom], mean(PC1)]
  if (micro_mean < macro_mean) scores[, c("PC1","PC2") := .(-PC1, -PC2)]

  scores[, PC1_z := scale(PC1)[, 1]]
  scores  <- merge(scores, class_tbl, by = "chrom")
  fusion_vec <- FUSION_CLASS[[sp]]
  scores[, fusion_class := fifelse(class == "micro_derived",
                                    fusion_vec[chrom], class)]
  scores[, species    := sp]
  scores[, pc1_var    := pca_res$var_explained[1]]

  results_list[[sp]] <- scores
}

## ── 7. Assemble combined tables ───────────────────────────────────────────────

all_results <- rbindlist(results_list)
all_trans   <- rbindlist(trans_list)
all_log2oe  <- rbindlist(log2oe_pairs)

fwrite(all_trans,  file.path(DATA_DIR, "hic_trans_fractions.tsv"), sep = "\t")
fwrite(all_log2oe, file.path(DATA_DIR, "hic_log2oe_pairs.tsv"),    sep = "\t")

## ── 8. Statistical tests ──────────────────────────────────────────────────────

stat_per_sp <- all_results[class %in% c("micro_derived","macro_derived"),
  .(p_value   = wilcox.test(PC1_z[class == "micro_derived"],
                            PC1_z[class == "macro_derived"],
                            alternative = "greater")$p.value,
    micro_mean = mean(PC1_z[class == "micro_derived"]),
    macro_mean = mean(PC1_z[class == "macro_derived"])),
  by = species]

pooled_wt <- wilcox.test(
  all_results[class == "micro_derived", PC1_z],
  all_results[class == "macro_derived", PC1_z],
  alternative = "greater"
)

trans_wt <- wilcox.test(
  all_trans[class == "micro_derived", trans_frac],
  all_trans[class == "macro_derived", trans_frac],
  alternative = "greater"
)

cat("Per-species Wilcoxon (PC1_z, micro > macro):\n")
print(stat_per_sp)
cat("\nPooled Wilcoxon (PC1_z): W =", pooled_wt$statistic,
    " p =", pooled_wt$p.value, "\n")
cat("Pooled Wilcoxon (trans_frac): W =", trans_wt$statistic,
    " p =", trans_wt$p.value, "\n")

## ── 9. Size regression — 5 models per species ─────────────────────────────────

reg_data  <- all_results[class %in% c("macro_derived","micro_derived")]
sizes_all <- all_trans[, .(species, chrom, size_mb, class)]
reg_data  <- merge(reg_data, sizes_all[, .(species, chrom, size_mb)],
                   by = c("species","chrom"), all.x = TRUE)
reg_data[, log10_size := log10(size_mb)]

fit_sigmoid <- function(x, y) {
  start <- list(a = max(y), d = min(y), c = median(x), b = -2)
  tryCatch(
    nls.lm(par = start, fn = function(p, x, y) {
      pred <- p$d + (p$a - p$d) / (1 + (x/p$c)^p$b)
      y - pred
    }, x = x, y = y)$par,
    error = function(e) NULL
  )
}

aic_table <- list()
for (sp in names(SPECIES)) {
  d  <- reg_data[species == sp]
  x  <- d$log10_size; y <- d$PC1_z

  m_lin  <- lm(y ~ x);              aic_lin  <- AIC(m_lin)
  m_quad <- lm(y ~ x + I(x^2));     aic_quad <- AIC(m_quad)
  sig_par <- fit_sigmoid(x, y)
  aic_sig <- if (!is.null(sig_par)) {
    pred <- sig_par$d + (sig_par$a - sig_par$d) / (1 + (x/sig_par$c)^sig_par$b)
    rss  <- sum((y - pred)^2)
    n    <- length(y); k <- 4
    n * log(rss/n) + 2*k
  } else NA_real_

  aic_table[[sp]] <- data.table(species = sp,
    AIC_linear = aic_lin, AIC_quadratic = aic_quad, AIC_sigmoid = aic_sig,
    dAIC_sigmoid_vs_linear = aic_sig - aic_lin)
}
aic_dt <- rbindlist(aic_table)
cat("\nAIC comparison (size models):\n"); print(aic_dt)

## ── 10. Compartment enrichment (class-pair log₂(O/E)) ────────────────────────

enrich_data  <- all_log2oe
class_lookup <- rbindlist(lapply(names(CHROM_CLASS), function(sp)
  cbind(species = sp, CHROM_CLASS[[sp]])))
enrich_data <- merge(enrich_data, class_lookup,
                     by.x = c("species","chrom1"), by.y = c("species","chrom"))
setnames(enrich_data, "class", "class1")
enrich_data <- merge(enrich_data, class_lookup,
                     by.x = c("species","chrom2"), by.y = c("species","chrom"))
setnames(enrich_data, "class", "class2")

enrich_data[, pair := paste(pmin(class1, class2), pmax(class1, class2), sep = "×")]

enrich_sum <- enrich_data[, .(mean_log2oe = mean(log2_oe, na.rm = TRUE)),
                            by = .(class1, class2, pair)]

## ── 11. Multiple regression + mediation ──────────────────────────────────────

reg_full <- fread(file.path(OUT_DIR, "hic_regression_dataset.tsv"))

m1 <- lm(PC1_z ~ is_micro, data = reg_full)
m2 <- lm(PC1_z ~ is_micro + log10_size, data = reg_full)
m3 <- lm(PC1_z ~ is_micro + log10_size + mean_GC + mean_gene_density, data = reg_full)
m4 <- lm(PC1_z ~ log10_size + mean_GC + mean_gene_density, data = reg_full)
m5 <- lmer(PC1_z ~ is_micro + log10_size + mean_GC + mean_gene_density + (1|species),
           data = reg_full, REML = FALSE)

m_gc   <- lm(mean_GC ~ is_micro + log10_size, data = reg_full)
a_path <- coef(m_gc)["is_micro"]
b_path <- coef(m3)["mean_GC"]
indirect  <- a_path * b_path
total_m2  <- coef(m2)["is_micro"]
prop_med  <- indirect / total_m2
cat("\nMediation through GC:\n",
    "  a =", round(a_path,4), " b =", round(b_path,2),
    " indirect =", round(indirect,4),
    " proportion =", round(prop_med,3), "\n")

## ── 12. Horseshoe arch deviation ─────────────────────────────────────────────

arch_dev_list <- list()
for (sp in names(SPECIES)) {
  d        <- all_results[species == sp]
  arch_fit <- lm(PC2 ~ PC1 + I(PC1^2), data = d)
  d[, arch_dev := PC2 - predict(arch_fit, newdata = d)]
  arch_dev_list[[sp]] <- d
}
all_arch <- rbindlist(arch_dev_list)

## ── 13. Outgroup analysis ─────────────────────────────────────────────────────

outgroup_results <- list()
for (sp in names(OUTGROUPS)) {
  message("Outgroup: ", sp)
  cool_data   <- read_cool(OUTGROUPS[[sp]])
  class_tbl   <- OUTGROUP_CLASS[[sp]]
  keep_chroms <- class_tbl$chrom

  mat     <- chrom_contact_matrix(cool_data, chroms_keep = keep_chroms)
  lo      <- log2oe(mat)
  pca_res <- run_pca(lo)
  scores  <- as.data.table(pca_res$scores)

  micro_mean <- scores[chrom %in% class_tbl[class == "micro", chrom], mean(PC1)]
  macro_mean <- scores[chrom %in% class_tbl[class == "macro", chrom], mean(PC1)]
  if (micro_mean < macro_mean) scores[, c("PC1","PC2") := .(-PC1, -PC2)]

  scores  <- merge(scores, class_tbl, by = "chrom")
  scores[, species := sp]

  tf       <- trans_fraction(mat)
  sizes_mb <- cool_data$bins[chrom %in% names(tf),
                               .(size_mb = sum(end - start) / 1e6), by = chrom]
  tf_dt    <- merge(data.table(chrom = names(tf), trans_frac = as.numeric(tf)),
                    class_tbl, by = "chrom")
  tf_dt    <- merge(tf_dt, sizes_mb, by = "chrom")
  tf_dt[, species := sp]

  wt <- wilcox.test(scores[class == "micro", PC1], scores[class == "macro", PC1],
                    alternative = "greater")
  cat(sp, ": PC1 var =", round(pca_res$var_explained[1]*100, 1), "%",
      " micro vs macro p =", signif(wt$p.value, 3), "\n")

  tf_dt[, log10_size := log10(size_mb)]
  m_out <- lm(trans_frac ~ class + log10_size, data = tf_dt)
  cat("  trans_frac ~ class + size: class(micro) p =",
      signif(summary(m_out)$coef["classmicro","Pr(>|t|)"], 3), "\n")

  outgroup_results[[sp]] <- list(scores = scores, trans = tf_dt,
                                  wt = wt, pca = pca_res)
}

## ── 14. Figures ───────────────────────────────────────────────────────────────

save_plot <- function(p, name, w = 10, h = 7) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h)
}

# 14a. Trans-fraction barplot
p_trans <- ggplot(
    all_trans[class != "Z" & class != "W"],
    aes(x = reorder(chrom, -size_mb), y = trans_frac, fill = class)) +
  geom_col(colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = CHR_COLORS) +
  facet_wrap(~species, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "Trans-contact fraction", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom")
save_plot(p_trans, "trans_fractions", w = 14, h = 5)

# 14b. Pooled PC1 z-score strip plot
p_pooled <- ggplot(
    all_results[class %in% c("macro_derived","micro_derived","true_micro")],
    aes(x = class, y = PC1_z, colour = class)) +
  geom_jitter(width = 0.15, size = 2.5, stroke = 0.4, shape = 21, fill = NA) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4,
               colour = "black", linewidth = 0.6) +
  scale_colour_manual(values = CHR_COLORS) +
  facet_wrap(~species, nrow = 1) +
  labs(x = NULL, y = "PC1 z-score", colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "none")
save_plot(p_pooled, "PCA_pooled", w = 12, h = 5)

# 14c. Per-species horseshoe (PC1 vs PC2)
p_horse_list <- lapply(names(SPECIES), function(sp) {
  d <- all_arch[species == sp]
  arc <- lm(PC2 ~ PC1 + I(PC1^2), data = d)
  x_seq    <- seq(min(d$PC1), max(d$PC1), length.out = 200)
  arc_line <- data.frame(PC1 = x_seq,
                          PC2 = predict(arc, newdata = data.frame(PC1 = x_seq)))
  ggplot(d, aes(x = PC1, y = PC2, colour = class)) +
    geom_line(data = arc_line, aes(x = PC1, y = PC2),
              colour = "grey60", linetype = "dashed", inherit.aes = FALSE) +
    geom_point(size = 3.5, stroke = 0.5, shape = 21, aes(fill = class, colour = "black")) +
    geom_text(aes(label = chrom), size = 2.2, vjust = -0.8, show.legend = FALSE) +
    scale_fill_manual(values = CHR_COLORS) +
    scale_colour_manual(values = CHR_COLORS) +
    labs(title = sp, x = "PC1", y = "PC2", fill = NULL) +
    theme_bw(base_size = 10) + theme(legend.position = "none")
})
p_horse <- wrap_plots(p_horse_list, nrow = 2)
save_plot(p_horse, "PCA12_horseshoe", w = 14, h = 10)

# 14d. PC1 vs size with sigmoidal fit
p_sig_list <- lapply(names(SPECIES), function(sp) {
  d  <- reg_data[species == sp]
  x  <- d$log10_size; y <- d$PC1_z
  sig_par <- fit_sigmoid(x, y)
  x_seq   <- seq(min(x), max(x), length.out = 200)
  fit_dt  <- if (!is.null(sig_par)) {
    data.table(log10_size = x_seq,
               PC1_z = sig_par$d + (sig_par$a - sig_par$d) /
                         (1 + (x_seq/sig_par$c)^sig_par$b))
  } else data.table(log10_size = x_seq, PC1_z = predict(lm(y~x), data.frame(x=x_seq)))

  ggplot(d, aes(x = log10_size, y = PC1_z, colour = class)) +
    geom_line(data = fit_dt, aes(x = log10_size, y = PC1_z),
              colour = "grey40", linetype = "solid", inherit.aes = FALSE) +
    geom_point(size = 3, stroke = 0.6, shape = 21, aes(fill = class)) +
    scale_fill_manual(values   = CHR_COLORS) +
    scale_colour_manual(values = CHR_COLORS) +
    labs(title = sp, x = "log10(chromosome size, Mb)", y = "PC1 z-score") +
    theme_bw(base_size = 10) + theme(legend.position = "none")
})
p_sig <- wrap_plots(p_sig_list, nrow = 2)
save_plot(p_sig, "PCA_size_sigmoidal", w = 12, h = 9)

# 14e. GC composition retention
reg_full[, class_label := fifelse(is_micro == 1, "micro_derived", "macro_derived")]
p_gc <- ggplot(reg_full, aes(x = log10_size, y = mean_GC, colour = class_label)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  scale_colour_manual(values = c(macro_derived = "#4A90D9", micro_derived = "#E8762B")) +
  labs(x = "log10(chromosome size, Mb)", y = "Mean GC content",
       colour = NULL, title = "GC composition retention") +
  theme_bw(base_size = 11)
save_plot(p_gc, "composition_GC_vs_size", w = 7, h = 5)

message("Done. All outputs written to ", OUT_DIR)
