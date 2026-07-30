#!/usr/bin/env Rscript
# ============================================================
# qPCR Analyzer - R / ggplot2 version
# ============================================================
# Same method as the Python engine (ddCt / Livak method), for users
# who prefer R output. Produces the same fold-change + stats + figure.
#
# Requires: tidyverse, rstatix, ggpubr
#   install.packages(c("tidyverse", "rstatix", "ggpubr"))
#
# Usage:
#   Rscript qpcr_analysis.R input.csv --reference GAPDH --control Control --output figure.pdf
#
# Input CSV format (long): Sample, Group, Gene, Ct   (see sample_data.csv)
# ============================================================

suppressMessages({
  library(tidyverse)
  library(rstatix)
  library(ggpubr)
})

# ---------------------------------------------------------------- CLI args
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  return(args[idx + 1])
}

input_file    <- args[1]
reference_gene <- get_arg("--reference", "GAPDH")
control_group  <- get_arg("--control", "Control")
output_file    <- get_arg("--output", "qpcr_figure.pdf")
alpha          <- as.numeric(get_arg("--alpha", "0.05"))

if (is.na(input_file) || is.null(input_file)) {
  stop("Usage: Rscript qpcr_analysis.R input.csv --reference GENE --control GROUP --output figure.pdf")
}

# ---------------------------------------------------------------- Load & ddCt
raw <- read_csv(input_file, show_col_types = FALSE)

ref_ct <- raw %>%
  filter(Gene == reference_gene) %>%
  select(Sample, Ct_ref = Ct)

target_genes <- setdiff(unique(raw$Gene), reference_gene)

ddct <- raw %>%
  filter(Gene %in% target_genes) %>%
  left_join(ref_ct, by = "Sample") %>%
  mutate(dCt = Ct - Ct_ref) %>%
  group_by(Gene) %>%
  mutate(control_mean_dct = mean(dCt[Group == control_group])) %>%
  ungroup() %>%
  mutate(
    ddCt = dCt - control_mean_dct,
    fold_change = 2^(-ddCt),
    log2_fold_change = log2(fold_change)
  )

write_csv(ddct, "qpcr_ddct_R.csv")

# ---------------------------------------------------------------- Stats per gene
stats_results <- list()
plots <- list()

for (g in unique(ddct$Gene)) {
  sub <- ddct %>% filter(Gene == g)
  n_groups <- n_distinct(sub$Group)

  if (n_groups == 2) {
    test_res <- sub %>% t_test(dCt ~ Group, var.equal = FALSE) %>%
      add_significance("p")
    test_res$test <- "Welch t-test"
  } else {
    anova_res <- sub %>% anova_test(dCt ~ Group)
    print(anova_res)
    test_res <- sub %>% tukey_hsd(dCt ~ Group) %>%
      add_significance("p.adj")
    test_res$test <- "Tukey HSD"
  }
  stats_results[[g]] <- test_res

  # ---------------- Plot: customize this block for your own formatting ----------------
  summary_df <- sub %>%
    group_by(Group) %>%
    summarise(mean_fc = mean(fold_change), sem_fc = sd(fold_change) / sqrt(n()), .groups = "drop")

  p <- ggplot(summary_df, aes(x = Group, y = mean_fc, fill = Group)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.4) +
    geom_errorbar(aes(ymin = mean_fc - sem_fc, ymax = mean_fc + sem_fc), width = 0.15) +
    geom_jitter(data = sub, aes(x = Group, y = fold_change),
                width = 0.08, size = 1.8, alpha = 0.7, color = "black", inherit.aes = FALSE) +
    scale_fill_manual(values = c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2")) +
    labs(title = g, x = NULL, y = "Relative expression (fold change)") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

  # add significance brackets if any comparisons are significant
  sig_pairs <- test_res %>% filter(if (n_groups == 2) p < alpha else p.adj < alpha)
  if (nrow(sig_pairs) > 0) {
    label_col <- if (n_groups == 2) "p.signif" else "p.adj.signif"
    p <- p + stat_pvalue_manual(
      sig_pairs %>% mutate(y.position = max(sub$fold_change) * seq(1.15, 1.6, length.out = nrow(sig_pairs))[1:nrow(sig_pairs)]),
      label = label_col
    )
  }

  plots[[g]] <- p
}

# ---------------------------------------------------------------- Save outputs
combined_plot <- ggarrange(plotlist = plots, ncol = length(plots))
ggsave(output_file, combined_plot, width = 5 * length(plots), height = 4.5, dpi = 300)

bind_rows(stats_results, .id = "Gene") %>% write_csv("qpcr_stats_R.csv")

cat("Done.\n")
cat("  Figure ->", output_file, "\n")
cat("  ddCt table -> qpcr_ddct_R.csv\n")
cat("  Stats table -> qpcr_stats_R.csv\n")
