# Human bulk RNA-seq analysis
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(fgsea)
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(xCell)
  library(MCPcounter)
  library(immunedeconv)
})

expr_file <- "Matrix.txt"
sample_group_file <- "human_sample_group.csv"
go_geneset_file <- "GO_human.xls"
kegg_geneset_file <- "KEGG_human.xls"

outdir <- "Human_bulk_RNAseq_results"
gsea_outdir <- file.path(outdir, "GO_cell_death_GSEA")
kegg_outdir <- file.path(outdir, "KEGG_pathway_GSEA")
immune_outdir <- file.path(outdir, "Immune_deconvolution")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(gsea_outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(kegg_outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(immune_outdir, showWarnings = FALSE, recursive = TRUE)

expr <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gene_col <- colnames(expr)[1]

expr_mat <- expr %>%
  distinct(.data[[gene_col]], .keep_all = TRUE) %>%
  column_to_rownames(gene_col)

expr_mat <- as.data.frame(expr_mat)
expr_mat[] <- lapply(expr_mat, function(x) as.numeric(as.character(x)))
expr_mat <- as.matrix(expr_mat)
expr_mat <- expr_mat[rowSums(is.na(expr_mat)) == 0, , drop = FALSE]
expr_mat <- expr_mat[rowSums(expr_mat) > 0, , drop = FALSE]

sample_info <- fread(sample_group_file) |> as.data.frame()

if (!all(c("sample", "group") %in% colnames(sample_info))) {
  stop("sample_group_file must contain columns named: sample and group")
}

sample_info$sample <- as.character(sample_info$sample)
sample_info$group <- trimws(as.character(sample_info$group))

missing_samples <- setdiff(sample_info$sample, colnames(expr_mat))
if (length(missing_samples) > 0) {
  stop("Samples in sample_group_file but not in expression matrix: ",
       paste(missing_samples, collapse = ", "))
}

extra_samples <- setdiff(colnames(expr_mat), sample_info$sample)
if (length(extra_samples) > 0) {
  stop("Samples in expression matrix but not in sample_group_file: ",
       paste(extra_samples, collapse = ", "))
}

if (!all(sample_info$group %in% c("Non-FILM", "FILM"))) {
  stop("The group column must contain only: Non-FILM and FILM")
}

expr_mat <- expr_mat[, sample_info$sample, drop = FALSE]

group <- sample_info$group
names(group) <- sample_info$sample

nonfilm_samples <- names(group)[group == "Non-FILM"]
film_samples <- names(group)[group == "FILM"]

expr_log <- log2(expr_mat + 1)

mean_nonfilm <- rowMeans(expr_log[, nonfilm_samples, drop = FALSE])
mean_film <- rowMeans(expr_log[, film_samples, drop = FALSE])

rank_metric <- mean_film - mean_nonfilm
names(rank_metric) <- rownames(expr_log)
ranks <- rank_metric[is.finite(rank_metric)]
ranks <- sort(ranks, decreasing = TRUE)

write.table(
  data.frame(Gene = names(ranks), RankMetric = as.numeric(ranks)),
  file = file.path(outdir, "Human_FILM_vs_NonFILM.rnk"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

read_local_genesets <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("Gene-set file not found: ", file_path)
  }

  sheet_names <- excel_sheets(file_path)
  pathways <- list()

  for (s in sheet_names) {
    gs <- read_excel(file_path, sheet = s, col_names = FALSE)[[1]]
    gs <- unique(na.omit(as.character(gs)))
    gs <- trimws(gs)
    gs <- gs[gs != ""]
    pathways[[s]] <- gs
  }

  pathways
}

run_gsea <- function(pathways, ranks, outdir, prefix) {
  overlap_info <- data.frame(
    pathway = names(pathways),
    geneset_size = sapply(pathways, length),
    overlap_with_matrix = sapply(pathways, function(x) sum(x %in% names(ranks)))
  )

  write.csv(
    overlap_info,
    file.path(outdir, paste0(prefix, "_geneset_overlap.csv")),
    row.names = FALSE
  )

  set.seed(123)

  res <- fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = 5,
    maxSize = 500,
    nperm = 10000
  )

  res <- res %>%
    arrange(padj) %>%
    mutate(
      direction = ifelse(
        NES > 0,
        "Enriched in FILM",
        "Enriched in Non-FILM"
      )
    )

  res_export <- res %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

  fwrite(
    res_export,
    file.path(outdir, paste0(prefix, "_results.csv"))
  )

  write.csv(
    res_export %>% select(pathway, leadingEdge),
    file.path(outdir, paste0(prefix, "_leading_edge_genes.csv")),
    row.names = FALSE
  )

  res
}

plot_dotplot <- function(res, outdir, prefix, title, top_n = 20) {
  plot_df <- res %>%
    arrange(padj) %>%
    slice_head(n = top_n) %>%
    mutate(
      pathway = factor(pathway, levels = rev(pathway)),
      padj_plot = ifelse(is.na(padj), 1, padj),
      neg_log10_padj = -log10(padj_plot + 1e-300)
    )

  p <- ggplot(plot_df, aes(x = NES, y = pathway)) +
    geom_point(aes(size = neg_log10_padj, color = NES)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_bw(base_size = 12) +
    labs(
      x = "Normalized enrichment score",
      y = NULL,
      size = "-log10(FDR)",
      color = "NES",
      title = title,
      subtitle = "FILM vs Non-FILM"
    )

  ggsave(file.path(outdir, paste0(prefix, "_dotplot.pdf")), p, width = 8, height = 6)
  ggsave(file.path(outdir, paste0(prefix, "_dotplot.png")), p, width = 8, height = 6, dpi = 300)

  p
}

plot_single_gsea <- function(pathways, res, ranks, outdir, selected_pathways, prefix) {
  single_dir <- file.path(outdir, "single_GSEA_plots")
  dir.create(single_dir, showWarnings = FALSE, recursive = TRUE)

  plot_list <- list()

  for (pw in selected_pathways) {
    if (!pw %in% names(pathways)) {
      next
    }

    res_pw <- res %>% filter(pathway == pw)

    if (nrow(res_pw) == 0) {
      next
    }

    label_text <- paste0(
      "ES=", round(res_pw$ES, 3), "\n",
      "NES=", round(res_pw$NES, 3), "\n",
      "FDR=", signif(res_pw$padj, 3)
    )

    p <- plotEnrichment(pathways[[pw]], ranks) +
      labs(title = pw, x = "Rank in ordered dataset", y = "Enrichment score") +
      annotate(
        "text",
        x = length(ranks) * 0.05,
        y = max(res_pw$ES, 0.05),
        label = label_text,
        hjust = 0,
        vjust = 1,
        size = 4
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        panel.grid.minor = element_blank()
      )

    ggsave(file.path(single_dir, paste0(make.names(pw), "_GSEA.pdf")), p, width = 6, height = 4)
    ggsave(file.path(single_dir, paste0(make.names(pw), "_GSEA.png")), p, width = 6, height = 4, dpi = 300)

    plot_list[[pw]] <- p +
      labs(x = NULL, y = "Enrichment score") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }

  if (length(plot_list) > 0) {
    combined_plot <- wrap_plots(plot_list, nrow = 1) +
      plot_annotation(title = prefix)

    ggsave(file.path(outdir, paste0(prefix, "_combined.pdf")), combined_plot, width = 15, height = 3.5)
    ggsave(file.path(outdir, paste0(prefix, "_combined.png")), combined_plot, width = 15, height = 3.5, dpi = 300)
  }
}

go_pathways <- read_local_genesets(go_geneset_file)

go_res <- run_gsea(
  pathways = go_pathways,
  ranks = ranks,
  outdir = gsea_outdir,
  prefix = "GO_cell_death_GSEA"
)

plot_dotplot(
  res = go_res,
  outdir = gsea_outdir,
  prefix = "GO_cell_death_GSEA",
  title = "GO cell death-related GSEA",
  top_n = nrow(go_res)
)

cell_death_order <- c(
  "AUTOPHAGY",
  "APOPTOSIS",
  "FERROPTOSIS",
  "NECROPTOTIC_PROCESS",
  "PYROPTOSIS"
)


cell_death_order <- unique(cell_death_order[cell_death_order %in% names(go_pathways)])

plot_single_gsea(
  pathways = go_pathways,
  res = go_res,
  ranks = ranks,
  outdir = gsea_outdir,
  selected_pathways = cell_death_order,
  prefix = "GO_cell_death_GSEA"
)

kegg_pathways <- read_local_genesets(kegg_geneset_file)

kegg_res <- run_gsea(
  pathways = kegg_pathways,
  ranks = ranks,
  outdir = kegg_outdir,
  prefix = "KEGG_pathway_GSEA"
)

plot_dotplot(
  res = kegg_res,
  outdir = kegg_outdir,
  prefix = "KEGG_pathway_GSEA",
  title = "KEGG pathway GSEA",
  top_n = 20
)

ppar_pathway <- names(kegg_pathways)[grepl("ppar", names(kegg_pathways), ignore.case = TRUE)]

if (length(ppar_pathway) > 0) {
  plot_single_gsea(
    pathways = kegg_pathways,
    res = kegg_res,
    ranks = ranks,
    outdir = kegg_outdir,
    selected_pathways = ppar_pathway[1],
    prefix = "PPAR_signaling_pathway"
  )
}

xcell_res <- tryCatch(
  {
    res <- xCellAnalysis(expr_mat)
    res <- as.data.frame(res)
    res$cell_type <- rownames(res)
    res <- res %>% relocate(cell_type)
    fwrite(res, file.path(immune_outdir, "xCell_scores.csv"))
    res
  },
  error = function(e) {
    message("xCell failed: ", e$message)
    NULL
  }
)

mcp_res <- tryCatch(
  {
    res <- MCPcounter.estimate(expr_mat, featuresType = "HUGO_symbols")
    res <- as.data.frame(res)
    res$cell_type <- rownames(res)
    res <- res %>% relocate(cell_type)
    fwrite(res, file.path(immune_outdir, "MCPcounter_scores.csv"))
    res
  },
  error = function(e) {
    message("MCPcounter failed: ", e$message)
    NULL
  }
)

quantiseq_res <- tryCatch(
  {
    res <- deconvolute(expr_mat, method = "quantiseq")
    res <- as.data.frame(res)
    colnames(res)[1] <- "cell_type"
    fwrite(res, file.path(immune_outdir, "quanTIseq_scores.csv"))
    res
  },
  error = function(e) {
    message("quanTIseq failed: ", e$message)
    NULL
  }
)

to_long <- function(res, method) {
  res %>%
    pivot_longer(
      cols = -cell_type,
      names_to = "sample",
      values_to = "score"
    ) %>%
    mutate(method = method)
}

deconv_list <- list()

if (!is.null(xcell_res)) {
  deconv_list[["xCell"]] <- to_long(xcell_res, "xCell")
}

if (!is.null(mcp_res)) {
  deconv_list[["MCPcounter"]] <- to_long(mcp_res, "MCPcounter")
}

if (!is.null(quantiseq_res)) {
  deconv_list[["quanTIseq"]] <- to_long(quantiseq_res, "quanTIseq")
}

if (length(deconv_list) > 0) {
  deconv_long <- bind_rows(deconv_list) %>%
    left_join(sample_info, by = "sample")

  fwrite(
    deconv_long,
    file.path(immune_outdir, "immune_deconvolution_long.csv")
  )

  deconv_summary <- deconv_long %>%
    group_by(method, cell_type, group) %>%
    summarize(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

  fwrite(
    deconv_summary,
    file.path(immune_outdir, "immune_deconvolution_group_mean.csv")
  )

  top_cells <- deconv_long %>%
    group_by(method, cell_type) %>%
    summarize(v = var(score, na.rm = TRUE), .groups = "drop") %>%
    group_by(method) %>%
    slice_max(v, n = 20, with_ties = FALSE) %>%
    ungroup()

  plot_df <- deconv_long %>%
    semi_join(top_cells, by = c("method", "cell_type"))

  p_immune <- ggplot(plot_df, aes(x = group, y = score, fill = group)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.25) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.7) +
    facet_grid(method ~ cell_type, scales = "free_y") +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.x = element_text(size = 5),
      strip.text.y = element_text(size = 7),
      legend.position = "none"
    ) +
    labs(
      x = NULL,
      y = "Estimated abundance",
      title = "Immune-cell deconvolution"
    )

  ggsave(
    file.path(immune_outdir, "immune_deconvolution_boxplot.pdf"),
    p_immune,
    width = 18,
    height = 8
  )

  ggsave(
    file.path(immune_outdir, "immune_deconvolution_boxplot.png"),
    p_immune,
    width = 18,
    height = 8,
    dpi = 300
  )
}

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "sessionInfo.txt")
)

message("Human bulk RNA-seq analysis completed.")
