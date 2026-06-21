############################################################
# Mouse bulk RNA-seq analysis
# FILM vs Non-FILM HCC tissues
# Analyses: focused GO cell-death GSEA + ImmunCellAI input
############################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(fgsea)
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

## =========================
## 1. Input files and output folders
## =========================

expr_file <- "Matrix.txt"      # mouse expression matrix: first column = gene symbol; columns = samples
geneset_file <- "GO.xls"       # mouse gene sets; each sheet = one pathway
outdir <- "Mouse_GSEA_GO_cell_death_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sample_group_file <- "mouse_sample_group.csv" # two columns required: sample, group

## =========================
## 2. Read expression matrix
## =========================

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

cat("Mouse matrix dimension:\n")
print(dim(expr_mat))
cat("Mouse samples:\n")
print(colnames(expr_mat))

sample_info <- fread(sample_group_file) |> as.data.frame()

required_cols <- c("sample", "group")
if (!all(required_cols %in% colnames(sample_info))) {
stop("sample_group_file must contain columns named: sample and group")
}

sample_info$sample <- as.character(sample_info$sample)
sample_info$group <- trimws(as.character(sample_info$group))

missing_samples <- setdiff(sample_info$sample, colnames(expr_mat))
if (length(missing_samples) > 0) {
stop("These samples are in sample_group_file but not in expression matrix: ",
paste(missing_samples, collapse = ", "))
}

extra_samples <- setdiff(colnames(expr_mat), sample_info$sample)
if (length(extra_samples) > 0) {
stop("These samples are in expression matrix but not in sample_group_file: ",
paste(extra_samples, collapse = ", "))
}

if (!all(sample_info$group %in% c("Non-FILM", "FILM"))) {
stop("The group column must contain only: Non-FILM and FILM")
}

expr_mat <- expr_mat[, sample_info$sample, drop = FALSE]

group <- sample_info$group
names(group) <- sample_info$sample

cat("Sample grouping used for analysis:\n")
print(table(group))
print(group)

nonfilm_samples <- names(group)[group == "Non-FILM"]
film_samples <- names(group)[group == "FILM"]

## If the matrix is TPM/FPKM/count-like values, log2(x + 1) is used here.
## If the matrix is already log-normalized, remove this transformation.
expr_log <- log2(expr_mat + 1)

mean_nonfilm <- rowMeans(expr_log[, nonfilm_samples, drop = FALSE])
mean_film <- rowMeans(expr_log[, film_samples, drop = FALSE])

## log2FC = FILM - Non-FILM
## NES > 0: enriched in FILM
## NES < 0: enriched in Non-FILM / suppressed in FILM
log2FC <- mean_film - mean_nonfilm
names(log2FC) <- rownames(expr_log)
ranks <- log2FC[is.finite(log2FC)]
ranks <- sort(ranks, decreasing = TRUE)

write.table(
  data.frame(Gene = names(ranks), RankMetric = as.numeric(ranks)),
  file = file.path(outdir, "Mouse_FILM_vs_NonFILM_log2FC.rnk"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

## =========================
## 3. Read gene sets
## =========================

sheet_names <- excel_sheets(geneset_file)
print(sheet_names)

pathways <- list()
for (s in sheet_names) {
  gs <- read_excel(geneset_file, sheet = s, col_names = FALSE)[[1]]
  gs <- unique(na.omit(as.character(gs)))
  gs <- gs[gs != ""]
  pathways[[s]] <- gs
}

cat("Gene set sizes:\n")
print(sapply(pathways, length))

overlap_info <- data.frame(
  pathway = names(pathways),
  geneset_size = sapply(pathways, length),
  overlap_with_matrix = sapply(pathways, function(x) sum(x %in% names(ranks)))
)
print(overlap_info)
write.csv(overlap_info, file.path(outdir, "GeneSet_overlap_with_matrix.csv"), row.names = FALSE)

## =========================
## 4. Run fgsea
## =========================

fgsea_res <- fgsea(
  pathways = pathways,
  stats = ranks,
  minSize = 5,
  maxSize = 500,
  nperm = 10000
)

fgsea_res <- fgsea_res %>%
  arrange(padj) %>%
  mutate(direction = ifelse(NES > 0, "Enriched in FILM", "Enriched in Non-FILM / suppressed in FILM"))

print(fgsea_res)
fwrite(fgsea_res, file.path(outdir, "Mouse_focused_GO_cell_death_GSEA_results.csv"))

leading_edge <- fgsea_res %>%
  select(pathway, leadingEdge) %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.csv(leading_edge, file.path(outdir, "Mouse_leading_edge_genes.csv"), row.names = FALSE)

## =========================
## 5. Summary dot plot
## =========================

plot_df <- fgsea_res %>%
  mutate(
    pathway = factor(pathway, levels = pathway[order(NES)]),
    neg_log10_padj = -log10(padj + 1e-300)
  )

p_dot <- ggplot(plot_df, aes(x = NES, y = pathway)) +
  geom_point(aes(size = neg_log10_padj, color = NES)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 12) +
  labs(
    x = "Normalized enrichment score (NES)",
    y = NULL,
    size = "-log10(FDR)",
    color = "NES",
    title = "Mouse focused GO cell death-related GSEA",
    subtitle = "FILM vs Non-FILM; NES < 0 indicates suppression in FILM"
  )

ggsave(file.path(outdir, "Mouse_focused_GO_cell_death_GSEA_dotplot.pdf"), p_dot, width = 7, height = 4)
ggsave(file.path(outdir, "Mouse_focused_GO_cell_death_GSEA_dotplot.png"), p_dot, width = 7, height = 4, dpi = 300)

## =========================
## 6. Single and combined enrichment plots
## =========================

single_plot_dir <- file.path(outdir, "single_GSEA_plots")
dir.create(single_plot_dir, showWarnings = FALSE)

pathway_label_map <- c(
  "Pyroptosis" = "PYROPTOSIS",
  "apoptosis" = "APOPTOSIS",
  "AUTOPHAGY" = "AUTOPHAGY",
  "NECROPTOTIC_PROCESS" = "NECROPTOTIC_PROCESS",
  "FERROPTOSIS" = "FERROPTOSIS"
)

plot_order <- c("AUTOPHAGY", "apoptosis", "FERROPTOSIS", "NECROPTOTIC_PROCESS", "Pyroptosis")
plot_list <- list()

for (pw in plot_order) {
  if (!pw %in% names(pathways)) {
    warning(paste("Pathway not found:", pw))
    next
  }
  
  res_pw <- fgsea_res %>% filter(pathway == pw)
  display_name <- ifelse(pw %in% names(pathway_label_map), pathway_label_map[[pw]], pw)
  
  label_text <- paste0(
    "ES=", round(res_pw$ES, 3), "\n",
    "NES=", round(res_pw$NES, 3), "\n",
    "FDR=", signif(res_pw$padj, 3)
  )
  
  p <- plotEnrichment(pathways[[pw]], ranks) +
    labs(title = display_name, x = "Rank in ordered dataset", y = "Enrichment score") +
    annotate("text", x = length(ranks) * 0.05, y = max(res_pw$ES, 0.05), label = label_text,
             hjust = 0, vjust = 1, size = 4) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          panel.grid.minor = element_blank())
  
  ggsave(file.path(single_plot_dir, paste0(display_name, "_Mouse_GSEA.pdf")), p, width = 6, height = 4)
  ggsave(file.path(single_plot_dir, paste0(display_name, "_Mouse_GSEA.png")), p, width = 6, height = 4, dpi = 300)
  
  plot_list[[display_name]] <- p +
    labs(x = NULL, y = "Enrichment Score") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

if (length(plot_list) > 0) {
  combined_plot <- wrap_plots(plot_list, nrow = 1) +
    plot_annotation(
      title = "HCC tissues from mouse",
      subtitle = "FILM vs Non-FILM; negative NES indicates suppression in FILM"
    )
  
  ggsave(file.path(outdir, "Mouse_focused_GO_cell_death_GSEA_combined.pdf"), combined_plot, width = 15, height = 3.5)
  ggsave(file.path(outdir, "Mouse_focused_GO_cell_death_GSEA_combined.png"), combined_plot, width = 15, height = 3.5, dpi = 300)
}

## =========================
## 7. Export mouse expression matrix for ImmunCellAI
## =========================

mouse_deconv_outdir <- "Mouse_immune_deconvolution_results"
dir.create(mouse_deconv_outdir, showWarnings = FALSE, recursive = TRUE)

mouse_anno <- data.frame(
  sample = colnames(expr_mat),
  group = group,
  row.names = colnames(expr_mat),
  check.names = FALSE
)
write.csv(mouse_anno, file.path(mouse_deconv_outdir, "mouse_sample_group.csv"), row.names = FALSE)

mouse_immunecellai_input <- data.frame(
  GeneSymbol = rownames(expr_mat),
  expr_mat,
  check.names = FALSE
)
write.csv(mouse_immunecellai_input, file.path(mouse_deconv_outdir, "mouse_expression_for_ImmunCellAI.csv"), row.names = FALSE)

## If ImmunCellAI results are generated externally, place the file below and re-run this section.
mouse_immunecellai_result_file <- file.path(mouse_deconv_outdir, "mouse_ImmunCellAI_results.csv")
if (file.exists(mouse_immunecellai_result_file)) {
  mouse_immunecellai <- read.csv(mouse_immunecellai_result_file, check.names = FALSE)
  write.csv(mouse_immunecellai, file.path(mouse_deconv_outdir, "mouse_ImmunCellAI_results_imported.csv"), row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_mouse_bulk_RNAseq_GSEA.txt"))
message("Mouse bulk RNA-seq GSEA analysis completed.")
