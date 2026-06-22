############################################################
# CD45+ immune-cell scRNA-seq analysis
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(DoubletFinder)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(Matrix)
  library(qs)
  library(future)
  library(pheatmap)
  library(scales)
})

# =========================
# 0. Paths, samples and parameters
# =========================
workdir <- "."
setwd(workdir)

outdir <- "results_FigS20_CD45_M1_M8"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots", "QC"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots", "CD45_branch"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots", "T_branch"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots", "TAM_branch"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots", "DC_branch"), showWarnings = FALSE, recursive = TRUE)

force_rebuild <- TRUE

raw_rds_all     <- file.path(outdir, "sce_raw_merged_M1_M8_geneSymbol.rds")
qc_rds_all      <- file.path(outdir, "sce_after_qc_M1_M8_geneSymbol.rds")
df_rds_all      <- file.path(outdir, "sce_after_DoubletFinder_M1_M8_geneSymbol.rds")
singlet_rds_all <- file.path(outdir, "sce_singlet_M1_M8_geneSymbol.rds")
cd45_rds_final  <- file.path(outdir, "sce_cd45_harmony_sample_annotated_M1_M8.qs")

# Optional sample metadata file; recommended for manuscript upload.
# Format: sample,folder,group,batch
# If absent, default group mapping below is used.
sample_info_file <- "sample_info_FigS20.csv"

if (!file.exists(sample_info_file)) {
  stop(
    "sample_info_FigS20.csv not found. ",
    "Please provide a sample metadata file with columns: sample, folder, group. ",
    "Optional column: batch."
  )
}

sample_info <- read.csv(sample_info_file, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("sample", "folder", "group")
if (!all(required_cols %in% colnames(sample_info))) {
  stop(
    "sample_info_FigS20.csv must contain columns: ",
    paste(required_cols, collapse = ", ")
  )
}

if (!"batch" %in% colnames(sample_info)) {
  sample_info$batch <- "batch1"
}

sample_dirs <- sample_info$folder
names(sample_dirs) <- sample_info$sample

sample_order <- sample_info$sample
cd45_samples <- sample_info$sample

sample_group_map <- setNames(sample_info$group, sample_info$sample)
sample_batch_map <- setNames(sample_info$batch, sample_info$sample)

group_order <- unique(unname(sample_group_map[sample_order]))

# QC parameters for mouse CD45+ TIIC scRNA-seq
min_features <- 200
max_features <- 7500
max_percent_mt <- 15

# Global CD45 clustering parameters
npcs_cd45 <- 30
dims_cd45 <- 1:15
res_cd45 <- 0.4
k_cd45 <- 20

# T-cell reclustering parameters
npcs_t <- 30
dims_t <- 1:10
res_t <- 0.15
k_t <- 20

# Myeloid/TAM and DC reclustering parameters
npcs_myeloid <- 30
dims_myeloid <- 1:15
res_myeloid <- 0.4
k_myeloid <- 20

npcs_tam <- 30
dims_tam <- 1:10
res_tam <- 0.1
k_tam <- 15

npcs_dc <- 20
dims_dc <- 1:10
res_dc <- 0.03
k_dc <- 10

# =========================
# 1. Utility functions
# =========================
read10x_auto_symbol <- function(folder_path) {
  feat_file <- file.path(folder_path, "features.tsv.gz")
  if (!file.exists(feat_file)) feat_file <- file.path(folder_path, "genes.tsv.gz")
  if (!file.exists(feat_file)) stop("Cannot find features.tsv.gz or genes.tsv.gz in: ", folder_path)
  feat <- read.delim(feat_file, header = FALSE, stringsAsFactors = FALSE)
  gene_col <- if (ncol(feat) >= 2) 2 else 1
  message("  using gene.column = ", gene_col)
  counts <- Read10X(data.dir = folder_path, gene.column = gene_col)
  rownames(counts) <- make.unique(rownames(counts))
  counts
}

join_if_needed <- function(obj, assay = "RNA") {
  DefaultAssay(obj) <- assay
  if (inherits(obj[[assay]], "Assay5")) {
    lyr_now <- Layers(obj[[assay]])
    need_join <- sum(grepl("^counts", lyr_now)) > 1 ||
      sum(grepl("^data", lyr_now)) > 1 ||
      sum(grepl("^scale.data", lyr_now)) > 1
    if (need_join) obj[[assay]] <- JoinLayers(obj[[assay]])
  }
  obj
}

get_counts_for_one_sample <- function(obj, sample_name, assay = "RNA") {
  DefaultAssay(obj) <- assay
  if (inherits(obj[[assay]], "Assay5")) {
    all_layers <- Layers(obj[[assay]])
    lyr <- grep("^counts", all_layers, value = TRUE)
    if (length(lyr) == 0) stop("No counts layer found for sample: ", sample_name)
    if (length(lyr) > 1) {
      lyr2 <- grep(paste0("^counts\\.", sample_name, "$"), lyr, value = TRUE)
      if (length(lyr2) == 1) {
        lyr <- lyr2
      } else if ("counts" %in% lyr) {
        lyr <- "counts"
      } else {
        stop("Multiple counts layers found for sample ", sample_name, ": ", paste(lyr, collapse = ", "))
      }
    }
    LayerData(obj[[assay]], layer = lyr)
  } else {
    GetAssayData(obj, assay = assay, slot = "counts")
  }
}

ensure_percent_mt <- function(obj, assay = "RNA") {
  DefaultAssay(obj) <- assay
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-|^Mt-")
  }
  obj$percent.mt[is.na(obj$percent.mt)] <- 0
  obj
}

save_plot <- function(plot, filename, width = 7, height = 6, dpi = 300) {
  ggsave(filename, plot, width = width, height = height, useDingbats = FALSE)
  if (grepl("\\.pdf$", filename)) {
    ggsave(sub("\\.pdf$", ".png", filename), plot, width = width, height = height, dpi = dpi)
  }
}

present_genes <- function(genes, obj) genes[genes %in% rownames(obj)]

add_module_scores_named <- function(obj, gene_sets, min_genes = 2) {
  gene_sets <- lapply(gene_sets, function(x) intersect(x, rownames(obj)))
  gene_sets <- gene_sets[lengths(gene_sets) >= min_genes]
  score_cols <- character(0)
  for (nm in names(gene_sets)) {
    old_cols <- colnames(obj@meta.data)
    obj <- AddModuleScore(obj, features = list(gene_sets[[nm]]), name = paste0(nm, "_Score"))
    new_col <- setdiff(colnames(obj@meta.data), old_cols)
    if (length(new_col) != 1) stop("Unexpected AddModuleScore columns for ", nm)
    colnames(obj@meta.data)[colnames(obj@meta.data) == new_col] <- nm
    score_cols <- c(score_cols, nm)
  }
  attr(obj, "module_score_cols") <- score_cols
  obj
}

plot_module_heatmap <- function(obj, group_col, score_cols, outfile, title = "Module scores") {
  score_df <- obj@meta.data %>%
    dplyr::filter(!is.na(.data[[group_col]])) %>%
    dplyr::select(dplyr::all_of(c(group_col, score_cols)))
  mean_score <- score_df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(score_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  write.csv(mean_score, sub("\\.pdf$", "_mean_scores.csv", outfile), row.names = FALSE)
  mean_mat <- as.data.frame(mean_score)
  rownames(mean_mat) <- mean_mat[[group_col]]
  mean_mat[[group_col]] <- NULL
  mat <- t(as.matrix(mean_mat))
  z_mat <- t(scale(t(mat)))
  z_mat[is.na(z_mat)] <- 0
  write.csv(mat, sub("\\.pdf$", "_raw_matrix.csv", outfile))
  write.csv(z_mat, sub("\\.pdf$", "_zscore_matrix.csv", outfile))
  pdf(outfile, width = 9, height = 6)
  pheatmap(z_mat, cluster_rows = FALSE, cluster_cols = FALSE,
           fontsize_row = 10, fontsize_col = 10, angle_col = 45,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           main = title)
  dev.off()
}

plot_composition <- function(meta, subtype_col, out_prefix, subtype_order = NULL, subtype_cols = NULL) {
  meta <- meta %>% dplyr::filter(!is.na(.data[[subtype_col]]))
  meta$sample <- factor(as.character(meta$sample), levels = sample_order)
  meta$sample_group <- unname(sample_group_map[as.character(meta$sample)])
  meta$sample_group <- factor(meta$sample_group, levels = group_order)
  if (!is.null(subtype_order)) meta[[subtype_col]] <- factor(meta[[subtype_col]], levels = subtype_order)
  
  # by sample
  df_sample <- meta %>%
    dplyr::group_by(sample, .data[[subtype_col]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
  colnames(df_sample)[2] <- "subtype"
  write.csv(df_sample, paste0(out_prefix, "_composition_by_sample.csv"), row.names = FALSE)
  p_sample <- ggplot(df_sample, aes(x = sample, y = prop, fill = subtype)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(x = "Sample", y = "Proportion", fill = subtype_col, title = paste0(subtype_col, " composition by sample")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
  if (!is.null(subtype_cols)) p_sample <- p_sample + scale_fill_manual(values = subtype_cols, drop = FALSE)
  save_plot(p_sample, paste0(out_prefix, "_composition_by_sample.pdf"), width = 10, height = 6)
  
  # by group
  df_group <- meta %>%
    dplyr::filter(!is.na(sample_group)) %>%
    dplyr::group_by(sample_group, .data[[subtype_col]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(sample_group) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
  colnames(df_group)[2] <- "subtype"
  write.csv(df_group, paste0(out_prefix, "_composition_by_group.csv"), row.names = FALSE)
  p_group <- ggplot(df_group, aes(x = sample_group, y = prop, fill = subtype)) +
    geom_bar(stat = "identity", width = 0.75) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(x = "Group", y = "Proportion", fill = subtype_col, title = paste0(subtype_col, " composition by group")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
  if (!is.null(subtype_cols)) p_group <- p_group + scale_fill_manual(values = subtype_cols, drop = FALSE)
  save_plot(p_group, paste0(out_prefix, "_composition_by_group.pdf"), width = 8, height = 6)
  
  p_facet <- ggplot(df_group, aes(x = sample_group, y = prop, fill = subtype)) +
    geom_bar(stat = "identity", width = 0.75) + facet_wrap(~ subtype, nrow = 1) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(x = "Group", y = "Proportion", fill = subtype_col, title = paste0(subtype_col, " by group")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), plot.title = element_text(hjust = 0.5), strip.background = element_rect(fill = "grey90")) +
    guides(fill = "none")
  if (!is.null(subtype_cols)) p_facet <- p_facet + scale_fill_manual(values = subtype_cols, drop = FALSE)
  save_plot(p_facet, paste0(out_prefix, "_composition_by_group_facet.pdf"), width = 13, height = 5)
}

# =========================
# 2. Read 10X matrices and build raw merged object
# =========================
if (!force_rebuild && file.exists(raw_rds_all)) {
  sce <- readRDS(raw_rds_all)
} else {
  sc_list <- list()
  for (sample_name in names(sample_dirs)) {
    message("Loading ", sample_name, " from ", sample_dirs[[sample_name]])
    counts <- read10x_auto_symbol(sample_dirs[[sample_name]])
    obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
    obj$sample <- sample_name
    obj$sample_group <- unname(sample_group_map[sample_name])
    obj$batch <- unname(sample_batch_map[sample_name])
    sc_list[[sample_name]] <- obj
  }
  sce <- merge(sc_list[[1]], y = sc_list[-1], add.cell.ids = names(sc_list))
  saveRDS(sce, raw_rds_all)
  rm(sc_list); gc()
}

message("Raw cells: ", ncol(sce), " | genes: ", nrow(sce))
message("ENSMUSG-like rownames: ", sum(grepl("^ENSMUSG", rownames(sce))))

# =========================
# 3. QC and DoubletFinder
# =========================
sce <- ensure_percent_mt(sce)

qc_summary_before <- sce@meta.data %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(n_cells = dplyr::n(), median_nFeature_RNA = median(nFeature_RNA), median_nCount_RNA = median(nCount_RNA), median_percent_mt = median(percent.mt), q95_nCount_RNA = quantile(nCount_RNA, 0.95), q99_nCount_RNA = quantile(nCount_RNA, 0.99), .groups = "drop")
write.csv(qc_summary_before, file.path(outdir, "QC_summary_before_filtering_by_sample.csv"), row.names = FALSE)

p_qc_before <- VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sample", pt.size = 0)
save_plot(p_qc_before, file.path(outdir, "plots", "QC", "QC_VlnPlot_before_filtering.pdf"), width = 14, height = 6)

sce <- subset(sce, subset = nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mt < max_percent_mt)
saveRDS(sce, qc_rds_all)

qc_summary_after <- sce@meta.data %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(n_cells = dplyr::n(), median_nFeature_RNA = median(nFeature_RNA), median_nCount_RNA = median(nCount_RNA), median_percent_mt = median(percent.mt), q95_nCount_RNA = quantile(nCount_RNA, 0.95), q99_nCount_RNA = quantile(nCount_RNA, 0.99), .groups = "drop")
write.csv(qc_summary_after, file.path(outdir, "QC_summary_after_filtering_by_sample.csv"), row.names = FALSE)

p_qc_after <- VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sample", pt.size = 0)
save_plot(p_qc_after, file.path(outdir, "plots", "QC", "QC_VlnPlot_after_filtering.pdf"), width = 14, height = 6)

# Per-sample DoubletFinder
sce_list <- SplitObject(sce, split.by = "sample")
df_list <- list()
for (i in names(sce_list)) {
  message("DoubletFinder for sample: ", i)
  tmp5 <- sce_list[[i]]
  counts_mat <- get_counts_for_one_sample(tmp5, sample_name = i, assay = "RNA")
  meta <- tmp5@meta.data[colnames(counts_mat), , drop = FALSE]
  tmp <- CreateSeuratObject(counts = counts_mat, meta.data = meta, project = i, min.cells = 0, min.features = 0)
  tmp <- NormalizeData(tmp, verbose = FALSE)
  tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  tmp <- ScaleData(tmp, verbose = FALSE)
  tmp <- RunPCA(tmp, npcs = 30, verbose = FALSE)
  df_dims <- 1:20
  tmp <- FindNeighbors(tmp, dims = df_dims, verbose = FALSE)
  tmp <- FindClusters(tmp, resolution = 0.4, verbose = FALSE)
  
  df_exports <- getNamespaceExports("DoubletFinder")
  paramSweep_fun <- if ("paramSweep_v3" %in% df_exports) DoubletFinder::paramSweep_v3 else DoubletFinder::paramSweep
  doublet_fun <- if ("doubletFinder_v3" %in% df_exports) DoubletFinder::doubletFinder_v3 else DoubletFinder::doubletFinder
  use_doubletFinder_v3 <- "doubletFinder_v3" %in% df_exports
  
  sweep.res.list <- paramSweep_fun(tmp, PCs = df_dims, sct = FALSE)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  bc_col <- if ("BCmetric" %in% colnames(bcmvn)) "BCmetric" else colnames(bcmvn)[grepl("BC", colnames(bcmvn))][1]
  pK_value <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn[[bc_col]])]))
  if (is.na(pK_value)) stop("pK_value is NA for sample ", i)
  
  doublet_rate <- min((ncol(tmp) / 1000) * 0.008, 0.15)
  nExp_poi <- max(round(doublet_rate * ncol(tmp)), 1)
  homotypic.prop <- DoubletFinder::modelHomotypic(as.character(Idents(tmp)))
  nExp_poi.adj <- max(round(nExp_poi * (1 - homotypic.prop)), 1)
  
  if (use_doubletFinder_v3) {
    tmp <- doublet_fun(tmp, PCs = df_dims, pN = 0.25, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  } else {
    tmp <- doublet_fun(tmp, PCs = df_dims, pN = 0.25, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
  }
  reuse_pann <- tail(grep("^pANN_", colnames(tmp@meta.data), value = TRUE), 1)
  if (use_doubletFinder_v3) {
    tmp <- doublet_fun(tmp, PCs = df_dims, pN = 0.25, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = reuse_pann, sct = FALSE)
  } else {
    tmp <- doublet_fun(tmp, PCs = df_dims, pN = 0.25, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = reuse_pann)
  }
  df_col <- tail(grep("DF.classifications", colnames(tmp@meta.data), value = TRUE), 1)
  tmp$Doublet_Class <- tmp@meta.data[[df_col]]
  df_list[[i]] <- tmp
}

sce_all <- merge(df_list[[1]], y = df_list[-1], add.cell.ids = names(df_list))
sce_all <- join_if_needed(sce_all, assay = "RNA")
write.csv(as.data.frame.matrix(table(sce_all$Doublet_Class, sce_all$sample)), file.path(outdir, "DoubletFinder_summary.csv"))
saveRDS(sce_all, df_rds_all)

sce_all <- subset(sce_all, subset = Doublet_Class == "Singlet")
sce_all$sample <- factor(as.character(sce_all$sample), levels = sample_order)
sce_all$sample_group <- unname(sample_group_map[as.character(sce_all$sample)])
sce_all$batch <- unname(sample_batch_map[as.character(sce_all$sample)])
saveRDS(sce_all, singlet_rds_all)
qs::qsave(sce_all, file.path(outdir, "sce_singlet_M1_M8_geneSymbol.qs"), nthreads = 8)

# =========================
# 4. CD45+ global immune-cell analysis
# =========================
sce_cd45 <- subset(sce_all, subset = sample %in% cd45_samples)
sce_cd45 <- join_if_needed(sce_cd45, assay = "RNA")
sce_cd45 <- ensure_percent_mt(sce_cd45)

sce_cd45 <- NormalizeData(sce_cd45, verbose = FALSE)
sce_cd45 <- FindVariableFeatures(sce_cd45, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sce_cd45 <- ScaleData(sce_cd45, features = VariableFeatures(sce_cd45), vars.to.regress = "percent.mt", verbose = FALSE)
sce_cd45 <- RunPCA(sce_cd45, npcs = npcs_cd45, verbose = FALSE)
pdf(file.path(outdir, "plots", "CD45_branch", "CD45_ElbowPlot_M1_M8.pdf"), width = 7, height = 5); print(ElbowPlot(sce_cd45, ndims = npcs_cd45)); dev.off()

sce_cd45 <- RunHarmony(sce_cd45, group.by.vars = "sample", reduction.use = "pca", dims.use = dims_cd45, reduction.save = "harmony_sample", verbose = FALSE)
sce_cd45 <- RunUMAP(sce_cd45, reduction = "harmony_sample", dims = dims_cd45, reduction.name = "umap_harmony_sample", umap.method = "uwot", verbose = FALSE)
sce_cd45 <- RunTSNE(sce_cd45, reduction = "harmony_sample", dims = dims_cd45, reduction.name = "tsne_harmony_sample", tsne.method = "Rtsne", check_duplicates = FALSE, verbose = FALSE)
sce_cd45 <- FindNeighbors(sce_cd45, reduction = "harmony_sample", dims = dims_cd45, k.param = k_cd45, graph.name = "RNA_nn_sample", verbose = FALSE)
sce_cd45 <- FindClusters(sce_cd45, graph.name = "RNA_nn_sample", resolution = res_cd45, cluster.name = "cluster_harmony_sample", verbose = FALSE)

p_cd45_sample <- DimPlot(sce_cd45, reduction = "tsne_harmony_sample", group.by = "sample", raster = FALSE, pt.size = 0.25) + ggtitle("CD45+ TIICs by sample")
p_cd45_cluster <- DimPlot(sce_cd45, reduction = "tsne_harmony_sample", group.by = "cluster_harmony_sample", label = TRUE, raster = FALSE, pt.size = 0.25) + ggtitle("CD45+ TIICs by cluster")
save_plot(p_cd45_sample | p_cd45_cluster, file.path(outdir, "plots", "CD45_branch", "FigS20B_CD45_tSNE_by_sample_and_cluster.pdf"), width = 12, height = 6)

# Major immune-cell annotation. Verify this map against marker dotplot if cluster numbers change.
celltype_major_map <- c(
  "0" = "T_NK", "1" = "Myeloid_TAM", "2" = "Myeloid_TAM", "3" = "Myeloid_TAM",
  "4" = "T_NK", "5" = "T_NK", "6" = "Neutrophil", "7" = "Myeloid_TAM",
  "8" = "Myeloid_TAM", "9" = "DC", "10" = "B", "11" = "Myeloid_TAM", "12" = NA
)
sce_cd45$celltype_major <- unname(celltype_major_map[as.character(sce_cd45$cluster_harmony_sample)])
write.csv(as.data.frame(table(sce_cd45$cluster_harmony_sample, sce_cd45$celltype_major, useNA = "ifany")), file.path(outdir, "CD45_cluster_to_major_celltype_table.csv"), row.names = FALSE)

p_major_tsne <- DimPlot(sce_cd45, reduction = "tsne_harmony_sample", group.by = "celltype_major", label = TRUE, raster = FALSE, pt.size = 0.25) + ggtitle("Major immune-cell lineages")
save_plot(p_major_tsne, file.path(outdir, "plots", "CD45_branch", "FigS20B_CD45_tSNE_by_major_lineage.pdf"), width = 8, height = 6)

marker_groups <- list(
  "CD45" = c("Ptprc"),
  "T" = c("Cd3d", "Cd3e", "Trbc2", "Il7r"),
  "NK" = c("Nkg7", "Klrd1", "Gzmb"),
  "B" = c("Ms4a1", "Cd79a", "Cd79b", "Cd19"),
  "TAM" = c("Csf1r", "Adgre1", "Apoe", "C1qa", "C1qb", "Trem2"),
  "Monocyte/Myeloid" = c("Lyz2", "Itgam", "Tyrobp", "Fcgr3", "Ccr2", "Plac8"),
  "Neutrophil" = c("S100a8", "S100a9", "Ly6g", "Retnlg"),
  "DC" = c("Xcr1", "Clec9a", "Itgax", "H2-Ab1", "Flt3", "Cd74"),
  "Mast" = c("Kit", "Ms4a2", "Cpa3", "Tpsab1")
)
marker_groups <- lapply(marker_groups, function(x) present_genes(x, sce_cd45))
p_major_dot <- DotPlot(sce_cd45, features = marker_groups, group.by = "celltype_major", cluster.idents = FALSE) + RotatedAxis() + ggtitle("Canonical markers of major immune-cell lineages")
save_plot(p_major_dot, file.path(outdir, "plots", "CD45_branch", "FigS20C_major_immune_marker_dotplot.pdf"), width = 14, height = 6)

Idents(sce_cd45) <- "cluster_harmony_sample"
markers_cd45 <- FindAllMarkers(sce_cd45, only.pos = TRUE, test.use = "wilcox", min.pct = 0.10, logfc.threshold = 0.25, return.thresh = 0.05, densify = TRUE, verbose = FALSE)
write.csv(markers_cd45, file.path(outdir, "CD45_global_FindAllMarkers.csv"), row.names = FALSE)

plot_composition(sce_cd45@meta.data, "celltype_major", file.path(outdir, "plots", "CD45_branch", "CD45_major_celltype"), subtype_order = c("T_NK", "Myeloid_TAM", "Neutrophil", "DC", "B"))
qs::qsave(sce_cd45, cd45_rds_final, nthreads = 8)

# =========================
# 5. T-cell subset analysis
# =========================
sce_tnk <- subset(sce_cd45, subset = celltype_major == "T_NK")
sce_tnk <- join_if_needed(sce_tnk, assay = "RNA")
sce_tnk <- ensure_percent_mt(sce_tnk)
sce_tnk <- NormalizeData(sce_tnk, verbose = FALSE)
sce_tnk <- FindVariableFeatures(sce_tnk, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sce_tnk <- ScaleData(sce_tnk, features = VariableFeatures(sce_tnk), vars.to.regress = "percent.mt", verbose = FALSE)
sce_tnk <- RunPCA(sce_tnk, npcs = 30, verbose = FALSE)
sce_tnk <- RunHarmony(sce_tnk, group.by.vars = "sample", reduction.use = "pca", dims.use = 1:15, reduction.save = "harmony_tnk", verbose = FALSE)
sce_tnk <- RunUMAP(sce_tnk, reduction = "harmony_tnk", dims = 1:15, reduction.name = "umap_tnk", verbose = FALSE)
sce_tnk <- FindNeighbors(sce_tnk, reduction = "harmony_tnk", dims = 1:15, k.param = 20, graph.name = c("tnk_nn", "tnk_snn"), verbose = FALSE)
sce_tnk <- FindClusters(sce_tnk, graph.name = "tnk_snn", resolution = 0.4, cluster.name = "cluster_tnk", verbose = FALSE)

# Keep T cells and exclude NK cluster according to marker validation.
tnk_lineage_map <- c("0" = "T", "1" = "T", "2" = "T", "3" = "T", "4" = "T", "5" = "T", "6" = "T", "7" = "NK", "8" = "T", "9" = "T")
sce_tnk$tnk_lineage <- unname(tnk_lineage_map[as.character(sce_tnk$cluster_tnk)])
sce_t <- subset(sce_tnk, subset = tnk_lineage == "T")
sce_t <- join_if_needed(sce_t, assay = "RNA")
sce_t <- ensure_percent_mt(sce_t)

sce_t <- NormalizeData(sce_t, verbose = FALSE)
sce_t <- FindVariableFeatures(sce_t, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sce_t <- ScaleData(sce_t, features = VariableFeatures(sce_t), vars.to.regress = "percent.mt", verbose = FALSE)
sce_t <- RunPCA(sce_t, npcs = npcs_t, verbose = FALSE)
pdf(file.path(outdir, "plots", "T_branch", "T_only_ElbowPlot.pdf"), width = 7, height = 5); print(ElbowPlot(sce_t, ndims = npcs_t)); dev.off()
sce_t <- RunHarmony(sce_t, group.by.vars = "sample", reduction.use = "pca", dims.use = dims_t, reduction.save = "harmony_t", verbose = FALSE)
sce_t <- RunUMAP(sce_t, reduction = "harmony_t", dims = dims_t, reduction.name = "umap_t", verbose = FALSE)
sce_t <- RunTSNE(sce_t, reduction = "harmony_t", dims = dims_t, reduction.name = "tsne_t", tsne.method = "Rtsne", check_duplicates = FALSE, verbose = FALSE)
sce_t <- FindNeighbors(sce_t, reduction = "harmony_t", dims = dims_t, k.param = k_t, graph.name = c("t_nn", "t_snn"), verbose = FALSE)
sce_t <- FindClusters(sce_t, graph.name = "t_snn", resolution = res_t, cluster.name = "cluster_t", verbose = FALSE)
sce_t$cluster_t <- factor(as.character(sce_t$cluster_t))

# T subtype annotation. Verify against marker dotplot if cluster IDs change.
t_subtype_map <- c(
  "0" = "Progenitor exhausted-like T",
  "1" = "Proliferating_Cytotoxic_like_T",
  "2" = "Resting_T",
  "3" = "Activated_Exhausted_like_T"
)
sce_t$t_subtype <- unname(t_subtype_map[as.character(sce_t$cluster_t)])
subtype_order_t <- c("Progenitor exhausted-like T", "Proliferating_Cytotoxic_like_T", "Resting_T", "Activated_Exhausted_like_T")
sce_t$t_subtype <- factor(sce_t$t_subtype, levels = subtype_order_t)
write.csv(as.data.frame(table(sce_t$cluster_t, sce_t$t_subtype, useNA = "ifany")), file.path(outdir, "T_cluster_to_subtype_table.csv"), row.names = FALSE)

p_t_tsne <- DimPlot(sce_t, reduction = "tsne_t", group.by = "t_subtype", label = FALSE, raster = FALSE, pt.size = 0.25) + ggtitle("T-cell subsets")
p_t_cluster <- DimPlot(sce_t, reduction = "tsne_t", group.by = "cluster_t", label = TRUE, raster = FALSE, pt.size = 0.25) + ggtitle("T-cell reclustering")
save_plot(p_t_cluster | p_t_tsne, file.path(outdir, "plots", "T_branch", "FigS20D_T_cell_reclustering_tSNE.pdf"), width = 12, height = 6)

t_marker_genes <- c("Tcf7", "Lef1", "Il7r", "Slamf6", "Ccr7", "Sell", "Klf2", "Ltb", "Mki67", "Top2a", "Birc5", "Stmn1", "Gzmb", "Prf1", "Cd8a", "Cd8b1", "Cxcr3", "Cd69", "Ifng", "Pdcd1", "Lag3", "Tigit", "Havcr2", "Ctla4", "Tox")
t_marker_genes <- present_genes(t_marker_genes, sce_t)
p_t_marker <- DotPlot(sce_t, features = t_marker_genes, group.by = "t_subtype", cluster.idents = FALSE) + RotatedAxis() + ggtitle("T-cell subset markers")
save_plot(p_t_marker, file.path(outdir, "plots", "T_branch", "FigS20E_T_cell_marker_dotplot.pdf"), width = 14, height = 6)

Idents(sce_t) <- "cluster_t"
markers_t <- FindAllMarkers(sce_t, only.pos = TRUE, test.use = "wilcox", min.pct = 0.10, logfc.threshold = 0.25, return.thresh = 0.05, densify = TRUE, verbose = FALSE)
write.csv(markers_t, file.path(outdir, "T_only_FindAllMarkers.csv"), row.names = FALSE)

# T functional module scores
module_sets_t <- list(
  Naive_Resting = c("Tcf7", "Lef1", "Il7r", "Sell", "Ccr7", "Klf2"),
  Progenitor_Exhaustion = c("Tcf7", "Slamf6", "Pdcd1", "Tox", "Cd27", "Cd28", "Id3", "Bcl6"),
  Terminal_Exhaustion = c("Pdcd1", "Lag3", "Havcr2", "Tigit", "Ctla4", "Tox", "Entpd1", "Prdm1", "Batf", "Cxcl13"),
  Cytotoxicity = c("Gzmb", "Gzma", "Prf1", "Nkg7", "Ifng", "Tnf", "Fasl"),
  Proliferation = c("Mki67", "Top2a", "Stmn1", "Birc5", "Cdk1", "Ccnb1", "Ccnb2", "Ube2c", "Hmgb2")
)
sce_t <- add_module_scores_named(sce_t, module_sets_t, min_genes = 2)
t_score_cols <- attr(sce_t, "module_score_cols")
plot_module_heatmap(sce_t, "t_subtype", t_score_cols, file.path(outdir, "plots", "T_branch", "FigS20E_T_functional_module_heatmap.pdf"), title = "Functional module scores across T-cell subsets")

t_subtype_cols <- c("Progenitor exhausted-like T" = "#7CAE00", "Proliferating_Cytotoxic_like_T" = "#00BFC4", "Resting_T" = "#C77CFF", "Activated_Exhausted_like_T" = "#F8766D")
plot_composition(sce_t@meta.data, "t_subtype", file.path(outdir, "plots", "T_branch", "FigS20F_T_subtype"), subtype_order = subtype_order_t, subtype_cols = t_subtype_cols)
qs::qsave(sce_t, file.path(outdir, "sce_T_cells_reclustered.qs"), nthreads = 8)

# =========================
# 6. TAM analysis
# =========================
sce_myeloid <- subset(sce_cd45, subset = celltype_major == "Myeloid_TAM")
sce_myeloid <- join_if_needed(sce_myeloid, assay = "RNA")
sce_myeloid <- ensure_percent_mt(sce_myeloid)
sce_myeloid <- NormalizeData(sce_myeloid, verbose = FALSE)
sce_myeloid <- FindVariableFeatures(sce_myeloid, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sce_myeloid <- ScaleData(sce_myeloid, features = VariableFeatures(sce_myeloid), vars.to.regress = "percent.mt", verbose = FALSE)
sce_myeloid <- RunPCA(sce_myeloid, npcs = npcs_myeloid, verbose = FALSE)
sce_myeloid <- RunHarmony(sce_myeloid, group.by.vars = "sample", reduction.use = "pca", dims.use = dims_myeloid, reduction.save = "harmony_myeloid", verbose = FALSE)
sce_myeloid <- RunUMAP(sce_myeloid, reduction = "harmony_myeloid", dims = dims_myeloid, reduction.name = "umap_myeloid", verbose = FALSE)
sce_myeloid <- FindNeighbors(sce_myeloid, reduction = "harmony_myeloid", dims = dims_myeloid, k.param = k_myeloid, graph.name = c("myeloid_nn", "myeloid_snn"), verbose = FALSE)
sce_myeloid <- FindClusters(sce_myeloid, graph.name = "myeloid_snn", resolution = res_myeloid, cluster.name = "cluster_myeloid", verbose = FALSE)

myeloid_lineage_map <- c("0" = "TAM", "1" = "TAM", "2" = "Monocyte_like", "3" = "Proliferating_TAM", "4" = "Non_myeloid_contam", "5" = "Neutrophil_like", "6" = "TAM", "7" = "TAM", "8" = "Fibroblast_like_contam", "9" = "Endothelial_like_contam")
sce_myeloid$myeloid_lineage <- unname(myeloid_lineage_map[as.character(sce_myeloid$cluster_myeloid)])
write.csv(as.data.frame(table(sce_myeloid$cluster_myeloid, sce_myeloid$myeloid_lineage, useNA = "ifany")), file.path(outdir, "Myeloid_cluster_to_lineage_table.csv"), row.names = FALSE)

# Include TAM clusters only, as in the submitted analysis.
sce_tam_raw <- subset(sce_myeloid, subset = cluster_myeloid %in% c("0", "1", "3", "6", "7"))
sce_tam_raw <- ensure_percent_mt(sce_tam_raw)

# TAM-specific QC
p_tam_qc <- VlnPlot(sce_tam_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "cluster_myeloid", pt.size = 0.02, ncol = 3)
save_plot(p_tam_qc, file.path(outdir, "plots", "TAM_branch", "TAM_raw_QC_vlnplot.pdf"), width = 14, height = 5)
sce_tam <- subset(sce_tam_raw, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 10)
sce_tam <- join_if_needed(sce_tam, assay = "RNA")

sce_tam <- NormalizeData(sce_tam, verbose = FALSE)
sce_tam <- FindVariableFeatures(sce_tam, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sce_tam <- ScaleData(sce_tam, features = VariableFeatures(sce_tam), vars.to.regress = "percent.mt", verbose = FALSE)
sce_tam <- RunPCA(sce_tam, npcs = npcs_tam, verbose = FALSE)
sce_tam <- RunHarmony(sce_tam, group.by.vars = "sample", reduction.use = "pca", dims.use = dims_tam, reduction.save = "harmony_tam", verbose = FALSE)
sce_tam <- RunUMAP(sce_tam, reduction = "harmony_tam", dims = dims_tam, reduction.name = "umap_tam", verbose = FALSE)
sce_tam <- RunTSNE(sce_tam, reduction = "harmony_tam", dims = dims_tam, reduction.name = "tsne_tam", tsne.method = "Rtsne", check_duplicates = FALSE, verbose = FALSE)
sce_tam <- FindNeighbors(sce_tam, reduction = "harmony_tam", dims = dims_tam, k.param = k_tam, graph.name = c("tam_nn", "tam_snn"), verbose = FALSE)
sce_tam <- FindClusters(sce_tam, graph.name = "tam_snn", resolution = res_tam, cluster.name = "cluster_tam", verbose = FALSE)
sce_tam$cluster_tam <- factor(as.character(sce_tam$cluster_tam))

# TAM subtype annotation. Verify against marker dotplot if cluster IDs change.
tam_subtype_map <- c(
  "0" = "C0_Malt1+ TAM",
  "1" = "C1_Mki67+ TAM",
  "2" = "C2_F13a1+ TAM",
  "3" = "C3_Arg1+ TAM",
  "4" = "C4_Marco+ TAM"
)
sce_tam$tam_subtype <- unname(tam_subtype_map[as.character(sce_tam$cluster_tam)])
tam_order <- c("C0_Malt1+ TAM", "C1_Mki67+ TAM", "C2_F13a1+ TAM", "C3_Arg1+ TAM", "C4_Marco+ TAM")
sce_tam$tam_subtype <- factor(sce_tam$tam_subtype, levels = tam_order)
write.csv(as.data.frame(table(sce_tam$cluster_tam, sce_tam$tam_subtype, useNA = "ifany")), file.path(outdir, "TAM_cluster_to_subtype_table.csv"), row.names = FALSE)

p_tam_tsne <- DimPlot(sce_tam, reduction = "tsne_tam", group.by = "tam_subtype", label = FALSE, raster = FALSE, pt.size = 0.25) + ggtitle("TAM subsets")
p_tam_cluster <- DimPlot(sce_tam, reduction = "tsne_tam", group.by = "cluster_tam", label = TRUE, raster = FALSE, pt.size = 0.25) + ggtitle("TAM reclustering")
save_plot(p_tam_cluster | p_tam_tsne, file.path(outdir, "plots", "TAM_branch", "FigS20G_TAM_reclustering_tSNE.pdf"), width = 12, height = 6)

tam_marker_features <- c(
  "Malt1", "Il1b", "Ccl3", "Ccl4", "Tnf", "Nfkbia",
  "Mki67", "Top2a", "Birc5", "Stmn1", "Cdk1", "Ccnb2",
  "F13a1", "Folr2", "Gas6", "Mrc1", "Lyve1", "Cd163",
  "Arg1", "Chil3", "Fn1", "Vegfa", "Spp1", "Mmp14",
  "Marco", "Clec4f", "Vsig4", "Timd4", "Cd5l", "C1qa"
)
tam_marker_features <- present_genes(tam_marker_features, sce_tam)
p_tam_marker <- DotPlot(sce_tam, features = tam_marker_features, group.by = "tam_subtype", cluster.idents = FALSE) + RotatedAxis() + ggtitle("TAM subset markers")
save_plot(p_tam_marker, file.path(outdir, "plots", "TAM_branch", "FigS20H_TAM_marker_dotplot.pdf"), width = 14, height = 6)

Idents(sce_tam) <- "cluster_tam"
markers_tam <- FindAllMarkers(sce_tam, only.pos = TRUE, test.use = "wilcox", min.pct = 0.10, logfc.threshold = 0.25, return.thresh = 0.05, densify = TRUE, verbose = FALSE)
write.csv(markers_tam, file.path(outdir, "TAM_FindAllMarkers.csv"), row.names = FALSE)

# TAM hallmark-like functional scores. Uses msigdbr if available; otherwise uses manual marker modules.
if (requireNamespace("msigdbr", quietly = TRUE)) {
  library(msigdbr)
  tam_hallmark_pathway_map <- c(
    Inflammatory_response = "HALLMARK_INFLAMMATORY_RESPONSE",
    TNFA_NFKB_signaling = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    IFN_alpha_response = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    IFN_gamma_response = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    IL6_JAK_STAT3_signaling = "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    Angiogenesis = "HALLMARK_ANGIOGENESIS",
    EMT_matrix_remodeling = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    ROS_pathway = "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    Glycolysis = "HALLMARK_GLYCOLYSIS",
    Oxidative_phosphorylation = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    Fatty_acid_metabolism = "HALLMARK_FATTY_ACID_METABOLISM",
    G2M_checkpoint = "HALLMARK_G2M_CHECKPOINT",
    E2F_targets = "HALLMARK_E2F_TARGETS"
  )
  msig_h <- tryCatch(msigdbr::msigdbr(db_species = "MM", species = "Mus musculus", collection = "H"), error = function(e) msigdbr::msigdbr(species = "Mus musculus", category = "H"))
  module_sets_tam <- lapply(tam_hallmark_pathway_map, function(gs) unique(msig_h$gene_symbol[msig_h$gs_name == gs]))
} else {
  module_sets_tam <- list(
    Inflammatory_response = c("Il1b", "Tnf", "Ccl3", "Ccl4", "Nfkbia", "Cxcl2"),
    IFN_response = c("Isg15", "Ifit1", "Ifit2", "Ifit3", "Irf7", "Stat1", "Cxcl10"),
    Proliferation = c("Mki67", "Top2a", "Stmn1", "Birc5", "Cdk1", "Ccnb1", "Ccnb2"),
    Angio_remodeling = c("Arg1", "Fn1", "Vegfa", "Spp1", "Mmp14", "Mmp9"),
    Resident_phagocytic = c("Marco", "Clec4f", "Vsig4", "Timd4", "Cd5l", "C1qa")
  )
}
sce_tam <- add_module_scores_named(sce_tam, module_sets_tam, min_genes = 5)
tam_score_cols <- attr(sce_tam, "module_score_cols")
plot_module_heatmap(sce_tam, "tam_subtype", tam_score_cols, file.path(outdir, "plots", "TAM_branch", "FigS20I_TAM_functional_module_heatmap.pdf"), title = "Functional pathway scores across TAM subsets")

tam_cols <- c("C0_Malt1+ TAM" = "#F8766D", "C1_Mki67+ TAM" = "#A3A500", "C2_F13a1+ TAM" = "#00BF7D", "C3_Arg1+ TAM" = "#00B0F6", "C4_Marco+ TAM" = "#E76BF3")
plot_composition(sce_tam@meta.data, "tam_subtype", file.path(outdir, "plots", "TAM_branch", "FigS20J_K_TAM_subtype"), subtype_order = tam_order, subtype_cols = tam_cols)
qs::qsave(sce_tam, file.path(outdir, "sce_TAM_reclustered.qs"), nthreads = 8)

# =========================
# 7. DC subset analysis
# =========================
sce_dc <- subset(sce_cd45, subset = celltype_major == "DC")
sce_dc <- join_if_needed(sce_dc, assay = "RNA")
sce_dc <- ensure_percent_mt(sce_dc)
sce_dc <- NormalizeData(sce_dc, verbose = FALSE)
sce_dc <- FindVariableFeatures(sce_dc, selection.method = "vst", nfeatures = 1500, verbose = FALSE)
sce_dc <- ScaleData(sce_dc, features = VariableFeatures(sce_dc), vars.to.regress = "percent.mt", verbose = FALSE)
sce_dc <- RunPCA(sce_dc, npcs = npcs_dc, verbose = FALSE)
sce_dc <- RunHarmony(sce_dc, group.by.vars = "sample", reduction.use = "pca", dims.use = dims_dc, reduction.save = "harmony_dc", verbose = FALSE)
sce_dc <- RunUMAP(sce_dc, reduction = "harmony_dc", dims = dims_dc, reduction.name = "umap_dc", verbose = FALSE)
sce_dc <- RunTSNE(sce_dc, reduction = "harmony_dc", dims = dims_dc, reduction.name = "tsne_dc", tsne.method = "Rtsne", check_duplicates = FALSE, verbose = FALSE)
sce_dc <- FindNeighbors(sce_dc, reduction = "harmony_dc", dims = dims_dc, k.param = k_dc, graph.name = c("dc_nn", "dc_snn"), verbose = FALSE)
sce_dc <- FindClusters(sce_dc, graph.name = "dc_snn", resolution = res_dc, cluster.name = "cluster_dc", verbose = FALSE)
sce_dc$cluster_dc <- factor(as.character(sce_dc$cluster_dc))

dc_subtype_map <- c(
  "0" = "Inflammatory_activated_DC",
  "1" = "cDC2_MHCII_DC",
  "2" = "cDC1_cross_presenting_DC",
  "3" = "IFN_response_DC"
)
sce_dc$dc_subtype <- unname(dc_subtype_map[as.character(sce_dc$cluster_dc)])
dc_order <- c("Inflammatory_activated_DC", "cDC1_cross_presenting_DC", "cDC2_MHCII_DC", "IFN_response_DC")
sce_dc$dc_subtype <- factor(sce_dc$dc_subtype, levels = dc_order)
write.csv(as.data.frame(table(sce_dc$cluster_dc, sce_dc$dc_subtype, useNA = "ifany")), file.path(outdir, "DC_cluster_to_subtype_table.csv"), row.names = FALSE)

p_dc_tsne <- DimPlot(sce_dc, reduction = "tsne_dc", group.by = "dc_subtype", label = FALSE, raster = FALSE, pt.size = 0.5) + ggtitle("DC subsets")
p_dc_cluster <- DimPlot(sce_dc, reduction = "tsne_dc", group.by = "cluster_dc", label = TRUE, raster = FALSE, pt.size = 0.5) + ggtitle("DC reclustering")
save_plot(p_dc_cluster | p_dc_tsne, file.path(outdir, "plots", "DC_branch", "FigS20L_DC_reclustering_tSNE.pdf"), width = 12, height = 6)

dc_marker_features <- c(
  "Ccr7", "Cd40", "Cd80", "Cd83", "Cd86", "Tnf", "Il1b", "Ccl17", "Ccl22", "Ccl19",
  "Sirpa", "Irf4", "Itgam", "Cd209a", "H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74",
  "Clec9a", "Xcr1", "Batf3", "Irf8", "Cadm1", "Wdfy4", "Tlr3",
  "Isg15", "Ifit1", "Ifit2", "Ifit3", "Irf7", "Stat1", "Cxcl10", "Oas1a", "Mx1"
)
dc_marker_features <- present_genes(dc_marker_features, sce_dc)
p_dc_marker <- DotPlot(sce_dc, features = dc_marker_features, group.by = "dc_subtype", cluster.idents = FALSE) + RotatedAxis() + ggtitle("DC subset markers")
save_plot(p_dc_marker, file.path(outdir, "plots", "DC_branch", "FigS20M_DC_marker_dotplot.pdf"), width = 14, height = 6)

Idents(sce_dc) <- "cluster_dc"
markers_dc <- FindAllMarkers(sce_dc, only.pos = TRUE, test.use = "wilcox", min.pct = 0.10, logfc.threshold = 0.25, return.thresh = 0.05, verbose = FALSE)
write.csv(markers_dc, file.path(outdir, "DC_FindAllMarkers.csv"), row.names = FALSE)

dc_cols <- c("Inflammatory_activated_DC" = "#F8766D", "cDC1_cross_presenting_DC" = "#00BFC4", "cDC2_MHCII_DC" = "#7CAE00", "IFN_response_DC" = "#C77CFF")
plot_composition(sce_dc@meta.data, "dc_subtype", file.path(outdir, "plots", "DC_branch", "FigS20N_DC_subtype"), subtype_order = dc_order, subtype_cols = dc_cols)
qs::qsave(sce_dc, file.path(outdir, "sce_DC_reclustered.qs"), nthreads = 8)

# =========================
# 8. Session information
# =========================
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("Figure S20 CD45+ scRNA-seq analysis finished. Results saved to: ", outdir)
