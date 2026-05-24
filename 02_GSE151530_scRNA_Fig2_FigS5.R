############################################################
# Figure S5 scRNA-seq analysis for GSE151530 HCC samples
# Revised version for manuscript upload
# Key changes from the original exploratory script:
#   1) Group names are prognosis-based: Favorable prognosis vs Poor prognosis
#   2) Added nCount_RNA filtering
#   3) Added per-sample doublet removal using DoubletFinder
#   4) Removed analyses not shown in Figure S5, e.g. CopyKAT/inferCNV
#   5) Added DEG-derived marker plot for Fig. S5B and T-cell DEG heatmap for Fig. S5D
#   6) Kept the original major Seurat parameters unless there was an obvious issue
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(2024)

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(DoubletFinder)
  library(dplyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(Matrix)
  library(patchwork)
  library(reshape2)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(TCellSI)
  library(pheatmap)
})

# ===================== 0. User-defined paths and parameters =====================
workdir <- "path/to/GSE151530_analysis"
setwd(workdir)

outdir <- "Figure_S5_scRNA_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "QC"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "Global"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "T_cells"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "HCC_cells"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "GSEA"), showWarnings = FALSE, recursive = TRUE)

favorable_dir <- "Favorable_prognosis"
poor_dir <- "Poor_prognosis"

# Fallback mapping if sample_group.csv is not provided:
favorable_dir_group <- "Favorable prognosis"
poor_dir_group <- "Poor prognosis"

sample_group_file <- "sample_group.csv"       # optional but strongly recommended
annot_file <- "GSE151530_Info.csv"            # columns should include cell/barcode and type/cell_type
pyro_gmt_file <- "GOBP_PYROPTOSIS.v2024.1.Hs.gmt"
kegg_gmt_file <- "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"

# QC parameters: original parameters retained, with nCount_RNA filter added.
min_features <- 200
max_features <- 5000
max_percent_mt <- 20
min_counts <- 500
max_counts <- Inf
# If your nCount_RNA distribution shows a clear high-count tail/doublet-enriched region,
# replace Inf with a dataset-specific cutoff, e.g. 50000.

# Seurat parameters retained from original script unless otherwise noted.
n_variable_features <- 2000
n_pcs_global <- 30
global_resolution <- 0.2
n_pcs_tcell <- 30
tcell_resolution <- 0.05

# Colors for prognosis groups.
group_cols <- c("Favorable prognosis" = "#4575b4", "Poor prognosis" = "#fc8d59")

# Marker/signature genes used for Figure S5-related validation.
epithelial_markers <- c(
  "KRT14", "KRT17", "KRT6A", "KRT5", "KRT19", "KRT8", "KRT16",
  "KRT18", "KRT6B", "KRT15", "KRT6C", "KRTCAP3", "SFN", "EPCAM"
)

canonical_markers <- c(
  "CD3D", "CD3E", "CD2", "CD8A", "CD8B", "NKG7", "GZMB", "PRF1", "IFNG",
  "CD79A", "MS4A1", "CD79B", "CD68", "CD14", "LYZ", "PECAM1", "ENG",
  "DCN", "COL1A1", "COL1A2", "COL3A1", "EPCAM", "KRT8", "KRT18", "KRT19", "AFP", "ALB"
)

# ===================== 1. Utility functions =====================
read_count_csv <- function(path) {
  x <- read.csv(path, row.names = 1, check.names = FALSE)
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  return(x)
}

extract_gsm <- function(fname) {
  str_split_fixed(basename(fname), "_", 2)[, 1]
}

make_seurat_object <- function(fname, input_dir, fallback_group, sample_group_tbl = NULL) {
  gsm <- extract_gsm(fname)
  cnts <- read_count_csv(file.path(input_dir, fname))

  obj <- CreateSeuratObject(
    counts = cnts,
    project = gsm,
    min.cells = 3,
    min.features = 0
  )
  obj$GSM <- gsm

  if (!is.null(sample_group_tbl)) {
    grp <- sample_group_tbl$prognosis_group[match(gsm, sample_group_tbl$GSM)]
    if (is.na(grp)) stop(paste0("No prognosis group found for sample: ", gsm))
    obj$prognosis_group <- grp
  } else {
    obj$prognosis_group <- fallback_group
  }

  obj$orig_folder <- input_dir
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.Ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  return(obj)
}

run_doublet_removal <- function(seurat_obj, expected_doublet_rate = 0.075) {
  # DoubletFinder is run per sample to avoid cross-sample artifacts.
  
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  
  n_pcs_df <- min(30, ncol(seurat_obj) - 1)
  seurat_obj <- RunPCA(seurat_obj, npcs = n_pcs_df, verbose = FALSE)
  
  sweep_res <- paramSweep(seurat_obj, PCs = 1:n_pcs_df, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  nExp <- round(expected_doublet_rate * ncol(seurat_obj))
  
  seurat_obj <- doubletFinder(
    seurat_obj,
    PCs = 1:n_pcs_df,
    pN = 0.25,
    pK = pK,
    nExp = nExp,
    reuse.pANN = FALSE,
    sct = FALSE
  )
  
  df_col <- grep("DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
  pANN_col <- grep("pANN", colnames(seurat_obj@meta.data), value = TRUE)
  
  seurat_obj$doublet_class <- seurat_obj@meta.data[[df_col[length(df_col)]]]
  seurat_obj$doublet_score <- seurat_obj@meta.data[[pANN_col[length(pANN_col)]]]
  
  seurat_obj <- subset(seurat_obj, subset = doublet_class == "Singlet")
  return(seurat_obj)
}

save_pdf <- function(filename, plot, width = 7, height = 6) {
  ggsave(filename = filename, plot = plot, width = width, height = height, useDingbats = FALSE)
}

get_present_genes <- function(genes, obj) {
  intersect(genes, rownames(obj))
}

safe_find_ppar_term <- function(gsea_obj) {
  res <- as.data.frame(gsea_obj)
  idx <- grep("PPAR", res$ID, ignore.case = TRUE)
  if (length(idx) == 0) idx <- grep("PPAR", res$Description, ignore.case = TRUE)
  if (length(idx) == 0) return(NA_character_)
  return(res$ID[idx[1]])
}

# ===================== 2. Read sample group information =====================
sample_group_tbl <- NULL
if (file.exists(sample_group_file)) {
  sample_group_tbl <- fread(sample_group_file) |> as.data.frame()
  colnames(sample_group_tbl) <- tolower(colnames(sample_group_tbl))
  if (!all(c("gsm", "prognosis_group") %in% colnames(sample_group_tbl))) {
    stop("sample_group.csv must contain columns: GSM, prognosis_group")
  }
  colnames(sample_group_tbl)[colnames(sample_group_tbl) == "gsm"] <- "GSM"
  valid_groups <- c("Favorable prognosis", "Poor prognosis")
  if (!all(sample_group_tbl$prognosis_group %in% valid_groups)) {
    stop("prognosis_group must be either 'Favorable prognosis' or 'Poor prognosis'.")
  }
} else {
  message("sample_group.csv not found. Using folder-based fallback mapping: Favorable_prognosis = Favorable prognosis; Poor_prognosis = Poor prognosis.")
}

# ===================== 3. Read count matrices and perform QC/doublet filtering =====================
filenames_favorable <- list.files(favorable_dir, full.names = FALSE)
filenames_poor <- list.files(poor_dir, full.names = FALSE)

obj_list_favorable <- lapply(
  filenames_favorable,
  make_seurat_object,
  input_dir = favorable_dir,
  fallback_group = favorable_dir_group,
  sample_group_tbl = sample_group_tbl
)

obj_list_poor <- lapply(
  filenames_poor,
  make_seurat_object,
  input_dir = poor_dir,
  fallback_group = poor_dir_group,
  sample_group_tbl = sample_group_tbl
)

obj_list <- c(obj_list_favorable, obj_list_poor)
names(obj_list) <- c(filenames_favorable, filenames_poor)


raw_meta <- do.call(rbind, lapply(obj_list, function(x) x@meta.data))
write.csv(raw_meta, file.path(outdir, "QC", "metadata_before_filtering.csv"), row.names = TRUE)

# Standard QC with added nCount_RNA thresholds.
obj_list_qc <- lapply(obj_list, function(x) {
  subset(
    x,
    subset = nFeature_RNA > min_features &
      nFeature_RNA < max_features &
      percent.mt < max_percent_mt &
      nCount_RNA > min_counts &
      nCount_RNA < max_counts
  )
})

qc_meta <- do.call(rbind, lapply(obj_list_qc, function(x) x@meta.data))
write.csv(qc_meta, file.path(outdir, "QC", "metadata_after_qc_before_doublet_removal.csv"), row.names = TRUE)

# Per-sample doublet removal.
obj_list_singlet <- lapply(obj_list_qc, run_doublet_removal)
singlet_meta <- do.call(rbind, lapply(obj_list_singlet, function(x) x@meta.data))
write.csv(singlet_meta, file.path(outdir, "QC", "metadata_after_doublet_removal.csv"), row.names = TRUE)

# Cell count summary.
cell_count_summary <- data.frame(
  GSM = unique(raw_meta$GSM)
) |>
  mutate(
    raw_count = as.numeric(table(raw_meta$GSM)[GSM]),
    after_qc_count = as.numeric(table(qc_meta$GSM)[GSM]),
    singlet_count = as.numeric(table(singlet_meta$GSM)[GSM])
  )
write.csv(cell_count_summary, file.path(outdir, "QC", "cell_count_by_sample.csv"), row.names = FALSE)

group_count_summary <- data.frame(
  group = names(table(raw_meta$prognosis_group)),
  raw_count = as.numeric(table(raw_meta$prognosis_group)),
  after_qc_count = as.numeric(table(qc_meta$prognosis_group)[names(table(raw_meta$prognosis_group))]),
  singlet_count = as.numeric(table(singlet_meta$prognosis_group)[names(table(raw_meta$prognosis_group))])
)
group_count_summary[is.na(group_count_summary)] <- 0
write.csv(group_count_summary, file.path(outdir, "QC", "cell_count_by_group.csv"), row.names = FALSE)

# Merge singlet cells.
sce <- merge(obj_list_singlet[[1]], y = obj_list_singlet[-1])
sce$prognosis_group <- factor(sce$prognosis_group, levels = c("Favorable prognosis", "Poor prognosis"))
sce$type <- sce$prognosis_group  # Compatibility with older plotting code.

# QC plots after filtering/doublet removal.
p_qc_vln <- VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    group.by = "prognosis_group", ncol = 3, pt.size = 0)
save_pdf(file.path(outdir, "QC", "QC_violin_after_filter_doublet.pdf"), p_qc_vln, width = 10, height = 4)

p_qc_scatter <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                               group.by = "prognosis_group")
save_pdf(file.path(outdir, "QC", "QC_nCount_vs_nFeature_after_filter_doublet.pdf"), p_qc_scatter, width = 6, height = 5)

# ===================== 4. Global Seurat workflow =====================
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = n_variable_features)
sce <- ScaleData(sce)
sce <- RunPCA(sce, features = VariableFeatures(sce))

p_elbow <- ElbowPlot(sce, ndims = 50)
save_pdf(file.path(outdir, "Global", "elbowplot_global.pdf"), p_elbow, width = 7, height = 5)

sce <- RunUMAP(sce, dims = 1:n_pcs_global)
sce <- RunTSNE(sce, dims = 1:n_pcs_global)
sce <- FindNeighbors(sce, dims = 1:n_pcs_global)
sce <- FindClusters(sce, resolution = global_resolution)

p_tsne_group <- DimPlot(sce, reduction = "tsne", group.by = "prognosis_group",
                        cols = group_cols, label = FALSE) +
  ggtitle("tSNE by prognosis group") + theme(plot.title = element_text(hjust = 0.5))
save_pdf(file.path(outdir, "Global", "FigS5_tsne_by_prognosis_group.pdf"), p_tsne_group, width = 7, height = 6)

p_umap_group <- DimPlot(sce, reduction = "umap", group.by = "prognosis_group",
                        cols = group_cols, label = FALSE) +
  ggtitle("UMAP by prognosis group") + theme(plot.title = element_text(hjust = 0.5))
save_pdf(file.path(outdir, "Global", "umap_by_prognosis_group.pdf"), p_umap_group, width = 7, height = 6)

p_tsne_cluster <- DimPlot(sce, reduction = "tsne", group.by = "seurat_clusters",
                          label = TRUE, repel = TRUE) +
  ggtitle("tSNE by cluster") + theme(plot.title = element_text(hjust = 0.5))
save_pdf(file.path(outdir, "Global", "tsne_by_cluster.pdf"), p_tsne_cluster, width = 7, height = 6)

# ===================== 5. Cell type annotation from GSE151530 metadata =====================
if (!file.exists(annot_file)) {
  stop("GSE151530_Info.csv not found. This revised Figure S5 script uses the reference metadata for cell annotation.")
}

annot_tbl <- fread(annot_file) |> as.data.frame()
colnames(annot_tbl) <- tolower(colnames(annot_tbl))
if (!"cell" %in% colnames(annot_tbl)) stop("Annotation file must contain a cell/barcode column named 'cell'.")
if (!"type" %in% colnames(annot_tbl) && !"cell_type" %in% colnames(annot_tbl)) {
  stop("Annotation file must contain a type or cell_type column.")
}
type_col <- ifelse("type" %in% colnames(annot_tbl), "type", "cell_type")

# Standardize barcodes for matching.
barcodes_sce <- gsub("\\.[0-9]+$", "", colnames(sce))
barcodes_sce <- gsub("-[0-9]+$", "", barcodes_sce)
barcodes_annot <- gsub("\\.[0-9]+$", "", as.character(annot_tbl$cell))
barcodes_annot <- gsub("-[0-9]+$", "", barcodes_annot)
match_idx <- match(barcodes_sce, barcodes_annot)
sce$cell_type <- annot_tbl[[type_col]][match_idx]

write.csv(table(is.na(sce$cell_type)), file.path(outdir, "Global", "cell_type_annotation_NA_table.csv"))
# Remove cells without reference annotation or labeled unclassified.
sce <- subset(sce, subset = !is.na(cell_type) & cell_type != "unclassified")

# Harmonize malignant cell naming for downstream analysis.
sce$cell_type <- as.character(sce$cell_type)
sce$cell_type[sce$cell_type %in% c("HCC cells", "HCC", "Tumor cells", "Malignant epithelial cells")] <- "Malignant cells"
sce$cell_type <- factor(sce$cell_type)

p_tsne_celltype <- DimPlot(sce, reduction = "tsne", group.by = "cell_type", label = TRUE, repel = TRUE) +
  ggtitle("tSNE by cell type") + theme(plot.title = element_text(hjust = 0.5))
save_pdf(file.path(outdir, "Global", "FigS5_tsne_by_cell_type.pdf"), p_tsne_celltype, width = 8, height = 7)

p_umap_celltype <- DimPlot(sce, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
  ggtitle("UMAP by cell type") + theme(plot.title = element_text(hjust = 0.5))
save_pdf(file.path(outdir, "Global", "umap_by_cell_type.pdf"), p_umap_celltype, width = 8, height = 7)

# Marker validation plot.
markers_use <- get_present_genes(canonical_markers, sce)
p_marker_dot <- DotPlot(sce, features = markers_use, group.by = "cell_type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Canonical marker expression by cell type")
save_pdf(file.path(outdir, "Global", "canonical_marker_dotplot_by_celltype.pdf"), p_marker_dot, width = 12, height = 6)

# Fig. S5B: marker genes identified by differential expression analysis for cell-type annotation.
# This is separated from the canonical marker dotplot above, because the figure legend describes DEG-derived markers.
DefaultAssay(sce) <- "RNA"
Idents(sce) <- sce$cell_type
celltype_markers <- FindAllMarkers(
  sce,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(celltype_markers, file.path(outdir, "Global", "FigS5B_celltype_DEG_markers.csv"), row.names = FALSE)

celltype_top_markers <- celltype_markers |>
  dplyr::filter(p_val_adj < 0.05) |>
  dplyr::group_by(cluster) |>
  dplyr::arrange(dplyr::desc(avg_log2FC), .by_group = TRUE) |>
  dplyr::slice_head(n = 5) |>
  dplyr::ungroup()

if (nrow(celltype_top_markers) > 0) {
  p_s5b_deg_marker <- DotPlot(sce, features = unique(celltype_top_markers$gene), group.by = "cell_type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("DEG-derived marker genes for cell-type annotation")
  save_pdf(file.path(outdir, "Global", "FigS5B_DEG_marker_dotplot_by_celltype.pdf"),
           p_s5b_deg_marker, width = 12, height = 6)
} else {
  warning("No significant cell-type DEG markers found for Fig. S5B under the current thresholds.")
}

# Epithelial score for malignant cell validation.
epi_use <- get_present_genes(epithelial_markers, sce)
sce <- AddModuleScore(sce, features = list(epi_use), name = "EpithelialScore")
p_epi_tsne <- FeaturePlot(sce, features = "EpithelialScore1", reduction = "tsne",
                          cols = c("grey90", "orange", "red"), pt.size = 0.5) +
  ggtitle("Epithelial score")
save_pdf(file.path(outdir, "Global", "FigS5_epithelial_score_tsne.pdf"), p_epi_tsne, width = 7, height = 6)

p_epi_vln <- VlnPlot(sce, features = "EpithelialScore1", group.by = "cell_type", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Epithelial score by cell type")
save_pdf(file.path(outdir, "Global", "epithelial_score_by_celltype.pdf"), p_epi_vln, width = 9, height = 5)

# Cell type composition by prognosis group.
celltype_counts <- table(sce$cell_type, sce$prognosis_group)
df_celltype <- as.data.frame(celltype_counts)
colnames(df_celltype) <- c("CellType", "Group", "Count")
df_celltype <- df_celltype |>
  group_by(Group) |>
  mutate(Percent = Count / sum(Count) * 100) |>
  ungroup()
write.csv(df_celltype, file.path(outdir, "Global", "celltype_composition_by_prognosis.csv"), row.names = FALSE)

p_celltype_stack <- ggplot(df_celltype, aes(x = Group, y = Percent, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  labs(title = "Cell type composition by prognosis group", x = NULL, y = "Cell percentage (%)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_pdf(file.path(outdir, "Global", "FigS5_celltype_composition_by_prognosis.pdf"), p_celltype_stack, width = 8, height = 6)

# ===================== 6. HCC/malignant cell analyses: pyroptosis, PPAR family, GSEA =====================
if (!"Malignant cells" %in% unique(sce$cell_type)) {
  stop("No 'Malignant cells' found after annotation harmonization. Please check cell_type labels in GSE151530_Info.csv.")
}

sce_hcc <- subset(sce, subset = cell_type == "Malignant cells")

# Pyroptosis score in HCC cells.
if (file.exists(pyro_gmt_file)) {
  pyro_gmt <- read.gmt(pyro_gmt_file)
  pyro_genes <- unique(pyro_gmt$gene)
  pyro_use <- intersect(pyro_genes, rownames(sce_hcc))
  if (length(pyro_use) >= 5) {
    sce_hcc <- AddModuleScore(sce_hcc, features = list(pyro_use), name = "PyroptosisScore")
    pyro_by_sample <- sce_hcc@meta.data |>
      group_by(GSM, prognosis_group) |>
      summarise(PyroptosisScore = mean(PyroptosisScore1, na.rm = TRUE), n_cell = n(), .groups = "drop")
    write.csv(pyro_by_sample, file.path(outdir, "HCC_cells", "pyroptosis_score_by_sample_HCC_cells.csv"), row.names = FALSE)

    p_pyro <- ggboxplot(pyro_by_sample, x = "prognosis_group", y = "PyroptosisScore",
                        fill = "prognosis_group", palette = group_cols, add = "jitter") +
      stat_compare_means(method = "wilcox.test") +
      labs(title = "Pyroptosis score in HCC cells", x = NULL, y = "Mean pyroptosis score per sample") +
      theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
    save_pdf(file.path(outdir, "HCC_cells", "FigS5_pyroptosis_score_HCC_cells_by_prognosis.pdf"), p_pyro, width = 5, height = 5)
  } else {
    warning("Too few pyroptosis genes found in the HCC cell expression matrix. Skipped pyroptosis scoring.")
  }
} else {
  warning("Pyroptosis GMT file not found. Skipped pyroptosis scoring.")
}

# PPARA/PPARG/PPARD expression across cell types and in HCC cells by prognosis.
ppar_genes <- intersect(c("PPARA", "PPARG", "PPARD"), rownames(sce))
if (length(ppar_genes) > 0) {
  p_ppar_dot <- DotPlot(sce, features = ppar_genes, group.by = "cell_type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("PPAR family expression across cell types")
  save_pdf(file.path(outdir, "Global", "FigS5_PPAR_family_dotplot_by_celltype.pdf"), p_ppar_dot, width = 8, height = 5)
}

for (gene in ppar_genes) {
  sce_hcc@meta.data[[paste0(gene, "_expression")]] <- FetchData(sce_hcc, vars = gene)[, 1]
  gene_by_sample <- sce_hcc@meta.data |>
    group_by(GSM, prognosis_group) |>
    summarise(expr_mean = mean(.data[[paste0(gene, "_expression")]], na.rm = TRUE), n_cell = n(), .groups = "drop")
  write.csv(gene_by_sample, file.path(outdir, "HCC_cells", paste0(gene, "_expression_by_sample_HCC_cells.csv")), row.names = FALSE)

  p_gene <- ggboxplot(gene_by_sample, x = "prognosis_group", y = "expr_mean",
                      fill = "prognosis_group", palette = group_cols, add = "jitter") +
    stat_compare_means(method = "wilcox.test") +
    labs(title = paste0(gene, " expression in HCC cells"), x = NULL, y = "Mean expression per sample") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
  save_pdf(file.path(outdir, "HCC_cells", paste0("FigS5_", gene, "_expression_HCC_cells_by_prognosis.pdf")), p_gene, width = 5, height = 5)
}

# Differential expression and GSEA in HCC cells: favorable versus poor prognosis.
# Pseudobulk DEG is used to reduce pseudoreplication.

counts_hcc <- GetAssayData(sce_hcc, assay = "RNA", slot = "counts")
meta_hcc <- sce_hcc@meta.data

hcc_cell_count_by_sample <- meta_hcc |>
  dplyr::count(GSM, prognosis_group, name = "n_hcc_cells")

write.csv(
  hcc_cell_count_by_sample,
  file.path(outdir, "HCC_cells", "HCC_cell_count_by_sample.csv"),
  row.names = FALSE
)

sample_ids <- unique(meta_hcc$GSM)

pb_counts <- sapply(sample_ids, function(sid) {
  cells_use <- rownames(meta_hcc)[meta_hcc$GSM == sid]
  Matrix::rowSums(counts_hcc[, cells_use, drop = FALSE])
})

pb_counts <- as.matrix(pb_counts)

pb_meta <- meta_hcc |>
  dplyr::distinct(GSM, prognosis_group) |>
  dplyr::arrange(match(GSM, colnames(pb_counts)))

rownames(pb_meta) <- pb_meta$GSM
pb_counts <- pb_counts[, rownames(pb_meta), drop = FALSE]

pb_meta$prognosis_group <- factor(
  pb_meta$prognosis_group,
  levels = c("Poor prognosis", "Favorable prognosis")
)

print(table(pb_meta$prognosis_group))

if (any(table(pb_meta$prognosis_group) < 2)) {
  warning("One prognosis group has fewer than 2 samples. Pseudobulk DEG may be underpowered.")
}

dds <- DESeqDataSetFromMatrix(
  countData = round(pb_counts),
  colData = pb_meta,
  design = ~ prognosis_group
)

dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
dds <- DESeq(dds)

res_hcc <- results(
  dds,
  contrast = c("prognosis_group", "Favorable prognosis", "Poor prognosis")
)

DEG_hcc <- as.data.frame(res_hcc)
DEG_hcc$gene <- rownames(DEG_hcc)
DEG_hcc <- DEG_hcc |> dplyr::arrange(padj)

write.csv(
  DEG_hcc,
  file.path(outdir, "HCC_cells", "Pseudobulk_DEG_HCC_cells_Favorable_vs_Poor_DESeq2.csv"),
  row.names = FALSE
)

geneList <- DEG_hcc$log2FoldChange
names(geneList) <- DEG_hcc$gene
geneList <- geneList[!is.na(geneList)]
geneList <- sort(geneList, decreasing = TRUE)

if (file.exists(kegg_gmt_file)) {
  kegg_gmt <- read.gmt(kegg_gmt_file)
  gsea_kegg <- GSEA(
    geneList,
    TERM2GENE = kegg_gmt,
    pvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 500,
    pAdjustMethod = "BH"
  )
  write.csv(as.data.frame(gsea_kegg), file.path(outdir, "GSEA", "GSEA_KEGG_HCC_cells_Favorable_vs_Poor.csv"), row.names = FALSE)

  p_kegg_dot <- dotplot(gsea_kegg, showCategory = 20) +
    ggtitle("KEGG GSEA in HCC cells: Favorable vs Poor prognosis")
  save_pdf(file.path(outdir, "GSEA", "FigS5_KEGG_GSEA_dotplot_HCC_cells.pdf"), p_kegg_dot, width = 9, height = 7)

  ppar_term <- safe_find_ppar_term(gsea_kegg)
  if (!is.na(ppar_term)) {
    res_df <- as.data.frame(gsea_kegg)
    idx <- match(ppar_term, res_df$ID)
    title_txt <- sprintf(
      "%s\nNES=%.3f, FDR=%.3g",
      res_df$Description[idx],
      res_df$NES[idx],
      res_df$p.adjust[idx]
    )
    p_ppar_gsea <- gseaplot2(gsea_kegg, geneSetID = ppar_term, title = title_txt)
    save_pdf(file.path(outdir, "GSEA", "FigS5_PPAR_pathway_GSEA_HCC_cells.pdf"), p_ppar_gsea, width = 6, height = 5)
  } else {
    warning("No PPAR-related pathway found in KEGG GSEA results.")
  }
} else {
  warning("KEGG GMT file not found. Skipped KEGG GSEA.")
}

# ===================== 7. T cell subpopulation analysis =====================
celltype_levels <- unique(as.character(sce$cell_type))

if ("T cells" %in% celltype_levels) {
  t_cell_labels <- "T cells"
} else {
  t_cell_labels <- celltype_levels[
    grepl("(^T cell|^T cells|CD4|CD8|Treg|Tconv|T/NK|NK/T)", celltype_levels, ignore.case = TRUE)
  ]
}

message("T-cell labels used for subclustering: ", paste(t_cell_labels, collapse = ", "))

if (length(t_cell_labels) == 0) {
  stop("No T-cell-related labels found. Please check cell_type labels in GSE151530_Info.csv.")
}

sce_t <- subset(sce, subset = cell_type %in% t_cell_labels)

sce_t <- NormalizeData(sce_t, normalization.method = "LogNormalize", scale.factor = 10000)
sce_t <- FindVariableFeatures(sce_t, selection.method = "vst", nfeatures = n_variable_features)
sce_t <- ScaleData(sce_t)
sce_t <- RunPCA(sce_t, features = VariableFeatures(sce_t))
p_elbow_t <- ElbowPlot(sce_t, ndims = 50)
save_pdf(file.path(outdir, "T_cells", "elbowplot_T_cells.pdf"), p_elbow_t, width = 7, height = 5)

sce_t <- RunUMAP(sce_t, dims = 1:n_pcs_tcell)
sce_t <- RunTSNE(sce_t, dims = 1:n_pcs_tcell)
sce_t <- FindNeighbors(sce_t, dims = 1:n_pcs_tcell)
sce_t <- FindClusters(sce_t, resolution = tcell_resolution)

p_t_tsne_cluster <- DimPlot(sce_t, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("tSNE by T cell cluster") + theme(plot.title = element_text(hjust = 0.5))
save_pdf(file.path(outdir, "T_cells", "FigS5_T_cell_tsne_by_cluster.pdf"), p_t_tsne_cluster, width = 7, height = 6)

# T cell marker genes.
tcell_marker_genes <- intersect(c(
  "CD3D", "CD3E", "CD4", "CD8A", "CD8B", "CCR7", "SELL", "TCF7", "IL7R",
  "NKG7", "GNLY", "GZMB", "PRF1", "IFNG", "PDCD1", "CTLA4", "LAG3", "TIGIT",
  "HAVCR2", "FOXP3", "IL2RA", "MKI67", "TOP2A", "TNFRSF9"
), rownames(sce_t))

p_t_marker_dot <- DotPlot(sce_t, features = tcell_marker_genes, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Representative marker genes for T cell clusters")
save_pdf(file.path(outdir, "T_cells", "FigS5_T_cell_marker_dotplot.pdf"), p_t_marker_dot, width = 11, height = 5)

p_t_marker_vln <- VlnPlot(sce_t, features = tcell_marker_genes, group.by = "seurat_clusters",
                          stack = TRUE, flip = TRUE, pt.size = 0) +
  ggtitle("T cell marker expression across clusters")
save_pdf(file.path(outdir, "T_cells", "T_cell_marker_vlnplot.pdf"), p_t_marker_vln, width = 10, height = 8)

# T Cell State Scores using TCellSI on average expression by T cell cluster.
av_t <- AverageExpression(sce_t, group.by = "seurat_clusters", assays = "RNA", slot = "data")
av_mat_t <- as.data.frame(av_t$RNA)
TCSS_scores <- TCSS_Calculate(av_mat_t)
write.csv(TCSS_scores, file.path(outdir, "T_cells", "FigS5_TCSS_scores_by_T_cell_cluster.csv"))

pdf(file.path(outdir, "T_cells", "FigS5_TCSS_heatmap_by_T_cell_cluster.pdf"), width = 7, height = 6)
pheatmap(TCSS_scores, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, main = "T Cell State Scores by cluster")
dev.off()

# T cell subset proportions by prognosis group.
tab_t <- table(sce_t$seurat_clusters, sce_t$prognosis_group)
tab_t_prop <- prop.table(tab_t, margin = 2)
write.csv(as.data.frame(tab_t_prop), file.path(outdir, "T_cells", "FigS5_T_cell_cluster_proportion_by_group.csv"), row.names = FALSE)

df_t_prop <- melt(tab_t_prop)
colnames(df_t_prop) <- c("T_cell_cluster", "Group", "Proportion")
df_t_prop$T_cell_cluster <- factor(df_t_prop$T_cell_cluster)

p_t_stack <- ggplot(df_t_prop, aes(x = Group, y = Proportion, fill = T_cell_cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  labs(title = "T cell subset proportions by prognosis group", x = NULL, y = "Proportion") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_pdf(file.path(outdir, "T_cells", "FigS5_T_cell_subset_proportions_by_prognosis.pdf"), p_t_stack, width = 6, height = 5)

# More conservative per-sample T cell proportion table for statistics.
t_per_sample <- sce_t@meta.data |>
  group_by(GSM, prognosis_group, seurat_clusters) |>
  summarise(count = n(), .groups = "drop") |>
  group_by(GSM) |>
  mutate(prop = count / sum(count)) |>
  ungroup()
write.csv(t_per_sample, file.path(outdir, "T_cells", "T_cell_cluster_proportion_per_sample.csv"), row.names = FALSE)

p_t_box <- ggplot(t_per_sample, aes(x = prognosis_group, y = prop, fill = prognosis_group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1) +
  facet_wrap(~ seurat_clusters, scales = "free_y") +
  scale_fill_manual(values = group_cols) +
  theme_bw() +
  labs(title = "T cell cluster proportions per sample", x = NULL, y = "Proportion") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
save_pdf(file.path(outdir, "T_cells", "T_cell_cluster_proportion_per_sample_boxplot.pdf"), p_t_box, width = 8, height = 6)

# Cluster marker genes for T cells.
tcell_cluster_markers <- FindAllMarkers(
  sce_t,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(tcell_cluster_markers, file.path(outdir, "T_cells", "T_cell_cluster_markers.csv"), row.names = FALSE)

# Fig. S5D: heatmap of differentially expressed genes among T cell subsets.
tcell_top_heatmap_markers <- tcell_cluster_markers |>
  dplyr::filter(p_val_adj < 0.05) |>
  dplyr::group_by(cluster) |>
  dplyr::arrange(dplyr::desc(avg_log2FC), .by_group = TRUE) |>
  dplyr::slice_head(n = 20) |>
  dplyr::ungroup()

if (nrow(tcell_top_heatmap_markers) > 0) {
  tcell_heatmap_genes <- unique(tcell_top_heatmap_markers$gene)
  avg_t_marker_expr <- AverageExpression(
    sce_t,
    features = tcell_heatmap_genes,
    group.by = "seurat_clusters",
    assays = "RNA",
    slot = "data"
  )$RNA
  avg_t_marker_expr <- as.matrix(avg_t_marker_expr)
  avg_t_marker_expr_z <- t(scale(t(avg_t_marker_expr)))
  avg_t_marker_expr_z[is.na(avg_t_marker_expr_z)] <- 0
  write.csv(avg_t_marker_expr_z, file.path(outdir, "T_cells", "FigS5D_T_cell_subset_DEG_heatmap_matrix_zscore.csv"))

  pdf(file.path(outdir, "T_cells", "FigS5D_T_cell_subset_DEG_heatmap.pdf"), width = 7, height = 10)
  pheatmap(
    avg_t_marker_expr_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 6,
    fontsize_col = 12,
    border_color = NA,
    main = "Differentially expressed genes among T cell subsets"
  )
  dev.off()
} else {
  warning("No significant T cell subset DEG markers found for Fig. S5D under the current thresholds.")
}

# Save final objects.
saveRDS(sce, file.path(outdir, "sce_global_FigureS5_final.rds"))
saveRDS(sce_hcc, file.path(outdir, "sce_HCC_cells_FigureS5_final.rds"))
saveRDS(sce_t, file.path(outdir, "sce_T_cells_FigureS5_final.rds"))

message("Figure S5 scRNA-seq analysis completed. Results saved to: ", outdir)
