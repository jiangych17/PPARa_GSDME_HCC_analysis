## =========================
## 1. 安装/加载需要的包
## =========================
library(readxl)
library(dplyr)
library(tibble)
library(fgsea)
library(ggplot2)
library(data.table)


## =========================
## 2. 设置文件路径
## =========================

expr_file <- "250508.txt"
geneset_file <- "GO.xls"

outdir <- "GSEA_GO_cell_death_results"
dir.create(outdir, showWarnings = FALSE)


## =========================
## 3. 读取表达矩阵
## =========================
## 假设 250508.txt 第一列是基因名，后面 9 列是样本
## 如果你的文件是 tab 分隔，用 sep = "\t"
## 如果是逗号分隔，改成 sep = ","

expr <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

## 看一下数据结构
head(expr)
dim(expr)

## 默认第一列为基因名
gene_col <- colnames(expr)[1]

expr_mat <- expr %>%
  distinct(.data[[gene_col]], .keep_all = TRUE) %>%
  column_to_rownames(gene_col)

## 确保全部是 numeric
expr_mat <- as.data.frame(expr_mat)
expr_mat[] <- lapply(expr_mat, function(x) as.numeric(as.character(x)))
expr_mat <- as.matrix(expr_mat)

## 去掉全是 NA 或无表达的基因
expr_mat <- expr_mat[rowSums(is.na(expr_mat)) == 0, ]
expr_mat <- expr_mat[rowSums(expr_mat) > 0, ]

## 检查样本数量
cat("Number of samples:", ncol(expr_mat), "\n")
print(colnames(expr_mat))


## =========================
## 4. 设置分组
## =========================
## 你说前三个样本为 Non-FILM，后六个为 FILM

group <- c(rep("Non-FILM", 3), rep("FILM", 6))
names(group) <- colnames(expr_mat)

print(group)


## =========================
## 5. 计算 ranking metric
## =========================
## 用 log2FC 作为 GSEA 排序指标
## 这里定义：
## log2FC = FILM - Non-FILM
## 所以：
##   NES > 0 表示该 gene set 在 FILM 富集
##   NES < 0 表示该 gene set 在 Non-FILM 富集，也就是 FILM 中下降

nonfilm_samples <- names(group)[group == "Non-FILM"]
film_samples <- names(group)[group == "FILM"]

## 如果你的矩阵是 raw count，建议先 log2(x + 1)
## 如果已经是 TPM/FPKM/log normalized，也可以保留。
## 这里统一用 log2(x + 1)，比较稳。
expr_log <- log2(expr_mat + 1)

mean_nonfilm <- rowMeans(expr_log[, nonfilm_samples, drop = FALSE])
mean_film <- rowMeans(expr_log[, film_samples, drop = FALSE])

log2FC <- mean_film - mean_nonfilm
names(log2FC) <- rownames(expr_log)

## 去掉 NA、重复、无限值
ranks <- log2FC[is.finite(log2FC)]
ranks <- sort(ranks, decreasing = TRUE)

head(ranks)
tail(ranks)


## =========================
## 6. 读取 GO.xls 里的 5 个 gene sets
## =========================

sheet_names <- excel_sheets(geneset_file)
print(sheet_names)

pathways <- list()

for (s in sheet_names) {
  gs <- read_excel(geneset_file, sheet = s, col_names = FALSE)[[1]]
  gs <- unique(na.omit(as.character(gs)))
  gs <- gs[gs != ""]
  
  ## 保持 mouse gene symbol 格式
  pathways[[s]] <- gs
}

## 查看每个 gene set 基因数
sapply(pathways, length)


## =========================
## 7. 检查 gene set 和表达矩阵基因的重叠
## =========================

overlap_info <- data.frame(
  pathway = names(pathways),
  geneset_size = sapply(pathways, length),
  overlap_with_matrix = sapply(pathways, function(x) sum(x %in% names(ranks)))
)

print(overlap_info)

write.csv(
  overlap_info,
  file = file.path(outdir, "GeneSet_overlap_with_matrix.csv"),
  row.names = FALSE
)


## =========================
## 8. 运行 fgsea
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
  mutate(
    direction = ifelse(NES > 0, "Enriched in FILM", "Enriched in Non-FILM / suppressed in FILM")
  )

print(fgsea_res)

## 保存完整结果
fwrite(
  fgsea_res,
  file = file.path(outdir, "Focused_GO_cell_death_GSEA_results.csv")
)


## =========================
## 9. 画 summary dot plot
## =========================

plot_df <- fgsea_res %>%
  mutate(
    pathway = factor(pathway, levels = pathway[order(NES)]),
    neg_log10_padj = -log10(padj + 1e-300)
  )

p <- ggplot(plot_df, aes(x = NES, y = pathway)) +
  geom_point(aes(size = neg_log10_padj, color = NES)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 12) +
  labs(
    x = "Normalized enrichment score (NES)",
    y = NULL,
    size = "-log10(FDR)",
    color = "NES",
    title = "Focused GO cell death-related GSEA",
    subtitle = "FILM vs Non-FILM; NES < 0 indicates suppression in FILM"
  )

print(p)

ggsave(
  filename = file.path(outdir, "Focused_GO_cell_death_GSEA_dotplot.pdf"),
  plot = p,
  width = 7,
  height = 4
)

ggsave(
  filename = file.path(outdir, "Focused_GO_cell_death_GSEA_dotplot.png"),
  plot = p,
  width = 7,
  height = 4,
  dpi = 300
)


## =========================
## 10. 单独输出每个通路的 enrichment plot
## =========================

for (pw in names(pathways)) {
  p_enrich <- plotEnrichment(pathways[[pw]], ranks) +
    labs(title = pw)
  
  ggsave(
    filename = file.path(outdir, paste0(pw, "_enrichment_plot.pdf")),
    plot = p_enrich,
    width = 6,
    height = 4
  )
  
  ggsave(
    filename = file.path(outdir, paste0(pw, "_enrichment_plot.png")),
    plot = p_enrich,
    width = 6,
    height = 4,
    dpi = 300
  )
}


## =========================
## 11. 输出 leading-edge genes
## =========================

leading_edge <- fgsea_res %>%
  select(pathway, leadingEdge) %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

write.csv(
  leading_edge,
  file = file.path(outdir, "Leading_edge_genes.csv"),
  row.names = FALSE
)

cat("Done! Results saved in:", outdir, "\n")


## =========================
## 单独输出每个 pathway 的 GSEA enrichment plot
## =========================

library(ggplot2)
library(fgsea)
library(dplyr)

single_plot_dir <- file.path(outdir, "single_GSEA_plots")
dir.create(single_plot_dir, showWarnings = FALSE)

## 如果 pathway 名字需要统一显示，可以在这里改
pathway_label_map <- c(
  "Pyroptosis" = "PYROPTOSIS",
  "apoptosis" = "APOPTOSIS",
  "AUTOPHAGY" = "AUTOPHAGY",
  "NECROPTOTIC_PROCESS" = "NECROPTOTIC_PROCESS",
  "FERROPTOSIS" = "FERROPTOSIS"
)

## 按你想展示的顺序排列
plot_order <- c(
  "AUTOPHAGY",
  "apoptosis",
  "FERROPTOSIS",
  "NECROPTOTIC_PROCESS",
  "Pyroptosis"
)

## 如果你的 pathway 名称和上面不完全一致，先看一下：
print(fgsea_res$pathway)

for (pw in plot_order) {
  
  if (!pw %in% names(pathways)) {
    warning(paste("Pathway not found:", pw))
    next
  }
  
  res_pw <- fgsea_res %>% filter(pathway == pw)
  
  display_name <- ifelse(
    pw %in% names(pathway_label_map),
    pathway_label_map[[pw]],
    pw
  )
  
  label_text <- paste0(
    "ES=", round(res_pw$ES, 3), "\n",
    "NES=", round(res_pw$NES, 3), "\n",
    "FDR=", signif(res_pw$padj, 3)
  )
  
  p <- plotEnrichment(pathways[[pw]], ranks) +
    labs(
      title = display_name,
      subtitle = "FILM vs Non-FILM; NES < 0 indicates suppression in FILM",
      x = "Rank in ordered dataset",
      y = "Enrichment score"
    ) +
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
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  ggsave(
    filename = file.path(single_plot_dir, paste0(display_name, "_GSEA.pdf")),
    plot = p,
    width = 6,
    height = 4
  )
  
  ggsave(
    filename = file.path(single_plot_dir, paste0(display_name, "_GSEA.png")),
    plot = p,
    width = 6,
    height = 4,
    dpi = 300
  )
}

## =========================
## 拼图输出：类似 Figure S2B
## =========================

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
library(patchwork)

plot_list <- list()

for (pw in plot_order) {
  
  if (!pw %in% names(pathways)) {
    next
  }
  
  res_pw <- fgsea_res %>% filter(pathway == pw)
  
  display_name <- ifelse(
    pw %in% names(pathway_label_map),
    pathway_label_map[[pw]],
    pw
  )
  
  label_text <- paste0(
    "ES=", round(res_pw$ES, 3), "\n",
    "NES=", round(res_pw$NES, 3), "\n",
    "FDR=", signif(res_pw$padj, 3)
  )
  
  p <- plotEnrichment(pathways[[pw]], ranks) +
    labs(
      title = display_name,
      x = NULL,
      y = "Enrichment Score"
    ) +
    annotate(
      "text",
      x = length(ranks) * 0.04,
      y = max(res_pw$ES, 0.05),
      label = label_text,
      hjust = 0,
      vjust = 1,
      size = 3.3
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  plot_list[[display_name]] <- p
}

combined_plot <- wrap_plots(plot_list, nrow = 1) +
  plot_annotation(
    title = "HCC tissues from mouse",
    subtitle = "FILM vs Non-FILM; negative NES indicates suppression in FILM"
  )

print(combined_plot)

ggsave(
  filename = file.path(outdir, "Focused_GO_cell_death_GSEA_combined.pdf"),
  plot = combined_plot,
  width = 15,
  height = 3.5
)

ggsave(
  filename = file.path(outdir, "Focused_GO_cell_death_GSEA_combined.png"),
  plot = combined_plot,
  width = 15,
  height = 3.5,
  dpi = 300
)



## =========================
## 12. 使用 clusterProfiler::GSEA 重新跑，供 gseaVis 作图
## =========================

packages2 <- c("clusterProfiler", "GseaVis", "cowplot", "ggplot2", "dplyr")

for (pkg in packages2) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("clusterProfiler")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else if (pkg == "GseaVis") {
      if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
      devtools::install_github("junjunlab/GseaVis")
    } else {
      install.packages(pkg)
    }
  }
}

library(clusterProfiler)
library(GseaVis)
library(cowplot)
library(ggplot2)
library(dplyr)

## 将 pathways list 转成 TERM2GENE 格式
term2gene <- do.call(
  rbind,
  lapply(names(pathways), function(pw) {
    data.frame(
      term = pw,
      gene = unique(pathways[[pw]]),
      stringsAsFactors = FALSE
    )
  })
)

## 确保 geneList 是 decreasing numeric vector
geneList <- ranks
geneList <- sort(geneList, decreasing = TRUE)

## clusterProfiler::GSEA
## pAdjustMethod = "BH" 是常规 FDR 校正方式
gsea_cp <- GSEA(
  geneList = geneList,
  TERM2GENE = term2gene,
  minGSSize = 5,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  eps = 0,
  verbose = FALSE,
  seed = TRUE
)

gsea_cp_df <- as.data.frame(gsea_cp)

## 保存结果
write.csv(
  gsea_cp_df,
  file = file.path(outdir, "clusterProfiler_GSEA_GO_cell_death_results.csv"),
  row.names = FALSE
)

print(gsea_cp_df)

## =========================
## 自定义 Nature 风格 GSEA 图：显示 ES / NES / FDR，不显示 P value
## =========================

make_gseaVis_plot <- function(gsea_object,
                              gsea_df,
                              geneSetID,
                              title_name,
                              label_size = 9) {
  
  ## 提取该通路结果
  res_pw <- gsea_df[gsea_df$ID == geneSetID, ]
  
  ## clusterProfiler 结果里通常是这些列：
  ## enrichmentScore, NES, p.adjust
  label_text <- paste0(
    "ES = ", round(res_pw$enrichmentScore, 3), "\n",
    "NES = ", round(res_pw$NES, 2), "\n",
    "FDR = ", signif(res_pw$p.adjust, 3)
  )
  
  ## 先画基础 GSEA 图
  ## 关键修改：
  ## 1. addPval = FALSE
  ## 2. subRatio 第三个值从 0.3 改成 0.15，rank list 高度减半
  base_p <- gseaNb(
    object = gsea_object,
    geneSetID = geneSetID,
    curveCol = "#5EABD6",
    htHeight = 1,
    htCol = c("#1450A3", "#8C1007"),
    htAlpha = 0.8,
    subRatio = c(0.58, 0.12, 0.15),
    rankCol = "#FFC436",
    addPval = FALSE,
    addZeroLine = TRUE
  ) +
    ggtitle(title_name) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 12
      ),
      text = element_text(
        face = "plain"
      )
    )
  
  ## 用 cowplot 叠加 ES/NES/FDR
  ## fontface = "plain" 确保不是斜体
  final_p <- cowplot::ggdraw(base_p) +
    cowplot::draw_label(
      label = label_text,
      x = 0.78,
      y = 0.82,
      hjust = 0,
      vjust = 1,
      size = label_size,
      fontface = "plain"
    )
  
  return(final_p)
}

## =========================
## 13. 使用 gseaVis n
#画 Nature 风格单图
## =========================

gseaVis_dir <- file.path(outdir, "GseaVis_Nature_style")
dir.create(gseaVis_dir, showWarnings = FALSE)

## 名称映射
label_map <- c(
  "Pyroptosis" = "PYROPTOSIS",
  "apoptosis" = "APOPTOSIS",
  "AUTOPHAGY" = "AUTOPHAGY",
  "NECROPTOTIC_PROCESS" = "NECROPTOTIC PROCESS",
  "FERROPTOSIS" = "FERROPTOSIS"
)

plot_order <- c(
  "Pyroptosis",
  "apoptosis",
  "AUTOPHAGY",
  "NECROPTOTIC_PROCESS",
  "FERROPTOSIS"
)

## 只保留实际存在于 gsea_cp 的 pathway
plot_order <- plot_order[plot_order %in% gsea_cp_df$ID]

gseaVis_dir <- file.path(outdir, "GseaVis_Nature_style")
dir.create(gseaVis_dir, showWarnings = FALSE)

for (pw in plot_order) {
  
  display_name <- ifelse(
    pw %in% names(label_map),
    label_map[[pw]],
    pw
  )
  
  p <- make_gseaVis_plot(
    gsea_object = gsea_cp,
    gsea_df = gsea_cp_df,
    geneSetID = pw,
    title_name = display_name,
    label_size = 10
  )
  
  ggsave(
    filename = file.path(gseaVis_dir, paste0(display_name, "_Nature_style_GSEA_ES_NES_FDR.pdf")),
    plot = p,
    width = 5.2,
    height = 4.0
  )
  
  ggsave(
    filename = file.path(gseaVis_dir, paste0(display_name, "_Nature_style_GSEA_ES_NES_FDR.png")),
    plot = p,
    width = 5.2,
    height = 4.0,
    dpi = 600
  )
}

## =========================
## 14. 合并 5 个 Nature 风格 GSEA 图
## =========================
## =========================
## 14. 合并 5 个 Nature 风格 GSEA 图
## =========================

gsea_list <- lapply(plot_order, function(pw) {
  
  display_name <- ifelse(
    pw %in% names(label_map),
    label_map[[pw]],
    pw
  )
  
  make_gseaVis_plot(
    gsea_object = gsea_cp,
    gsea_df = gsea_cp_df,
    geneSetID = pw,
    title_name = display_name,
    label_size = 8
  )
})

names(gsea_list) <- plot_order

combined_gsea <- cowplot::plot_grid(
  plotlist = gsea_list,
  ncol = 5,
  align = "hv"
)

ggsave(
  filename = file.path(gseaVis_dir, "Combined_5_GO_death_GSEA_Nature_style_ES_NES_FDR.pdf"),
  plot = combined_gsea,
  width = 22,
  height = 4.0
)

ggsave(
  filename = file.path(gseaVis_dir, "Combined_5_GO_death_GSEA_Nature_style_ES_NES_FDR.png"),
  plot = combined_gsea,
  width = 22,
  height = 4.0,
  dpi = 600
)

print(combined_gsea)
## =========================
## 15. 美化整体 summary 图
## =========================


summary_df <- gsea_cp_df %>%
  mutate(
    display_name = dplyr::case_when(
      ID %in% c("Pyroptosis", "PYROPTOSIS", "GOBP_PYROPTOSIS") ~ "Pyroptosis",
      ID %in% c("apoptosis", "APOPTOSIS", "GOBP_EXECUTION_PHASE_OF_APOPTOSIS") ~ "Apoptosis",
      ID %in% c("AUTOPHAGY", "Autophagy", "GOBP_REGULATION_OF_AUTOPHAGY") ~ "Autophagy",
      ID %in% c("NECROPTOTIC_PROCESS", "Necroptotic process", "GOBP_NECROPTOTIC_PROCESS") ~ "Necroptotic process",
      ID %in% c("FERROPTOSIS", "Ferroptosis", "GOBP_FERROPTOSIS") ~ "Ferroptosis",
      TRUE ~ ID
    ),
    abs_NES = abs(NES),
    neg_log10_fdr = -log10(p.adjust + 1e-300)
  ) %>%
  arrange(desc(abs_NES))

## ggplot 的 y 轴 factor 顺序通常第一个在最下面
## 所以这里用 rev()，让 abs_NES 最大的显示在最上面
summary_df$display_name <- factor(
  summary_df$display_name,
  levels = rev(summary_df$display_name)
)

p_summary <- ggplot(summary_df, aes(x = abs_NES, y = display_name)) +
  geom_segment(
    aes(x = 0, xend = abs_NES, y = display_name, yend = display_name),
    linewidth = 2.5,
    lineend = "round",
    color = "#7FA9D6"
  ) +
  geom_point(
    aes(size = neg_log10_fdr),
    shape = 21,
    fill = "#7FA9D6",
    color = "#5C88BE",
    stroke = 0.3
  ) +
  scale_x_continuous(
    breaks = seq(
      0,
      ceiling(max(summary_df$abs_NES, na.rm = TRUE)),
      by = 0.5
    ),
    limits = c(
      0,
      ceiling(max(summary_df$abs_NES, na.rm = TRUE)) + 0.2
    )
  ) +
  scale_size_continuous(
    range = c(3, 8),
    breaks = c(0.5, 1.0, 1.5, 2.0),
    name = "-log10(FDR)"
  ) +
  theme_classic(base_size = 13) +
  labs(
    x = "|Normalized enrichment score|",
    y = NULL,
    title = "Focused GO cell death-related GSEA",
    subtitle = "FILM vs Non-FILM; bar length shows absolute NES"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11),
    legend.position = "right"
  )

print(p_summary)

ggsave(
  filename = file.path(gseaVis_dir, "Summary_GO_death_GSEA_absNES_lollipop.pdf"),
  plot = p_summary,
  width = 7,
  height = 4.5
)

ggsave(
  filename = file.path(gseaVis_dir, "Summary_GO_death_GSEA_absNES_lollipop.png"),
  plot = p_summary,
  width = 7,
  height = 4.5,
  dpi = 600
)

#------------------------------HUman-----------------------------
#------------------------------HUman-----------------------------
#------------------------------HUman-----------------------------
## =========================
## Human 250514 focused GO cell death GSEA
## 前9个样本 = Non-FILM
## 后9个样本 = FILM
## Gene set file = GOH.xls
## =========================

## =========================
## 1. 安装/加载需要的包
## =========================


library(readxl)
library(dplyr)
library(tibble)
library(fgsea)
library(ggplot2)
library(data.table)
library(clusterProfiler)
library(GseaVis)
library(cowplot)


## =========================
## 2. 设置文件路径
## =========================

expr_file <- "260514.txt"
geneset_file <- "GOH.xls"

outdir <- "GSEA_GOH_cell_death_results_human"
dir.create(outdir, showWarnings = FALSE)


## =========================
## 3. 读取表达矩阵
## =========================

expr <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Raw matrix dimension:\n")
print(dim(expr))
cat("First few columns:\n")
print(head(colnames(expr)))

## 默认第一列是 gene symbol
gene_col <- colnames(expr)[1]

expr_mat <- expr %>%
  distinct(.data[[gene_col]], .keep_all = TRUE) %>%
  column_to_rownames(gene_col)

expr_mat <- as.data.frame(expr_mat)
expr_mat[] <- lapply(expr_mat, function(x) as.numeric(as.character(x)))
expr_mat <- as.matrix(expr_mat)

## 去掉 NA 和全 0 基因
expr_mat <- expr_mat[rowSums(is.na(expr_mat)) == 0, ]
expr_mat <- expr_mat[rowSums(expr_mat) > 0, ]

cat("Clean matrix dimension:\n")
print(dim(expr_mat))
cat("Samples:\n")
print(colnames(expr_mat))


## =========================
## 4. 设置分组
## =========================
## 前9个 Non-FILM，后9个 FILM

stopifnot(ncol(expr_mat) == 18)

group <- c(rep("Non-FILM", 9), rep("FILM", 9))
names(group) <- colnames(expr_mat)

print(group)


## =========================
## 5. 计算 ranking metric
## =========================
## log2FC = FILM - Non-FILM
## NES < 0 表示该通路在 FILM 中被抑制

nonfilm_samples <- names(group)[group == "Non-FILM"]
film_samples <- names(group)[group == "FILM"]

expr_log <- log2(expr_mat + 1)

mean_nonfilm <- rowMeans(expr_log[, nonfilm_samples, drop = FALSE])
mean_film <- rowMeans(expr_log[, film_samples, drop = FALSE])

log2FC <- mean_film - mean_nonfilm
names(log2FC) <- rownames(expr_log)

ranks <- log2FC[is.finite(log2FC)]
ranks <- sort(ranks, decreasing = TRUE)

cat("Ranking metric summary:\n")
print(summary(ranks))
cat("Top genes enriched in FILM:\n")
print(head(ranks, 10))
cat("Top genes enriched in Non-FILM / suppressed in FILM:\n")
print(tail(ranks, 10))


## 输出 rnk 文件，方便用 GSEAPreranked 复现
rnk_df <- data.frame(
  Gene = names(ranks),
  RankMetric = as.numeric(ranks)
)

write.table(
  rnk_df,
  file = file.path(outdir, "Human_FILM_vs_NonFILM_log2FC.rnk"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


## =========================
## 6. 读取 GOH.xls 里的 gene sets
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


## 输出 GMT 文件，方便用 GSEAPreranked 复现
gmt_file <- file.path(outdir, "GOH_cell_death_5sets.gmt")
con <- file(gmt_file, open = "wt")

for (pw in names(pathways)) {
  genes <- unique(pathways[[pw]])
  genes <- genes[!is.na(genes) & genes != ""]
  line <- paste(c(pw, "na", genes), collapse = "\t")
  writeLines(line, con = con)
}

close(con)


## =========================
## 7. 检查 gene set 和表达矩阵基因重叠
## =========================

overlap_info <- data.frame(
  pathway = names(pathways),
  geneset_size = sapply(pathways, length),
  overlap_with_matrix = sapply(pathways, function(x) sum(x %in% names(ranks)))
)

print(overlap_info)

write.csv(
  overlap_info,
  file = file.path(outdir, "GeneSet_overlap_with_matrix.csv"),
  row.names = FALSE
)


## =========================
## 8. fgsea 分析
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
  mutate(
    direction = ifelse(
      NES > 0,
      "Enriched in FILM",
      "Enriched in Non-FILM / suppressed in FILM"
    )
  )

print(fgsea_res)

fwrite(
  fgsea_res,
  file = file.path(outdir, "fgsea_Focused_GOH_cell_death_results.csv")
)

leading_edge <- fgsea_res %>%
  select(pathway, leadingEdge) %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

write.csv(
  leading_edge,
  file = file.path(outdir, "fgsea_Leading_edge_genes.csv"),
  row.names = FALSE
)


## =========================
## 9. clusterProfiler::GSEA
## 供 GseaVis 画 Nature 风格图
## =========================

term2gene <- do.call(
  rbind,
  lapply(names(pathways), function(pw) {
    data.frame(
      term = pw,
      gene = unique(pathways[[pw]]),
      stringsAsFactors = FALSE
    )
  })
)

geneList <- ranks
geneList <- sort(geneList, decreasing = TRUE)

gsea_cp <- GSEA(
  geneList = geneList,
  TERM2GENE = term2gene,
  minGSSize = 5,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  eps = 0,
  verbose = FALSE,
  seed = TRUE
)

gsea_cp_df <- as.data.frame(gsea_cp)

write.csv(
  gsea_cp_df,
  file = file.path(outdir, "clusterProfiler_GSEA_GOH_cell_death_results.csv"),
  row.names = FALSE
)

print(gsea_cp_df)


## =========================
## 10. Nature 风格 GSEA 图函数
## 显示 ES / NES / FDR
## =========================

make_gseaVis_plot <- function(gsea_object,
                              gsea_df,
                              geneSetID,
                              title_name,
                              label_size = 9) {
  
  res_pw <- gsea_df[gsea_df$ID == geneSetID, ]
  
  label_text <- paste0(
    "ES = ", round(res_pw$enrichmentScore, 3), "\n",
    "NES = ", round(res_pw$NES, 2), "\n",
    "FDR = ", signif(res_pw$p.adjust, 3)
  )
  
  base_p <- gseaNb(
    object = gsea_object,
    geneSetID = geneSetID,
    curveCol = "#5EABD6",
    htHeight = 1,
    htCol = c("#1450A3", "#8C1007"),
    htAlpha = 0.8,
    subRatio = c(0.58, 0.12, 0.15),
    rankCol = "#FFC436",
    addPval = FALSE,
    addZeroLine = TRUE
  ) +
    ggtitle(title_name) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 12
      ),
      text = element_text(face = "plain")
    )
  
  final_p <- cowplot::ggdraw(base_p) +
    cowplot::draw_label(
      label = label_text,
      x = 0.78,
      y = 0.82,
      hjust = 0,
      vjust = 1,
      size = label_size,
      fontface = "plain"
    )
  
  return(final_p)
}


## =========================
## 11. 输出单个 Nature 风格 GSEA 图
## =========================
## =========================
## 11. 输出单个 Nature 风格 GSEA 图
## =========================

gseaVis_dir <- file.path(outdir, "GseaVis_Nature_style")
dir.create(gseaVis_dir, showWarnings = FALSE)

## 使用 GOH.xls 的真实 pathway ID
label_map <- c(
  "GOBP_PYROPTOTIC_INFLAMMATORY_RE" = "PYROPTOSIS",
  "GOBP_EXECUTION_PHASE_OF_APOPTOS" = "APOPTOSIS",
  "GOBP_REGULATION_OF_AUTOPHAGY" = "AUTOPHAGY",
  "GOBP_NECROPTOTIC_PROCESS" = "NECROPTOTIC PROCESS",
  "FERROPTOSIS" = "FERROPTOSIS"
)

plot_order <- c(
  "GOBP_PYROPTOTIC_INFLAMMATORY_RE",
  "GOBP_EXECUTION_PHASE_OF_APOPTOS",
  "GOBP_REGULATION_OF_AUTOPHAGY",
  "GOBP_NECROPTOTIC_PROCESS",
  "FERROPTOSIS"
)

plot_order <- plot_order[plot_order %in% gsea_cp_df$ID]

cat("Pathways in gsea_cp_df:\n")
print(gsea_cp_df$ID)

cat("Pathways selected for plotting:\n")
print(plot_order)

cat("Missing pathways:\n")
print(setdiff(names(label_map), gsea_cp_df$ID))


## 单独输出每个通路
for (pw in plot_order) {
  
  display_name <- label_map[[pw]]
  
  p <- make_gseaVis_plot(
    gsea_object = gsea_cp,
    gsea_df = gsea_cp_df,
    geneSetID = pw,
    title_name = display_name,
    label_size = 8
  )
  
  ggsave(
    filename = file.path(gseaVis_dir, paste0(display_name, "_Human_Nature_style_GSEA_ES_NES_FDR.pdf")),
    plot = p,
    width = 5.2,
    height = 4.0
  )
  
  ggsave(
    filename = file.path(gseaVis_dir, paste0(display_name, "_Human_Nature_style_GSEA_ES_NES_FDR.png")),
    plot = p,
    width = 5.2,
    height = 4.0,
    dpi = 600
  )
}


## =========================
## 12. 合并 Nature 风格 GSEA 图
## =========================

gsea_list <- lapply(plot_order, function(pw) {
  
  display_name <- label_map[[pw]]
  
  make_gseaVis_plot(
    gsea_object = gsea_cp,
    gsea_df = gsea_cp_df,
    geneSetID = pw,
    title_name = display_name,
    label_size = 6
  )
})

names(gsea_list) <- label_map[plot_order]

combined_gsea <- cowplot::plot_grid(
  plotlist = gsea_list,
  ncol = length(gsea_list),
  align = "hv"
)

ggsave(
  filename = file.path(gseaVis_dir, "Combined_Human_GOH_death_GSEA_Nature_style_ES_NES_FDR_fixed.pdf"),
  plot = combined_gsea,
  width = 5.0 * length(gsea_list),
  height = 4.0
)

ggsave(
  filename = file.path(gseaVis_dir, "Combined_Human_GOH_death_GSEA_Nature_style_ES_NES_FDR_fixed.png"),
  plot = combined_gsea,
  width = 5.0 * length(gsea_list),
  height = 4.0,
  dpi = 600
)

print(combined_gsea)
## =========================
## 13. 美化 summary lollipop 图
## =========================
## 适配 GOH 人基因集真实名称
## =========================
summary_df <- gsea_cp_df %>%
  mutate(
    display_name = dplyr::case_when(
      ID == "GOBP_PYROPTOTIC_INFLAMMATORY_RE" ~ "Pyroptosis",
      ID == "GOBP_EXECUTION_PHASE_OF_APOPTOS" ~ "Apoptosis",
      ID == "GOBP_REGULATION_OF_AUTOPHAGY" ~ "Autophagy",
      ID == "GOBP_NECROPTOTIC_PROCESS" ~ "Necroptotic process",
      ID == "FERROPTOSIS" ~ "Ferroptosis",
      TRUE ~ ID
    ),
    direction = ifelse(NES > 0, "Enriched in FILM", "Suppressed in FILM"),
    abs_NES = abs(NES),          # 新增：NES 绝对值
    neg_log10_fdr = -log10(p.adjust + 1e-300)
  )

summary_df <- summary_df %>%
  arrange(desc(abs_NES))

summary_df$display_name <- factor(
  summary_df$display_name,
  levels = rev(summary_df$display_name)
)

print(summary_df[, c("ID", "display_name", "NES", "abs_NES", "pvalue", "p.adjust", "neg_log10_fdr")])

x_min <- 0
x_max <- ceiling(max(summary_df$abs_NES, na.rm = TRUE) * 10) / 10 + 0.2

p_summary <- ggplot(summary_df, aes(x = abs_NES, y = display_name)) +
  geom_segment(
    aes(x = 0, xend = abs_NES, y = display_name, yend = display_name),
    linewidth = 2.5,
    lineend = "round",
    color = "#7FA9D6"
  ) +
  geom_point(
    aes(size = neg_log10_fdr),
    shape = 21,
    fill = "#7FA9D6",
    color = "#5C88BE",
    stroke = 0.3
  ) +
  scale_x_continuous(
    breaks = seq(0, ceiling(x_max), by = 0.5),
    limits = c(x_min, x_max)
  ) +
  scale_size_continuous(
    range = c(3, 8),
    breaks = c(0.5, 1.0, 1.5, 2.0),
    name = "-log10(FDR)"
  ) +
  theme_classic(base_size = 13) +
  labs(
    x = "|Normalized enrichment score| (absolute NES)",
    y = NULL,
    title = "Focused GO cell death-related GSEA",
    subtitle = "Human HCC tissues: FILM vs Non-FILM; bar length shows absolute NES"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11),
    legend.position = "right"
  )

print(p_summary)

############################################################
# Immune deconvolution analysis
# Human tissues: xCell, MCPcounter, Quantiseq, ImmunCellAI
# Mouse tissues: ImmunCellAI only
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

# =========================
# 1. Helper functions
# =========================

read_expression_matrix <- function(expr_file) {
  expr <- fread(expr_file) |> as.data.frame()
  rownames(expr) <- expr[[1]]
  expr[[1]] <- NULL
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  expr <- expr[rowSums(is.na(expr)) == 0, , drop = FALSE]
  expr <- expr[rowSums(expr) > 0, , drop = FALSE]
  return(expr)
}

make_group_annotation <- function(samples, group_vector) {
  anno <- data.frame(
    sample = samples,
    group = group_vector,
    row.names = samples,
    check.names = FALSE
  )
  return(anno)
}

plot_deconv_heatmap <- function(deconv_mat, annotation_col, outfile) {
  common_samples <- intersect(colnames(deconv_mat), rownames(annotation_col))
  deconv_mat <- deconv_mat[, common_samples, drop = FALSE]
  annotation_col <- annotation_col[common_samples, , drop = FALSE]
  
  pdf(outfile, width = 10, height = 8)
  pheatmap(
    deconv_mat,
    scale = "row",
    annotation_col = annotation_col,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 8,
    main = basename(outfile)
  )
  dev.off()
}

# =========================
# 2. Human immune deconvolution
# =========================

human_expr_file <- "260514.txt"   # rows = genes, columns = 18 human samples
human_outdir <- "Human_immune_deconvolution_results"
dir.create(human_outdir, showWarnings = FALSE, recursive = TRUE)

human_expr <- read_expression_matrix(human_expr_file)

# 前9个 Non-FILM，后9个 FILM；如果你的样本顺序不同，请在这里修改
human_group <- c(rep("Non-FILM", 9), rep("FILM", 9))
names(human_group) <- colnames(human_expr)

human_anno <- make_group_annotation(colnames(human_expr), human_group)
write.csv(human_anno, file.path(human_outdir, "human_sample_group.csv"), row.names = FALSE)

# immunedeconv supports xCell, MCPcounter, Quantiseq and other methods
# Install if needed:
# install.packages("remotes")
# remotes::install_github("omnideconv/immunedeconv")

if (requireNamespace("immunedeconv", quietly = TRUE)) {
  library(immunedeconv)
  
  # xCell
  human_xcell <- deconvolute(human_expr, method = "xcell")
  write.csv(human_xcell, file.path(human_outdir, "human_xCell_results.csv"), row.names = FALSE)
  
  # MCPcounter
  human_mcp <- deconvolute(human_expr, method = "mcp_counter")
  write.csv(human_mcp, file.path(human_outdir, "human_MCPcounter_results.csv"), row.names = FALSE)
  
  # Quantiseq
  human_quantiseq <- deconvolute(human_expr, method = "quantiseq")
  write.csv(human_quantiseq, file.path(human_outdir, "human_Quantiseq_results.csv"), row.names = FALSE)
  
  # Convert immunedeconv output to matrix for heatmap
  format_deconv_matrix <- function(x) {
    x <- as.data.frame(x)
    rownames(x) <- x$cell_type
    x$cell_type <- NULL
    x <- as.matrix(x)
    storage.mode(x) <- "numeric"
    return(x)
  }
  
  plot_deconv_heatmap(
    format_deconv_matrix(human_xcell),
    human_anno,
    file.path(human_outdir, "human_xCell_heatmap.pdf")
  )
  
  plot_deconv_heatmap(
    format_deconv_matrix(human_mcp),
    human_anno,
    file.path(human_outdir, "human_MCPcounter_heatmap.pdf")
  )
  
  plot_deconv_heatmap(
    format_deconv_matrix(human_quantiseq),
    human_anno,
    file.path(human_outdir, "human_Quantiseq_heatmap.pdf")
  )
  
} else {
  warning("Package immunedeconv is not installed. xCell/MCPcounter/Quantiseq were skipped.")
}

# =========================
# 3. Prepare human expression matrix for ImmunCellAI
# =========================
# If ImmunCellAI was run through its web server or official software,
# use this exported matrix as input and save returned results manually.

human_immunecellai_input <- data.frame(
  GeneSymbol = rownames(human_expr),
  human_expr,
  check.names = FALSE
)

write.csv(
  human_immunecellai_input,
  file.path(human_outdir, "human_expression_for_ImmunCellAI.csv"),
  row.names = FALSE
)

# If you already have ImmunCellAI result file, put it in this folder and read it:
human_immunecellai_result_file <- file.path(human_outdir, "human_ImmunCellAI_results.csv")

if (file.exists(human_immunecellai_result_file)) {
  human_immunecellai <- read.csv(human_immunecellai_result_file, check.names = FALSE)
  write.csv(human_immunecellai, file.path(human_outdir, "human_ImmunCellAI_results_imported.csv"), row.names = FALSE)
}

# =========================
# 4. Mouse ImmunCellAI only
# =========================

mouse_expr_file <- "250508.txt"   # rows = genes, columns = 9 mouse samples
mouse_outdir <- "Mouse_immune_deconvolution_results"
dir.create(mouse_outdir, showWarnings = FALSE, recursive = TRUE)

mouse_expr <- read_expression_matrix(mouse_expr_file)

# 前3个 Non-FILM，后6个 FILM；如果你的样本顺序不同，请在这里修改
mouse_group <- c(rep("Non-FILM", 3), rep("FILM", 6))
names(mouse_group) <- colnames(mouse_expr)

mouse_anno <- make_group_annotation(colnames(mouse_expr), mouse_group)
write.csv(mouse_anno, file.path(mouse_outdir, "mouse_sample_group.csv"), row.names = FALSE)

# Prepare mouse expression matrix for ImmunCellAI
mouse_immunecellai_input <- data.frame(
  GeneSymbol = rownames(mouse_expr),
  mouse_expr,
  check.names = FALSE
)

write.csv(
  mouse_immunecellai_input,
  file.path(mouse_outdir, "mouse_expression_for_ImmunCellAI.csv"),
  row.names = FALSE
)

# If you already have mouse ImmunCellAI result file, put it in this folder and read it:
mouse_immunecellai_result_file <- file.path(mouse_outdir, "mouse_ImmunCellAI_results.csv")

if (file.exists(mouse_immunecellai_result_file)) {
  mouse_immunecellai <- read.csv(mouse_immunecellai_result_file, check.names = FALSE)
  write.csv(mouse_immunecellai, file.path(mouse_outdir, "mouse_ImmunCellAI_results_imported.csv"), row.names = FALSE)
}

# =========================
# 5. Save session information
# =========================

writeLines(
  capture.output(sessionInfo()),
  file.path("sessionInfo_bulk_RNAseq_GSEA_immune_deconvolution.txt")
)

message("Bulk RNA-seq immune deconvolution analysis completed.")