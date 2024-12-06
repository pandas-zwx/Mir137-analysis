rm(list=ls())
library(Seurat)
library(ggplot2)
setwd("/home/u24211510018/workspace/BI_project/Sc_RNAseq_cKO_Mir137/GSE222858_CellRanger_out_matrix")


{
  # 读取 cWT 样本
  cWT1 <- Read10X(data.dir = "./cWT-1/")
  cWT2 <- Read10X(data.dir = "./cWT-2/")
  
  # 创建 Seurat 对象
  cWT1_seurat <- CreateSeuratObject(counts = cWT1, project = "cWT1",min.cells = 3, min.features = 200)
  cWT2_seurat <- CreateSeuratObject(counts = cWT2, project = "cWT2",min.cells = 3, min.features = 200)
  
  # 读取 cKO 样本
  cKO1 <- Read10X(data.dir = "./cKO-1/")
  cKO2 <- Read10X(data.dir = "./cKO-2/")
  
  # 创建 Seurat 对象
  cKO1_seurat <- CreateSeuratObject(counts = cKO1, project = "cKO1",min.cells = 3, min.features = 200)
  cKO2_seurat <- CreateSeuratObject(counts = cKO2, project = "cKO2",min.cells = 3, min.features = 200)
  
  # 合并两个 cWT 样本
  cWT_combined <- merge(cWT1_seurat, y = cWT2_seurat, 
                        add.cell.ids = c("cWT1", "cWT2"), 
                        project = "cWT")
  
  # 合并两个 cKO 样本
  cKO_combined <- merge(cKO1_seurat, y = cKO2_seurat, 
                        add.cell.ids = c("cKO1", "cKO2"), 
                        project = "cKO")
  
  # 合并所有样本
  mir137_combined <- merge(cWT_combined, y = cKO_combined, 
                           add.cell.ids = c("cWT", "cKO"), 
                           project = "Mir137_All")
}

# QC
{
  Idents(mir137_combined) <- "Total cells"
  mir137_combined[["percent.mt"]] <- PercentageFeatureSet(mir137_combined, pattern = "(?i)^MT-")
  mir137_combined[["percent.ribo"]] <- PercentageFeatureSet(object = mir137_combined, pattern = "Rps|Rpl|Mrpl|Mrps")
  pdf("Quality control1.pdf",width = 10,height = 5)
  VlnPlot(mir137_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4,group.by = NULL)
  dev.off()
  pdf("Quality control2.pdf",width = 12,height = 5)
  plot1 <- FeatureScatter(mir137_combined, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = NULL)
  plot2 <- FeatureScatter(mir137_combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = NULL)
  plot3 <- FeatureScatter(mir137_combined, feature1 = "nCount_RNA", feature2 = "percent.ribo",group.by = NULL)
  plot1 + plot2 + plot3
  dev.off()
}

# filter condition: nCount: > 2000 and up to 20,000, nFeature: > 800 and up to 6,000, and percent.mt < 15%
{
  mir137_combined <- subset(mir137_combined, subset = nFeature_RNA > 800 & nFeature_RNA < 6000 & percent.mt < 15 & nCount_RNA >2000 & nCount_RNA < 20000)
}

# normalize data
{
  mir137_combined <- NormalizeData(mir137_combined, normalization.method = "LogNormalize", scale.factor = 10000)
  head(mir137_combined[["RNA"]]$data.cWT1.cWT)
}

# findvariable genes
{
  mir137_combined <- FindVariableFeatures(mir137_combined, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(mir137_combined), 10)
}

# Scaling the data 均值为0，标准差为1
{
  plan("multicore", workers = 8)
  all.genes <- rownames(mir137_combined)
  mir137_combined <- ScaleData(mir137_combined, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"),verbose = TRUE)
  plan("multicore", workers = 1)
}

# RunPCA
{
  # 之前选择的2k个高变异基因
  # mir137_combined <- RunPCA(mir137_combined, features = VariableFeatures(object = mir137_combined))
  mir137_combined <- RunPCA(mir137_combined, features = VariableFeatures(object = mir137_combined), ndims.print = 1:2,npcs = 100)
  ElbowPlot(object = mir137_combined, ndims = 100)
  # determine PC value
  {
    library(ggplot2)
    
    # 创建 ElbowPlot 并直接修改 Y 轴刻度
    ElbowPlot(object = mir137_combined, ndims = 100) +
      scale_y_continuous(breaks = seq(0, 5, by = 1)) +
      scale_x_continuous(breaks = seq(0, 100, by = 5)) # 修改刻度间隔
  }
  print(mir137_combined[["pca"]], dims = 1:5)
  DimHeatmap(object = mir137_combined,dims = 1:6)
  DimPlot(object = mir137_combined,dims = c(1,2))
  VizDimLoadings(mir137_combined, dims = 1:2, reduction = "pca")
}

# 1. ----------------------------------------------------------------------
if (F) {
  saveRDS(mir137_combined, file = "./mir137_combined.1")
  mir137_combined <- readRDS(file = "./mir137_combined.1")
}

# 2.-------------------------------------------------------------------------

# Cluster the cells
{
  mir137_combined <- FindNeighbors(mir137_combined, dims = 1:20)
  mir137_combined <- FindClusters(mir137_combined, resolution = c(1.0))
  mir137_combined <- RunUMAP(mir137_combined, dims = 1:20)
  pdf("output_plot.pdf")
  DimPlot(mir137_combined, reduction = "umap")
  dev.off()
  DimPlot(mir137_combined, reduction = "umap")# r02 <- DimPlot(mir137_combined, label = T, reduction = "umap", group.by = "RNA_snn_res.0.2")+labs(title="resolution 0.2")
  # r04 <- DimPlot(mir137_combined, label = T, reduction = "umap", group.by = "RNA_snn_res.0.4")+labs(title="resolution 0.4")
  # r06 <- DimPlot(mir137_combined, label = T, reduction = "umap", group.by = "RNA_snn_res.0.6")+labs(title="resolution 0.6")
  # r08 <- DimPlot(mir137_combined, label = T, reduction = "umap", group.by = "RNA_snn_res.0.8")+labs(title="resolution 0.8")
  # CombinePlots(plots=list(r02,r04,r06,r08), ncol = 2)
  # FeaturePlot(mir137_combined, c("Ifit3", "Pik3ip1",'Malat1','Olfml3'))
  FeaturePlot(mir137_combined, c("Ifit3", "Ifit2",'Isg15','Ifit3b','Tmem119'))
}



# Filter doublets and multiplets ----------------------------------------
# 每个细胞表达至少5个该cluster的top_genes,然后还有其他其他细胞的top_genes
# 不能是该细胞的top_genes为什么呢？有可能是管家基因，有可能markergenes低表达，但是其他cluster不表达
{
  mir137_combined <- JoinLayers(mir137_combined)
  cluster_markers <- FindAllMarkers(object = mir137_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,)
  top10_genes <- cluster_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))

  scale_data <- GetAssayData(mir137_combined, slot = "scale.data")
  # scale_data <- GetAssayData(mir137_combined, slot = "count")
  # scale_data <- GetAssayData(mir137_combined, slot = "data")
  
  # 定义函数：提取单个细胞的 top_genes
  get_top_genes <- function(cell_expression, n = 10) {
    ranked_genes <- names(sort(cell_expression, decreasing = TRUE))  # 基因按表达值降序排序
    return(ranked_genes[1:n])  # 返回前 n 个基因
  }
  
  # 使用 lapply() 生成命名列表
  cell_top_genes <- lapply(colnames(scale_data), function(cell) {
    get_top_genes(scale_data[, cell])
  })
  names(cell_top_genes) <- colnames(scale_data)  # 设置名称为细胞名
  
  # 检查一个细胞是否包含 5 个其他簇的 top_genes
  is_multiplet <- function(cell_genes, top10_genes, current_cluster) {
    # 过滤出其他簇的特异性基因
    other_clusters_genes <- top10_genes %>%
      filter(cluster != current_cluster) %>%
      pull(gene)
    
    # 计算当前细胞的 top_genes 与其他簇 top_genes 的交集数量
    overlap <- sum(cell_genes %in% other_clusters_genes)
    return(overlap >= 0)  # 返回是否满足条件
  }
  
  # 对每个细胞进行判断
  multiplet_cells <- sapply(colnames(scale_data), function(cell) {
    current_cluster <- Idents(mir137_combined)[cell]  # 获取该细胞所属的簇
    cell_genes <- cell_top_genes[[cell]]  # 获取该细胞的 top_genes
    is_multiplet(cell_genes, top10_genes, current_cluster)
  })
  
  # 将每个细胞是否为 multiplet 结果加入 Seurat 对象元数据
  mir137_combined <- AddMetaData(mir137_combined, metadata = multiplet_cells, col.name = "is_multiplet")
  
  # 计算每个簇中 multiplet 的比例
  multiplet_rate_per_cluster <- mir137_combined@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(multiplet_rate = mean(is_multiplet))
  
  # 标记需要移除的簇
  clusters_to_remove <- multiplet_rate_per_cluster %>%
    filter(multiplet_rate > 0.3) %>%
    pull(seurat_clusters)
  
  # 移除双细胞比例高的簇
  mir137_combined <- subset(mir137_combined, idents = setdiff(Idents(mir137_combined), clusters_to_remove))
  
}

{
  scale_data <- GetAssayData(mir137_combined, slot = "data")
  mir137_combined <- JoinLayers(mir137_combined)
  cluster_markers <- FindAllMarkers(object = mir137_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,)
  # 获取 FindAllMarkers() 结果的 top10 marker genes
  top10_genes <- cluster_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))
  
  # 定义函数：检查某细胞是否为 multiplet
  is_multiplet <- function(cell_name, top10_genes, current_cluster, scale_data) {
    # 获取当前簇的 marker genes
    current_cluster_genes <- top10_genes %>%
      filter(cluster == current_cluster) %>%
      pull(gene)
    
    # 获取其他簇的 marker genes
    other_clusters_genes <- top10_genes %>%
      filter(cluster != current_cluster) %>%
      pull(gene)
    
    # 获取当前细胞的表达值
    cell_expression <- scale_data[, cell_name]
    
    # 检查是否表达了当前簇的至少 5 个 marker genes
    expresses_current_genes <- sum(cell_expression[current_cluster_genes] > 0) >= 5
    
    # 检查是否表达了其他簇的 marker genes
    expresses_other_genes <- sum(cell_expression[other_clusters_genes] > 0) >= 5
    
    # 如果满足两个条件，则标记为 multiplet
    return(expresses_current_genes && expresses_other_genes)
  }
  
  # 对每个细胞进行判断
  multiplet_cells <- sapply(colnames(scale_data), function(cell_name) {
    current_cluster <- as.character(Idents(mir137_combined)[cell_name])  # 获取细胞所属簇
    is_multiplet(cell_name, top10_genes, current_cluster, scale_data)
  })
  
  # 将结果添加到 Seurat 对象的元数据
  mir137_combined <- AddMetaData(mir137_combined, metadata = multiplet_cells, col.name = "is_multiplet")
  
  # 计算每个簇中 multiplet 的比例
  multiplet_rate_per_cluster <- mir137_combined@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(multiplet_rate = mean(is_multiplet))
  
  # 标记需要移除的簇
  clusters_to_remove <- multiplet_rate_per_cluster %>%
    filter(multiplet_rate > 0.3) %>%
    pull(seurat_clusters)
  print(clusters_to_remove)
  
  # 移除双细胞比例高的簇
  mir137_combined <- subset(mir137_combined, idents = setdiff(Idents(mir137_combined), clusters_to_remove)) 
}

# OK:0/1/
# no_OK：6/10/
mir137_combined <- JoinLayers(mir137_combined)
plan("multicore", workers = 4)
cluster_markers <- FindAllMarkers(object = mir137_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,)
plan("multicore", workers = 1)
# 获取 FindAllMarkers() 结果的 top10 marker genes
library(dplyr)
top10_genes <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(top10_genes, file = "top10_genes.csv", row.names = FALSE)
# 选择所需列
selected_genes <- top10_genes %>%
  select(p_val, avg_log2FC, cluster, gene)
write.csv(selected_genes, file = "selected_genes.csv", row.names = FALSE)


# filter 6/10 -------------------------------------------------------------

mir137_combined <- subset(mir137_combined, idents = c(0:28)[-c(7,11)])
pdf("output_plot.pdf")
DimPlot(mir137_combined, label = T, reduction = "umap")
dev.off()


mir137_combined <- NormalizeData(mir137_combined, normalization.method = "LogNormalize", scale.factor = 10000)
mir137_combined <- FindVariableFeatures(mir137_combined, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mir137_combined)
plan("multicore", workers = 8)
mir137_combined <- ScaleData(mir137_combined, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"))
plan("multicore", workers = 1)

# 2. ----------------------------------------------------------------------
if (F) {
  saveRDS(mir137_combined, file = "./mir137_combined.2")
  mir137_combined <- readRDS(file = "./mir137_combined.2")
}

# -------------------------------------------------------------------------


mir137_combined <- RunPCA (mir137_combined, features = VariableFeatures(object = mir137_combined))
mir137_combined <- FindNeighbors(mir137_combined, dims = 1:20)
mir137_combined <- FindClusters(mir137_combined, resolution = c(0.4))
mir137_combined <- RunUMAP(mir137_combined, dims = 1:20)

pdf("output_plot1.pdf")
DimPlot(mir137_combined, label = T, reduction = "umap")
dev.off()

genes <- c(
  "Cx3cr1", "Hexb", "P2ry12", "Tmem119", "Sparc", "Cst3", "Olfml3", "Fcrls", 
  "Malat1", "Slc1a3", "Dleu2", "Fchsd2", "Smad7", "Pik3ip1", 
  "Ccl4", "Ccl3", "B4galt1", "Fkbp5", "Ifit2", "Ifit3", "Isg15", "Usp18", "Phf11d",
  "Mrc1", "Pf4", "F13a1", "Lyve1", "Dab2", "S100a9", "S100a8", "Retnlg", 
  "Ifitm1", "Wfdc17", "Chil3", "Plac8", "S100a4", "Napsa", "Smpdl3a", 
  "H2-Aa", "H2-Eb1", "H2-Ab1", "Ciita", "Axl", "Nkg7", "Ms4a4b", 
  "Gzmb", "Klrk1", "Prf1"
)
pdf("dot_plot.pdf", width = 15, height = 8)  # 设置宽度为 10 英寸，高度为 8 英寸
DotPlot(mir137_combined, features = genes) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
# 排除：15 10 8 9 5 3
# 留下：0/1/2/4/12/8

# filter2 -----------------------------------------------------------------


# 子集保留指定的 cluster
mir137_combined <- subset(mir137_combined, idents = c(0, 1, 2, 4, 12, 8))
pdf("output_plot.pdf")
DimPlot(mir137_combined, label = T, reduction = "umap")
dev.off()


mir137_combined <- NormalizeData(mir137_combined, normalization.method = "LogNormalize", scale.factor = 10000)
mir137_combined <- FindVariableFeatures(mir137_combined, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mir137_combined)
plan("multicore", workers = 8)
mir137_combined <- ScaleData(mir137_combined, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"))
plan("multicore", workers = 1)
if (F) {
  saveRDS(mir137_combined, file = "./mir137_combined.3")
  mir137_combined <- readRDS(file = "./mir137_combined.3")
}

# -------------------------------------------------------------------------
mir137_combined <- RunPCA (mir137_combined, features = VariableFeatures(object = mir137_combined))
mir137_combined <- FindNeighbors(mir137_combined, dims = 1:20)
mir137_combined <- FindClusters(mir137_combined, resolution = c(0.6))
mir137_combined <- RunUMAP(mir137_combined, dims = 1:20)

pdf("output_plot2.pdf")
DimPlot(mir137_combined, label = T, reduction = "umap")
dev.off()

genes <- c(
  "Cx3cr1", "Hexb", "P2ry12", "Tmem119", "Sparc", "Cst3", "Olfml3", "Fcrls", 
  "Malat1", "Slc1a3", "Dleu2", "Fchsd2", "Smad7", "Pik3ip1", 
  "Ccl4", "Ccl3", "B4galt1", "Fkbp5", "Ifit2", "Ifit3", "Isg15", "Usp18", "Phf11d",
  "Mrc1", "Pf4", "F13a1", "Lyve1", "Dab2", "S100a9", "S100a8", "Retnlg", 
  "Ifitm1", "Wfdc17", "Chil3", "Plac8", "S100a4", "Napsa", "Smpdl3a", 
  "H2-Aa", "H2-Eb1", "H2-Ab1", "Ciita", "Axl", "Nkg7", "Ms4a4b", 
  "Gzmb", "Klrk1", "Prf1"
)
pdf("dot_plot1.pdf", width = 15, height = 8)  # 设置宽度为 10 英寸，高度为 8 英寸
DotPlot(mir137_combined, features = genes) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
# 0.6

# final filter ------------------------------------------------------------
# filter cluster 7/9
mir137_combined <- subset(mir137_combined, idents = c(0, 1, 2, 3, 4, 5, 6, 8, 10))
pdf("output_plot3.pdf")
DimPlot(mir137_combined, label = T, reduction = "umap")
dev.off()
# 将所有剩余 clusters 的标识符改为 "Microglia"
Idents(mir137_combined) <- "Microglia"

pdf("cWT VS cKO.pdf")
DimPlot(mir137_combined, group.by = "cGroup", reduction = "umap", label = TRUE) +
  labs(title = "UMAP Plot with cWT and cKO Labels")
dev.off()

# 创建一个新的分组列
mir137_combined$cGroup <- ifelse(mir137_combined$orig.ident %in% c("cWT1", "cWT2"), "cWT", "cKO")
pdf("Jdp2_feature_plot3.pdf")
FeaturePlot(mir137_combined, features = "Jdp2")
dev.off()
pdf("Jdp2_Vlnplot3.pdf")
VlnPlot(mir137_combined, features = "Jdp2", group.by = "cGroup", pt.size = 0.5)
dev.off()

# final RDS ---------------------------------------------------------------
if (F) {
  saveRDS(mir137_combined, file = "./mir137_combined_final")
  mir137_combined <- readRDS(file = "./mir137_combined_final")
}

# count cell per cluster
{
  table(mir137_combined$orig.ident)
  table(mir137_combined$cGroup)

}
# find DEG between MG(cWT VS cKO)-----------------------------------------------------
Idents(mir137_combined) <- "cGroup"
pdf("cKO_VS_cWT_MG.pdf")
DimPlot(mir137_combined, reduction = "umap",group.by = "cGroup")
dev.off()

markers <- FindMarkers(mir137_combined, ident.1 = "cKO", ident.2 = "cWT", test.use = "wilcox")
write.csv(markers, file = "FindMarkers_results.csv", row.names = TRUE)
# 筛选显著差异基因：调整后 P 值 < 0.05，且 log2 倍数变化 > 0.5
significant_markers <- markers %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > log2(1.5))

# 查看显著差异基因
head(significant_markers)
write.csv(significant_markers, "cKO_vs_cWT_differential_genes.csv")



# volcano plot ------------------------------------------------------------
{library(ggplot2)
  library(dplyr)
  
  # 添加显著性标记
  markers$significance <- ifelse(markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) > log2(1.5), "Significant", "Not Significant")
  
  # 获取显著上调的前五个基因
  top_upregulated <- markers %>%
    filter(significance == "Significant" & avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC)) %>%
    head(5)
  
  # 获取显著下调的前五个基因
  top_downregulated <- markers %>%
    filter(significance == "Significant" & avg_log2FC < 0) %>%
    arrange(avg_log2FC) %>%
    head(5)
  
  # 合并上调和下调的前五个显著基因
  top_genes <- rbind(top_upregulated, top_downregulated)
  
  # 获取前五个显著基因的名称和 p 值
  top_genes_info <- data.frame(
    gene_name = rownames(top_genes),
    p_value = top_genes$p_val_adj,
    avg_log2FC = top_genes$avg_log2FC
  )
  
  # 保存火山图
  pdf("cKO VS cWT volcano with top genes.pdf", width = 15, height = 10)
  ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
    theme_minimal() +
    # 在火山图中标记显著的前五个上调和下调基因
    geom_text(data = top_genes_info, aes(x = avg_log2FC, y = -log10(p_value), label = gene_name), 
              vjust = 1, color = "black", size = 3)
  dev.off()
}
# library(ggplot2)
# 
# # 添加显著性标记
# markers$significance <- ifelse(markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) > log2(1.5), "Significant", "Not Significant")
# 
# pdf("cKO VS cWT volcano.pdf",width = 15,height = 10)
# # 火山图
# ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
#   geom_point(alpha = 0.8) +
#   scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
#   labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
#   theme_minimal()
# dev.off()

pdf("cWT VS cKO heatmeap.pdf")
top_genes <- rownames(head(significant_markers[order(significant_markers$p_val_adj), ], 40))


DoHeatmap(mir137_combined, features = top_genes, group.by = "cGroup") +
  labs(title = "Top Differentially Expressed Genes")
dev.off()

pdf("cWT VS cKO Jdp2 heatmeap.pdf")
DoHeatmap(mir137_combined, features = "Jdp2", group.by = "cGroup") +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +  # 设置渐变颜色
  labs(title = "Heatmap of Jdp2 Expression", x = "Samples", y = "Genes") +
  theme(axis.text.y = element_text(size = 10))  # 调整基因名的字体大小
dev.off()


# bulk_RNA-seq/Sc/predicted -----------------------------------------------
library(openxlsx)


Bulk_DEG = read.table(file = "/home/u24211510018/workspace/BI_project/Bulk-RNA-seq/4.Differential_Expression_Genes/DiffGenes_MiR137/Mir_137_DEG.xls",sep="\t",header = T)
sc_DEG = read.table(file = "cKO_vs_cWT_differential_genes.csv", sep = ",", header = TRUE, quote = "\"")
predicted <- read.xlsx("/home/u24211510018/workspace/BI_project/Bulk-RNA-seq/4.Differential_Expression_Genes/DiffGenes_MiR137/TargetScan7.1__miR-137-3p.predicted_targets.xlsx")

Bulk_DEG1 = Bulk_DEG[abs(Bulk_DEG$foldchange)>log2(2),c(1)]
sc_DEG1 = sc_DEG[abs(sc_DEG$avg_log2FC)>log2(2),1]
predicted1 = predicted[,1]

# 求三个向量的交集
common_genes <- Reduce(intersect, list(Bulk_DEG1, sc_DEG1, predicted1))

# 查看结果
print(common_genes)

Bulk_DEG2 = Bulk_DEG[Bulk_DEG$foldchange>log2(2),c(1)]
sc_DEG2 = sc_DEG[sc_DEG$avg_log2FC>log2(2),1]
predicted2 = predicted[,1]

# 求三个向量的交集
common_genes2 <- Reduce(intersect, list(Bulk_DEG2, sc_DEG2, predicted2))

# 查看结果
print(common_genes2)





# venn --------------------------------------------------------------------

library("ggvenn")

# 将三个基因列表准备为命名列表
gene_lists <- list(
  Bulk_RNA_seq = Bulk_DEG2,
  scRNA = sc_DEG2,
  Predicted_Mir137_Targets = predicted2
)
# 绘制 Venn 图
pdf("Up_DEG_Venn.pdf")
ggvenn(gene_lists,show_percentage = F)
dev.off()



