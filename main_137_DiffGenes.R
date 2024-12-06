# create: 21 Nov, 2023
# Update: 10 Sep, 2024
# by Feizhen Wu 

setwd("/home/u24211510018/workspace/BI_project/Bulk-RNA-seq/4.Differential_Expression_Genes/DiffGenes_MiR137/")
rm(list=ls())

library(ggplot2)

#Data preparation
{
  DEG=read.table("gene_exp.diff",header = T)
  DEG=DEG[,c(3,10,12)]
  DEG=DEG[is.finite(DEG$log2.fold_change.),]
  DEG=DEG[abs(DEG$log2.fold_change.)>log2(1.5) & DEG$p_value<0.05,]
  names(DEG)=c("genes","foldchange","pvalue")
  DEG$regulation="up"
  DEG$regulation[DEG$foldchange<0]="down"
}

#bar-plot
{
  tab = as.data.frame(table(DEG$regulation))
  tab$Var1=factor(tab$Var1,levels=c("up","down"))
  p=ggplot(tab,aes(x=Var1,y=Freq,label = Freq,fill=Var1))+geom_bar(stat = "identity")
  p=p+geom_text(position = position_dodge(0.9),vjust = 0,size=3)+ylim(0,max(tab$Freq)*1.1)
  p=p+theme_classic(8)+xlab("differential expression")+ylab("Number of genes")
  p=p+ggtitle("Zika-diffgenes")+theme(legend.position = "none")
  p=p+theme(plot.title = element_text(hjust = 0.5))
  p
  ggsave(p,filename = "MiR137_diffgene_number_barplot.pdf",width = 2.2,height = 2.2)
}

#heatmap
{
  library(scales)
  library(pheatmap)
  # dd是所有基因在所有样品中的表达量，没有任何的过滤
  dd=read.table("../Expression_TPM_TPM.xls",header = T)
  DEG1=DEG[order(abs(DEG$foldchange),decreasing = T),]
  
  DEG1=DEG1[1:40,]
  dd1=dd[dd$gene_id %in% DEG1$genes,]
  row.names(dd1)=dd1$gene_id
  dd1$gene_id=NULL
  names(dd1)=c("cWT1","cWT2","cWT3","cKO1","cKO2","cKO3")
  dd2=t(apply(dd1,1,rescale))
  pdf(file = "top_40gene_MiR137.pdf",width = 3,height = 5.8)
  pheatmap(dd2,cutree_rows = 2,cutree_cols = 2,fontsize_row =8)
  dev.off()
}

#enrichment analysis
{
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(cowplot)
  gene <- bitr(DEG$genes, fromType ="SYMBOL",
                  toType =  "ENTREZID",
                  OrgDb = org.Mm.eg.db)
  
  geneList <- bitr(dd$gene_id, fromType ="SYMBOL",
                          toType =  "ENTREZID",
                          OrgDb = org.Mm.eg.db)
  
  ego <- enrichGO(gene          = gene$ENTREZID,
                  universe      = names(geneList$ENTREZID),
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  kk <- enrichKEGG(gene         = gene$ENTREZID,
                   organism     = 'mmu',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05)
  
  p1 <- dotplot(ego, showCategory=5,orderBy = "x") + ggtitle("dotplot for GOBP")
  p2 <- dotplot(kk, showCategory=5,orderBy = "x") + ggtitle("dotplot for KEGG")
  pp=plot_grid(p1, p2, ncol=2)
  pp
  ggsave(pp,filename = "Mir137_DEG_enrichment.pdf",width = 12,height = 3.8)

  write.table(DEG,file = "Mir_137_DEG.xls",sep="\t",quote = F,row.names = F)
}


# predicted target of Mir137 ----------------------------------------------
{
  library(openxlsx)
  
  # 读取 Excel 文件
  predicted <- read.xlsx("./TargetScan7.1__miR-137-3p.predicted_targets.xlsx")
  head(predicted)
  predicted[predicted$Target.gene=="Jdp2",]
  predicted1 = predicted[,1]
  head(predicted1)
  length(predicted1)
  DEG_gene = DEG[,1]
  length(DEG_gene)
  length(DEG$genes)
  head(DEG$genes)
  
  # 交集
  common_genes <- intersect(predicted1, DEG_gene)
  
  # 查看交集
  print(common_genes)
}


# vocanol plot ------------------------------------------------------------
{
  
  #Data preparation
  {
    DEG1=read.table("gene_exp.diff",header = T)
    DEG1=DEG1[,c(3,10,12)]
    DEG1=DEG1[is.finite(DEG1$log2.fold_change.),]
    DEG1=DEG1[abs(DEG1$log2.fold_change.)>log2(2) & DEG1$p_value<0.05,]
    names(DEG1)=c("genes","foldchange","pvalue")
    DEG1$regulation="up"
    DEG1$regulation[DEG1$foldchange<0]="down"
  }
  library(ggplot2)
  library(ggrepel)
  
  ALL_DEG=read.table("gene_exp.diff",header = T)
  ALL_DEG=ALL_DEG[,c(3,10,12)]
  ALL_DEG=ALL_DEG[is.finite(ALL_DEG$log2.fold_change.),]
  names(ALL_DEG)=c("genes","foldchange","pvalue")
  ALL_DEG$significance <- ifelse(ALL_DEG$pvalue < 0.05, "Significant", "Not Significant")
  ALL_DEG = ALL_DEG[ALL_DEG$foldchange != 0 & !is.na(ALL_DEG$foldchange), ]
  
  
  # 假设 DEG 数据框包含了 pvalue 列，你可以根据 pvalue 创建 significance 列
  DEG1$significance <- ifelse(DEG1$pvalue < 0.05, "Significant", "Not Significant")
  
  # 按照 foldchange 排序并取前 20 和后 20 个基因
  DEG_top = DEG1[order(DEG1$foldchange, decreasing = TRUE), ][1:20, -c(4)]
  DEG_bottom = DEG1[order(DEG1$foldchange, decreasing = FALSE), ][1:20,-c(4) ]
  
  # 设置 y 轴的最大值
  max_y <- max(-log(ALL_DEG$pvalue), na.rm = TRUE) + 1
  # max_y <- 100
  pdf(file = "volcanol_MiR137.pdf",width = 7.84,height = 4.76)
  ggplot(ALL_DEG, aes(foldchange, -log(pvalue))) +
    geom_point(aes(color = significance), size = 1) +  # 根据 significance 列着色
    scale_y_continuous(limits = c(-1, max_y)) +  # 设置 y 轴的限制
    scale_color_manual(values = c("black", "red")) +  # 显著性为红色，不显著为黑色
    ggrepel::geom_text_repel(
      data = DEG_top,  # 使用 DEG_top 数据来显示前 20 个基因的标签
      aes(label = genes),  # 使用 genes 列作为标签
      size = 1.5,
      box.padding = 0.35, max.overlaps = Inf  # 避免标签重叠
    ) +
    ggrepel::geom_text_repel(
      data = DEG_bottom,  # 使用 DEG_bottom 数据来显示后 20 个基因的标签
      aes(label = genes),  # 使用 genes 列作为标签
      color = "#619CFF",  # 设置标签颜色
      size = 1.5,
      box.padding = 0.35, max.overlaps = Inf  # 避免标签重叠
    ) +
    labs(x = "Log2FoldChange", y = "-(Log normalized p-value)") +  # 设置坐标轴标签
    geom_vline(xintercept = 0, linetype = "dotted") +  # 添加垂直参考线
    theme_minimal()  # 使用简洁的主题
  dev.off()
  
}
