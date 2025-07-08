
#0.Set Work Dirtory----

# setwd("D:/Work/Spatial_Project/空间转录组点对点/Spatial_project")
getwd()

#1.Load R Packages----

library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(scRNAtoolVis)
library(openxlsx)
library(pheatmap)

#2.Load R data----

# maize_se <- readRDS("D:/Work/Spatial_Project/空间转录组点对点/空转代点对点DZOE2025040772-b1/RDS/maize_merge_new.rds")
maize_se <- readRDS(file.path(getwd(),'Git_data','maize_merge_new.rds'))
#3.Find Markers Gene-----

#3.1.Find All Marker-----

maize_se <- PrepSCTFindMarkers(maize_se)

Markers_Gene <- FindAllMarkers(maize_se,
               assay = "SCT",
               logfc.threshold = 0.25,
               min.diff.pct = 0.2,
               only.pos = TRUE)

Markers_Data <- as.data.frame(Markers_Gene)

#3.2.Marker Gene data Save-----

write.xlsx(Markers_Data,file.path(getwd(),"Git_data","Markers_List.xlsx"))

#4.Extract Markers Gene Expression-----

Markers_Gene <- read.xlsx(file.path(getwd(),"Git_data","Markers_List.xlsx"))

Markers_Gene_List <- Markers_Gene$gene

Markers_Gene_Ex <- AverageExpression(maize_se,
                                     assays = 'SCT',
                                     features = Markers_Gene_List,
                                     group.by = "newclusters")

Markers_Gene_Data <- data.frame(Markers_Gene_Ex$SCT)

#5.Draw Markers Gene HeatMap-----

colnames(Markers_Gene_Data)
Markers_Gene_Data <- Markers_Gene_Data %>% select("Xyl","Phl","Scl","TB.EP","LDGT","EDGT","Hypo","Epi","Unk1","Unk2","Unk3")

heatmap_mrkGene <- pheatmap(as.matrix(Markers_Gene_Data), 
               color = colorRampPalette(c("navy",'#1F78B4', "#FB9A99" ,'#E31A1C'))(50),  # 定义颜色梯度
               # color = colorRampPalette(c("navy", "purple4" ,"#FFFF99","#E31A1C"))(50),  # 定义颜色梯度
               cluster_rows = F,  # 对行（GeneID）进行聚类
               cluster_cols = F,  # 对列（Cluster）进行聚类
               scale = "row",  # 对行进行标准化
               main = "Markers Gene Heatmap",
               treeheight_row = 0,
               fontsize_col = 10,
               angle_col = 45,
               show_rownames = FALSE)  # 隐藏行名称（Y轴标签）

#6.Save plot----

MarkerHeatmap_save <- "Figures"

if (!dir.exists(MarkerHeatmap_save)) {
  dir.create(file.path(getwd(),MarkerHeatmap_save))
}
# 保存为 SVG 矢量图
ggsave(plot=heatmap_mrkGene, file.path(getwd(),MarkerHeatmap_save,"MarkerHeatmap_plot.png"), width = 4, height = 12)  
# heatmap_mrkGene
# dev.off()


