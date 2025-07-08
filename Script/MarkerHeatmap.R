
#0.Set Work Dirtory----

setwd("D:/Work/Spatial_Project/空间转录组点对点/Spatial_project")

#1.Load R Packages----

library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(scRNAtoolVis)
library(openxlsx)
library(pheatmap)

#2.Load R data----

maize_se <- readRDS("D:/Work/Spatial_Project/空间转录组点对点/空转代点对点DZOE2025040772-b1/RDS/maize_merge_new.rds")

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

Markers_Data_Path <- "D:/Work/Spatial_Project/ST_SM_MazieStem/Git_data/"

write.xlsx(Markers_Data,paste0(Markers_Data_Path,"Markers_List.xlsx"))

#4.Extract Markers Gene Expression-----

Markers_Gene <- read.xlsx(paste0(Markers_Data_Path,"Markers_List.xlsx"))

Markers_Gene_List <- Markers_Gene$gene

Markers_Gene_Ex <- AverageExpression(maize_se,
                                     assays = 'SCT',
                                     features = Markers_Gene_List,
                                     group.by = "newclusters")

Markers_Gene_Data <- data.frame(Markers_Gene_Ex$SCT)

#5.Draw Markers Gene HeatMap-----

colnames(Markers_Gene_Data)
Markers_Gene_Data <- Markers_Gene_Data %>% select("Xyl","Phl","Scl","TB.EP","LDGT","EDGT","Hypo","Epi","Unk1","Unk2","Unk3")

p1 <- pheatmap(as.matrix(Markers_Gene_Data), 
               color = colorRampPalette(c("#2e77b4", "#f7f1ee", "#940e26"))(100),  # 定义颜色梯度
               cluster_rows = F,  # 对行（GeneID）进行聚类
               cluster_cols = F,  # 对列（Cluster）进行聚类
               scale = "row",  # 对行进行标准化
               main = "Markers Gene Heatmap",
               treeheight_row = 0,
               fontsize_col = 10,
               angle_col = 45,
               show_rownames = FALSE)  # 隐藏行名称（Y轴标签）

#6.Save plot----

MarkerHeatmap_save <- "D:/Work/Spatial_Project/空间转录组点对点/Spatial_project/MarkerHeatmap/"

if(!dir.exists(MarkerHeatmap_save)){
  dir.create(MarkerHeatmap_save)
}

# 保存为 SVG 矢量图
svg(paste0(MarkerHeatmap_save,"MarkerHeatmap_plot.svg"), width = 10, height = 8)  
p1
dev.off()






