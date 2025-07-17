
#0.Set Work Dirtory----

# setwd("D:/Work/Spatial_Project/空间转录组点对点/Spatial_project")

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

B73 <- subset(maize_se,orig.ident %in% c("B3","B7"))
Teo <- subset(maize_se,orig.ident %in% c("T3","T1"))

B73 <- PrepSCTFindMarkers(B73)
Teo <- PrepSCTFindMarkers(Teo)

Marker_gene <- function(ST_Data,threshold = 1){
  Markers_Gene <- FindAllMarkers(ST_Data,
                                 assay = "SCT",
                                 logfc.threshold = threshold,
                                 min.diff.pct = 0.2,
                                 only.pos = TRUE)
  Markers_Gene <- as.data.frame(Markers_Gene)
  return(Markers_Gene)
}

B73_Marker <- Marker_gene(B73,1)
Teo_Marker <- Marker_gene(Teo,1)

#3.2.Marker Gene data Save-----

write.xlsx(B73_Marker,file.path(getwd(),"Git_data","B73_Markers_List.xlsx"))

write.xlsx(Teo_Marker,file.path(getwd(),"Git_data","Teo_Markers_List.xlsx"))

#4.Extract Markers Gene Expression-----

Extract_Marker_GeneExpression <- function(Markers_List,ST_Data){
  Markers_Gene <- read.xlsx(file.path(getwd(),"Git_data",paste0(Markers_List,".xlsx")))
  
  Markers_Gene_List <- Markers_Gene$gene
  
  Markers_Gene_Ex <- AverageExpression(ST_Data,
                                       assays = 'SCT',
                                       features = Markers_Gene_List,
                                       group.by = "newclusters")
  
  Markers_Gene_Data <- data.frame(Markers_Gene_Ex$SCT)
  return(Markers_Gene_Data)
}

B73_Marker_Gene_Data <- Extract_Marker_GeneExpression("B73_Markers_List",B73)
# colnames(B73_Marker_Gene_Data) <- paste0("B73_",colnames(B73_Marker_Gene_Data))
Teo_Marker_Gene_Data <- Extract_Marker_GeneExpression("Teo_Markers_List",Teo)
# colnames(Teo_Marker_Gene_Data) <- paste0("Teo_",colnames(Teo_Marker_Gene_Data))

B73_Marker_Gene_Data <- B73_Marker_Gene_Data %>% select("Xyl","Phl","Scl","TB.EP","LDGT","EDGT","Hypo","Epi","Unk1","Unk2","Unk3")

Teo_Marker_Gene_Data <- Teo_Marker_Gene_Data %>% select("Xyl","Phl","Scl","TB.EP","LDGT","EDGT","Hypo","Epi","Unk1","Unk2","Unk3")

#5.Draw Markers Gene HeatMap-----

n_cols <- ncol(Teo_Marker_Gene_Data)
col_gaps <- seq(1, n_cols-1)

B73_Heatmap <- pheatmap(as.matrix(B73_Marker_Gene_Data), 
         color = colorRampPalette(c("#2e77b4", "#f7f1ee", "#940e26"))(50),  # 定义颜色梯度
         cluster_rows = F,  # 对行（GeneID）进行聚类
         cluster_cols = F,  # 对列（Cluster）进行聚类
         scale = "row",  # 对行进行标准化
         main = "B73 Markers Gene Heatmap",
         treeheight_row = 0,
         fontsize_col = 10,
         # angle_col = 45,
         show_rownames = FALSE,
         gaps_col = col_gaps,      # 指定分隔位置
         border_color = "black",  # 分隔线颜色
         silent = TRUE   
         )  # 隐藏行名称（Y轴标签）

 

Teo_Heatmap <- pheatmap(as.matrix(Teo_Marker_Gene_Data), 
         color = colorRampPalette(c("#2e77b4", "#f7f1ee", "#940e26"))(50),  # 定义颜色梯度
         cluster_rows = F,  # 对行（GeneID）进行聚类
         cluster_cols = F,  # 对列（Cluster）进行聚类
         scale = "row",  # 对行进行标准化
         main = "Teo Markers Gene Heatmap",
         use_raster = FALSE,
         treeheight_row = 0,
         fontsize_col = 10,
         # angle_col = 45,
         show_rownames = FALSE,
         gaps_col = col_gaps,      # 指定分隔位置
         border_color = "black",  # 分隔线颜色
         silent = TRUE             # 关闭提示消息 # 分隔线颜色  # 隐藏行名称（Y轴标签）
)
#6.Save plot----

MarkerHeatmap_save <- "Figures"

if (!dir.exists(MarkerHeatmap_save)) {
  dir.create(file.path(getwd(),MarkerHeatmap_save))
}

# 保存为 SVG 矢量图
svg(paste0(MarkerHeatmap_save,"/B73_Heatmap.svg"), width = 10, height = 8)  
B73_Heatmap
dev.off()
svg(paste0(MarkerHeatmap_save,"/Teo_Heatmap.svg"), width = 10, height = 8)  
Teo_Heatmap
dev.off()
# 保存为 PNG
ggsave(plot=heatmap_mrkGene, file.path(getwd(),MarkerHeatmap_save,"MarkerHeatmap_plot.png"), width = 4, height = 12)  


