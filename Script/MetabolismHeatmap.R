
#0.Set Work Dirtory----

setwd("D:/Work/Spatial_Project/ST_SM_MazieStem")

#1.Load R Packages----

library(Seurat)
library(tidyr)
library(dplyr)
library(openxlsx)
library(R.utils)
library(ggplot2)
library(SCP)
library(circlize)
library(patchwork)
library(ComplexHeatmap)
library(officer)
library(rvg)
library(devEMF)

#2.Load R data----

#2.1.Read RDS Data----

maize_se <- readRDS(file.path(getwd(),'Git_data','maize_merge_new.rds'))

#2.2.Read Metabolism List----

Ac_neg <- read.xlsx(file.path(getwd(),'Git_data','Amino_acid_neg.xlsx'))
Ac_pos <- read.xlsx(file.path(getwd(),'Git_data','Amino_acid_pos.xlsx'))

Ac_neg_list <- Ac_neg %>% select(mz,Metabolites)
Ac_pos_list <- Ac_pos %>% select(mz,Metabolites)

Ac_neg_list <- unique(Ac_neg_list)
Ac_pos_list <- unique(Ac_pos_list)

features_pos <- paste0(Ac_neg_list$mz)
features_neg <- paste0(Ac_pos_list$mz)
features_all <- c(features_pos, features_neg)

Ac_list <- rbind(Ac_neg_list,Ac_pos_list)

#3.Draw GroupHeatmap-----

#3.1.Merge Assay----

pos_mat <- GetAssayData(maize_se, assay = "Metabolism_pos", slot = "data")
neg_mat <- GetAssayData(maize_se, assay = "Metabolism_neg", slot = "data")

common_cells <- intersect(colnames(pos_mat), colnames(neg_mat))
pos_mat <- pos_mat[, common_cells, drop = FALSE]
neg_mat <- neg_mat[, common_cells, drop = FALSE]

# 如果行名重复，可加前缀
if(length(intersect(rownames(pos_mat), rownames(neg_mat))) > 0){
  rownames(pos_mat) <- paste0("pos_", rownames(pos_mat))
  rownames(neg_mat) <- paste0("neg_", rownames(neg_mat))
}

merged_mat <- rbind(pos_mat, neg_mat)

new_assay <- CreateAssayObject(counts = merged_mat)
maize_se[["Metabolism_all"]] <- new_assay

#3.2.Extract Metabolism Cell----

# 获取 Metabolism_neg 中有数据的细胞
valid_cells <- colnames(maize_se[["Metabolism_all"]])

# 在 Seurat 对象中保留这些细胞
maize_se_sub <- subset(maize_se, cells = valid_cells)

# 传入features参数
features_use <- features_all[features_all %in% rownames(maize_se_sub[["Metabolism_all"]])]

# 确保features_use不为空
length(features_use)

# 确保是 character 类型
features_use <- as.character(features_use)
Ac_list$mz <- as.character(Ac_list$mz)

# 替换 mz 为代谢物名（保留顺序）
features_label <- mz_to_name[features_use]
names(features_label) <- features_use  # 确保 names 对应 rownames

# 创建映射向量
mz_to_name <- setNames(Ac_list$Metabolites, Ac_list$mz)

#3.3.Draw Goup HeatMap----

maize_se_sub@meta.data$Clusters <- maize_se_sub@meta.data$newclusters

# 1. 指定 Cluster 顺序
desired_order <- c("Xyl", "Phl", "Scl","TB/EP", 
                   "LDGT", "EDGT",
                   "Hypo", "Epi", 
                   "Unk1","Unk2", "Unk3")

maize_se_sub@meta.data$Clusters <- factor(maize_se_sub@meta.data$Clusters, levels = desired_order)


# 2. 给每个 Cluster 定义颜色（你可按需求换色）
cluster_colors <- c(
  "Xyl" = "#1f77b4",
  "Phl" = "#ff7f0e",
  "Scl" = "#2ca02c",
  "TB/EP" = "#d62728",
  "LDGT" = "#9467bd",
  "EDGT" = "#8c564b",
  "Hypo" = "#e377c2",
  "Epi" = "#7f7f7f",
  "Unk1" = "#bcbd22",
  "Unk2" = "#17becf",
  "Unk3" = "#aec7e8"
)

variety_colors <- c("B73" = "#66c2a5", "Teo" = "#fc8d62")

# 3. 调用 GroupHeatmap，传入颜色列表
ht <- GroupHeatmap(
  srt = maize_se_sub,
  features = features_use,
  # features_label = features_use,
  group.by = "Clusters",                 # 用你的分组列
  heatmap_palette = "YlOrRd",
  cell_annotation = "variety",         
  cell_annotation_palcolor  = list(variety = variety_colors),    # 这个参数用于默认调色板，可保持或删
  group_palcolor = list(Clusters = cluster_colors),
  heatmap_palcolor = c("#2e77b4", "#f7f1ee", "#940e26"),
  show_row_names = FALSE,
  row_names_side = "left",
  add_dot = F,
  show_column_names = T ,
  add_reticle = FALSE,
  assay = "Metabolism_all",
  slot = "data",
  label_color = "black",
  label_size = 12 ,
  cluster_columns = T ,
  cluster_rows = T,
  column_title = "Cluster",
  row_title = "Metabolite",
  row_title_rot = 90,
  ht_params = list(
    show_row_dend = FALSE,
    show_column_dend = FALSE),  # 传递给 ComplexHeatmap::Heatmap，用于隐藏聚类线条
  exp_cutoff = 1
  
)

print(ht$plot)

heatmap_obj <- ht$matrix_list[[1]]

rownames(heatmap_obj)

rownames(heatmap_obj) <- Ac_list$Metabolites[match(rownames(heatmap_obj), Ac_list$mz)]

rownames(heatmap_obj)[is.na(rownames(heatmap_obj))] <- "Isoleucine"

row_anno_labels <- rownames(heatmap_obj)

ht <- Heatmap(
  heatmap_obj,
  # name = "Expression",
  name = " ", #隐藏图例名称
  col = colorRamp2(c(-4, 0, 4), c("#2e77b4", "#f7f1ee", "#940e26")), 
  rect_gp = gpar(col = "gray80", lwd = 0.5),
  border = "gray50",
  heatmap_legend_param = list(
    title = NULL, #不显示标题
    border = "black",  # 设置图例边框颜色
    at = c(-4, -2,0,2, 4), #图例刻度
    direction = "vertical",#垂直图例
    title_position = "topcenter" ,#标题位置
    legend_height = unit(3, "cm")  # 图例高度
  ),
  # border_gp = gpar(col = "black", lty = 2, lwd = 2) , # 虚线，宽度2
  cluster_rows = TRUE,  # 保持行聚类
  cluster_columns = TRUE,  # 保持列聚类
  show_row_dend = FALSE,  # 隐藏行聚类树
  show_column_dend = FALSE,  # 隐藏列聚类树
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  row_labels = row_anno_labels,
  column_title = "Cluster",
  # gap = unit(2, "mm"),
  jitter = F,
  row_title = "Metabolite"
)
class(ht)

#4.Save Plot----

svg(file.path(getwd(),'Figures','MetabolismHeatmap.svg'), width = 10, height = 7)
draw(ht, heatmap_legend_side = "right")
dev.off()






















