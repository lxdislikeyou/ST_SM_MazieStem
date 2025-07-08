
#0.Set Work Dirtory----

setwd("D:/Work/Spatial_Project/空间转录组点对点/Spatial_project")

#1.Load R Packages----

library(Seurat)
library(tidyr)
library(dplyr)
library(openxlsx)
library(R.utils)
library(ggplot2)
library(SCP)
library(patchwork)
library(cetcolor)
library(officer)
library(rvg)

#2.Load R data----

maize_se <- readRDS("D:/Work/Spatial_Project/空间转录组点对点/空转代点对点DZOE2025040772-b1/RDS/maize_merge_new.rds")

#3.Draw Heat Umap----

DefaultAssay(maize_se) <- "SCT"

B73 <- subset(maize_se,orig.ident == c("B3","B7"))
Teo <- subset(maize_se,orig.ident == c("T1","T3"))

DefaultAssay(B73) <- "SCT"
DefaultAssay(Teo) <- "SCT"

#3.1.UMAP Plot----

cluster_levels <- levels(B73@meta.data$newclusters)  

default_colors <- scales::hue_pal()(length(cluster_levels))

stopifnot(length(default_colors) == length(cluster_levels))

named_colors <- setNames(default_colors, cluster_levels)

B73_plot <- CellDimPlot(B73, group.by = "newclusters", reduction = "UMAP", theme_use = "theme_blank",
                        label = T,
                        pt.alpha = 0.7,
                        label.bg = "black",
                        label.bg.r = 0.1,
                        label_insitu = FALSE,
                        label_repel = FALSE,
                        label_segment_color = "black",
                        lineages_trim = c(0.01, 0.99),
                        streamline_L = 5,
                        add_mark = F,
                        palcolor = named_colors,
                        mark_type = "hull",  # 或 "ellipse"
                        mark_alpha = 0.1,    # 设置透明度
                        mark_linetype = 4,
                        mark_expand = unit(1, "mm"),  # 设置轮廓扩展范围
                        add_density =T ,
                        density_color = "grey80",
                        # density_filled = T,
                        # density_filled_palette = "Greys",
                        pt.size = 1.5)

Teo_plot <- CellDimPlot(Teo, group.by = "newclusters", reduction = "UMAP", theme_use = "theme_blank",
                        label = T,
                        pt.alpha = 0.7,
                        label.bg = "black",
                        label.bg.r = 0.1,
                        label_insitu = FALSE,
                        label_repel = FALSE,
                        label_segment_color = "black",
                        lineages_trim = c(0.01, 0.99),
                        streamline_L = 5,
                        add_mark = F,
                        palcolor = named_colors,
                        mark_type = "hull",  # 或 "ellipse"
                        mark_alpha = 0.1,    # 设置透明度
                        mark_linetype = 4,
                        mark_expand = unit(1, "mm"),  # 设置轮廓扩展范围
                        add_density =T ,
                        density_color = "grey80",
                        # density_filled = T,
                        # density_filled_palette = "Greys",
                        pt.size = 1.5)



BT_plot <- B73_plot|Teo_plot
BT_plot

umap_save <- "./UMAP_plot/"

if (!dir.exists(umap_save)) {
  dir.create(umap_save)
}

ggsave(plot = BT_plot,paste0(umap_save,"BT_umap.png"),dpi = 600,width = 13,height = 5)

#3.1.UMAPplot----

B73_nUmap <- CellDimPlot(B73, group.by = "newclusters", reduction = "UMAP")
Teo_nUmap <- CellDimPlot(Teo, group.by = "newclusters", reduction = "UMAP")

BT_nplot <- B73_nUmap|Teo_nUmap

ggsave(plot = BT_nplot,paste0(umap_save,"BT_numap.png"),dpi = 600,width = 13,height = 5)

# For the first plot (BT_plot)
ggsave(plot = BT_plot, paste0(umap_save, "BT_umap.pdf"), device = "pdf", 
       width = 13, height = 5)
# Or as SVG
ggsave(plot = BT_plot, paste0(umap_save, "BT_umap.svg"), device = "svg", 
       width = 13, height = 5)

# For the second plot (BT_nplot)
ggsave(plot = BT_nplot, paste0(umap_save, "BT_numap.pdf"), device = "pdf", 
       width = 13, height = 5)
# Or as SVG
ggsave(plot = BT_nplot, paste0(umap_save, "BT_numap.svg"), device = "svg", 
       width = 13, height = 5)


ppt <- read_pptx()  

# 添加第一张图
ppt <- ppt %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(
    dml(ggobj = BT_plot),
    location = ph_location_fullsize()
  )

# 添加第二张图
ppt <- ppt %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(
    dml(ggobj = BT_nplot),
    location = ph_location_fullsize()
  )

# 保存合并后的 PPT
print(ppt, target = paste0(umap_save, "BT_UMAP_combined.pptx"))

