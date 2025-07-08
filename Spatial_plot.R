
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

B3 <- subset(maize_se,orig.ident == c("B3"))
B7 <- subset(maize_se,orig.ident == c("B7"))

T1 <- subset(maize_se,orig.ident == c("T1"))
T3 <- subset(maize_se,orig.ident == c("T3"))


DefaultAssay(B73) <- "SCT"
DefaultAssay(Teo) <- "SCT"

DefaultAssay(B3) <- "SCT"
DefaultAssay(B7) <- "SCT"

DefaultAssay(T1) <- "SCT"
DefaultAssay(T3) <- "SCT"

#4.Feature Plot----

#4.1.Set Colors-----

cluster_levels <- levels(B73@meta.data$newclusters)  

default_colors <- scales::hue_pal()(length(cluster_levels))

stopifnot(length(default_colors) == length(cluster_levels))

named_colors <- setNames(default_colors, cluster_levels)

#4.2.Spatial Plot----

B3_plot <- SpatialPlot(B3,
                        group.by = "newclusters",
                        cols = named_colors,
                        slot = "SCT",
                        # label = T,
                        # label.box = T,
                        # repel = T,
                        # alpha = c(0.5, 0.2),
                        image.alpha = 0.6,
                        # combine = F,
                        pt.size.factor = 2,
                        stroke = NA)

B7_plot <- SpatialPlot(B7,
                        group.by = "newclusters",
                        cols = named_colors,
                        slot = "SCT",
                        # label = T,
                        # label.box = T,
                        # repel = T,
                        # alpha = c(0.5, 0.2),
                        image.alpha = 0.6,
                        # combine = F,
                        pt.size.factor = 1.8,
                        stroke = NA)


T1_plot <- SpatialPlot(T1,
                       group.by = "newclusters",
                       cols = named_colors,
                       slot = "SCT",
                       # label = T,
                       # label.box = T,
                       # repel = T,
                       # alpha = c(0.5, 0.2),
                       image.alpha = 0.6,
                       combine = F,
                       pt.size.factor = 2,
                       stroke = NA)

T3_plot <- SpatialPlot(T3,
                       group.by = "newclusters",
                       cols = named_colors,
                       slot = "SCT",
                       # label = T,
                       # label.box = T,
                       # repel = T,
                       # alpha = c(0.5, 0.2),
                       image.alpha = 0.6,
                       combine = F,
                       pt.size.factor = 2,
                       stroke = NA)

B73_plot <- B3_plot|B7_plot
Teo_plot <- B3_plot|B7_plot

Spatialplot_save <- "./Spatial_plot/"

if (!dir.exists(Spatialplot_save)) {
  dir.create(Spatialplot_save)
}

#5.Save Plot----

#5.1.SVG ----
  
ggsave(plot = B73_plot, paste0(Spatialplot_save, "B73_SpatialPlot.svg"), device = "svg", 
       width = 13, height = 5)
ggsave(plot = Teo_plot, paste0(Spatialplot_save, "Teo_SpatialPlot.svg"), device = "svg", 
       width = 13, height = 5)

#5.2.PPT----

ppt <- read_pptx()  

# 添加第一张图
ppt <- ppt %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(
    dml(ggobj = B73_plot),
    location = ph_location_fullsize()
  )

# 添加第二张图
ppt <- ppt %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(
    dml(ggobj = Teo_plot),
    location = ph_location_fullsize()
  )

# 保存合并后的 PPT
print(ppt, target = paste0(Spatialplot_save, "BT_SpatialPlot_combined.pptx"))













