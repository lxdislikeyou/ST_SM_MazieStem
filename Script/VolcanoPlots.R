#1.Load R Packages----

library(Seurat)
library(tidyr)
library(dplyr)
library(openxlsx)
library(R.utils)
library(ggplot2)
library(SCP)
library(officer)
library(rvg)
library(patchwork)
library(svglite)

#2.Load R data----

maize_se <- readRDS(file.path(getwd(),'Git_data','maize_merge_new.rds'))

#3.Different expression analysis----

#3.1.Set DefaultAssay----
maize_se <- PrepSCTFindMarkers(maize_se) 
DefaultAssay(maize_se) <- "SCT"        

#3.2.Create new Group----
maize_se$cluster_variety <- paste(maize_se$newclusters, maize_se$variety, sep = "_")

#3.3.To Loop Create VolcanoPlot----
Cluster_names <- names(table(maize_se$cluster_variety))

Cluster_name <- data.frame(Cluster_B73 = Cluster_names[grep("B73",Cluster_names)],
                           Cluster_Teo = Cluster_names[grep("TEO",Cluster_names)])


#3.4.Create ppt----
ppt <- read_pptx()

for (i in 1:nrow(Cluster_name)) {
  
  maize_sub <- RunDEtest(
    srt = maize_se,
    group_by = "cluster_variety",
    group1 = Cluster_name$Cluster_Teo[i],
    group2 = Cluster_name$Cluster_B73[i],
    fc.threshold = 1,
    min.cells.group = 3, 
    BPPARAM = BiocParallel::SerialParam(),
    only.pos = F
  )
  
  maize_sub@tools$DEtest_custom$AllMarkers_wilcox$group1 <- paste0(Cluster_name$Cluster_Teo[i],"_vs_",Cluster_name$Cluster_B73[i])
  maize_sub@tools$DEtest_custom$AllMarkers_wilcox$group1 <- as.factor(maize_sub@tools$DEtest_custom$AllMarkers_wilcox$group1)
  
  plot <- VolcanoPlot(
    srt = maize_sub,                       
    DE_threshold = "avg_log2FC > 1 & p_val_adj < 0.05")
  
  ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")
  ppt <- ph_with(ppt, value = dml(ggobj = plot), location = ph_location_fullsize())
  
}

#3.5.Save PPT Plot----

VolcanoPlotsave <- "Figures"

if (!dir.exists(VolcanoPlotsave)) {
  dir.create(file.path(getwd(),VolcanoPlotsave))
}

print(ppt, target = paste0(VolcanoPlotsave,"/","VolcanoPlots_DEcomparison.pptx"))

#3.6.拼接全部图

volcano_list <- list()  # 用于存储所有火山图

for (i in 1:nrow(Cluster_name)) {
  
  # 差异分析
  maize_sub <- RunDEtest(
    srt = maize_se,
    group_by = "cluster_variety",
    group1 = Cluster_name$Cluster_Teo[i],
    group2 = Cluster_name$Cluster_B73[i],
    fc.threshold = 1,
    min.cells.group = 3, 
    BPPARAM = BiocParallel::SerialParam(),
    only.pos = FALSE
  )
  
  comp_name <- paste0(Cluster_name$Cluster_Teo[i], "_vs_", Cluster_name$Cluster_B73[i])
  maize_sub@tools$DEtest_custom$AllMarkers_wilcox$group1 <- comp_name
  maize_sub@tools$DEtest_custom$AllMarkers_wilcox$group1 <- as.factor(comp_name)
  
  # 火山图
  plot <- VolcanoPlot(
    srt = maize_sub,
    DE_threshold = "avg_log2FC > 1 & p_val_adj < 0.05"
  )
  
  volcano_list[[i]] <- plot
}

# 拼接图形（按行排列，每行3个）
combined_plot <- wrap_plots(volcano_list, ncol = 3)

# 保存为 PDF（矢量图）
ggsave(paste0(VolcanoPlotsave,"/","All_VolcanoPlots.pdf"), combined_plot, width = 16, height = ceiling(length(volcano_list)/3) * 5)

# 或者保存为 SVG（也是矢量图）
ggsave(paste0(VolcanoPlotsave,"/","All_VolcanoPlots.svg"), combined_plot, width = 16, height = ceiling(length(volcano_list)/3) * 5)


#4.Extract Different Gene----

wb <- createWorkbook()

#4.1.提取结果表----

for (i in 1:nrow(Cluster_name)) {
  
  maize_sub <- RunDEtest(
    srt = maize_se,
    group_by = "cluster_variety",
    group1 = Cluster_name$Cluster_Teo[i],
    group2 = Cluster_name$Cluster_B73[i],
    fc.threshold = 1,
    min.cells.group = 3, 
    BPPARAM = BiocParallel::SerialParam(),
    only.pos = FALSE
  )
  
  comp_name <- paste0(Cluster_name$Cluster_Teo[i], "_vs_", Cluster_name$Cluster_B73[i])
  comp_name <- gsub("[^A-Za-z0-9_]", ".", comp_name)  
  maize_sub@tools$DEtest_custom$AllMarkers_wilcox$test_group <- comp_name
  maize_sub@tools$DEtest_custom$AllMarkers_wilcox$test_group <- as.factor(comp_name)
  
  de_df <- maize_sub@tools$DEtest_custom$AllMarkers_wilcox
  
  #4.2.设置阈值（和火山图一致）----
  
  de_df_filter <- de_df[de_df$avg_log2FC > 1 & de_df$p_val_adj < 0.05, ]
  
  #4.3.Add sheet names----
  
  addWorksheet(wb, sheetName = comp_name)
  writeData(wb, sheet = comp_name, x = de_df_filter)
  
}

#4.4.Save Excel File----

saveWorkbook(wb, file =paste0("Git_data/","DEG_AllComparisons.xlsx"), overwrite = TRUE)




