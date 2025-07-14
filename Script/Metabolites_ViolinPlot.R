#1.Load R Packages----

library(tidyr)
library(dplyr)
library(openxlsx)
library(R.utils)
library(ggplot2)
library(patchwork)
library(ggh4x)
library(ggpubr)

#2.Load R data----

#2.1.Read Metabolism Levels----

Save_path <- "./Metabolism_Levels"

combined_GEM <- readRDS(file.path(getwd(),"Git_data/BT_combined_MEM_objects.rds"))

B73_neg <- combined_GEM$B73_MEM_neg
B73_pos <- combined_GEM$B73_MEM_pos
Teo_neg <- combined_GEM$Teo_MEM_neg
Teo_pos <- combined_GEM$Teo_MEM_pos

#2.2.Read Metabolism List----

Ac_neg <- read.xlsx(file.path(getwd(),"Git_data/Amino_acid_neg.xlsx"))
Ac_pos <- read.xlsx(file.path(getwd(),"Git_data/Amino_acid_pos.xlsx"))

Ac_neg_list <- Ac_neg %>% select(mz,Metabolites)
Ac_pos_list <- Ac_pos %>% select(mz,Metabolites)

Ac_list <- rbind(Ac_neg_list,Ac_pos_list)

#3.Extract Ac Matrix----

Extract_MetaEx <- function(Sample_MEM,metalist){
  selected_columns <- c("Cluster", metalist$mz[metalist$mz %in% colnames(Sample_MEM)])
  Sample_MEM_filtered <- Sample_MEM[, selected_columns, drop = FALSE]
  return(Sample_MEM_filtered)
}

B73_meta_neg <- Extract_MetaEx(B73_neg,Ac_neg_list)
B73_meta_pos <- Extract_MetaEx(B73_pos,Ac_pos_list)
Teo_meta_neg <- Extract_MetaEx(Teo_neg,Ac_neg_list)
Teo_meta_pos <- Extract_MetaEx(Teo_pos,Ac_pos_list)

#3.2.Filter Cluster----

Cluster_dx <- c("Phl","Xyl","Scl","TB/EP")
B73_Cluster_neg <- B73_meta_neg[B73_meta_neg$Cluster %in% Cluster_dx, ]
B73_Cluster_pos <- B73_meta_pos[B73_meta_pos$Cluster %in% Cluster_dx, ]
Teo_Cluster_neg <- Teo_meta_neg[Teo_meta_neg$Cluster %in% Cluster_dx, ]
Teo_Cluster_pos <- Teo_meta_pos[Teo_meta_pos$Cluster %in% Cluster_dx, ]

B73_Cluster_neg$Sample <- "B73"
B73_Cluster_pos$Sample <- "B73"
Teo_Cluster_neg$Sample <- "Teo"
Teo_Cluster_pos$Sample <- "Teo"

B73_Cluster_neg <- B73_Cluster_neg %>% select(Sample,everything())
B73_Cluster_pos <- B73_Cluster_pos %>% select(Sample,everything())
Teo_Cluster_neg <- Teo_Cluster_neg %>% select(Sample,everything())
Teo_Cluster_pos <- Teo_Cluster_pos %>% select(Sample,everything())

B73_Cluster_neg_long <- B73_Cluster_neg %>%
  pivot_longer(
    cols = -c(Sample, Cluster),,  # 所有基因表达列
    names_to = "Metabolite",
    values_to = "Level"
  )
B73_Cluster_pos_long <- B73_Cluster_pos %>%
  pivot_longer(
    cols = -c(Sample, Cluster),,  # 所有基因表达列
    names_to = "Metabolite",
    values_to = "Level"
  )

B73_Cluster_long <- rbind(B73_Cluster_neg_long,B73_Cluster_pos_long)


Teo_Cluster_neg_long <- Teo_Cluster_neg %>%
  pivot_longer(
    cols = -c(Sample, Cluster),,  # 所有基因表达列
    names_to = "Metabolite",
    values_to = "Level"
  )
Teo_Cluster_pos_long <- Teo_Cluster_pos %>%
  pivot_longer(
    cols = -c(Sample, Cluster),,  # 所有基因表达列
    names_to = "Metabolite",
    values_to = "Level"
  )

Teo_Cluster_long <- rbind(Teo_Cluster_neg_long,Teo_Cluster_pos_long)

BT_Cluster_long <- rbind(B73_Cluster_long,Teo_Cluster_long)

BT_Cluster_long$Level <- log2(BT_Cluster_long$Level + 1)


#4.Plot Violin----

cluster_colors <- c(
  "Phl" = "#fbb4ae",
  "Xyl" = "#b3cde3",
  "Scl" = "#ccebc5",
  "TB/EP" = "#decbe4"
)


Ac_list <- unique(Ac_list)
for (i in 1:nrow(Ac_list)) {
  i=3
  mz_of_interest <- Ac_list$mz[i]
  metaid <- Ac_list$Metabolites[i]
  
  print(paste0("目前正在绘制第",i,"个:",metaid))
  
  p1 <- ggplot(
    BT_Cluster_long %>% filter(Metabolite == mz_of_interest),
    aes(x = Sample, y = Level, fill = Sample)
  ) +
    geom_violin(trim = TRUE, width = 0.7, alpha = 0.7) +  # 调整透明度避免遮挡散点
    geom_jitter(
      width = 0.2 ,     # 控制散点的水平抖动范围（避免重叠）
      size = 0.5,         # 点的大小
      alpha = 0.3,      # 点的透明度
      shape = 21,       # 21是带填充色的点
      color = "black"   # 点的边框颜色
    ) +
    stat_compare_means(
      method = "wilcox.test", 
      label = "p.signif",
      comparisons = list(c("B73", "Teo")),
      label.y = max(BT_Cluster_long$Level, na.rm = TRUE) * 1
    ) +
    facet_wrap2(
      ~ Cluster,
      strip = strip_themed(
        background_x = elem_list_rect(fill = cluster_colors)
      ), scales = "free_x"
    ) +
    scale_fill_manual(values = c("B73" = "#008080", "Teo" = "#3CB371")) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 5)  # 自动生成5个合理的Y轴刻度
    ) +
    ggtitle(metaid) +
    theme_bw() +
    theme(
      panel.spacing = unit(0.2, "lines"),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(size = 1.2),
      axis.text.y = element_text(family = "SimHei", size = 10, color = "black"),
      text = element_text(size = 16),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      strip.text = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      plot.title = element_text(
        size = 20,
        hjust = 0.5
      )
    ) +
    coord_cartesian(ylim = c(min(BT_Cluster_long$Level, na.rm = TRUE), max(BT_Cluster_long$Level, na.rm = TRUE)*1.1))
  
  test_result <- wilcox.test(
    Level ~ Sample, 
    data = BT_Cluster_long %>% filter(Metabolite == mz_of_interest & Sample %in% c("B73", "Teo")))
  print(test_result)
  
  BT_Cluster_long %>%
    filter(Metabolite == mz_of_interest) %>%
    group_by(Sample) %>%
    summarise(Median = median(Level, na.rm = TRUE))
  
  #4.2.Save Plot----
  
  ViolinPlot_dir <- file.path(getwd(),"Figures/ViolinPlot/Metabolism")
  if (!dir.exists(ViolinPlot_dir)) {
    dir.create(ViolinPlot_dir)
  }
  ggsave(plot = p1,paste0(ViolinPlot_dir,'/',metaid,".png"),dpi = 600,width = 5,height = 5)
  
}






