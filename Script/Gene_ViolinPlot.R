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

#2.1.Read Gene Expression----

combined_GEM <- readRDS(file.path(getwd(),"Git_data/combined_GEM_objects.rds"))

B3_GEM <- combined_GEM$B3_GEM
B7_GEM <- combined_GEM$B7_GEM
T1_GEM <- combined_GEM$T1_GEM
T3_GEM <- combined_GEM$T3_GEM

#2.2.Read Gene List----

Save_Flode <- "Violin"

Genelist <- read.xlsx(file.path(getwd(),"Git_data/Genelist.xlsx"),colNames = F)

colnames(Genelist)[1] <- "geneID"
Genelist <- unique(Genelist)

#2.3.Read Gene Annotation----

gene_annotation <- read.csv(file.path(getwd(),"Git_data/gene_annotation.csv"))

Genelist_annotation <- left_join(Genelist,gene_annotation,by = "geneID")

Genelist_annotation <- Genelist_annotation %>% select(geneID,maizeIDv5,maizeIDv4)

Genelist_annotation$maizeIDv5 <- ifelse(
  is.na(Genelist_annotation$maizeIDv5),  
  Genelist_annotation$maizeIDv4,            
  Genelist_annotation$maizeIDv5          
)

#3.Extract Gene Exprssion----

#3.1.Extract GeneList Expression----

B3_GEM <- B3_GEM %>% select(Cluster,everything())
B7_GEM <- B7_GEM %>% select(Cluster,everything())

T1_GEM <- T1_GEM %>% select(Cluster,everything())
T3_GEM <- T3_GEM %>% select(Cluster,everything())

Extract_GeneEx <- function(Sample_GEM,genelist){
  selected_columns <- c("Cluster", genelist$geneID[genelist$geneID %in% colnames(Sample_GEM)])
  Sample_GEM_filtered <- Sample_GEM[, selected_columns, drop = FALSE]
  return(Sample_GEM_filtered)
}

B3_ex <- Extract_GeneEx(B3_GEM,Genelist)
B7_ex <- Extract_GeneEx(B7_GEM,Genelist)
T1_ex <- Extract_GeneEx(T1_GEM,Genelist)
T3_ex <- Extract_GeneEx(T3_GEM,Genelist)


#3.2.Filter Cluster----

Cluster_dx <- c("Phl","Xyl","Scl","TB/EP")
B3_Cluster_ex <- B3_ex[B3_ex$Cluster %in% Cluster_dx, ]
B7_Cluster_ex <- B7_ex[B7_ex$Cluster %in% Cluster_dx, ]
T1_Cluster_ex <- T1_ex[T1_ex$Cluster %in% Cluster_dx, ]
T3_Cluster_ex <- T3_ex[T3_ex$Cluster %in% Cluster_dx, ]

B73_ex <- rbind(B3_Cluster_ex,B7_Cluster_ex)
Teo_ex <- rbind(T1_Cluster_ex,T3_Cluster_ex)
B73_ex$Sample <- "B73"
Teo_ex$Sample <- "Teo"

df1 <- rbind(B73_ex,Teo_ex)
df1 <- df1 %>% select(Sample,everything())

df_long <- df1 %>%
  pivot_longer(
    cols = starts_with("LOC"),  # 所有基因表达列
    names_to = "Gene",
    values_to = "Expression"
  )

#4.Plot Violin----

cluster_colors <- c(
  "Phl" = "#fbb4ae",
  "Xyl" = "#b3cde3",
  "Scl" = "#ccebc5",
  "TB/EP" = "#decbe4"
)

for (i in 1:nrow(Genelist_annotation)) {
  gene_of_interest <- Genelist_annotation$geneID[i]
  
  maizeid <- Genelist_annotation$maizeIDv5[i]
  
  print(paste0("目前正在绘制第",i,"个:",gene_of_interest))
  
  
  #4.1.绘图----
  p1 <-  ggplot(
    df_long %>% filter(Gene == gene_of_interest),
    aes(x = Sample, y = Expression, fill = Sample)
  ) +
    geom_violin(trim = T, width = 0.7) +
    stat_compare_means(
      method = "wilcox.test", 
      label = "p.signif",
      comparisons = list(c("B73", "Teo")),
      label.y = max(df_long$Expression, na.rm = TRUE) * 1
      # geom = "text",
      # tip.length = 0
      # ) + stat_compare_means(
      #   method = "wilcox.test",
      #   comparisons = list(c("B73", "Teo")),
      #   label = "p.format",
      #   label.y = max(df_long$Expression, na.rm = TRUE) * 0.52,  # 标签低一点
      #   # tip.length = 0,                                       # 显示连接线
      #   # geom = "point",
      #   # bracket.size = 0.5
    )+
    facet_wrap2(
      ~ Cluster,
      strip = strip_themed(
        background_x = elem_list_rect(fill = cluster_colors)
      ), scales = "free_x"
    ) +
    scale_fill_manual(values = c("B73" = "#008080", "Teo" = "#3CB371")) +
    ggtitle(maizeid) +
    theme_bw() +
    theme(
      panel.spacing = unit(0.2, "lines"),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(size = 1.2),
      axis.text.y = element_text(family = "SimHei", size = 10, color = "black"),
      text = element_text(size = 16),               # 整体文字大小
      axis.text = element_text(size = 16),          # 坐标轴刻度文字
      axis.title = element_text(size = 18),         # 坐标轴标题
      strip.text = element_text(size = 16),         # facet标签文字
      legend.text = element_text(size = 14),        # 图例文字
      legend.title = element_text(size = 16),        # 图例标题
      plot.title = element_text(
        size = 20,            # 字号
        # face = "bold",        # 加粗
        hjust = 0.5           # 居中
      )
    ) +
    coord_cartesian(ylim = c(min(df_long$Expression, na.rm = TRUE), max(df_long$Expression, na.rm = TRUE)*1.1))
  
  
  #4.2.Save Plot----
  
  ViolinPlot_dir <- file.path(getwd(),"Figures/ViolinPlot/Gene")
  if (!dir.exists(ViolinPlot_dir)) {
    dir.create(ViolinPlot_dir)
  }
  
  ggsave(plot = p1,paste0(ViolinPlot_dir,'/',maizeid,".png"),dpi = 600,width = 5,height = 5)
}
















