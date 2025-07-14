#1.Load R Packages----

library(Seurat)
library(tidyr)
library(dplyr)
library(openxlsx)
library(R.utils)
library(ggplot2)
library(SCP)
library(patchwork)

#2.Load R data----

maize_se <- readRDS(file.path(getwd(),"Git_data/maize_merge_new.rds"))

maize_se <- PrepSCTFindMarkers(maize_se,assay = "SCT", verbose = TRUE)
DefaultAssay(maize_se) <- "SCT"

#2.1.Extract B73 and Teo Sample----

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

#3.Read Gene list -----

Save_Flode <- "FeaturePlot"

Genelist <- read.xlsx(file.path(getwd(),"Git_data/Genelist.xlsx"),colNames = F)

colnames(Genelist)[1] <- "geneID"

#3.1.Read Gene Annotation----

gene_annotation <- read.csv(file.path(getwd(),"Git_data/gene_annotation.csv"))

Genelist_annotation <- left_join(Genelist,gene_annotation,by = "geneID")

Genelist_annotation <- Genelist_annotation %>% select(geneID,maizeIDv5,maizeIDv4)

#3.2.Read Teo Gene List----

Teo_list <- read.xlsx(file.path(getwd(),"Git_data/B73_Ames21814_OrthoFinder.xlsx"))

Genelist_Teo <- left_join(Genelist_annotation, Teo_list, by = c("maizeIDv5" = "B73ID"))
Genelist_Teo <- Genelist_Teo %>% select(geneID,maizeIDv5,maizeIDv4,TeoID)

Genelist_Teo <- unique(Genelist_Teo)

Genelist_Teo$maizeIDv5 <- ifelse(
  is.na(Genelist_Teo$maizeIDv5),  
  Genelist_Teo$geneID,            
  Genelist_Teo$maizeIDv5          
)

#4.Save Path----

Save_path <- file.path(getwd(),"Figures",Save_Flode)

if (!dir.exists(Save_path)) {
  dir.create(Save_path)
}

#5.Plot FeaturePlot----

#5.1.Plot FeatureDimPlot----

for (i in 1:nrow(Genelist_Teo)) {
  
  GeneID <- Genelist_Teo$geneID[i]
  maizeIDv5 <- Genelist_Teo$maizeIDv5[i]
  TeoID <- Genelist_Teo$TeoID[i]
  
  print(paste0("正在进行第",i,"个:",maizeIDv5,"_",TeoID))
  
  B73_plot <- SCP::FeatureDimPlot(
    srt = B73,
    features = GeneID,
    reduction = "UMAP",
    theme_use = "theme_blank",
    assay = "SCT",
    pt.size = 2,
    pt.alpha = 0.8,
    # label = T,
    # label.size = 4,
    # label.fg = "white",
    # label.bg = "black",
    # label_repel = T,
    # label_repulsion = 100,
    # label_insitu = T,
    # label_point_size = 1,
    theme_args = list(strip.text = ggplot2::element_blank())
  )+ 
    ggtitle(paste0(maizeIDv5)) +                     
    theme(plot.title = element_text(hjust = 0.5, size = 16)) ; B73_plot
  
  Teo_plot <- SCP::FeatureDimPlot(
    srt = Teo,
    features = GeneID,
    reduction = "UMAP",
    theme_use = "theme_blank",
    assay = "SCT",
    pt.size = 2,
    pt.alpha = 0.8,
    # label = T,
    # label.size = 4,
    # label.fg = "white",
    # label.bg = "black",
    # label_repel = T,
    # label_repulsion = 100,
    # label_insitu = T,
    # label_point_size = 1,
    theme_args = list(strip.text = ggplot2::element_blank())
  )+ 
    ggtitle(paste0(TeoID)) +                     
    theme(plot.title = element_text(hjust = 0.5, size = 16)) ; Teo_plot

  #5.2.Combind Plot----
  combined_plot <- B73_plot|Teo_plot
  
  #5.3.Save Combind_plot----
  
  UMAP_dir <- paste0(Save_path,"/UMAP/")
  
  if (!dir.exists(UMAP_dir)) {
    dir.create(UMAP_dir)
  }
  
  ggsave(plot = combined_plot,paste0(UMAP_dir,maizeIDv5,"_",TeoID,".png"),dpi = 600,width = 8,height = 5)
}

#5.4.Save Gene list----

write.xlsx(Genelist_Teo,paste0(Save_path,'/',Save_Flode,".xlsx"))

#6.SpatialFeaturePlot----

SpatialPlot_Function <- function(Sample,ID,titleID){
  SFP <- Seurat::SpatialFeaturePlot(Sample,
                                    features = ID,
                                    crop = T ,
                                    keep.scale = "feature",
                                    pt.size.factor = 2.1,
                                    image.alpha = 1,
                                    alpha = c(0.6, 1),
                                    stroke = NA,
                                    shape = 21
  )  + scale_fill_gradient(
    low = "white",
    high = "darkred",
    guide = guide_colourbar(  # 控制颜色图例
      frame.colour = "black",  # 色条边框颜色
      frame.linewidth = 1,     # 边框线宽
      ticks.colour = "black",  # 刻度线颜色（可选）
      ticks.linewidth = 0.5    # 刻度线粗细（可选）
    )
  ) + ggtitle(titleID) +
    theme(
      plot.title = element_text(
        size = 16,           # 字体大小
        # face = "bold",       # 粗体
        hjust = 0.5,         # 水平居中 (0.5)
        color = "black" ),   # 颜色
      legend.position = "right",     # 图例位置
      legend.title = element_blank() # 移除图例标题（可选）
    )
}

for (i in 1:nrow(Genelist_Teo)) {
  
  GeneID <- Genelist_Teo$geneID[i]
  maizeIDv5 <- Genelist_Teo$maizeIDv5[i]
  TeoID <- Genelist_Teo$TeoID[i]
  
  print(paste0("正在进行第",i,"个:",maizeIDv5,"_",TeoID))
  
  #6.1.Plot Feature Plot----
  
  B3_SPF <- SpatialPlot_Function(B3,GeneID,maizeIDv5)
  
  B7_SPF <- SpatialPlot_Function(B7,GeneID,maizeIDv5)
  
  T1_SPF <- SpatialPlot_Function(T1,GeneID,TeoID)
  
  T3_SPF <- SpatialPlot_Function(T3,GeneID,TeoID)
  
  combind_B3T1 <- B3_SPF|T1_SPF
  
  #6.2.Save Feature Plot----
  
  SFP_dir <- paste0(Save_path,"/SFP/")
  
  if (!dir.exists(SFP_dir)) {
    dir.create(SFP_dir)
  }
  
  gsf <- function(plot,id){
    ggsave(plot = plot,paste0(SFP_dir,id,".png"),dpi = 600,width = 5 ,height = 5)
  }
  
  gsf(B3_SPF,paste0("B3_",maizeIDv5))
  gsf(B7_SPF,paste0("B7_",maizeIDv5))
  gsf(T1_SPF,paste0("T1_",maizeIDv5))
  gsf(T3_SPF,paste0("T3_",maizeIDv5))
  gsf(combind_B3T1,paste0("B3T1_",maizeIDv5))
  
}



