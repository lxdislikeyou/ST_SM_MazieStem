#1.Load R Packages-----

library(Seurat)
library(ggplot2)
library(R.utils)
library(openxlsx)
library(patchwork)
library(officer)
library(rvg)

#2.Load Maize RDS File-----

maize <- readRDS(file.path(getwd(),"Git_data/maize_merge_new.rds"))

#3.Draw mz Spatial Plot----

#3.1.Read mz File

neg <- read.xlsx(file.path(getwd(),"Git_data/Amino_acid_neg.xlsx"))
pos <- read.xlsx(file.path(getwd(),"Git_data/Amino_acid_pos.xlsx"))

neg <- neg %>% select(mz,Metabolites)
neg_list <- unique(neg)
neg_mz <- neg_list$mz
neg_names <- neg_list$Metabolites

pos <- pos %>% select(mz,Metabolites)
pos_list <- unique(pos)
pos_mz <- pos_list$mz
pos_names <- pos_list$Metabolites

#4.3.Spatial Plot

#4.3.1.Neg Spatial Plot

DefaultAssay(maize) <- "Metabolism_neg"

neg_path <- file.path(getwd(),"Figures/Metabolism_SpatialPlot")

if (!dir.exists(neg_path)) {
  dir.create(neg_path)
}

doc <- read_pptx()

for (i in 1:length(neg_mz)) {
  
  p_list <- SpatialFeaturePlot(maize,
                               features = neg_mz[i],
                               pt.size.factor = 2,
                               combine = FALSE,
                               min.cutoff = 2,
                               alpha = c(1, 0.6),
                               # min.cutoff = 1,
                               crop = FALSE)
  
  # 修改颜色（每个图都应用）
  p_list <- lapply(p_list, function(p) {
    p + 
      scale_fill_gradientn(colors = c("#2e77b4", "#f7f1ee", "#940e26")) +
      theme(legend.position = "none")
  })
  
  combined_plot <- wrap_plots(p_list, ncol = 2)
  
  neg_name <- neg_names[i]
  
  for (x in 1:4) {
    
    name_list <- c("B3","B7","T1","T3")
    name <- name_list[x]

    ggsave(plot = p_list[[x]] , paste0(neg_path,'/',name,"_",neg_name,".png"),
           dpi = 300,width = 5,height = 5)
  }
  
  doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
  doc <- ph_with(doc, value = neg_name , location = ph_location_type(type = "title"))
  doc <- ph_with(doc,
                 dml(ggobj = combined_plot),
                 location = ph_location_type(type = "body"))
}
print(doc, target = paste0(neg_path,"/neg_spatial.pptx"))


#4.3.2.pos Spatial Plot

DefaultAssay(maize) <- "Metabolism_pos"
pos_path <- neg_path  

for (i in 1:length(pos_mz)) {
  
  p_list <- SpatialFeaturePlot(maize,
                               features = pos_mz[i],
                               pt.size.factor = 2,
                               combine = FALSE,
                               min.cutoff = 2,
                               alpha = c(1, 0.6),
                               # min.cutoff = 1,
                               crop = FALSE)
  
  # 修改颜色（每个图都应用）
  p_list <- lapply(p_list, function(p) {
    p + 
      scale_fill_gradientn(colors = c("#F8F8FF","#F5F5DC","#FA8072","#800000")) +
      theme(legend.position = "none")
  })
  
  combined_plot <- wrap_plots(p_list, ncol = 2)
  
  pos_name <- pos_names[i]
  
  for (x in 1:4) {
    
    name_list <- c("B3","B7","T1","T3")
    name <- name_list[x]
    
    ggsave(plot = p_list[[x]] , 
           paste0(pos_path,'/',name,"_",pos_name,".png"),
           dpi = 300,width = 5,height = 5)
  }
  
  doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
  doc <- ph_with(doc, value = pos_name , location = ph_location_type(type = "title"))
  doc <- ph_with(doc,
                 dml(ggobj = combined_plot),
                 location = ph_location_type(type = "body"))
}

print(doc, target =paste0(pos_path,"/pos_spatial.pptx"))







