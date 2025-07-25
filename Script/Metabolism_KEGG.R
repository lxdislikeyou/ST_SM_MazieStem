#1.Load R Packages----

library(patchwork)
library(ggplot2)
library(cowplot)
library(scRNAtoolVis)
library(openxlsx)
library(AnnotationHub)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(openxlsx)
library(KEGGREST)
library(FELLA)
library(ggraph)
library(magick)

#2.Read Metabolism File----

Mt_neg <- read.xlsx(file.path(getwd(),"Git_data/Qualitative.xlsx"),sheet = 2)
Mt_pos <- read.xlsx(file.path(getwd(),"Git_data/Qualitative.xlsx"),sheet = 4)

Mt_neg_list <- Mt_neg %>% select(Metabolites,mz,cid,KEGG)
Mt_pos_list <- Mt_pos %>% select(Metabolites,mz,cid,KEGG)

Mt_list <- rbind(Mt_neg_list,Mt_pos_list)
Mt_list <- unique(Mt_list)
Mt_list <- Mt_list[!is.na(Mt_list$KEGG),]

#3.Run Kegg Enrichment----

#3.1.Build Graph from Kegg====

# kegg_graph <- buildGraphFromKEGGREST(organism = "zma")

#3.2.Build fella_zma_db Data from Graph====
# 
# buildDataFromGraph(
#   keggdata.graph = kegg_graph,
#   databaseDir = "fella_zma_db_v2",   # 可自定义名称
#   internalDir = T ,            # 保存到工作目录而非包内部
#   matrices = c("hypergeom", "diffusion", "pagerank"),  # 可选哪些矩阵
#   normality = c("diffusion", "pagerank"),              # 哪些方法做标准化
#   dampingFactor = 0.85,
#   niter = 500
# )

#3.3.Load Kegg Data====

fella.data <- loadKEGGdata(
  databaseDir = "fella_zma_db_v2",  
  internalDir = T,           
  loadMatrix = c("diffusion", "pagerank","hypergeom") 
)

#3.4.Read Different Metabolites====

sheet_names <- c("EDGT",'Epi','Hypo','LDGT','Phl','Scl','TB.EP','Unk1','Unk2','Unk3','Xyl')

cluster_colors <- c("EDGT" = "#33A02C",
                    "Epi" = "#FF7F00",
                    "Hypo" = "#FDBF6F",
                    "LDGT" = "#1F78B4",
                    "Phl" = "#B2DF8A", 
                    "Scl" = "#FB9A99",
                    "TB/EP" = "#E31A1C",
                    "Unk1" = "#CAB2D6",
                    "Unk2" = "#6A3D9A",
                    "Unk3" = "#FFFF99",
                    "Xyl" = "#A6CEE3"
)

KEGG_Photo <- file.path(getwd(),"Figures/Kegg")

if (!dir.exists(KEGG_Photo)) {
  dir.create(KEGG_Photo)
}

for (i in 1:length(sheet_names)) {
  
print(paste0("目前正在绘制第",i,"个组织:",sheet_names[i],"的KEGG"))

sig_list <- read.xlsx(file.path(getwd(),"Git_data/KEGG/Diff_metabolites.xlsx"),
                      sheet = i ,colNames = F)
colnames(sig_list)[1] <- "mz"
Metabolites_list <- sig_list$mz
Tissue_mt <- Mt_list[Mt_list$mz %in% Metabolites_list,]
kegg_ids <- (Tissue_mt$KEGG)
kegg_ids <- unique(kegg_ids) 

#3.5.Run Kegg Analysis====

myAnalysis <- enrich(
  compounds = kegg_ids,
  method = listMethods(),      #===parameter : "hypergeom", "diffusion", "pagerank" or listMethods()
  approx = "normality",     #===parameter : "normality","simulation"
  data = fella.data          
)


#3.6.Generate Results Table====

result_table <- generateResultsTable(
  object = myAnalysis,
  method = "hypergeom",
  data = fella.data,
  threshold = 0.05,
  plimit = 15,
  nlimit = 250,
  LabelLengthAtPlot = 40 
)

#3.7.Generate Results Graph File====

myGraph1 <- generateResultsGraph(
  object = myAnalysis, 
  method = "hypergeom", 
  threshold = 0.05, 
  data = fella.data)

igraph::V(myGraph1)$mycolor <- ifelse(igraph::V(myGraph1)$type, cluster_colors[i], "coral1")

p1 <- local({
  current_graph <- myGraph1  # 复制当前图的副本
  ggraph(current_graph) +  
    geom_edge_link(alpha = 0.3) +
    geom_node_point(
      aes(
        color = I(igraph::V(current_graph)$mycolor),  # 引用当前图的属性
        size = factor(igraph::V(current_graph)$type)
      )
    ) +  
    scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4)) +
    geom_node_text(aes(label = igraph::V(current_graph)$label), repel = TRUE, size = 2) +
    theme_void() +
    theme(legend.position = "none")
})

ggsave(plot = p1,paste0(KEGG_Photo,'/',sheet_names[i],"_Kegg.png"),width = 5,height = 5,dpi = 600)

}

######==================#######

# 设置图片目录
KEGG_Photo <- file.path(getwd(), "Figures/Kegg")

# 定义组织名称顺序（必须与保存图片时的名称一致）
sheet_names <- c("EDGT", "Epi", "Hypo", "LDGT", "Phl", "Scl", "TB.EP", "Unk1", "Unk2", "Unk3", "Xyl")

# 获取所有组织对应的图片文件
image_files <- list.files(KEGG_Photo, pattern = "_Kegg\\.png$", full.names = TRUE)

# 提取基础名称（不带路径和扩展名）
base_names <- gsub("_Kegg\\.png$", "", basename(image_files))

# 确保顺序与sheet_names一致
image_files <- image_files[match(sheet_names, base_names)]

# 读取所有图片并转换为ggplot对象
plot_list <- lapply(image_files, function(file) {
  img <- image_read(file)
  ggplot() + 
    annotation_custom(grid::rasterGrob(img, interpolate = TRUE)) + 
    theme_void()
})

# 使用patchwork拼接（3列布局）
combined_plot <- wrap_plots(plot_list, ncol = 3)  # 可以调整列数

# 保存拼接后的图
ggsave(paste0(KEGG_Photo, "/combined_KEGG.png"), 
       combined_plot, 
       width = 15,  # 宽度根据列数调整
       height = ceiling(length(plot_list)/3) * 5,  # 高度根据行数调整
       dpi = 600)

