#1.Load R Packages----

library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(scRNAtoolVis)
library(openxlsx)
library(pheatmap)
library(AnnotationHub)
library(AnnotationDbi)
library(clusterProfiler)
library(ggsankeyfier)
library(networkD3)
library(DOSE)
library(enrichplot)

#2.Load R data----

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

colnames(B73_Marker_Gene_Data) <- paste0("B73_",colnames(B73_Marker_Gene_Data))

B73_Heatmap <- pheatmap(t(as.matrix(B73_Marker_Gene_Data)), 
                        color = colorRampPalette(c("#2e77b4", "#f7f1ee", "#940e26"))(50),  # 定义颜色梯度
                        cluster_rows = F,  # 对行（Cluster）进行聚类
                        cluster_cols = F,  # 对列（GeneID）进行聚类
                        scale = "column",  # 因为转置后GeneID在列，对列标准化
                        main = "B73 Markers Gene Heatmap (Transposed)",
                        treeheight_row = 0,
                        fontsize_col = 10,
                        # angle_col = 45,
                        show_rownames = T,
                        show_colnames = F,
                        # gaps_row = col_gaps,      # 转置后用gaps_row
                        border_color = "black",  # 分隔线颜色
                        silent = TRUE   
)

svg(file.path(getwd(),"Figures/B73_GeneHeatmap.svg")), width = 6, height = 5)

print(B73_Heatmap)

dev.off()

colnames(Teo_Marker_Gene_Data) <- paste0("Teo_",colnames(Teo_Marker_Gene_Data))

Teo_Heatmap <- pheatmap(t(as.matrix(Teo_Marker_Gene_Data)), 
                        color = colorRampPalette(c("#2e77b4", "#f7f1ee", "#940e26"))(50),  # 定义颜色梯度
                        cluster_rows = F,  # 对行（Cluster）进行聚类
                        cluster_cols = F,  # 对列（GeneID）进行聚类
                        scale = "column",  # 因为转置后GeneID在列，对列标准化
                        main = "Teo Markers Gene Heatmap (Transposed)",
                        treeheight_row = 0,
                        fontsize_col = 10,
                        # angle_col = 45,
                        show_rownames = T,
                        show_colnames = F,
                        # gaps_row = col_gaps,      # 转置后用gaps_row
                        border_color = "black",  # 分隔线颜色
                        silent = TRUE   
)

svg(file.path(getwd(),"Figures/Teo_GeneHeatmap.svg")), width = 6, height = 5)

print(Teo_Heatmap)

dev.off()


#B GO Enrichment Plot----

B73_Markers_Gene <- FindAllMarkers(B73,
                                   assay = "SCT",
                                   logfc.threshold = 1,
                                   min.diff.pct = 0.2,
                                   group.by = "newclusters",
                                   only.pos = TRUE)

B73_allmarker_list <- B73_Markers_Gene |> group_by(cluster) |>
  filter(p_val_adj<0.001)|>ungroup()


Teo_Markers_Gene <- FindAllMarkers(Teo,
                                   assay = "SCT",
                                   logfc.threshold = 1,
                                   min.diff.pct = 0.2,
                                   group.by = "newclusters",
                                   only.pos = TRUE)

Teo_allmarker_list <- Teo_Markers_Gene |> group_by(cluster) |>
  filter(p_val_adj<0.001)|>ungroup()

# Run Go Enrichment-----

maize <- loadDb(file.path(getwd(),"Git_data/maize.OrgDb"))

#5.2.Extract ENTREZID----

B73_allmarker_list$gene <- sub("^LOC", "", B73_allmarker_list$gene)

Teo_allmarker_list$gene <- sub("^LOC", "", Teo_allmarker_list$gene)

#5.3.Renames Colnames----

colnames(B73_allmarker_list)[7] <- "ENTREZID"
colnames(Teo_allmarker_list)[7] <- "ENTREZID"

#5.4.Compare Cluster-----

B73_GO_enrich = compareCluster(ENTREZID ~ cluster,
                               data = B73_allmarker_list,
                               fun = "enrichGO",
                               pvalueCutoff=0.05,
                               ont = "BP", 
                               OrgDb= maize )

B73_GO_Enrich_df <- data.frame(B73_GO_enrich)

Teo_GO_enrich = compareCluster(ENTREZID ~ cluster,
                               data = Teo_allmarker_list,
                               fun = "enrichGO",
                               pvalueCutoff=0.05,
                               ont = "BP", 
                               OrgDb= maize )

Teo_GO_Enrich_df <- data.frame(Teo_GO_enrich)

GO_df_modif1 <- function(Go_df,Accession){
  GO_Enrich_df1 <- Go_df %>% select(-cluster)
  GO_Enrich_df1$Accession <- Accession
  GO_Enrich_df2 <- GO_Enrich_df1 %>% select("Accession",everything())
  GO_Enrich_df3 <- GO_Enrich_df2 %>% group_by(Cluster) %>% arrange(p.adjust) %>% slice_min(p.adjust,n=3)
  return(GO_Enrich_df3)
}
Teo_GO_Enrich_df1 <- GO_df_modif1(Teo_GO_Enrich_df,"Teo")
B73_GO_Enrich_df1 <- GO_df_modif1(B73_GO_Enrich_df,"B73")
Go_Enrich <- rbind(Teo_GO_Enrich_df1,B73_GO_Enrich_df1)

##

cnetplot(Teo_GO_enrich, 
         layout = igraph::layout_with_fr,
         showCategory = 3,  # 显示每个cluster的前5个通路
         foldChange = Teo_allmarker_list$log2FoldChange, # 可选：添加logFC颜色
         color_item = "#B3B3B3",
         color_category = "#E5C494",
         node_label = "category",
         circular = T,   # 环形布局
         colorEdge = TRUE) +
  ggtitle("Gene-Pathway Network") +
  theme(legend.position = "right")


#5.5.Draw GO Enrichment Plot----

Go_Enrich <- Go_Enrich %>%
  mutate(Description = fct_reorder(Description, Count, .desc = TRUE))

p0 <- ggplot(data = Go_Enrich, aes(
  x = Description,
  y = Count,
  fill = Cluster
)) +
  geom_bar(
    stat = "identity",
    width = 0.8
  ) +
  scale_fill_manual(
    values = c("Xyl" = "#A6CEE3", 
               "Phl" = "#B2DF8A", 
               "Scl" = "#FB9A99",
               "TB/EP" = "#E31A1C",
               "LDGT" = "#1F78B4",
               "EDGT" = "#33A02C",
               "Hypo" = "#FDBF6F",
               "Epi" = "#FF7F00",
               "Unk1" = "#CAB2D6",
               "Unk2" = "#6A3D9A",
               "Unk3" = "#FFFF99")
  ) +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Gene Count",
    # title = "GO Enrichment"
  ) +
  theme_minimal() +
  theme(
    axis.text.y  = element_text(size = 10),
    legend.position = "right",
    
    # 🔥 取消网格
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # ✅ y轴保留完整线
    axis.line.x = element_line(color = "black"),
    
    # ✅ x轴只保留刻度，不要轴线
    axis.line.y = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black")
  )

ggsave(
  filename = file.path(getwd(),"Figures/Gene_GoEnrichment.svg"),  # 文件名
  plot = p0,              # 图形对象
  device = "svg",        # 设备类型
  width = 8,             # 宽度（英寸）
  height = 10,            # 高度（英寸）
  dpi = 600              # 分辨率（影响文本和点的大小）
)








