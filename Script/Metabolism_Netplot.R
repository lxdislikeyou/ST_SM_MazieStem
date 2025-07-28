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

#3.5.Set Cluster Colors---- 

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

#3.6.Set Save Path----

KEGG_Photo <- file.path(getwd(),"Figures/Kegg")
KEGG_list <- file.path(getwd(),"Git_data/KEGG")

if (!dir.exists(KEGG_Photo)) {
  dir.create(KEGG_Photo)
}
if (!dir.exists(KEGG_list)) {
  dir.create(KEGG_list)
}

#3.7.Read File ----

sig_list <- read.xlsx(file.path(getwd(), "Git_data/KEGG/Diff_metabolites_with_direction.xlsx"))

sig_list <- sig_list %>%
  group_by(Cluster) %>%
  arrange(p.value) %>%     # 从最小开始排序（越负越小）
  slice_head(n = 100) %>%      # 每组取前100个
  ungroup()

#4.Run KEGG Enrichment----

#4.1.Create an Empty List====

Grap_list <- list()

clusters <- unique(sig_list$Cluster)

#4.2.Run Loop====

for (clust in clusters) {
  message("Analyzing the cluster: ", clust)
  
    # 4.2.1.Significant Metabolites in The Current Cluster
  sub_sig <- sig_list %>% filter(Cluster == clust)
  Metabolites_list <- sub_sig$mz
  
  # 4.2.2.Filter Metabolites Matching KEGG IDs
  Tissue_mt <- Mt_list[Mt_list$mz %in% Metabolites_list, ]
  kegg_ids <- unique(Tissue_mt$KEGG)
  
  # 4.2.3.Skip unmatched cases
  if (length(kegg_ids) == 0 || all(is.na(kegg_ids))) {
    message("Skip Cluster ", clust, "，Invalid KEGG ID")
    next
  }
  
  # 4.3 Run FELLA Analysis====
  myAnalysis <- enrich(
    compounds = kegg_ids,
    method = "hypergeom",
    approx = "normality",
    data = fella.data
  )
  
  # 4.4 Generate Results Graph====
  myGraph <- generateResultsGraph(
    object = myAnalysis,
    method = "hypergeom",
    threshold = 0.05,
    data = fella.data
  )
  
  # 4.5.Check for Plot Generation====
  if (is.null(myGraph)) {
    message("Skip Cluster ", clust, "，The plot is empty")
    next
  }
  
  #4.5 Change to Data.frame，Add Cluster Information====
  g_df <- igraph::as_long_data_frame(myGraph)
  g_df$Cluster <- clust
  
  # 4.6 Save Grap_list====
  Grap_list[[clust]] <- g_df
}

#4.7.Merge All Grap List====

Grap_df_all <- bind_rows(Grap_list)
  
#5.1.Ensure Both 'from' and 'to' are Character Type====

Grap_df_all$from <- as.character(Grap_df_all$from)
Grap_df_all$to   <- as.character(Grap_df_all$to)

#5.2.Add cluster prefix to ensure uniqueness====

Grap_df_all <- Grap_df_all %>%
  mutate(
    from = paste0(Cluster, "_", from),
    to   = paste0(Cluster, "_", to)
  )

#5.3.Identify hub nodes (all nodes appearing in 'to' column)====

hub_nodes <- unique(Grap_df_all$to)

#5.4.Construct igraph object====

g_all <- igraph::graph_from_data_frame(
  d = Grap_df_all[, c("from", "to")],
  directed = FALSE
)

#5.5. Assign cluster information to nodes====

# (first try matching with 'from', fall back to 'to' if no match)
cluster_from <- Grap_df_all[, c("from", "Cluster")] %>% distinct()
cluster_to <- Grap_df_all[, c("to", "Cluster")] %>% distinct()

cluster_map_from <- setNames(cluster_from$Cluster, cluster_from$from)
cluster_map_to <- setNames(cluster_to$Cluster, cluster_to$to)

all_nodes <- igraph::V(g_all)$name
cluster_values <- sapply(all_nodes, function(n) {
  if (!is.na(cluster_map_from[n])) {
    cluster_map_from[n]
  } else if (!is.na(cluster_map_to[n])) {
    cluster_map_to[n]
  } else {
    NA
  }
})
igraph::V(g_all)$Cluster <- cluster_values

#5.6.Label hub nodes====

igraph::V(g_all)$Type <- ifelse(all_nodes %in% hub_nodes, "hub", "node")

#5.7.Set node styles====

vertex_sizes <- ifelse(igraph::V(g_all)$Type == "hub", 4, 3)
vertex_shapes <- ifelse(igraph::V(g_all)$Type == "hub", "circle", "circle")

#5.8.Custom color mapping vector for Clusters====

cluster_colors <- c(
  "EDGT" = "#33A02C",
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

#6. Assign colors: hub nodes use Cluster colors, regular nodes use gray====
vertex_colors <- ifelse(
  igraph::V(g_all)$Type == "hub",
  sapply(igraph::V(g_all)$Cluster, function(cl) {
    if (!is.na(cl) && cl %in% names(cluster_colors)) cluster_colors[[cl]] else "grey80"
  }),
  "grey80"
)

#7. Set node border colors: red for hubs, none (NA) for regular nodes====
vertex_frame_colors <- ifelse(igraph::V(g_all)$Type == "hub", "red", NA)

#7.1—— New section: Create to_label mapping and assign to node label attribute —— #
to_label_map <- Grap_df_all %>%
  dplyr::select(to, to_label) %>%
  dplyr::distinct()

to_label_vec <- setNames(to_label_map$to_label, to_label_map$to)

#7.2.Default node labels use node names====
labels <- rep("", length(all_nodes))
names(labels) <- all_nodes

#7.3.Only for hub nodes that have to_label available, replace labels with to_label====
labels[all_nodes %in% names(to_label_vec) & igraph::V(g_all)$Type == "hub"] <-
  to_label_vec[all_nodes[all_nodes %in% names(to_label_vec) & igraph::V(g_all)$Type == "hub"]]

igraph::V(g_all)$label <- labels


#8.Draw NetPlot----

#8.1.Netplot layout====

set.seed(123)
layout <- igraph::layout_with_graphopt(g_all, charge = 0.03, niter = 1500)

#8.2.Draw Netplot====
igraph::plot.igraph(
  g_all,
  layout = layout,
  vertex.size = vertex_sizes,
  vertex.color = vertex_colors,
  vertex.shape = vertex_shapes,
  vertex.label = igraph::V(g_all)$label,
  vertex.label.cex = 0.4,
  vertex.label.color = "black",
  vertex.frame.color = vertex_frame_colors,
  edge.color = adjustcolor("gray30", alpha.f = 0.3),
  edge.width = 1,
  rescale = T,
  margin = 0
)
  
  
  
  














