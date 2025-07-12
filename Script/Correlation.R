#1.Load R Packages----

library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)
library(grid)
library(gridExtra)
library(officer)
library(rvg)
library(patchwork)
library(networkD3)
library(tidyverse)
library(htmlwidgets)
library(webshot2)
library(ggrepel)

#2.Read RDS data-----

#2.1.Read Gene Expression----

BT_Ex <- readRDS(file.path(getwd(),"Git_data/BT_combined_GEM_objects.rds"))

B73_Ex <- BT_Ex[[1]]
Teo_Ex <- BT_Ex[[2]]

#2.2.Read Metabolism Levels----

SM_BT <- readRDS(file.path(getwd(),"Git_data/BT_combined_MEM_objects.rds"))

B73_SM_neg <- SM_BT[[1]]
B73_SM_pos <- SM_BT[[2]]
Teo_SM_neg <- SM_BT[[3]]
Teo_SM_pos <- SM_BT[[4]]

#3.Read list-----

#3.1. AA transport ====

Transport_protein <- read.xlsx("D:/Work/Spatial_Project/空间转录组点对点/Spatial_project/Transport_Protein/Transport_protein.xlsx")

AA_Tr <- c("AAP","UMAMIT")

AA_list <- Transport_protein[Transport_protein$Family %in% AA_Tr,]
AA_list <- AA_list[order(AA_list$Family, decreasing = TRUE), ]
AA_list <- AA_list[!is.na(AA_list$geneID),]

#3.2.Read AA list====

Ac_neg <- read.xlsx(file.path(getwd(),'Git_data','Amino_acid_neg.xlsx'))
Ac_pos <- read.xlsx(file.path(getwd(),'Git_data','Amino_acid_pos.xlsx'))
Ac_neg <- Ac_neg %>% select("mz","Metabolites")
Ac_pos <- Ac_pos %>% select("mz","Metabolites")
Ac_neg <- unique(Ac_neg)
Ac_pos <- unique(Ac_pos)
Ac_all <- rbind(Ac_neg,Ac_pos)


#4.Extract Gene Expression and Metabolites Levels----
AA_df1 <- AA_list %>% select(geneID,maizeIDv5,Family)

B73_df1 <- B73_Ex %>% select(Cluster,any_of(AA_df1$geneID))
Teo_df1 <- Teo_Ex %>% select(Cluster,any_of(AA_df1$geneID))

B73_aa_sum <- B73_df1 %>% 
  mutate(Total_ex = rowSums(across(where(is.numeric)), na.rm = TRUE))
B73_aa_sum$Cell <- rownames(B73_aa_sum)
# B73_aa_filter <- B73_aa_sum %>% select(Cell,Cluster,Total_ex)
B73_aa_filter <- B73_aa_sum 
Teo_aa_sum <- Teo_df1 %>% 
  mutate(Total_ex = rowSums(across(where(is.numeric)), na.rm = TRUE))
Teo_aa_sum$Cell <- rownames(Teo_aa_sum)
# Teo_aa_filter <- Teo_aa_sum %>% select(Cell,Cluster,Total_ex)
Teo_aa_filter <- Teo_aa_sum 

#5.Merge ST and SM ----

B73_aa_neg <- B73_SM_neg[,colnames(B73_SM_neg) %in% Ac_neg$mz]
B73_aa_neg$Cluster <- B73_SM_neg$Cluster
B73_aa_neg$Cell <- rownames(B73_aa_neg)
Teo_aa_neg <- Teo_SM_neg[,colnames(Teo_SM_neg) %in% Ac_neg$mz]
Teo_aa_neg$Cluster <- Teo_SM_neg$Cluster
Teo_aa_neg$Cell <- rownames(Teo_aa_neg)
B73_aa_pos <- B73_SM_pos[,colnames(B73_SM_pos) %in% Ac_pos$mz]
B73_aa_pos$Cluster <- B73_SM_pos$Cluster
B73_aa_pos$Cell <- rownames(B73_aa_pos)
Teo_aa_pos <- Teo_SM_pos[,colnames(Teo_SM_pos) %in% Ac_pos$mz]
Teo_aa_pos$Cluster <- Teo_SM_pos$Cluster
Teo_aa_pos$Cell <- rownames(Teo_aa_pos)

B73_BT <- inner_join(B73_aa_filter,B73_aa_neg,by = c('Cell','Cluster'))
B73_BT <- inner_join(B73_BT,B73_aa_pos,by = c('Cell','Cluster'))

Teo_BT <- inner_join(Teo_aa_filter,Teo_aa_neg,by = c('Cell','Cluster'))
Teo_BT <- inner_join(Teo_BT,Teo_aa_pos,by = c('Cell','Cluster'))

#6.Run Pearson Correlation----

AA_df1 <- unique(AA_df1)
#6.1.Create metabolite_col and Gene_cols====
Teo_metabolite_cols <- setdiff(colnames(Teo_BT), c(colnames(Teo_aa_filter)))
gene_cols <- setdiff(colnames(Teo_BT), c("Cell", "Cluster", Teo_metabolite_cols))
B73_metabolite_cols <- setdiff(colnames(B73_BT), c(colnames(B73_aa_filter)))
gene_cols <- setdiff(colnames(B73_BT), c("Cell", "Cluster", B73_metabolite_cols))
#6.2.Run Pearson Anaylsis====
Teo_cor_results <- Teo_BT %>%
  group_by(Cluster) %>%
  summarise(
    results = list({
      df <- cur_data()
      expand_grid(Gene = gene_cols, Metabolite = Teo_metabolite_cols) %>%
        mutate(
          cor_result = map2(Gene, Metabolite, function(g, m) {
            cor_test <- cor.test(
              x = df[[g]],
              y = df[[m]],
              method = "pearson"
            )
            tibble(
              Pearson_r = cor_test$estimate,
              P_value = cor_test$p.value
            )
          })
        ) %>%
        unnest(c(cor_result))
    }),
    .groups = "drop"
  ) %>%
  unnest(results)

B73_cor_results <- B73_BT %>%
  group_by(Cluster) %>%
  summarise(
    results = list({
      df <- cur_data()
      expand_grid(Gene = gene_cols, Metabolite = B73_metabolite_cols) %>%
        mutate(
          cor_result = map2(Gene, Metabolite, function(g, m) {
            cor_test <- cor.test(
              x = df[[g]],
              y = df[[m]],
              method = "pearson"
            )
            tibble(
              Pearson_r = cor_test$estimate,
              P_value = cor_test$p.value
            )
          })
        ) %>%
        unnest(c(cor_result))
    }),
    .groups = "drop"
  ) %>%
  unnest(results)

#6.3.Change Data format====

Teo_cor_results_df <- as.data.frame(Teo_cor_results)
Teo_cor_results_df$Accession <- 'Teo'
B73_cor_results_df <- as.data.frame(B73_cor_results)
B73_cor_results_df$Accession <- 'B73'

#6.4.Modif mz to Metabolites and geneID to MaizeIDv5====
Teo_cor_results_df <- Teo_cor_results_df %>%
  left_join(Ac_all %>% select(mz, Metabolites), 
            by = c("Metabolite" = "mz")) %>%
  mutate(Metabolite = Metabolites) %>%  # 替换原有列
  select(-Metabolites) 
Teo_cor_results_df <- Teo_cor_results_df %>%
  left_join(AA_df1 %>% select(geneID, maizeIDv5), 
            by = c("Gene" = "geneID")) %>%
  mutate(Gene = maizeIDv5) %>%  # 替换原有列
  select(-Gene) 
Teo_cor_results_df$maizeIDv5[is.na(Teo_cor_results_df$maizeIDv5)] <- "Total" 

B73_cor_results_df <- B73_cor_results_df %>%
  left_join(Ac_all %>% select(mz, Metabolites), 
            by = c("Metabolite" = "mz")) %>%
  mutate(Metabolite = Metabolites) %>%  # 替换原有列
  select(-Metabolites) 
B73_cor_results_df <- B73_cor_results_df %>%
  left_join(AA_df1 %>% select(geneID, maizeIDv5), 
            by = c("Gene" = "geneID")) %>%
  mutate(Gene = maizeIDv5) %>%  # 替换原有列
  select(-Gene) 
B73_cor_results_df$maizeIDv5[is.na(B73_cor_results_df$maizeIDv5)] <- "Total" 

#6.5.Merge Data B73 and Teo Data====

merged_results <- B73_cor_results_df %>%
  rename(B73_r = Pearson_r, B73_p = P_value) %>%
  inner_join(
    Teo_cor_results_df %>% rename(Teo_r = Pearson_r, Teo_p = P_value),
    by = c("Cluster", "Metabolite","maizeIDv5")
  ) %>%
  filter(B73_p <= 0.05 | Teo_p <= 0.05) # 至少一个群体显著

merged_results <- merged_results %>% select("Cluster","Metabolite" , "B73_r","B73_p",
                                            "Teo_r"   , "Teo_p" ,"maizeIDv5")

#7.Draw Corraltion Plot----

All_Genes <- unique(merged_results$maizeIDv5)
cluster_colors <- c(
  "Xyl" = "#1f77b4",
  "Phl" = "#ff7f0e",
  "Scl" = "#2ca02c",
  "TB/EP" = "#d62728",
  "LDGT" = "#9467bd",
  "EDGT" = "#8c564b",
  "Hypo" = "#e377c2",
  "Epi" = "#7f7f7f",
  "Unk1" = "#bcbd22",
  "Unk2" = "#17becf",
  "Unk3" = "#aec7e8"
)

# 创建空PPT对象
ppt <- read_pptx()

for (i in seq_along(All_Genes)) {
  GeneID <- All_Genes[i]
  print(paste0(GeneID, "_", i))
  
  Gene_df1 <- merged_results %>% filter(maizeIDv5 == GeneID)
  Gene_df1_clean <- Gene_df1 %>%
    filter(!is.na(Teo_r), !is.na(B73_r))
  
  if (nrow(Gene_df1_clean) == 0) {
    message("❌ Gene ", GeneID, " 在 B73 和 Teo 都没有有效相关性，跳过绘图")
    next
  }
  
  label_data <- Gene_df1_clean %>%
    arrange(pmin(Teo_p, B73_p, na.rm = TRUE)) %>%
    slice_head(n = 10)
  
  p1 <- ggplot(Gene_df1_clean, aes(
    x = Teo_r,
    y = B73_r,
    color = Cluster,
    size = -log10(pmin(Teo_p, B73_p, na.rm = TRUE))
  )) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_text_repel(
      data = label_data,
      aes(x = Teo_r, y = B73_r, label = Metabolite),
      color = "black",
      size = 4,
      box.padding = 0.5,
      max.overlaps = 100,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = cluster_colors) +
    scale_size_continuous(name = "-log10(p)", range = c(2, 6)) +
    theme_classic(base_size = 14) +  # 🎯 去掉背景网格
    theme(
      axis.line = element_line(color = "black", size = 0.8),  # 🎯 加粗坐标轴
      axis.ticks = element_line(color = "black"),             # 🎯 保留刻度线
      axis.text = element_text(color = "black"),              # 🎯 坐标刻度文字黑色
      axis.title = element_text(face = "bold")                # 🎯 坐标标题加粗
    ) +
    labs(
      title = paste("B73 vs Teo 相关性 (Gene:", GeneID, ")"),
      x = "Teo Pearson r", y = "B73 Pearson r"
    )
  
  
  # 把plot转成pptx可以识别的图形对象
  # dml_plot <- dml(ggobj = p1)
  dml_plot <- rvg::dml(ggobj = p1)
  
  # 新建一页 并添加图形
  ppt <- add_slide(ppt, layout = "Title and Content", master = "Office Theme")
  ppt <- ph_with(ppt, value = dml_plot, location = ph_location_fullsize())
  
}

# 保存ppt文件
out_pptx <- file.path(getwd(),"Figures/Correlation_plot/")
if(!dir.exists(out_pptx)){
  dir.create(out_pptx)
}
print(ppt, target = paste0(out_pptx,"/Correlation_plot.pptx"))
message("PPT已保存: ", out_pptx)




