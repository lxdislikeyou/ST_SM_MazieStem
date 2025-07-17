##########################################################################################################################################################
#### Metabolites (B73 vs TEO)--> cluster
###########################################################################################################################################################

SM_BT <- readRDS(file.path(getwd(),"Git_data/BT_combined_MEM_objects.rds"))

B73_SM_neg <- SM_BT[[1]]
B73_SM_pos <- SM_BT[[2]]
Teo_SM_neg <- SM_BT[[3]]
Teo_SM_pos <- SM_BT[[4]]

# library(dplyr)
library(tibble)   # for rownames_to_column()

# B73_Ex      <- B73_Ex      %>% rownames_to_column(var = "label")
B73_SM_neg  <- B73_SM_neg  %>% rownames_to_column(var = "label")
B73_SM_pos  <- B73_SM_pos  %>% rownames_to_column(var = "label")
Teo_SM_neg  <- Teo_SM_neg  %>% rownames_to_column(var = "label")
Teo_SM_pos  <- Teo_SM_pos  %>% rownames_to_column(var = "label")

# 加载必要包
library(dplyr)
library(purrr)
library(tidyr)
library(rstatix)
library(broom)   
library(openxlsx)

#Read AA list

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


wide_to_long <- function(Sm_df,mzid){
  Sm_df <- Sm_df %>% select(Sample,Cluster,mzid)
  
  Sm_df_long <- Sm_df %>% pivot_longer(
    cols = where(is.numeric),
    names_to = "mz",
    values_to = "Value"
  )
  return(Sm_df_long)
}

B73_neg_long <- wide_to_long(B73_SM_neg,neg_mz)
B73_pos_long <- wide_to_long(B73_SM_pos,pos_mz)
Teo_neg_long <- wide_to_long(Teo_SM_neg,neg_mz)
Teo_pos_long <- wide_to_long(Teo_SM_pos,pos_mz)

B73_sm_long <- rbind(B73_neg_long,B73_pos_long)
Teo_sm_long <- rbind(Teo_neg_long,Teo_pos_long)

long_to_wide <- function(long_df) {
  wide_df <- long_df %>% 
    pivot_wider(
      names_from = "mz",     
      values_from = "Value"   
    )
  return(wide_df)
}

B73_sm_wide <- long_to_wide(B73_sm_long) %>% unnest(cols = everything())
Teo_sm_wide <- long_to_wide(Teo_sm_long) %>% unnest(cols = everything())



Pvalue_function <- function(df){
  metabolites <- names(df)[3:ncol(df)]
  
  results <- map_dfr(metabolites, function(met) {
    fml <- as.formula(sprintf("`%s` ~ Cluster", met))
    
    kr <- kruskal_test(df, formula = fml)
    
    means <- df %>% 
      group_by(Cluster) %>% 
      summarise(mean_val = mean(.data[[met]], na.rm = TRUE), .groups = "drop")
    top_cluster <- means %>% 
      slice_max(mean_val, n = 1) %>% 
      pull(Cluster) %>% 
      as.character()
    
    pw <- pairwise_wilcox_test(df, formula = fml, p.adjust.method = "BH")
    
    sig_pairs <- pw %>% 
      filter(p.adj < 0.05) %>% 
      transmute(pair = paste(group1, "vs", group2)) %>% 
      pull(pair) %>% 
      { if(length(.)>0) paste(., collapse = "; ") else NA_character_ }
    
    tibble(
      metabolite      = met,
      kruskal_p       = kr$p,                       # Kruskal–Wallis p
      kruskal_p.adj   = p.adjust(kr$p, method = "BH"),  # BH 校正
      top_cluster     = top_cluster,
      sig_comparisons = sig_pairs
    )
  })
  
}

B73_P <- Pvalue_function(B73_sm_wide)
Teo_P <- Pvalue_function(Teo_sm_wide)

