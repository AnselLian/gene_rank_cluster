library(tidyverse)
library(janitor)
library(openxlsx)
library(smoother)
library(dplyr)
library(tidyr)
library(plotrix)

setwd("~/src/gene_count")
cluster_cortex = read_csv('data/assignments Features.csv') %>% clean_names()
cluster_2 = read_csv('data/assignments Features_2.csv') %>% clean_names()
cluster_cerebellum = read_csv('data/assignments_cerebellum.csv') %>% clean_names()
cluster_4 = read_csv('data/de.csv') %>% clean_names()



## filter on differential expression
filter_cortex <- cluster_cortex %>% 
  filter(
         astrocyte_log2_fold_change>0,
         inhibitory_log2_fold_change>0,
         excitatory_log2_fold_change>0,
         oligodendrocyte_log2_fold_change<0,
         microglia_log2_fold_change<0,
         opc_log2_fold_change<0,
         pericyte_log2_fold_change<0,
         endothelial_stalk_log2_fold_change<0,
         ependymal_log2_fold_change<0
         
         
         
         )



filter_2 <- cluster_2 %>% 
  filter(
    astrocyte_log2_fold_change>0,
    inhibitory_log2_fold_change>0,
    excitatory_log2_fold_change>0)
    

filter_cerebellum <- cluster_cerebellum %>% 
  filter(
    purkinje_log2_fold_change>0,
    granule_log2_fold_change>0,
    mli_log2_fold_change>0,
    golgi_log2_fold_change>0,
    pli_log2_fold_change>0,
    ubc_log2_fold_change>0,
    astrocyte_log2_fold_change>0,
    bergmann_glia_log2_fold_change<0,
    oligodendrocyte_log2_fold_change<0,
    opc_log2_fold_change<0,
    microglia_log2_fold_change<0,
    endothelial_stalk_log2_fold_change<0,
    fibroblast_log2_fold_change<0,
    pericyte_log2_fold_change<0
    
  )

cluster_cerebellum$purkinje_rank <- rank(-cluster_cerebellum$purkinje_log2_fold_change)
cluster_cerebellum$granule_rank <- rank(-cluster_cerebellum$granule_log2_fold_change)
cluster_cerebellum$mli_rank <- rank(-cluster_cerebellum$mli_log2_fold_change)
cluster_cerebellum$golgi_rank <- rank(-cluster_cerebellum$golgi_log2_fold_change)
cluster_cerebellum$pli_rank <- rank(-cluster_cerebellum$pli_log2_fold_change)
cluster_cerebellum$ubc_rank <- rank(-cluster_cerebellum$ubc_log2_fold_change)
cluster_cerebellum$astrocyte_rank <- rank(-cluster_cerebellum$astrocyte_log2_fold_change)
cluster_cerebellum$bergmann_glia_rank <- rank(-cluster_cerebellum$bergmann_glia_log2_fold_change)
cluster_cerebellum$oligodendrocyte_rank <- rank(-cluster_cerebellum$oligodendrocyte_log2_fold_change)
cluster_cerebellum$microglia_rank <- rank(-cluster_cerebellum$microglia_log2_fold_change)
cluster_cerebellum$endothelial_stalk_rank <- rank(-cluster_cerebellum$endothelial_stalk_log2_fold_change)
cluster_cerebellum$fibroblast_rank <- rank(-cluster_cerebellum$fibroblast_log2_fold_change)
cluster_cerebellum$pericyte_rank <- rank(-cluster_cerebellum$pericyte_log2_fold_change)




# Aggregate rankings to create a combined score
cluster_cerebellum$combined_rank <- cluster_cerebellum$purkinje_rank +
  cluster_cerebellum$granule_rank +
  cluster_cerebellum$mli_rank +
  cluster_cerebellum$golgi_rank +
  cluster_cerebellum$pli_rank +
  cluster_cerebellum$ubc_rank +
  cluster_cerebellum$astrocyte_rank- 
  cluster_cerebellum$bergmann_glia_rank-
  cluster_cerebellum$oligodendrocyte_rank-
  cluster_cerebellum$microglia_rank-
  cluster_cerebellum$endothelial_stalk_rank-
  cluster_cerebellum$fibroblast_rank-
  cluster_cerebellum$pericyte_rank
  


top_genes <- cluster_cerebellum[order(cluster_cerebellum$combined_rank), ]




## filter on gene count from Terra
gene_count = read_tsv('data/gene_counts_smry.tsv')  %>%
  mutate(pos=rank(celltype))


gene_count_summary <- gene_count %>%
  group_by(gene, celltype) %>%
  summarize(rpm = mean(rpm, na.rm = TRUE), 
            .groups = 'keep') %>%
  pivot_wider(names_from = celltype, values_from = rpm) %>%
  column_to_rownames(var = "gene")


  
neuron = c("Golgi", "MLI", "UBC", "PLI", "Purkinje", "granule", "astrocyte")
other = c("Bergmann glia", "endothelial stalk", "OPC", "fibroblast", "microglia", "oligodendrocyte", "pericyte")

filter_neuron <- gene_count_summary %>%
  filter(if_all(all_of(neuron), ~ . > 100)) %>%
  filter(if_all(all_of(other), ~ . < 100))



plot(gene_count$pos,gene_count$rpm, xlab = "cell_type", ylab = "rpm")









