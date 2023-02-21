##Phylogeny Tree of S. aureus isolates from cutaneous leishmaniasis patients
##Victoria Lovins
##October 2021

##Installing ggtree for visualization of the phylogeny tree ####

#setwd('/Users/vlovins/Desktop/Phylogeny_saureus')
library('ape')
library(gplots)
#library('phyloseq')
#library('magrittr')
#library('circlize')
#library('phylogram')
library('ggtree')
library(patchwork)
#library('RColorBrewer')
#library('ggnewscale')
library(tidyverse) # Added by Camila


##Build the Phylogeny tree using ggtree ####

##Read the Newick file in

#tree <- read.tree('saureus_tree_newick')

# Camila Note: Importing "tree"object from Tori's R environment:
load("tree")
class(tree)

###Remove underscores from tip labels
tree$tip.label <- tree$tip.label %>% str_replace_all('_', ' ')
tree$tip.label <- tree$tip.label %>% str_replace_all('s epidermidis', 'S. epidermidis')

## Rename the tips??
#tree$tip.label %>% c("")

##Build the image tree
pre <- tree %>%
  ggtree(branch.length = "none") +
  geom_tiplab(label=tree$tip.label) +
  ## geom_text(aes(label = node)) +
  theme(legend.position = "none") +
  xlim(0,40)

pre

# From here, Camila is working on the tree ----
# Load color scheme ----
load("../../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/Dark24")

# Ultrametric tree ----
# Read methods for phylogenetic trees: https://pedrohbraga.github.io/PhyloCompMethods-in-R-workshop/PhyloCompMethodsMaterial.html
# To cut tree it's necessary to ultrametric-it. But I chose not to cut it since there are few isolates in each cluster. I did it manually below.
is.ultrametric(tree)
tree_dendrogram <- chronos(tree)
is.ultrametric(tree_dendrogram)

# Convert to hclust object:
tree_hclust <- as.hclust.phylo(tree_dendrogram)

# Effect of chronalization:
#pre
#post

# Cut tree
clinical_module <- cutree(tree_hclust, k=2)
clinical_module_fac <- factor(clinical_module)

# IMPORTANT: The cut tree will be manual, so Clinical Isolates will be analyzed independently of References.

# TREE VIZ ----
'%ni%' <- Negate('%in%')
load("../../Dataset_Integration/outputs/phenotype")
phenotype <- rownames_to_column(phenotype, "label2")
phenotype$label <- str_replace_all(phenotype$label,"^29","SA29")

tree_tb <- as_tibble(tree)
tree_tb <- tree_tb %>%
  mutate(label2 = case_when(
    label %in% ClinicalIsolates ~ "CL S. aureus isolates",
    label %in% "S. epidermidis" ~ "S. epidermidis reference",
    TRUE ~ "S. aureus references"
  )) %>%
  mutate(label = case_when(
    label %in% "SA29417" ~ "CLSA62",
    label %in% "SA29151" ~ "CLSA38",
    label %in% "SA29271" ~ "CLSA19",
    label %in% "SA29150" ~ "CLSA2",
    label %in% "SA29303" ~ "CLSA51",
    label %in% "SA29256" ~ "CLSA16",
    label %in% "SA29137" ~ "CLSA41",
    label %in% "SA29434" ~ "CLSA57",
    label %in% "SA29253" ~ "CLSA14",
    label %in% "SA29297" ~ "CLSA50",
    label %in% "SA29412" ~ "CLSA54",
    label %in% "SA29289" ~ "CLSA48",
    label %in% "SA29258" ~ "CLSA17",
    TRUE ~ label
  )) %>%
  left_join(phenotype)

# For Camila's paper:
as.phylo(tree_tb) %>%
  #ggtree(branch.length = "none") +
  ggtree(branch.length = "none", layout = "fan", open.angle = 190) +
  #ggtree(branch.length = "none", layout = "circular") +
  geom_hilight(node=36, fill="steelblue", alpha=.6) +
  geom_vline(xintercept = 7.5, color=Dark24[7], linetype = 1) +
  geom_tiplab(aes(color=tree_tb$label2)) +
  #geom_strip("20","23", barsize = 2,label = "try", offset.text = .1) +
  theme(legend.position = "right") +
  scale_color_manual(values = c(Dark24[1],"#2B2C39", Dark24[8]), name = 'Key') +
  #annotate("text", x=0, y=0, label="blabal", size=5, color=Dark24[2], fontface="bold") +
  xlim(0,20)






tree_tb

# Export ID of Clinical Isolates ----
ClinicalIsolates
Cluster_order <- c(6,6,6,3,3,3,3,1,2,4,5,5,5) # Not ideal to do by eye

Saureus_clusters <- as.data.frame(cbind(ClinicalIsolates,
                                        Cluster_order))
colnames(Saureus_clusters)[1] <- "label"
#Saureus_clusters <- column_to_rownames(Saureus_clusters, "LTCP_patient_ID")
#rownames(Saureus_clusters) <- str_replace_all(rownames(Saureus_clusters), "SA","CL")
#save(Saureus_clusters, file = "Saureus_clusters")

# Associating Clinical Isolates with Clinical metadata and other features ----
# Import parasite load
load("../../3rd_LesionRNAseq_Dataset/qPCRs/Lbraziliensis_qPCR_ImportingDataIntoR/parasiteload_3rdDataset")
parasiteload_3rdDataset <- parasiteload_3rdDataset %>%
  select(sample_RNAseq, Lbraziliensis_mg) %>%
  rename(label = sample_RNAseq)
parasiteload_3rdDataset$label <- str_replace_all(parasiteload_3rdDataset$label, "^CL","SA")

# Import S. aureus CPM
load("../../3rd_LesionRNAseq_Dataset/mapping_Saureus_pangenome/outputs/targets_info")
targets_info <- targets_info %>%
  select(LTCP_patient_ID, cpm_saureus) %>%
  rename(label = LTCP_patient_ID)
targets_info$label <- str_replace_all(targets_info$label, "^29","SA29")

# Import 16S qPCR from biopsies (bacterial burden)
load("../../3rd_LesionRNAseq_Dataset/qPCRs/Saureus16S_qPCR_ImportingDataIntoR/qPCR_23")
qPCR_23 <- qPCR_23 %>%
  select(LTCP_patient_ID, s16_mg) %>%
  rename(label = LTCP_patient_ID)
qPCR_23$label <- str_replace_all(qPCR_23$label, "^29","SA29")

clinc <- phenotype %>%
  filter(label %in% ClinicalIsolates) %>%
  left_join(Saureus_clusters) %>%
  left_join(qPCR_23) %>%
  left_join(targets_info) %>%
  left_join(parasiteload_3rdDataset) %>%
  arrange(Cluster_order) %>%
  select(Cluster_order, label, age:Clinical_Outcome, microbiome_cluster, s16_mg, cpm_saureus, Lbraziliensis_mg) %>%
  rename(Cluster = Cluster_order,
         "CL S. aureus" = label, "Bacterial load" = s16_mg, "S. aureus (CPM)" = cpm_saureus, "Parasite load" = Lbraziliensis_mg,
         Age = age, Sex = sex, "Lesion size (mm2)" = size_lesion_mm2, Lymphadenopathy = lymphadenopathy,
         "Healing time (days)" = healing_time_days, "Clinical Outcome" = Clinical_Outcome, "Main dysbiosis" = microbiome_cluster)
clinc
#write_tsv(clinc, "../../../../../Desktop/clinical.txt")


# Associating Clinical Isolates with Variable up gene expression (Inflammation) heatmap ----
load("../../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/log2.cpm.filtered.norm")
load("../outputs/orderPC2_v3_pos") # From data dimensionality reduction
load("../../3rd_LesionRNAseq_Dataset/mapping_Saureus_pangenome/outputs/sig_genes_up8") # up DGE S. aureus detect (from pangenome mapping)
load("../../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/myTopHits_treat_var_imp") # Biomarkers for treatment outcome

#micro_colors <- function(micro_module_mean_fac){ if (micro_module_mean_fac==1) Dark24[1] else if(micro_module_mean_fac==2) Dark24[2] else if (micro_module_mean_fac==3) Dark24[3] else if (micro_module_mean_fac==4) "gray" else if (micro_module_mean_fac==5) Dark24[4]}
#micro_colors <- unlist(lapply(micro_module_mean_fac, micro_colors))

heatmap.2(log2.cpm.filtered.norm[myTopHits_treat_var_imp,
                                 colnames(log2.cpm.filtered.norm) %in% str_replace_all(clinc$`CL S. aureus`,"SA","CL")],
          #Colv=as.dendrogram(microclust),
          #Rowv = as.dendrogram(impgenes),
          dendrogram = "both",
          #RowSideColors=module.color,
          col=colorRampPalette(colors=c("blue","white","red"))(100),
          #ColSideColors = micro_colors_mean,
          scale='row',
          density.info="none", trace="none",
          cexRow=1.1, cexCol=1.1,margins = c(5,15))

mat2 <- log2.cpm.filtered.norm[myTopHits_treat_var_imp,
                       colnames(log2.cpm.filtered.norm) %in% str_replace_all(clinc$`CL S. aureus`,"SA","CL")]
heatmap.2(mat2,
          #Colv=as.dendrogram(microclust),
          #Rowv = as.dendrogram(impgenes),
          #Colv = F,
          dendrogram = "both",
          #RowSideColors=module.color,
          col=colorRampPalette(colors=c("blue","white","red"))(100),
          #ColSideColors = micro_colors_mean,
          scale='row',
          density.info="none", trace="none",
          cexRow=1.1, cexCol=1.1,margins = c(5,15))


class(log2.cpm.filtered.norm)

# ssGSEA for Failure gene signature ----
library(GSVA)
library(GSEABase)
library(ggforce)
library(ggpubr)
library(ggrepel)
Failure_geneSig <- getGmt("../../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/myTopHits_treat_var_imp.gmt", geneIdType=SymbolIdentifier())

# GSVA Amorim ----
GSVA.res_Failure <- gsva(log2.cpm.filtered.norm, #your data
                         Failure_geneSig, # your gene set collection
                         #min.sz=5, max.sz=500, #criteria for filtering gene sets
                         mx.diff=FALSE,
                         method="ssgsea")
my_comparisons <- list(c("CL Cure", "CL Failure"))

FailureES_df <- rownames_to_column(as.data.frame(GSVA.res_Failure),"feature") %>%
  pivot_longer(!feature, names_to = "label", values_to = "ES") %>%
  dplyr::select(-feature) %>%
  mutate(label = str_replace_all(label,"CL","SA")) %>%
  left_join(phenotype) %>%
  mutate(withBiopsy = case_when(
    label %in% ClinicalIsolates ~ "Clinical_Isolates",
    TRUE ~ "other")) %>%
  filter(Clinical_Outcome != "NA" | label %in% str_subset(label,"^HS")) %>%
  mutate(group = case_when(
    label %in% str_subset(label,"^SA") & Clinical_Outcome == "Cure" ~ "CL Cure",
    label %in% str_subset(label,"^SA") & Clinical_Outcome == "Failure" ~ "CL Failure",
    TRUE ~ "HS")) %>%
  mutate(group2 = case_when(
    group %in% c("CL Cure","CL Failure") ~ "CL",
    TRUE ~ "HS")) %>%
  mutate(group = factor(group, levels = c("HS","CL Cure","CL Failure"))) %>%
  mutate(group2 = factor(group2, levels = c("HS","CL"))) %>%
  left_join(Saureus_clusters %>% mutate(Clinical_IsolateID = label))

FailureES_df %>%
  ggplot(.,aes(x=group, y=ES, fill=group)) +
  #geom_boxplot(fill=NA) +
  geom_jitter(position = position_jitter(0.1),size=3, shape=21) + 
  #geom_violin(size = 0.7, trim = F) +
  #geom_sina(size=2.3) +
  theme_classic() + scale_fill_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="right") +
  #geom_text_repel(aes(label = Clinical_IsolateID), size = 3, fontface="bold", nudge_x = 0.5) +
  stat_summary(fun = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test", label.y = .50,
                     aes(label = ..p.signif..),
                     #aes(label = ..p.value..),
                     #label.x = 1.5,
                     size = 5, color="black") +
  xlab("") + ylab(paste("Failure gene signature ES")) +
  ylim(-.60,.60)

FailureES_df %>%
  filter(group != "HS") %>%
  ggplot(.,aes(x=group2, y=ES)) +
    geom_point(aes(fill=group), size=3, shape=21) + 
    geom_boxplot(fill=NA, width=0.3, color=Dark24[20], size=0.5) +
    annotate("text", x=0.63, y=0.4, label="Percentile:", color="black", size=3, fontface=2) +
    annotate("text", x=0.63, y=0.3, label=">75th", color=Dark24[1], size=3) +
    annotate("text", x=0.63, y=0.23, label="50-75th", color="#CF7487", size=3) +
    annotate("text", x=0.63, y=0.15, label="25-50th", color="dark gray", size=3) +
    annotate("text", x=0.63, y=0.05, label="0-25th", color="gray", size=3) +
    #geom_violin(size = 0.7, trim = F) +
    #geom_sina(size=2.3) +
    theme_classic() + scale_fill_manual(values = Dark24[-1]) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 13), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
          legend.position="right") +
    geom_text_repel(aes(label = Clinical_IsolateID), size = 3, fontface="bold", nudge_x = 0.5) +
    xlab("") + ylab(paste("Failure gene signature ES")) +
    ylim(-.50,.40)
















