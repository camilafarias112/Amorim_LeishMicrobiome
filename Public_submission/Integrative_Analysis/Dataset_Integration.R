# In this script, I analyze Integrative data, from different workflows.
# I have used rexposome and Maaslin for dataset integration.
# The inputs for those pipelines were reduced datasets:
# 1) A list of reduced genes from RNA-seq, and script can be found: /Users/amorimc/Dropbox/CamilaAnalysis/Lesion3rd_dataset_Microbiome/3rd_LesionRNAseq_Dataset/Unsup_geneExpr_DataReduction
# 2) Microbiome clusters reduced from 16S-seq, and some diversity metrics included. Find in: /Users/amorimc/Dropbox/CamilaAnalysis/Lesion3rd_dataset_Microbiome/microbiome_Lesion16S-seqDataset
# 3) Clinical metadata.

# Here I also have the inventory information for this study.

# Libraries ----
library(gplots)
library(gt)
library(ggdensity)
library(ggrepel)
library(ggpubr)
library(ggforce)
library(reshape)
library(tidyverse)
library(RColorBrewer)


# Color scheme ----
load("../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/Dark24")

# Overview of all datasets and number of samples used for analysis INVENTORY ----
load("../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/targets_all") # bring main study design
LeishOmics <- targets_all %>%
  select(LTCP_patient_ID,
         matching_biopsies_for_RNAseq,
         age,
         sex,
         illness_duration_days, 
         local_of_lesion,
         size_lesion_mm2,
         n_lesions,
         lymphadenopathy,
         DTH_mm2,
         treatment_outcome,
         healing_time_days,
         treatment_other_drug,
         group_RNAseq
         ) 

load("../microbiome_Lesion16S-seqDataset/Qiime2/metadata_mod")
metadata_mod <- rownames_to_column(metadata_mod, "LTCP_patient_ID")
metadata_mod <- metadata_mod[,-11]
LeishOmics <- LeishOmics %>%
  left_join(metadata_mod)

#write_tsv(LeishOmics %>%
#            select(LTCP_patient_ID, age, sex, size_lesion_mm2, n_lesions,
#                   lymphadenopathy, DTH_mm2, treatment_outcome, healing_time_days) %>%
#            filter(LTCP_patient_ID %in% c("29271","29213","29256","29272","29513",
#                                          "29159","29231","29151")), "../../../../Desktop/LeishOmics_Arca.txt")


LeishOmics_heat <- LeishOmics %>%
  filter(LTCP_patient_ID != "NA") %>%
  mutate(day0_biopsy = case_when(
    matching_biopsies_for_RNAseq == "yes" ~ 1,
    matching_biopsies_for_RNAseq == "no" ~ 0)) %>%
  select(LTCP_patient_ID,day0_biopsy,day0,day30,day60,day90,day120,day150,day180,day210,day240) %>%
  arrange(desc(day0_biopsy), desc(day0)) %>%
  replace(is.na(.), 0) %>%
  filter(day0_biopsy == 1 | day0 == 1) %>%
  column_to_rownames("LTCP_patient_ID")

# matrix:
colorheat <- colorRampPalette(colors=c("white",Dark24[1]))(2)
heatInventory <- heatmap.2(as.matrix(LeishOmics_heat), 
          dendrogram = "none",
          Colv=F,
          col=colorheat, key = F,
          density.info="none", trace="none",  
          cexCol=1, cexRow = .7,
          colsep=1:nrow(LeishOmics_heat),
          rowsep=1:nrow(LeishOmics_heat),
          sepcolor = "black", margins = c(0.1,10),
          NULL)

heat_sampleorder <- rev(rownames(as.matrix(LeishOmics_heat))[heatInventory$rowInd])

LeishOmics_gt <- LeishOmics %>%
  filter(LTCP_patient_ID %in% heat_sampleorder) %>%
  column_to_rownames("LTCP_patient_ID")
LeishOmics_gt <- LeishOmics_gt[heat_sampleorder,] %>%
  rownames_to_column("LTCP_patient_ID")

LeishOmics_gt <- LeishOmics_gt %>%
  select(LTCP_patient_ID,
         age, sex, size_lesion_mm2, lymphadenopathy, DTH_mm2, treatment_outcome,
         healing_time_days, treatment_other_drug) %>%
  mutate(treatment_other_drug = case_when(
    treatment_other_drug == "sbv" ~ "antimony",
    treatment_other_drug == "miltefosine" ~ "DOT")) %>%
  mutate(treatment_outcome = case_when(
    treatment_outcome %in% "failure" ~ "Failure",
    treatment_outcome %in% "cure" ~ "Cure",
    treatment_outcome %in% c("cure_miltefosine","failure_miltefosine") ~ "DOT"))
  
LeishOmics_gt <- rownames_to_column(LeishOmics_gt, "PatientID")
LeishOmics_gt$PatientID <- str_replace_all(LeishOmics_gt$PatientID, "^", "Patient")
write_tsv(LeishOmics_gt %>%
            select(-LTCP_patient_ID), "outputs/LeishOmics_inventoryV2.txt")

phenotype_quick <- phenotype
phenotype_quick$sample_RNAseq <- str_replace_all(phenotype_quick$sample_RNAseq, "CL29", "29")
phenotype_quick <- phenotype_quick %>% mutate(LTCP_patient_ID = sample_RNAseq) %>%
  select(LTCP_patient_ID, microbiome_cluster)
                             
# SUPPLEMENTAL TABLE 1 ----              
LeishOmics_gt3 <- LeishOmics_gt %>%
  left_join(phenotype_quick) 
write_tsv(LeishOmics_gt3, "../../../myPapers_submissions/Microbiome/SupplementalMaterial/SupplemmentalTable1_LeishOmics_inventory.txt")

# Add total bacterial load
load("../3rd_LesionRNAseq_Dataset/qPCRs/Saureus16S_qPCR_ImportingDataIntoR/qPCR_23")
LeishOmics_gt3 <- LeishOmics_gt3 %>%
  left_join(qPCR_23 %>% select(LTCP_patient_ID, s16_mg))
LeishOmics_gt3 <- LeishOmics_gt3 %>% rename("Bacterial load (qPCR)" = s16_mg)

# Add S. aureus transcript abundances
load("../3rd_LesionRNAseq_Dataset/mapping_Saureus_pangenome/outputs/targets_CL")
LeishOmics_gt3 <- LeishOmics_gt3 %>%
  left_join(targets_CL %>% select(LTCP_patient_ID, cpm_saureus, cpm_saureus_cat))
LeishOmics_gt3 <- LeishOmics_gt3 %>% rename("S. aureus abundance (CPM)" = cpm_saureus,
                                            "S. aureus abundance - categorical" = cpm_saureus_cat)

# Add Parasite Loads
load("../3rd_LesionRNAseq_Dataset/qPCRs/Lbraziliensis_qPCR_ImportingDataIntoR/parasiteload_3rdDataset")
parasiteload_3rdDataset$sample_RNAseq <- str_replace_all(parasiteload_3rdDataset$sample_RNAseq, "^CL", "")
parasiteload_3rdDataset <- parasiteload_3rdDataset %>% rename(LTCP_patient_ID = sample_RNAseq)
LeishOmics_gt3 <- LeishOmics_gt3 %>%
  left_join(parasiteload_3rdDataset %>% select(LTCP_patient_ID, Lbraziliensis_mg))
LeishOmics_gt3 <- LeishOmics_gt3 %>% rename("L.braziliensis load (qPCR)" = Lbraziliensis_mg)

write_tsv(LeishOmics_gt3, "../../../myPapers_submissions/Microbiome/SupplementalMaterial/SupplemmentalTable1_LeishOmics_inventory.txt")

# Change name of the clusters:
load("../3rd_LesionRNAseq_Dataset/qPCRs/Lbraziliensis_qPCR_ImportingDataIntoR/parasiteload_3rdDataset")

LeishOmics_gt3 <- LeishOmics_gt3 %>%
  mutate(microbiome_cluster = case_when(
    microbiome_cluster == "Staphylococcus" ~ "M6",
    microbiome_cluster == "Arcanobacterium" ~ "M7",
    microbiome_cluster == "Staphylococcus_Streptococcus" ~ "M3",
    microbiome_cluster == "Streptococcus" ~ "M4",
    microbiome_cluster == "Heterogeneous" ~ "M0",
    microbiome_cluster == "Corynebacterium" ~ "M5",
    microbiome_cluster == "Finegoldia" ~ "M2",
    microbiome_cluster == "Lactobacillus" ~ "M1")) 

write_tsv(LeishOmics_gt3, "../../../myPapers_submissions/Microbiome/SupplementalMaterial/SupplemmentalTable1_LeishOmics_inventory.txt")

# Export for GEO:
write_tsv(LeishOmics_gt3, "../GEO_submission/Amorim2022_SupplemmentalTable1_LeishOmics_StudyDesign.txt")

# Continuation of the previous chunk ----
biopsy_samples <- LeishOmics_gt %>%
  filter(matching_biopsies_for_RNAseq == "yes")
biopsy_samples <- biopsy_samples$LTCP_patient_ID

swab0_samples <- LeishOmics_gt %>%
  filter(day0 == 1)
swab0_samples <- swab0_samples$LTCP_patient_ID

swabfollow_samples <- LeishOmics_gt %>%
  filter(day0 == 1) %>%
  filter(day60 == 1 | day90 == 1)
swabfollow_samples <- swabfollow_samples$LTCP_patient_ID

venn_samples1 <- venn(list(Biopsy = biopsy_samples,
                          "Swab (day 0)" = swab0_samples))

venn_samples2 <- venn(list("Swab (day 0)" = swab0_samples,
                           "Swab follow-up" = swabfollow_samples))

# Create the phenotype data for rexposome ----
load("../microbiome_Lesion16S-seqDataset/outputs/phenotype_pre")

#'%ni%' <- Negate('%in%')
phenotype <- LeishOmics_gt %>%
  left_join(phenotype_pre) %>%
  mutate(Clinical_Outcome = case_when(
    treatment_outcome %in% "Cure" ~ "1 Sbv",
    treatment_outcome %in% "Failure" ~ ">1 Sbv",
    TRUE ~ "DOT")) %>%
  #mutate_at(vars(Clinical_Outcome), na_if, "miltefosine") %>%
  select(LTCP_patient_ID,
         microbiome_cluster,
         shannon_entropy,
         observed_features,
         age,
         sex,
         size_lesion_mm2,
         lymphadenopathy,
         healing_time_days,
         Clinical_Outcome) %>%
  column_to_rownames("LTCP_patient_ID")

healing_v <- phenotype$Clinical_Outcome %in% NA
phenotype$healing_time_days[healing_v] <- NA
save(phenotype, file = "outputs/phenotype")

phenotype <- rownames_to_column(phenotype, "sample_RNAseq")
phenotype$sample_RNAseq <- str_replace_all(phenotype$sample_RNAseq, "^29", "CL29")

phenotype$microbiome_cluster <- factor(phenotype$microbiome_cluster,
                                       levels = c("Staphylococcus",
                                                  "Arcanobacterium",
                                                  "Corynebacterium",
                                                  "Streptococcus",
                                                  "Staphylococcus_Streptococcus",
                                                  "Heterogeneous",
                                                  "Lactobacillus",
                                                  "Finegoldia"))


# Imported rexposome outputs (ONLY UP VARIABLE GENES) ----
# Samples with missing info in the phenotype file:
load("inputs_rexposome_maaslin/missing_info")
missing_info

#PCA samples
load("inputs_rexposome_maaslin/pcaplot_microbiome_up")
pcaplot_microbiome_up[["data"]][["Label"]] <- str_replace_all(pcaplot_microbiome_up[["data"]][["Label"]], "^29", "CL29")

pcaplot_microbiome_up[["data"]] %>% 
  rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  mutate(microbiome_cluster2 = case_when(
    microbiome_cluster == "Staphylococcus" ~ "M6",
    microbiome_cluster == "Arcanobacterium" ~ "M7",
    microbiome_cluster == "Streptococcus" ~ "M4",
    microbiome_cluster == "Heterogeneous" ~ "M0",
    microbiome_cluster == "Corynebacterium" ~ "M5",
    microbiome_cluster == "Finegoldia" ~ "M2",
    microbiome_cluster == "Lactobacillus" ~ "M1",
    TRUE ~ "HS")) %>%
  ggplot(., aes(x=Dim.1, y=Dim.2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 4) +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),legend.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17)) +
  scale_color_manual(values = c("gray",
                                Dark24[5],
                                Dark24[6],
                                Dark24[3],
                                Dark24[2],
                                Dark24[1],
                                Dark24[4])) +
  #geom_text_repel(aes(label = phenotype), size = 2, fontface="bold",
  #                color="black", max.overlaps = 20) +
  #stat_ellipse(level = 0.95) +
  xlab(pcaplot_loadings_up[["labels"]][["x"]]) + ylab(pcaplot_loadings_up[["labels"]][["y"]]) + 
  coord_fixed()

pcaplot_microbiome_up[["data"]] %>% 
  rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  mutate(microbiome_cluster2 = case_when(
    microbiome_cluster == "Staphylococcus" ~ "M6",
    microbiome_cluster == "Arcanobacterium" ~ "M7",
    microbiome_cluster == "Streptococcus" ~ "M4",
    microbiome_cluster == "Heterogeneous" ~ "M0",
    microbiome_cluster == "Corynebacterium" ~ "M5",
    microbiome_cluster == "Finegoldia" ~ "M2",
    microbiome_cluster == "Lactobacillus" ~ "M1",
    TRUE ~ "HS")) %>%
ggplot(., aes(x=Dim.1, y=Dim.2, color=microbiome_cluster2)) +
  geom_mark_hull(aes(fill = microbiome_cluster2), alpha = 0.1, show.legend = F, expand = unit(2.5, "mm")) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 4) +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),legend.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17)) +
  scale_color_manual(values = c("gray",
                                  Dark24[5],
                                  Dark24[6],
                                  Dark24[3],
                                  Dark24[2],
                                  Dark24[1],
                                  Dark24[4])) +
  #geom_text_repel(aes(label = phenotype), size = 2, fontface="bold",
  #                color="black", max.overlaps = 20) +
  #stat_ellipse(level = 0.95) +
  xlab(pcaplot_loadings_up[["labels"]][["x"]]) + ylab(pcaplot_loadings_up[["labels"]][["y"]]) + 
  coord_fixed()

pcaplot_microbiome_up[["data"]] %>% 
  rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  mutate(microbiome_cluster2 = case_when(
    microbiome_cluster == "Staphylococcus" ~ "M6",
    microbiome_cluster == "Arcanobacterium" ~ "M7",
    microbiome_cluster == "Streptococcus" ~ "M4",
    microbiome_cluster == "Heterogeneous" ~ "M0",
    microbiome_cluster == "Corynebacterium" ~ "M5",
    microbiome_cluster == "Finegoldia" ~ "M2",
    microbiome_cluster == "Lactobacillus" ~ "M1",
    TRUE ~ "HS")) %>%
  ggplot(., aes(x=microbiome_cluster2, y=Dim.2, color=microbiome_cluster2)) +
  geom_hline(yintercept = 0) +
  geom_violin(size = 0.7, trim = T) +
  geom_point(size = 2.5) +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17),
        axis.text.x = element_blank(), axis.ticks = element_blank()) +
  scale_color_manual(values = c("gray",
                                Dark24[5],
                                Dark24[6],
                                Dark24[3],
                                Dark24[2],
                                Dark24[1],
                                Dark24[4])) +
  ylab(pcaplot_loadings_up[["labels"]][["y"]]) + xlab("")

pcaplot_microbiome_up[["data"]] %>% 
  rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  mutate(microbiome_cluster2 = case_when(
    microbiome_cluster == "Staphylococcus" ~ "M6",
    microbiome_cluster == "Arcanobacterium" ~ "M7",
    microbiome_cluster == "Streptococcus" ~ "M4",
    microbiome_cluster == "Heterogeneous" ~ "M0",
    microbiome_cluster == "Corynebacterium" ~ "M5",
    microbiome_cluster == "Finegoldia" ~ "M2",
    microbiome_cluster == "Lactobacillus" ~ "M1",
    TRUE ~ "HS")) %>%
  filter(Clinical_Outcome != "DOT") %>%
  ggplot(., aes(x=observed_features, y=Dim.2, color=observed_features)) +
  geom_hline(yintercept = 0) +
  geom_smooth(method=lm, color="#000080") +
  geom_point(size = 2.5, color="black") +
  theme_classic() +
  scale_color_gradient(high = "#750D86", low = "gray") +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17)) +
  ylab(pcaplot_loadings_up[["labels"]][["y"]]) + xlab("# observed OTUs") +
  ylim(-15,15)

pcaplot_microbiome_up[["data"]] %>% 
  rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  mutate(microbiome_cluster2 = case_when(
    microbiome_cluster == "Staphylococcus" ~ "M6",
    microbiome_cluster == "Arcanobacterium" ~ "M7",
    microbiome_cluster == "Streptococcus" ~ "M4",
    microbiome_cluster == "Heterogeneous" ~ "M0",
    microbiome_cluster == "Corynebacterium" ~ "M5",
    microbiome_cluster == "Finegoldia" ~ "M2",
    microbiome_cluster == "Lactobacillus" ~ "M1",
    TRUE ~ "HS")) %>%
  filter(Clinical_Outcome != "DOT") %>%
  ggplot(., aes(x=fct_rev(Clinical_Outcome), y=Dim.2, color=Clinical_Outcome)) +
  geom_hline(yintercept = 0) +
  geom_violin(size = 0.7, trim = T) +
  geom_point(size = 2.5) +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17)) +
  scale_color_manual(values = c("#000080",
                                Dark24[9])) +
  ylab(pcaplot_loadings_up[["labels"]][["y"]]) + xlab("")

pcaplot_microbiome_up[["data"]] %>% 
  rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  mutate(microbiome_cluster2 = case_when(
    microbiome_cluster == "Staphylococcus" ~ "M6",
    microbiome_cluster == "Arcanobacterium" ~ "M7",
    microbiome_cluster == "Streptococcus" ~ "M4",
    microbiome_cluster == "Heterogeneous" ~ "M0",
    microbiome_cluster == "Corynebacterium" ~ "M5",
    microbiome_cluster == "Finegoldia" ~ "M2",
    microbiome_cluster == "Lactobacillus" ~ "M1",
    TRUE ~ "HS")) %>%
  filter(Clinical_Outcome != "DOT") %>%
  ggplot(., aes(x=healing_time_days, y=Dim.2, color=healing_time_days)) +
  geom_hline(yintercept = 0) +
  geom_smooth(method=lm, color="#000080") +
  geom_point(size = 2.5, color="black") +
  theme_classic() +
  scale_color_gradient(high = "#750D86", low = "gray") +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17)) +
  ylab(pcaplot_loadings_up[["labels"]][["y"]]) + xlab("Healing time (days)")


orderPC2s <- pcaplot_microbiome_up[["data"]] %>%
  arrange(Dim.2)
orderPC2s <- orderPC2s$Label

pcaplot_microbiome_up[["data"]] %>%
  arrange(Dim.2) %>%
  dplyr::rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  ggplot(., aes(x=factor(sample_RNAseq, levels = rev(orderPC2s)), y=Dim.2, fill=phenotype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(Dark24[4],Dark24[2],"gray",Dark24[12],Dark24[1],Dark24[3])) +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  #geom_text(aes(y=0, label=Label), angle = 90,
  #          fontface = 2, size=1, vjust = 0.5, hjust=1.1) +
  geom_hline(yintercept = 0) +
  ylab(pcaplot_microbiome_up[["labels"]][["y"]])

pcaplot_microbiome_up[["data"]] %>%
  arrange(Dim.2) %>%
  dplyr::rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  ggplot(., aes(x=factor(sample_RNAseq, levels = rev(orderPC2s)), y=Dim.2, fill=Clinical_Outcome)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(Dark24[2],Dark24[1])) +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  #geom_text(aes(y=0, label=Label), angle = 90,
  #          fontface = 2, size=1, vjust = 0.5, hjust=1.1) +
  geom_hline(yintercept = 0) +
  ylab(pcaplot_microbiome_up[["labels"]][["y"]])

pcaplot_microbiome_up[["data"]] %>%
  arrange(Dim.2) %>%
  dplyr::rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  ggplot(., aes(x=factor(sample_RNAseq, levels = rev(orderPC2s)), y=Dim.2, fill=healing_time_days)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "gray",
                      high = "purple",
                      space = "Lab",
                      na.value = "white") +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        #axis.text.x = element_blank(), 
        axis.text.x = element_text(angle = 90),
        axis.ticks.x = element_blank()) +
  #geom_text(aes(y=0, label=Label), angle = 90,
  #          fontface = 2, size=1, vjust = 0.5, hjust=1.1) +
  geom_hline(yintercept = 0) +
  ylab(pcaplot_microbiome_up[["labels"]][["y"]])

pcaplot_microbiome_up[["data"]] %>%
  arrange(Dim.2) %>%
  dplyr::rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  ggplot(., aes(x=factor(sample_RNAseq, levels = rev(orderPC2s)), y=Dim.2, fill=observed_features)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "gray",
    high = "purple",
    space = "Lab",
    na.value = "white") +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  #geom_text(aes(y=0, label=Label), angle = 90,
  #          fontface = 2, size=1, vjust = 0.5, hjust=1.1) +
  geom_hline(yintercept = 0) +
  ylab(pcaplot_microbiome_up[["labels"]][["y"]])

pcaplot_microbiome_up[["data"]] %>%
  arrange(Dim.2) %>%
  dplyr::rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  ggplot(., aes(x=factor(sample_RNAseq, levels = rev(orderPC2s)), y=Dim.2, fill=shannon_entropy)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "gray",
                      high = "purple",
                      space = "Lab",
                      na.value = "white") +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  #geom_text(aes(y=0, label=Label), angle = 90,
  #          fontface = 2, size=1, vjust = 0.5, hjust=1.1) +
  geom_hline(yintercept = 0) +
  ylab(pcaplot_microbiome_up[["labels"]][["y"]])

pcaplot_microbiome_up[["data"]] %>%
  arrange(Dim.2) %>%
  dplyr::rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  ggplot(., aes(x=factor(sample_RNAseq, levels = rev(orderPC2s)), y=Dim.2, fill=size_lesion_mm2)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "gray",
                      high = "purple",
                      space = "Lab",
                      na.value = "white") +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  #geom_text(aes(y=0, label=Label), angle = 90,
  #          fontface = 2, size=1, vjust = 0.5, hjust=1.1) +
  geom_hline(yintercept = 0) +
  ylab(pcaplot_microbiome_up[["labels"]][["y"]])

pcaplot_microbiome_up[["data"]] %>%
  arrange(Dim.2) %>%
  dplyr::rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  ggplot(., aes(x=factor(sample_RNAseq, levels = rev(orderPC2s)), y=Dim.2, fill=age)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "gray",
                      high = "purple",
                      space = "Lab",
                      na.value = "white") +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  #geom_text(aes(y=0, label=Label), angle = 90,
  #          fontface = 2, size=1, vjust = 0.5, hjust=1.1) +
  geom_hline(yintercept = 0) +
  ylab(pcaplot_microbiome_up[["labels"]][["y"]])

#PCA PHE:
load("inputs_rexposome_maaslin/plotPHE_up")

plotPHE_up
  
plotPHE_up[["data"]] %>%
  filter(Dim == "PC 02")

plotPHE_up[["data"]] %>%
  filter(Dim == "PC 09")

#PCA loadings:
load("inputs_rexposome_maaslin/pcaplot_loadings_up")

pcaplot_loadings_up[["data"]] %>%
  ggplot(., aes(x=Dim.1, y=Dim.2)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #stat_ellipse(level = 0.95) +
  xlab(pcaplot_loadings_up[["labels"]][["x"]]) + ylab(pcaplot_loadings_up[["labels"]][["y"]]) + 
  coord_fixed()

pcaplot_loadings_up[["data"]] %>%
  ggplot(., aes(x=Dim.1, y=Dim.2)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        panel.background = element_rect_wiggle(sides = c("lb"),
                                               colour = "black")) +
  geom_hline(yintercept = 0, color = Dark24[1], size=2) +
  #geom_vline(xintercept = 0, color = Dark24[4]) +
  #stat_ellipse(level = 0.95) +
  xlab(pcaplot_loadings_up[["labels"]][["x"]]) + ylab(pcaplot_loadings_up[["labels"]][["y"]]) + 
  coord_fixed()

pcaplot_loadings_up[["data"]] %>%
  mutate(Label_IGs = case_when(
    startsWith(Label, 'IG') ~ "Immunoglobulin encoding genes",
    TRUE ~ "other")) %>%
  ggplot(., aes(x=Dim.1, y=Dim.2, color=Label_IGs)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  scale_color_manual(values = c(Dark24[18],"black")) +
  geom_text_repel(aes(label = Label), size = 2.5, fontface="bold", #color="black",
                  max.overlaps = 50) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #stat_ellipse(level = 0.95) +
  xlab(pcaplot_loadings_up[["labels"]][["x"]]) + ylab(pcaplot_loadings_up[["labels"]][["y"]]) + 
  coord_fixed()

orderPC2 <- pcaplot_loadings_up[["data"]] %>%
  arrange(Dim.2)
orderPC2 <- orderPC2$Label
#write_tsv(as_tibble(orderPC2), "outputs_inputs/orderPC2.txt")

pcaplot_loadings_up[["data"]] %>%
  arrange(Dim.2) %>%
  ggplot(., aes(x=factor(Label, levels = rev(orderPC2)), y=Dim.2)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(legend.position="none",legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_text(aes(y=0, label=Label), angle = 90,
            fontface = 2, size=1, vjust = 0.5, hjust=1.1) +
  geom_hline(yintercept = 0) +
  ylab(pcaplot_loadings_up[["labels"]][["y"]])

ggsave("outputs_inputs/PCA2_loadings.tiff",
       width = 10, height = 2.5)
# Import plots from rexposome for better aesthetics, stats and modeling (ONLY DOWN VARIABLE GENES) ----
# Samples with missing info in the phenotype file:
#PCA PHE:
load("inputs_rexposome_maaslin/plotPHE_down")

plotPHE_down
plotPHE_down[["data"]] %>%
  filter(Dim %in% c("PC 01","PC 02"))

#PCA samples
load("inputs_rexposome_maaslin/pcaplot_microbiome_down")
pcaplot_microbiome_down[["data"]][["Label"]] <- str_replace_all(pcaplot_microbiome_down[["data"]][["Label"]], "^29", "CL29")

pcaplot_microbiome_down[["data"]] %>% 
  ggplot(., aes(x=Dim.1, y=Dim.2,
                #color=phenotype,
                NULL)) +
  geom_point(size = 4) +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  #scale_color_manual(values = c(Dark24[4],Dark24[2],"gray",Dark24[12],Dark24[1],Dark24[3])) +
  #geom_text_repel(aes(label = phenotype), size = 2, fontface="bold",
  #                color="black", max.overlaps = 20) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #stat_ellipse(level = 0.95) +
  xlab(pcaplot_microbiome_down[["labels"]][["x"]]) + ylab(pcaplot_microbiome_down[["labels"]][["y"]]) + 
  coord_fixed()

pcaplot_microbiome_down[["data"]] %>% 
  ggplot(., aes(x=Dim.1, y=Dim.2, color=phenotype)) +
  geom_point(size = 4) +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  scale_color_manual(values = c(Dark24[4],Dark24[2],"gray",Dark24[12],Dark24[1],Dark24[3])) +
  #geom_text_repel(aes(label = phenotype), size = 2, fontface="bold",
  #                color="black", max.overlaps = 20) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #stat_ellipse(level = 0.95) +
  xlab(pcaplot_microbiome_down[["labels"]][["x"]]) + ylab(pcaplot_microbiome_down[["labels"]][["y"]]) + 
  coord_fixed()

pcaplot_microbiome_down[["data"]] %>% 
  rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  ggplot(., aes(x=Dim.1, y=Dim.2, color=phenotype, shape=Clinical_Outcome)) +
  geom_point(size = 4) +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  scale_color_manual(values = c(Dark24[4],Dark24[2],"gray",Dark24[12],Dark24[1],Dark24[3])) +
  #geom_text_repel(aes(label = phenotype), size = 2, fontface="bold",
  #                color="black", max.overlaps = 20) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #stat_ellipse(level = 0.95) +
  xlab(pcaplot_microbiome_down[["labels"]][["x"]]) + ylab(pcaplot_microbiome_down[["labels"]][["y"]]) + 
  coord_fixed()

pcaplot_microbiome_down[["data"]] %>% 
  dplyr::rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  ggplot(., aes(x=Dim.1, y=Dim.2, color=Clinical_Outcome)) +
  geom_point(size = 4) +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  scale_color_manual(values = c(Dark24[2],Dark24[1])) +
  #geom_text_repel(aes(label = phenotype), size = 2, fontface="bold",
  #                color="black", max.overlaps = 20) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #stat_ellipse(level = 0.95) +
  xlab(pcaplot_loadings_down[["labels"]][["x"]]) + ylab(pcaplot_loadings_down[["labels"]][["y"]]) + 
  coord_fixed()

orderPC1s_down <- pcaplot_microbiome_down[["data"]] %>%
  arrange(Dim.1)
orderPC1s_down <- orderPC1s_down$Label

orderPC2s_down <- pcaplot_microbiome_down[["data"]] %>%
  arrange(Dim.2)
orderPC2s_down <- orderPC2s_down$Label

pcaplot_microbiome_down[["data"]] %>%
  arrange(Dim.1) %>%
  dplyr::rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  ggplot(., aes(x=factor(sample_RNAseq, levels = rev(orderPC1s_down)), y=Dim.1, fill=phenotype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(Dark24[4],Dark24[2],"gray",Dark24[12],Dark24[1],Dark24[3])) +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  #geom_text(aes(y=0, label=Label), angle = 90,
  #          fontface = 2, size=1, vjust = 0.5, hjust=1.1) +
  geom_hline(yintercept = 0) +
  ylab(pcaplot_microbiome_down[["labels"]][["x"]])

pcaplot_microbiome_down[["data"]] %>%
  arrange(Dim.1) %>%
  dplyr::rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  ggplot(., aes(x=factor(sample_RNAseq, levels = rev(orderPC1s_down)), y=Dim.1, fill=observed_features)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "gray",
                      high = "purple",
                      space = "Lab",
                      na.value = "white") +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  #geom_text(aes(y=0, label=Label), angle = 90,
  #          fontface = 2, size=1, vjust = 0.5, hjust=1.1) +
  geom_hline(yintercept = 0) +
  ylab(pcaplot_microbiome_down[["labels"]][["x"]])

pcaplot_microbiome_down[["data"]] %>%
  arrange(Dim.2) %>%
  dplyr::rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  ggplot(., aes(x=factor(sample_RNAseq, levels = rev(orderPC2s_down)), y=Dim.2, fill=Clinical_Outcome)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(Dark24[2],Dark24[1])) +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  #geom_text(aes(y=0, label=Label), angle = 90,
  #          fontface = 2, size=1, vjust = 0.5, hjust=1.1) +
  geom_hline(yintercept = 0) +
  ylab(pcaplot_microbiome_down[["labels"]][["y"]])

pcaplot_microbiome_down[["data"]] %>%
  arrange(Dim.2) %>%
  dplyr::rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  ggplot(., aes(x=factor(sample_RNAseq, levels = rev(orderPC2s_down)), y=Dim.2, fill=healing_time_days)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "gray",
                      high = "purple",
                      space = "Lab",
                      na.value = "white") +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        #axis.text.x = element_blank(), 
        axis.text.x = element_text(angle = 90),
        axis.ticks.x = element_blank()) +
  #geom_text(aes(y=0, label=Label), angle = 90,
  #          fontface = 2, size=1, vjust = 0.5, hjust=1.1) +
  geom_hline(yintercept = 0) +
  ylab(pcaplot_microbiome_down[["labels"]][["y"]])

pcaplot_microbiome_down[["data"]] %>%
  arrange(Dim.2) %>%
  dplyr::rename(sample_RNAseq = Label) %>%
  left_join(phenotype) %>%
  ggplot(., aes(x=factor(sample_RNAseq, levels = rev(orderPC2s_down)), y=Dim.2, fill=size_lesion_mm2)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "gray",
                      high = "purple",
                      space = "Lab",
                      na.value = "white") +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        #axis.text.x = element_blank(), 
        axis.text.x = element_text(angle = 90),
        axis.ticks.x = element_blank()) +
  #geom_text(aes(y=0, label=Label), angle = 90,
  #          fontface = 2, size=1, vjust = 0.5, hjust=1.1) +
  geom_hline(yintercept = 0) +
  ylab(pcaplot_microbiome_down[["labels"]][["y"]])

#PCA loadings:
load("inputs_rexposome_maaslin/pcaplot_loadings_down")

pcaplot_loadings_down[["data"]] %>%
  ggplot(., aes(x=Dim.1, y=Dim.2)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #stat_ellipse(level = 0.95) +
  xlab(pcaplot_loadings_down[["labels"]][["x"]]) + ylab(pcaplot_loadings_down[["labels"]][["y"]]) + 
  coord_fixed()

pcaplot_loadings_down[["data"]] %>%
  ggplot(., aes(x=Dim.1, y=Dim.2)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  #scale_color_manual(values = c(Salvador[5],Salvador[3],Salvador[4],Salvador[6],"gray")) +
  geom_text_repel(aes(label = Label), size = 2, fontface="bold",
                  color="black", max.overlaps = 20) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #stat_ellipse(level = 0.95) +
  xlab(pcaplot_loadings_down[["labels"]][["x"]]) + ylab(pcaplot_loadings_down[["labels"]][["y"]]) + 
  coord_fixed()

orderPC1s_down2 <- pcaplot_loadings_down[["data"]] %>%
  arrange(Dim.1)
orderPC1s_down2 <- orderPC1s_down2$Label
#write_tsv(as_tibble(orderPC1s_down2), "outputs_inputs/orderPC1s_down2.txt")

orderPC2s_down2 <- pcaplot_loadings_down[["data"]] %>%
  arrange(Dim.2)
orderPC2s_down2 <- orderPC2s_down2$Label
#write_tsv(as_tibble(orderPC2s_down2), "outputs_inputs/orderPC2s_down2.txt")

# Import Maaslin outputs NOT UPDATED BECAUSE I DON'T THINK THIS IS IDEAL ----
all_results <- read_delim("inputs_rexposome_maaslin/Maaslin_output_MicroClustSamples/all_results.tsv", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)


# Association Microbiome clusters and gene expression (Maaslin)
all_results_microbiome <- all_results %>%
  filter(metadata == "microbiome_cluster") %>%
  filter(value != "Other") #%>%
  #filter(pval < 0.05) %>% filter(coef > 1 | coef < -1)

# Colunm of significant associations
sig_asso_up <- all_results_microbiome %>%
  filter(pval < 0.05) %>% filter(coef > 1)
sig_asso_down <- all_results_microbiome %>%
  filter(pval < 0.05) %>% filter(coef < -1)
sig_asso_up <- sig_asso_up$feature
sig_asso_down <- sig_asso_down$feature
#write_tsv(as_tibble(sig_asso_up), "up.txt") # to be used in GO
#write_tsv(as_tibble(sig_asso_down), "down.txt") # to be used in GO
#save(sig_asso_up, file = "RStudio_outputs/data/sig_asso_up")
#save(sig_asso_down, file = "RStudio_outputs/data/sig_asso_down")

all_results_microbiome$col0 <- all_results_microbiome$feature
col0_s <- all_results_microbiome$col0 %in% sig_asso_up
col0_s2 <- all_results_microbiome$col0 %in% sig_asso_down
all_results_microbiome$col0[!col0_s & !col0_s2] <- NA

all_results_microbiome$col1 <- all_results_microbiome$feature
col1_s <- all_results_microbiome$col1 %in% sig_asso_up
col1_s2 <- all_results_microbiome$col1 %in% sig_asso_down
all_results_microbiome$col1[col1_s] <- "sig_up"
all_results_microbiome$col1[col1_s2] <- "sig_down"
all_results_microbiome$col1[!col1_s & !col1_s2] <- "notsig"

all_results_microbiome$col2 <- all_results_microbiome$col0
col2_s <- all_results_microbiome$col2 %in% NA
col2_s2 <- all_results_microbiome$col2 %in% sig_asso_up
col2_s3 <- all_results_microbiome$col2 %in% sig_asso_down
all_results_microbiome$col2[col2_s] <- "notsig"
all_results_microbiome$col2[col2_s2 | col2_s3] <- "sig"

# Staphylococcus correlations
View(all_results_microbiome %>%
  filter(value == "Staphylococcus") %>%
    filter(pval <.05) %>%
    arrange(desc(coef)))

View(all_results_microbiome %>%
       filter(value == "Streptococcus") %>%
       filter(pval <.05) %>%
       arrange(desc(coef)))

View(all_results_microbiome %>%
       filter(value == "Heterogenous") %>%
       filter(pval <.05) %>%
       arrange(desc(coef)))

View(all_results_microbiome %>%
       filter(value == "Corynebacterium") %>%
       filter(pval <.05) %>%
       arrange(desc(coef)))

View(all_results_microbiome %>%
       filter(value == "Arcanobacterium") %>%
       filter(pval <.05) %>%
       arrange(desc(coef)))

  
ggplot(., aes(y=-log10(pval), x=coef,size=col1)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("dark gray", Dark24[2], Dark24[1])) +
  scale_size_manual(values = c(0.5,1,1)) +
  geom_text_repel(aes(label = feature), size = 0.5, fontface="bold",colour="black") +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  geom_hline(yintercept = -log10(0.05)) + geom_vline(xintercept = 1) + geom_vline(xintercept = -1) +
  #annotate("text", x=-5, y=32, label=paste(length(sig_genes_down)), size=6, color=Dark24[2], fontface="bold") +
  #annotate("text", x=5, y=32, label=paste(length(sig_genes_up)), size=6, color=Dark24[1], fontface="bold") +
  xlab("logFC CL vs. HS")
  


# Heatmap of microbiome clusters and gene expression OVERrepresented genes ----
load("../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/log2.cpm.filtered.norm")
load("../3rd_LesionRNAseq_Dataset/Unsup_geneExpr_DataReduction/viz_var")
load("../microbiome_Lesion16S-seqDataset/outputs/clusters_che")
load("../microbiome_Lesion16S-seqDataset/outputs/micro_colors")

load("../microbiome_Lesion16S-seqDataset/outputs/Arca_samples")
load("../microbiome_Lesion16S-seqDataset/outputs/Staphy_samples")
load("../microbiome_Lesion16S-seqDataset/outputs/Strep_samples")
load("../microbiome_Lesion16S-seqDataset/outputs/Hetero_samples")
load("../microbiome_Lesion16S-seqDataset/outputs/Coryne_samples")


# Modeling gene expression mt:
orderPC2_v2 <- str_replace_all(orderPC2, "_", "-")
# Only top PC2 loadings:
orderPC2_v3_pos <- rev(orderPC2_v2)[1:50]
orderPC2_v3_neg <- orderPC2_v2[1:50]
save(orderPC2_v3_pos, file = "outputs/orderPC2_v3_pos")

onlymicroclust_samples <- c(Arca_samples,Staphy_samples,
                            Strep_samples,Coryne_samples,
                            Hetero_samples)

log2.cpm.filtered.norm_mtx <- as.data.frame(log2.cpm.filtered.norm[orderPC2_v2,
                                                                   colnames(log2.cpm.filtered.norm) %in% onlymicroclust_samples])

log2.cpm.filtered.norm_mtx <- as.data.frame(log2.cpm.filtered.norm[c(orderPC2_v3_neg,orderPC2_v3_pos),
                                                                   colnames(log2.cpm.filtered.norm) %in% onlymicroclust_samples])

log2.cpm.filtered.norm_mtx <- as.data.frame(log2.cpm.filtered.norm[orderPC2_v3_pos,
                                                                   colnames(log2.cpm.filtered.norm) %in% onlymicroclust_samples])


heatmap.2(as.matrix(log2.cpm.filtered.norm_mtx),
          #Colv=as.dendrogram(microclust),
          #Rowv = as.dendrogram(impgenes),
          dendrogram = "both",
          #RowSideColors=module.color,
          col=colorRampPalette(colors=c("blue","white","red"))(100),
          #ColSideColors = micro_colors_mean,
          scale='row',
          density.info="none", trace="none",
          cexRow=0.6, cexCol=1,margins = c(5,15))

heatmap.2(as.matrix(log2.cpm.filtered.norm[orderPC2_v3_pos,c("CL29180","CL29151","CL29271","CL29256")]),
          cellnote=round(as.matrix(log2.cpm.filtered.norm[orderPC2_v3_pos,c("CL29180","CL29151","CL29271","CL29256")]), 2),
          #Colv=as.dendrogram(microclust),
          Rowv = as.dendrogram(impgenes),
          dendrogram = "both",
          #RowSideColors=module.color,
          col=myheatcol,
          #ColSideColors = micro_colors_mean,
          scale='row',
          density.info="none", trace="none",
          cexRow=1, cexCol=1,margins = c(5,15))


log2.cpm.filtered.norm_mtx2 <- log2.cpm.filtered.norm_mtx %>%
  mutate(Arcanobacterium = rowMeans(.[,Arca_samples])) %>%
  mutate(Streptococcus = rowMeans(.[,Strep_samples])) %>%
  mutate(Heterogenous = rowMeans(.[,colnames(log2.cpm.filtered.norm_mtx) %in% Hetero_samples])) %>%
  mutate(Corynebacterium = rowMeans(.[,colnames(log2.cpm.filtered.norm_mtx) %in% Coryne_samples])) %>%
  mutate(Staphylococcus = rowMeans(.[,colnames(log2.cpm.filtered.norm_mtx) %in% Staphy_samples])) %>%
  select(Staphylococcus,Corynebacterium,
         Streptococcus,Heterogenous,Arcanobacterium)

  
myheatcol <- colorRampPalette(colors=c(Dark24[8],"white","#000080"))(50)

dist_microclust <- dist(t(log2.cpm.filtered.norm_mtx2), method = "canberra") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
microclust <- hclust(dist_microclust, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(microclust)

micro_module_mean <- cutree(microclust, k=5)
micro_module_mean_fac <- factor(micro_module_mean)
micro_colors_mean <- function(micro_module_mean_fac){ if (micro_module_mean_fac==1) Dark24[1] else if(micro_module_mean_fac==2) Dark24[2] else if (micro_module_mean_fac==3) Dark24[3] else if (micro_module_mean_fac==4) "gray" else if (micro_module_mean_fac==5) Dark24[4]}
micro_colors_mean <- unlist(lapply(micro_module_mean_fac, micro_colors_mean))

# Calculate from the mean:
dist_impgenes <- dist(log2.cpm.filtered.norm_mtx2, method = "euclidean") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
impgenes <- hclust(dist_impgenes, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(impgenes)

# Calculate from the gene expression
dist_impgenes <- dist(log2.cpm.filtered.norm_mtx, method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
impgenes <- hclust(dist_impgenes, method = "ward.D2") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(impgenes)

heatmap.2(as.matrix(log2.cpm.filtered.norm_mtx2),
          Colv=as.dendrogram(microclust),
          Rowv = as.dendrogram(impgenes),
          dendrogram = "both",
          #RowSideColors=module.color,
          col=myheatcol,
          ColSideColors = micro_colors_mean,
          scale='row',
          labCol = NA, key = T,
          density.info="none", trace="none",
          cexRow=1.05, cexCol=0.5,margins = c(1,15))

heatmap.2(as.matrix(log2.cpm.filtered.norm_mtx2),
          Colv=as.dendrogram(microclust),
          Rowv = as.dendrogram(impgenes),
          dendrogram = "both",
          #RowSideColors=module.color,
          col=myheatcol,
          ColSideColors = micro_colors_mean,
          scale='row',
          labCol = NA,
          density.info="none", trace="none",
          cexRow=1, cexCol=0.5,margins = c(15,15))

# Heatmap from K-means ----
# Another way to cluster genes would be by k-means (instead what I did with manual mean calculation).
# Read: https://statsandr.com/blog/clustering-analysis-k-means-and-hierarchical-clustering-by-hand-and-in-r/

kmeansclust <- kmeans(log2.cpm.filtered.norm_mtx, centers = 5)
# I just didn't continue...

heatmap.2(as.matrix(log2.cpm.filtered.norm[names(sort(kmeansclust[["cluster"]])),]),
          Colv=F,
          Rowv = F,
          dendrogram = "none",
          #RowSideColors=module.color,
          col=myheatcol,
          #ColSideColors = micro_colors_mean,
          scale='row',
          labCol = NA,
          density.info="none", trace="none",
          cexRow=0.5, cexCol=0.5,margins = c(1,15))
      
# Heatmap of microbiome clusters and gene expression UNDERrepresented genes ----
# Modeling gene expression mt:
orderPC1s_down22 <- str_replace_all(orderPC1s_down2, "_", "-")

# Only top PC1 loadings:
orderPC1s_down22_pos <- rev(orderPC1s_down22)[1:50]

log2.cpm.filtered.norm_mtx_down <- as.data.frame(log2.cpm.filtered.norm[orderPC1s_down22_pos,
                                                                   colnames(log2.cpm.filtered.norm) %in% onlymicroclust_samples])


heatmap.2(as.matrix(log2.cpm.filtered.norm_mtx_down),
          dendrogram = "both",
          #RowSideColors=module.color,
          col=myheatcol,
          #ColSideColors = micro_colors_mean,
          scale='row',
          density.info="none", trace="none",
          cexRow=1, cexCol=0.6,margins = c(5,15))

heatmap.2(as.matrix(log2.cpm.filtered.norm[orderPC2_v3_pos,c("CL29180","CL29151","CL29271","CL29256")]),
          cellnote=round(as.matrix(log2.cpm.filtered.norm[orderPC2_v3_pos,c("CL29180","CL29151","CL29271","CL29256")]), 2),
          #Colv=as.dendrogram(microclust),
          Rowv = as.dendrogram(impgenes),
          dendrogram = "both",
          #RowSideColors=module.color,
          col=myheatcol,
          #ColSideColors = micro_colors_mean,
          scale='row',
          density.info="none", trace="none",
          cexRow=1, cexCol=1,margins = c(5,15))


log2.cpm.filtered.norm_mtx2 <- log2.cpm.filtered.norm_mtx %>%
  mutate(Arcanobacterium = rowMeans(.[,Arca_samples])) %>%
  mutate(Streptococcus = rowMeans(.[,Strep_samples])) %>%
  mutate(Heterogenous = rowMeans(.[,colnames(log2.cpm.filtered.norm_mtx) %in% Hetero_samples])) %>%
  mutate(Corynebacterium = rowMeans(.[,colnames(log2.cpm.filtered.norm_mtx) %in% Coryne_samples])) %>%
  mutate(Staphylococcus = rowMeans(.[,colnames(log2.cpm.filtered.norm_mtx) %in% Staphy_samples])) %>%
  select(Staphylococcus,Corynebacterium,
         Streptococcus,Heterogenous,Arcanobacterium)

myheatcol <- colorRampPalette(colors=c(Dark24[8],"white","purple"))(100)

dist_microclust <- dist(t(log2.cpm.filtered.norm_mtx2), method = "canberra") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
microclust <- hclust(dist_microclust, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(microclust)

micro_module_mean <- cutree(microclust, k=5)
micro_module_mean_fac <- factor(micro_module_mean)
micro_colors_mean <- function(micro_module_mean_fac){ if (micro_module_mean_fac==1) Dark24[1] else if(micro_module_mean_fac==2) Dark24[2] else if (micro_module_mean_fac==3) Dark24[3] else if (micro_module_mean_fac==4) "gray" else if (micro_module_mean_fac==5) Dark24[4]}
micro_colors_mean <- unlist(lapply(micro_module_mean_fac, micro_colors_mean))

# Calculate from the mean:
dist_impgenes <- dist(log2.cpm.filtered.norm_mtx2, method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
impgenes <- hclust(dist_impgenes, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(impgenes)

# Calculate from the gene expression
dist_impgenes <- dist(log2.cpm.filtered.norm_mtx, method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
impgenes <- hclust(dist_impgenes, method = "ward.D2") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(impgenes)

heatmap.2(as.matrix(log2.cpm.filtered.norm_mtx2),
          Colv=as.dendrogram(microclust),
          Rowv = as.dendrogram(impgenes),
          dendrogram = "both",
          #RowSideColors=module.color,
          col=myheatcol,
          ColSideColors = micro_colors_mean,
          scale='row',
          labCol = NA,
          density.info="none", trace="none",
          cexRow=1, cexCol=0.5,margins = c(1,15))

# Visualizing individual gene and microbiome clusters ----
load("../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/CL_samples")
gene_viz <- log2.cpm.filtered.norm[c("LTB4R","SLPI","PRF1","GNLY","ISG15",
                                     "GZMB","HEPHL1","IL1B","APOBEC3A",
                                     "UPP1","CXCL8","OSM","IL1A"),CL_samples]

gene_viz <- log2.cpm.filtered.norm[c("GBP5", "TNFRSF6B", "TNFRSF18", "CCL4L2", "SLC11A1",
                                     "IL24", "ACOD1", "CXCL1", "IL27", "FPR2", "CXCL10",
                                     "CXCL11", "CCL8", "CCL7", "IL1B", "TNIP3", "CCL4", "CCL3", "CCL2"),CL_samples]

gene_viz <- log2.cpm.filtered.norm[c("LTB4R","PTGES2","IL1A","ROS1"),CL_samples]
gene_viz <- log2.cpm.filtered.norm[c("CCL20","RORA","RORC","RBPJ","IL23R"),CL_samples]

gene_viz <- log2.cpm.filtered.norm[c("IL7","IL7R"),CL_samples]

melt_viz <- as.data.frame(melt(gene_viz))
melt_viz$X2 <- as.character(melt_viz$X2)
colnames(melt_viz)[2] <- "sample_RNAseq"
melt_viz <- left_join(melt_viz, phenotype, by="sample_RNAseq")
#colnames(melt_viz)


melt_viz %>%
  filter(microbiome_cluster %ni% c(NA, "Other")) %>%
ggplot(.,
       aes(x=microbiome_cluster, y=value, color=microbiome_cluster)) +
  geom_jitter(position=position_jitter(0.1), size=1) +
  theme_classic() + #scale_color_manual(values = c("gray",Dark24[3],Dark24[2],Dark24[1],Dark24[4],Dark24[5],Dark24[6])) +
  scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="black") +
  stat_summary(fun = mean,
               geom = "crossbar", width = 0.3, color="black") +
  xlab("") + ylab(paste("CPM (log2)")) + 
  #  ylim(0,5) +
  # facet_grid(. ~ Var1) +
  facet_wrap(. ~ X1, scales = "free", nrow = 3)

melt_viz %>%
  filter(microbiome_cluster %in% c("Staphylococcus", "Heterogenous")) %>%
  ggplot(.,
         aes(x=microbiome_cluster, y=value, color=microbiome_cluster)) +
  geom_violin(size = 0.7, trim = F) +
  theme_classic() + scale_color_manual(values = c("gray",Dark24[1])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="black") +
  #stat_summary(fun = mean,
  #             geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "t.test", #label.y = 0,
                     aes(label = ..p.signif..),
                     label.x = 1, size = 3.5, color="black") +
  xlab("") + ylab(paste("CPM (log2)")) + 
  #  ylim(0,5) +
  # facet_grid(. ~ Var1) +
  facet_wrap(. ~ X1, scales = "free", nrow = 3)

# xCell and MCP counter for Microbiome clusters and Differential cell analysis ----
load("../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/res_xcell")
load("../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/res_mcp_counter")
res_xcell_mtx <- column_to_rownames(res_xcell,"cell_type")
res_xcell_mtx <- as.data.frame(res_xcell_mtx[,colnames(res_xcell_mtx) %in% phenotype_bothData$sample_RNAseq]) # from MergingMicrobiomeProfile...R

res_mcp_mtx <- column_to_rownames(res_mcp_counter,"cell_type")
res_mcp_mtx <- as.data.frame(res_mcp_mtx[,colnames(res_mcp_mtx) %in% phenotype_bothData$sample_RNAseq]) # from MergingMicrobiomeProfile...R
res_mcp_mtx <- res_mcp_mtx[-11,]
rownames(res_mcp_mtx)[3] <- "CTLs"

# Filtering for cell subtypes present in at least certain number of samples:
#keepers <- rowSums(res_xcell_mtx>0)>=4
#dim(res_xcell_mtx)
#res_xcell_mtx <- res_xcell_mtx[keepers,]
#dim(res_xcell_mtx)

# Log2 scale?
res_xcell_mtx <- log2(res_xcell_mtx)

res_xcell_mtx2 <- res_xcell_mtx %>%
  mutate(Arcanobacterium = rowMeans(.[,Arca_samples])) %>%
  mutate(Streptococcus = rowMeans(.[,Strep_samples])) %>%
  mutate(Heterogenous = rowMeans(.[,colnames(res_xcell_mtx) %in% Hetero_samples])) %>%
  mutate(Corynebacterium = rowMeans(.[,colnames(res_xcell_mtx) %in% Coryne_samples])) %>%
  mutate(Staphylococcus = rowMeans(.[,colnames(res_xcell_mtx) %in% Staphy_samples])) %>%
  select(Staphylococcus,Corynebacterium,
         Streptococcus,Heterogenous,Arcanobacterium)

res_xcell_mtx3 <- res_xcell_mtx %>%
  mutate(Arcanobacterium = rowMeans(.[,Arca_samples])) %>%
  mutate(Streptococcus = rowMeans(.[,Strep_samples])) %>%
  mutate(Heterogenous = rowMeans(.[,colnames(res_xcell_mtx) %in% Hetero_samples])) %>%
  mutate(Corynebacterium = rowMeans(.[,colnames(res_xcell_mtx) %in% Coryne_samples])) %>%
  mutate(Staphylococcus = rowMeans(.[,colnames(res_xcell_mtx) %in% Staphy_samples])) %>%
  select(Staphylococcus,Heterogenous)

# Calculate from mean:
dist_impcells <- dist(res_xcell_mtx2, method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
impcells <- hclust(dist_impcells, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(impcells)

# Calculate from total scores:
dist_impcells <- dist(res_xcell_mtx2, method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
impcells <- hclust(dist_impcells, method = "ward.D2") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(impcells)

heatmap.2(as.matrix(res_xcell_mtx2),
          Colv= as.dendrogram(microclust),
          Rowv = as.dendrogram(impcells),
          dendrogram = "both",
          #RowSideColors=module.color,
          col=myheatcol,
          ColSideColors = micro_colors_mean,
          scale='row',
          labCol = NA,
          density.info="none", trace="none",
          cexRow=0.7, cexCol=0.5,margins = c(1,15))

# Only Hetero vs Staphy
heatmap.2(as.matrix(res_xcell_mtx3),
          #Colv= as.dendrogram(microclust),
          #Rowv = as.dendrogram(impcells),
          dendrogram = "none",
          Colv = F,
          Rowv = F,
          #RowSideColors=module.color,
          col=myheatcol,
          #ColSideColors = micro_colors_mean,
          scale='row',
          labCol = T,
          density.info="none", trace="none",
          cexRow=0.7, cexCol=0.5,margins = c(5,15))

res_xcell_mod <- res_xcell %>%
  pivot_longer(!cell_type, names_to = "sample_RNAseq", values_to = "xCell score") %>%
  left_join(phenotype)

res_xcell_mod %>%
  filter(microbiome_cluster %ni% c(NA, "Other")) %>%
  ggplot(.,
         aes(x=microbiome_cluster, y=`xCell score`, color=microbiome_cluster)) +
  geom_jitter(position=position_jitter(0.1), size=1) +
  theme_classic() + scale_color_manual(values = c("gray",Dark24[3],Dark24[2],Dark24[1],Dark24[4])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="black") +
  stat_summary(fun = median,
               geom = "crossbar", width = 0.3, color="black") +
  xlab("") + ylab(paste("xCell score")) + 
  #  ylim(0,5) +
  # facet_grid(. ~ Var1) +
  facet_wrap(. ~ cell_type, scales = "free")

# Just Staphy and Hetero
res_xcell_mod %>%
  filter(microbiome_cluster %in% c("Staphylococcus", "Heterogenous")) %>%
  ggplot(.,
         aes(x=microbiome_cluster, y=`xCell score`, color=microbiome_cluster)) +
  geom_violin(size = 0.7, trim = F) +
  theme_classic() + scale_color_manual(values = c("gray",Dark24[1])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="black") +
  stat_summary(fun = mean,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "t.test", #label.y = 0,
                     aes(label = ..p.signif..),
                     label.x = 1, size = 3.5, color="black") +
  xlab("") + ylab(paste("xCell score")) + 
  #  ylim(0,5) +
  # facet_grid(. ~ Var1) +
  facet_wrap(. ~ cell_type, scales = "free")

# Differential cell analysis:
group_microclust <- factor(phenotype_bothData$microbiome_cluster)
design_microclust <- model.matrix(~0 + group_microclust)
colnames(design_microclust) <- levels(group_microclust)

v.DEGList.cell <- voom(res_mcp_mtx, design_microclust, plot = FALSE)
fit_micro_cell <- lmFit(v.DEGList.cell, design_microclust)
contrast.matrix_micro_cell <- makeContrasts(Staphylococcus = Staphylococcus - Heterogeneous,
                                       Arcanobacterium = Arcanobacterium - Heterogeneous,
                                       Corynebacterium = Corynebacterium - Heterogeneous,
                                       Streptococcus = Streptococcus - Heterogeneous,
                                       levels=design_microclust)

fits_micro_cell <- contrasts.fit(fit_micro_cell, contrast.matrix_micro_cell)
ebFit_micro_cell <- eBayes(fits_micro_cell)
myTopHitsmicro_cell_Staph <- topTable(ebFit_micro_cell, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHitsmicro_cell_Staph <- myTopHitsmicro_cell_Staph %>%
  add_column(microbiome_cluster = "M6") %>%
  rownames_to_column("Cell type")
myTopHitsmicro_cell_Arcano <- topTable(ebFit_micro_cell, adjust ="BH", coef=2, number=40000, sort.by="logFC")
myTopHitsmicro_cell_Arcano <- myTopHitsmicro_cell_Arcano %>%
  add_column(microbiome_cluster = "M7") %>%
  rownames_to_column("Cell type")
myTopHitsmicro_cell_Coryne <- topTable(ebFit_micro_cell, adjust ="BH", coef=3, number=40000, sort.by="logFC")
myTopHitsmicro_cell_Coryne <- myTopHitsmicro_cell_Coryne %>%
  add_column(microbiome_cluster = "M5") %>%
  rownames_to_column("Cell type")
myTopHitsmicro_cell_Strep <- topTable(ebFit_micro_cell, adjust ="BH", coef=4, number=40000, sort.by="logFC")
myTopHitsmicro_cell_Strep <- myTopHitsmicro_cell_Strep %>%
  add_column(microbiome_cluster = "M4") %>%
  rownames_to_column("Cell type")

myTopHitsmicro_cell_all <- rbind(myTopHitsmicro_cell_Staph,
                                 myTopHitsmicro_cell_Arcano,
                                 myTopHitsmicro_cell_Coryne,
                                 myTopHitsmicro_cell_Strep)

myTopHitsmicro_cell_all %>%
  mutate(`Cell type` = as.factor(`Cell type`)) %>%
  mutate(FC = 2^logFC) %>%
  mutate(P.value = case_when(
    P.Value < 0.001 ~ "***",
    P.Value < 0.01 ~ "**",
    P.Value < 0.05 ~ "*",
    TRUE ~ "NA")) %>%
  mutate(P.value = na_if(P.value, "NA")) %>%
  ggplot(., aes(microbiome_cluster,`Cell type`, color = P.Value, size=FC)) +
  geom_point() +
  theme_classic() +
  scale_color_gradient(low = Dark24[15],
                       high = "white") +
  #geom_text(aes(label = P.value), color = "white", vjust=0.8, size=5) +
  theme(axis.text = element_text(size = 13), legend.text=element_text(size=13),
        axis.title = element_blank(), legend.title = element_text(size = 13))
  




# Parasite Load and Microbiome clusters ----
load("../3rd_LesionRNAseq_Dataset/qPCRs/Lbraziliensis_qPCR_ImportingDataIntoR/parasiteload_3rdDataset")
parasiteload_3rdDataset <- parasiteload_3rdDataset %>%
  left_join(phenotype)

parasiteload_3rdDataset %>%
  filter(microbiome_cluster %ni% NA) %>%
  ggplot(.,
         aes(x=factor(microbiome_cluster, levels = c("Staphylococcus","Arcanobacterium","Corynebacterium",
                                                       "Streptococcus","Heterogenous","Lactobacillus","Finegoldia")),
                      y=Lbraziliensis_mg/1000, color=microbiome_cluster)) +
  geom_jitter(position=position_jitter(0.1), size=3) +
  theme_classic() + scale_color_manual(values = c("gray",Dark24[3],Dark24[2],Dark24[1],Dark24[4],Dark24[5],Dark24[6])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=11, angle = 90, hjust = 1, vjust = 0.5), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none") +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="black") +
  stat_summary(fun = median,
               geom = "crossbar", width = 0.3, color="black") +
  xlab("") + ylab(paste("number of L. braziliensis/mg biopsy (/1000)"))

parasiteload_3rdDataset %>%
  filter(microbiome_cluster %ni% c(NA, "Finegoldia", "Lactobacillus")) %>%
  ggplot(.,
         aes(x=factor(microbiome_cluster, levels = c("Heterogenous","Streptococcus","Corynebacterium",
                                                     "Staphylococcus","Arcanobacterium")),
             y=Lbraziliensis_mg/1000, color=microbiome_cluster)) +
  geom_jitter(position=position_jitter(0.1), size=3) +
  theme_classic() + scale_color_manual(values = c("gray",Dark24[3],Dark24[2],Dark24[1],Dark24[4],Dark24[5],Dark24[6])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="right") +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="black") +
  stat_summary(fun = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",# label.y = 2,
  #                   aes(label = ..p.signif..),
  aes(label = ..p.value..),
  #label.x = 1.5,
                     size = 5, color="black") +
  xlab("") + ylab(paste("number of L. braziliensis/mg biopsy (/1000)"))

parasiteload_3rdDataset %>%
  filter(microbiome_cluster %in% c("Staphylococcus", "Heterogenous")) %>%
  filter(Lbraziliensis_mg != "NA") %>%
  ggplot(.,
         aes(x=factor(microbiome_cluster, levels = c("Heterogenous","Staphylococcus")),
             y=Lbraziliensis_mg/1000, color=microbiome_cluster)) +
  geom_jitter(position=position_jitter(0.1), size=3) +
  theme_classic() + scale_color_manual(values = c("gray",Dark24[1])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="black") +
  stat_summary(fun = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "wilcox.test", label.y = 30,
                     aes(label = ..p.signif..),
                     #aes(label = ..p.value..),
                     label.x = 1.3,
                     size = 8, color="black") +
  xlab("") + ylab(paste("number of L. braziliensis/mg biopsy (/1000)"))

# or another way to plot this, and match with the RA relative plots from the 16S sequencing dataset:
load("../microbiome_Lesion16S-seqDataset/outputs/patient_order2")
patient_order2 <- str_replace_all(patient_order2, "^", "CL")
parasiteload_3rdDataset %>%
  filter(microbiome_cluster != "NA") %>%
ggplot(.,
       aes(x=factor(sample_RNAseq, levels = patient_order2), y=Lbraziliensis_mg/10000)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5, angle = 90), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="right") +
  xlab("") + ylab(paste("number of parasites")) +
  #geom_hline(yintercept = 10, color = Dark24[3], linetype=2) +
  #ylim(0, 400) +
  NULL

# Biomarker for treatment outcome gene signature and Microbiome clusters ----
# Because I've found that Arcanobacterium samples expressed more genes that we know are associated with treatment outcome,
# I'll make this information official and run ssGSEA with the genes that are biomarkers for TA in at least 3 RNA-seq datasets.
# and verify funcional enrichment in all samples.
library(GSVA)
library(GSEABase)
Failure_geneSig <- getGmt("../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/myTopHits_treat_var_imp.gmt", geneIdType=SymbolIdentifier())

# GSVA Amorim ----
GSVA.res_Failure <- gsva(log2.cpm.filtered.norm[,CL_samples], #your data
                 Failure_geneSig, # your gene set collection
                 #min.sz=5, max.sz=500, #criteria for filtering gene sets
                 mx.diff=FALSE,
                 method="ssgsea")

GSVA.res_Failure <- rownames_to_column(as.data.frame(GSVA.res_Failure), "Signature")

my_comparisons <- list(#c("Staphylococcus", "Arcanobacterium"),
                       c("Heterogeneous", "Streptococcus"),
                       c("Heterogeneous", "Corynebacterium"),
                       c("Heterogeneous", "Staphylococcus"),
                       c("Heterogeneous", "Arcanobacterium")
                       )

my_comparisonsB <- list(c("Heterogeneous", "Staphylococcus"))
GSVA.res_Failure %>%
  pivot_longer(!Signature, names_to = "sample_RNAseq", values_to = "Failure_geneSig") %>%
  dplyr::select(-Signature) %>%
  left_join(phenotype) %>%
  filter(microbiome_cluster %ni% c(NA, "Finegoldia","Lactobacillus")) %>%
  ggplot(.,
         aes(x=microbiome_cluster, y=Failure_geneSig, color=microbiome_cluster)) +
  geom_violin(size = 0.7, trim = F) +
  geom_sina(size=2.3) +
  theme_classic() + scale_color_manual(values = c("gray",Dark24[3],Dark24[2],Dark24[1],Dark24[4],Dark24[3],"black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="right") +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="black") +
  stat_summary(fun = median,
               geom = "crossbar", width = 0.3, color="black") +
  #stat_compare_means(comparisons = my_comparisons,
  #                   method = "wilcox.test",# label.y = 2,
  #                   aes(label = ..p.signif..),
                     #aes(label = ..p.value..),
                     #label.x = 1.5,
  #                   size = 5, color="black") +
  xlab("") + ylab(paste("Failure gene signature ES"))

GSVA.res_Failure %>%
  pivot_longer(!Signature, names_to = "sample_RNAseq", values_to = "Failure_geneSig") %>%
  dplyr::select(-Signature) %>%
  left_join(phenotype) %>%
  filter(microbiome_cluster %in% c("Staphylococcus", "Heterogeneous")) %>%
  ggplot(.,
         aes(x=factor(microbiome_cluster, levels = c("Heterogeneous","Staphylococcus")), y=Failure_geneSig, color=microbiome_cluster)) +
  geom_violin(size = 0.7, trim = F) +
  geom_sina(size=2.3) +
  theme_classic() + scale_color_manual(values = c(Dark24[1],"gray")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="black") +
  stat_summary(fun = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "wilcox.test",# label.y = 2,
                     aes(label = ..p.signif..),
                     #aes(label = ..p.value..),
                     label.x = 1.4,
                     size = 8, color="black") +
  xlab("") + ylab(paste("Failure gene signature ES"))

# IL17 gene signature and Microbiome clusters ----
# Because we know that IL-17 is an important signaling associated with bacterial infections/colonization
# I took the gene list created by Tej in:
# ~/Dropbox/Analysis_OtherLabs_People/Tej/IL17_pathway_for_human_CamilaMod2.txt
# ~/Dropbox/Analysis_OtherLabs_People/Tej/ILCs_project.R
# I will run ssGSEA with this IL17 signature and verify funcional enrichment in all samples.
library(GSVA)
library(GSEABase)
IL17_geneSig <- getGmt("../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/tejsIL17geneset.gmt", geneIdType=SymbolIdentifier())

GSVA.res_IL17 <- gsva(log2.cpm.filtered.norm[,CL_samples], #your data
                         IL17_geneSig, # your gene set collection
                         #min.sz=5, max.sz=500, #criteria for filtering gene sets
                         mx.diff=FALSE,
                         method="ssgsea")

GSVA.res_IL17 <- rownames_to_column(as.data.frame(GSVA.res_IL17), "Signature")

GSVA.res_IL17 %>%
  pivot_longer(!Signature, names_to = "sample_RNAseq", values_to = "IL17_geneSig") %>%
  dplyr::select(-Signature) %>%
  left_join(phenotype) %>%
  filter(microbiome_cluster %ni% c(NA, "Other")) %>%
  ggplot(.,
         aes(x=microbiome_cluster, y=IL17_geneSig, color=microbiome_cluster)) +
  geom_jitter(position=position_jitter(0.1), size=3) +
  theme_classic() + scale_color_manual(values = c("gray",Dark24[3],Dark24[2],Dark24[1],Dark24[4])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="black") +
  stat_summary(fun = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test",# label.y = 2,
                     aes(label = ..p.signif..),
                     #aes(label = ..p.value..),
                     #label.x = 1.5,
                     size = 3, color="black") +
  xlab("") + ylab(paste("Tej's IL17 gene signature ES"))

# Arcanobacteirum samples info ----
mean(LeishOmics$age, na.rm = TRUE)
median(LeishOmics$age, na.rm = TRUE)

mean(LeishOmics$DTH_mm2, na.rm = TRUE)
median(LeishOmics$DTH_mm2, na.rm = TRUE)

mean(LeishOmics$size_lesion_mm2, na.rm = TRUE)
median(LeishOmics$size_lesion_mm2, na.rm = TRUE)

Arca_factor <- factor(c("29271","29213","29256","29272","29513","29159","29231","29151"),
                      levels = c("29271","29213","29256","29272","29513","29159","29231","29151"))
Arca_meta <- LeishOmics %>%
  filter(LTCP_patient_ID %in% Arca_factor) %>%
  mutate(LTCP_patient_ID = fct_relevel(LTCP_patient_ID, levels(Arca_factor))) %>%
  arrange(LTCP_patient_ID) %>%
  select(LTCP_patient_ID,age,sex,local_of_lesion,size_lesion_mm2,n_lesions,lymphadenopathy,
         DTH_mm2,treatment_outcome,healing_time_days)
Arca_meta

heatmap.2(as.matrix(LeishOmics_heat[levels(Arca_factor),]), 
                           dendrogram = "none",
                           Colv=F, Rowv = F,
                           col=colorheat, key = F,
                           density.info="none", trace="none",  
                           cexCol=1, cexRow = 1.5,
                           colsep=1:nrow(LeishOmics_heat),
                           rowsep=1:nrow(LeishOmics_heat),
                           sepcolor = "black", margins = c(10,10))

# DGE between micro clusters, intercepting for ~Hetereogeneous group, HC and Networks from gene modules ----
library(limma)
library(edgeR) 

# Get myDGEList.filtered.norm:
v.DEGList.filtered.norm_micro <- voom(myDGEList.filtered.norm[,phenotype_bothData$sample_RNAseq],
                                      design_microclust, plot = FALSE)
fit_micro <- lmFit(v.DEGList.filtered.norm_micro, design_microclust)
contrast.matrix_micro <- makeContrasts(Staphylococcus = Staphylococcus - Heterogeneous,
                                 Arcanobacterium = Arcanobacterium - Heterogeneous,
                                 Corynebacterium = Corynebacterium - Heterogeneous,
                                 Streptococcus = Streptococcus - Heterogeneous,
                                 levels=design_microclust)

fits_micro <- contrasts.fit(fit_micro, contrast.matrix_micro)
ebFit_micro <- eBayes(fits_micro)
myTopHitsmicro <- topTable(ebFit_micro, adjust ="BH", coef=1, number=40000, sort.by="logFC")

write_tsv(rownames_to_column(myTopHitsmicro, "geneSymbol"), "../../../../Desktop/SupplementalTable5_DEGs_microbiomeClusters.txt")


results_micro <- decideTests(ebFit_micro, method="global", adjust.method="BH", p.value=0.76, lfc=0.59)
colnames(v.DEGList.filtered.norm_micro$E) <- phenotype_bothData$sample_RNAseq
diffGenes_micro <- v.DEGList.filtered.norm_micro$E[results_micro[,1] !=0,]
clustRows_micro <- hclust(as.dist(1-cor(t(diffGenes_micro), method="pearson")), method="complete") 
clustColumns_micro <- hclust(as.dist(1-cor(diffGenes_micro, method="spearman")), method="ward.D2") #cluster columns by spearman correlation

micro_colors2 <- function(group_microclust){ if (group_microclust=="Staphylococcus") Dark24[1] else if (group_microclust=="Arcanobacterium") Dark24[4] else if (group_microclust=="Corynebacterium") Dark24[2] else if (group_microclust=="Streptococcus") Dark24[3] else if (group_microclust=="Heterogeneous") "gray" else if (group_microclust=="Lactobacillus") Dark24[5] else Dark24[6]}
micro_colors2 <- unlist(lapply(group_microclust, micro_colors2))

clust.assign <- cutree(clustRows_micro, k=2)
module.color <- c(Dark24[14],Dark24[22])
module.color <- module.color[as.vector(clust.assign)] 

myheatcol3 <- colorRampPalette(colors=c(Dark24[8],"white","#000080"))(50)
library(RColorBrewer)
myheatcol3 <- brewer.pal(name="BrBG", n=20)
myheatcol3 <- myheatcol3[c(-5,-6,-7)]
heatmap.2(diffGenes_micro, 
          Rowv=as.dendrogram(clustRows_micro), 
          #Colv=as.dendrogram(clustColumns_micro),
          RowSideColors=module.color,
          ColSideColors=micro_colors2,
          col=rev(brewer.pal(name="BrBG", n=20)[-6]),
          scale='row', labRow=NA, #labCol = NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=0.5, margins=c(8,20)) 

patient_orderT <- heatmap.2(diffGenes_micro, 
          Rowv=as.dendrogram(clustRows_micro), 
          #Colv=as.dendrogram(clustColumns_micro),
          RowSideColors=module.color,
          ColSideColors=micro_colors2,
          col=rev(brewer.pal(name="BrBG", n=20)[-6]),
          scale='row', labRow=NA, #labCol = NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=0.5, margins=c(8,20)) 

# Picking module:
# Module 1:
myModule1 <- diffGenes_micro[names(clust.assign[clust.assign == 1]),] 
hrsub1 <- hclust(as.dist(1-cor(t(myModule1), method="pearson")), method="complete") 

clust.assign1 <- cutree(hrsub1, k=5)
module.color1 <- Dark24[6:20]
module.color1 <- module.color1[as.vector(clust.assign1)] 

# Create heatmap for chosen sub-cluster.
heatmap.2(myModule1, 
          Rowv=as.dendrogram(hrsub1), 
          Colv = as.dendrogram(clustColumns_micro),
          ColSideColors=micro_colors2,
          RowSideColors=module.color1,

          scale='row', labRow=NA, labCol = NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20))

rownames(myModule1[names(clust.assign1[clust.assign1 == 5]),])
write_tsv(as.data.frame(rownames(myModule1)), "outputs/myModule1.txt")

# GO BP on this module:
Module1_GO <- read_delim("outputs/myModule1_GOBP.txt", 
                                                  delim = "\t", escape_double = FALSE, 
                                                  trim_ws = TRUE)
Module1_GO %>%
  slice_head(n = 20) %>%
  mutate(Term = str_remove(Term,".*~")) %>%
  mutate(Term = replace(Term, Term == "antigen processing and presentation, endogenous lipid antigen via MHC class Ib", "APC, endogenous lipid antigen via MHC class Ib")) %>%
  mutate(Term = replace(Term, Term == "antigen processing and presentation, exogenous lipid antigen via MHC class Ib", "APC, exogenous lipid antigen via MHC class Ib")) %>%
  arrange(desc(Benjamini)) %>%
  filter(Benjamini < 0.05) %>%
  ggplot(., aes(x=Count, y=fct_inorder(Term), fill=-log10(Benjamini))) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradient(
    low = "#DC91A9",
    high = Dark24[1]) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size = 15)) +
  xlab("Gene Count")


# Module 2:
myModule2 <- diffGenes_micro[names(clust.assign[clust.assign == 2]),] 
hrsub2 <- hclust(as.dist(1-cor(t(myModule2), method="pearson")), method="complete") 

clust.assign2 <- cutree(hrsub2, k=5)
module.color2 <- Dark24[6:20]
module.color2 <- module.color2[as.vector(clust.assign2)] 

# Create heatmap for chosen sub-cluster.
heatmap.2(myModule2, 
          Rowv=as.dendrogram(hrsub2), 
          Colv = as.dendrogram(clustColumns_micro),
          ColSideColors=micro_colors2,
          RowSideColors=module.color2,

          scale='row', labRow=NA, labCol = NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20))

rownames(myModule2[names(clust.assign2[clust.assign2 == 3]),])

write_tsv(as.data.frame(rownames(myModule2)), "outputs/myModule2.txt")

# GO BP on this module:
Module2_GO <- read_delim("outputs/myModule2_GOBP.txt", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)


Module2_GO %>%
  slice_head(n = 20) %>%
  mutate(Term = str_remove(Term,".*~")) %>%
  #mutate(Term = replace(Term, Term == "homophilic cell adhesion via plasma membrane adhesion molecules", "cellular adhesion")) %>%
  arrange(desc(Benjamini)) %>%
  filter(Benjamini < 0.05) %>%
  ggplot(., aes(x=Count, y=fct_inorder(Term), fill=-log10(Benjamini))) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradient(
    low = "#DC91A9",
    high = Dark24[1]) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size = 15)) +
  xlab("Gene Count")


# Networks from these 2 modules:
library(broom)
library(patchwork)
library(corrplot)
library(Hmisc)
library(tidygraph)
library(ggraph)
library(corrr)

corr1 <- t(myModule1) %>% correlate(method = "pearson") %>%
  shave(upper = TRUE) %>%
  stretch(na.rm = TRUE) %>%
  arrange(desc(r))# %>%
#top_n(1:25)

# I can't run this network in my laptop. So I'm moving this object to the server:
save(corr1, file = "outputs/corr1")

# Association Healing time and sample clustering from master heatmap of DGEs above ----
#extract sample order from dendrogram:
patient_orderT <- as.data.frame(colnames(diffGenes_micro)[patient_orderT$colInd])[,1]

phenotype %>%
  filter(sample_RNAseq %in% patient_orderT) %>%
  ggplot(., aes(x=factor(sample_RNAseq, levels = patient_orderT),
                y=healing_time_days)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text = element_text(angle = 90)) +
  geom_hline(yintercept = 90) + geom_hline(yintercept = 180)
  








# MyTopHits microclusters top 10 genes ----
StaphGenes <- topTable(ebFit_micro, adjust ="BH", coef=1, number=40000, sort.by="logFC")
StaphGenes <- StaphGenes %>%
  rownames_to_column("geneSymbol") %>%
  filter(P.Value < 0.01) %>%
  arrange(desc(logFC))
StaphGenes <- c(StaphGenes$geneSymbol[1:10],rev(StaphGenes$geneSymbol[1:10]))

ArcanoGenes <- topTable(ebFit_micro, adjust ="BH", coef=2, number=40000, sort.by="logFC")
ArcanoGenes <- ArcanoGenes %>%
  rownames_to_column("geneSymbol") %>%
  filter(P.Value < 0.01) %>%
  arrange(desc(logFC))
ArcanoGenes <- c(ArcanoGenes$geneSymbol[1:10],rev(ArcanoGenes$geneSymbol[1:10]))

CoryneGenes <- topTable(ebFit_micro, adjust ="BH", coef=3, number=40000, sort.by="logFC")
CoryneGenes <- CoryneGenes %>%
  rownames_to_column("geneSymbol") %>%
  filter(P.Value < 0.01) %>%
  arrange(desc(logFC))
CoryneGenes <- c(CoryneGenes$geneSymbol[1:10],rev(CoryneGenes$geneSymbol[1:10]))

StrepGenes <- topTable(ebFit_micro, adjust ="BH", coef=4, number=40000, sort.by="logFC")
StrepGenes <- StrepGenes %>%
  rownames_to_column("geneSymbol") %>%
  filter(P.Value < 0.01) %>%
  arrange(desc(logFC))
StrepGenes <- c(StrepGenes$geneSymbol[1:10],rev(StrepGenes$geneSymbol[1:10]))

MicroclusterGenes <- unique(c(StaphGenes, ArcanoGenes, CoryneGenes, StrepGenes))

topdiffGenes_micro <- diffGenes_micro[rownames(diffGenes_micro) %in% MicroclusterGenes,]
topclustRows_micro <- hclust(as.dist(1-cor(t(topdiffGenes_micro), method="pearson")), method="complete") 
topclustColumns_micro <- hclust(as.dist(1-cor(topdiffGenes_micro, method="spearman")), method="ward.D2") #cluster columns by spearman correlation


heatmap.2(topdiffGenes_micro, 
          Rowv=as.dendrogram(topclustRows_micro), 
          #Colv=as.dendrogram(topclustColumns_micro),
          ColSideColors=micro_colors2,
          col=rev(brewer.pal(name="BrBG", n=20)[-6]),
          scale='row', labCol = NA,
          density.info="none", trace="none",  
          cexRow=1.4, cexCol=1, margins=c(8,20))








# Multivariate linear regression analysis Microbiome clusters vs. clinical metadata ----
# Build data frame of associations:

mvlrm <- LeishOmics_gt3 %>%
  # Select variables:
  select(microbiome_cluster, age, sex, size_lesion_mm2, DTH_mm2, lymphadenopathy) %>%
  filter(microbiome_cluster != "NA") %>%
  
  # Modeling categorical variables:
  mutate(sex = case_when(sex %in% "female" ~ 1, TRUE ~ 2),
         lymphadenopathy = case_when(lymphadenopathy %in% "yes" ~ 1, TRUE ~ 2),
         microbiome_cluster = str_remove(microbiome_cluster, "M"))
mvlrm

# Build model:
model <- lm(microbiome_cluster ~ age + sex + size_lesion_mm2 + DTH_mm2 + lymphadenopathy, data = mvlrm)
summary(model)
write_tsv(as.data.frame(model$coefficients), "outputs/model.txt")


