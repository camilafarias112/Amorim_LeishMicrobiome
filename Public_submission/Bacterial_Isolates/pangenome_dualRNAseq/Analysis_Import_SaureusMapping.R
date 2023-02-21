# In this script, I bring the S. aureus pangenome mapping into my analysis to investigate differences between lesions
# with increased levels of reads mapped to the bacteria.

# Libraries ----
library(reshape2)
library(gplots)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(ggforce)
#library(limma)
library(ggbreak)
library(ggExtra)
library(tidyverse)


# Load color palette ----
#load("../RStudio_outputs/data/Dark24")

# Importing multiQC reports from the bowtie2 logs ----
mapping_afterkneaddata_saureus <- read_csv("mapping_files/QC/bowtie2_se_plot_multiQCreport.csv")

# dplyR data:
mapping_afterkneaddata_saureus <- mapping_afterkneaddata_saureus %>%
  mutate(total_reads = reads_mapped + `SE not aligned`,
         cpm_saureus = reads_mapped / total_reads * 1000000)

# Load study design, gene expression data and other vectors from main directory ----
load("../../microbiome_Lesion16S-seqDataset/outputs/phenotype_pre")
load("../RStudio_outputs/data/targets_all")
targets <- left_join(targets_all, mapping_afterkneaddata_saureus[,c("cpm_saureus","sample_RNAseq")],
                     by="sample_RNAseq")
targets <- left_join(targets, phenotype_pre)
#targets_info <- targets %>%
#  select(LTCP_patient_ID,cpm_saureus, microbiome_cluster)

write_tsv(targets %>%
            filter(treatment_other_drug == "sbv"), "outputs/tagerts_withSaureusCPM.txt")
#save(targets_info, file = "outputs/targets_info")
load("../RStudio_outputs/data/myDGEList.filtered.norm")
load("../RStudio_outputs/data/CL_samples")

# S. aureus CPM in all samples ----
# In log10:
targets %>%
  mutate(cpm_saureus = case_when(cpm_saureus == 0 ~ (cpm_saureus + 0.0000001), TRUE ~ cpm_saureus)) %>%
  mutate(cpm_saureus_log = log10(cpm_saureus)) %>%
  mutate(Clinical_Outcome = case_when(
    treatment_outcome %in% "cure" ~ "1 Sbv",
    treatment_outcome %in% "failure" ~ ">1 Sbv",
    TRUE ~ "DOT")) %>%
  filter(group_RNAseq != "NA") %>%
ggplot(., aes(x=factor(group_RNAseq, levels = c("HS", "CL")), y=cpm_saureus_log,
              color=group_RNAseq)) +
  geom_violin(size = 0.7, trim = T, adjust = 0.4) +
  geom_sina(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 17, vjust = 0.5), 
        axis.text.y = element_text(size = 17), axis.title = element_text(size = 17),
        legend.position="none", legend.text = element_text(size = 14), legend.title = element_blank()) +
  #geom_text_repel(aes(label = LTCP_patient_ID), size = 4, fontface="bold",
  #                max.overlaps = 7) +
  #geom_text_repel(aes(label = treatment_outcome), size = 2.5, fontface="bold",colour="black") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  #stat_compare_means(method = "wilcox.text", #label.y = 8.2,
  #                   aes(label = ..p.signif..),
  #                   label.x = 1.5, size = 5, color="black") +
  xlab("") + ylab(paste("S. aureus transcripts (log2CPM)")) +
  #geom_hline(yintercept = 1, color = Dark24[3], linetype=2, size=1) +
  #geom_hline(yintercept = 2, color = Dark24[14], linetype=2, size=1) +
  scale_y_break(c(-6.75, -1)) + 
  ylim(-7.5, 3.5)

targets %>%
  mutate(cpm_saureus = case_when(cpm_saureus == 0 ~ (cpm_saureus + 0.0000001), TRUE ~ cpm_saureus)) %>%
  mutate(cpm_saureus_log = log10(cpm_saureus)) %>%
  mutate(Clinical_Outcome = case_when(
    treatment_outcome %in% "cure" ~ "1 Sbv",
    treatment_outcome %in% "failure" ~ ">1 Sbv",
    TRUE ~ "DOT")) %>%
  filter(group_RNAseq != "NA") %>%
  filter(group_RNAseq == "CL") %>%
  ggplot(., aes(y=cpm_saureus_log,
                color=group_RNAseq)) +
  geom_line(stat = "density", adjust = .35, size=1.6) +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 17, vjust = 0.5), 
        axis.text.y = element_text(size = 17), axis.title = element_text(size = 17),
        legend.position="none", legend.text = element_text(size = 14), legend.title = element_blank()) +
  xlab("") + ylab(paste("S. aureus transcripts (log2CPM)")) +
  #geom_hline(yintercept = 1, color = Dark24[3], linetype=2, size=1) +
  #geom_hline(yintercept = 2, color = Dark24[14], linetype=2, size=1) +
  scale_y_break(c(-6.75, -1)) + 
  ylim(-7.5, 3.5)


targets %>%
  mutate(cpm_saureus = case_when(cpm_saureus == 0 ~ (cpm_saureus + 0.0000001), TRUE ~ cpm_saureus)) %>% # An adjustment was done
  #mutate(cpm_saureus_log = log10(cpm_saureus)) %>%
  mutate(Clinical_Outcome = case_when(
    treatment_outcome %in% "cure" ~ "1 Sbv",
    treatment_outcome %in% "failure" ~ ">1 Sbv",
    TRUE ~ "DOT")) %>%
  filter(group_RNAseq != "NA") %>%
  ggplot(., aes(x=factor(group_RNAseq, levels = c("HS", "CL")), y=cpm_saureus,
                color=group_RNAseq)) +
  geom_violin(size = 0.7, trim = F) +
  geom_sina(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 17, vjust = 0.5), 
        axis.text.y = element_text(size = 17), axis.title = element_text(size = 17),
        legend.position="none", legend.text = element_text(size = 14), legend.title = element_blank()) +
  #geom_text_repel(aes(label = LTCP_patient_ID), size = 4, fontface="bold",
  #                max.overlaps = 7) +
  #geom_text_repel(aes(label = treatment_outcome), size = 2.5, fontface="bold",colour="black") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  #stat_compare_means(method = "wilcox.text", #label.y = 8.2,
  #                   aes(label = ..p.signif..),
  #                   label.x = 1.5, size = 5, color="black") +
  xlab("") + ylab(paste("S. aureus transcripts (log2CPM)")) +
  geom_hline(yintercept = 1, color = Dark24[3], linetype=2, size=1) +
  geom_hline(yintercept = 2, color = Dark24[14], linetype=2, size=1) +
  #ylim(1, 5) +
  NULL
write_tsv(targets %>%
            mutate(cpm_saureus = case_when(cpm_saureus == 0 ~ (cpm_saureus + 0.0000001), TRUE ~ cpm_saureus)) %>%
            select(LTCP_patient_ID, group_RNAseq, cpm_saureus),
          "../../../../../Desktop/targets_for_SaureusCPM.txt")


targets %>%
  mutate(cpm_saureus = case_when(cpm_saureus == 0 ~ (cpm_saureus + 0.0000001), TRUE ~ cpm_saureus)) %>%
  mutate(cpm_saureus_log = log10(cpm_saureus)) %>%
  mutate(Clinical_Outcome = case_when(
    treatment_outcome %in% "cure" ~ "1 Sbv",
    treatment_outcome %in% "failure" ~ ">1 Sbv",
    TRUE ~ "DOT")) %>%
  filter(group_RNAseq != "NA") %>%
  ggplot(., aes(x=factor(group_RNAseq, levels = c("HS", "CL")), y=cpm_saureus_log,
                color=group_RNAseq)) +
  geom_violin(size = 0.7, trim = F) +
  geom_sina(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 17, vjust = 0.5), 
        axis.text.y = element_text(size = 17), axis.title = element_text(size = 17),
        legend.position="none", legend.text = element_text(size = 14), legend.title = element_blank()) +
  #geom_text_repel(aes(label = LTCP_patient_ID), size = 4, fontface="bold",
  #                max.overlaps = 7) +
  #geom_text_repel(aes(label = treatment_outcome), size = 2.5, fontface="bold",colour="black") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  #stat_compare_means(method = "wilcox.text", #label.y = 8.2,
  #                   aes(label = ..p.signif..),
  #                   label.x = 1.5, size = 5, color="black") +
  xlab("") + ylab(paste("S. aureus transcripts (log2CPM)")) +
  geom_hline(yintercept = 1, color = Dark24[3], linetype=2, size=1) +
  geom_hline(yintercept = 2, color = Dark24[14], linetype=2, size=1) +
  ylim(-8, 5)


# Microbiome clusters:
my_comparisons <- list(M4 = c("M0","M4"),
                       #M5 = c("M0","M5"),
                       M6 = c("M0","M6"),
                       M7 = c("M0","M7"))
targets %>%
  mutate(cpm_saureus = case_when(cpm_saureus == 0 ~ (cpm_saureus + 0.0000001), TRUE ~ cpm_saureus)) %>%
  mutate(cpm_saureus_log = log10(cpm_saureus)) %>%
  mutate(Clinical_Outcome = case_when(
    treatment_outcome %in% "cure" ~ "1 Sbv",
    treatment_outcome %in% "failure" ~ ">1 Sbv",
    TRUE ~ "DOT")) %>%
  filter(microbiome_cluster != "NA") %>%
  filter(group_RNAseq == "CL") %>%
  mutate(microbiome_cluster2 = case_when(
    microbiome_cluster == "Staphylococcus" ~ "M6",
    microbiome_cluster == "Arcanobacterium" ~ "M7",
    microbiome_cluster == "Streptococcus" ~ "M4",
    microbiome_cluster == "Heterogeneous" ~ "M0",
    microbiome_cluster == "Corynebacterium" ~ "M5",
    microbiome_cluster == "Finegoldia" ~ "M2",
    microbiome_cluster == "Lactobacillus" ~ "M1",
    TRUE ~ "HS")) %>%
  ggplot(., aes(x=microbiome_cluster2, y=cpm_saureus_log,
                color=microbiome_cluster2)) +
  geom_violin(size = 0.7, trim = T, adjust = 0.4) +
  geom_sina(size=3) +
  theme_classic() +
  scale_color_manual(values = c("gray",
                                Dark24[5],
                                Dark24[6],
                                Dark24[3],
                                Dark24[2],
                                Dark24[1],
                                Dark24[4])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 17), axis.title = element_text(size = 17),
        legend.position="right", legend.text = element_text(size = 17), legend.title = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  xlab("") + ylab(paste("S. aureus transcripts (log2CPM)")) +
  scale_y_break(c(-6.75, -1)) + 
  ylim(-7.5, 3.5)


export_saureus_microclusters <- targets %>%
  mutate(cpm_saureus = case_when(cpm_saureus == 0 ~ (cpm_saureus + 0.0000001), TRUE ~ cpm_saureus)) %>%
  mutate(cpm_saureus_log = log10(cpm_saureus)) %>%
  mutate(Clinical_Outcome = case_when(
    treatment_outcome %in% "cure" ~ "1 Sbv",
    treatment_outcome %in% "failure" ~ ">1 Sbv",
    TRUE ~ "DOT")) %>%
  filter(group_RNAseq != "NA",
         treatment_other_drug == "sbv") %>%
  select(LTCP_patient_ID, healing_time_days, treatment_outcome, cpm_saureus,cpm_saureus_log, microbiome_cluster)

write_tsv(export_saureus_microclusters, "saureus_microclusters.txt")


# Percentiles:
targets %>%
  filter(group_RNAseq == "CL") %>%
  ggplot(.,aes(x=group_RNAseq, y=cpm_saureus)) +
  geom_point(aes(fill=group_RNAseq), size=3, shape=21) + 
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
  geom_text_repel(aes(label = LTCP_patient_ID), size = 2, nudge_x = 0.5, max.overlaps = 100) +
  xlab("") + ylab(paste("Failure gene signature ES")) +
  #ylim(0,10) +
  NULL



high_Saureus <- 
targets %>%
  #left_join(phenotype, by= "LTCP_patient_ID") %>%
  #filter(microbiome_cluster != "NA") %>%
  filter(cpm_saureus >10) %>% 
  arrange(desc(cpm_saureus))
write_tsv(high_Saureus, "../../../../../Desktop/high_Saureus.txt")

# Export S. aureus transcript values ----
save(targets_CL, file = "outputs/targets_CL")

# Add categorical vector to parasite and bacteria CPMs ----
# S. aureus
#targets_w_Gricedata <- targets_w_Gricedata %>%
#  mutate(cpm_saureus_cat = case_when(
#    cpm_saureus > .05*10000 ~ "Saureus_detec",
#    cpm_saureus < .05*10000 ~ "Saureus_NOTdetec"))
#View(targets_w_Gricedata %>% select(LTCP_patient_ID, cpm_saureus_cat, cpm_saureus))

targets <- targets %>%
  mutate(cpm_saureus_cat = case_when(
    cpm_saureus > 10 ~ "Saureus_detec",
    cpm_saureus < 10 ~ "Saureus_NOTdetec"))

Saureus_detec_samples <- targets %>% filter(cpm_saureus_cat == "Saureus_detec")
Saureus_detec_samples <- Saureus_detec_samples$sample_RNAseq

Saureus_nondetec_samples <- targets %>% filter(cpm_saureus_cat == "Saureus_NOTdetec")
Saureus_nondetec_samples <- Saureus_nondetec_samples$sample_RNAseq

#View(targets %>% select(sample_RNAseq, cpm_saureus_cat, cpm_saureus))

# Association between parasite and bacteria CPM and clinical metadata ----
targets %>%
  mutate(cpm_saureus = case_when(cpm_saureus == 0 ~ (cpm_saureus + 0.0000001), TRUE ~ cpm_saureus)) %>%
  mutate(cpm_saureus_log = log10(cpm_saureus)) %>%
  mutate(Clinical_Outcome = case_when(
    treatment_outcome %in% "cure" ~ "1 Sbv",
    treatment_outcome %in% "failure" ~ ">1 Sbv",
    TRUE ~ "DOT")) %>%
  filter(group_RNAseq != "NA",
         Clinical_Outcome != "DOT") %>%
  ggplot(., aes(x=factor(Clinical_Outcome, levels = c("1 Sbv", ">1 Sbv")), y=cpm_saureus_log,
                color=Clinical_Outcome)) +
  geom_violin(size = 0.7, trim = F) +
  geom_sina(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 17, vjust = 0.5), 
        axis.text.y = element_text(size = 17), axis.title = element_text(size = 17),
        legend.position="none", legend.text = element_text(size = 14), legend.title = element_blank()) +
  #geom_text_repel(aes(label = LTCP_patient_ID), size = 4, fontface="bold",
  #                max.overlaps = 7) +
  #geom_text_repel(aes(label = treatment_outcome), size = 2.5, fontface="bold",colour="black") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "wilcox.test", #label.y = 8.2,
                     aes(label = ..p.signif..),
                     label.x = 1.5, size = 5, color="black") +
  xlab("") + ylab(paste("S. aureus transcripts (log2CPM)")) +
  geom_hline(yintercept = 1, color = Dark24[3], linetype=2) +
  ylim(-8, 5)

targets %>%
  mutate(cpm_saureus = case_when(cpm_saureus == 0 ~ (cpm_saureus + 0.0000001), TRUE ~ cpm_saureus)) %>%
  mutate(cpm_saureus_log = log10(cpm_saureus)) %>%
  mutate(Clinical_Outcome = case_when(
    treatment_outcome %in% "cure" ~ "1 Sbv",
    treatment_outcome %in% "failure" ~ ">1 Sbv",
    TRUE ~ "DOT")) %>%
  filter(group_RNAseq != "NA",
         Clinical_Outcome != "DOT") %>%
  ggplot(., aes(x=healing_time_days, y=cpm_saureus_log,
                color=Clinical_Outcome)) +
  geom_smooth(method=lm, color=Dark24[1]) +
  geom_point(size=2, color= "#424242") +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 17, vjust = 0.5), 
        axis.text.y = element_text(size = 17), axis.title = element_text(size = 17),
        legend.position="none", legend.text = element_text(size = 14), legend.title = element_blank()) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE, color=Dark24[1]) + 
  xlab("") + ylab(paste("S. aureus transcripts (log2CPM)")) +
  geom_hline(yintercept = 1, color = Dark24[3], linetype=2) +
  ylim(-8, 5)


# Parasite loads and S. aureus CPM ----
load("../../3rd_LesionRNAseq_Dataset/qPCRs/Lbraziliensis_qPCR_ImportingDataIntoR/parasiteload_3rdDataset")
targets %>%
  left_join(parasiteload_3rdDataset, by = "sample_RNAseq") %>%
  filter(cpm_saureus_cat == "Saureus_detec") %>%
  ggplot(., aes(x=cpm_saureus, y=Lbraziliensis_mg)) +
  geom_point(size=2, color= "#424242") +
  theme_classic() + #scale_color_manual(values = PalC) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11),
        legend.position="right", strip.background = element_blank(), strip.text = element_text(size = 10, face="bold"),
        axis.text = element_blank(), axis.ticks = element_blank()) +
  stat_cor(method = "pearson",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE, color=Dark24[1])



# DE S. aureus detec vs S. aureus NONdetec ----
targets_CL <- targets %>% filter(group_RNAseq == "CL")
table(targets_CL$cpm_saureus_cat) # number of samples for this comparison
design8 <- model.matrix(~0 + factor(targets_CL$cpm_saureus_cat)) 
colnames(design8) <- levels(factor(targets_CL$cpm_saureus_cat))

# A dimentionally reduction happens here: Only working with VITALS!
load("../Unsup_geneExpr_DataReduction/viz_var_up")

v.DEGList.filtered.norm8 <- voom(myDGEList.filtered.norm[viz_var_up,
                                                         CL_samples], design8)
fit8 <- lmFit(v.DEGList.filtered.norm8, design8)
contrast.matrix8 <- makeContrasts(comparison = Saureus_detec - Saureus_NOTdetec,
                                  levels=design8)
fits8 <- contrasts.fit(fit8, contrast.matrix8)
ebFit8 <- eBayes(fits8)
myTopHits8 <- topTable(ebFit8, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits8 <- as_tibble(myTopHits8, rownames = "geneSymbol")

myTopHits8 <- myTopHits8 %>%
  mutate(FC = 2^logFC)
sig_genes_up8 <- myTopHits8 %>%
  filter(P.Value < 0.05) %>% filter(FC > 1.5)
sig_genes_down8 <- myTopHits8 %>%
  filter(P.Value < 0.05) %>% filter(FC < 1.5)
sig_genes_up8 <- sig_genes_up8$geneSymbol
sig_genes_down8 <- sig_genes_down8$geneSymbol

save(sig_genes_up8, file = "outputs/sig_genes_up8")
write_tsv(as_tibble(sig_genes_up8), "../RStudio_outputs/data/Staph_detect_sig_genes_up8_FC1.5.txt")
write_tsv(as_tibble(sig_genes_down8), "../RStudio_outputs/data/Staph_detect_sig_genes_down8_FC1.5.txt")


write_tsv(myTopHits8, "../../../../../Desktop/SupplementalTable9_DEGs_Saureus.txt")
#sig_genes_up_cpm_mapping <- myTopHits8 %>%
#  filter(P.Value < 0.05) %>% filter(logFC > 1)
#sig_genes_down_cpm_mapping <- myTopHits8 %>%
#  filter(P.Value < 0.05) %>% filter(logFC < -1)
#sig_genes_up_cpm_mapping <- sig_genes_up_cpm_mapping$geneSymbol
#sig_genes_down_cpm_mapping <- sig_genes_down_cpm_mapping$geneSymbol
#save(sig_genes_up_cpm_mapping, file = "RStudio_outputs/data/sig_genes_up_cpm_mapping")

myTopHits8$col0 <- myTopHits8$geneSymbol
col0_v8 <- myTopHits8$col0 %in% sig_genes_up8
col0_v88 <- myTopHits8$col0 %in% sig_genes_down8
myTopHits8$col0[!col0_v8 & !col0_v88] <- NA

myTopHits8$col00 <- myTopHits8$geneSymbol
myTopHits8$col00[col0_v8] <- "up"
myTopHits8$col00[col0_v88] <- "down"
myTopHits8$col00[!col0_v8 & !col0_v88] <- "other"

# Top genes ----
sig_genes_up_top <- myTopHits8 %>%
  filter(P.Value < 0.05) %>%
  arrange(desc(FC), P.Value) %>%
  top_n(7, FC)

sig_genes_up_top2 <- myTopHits8 %>%
  filter(P.Value < 0.05 & logFC > 0) %>%
  arrange(P.Value, desc(FC)) %>%
  top_n(5, desc(P.Value))
  
sig_genes_down_top <- myTopHits8 %>%
  filter(P.Value < 0.05) %>%
  arrange(FC, P.Value) %>%
  top_n(-10, FC)

sig_genes_up_top3 <- myTopHits8 %>%
  filter(P.Value < 0.05 & FC > 1.5) %>%
  arrange(P.Value, desc(FC))

sig_genes_up_top <- sig_genes_up_top$geneSymbol
sig_genes_down_top <- sig_genes_down_top$geneSymbol
sig_genes_up_top2 <- sig_genes_up_top2$geneSymbol
sig_genes_up_top3 <- sig_genes_up_top3$geneSymbol

myTopHits8$col10 <- myTopHits8$geneSymbol
col10_v8 <- myTopHits8$col10 %in% sig_genes_up_top
col10_v88 <- myTopHits8$col10 %in% sig_genes_down_top
myTopHits8$col10[!col10_v8 & !col10_v88] <- NA

myTopHits8$col11 <- myTopHits8$geneSymbol
myTopHits8$col11[col10_v8] <- "up"
myTopHits8$col11[col10_v88] <- "down"
myTopHits8$col11[!col10_v8 & !col10_v88] <- "other"

myTopHits8$col12 <- myTopHits8$geneSymbol
col10_v82 <- myTopHits8$col12 %in% c(sig_genes_up_top,sig_genes_up_top2)
myTopHits8$col12[!col10_v82] <- NA

myTopHits8$col13 <- myTopHits8$geneSymbol
myTopHits8$col13[col10_v82] <- "tops"
myTopHits8$col13[!col10_v82] <- "other"

myTopHits8$col14 <- myTopHits8$geneSymbol
col10_v824 <- myTopHits8$col14 %in% sig_genes_up_top3
myTopHits8$col14[!col10_v824] <- NA

myTopHits8$col15 <- myTopHits8$geneSymbol
myTopHits8$col15[col10_v824] <- "tops"
myTopHits8$col15[!col10_v824] <- "other"

ggplot(myTopHits8 %>%
         arrange(desc(col10)), aes(y=-log10(P.Value), x=FC, color=col11, size=col11)) +
  geom_point() +
  theme_classic() + scale_color_manual(values = c(Dark24[2],"gray",Dark24[1])) +
  scale_size_manual(values = c(2, 0.05,2)) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  geom_text_repel(aes(label = col10, color=col11), size = 5, fontface="bold",
                  max.overlaps = 50) +

  xlab("logFC S. aureus + vs. -")

ggplot(myTopHits8 %>%
         arrange(desc(col10)), aes(y=-log10(P.Value), x=FC)) +
  geom_point(color="gray") +
  theme_classic() +
  scale_color_manual(values = c("dark gray", "dark gray", Dark24[1])) +
  scale_size_manual(values = c(0.5,1,1)) +
  geom_text_repel(aes(label = geneSymbol), size = 3.2, fontface="bold",color="black",
                  max.overlaps = 14) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13), panel.border = element_rect(fill = NA, color = Dark24[1], size = 2)) +
  xlim(1.5, 3.6) + ylim(-log10(0.05),6) +
  xlab("FC S. aureus + vs. -")

ggplot(myTopHits8, aes(y=-log10(P.Value), x=FC,
                      color=col00, size=col15, alpha=col15)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("dark gray", "dark gray", Dark24[21])) +
  scale_size_manual(values = c(2,3)) +
  scale_alpha_manual(values = c(0.7,1)) +
  geom_text_repel(aes(label = col14), size = 3.8, fontface=4, color="black", max.overlaps = 50) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = -log10(0.05)) + geom_vline(xintercept = 1.5) +
  #annotate("text", x=3, y=4.5, label=paste(length(sig_genes_up8)), size=5, color=Dark24[1], fontface="bold") +
  xlab("FC S. aureus + vs. -")

ggplot(myTopHits8, aes(y=-log10(P.Value), x=FC,
                       color=col00, size=col15, alpha=col15)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("dark gray", "dark gray", Dark24[21])) +
  scale_size_manual(values = c(2,4)) +
  scale_alpha_manual(values = c(0.7,1)) +
  #geom_text_repel(aes(label = col14), size = 4, fontface=4, color="black") +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = -log10(0.05)) + geom_vline(xintercept = 1.5) +
  #annotate("text", x=3, y=4.5, label=paste(length(sig_genes_up8)), size=5, color=Dark24[1], fontface="bold") +
  xlab("FC S. aureus + vs. -")

ggplot(myTopHits8, aes(y=-log10(P.Value), x=FC,
                       color=col00, size=col15, alpha=col15)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("dark gray", "dark gray", Dark24[21])) +
  scale_size_manual(values = c(2,4)) +
  scale_alpha_manual(values = c(0.7,1)) +
  geom_text_repel(aes(label = col0), size = 4, fontface=4, color="black",
                  max.overlaps = 50) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = -log10(0.05)) + geom_vline(xintercept = 1.5) +
  #annotate("text", x=3, y=4.5, label=paste(length(sig_genes_up8)), size=5, color=Dark24[1], fontface="bold") +
  xlab("FC S. aureus + vs. -") + xlim(1, 4)


# Add to volcano the GO genesets  (HAVE TO UPDATE, not been used) ----
#GO_sig_up_SaureusDetec <- read_delim("../GSEA_GO_outputs/GO/GO_StaphyDetect_sig_genes_up8_FC1.5_modToImport_pangenomeSaureus.txt", 
#                                     "\t", escape_double = FALSE,
#                                     trim_ws = TRUE)
#geneset_clus1 <- GO_sig_up_SaureusDetec$cluster1
#geneset_clus1 <- geneset_clus1[!duplicated(geneset_clus1)]
#geneset_clus1 <- na.omit(geneset_clus1)

#geneset_clus2 <- GO_sig_up_SaureusDetec$cluster2 
#geneset_clus2 <- geneset_clus2[!duplicated(geneset_clus2)]
#geneset_clus2 <- na.omit(geneset_clus2)

#geneset_clus3 <- GO_sig_up_SaureusDetec$cluster3 
#geneset_clus3 <- geneset_clus3[!duplicated(geneset_clus3)]
#geneset_clus3 <- na.omit(geneset_clus3)

#geneset_clus4 <- GO_sig_up_SaureusDetec$cluster4
#geneset_clus4 <- geneset_clus4[!duplicated(geneset_clus4)]
#geneset_clus4 <- na.omit(geneset_clus4)

#allgeneset_clus <- c(geneset_clus1,geneset_clus2,geneset_clus3,geneset_clus4)
#allgeneset_clus <- allgeneset_clus[!duplicated(allgeneset_clus)]

# Modeling GO dataframe:
#GO_sig_up_SaureusDetec

#data.frame(GeneCount = c(length(geneset_clus1),
#                         length(geneset_clus2),
#                         length(geneset_clus3),
#                         length(geneset_clus4)),
#           AnnotationCluster = c("Chemotaxis",
#                                 "Lymphoid and Non-lymphoid cell Interactions",
#                                 "Infectious and auto-immune diseases",
#                                 "Epidermis development, cellular adhesion and ECM")) %>%
#  ggplot(., aes(x=AnnotationCluster, y=GeneCount, fill=AnnotationCluster)) +
#  geom_bar(width = 1, stat = "identity") +
#  theme_classic() +
#  scale_fill_manual(values = c(Dark24[4],Dark24[1],Dark24[2],Dark24[6],"gray",Dark24[8])) +
#  theme(legend.position = "none", axis.title =element_text(size=20),
#        axis.title.x = element_blank(), axis.text = element_text(size=18),
#        axis.ticks.length = unit(0, "cm")) +
#  coord_flip() +
#  NULL

# Annotations were given according to the main pathway in the cluster - This CLUSTERING HERE DOESN'T MAKE SENSE
#venn_sigup8 <- venn(list(Cluster1 = geneset_clus1,
#                         Cluster2 = geneset_clus2,
#                         Cluster3 = geneset_clus3,
#                         Cluster4 = geneset_clus4))

#myTopHits8$col3 <- myTopHits8$geneSymbol
#col3_v8_2 <- myTopHits8$col3 %in% attr(venn_sigup8,"intersections")$Cluster1
#col3_v8_22 <- myTopHits8$col3 %in% attr(venn_sigup8,"intersections")$Cluster2
#col3_v8_222 <- myTopHits8$col3 %in% attr(venn_sigup8,"intersections")$Cluster3
#col3_v8_2222 <- myTopHits8$col3 %in% attr(venn_sigup8,"intersections")$Cluster4
#col3_v8_22222 <- myTopHits8$col3 %in% c(attr(venn_sigup8,"intersections")$`Cluster1:Cluster2`,
#                                        attr(venn_sigup8,"intersections")$`Cluster1:Cluster3`,
#                                        attr(venn_sigup8,"intersections")$`Cluster1:Cluster4`,
#                                        attr(venn_sigup8,"intersections")$`Cluster1:Cluster2:Cluster4`)


#myTopHits8$col3[col3_v8_2] <- "Chemotaxis"
#myTopHits8$col3[col3_v8_22] <- "Lymphoid and a non-Lymphoid cell interactions"
#myTopHits8$col3[col3_v8_222] <- "Infectious and auto-immune diseases"
#myTopHits8$col3[col3_v8_2222] <- "Epidermis development, cellular adhesion and ECM"
#myTopHits8$col3[col3_v8_22222] <- "Shared"
#myTopHits8$col3[!col3_v8_2 & !col3_v8_22 & !col3_v8_222 & !col3_v8_2222 & !col3_v8_22222] <- "others"

#myTopHits8$col4 <- myTopHits8$geneSymbol
#myTopHits8$col4[!col3_v8_2 & !col3_v8_22 & !col3_v8_222 & !col3_v8_2222 & !col3_v8_22222] <- NA

# Picking only selected GO genesets:
#myTopHits8$col5 <- myTopHits8$geneSymbol
#myTopHits8$col5[!col3_v8_2 & !col3_v8_22222] <- NA

# Only clusters:
#`Cluster1:Cluster3:Cluster4` and attr(,"intersections")$Cluster3

#myTopHits8$col6 <- myTopHits8$geneSymbol
#col6_v <- myTopHits8$col6 %in% attr(venn_sigup8,"intersections")$`Cluster1:Cluster3:Cluster4`
#col6_v2 <- myTopHits8$col6 %in% attr(venn_sigup8,"intersections")$Cluster3
#myTopHits8$col6[col6_v] <- "Cluster 1/3/4"
#myTopHits8$col6[col6_v2] <- "Cluster 3"
#myTopHits8$col6[!col6_v & !col6_v2] <- "Others"

#myTopHits8$col7 <- myTopHits8$geneSymbol
#myTopHits8$col7[!col6_v & !col6_v2] <- NA

# Only Cluster 3:
#myTopHits8$col8 <- myTopHits8$geneSymbol
#col8_v <- myTopHits8$col8 %in% geneset_clus3
#myTopHits8$col8[col8_v] <- "Cluster 3"
#myTopHits8$col8[!col8_v] <- "Others"

#myTopHits8$col9 <- myTopHits8$geneSymbol
#myTopHits8$col9[!col8_v] <- NA

# Specific genes to show ----
genes_toshow1 <- c("IL1A","IL1B","CXCL5","CXCL8", "GNLY","GZMB","PRF1")

genes_toshow2 <- c("CXCL8","CXCL5","CXCL1","CXCL2",
                   "CXCR1")

genes_toshow3 <- c("CCL3L3","CCL4L2","CCL3","CCL4","CCL7","CCL8")

load("../RStudio_outputs/data/myTopHits_treat_var_imp")
genes_toshow3 <- myTopHits_treat_var_imp

myTopHits8$col6 <- myTopHits8$geneSymbol
col6_v8_1 <- myTopHits8$col6 %in% genes_toshow1
col6_v8_2 <- myTopHits8$col6 %in% genes_toshow2
col6_v8_3 <- myTopHits8$col6 %in% genes_toshow3
myTopHits8$col6[!col6_v8_1 & !col6_v8_2 & !col6_v8_3] <- NA

myTopHits8$col6 <- myTopHits8$geneSymbol
col6_v8_3 <- myTopHits8$col6 %in% genes_toshow1
myTopHits8$col6[!col6_v8_3] <- NA

myTopHits8$col8 <- myTopHits8$geneSymbol
myTopHits8$col8[!col6_v8_3] <- "notshow"
myTopHits8$col8[col6_v8_3] <- "show"

myTopHits8$col7 <- myTopHits8$geneSymbol
myTopHits8$col7[col6_v8_1] <- "IL1s"
myTopHits8$col7[col6_v8_2] <- "CXCs"
myTopHits8$col7[col6_v8_3] <- "CCLs"
myTopHits8$col7[!col6_v8_1 & !col6_v8_2 & !col6_v8_3] <- "others"

ggplot(myTopHits8 %>%
         arrange(col8), aes(y=-log10(P.Value), x=FC,
                       color=col8, size=col8, alpha=col8)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("gray", Dark24[1])) +
  scale_size_manual(values = c(1,2)) +
  scale_alpha_manual(values = c(0.7,1)) +
  geom_text_repel(aes(label = col6), size = 4.2, fontface=4,
                  max.overlaps = 50) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13)) +
  #geom_hline(yintercept = -log10(0.05)) + geom_vline(xintercept = 1.5) +
  #annotate("text", x=3, y=4.5, label=paste(length(sig_genes_up8)), size=5, color=Dark24[1], fontface="bold") +
  xlab("FC S. aureus + vs. -")

# Show individual gene expression Violins ----
genespick <- c("IL1A","IL1B","IFNG")
gene_viz <- log2.cpm.filtered.norm[genespick,CL_samples]

melt_viz <- as.data.frame(melt(gene_viz))
melt_viz$Var2 <- as.character(melt_viz$Var2)
colnames(melt_viz)[2] <- "sample_RNAseq"

#my_comparisons1 <- list(High = c("High","S. aureus -"),
#                       #M5 = c("M0","M5"),
#                       Intermediate = c("Intermediate","S. aureus -"),
#                       High_Inter = c("High","Intermediate"))

melt_viz %>%
  left_join(targets %>%
              mutate(cpm_saureus = case_when(cpm_saureus == 0 ~ (cpm_saureus + 0.0000001), TRUE ~ cpm_saureus)) %>%
              mutate(cpm_saureus_log = log10(cpm_saureus)) %>%
              mutate(Clinical_Outcome = case_when(
                treatment_outcome %in% "cure" ~ "1 Sbv",
                treatment_outcome %in% "failure" ~ ">1 Sbv",
                TRUE ~ "DOT")) %>%
              filter(group_RNAseq == "CL"), by="sample_RNAseq") %>%
  mutate(cpm_saureus_log_cat = case_when(
    cpm_saureus_log > 2 ~ "High",
    cpm_saureus_log > 1 ~ "Intermediate",
    cpm_saureus_log < 1 ~ "S. aureus -"
  )) %>%
  select(Var1, sample_RNAseq, value, cpm_saureus_log_cat, cpm_saureus_log) %>%
  ggplot(., aes(x=factor(cpm_saureus_log_cat, levels = c("S. aureus -","Intermediate","High")),
                y=value,
                color=cpm_saureus_log_cat)) +
  geom_violin(size = 0.7, trim = F) +
  geom_sina(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 17, vjust = 1, hjust = 1,angle = 45), 
        axis.text.y = element_text(size = 17), axis.title = element_text(size = 17),
        legend.position="none", legend.text = element_text(size = 14), legend.title = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 14)) +
  #geom_text_repel(aes(label = LTCP_patient_ID), size = 4, fontface="bold",
  #                max.overlaps = 7) +
  #geom_text_repel(aes(label = treatment_outcome), size = 2.5, fontface="bold",colour="black") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "wilcox.test", #label.y = 8.2,
                     aes(label = ..p.signif..), comparisons = my_comparisons1,
                     label.x = 1.5, size = 5, color="black") +
  xlab("") + ylab(paste("")) +
  facet_wrap(. ~ Var1, scales = "free") +
  scale_y_break(c(-6.75, -1))


test <- melt_viz %>%
  left_join(targets %>%
              mutate(cpm_saureus = case_when(cpm_saureus == 0 ~ (cpm_saureus + 0.0000001), TRUE ~ cpm_saureus)) %>%
              mutate(cpm_saureus_log = log10(cpm_saureus)) %>%
              mutate(Clinical_Outcome = case_when(
                treatment_outcome %in% "cure" ~ "1 Sbv",
                treatment_outcome %in% "failure" ~ ">1 Sbv",
                TRUE ~ "DOT")) %>%
              filter(group_RNAseq == "CL"), by="sample_RNAseq") %>%
  mutate(cpm_saureus_log_cat = case_when(
    cpm_saureus_log > 2 ~ "High",
    cpm_saureus_log > 1 ~ "Intermediate",
    cpm_saureus_log < 1 ~ "S. aureus -"
  )) %>%
  select(Var1, sample_RNAseq, value, cpm_saureus_log_cat, cpm_saureus_log, microbiome_cluster, treatment_outcome,)
write_tsv(test, "../../../../../Desktop/text.txt")


# Fold change pick genes ----
myTopHits8 %>%
  filter(geneSymbol %in% genespick) %>%
  ggplot(., aes(x=2^logFC, y=geneSymbol)) +
  geom_bar(stat = "identity")

# Correlation networks for healing time and gene expression ----
targets_sbv <- targets_CL %>% filter(treatment_other_drug == "sbv")

library(broom)
library(patchwork)
library(corrplot)
library(Hmisc)
library(tidygraph)
library(ggraph)
library(corrr)

load("../RStudio_outputs/data/log2.cpm.filtered.norm")

healing_variable <- "Healing time"
healing_corr_mat <- rbind(log2.cpm.filtered.norm[sig_genes_up8,targets_sbv$sample_RNAseq],targets_sbv$healing_time_days)
rownames(healing_corr_mat)[dim(healing_corr_mat)[1]] <- healing_variable

res.cor_heal <- t(healing_corr_mat) %>% correlate(method = "pearson") %>%
  shave(upper = TRUE) %>%
  stretch(na.rm = TRUE) %>%
  filter(y == healing_variable) %>%
  arrange(desc(r))# %>%
  #top_n(1:25)

cor.graph_heal <- as_tbl_graph(res.cor_heal, directed = FALSE)

cor.graph_heal <- cor.graph_heal %>% 
  activate(nodes) %>% 
  rename(Feature = name) %>%
  #left_join(description, by = "Gene") %>%
  mutate(Feature2 = case_when(
    Feature %in% healing_variable ~ "healing",
    Feature %ni% healing_variable ~ "genes")) %>%
  activate(edges) %>%
  arrange(desc(r)) #%>%
#mutate(direction = case_when( #there aren't negative correlations with spearman < -.70
#  r > 0 ~ "positive",
#  r < 0 ~ "negative"))

cor.graph_heal
cor.graph_heal %>%
  ggraph(., layout = "kk") +
  #geom_edge_density(aes(fill = r)) + 
  geom_edge_link(aes(width = r,edge_colour = r)) +
  #scale_edge_color_manual(values = c("blue", "red")) +
  scale_edge_color_gradient(low = "white", high = Dark24[1]) +
  scale_edge_width(range = c(0.2, 2)) +
  #geom_node_point(aes(color = family, shape = impgenes), size = 3) +
  geom_node_point(aes(color=Feature2), size = 2) +
  #scale_shape_manual(values = c(15,16)) +
  scale_color_manual(values = c("black",Dark24[15])) +
  geom_node_text(aes(label = Feature, color=Feature2), size = 4, fontface="bold",
                 max.overlaps = 20, repel = TRUE) +
  theme_graph()

# Individual correlations:
corr_transp_heal <- rownames_to_column(as.data.frame(t(healing_corr_mat[c(res.cor_heal$x[1:25],healing_variable),])), "sample_RNAseq")

corr_transp_heal %>%
  pivot_longer(!sample_RNAseq & !`Healing time`, names_to = "feature", values_to = "value") %>%
  ggplot(., aes(x=value, y=`Healing time`)) +
  geom_smooth(method=lm, color=Dark24[1]) +
  geom_point(size=2, color= "#424242") +
  theme_classic() + #scale_color_manual(values = PalC) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11),
        legend.position="right", strip.background = element_blank(), strip.text = element_text(size = 10, face="bold"),
        axis.text = element_blank(), axis.ticks = element_blank()) +
  stat_cor(method = "pearson",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE, color=Dark24[1]) + 
  xlab("Gene expression (log2CPM)") + ylab("Healing time") +
  facet_wrap(. ~ factor(feature, levels = res.cor_heal$x[1:30]), nrow = 3, scales = "free")

# corrplot:
genehealing_stats <- rcorr(as.matrix(t(healing_corr_mat)), type = "pearson")

genehealing_stats_tidy <- tidy(genehealing_stats)
genehealing_stats_tidy <- genehealing_stats_tidy %>%
  filter(column1 == healing_variable)

library(paletteer)
colormat <- paletteer_c("ggthemes::Classic Red-Blue", 30)
corrplot(genehealing_stats$r, type="upper", #order="hclust", 
         p.mat = genehealing_stats$P, sig.level = 0.05, insig = "blank",
         col = rev(colormat),tl.col = "black",tl.cex = 0.7,
         pch.cex = 1, pch.col = "white")

# Volcano S. aureus detect with biomarkers for treatment outcome in 3rd lesion dataset ----
load("../RStudio_outputs/data/biomarkers_3rd")
myTopHits8$col14 <- myTopHits8$geneSymbol
col14_v <- myTopHits8$col14 %in% biomarkers_3rd
myTopHits8$col14[!col14_v] <- NA

myTopHits8$col15 <- myTopHits8$geneSymbol
myTopHits8$col15[col14_v] <- "associated_failure"
myTopHits8$col15[!col14_v] <- "other"

ggplot(myTopHits8 %>%
         arrange(desc(col10)), aes(y=-log10(P.Value), x=FC, color=col15)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c(Dark24[1],"dark gray")) +
  geom_text_repel(aes(label = geneSymbol), size = 3.2, fontface="bold",
                  max.overlaps = 14) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13), panel.border = element_rect(fill = NA, color = Dark24[1], size = 2)) +
  xlim(1.5, 3.6) + ylim(-log10(0.05),6) +
  xlab("FC S. aureus + vs. -")

# Correlation networks for cpm_aureus and gene expression ----
Saureus_variable <- "S.aureus CPM abundance"
Saureus_corr_mat <- rbind(log2.cpm.filtered.norm[sig_genes_up8,targets_CL_Saureus$sample_RNAseq],targets_CL_Saureus$cpm_saureus)
rownames(Saureus_corr_mat)[dim(Saureus_corr_mat)[1]] <- Saureus_variable

res.cor <- t(Saureus_corr_mat) %>% correlate(method = "pearson") %>%
  shave(upper = TRUE) %>%
  stretch(na.rm = TRUE) %>%
  filter(y == Saureus_variable) %>%
  arrange(desc(r)) %>%
  top_n(1:30)

cor.graph <- as_tbl_graph(res.cor, directed = FALSE)

cor.graph <- cor.graph %>% 
  activate(nodes) %>% 
  rename(Feature = name) %>%
  #left_join(description, by = "Gene") %>%
  mutate(Feature2 = case_when(
    Feature %in% Saureus_variable ~ "main",
    Feature %ni% Saureus_variable ~ "genes")) %>%
  activate(edges) %>%
  arrange(desc(r)) #%>%
#mutate(direction = case_when( #there aren't negative correlations with spearman < -.70
#  r > 0 ~ "positive",
#  r < 0 ~ "negative"))

cor.graph
cor.graph %>%
  ggraph(., layout = "kk") +
  #geom_edge_density(aes(fill = r)) + 
  geom_edge_link(aes(width = r,edge_colour = r)) +
  #scale_edge_color_manual(values = c("blue", "red")) +
  scale_edge_color_gradient(low = "white", high = Dark24[1]) +
  scale_edge_width(range = c(0.2, 2)) +
  #geom_node_point(aes(color = family, shape = impgenes), size = 3) +
  geom_node_point(aes(color=Feature2), size = 2) +
  #scale_shape_manual(values = c(15,16)) +
  scale_color_manual(values = c("black",Dark24[15])) +
  geom_node_text(aes(label = Feature, color=Feature2), size = 4, fontface="bold",
                 max.overlaps = 20, repel = TRUE) +
  theme_graph()

# Individual correlations:
corr_transp <- rownames_to_column(as.data.frame(t(Saureus_corr_mat[c(res.cor$x[1:5],Saureus_variable),])), "sample_RNAseq")

corr_transp %>%
  pivot_longer(!sample_RNAseq & !`S.aureus CPM abundance`, names_to = "feature", values_to = "value") %>%
  ggplot(., aes(x=value, y=`S.aureus CPM abundance`)) +
  geom_smooth(method=lm, color=Dark24[1]) +
  geom_point(size=2, color= "#424242") +
  theme_classic() + #scale_color_manual(values = PalC) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11),
        legend.position="right", strip.background = element_blank(), strip.text = element_text(size = 10, face="bold"),
        axis.text = element_blank(), axis.ticks = element_blank()) +
  stat_cor(method = "pearson",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE, color=Dark24[1]) + 
  xlab("Gene expression (log2CPM)") + ylab("S. aureus CPM abundance") +
  facet_wrap(. ~ feature, nrow = 1, scales = "free")

# Biomart gene description ----
library(biomaRt)
anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
Tx.other <- getBM(attributes=c('external_gene_name',
                               'entrezgene_description'),
                  mart = anno)

gene_Fullname <- res.cor %>%
  mutate(external_gene_name = x) %>%
  left_join(Tx.other) %>%
  dplyr::select(external_gene_name, entrezgene_description ,r)

write_tsv(gene_Fullname, "gene_Fullname.txt")

# xCell and S. aureus detect. non-detect ----
library(immunedeconv)
library(broom)
library(corrplot)
library(patchwork)
library(ggpubr)
library(reshape2)
library(gplots)
library(Hmisc)

library(tidyverse)
deconvolution_methods
matrix <- 2^log2.cpm.filtered.norm[,targets_CL_Saureus$sample_RNAseq]

# `xcell`
res_cell <- deconvolute(matrix, "xcell")
#res_cell <- deconvolute(matrix, "mcp_counter")
res_cell <- as.matrix(column_to_rownames(res_cell, "cell_type"))
res_cell_mat <- rbind(res_cell, targets_CL_Saureus$cpm_saureus)
rownames(res_cell_mat)[40] <- "S. aureus abundances"

# the correlation:
cellgene_stats <- rcorr(as.matrix(t(res_cell_mat)), type = "pearson")
#cellgene_cormat <- format(as.matrix(cellgene_stats$r), digits = 3)
#cellgene_pmat <- format(as.matrix(cellgene_stats$P), digits = 3)

cellgene_stats_tidy <- tidy(cellgene_stats)
cellgene_stats_tidy <- cellgene_stats_tidy%>%
  filter(estimate > 0.5 & p.value < 0.05)
#View(cellgene_stats_tidy)

# So, correlations were not informative.
colorheat <- colorRampPalette(colors=c(Dark24[2],"white",Dark24[1]))(20)
corrplot(cellgene_stats$r, method = "color", #order="hclust", 
         p.mat = cellgene_stats$P, #sig.level = 0.05, insig = "blank",
         col = colorheat,tl.col = "black",tl.cex = 0.6, pch.col = "white")

# Categorically then:
matrix2 <- 2^log2.cpm.filtered.norm[,CL_samples]

# `xcell`
res_cell2 <- deconvolute(matrix2, "xcell")
#res_cell2 <- deconvolute(matrix2, "mcp_counter")

#all:
View(rownames_to_column(as.data.frame(t(column_to_rownames(res_cell2, "cell_type"))), "sample_RNAseq") )
rownames_to_column(as.data.frame(t(column_to_rownames(res_cell2, "cell_type"))), "sample_RNAseq") %>%
  pivot_longer(!sample_RNAseq, names_to = "Cell_type", values_to = "Cell_score") %>%
  left_join(targets_CL, by="sample_RNAseq") %>%
  mutate(cpm_saureus_cat = factor(cpm_saureus_cat, levels = c("Saureus_NOTdetec","Saureus_detec"))) %>%
  ggplot(.,
         aes(x=cpm_saureus_cat, y=Cell_score, color=cpm_saureus_cat)) +
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
  facet_wrap(. ~ Cell_type, scales = "free")

# Network in the upgenes in Saureus_detect NOT UPDATEDDDDDD ----
library(corrplot)
library(corrr)
library(tidygraph)
library(ggraph)

clusters_genes <- hclust(as.dist(1-cor(t(myDGEList.filtered.norm$counts[c(sig_genes_up8),targets_CL$sample_RNAseq]), method="pearson")), method="ward.D2") 
plot(clusters_genes)

as_tbl_graph(as.dendrogram(clusters_genes)) %>%
  ggraph(., layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal()+
  geom_node_point() +
  geom_node_text(aes(label=label, filter=leaf), size = 4, fontface="bold",
                 color="black", max.overlaps = 20, repel = TRUE)



res.cor <- t(myDGEList.filtered.norm$counts[c(sig_genes_up8),targets_CL$sample_RNAseq]) %>% correlate(method = "pearson") %>%
  shave(upper = TRUE) %>%
  stretch(na.rm = TRUE)# %>%
  #filter(r >= 0.75 | r <= -0.75)

cor.graph <- as_tbl_graph(res.cor, directed = FALSE)

#'%ni%' <- Negate('%in%')
cor.graph <- cor.graph %>% 
  activate(nodes) %>% 
  rename(Gene = name) %>%
  #left_join(description, by = "Gene") %>%
  #mutate(impgenes = case_when(
  #  Gene %in% impgenes ~ "impgenes",
  #  Gene %ni% impgenes ~ "notimpgenes")) %>%
  activate(edges) %>%
  arrange(desc(r)) #%>%
#mutate(direction = case_when( #there aren't negative correlations with spearman < -.70
#  r > 0 ~ "positive",
#  r < 0 ~ "negative"))

cor.graph
cor.graph %>%
  ggraph(., layout = "treemap") +
  #geom_node_circle(aes(fill = "r")) +
  #geom_edge_density(aes(fill = r)) + 
  #geom_edge_link(aes(width = r,edge_colour = r)) +
  #scale_edge_color_manual(values = c("blue", "red")) +
  #scale_edge_color_gradient(low = "gray", high = Dark24) +
  #scale_edge_width(range = c(0.2, 1.3)) +
  #geom_node_point(aes(color = family, shape = impgenes), size = 3) +
  #geom_node_point(aes(color = Family), size = 3) +
  #scale_shape_manual(values = c(15,16)) +
  #scale_color_manual(values = Salvador) +
  geom_node_text(aes(label = Gene), size = 4, fontface="bold",
                 color="black", max.overlaps = 20, repel = TRUE) +
  theme_graph()

ggsave("outputs_inputs/Network Chemokines and Receptors.tiff",
       width = 10, height = 6)


# Venn in the upgenes in Saureus_detect with GO BP (NOT UPDATED since working with VITALS) ----
# The list of 78 genes was submitted to DAVID GO analysis. Reactome and BP were accessed and 46 unique genes were associated
# With at least one geneset from 11 GO annotation clusters (genesets P value < 0.05)
GO_StaphDetec_venn <- venn(list(Cluster1 = c(unique(c("IL11","ANXA1","CXCL8","CSF2","MMP1","CCL3L3","IL24","OSM","LILRA3","LILRA5","IL1A","IL1B","CCL4","CCL3","IL21","IL11","ANXA1","CXCL8","CSF2","MMP1","CCL3L3","IL24","OSM","IL1A","IFNG","IL1B","CCL4","CCL3","S100A12","IL1A","CXCL8","CSF2","IL1B","CCL3L3","CCL4","CCL3","IL21","IL11","TNFRSF6B","ANXA1","CXCL8","CSF2","MMP1","CCL3L3","IL24","OSM","ISG15","IL1A","IFNG","IL1B","CCL4","CCL3","S100A12","IL1A","CXCL8","ANXA1","MMP1","IL1B","OSM","IL1A","IL1B","OSM","IL1A","CXCL8","IL1B","CCL3L3","IL24","OSM"),
                                                    fromLast = FALSE)),
                                Cluster2 = c(unique(c("DEFB103B","KLRC2","TNFAIP6","ISG15","FPR2","TREM1","LILRA3","ADGRG3","CEACAM3","KLRK1","CLEC4D","SLPI","GNLY","IL1B","CLEC6A","S100A12","CD300E","FOLR3","SLC27A2","ADGRG3","CEACAM3","CLEC4D","TNFAIP6","SLPI","S100A12","FOLR3","FPR2","LILRA3","SLC27A2","ADGRG3","CEACAM3","CLEC4D","TNFAIP6","SLPI","S100A12","FOLR3","FPR2","LILRA3","SLC27A2"),
                                                    fromLast = FALSE)),
                                Cluster3 = c(unique(c("CCL8","CXCL8","CCL7","IL1B","CCL3L3","CCL4","CCL3","S100A12","TREM1","CXCL5","ANXA1","CXCL8","TNFAIP6","CCL3L3","ACOD1","FPR2","CXCL5","IL1A","CCL8","CCL7","IL1B","CCL4","CCL3","S100A12","CCL8","ANXA1","CCL7","CCL3L3","CCL4","CCL3","S100A12","CCL8","CXCL8","CCL7","CCL3L3","CCL4","CCL3","CXCL5","CCL8","CXCL8","CCL7","CCL3L3","CCL4","CCL3","ACOD1","CCL8","CCL7","CCL3L3","CCL4","CCL3","CXCL8","CCL7","CCL3L3","CCL4","CCL3","CXCL5","CCL8","CXCL8","CCL7","CCL3L3","CCL4","CCL3","ACOD1","CCL8","CCL7","CCL3L3","CCL4","CCL3","ACOD1","CCL8","CCL7","CCL4","CCL3","CXCL8","ANXA1","CCL7","CCL3L3","CCL4","CCL3","FPR2","CXCL5","CCL8","CXCL8","CCL7","CCL3","FPR2","CXCL5","IL1A","CCL8","CCL7","CCL3L3","CCL4","CCL3","FPR2","CCL8","CCL7","TNFAIP6","IL1B","CCL4","CCL3","CXCL5","CCL7","CCL4","CCL3","CXCL8","ANXA1","CCL7","CCL3L3","CCL4","CCL3","FFAR2","FPR2","CXCL5","CXCL8","ANXA1","CCL7","CCL3L3","CCL4","CCL3","FFAR2","FPR2","CXCL5","CXCL8","ANXA1","CCL7","CCL4L2","CCL3L3","CCL4","CCL3","FFAR2","FPR2","CXCL5","ADGRG3","CCL8","CXCL8","ANXA1","CCL7","CCL3L3","CCL4","CCL3","FFAR2","FPR2","CXCL5","CXCL8","ANXA1","CCL4L2","CCL4","FPR2","CXCL5","CCL8","CCL7","CCL3L3","CCL4","CCL3","CXCL8","CCL3","FPR2"),
                                                    fromLast = FALSE)),
                                Cluster4 = c(unique(c("IL1A","CXCL8","CSF2","IL1B","CCL3L3","CCL4","CCL3","IL1A","KLRK1","CXCL8","CSF2","IL1B","IL24","ACOD1","CXCL5","IL1A","CSF2","IL1B","IL1A","CXCL8","IL1B","CCL3L3","IL24","OSM","IL1A","CXCL8","CSF2","IFNG","IL1B","CCL3"),
                                                    fromLast = FALSE)),
                                Cluster5 = c(unique(c("IL21","IFNG","IL1B","OSM","CCL3","S100A12","LILRA5","IL1A","IFNG","IL1B","LILRA5","IL1A","IFNG","CCL3","LILRA5","IFNG","CCL3","LILRA5","IL1A","CXCL8","CSF2","IFNG","IL1B","CCL3"),
                                                    fromLast = FALSE)),
                                Cluster6 = c(unique(c("MT1L","MT1G","MT1H","MT1L","MT1G","MT1H","MT1L","MT1G","MT1H","MT1L","MT1G","MT1H","MT1L","MT1G","MT1H","MT1L","MT1G","MT1H"),
                                                    fromLast = FALSE)),
                                Cluster7 = c(unique(c("KLRK1","CD300E","LILRA3","TREM1","KIR2DL3","LILRA5","KLRK1","KLRC2","CD300E","TREM1","KLRK1","CD300E","IGHV3-48","TREM1","KIR2DL3"),
                                                    fromLast = FALSE)),
                                Cluster8 = c(unique(c("IFNG","IL1B","FPR2","IFNG","IL1B","FPR2"),
                                                    fromLast = FALSE)),
                                Cluster9 = c(unique(c("IL1A","CLEC4D","IL1B","CLEC6A","S100A12"),
                                                    fromLast = FALSE)),
                                Cluster10 = c(unique(c("IL1A","IL1B","GZMB","IL1A","IL1B","GZMB"),
                                                    fromLast = FALSE)),
                                Cluster11 = c(unique(c("CEACAM3","MMP1","IGHV3-48","TREM1"),
                                                    fromLast = FALSE))))

BP_genes <- unique(c(unique(c("IL11","ANXA1","CXCL8","CSF2","MMP1","CCL3L3","IL24","OSM","LILRA3","LILRA5","IL1A","IL1B","CCL4","CCL3","IL21","IL11","ANXA1","CXCL8","CSF2","MMP1","CCL3L3","IL24","OSM","IL1A","IFNG","IL1B","CCL4","CCL3","S100A12","IL1A","CXCL8","CSF2","IL1B","CCL3L3","CCL4","CCL3","IL21","IL11","TNFRSF6B","ANXA1","CXCL8","CSF2","MMP1","CCL3L3","IL24","OSM","ISG15","IL1A","IFNG","IL1B","CCL4","CCL3","S100A12","IL1A","CXCL8","ANXA1","MMP1","IL1B","OSM","IL1A","IL1B","OSM","IL1A","CXCL8","IL1B","CCL3L3","IL24","OSM"),
                fromLast = FALSE),
         unique(c("DEFB103B","KLRC2","TNFAIP6","ISG15","FPR2","TREM1","LILRA3","ADGRG3","CEACAM3","KLRK1","CLEC4D","SLPI","GNLY","IL1B","CLEC6A","S100A12","CD300E","FOLR3","SLC27A2","ADGRG3","CEACAM3","CLEC4D","TNFAIP6","SLPI","S100A12","FOLR3","FPR2","LILRA3","SLC27A2","ADGRG3","CEACAM3","CLEC4D","TNFAIP6","SLPI","S100A12","FOLR3","FPR2","LILRA3","SLC27A2"),
                fromLast = FALSE),
         unique(c("CCL8","CXCL8","CCL7","IL1B","CCL3L3","CCL4","CCL3","S100A12","TREM1","CXCL5","ANXA1","CXCL8","TNFAIP6","CCL3L3","ACOD1","FPR2","CXCL5","IL1A","CCL8","CCL7","IL1B","CCL4","CCL3","S100A12","CCL8","ANXA1","CCL7","CCL3L3","CCL4","CCL3","S100A12","CCL8","CXCL8","CCL7","CCL3L3","CCL4","CCL3","CXCL5","CCL8","CXCL8","CCL7","CCL3L3","CCL4","CCL3","ACOD1","CCL8","CCL7","CCL3L3","CCL4","CCL3","CXCL8","CCL7","CCL3L3","CCL4","CCL3","CXCL5","CCL8","CXCL8","CCL7","CCL3L3","CCL4","CCL3","ACOD1","CCL8","CCL7","CCL3L3","CCL4","CCL3","ACOD1","CCL8","CCL7","CCL4","CCL3","CXCL8","ANXA1","CCL7","CCL3L3","CCL4","CCL3","FPR2","CXCL5","CCL8","CXCL8","CCL7","CCL3","FPR2","CXCL5","IL1A","CCL8","CCL7","CCL3L3","CCL4","CCL3","FPR2","CCL8","CCL7","TNFAIP6","IL1B","CCL4","CCL3","CXCL5","CCL7","CCL4","CCL3","CXCL8","ANXA1","CCL7","CCL3L3","CCL4","CCL3","FFAR2","FPR2","CXCL5","CXCL8","ANXA1","CCL7","CCL3L3","CCL4","CCL3","FFAR2","FPR2","CXCL5","CXCL8","ANXA1","CCL7","CCL4L2","CCL3L3","CCL4","CCL3","FFAR2","FPR2","CXCL5","ADGRG3","CCL8","CXCL8","ANXA1","CCL7","CCL3L3","CCL4","CCL3","FFAR2","FPR2","CXCL5","CXCL8","ANXA1","CCL4L2","CCL4","FPR2","CXCL5","CCL8","CCL7","CCL3L3","CCL4","CCL3","CXCL8","CCL3","FPR2"),
                fromLast = FALSE),
         unique(c("IL1A","CXCL8","CSF2","IL1B","CCL3L3","CCL4","CCL3","IL1A","KLRK1","CXCL8","CSF2","IL1B","IL24","ACOD1","CXCL5","IL1A","CSF2","IL1B","IL1A","CXCL8","IL1B","CCL3L3","IL24","OSM","IL1A","CXCL8","CSF2","IFNG","IL1B","CCL3"),
                fromLast = FALSE),
         unique(c("IL21","IFNG","IL1B","OSM","CCL3","S100A12","LILRA5","IL1A","IFNG","IL1B","LILRA5","IL1A","IFNG","CCL3","LILRA5","IFNG","CCL3","LILRA5","IL1A","CXCL8","CSF2","IFNG","IL1B","CCL3"),
                fromLast = FALSE),
         unique(c("MT1L","MT1G","MT1H","MT1L","MT1G","MT1H","MT1L","MT1G","MT1H","MT1L","MT1G","MT1H","MT1L","MT1G","MT1H","MT1L","MT1G","MT1H"),
                fromLast = FALSE),
         unique(c("KLRK1","CD300E","LILRA3","TREM1","KIR2DL3","LILRA5","KLRK1","KLRC2","CD300E","TREM1","KLRK1","CD300E","IGHV3-48","TREM1","KIR2DL3"),
                fromLast = FALSE),
         unique(c("IFNG","IL1B","FPR2","IFNG","IL1B","FPR2"),
                fromLast = FALSE),
         unique(c("IL1A","CLEC4D","IL1B","CLEC6A","S100A12"),
                fromLast = FALSE),
         unique(c("IL1A","IL1B","GZMB","IL1A","IL1B","GZMB"),
                fromLast = FALSE),
         unique(c("CEACAM3","MMP1","IGHV3-48","TREM1"),
                fromLast = FALSE)), fromLast = FALSE)

# Wordcloud ----
library("tm")
library("SnowballC")
library(ggwordcloud)

# Load the text:
geneset <- readLines(file.choose())

# Load the data as a corpus
geneset <- Corpus(VectorSource(geneset))

# Text transformation:
geneset <- tm_map(geneset, removePunctuation) # Remove punctuation

# Inspect:
inspect(geneset)

# Build a term-document matrix
dtm <- TermDocumentMatrix(geneset)
m <- as.matrix(dtm)
rownames(m) <- str_to_upper(rownames(m))
v <- sort(rowSums(m),decreasing=TRUE)
geneset <- data.frame(word = names(v),freq=v)
geneset$freq <- as.numeric(geneset$freq)
geneset <- arrange(geneset, desc(freq))

geneset$word <- as.factor(geneset$word)
#View(geneset8)

# Highlight biomarkers for treatment outcome ----
#load("../RStudio_outputs/data/myTopHits_treat_var_imp")
geneset <- geneset %>%
  mutate(Failure_biomarker = case_when(
    word %in% myTopHits_treat_var_imp ~ "biomarker",
    TRUE ~ "not_biomarker"
  ))


#barplot(geneset$freq, las = 2, names.arg = geneset$word,
#        col ="pink", main = "blabla",
#        ylab = "Gene frequency in the geneset")


ggplot(geneset, aes(label = word, size = freq, color=word)) +
  scale_radius(range = c(1, 10), limits = c(0, NA)) +
  geom_text_wordcloud(area_corr = TRUE) +
  scale_size_area(max_size = 24) +
  scale_color_manual(values = Dark24) +
  theme_minimal()

# Top 3 clusters:
geneset1_plot +
geneset2_plot +
geneset3_plot




