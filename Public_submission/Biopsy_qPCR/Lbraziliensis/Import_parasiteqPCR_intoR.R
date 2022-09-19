# Here I import the results from Ba's parasite qPCR.
# See ~Dropbox/CamilaAnalysis/Lesion3rd_dataset_Microbiome/3rd_LesionRNAseq_Dataset/qPCRs/Notes_qPCR to see how I processed the calculations.

# Load R objects ----
load("../../RStudio_outputs/data/Dark24")
load("../../RStudio_outputs/data/targets")
load("../../RStudio_outputs/data/myDGEList.filtered.norm")
load("../../RStudio_outputs/data/CL_samples")

# Libraries ----
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(patchwork)

# Processing and normalizing the interpolated values ----
# Prism's export:
qPCR_1 <- read_delim("Nonlin fit of Lbraziliensis_removeOutlier.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
qPCR_1 <- qPCR_1 %>%
  select(1,3) %>%
  rename(interpolated_1 = 2)

# Tissue weight:
tw <- read_delim("~/Dropbox/CamilaAnalysis/Lesion3rd_dataset_Microbiome/3rd_LesionRNAseq_Dataset/RStudio_inputs/TissueWeight_FirstHalfBiopsy.txt", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)

tw <- tw %>%
  mutate(TissueWeight_mg = TissueWeight_g*1000)

# Normalize Bacteria data to tissue weight, and merge everything:
qPCR_1 <- left_join(targets, qPCR_1, by="sample_RNAseq")
qPCR_1 <- left_join(qPCR_1, tw, by="sample_RNAseq")
colnames(qPCR_1)
qPCR_1 %>%
  select(sample_RNAseq, interpolated_1, TissueWeight_mg)

qPCR_1 <- qPCR_1 %>%
  mutate(Lbraziliensis_mg = interpolated_1/TissueWeight_mg)

# Visualizing ----
ggplot(qPCR_1,
       aes(x=factor(group_RNAseq, levels=c("HS","CL")),
           y=Lbraziliensis_mg, color=group_RNAseq)) +
  geom_jitter(position=position_jitter(0.1), size=1) +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="black") +
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1), color="black") +
  stat_compare_means(method = "t.test", #label.y = 30000,
                    aes(label = ..p.signif..),
                     label.x = 1.5, size = 5, color="red") +
  xlab("") + ylab(paste("number of L. braziliensis/mg biopsy)")) +

# cpm vs. pg:
ggplot(qPCR_1 %>%
         filter(group_RNAseq == "CL"),
       aes(x=cpm_leishmania/10000, y=Lbraziliensis_mg)) +
  geom_smooth(method=lm, color=Dark24[1]) +
  geom_point(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2.5, fontface="bold",colour="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) +
  xlab("L. braziliensis (CPM/10000)") + ylab(paste("number of L. braziliensis/mg biopsy)")) +

# Failure vs. Cure:
ggplot(qPCR_1 %>%
         filter(group_RNAseq == "CL") %>%
         filter(treatment_other_drug == "sbv"),
       aes(x=treatment_outcome, y=Lbraziliensis_mg,
           color=treatment_outcome)) +
  geom_jitter(position=position_jitter(0.1), size=2) +
  #annotate("text", x=1.5,y=3500, label=paste("2000 parasite/mg"),
  #         size=2, fontface="bold",color="black") +
  #geom_hline(yintercept = 2000) +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, hjust = 1, angle = 45), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  stat_compare_means(method = "t.test", #label.y = 30000,
                     aes(label = ..p.signif..),
                     label.x = 1.5, size = 5, color="red") +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2.5, fontface="bold",colour="black") +
  xlab("") + ylab(paste("number of L. braziliensis/mg biopsy)")) +
  NULL

# Percentiles:
qPCR_1 %>%
  filter(group_RNAseq != "HS") %>%
  ggplot(.,aes(x=group_RNAseq, y=Lbraziliensis_mg)) +
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
  geom_text_repel(aes(label = sample_RNAseq), size = 3, fontface="bold", nudge_x = 0.5, max.overlaps = 100) +
  xlab("") + ylab(paste("Failure gene signature ES")) +
  #ylim(0,.40) +
  NULL

# DE qPCR parasite biopsy, absolute (continuous) ----
design9 <- model.matrix(~qPCR_1$Lbraziliensis_mg[1:51]) # pay attention here, it's not automatic! Only selecting values from patients.
colnames(design9)[2] <- "qPCR_2_absolute"
head(design9)
dim(design9)
v.DEGList.filtered.norm9 <- voom(myDGEList.filtered.norm[,CL_samples[-47]], design9)
fit9 <- lmFit(v.DEGList.filtered.norm9, design9)
fits9 <- contrasts.fit(fit9, coefficients = 2)
ebFit9 <- eBayes(fits9)
myTopHits9 <- topTable(ebFit9, adjust ="BH", number=40000, sort.by="logFC")
myTopHits9 <- as_tibble(myTopHits9, rownames = "geneSymbol")

sig_genes_up9 <- myTopHits9 %>%
  filter(P.Value < 0.05, logFC > 0.00002)
sig_genes_down9 <- myTopHits9 %>%
  filter(P.Value < 0.05, logFC < -0.00002)
sig_genes_up9 <- sig_genes_up9$geneSymbol
sig_genes_down9 <- sig_genes_down9$geneSymbol
#write_tsv(as_tibble(sig_genes_up9), "../../../../../../Desktop/sig_genes_up_pval05Continuous.txt")

myTopHits9$col0 <- myTopHits9$geneSymbol
col0_v9 <- myTopHits9$col0 %in% sig_genes_up9
col0_v99 <- myTopHits9$col0 %in% sig_genes_down9
myTopHits9$col0[!col0_v9 & !col0_v99] <- NA

#volcano:
ggplot(myTopHits9, aes(y=-log10(P.Value), x=logFC)) +
  geom_point(size=1, color="gray") +
  theme_classic() +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0.00002) + geom_vline(xintercept = -0.00002) +
  #geom_text_repel(aes(label = col0), size = 2.2, fontface="bold",color="black") +
  xlab("logFC")

# Export parasite numbers data ----
write_tsv(qPCR_1 %>%
          select(sample_RNAseq, cpm_leishmania, Lbraziliensis_mg),
          "../../RStudio_outputs/data/parasiteload_3rdDataset.txt")


parasiteload_3rdDataset <- qPCR_1 %>%
  select(sample_RNAseq, cpm_leishmania, Lbraziliensis_mg)
save(parasiteload_3rdDataset, file = "parasiteload_3rdDataset")
