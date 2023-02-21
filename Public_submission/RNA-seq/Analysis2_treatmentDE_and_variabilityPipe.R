# In the previous manuscript, we found biomarkers for treatment outcome in two independent cohorts
# Now we are confirming in this current dataset

# Making vectors and groups Failure vs. Cure ----
# Selecting patients with treatment outcome:
cured <- targets %>% filter(treatment_outcome == "cure") %>%
  filter(treatment_other_drug == "sbv")
cured <- cured$sample_RNAseq

failed <- targets %>% filter(treatment_outcome == "failure") %>%
  filter(treatment_other_drug == "sbv")
failed <- failed$sample_RNAseq

targets_treat <- targets %>% filter(group_RNAseq == "CL") %>% # targets_treat excludes the patients treated with miltefosine
  filter(treatment_other_drug == "sbv")
groups1_treat <- targets_treat$treatment_outcome
groups1_treat <- factor(groups1_treat)

# Processing ----
# Filtering
min_sample2 = as.vector(table(groups1_treat))
num_classes2 = length(min_sample2)
min_size2 = min_sample2[order(min_sample2,decreasing=FALSE)[1]]

cpm_treat <- cpm(myDGEList[,targets_treat$sample_RNAseq])
keepers_treat <- rowSums(cpm_treat>1)>=min_size2
myDGEList.filtered_treat <- myDGEList[keepers_treat,targets_treat$sample_RNAseq]

dim(myDGEList[,targets_treat$sample_RNAseq])
dim(myDGEList.filtered_treat)

# Normalization
myDGEList.filtered.norm_treat <- calcNormFactors(myDGEList.filtered_treat, method = "TMM")
log2.cpm.filtered.norm_treat <- cpm(myDGEList.filtered.norm_treat, log=TRUE)

# PCA ----
pca.res_treat <- prcomp(t(log2.cpm.filtered.norm_treat), scale.=F, retx=T)
pc.var_treat<-pca.res_treat$sdev^2
pc.per_treat<-round(pc.var_treat/sum(pc.var_treat)*100, 1)

pca.res.df_treat <- as_tibble(pca.res_treat$x)
ggplot(pca.res.df_treat, aes(x=PC1, y=PC2,
                             color=factor(targets_treat$treatment_outcome, levels = c("failure","cure")))) +
  theme_classic() + scale_color_manual(values = Dark24) +
  geom_point(size=2.5) +
#  geom_text_repel(aes(label = targets_treat$sample_RNAseq), size = 2.5, fontface="bold",colour="black") +
  theme(legend.position="top",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  stat_ellipse(level = 0.95) +
  xlab(paste("PC1 -",pc.per_treat[1],"%")) +
  ylab(paste("PC2 -",pc.per_treat[2],"%"))


dist_treat <- vegdist(t(2^log2.cpm.filtered.norm_treat),
                             method = "bray")
vegan_treat <- adonis2(dist_treat~targets_treat$treatment_outcome,
                              permutations = 999, method="bray")

# DE analysis Failure vs Cure ----
# design matrix:
design_treat <- model.matrix(~0 + groups1_treat)
colnames(design_treat) <- levels(groups1_treat)

# Use VOOM function from Limma package to model the mean-variance relationship
v.DEGList.filtered.norm_treat <- voom(myDGEList.filtered.norm_treat, design_treat)

# fit a linear model to your data
fit_treat <- lmFit(v.DEGList.filtered.norm_treat, design_treat)

# Contrast matrix
contrast.matrix_treat <- makeContrasts(treatment = failure - cure,
                                 levels=design_treat)

# extract the linear model fit
fits_treat <- contrasts.fit(fit_treat, contrast.matrix_treat)

#get bayesian stats for your linear model fit
ebFit_treat <- eBayes(fits_treat)

# MyTopHits and Volcano ----
myTopHits_treat <- topTable(ebFit_treat, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits_treat <- as_tibble(myTopHits_treat, rownames = "geneSymbol")
save(myTopHits_treat, file = "RStudio_outputs/data/myTopHits_treat")

# Variability workflow ----
# Filtering all the genes that were up-upregulated CL lesions vs. HS:
# getting the expression levels of these genes for each CL sample:
upCL_exp <- 2^notfilt_norm_log2cpm[sig_genes_up,targets_treat$sample_RNAseq]

# calculate sd, mean, cv for regulated genes in CL lesions
mean_genes <- rowMeans(upCL_exp)
sd_genes <- apply(upCL_exp, 1, sd)
cv_genes <- coefficient.variation(sd_genes, mean_genes)

# making a dataframe with myTopHits and cv:
upCL_myTopHits <- myTopHits %>%
  filter(logFC > 1) %>%
  filter(adj.P.Val < 0.05)
allgenes_stats <- as_tibble(cbind(upCL_myTopHits,
                                  cv_genes))


View(as_tibble(cbind(upCL_myTopHits,
                cv_genes,sd_genes)) %>%
       arrange(cv_genes,sd_genes))

# visualize names:
viz_var <- allgenes_stats %>%
  filter(cv_genes > 0.6506)
viz_var_names <- viz_var$geneSymbol
length(viz_var_names)
#write_tsv(as_tibble(viz_var_names), "viz_var_names.txt")

# Significant genes:
allgenes_stats$colV <- allgenes_stats$geneSymbol
viz_var_f <- allgenes_stats$colV %in% viz_var_names
allgenes_stats$colV[!viz_var_f] <- NA

allgenes_stats$col_show <- allgenes_stats$geneSymbol
col_show_v <- allgenes_stats$col_show %in% c("GZMB","GNLY","PRF1","IL1B")
allgenes_stats$col_show[!col_show_v] <- NA

# getting FC between CL vs HS, using limma DE analysis:
ggExtra::ggMarginal(
  ggplot(allgenes_stats, aes(y=cv_genes, x=logFC)) +
    geom_point(size=1, color="gray") +
    theme_classic() + scale_color_few() +
    #geom_text_repel(aes(label = col_show), size = 3.5, fontface="bold",colour="black",nudge_y = 0.5) +
    geom_hline(yintercept = 0.6506, color=Dark24[1]) +
#    geom_hline(yintercept = 0.5, color=Dark24[2]) +
    theme(axis.text = element_text(size=13), axis.title = element_text(size=13)) +
    annotate("text", x=9, y=0.6506+0.2, label=paste("Cv=0.65"), fontface=4, size=4, color=Dark24[1]) +
#    annotate("text", x=9, y=0.5+0.08, label=paste("Cv=0.5"), fontface=4, size=4, color=Dark24[2]) +
    xlab("logFC CL relative to HS") +
    ylab("Coefficient of Variation (CV)"),
  type = 'density', margins = 'both', size = 10, col = 'grey', fill = 'grey')

myTopHits_treat_var <- myTopHits_treat[myTopHits_treat$geneSymbol %in% viz_var_names,]
sigTreat <- myTopHits_treat_var %>% dplyr::arrange(P.Value) %>%
  filter(P.Value < 0.05)
sigTreat$geneSymbol

# Volcano Failure vs. Cure ----
myTopHits_treat_var$colV <- myTopHits_treat_var$geneSymbol
viz_var_f_treat <- myTopHits_treat_var$colV %in% sigTreat$geneSymbol
myTopHits_treat_var$colV[!viz_var_f_treat] <- NA

myTopHits_treat_var$col_show2 <- myTopHits_treat_var$geneSymbol
col_show2_v <- myTopHits_treat_var$col_show2 %in% c("GZMB","GNLY","PRF1","IL1B")
myTopHits_treat_var$col_show2[!col_show2_v] <- NA

ggplot(myTopHits_treat_var, aes(y=-log10(P.Value), x=logFC)) +
  annotate("text", x=0, y=5, label=paste("ngenes P<0.05 =", length(sigTreat$geneSymbol)), fontface=4, size=4, color=Dark24[1]) +
  annotate("text", x=0, y=-log10(0.05)+0.1, label=paste("P<0.05"), fontface=4, size=4, color=Dark24[10]) +
  geom_hline(yintercept = -log10(0.05), color=Dark24[10]) + geom_vline(xintercept = 1) +
  geom_point(size=1) +
  theme_classic() +
#  scale_color_manual(values = c("dark gray", Dark24[2], Dark24[1])) +
  geom_text_repel(aes(label = colV), size = 3.5, fontface="bold",colour="black") +
  geom_text_repel(aes(label = col_show2), size = 3.5, fontface="bold",colour=Dark24[4],nudge_x = 0.5) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +

  xlab("logFC Failure vs. Cure")

# Venn diagram with biomarkers for treatment outcome in 3 datasets ----
load("RStudio_inputs/biomarkers_1st") # biomarkers from Mosser's dataset
load("RStudio_inputs/biomarkers_2nd") # biomarkers from Amorim 2019 dataset
biomarkers_3rd <- sigTreat$geneSymbol # biomarkers for this dataset
save(biomarkers_3rd, file = "RStudio_outputs/data/biomarkers_3rd")

length(biomarkers_1st)
length(biomarkers_2nd)
length(biomarkers_3rd)


venn_biomarkers <- venn(list(Mosser2016 = biomarkers_1st,
                             Amorim2019 = biomarkers_2nd,
                             Current = biomarkers_3rd))
attr(venn_biomarkers,"intersections")$`Amorim2019:Current`
attr(venn_biomarkers,"intersections")$`Mosser2016:Current`
attr(venn_biomarkers,"intersections")$`Mosser2016:Amorim2019`
attr(venn_biomarkers,"intersections")$`Mosser2016:Amorim2019:Current`

# modifying the mytophits_treat a little bit more to include shared biomarkers for treatment outcome in the plot ----
myTopHits_treat_var$col_col <- myTopHits_treat_var$geneSymbol
uniquebio <- myTopHits_treat_var$col_col %in% c(attr(venn_biomarkers,"intersections")$Current,
                                                attr(venn_biomarkers,"intersections")$Mosser2016,
                                                attr(venn_biomarkers,"intersections")$Amorim2019)
shared2 <- myTopHits_treat_var$col_col %in% c(attr(venn_biomarkers,"intersections")$`Amorim2019:Current`,
                                              attr(venn_biomarkers,"intersections")$`Mosser2016:Current`,
                                              attr(venn_biomarkers,"intersections")$`Mosser2016:Amorim2019`)
shared3 <- myTopHits_treat_var$col_col %in% c(attr(venn_biomarkers,"intersections")$`Mosser2016:Amorim2019:Current`)

save(uniquebio, file = "../../CD8_dependent_pathology/Robjects/uniquebio")
save(shared2, file = "../../CD8_dependent_pathology/Robjects/shared2")
save(shared3, file = "../../CD8_dependent_pathology/Robjects/shared3")

myTopHits_treat_var$col_col[!uniquebio & !shared2 & !shared3] <- NA

myTopHits_treat_var$col_col2 <- myTopHits_treat_var$geneSymbol
myTopHits_treat_var$col_col2[uniquebio] <- "unique"
myTopHits_treat_var$col_col2[shared2] <- "shared in 2"
myTopHits_treat_var$col_col2[shared3] <- "shared in 3"
myTopHits_treat_var$col_col2[!uniquebio & !shared2 & !shared3] <- "not biomarker"

ggplot(myTopHits_treat_var, aes(y=-log10(P.Value), x=logFC,
                                color=factor(col_col2, levels = c("unique","shared in 2",
                                                                  "shared in 3", "not biomarker")))) +
  annotate("text", x=-1.5, y=-log10(0.05)+0.1, label=paste("P<0.05"), fontface=4, size=4) +
  geom_hline(yintercept = -log10(0.05)) + #geom_vline(xintercept = 1) +
  geom_point(size=1) +
  theme_classic() +
  scale_color_manual(values = c(Dark24[8],Dark24[2],Dark24[1],"gray")) +
  geom_text_repel(aes(label = col_col, color=col_col2), size = 3.5, fontface=2,
                  max.overlaps = 500) +
  theme(legend.position="right",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  xlab("logFC Failure vs. Cure")

# Export gene names biomarkers for treatmento outcome in at least 2 datasets:
myTopHits_treat_var_imp <- unique(c(attr(venn_biomarkers,"intersections")$`Amorim2019:Current`,
                                    attr(venn_biomarkers,"intersections")$`Mosser2016:Current`,
                                    attr(venn_biomarkers,"intersections")$`Mosser2016:Amorim2019`,
                                    attr(venn_biomarkers,"intersections")$`Mosser2016:Amorim2019:Current`))

save(myTopHits_treat_var_imp, file = "RStudio_outputs/data/myTopHits_treat_var_imp")
write_tsv(as.data.frame(myTopHits_treat_var_imp), "RStudio_outputs/data/myTopHits_treat_var_imp.txt")

# Take topHits for treatment outcome from 3 datasets and see similar ----
load("RStudio_inputs/myTopHits.treat") # From Amorim et al. 2019
load("RStudio_inputs/myTopHits.OP2") # From Christensen et al. 2016
head(myTopHits_treat) # Current Microbiome project

tophitsTreat_amorim <- myTopHits.treat %>%
  filter(P.Value < .05) %>%
  arrange(desc(logFC))
tophitsTreat_amorim <- rbind(head(tophitsTreat_amorim,100), tail(tophitsTreat_amorim, 100))
tophitsTreat_amorim <- head(tophitsTreat_amorim,500)

tophitsTreat_Chris <- myTopHits.OP2 %>%
  filter(P.Value < .05) %>%
  arrange(desc(logFC))
tophitsTreat_Chris <- rbind(head(tophitsTreat_Chris,100), tail(tophitsTreat_Chris, 100))
tophitsTreat_Chris <- head(tophitsTreat_Chris,500)

tophitsTreat_Current <- myTopHits_treat %>%
  filter(P.Value < .05) %>%
  arrange(desc(logFC))
tophitsTreat_Current <- rbind(head(tophitsTreat_Current,100), tail(tophitsTreat_Current, 100))
tophitsTreat_Current <- head(tophitsTreat_Current,500)

venn_biomarkers <- venn(list(Mosser2016_topHits = tophitsTreat_amorim$geneID,
                             Amorim2019_topHits = tophitsTreat_Chris$geneID,
                             Current_topHits = tophitsTreat_Current$geneSymbol))
venn_biomarkers
