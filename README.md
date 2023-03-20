# Farias Amorim et al. 2023

> This is the R code associated with the manuscript:
<p align="center"><strong>Multi-omic profiling of cutaneous leishmaniasis infections reveals microbiota-driven mechanisms underlying disease severity</strong></p>
<p align="center"><em>Camila Farias Amorim, Victoria M. Lovins, Tej Pratap Singh, Fernanda O. Novais, Jordan C. Harris, Alexsandro S. Lago, Lucas P. Carvalho, Edgar M. Carvalho, Daniel P. Beiting, Phillip Scott*, Elizabeth A. Grice*</em></p>

Code for my publication "Multi-omic profiling of cutaneous leishmaniasis infections reveals microbiota-driven mechanisms underlying disease severity", 2023

<p>[medRxiv](https://www.medrxiv.org/content/10.1101/2023.02.02.23285247v2)</p>

## Abstract
<i>Leishmania braziliensis</i> infection results in inflammation and skin injury, with highly variable and unpredictable clinical outcomes. Here, we investigated the potential impact of microbiota on infection-induced inflammatory responses and disease resolution by conducting an integrated analysis of the skin microbiome and host transcriptome on a cohort of 62 <i>L. braziliensis</i>-infected patients. We found that overall bacterial burden and microbiome configurations dominated with <i>Staphylococcus spp.</i> were associated with delayed healing and enhanced inflammatory responses, especially by IL-1 family members. Dual RNA-seq of human lesions revealed that high lesional <i>S. aureus</i> transcript abundance was associated with delayed healing and increased expression of IL-1β.  This cytokine was critical for modulating disease outcome in <i>L. braziliensis</i>-infected mice colonized with <i>S. aureus</i>, as its neutralization reduced pathology and inflammation. These results implicate the microbiome in cutaneous leishmaniasis disease outcomes in humans and suggest host-directed therapies to mitigate the inflammatory consequences. 

<img align="center" width="1000" height="420" src="/Public_submission/MultipliedFactor.png">

> The locations of the core components of this repo are outlined in the file system map below. In short, there are the following main directories:

 - Public_submission - contains the main .Rscripts, main Robjects, tables, raw data associated this manuscript (with exception of raw sequencing files). The subdirectories included here are mostly divided per datasets or experiment (RNA-seq, 16S-seq, S. aureus isolates data, mice-experimental data, ...).

 - Amorim2022_SupplemmentalTable1_LeishOmics_StudyDesign.txt - Clinical metadata associated with the patients included in this study (n=62).
 - /Bacterial_isolates/ - Contains data about the live bacterial isolates collected from CL lesions. See isolation process and methodology in the official manuscript.
 - /Biopsy_qPCR/ - contains the raw data from the 2 qPCR experiments included in this manuscript: L. braziliensis's 18S and total bacteria's 16S ribosomal subunits. Ps.: In some parts there are still included the code for a failed qPCR trying to quantify specifically S. aureus. This experiment was never included in the manuscript.
 - /Integrative_Analysis/ - integration pipeline using the rexposome R package
 - /16S-seq/ - Data regarding the 16S-seq analyses performed with Qiime2 and final modeling in R.
 - /RNA-seq/ - Data regarding the RNA-seq dataset: from pre-processing, mapping with kallisto pseudoaligner, gene annotation and final modeling in R.
 - /Mice_experiments/ - Data regarding the experiments performed by Tej Singh. Here with blocked IL1-signaling in mice in the context of <i>S. aureus</i> co-infection.



```
