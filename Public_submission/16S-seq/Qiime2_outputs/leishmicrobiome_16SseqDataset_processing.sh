# This is a shell script reporting the processing of going over the 16S sequencing dataset from the Leish-Microbiome project 2020.
# Camila Farias Amorim

# Ps. PATHS are not updated.

# Intro here ----
# These are the steps for the processing of the 16S-seq dataset associated with the 3rd lesion dataset 2019-2020.
# We used the CHMI server in november, 2020-21.
# Update: some samples were removed, and a re-processing with samples removed directly from the metadata file happenned in April 2021.

# We followed the workflow described in the CHMI website: https://www.notion.so/Analysis-of-16S-data-using-QIIME2-f225c95e34784cc3b726cff2d4d137cc
# This workflow was written by Alex Berry 2019.
# As well as parts of the Qiime1-2 tutorials.

# We had in our hands:
# 2 index files named: Undetermined_S0_L001_I1_001.fastq.gz and Undetermined_S0_L001_I2_001.fastq.gz
# 2 Multiplexed files named: Undetermined_S0_L001_R1_001.fastq.gz and Undetermined_S0_L001_R2_001.fastq.gz
# 1 metadata file named: MiSeqV1V3_38_humanLesion_metadata.tsv - this file was modified MANUALLY by Qi from the Grice lab to add the colunms of interest: type, day of collection, patient ID.
# Updated metadata file named: MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv 
# I (Camila) did not check the accuracy of this metadata file, since is too big. I will accept that it is correct.
# Associations with the metadata will only be important in really downstream analysis, so the only columns that should be correct is patient ID, and that looks correct.

# Important ps given by Alex Berry ----
# qiime 1 and 2 store temp files in the conda environment, and since ours in the CHMI linux is small we can set up a temp folder anywhere else.
# See his example, and run this everytime before working with qiime:
#export TMPDIR=/publicData/amorimc/tempFiles_Camila
#echo $TMPDIR

# Quick change because the server was updated:
export TMPDIR=/rome/amorimc/tempFiles_Camila
echo $TMPDIR

# Command lines in the CHMI server QIIME1 PART1 ----
# Initial exploration:
source activate qiime2
qiime tools inspect-metadata OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv

# Generate the .qzv file, which is the right format for Qiime
qiime metadata tabulate --m-input-file OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv --o-visualization qiimeProcessingFolder/tabulated-metadata.qzv

# Extract barcodes from Index files:
# Based on Qiime1: http://qiime.org/tutorials/processing_illumina_data.html#two-index-barcode-reads-and-two-fastq-reads
# The extract_barcodes.py is from Qiime1
# Extracting barcodes happenned with the Index files because of the format of the library preparation, which is a bit different from the CHMI tutorial.
source activate qiime1
extract_barcodes.py --input_type barcode_paired_end -f OriginalFiles/Undetermined_S0_L001_I1_001.fastq.gz -r OriginalFiles/Undetermined_S0_L001_I2_001.fastq.gz --bc1_len 12 --bc2_len 12 -o parsed_barcodes/
# --bc1/2_len of 12 because this is the number of bases in each foward and reverse included in the metadata.

# Move the barcodes.fastq to the folder that we will be working, and gzip it:
cp parsed_barcodes/barcodes.fastq qiimeProcessingFolder/barcodes.fastq
gzip qiimeProcessingFolder/barcodes.fastq
chmod u+x barcodes.fastq.gz # Not sure if it was necessary to change permissions, but I did it anyway.

# Also move the R1 and R2 original files to the same working directory and rename them as forward/reverse.fastq.gz because the script only accepts this way.
cp OriginalFiles/Undetermined_S0_L001_R1_001.fastq.gz qiimeProcessingFolder/forward.fastq.gz
cp OriginalFiles/Undetermined_S0_L001_R2_001.fastq.gz qiimeProcessingFolder/reverse.fastq.gz

# Making individual fastq per sample (FOR SRA SUBMISSION):
# Forward reads:
split_libraries_fastq.py \
  -i qiimeProcessingFolder/forward.fastq.gz \
  -o qiimeProcessingFolder/sra_submission \
  -m OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv \
  -b qiimeProcessingFolder/barcodes.fastq.gz \
  --barcode_type 16 \
  --store_demultiplexed_fastq #Without this flag the output is only a demultiplexed fasta file

# renaming stuff
mv qiimeProcessingFolder/sra_submission/seqs.fastq qiimeProcessingFolder/sra_submission/LeishMicro_forward_seqs.fastq
mv qiimeProcessingFolder/sra_submission/seqs.fna qiimeProcessingFolder/sra_submission/LeishMicro_forward.fna
mv qiimeProcessingFolder/sra_submission/histograms.txt qiimeProcessingFolder/sra_submission/LeishMicro_forward_histograms.txt
mv qiimeProcessingFolder/sra_submission/split_library_log.txt qiimeProcessingFolder/sra_submission/LeishMicro_forward_split_library_log.txt

split_sequence_file_on_sample_ids.py \
  -i qiimeProcessingFolder/sra_submission/LeishMicro_forward_seqs.fastq \
  -o qiimeProcessingFolder/sra_submission/forward_fastqs \
  --file_type fastq

cd qiimeProcessingFolder/sra_submission/forward_fastqs
for file in *.fastq; do mv "$file" "${file%.fastq_F.fastq}_F.fastq"; done
cd /rome/amorimc/LeishMicrobiome16SdatasetProcessing2020/JordansFiles

# Reverse now:
split_libraries_fastq.py \
  -i qiimeProcessingFolder/reverse.fastq.gz \
  -o qiimeProcessingFolder/sra_submission \
  -m OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv \
  -b qiimeProcessingFolder/barcodes.fastq.gz \
  --barcode_type 16 \
  --store_demultiplexed_fastq #Without this flag the output is only a demultiplexed fasta file

# renaming stuff
mv qiimeProcessingFolder/sra_submission/seqs.fastq qiimeProcessingFolder/sra_submission/LeishMicro_reverse_seqs.fastq
mv qiimeProcessingFolder/sra_submission/seqs.fna qiimeProcessingFolder/sra_submission/LeishMicro_reverse.fna
mv qiimeProcessingFolder/sra_submission/histograms.txt qiimeProcessingFolder/sra_submission/LeishMicro_reverse_histograms.txt
mv qiimeProcessingFolder/sra_submission/split_library_log.txt qiimeProcessingFolder/sra_submission/LeishMicro_reverse_split_library_log.txt

split_sequence_file_on_sample_ids.py \
  -i qiimeProcessingFolder/sra_submission/LeishMicro_reverse_seqs.fastq \
  -o qiimeProcessingFolder/sra_submission/reverse_fastqs \
  --file_type fastq

cd qiimeProcessingFolder/sra_submission/reverse_fastqs
for file in *.fastq; do mv "$file" "${file%.fastq_R.fastq}_R.fastq"; done
cd /rome/amorimc/LeishMicrobiome16SdatasetProcessing2020/JordansFiles


# Create an artifact of our data with QIIME2 PART2 ANALYSIS ----
# In this folder, have only the forward/reverse.fastq.gz and barcodes.fastq.gz
source activate qiime2
qiime tools import \
  --type EMPPairedEndSequences \
  --input-path /home/amorimc/LeishMicrobiome16SdatasetProcessing2020/JordansFiles/qiimeProcessingFolder/ \
  --output-path /home/amorimc/LeishMicrobiome16SdatasetProcessing2020/JordansFiles/qiimeProcessingFolder/emp-paired-end-sequences.qza

# Check if qiime2 recognizes your imports:
qiime tools peek qiimeProcessingFolder/emp-paired-end-sequences.qza

#I got this output:
#UUID:        bc78eef3-4e53-447d-b404-82d56d98cf4c
#Type:        EMPPairedEndSequences
#Data format: EMPPairedEndDirFmt

# Demultiplexing ----
qiime demux emp-paired \
  --i-seqs qiimeProcessingFolder/emp-paired-end-sequences.qza \
  --m-barcodes-file OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv \
  --m-barcodes-column BarcodeSequence \
  --o-error-correction-details qiimeProcessingFolder/demux-details.qza  \
  --o-per-sample-sequences qiimeProcessingFolder/demux.qza \
  --p-no-golay-error-correction

# Visualizing output:
qiime demux summarize \
  --i-data qiimeProcessingFolder/demux.qza \
  --o-visualization qiimeProcessingFolder/demux.qzv


# Denoising and QC filtering ----
# These cutoffs were not super restringent.
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiimeProcessingFolder/demux.qza \
  --p-trunc-len-f 285 \
  --p-trunc-len-r 247 \
  --o-representative-sequences qiimeProcessingFolder/rep-seqs-dada2.qza \
  --o-table qiimeProcessingFolder/table.qza \
  --o-denoising-stats qiimeProcessingFolder/denoising-stats.qza \
  --p-n-threads 20

# Restringent ones Having bottom of the Box = 20 QC from the demux.qzv
# I didn't go with this cutoff because too many samples didn't pass it >20 
#qiime dada2 denoise-paired \
#  --i-demultiplexed-seqs qiimeProcessingFolder/demux.qza \
#  --p-trunc-len-f 260 \
#  --p-trunc-len-r 214 \
#  --o-representative-sequences qiimeProcessingFolder/rep-seqs-dada2.qza \
#  --o-table qiimeProcessingFolder/pet-table.qza \
#  --o-denoising-stats qiimeProcessingFolder/denoising-stats.qza \
#  --p-n-threads 20

#  Summarize filtered/denoised data:
## Feature table:
qiime feature-table summarize \
  --i-table qiimeProcessingFolder/pet-table.qza \
  --o-visualization qiimeProcessingFolder/table.qzv \
  --m-sample-metadata-file OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv

## Feature sequences:
qiime feature-table tabulate-seqs \
  --i-data qiimeProcessingFolder/rep-seqs-dada2.qza \
  --o-visualization qiimeProcessingFolder/rep-seqs.qzv

# Summary of denoising stats:
qiime metadata tabulate \
  --m-input-file qiimeProcessingFolder/denoising-stats.qza \
  --o-visualization qiimeProcessingFolder/denoising-stats.qzv

# Build a phylogenetic tree ----
#carry out a multiple sequence alignment using Mafft
qiime alignment mafft \
  --i-sequences qiimeProcessingFolder/rep-seqs-dada2.qza \
  --o-alignment qiimeProcessingFolder/aligned-rep-seqs.qza

#mask (or filter) the alignment to remove positions that are highly variable. These positions are generally considered to add noise to a resulting phylogenetic tree.
qiime alignment mask \
  --i-alignment qiimeProcessingFolder/aligned-rep-seqs.qza \
  --o-masked-alignment qiimeProcessingFolder/masked-aligned-rep-seqs.qza

#create the tree using the Fasttree program
qiime phylogeny fasttree \
  --i-alignment qiimeProcessingFolder/masked-aligned-rep-seqs.qza \
  --o-tree qiimeProcessingFolder/unrooted-tree.qza

#root the tree using the longest root
qiime phylogeny midpoint-root \
  --i-tree qiimeProcessingFolder/unrooted-tree.qza \
  --o-rooted-tree qiimeProcessingFolder/rooted-tree.qza

# Alpha rarefaction
# max-depth of 17761 because this was the sample with the highest level of feature count.
qiime diversity alpha-rarefaction \
--i-table qiimeProcessingFolder/pet-table.qza \
--i-phylogeny qiimeProcessingFolder/rooted-tree.qza \
--p-max-depth 17761 \
--m-metadata-file OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv \
--o-visualization qiimeProcessingFolder/alpha-rarefaction.qzv

# Calculate and explore diversity metrics (all at once)
## Setting up sampling-depth to 896 because this was the samples with values above the environmental control.
## There were samples with poor sequencing, those will be excluded from these calculations.
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny qiimeProcessingFolder/rooted-tree.qza \
  --i-table qiimeProcessingFolder/pet-table.qza \
  --p-sampling-depth 896 \
  --m-metadata-file OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv \
  --output-dir qiimeProcessingFolder/core-metrics-results

# and vizualization of diversity metrics (alpha):
## I didn't run this one cause I don't know what is faith_pd:
qiime diversity alpha-group-significance \
  --i-alpha-diversity qiimeProcessingFolder/core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv \
  --o-visualization qiimeProcessingFolder/core-metrics-results/viz/faith-pd-group-significance.qzv

## I didn't run this one cause I don't know what is evenness:
qiime diversity alpha-group-significance \
  --i-alpha-diversity qiimeProcessingFolder/core-metrics-results/evenness_vector.qza \
  --m-metadata-file OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv \
  --o-visualization qiimeProcessingFolder/core-metrics-results/viz/evenness-group-significance.qzv

## I ran shannon:
qiime diversity alpha-group-significance \
  --i-alpha-diversity qiimeProcessingFolder/core-metrics-results/shannon_vector.qza \
  --m-metadata-file OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv \
  --o-visualization qiimeProcessingFolder/core-metrics-results/viz/shannon_group-significance.qzv

# and vizualization of diversity metrics (beta): --> Didnt do this yet
#qiime diversity beta-group-significance \
#  --i-distance-matrix qiimeProcessingFolder/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
#  --m-metadata-file OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv \
#  --m-metadata-column SubjectID \
#  --o-visualization qiimeProcessingFolder/core-metrics-results/viz/unweighted-unifrac-body-site-significance.qzv \
#  --p-pairwise

# PCoA to explore beta diversity metric:
#first, use the unweighted unifrac data as input
qiime emperor plot \
  --i-pcoa qiimeProcessingFolder/core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv \
  --o-visualization qiimeProcessingFolder/core-metrics-results/viz/unweighted-unifrac-emperor.qzv

qiime tools export \
  --input-path qiimeProcessingFolder/core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --output-path qiimeProcessingFolder/core-metrics-results/viz/exported-feature-table

# Assign taxonomy - GreenGenes:
# I had problems with downloading the taxonomy GreenGenes classifier. Anyway, my intention was to use Silva's newest classifier.
#wget -O "qiimeProcessingFolder/gg-13-8-99-515-806-nb-classifier.qza" "https://data.qiime2.org/2018.2/common/gg-13-8-99-515-806-nb-classifier.qza"
#wget -O "qiimeProcessingFolder/gg_13_8_otus.tar.gz" "ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz"

# I downloaded Silva's post-processed classifier from this page: https://docs.qiime2.org/2020.11/data-resources/
qiime tools peek qiimeProcessingFolder/silva-138-99-nb-classifier.qza

# this part took forever, consider next time to add thread argument
qiime feature-classifier classify-sklearn \
  --i-classifier qiimeProcessingFolder/silva-138-99-nb-classifier.qza \
  --i-reads qiimeProcessingFolder/rep-seqs-dada2.qza \
  --o-classification qiimeProcessingFolder/taxonomy.qza

qiime metadata tabulate \
  --m-input-file qiimeProcessingFolder/taxonomy.qza \
  --o-visualization qiimeProcessingFolder/taxonomy.qzv

qiime taxa barplot \
  --i-table qiimeProcessingFolder/pet-table.qza \
  --i-taxonomy qiimeProcessingFolder/taxonomy.qza \
  --m-metadata-file OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv \
  --o-visualization qiimeProcessingFolder/taxa-bar-plots.qzv

# From here, the .qzv and .qza objects were imported into R.

# Nmit for patients (copied from https://docs.qiime2.org/2018.11/tutorials/longitudinal/)

# Filter only Staphylococcus dysbiosis samples (more than 2 time points):
qiime feature-table filter-samples \
  --i-table qiimeProcessingFolder/pet-table.qza \
  --m-metadata-file OriginalFiles/StaphyDysbiosis_metadata.tsv \
  --p-where "StaphyDysbiosisincludenmit='yes'" \
  --o-filtered-table qiimeProcessingFolder/StaphyDysbiosis_pet-table.qza

# Just to vizualize:
qiime feature-table summarize \
  --i-table qiimeProcessingFolder/StaphyDysbiosis_pet-table.qza \
  --o-visualization qiimeProcessingFolder/StaphyDysbiosis_pet-table.qzv

# We use the taxa table collapsed to the genus level
qiime taxa collapse \
  --i-table qiimeProcessingFolder/StaphyDysbiosis_pet-table.qza \
  --i-taxonomy qiimeProcessingFolder/taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table qiimeProcessingFolder/StaphyDysbiosis_genus_table.qza

# Relative Feature table:
qiime feature-table relative-frequency \
  --i-table qiimeProcessingFolder/StaphyDysbiosis_genus_table.qza \
  --o-relative-frequency-table qiimeProcessingFolder/StaphyDysbiosis_relative_genus_table.qza

#Perform the NMIT
qiime longitudinal nmit \
  --i-table qiimeProcessingFolder/StaphyDysbiosis_relative_genus_table.qza \
  --m-metadata-file OriginalFiles/StaphyDysbiosis_metadata.tsv \
  --p-individual-id-column SubjectID \
  --p-corr-method pearson \
  --o-distance-matrix qiimeProcessingFolder/StaphySamples_nmit-dm.qza

# Include time point beta distance comparison
qiime diversity beta-group-significance \
  --i-distance-matrix qiimeProcessingFolder/StaphySamples_nmit-dm.qza \
  --m-metadata-file OriginalFiles/StaphyDysbiosis_metadata.tsv \
  --m-metadata-column Time_point \
  --o-visualization qiimeProcessingFolder/StaphySamples_nmit.qzv

#-----
# Filter all samples with more than 2 time points:
qiime feature-table filter-samples \
  --i-table qiimeProcessingFolder/pet-table.qza \
  --m-metadata-file OriginalFiles/sampDysbiosis_metadata.tsv \
  --p-where "sampDysbiosisincludenmit='yes'" \
  --o-filtered-table qiimeProcessingFolder/sampDysbiosis_pet-table.qza

# Just to vizualize:
qiime feature-table summarize \
  --i-table qiimeProcessingFolder/sampDysbiosis_pet-table.qza \
  --o-visualization qiimeProcessingFolder/sampDysbiosis_pet-table.qzv

# We use the taxa table collapsed to the genus level
qiime taxa collapse \
  --i-table qiimeProcessingFolder/sampDysbiosis_pet-table.qza \
  --i-taxonomy qiimeProcessingFolder/taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table qiimeProcessingFolder/sampDysbiosis_genus_table.qza

# Relative Feature table:
qiime feature-table relative-frequency \
  --i-table qiimeProcessingFolder/sampDysbiosis_genus_table.qza \
  --o-relative-frequency-table qiimeProcessingFolder/sampDysbiosis_relative_genus_table.qza

#Perform the NMIT
qiime longitudinal nmit \
  --i-table qiimeProcessingFolder/sampDysbiosis_relative_genus_table.qza \
  --m-metadata-file OriginalFiles/sampDysbiosis_metadata.tsv \
  --p-individual-id-column SubjectID \
  --p-corr-method pearson \
  --o-distance-matrix qiimeProcessingFolder/sampSamples_nmit-dm.qza

# Include time point beta distance comparison
qiime diversity beta-group-significance \
  --i-distance-matrix qiimeProcessingFolder/sampSamples_nmit-dm.qza \
  --m-metadata-file OriginalFiles/sampDysbiosis_metadata.tsv \
  --m-metadata-column Time_point \
  --o-visualization qiimeProcessingFolder/sampSamples_nmit.qzv

# The Pcoa object
qiime diversity pcoa \
  --i-distance-matrix qiimeProcessingFolder/sampSamples_nmit-dm.qza \
  --o-pcoa qiimeProcessingFolder/sampSamples_nmit-pc.qza


# This didnt work and I got the following error message:
#Plugin error from diversity:
#All values in the grouping vector are the same. This method cannot operate on a grouping vector with only a single group of objects (e.g., there are no 'between' distances because there is only a single group).
#Debug info has been saved to /publicData/amorimc/tempFiles_Camila/qiime2-q2cli-err-0s4wx3a7.log


# The Pcoa object
qiime diversity pcoa \
  --i-distance-matrix qiimeProcessingFolder/StaphySamples_nmit-dm.qza \
  --o-pcoa qiimeProcessingFolder/StaphySamples_nmit-pc.qza


# Every new session chunk ----
cd /publicData/amorimc/LeishMicrobiome16SdatasetProcessing2020/JordansFiles
export TMPDIR=/publicData/amorimc/tempFiles_Camila
echo $TMPDIR
source activate qiime2


# Convert demultiplexed objects to individual fastq for SRA submission:
# PATHS ARE CHANGED!!
split_libraries_fastq.py \
  -i forward.fastq.gz \
  -o sra_submission \
  -m OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv \
  -b barcodes.fastq.gz \
  --barcode_type 16 \
  --store_demultiplexed_fastq #Without this flag the output is only a demultiplexed fasta file

mv seqs.fastq LeishMicro_forward_seqs.fastq
mv seqs.fna LeishMicro_forward.fna
mv histograms.txt LeishMicro_orward_histograms.txt
mv split_library_log.txt LeishMicro_forward_split_library_log.txt

split_sequence_file_on_sample_ids.py \
  -i LeishMicro_forward_seqs.fastq \
  -o forward_fastqs \
  --file_type fastq


for file in *.fastq; do mv "$file" "${file%.fastq_F.fastq}_F.fastq"; done



###   Reverse   #####

split_libraries_fastq.py \
  -i reverse.fastq.gz \
  -o sra_submission \
  -m OriginalFiles/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv \
  -b barcodes.fastq.gz \
  --barcode_type 16 \
  --store_demultiplexed_fastq #Without this flag the output is only a demultiplexed fasta file

mv seqs.fastq LeishMicro_reverse_seqs.fastq
mv seqs.fna LeishMicro_reverse.fna
mv histograms.txt rLeishMicro_everse_histograms.txt
mv split_library_log.txt LeishMicro_reverse_split_library_log.txt

split_sequence_file_on_sample_ids.py \
  -i rLeishMicro_reverse_seqs.fastq \
  -o reverse_fastqs \
  --file_type fastq







