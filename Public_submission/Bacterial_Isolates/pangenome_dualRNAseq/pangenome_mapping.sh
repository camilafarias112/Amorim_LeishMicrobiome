# Mapping reads from the 3rd lesion CL RNA-seq dataset to a Staphylococcus aureus pan-genome reference.
# The pan-genome was put together by Tori Lovins using the Roary (https://sanger-pathogens.github.io/Roary/) workflow. This pan-genome includes genomes from our isolates and known available S. aureus strains.
# First I filtered the non-human reads with kneadData and used those.
# Secondly, performed the pangenome mapping with Bowtie2.

# Temporary folder
#export TMPDIR=/publicData/amorimc/tempFiles_Camila
#echo $TMPDIR

#bowtie2-build --threads 20 -f pan_genome_reference.fa Saureus_pangenome_index

# Mapping:
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29075_kneaddata.fastq -S Saureus_pangenome_outputs/CL29075_Saureus.sam >& Saureus_pangenome_outputs/CL29075_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29107_kneaddata.fastq -S Saureus_pangenome_outputs/CL29107_Saureus.sam >& Saureus_pangenome_outputs/CL29107_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29077_kneaddata.fastq -S Saureus_pangenome_outputs/CL29077_Saureus.sam >& Saureus_pangenome_outputs/CL29077_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29151_kneaddata.fastq -S Saureus_pangenome_outputs/CL29151_Saureus.sam >& Saureus_pangenome_outputs/CL29151_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29148_kneaddata.fastq -S Saureus_pangenome_outputs/CL29148_Saureus.sam >& Saureus_pangenome_outputs/CL29148_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29150_kneaddata.fastq -S Saureus_pangenome_outputs/CL29150_Saureus.sam >& Saureus_pangenome_outputs/CL29150_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29085_kneaddata.fastq -S Saureus_pangenome_outputs/CL29085_Saureus.sam >& Saureus_pangenome_outputs/CL29085_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29118_kneaddata.fastq -S Saureus_pangenome_outputs/CL29118_Saureus.sam >& Saureus_pangenome_outputs/CL29118_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29137_kneaddata.fastq -S Saureus_pangenome_outputs/CL29137_Saureus.sam >& Saureus_pangenome_outputs/CL29137_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29161_kneaddata.fastq -S Saureus_pangenome_outputs/CL29161_Saureus.sam >& Saureus_pangenome_outputs/CL29161_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29170_kneaddata.fastq -S Saureus_pangenome_outputs/CL29170_Saureus.sam >& Saureus_pangenome_outputs/CL29170_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29180_kneaddata.fastq -S Saureus_pangenome_outputs/CL29180_Saureus.sam >& Saureus_pangenome_outputs/CL29180_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29159_kneaddata.fastq -S Saureus_pangenome_outputs/CL29159_Saureus.sam >& Saureus_pangenome_outputs/CL29159_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29168_kneaddata.fastq -S Saureus_pangenome_outputs/CL29168_Saureus.sam >& Saureus_pangenome_outputs/CL29168_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29213_kneaddata.fastq -S Saureus_pangenome_outputs/CL29213_Saureus.sam >& Saureus_pangenome_outputs/CL29213_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29231_kneaddata.fastq -S Saureus_pangenome_outputs/CL29231_Saureus.sam >& Saureus_pangenome_outputs/CL29231_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29201_kneaddata.fastq -S Saureus_pangenome_outputs/CL29201_Saureus.sam >& Saureus_pangenome_outputs/CL29201_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29219_kneaddata.fastq -S Saureus_pangenome_outputs/CL29219_Saureus.sam >& Saureus_pangenome_outputs/CL29219_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29234_kneaddata.fastq -S Saureus_pangenome_outputs/CL29234_Saureus.sam >& Saureus_pangenome_outputs/CL29234_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29246_kneaddata.fastq -S Saureus_pangenome_outputs/CL29246_Saureus.sam >& Saureus_pangenome_outputs/CL29246_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29223_kneaddata.fastq -S Saureus_pangenome_outputs/CL29223_Saureus.sam >& Saureus_pangenome_outputs/CL29223_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29218_kneaddata.fastq -S Saureus_pangenome_outputs/CL29218_Saureus.sam >& Saureus_pangenome_outputs/CL29218_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29253_kneaddata.fastq -S Saureus_pangenome_outputs/CL29253_Saureus.sam >& Saureus_pangenome_outputs/CL29253_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29254_kneaddata.fastq -S Saureus_pangenome_outputs/CL29254_Saureus.sam >& Saureus_pangenome_outputs/CL29254_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29256_kneaddata.fastq -S Saureus_pangenome_outputs/CL29256_Saureus.sam >& Saureus_pangenome_outputs/CL29256_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29257_kneaddata.fastq -S Saureus_pangenome_outputs/CL29257_Saureus.sam >& Saureus_pangenome_outputs/CL29257_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29258_kneaddata.fastq -S Saureus_pangenome_outputs/CL29258_Saureus.sam >& Saureus_pangenome_outputs/CL29258_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29264_kneaddata.fastq -S Saureus_pangenome_outputs/CL29264_Saureus.sam >& Saureus_pangenome_outputs/CL29264_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29270_kneaddata.fastq -S Saureus_pangenome_outputs/CL29270_Saureus.sam >& Saureus_pangenome_outputs/CL29270_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29271_kneaddata.fastq -S Saureus_pangenome_outputs/CL29271_Saureus.sam >& Saureus_pangenome_outputs/CL29271_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29272_kneaddata.fastq -S Saureus_pangenome_outputs/CL29272_Saureus.sam >& Saureus_pangenome_outputs/CL29272_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29274_kneaddata.fastq -S Saureus_pangenome_outputs/CL29274_Saureus.sam >& Saureus_pangenome_outputs/CL29274_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29276_kneaddata.fastq -S Saureus_pangenome_outputs/CL29276_Saureus.sam >& Saureus_pangenome_outputs/CL29276_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29289_kneaddata.fastq -S Saureus_pangenome_outputs/CL29289_Saureus.sam >& Saureus_pangenome_outputs/CL29289_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29291_kneaddata.fastq -S Saureus_pangenome_outputs/CL29291_Saureus.sam >& Saureus_pangenome_outputs/CL29291_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29297_kneaddata.fastq -S Saureus_pangenome_outputs/CL29297_Saureus.sam >& Saureus_pangenome_outputs/CL29297_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29298_kneaddata.fastq -S Saureus_pangenome_outputs/CL29298_Saureus.sam >& Saureus_pangenome_outputs/CL29298_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29303_kneaddata.fastq -S Saureus_pangenome_outputs/CL29303_Saureus.sam >& Saureus_pangenome_outputs/CL29303_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29491_kneaddata.fastq -S Saureus_pangenome_outputs/CL29491_Saureus.sam >& Saureus_pangenome_outputs/CL29491_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29498_kneaddata.fastq -S Saureus_pangenome_outputs/CL29498_Saureus.sam >& Saureus_pangenome_outputs/CL29498_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29513_kneaddata.fastq -S Saureus_pangenome_outputs/CL29513_Saureus.sam >& Saureus_pangenome_outputs/CL29513_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29519_kneaddata.fastq -S Saureus_pangenome_outputs/CL29519_Saureus.sam >& Saureus_pangenome_outputs/CL29519_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29529_kneaddata.fastq -S Saureus_pangenome_outputs/CL29529_Saureus.sam >& Saureus_pangenome_outputs/CL29529_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29539_kneaddata.fastq -S Saureus_pangenome_outputs/CL29539_Saureus.sam >& Saureus_pangenome_outputs/CL29539_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29558_kneaddata.fastq -S Saureus_pangenome_outputs/CL29558_Saureus.sam >& Saureus_pangenome_outputs/CL29558_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29663_kneaddata.fastq -S Saureus_pangenome_outputs/CL29663_Saureus.sam >& Saureus_pangenome_outputs/CL29663_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29680_kneaddata.fastq -S Saureus_pangenome_outputs/CL29680_Saureus.sam >& Saureus_pangenome_outputs/CL29680_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29687_kneaddata.fastq -S Saureus_pangenome_outputs/CL29687_Saureus.sam >& Saureus_pangenome_outputs/CL29687_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29727_kneaddata.fastq -S Saureus_pangenome_outputs/CL29727_Saureus.sam >& Saureus_pangenome_outputs/CL29727_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29770_kneaddata.fastq -S Saureus_pangenome_outputs/CL29770_Saureus.sam >& Saureus_pangenome_outputs/CL29770_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/CL29776_kneaddata.fastq -S Saureus_pangenome_outputs/CL29776_Saureus.sam >& Saureus_pangenome_outputs/CL29776_Saureus.txt

bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/HS1_kneaddata.fastq -S Saureus_pangenome_outputs/HS1_Saureus.sam >& Saureus_pangenome_outputs/HS1_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/HS2_kneaddata.fastq -S Saureus_pangenome_outputs/HS2_Saureus.sam >& Saureus_pangenome_outputs/HS2_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/HS3_kneaddata.fastq -S Saureus_pangenome_outputs/HS3_Saureus.sam >& Saureus_pangenome_outputs/HS3_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/HS4_kneaddata.fastq -S Saureus_pangenome_outputs/HS4_Saureus.sam >& Saureus_pangenome_outputs/HS4_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/HS5_kneaddata.fastq -S Saureus_pangenome_outputs/HS5_Saureus.sam >& Saureus_pangenome_outputs/HS5_Saureus.txt
bowtie2 -x Saureus_pangenome_index -U /publicData/amorimc/3rdDataset/kneaddata_nonhuman_fastq/HS6_kneaddata.fastq -S Saureus_pangenome_outputs/HS6_Saureus.sam >& Saureus_pangenome_outputs/HS6_Saureus.txt

# MultiQC:
multiqc -d Saureus_pangenome_outputs/.

echo "Finished"



# Make a bam file:
# samtools view -bS CL29075_Saureus.sam > CL29075_Saureus.bam
# samtools sort CL29075_Saureus.bam -o CL29075_Saureus.sorted.bam


