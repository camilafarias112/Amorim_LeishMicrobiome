# In this shell I performed the mapping with kallisto on the samples from the Microbiome project.
# These samples are biopsies from CL patients an biopsies from health individuals from the hospital here in Philly.

# Cat Runs:
#1
#cat CL29075-batch1_S57_L001_R1_001.fastq.gz CL29075-batch1_S57_L002_R1_001.fastq.gz CL29075-batch1_S57_L003_R1_001.fastq.gz CL29075-batch1_S57_L004_R1_001.fastq.gz > CL29075_1.fastq.gz 
#cat CL29107-batch1_S54_L001_R1_001.fastq.gz CL29107-batch1_S54_L002_R1_001.fastq.gz CL29107-batch1_S54_L003_R1_001.fastq.gz CL29107-batch1_S54_L004_R1_001.fastq.gz > CL29107_1.fastq.gz 
#cat CL29077-batch1_S56_L001_R1_001.fastq.gz CL29077-batch1_S56_L002_R1_001.fastq.gz CL29077-batch1_S56_L003_R1_001.fastq.gz CL29077-batch1_S56_L004_R1_001.fastq.gz > CL29077_1.fastq.gz 
#cat CL29151-batch1_S49_L001_R1_001.fastq.gz CL29151-batch1_S49_L002_R1_001.fastq.gz CL29151-batch1_S49_L003_R1_001.fastq.gz CL29151-batch1_S49_L004_R1_001.fastq.gz > CL29151_1.fastq.gz 
#cat CL29148-batch1_S51_L001_R1_001.fastq.gz CL29148-batch1_S51_L002_R1_001.fastq.gz CL29148-batch1_S51_L003_R1_001.fastq.gz CL29148-batch1_S51_L004_R1_001.fastq.gz > CL29148_1.fastq.gz 
#cat CL29150-batch1_S50_L001_R1_001.fastq.gz CL29150-batch1_S50_L002_R1_001.fastq.gz CL29150-batch1_S50_L003_R1_001.fastq.gz CL29150-batch1_S50_L004_R1_001.fastq.gz > CL29150_1.fastq.gz 
#cat CL29085-batch2_S55_L001_R1_001.fastq.gz CL29085-batch2_S55_L002_R1_001.fastq.gz CL29085-batch2_S55_L003_R1_001.fastq.gz CL29085-batch2_S55_L004_R1_001.fastq.gz > CL29085_1.fastq.gz 
#cat CL29118-batch1_S53_L001_R1_001.fastq.gz CL29118-batch1_S53_L002_R1_001.fastq.gz CL29118-batch1_S53_L003_R1_001.fastq.gz CL29118-batch1_S53_L004_R1_001.fastq.gz > CL29118_1.fastq.gz 
#cat CL29137-batch1_S52_L001_R1_001.fastq.gz CL29137-batch1_S52_L002_R1_001.fastq.gz CL29137-batch1_S52_L003_R1_001.fastq.gz CL29137-batch1_S52_L004_R1_001.fastq.gz > CL29137_1.fastq.gz 
#cat CL29161-batch1_S47_L001_R1_001.fastq.gz CL29161-batch1_S47_L002_R1_001.fastq.gz CL29161-batch1_S47_L003_R1_001.fastq.gz CL29161-batch1_S47_L004_R1_001.fastq.gz > CL29161_1.fastq.gz 
#cat CL29170-batch2_S45_L001_R1_001.fastq.gz CL29170-batch2_S45_L002_R1_001.fastq.gz CL29170-batch2_S45_L003_R1_001.fastq.gz CL29170-batch2_S45_L004_R1_001.fastq.gz > CL29170_1.fastq.gz 
#cat CL29180-batch2_S44_L001_R1_001.fastq.gz CL29180-batch2_S44_L002_R1_001.fastq.gz CL29180-batch2_S44_L003_R1_001.fastq.gz CL29180-batch2_S44_L004_R1_001.fastq.gz > CL29180_1.fastq.gz 
#cat CL29159-batch1_S48_L001_R1_001.fastq.gz CL29159-batch1_S48_L002_R1_001.fastq.gz CL29159-batch1_S48_L003_R1_001.fastq.gz CL29159-batch1_S48_L004_R1_001.fastq.gz > CL29159_1.fastq.gz 
#cat CL29168-batch2_S46_L001_R1_001.fastq.gz CL29168-batch2_S46_L002_R1_001.fastq.gz CL29168-batch2_S46_L003_R1_001.fastq.gz CL29168-batch2_S46_L004_R1_001.fastq.gz > CL29168_1.fastq.gz 
#cat CL29213-batch2_S42_L001_R1_001.fastq.gz CL29213-batch2_S42_L002_R1_001.fastq.gz CL29213-batch2_S42_L003_R1_001.fastq.gz CL29213-batch2_S42_L004_R1_001.fastq.gz > CL29213_1.fastq.gz 
#cat CL29231-batch2_S38_L001_R1_001.fastq.gz CL29231-batch2_S38_L002_R1_001.fastq.gz CL29231-batch2_S38_L003_R1_001.fastq.gz CL29231-batch2_S38_L004_R1_001.fastq.gz > CL29231_1.fastq.gz 
#cat CL29201-batch2_S43_L001_R1_001.fastq.gz CL29201-batch2_S43_L002_R1_001.fastq.gz CL29201-batch2_S43_L003_R1_001.fastq.gz CL29201-batch2_S43_L004_R1_001.fastq.gz > CL29201_1.fastq.gz 
#cat CL29219-batch2_S40_L001_R1_001.fastq.gz CL29219-batch2_S40_L002_R1_001.fastq.gz CL29219-batch2_S40_L003_R1_001.fastq.gz CL29219-batch2_S40_L004_R1_001.fastq.gz > CL29219_1.fastq.gz 
#cat CL29234-batch2_S37_L001_R1_001.fastq.gz CL29234-batch2_S37_L002_R1_001.fastq.gz CL29234-batch2_S37_L003_R1_001.fastq.gz CL29234-batch2_S37_L004_R1_001.fastq.gz > CL29234_1.fastq.gz 
#cat CL29246-batch2_S36_L001_R1_001.fastq.gz CL29246-batch2_S36_L002_R1_001.fastq.gz CL29246-batch2_S36_L003_R1_001.fastq.gz CL29246-batch2_S36_L004_R1_001.fastq.gz > CL29246_1.fastq.gz 
#cat CL29223-batch2_S39_L001_R1_001.fastq.gz CL29223-batch2_S39_L002_R1_001.fastq.gz CL29223-batch2_S39_L003_R1_001.fastq.gz CL29223-batch2_S39_L004_R1_001.fastq.gz > CL29223_1.fastq.gz 
#cat CL29218-batch2_S41_L001_R1_001.fastq.gz CL29218-batch2_S41_L002_R1_001.fastq.gz CL29218-batch2_S41_L003_R1_001.fastq.gz CL29218-batch2_S41_L004_R1_001.fastq.gz > CL29218_1.fastq.gz 
#cat CL29253-batch2_S35_L001_R1_001.fastq.gz CL29253-batch2_S35_L002_R1_001.fastq.gz CL29253-batch2_S35_L003_R1_001.fastq.gz CL29253-batch2_S35_L004_R1_001.fastq.gz > CL29253_1.fastq.gz 
#cat CL29254-batch2_S34_L001_R1_001.fastq.gz CL29254-batch2_S34_L002_R1_001.fastq.gz CL29254-batch2_S34_L003_R1_001.fastq.gz CL29254-batch2_S34_L004_R1_001.fastq.gz > CL29254_1.fastq.gz 
#cat CL29256-batch2_S33_L001_R1_001.fastq.gz CL29256-batch2_S33_L002_R1_001.fastq.gz CL29256-batch2_S33_L003_R1_001.fastq.gz CL29256-batch2_S33_L004_R1_001.fastq.gz > CL29256_1.fastq.gz 
#cat CL29257-batch2_S32_L001_R1_001.fastq.gz CL29257-batch2_S32_L002_R1_001.fastq.gz CL29257-batch2_S32_L003_R1_001.fastq.gz CL29257-batch2_S32_L004_R1_001.fastq.gz > CL29257_1.fastq.gz 
#cat CL29258-batch2_S31_L001_R1_001.fastq.gz CL29258-batch2_S31_L002_R1_001.fastq.gz CL29258-batch2_S31_L003_R1_001.fastq.gz CL29258-batch2_S31_L004_R1_001.fastq.gz > CL29258_1.fastq.gz 
#cat CL29264-batch2_S30_L001_R1_001.fastq.gz CL29264-batch2_S30_L002_R1_001.fastq.gz CL29264-batch2_S30_L003_R1_001.fastq.gz CL29264-batch2_S30_L004_R1_001.fastq.gz > CL29264_1.fastq.gz 
#cat CL29270-batch2_S29_L001_R1_001.fastq.gz CL29270-batch2_S29_L002_R1_001.fastq.gz CL29270-batch2_S29_L003_R1_001.fastq.gz CL29270-batch2_S29_L004_R1_001.fastq.gz > CL29270_1.fastq.gz 
#cat CL29271-batch3_S28_L001_R1_001.fastq.gz CL29271-batch3_S28_L002_R1_001.fastq.gz CL29271-batch3_S28_L003_R1_001.fastq.gz CL29271-batch3_S28_L004_R1_001.fastq.gz > CL29271_1.fastq.gz 
#cat CL29272-batch2_S27_L001_R1_001.fastq.gz CL29272-batch2_S27_L002_R1_001.fastq.gz CL29272-batch2_S27_L003_R1_001.fastq.gz CL29272-batch2_S27_L004_R1_001.fastq.gz > CL29272_1.fastq.gz 
#cat CL29274-batch3_S26_L001_R1_001.fastq.gz CL29274-batch3_S26_L002_R1_001.fastq.gz CL29274-batch3_S26_L003_R1_001.fastq.gz CL29274-batch3_S26_L004_R1_001.fastq.gz > CL29274_1.fastq.gz 
#cat CL29276-batch3_S25_L001_R1_001.fastq.gz CL29276-batch3_S25_L002_R1_001.fastq.gz CL29276-batch3_S25_L003_R1_001.fastq.gz CL29276-batch3_S25_L004_R1_001.fastq.gz > CL29276_1.fastq.gz 
#cat CL29289-batch3_S24_L001_R1_001.fastq.gz CL29289-batch3_S24_L002_R1_001.fastq.gz CL29289-batch3_S24_L003_R1_001.fastq.gz CL29289-batch3_S24_L004_R1_001.fastq.gz > CL29289_1.fastq.gz 
#cat CL29291-batch3_S23_L001_R1_001.fastq.gz CL29291-batch3_S23_L002_R1_001.fastq.gz CL29291-batch3_S23_L003_R1_001.fastq.gz CL29291-batch3_S23_L004_R1_001.fastq.gz > CL29291_1.fastq.gz 
#cat CL29297-batch3_S22_L001_R1_001.fastq.gz CL29297-batch3_S22_L002_R1_001.fastq.gz CL29297-batch3_S22_L003_R1_001.fastq.gz CL29297-batch3_S22_L004_R1_001.fastq.gz > CL29297_1.fastq.gz 
#cat CL29298-batch3_S21_L001_R1_001.fastq.gz CL29298-batch3_S21_L002_R1_001.fastq.gz CL29298-batch3_S21_L003_R1_001.fastq.gz CL29298-batch3_S21_L004_R1_001.fastq.gz > CL29298_1.fastq.gz 
#cat CL29303-batch3_S20_L001_R1_001.fastq.gz CL29303-batch3_S20_L002_R1_001.fastq.gz CL29303-batch3_S20_L003_R1_001.fastq.gz CL29303-batch3_S20_L004_R1_001.fastq.gz > CL29303_1.fastq.gz 
#cat CL29491-batch3_S19_L001_R1_001.fastq.gz CL29491-batch3_S19_L002_R1_001.fastq.gz CL29491-batch3_S19_L003_R1_001.fastq.gz CL29491-batch3_S19_L004_R1_001.fastq.gz > CL29491_1.fastq.gz 
#cat CL29498-batch3_S18_L001_R1_001.fastq.gz CL29498-batch3_S18_L002_R1_001.fastq.gz CL29498-batch3_S18_L003_R1_001.fastq.gz CL29498-batch3_S18_L004_R1_001.fastq.gz > CL29498_1.fastq.gz 
#cat CL29513-batch3_S17_L001_R1_001.fastq.gz CL29513-batch3_S17_L002_R1_001.fastq.gz CL29513-batch3_S17_L003_R1_001.fastq.gz CL29513-batch3_S17_L004_R1_001.fastq.gz > CL29513_1.fastq.gz 
#cat CL29519-batch3_S16_L001_R1_001.fastq.gz CL29519-batch3_S16_L002_R1_001.fastq.gz CL29519-batch3_S16_L003_R1_001.fastq.gz CL29519-batch3_S16_L004_R1_001.fastq.gz > CL29519_1.fastq.gz 
#cat CL29529-batch3_S15_L001_R1_001.fastq.gz CL29529-batch3_S15_L002_R1_001.fastq.gz CL29529-batch3_S15_L003_R1_001.fastq.gz CL29529-batch3_S15_L004_R1_001.fastq.gz > CL29529_1.fastq.gz 
#cat CL29539-batch3_S14_L001_R1_001.fastq.gz CL29539-batch3_S14_L002_R1_001.fastq.gz CL29539-batch3_S14_L003_R1_001.fastq.gz CL29539-batch3_S14_L004_R1_001.fastq.gz > CL29539_1.fastq.gz 
#cat CL29558-batch3_S13_L001_R1_001.fastq.gz CL29558-batch3_S13_L002_R1_001.fastq.gz CL29558-batch3_S13_L003_R1_001.fastq.gz CL29558-batch3_S13_L004_R1_001.fastq.gz > CL29558_1.fastq.gz 
#cat CL29663-batch3_S12_L001_R1_001.fastq.gz CL29663-batch3_S12_L002_R1_001.fastq.gz CL29663-batch3_S12_L003_R1_001.fastq.gz CL29663-batch3_S12_L004_R1_001.fastq.gz > CL29663_1.fastq.gz 
#cat CL29680-batch3_S11_L001_R1_001.fastq.gz CL29680-batch3_S11_L002_R1_001.fastq.gz CL29680-batch3_S11_L003_R1_001.fastq.gz CL29680-batch3_S11_L004_R1_001.fastq.gz > CL29680_1.fastq.gz 
#cat CL29687-batch3_S10_L001_R1_001.fastq.gz CL29687-batch3_S10_L002_R1_001.fastq.gz CL29687-batch3_S10_L003_R1_001.fastq.gz CL29687-batch3_S10_L004_R1_001.fastq.gz > CL29687_1.fastq.gz 
#cat CL29727-batch3_S9_L001_R1_001.fastq.gz CL29727-batch3_S9_L002_R1_001.fastq.gz CL29727-batch3_S9_L003_R1_001.fastq.gz CL29727-batch3_S9_L004_R1_001.fastq.gz > CL29727_1.fastq.gz 
#cat CL29770-batch3_S8_L001_R1_001.fastq.gz CL29770-batch3_S8_L002_R1_001.fastq.gz CL29770-batch3_S8_L003_R1_001.fastq.gz CL29770-batch3_S8_L004_R1_001.fastq.gz > CL29770_1.fastq.gz 
#cat CL29776-batch3_S7_L001_R1_001.fastq.gz CL29776-batch3_S7_L002_R1_001.fastq.gz CL29776-batch3_S7_L003_R1_001.fastq.gz CL29776-batch3_S7_L004_R1_001.fastq.gz > CL29776_1.fastq.gz 
#cat HS1-batch1_S6_L001_R1_001.fastq.gz HS1-batch1_S6_L002_R1_001.fastq.gz HS1-batch1_S6_L003_R1_001.fastq.gz HS1-batch1_S6_L004_R1_001.fastq.gz > HS1_1.fastq.gz 
#cat HS2-batch2_S5_L001_R1_001.fastq.gz HS2-batch2_S5_L002_R1_001.fastq.gz HS2-batch2_S5_L003_R1_001.fastq.gz HS2-batch2_S5_L004_R1_001.fastq.gz > HS2_1.fastq.gz 
#cat HS3-batch2_S4_L001_R1_001.fastq.gz HS3-batch2_S4_L002_R1_001.fastq.gz HS3-batch2_S4_L003_R1_001.fastq.gz HS3-batch2_S4_L004_R1_001.fastq.gz > HS3_1.fastq.gz 
#cat HS4-batch3_S3_L001_R1_001.fastq.gz HS4-batch3_S3_L002_R1_001.fastq.gz HS4-batch3_S3_L003_R1_001.fastq.gz HS4-batch3_S3_L004_R1_001.fastq.gz > HS4_1.fastq.gz 
#cat HS5-batch3_S2_L001_R1_001.fastq.gz HS5-batch3_S2_L002_R1_001.fastq.gz HS5-batch3_S2_L003_R1_001.fastq.gz HS5-batch3_S2_L004_R1_001.fastq.gz > HS5_1.fastq.gz 
#cat HS6-batch3_S1_L001_R1_001.fastq.gz HS6-batch3_S1_L002_R1_001.fastq.gz HS6-batch3_S1_L003_R1_001.fastq.gz HS6-batch3_S1_L004_R1_001.fastq.gz > HS6_1.fastq.gz #

##2
#cat CL29075-batch1_S57_L001_R1_001.fastq.gz CL29075-batch1_S57_L002_R1_001.fastq.gz CL29075-batch1_S57_L003_R1_001.fastq.gz CL29075-batch1_S57_L004_R1_001.fastq.gz > CL29075_2.fastq.gz 
#cat CL29107-batch1_S54_L001_R1_001.fastq.gz CL29107-batch1_S54_L002_R1_001.fastq.gz CL29107-batch1_S54_L003_R1_001.fastq.gz CL29107-batch1_S54_L004_R1_001.fastq.gz > CL29107_2.fastq.gz 
#cat CL29077-batch1_S56_L001_R1_001.fastq.gz CL29077-batch1_S56_L002_R1_001.fastq.gz CL29077-batch1_S56_L003_R1_001.fastq.gz CL29077-batch1_S56_L004_R1_001.fastq.gz > CL29077_2.fastq.gz 
#cat CL29151-batch1_S49_L001_R1_001.fastq.gz CL29151-batch1_S49_L002_R1_001.fastq.gz CL29151-batch1_S49_L003_R1_001.fastq.gz CL29151-batch1_S49_L004_R1_001.fastq.gz > CL29151_2.fastq.gz 
#cat CL29148-batch1_S51_L001_R1_001.fastq.gz CL29148-batch1_S51_L002_R1_001.fastq.gz CL29148-batch1_S51_L003_R1_001.fastq.gz CL29148-batch1_S51_L004_R1_001.fastq.gz > CL29148_2.fastq.gz 
#cat CL29150-batch1_S50_L001_R1_001.fastq.gz CL29150-batch1_S50_L002_R1_001.fastq.gz CL29150-batch1_S50_L003_R1_001.fastq.gz CL29150-batch1_S50_L004_R1_001.fastq.gz > CL29150_2.fastq.gz 
#cat CL29085-batch2_S55_L001_R1_001.fastq.gz CL29085-batch2_S55_L002_R1_001.fastq.gz CL29085-batch2_S55_L003_R1_001.fastq.gz CL29085-batch2_S55_L004_R1_001.fastq.gz > CL29085_2.fastq.gz 
#cat CL29118-batch1_S53_L001_R1_001.fastq.gz CL29118-batch1_S53_L002_R1_001.fastq.gz CL29118-batch1_S53_L003_R1_001.fastq.gz CL29118-batch1_S53_L004_R1_001.fastq.gz > CL29118_2.fastq.gz 
#cat CL29137-batch1_S52_L001_R1_001.fastq.gz CL29137-batch1_S52_L002_R1_001.fastq.gz CL29137-batch1_S52_L003_R1_001.fastq.gz CL29137-batch1_S52_L004_R1_001.fastq.gz > CL29137_2.fastq.gz 
#cat CL29161-batch1_S47_L001_R1_001.fastq.gz CL29161-batch1_S47_L002_R1_001.fastq.gz CL29161-batch1_S47_L003_R1_001.fastq.gz CL29161-batch1_S47_L004_R1_001.fastq.gz > CL29161_2.fastq.gz 
#cat CL29170-batch2_S45_L001_R1_001.fastq.gz CL29170-batch2_S45_L002_R1_001.fastq.gz CL29170-batch2_S45_L003_R1_001.fastq.gz CL29170-batch2_S45_L004_R1_001.fastq.gz > CL29170_2.fastq.gz 
#cat CL29180-batch2_S44_L001_R1_001.fastq.gz CL29180-batch2_S44_L002_R1_001.fastq.gz CL29180-batch2_S44_L003_R1_001.fastq.gz CL29180-batch2_S44_L004_R1_001.fastq.gz > CL29180_2.fastq.gz 
#cat CL29159-batch1_S48_L001_R1_001.fastq.gz CL29159-batch1_S48_L002_R1_001.fastq.gz CL29159-batch1_S48_L003_R1_001.fastq.gz CL29159-batch1_S48_L004_R1_001.fastq.gz > CL29159_2.fastq.gz 
#cat CL29168-batch2_S46_L001_R1_001.fastq.gz CL29168-batch2_S46_L002_R1_001.fastq.gz CL29168-batch2_S46_L003_R1_001.fastq.gz CL29168-batch2_S46_L004_R1_001.fastq.gz > CL29168_2.fastq.gz 
#cat CL29213-batch2_S42_L001_R1_001.fastq.gz CL29213-batch2_S42_L002_R1_001.fastq.gz CL29213-batch2_S42_L003_R1_001.fastq.gz CL29213-batch2_S42_L004_R1_001.fastq.gz > CL29213_2.fastq.gz 
#cat CL29231-batch2_S38_L001_R1_001.fastq.gz CL29231-batch2_S38_L002_R1_001.fastq.gz CL29231-batch2_S38_L003_R1_001.fastq.gz CL29231-batch2_S38_L004_R1_001.fastq.gz > CL29231_2.fastq.gz 
#cat CL29201-batch2_S43_L001_R1_001.fastq.gz CL29201-batch2_S43_L002_R1_001.fastq.gz CL29201-batch2_S43_L003_R1_001.fastq.gz CL29201-batch2_S43_L004_R1_001.fastq.gz > CL29201_2.fastq.gz 
#cat CL29219-batch2_S40_L001_R1_001.fastq.gz CL29219-batch2_S40_L002_R1_001.fastq.gz CL29219-batch2_S40_L003_R1_001.fastq.gz CL29219-batch2_S40_L004_R1_001.fastq.gz > CL29219_2.fastq.gz 
#cat CL29234-batch2_S37_L001_R1_001.fastq.gz CL29234-batch2_S37_L002_R1_001.fastq.gz CL29234-batch2_S37_L003_R1_001.fastq.gz CL29234-batch2_S37_L004_R1_001.fastq.gz > CL29234_2.fastq.gz 
#cat CL29246-batch2_S36_L001_R1_001.fastq.gz CL29246-batch2_S36_L002_R1_001.fastq.gz CL29246-batch2_S36_L003_R1_001.fastq.gz CL29246-batch2_S36_L004_R1_001.fastq.gz > CL29246_2.fastq.gz 
#cat CL29223-batch2_S39_L001_R1_001.fastq.gz CL29223-batch2_S39_L002_R1_001.fastq.gz CL29223-batch2_S39_L003_R1_001.fastq.gz CL29223-batch2_S39_L004_R1_001.fastq.gz > CL29223_2.fastq.gz 
#cat CL29218-batch2_S41_L001_R1_001.fastq.gz CL29218-batch2_S41_L002_R1_001.fastq.gz CL29218-batch2_S41_L003_R1_001.fastq.gz CL29218-batch2_S41_L004_R1_001.fastq.gz > CL29218_2.fastq.gz 
#cat CL29253-batch2_S35_L001_R1_001.fastq.gz CL29253-batch2_S35_L002_R1_001.fastq.gz CL29253-batch2_S35_L003_R1_001.fastq.gz CL29253-batch2_S35_L004_R1_001.fastq.gz > CL29253_2.fastq.gz 
#cat CL29254-batch2_S34_L001_R1_001.fastq.gz CL29254-batch2_S34_L002_R1_001.fastq.gz CL29254-batch2_S34_L003_R1_001.fastq.gz CL29254-batch2_S34_L004_R1_001.fastq.gz > CL29254_2.fastq.gz 
#cat CL29256-batch2_S33_L001_R1_001.fastq.gz CL29256-batch2_S33_L002_R1_001.fastq.gz CL29256-batch2_S33_L003_R1_001.fastq.gz CL29256-batch2_S33_L004_R1_001.fastq.gz > CL29256_2.fastq.gz 
#cat CL29257-batch2_S32_L001_R1_001.fastq.gz CL29257-batch2_S32_L002_R1_001.fastq.gz CL29257-batch2_S32_L003_R1_001.fastq.gz CL29257-batch2_S32_L004_R1_001.fastq.gz > CL29257_2.fastq.gz 
#cat CL29258-batch2_S31_L001_R1_001.fastq.gz CL29258-batch2_S31_L002_R1_001.fastq.gz CL29258-batch2_S31_L003_R1_001.fastq.gz CL29258-batch2_S31_L004_R1_001.fastq.gz > CL29258_2.fastq.gz 
#cat CL29264-batch2_S30_L001_R1_001.fastq.gz CL29264-batch2_S30_L002_R1_001.fastq.gz CL29264-batch2_S30_L003_R1_001.fastq.gz CL29264-batch2_S30_L004_R1_001.fastq.gz > CL29264_2.fastq.gz 
#cat CL29270-batch2_S29_L001_R1_001.fastq.gz CL29270-batch2_S29_L002_R1_001.fastq.gz CL29270-batch2_S29_L003_R1_001.fastq.gz CL29270-batch2_S29_L004_R1_001.fastq.gz > CL29270_2.fastq.gz 
#cat CL29271-batch3_S28_L001_R1_001.fastq.gz CL29271-batch3_S28_L002_R1_001.fastq.gz CL29271-batch3_S28_L003_R1_001.fastq.gz CL29271-batch3_S28_L004_R1_001.fastq.gz > CL29271_2.fastq.gz 
#cat CL29272-batch2_S27_L001_R1_001.fastq.gz CL29272-batch2_S27_L002_R1_001.fastq.gz CL29272-batch2_S27_L003_R1_001.fastq.gz CL29272-batch2_S27_L004_R1_001.fastq.gz > CL29272_2.fastq.gz 
#cat CL29274-batch3_S26_L001_R1_001.fastq.gz CL29274-batch3_S26_L002_R1_001.fastq.gz CL29274-batch3_S26_L003_R1_001.fastq.gz CL29274-batch3_S26_L004_R1_001.fastq.gz > CL29274_2.fastq.gz 
#cat CL29276-batch3_S25_L001_R1_001.fastq.gz CL29276-batch3_S25_L002_R1_001.fastq.gz CL29276-batch3_S25_L003_R1_001.fastq.gz CL29276-batch3_S25_L004_R1_001.fastq.gz > CL29276_2.fastq.gz 
#cat CL29289-batch3_S24_L001_R1_001.fastq.gz CL29289-batch3_S24_L002_R1_001.fastq.gz CL29289-batch3_S24_L003_R1_001.fastq.gz CL29289-batch3_S24_L004_R1_001.fastq.gz > CL29289_2.fastq.gz 
#cat CL29291-batch3_S23_L001_R1_001.fastq.gz CL29291-batch3_S23_L002_R1_001.fastq.gz CL29291-batch3_S23_L003_R1_001.fastq.gz CL29291-batch3_S23_L004_R1_001.fastq.gz > CL29291_2.fastq.gz 
#cat CL29297-batch3_S22_L001_R1_001.fastq.gz CL29297-batch3_S22_L002_R1_001.fastq.gz CL29297-batch3_S22_L003_R1_001.fastq.gz CL29297-batch3_S22_L004_R1_001.fastq.gz > CL29297_2.fastq.gz 
#cat CL29298-batch3_S21_L001_R1_001.fastq.gz CL29298-batch3_S21_L002_R1_001.fastq.gz CL29298-batch3_S21_L003_R1_001.fastq.gz CL29298-batch3_S21_L004_R1_001.fastq.gz > CL29298_2.fastq.gz 
#cat CL29303-batch3_S20_L001_R1_001.fastq.gz CL29303-batch3_S20_L002_R1_001.fastq.gz CL29303-batch3_S20_L003_R1_001.fastq.gz CL29303-batch3_S20_L004_R1_001.fastq.gz > CL29303_2.fastq.gz 
#cat CL29491-batch3_S19_L001_R1_001.fastq.gz CL29491-batch3_S19_L002_R1_001.fastq.gz CL29491-batch3_S19_L003_R1_001.fastq.gz CL29491-batch3_S19_L004_R1_001.fastq.gz > CL29491_2.fastq.gz 
#cat CL29498-batch3_S18_L001_R1_001.fastq.gz CL29498-batch3_S18_L002_R1_001.fastq.gz CL29498-batch3_S18_L003_R1_001.fastq.gz CL29498-batch3_S18_L004_R1_001.fastq.gz > CL29498_2.fastq.gz 
#cat CL29513-batch3_S17_L001_R1_001.fastq.gz CL29513-batch3_S17_L002_R1_001.fastq.gz CL29513-batch3_S17_L003_R1_001.fastq.gz CL29513-batch3_S17_L004_R1_001.fastq.gz > CL29513_2.fastq.gz 
#cat CL29519-batch3_S16_L001_R1_001.fastq.gz CL29519-batch3_S16_L002_R1_001.fastq.gz CL29519-batch3_S16_L003_R1_001.fastq.gz CL29519-batch3_S16_L004_R1_001.fastq.gz > CL29519_2.fastq.gz 
#cat CL29529-batch3_S15_L001_R1_001.fastq.gz CL29529-batch3_S15_L002_R1_001.fastq.gz CL29529-batch3_S15_L003_R1_001.fastq.gz CL29529-batch3_S15_L004_R1_001.fastq.gz > CL29529_2.fastq.gz 
#cat CL29539-batch3_S14_L001_R1_001.fastq.gz CL29539-batch3_S14_L002_R1_001.fastq.gz CL29539-batch3_S14_L003_R1_001.fastq.gz CL29539-batch3_S14_L004_R1_001.fastq.gz > CL29539_2.fastq.gz 
#cat CL29558-batch3_S13_L001_R1_001.fastq.gz CL29558-batch3_S13_L002_R1_001.fastq.gz CL29558-batch3_S13_L003_R1_001.fastq.gz CL29558-batch3_S13_L004_R1_001.fastq.gz > CL29558_2.fastq.gz 
#cat CL29663-batch3_S12_L001_R1_001.fastq.gz CL29663-batch3_S12_L002_R1_001.fastq.gz CL29663-batch3_S12_L003_R1_001.fastq.gz CL29663-batch3_S12_L004_R1_001.fastq.gz > CL29663_2.fastq.gz 
#cat CL29680-batch3_S11_L001_R1_001.fastq.gz CL29680-batch3_S11_L002_R1_001.fastq.gz CL29680-batch3_S11_L003_R1_001.fastq.gz CL29680-batch3_S11_L004_R1_001.fastq.gz > CL29680_2.fastq.gz 
#cat CL29687-batch3_S10_L001_R1_001.fastq.gz CL29687-batch3_S10_L002_R1_001.fastq.gz CL29687-batch3_S10_L003_R1_001.fastq.gz CL29687-batch3_S10_L004_R1_001.fastq.gz > CL29687_2.fastq.gz 
#cat CL29727-batch3_S9_L001_R1_001.fastq.gz CL29727-batch3_S9_L002_R1_001.fastq.gz CL29727-batch3_S9_L003_R1_001.fastq.gz CL29727-batch3_S9_L004_R1_001.fastq.gz > CL29727_2.fastq.gz 
#cat CL29770-batch3_S8_L001_R1_001.fastq.gz CL29770-batch3_S8_L002_R1_001.fastq.gz CL29770-batch3_S8_L003_R1_001.fastq.gz CL29770-batch3_S8_L004_R1_001.fastq.gz > CL29770_2.fastq.gz 
#cat CL29776-batch3_S7_L001_R1_001.fastq.gz CL29776-batch3_S7_L002_R1_001.fastq.gz CL29776-batch3_S7_L003_R1_001.fastq.gz CL29776-batch3_S7_L004_R1_001.fastq.gz > CL29776_2.fastq.gz 
#cat HS1-batch1_S6_L001_R1_001.fastq.gz HS1-batch1_S6_L002_R1_001.fastq.gz HS1-batch1_S6_L003_R1_001.fastq.gz HS1-batch1_S6_L004_R1_001.fastq.gz > HS1_2.fastq.gz 
#cat HS2-batch2_S5_L001_R1_001.fastq.gz HS2-batch2_S5_L002_R1_001.fastq.gz HS2-batch2_S5_L003_R1_001.fastq.gz HS2-batch2_S5_L004_R1_001.fastq.gz > HS2_2.fastq.gz 
#cat HS3-batch2_S4_L001_R1_001.fastq.gz HS3-batch2_S4_L002_R1_001.fastq.gz HS3-batch2_S4_L003_R1_001.fastq.gz HS3-batch2_S4_L004_R1_001.fastq.gz > HS3_2.fastq.gz 
#cat HS4-batch3_S3_L001_R1_001.fastq.gz HS4-batch3_S3_L002_R1_001.fastq.gz HS4-batch3_S3_L003_R1_001.fastq.gz HS4-batch3_S3_L004_R1_001.fastq.gz > HS4_2.fastq.gz 
#cat HS5-batch3_S2_L001_R1_001.fastq.gz HS5-batch3_S2_L002_R1_001.fastq.gz HS5-batch3_S2_L003_R1_001.fastq.gz HS5-batch3_S2_L004_R1_001.fastq.gz > HS5_2.fastq.gz 
#cat HS6-batch3_S1_L001_R1_001.fastq.gz HS6-batch3_S1_L002_R1_001.fastq.gz HS6-batch3_S1_L003_R1_001.fastq.gz HS6-batch3_S1_L004_R1_001.fastq.gz > HS6_2.fastq.gz #

##3
#cat CL29075-batch1_S57_L001_R1_001.fastq.gz CL29075-batch1_S57_L002_R1_001.fastq.gz CL29075-batch1_S57_L003_R1_001.fastq.gz CL29075-batch1_S57_L004_R1_001.fastq.gz > CL29075_3.fastq.gz 
#cat CL29107-batch1_S54_L001_R1_001.fastq.gz CL29107-batch1_S54_L002_R1_001.fastq.gz CL29107-batch1_S54_L003_R1_001.fastq.gz CL29107-batch1_S54_L004_R1_001.fastq.gz > CL29107_3.fastq.gz 
#cat CL29077-batch1_S56_L001_R1_001.fastq.gz CL29077-batch1_S56_L002_R1_001.fastq.gz CL29077-batch1_S56_L003_R1_001.fastq.gz CL29077-batch1_S56_L004_R1_001.fastq.gz > CL29077_3.fastq.gz 
#cat CL29151-batch1_S49_L001_R1_001.fastq.gz CL29151-batch1_S49_L002_R1_001.fastq.gz CL29151-batch1_S49_L003_R1_001.fastq.gz CL29151-batch1_S49_L004_R1_001.fastq.gz > CL29151_3.fastq.gz 
#cat CL29148-batch1_S51_L001_R1_001.fastq.gz CL29148-batch1_S51_L002_R1_001.fastq.gz CL29148-batch1_S51_L003_R1_001.fastq.gz CL29148-batch1_S51_L004_R1_001.fastq.gz > CL29148_3.fastq.gz 
#cat CL29150-batch1_S50_L001_R1_001.fastq.gz CL29150-batch1_S50_L002_R1_001.fastq.gz CL29150-batch1_S50_L003_R1_001.fastq.gz CL29150-batch1_S50_L004_R1_001.fastq.gz > CL29150_3.fastq.gz 
#cat CL29085-batch2_S55_L001_R1_001.fastq.gz CL29085-batch2_S55_L002_R1_001.fastq.gz CL29085-batch2_S55_L003_R1_001.fastq.gz CL29085-batch2_S55_L004_R1_001.fastq.gz > CL29085_3.fastq.gz 
#cat CL29118-batch1_S53_L001_R1_001.fastq.gz CL29118-batch1_S53_L002_R1_001.fastq.gz CL29118-batch1_S53_L003_R1_001.fastq.gz CL29118-batch1_S53_L004_R1_001.fastq.gz > CL29118_3.fastq.gz 
#cat CL29137-batch1_S52_L001_R1_001.fastq.gz CL29137-batch1_S52_L002_R1_001.fastq.gz CL29137-batch1_S52_L003_R1_001.fastq.gz CL29137-batch1_S52_L004_R1_001.fastq.gz > CL29137_3.fastq.gz 
#cat CL29161-batch1_S47_L001_R1_001.fastq.gz CL29161-batch1_S47_L002_R1_001.fastq.gz CL29161-batch1_S47_L003_R1_001.fastq.gz CL29161-batch1_S47_L004_R1_001.fastq.gz > CL29161_3.fastq.gz 
#cat CL29170-batch2_S45_L001_R1_001.fastq.gz CL29170-batch2_S45_L002_R1_001.fastq.gz CL29170-batch2_S45_L003_R1_001.fastq.gz CL29170-batch2_S45_L004_R1_001.fastq.gz > CL29170_3.fastq.gz 
#cat CL29180-batch2_S44_L001_R1_001.fastq.gz CL29180-batch2_S44_L002_R1_001.fastq.gz CL29180-batch2_S44_L003_R1_001.fastq.gz CL29180-batch2_S44_L004_R1_001.fastq.gz > CL29180_3.fastq.gz 
#cat CL29159-batch1_S48_L001_R1_001.fastq.gz CL29159-batch1_S48_L002_R1_001.fastq.gz CL29159-batch1_S48_L003_R1_001.fastq.gz CL29159-batch1_S48_L004_R1_001.fastq.gz > CL29159_3.fastq.gz 
#cat CL29168-batch2_S46_L001_R1_001.fastq.gz CL29168-batch2_S46_L002_R1_001.fastq.gz CL29168-batch2_S46_L003_R1_001.fastq.gz CL29168-batch2_S46_L004_R1_001.fastq.gz > CL29168_3.fastq.gz 
#cat CL29213-batch2_S42_L001_R1_001.fastq.gz CL29213-batch2_S42_L002_R1_001.fastq.gz CL29213-batch2_S42_L003_R1_001.fastq.gz CL29213-batch2_S42_L004_R1_001.fastq.gz > CL29213_3.fastq.gz 
#cat CL29231-batch2_S38_L001_R1_001.fastq.gz CL29231-batch2_S38_L002_R1_001.fastq.gz CL29231-batch2_S38_L003_R1_001.fastq.gz CL29231-batch2_S38_L004_R1_001.fastq.gz > CL29231_3.fastq.gz 
#cat CL29201-batch2_S43_L001_R1_001.fastq.gz CL29201-batch2_S43_L002_R1_001.fastq.gz CL29201-batch2_S43_L003_R1_001.fastq.gz CL29201-batch2_S43_L004_R1_001.fastq.gz > CL29201_3.fastq.gz 
#cat CL29219-batch2_S40_L001_R1_001.fastq.gz CL29219-batch2_S40_L002_R1_001.fastq.gz CL29219-batch2_S40_L003_R1_001.fastq.gz CL29219-batch2_S40_L004_R1_001.fastq.gz > CL29219_3.fastq.gz 
#cat CL29234-batch2_S37_L001_R1_001.fastq.gz CL29234-batch2_S37_L002_R1_001.fastq.gz CL29234-batch2_S37_L003_R1_001.fastq.gz CL29234-batch2_S37_L004_R1_001.fastq.gz > CL29234_3.fastq.gz 
#cat CL29246-batch2_S36_L001_R1_001.fastq.gz CL29246-batch2_S36_L002_R1_001.fastq.gz CL29246-batch2_S36_L003_R1_001.fastq.gz CL29246-batch2_S36_L004_R1_001.fastq.gz > CL29246_3.fastq.gz 
#cat CL29223-batch2_S39_L001_R1_001.fastq.gz CL29223-batch2_S39_L002_R1_001.fastq.gz CL29223-batch2_S39_L003_R1_001.fastq.gz CL29223-batch2_S39_L004_R1_001.fastq.gz > CL29223_3.fastq.gz 
#cat CL29218-batch2_S41_L001_R1_001.fastq.gz CL29218-batch2_S41_L002_R1_001.fastq.gz CL29218-batch2_S41_L003_R1_001.fastq.gz CL29218-batch2_S41_L004_R1_001.fastq.gz > CL29218_3.fastq.gz 
#cat CL29253-batch2_S35_L001_R1_001.fastq.gz CL29253-batch2_S35_L002_R1_001.fastq.gz CL29253-batch2_S35_L003_R1_001.fastq.gz CL29253-batch2_S35_L004_R1_001.fastq.gz > CL29253_3.fastq.gz 
#cat CL29254-batch2_S34_L001_R1_001.fastq.gz CL29254-batch2_S34_L002_R1_001.fastq.gz CL29254-batch2_S34_L003_R1_001.fastq.gz CL29254-batch2_S34_L004_R1_001.fastq.gz > CL29254_3.fastq.gz 
#cat CL29256-batch2_S33_L001_R1_001.fastq.gz CL29256-batch2_S33_L002_R1_001.fastq.gz CL29256-batch2_S33_L003_R1_001.fastq.gz CL29256-batch2_S33_L004_R1_001.fastq.gz > CL29256_3.fastq.gz 
#cat CL29257-batch2_S32_L001_R1_001.fastq.gz CL29257-batch2_S32_L002_R1_001.fastq.gz CL29257-batch2_S32_L003_R1_001.fastq.gz CL29257-batch2_S32_L004_R1_001.fastq.gz > CL29257_3.fastq.gz 
#cat CL29258-batch2_S31_L001_R1_001.fastq.gz CL29258-batch2_S31_L002_R1_001.fastq.gz CL29258-batch2_S31_L003_R1_001.fastq.gz CL29258-batch2_S31_L004_R1_001.fastq.gz > CL29258_3.fastq.gz 
#cat CL29264-batch2_S30_L001_R1_001.fastq.gz CL29264-batch2_S30_L002_R1_001.fastq.gz CL29264-batch2_S30_L003_R1_001.fastq.gz CL29264-batch2_S30_L004_R1_001.fastq.gz > CL29264_3.fastq.gz 
#cat CL29270-batch2_S29_L001_R1_001.fastq.gz CL29270-batch2_S29_L002_R1_001.fastq.gz CL29270-batch2_S29_L003_R1_001.fastq.gz CL29270-batch2_S29_L004_R1_001.fastq.gz > CL29270_3.fastq.gz 
#cat CL29271-batch3_S28_L001_R1_001.fastq.gz CL29271-batch3_S28_L002_R1_001.fastq.gz CL29271-batch3_S28_L003_R1_001.fastq.gz CL29271-batch3_S28_L004_R1_001.fastq.gz > CL29271_3.fastq.gz 
#cat CL29272-batch2_S27_L001_R1_001.fastq.gz CL29272-batch2_S27_L002_R1_001.fastq.gz CL29272-batch2_S27_L003_R1_001.fastq.gz CL29272-batch2_S27_L004_R1_001.fastq.gz > CL29272_3.fastq.gz 
#cat CL29274-batch3_S26_L001_R1_001.fastq.gz CL29274-batch3_S26_L002_R1_001.fastq.gz CL29274-batch3_S26_L003_R1_001.fastq.gz CL29274-batch3_S26_L004_R1_001.fastq.gz > CL29274_3.fastq.gz 
#cat CL29276-batch3_S25_L001_R1_001.fastq.gz CL29276-batch3_S25_L002_R1_001.fastq.gz CL29276-batch3_S25_L003_R1_001.fastq.gz CL29276-batch3_S25_L004_R1_001.fastq.gz > CL29276_3.fastq.gz 
#cat CL29289-batch3_S24_L001_R1_001.fastq.gz CL29289-batch3_S24_L002_R1_001.fastq.gz CL29289-batch3_S24_L003_R1_001.fastq.gz CL29289-batch3_S24_L004_R1_001.fastq.gz > CL29289_3.fastq.gz 
#cat CL29291-batch3_S23_L001_R1_001.fastq.gz CL29291-batch3_S23_L002_R1_001.fastq.gz CL29291-batch3_S23_L003_R1_001.fastq.gz CL29291-batch3_S23_L004_R1_001.fastq.gz > CL29291_3.fastq.gz 
#cat CL29297-batch3_S22_L001_R1_001.fastq.gz CL29297-batch3_S22_L002_R1_001.fastq.gz CL29297-batch3_S22_L003_R1_001.fastq.gz CL29297-batch3_S22_L004_R1_001.fastq.gz > CL29297_3.fastq.gz 
#cat CL29298-batch3_S21_L001_R1_001.fastq.gz CL29298-batch3_S21_L002_R1_001.fastq.gz CL29298-batch3_S21_L003_R1_001.fastq.gz CL29298-batch3_S21_L004_R1_001.fastq.gz > CL29298_3.fastq.gz 
#cat CL29303-batch3_S20_L001_R1_001.fastq.gz CL29303-batch3_S20_L002_R1_001.fastq.gz CL29303-batch3_S20_L003_R1_001.fastq.gz CL29303-batch3_S20_L004_R1_001.fastq.gz > CL29303_3.fastq.gz 
#cat CL29491-batch3_S19_L001_R1_001.fastq.gz CL29491-batch3_S19_L002_R1_001.fastq.gz CL29491-batch3_S19_L003_R1_001.fastq.gz CL29491-batch3_S19_L004_R1_001.fastq.gz > CL29491_3.fastq.gz 
#cat CL29498-batch3_S18_L001_R1_001.fastq.gz CL29498-batch3_S18_L002_R1_001.fastq.gz CL29498-batch3_S18_L003_R1_001.fastq.gz CL29498-batch3_S18_L004_R1_001.fastq.gz > CL29498_3.fastq.gz 
#cat CL29513-batch3_S17_L001_R1_001.fastq.gz CL29513-batch3_S17_L002_R1_001.fastq.gz CL29513-batch3_S17_L003_R1_001.fastq.gz CL29513-batch3_S17_L004_R1_001.fastq.gz > CL29513_3.fastq.gz 
#cat CL29519-batch3_S16_L001_R1_001.fastq.gz CL29519-batch3_S16_L002_R1_001.fastq.gz CL29519-batch3_S16_L003_R1_001.fastq.gz CL29519-batch3_S16_L004_R1_001.fastq.gz > CL29519_3.fastq.gz 
#cat CL29529-batch3_S15_L001_R1_001.fastq.gz CL29529-batch3_S15_L002_R1_001.fastq.gz CL29529-batch3_S15_L003_R1_001.fastq.gz CL29529-batch3_S15_L004_R1_001.fastq.gz > CL29529_3.fastq.gz 
#cat CL29539-batch3_S14_L001_R1_001.fastq.gz CL29539-batch3_S14_L002_R1_001.fastq.gz CL29539-batch3_S14_L003_R1_001.fastq.gz CL29539-batch3_S14_L004_R1_001.fastq.gz > CL29539_3.fastq.gz 
#cat CL29558-batch3_S13_L001_R1_001.fastq.gz CL29558-batch3_S13_L002_R1_001.fastq.gz CL29558-batch3_S13_L003_R1_001.fastq.gz CL29558-batch3_S13_L004_R1_001.fastq.gz > CL29558_3.fastq.gz 
#cat CL29663-batch3_S12_L001_R1_001.fastq.gz CL29663-batch3_S12_L002_R1_001.fastq.gz CL29663-batch3_S12_L003_R1_001.fastq.gz CL29663-batch3_S12_L004_R1_001.fastq.gz > CL29663_3.fastq.gz 
#cat CL29680-batch3_S11_L001_R1_001.fastq.gz CL29680-batch3_S11_L002_R1_001.fastq.gz CL29680-batch3_S11_L003_R1_001.fastq.gz CL29680-batch3_S11_L004_R1_001.fastq.gz > CL29680_3.fastq.gz 
#cat CL29687-batch3_S10_L001_R1_001.fastq.gz CL29687-batch3_S10_L002_R1_001.fastq.gz CL29687-batch3_S10_L003_R1_001.fastq.gz CL29687-batch3_S10_L004_R1_001.fastq.gz > CL29687_3.fastq.gz 
#cat CL29727-batch3_S9_L001_R1_001.fastq.gz CL29727-batch3_S9_L002_R1_001.fastq.gz CL29727-batch3_S9_L003_R1_001.fastq.gz CL29727-batch3_S9_L004_R1_001.fastq.gz > CL29727_3.fastq.gz 
#cat CL29770-batch3_S8_L001_R1_001.fastq.gz CL29770-batch3_S8_L002_R1_001.fastq.gz CL29770-batch3_S8_L003_R1_001.fastq.gz CL29770-batch3_S8_L004_R1_001.fastq.gz > CL29770_3.fastq.gz 
#cat CL29776-batch3_S7_L001_R1_001.fastq.gz CL29776-batch3_S7_L002_R1_001.fastq.gz CL29776-batch3_S7_L003_R1_001.fastq.gz CL29776-batch3_S7_L004_R1_001.fastq.gz > CL29776_3.fastq.gz 
#cat HS1-batch1_S6_L001_R1_001.fastq.gz HS1-batch1_S6_L002_R1_001.fastq.gz HS1-batch1_S6_L003_R1_001.fastq.gz HS1-batch1_S6_L004_R1_001.fastq.gz > HS1_3.fastq.gz 
#cat HS2-batch2_S5_L001_R1_001.fastq.gz HS2-batch2_S5_L002_R1_001.fastq.gz HS2-batch2_S5_L003_R1_001.fastq.gz HS2-batch2_S5_L004_R1_001.fastq.gz > HS2_3.fastq.gz 
#cat HS3-batch2_S4_L001_R1_001.fastq.gz HS3-batch2_S4_L002_R1_001.fastq.gz HS3-batch2_S4_L003_R1_001.fastq.gz HS3-batch2_S4_L004_R1_001.fastq.gz > HS3_3.fastq.gz 
#cat HS4-batch3_S3_L001_R1_001.fastq.gz HS4-batch3_S3_L002_R1_001.fastq.gz HS4-batch3_S3_L003_R1_001.fastq.gz HS4-batch3_S3_L004_R1_001.fastq.gz > HS4_3.fastq.gz 
#cat HS5-batch3_S2_L001_R1_001.fastq.gz HS5-batch3_S2_L002_R1_001.fastq.gz HS5-batch3_S2_L003_R1_001.fastq.gz HS5-batch3_S2_L004_R1_001.fastq.gz > HS5_3.fastq.gz 
#cat HS6-batch3_S1_L001_R1_001.fastq.gz HS6-batch3_S1_L002_R1_001.fastq.gz HS6-batch3_S1_L003_R1_001.fastq.gz HS6-batch3_S1_L004_R1_001.fastq.gz > HS6_3.fastq.gz #

##4
#cat CL29075-batch1_S57_L001_R1_001.fastq.gz CL29075-batch1_S57_L002_R1_001.fastq.gz CL29075-batch1_S57_L003_R1_001.fastq.gz CL29075-batch1_S57_L004_R1_001.fastq.gz > CL29075_4.fastq.gz 
#cat CL29107-batch1_S54_L001_R1_001.fastq.gz CL29107-batch1_S54_L002_R1_001.fastq.gz CL29107-batch1_S54_L003_R1_001.fastq.gz CL29107-batch1_S54_L004_R1_001.fastq.gz > CL29107_4.fastq.gz 
#cat CL29077-batch1_S56_L001_R1_001.fastq.gz CL29077-batch1_S56_L002_R1_001.fastq.gz CL29077-batch1_S56_L003_R1_001.fastq.gz CL29077-batch1_S56_L004_R1_001.fastq.gz > CL29077_4.fastq.gz 
#cat CL29151-batch1_S49_L001_R1_001.fastq.gz CL29151-batch1_S49_L002_R1_001.fastq.gz CL29151-batch1_S49_L003_R1_001.fastq.gz CL29151-batch1_S49_L004_R1_001.fastq.gz > CL29151_4.fastq.gz 
#cat CL29148-batch1_S51_L001_R1_001.fastq.gz CL29148-batch1_S51_L002_R1_001.fastq.gz CL29148-batch1_S51_L003_R1_001.fastq.gz CL29148-batch1_S51_L004_R1_001.fastq.gz > CL29148_4.fastq.gz 
#cat CL29150-batch1_S50_L001_R1_001.fastq.gz CL29150-batch1_S50_L002_R1_001.fastq.gz CL29150-batch1_S50_L003_R1_001.fastq.gz CL29150-batch1_S50_L004_R1_001.fastq.gz > CL29150_4.fastq.gz 
#cat CL29085-batch2_S55_L001_R1_001.fastq.gz CL29085-batch2_S55_L002_R1_001.fastq.gz CL29085-batch2_S55_L003_R1_001.fastq.gz CL29085-batch2_S55_L004_R1_001.fastq.gz > CL29085_4.fastq.gz 
#cat CL29118-batch1_S53_L001_R1_001.fastq.gz CL29118-batch1_S53_L002_R1_001.fastq.gz CL29118-batch1_S53_L003_R1_001.fastq.gz CL29118-batch1_S53_L004_R1_001.fastq.gz > CL29118_4.fastq.gz 
#cat CL29137-batch1_S52_L001_R1_001.fastq.gz CL29137-batch1_S52_L002_R1_001.fastq.gz CL29137-batch1_S52_L003_R1_001.fastq.gz CL29137-batch1_S52_L004_R1_001.fastq.gz > CL29137_4.fastq.gz 
#cat CL29161-batch1_S47_L001_R1_001.fastq.gz CL29161-batch1_S47_L002_R1_001.fastq.gz CL29161-batch1_S47_L003_R1_001.fastq.gz CL29161-batch1_S47_L004_R1_001.fastq.gz > CL29161_4.fastq.gz 
#cat CL29170-batch2_S45_L001_R1_001.fastq.gz CL29170-batch2_S45_L002_R1_001.fastq.gz CL29170-batch2_S45_L003_R1_001.fastq.gz CL29170-batch2_S45_L004_R1_001.fastq.gz > CL29170_4.fastq.gz 
#cat CL29180-batch2_S44_L001_R1_001.fastq.gz CL29180-batch2_S44_L002_R1_001.fastq.gz CL29180-batch2_S44_L003_R1_001.fastq.gz CL29180-batch2_S44_L004_R1_001.fastq.gz > CL29180_4.fastq.gz 
#cat CL29159-batch1_S48_L001_R1_001.fastq.gz CL29159-batch1_S48_L002_R1_001.fastq.gz CL29159-batch1_S48_L003_R1_001.fastq.gz CL29159-batch1_S48_L004_R1_001.fastq.gz > CL29159_4.fastq.gz 
#cat CL29168-batch2_S46_L001_R1_001.fastq.gz CL29168-batch2_S46_L002_R1_001.fastq.gz CL29168-batch2_S46_L003_R1_001.fastq.gz CL29168-batch2_S46_L004_R1_001.fastq.gz > CL29168_4.fastq.gz 
#cat CL29213-batch2_S42_L001_R1_001.fastq.gz CL29213-batch2_S42_L002_R1_001.fastq.gz CL29213-batch2_S42_L003_R1_001.fastq.gz CL29213-batch2_S42_L004_R1_001.fastq.gz > CL29213_4.fastq.gz 
#cat CL29231-batch2_S38_L001_R1_001.fastq.gz CL29231-batch2_S38_L002_R1_001.fastq.gz CL29231-batch2_S38_L003_R1_001.fastq.gz CL29231-batch2_S38_L004_R1_001.fastq.gz > CL29231_4.fastq.gz 
#cat CL29201-batch2_S43_L001_R1_001.fastq.gz CL29201-batch2_S43_L002_R1_001.fastq.gz CL29201-batch2_S43_L003_R1_001.fastq.gz CL29201-batch2_S43_L004_R1_001.fastq.gz > CL29201_4.fastq.gz 
#cat CL29219-batch2_S40_L001_R1_001.fastq.gz CL29219-batch2_S40_L002_R1_001.fastq.gz CL29219-batch2_S40_L003_R1_001.fastq.gz CL29219-batch2_S40_L004_R1_001.fastq.gz > CL29219_4.fastq.gz 
#cat CL29234-batch2_S37_L001_R1_001.fastq.gz CL29234-batch2_S37_L002_R1_001.fastq.gz CL29234-batch2_S37_L003_R1_001.fastq.gz CL29234-batch2_S37_L004_R1_001.fastq.gz > CL29234_4.fastq.gz 
#cat CL29246-batch2_S36_L001_R1_001.fastq.gz CL29246-batch2_S36_L002_R1_001.fastq.gz CL29246-batch2_S36_L003_R1_001.fastq.gz CL29246-batch2_S36_L004_R1_001.fastq.gz > CL29246_4.fastq.gz 
#cat CL29223-batch2_S39_L001_R1_001.fastq.gz CL29223-batch2_S39_L002_R1_001.fastq.gz CL29223-batch2_S39_L003_R1_001.fastq.gz CL29223-batch2_S39_L004_R1_001.fastq.gz > CL29223_4.fastq.gz 
#cat CL29218-batch2_S41_L001_R1_001.fastq.gz CL29218-batch2_S41_L002_R1_001.fastq.gz CL29218-batch2_S41_L003_R1_001.fastq.gz CL29218-batch2_S41_L004_R1_001.fastq.gz > CL29218_4.fastq.gz 
#cat CL29253-batch2_S35_L001_R1_001.fastq.gz CL29253-batch2_S35_L002_R1_001.fastq.gz CL29253-batch2_S35_L003_R1_001.fastq.gz CL29253-batch2_S35_L004_R1_001.fastq.gz > CL29253_4.fastq.gz 
#cat CL29254-batch2_S34_L001_R1_001.fastq.gz CL29254-batch2_S34_L002_R1_001.fastq.gz CL29254-batch2_S34_L003_R1_001.fastq.gz CL29254-batch2_S34_L004_R1_001.fastq.gz > CL29254_4.fastq.gz 
#cat CL29256-batch2_S33_L001_R1_001.fastq.gz CL29256-batch2_S33_L002_R1_001.fastq.gz CL29256-batch2_S33_L003_R1_001.fastq.gz CL29256-batch2_S33_L004_R1_001.fastq.gz > CL29256_4.fastq.gz 
#cat CL29257-batch2_S32_L001_R1_001.fastq.gz CL29257-batch2_S32_L002_R1_001.fastq.gz CL29257-batch2_S32_L003_R1_001.fastq.gz CL29257-batch2_S32_L004_R1_001.fastq.gz > CL29257_4.fastq.gz 
#cat CL29258-batch2_S31_L001_R1_001.fastq.gz CL29258-batch2_S31_L002_R1_001.fastq.gz CL29258-batch2_S31_L003_R1_001.fastq.gz CL29258-batch2_S31_L004_R1_001.fastq.gz > CL29258_4.fastq.gz 
#cat CL29264-batch2_S30_L001_R1_001.fastq.gz CL29264-batch2_S30_L002_R1_001.fastq.gz CL29264-batch2_S30_L003_R1_001.fastq.gz CL29264-batch2_S30_L004_R1_001.fastq.gz > CL29264_4.fastq.gz 
#cat CL29270-batch2_S29_L001_R1_001.fastq.gz CL29270-batch2_S29_L002_R1_001.fastq.gz CL29270-batch2_S29_L003_R1_001.fastq.gz CL29270-batch2_S29_L004_R1_001.fastq.gz > CL29270_4.fastq.gz 
#cat CL29271-batch3_S28_L001_R1_001.fastq.gz CL29271-batch3_S28_L002_R1_001.fastq.gz CL29271-batch3_S28_L003_R1_001.fastq.gz CL29271-batch3_S28_L004_R1_001.fastq.gz > CL29271_4.fastq.gz 
#cat CL29272-batch2_S27_L001_R1_001.fastq.gz CL29272-batch2_S27_L002_R1_001.fastq.gz CL29272-batch2_S27_L003_R1_001.fastq.gz CL29272-batch2_S27_L004_R1_001.fastq.gz > CL29272_4.fastq.gz 
#cat CL29274-batch3_S26_L001_R1_001.fastq.gz CL29274-batch3_S26_L002_R1_001.fastq.gz CL29274-batch3_S26_L003_R1_001.fastq.gz CL29274-batch3_S26_L004_R1_001.fastq.gz > CL29274_4.fastq.gz 
#cat CL29276-batch3_S25_L001_R1_001.fastq.gz CL29276-batch3_S25_L002_R1_001.fastq.gz CL29276-batch3_S25_L003_R1_001.fastq.gz CL29276-batch3_S25_L004_R1_001.fastq.gz > CL29276_4.fastq.gz 
#cat CL29289-batch3_S24_L001_R1_001.fastq.gz CL29289-batch3_S24_L002_R1_001.fastq.gz CL29289-batch3_S24_L003_R1_001.fastq.gz CL29289-batch3_S24_L004_R1_001.fastq.gz > CL29289_4.fastq.gz 
#cat CL29291-batch3_S23_L001_R1_001.fastq.gz CL29291-batch3_S23_L002_R1_001.fastq.gz CL29291-batch3_S23_L003_R1_001.fastq.gz CL29291-batch3_S23_L004_R1_001.fastq.gz > CL29291_4.fastq.gz 
#cat CL29297-batch3_S22_L001_R1_001.fastq.gz CL29297-batch3_S22_L002_R1_001.fastq.gz CL29297-batch3_S22_L003_R1_001.fastq.gz CL29297-batch3_S22_L004_R1_001.fastq.gz > CL29297_4.fastq.gz 
#cat CL29298-batch3_S21_L001_R1_001.fastq.gz CL29298-batch3_S21_L002_R1_001.fastq.gz CL29298-batch3_S21_L003_R1_001.fastq.gz CL29298-batch3_S21_L004_R1_001.fastq.gz > CL29298_4.fastq.gz 
#cat CL29303-batch3_S20_L001_R1_001.fastq.gz CL29303-batch3_S20_L002_R1_001.fastq.gz CL29303-batch3_S20_L003_R1_001.fastq.gz CL29303-batch3_S20_L004_R1_001.fastq.gz > CL29303_4.fastq.gz 
#cat CL29491-batch3_S19_L001_R1_001.fastq.gz CL29491-batch3_S19_L002_R1_001.fastq.gz CL29491-batch3_S19_L003_R1_001.fastq.gz CL29491-batch3_S19_L004_R1_001.fastq.gz > CL29491_4.fastq.gz 
#cat CL29498-batch3_S18_L001_R1_001.fastq.gz CL29498-batch3_S18_L002_R1_001.fastq.gz CL29498-batch3_S18_L003_R1_001.fastq.gz CL29498-batch3_S18_L004_R1_001.fastq.gz > CL29498_4.fastq.gz 
#cat CL29513-batch3_S17_L001_R1_001.fastq.gz CL29513-batch3_S17_L002_R1_001.fastq.gz CL29513-batch3_S17_L003_R1_001.fastq.gz CL29513-batch3_S17_L004_R1_001.fastq.gz > CL29513_4.fastq.gz 
#cat CL29519-batch3_S16_L001_R1_001.fastq.gz CL29519-batch3_S16_L002_R1_001.fastq.gz CL29519-batch3_S16_L003_R1_001.fastq.gz CL29519-batch3_S16_L004_R1_001.fastq.gz > CL29519_4.fastq.gz 
#cat CL29529-batch3_S15_L001_R1_001.fastq.gz CL29529-batch3_S15_L002_R1_001.fastq.gz CL29529-batch3_S15_L003_R1_001.fastq.gz CL29529-batch3_S15_L004_R1_001.fastq.gz > CL29529_4.fastq.gz 
#cat CL29539-batch3_S14_L001_R1_001.fastq.gz CL29539-batch3_S14_L002_R1_001.fastq.gz CL29539-batch3_S14_L003_R1_001.fastq.gz CL29539-batch3_S14_L004_R1_001.fastq.gz > CL29539_4.fastq.gz 
#cat CL29558-batch3_S13_L001_R1_001.fastq.gz CL29558-batch3_S13_L002_R1_001.fastq.gz CL29558-batch3_S13_L003_R1_001.fastq.gz CL29558-batch3_S13_L004_R1_001.fastq.gz > CL29558_4.fastq.gz 
#cat CL29663-batch3_S12_L001_R1_001.fastq.gz CL29663-batch3_S12_L002_R1_001.fastq.gz CL29663-batch3_S12_L003_R1_001.fastq.gz CL29663-batch3_S12_L004_R1_001.fastq.gz > CL29663_4.fastq.gz 
#cat CL29680-batch3_S11_L001_R1_001.fastq.gz CL29680-batch3_S11_L002_R1_001.fastq.gz CL29680-batch3_S11_L003_R1_001.fastq.gz CL29680-batch3_S11_L004_R1_001.fastq.gz > CL29680_4.fastq.gz 
#cat CL29687-batch3_S10_L001_R1_001.fastq.gz CL29687-batch3_S10_L002_R1_001.fastq.gz CL29687-batch3_S10_L003_R1_001.fastq.gz CL29687-batch3_S10_L004_R1_001.fastq.gz > CL29687_4.fastq.gz 
#cat CL29727-batch3_S9_L001_R1_001.fastq.gz CL29727-batch3_S9_L002_R1_001.fastq.gz CL29727-batch3_S9_L003_R1_001.fastq.gz CL29727-batch3_S9_L004_R1_001.fastq.gz > CL29727_4.fastq.gz 
#cat CL29770-batch3_S8_L001_R1_001.fastq.gz CL29770-batch3_S8_L002_R1_001.fastq.gz CL29770-batch3_S8_L003_R1_001.fastq.gz CL29770-batch3_S8_L004_R1_001.fastq.gz > CL29770_4.fastq.gz 
#cat CL29776-batch3_S7_L001_R1_001.fastq.gz CL29776-batch3_S7_L002_R1_001.fastq.gz CL29776-batch3_S7_L003_R1_001.fastq.gz CL29776-batch3_S7_L004_R1_001.fastq.gz > CL29776_4.fastq.gz 
#cat HS1-batch1_S6_L001_R1_001.fastq.gz HS1-batch1_S6_L002_R1_001.fastq.gz HS1-batch1_S6_L003_R1_001.fastq.gz HS1-batch1_S6_L004_R1_001.fastq.gz > HS1_4.fastq.gz 
#cat HS2-batch2_S5_L001_R1_001.fastq.gz HS2-batch2_S5_L002_R1_001.fastq.gz HS2-batch2_S5_L003_R1_001.fastq.gz HS2-batch2_S5_L004_R1_001.fastq.gz > HS2_4.fastq.gz 
#cat HS3-batch2_S4_L001_R1_001.fastq.gz HS3-batch2_S4_L002_R1_001.fastq.gz HS3-batch2_S4_L003_R1_001.fastq.gz HS3-batch2_S4_L004_R1_001.fastq.gz > HS3_4.fastq.gz 
#cat HS4-batch3_S3_L001_R1_001.fastq.gz HS4-batch3_S3_L002_R1_001.fastq.gz HS4-batch3_S3_L003_R1_001.fastq.gz HS4-batch3_S3_L004_R1_001.fastq.gz > HS4_4.fastq.gz 
#cat HS5-batch3_S2_L001_R1_001.fastq.gz HS5-batch3_S2_L002_R1_001.fastq.gz HS5-batch3_S2_L003_R1_001.fastq.gz HS5-batch3_S2_L004_R1_001.fastq.gz > HS5_4.fastq.gz 
#cat HS6-batch3_S1_L001_R1_001.fastq.gz HS6-batch3_S1_L002_R1_001.fastq.gz HS6-batch3_S1_L003_R1_001.fastq.gz HS6-batch3_S1_L004_R1_001.fastq.gz > HS6_4.fastq.gz #

##5
#cat CL29075-batch1_S57_L001_R1_001.fastq.gz CL29075-batch1_S57_L002_R1_001.fastq.gz CL29075-batch1_S57_L003_R1_001.fastq.gz CL29075-batch1_S57_L004_R1_001.fastq.gz > CL29075_5.fastq.gz 
#cat CL29107-batch1_S54_L001_R1_001.fastq.gz CL29107-batch1_S54_L002_R1_001.fastq.gz CL29107-batch1_S54_L003_R1_001.fastq.gz CL29107-batch1_S54_L004_R1_001.fastq.gz > CL29107_5.fastq.gz 
#cat CL29077-batch1_S56_L001_R1_001.fastq.gz CL29077-batch1_S56_L002_R1_001.fastq.gz CL29077-batch1_S56_L003_R1_001.fastq.gz CL29077-batch1_S56_L004_R1_001.fastq.gz > CL29077_5.fastq.gz 
#cat CL29151-batch1_S49_L001_R1_001.fastq.gz CL29151-batch1_S49_L002_R1_001.fastq.gz CL29151-batch1_S49_L003_R1_001.fastq.gz CL29151-batch1_S49_L004_R1_001.fastq.gz > CL29151_5.fastq.gz 
#cat CL29148-batch1_S51_L001_R1_001.fastq.gz CL29148-batch1_S51_L002_R1_001.fastq.gz CL29148-batch1_S51_L003_R1_001.fastq.gz CL29148-batch1_S51_L004_R1_001.fastq.gz > CL29148_5.fastq.gz 
#cat CL29150-batch1_S50_L001_R1_001.fastq.gz CL29150-batch1_S50_L002_R1_001.fastq.gz CL29150-batch1_S50_L003_R1_001.fastq.gz CL29150-batch1_S50_L004_R1_001.fastq.gz > CL29150_5.fastq.gz 
#cat CL29085-batch2_S55_L001_R1_001.fastq.gz CL29085-batch2_S55_L002_R1_001.fastq.gz CL29085-batch2_S55_L003_R1_001.fastq.gz CL29085-batch2_S55_L004_R1_001.fastq.gz > CL29085_5.fastq.gz 
#cat CL29118-batch1_S53_L001_R1_001.fastq.gz CL29118-batch1_S53_L002_R1_001.fastq.gz CL29118-batch1_S53_L003_R1_001.fastq.gz CL29118-batch1_S53_L004_R1_001.fastq.gz > CL29118_5.fastq.gz 
#cat CL29137-batch1_S52_L001_R1_001.fastq.gz CL29137-batch1_S52_L002_R1_001.fastq.gz CL29137-batch1_S52_L003_R1_001.fastq.gz CL29137-batch1_S52_L004_R1_001.fastq.gz > CL29137_5.fastq.gz 
#cat CL29161-batch1_S47_L001_R1_001.fastq.gz CL29161-batch1_S47_L002_R1_001.fastq.gz CL29161-batch1_S47_L003_R1_001.fastq.gz CL29161-batch1_S47_L004_R1_001.fastq.gz > CL29161_5.fastq.gz 
#cat CL29170-batch2_S45_L001_R1_001.fastq.gz CL29170-batch2_S45_L002_R1_001.fastq.gz CL29170-batch2_S45_L003_R1_001.fastq.gz CL29170-batch2_S45_L004_R1_001.fastq.gz > CL29170_5.fastq.gz 
#cat CL29180-batch2_S44_L001_R1_001.fastq.gz CL29180-batch2_S44_L002_R1_001.fastq.gz CL29180-batch2_S44_L003_R1_001.fastq.gz CL29180-batch2_S44_L004_R1_001.fastq.gz > CL29180_5.fastq.gz 
#cat CL29159-batch1_S48_L001_R1_001.fastq.gz CL29159-batch1_S48_L002_R1_001.fastq.gz CL29159-batch1_S48_L003_R1_001.fastq.gz CL29159-batch1_S48_L004_R1_001.fastq.gz > CL29159_5.fastq.gz 
#cat CL29168-batch2_S46_L001_R1_001.fastq.gz CL29168-batch2_S46_L002_R1_001.fastq.gz CL29168-batch2_S46_L003_R1_001.fastq.gz CL29168-batch2_S46_L004_R1_001.fastq.gz > CL29168_5.fastq.gz 
#cat CL29213-batch2_S42_L001_R1_001.fastq.gz CL29213-batch2_S42_L002_R1_001.fastq.gz CL29213-batch2_S42_L003_R1_001.fastq.gz CL29213-batch2_S42_L004_R1_001.fastq.gz > CL29213_5.fastq.gz 
#cat CL29231-batch2_S38_L001_R1_001.fastq.gz CL29231-batch2_S38_L002_R1_001.fastq.gz CL29231-batch2_S38_L003_R1_001.fastq.gz CL29231-batch2_S38_L004_R1_001.fastq.gz > CL29231_5.fastq.gz 
#cat CL29201-batch2_S43_L001_R1_001.fastq.gz CL29201-batch2_S43_L002_R1_001.fastq.gz CL29201-batch2_S43_L003_R1_001.fastq.gz CL29201-batch2_S43_L004_R1_001.fastq.gz > CL29201_5.fastq.gz 
#cat CL29219-batch2_S40_L001_R1_001.fastq.gz CL29219-batch2_S40_L002_R1_001.fastq.gz CL29219-batch2_S40_L003_R1_001.fastq.gz CL29219-batch2_S40_L004_R1_001.fastq.gz > CL29219_5.fastq.gz 
#cat CL29234-batch2_S37_L001_R1_001.fastq.gz CL29234-batch2_S37_L002_R1_001.fastq.gz CL29234-batch2_S37_L003_R1_001.fastq.gz CL29234-batch2_S37_L004_R1_001.fastq.gz > CL29234_5.fastq.gz 
#cat CL29246-batch2_S36_L001_R1_001.fastq.gz CL29246-batch2_S36_L002_R1_001.fastq.gz CL29246-batch2_S36_L003_R1_001.fastq.gz CL29246-batch2_S36_L004_R1_001.fastq.gz > CL29246_5.fastq.gz 
#cat CL29223-batch2_S39_L001_R1_001.fastq.gz CL29223-batch2_S39_L002_R1_001.fastq.gz CL29223-batch2_S39_L003_R1_001.fastq.gz CL29223-batch2_S39_L004_R1_001.fastq.gz > CL29223_5.fastq.gz 
#cat CL29218-batch2_S41_L001_R1_001.fastq.gz CL29218-batch2_S41_L002_R1_001.fastq.gz CL29218-batch2_S41_L003_R1_001.fastq.gz CL29218-batch2_S41_L004_R1_001.fastq.gz > CL29218_5.fastq.gz 
#cat CL29253-batch2_S35_L001_R1_001.fastq.gz CL29253-batch2_S35_L002_R1_001.fastq.gz CL29253-batch2_S35_L003_R1_001.fastq.gz CL29253-batch2_S35_L004_R1_001.fastq.gz > CL29253_5.fastq.gz 
#cat CL29254-batch2_S34_L001_R1_001.fastq.gz CL29254-batch2_S34_L002_R1_001.fastq.gz CL29254-batch2_S34_L003_R1_001.fastq.gz CL29254-batch2_S34_L004_R1_001.fastq.gz > CL29254_5.fastq.gz 
#cat CL29256-batch2_S33_L001_R1_001.fastq.gz CL29256-batch2_S33_L002_R1_001.fastq.gz CL29256-batch2_S33_L003_R1_001.fastq.gz CL29256-batch2_S33_L004_R1_001.fastq.gz > CL29256_5.fastq.gz 
#cat CL29257-batch2_S32_L001_R1_001.fastq.gz CL29257-batch2_S32_L002_R1_001.fastq.gz CL29257-batch2_S32_L003_R1_001.fastq.gz CL29257-batch2_S32_L004_R1_001.fastq.gz > CL29257_5.fastq.gz 
#cat CL29258-batch2_S31_L001_R1_001.fastq.gz CL29258-batch2_S31_L002_R1_001.fastq.gz CL29258-batch2_S31_L003_R1_001.fastq.gz CL29258-batch2_S31_L004_R1_001.fastq.gz > CL29258_5.fastq.gz 
#cat CL29264-batch2_S30_L001_R1_001.fastq.gz CL29264-batch2_S30_L002_R1_001.fastq.gz CL29264-batch2_S30_L003_R1_001.fastq.gz CL29264-batch2_S30_L004_R1_001.fastq.gz > CL29264_5.fastq.gz 
#cat CL29270-batch2_S29_L001_R1_001.fastq.gz CL29270-batch2_S29_L002_R1_001.fastq.gz CL29270-batch2_S29_L003_R1_001.fastq.gz CL29270-batch2_S29_L004_R1_001.fastq.gz > CL29270_5.fastq.gz 
#cat CL29271-batch3_S28_L001_R1_001.fastq.gz CL29271-batch3_S28_L002_R1_001.fastq.gz CL29271-batch3_S28_L003_R1_001.fastq.gz CL29271-batch3_S28_L004_R1_001.fastq.gz > CL29271_5.fastq.gz 
#cat CL29272-batch2_S27_L001_R1_001.fastq.gz CL29272-batch2_S27_L002_R1_001.fastq.gz CL29272-batch2_S27_L003_R1_001.fastq.gz CL29272-batch2_S27_L004_R1_001.fastq.gz > CL29272_5.fastq.gz 
#cat CL29274-batch3_S26_L001_R1_001.fastq.gz CL29274-batch3_S26_L002_R1_001.fastq.gz CL29274-batch3_S26_L003_R1_001.fastq.gz CL29274-batch3_S26_L004_R1_001.fastq.gz > CL29274_5.fastq.gz 
#cat CL29276-batch3_S25_L001_R1_001.fastq.gz CL29276-batch3_S25_L002_R1_001.fastq.gz CL29276-batch3_S25_L003_R1_001.fastq.gz CL29276-batch3_S25_L004_R1_001.fastq.gz > CL29276_5.fastq.gz 
#cat CL29289-batch3_S24_L001_R1_001.fastq.gz CL29289-batch3_S24_L002_R1_001.fastq.gz CL29289-batch3_S24_L003_R1_001.fastq.gz CL29289-batch3_S24_L004_R1_001.fastq.gz > CL29289_5.fastq.gz 
#cat CL29291-batch3_S23_L001_R1_001.fastq.gz CL29291-batch3_S23_L002_R1_001.fastq.gz CL29291-batch3_S23_L003_R1_001.fastq.gz CL29291-batch3_S23_L004_R1_001.fastq.gz > CL29291_5.fastq.gz 
#cat CL29297-batch3_S22_L001_R1_001.fastq.gz CL29297-batch3_S22_L002_R1_001.fastq.gz CL29297-batch3_S22_L003_R1_001.fastq.gz CL29297-batch3_S22_L004_R1_001.fastq.gz > CL29297_5.fastq.gz 
#cat CL29298-batch3_S21_L001_R1_001.fastq.gz CL29298-batch3_S21_L002_R1_001.fastq.gz CL29298-batch3_S21_L003_R1_001.fastq.gz CL29298-batch3_S21_L004_R1_001.fastq.gz > CL29298_5.fastq.gz 
#cat CL29303-batch3_S20_L001_R1_001.fastq.gz CL29303-batch3_S20_L002_R1_001.fastq.gz CL29303-batch3_S20_L003_R1_001.fastq.gz CL29303-batch3_S20_L004_R1_001.fastq.gz > CL29303_5.fastq.gz 
#cat CL29491-batch3_S19_L001_R1_001.fastq.gz CL29491-batch3_S19_L002_R1_001.fastq.gz CL29491-batch3_S19_L003_R1_001.fastq.gz CL29491-batch3_S19_L004_R1_001.fastq.gz > CL29491_5.fastq.gz 
#cat CL29498-batch3_S18_L001_R1_001.fastq.gz CL29498-batch3_S18_L002_R1_001.fastq.gz CL29498-batch3_S18_L003_R1_001.fastq.gz CL29498-batch3_S18_L004_R1_001.fastq.gz > CL29498_5.fastq.gz 
#cat CL29513-batch3_S17_L001_R1_001.fastq.gz CL29513-batch3_S17_L002_R1_001.fastq.gz CL29513-batch3_S17_L003_R1_001.fastq.gz CL29513-batch3_S17_L004_R1_001.fastq.gz > CL29513_5.fastq.gz 
#cat CL29519-batch3_S16_L001_R1_001.fastq.gz CL29519-batch3_S16_L002_R1_001.fastq.gz CL29519-batch3_S16_L003_R1_001.fastq.gz CL29519-batch3_S16_L004_R1_001.fastq.gz > CL29519_5.fastq.gz 
#cat CL29529-batch3_S15_L001_R1_001.fastq.gz CL29529-batch3_S15_L002_R1_001.fastq.gz CL29529-batch3_S15_L003_R1_001.fastq.gz CL29529-batch3_S15_L004_R1_001.fastq.gz > CL29529_5.fastq.gz 
#cat CL29539-batch3_S14_L001_R1_001.fastq.gz CL29539-batch3_S14_L002_R1_001.fastq.gz CL29539-batch3_S14_L003_R1_001.fastq.gz CL29539-batch3_S14_L004_R1_001.fastq.gz > CL29539_5.fastq.gz 
#cat CL29558-batch3_S13_L001_R1_001.fastq.gz CL29558-batch3_S13_L002_R1_001.fastq.gz CL29558-batch3_S13_L003_R1_001.fastq.gz CL29558-batch3_S13_L004_R1_001.fastq.gz > CL29558_5.fastq.gz 
#cat CL29663-batch3_S12_L001_R1_001.fastq.gz CL29663-batch3_S12_L002_R1_001.fastq.gz CL29663-batch3_S12_L003_R1_001.fastq.gz CL29663-batch3_S12_L004_R1_001.fastq.gz > CL29663_5.fastq.gz 
#cat CL29680-batch3_S11_L001_R1_001.fastq.gz CL29680-batch3_S11_L002_R1_001.fastq.gz CL29680-batch3_S11_L003_R1_001.fastq.gz CL29680-batch3_S11_L004_R1_001.fastq.gz > CL29680_5.fastq.gz 
#cat CL29687-batch3_S10_L001_R1_001.fastq.gz CL29687-batch3_S10_L002_R1_001.fastq.gz CL29687-batch3_S10_L003_R1_001.fastq.gz CL29687-batch3_S10_L004_R1_001.fastq.gz > CL29687_5.fastq.gz 
#cat CL29727-batch3_S9_L001_R1_001.fastq.gz CL29727-batch3_S9_L002_R1_001.fastq.gz CL29727-batch3_S9_L003_R1_001.fastq.gz CL29727-batch3_S9_L004_R1_001.fastq.gz > CL29727_5.fastq.gz 
#cat CL29770-batch3_S8_L001_R1_001.fastq.gz CL29770-batch3_S8_L002_R1_001.fastq.gz CL29770-batch3_S8_L003_R1_001.fastq.gz CL29770-batch3_S8_L004_R1_001.fastq.gz > CL29770_5.fastq.gz 
#cat CL29776-batch3_S7_L001_R1_001.fastq.gz CL29776-batch3_S7_L002_R1_001.fastq.gz CL29776-batch3_S7_L003_R1_001.fastq.gz CL29776-batch3_S7_L004_R1_001.fastq.gz > CL29776_5.fastq.gz 
#cat HS1-batch1_S6_L001_R1_001.fastq.gz HS1-batch1_S6_L002_R1_001.fastq.gz HS1-batch1_S6_L003_R1_001.fastq.gz HS1-batch1_S6_L004_R1_001.fastq.gz > HS1_5.fastq.gz 
#cat HS2-batch2_S5_L001_R1_001.fastq.gz HS2-batch2_S5_L002_R1_001.fastq.gz HS2-batch2_S5_L003_R1_001.fastq.gz HS2-batch2_S5_L004_R1_001.fastq.gz > HS2_5.fastq.gz 
#cat HS3-batch2_S4_L001_R1_001.fastq.gz HS3-batch2_S4_L002_R1_001.fastq.gz HS3-batch2_S4_L003_R1_001.fastq.gz HS3-batch2_S4_L004_R1_001.fastq.gz > HS3_5.fastq.gz 
#cat HS4-batch3_S3_L001_R1_001.fastq.gz HS4-batch3_S3_L002_R1_001.fastq.gz HS4-batch3_S3_L003_R1_001.fastq.gz HS4-batch3_S3_L004_R1_001.fastq.gz > HS4_5.fastq.gz 
#cat HS5-batch3_S2_L001_R1_001.fastq.gz HS5-batch3_S2_L002_R1_001.fastq.gz HS5-batch3_S2_L003_R1_001.fastq.gz HS5-batch3_S2_L004_R1_001.fastq.gz > HS5_5.fastq.gz 
#cat HS6-batch3_S1_L001_R1_001.fastq.gz HS6-batch3_S1_L002_R1_001.fastq.gz HS6-batch3_S1_L003_R1_001.fastq.gz HS6-batch3_S1_L004_R1_001.fastq.gz > HS6_5.fastq.gz #

## Cat 5 runs
#cat CL29075*fastq.gz > CL29075.fastq.gz 
#cat CL29107*fastq.gz > CL29107.fastq.gz 
#cat CL29077*fastq.gz > CL29077.fastq.gz 
#cat CL29151*fastq.gz > CL29151.fastq.gz 
#cat CL29148*fastq.gz > CL29148.fastq.gz 
#cat CL29150*fastq.gz > CL29150.fastq.gz 
#cat CL29085*fastq.gz > CL29085.fastq.gz 
#cat CL29118*fastq.gz > CL29118.fastq.gz 
#cat CL29137*fastq.gz > CL29137.fastq.gz 
#cat CL29161*fastq.gz > CL29161.fastq.gz 
#cat CL29170*fastq.gz > CL29170.fastq.gz 
#cat CL29180*fastq.gz > CL29180.fastq.gz 
#cat CL29159*fastq.gz > CL29159.fastq.gz 
#cat CL29168*fastq.gz > CL29168.fastq.gz 
#cat CL29213*fastq.gz > CL29213.fastq.gz 
#cat CL29231*fastq.gz > CL29231.fastq.gz 
#cat CL29201*fastq.gz > CL29201.fastq.gz 
#cat CL29219*fastq.gz > CL29219.fastq.gz 
#cat CL29234*fastq.gz > CL29234.fastq.gz 
#cat CL29246*fastq.gz > CL29246.fastq.gz 
#cat CL29223*fastq.gz > CL29223.fastq.gz 
#cat CL29218*fastq.gz > CL29218.fastq.gz 
#cat CL29253*fastq.gz > CL29253.fastq.gz 
#cat CL29254*fastq.gz > CL29254.fastq.gz 
#cat CL29256*fastq.gz > CL29256.fastq.gz 
#cat CL29257*fastq.gz > CL29257.fastq.gz 
#cat CL29258*fastq.gz > CL29258.fastq.gz 
#cat CL29264*fastq.gz > CL29264.fastq.gz 
#cat CL29270*fastq.gz > CL29270.fastq.gz 
#cat CL29271*fastq.gz > CL29271.fastq.gz 
#cat CL29272*fastq.gz > CL29272.fastq.gz 
#cat CL29274*fastq.gz > CL29274.fastq.gz 
#cat CL29276*fastq.gz > CL29276.fastq.gz 
#cat CL29289*fastq.gz > CL29289.fastq.gz 
#cat CL29291*fastq.gz > CL29291.fastq.gz 
#cat CL29297*fastq.gz > CL29297.fastq.gz 
#cat CL29298*fastq.gz > CL29298.fastq.gz 
#cat CL29303*fastq.gz > CL29303.fastq.gz 
#cat CL29491*fastq.gz > CL29491.fastq.gz 
#cat CL29498*fastq.gz > CL29498.fastq.gz 
#cat CL29513*fastq.gz > CL29513.fastq.gz 
#cat CL29519*fastq.gz > CL29519.fastq.gz 
#cat CL29529*fastq.gz > CL29529.fastq.gz 
#cat CL29539*fastq.gz > CL29539.fastq.gz 
#cat CL29558*fastq.gz > CL29558.fastq.gz 
#cat CL29663*fastq.gz > CL29663.fastq.gz 
#cat CL29680*fastq.gz > CL29680.fastq.gz 
#cat CL29687*fastq.gz > CL29687.fastq.gz 
#cat CL29727*fastq.gz > CL29727.fastq.gz 
#cat CL29770*fastq.gz > CL29770.fastq.gz 
#cat CL29776*fastq.gz > CL29776.fastq.gz 
#cat HS1*fastq.gz > HS1.fastq.gz 
#cat HS2*fastq.gz > HS2.fastq.gz 
#cat HS3*fastq.gz > HS3.fastq.gz 
#cat HS4*fastq.gz > HS4.fastq.gz 
#cat HS5*fastq.gz > HS5.fastq.gz 
#cat HS6*fastq.gz > HS6.fastq.gz 

# 2) mapping with kallisto:
# index
# reference from 3/5/20, 12:54:00 AM
kallisto index -i myHumanIndex Homo_sapiens.GRCh38.cdna.all.fa.gz

# mapping
## patients
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29075 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29075.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29075.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29107 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29107.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29107.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29077 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29077.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29077.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29151 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29151.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29151.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29148 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29148.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29148.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29150 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29150.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29150.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29085 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29085.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29085.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29118 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29118.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29118.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29137 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29137.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29137.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29161 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29161.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29161.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29170 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29170.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29170.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29180 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29180.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29180.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29159 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29159.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29159.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29168 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29168.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29168.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29213 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29213.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29213.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29231 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29231.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29231.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29201 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29201.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29201.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29219 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29219.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29219.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29234 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29234.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29234.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29246 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29246.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29246.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29223 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29223.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29223.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29218 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29218.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29218.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29253 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29253.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29253.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29254 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29254.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29254.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29256 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29256.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29256.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29257 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29257.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29257.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29258 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29258.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29258.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29264 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29264.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29264.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29270 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29270.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29270.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29271 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29271.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29271.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29272 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29272.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29272.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29274 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29274.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29274.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29276 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29276.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29276.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29289 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29289.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29289.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29291 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29291.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29291.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29297 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29297.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29297.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29298 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29298.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29298.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29303 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29303.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29303.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29491 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29491.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29491.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29498 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29498.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29498.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29513 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29513.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29513.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29519 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29519.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29519.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29529 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29529.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29529.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29539 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29539.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29539.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29558 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29558.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29558.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29663 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29663.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29663.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29680 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29680.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29680.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29687 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29687.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29687.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29727 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29727.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29727.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29770 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29770.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29770.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/CL29776 -b 60 -t 22 --single -l 250 -s 30 /home/amorimc/3rdDataset/merged_samples/CL29776.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/CL29776.log

## controls
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/HS1 -b 60 -t 22 -l 250 -s 30 --single /home/amorimc/3rdDataset/merged_samples/HS1.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/HS1.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/HS2 -b 60 -t 22 -l 250 -s 30 --single /home/amorimc/3rdDataset/merged_samples/HS2.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/HS2.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/HS3 -b 60 -t 22 -l 250 -s 30 --single /home/amorimc/3rdDataset/merged_samples/HS3.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/HS3.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/HS4 -b 60 -t 22 -l 250 -s 30 --single /home/amorimc/3rdDataset/merged_samples/HS4.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/HS4.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/HS5 -b 60 -t 22 -l 250 -s 30 --single /home/amorimc/3rdDataset/merged_samples/HS5.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/HS5.log
kallisto quant -i myHumanIndex -o /home/amorimc/3rdDataset/kallisto_outputs/human/HS6 -b 60 -t 22 -l 250 -s 30 --single /home/amorimc/3rdDataset/merged_samples/HS6.fastq.gz &> /home/amorimc/3rdDataset/kallisto_outputs/human/HS6.log

# 3) Checking quality:

#mkdir fastqc
#fastqc *.gz -t 24 -o fastqc
#multiqc -d .

echo "Finished"