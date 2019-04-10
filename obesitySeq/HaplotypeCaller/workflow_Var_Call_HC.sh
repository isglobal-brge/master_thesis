#!/bin/bash -x
# Variant Calling GATK Haplotype caller 
# Isaac De la Hoz

##############################################################################################
#                              DIRECTORIES                                                   #
##############################################################################################
DWD="$(pwd)"
DIRBAM=$DWD/Data
DIRVCF=$DWD/VCF_HapCaller
GREF=$DWD/hg38.fa

##############################################################################################
#                              EXECUTABLES                                                   #
##############################################################################################

export PATH=$PATH:/home/isglobal.lan/idelahoz/data/software/gatk/gatk-4.1.1.0
gatk = gatk

##############################################################################################
#                              SAMPLE NAME AND MAIN PARAMETERS                               #
##############################################################################################

printf "Samples that will be analized: \n\n`ls $DIRBAM| grep "bam$"|grep -v "2F759.bam"`\n\n"

################################################################################
# 1. Haplotype Caller                                                          #
################################################################################

#Changing from chrM to chrMT 
sed 's/chrM/chrMT/g' hg38.fa 
#Reference dictionary creation
gatk CreateSequenceDictionary -R hg38.fa -O hg38.dict
#Creating index file
samtools faidx hg38.fa 

##HaplotypeCaller 
for BAM in `ls $DIRBAM| grep "bam$"|grep -v "2F759.bam"`
do
    gatk --java-options "-Xmx4g" HaplotypeCaller -R $GREF -I $BAM -O $DIRVCF/${BAM%".bam"}.raw.snps.indels.g.vcf -ERC GVCF &
done
wait
################################################################################
# 2. VAriant validation                                                        #
################################################################################

for FILE in `find $DIRVCF -name "*.raw.snps.indels.g.vcf"`
do
    gatk ValidateVariants -R $GREF -V $FILE &
done
wait
################################################################################
# 3. GVCF combining                                                            #
################################################################################

find $DIRVCF -name "*.raw.snps.indels.g.vcf" > $DWD/input.list
gatk CombineGVCFs -R $GREF -V $DWD/input.list -O $DIRVCF/RawVariants.vcf 

################################################################################
# 4. GVCF Genotyping                                                           #
################################################################################

mkdir $DIRVCF/finalVCF

gatk --java-options "-Xmx4g" GenotypeGVCFs -R $GREF -V $DIRVCF/RawVariants.vcf -O $DIRVCF/finalVCF/variants.vcf

################################################################################
# 4. Variant filtering via Hard filtering                                      #
################################################################################

gatk SelectVariants -V $DIRVCF/finalVCF/variants.vcf -select-type SNP -O $DIRVCF/finalVCF/variants.snps.vcf &
gatk SelectVariants -V $DIRVCF/finalVCF/variants.vcf -select-type INDEL -O $DIRVCF/finalVCF/variants.indels.vcf
wait
##SNPs filtration
gatk VariantFiltration -V $DIRVCF/finalVCF/variants.snps.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O $DIRVCF/finalVCF/variants.snps_filtered.vcf

##Indels filtration

gatk VariantFiltration -V $DIRVCF/finalVCF/variants.indels.vcf -filter "QD < 2.0" \
--filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O $DIRVCF/finalVCF/variants.indels_filtered.vcf

################################################################################
# 5. Merging filtered SNPs and indels                                      #
################################################################################

gatk MergeVcfs \
          -I $DIRVCF/finalVCF/variants.snps_filtered.vcf \
          -I $DIRVCF/finalVCF/variants.indels_filtered.vcf \
          -O $DIRVCF/finalVCF/variants_filtered.vcf

# Selection of unfiltered variants
gatk SelectVariants -V $DIRVCF/finalVCF/variants_filtered.vcf -exclude-filtered true -O $DIRVCF/finalVCF/variant_filtered.vcf

bgzip -c $DIRVCF/finalVCF/variants_filtered.vcf > $DIRVCF/finalVCF/variants_filtered.vcf.gz
tabix -f -p vcf $DIRVCF/finalVCF/variants_filtered.vcf.gz
################################################################################
# 6. Variants that have pass the filter are located to a table                 #
################################################################################

gatk VariantsToTable \
    -V $DIRVCF/finalVCF/variants_filtered.vcf \
    -F CHROM -F POS -F REF -F ALT -F TYPE -F AF \
    -O $DIRVCF/finalVCF/Variants.table
exit 0