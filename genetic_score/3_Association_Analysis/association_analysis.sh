## Association Analysis

cd PRS_analysis/association

# In this script we will perform the association analysis (3 different ways) between the SNPs and a binary trait (case-control)
# And with the outputs we will prepare the files to the next analysis, the PRS calculation


# Direct Association

plink --bfile ../QC/SNP_file_10 --assoc --out assoc_results

# output: assoc_results.assoc, this file contains for each SNP different data including the p-value the Odd ratio and 
# specific information of this kind of associacion


# Correcting for Population Stratification (10 PC) Association

plink --bfile ../QC/SNP_file_10 --extract ../QC/indepSNP.prune.in --genome --out SNP_file_10
plink --bfile ../QC/SNP_file_10 --read-genome SNP_file_10.genome --cluster --mds-plot 10 --out SNP_file_10_mds

awk '{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' SNP_file_10_mds.mds > covar_mds.txt

plink --bfile ../QC/SNP_file_10 --covar covar_mds.txt --logistic --hide-covar --out logistic_results

# output: logistic_results.assoc.logistic, this file contains for each SNP different data including the p-value the Odd ratio and 
# specific information of this kind of associacion

# Correcting for trio data (family datasets) Association

plink --bfile ../QC/SNP_file_10 --tdt --out tdt_analysis

# output: tdt_analysis.tdt, this file contains for each SNP different data including the p-value the Odd ratio and 
# specific information of this kind of associacion

