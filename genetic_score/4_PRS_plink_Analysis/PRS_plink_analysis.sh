## PRS Plink Analysis

cd analysis/plink_scores

# To calculate a polygenic risk score (PRS) with Plink we have to prepare several files first with the data the program will take
# to compute the score. The files are listed below:

# OR_to_score.txt: SNP_ID, ref_allele and OR
# LOGOR_to_score.txt: SNP_ID, ref_allele and log(OR), it is obtained from OR_to_score.txt
# q_score_file.txt: SNP_ID and p-val
# ranges.txt: range_name lower_threshold upper_threshold (no header)


plink --bfile ../QC/SNP_file_10 --score LOGOR_to_score.txt --q-score-file q_score_file.txt --q-score-range ranges.txt


# According to the kind of association you performed the different files are obtained:

# To create ranges.txt we recommend to explore the range of p-val of your data and based on that choose 
# different thresholds to compute several scores, from more stringent thresholds to less stringent ones

# Direct Association

awk '{print $2,$9}' ../association/assoc_results.assoc > q_score_file.txt
awk '{print $2,$4,$10}' ../association/assoc_results.assoc > OR_to_score.txt
awk '$3=log($3)' OR_to_score.txt > LOGOR_to_score.txt # Edit the header using a text editor like nano or vim


# Correcting for Population Stratification (10 PC) Association

awk '{print $2,$9}' ../association/logistic_results.assoc.logistic > q_score_file.txt 
awk '{print $2,$4,$7}' ../association/logistic_results.assoc.logistic > OR_to_score.txt
awk '$3=log($3)' OR_to_score.txt > LOGOR_to_score.txt # Edit the header using a text editor like nano or vim


# Correcting for trio data (family datasets) Association

awk '{print $2,$10}' ../association/tdt_analysis.tdt > q_score_file.txt
awk '{print $2,$4,$8}' ../association/tdt_analysis.tdt > OR_to_score.txt 
awk '$3=log($3)' OR_to_score.txt > LOGOR_to_score.txt # Edit the header using a text editor like nano or vim


# Once the different scores are computed the outputs are files with the name plink.(range_name).profile.

# For each of the ranges in ranges.txt we listed the SNPs taken into account (snp_score_[range_name].txt) and 
# we perform a generalized linear model to compute the rsquared of the model (model_[range_name].txt)

for row in $(seq 1 $(wc -l ranges.txt | awk '{print $1}'));
do
echo $p
range=$(sed $row'q;d' ranges.txt | awk '{print $1}')
echo $range

p=$(sed $row'q;d' ranges.txt | awk '{print $3}')

awk -v pval="$p" '{ if ($2 <= pval) { print } }' q_score_file.txt > snp_score_${range}.txt

Rscript model_score.R ${range}
sed -i 's/"//g' snp_score_${range}.txt

echo "Number of SNP in Score:" >> model_${range}.txt 
echo $(wc -l snp_score_${range}.txt) >> model_${range}.txt

done

# We make a barplot with the rsquared of the models and save them in the file rsquared_scores.txt
Rscript barplot_rsquared.R
