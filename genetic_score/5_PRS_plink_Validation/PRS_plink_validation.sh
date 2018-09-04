## PRS Plink Validation 

cd /PRS/analysis/plink_scores

# First we have to create a q_score_file with false p-values according to each score of the analysis data
Rscript q_ranges_validation.R

sed -i 's/"// g' q_score_to_new_data.txt
sed -i 's/"// g' ranges_to_new_data.txt

cp q_score_to_new_data.txt ranges_to_new_data.txt LOGOR_to_score.txt /PRS/validation/plink_scores

cd /PRS/validation/plink_scores

#log(or) file have to be check to change the coding of the alleles in case the analysis data was coded by numbers
# and the validation data was coded by letters, or viceversa
Rscript logor_score_validation.R

sed -i 's/"// g' logor_new_data.txt

# Then we can compute the scores
plink --bfile ../QC/SNP_file_validation_10 --score logor_new_data.txt --q-score-file q_score_to_new_data.txt --q-score-range ranges_to_new_data.txt


# Once the different scores are computed the outputs are files with the name plink.(range_name).profile.

# For each of the ranges in ranges.txt we listed the SNPs taken into account (snp_score_[range_name].txt) because
# they could be less than the ones used in the original score and we perform a generalized linear model to compute 
# the rsquared of the model (model_[range_name].txt) to compare with the previous ones

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

# We make a barplot with the rsquared of the models and save them in the file rsquared_scores.txt to compare visually
Rscript barplot_rsquared.R