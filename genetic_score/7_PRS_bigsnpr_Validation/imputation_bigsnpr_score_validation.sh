# Imputation and Bigsnpr scores analysis (machine-learning)

cd PRS/validation/imputation

bash imputation.sh

cd PRS/analysis/imputation/imputed/scores

cp bigsnpr_pvals_to_new_data.txt bigsnpr_ranges_to_new_data.txt betas.txt PRS/validation/imputation/imputed/scores

cd PRS/validation/imputation/imputed/scores

Rscript bigsnpr_score_validation.R

