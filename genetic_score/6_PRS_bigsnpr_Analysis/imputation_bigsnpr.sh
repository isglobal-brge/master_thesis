# Imputation and Bigsnpr scores analysis (machine-learning)

cd PRS/analysis/imputation

bash imputation.sh

cd imputed

Rscript bigsnpr_score.R

cd scores

for i in *.txt
do
sed -i 's/"//g' $i
done
