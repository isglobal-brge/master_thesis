#!/bin/sh

plink='/path/to/plink'

for i in `seq 1 22`
do
	$plink --bfile ../QC/SNP_file_10 --chr $i --make-bed --out SNP_file_chr${i}
done

echo files_to_impute > files_to_impute.txt
for i in *chr*.bim
do  
	lines=`wc -l $i | cut -f1 -d' '`
	if [ $lines -gt 30000 ]
	then 
		echo ${i%.*}
		name=${i%.*} 
		snp_init1=$(sed -n 1p $i | cut -f2)
		snp_end1=$(sed -n $(($lines/2))p $i | cut -f2)
		snp_init2=$(sed -n $((($lines/2)+1))p $i | cut -f2)
		snp_end2=$(sed -n ${lines}p $i | cut -f2)

$plink --noweb --bfile ${name} --from $snp_init1 --to $snp_end1 --make-bed --out ${name}_chunk1
$plink --noweb --bfile ${name} --from $snp_init2 --to $snp_end2 --make-bed --out ${name}_chunk2

echo ${name}_chunk1 >> files_to_impute.txt
echo ${name}_chunk2 >> files_to_impute.txt

lines1=`wc -l ${name}_chunk1.bim | cut -f1 -d' '`
lines2=`wc -l ${name}_chunk2.bim | cut -f1 -d' '`
total_lines=$((lines1+lines2))
echo 'lines chunks' $total_lines 'and lines original:' $lines
	else
		echo ${i%.*} >> files_to_impute.txt
	fi
done


## Then impute them with bigsnpr_imputation

Rscript bigsnpr_imputation.R

## Use this command if the files' names content spaces in them
## for i in *; do mv "${i}" "$(echo ${i} | tr -d ' ')"; done 

cd imputed

ls -1 *bed | sed -e 's/\..*$//' > files_to_merge.txt

file=$(head -n 1 files_to_merge.txt)

$plink --bfile $file --merge-list files_to_merge.txt --allow-no-sex --make-bed --out SNP_file_Imputed
