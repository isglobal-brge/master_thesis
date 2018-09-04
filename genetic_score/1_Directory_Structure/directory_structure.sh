#!/bin/sh

## Making directory structure

# PRS
#	- analysis
#		 - Original_Data
#		 - QC
#		 - plink_scores
#		 - association
#		 - imputation
#		 	- imputed
#		 		- scores
#		 - (database_score) # diseases with gene databases will need this directory to perform scores based on SNPassoc package
#	- validation
#		 - Original_Data
#		 - QC
#		 - plink_scores
#		 - imputation
#		 	- imputed
#		 		- scores
#		 - (database_score) # diseases with gene databases will need this directory to perform scores based on SNPassoc package


mkdir PRS

cd PRS

mkdir validation analysis

cd analysis

mkdir Original_Data QC plink_scores association imputation database_scores

cd imputation

mkdir imputed

cd imputed 

mkdir scores

cd ../../../validation

mkdir Original_Data QC plink_scores imputation database_scores

cd imputation

mkdir imputed

cd imputed 

mkdir scores