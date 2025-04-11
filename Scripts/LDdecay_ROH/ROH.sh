#!/bin/bash

######## GENO Filtering 
plink --vcf M_aethio_MAN_HAM_LIN_DUN_IRE_bialSNP.vcf --recode vcf-iid --geno 0.2 --out Maethio_MAF_geno  --allow-extra-chr --double-id --set-missing-var-ids @:#
  
####### Run plink ########  
plink --vcf Maethio_geno.vcf --homozyg --homozyg-window-snp 50 --homozyg-density 40 --homozyg-window-het 1 --homozyg-gap 100 --homozyg-kb 20 --homozyg-snp 50 --out ROH_results --allow-extra-chr





