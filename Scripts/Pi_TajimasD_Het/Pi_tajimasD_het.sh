####### Sample preparation ################
sed -i 's/##fileformat=VCFv4.3/##fileformat=VCFv4.2/' M_aethio_MAN_HAM_LIN_DUN_IRE_bialSNP_MAF_geno.vcf 

vcftools --vcf M_aethio_MAN_HAM_LIN_DUN_IRE_bialSNP_MAF_geno.vcf --keep ../00_pop_files/HAM_pop.txt --out HAM --recode 
vcftools --vcf M_aethio_MAN_HAM_LIN_DUN_IRE_bialSNP_MAF_geno.vcf --keep ../00_pop_files/INV_pop.txt --out DUN --recode 
vcftools --vcf M_aethio_MAN_HAM_LIN_DUN_IRE_bialSNP_MAF_geno.vcf --keep ../00_pop_files/NOR_pop.txt --out MAN --recode 
vcftools --vcf M_aethio_MAN_HAM_LIN_DUN_IRE_bialSNP_MAF_geno.vcf --keep ../00_pop_files/LIN_pop.txt --out LIN --recode 
vcftools --vcf M_aethio_MAN_HAM_LIN_DUN_IRE_bialSNP_MAF_geno.vcf --keep ../00_pop_files/IRE_pop.txt --out IRE --recode 

####### Pi Estimation ############
vcftools --vcf HAM.recode.vcf --window-pi 100000 --out HAM_pi_100kb
vcftools --vcf MAN.recode.vcf --window-pi 100000 --out MAN_pi_100kb
vcftools --vcf IRE.recode.vcf --window-pi 100000 --out IRE_pi_100kb
vcftools --vcf LIN.recode.vcf --window-pi 100000 --out LIN_pi_100kb
vcftools --vcf DUN.recode.vcf --window-pi 100000 --out DUN_pi_100kb


####### TajimasD Estimation ############
vcftools --vcf HAM.recode.vcf --TajimaD 100000 --out HAM_tjd_100kb
vcftools --vcf DUN.recode.vcf --TajimaD 100000 --out DUN_tjd_100kb
vcftools --vcf MAN.recode.vcf --TajimaD 100000 --out MAN_tjd_100kb
vcftools --vcf LIN.recode.vcf --TajimaD 100000 --out LIN_tjd_100kb
vcftools --vcf IRE.recode.vcf --TajimaD 100000 --out IRE_tjd_100kb

####### Individual Heterozygosity ##########

vcftools --vcf M_aethio_MAN_HAM_LIN_DUN_IRE_bialSNP_MAF_geno.vcf --het --out All_samplw_het

####### Convert obs(HOM) to obs(HET)
import pandas as pd
import sys

# Check if input file is provided
if len(sys.argv) != 2:
    print("Usage: python script.py <input_file.tsv>")
    sys.exit(1)

# Input file from command-line
input_file = sys.argv[1]

# Load TSV data into a DataFrame
df = pd.read_csv(input_file, sep="\t")

# Calculate Observed Heterozygosity (Obs_Het) and Expected Heterozygosity (Exp_Het)
df["Obs_Het"] = 1 - (df["O(HOM)"] / df["N_SITES"])
df["Exp_Het"] = 1 - (df["E(HOM)"] / df["N_SITES"])
df["F_inbreeding"] = 1 - (df["Obs_Het"] / df["Exp_Het"])

# Output file name
output_file = "heterozygosity_results.csv"
df.to_csv(output_file, index=False)

print(f"Heterozygosity results saved to {output_file}")