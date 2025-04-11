#!/bin/bash -e
#SBATCH --account=uow03744
#SBATCH --partition=milan
#SBATCH --job-name=BCFtools
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=30G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output bcftools_%j.out    # save the output into a file
#SBATCH --error bcftools_%j.err     # save the error output into a file

# purge all other modules that may be loaded, and might interfare
module purge

## load BCF

module load BCFtools/1.19-GCC-11.3.0 VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1 Miniconda3/23.10.0-1
module load Python/3.11.6-foss-2023a
source /nesi/nobackup/uow03744/softwares/PIXY/bin/activate
### BCF run

bcftools mpileup --threads 12 -f /nesi/nobackup/uow03744/Population_genomics_Maethio/02_Popgen_NZ/02_reference/03_maethiopoides_NZI/Maethio_NZI.fasta -b list_bam_VC.txt --annotate FORMAT/DP | bcftools call --threads 12 -m -Oz -o Maethio_MAN_HAM_LIN_DUN_IRE.vcf.gz

##### Filter VCF
    
vcftools --gzvcf Maethio_MAN_HAM_LIN_DUN_IRE.vcf.gz --remove-indels --max-missing 0.8 --min-meanDP 20 --max-meanDP 500 --recode --stdout | gzip -c > Maethio_MAN_HAM_LIN_DUN_IRE_fil.vcf.gz

### Invariant

vcftools --gzvcf Maethio_MAN_HAM_LIN_DUN_IRE_fil.vcf.gz --max-maf 0 --recode --stdout | bgzip -c > Maethio_invariantsites_fil.vcf.gz


#### Variant sites

vcftools --gzvcf Maethio_MAN_HAM_LIN_DUN_IRE_fil.vcf.gz --mac 1 --maf 0.05 --recode --stdout | bgzip -c > Maethio_variantsites_fil.vcf.gz


### tabix
tabix Maethio_invariantsites_fil.vcf.gz
tabix Maethio_variantsites_fil.vcf.gz

### BCF concatenate
bcftools concat --allow-overlaps Maethio_variantsites_fil.vcf.gz Maethio_invariantsites_fil.vcf.gz -O z -o merged_fil_var_invar.vcf.gz

bcftools sort merged_fil_var_invar.vcf.gz -Oz -o merged_fil_var_invar_sorted.vcf.gz

tabix merged_fil_var_invar_sorted.vcf.gz


#RUN pixy

pixy --stats pi fst dxy --vcf merged_fil_var_invar_sorted.vcf.gz --populations pop_map.txt --window_size 100000 --n_cores 24 --output_prefix pixy_100kb_sort

