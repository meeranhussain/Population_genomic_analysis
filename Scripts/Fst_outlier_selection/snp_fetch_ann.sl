#!/bin/bash -e
#SBATCH --account=uow03744
#SBATCH --partition=milan
#SBATCH --job-name=snp_fetch_ann
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output snp_fetch_ann_%j.out    # save the output into a file
#SBATCH --error snp_fetch_ann_%j.err     # save the error output into a file

# purge all other modules that may be loaded, and might interfare
module purge

ml snpEff/5.0-Java-11.0.4 BCFtools/1.19-GCC-11.3.0 R/4.3.2-foss-2023a


########## Task run ##########################

for i in HAM_IRE HAM_LIN HAM_MAN IRE_DUN IRE_LIN IRE_MAN LIN_DUN MAN_DUN MAN_LIN; 
do 

########################## Fetch regions ###############################################
#echo "Running Fetch regions....."
#awk 'NR > 1 {print $3"\t"$4"\t"$5}' ${i}_high_fst_windows_tp_1.txt > ${i}/01_fetch_snps/regions_list.txt;
#echo "Completed Fetch regions"

########## Select sample Population ####################################################
#echo "Selecting Population ....."

#Split the variable into two parts
  #part1=${i%_*}  # Extract everything before the underscore
  #part2=${i#*_}  # Extract everything after the underscore
   # Optional: print for debugging
  #echo "Processing $i -> Part1: $part1, Part2: $part2"
  # Use the parts in the awk command
#awk -v p1="$part1" -v p2="$part2" '$2 == p1 || $2 == p2 {print $1}' pop_map.txt >  ${i}/01_fetch_snps/selected_samples.txt
  

############# Fetch snps ##############################################################
#echo "Fetching snps associated with fetched regions....."

#bcftools view -v snps -m2 -M2 -R ${i}/01_fetch_snps/regions_list.txt --samples-file ${i}/01_fetch_snps/selected_samples.txt -Oz -o  ${i}/01_fetch_snps/${i}_hgfst_snps_regions.vcf.gz merged_fil_var_invar.vcf.gz;

############## Annotation ##############################################
#echo "Running annotation....."

#cd ${i}/01_fetch_snps

#java -Xmx4g -jar $EBROOTSNPEFF/snpEff.jar ann -c my_snpEff.config -v Maethio_gnm ${i}_hgfst_snps_regions.vcf.gz > ${i}_hgfst_snps_ann.vcf;

#cd ../../
#echo "Completed $i annotation....."

################## Fetch GO terms #############################
#cat  ${i}/01_fetch_snps/snpEff_genes.txt | grep -v "#" | cut -f3 | uniq >  ${i}/02_fetch_GOs/${i}_trans_idssnpeff.txt


######################### GO enrichment ##########################

Rscript GO_enrichment.r ${i};

done