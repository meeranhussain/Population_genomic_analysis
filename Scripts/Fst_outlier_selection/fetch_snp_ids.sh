#!/bin/bash


# Loop through each VCF file in the directory
for i in IRE_vsManLinDun HAM_vsManLinDun; do

    ######## Fetch regions #######################  
  awk 'NR > 1 {print $3"\t"$4"\t"$5}' ${i}/${i}_high_fst_windows_tp_1.txt > ${i}/${i}_regions_list.txt

  ########## Select sample Population ####################################################
  # Split the variable into two parts
  part1=${i%_*}  # Extract everything before the underscore
  part2=${i#*_}  # Extract everything after the underscore

  # Optional: print for debugging
  echo "Processing $i -> Part1: $part1, Part2: $part2"

  # Use awk with -v to pass shell variables
  awk -v p1="$part1" -v p2="$part2" '$2 == p1 || $2 == p2 {print $1}' pop_map.txt > ${i}/selected_samples.txt
  ############# Fetch snps from vcf ##############################################################
  bcftools view  -v snps -m2 -M2 -R ${i}/${i}_regions_list.txt --samples-file ${i}/selected_samples.txt -Oz -o ${i}/${i}_hgfst_snps_regions.vcf.gz ./merged_fil_var_invar.vcf.gz;
  
  ################# Fetch snp IDS ##################################################################
  # Define the output file name
   output_file="${i}/${i}_Fst_SNPIDs.txt"

   # Process the VCF file to extract CHROM and POS and format as "CHROM:POS"
   zcat "${i}/${i}_hgfst_snps_regions.vcf.gz" | awk -F'\t' 'BEGIN {OFS=":"} !/^#/ {print $1, $2}' > "$output_file"
   echo "All files processed. Results saved in $output_file."

done

