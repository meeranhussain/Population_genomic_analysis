import os
import subprocess
import pandas as pd

from os.path import join
#from pyfiglet import Figlet
from pathlib import Path
from glob import glob



#Reference path
fasta_path = config["fasta_path"]

# Select data files
SAMPLE_DIR = "01_Data"
SAMPLES, = glob_wildcards(SAMPLE_DIR + "/{sample}_r1.fq.gz")
R1 = '{sample}_r1.fq.gz'
R2 = '{sample}_r2.fq.gz'

#################### VCF name ############################
vcf_name_global = config['vcf_name']


################################################# Defined Function #############################################

################### mapped file listing funtion ##################
def list_bam_files(root_directory, output_file):
    try:
        with open(output_file, 'a') as file:  # Use 'a' for append mode

            # Walk through the directory tree
            for root, dirs, files in os.walk(root_directory):
                for file_name in files:
                    # Check if the file has the specified extension
                    if file_name.endswith("mapped.bam"):
                        base_name = os.path.splitext(file_name)[0]  # Remove the ".bam" extension
                        base_name_without_mapped = base_name.replace("mapped", "")  # Remove "mapped"
                        file_path = os.path.join(root, file_name)
                        file.write(f"{base_name_without_mapped}\t{file_path}\n")

        print(f"BAM file names (without extension) and paths have been appended to '{output_file}'.")

    except FileNotFoundError:
        print(f"Error: The specified root directory '{root_directory}' was not found.")


################### mapped file listing funtion for Variat calling (VC) ##################
def list_bam_files_sorted(root_directory, output_file):
    try:
        bam_paths = []  # Create an empty list to store BAM file paths

        # Walk through the directory tree
        for root, dirs, files in os.walk(root_directory):
            for file_name in files:
                # Check if the file has the specified extension
                if file_name.endswith("mapped.bam"):
                    file_path = os.path.join(root, file_name)
                    bam_paths.append(file_path)

        # Sort the list of BAM file paths based on file names
        bam_paths_sorted = sorted(bam_paths, key=lambda x: os.path.basename(x))

        # Write sorted paths to the output file
        with open(output_file, 'w') as file:  # Use 'w' to overwrite existing content
            for path in bam_paths_sorted:
                file.write(f"{path}\n")

        print(f"BAM file paths have been sorted and written to '{output_file}'.")

    except FileNotFoundError:
        print(f"Error: The specified root directory '{root_directory}' was not found.")




############################# Pop Map file creation for Stacks #############################
# Open the input file for reading
def pop_map(input_file, output_dir):
    try:
        # Open the input file for reading
        with open(input_file, "r") as infile:
            # Read lines from the file
            lines = infile.readlines()

        # Create a list to store the modified lines
        modified_lines = []

        # Process each line
        for line in lines:
            # Split the line into two columns
            columns = line.strip().split("\t")

            # Remove "_" and numeric characters from the first column
            first_column = columns[0].rstrip("0123456789_")

            # Remove trailing "_" and numeric characters from the original first column
            original_first_column = columns[0].rstrip("_")

            # Add the modified line to the list
            modified_lines.append([original_first_column, first_column])

        # Convert the list to a DataFrame
        df = pd.DataFrame(modified_lines, columns=['modified_first_column', 'original_first_column'])

        # Sort the DataFrame by the first column (modified_first_column)
        df_sorted = df.sort_values(by='modified_first_column')

        # Write the modified lines to an output file
        df_sorted.to_csv(output_dir, sep="\t", index=False, header=False)
        
    except FileNotFoundError:
        print(f"Error: '{output_dir}' was not created.")

##################################################################################################


############################ Rules ###############################
rule all:
    input: 
        expand("01_Data/galore/{sample}_r1_val_1.fq.gz", sample=SAMPLES),
        expand("01_Data/galore/{sample}_r2_val_2.fq.gz", sample=SAMPLES),
        expand("03_Analysis/01_BWA/{sample}/{sample}.sam", sample=SAMPLES),
#        expand("03_Analysis/01_BWA/{sample}/{sample}_Sorted.bam", sample=SAMPLES)
        expand("03_Analysis/01_BWA/{sample}/{sample}_markdup.bam", sample=SAMPLES),
        expand("03_Analysis/01_BWA/{sample}/{sample}_markdup.bam.bai", sample=SAMPLES),
        expand("03_Analysis/01_BWA/{sample}/{sample}_markdup.depth", sample=SAMPLES),
        expand("03_Analysis/01_BWA/{sample}/{sample}_markdup.stats", sample=SAMPLES),
        expand("03_Analysis/01_BWA/{sample}/{sample}_mapped.bam", sample=SAMPLES),
        '03_Analysis/01_BWA/list_bam.txt',
        '03_Analysis/01_BWA/list_bam_VC.txt',
        '03_Analysis/02_variant_calling',
        expand("03_Analysis/03_filtered_vcf/{vcf_name_global}_bialSNP_MAF.vcf", vcf_name_global=config['vcf_name']),
        expand("03_Analysis/04_plink/{vcf_name_global}_bialSNP_MAF_geno.vcf", vcf_name_global=config['vcf_name']),
        '03_Analysis/04_plink/pop_map.txt',
        '03_Analysis/05_stats'

        


rule trim:
    input:
        r1 = join(SAMPLE_DIR, R1),
        r2 = join(SAMPLE_DIR, R2)
    output:
        "01_Data/galore/{sample}_r1_val_1.fq.gz",
        "01_Data/galore/{sample}_r2_val_2.fq.gz"  
    shell:"""
        mkdir -p 01_Data/galore;
        mkdir -p 01_Data/fastqc;
        trim_galore --gzip --fastqc --fastqc_args "--outdir 01_Data/fastqc" -o 01_Data/galore --paired {input.r1} {input.r2} --cores 8"""

rule alignment:
    input: 
        r1 = "01_Data/galore/{sample}_r1_val_1.fq.gz",
        r2 = "01_Data/galore/{sample}_r2_val_2.fq.gz"
    params:
        fasta = fasta_path,
        PL = config['PL'],
        threads = config['THREADS'],
        PM = config['PM']
    output: 
        "03_Analysis/01_BWA/{sample}/{sample}.sam"
    message:
        "--- Mapping using BWA---"
    shell:"""
        mkdir -p 03_Analysis;
        mkdir -p 03_Analysis/01_BWA;
        mkdir -p 03_Analysis/01_BWA/{wildcards.sample};
        bwa mem -t {params.threads} -M -R "@RG\\tID:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:{params.PL}\\tPM:{params.PM}\\tSM:{wildcards.sample}" {params.fasta} {input.r1} {input.r2} > 03_Analysis/01_BWA/{wildcards.sample}/{wildcards.sample}.sam"""

#rule Sorting:
#    input:
#        "03_Analysis/01_BWA/{sample}/{sample}.sam"
#    output:
#        "03_Analysis/01_BWA/{sample}/{sample}_Sorted.bam"
#    threads: 16
#    message:
#        "---Sorting---"
#    shell:"""
#        samtools sort -n -@ 16 03_Analysis/01_BWA/{wildcards.sample}/{wildcards.sample}.sam -o 03_Analysis/01_BWA/{wildcards.sample}/{wildcards.sample}_Sorted.bam"""

rule rm_dup:
    input:
        "03_Analysis/01_BWA/{sample}/{sample}.sam"
    output:
        "03_Analysis/01_BWA/{sample}/{sample}_markdup.bam"
    params:
        threads = config['THREADS']
    shell:"""
        samtools fixmate -@ {params.threads} -m 03_Analysis/01_BWA/{wildcards.sample}/{wildcards.sample}.sam - | samtools sort -@ {params.threads} -O BAM | samtools markdup -@ {params.threads} - 03_Analysis/01_BWA/{wildcards.sample}/{wildcards.sample}_markdup.bam""" 

rule index_bam:
    input:
        "03_Analysis/01_BWA/{sample}/{sample}_markdup.bam"
    output:
        "03_Analysis/01_BWA/{sample}/{sample}_markdup.bam.bai"
    params:
        threads = config['THREADS']
    shell:"""
        samtools index -@ {params.threads} 03_Analysis/01_BWA/{wildcards.sample}/{wildcards.sample}_markdup.bam """
        
 
rule depth:
    input: 
        "03_Analysis/01_BWA/{sample}/{sample}_markdup.bam"
    output: 
        "03_Analysis/01_BWA/{sample}/{sample}_markdup.depth"
    params:
        threads = config['THREADS']
    message:
        "--- Depth check on all BAM files --- "
    shell:"""
        samtools depth -@ {params.threads} 03_Analysis/01_BWA/{wildcards.sample}/{wildcards.sample}_markdup.bam > 03_Analysis/01_BWA/{wildcards.sample}/{wildcards.sample}_markdup.depth"""

rule stats:
    input: 
        "03_Analysis/01_BWA/{sample}/{sample}_markdup.bam"
    output: 
        "03_Analysis/01_BWA/{sample}/{sample}_markdup.stats"
    params:
        threads = config['THREADS']
    message:
        "--- Stats check on all BAM files --- "
    shell:"""
        samtools flagstat -@ {params.threads} 03_Analysis/01_BWA/{wildcards.sample}/{wildcards.sample}_markdup.bam > 03_Analysis/01_BWA/{wildcards.sample}/{wildcards.sample}_markdup.stats"""

rule view:
    input:
        "03_Analysis/01_BWA/{sample}/{sample}_markdup.bam"
    output:
        "03_Analysis/01_BWA/{sample}/{sample}_mapped.bam"
    params:
        threads = config['THREADS']
    message:
        "--- Stats check on all BAM files --- "
    shell:"""
        samtools view -F 4 -b 03_Analysis/01_BWA/{wildcards.sample}/{wildcards.sample}_markdup.bam > 03_Analysis/01_BWA/{wildcards.sample}/{wildcards.sample}_mapped.bam"""


rule print_map:
    input:
        "03_Analysis/01_BWA"
    output:
        "03_Analysis/01_BWA/list_bam.txt"
    message:
        "--- Make file list for QC check of Mapped BAM files ---"
    run:
        list_bam_files(input[0], output[0])

rule print_BAM_VC:
    input:
        "03_Analysis/01_BWA"
    output:
        "03_Analysis/01_BWA/list_bam_VC.txt"
    message:
        "--- Make file list for Variant calling ---"
    run:
        list_bam_files_sorted(input[0], output[0])


rule bcf_mpileup:
    input: 
        "03_Analysis/01_BWA/list_bam_VC.txt"
    output:
        directory("03_Analysis/02_variant_calling")
    params:
        fasta = fasta_path,
        min_MQ = config['min_MQ'],
        threads = config['THREADS'],
        min_BQ = config['min_BQ'],
        vcf_name = config['vcf_name']
    message:
        "---- Running Variant calling using BCFtools mpileup----"
    shell:"""
        mkdir -p 03_Analysis/02_variant_calling;
        bcftools mpileup -Ou -f {params.fasta} --bam-list 03_Analysis/01_BWA/list_bam_VC.txt --threads {params.threads} --min-MQ {params.min_MQ} --min-BQ {params.min_BQ} | bcftools call -mv -Ob -o 03_Analysis/02_variant_calling/{params.vcf_name}.vcf --threads {params.threads}"""

rule bcf_filter:
    input: 
        "03_Analysis/02_variant_calling"
    output:
        "03_Analysis/03_filtered_vcf/{vcf_name_global}_bialSNP_MAF.vcf"
    params:
        threads = config['THREADS'],
        vcf_name = config['vcf_name'],
        MAF = config['MAF']
    message:
        "---- Running Filtration of VCF file----"
    shell:"""
        bcftools view -m2 -M2 -v snps -O b -o 03_Analysis/03_filtered_vcf/{params.vcf_name}_bialSNP.vcf 03_Analysis/02_variant_calling/{params.vcf_name}.vcf;
        bcftools filter -i '{params.MAF}' -O v -o 03_Analysis/03_filtered_vcf/{params.vcf_name}_bialSNP_MAF.vcf 03_Analysis/03_filtered_vcf/{params.vcf_name}_bialSNP.vcf;"""


rule plink:
    input: 
        "03_Analysis/03_filtered_vcf/{vcf_name_global}_bialSNP_MAF.vcf"
    output:
        "03_Analysis/04_plink/{vcf_name_global}_bialSNP_MAF_geno.vcf"
    params:
        GENO = config['GENO'],
        vcf_name = config['vcf_name']
    message:
        "---- Running Plink----"
    shell:"""
        plink2 --vcf 03_Analysis/03_filtered_vcf/{params.vcf_name}_bialSNP_MAF.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --geno {params.GENO} --recode vcf-iid --out 03_Analysis/04_plink/{params.vcf_name}_bialSNP_MAF_geno;
        plink2 --vcf 03_Analysis/04_plink/{params.vcf_name}_bialSNP_MAF_geno.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --indep-pairwise 50 5 0.2 --out 03_Analysis/04_plink/plink;
        plink2 --vcf 03_Analysis/04_plink/{params.vcf_name}_bialSNP_MAF_geno.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --exclude 03_Analysis/04_plink/plink.prune.out --recode vcf-iid --out 03_Analysis/04_plink/{params.vcf_name}_bialSNP_MAF_geno_LD;
"""

rule pop_map:
    input:
        "03_Analysis/01_BWA/list_bam.txt"
    output:
        "03_Analysis/04_plink/pop_map.txt"

    message:
        "--- Make file list for Variant calling ---"
    run:
        pop_map(input[0], output[0])

rule Rstats:
    input:
        "03_Analysis/04_plink/pop_map.txt"
    output:
        directory("03_Analysis/05_stats")
    params:
        vcf_name = config['vcf_name'],
        pop_map = "03_Analysis/04_plink/pop_map.txt"
    message:
        "-----Running Stats.r for statistical work--------"
    shell:"""
        mkdir -p 03_Analysis/05_stats;
        Rscript Stats.r 03_Analysis/04_plink/{params.vcf_name}_bialSNP_MAF_geno_LD.vcf {params.pop_map} {params.vcf_name}
    """

