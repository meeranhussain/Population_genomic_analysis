# Population genomics
Steps associated with population genomic analysis compiled to a snake make workflow language
![image](https://github.com/meeranhussain/Population_genomics_analysis/assets/40800675/26de3f9c-b8c4-4b46-95a1-fd495d22d0cf)

# FOR SNAKEMAKE RUN
## Step 1: Project Folder
Create a project folder and give it a meaningful Project_ID.
```bash
mkdir <project-id>
```
## Step 2: Copy Files into Project Folder
Copy the following files into the project folder:
- `Snakefile`
- `Stats.R`
- `config.yaml`

## Step 3: Create a Sub-folder "01_Data"
Inside the project folder, create a sub-folder named `01_Data`.
```bash
mkdir <project-id>/01_Data
```
## Step 4: Copy Sample Files to 01_Data
Copy the sample files into the `01_Data` folder. 
```bash
cp *.fq <project-id>/01_Data
```
Ensure that the fastq files are named according to this pattern. **Ex: Featherston_01.fastq, Featherston_02.fastq, Mosburn_01.fastq**
#### Explanation of Example
To help you understand how to label the files correctly:
- `Each file should have a name followed by an underscore and a two-digit number.`
- `The name represents a specific population or sample, such as "Featherston" or "Mosburn".`
- `The two-digit number distinguishes different files from the same population or sample.`

## Step 5: Use Config File to Add Additional Information
Utilize the `config.yaml` file to add any additional information required for the workflow.
#### config.yaml content for snakemake workflow (Example file)

```yaml
##### Sequencing platform info (mostly keep this constant)
PL: "Illumina"
PM: "HISEQ"

##### Assign threads 
THREADS: 16

##### Provide path to reference file (Ensure reference is indexed using BWA index command and available in path provided)
fasta_path: "/nesi/nobackup/uow03744/Biocontrol_project/03_snakemake/02_reference/01_maethiopoides_IR/Maeth_IR_genomic.fasta"

######## Variant calling filter parameters ########
min_MQ: 20
min_BQ: 20
vcf_name: 'M_aethio_MOSS_LIN_FEA' #Used to assign names to output files generated in most steps
MAF: 'MAF > 0.05'

######## PLINK parameters #################
GENO: 0.1
```
## Step 6: Open Terminal in Project Folder
Navigate to the project folder in your terminal.

## Step 7: Run Snakemake
Type the following command in the terminal:
```bash
snakemake --configfile=config.yaml --cores 8 
```
