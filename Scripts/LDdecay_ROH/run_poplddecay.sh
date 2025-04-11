#!/bin/bash

# Paths and Input Files
VCF_FILE="M_aethio_MAN_HAM_LIN_DUN_IRE_bialSNP_MAF_geno.vcf"
POP_LD_DECAY="/nesi/nobackup/uow03744/softwares/PopLDdecay/PopLDdecay"
PLOT_SCRIPT="/nesi/nobackup/uow03744/softwares/PopLDdecay/bin/Plot_MultiPop.pl"

# Population sample prefixes
declare -A POPULATIONS=( 
    ["IRE"]="IRELAND" 
    ["MAN"]="MANGONUI" 
    ["HAM"]="HAMILTON" 
    ["LIN"]="LINCOLN" 
    ["DUN"]="DUNEDIN" 
)

# Number of replicates
REPLICATES=1

# Create main output directory
MAIN_OUTPUT_DIR="PopLDdecay_Replicates_1"
mkdir -p $MAIN_OUTPUT_DIR

# Function to extract random samples
extract_samples() {
    POP=$1
    OUTFILE=$2

    # Extract sample names for the population
    bcftools query -l $VCF_FILE | grep "${POP}_" > tmp_${POP}_samples.txt

    # Randomly select 5 samples if more than 5 exist
    if [ "$POP" == "MAN" ]; then
        cp tmp_${POP}_samples.txt $OUTFILE  # Keep all MAN individuals
    else
        shuf -n 5 tmp_${POP}_samples.txt > $OUTFILE  # Randomly select 5 individuals
    fi

    rm tmp_${POP}_samples.txt  # Cleanup
}

# Loop through replicates
for ((rep=1; rep<=REPLICATES; rep++)); do
    echo "Running replicate $rep..."

    # Create a folder for this replicate
    REP_DIR="${MAIN_OUTPUT_DIR}/Replicate_${rep}"
    mkdir -p $REP_DIR

    PLOT_LIST="${REP_DIR}/all_plot_rep${rep}.txt"
    > $PLOT_LIST  # Clear previous plot list

    # Loop through populations and prepare subpopulation files
    for POP in "${!POPULATIONS[@]}"; do
        SUBPOP_FILE="${REP_DIR}/${POP}_rep${rep}.txt"

        echo "Selecting random samples for $POP (Replicate $rep)..."
        extract_samples "${POP}" "$SUBPOP_FILE"

        # Run PopLDdecay for each replicate and population
        echo "Running PopLDdecay for $POP (Replicate $rep)..."
        STAT_FILE="${REP_DIR}/${POP,,}_rep${rep}.stat.gz"
        $POP_LD_DECAY -InVCF $VCF_FILE -OutStat $STAT_FILE -SubPop $SUBPOP_FILE -MaxDist 10 -MAF 0.05 -Miss 0.2

        # Add output in two-column format: <stat_file> <PopulationName>
        echo "$STAT_FILE ${POPULATIONS[$POP]}" >> $PLOT_LIST
    done

    # Plotting LD decay for all populations in this replicate
    echo "Plotting LD decay for replicate $rep..."
    perl $PLOT_SCRIPT -inList $PLOT_LIST -output ${REP_DIR}/plot_all_rep${rep} -keepR -maxX 10
done

echo "Multiple sampling and analysis complete."
