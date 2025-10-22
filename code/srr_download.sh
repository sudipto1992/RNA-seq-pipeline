#!/bin/bash

# Output directories
RAW_SRA_DIR="sra_raw"
FASTQ_DIR="fastq_input"

# Make folders
mkdir -p $RAW_SRA_DIR
mkdir -p $FASTQ_DIR

# Input list of SRR IDs
SRR_LIST="SRR_Acc_List.txt"

# Loop through each SRR ID
while read SRR; do
    echo "### Processing $SRR ###"

    # Step 1: Download SRA
    prefetch $SRR --output-directory $RAW_SRA_DIR

    # Step 2: Convert to FASTQ (paired-end or single-end auto-detect)
    fasterq-dump $RAW_SRA_DIR/$SRR/$SRR.sra --threads 8 --outdir $FASTQ_DIR

    # Step4: Remove SRR
    rm -r ${RAW_SRA_DIR}/${SRR}

    echo "### Completed $SRR ###"
    echo
done < $SRR_LIST
