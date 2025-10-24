#!/bin/bash

# Extract sample names from config
samples=$(yq e '.samples[]' config.yaml)

for sample in $samples
do
    echo "====================================="
    echo " Running Snakemake for SAMPLE: $sample"
    echo "====================================="

    # Create temporary single-sample config
    yq e ".samples = [\"$sample\"]" config.yaml > config_single.yaml

    # Run Snakemake for that sample
    snakemake --cores 8 --configfile config_single.yaml --use-conda --scheduler greedy

    # (OPTIONAL) Cleanup intermediate files
    # snakemake --cores 8 --configfile config_single.yaml --delete-temp-output

    echo "Finished $sample"
    echo ""
done

echo "ALL SAMPLES DONE"
