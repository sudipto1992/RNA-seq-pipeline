
Guide to run this RNASeq pipeline.

Download SRR from NCBI SRA archive: 
    1. To enable, make sure you have the list of SRR in the config file.
    2. Flag downloadSRR: True
    3. If you fetching input files from SRR, flag iscompressed: False, because input srr files are automatically converted to fastq.

Current Pipeline:

1.Quality checking

|--Trimming
    -Adapter trimming
|--Quality checking
    -FastQC
    -MultQC
    
2.Buidling Genome files using hisat
    -Generating splice sites
    -Generating exons
    -Generating genome index files
    
3.Alignment
    -Using hisat
    -Output as bam
    
4.Gene count
    -Using featurecount of subread package

Input files required:

1. Create a folder named fastq_input and place the fastq samples in them
2. Create a folder named annotation and put the genome annotation file (gtf)
3. Create a folder name genome and put the genome sequence file (.fna)

=======
1. Create a folder fastq_input and place the fastq or fastq.gz samples in them. 
2. Put the genome sequence file (.fna or .fa or .fasta) in a folder named as genome.
3. Create a folder named annotation and place the genome annotation file (gtf) in it.


Make sure you change the following flag in the config.yaml

1. Put the sample numbers 
2. Whether your input samples are gz or fastq
3. Mention the gtf and genome file.
4. Flag the appropriate steps which you need to execute as True. Rest keep as False.


Future integrations:

1. Use STAR aligner as an alternative to hisat
2. Using htseq-count as an alternative subread package
3. Integrate DESeq2 


Command:

1. Create a virtual conda environment with snakemake installed
2. Activate the envrionment and run the below command
    snakemake --use-conda --cores <number> --scheduler greedy
    snakemake --use-conda --cores 4 --scheduler greedy 
