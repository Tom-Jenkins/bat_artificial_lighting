#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=24:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-T121362 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
##SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=t.l.jenkins@exeter.ac.uk # email address
##SBATCH -p highmem

# Variables
reads=/lustre/projects/Research_Project-T121362/nerc_lighting_rna/00_raw_reads
adapters=~/nextflow-pipelines/misc/adapters.fasta
trimmed_reads=/lustre/projects/Research_Project-T121362/nerc_lighting_rna/01_trimmed_reads
outdir=/lustre/projects/Research_Project-T121362/nerc_lighting_rna/02_expression
genome=/lustre/projects/Research_Project-T121362/Rhinolophus_hipposideros_genome/GCA_964194185.1_mRhiHip2.hap1.1_genomic.fna
transcriptome=/lustre/projects/Research_Project-T121362/Rhinolophus_hipposideros_genome/HLrhiHip1A_CDS.fasta
cpus=16

## 1. Trim raw reads

# # Activate conda environment
# source activate fastp

# # Run fastp pipeline
# nextflow run ~/nextflow-pipelines/src/fastp.nf \
#     --reads ${reads} \
#     --suffix "_{R1,R2}.fastq.gz" \
#     --adapters ${adapters} \
#     --filter "--qualified_quality_phred 30 --length_required 100 --trim_poly_g --trim_poly_x --dont_eval_duplication" \
#     --outdir ${trimmed_reads} \
#     --cpus ${cpus}

# # Move all output files into ${trimmed_reads} directory
# mv ${trimmed_reads}/trimmed_reads/* ${trimmed_reads}
# rm -r ${trimmed_reads}/trimmed_reads/

# # Re-run fastqc
# fastqc ${trimmed_reads}/*.fq.gz --outdir ${trimmed_reads} --threads ${cpus}

# # Run Multiqc
# multiqc ${trimmed_reads}


## 2. Quantify expression

# # Activate conda environment
# source activate quantifyexpression

# # Run pipeline
# nextflow run ~/nextflow-pipelines/src/quantifyexpression.nf \
#     --reads ${trimmed_reads} \
#     --suffix "_{1,2}.fp.fq.gz" \
#     --transcriptome ${transcriptome} \
#     --salmon \
#     --salmonParams "--libType A --gcBias --seqBias" \
#     --kallisto \
#     --outdir ${outdir} \
#     --cpus ${cpus}

## 3. Merge files using custom Python script
python ~/nextflow-pipelines/misc/merge_quant_files.py --salmon salmon_output/ salmon_quant_results.csv
python ~/nextflow-pipelines/misc/merge_quant_files.py --kallisto kallisto_output/ kallisto_quant_results.csv