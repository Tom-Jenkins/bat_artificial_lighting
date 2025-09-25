#!/bin/bash

#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=72:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-T121362 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
##SBATCH -p highmem # High Memory
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=t.l.jenkins@exeter.ac.uk # email address

# Global variables
cpus=16
raw_reads=/lustre/projects/Research_Project-T121362/nerc_lighting_dna/00_raw_reads
trimmed_reads=/lustre/projects/Research_Project-T121362/nerc_lighting_dna/01_trimmed_reads
variantDIR=/lustre/projects/Research_Project-T121362/nerc_lighting_dna/02_variantcalling
genome=/lustre/projects/Research_Project-T121362/Rhinolophus_hipposideros_genome/GCA_964194185.1_mRhiHip2.hap1.1_genomic.fna

# ---------- #
# Trim reads using Fastp
# ---------- #

# Nextflow pipeline
# https://github.com/Tom-Jenkins/nextflow-pipelines/tree/main
# https://github.com/Tom-Jenkins/nextflow-pipelines/blob/main/docs/trim-paired-reads.md

# # Activate conda environment
# source activate fastp

# # Run Fastp
# nextflow run ~/nextflow-pipelines/src/fastp.nf \
#     --cpus ${cpus} \
#     --reads ${raw_reads} \
#     --suffix "*_{R1,R2}_001.fastq.gz" \
#     --adapters ~/nextflow-pipelines/misc/adapters.fasta \
#     --outdir ${trimmed_reads} \
#     --filter "--qualified_quality_phred 30 --length_required 100 --trim_poly_g"

# ---------- #
# Count number of reads
# ---------- #

# python ~/nextflow-pipelines/misc/count_reads_in_fastq.py ${raw_reads} raw_reads_counts.csv
# python ~/nextflow-pipelines/misc/count_reads_in_fastq.py ${trimmed_reads} trimmed_reads_counts.csv

# ---------- #
# Variant calling
# ---------- #

# Nextflow pipeline
# https://github.com/Tom-Jenkins/nextflow-pipelines/tree/main
# https://github.com/Tom-Jenkins/nextflow-pipelines/blob/main/docs/variant-calling.md

Activate conda environment
source activate variantcalling

# Index reference genome
# bowtie2-build ${genome} ${genome}

# Align reads to reference genome 
# nextflow run ~/nextflow-pipelines/src/variantcalling.nf \
#     --sampleSheet ${variantDIR}/sample_sheet_lesser_horseshoe.csv \
#     --genome ${genome} \
#     --outdir ${variantDIR} \
#     --alignOnly \
#     --cpus ${cpus}

# Create BAM list text file
# ls -d -1 ${variantDIR}/processed_bams/*.bam > bamlist.txt

# Create .fai index
# samtools faidx ${genome}

# Create BAM indexes
# samtools index -@ ${cpus} -M ${variantDIR}/processed_bams/*.bam

# Call variants
# freebayes-parallel <(fasta_generate_regions.py ${genome}.fai 100000) ${cpus} \
#         --fasta-reference ${genome} \
#         --bam-list ${variantDIR}/bamlist.txt \
#         --ploidy 2 --min-base-quality 20 --min-mapping-quality 20 --genotype-qualities -g 1000 \
#         > ${variantDIR}/variants_freebayes.vcf
# gzip ${variantDIR}/variants_freebayes.vcf

# Filter variants #1
# vcftools --gzvcf ${variantDIR}/variants_freebayes.vcf.gz --max-missing 0.50 --mac 3 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out filter1

# Assess missing data
# vcftools --vcf filter1.recode.vcf --missing-indv
# cat out.imiss | sort -n -k 4

# Filter samples with high missing data
# bcftools view \
#     --threads ${cpus} \
#     --samples ^BRO10,HID11 \
#     --output-type z \
#     --output "filter2.vcf.gz" \
#     ${variantDIR}/filter1.recode.vcf

# Count number of samples left
# bcftools query -l filter2.vcf.gz | wc -l

# Final filter(s)
# bcftools view \
#     --threads ${cpus} \
#     --output-type z \
#     -v snps \
#     --include "F_MISSING <= 0.10 & FORMAT/GQ >= 15 & MAF >= 0.03" \
#     --output "Rhipposideros_0.10missing_0.03maf.vcf.gz" \
#    ${variantDIR}/filter2.vcf.gz
# bcftools view \
#     --threads ${cpus} \
#     --output-type z \
#     -v snps \
#     --include "F_MISSING <= 0.15 & FORMAT/GQ >= 15 & MAF >= 0.03" \
#     --output "Rhipposideros_0.15missing_0.03maf.vcf.gz" \
#    ${variantDIR}/filter2.vcf.gz
# bcftools view \
#     --threads ${cpus} \
#     --output-type z \
#     -v snps \
#     --include "F_MISSING <= 0.20 & FORMAT/GQ >= 15 & MAF >= 0.03" \
#     --output "Rhipposideros_0.20missing_0.03maf.vcf.gz" \
#    ${variantDIR}/filter2.vcf.gz
# bcftools view \
#     --threads ${cpus} \
#     --output-type z \
#     -v snps \
#     --include "F_MISSING <= 0.10 & FORMAT/GQ >= 15 & MAC >= 5" \
#     --output "Rhipposideros_0.10missing_5mac.vcf.gz" \
#    ${variantDIR}/filter2.vcf.gz
# bcftools view \
#     --threads ${cpus} \
#     --output-type z \
#     -v snps \
#     --include "F_MISSING <= 0.10 & FORMAT/GQ >= 15 & MAC >= 4" \
#     --output "Rhipposideros_0.10missing_4mac.vcf.gz" \
#    ${variantDIR}/filter2.vcf.gz
# bcftools view \
#     --threads ${cpus} \
#     --output-type z \
#     -v snps \
#     --include "F_MISSING <= 0.10 & FORMAT/GQ >= 15 & MAC >= 3" \
#     --output "Rhipposideros_0.10missing_3mac.vcf.gz" \
#     ${variantDIR}/filter2.vcf.gz

# ---------- #
# ROHan: Runs of Homozygosity
# ---------- #

# Change to ROH directory
# cd /lustre/home/tj311/bats/nerc_lighting_dna/03_ROH

# # Extract stats after computations complete
# cat */*.summary.txt | grep -o -P "(?<=--out  ).{7}"
# cat */*.summary.txt | grep -o -P "(?<=Genome-wide theta outside ROH:\t).{9}"
# cat */*.summary.txt | grep -o -P "(?<=Genome-wide theta inc. ROH:\t).{9}"
# cat */*.summary.txt | grep -o -P "Segments in ROH\(\%\)\s+:\s+[0-9].[0-9]+" | cut -f 2
# cat */*.summary.txt | grep "Avg. length of ROH" | awk '{print $6}'

# # Save stats to text file
# (echo -e "Ind\tHet_outside_ROH\tHet_inc_ROH\tSegments_in_ROH_%\tAvg_length_of_ROH" && paste \
#     <(cat */*.summary.txt | grep -o -P "(?<=--out  ).{7}") \
#     <(cat */*.summary.txt | grep -o -P "(?<=Genome-wide theta outside ROH:\t).{9}") \
#     <(cat */*.summary.txt | grep -o -P "(?<=Genome-wide theta inc. ROH:\t).{9}") \
#     <(cat */*.summary.txt | grep -o -P "Segments in ROH\(\%\)\s+:\s+[0-9].[0-9]+" | cut -f 2) \
#     <(cat */*.summary.txt | grep "Avg. length of ROH" | awk '{print $6}') \
# ) > heterozygosity_ROHan.tsv
    
# Run for a single sample (~5 hours)
# mkdir BAT01 && cd BAT01 && rohan --rohmu 2e-5 -t ${cpus} --out BAT01 ${genome} ${variantDIR}/processed_bams/BAT01.rg.md.bam && cd ..

# Run for all samples (did not have enough time to automate so ran in batches per site)
# mkdir BAT02 && cd BAT02 && rohan --rohmu 2e-5 -t ${cpus} --out BAT02 ${genome} ${variantDIR}/processed_bams/BAT02.rg.md.bam && cd ..
# mkdir BAT03 && cd BAT03 && rohan --rohmu 2e-5 -t ${cpus} --out BAT03 ${genome} ${variantDIR}/processed_bams/BAT03.rg.md.bam && cd ..
# mkdir BAT04 && cd BAT04 && rohan --rohmu 2e-5 -t ${cpus} --out BAT04 ${genome} ${variantDIR}/processed_bams/BAT04.rg.md.bam && cd ..
# mkdir BAT05 && cd BAT05 && rohan --rohmu 2e-5 -t ${cpus} --out BAT05 ${genome} ${variantDIR}/processed_bams/BAT05.rg.md.bam && cd ..
# mkdir BAT06 && cd BAT06 && rohan --rohmu 2e-5 -t ${cpus} --out BAT06 ${genome} ${variantDIR}/processed_bams/BAT06.rg.md.bam && cd ..
# mkdir BAT07 && cd BAT07 && rohan --rohmu 2e-5 -t ${cpus} --out BAT07 ${genome} ${variantDIR}/processed_bams/BAT07.rg.md.bam && cd ..
# mkdir BAT08 && cd BAT08 && rohan --rohmu 2e-5 -t ${cpus} --out BAT08 ${genome} ${variantDIR}/processed_bams/BAT08.rg.md.bam && cd ..
# mkdir BAT09 && cd BAT09 && rohan --rohmu 2e-5 -t ${cpus} --out BAT09 ${genome} ${variantDIR}/processed_bams/BAT09.rg.md.bam && cd ..
# mkdir BAT10 && cd BAT10 && rohan --rohmu 2e-5 -t ${cpus} --out BAT10 ${genome} ${variantDIR}/processed_bams/BAT10.rg.md.bam && cd ..

# mkdir ARL01 && cd ARL01 && rohan --rohmu 2e-5 -t ${cpus} --out ARL01 ${genome} ${variantDIR}/processed_bams/ARL01.rg.md.bam && cd ..
# mkdir ARL02 && cd ARL02 && rohan --rohmu 2e-5 -t ${cpus} --out ARL02 ${genome} ${variantDIR}/processed_bams/ARL02.rg.md.bam && cd ..
# mkdir ARL03 && cd ARL03 && rohan --rohmu 2e-5 -t ${cpus} --out ARL03 ${genome} ${variantDIR}/processed_bams/ARL03.rg.md.bam && cd ..
# mkdir ARL04 && cd ARL04 && rohan --rohmu 2e-5 -t ${cpus} --out ARL04 ${genome} ${variantDIR}/processed_bams/ARL04.rg.md.bam && cd ..
# mkdir ARL05 && cd ARL05 && rohan --rohmu 2e-5 -t ${cpus} --out ARL05 ${genome} ${variantDIR}/processed_bams/ARL05.rg.md.bam && cd ..
# mkdir ARL06 && cd ARL06 && rohan --rohmu 2e-5 -t ${cpus} --out ARL06 ${genome} ${variantDIR}/processed_bams/ARL06.rg.md.bam && cd ..
# mkdir ARL07 && cd ARL07 && rohan --rohmu 2e-5 -t ${cpus} --out ARL07 ${genome} ${variantDIR}/processed_bams/ARL07.rg.md.bam && cd ..
# mkdir ARL08 && cd ARL08 && rohan --rohmu 2e-5 -t ${cpus} --out ARL08 ${genome} ${variantDIR}/processed_bams/ARL08.rg.md.bam && cd ..
# mkdir ARL09 && cd ARL09 && rohan --rohmu 2e-5 -t ${cpus} --out ARL09 ${genome} ${variantDIR}/processed_bams/ARL09.rg.md.bam && cd ..

# mkdir ARN01 && cd ARN01 && rohan --rohmu 2e-5 -t ${cpus} --out ARN01 ${genome} ${variantDIR}/processed_bams/ARN01.rg.md.bam && cd ..
# mkdir ARN02 && cd ARN02 && rohan --rohmu 2e-5 -t ${cpus} --out ARN02 ${genome} ${variantDIR}/processed_bams/ARN02.rg.md.bam && cd ..
# mkdir ARN03 && cd ARN03 && rohan --rohmu 2e-5 -t ${cpus} --out ARN03 ${genome} ${variantDIR}/processed_bams/ARN03.rg.md.bam && cd ..
# mkdir ARN04 && cd ARN04 && rohan --rohmu 2e-5 -t ${cpus} --out ARN04 ${genome} ${variantDIR}/processed_bams/ARN04.rg.md.bam && cd ..
# mkdir ARN05 && cd ARN05 && rohan --rohmu 2e-5 -t ${cpus} --out ARN05 ${genome} ${variantDIR}/processed_bams/ARN05.rg.md.bam && cd ..
# mkdir ARN06 && cd ARN06 && rohan --rohmu 2e-5 -t ${cpus} --out ARN06 ${genome} ${variantDIR}/processed_bams/ARN06.rg.md.bam && cd ..
# mkdir ARN07 && cd ARN07 && rohan --rohmu 2e-5 -t ${cpus} --out ARN07 ${genome} ${variantDIR}/processed_bams/ARN07.rg.md.bam && cd ..
# mkdir ARN08 && cd ARN08 && rohan --rohmu 2e-5 -t ${cpus} --out ARN08 ${genome} ${variantDIR}/processed_bams/ARN08.rg.md.bam && cd ..
# mkdir ARN09 && cd ARN09 && rohan --rohmu 2e-5 -t ${cpus} --out ARN09 ${genome} ${variantDIR}/processed_bams/ARN09.rg.md.bam && cd ..
# mkdir ARN10 && cd ARN10 && rohan --rohmu 2e-5 -t ${cpus} --out ARN10 ${genome} ${variantDIR}/processed_bams/ARN10.rg.md.bam && cd ..
# mkdir ARN11 && cd ARN11 && rohan --rohmu 2e-5 -t ${cpus} --out ARN11 ${genome} ${variantDIR}/processed_bams/ARN11.rg.md.bam && cd ..
# mkdir ARN12 && cd ARN12 && rohan --rohmu 2e-5 -t ${cpus} --out ARN12 ${genome} ${variantDIR}/processed_bams/ARN12.rg.md.bam && cd ..

# mkdir BRO01 && cd BRO01 && rohan --rohmu 2e-5 -t ${cpus} --out BRO01 ${genome} ${variantDIR}/processed_bams/BRO01.rg.md.bam && cd ..
# mkdir BRO02 && cd BRO02 && rohan --rohmu 2e-5 -t ${cpus} --out BRO02 ${genome} ${variantDIR}/processed_bams/BRO02.rg.md.bam && cd ..
# mkdir BRO03 && cd BRO03 && rohan --rohmu 2e-5 -t ${cpus} --out BRO03 ${genome} ${variantDIR}/processed_bams/BRO03.rg.md.bam && cd ..
# mkdir BRO04 && cd BRO04 && rohan --rohmu 2e-5 -t ${cpus} --out BRO04 ${genome} ${variantDIR}/processed_bams/BRO04.rg.md.bam && cd ..
# mkdir BRO05 && cd BRO05 && rohan --rohmu 2e-5 -t ${cpus} --out BRO05 ${genome} ${variantDIR}/processed_bams/BRO05.rg.md.bam && cd ..
# mkdir BRO06 && cd BRO06 && rohan --rohmu 2e-5 -t ${cpus} --out BRO06 ${genome} ${variantDIR}/processed_bams/BRO06.rg.md.bam && cd ..
# mkdir BRO07 && cd BRO07 && rohan --rohmu 2e-5 -t ${cpus} --out BRO07 ${genome} ${variantDIR}/processed_bams/BRO07.rg.md.bam && cd ..
# mkdir BRO08 && cd BRO08 && rohan --rohmu 2e-5 -t ${cpus} --out BRO08 ${genome} ${variantDIR}/processed_bams/BRO08.rg.md.bam && cd ..
# mkdir BRO10 && cd BRO10 && rohan --rohmu 2e-5 -t ${cpus} --out BRO10 ${genome} ${variantDIR}/processed_bams/BRO10.rg.md.bam && cd ..
# mkdir BRO11 && cd BRO11 && rohan --rohmu 2e-5 -t ${cpus} --out BRO11 ${genome} ${variantDIR}/processed_bams/BRO11.rg.md.bam && cd ..
# mkdir BRO12 && cd BRO12 && rohan --rohmu 2e-5 -t ${cpus} --out BRO12 ${genome} ${variantDIR}/processed_bams/BRO12.rg.md.bam && cd ..

# mkdir CAN01 && cd CAN01 && rohan --rohmu 2e-5 -t ${cpus} --out CAD01 ${genome} ${variantDIR}/processed_bams/CAD01.rg.md.bam && cd ..
# mkdir CAD02 && cd CAD02 && rohan --rohmu 2e-5 -t ${cpus} --out CAD02 ${genome} ${variantDIR}/processed_bams/CAD02.rg.md.bam && cd ..
# mkdir CAD03 && cd CAD03 && rohan --rohmu 2e-5 -t ${cpus} --out CAD03 ${genome} ${variantDIR}/processed_bams/CAD03.rg.md.bam && cd ..
# mkdir CAD04 && cd CAD04 && rohan --rohmu 2e-5 -t ${cpus} --out CAD04 ${genome} ${variantDIR}/processed_bams/CAD04.rg.md.bam && cd ..
# mkdir CAD05 && cd CAD05 && rohan --rohmu 2e-5 -t ${cpus} --out CAD05 ${genome} ${variantDIR}/processed_bams/CAD05.rg.md.bam && cd ..
# mkdir CAD06 && cd CAD06 && rohan --rohmu 2e-5 -t ${cpus} --out CAD06 ${genome} ${variantDIR}/processed_bams/CAD06.rg.md.bam && cd ..
# mkdir CAD07 && cd CAD07 && rohan --rohmu 2e-5 -t ${cpus} --out CAD07 ${genome} ${variantDIR}/processed_bams/CAD07.rg.md.bam && cd ..
# mkdir CAD08 && cd CAD08 && rohan --rohmu 2e-5 -t ${cpus} --out CAD08 ${genome} ${variantDIR}/processed_bams/CAD08.rg.md.bam && cd ..
# mkdir CAD09 && cd CAD09 && rohan --rohmu 2e-5 -t ${cpus} --out CAD09 ${genome} ${variantDIR}/processed_bams/CAD09.rg.md.bam && cd ..
# mkdir CAD10 && cd CAD10 && rohan --rohmu 2e-5 -t ${cpus} --out CAD10 ${genome} ${variantDIR}/processed_bams/CAD10.rg.md.bam && cd ..
# mkdir CAD11 && cd CAD11 && rohan --rohmu 2e-5 -t ${cpus} --out CAD11 ${genome} ${variantDIR}/processed_bams/CAD11.rg.md.bam && cd ..
# mkdir CAD12 && cd CAD12 && rohan --rohmu 2e-5 -t ${cpus} --out CAD12 ${genome} ${variantDIR}/processed_bams/CAD12.rg.md.bam && cd ..

# mkdir CAN01 && cd CAN01 && rohan --rohmu 2e-5 -t ${cpus} --out CAN01 ${genome} ${variantDIR}/processed_bams/CAN01.rg.md.bam && cd ..
# mkdir CAN03 && cd CAN03 && rohan --rohmu 2e-5 -t ${cpus} --out CAN03 ${genome} ${variantDIR}/processed_bams/CAN03.rg.md.bam && cd ..
# mkdir CAN05 && cd CAN05 && rohan --rohmu 2e-5 -t ${cpus} --out CAN05 ${genome} ${variantDIR}/processed_bams/CAN05.rg.md.bam && cd ..
# mkdir CAN06 && cd CAN06 && rohan --rohmu 2e-5 -t ${cpus} --out CAN06 ${genome} ${variantDIR}/processed_bams/CAN06.rg.md.bam && cd ..
# mkdir CAN07 && cd CAN07 && rohan --rohmu 2e-5 -t ${cpus} --out CAN07 ${genome} ${variantDIR}/processed_bams/CAN07.rg.md.bam && cd ..
# mkdir CAN08 && cd CAN08 && rohan --rohmu 2e-5 -t ${cpus} --out CAN08 ${genome} ${variantDIR}/processed_bams/CAN08.rg.md.bam && cd ..
# mkdir CAN11 && cd CAN11 && rohan --rohmu 2e-5 -t ${cpus} --out CAN11 ${genome} ${variantDIR}/processed_bams/CAN11.rg.md.bam && cd ..
# mkdir CAN12 && cd CAN12 && rohan --rohmu 2e-5 -t ${cpus} --out CAN12 ${genome} ${variantDIR}/processed_bams/CAN12.rg.md.bam && cd ..

# mkdir CAS01 && cd CAS01 && rohan --rohmu 2e-5 -t ${cpus} --out CAS01 ${genome} ${variantDIR}/processed_bams/CAS01.rg.md.bam && cd ..
# mkdir CAS02 && cd CAS02 && rohan --rohmu 2e-5 -t ${cpus} --out CAS02 ${genome} ${variantDIR}/processed_bams/CAS02.rg.md.bam && cd ..
# mkdir CAS03 && cd CAS03 && rohan --rohmu 2e-5 -t ${cpus} --out CAS03 ${genome} ${variantDIR}/processed_bams/CAS03.rg.md.bam && cd ..
# mkdir CAS05 && cd CAS05 && rohan --rohmu 2e-5 -t ${cpus} --out CAS05 ${genome} ${variantDIR}/processed_bams/CAS05.rg.md.bam && cd ..
# mkdir CAS06 && cd CAS06 && rohan --rohmu 2e-5 -t ${cpus} --out CAS06 ${genome} ${variantDIR}/processed_bams/CAS06.rg.md.bam && cd ..
# mkdir CAS07 && cd CAS07 && rohan --rohmu 2e-5 -t ${cpus} --out CAS07 ${genome} ${variantDIR}/processed_bams/CAS07.rg.md.bam && cd ..
# mkdir CAS08 && cd CAS08 && rohan --rohmu 2e-5 -t ${cpus} --out CAS08 ${genome} ${variantDIR}/processed_bams/CAS08.rg.md.bam && cd ..
# mkdir CAS09 && cd CAS09 && rohan --rohmu 2e-5 -t ${cpus} --out CAS09 ${genome} ${variantDIR}/processed_bams/CAS09.rg.md.bam && cd ..
# mkdir CAS10 && cd CAS10 && rohan --rohmu 2e-5 -t ${cpus} --out CAS10 ${genome} ${variantDIR}/processed_bams/CAS10.rg.md.bam && cd ..
# mkdir CAS11 && cd CAS11 && rohan --rohmu 2e-5 -t ${cpus} --out CAS11 ${genome} ${variantDIR}/processed_bams/CAS11.rg.md.bam && cd ..
# mkdir CAS12 && cd CAS12 && rohan --rohmu 2e-5 -t ${cpus} --out CAS12 ${genome} ${variantDIR}/processed_bams/CAS12.rg.md.bam && cd ..

# mkdir CLA01 && cd CLA01 && rohan --rohmu 2e-5 -t ${cpus} --out CLA01 ${genome} ${variantDIR}/processed_bams/CLA01.rg.md.bam && cd ..
# mkdir CLA02 && cd CLA02 && rohan --rohmu 2e-5 -t ${cpus} --out CLA02 ${genome} ${variantDIR}/processed_bams/CLA02.rg.md.bam && cd ..
# mkdir CLA03 && cd CLA03 && rohan --rohmu 2e-5 -t ${cpus} --out CLA03 ${genome} ${variantDIR}/processed_bams/CLA03.rg.md.bam && cd ..
# mkdir CLA04 && cd CLA04 && rohan --rohmu 2e-5 -t ${cpus} --out CLA04 ${genome} ${variantDIR}/processed_bams/CLA04.rg.md.bam && cd ..
# mkdir CLA05 && cd CLA05 && rohan --rohmu 2e-5 -t ${cpus} --out CLA05 ${genome} ${variantDIR}/processed_bams/CLA05.rg.md.bam && cd ..
# mkdir CLA06 && cd CLA06 && rohan --rohmu 2e-5 -t ${cpus} --out CLA06 ${genome} ${variantDIR}/processed_bams/CLA06.rg.md.bam && cd ..
# mkdir CLA07 && cd CLA07 && rohan --rohmu 2e-5 -t ${cpus} --out CLA07 ${genome} ${variantDIR}/processed_bams/CLA07.rg.md.bam && cd ..
# mkdir CLA08 && cd CLA08 && rohan --rohmu 2e-5 -t ${cpus} --out CLA08 ${genome} ${variantDIR}/processed_bams/CLA08.rg.md.bam && cd ..
# mkdir CLA09 && cd CLA09 && rohan --rohmu 2e-5 -t ${cpus} --out CLA09 ${genome} ${variantDIR}/processed_bams/CLA09.rg.md.bam && cd ..
# mkdir CLA10 && cd CLA10 && rohan --rohmu 2e-5 -t ${cpus} --out CLA10 ${genome} ${variantDIR}/processed_bams/CLA10.rg.md.bam && cd ..
# mkdir CLA11 && cd CLA11 && rohan --rohmu 2e-5 -t ${cpus} --out CLA11 ${genome} ${variantDIR}/processed_bams/CLA11.rg.md.bam && cd ..
# mkdir CLA12 && cd CLA12 && rohan --rohmu 2e-5 -t ${cpus} --out CLA12 ${genome} ${variantDIR}/processed_bams/CLA12.rg.md.bam && cd ..

# mkdir CRA01 && cd CRA01 && rohan --rohmu 2e-5 -t ${cpus} --out CRA01 ${genome} ${variantDIR}/processed_bams/CRA01.rg.md.bam && cd ..
# mkdir CRA02 && cd CRA02 && rohan --rohmu 2e-5 -t ${cpus} --out CRA02 ${genome} ${variantDIR}/processed_bams/CRA02.rg.md.bam && cd ..
# mkdir CRA03 && cd CRA03 && rohan --rohmu 2e-5 -t ${cpus} --out CRA03 ${genome} ${variantDIR}/processed_bams/CRA03.rg.md.bam && cd ..
# mkdir CRA04 && cd CRA04 && rohan --rohmu 2e-5 -t ${cpus} --out CRA04 ${genome} ${variantDIR}/processed_bams/CRA04.rg.md.bam && cd ..
# mkdir CRA05 && cd CRA05 && rohan --rohmu 2e-5 -t ${cpus} --out CRA05 ${genome} ${variantDIR}/processed_bams/CRA05.rg.md.bam && cd ..
# mkdir CRA06 && cd CRA06 && rohan --rohmu 2e-5 -t ${cpus} --out CRA06 ${genome} ${variantDIR}/processed_bams/CRA06.rg.md.bam && cd ..
# mkdir CRA07 && cd CRA07 && rohan --rohmu 2e-5 -t ${cpus} --out CRA07 ${genome} ${variantDIR}/processed_bams/CRA07.rg.md.bam && cd ..
# mkdir CRA08 && cd CRA08 && rohan --rohmu 2e-5 -t ${cpus} --out CRA08 ${genome} ${variantDIR}/processed_bams/CRA08.rg.md.bam && cd ..
# mkdir CRA09 && cd CRA09 && rohan --rohmu 2e-5 -t ${cpus} --out CRA09 ${genome} ${variantDIR}/processed_bams/CRA09.rg.md.bam && cd ..

# mkdir HEN01 && cd HEN01 && rohan --rohmu 2e-5 -t ${cpus} --out HEN01 ${genome} ${variantDIR}/processed_bams/HEN01.rg.md.bam && cd ..
# mkdir HEN02 && cd HEN02 && rohan --rohmu 2e-5 -t ${cpus} --out HEN02 ${genome} ${variantDIR}/processed_bams/HEN02.rg.md.bam && cd ..
# mkdir HEN03 && cd HEN03 && rohan --rohmu 2e-5 -t ${cpus} --out HEN03 ${genome} ${variantDIR}/processed_bams/HEN03.rg.md.bam && cd ..
# mkdir HEN04 && cd HEN04 && rohan --rohmu 2e-5 -t ${cpus} --out HEN04 ${genome} ${variantDIR}/processed_bams/HEN04.rg.md.bam && cd ..
# mkdir HEN05 && cd HEN05 && rohan --rohmu 2e-5 -t ${cpus} --out HEN05 ${genome} ${variantDIR}/processed_bams/HEN05.rg.md.bam && cd ..
# mkdir HEN06 && cd HEN06 && rohan --rohmu 2e-5 -t ${cpus} --out HEN06 ${genome} ${variantDIR}/processed_bams/HEN06.rg.md.bam && cd ..
# mkdir HEN07 && cd HEN07 && rohan --rohmu 2e-5 -t ${cpus} --out HEN07 ${genome} ${variantDIR}/processed_bams/HEN07.rg.md.bam && cd ..
# mkdir HEN08 && cd HEN08 && rohan --rohmu 2e-5 -t ${cpus} --out HEN08 ${genome} ${variantDIR}/processed_bams/HEN08.rg.md.bam && cd ..
# mkdir HEN09 && cd HEN09 && rohan --rohmu 2e-5 -t ${cpus} --out HEN09 ${genome} ${variantDIR}/processed_bams/HEN09.rg.md.bam && cd ..
# mkdir HEN10 && cd HEN10 && rohan --rohmu 2e-5 -t ${cpus} --out HEN10 ${genome} ${variantDIR}/processed_bams/HEN10.rg.md.bam && cd ..
# mkdir HEN11 && cd HEN11 && rohan --rohmu 2e-5 -t ${cpus} --out HEN11 ${genome} ${variantDIR}/processed_bams/HEN11.rg.md.bam && cd ..
# mkdir HEN12 && cd HEN12 && rohan --rohmu 2e-5 -t ${cpus} --out HEN12 ${genome} ${variantDIR}/processed_bams/HEN12.rg.md.bam && cd ..

# mkdir HID01 && cd HID01 && rohan --rohmu 2e-5 -t ${cpus} --out HID01 ${genome} ${variantDIR}/processed_bams/HID01.rg.md.bam && cd ..
# mkdir HID02 && cd HID02 && rohan --rohmu 2e-5 -t ${cpus} --out HID02 ${genome} ${variantDIR}/processed_bams/HID02.rg.md.bam && cd ..
# mkdir HID03 && cd HID03 && rohan --rohmu 2e-5 -t ${cpus} --out HID03 ${genome} ${variantDIR}/processed_bams/HID03.rg.md.bam && cd ..
# mkdir HID04 && cd HID04 && rohan --rohmu 2e-5 -t ${cpus} --out HID04 ${genome} ${variantDIR}/processed_bams/HID04.rg.md.bam && cd ..
# mkdir HID05 && cd HID05 && rohan --rohmu 2e-5 -t ${cpus} --out HID05 ${genome} ${variantDIR}/processed_bams/HID05.rg.md.bam && cd ..
# mkdir HID06 && cd HID06 && rohan --rohmu 2e-5 -t ${cpus} --out HID06 ${genome} ${variantDIR}/processed_bams/HID06.rg.md.bam && cd ..
# mkdir HID07 && cd HID07 && rohan --rohmu 2e-5 -t ${cpus} --out HID07 ${genome} ${variantDIR}/processed_bams/HID07.rg.md.bam && cd ..
# mkdir HID08 && cd HID08 && rohan --rohmu 2e-5 -t ${cpus} --out HID08 ${genome} ${variantDIR}/processed_bams/HID08.rg.md.bam && cd ..
# mkdir HID09 && cd HID09 && rohan --rohmu 2e-5 -t ${cpus} --out HID09 ${genome} ${variantDIR}/processed_bams/HID09.rg.md.bam && cd ..
# mkdir HID10 && cd HID10 && rohan --rohmu 2e-5 -t ${cpus} --out HID10 ${genome} ${variantDIR}/processed_bams/HID10.rg.md.bam && cd ..
# mkdir HID11 && cd HID11 && rohan --rohmu 2e-5 -t ${cpus} --out HID11 ${genome} ${variantDIR}/processed_bams/HID11.rg.md.bam && cd ..
# mkdir HID12 && cd HID12 && rohan --rohmu 2e-5 -t ${cpus} --out HID12 ${genome} ${variantDIR}/processed_bams/HID12.rg.md.bam && cd ..

# mkdir LAN01 && cd LAN01 && rohan --rohmu 2e-5 -t ${cpus} --out LAN01 ${genome} ${variantDIR}/processed_bams/LAN01.rg.md.bam && cd ..
# mkdir LAN02 && cd LAN02 && rohan --rohmu 2e-5 -t ${cpus} --out LAN02 ${genome} ${variantDIR}/processed_bams/LAN02.rg.md.bam && cd ..
# mkdir LAN03 && cd LAN03 && rohan --rohmu 2e-5 -t ${cpus} --out LAN03 ${genome} ${variantDIR}/processed_bams/LAN03.rg.md.bam && cd ..
# mkdir LAN04 && cd LAN04 && rohan --rohmu 2e-5 -t ${cpus} --out LAN04 ${genome} ${variantDIR}/processed_bams/LAN04.rg.md.bam && cd ..
# mkdir LAN05 && cd LAN05 && rohan --rohmu 2e-5 -t ${cpus} --out LAN05 ${genome} ${variantDIR}/processed_bams/LAN05.rg.md.bam && cd ..
# mkdir LAN06 && cd LAN06 && rohan --rohmu 2e-5 -t ${cpus} --out LAN06 ${genome} ${variantDIR}/processed_bams/LAN06.rg.md.bam && cd ..
# mkdir LAN07 && cd LAN07 && rohan --rohmu 2e-5 -t ${cpus} --out LAN07 ${genome} ${variantDIR}/processed_bams/LAN07.rg.md.bam && cd ..
# mkdir LAN08 && cd LAN08 && rohan --rohmu 2e-5 -t ${cpus} --out LAN08 ${genome} ${variantDIR}/processed_bams/LAN08.rg.md.bam && cd ..
# mkdir LAN09 && cd LAN09 && rohan --rohmu 2e-5 -t ${cpus} --out LAN09 ${genome} ${variantDIR}/processed_bams/LAN09.rg.md.bam && cd ..
# mkdir LAN10 && cd LAN10 && rohan --rohmu 2e-5 -t ${cpus} --out LAN10 ${genome} ${variantDIR}/processed_bams/LAN10.rg.md.bam && cd ..
# mkdir LAN11 && cd LAN11 && rohan --rohmu 2e-5 -t ${cpus} --out LAN11 ${genome} ${variantDIR}/processed_bams/LAN11.rg.md.bam && cd ..
# mkdir LAN12 && cd LAN12 && rohan --rohmu 2e-5 -t ${cpus} --out LAN12 ${genome} ${variantDIR}/processed_bams/LAN12.rg.md.bam && cd ..

# mkdir PAR01 && cd PAR01 && rohan --rohmu 2e-5 -t ${cpus} --out PAR01 ${genome} ${variantDIR}/processed_bams/PAR01.rg.md.bam && cd ..
# mkdir PAR02 && cd PAR02 && rohan --rohmu 2e-5 -t ${cpus} --out PAR02 ${genome} ${variantDIR}/processed_bams/PAR02.rg.md.bam && cd ..
# mkdir PAR03 && cd PAR03 && rohan --rohmu 2e-5 -t ${cpus} --out PAR03 ${genome} ${variantDIR}/processed_bams/PAR03.rg.md.bam && cd ..
# mkdir PAR04 && cd PAR04 && rohan --rohmu 2e-5 -t ${cpus} --out PAR04 ${genome} ${variantDIR}/processed_bams/PAR04.rg.md.bam && cd ..
# mkdir PAR05 && cd PAR05 && rohan --rohmu 2e-5 -t ${cpus} --out PAR05 ${genome} ${variantDIR}/processed_bams/PAR05.rg.md.bam && cd ..
# mkdir PAR06 && cd PAR06 && rohan --rohmu 2e-5 -t ${cpus} --out PAR06 ${genome} ${variantDIR}/processed_bams/PAR06.rg.md.bam && cd ..
# mkdir PAR07 && cd PAR07 && rohan --rohmu 2e-5 -t ${cpus} --out PAR07 ${genome} ${variantDIR}/processed_bams/PAR07.rg.md.bam && cd ..
# mkdir PAR08 && cd PAR08 && rohan --rohmu 2e-5 -t ${cpus} --out PAR08 ${genome} ${variantDIR}/processed_bams/PAR08.rg.md.bam && cd ..
# mkdir PAR09 && cd PAR09 && rohan --rohmu 2e-5 -t ${cpus} --out PAR09 ${genome} ${variantDIR}/processed_bams/PAR09.rg.md.bam && cd ..
# mkdir PAR10 && cd PAR10 && rohan --rohmu 2e-5 -t ${cpus} --out PAR10 ${genome} ${variantDIR}/processed_bams/PAR10.rg.md.bam && cd ..
# mkdir PAR11 && cd PAR11 && rohan --rohmu 2e-5 -t ${cpus} --out PAR11 ${genome} ${variantDIR}/processed_bams/PAR11.rg.md.bam && cd ..
# mkdir PAR12 && cd PAR12 && rohan --rohmu 2e-5 -t ${cpus} --out PAR12 ${genome} ${variantDIR}/processed_bams/PAR12.rg.md.bam && cd ..

# mkdir PEN01 && cd PEN01 && rohan --rohmu 2e-5 -t ${cpus} --out PEN01 ${genome} ${variantDIR}/processed_bams/PEN01.rg.md.bam && cd ..
# mkdir PEN02 && cd PEN02 && rohan --rohmu 2e-5 -t ${cpus} --out PEN02 ${genome} ${variantDIR}/processed_bams/PEN02.rg.md.bam && cd ..
# mkdir PEN03 && cd PEN03 && rohan --rohmu 2e-5 -t ${cpus} --out PEN03 ${genome} ${variantDIR}/processed_bams/PEN03.rg.md.bam && cd ..
# mkdir PEN04 && cd PEN04 && rohan --rohmu 2e-5 -t ${cpus} --out PEN04 ${genome} ${variantDIR}/processed_bams/PEN04.rg.md.bam && cd ..
# mkdir PEN05 && cd PEN05 && rohan --rohmu 2e-5 -t ${cpus} --out PEN05 ${genome} ${variantDIR}/processed_bams/PEN05.rg.md.bam && cd ..
# mkdir PEN06 && cd PEN06 && rohan --rohmu 2e-5 -t ${cpus} --out PEN06 ${genome} ${variantDIR}/processed_bams/PEN06.rg.md.bam && cd ..
# mkdir PEN07 && cd PEN07 && rohan --rohmu 2e-5 -t ${cpus} --out PEN07 ${genome} ${variantDIR}/processed_bams/PEN07.rg.md.bam && cd ..
# mkdir PEN08 && cd PEN08 && rohan --rohmu 2e-5 -t ${cpus} --out PEN08 ${genome} ${variantDIR}/processed_bams/PEN08.rg.md.bam && cd ..
# mkdir PEN09 && cd PEN09 && rohan --rohmu 2e-5 -t ${cpus} --out PEN09 ${genome} ${variantDIR}/processed_bams/PEN09.rg.md.bam && cd ..
# mkdir PEN10 && cd PEN10 && rohan --rohmu 2e-5 -t ${cpus} --out PEN10 ${genome} ${variantDIR}/processed_bams/PEN10.rg.md.bam && cd ..
# mkdir PEN11 && cd PEN11 && rohan --rohmu 2e-5 -t ${cpus} --out PEN11 ${genome} ${variantDIR}/processed_bams/PEN11.rg.md.bam && cd ..
# mkdir PEN12 && cd PEN12 && rohan --rohmu 2e-5 -t ${cpus} --out PEN12 ${genome} ${variantDIR}/processed_bams/PEN12.rg.md.bam && cd ..

# mkdir PYS01 && cd PYS01 && rohan --rohmu 2e-5 -t ${cpus} --out PYS01 ${genome} ${variantDIR}/processed_bams/PYS01.rg.md.bam && cd ..
# mkdir PYS02 && cd PYS02 && rohan --rohmu 2e-5 -t ${cpus} --out PYS02 ${genome} ${variantDIR}/processed_bams/PYS02.rg.md.bam && cd ..
# mkdir PYS03 && cd PYS03 && rohan --rohmu 2e-5 -t ${cpus} --out PYS03 ${genome} ${variantDIR}/processed_bams/PYS03.rg.md.bam && cd ..
# mkdir PYS04 && cd PYS04 && rohan --rohmu 2e-5 -t ${cpus} --out PYS04 ${genome} ${variantDIR}/processed_bams/PYS04.rg.md.bam && cd ..
# mkdir PYS05 && cd PYS05 && rohan --rohmu 2e-5 -t ${cpus} --out PYS05 ${genome} ${variantDIR}/processed_bams/PYS05.rg.md.bam && cd ..
# mkdir PYS06 && cd PYS06 && rohan --rohmu 2e-5 -t ${cpus} --out PYS06 ${genome} ${variantDIR}/processed_bams/PYS06.rg.md.bam && cd ..

# mkdir SHE01 && cd SHE01 && rohan --rohmu 2e-5 -t ${cpus} --out SHE01 ${genome} ${variantDIR}/processed_bams/SHE01.rg.md.bam && cd ..
# mkdir SHE02 && cd SHE02 && rohan --rohmu 2e-5 -t ${cpus} --out SHE02 ${genome} ${variantDIR}/processed_bams/SHE02.rg.md.bam && cd ..
# mkdir SHE03 && cd SHE03 && rohan --rohmu 2e-5 -t ${cpus} --out SHE03 ${genome} ${variantDIR}/processed_bams/SHE03.rg.md.bam && cd ..
# mkdir SHE07 && cd SHE07 && rohan --rohmu 2e-5 -t ${cpus} --out SHE07 ${genome} ${variantDIR}/processed_bams/SHE07.rg.md.bam && cd ..
# mkdir SHE08 && cd SHE08 && rohan --rohmu 2e-5 -t ${cpus} --out SHE08 ${genome} ${variantDIR}/processed_bams/SHE08.rg.md.bam && cd ..
# mkdir SHE09 && cd SHE09 && rohan --rohmu 2e-5 -t ${cpus} --out SHE09 ${genome} ${variantDIR}/processed_bams/SHE09.rg.md.bam && cd ..
# mkdir SHE10 && cd SHE10 && rohan --rohmu 2e-5 -t ${cpus} --out SHE10 ${genome} ${variantDIR}/processed_bams/SHE10.rg.md.bam && cd ..
# mkdir SHE11 && cd SHE11 && rohan --rohmu 2e-5 -t ${cpus} --out SHE11 ${genome} ${variantDIR}/processed_bams/SHE11.rg.md.bam && cd ..
# mkdir SHE12 && cd SHE12 && rohan --rohmu 2e-5 -t ${cpus} --out SHE12 ${genome} ${variantDIR}/processed_bams/SHE12.rg.md.bam && cd ..

# mkdir STO01 && cd STO01 && rohan --rohmu 2e-5 -t ${cpus} --out STO01 ${genome} ${variantDIR}/processed_bams/STO01.rg.md.bam && cd ..
# mkdir STO02 && cd STO02 && rohan --rohmu 2e-5 -t ${cpus} --out STO02 ${genome} ${variantDIR}/processed_bams/STO02.rg.md.bam && cd ..
# mkdir STO03 && cd STO03 && rohan --rohmu 2e-5 -t ${cpus} --out STO03 ${genome} ${variantDIR}/processed_bams/STO03.rg.md.bam && cd ..
# mkdir STO04 && cd STO04 && rohan --rohmu 2e-5 -t ${cpus} --out STO04 ${genome} ${variantDIR}/processed_bams/STO04.rg.md.bam && cd ..
# mkdir STO05 && cd STO05 && rohan --rohmu 2e-5 -t ${cpus} --out STO05 ${genome} ${variantDIR}/processed_bams/STO05.rg.md.bam && cd ..
# mkdir STO06 && cd STO06 && rohan --rohmu 2e-5 -t ${cpus} --out STO06 ${genome} ${variantDIR}/processed_bams/STO06.rg.md.bam && cd ..
# mkdir STO07 && cd STO07 && rohan --rohmu 2e-5 -t ${cpus} --out STO07 ${genome} ${variantDIR}/processed_bams/STO07.rg.md.bam && cd ..
# mkdir STO08 && cd STO08 && rohan --rohmu 2e-5 -t ${cpus} --out STO08 ${genome} ${variantDIR}/processed_bams/STO08.rg.md.bam && cd ..

# mkdir THE01 && cd THE01 && rohan --rohmu 2e-5 -t ${cpus} --out THE01 ${genome} ${variantDIR}/processed_bams/THE01.rg.md.bam && cd ..
# mkdir THE02 && cd THE02 && rohan --rohmu 2e-5 -t ${cpus} --out THE02 ${genome} ${variantDIR}/processed_bams/THE02.rg.md.bam && cd ..
# mkdir THE03 && cd THE03 && rohan --rohmu 2e-5 -t ${cpus} --out THE03 ${genome} ${variantDIR}/processed_bams/THE03.rg.md.bam && cd ..
# mkdir THE05 && cd THE05 && rohan --rohmu 2e-5 -t ${cpus} --out THE05 ${genome} ${variantDIR}/processed_bams/THE05.rg.md.bam && cd ..
# mkdir THE06 && cd THE06 && rohan --rohmu 2e-5 -t ${cpus} --out THE06 ${genome} ${variantDIR}/processed_bams/THE06.rg.md.bam && cd ..
# mkdir THE07 && cd THE07 && rohan --rohmu 2e-5 -t ${cpus} --out THE07 ${genome} ${variantDIR}/processed_bams/THE07.rg.md.bam && cd ..
# mkdir THE08 && cd THE08 && rohan --rohmu 2e-5 -t ${cpus} --out THE08 ${genome} ${variantDIR}/processed_bams/THE08.rg.md.bam && cd ..
# mkdir THE09 && cd THE09 && rohan --rohmu 2e-5 -t ${cpus} --out THE09 ${genome} ${variantDIR}/processed_bams/THE09.rg.md.bam && cd ..
# mkdir THE10 && cd THE10 && rohan --rohmu 2e-5 -t ${cpus} --out THE10 ${genome} ${variantDIR}/processed_bams/THE10.rg.md.bam && cd ..
# mkdir THE11 && cd THE11 && rohan --rohmu 2e-5 -t ${cpus} --out THE11 ${genome} ${variantDIR}/processed_bams/THE11.rg.md.bam && cd ..
# mkdir THE12 && cd THE12 && rohan --rohmu 2e-5 -t ${cpus} --out THE12 ${genome} ${variantDIR}/processed_bams/THE12.rg.md.bam && cd ..

# mkdir WES01 && cd WES01 && rohan --rohmu 2e-5 -t ${cpus} --out WES01 ${genome} ${variantDIR}/processed_bams/WES01.rg.md.bam && cd ..
# mkdir WES02 && cd WES02 && rohan --rohmu 2e-5 -t ${cpus} --out WES02 ${genome} ${variantDIR}/processed_bams/WES02.rg.md.bam && cd ..
# mkdir WES03 && cd WES03 && rohan --rohmu 2e-5 -t ${cpus} --out WES03 ${genome} ${variantDIR}/processed_bams/WES03.rg.md.bam && cd ..
# mkdir WES04 && cd WES04 && rohan --rohmu 2e-5 -t ${cpus} --out WES04 ${genome} ${variantDIR}/processed_bams/WES04.rg.md.bam && cd ..
# mkdir WES05 && cd WES05 && rohan --rohmu 2e-5 -t ${cpus} --out WES05 ${genome} ${variantDIR}/processed_bams/WES05.rg.md.bam && cd ..
# mkdir WES06 && cd WES06 && rohan --rohmu 2e-5 -t ${cpus} --out WES06 ${genome} ${variantDIR}/processed_bams/WES06.rg.md.bam && cd ..
# mkdir WES07 && cd WES07 && rohan --rohmu 2e-5 -t ${cpus} --out WES07 ${genome} ${variantDIR}/processed_bams/WES07.rg.md.bam && cd ..
# mkdir WES08 && cd WES08 && rohan --rohmu 2e-5 -t ${cpus} --out WES08 ${genome} ${variantDIR}/processed_bams/WES08.rg.md.bam && cd ..
# mkdir WES09 && cd WES09 && rohan --rohmu 2e-5 -t ${cpus} --out WES09 ${genome} ${variantDIR}/processed_bams/WES09.rg.md.bam && cd ..
# mkdir WES10 && cd WES10 && rohan --rohmu 2e-5 -t ${cpus} --out WES10 ${genome} ${variantDIR}/processed_bams/WES10.rg.md.bam && cd ..
# mkdir WES11 && cd WES11 && rohan --rohmu 2e-5 -t ${cpus} --out WES11 ${genome} ${variantDIR}/processed_bams/WES11.rg.md.bam && cd ..
# mkdir WES12 && cd WES12 && rohan --rohmu 2e-5 -t ${cpus} --out WES12 ${genome} ${variantDIR}/processed_bams/WES12.rg.md.bam && cd ..
