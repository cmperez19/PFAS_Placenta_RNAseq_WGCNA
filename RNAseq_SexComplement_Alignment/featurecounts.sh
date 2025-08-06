#!/bin/bash
#SBATCH --partition=week-long-cpu
#SBATCH --time=02-00:00:00
#SBATCH --job-name="readcounts"
#SBATCH --mem=20GB
#SBATCH --output=readcounts.out
#SBATCH --error=readcounts.err
#SBATCH --mail-user=cmpere2@emory.edu
#SBATCH --mail-type=BEGIN,END

# Activate conda environment
#source ~/.bashrc 
#conda init --all
#conda activate ~/home/cmpere2/miniconda3/envs/featurecounts
output_dir="/projects/eversongroup/glowing/transcriptcounts"
for bam_file in "$output_dir"/G-*.sorted.markdup.bam; do
    sample_id=$(basename "$bam_file" | cut -d. -f1)
    /home/cmpere2/miniconda3/pkgs/subread-2.0.6-he4a0461_2/bin/featureCounts -T 8 --primary -s 0 -t gene -g gene_name \
             -a /projects/eversongroup/glowing/Homo_sapiens.GRCh38.92.gtf \
              -o "${output_dir}/${sample_id}_FC.txt" \
               "${output_dir}/${sample_id}.sorted.markdup.bam"
done
