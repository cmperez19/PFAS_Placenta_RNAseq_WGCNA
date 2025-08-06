#!/bin/bash
#SBATCH --partition=week-long-cpu
#SBATCH --time=02-00:00:00
#SBATCH --job-name="bam_stuff_cp"
#SBATCH --mem=20GB
#SBATCH --output=bam_stuff.out
#SBATCH --error=bam_stuff.err
#SBATCH --mail-user=cmpere2@emory.edu
#SBATCH --mail-type=BEGIN,END

# Activate conda environment
#conda activate featurecounts

# Set PICARD tool path
PICARD=/home/cmpere2/picard-3.1.1/picard.jar

# Define the directory containing the BAM files and output directory
bam_directory="/projects/eversongroup/glowing/alignment_output"
output_dir="/projects/eversongroup/glowing/transcriptcounts"

# Loop through each .bam file in the directory
for bam_file in "$bam_directory"/G-*.bam; do
    # Check if the file is a BAM file and starts with the required prefix
    if [ -f "$bam_file" ]; then
        filename=$(basename -- "$bam_file")
        
        # Extract the numerical part after "G-" and before the first underscore
        prefix=$(echo "$filename" | grep -oP '^G-\d+')

        # Check if the prefix starts with G-250 or higher
       # if [[ $prefix =~ ^G-[2-9][5-9][0-9] ]]; then
	#if [[ $prefix -lt 250 ]]; then
            # Extract the parts needed for the sample ID
            part1=$(echo "$filename" | cut -d'_' -f1)
            part2=$(echo "$filename" | cut -d'_' -f5 | cut -d'.' -f1)
            
            # Construct the sample ID
            sample_id="${part1}_${part2}"
            
            # Sort the BAM file using bamtools
            bamtools sort -in "$bam_file" -out "${output_dir}/${sample_id}.sorted.bam"
            
            # Mark duplicates using Picard
            java -Xmx8g -jar $PICARD MarkDuplicates \
                INPUT="${output_dir}/${sample_id}.sorted.bam" \
                OUTPUT="${output_dir}/${sample_id}.sorted.markdup.bam" \
                METRICS_FILE="${output_dir}/${sample_id}.markdup.picardMetrics.txt" \
                REMOVE_DUPLICATES=false \
                ASSUME_SORTED=true \
                VALIDATION_STRINGENCY=LENIENT
            
            # Add or replace read groups (commented out)
            # java -Xmx8g -jar $PICARD AddOrReplaceReadGroups \
            #     INPUT="${output_dir}/${sample_id}.sorted.markdup.bam" \
            #     OUTPUT="${output_dir}/${sample_id}.sorted.markdup.addReadGr.bam" \
            #     RGLB="${sample_id}" \
            #     RGPL=ILLUMINA \
            #     RGPU="laneUsed" \
            #     RGSM="sampleName" \
            #     RGCN="location" \
            #     RGDS="species" \
            #     VALIDATION_STRINGENCY=LENIENT
            
            # Index BAM files
            bamtools index -in "${output_dir}/${sample_id}.sorted.markdup.bam"
            
            # Generate stats on final processed bam files
            bamtools stats -in "${output_dir}/${sample_id}.sorted.markdup.bam"
            
            # Get raw transcript counts using featureCounts
            #featureCounts -T 8 --primary -p -s 0 -t gene -g gene_name \
             #   -a /projects/eversongroup/glowing/Homo_sapiens.GRCh38.92.gtf \
              #  -o "${output_dir}/${sample_id}_FC.txt" \
               # "${output_dir}/${sample_id}.sorted.markdup.bam"
        #fi
    fi
done

