#!/bin/bash
#SBATCH --partition=week-long-cpu
#SBATCH --time=03-00:00:00
#SBATCH --job-name="trim_fastq"
#SBATCH --mem=31GB
#SBATCH --output=run_trim.out
#SBATCH --error=run_trim.err
#SBATCH --mail-user=cmpere2@emory.edu
#SBATCH --mail-type=BEGIN,END


# Directory containing FASTQ files
input_dir="fastq_files" 

# Directory to store trimmed FASTQ files
output_dir="trimmed_fastq"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each FASTQ file in the input directory
for file in "$input_dir"/*.fastq.gz; do
    if [ -e "$file" ]; then
        # Extract filename without extension
        filename=$(basename "$file" .fastq.gz)

        # Extract sample ID
        sample_id=$(echo "$filename" | cut -d'_' -f1)

        # Trim using Trim Galore!
        echo "Trimming $filename..."
        trim_galore --quality 30 --phred33 --stringency 5 --length 20 "$file" -o "$output_dir"

        # Renaming output files to reflect sample ID
        mv "$output_dir/${filename}_trimmed.fq.gz" "$output_dir/${sample_id}_trimmed.fq.gz"

        # FastQC on trimmed files
        echo "Running FastQC on trimmed $sample_id..."
        fastqc "$output_dir/${sample_id}_trimmed.fq.gz" -o "$output_dir"
    fi
done

echo "Processing complete!"

