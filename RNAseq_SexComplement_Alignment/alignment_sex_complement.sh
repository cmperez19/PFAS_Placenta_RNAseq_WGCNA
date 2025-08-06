#!/bin/bash
#SBATCH --partition=week-long-cpu
#SBATCH --time=02-00:00:00
#SBATCH --job-name="alignment_cp"
#SBATCH --mem=31GB
#SBATCH --output=run_align.out
#SBATCH --error=run_align.err
#SBATCH --mail-user=cmpere2@emory.edu
#SBATCH --mail-type=BEGIN,END


cd /projects/eversongroup/glowing/trimmed_fastq
output_dir='/projects/eversongroup/glowing/alignment_output'
# Read glowing_ID_sex.txt file and iterate through each line
while IFS=$'\t' read -r participant_id child_sex || [ -n "$participant_id" ]; do
    # Extract the sample ID from the file name
    file_name="${participant_id}_trimmed.fq.gz"
    echo "Debug: participant_id=$participant_id, child_sex=$child_sex, file_name=$file_name"
    # Check if the file exists
    if [ -e "$file_name" ]; then
        # Execute commands for the file
        echo "Processing file: $file_name"
        if [[ "$child_sex" == *"Male"* ]]; then
	    echo "Ypars:$child_sex"
            STAR --genomeDir /projects/eversongroup/glowing/YPARsMasked --sjdbGTFfile /projects/eversongroup/glowing/YPARsMasked/Homo_sapiens.GRCh38.92.gtf --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesIn "$file_name" --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix ${output_dir}/${participant_id}_STAR_M_std_YPARsMasked. --runThreadN 4
        else
           echo "Female"
           STAR --genomeDir /projects/eversongroup/glowing/Ymasked --sjdbGTFfile /projects/eversongroup/glowing/Ymasked/Homo_sapiens.GRCh38.92.gtf --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesIn "$file_name" --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix ${output_dir}/${participant_id}_STAR_F_std_Ymasked. --runThreadN 4

	   STAR --genomeDir /projects/eversongroup/glowing/YPARsMasked --sjdbGTFfile /projects/eversongroup/glowing/YPARsMasked/Homo_sapiens.GRCh38.92.gtf --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesIn "$file_name" --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix ${output_dir}/${participant_id}_STAR_F_std_YPARsMasked. --runThreadN 4
        fi
    else
       echo "File not found:$file_name"
    fi
done < glowing_id_sex.txt



