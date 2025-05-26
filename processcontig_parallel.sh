#!/bin/bash

module load pythonlib samtools parallel htslib

set -ex

cd ${PBS_JOBFS}

# Sort and uniq the contigs file, remove .tmp suffix
sort "${contigs}" | uniq > "$(basename "$contigs").tmp"
mv "$(basename "$contigs").tmp" "$(basename "$contigs" .tmp)"
contigs="$(basename "$contigs" .tmp)"
chr="$(basename "$contigs" .contigs)"

# Clear existing contig files
for file in ${chr}.asm.*.contigs; do
  if [ -e "$file" ]; then
    rm "$file"
  fi
done

# Read contigs and prepare output files
while read -r line; do
  assembly=$(echo "$line" | cut -d'#' -f1-2)
  output_file="${chr}.asm.${assembly}.contigs"
  if [[ ! " ${output_files[@]} " =~ " ${output_file} " ]]; then
    output_files+=("$output_file")
  fi
  echo "$line" >> "$output_file"
done < "$contigs"

# Use GNU parallel to process the files
export metadata chr PBS_NCPUS  # Export variables to be used in parallel jobs
parallel --env metadata --env chr --env PBS_NCPUS --jobs ${PBS_NCPUS} '
  file={}
  sample=$(echo "$file" | awk -F"[.]" '"'"'{print $3}'"'"')
  fasta=$(grep "$sample" "$metadata" | awk -F "," -v donor="$sample" '"'"'$1 == donor {print $2}'"'"')
  samtools faidx ${fasta} -r ${file} -o "${chr}.asm.${sample}.fasta"
  bgzip -f -@ ${PBS_NCPUS} "${chr}.asm.${sample}.fasta"
' ::: "${output_files[@]}"

# Remove old combined fasta if exists
if [ -e "${chr}.fasta.gz" ]; then
    rm "${chr}.fasta.gz"
fi

# Concatenate all gzipped fasta files
for file in ${chr}.asm.*.fasta.gz; do
  cat ${file} >> ${chr}.fasta.gz
done

# Index the combined fasta file
samtools faidx ${chr}.fasta.gz

mv ${chr}.fasta.gz* ${outdir}
