#!/bin/bash

# Vector-Genome Fusion Reads Analysis Pipeline
# This script automates the three rounds of mapping to identify vector-genome fusion reads

# Exit on error
set -e

# Function to check if a command exists
check_command() {
    if ! command -v "$1" &> /dev/null; then
        echo "Error: $1 is required but not installed."
        exit 1
    fi
}

# Define paths
BOWTIE="/home/wz1424/.conda/envs/piRNA_fusion/bin/bowtie"
SAMTOOLS="/home/wz1424/.conda/envs/piRNA_fusion/bin/samtools"

# Check if executables exist at specified paths
if [ ! -x "$BOWTIE" ]; then
    echo "Error: bowtie not found at $BOWTIE or not executable"
    exit 1
fi

if [ ! -x "$SAMTOOLS" ]; then
    echo "Error: samtools not found at $SAMTOOLS or not executable"
    exit 1
fi

# Check required commands
check_command "python"

# Input validation: now require 5 parameters
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <input_fastq> <vector_index> <genome_index> <output_prefix> <extract_bp>"
    echo "Example: $0 input.fastq Original_vector dm6 output 20"
    exit 1
fi

# Input parameters
INPUT_FASTQ=$1
VECTOR_INDEX=$2
GENOME_INDEX=$3
OUTPUT_PREFIX=$4
EXTRACT_BP=$5

# Create output directory
mkdir -p "${OUTPUT_PREFIX}_output"
cd "${OUTPUT_PREFIX}_output"

echo "=== Starting Vector-Genome Fusion Reads Analysis ==="

# Round 1: Get reads containing vector sequences
echo "Round 1: Identifying reads with vector sequences..."

# Extract EXTRACT_BP bp from 5' end
echo "Extracting ${EXTRACT_BP}bp from 5' end..."
python ~/piRNA_fusion_read/scripts/trimfastq_python3.py "$INPUT_FASTQ" "$EXTRACT_BP" -stdout > "${EXTRACT_BP}fiveprime.fastq"

# Map to vector and extract mapped reads
echo "Mapping to vector..."
"$BOWTIE" -x "$VECTOR_INDEX" -p 8 -v 2 -k 1 -t -q \
    --al alltrimmedVectorfastq "${EXTRACT_BP}fiveprime.fastq" -S "${EXTRACT_BP}fiveprime.sam"

# Retrieve original full length reads
echo "Retrieving full length reads..."
grep -f <(cat alltrimmedVectorfastq | paste - - - - | cut -f 1) \
    <(cat "$INPUT_FASTQ" | sed 's/ /_/g' | paste - - - - ) | \
    tr "\t" "\n" > alltrimmedVectorFullLengthfastq

# Round 2: Filter out pure vector reads and vector-vector reads
echo "Round 2: Filtering pure vector reads and vector-vector reads..."

# Extract EXTRACT_BP bp from 3' end
echo "Extracting ${EXTRACT_BP}bp from 3' end..."
paste <(cat alltrimmedVectorFullLengthfastq | paste - - - - | cut -f1 | tr "\t" "\n") \
      <(cat alltrimmedVectorFullLengthfastq | paste - - - - | cut -f2 | tr "\t" "\n" | grep -o ".\{${EXTRACT_BP}\}$") \
      <(cat alltrimmedVectorFullLengthfastq | paste - - - - | cut -f3 | tr "\t" "\n") \
      <(cat alltrimmedVectorFullLengthfastq | paste - - - - | cut -f4 | tr "\t" "\n" | grep -o ".\{${EXTRACT_BP}\}$") | tr "\t" "\n" > "${EXTRACT_BP}threeprime.fastq"

# Map to vector index and get unmapped reads
"$BOWTIE" -x "$VECTOR_INDEX" -p 8 -v 2 -k 1 -t -q \
    --un 3primeVectorUnmapfastq "${EXTRACT_BP}threeprime.fastq" -S /dev/null

# Round 3: Identify genome sequences
echo "Round 3: Identifying genome sequences..."

# Map to genome and create sorted BAM
echo "Mapping to genome and creating BAM file..."
"$BOWTIE" -x "$GENOME_INDEX" -p 8 -v 0 -k 1 -t -q 3primeVectorUnmapfastq -S 3primeVectorUnmap.sam

# Identify mapped transposon vector and position
gawk '($1 !~ /^@/ && $4 > 0) {print $1 "\t" $3}' "${EXTRACT_BP}fiveprime.sam" > "${EXTRACT_BP}fiveprimeMap.txt"

# Link the mapping vector position to genome position
awk -v mapfile="${EXTRACT_BP}fiveprimeMap.txt" 'BEGIN {
    while(getline < mapfile) {
        vec[$1] = $2
    }
}
/^@/ { print; next }
{
    if($1 in vec) {
        $1 = $1 ":" vec[$1]
    }
    print
}' 3primeVectorUnmap.sam > 3primeVectorUnmapAnnot.sam


"$SAMTOOLS" view -bT "$GENOME_INDEX.fa" 3primeVectorUnmapAnnot.sam | "$SAMTOOLS" sort -o "${OUTPUT_PREFIX}.3primeVectorUnmapAnnot.dm6.bam"

# Index BAM file
echo "Indexing BAM file..."
"$SAMTOOLS" index "${OUTPUT_PREFIX}.3primeVectorUnmapAnnot.dm6.bam"

# Identify mapped transposon vector and position
gawk '($1 !~ /^@/ && $4 > 0) {print $1 "\t" $3 "\t" $4}' "${EXTRACT_BP}fiveprime.sam" > "${EXTRACT_BP}fiveprimeMap.txt"
gawk '($1 !~ /^@/ && $4 > 0) {print $1}' 3primeVectorUnmap.sam > 3primeVectorUnmapMap.txt
sort -k1,1 "${EXTRACT_BP}fiveprimeMap.txt" > "${EXTRACT_BP}fiveprimeMap.sorted.txt"
sort -k1,1 3primeVectorUnmapMap.txt > 3primeVectorUnmapMap.sorted.txt
join -t $'\t' -1 1 -2 1 -o 1.1,1.2,1.3 "${EXTRACT_BP}fiveprimeMap.sorted.txt" 3primeVectorUnmapMap.sorted.txt > "${OUTPUT_PREFIX}.vectorPos.txt"

# Generate summary statistics
echo "Generating summary statistics..."
echo "Initial reads: $(expr $(cat "$INPUT_FASTQ" | wc -l) / 4)" > "${OUTPUT_PREFIX}_statistics.txt"
echo "Reads with vector (5'): $(expr $(cat alltrimmedVectorfastq | wc -l) / 4)" >> "${OUTPUT_PREFIX}_statistics.txt"
echo "Reads after vector filtering: $(expr $(cat 3primeVectorUnmapfastq | wc -l) / 4)" >> "${OUTPUT_PREFIX}_statistics.txt"
echo "Final fusion reads: $("$SAMTOOLS" view -c "${OUTPUT_PREFIX}.3primeVectorUnmapAnnot.dm6.bam")" >> "${OUTPUT_PREFIX}_statistics.txt"

echo "=== Analysis Complete ==="
echo "Results can be found in ${OUTPUT_PREFIX}_output directory"
echo "BAM file: ${OUTPUT_PREFIX}.3primeVectorUnmapAnnot.dm6.bam"
echo "Statistics file: ${OUTPUT_PREFIX}_statistics.txt"
echo "Vector position file: ${OUTPUT_PREFIX}.vectorPos.txt"
