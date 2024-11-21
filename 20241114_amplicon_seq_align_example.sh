#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=align
#SBATCH --ntasks=1 --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1:00:00
#SBATCH --error=align.%A_%a.err
#SBATCH --output=align.%A_%a.out
#SBATCH --mail-type=ALL

module load SAMtools/1.21-GCC-12.2.0
module load BWA/0.7.17-GCCcore-12.2.0

# Define file paths (update these paths as needed)
REFERENCE="/home/df555/palmer_scratch/amplicon_seq/reference_Cell_SEC62.fasta"   # Path to the reference FASTA file, generated with R code
FASTQ_R1="/home/df555/palmer_scratch/amplicon_seq/U2-SEC62_S1_L001_R1_001.fastq.gz"  # Path to the first FASTQ file
FASTQ_R2="/home/df555/palmer_scratch/amplicon_seq/U2-SEC62_S1_L001_R2_001.fastq.gz"  # Path to the second FASTQ file
OUTPUT_BAM="/home/df555/palmer_scratch/amplicon_seq/IVT_2A.bam"                 # Output BAM file - just an intermediate file
SORTED_BAM="/home/df555/palmer_scratch/amplicon_seq/sorted_output2A.bam"              #Sorted BAM file - this is what you will put back into the R code 

# Check if bwa is installed
if ! command -v bwa &> /dev/null
then
    echo "bwa could not be found. Please install bwa."
    exit 1
fi

# Check if samtools is installed
if ! command -v samtools &> /dev/null
then
    echo "samtools could not be found. Please install samtools."
    exit 1
fi

# Step 1: Index the reference sequence (only needs to be done once)
echo "Indexing reference sequence..."
bwa index "$REFERENCE"

# Step 2: Align FASTQ files to the reference and convert to BAM format
echo "Aligning FASTQ files to reference..."
bwa mem "$REFERENCE" "$FASTQ_R1" "$FASTQ_R2" | samtools view -b -o "$OUTPUT_BAM"

# Step 3: Sort the BAM file
echo "Sorting BAM file..."
samtools sort "$OUTPUT_BAM" -o "$SORTED_BAM"

# Step 4: Index the sorted BAM file
echo "Indexing sorted BAM file..."
samtools index "$SORTED_BAM"

echo "Alignment and processing complete. Output BAM: $SORTED_BAM"