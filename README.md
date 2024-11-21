Most amplicon sequencing returns fastq.gz raw files. Steps for analysis will be:
1) Generate a .fasta file with your reference sequence of interest with the .Rmd notebook for amplicon sequencing alignment. This reference will be used to turn fastq.gz files into aligned bam files - it is important this is accurate in order for bwa to call mutations/insertions/deletions etc.
2) Align your fastq.gz data using bwa with the shell script provided, replacing file names as eneded. 
3) Return to .Rmd notebook for amplicon sequencing alignment, replacing file names and sequences used for filtering as needed.
4) Run the amplicon sequencing analysis code using the final filtered table generated from the alignment code, or with the abundance table provided by Azenta directly.
