#!/bin/bash

# carrierseq.sh
# Angel Mojarro <mojarro at mit dot edu>
# angelmojarro.com

# Find and replace <reference_1> with your reference genome

# your working directory
DataFolder="your/working/directory" # the DataFolder contains the fastq, python, and reference folders containing your all_reads.fast file, the included quality filter script, and *.fa reference genomes
FastQ="$DataFolder/fastq" # fastq to be analyzed - "all_reads.fastq"
Reference="$DataFolder/reference" # reference genome(s). example - "lambda/lambda.fa"
PythonScript="$DataFolder/python" # fastq quality filter python script by Michael Micorescu <Michael dot Micorescu at nanoporetech dot com> and edited by Angel Mojarro

# make output directories
mkdir -p $DataFolder/00_bwa # map all reads to carrier reference genome
mkdir -p $DataFolder/01_samtools # extract unmapped sam file
mkdir -p $DataFolder/02_seqtk # extract unmapped reads
mkdir -p $DataFolder/03_fastq9 # discard low-quality reads < q9
mkdir -p $DataFolder/03_01_low_quality_reads # save low-quality reads
mkdir -p $DataFolder/04_fqtrim_dusted # discard low complexity reads
mkdir -p $DataFolder/04_01_low_complexity_reads # save low-complexity reads
mkdir -p $DataFolder/05_target_reads # final output reads to be analyzed if target is unknown
mkdir -p $DataFolder/06_poisson_calculation # calculations for sorting "real" reads versus "possible noise"

# -01 bwa - index reference genome #1
Cmd="bwa index"
$Cmd $DataFolder/reference/<reference_1>/<reference_1>.fa

# 00 bwa - map all_reads.fastq the reference genome 1/2
Cmd="bwa mem -x ont2d -t 6" # -t, --threads check your cpu, -t 6 default on my 4 GHz Intel Core i7
$Cmd $DataFolder/reference/<reference_1>/<reference_1>.fa $DataFolder/fastq/all_reads.fastq > $DataFolder/00_bwa/bwa_mapped.sam
 
# 01 samtools - extract unmapped reads as sam file
Cmd="samtools view -S -f4"
$Cmd $DataFolder/00_bwa/bwa_mapped.sam > $DataFolder/01_samtools/bwa_unmapped.sam

# 01.1 samtools - identify unmapped reads
cut -f1 $DataFolder/01_samtools/bwa_unmapped.sam | sort | uniq > $DataFolder/01_samtools/bwa_unmapped_reads.lst

# 02 seqtk - extract unmapped reads as fastq file
Cmd="seqtk subseq"
$Cmd $DataFolder/fastq/all_reads.fastq $DataFolder/01_samtools/bwa_unmapped_reads.lst > $DataFolder/02_seqtk/unmapped_reads.fastq

# 02.1 seqtk - make fasta file
Cmd="seqtk seq -a"
$Cmd $DataFolder/02_seqtk/unmapped_reads.fastq > $DataFolder/02_seqtk/unmapped_reads.fasta

# 02.2 grep - count reads from fasta file
grep -c ">" $DataFolder/02_seqtk/unmapped_reads.fasta > $DataFolder/02_seqtk/unmapped_reads.txt

# 03 MM_filter_ont_1.py - discard low-quality reads
python $DataFolder/python/MM_filter_ont_1_AM.py

# 03.1 grep - count "high-qualit" reads
grep -c ">" $DataFolder/03_fastq9/unmapped_reads_q9.fa > $DataFolder/03_fastq9/unmapped_reads_q9.txt

# 03.01 grep - identify discarded low-quality reads
grep -e ">" $DataFolder/03_fastq9/unmapped_reads_q9.fa | awk '{print $1}' | sed 's/^.//' > $DataFolder/03_fastq9/unmapped_reads_q9.lst
grep -Fxvf $DataFolder/03_fastq9/unmapped_reads_q9.lst $DataFolder/01_samtools/bwa_unmapped_reads.lst > $DataFolder/03_01_low_quality_reads/low_quality_unmapped_reads.lst

# 03.01.1 seqtk - save discarded low-quality read
Cmd="seqtk subseq"
$Cmd $DataFolder/02_seqtk/unmapped_reads.fastq $DataFolder/03_01_low_quality_reads/low_quality_unmapped_reads.lst > $DataFolder/03_01_low_quality_reads/low_quality_unmapped_reads.fastq

# 03.01.2 seqtk - make fasta file
Cmd="seqtk seq -a"
$Cmd $DataFolder/03_01_low_quality_reads/low_quality_unmapped_reads.fastq > $DataFolder/03_01_low_quality_reads/low_quality_unmapped_reads.fasta

# 03.01.3 grep - count low-quality reads from fasta file
grep -c ">" $DataFolder/03_01_low_quality_reads/low_quality_unmapped_reads.fasta > $DataFolder/03_01_low_quality_reads/low_quality_unmapped_reads.txt

# 04 fqtrim - discard low complexity reads
Cmd="fqtrim -D"
$Cmd $DataFolder/03_fastq9/unmapped_reads_q9.fq > $DataFolder/04_fqtrim_dusted/unmapped_reads_q9_dusted.fastq

# 04.1 seqtk - make fasta file
Cmd="seqtk seq -a"
$Cmd $DataFolder/04_fqtrim_dusted/unmapped_reads_q9_dusted.fastq > $DataFolder/04_fqtrim_dusted/unmapped_reads_q9_dusted.fasta

# 04.2 grep - count dusted reads from fasta file
grep -c ">" $DataFolder/04_fqtrim_dusted/unmapped_reads_q9_dusted.fasta > $DataFolder/04_fqtrim_dusted/unmapped_reads_q9_dusted.txt

# 04.01 - identify discarded low-complexity reads
grep -e ">" $DataFolder/04_fqtrim_dusted/unmapped_reads_q9_dusted.fasta | awk '{print $1}' | sed 's/^.//' > $DataFolder/04_fqtrim_dusted/unmapped_reads_q9_dusted.lst
grep -Fxvf $DataFolder/04_fqtrim_dusted/unmapped_reads_q9_dusted.lst $DataFolder/03_fastq9/unmapped_reads_q9.lst > $DataFolder/04_01_low_complexity_reads/low_complexity_reads_q9.lst

# 04.01.1 seqtk - save discarded low-complexity reads
Cmd="seqtk subseq"
$Cmd $DataFolder/03_fastq9/unmapped_reads_q9.fq $DataFolder/04_01_low_complexity_reads/low_complexity_reads_q9.lst > $DataFolder/04_01_low_complexity_reads/low_complexity_reads_q9.fastq

# 04.01.2 seqtk - make fasta file
Cmd="seqtk seq -a"
$Cmd $DataFolder/04_01_low_complexity_reads/low_complexity_reads_q9.fastq > $DataFolder/04_01_low_complexity_reads/low_complexity_reads_q9.fasta

# 04.01.3 grep - count low-complexity reads from fasta file
grep -c ">" $DataFolder/04_01_low_complexity_reads/low_complexity_reads_q9.fasta > $DataFolder/04_01_low_complexity_reads/low_complexity_reads_q9.txt

# 05 copy "high-quality reads" q9 and higher + complex reads to 05_target_reads for further analysis
cp $DataFolder/04_fqtrim_dusted/unmapped_reads_q9_dusted.fastq $DataFolder/05_target_reads/carrierseq_out.fastq
cp $DataFolder/04_fqtrim_dusted/unmapped_reads_q9_dusted.fasta $DataFolder/05_target_reads/carrierseq_out.fasta 

# 05.1 grep - count target reads from fasta file
grep -c ">" $DataFolder/05_target_reads/carrierseq_out.fasta > $DataFolder/05_target_reads/carrierseq_out.txt

# 06 grep - extract all channels used, delete duplicates to count unique (n/512) channels used
grep -Eio "_ch[0-9]+_" $DataFolder/fastq/all_reads.fastq | awk '!seen[$0]++' > $DataFolder/06_poisson_calculation/channels_used.lst

# 06.01 - count unique channels (n/512)
grep -c "ch" $DataFolder/06_poisson_calculation/channels_used.lst > $DataFolder/06_poisson_calculation/channels_in_use.txt

# 06.02 python - calculate lambda for poisson calculation
python $DataFolder/python/calculate_lambda.py > $DataFolder/06_poisson_calculation/lambda_value.txt

# 06.02.1 python - calculate x_critical
python $DataFolder/python/xcrit.py > $DataFolder/06_poisson_calculation/read_channel_threshold.txt