#!/bin/bash

# carrierseqxl.sh
# Angel Mojarro <mojarro at mit dot edu>
# angelmojarro.com

# Find and replace <reference_1> with your reference genome
# Find and replace <reference_2> with your reference genome

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
mkdir -p $DataFolder/06_bwa_<reference_2> # map target reads to 2nd reference genome
mkdir -p $DataFolder/07_samtools_<reference_2> # extract mapped sam file
mkdir -p $DataFolder/08_samtools_unknown_reads # extract unmapped sam file
mkdir -p $DataFolder/09_seqtk_<reference_2> # extract mapped reads
mkdir -p $DataFolder/10_seqtk_unknown_reads # extract unmapped reads
mkdir -p $DataFolder/11_final # write final "target" reads to be analyzed

# -02 bwa - index reference genome #2
Cmd="bwa index"
$Cmd $DataFolder/reference/<reference_2>/<reference_2>.fa

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

# 06 Map carrierseq_out.fastq to reference genome 2/2
Cmd="bwa mem -x ont2d -t 6" # -t, --threads check your cpu, -t 6 default on my 4 GHz Intel Core i7
$Cmd $DataFolder/reference/<reference_2>/<reference_2>.fa $DataFolder/05_target_reads/carrierseq_out.fastq > $DataFolder/06_bwa_<reference_2>/bwa_<reference_2>.sam

# 07 samtools - extract **mapped** reads as sam file
Cmd="samtools view -S -F4"
$Cmd $DataFolder/06_bwa_<reference_2>/bwa_<reference_2>.sam > $DataFolder/07_samtools_<reference_2>/bwa_<reference_2>_mapped.sam

# 07.1 samtools - identify mapped reads
cut -f1 $DataFolder/07_samtools_<reference_2>/bwa_<reference_2>_mapped.sam | sort | uniq > $DataFolder/07_samtools_<reference_2>/bwa_<reference_2>_mapped_reads.lst

# 08 samtools - extract **un-mapped** reads as sam file
Cmd="samtools view -S -f4"
$Cmd $DataFolder/06_bwa_<reference_2>/bwa_<reference_2>.sam > $DataFolder/08_samtools_unknown_reads/bwa_unknown_unmapped.sam

# 08.1 samtools - identify unmapped reads
cut -f1 $DataFolder/08_samtools_unknown_reads/bwa_unknown_unmapped.sam | sort | uniq > $DataFolder/08_samtools_unknown_reads/bwa_unknown_unmapped_reads.lst

# 9 seqtk - extract **mapped** reads as fastq file
Cmd="seqtk subseq"
$Cmd $DataFolder/05_target_reads/carrierseq_out.fastq $DataFolder/07_samtools_<reference_2>/bwa_<reference_2>_mapped_reads.lst > $DataFolder/09_seqtk_<reference_2>/<reference_2>_reads.fastq

# 9.1 seqtk - make fasta file from mapped fastq
Cmd="seqtk seq -a"
$Cmd $DataFolder/09_seqtk_<reference_2>/<reference_2>_reads.fastq > $DataFolder/09_seqtk_<reference_2>/<reference_2>_reads.fasta

# 9.2 grep - count mapped reads from fasta file
grep -c ">" $DataFolder/09_seqtk_<reference_2>/<reference_2>_reads.fasta > $DataFolder/09_seqtk_<reference_2>/<reference_2>_reads.txt

# 10 seqtk - extract **un-mapped** reads as fastq file
Cmd="seqtk subseq"
$Cmd $DataFolder/05_target_reads/carrierseq_out.fastq $DataFolder/08_samtools_unknown_reads/bwa_unknown_unmapped_reads.lst > $DataFolder/10_seqtk_unknown_reads/unknown_reads.fastq

# 10.1 seqtk - make fasta file from unmapped fastq
Cmd="seqtk seq -a"
$Cmd $DataFolder/10_seqtk_unknown_reads/unknown_reads.fastq > $DataFolder/10_seqtk_unknown_reads/unknown_reads.fasta

# 10.2 grep - count unmapped reads from fasta file
grep -c ">" $DataFolder/10_seqtk_unknown_reads/unknown_reads.fasta > $DataFolder/10_seqtk_unknown_reads/unknown_reads.txt

# 11 copy <reference_2> reads and unknown reads to 10_final folder for further analysis
cp $DataFolder/09_seqtk_<reference_2>/<reference_2>_reads.fastq $DataFolder/11_final/carrierseq_<reference_2>.fastq
cp $DataFolder/10_seqtk_unknown_reads/unknown_reads.fastq $DataFolder/11_final/carrierseq_unknown.fastq

# 11.1 Final seqtk - make fasta file from fastq
Cmd="seqtk seq -a"
$Cmd $DataFolder/11_final/carrierseq_<reference_2>.fastq > $DataFolder/11_final/carrierseq_<reference_2>.fasta
$Cmd $DataFolder/11_final/carrierseq_unknown.fastq > $DataFolder/11_final/carrierseq_unknown.fasta

