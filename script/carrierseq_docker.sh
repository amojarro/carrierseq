#!/bin/bash

# carrierseq.sh
# Angel Mojarro <mojarro at mit dot edu>
# angelmojarro.com

# Usage info
show_help() {
cat << EOF
Usage: ${0##*/} [-i INPUT] [-r REFERENCE] [-o OUTPUT]...
CarrierSeq requires bwa, samtools, seqtk, and fqtrim. 
Reads to be analyzed must be compiled into a single fastq file and the reference genome must be in fasta format.
     -i          All reads to be analyzed *.fastq
     -r          Carrier reference genome *.fasta
     -t          Number of threads used for BWA mapping (default = 1)
     -q          User-defined quality (phred) score (default = 9)
     -p          User-defined p-value 
                 (default = 0.0001 or 0.05/512 active channels)
     -o          Output directory
EOF
}

# Getops variables
all_reads=""
reference_genome=""
bwa_threads="1"
q_score="9"
p_value="0.0001"
output_folder=""

OPTIND=1   

while getopts "h?i:r:t:q:p:o:v:" opt; do
    case "$opt" in
    h)
        show_help
        exit 0
        ;;
    i)  all_reads=$OPTARG
        ;;
    r)  reference_genome=$OPTARG
        ;;
    t)  bwa_threads=$OPTARG
        ;;
    q)  q_score=$OPTARG
        ;;
    p)  p_value=$OPTARG
        ;;
    o)  output_folder=$OPTARG
        ;;       
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

echo "all_reads='$all_reads', reference_genome='$reference_genome', bwa_threads='$bwa_threads', q_score='$q_score', p_value='$p_value', output_folder='$output_folder'"

# Make output directories
mkdir -p $output_folder/fastq_tmp # Temporarily copies reads into docker working folder
mkdir -p $output_folder/reference_tmp # Temporarily copies reference genome into docker working folder
mkdir -p $output_folder/00_bwa # map all reads to carrier reference genome
mkdir -p $output_folder/01_samtools # extract unmapped sam file
mkdir -p $output_folder/02_seqtk # extract unmapped reads
mkdir -p $output_folder/03_fastqc # discard low-quality reads < qx
mkdir -p $output_folder/03_01_low_quality_reads # save low-quality reads
mkdir -p $output_folder/04_fqtrim_dusted # discard low complexity reads
mkdir -p $output_folder/04_01_low_complexity_reads # save low-complexity reads
mkdir -p $output_folder/05_reads_of_interest # 
mkdir -p $output_folder/06_poisson_calculation # calculations for sorting "real" reads versus "possible noise"
mkdir -p $output_folder/07_hqnrs #
mkdir -p $output_folder/08_target_reads # final output reads to be analyzed if target is unknown

# Docker settings
CarrierSeq="mojarro/carrierseq:latest" 
DockerPath="/carrierseq"
DockerOptions="-v $output_folder:$DockerPath $CarrierSeq"

# Copy reads and reference genome into temporary folders for Docker
echo Copying $all_reads & $reference_genome into temporary docker container
cp $all_reads $output_folder/fastq_tmp/all_reads.fastq
cp $reference_genome $output_folder/reference_tmp/reference_genome.fasta

# -01 bwa - Index carrier reference genome
echo Indexing reference genome...
Cmd="$DockerOptions bwa index"
docker run $Cmd $DockerPath/reference_tmp/reference_genome.fasta

# 00 bwa - map $all_reads to the $reference_genome
echo Mapping all reads to reference genome...
Cmd="$DockerOptions bwa mem -x ont2d -t $bwa_threads"
docker run $Cmd $DockerPath/reference_tmp/reference_genome.fasta $DockerPath/fastq_tmp/all_reads.fastq > $output_folder/00_bwa/bwa_mapped.sam

# 01 samtools - extract unmapped reads as sam file
echo Extracting unmapped reads...
Cmd="$DockerOptions samtools view -S -f4"
docker run $Cmd $DockerPath/00_bwa/bwa_mapped.sam > $output_folder/01_samtools/bwa_unmapped.sam

# 01.1 samtools - identify unmapped reads
cut -f1 $output_folder/01_samtools/bwa_unmapped.sam | sort | uniq > $output_folder/01_samtools/bwa_unmapped_reads.lst

# 02 seqtk - extract unmapped reads as fastq file
Cmd="$DockerOptions seqtk subseq"
docker run $Cmd $DockerPath/fastq_tmp/all_reads.fastq $DockerPath/01_samtools/bwa_unmapped_reads.lst > $output_folder/02_seqtk/unmapped_reads.fastq

# 02.1 seqtk - make fasta file
Cmd="$DockerOptions seqtk seq -a"
docker run $Cmd $DockerPath/02_seqtk/unmapped_reads.fastq > $output_folder/02_seqtk/unmapped_reads.fasta

# 02.2 grep - count reads from fasta file
grep -c ">" $output_folder/02_seqtk/unmapped_reads.fasta > $output_folder/02_seqtk/unmapped_reads.txt
echo Total unmapped reads:
cat $output_folder/02_seqtk/unmapped_reads.txt

# 03 quality_score_filter.py - discard low-quality reads
echo Applying quality filter...
python python/quality_score_filter.py $output_folder/02_seqtk/unmapped_reads.fastq $output_folder/03_fastqc/unmapped_reads_qc $q_score

# 03.1 grep - count "high-quality" reads
grep -c ">" $output_folder/03_fastqc/unmapped_reads_qc.fa > $output_folder/03_fastqc/unmapped_reads_qc.txt
echo Reads ≥ $q_score quality score:
cat $output_folder/03_fastqc/unmapped_reads_qc.txt

# 03.01 grep - identify discarded low-quality reads
grep -e ">" $output_folder/03_fastqc/unmapped_reads_qc.fa | awk '{print $1}' | sed 's/^.//' > $output_folder/03_fastqc/unmapped_reads_qc.lst
grep -Fxvf $output_folder/03_fastqc/unmapped_reads_qc.lst $output_folder/01_samtools/bwa_unmapped_reads.lst > $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.lst

# 03.01.1 seqtk - save discarded low-quality read
Cmd="$DockerOptions seqtk subseq"
docker run $Cmd $DockerPath/02_seqtk/unmapped_reads.fastq $DockerPath/03_01_low_quality_reads/low_quality_unmapped_reads.lst > $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.fastq

# 03.01.2 seqtk - make fasta file
Cmd="$DockerOptions seqtk seq -a"
docker run $Cmd $DockerPath/03_01_low_quality_reads/low_quality_unmapped_reads.fastq > $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.fasta

# 03.01.3 grep - count low-quality reads from fasta file
grep -c ">" $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.fasta > $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.txt

# 04 fqtrim - discard low complexity reads
echo Applying DUST filter...
Cmd="$DockerOptions fqtrim -D"
docker run $Cmd $DockerPath/03_fastqc/unmapped_reads_qc.fq > $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.fastq

# 04.1 seqtk - make fasta file
Cmd="$DockerOptions seqtk seq -a"
docker run $Cmd $DockerPath/04_fqtrim_dusted/unmapped_reads_qc_dusted.fastq > $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.fasta

# 04.2 grep - count dusted reads from fasta file
grep -c ">" $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.fasta > $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.txt
echo Reads after DUST filter:
cat $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.txt

# 04.01 - identify discarded low-complexity reads
grep -e ">" $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.fasta | awk '{print $1}' | sed 's/^.//' > $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.lst
grep -Fxvf $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.lst $output_folder/03_fastqc/unmapped_reads_qc.lst > $output_folder/04_01_low_complexity_reads/low_complexity_reads_qc.lst

# 04.01.1 seqtk - save discarded low-complexity reads
Cmd="$DockerOptions seqtk subseq"
docker run $Cmd $DockerPath/03_fastqc/unmapped_reads_qc.fq $DockerPath/04_01_low_complexity_reads/low_complexity_reads_qc.lst > $output_folder/04_01_low_complexity_reads/low_complexity_reads_qc.fastq

# 04.01.2 seqtk - make fasta file
Cmd="$DockerOptions seqtk seq -a"
docker run $Cmd $DockerPath/04_01_low_complexity_reads/low_complexity_reads_qc.fastq > $output_folder/04_01_low_complexity_reads/low_complexity_reads_qc.fasta

# 04.01.3 grep - count low-complexity reads from fasta file
grep -c ">" $output_folder/04_01_low_complexity_reads/low_complexity_reads_qc.fasta > $output_folder/04_01_low_complexity_reads/low_complexity_reads_qc.txt

# 05 copy "high-quality reads" qc and higher + complex reads to 05_target_reads for further analysis
echo Saving Reads of Interest...
cp $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.fastq $output_folder/05_reads_of_interest/carrierseq_roi.fastq
cp $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.fasta $output_folder/05_reads_of_interest/carrierseq_roi.fasta 
echo Done!

# 05.1 grep - count target reads from fasta file
grep -c ">" $output_folder/05_reads_of_interest/carrierseq_roi.fasta > $output_folder/05_reads_of_interest/carrierseq_roi.txt
echo Starting Poisson Calculation...

ChannelsInUse="$output_folder/06_poisson_calculation/channels_in_use.txt"
TotalROIs="$output_folder/05_reads_of_interest/carrierseq_roi.txt"
LambdaValue="$output_folder/06_poisson_calculation/lambda_value.txt"

##### FOR OLD FASTQ HEADER #####
# 06 grep - extract all channels used, delete duplicates to count unique (n/512) channels used
grep -Eio "_ch[0-9]+_" $all_reads | awk '!seen[$0]++' > $output_folder/06_poisson_calculation/channels_used.lst

##### FOR NEW FASTQ HEADER #####
# grep -Eio "ch=[0-9]+" $all_reads | awk '!seen[$0]++' > $output_folder/06_poisson_calculation/channels_used.lst

# 06.01 - count unique channels (n/512)
grep -c "ch" $output_folder/06_poisson_calculation/channels_used.lst > $output_folder/06_poisson_calculation/channels_in_use.txt
echo Channels in use:
cat $output_folder/06_poisson_calculation/channels_in_use.txt

# 06.02 python - calculate lambda for poisson calculation
echo Calculating lambda value and x_crit...
python python/calculate_lambda.py $TotalROIs $ChannelsInUse > $output_folder/06_poisson_calculation/lambda_value.txt

# 06.02.1 python - calculate x_critical
python python/xcrit.py $LambdaValue $p_value > $output_folder/06_poisson_calculation/read_channel_threshold.txt
cat $output_folder/06_poisson_calculation/read_channel_threshold.txt

# Cleaning up
echo Deleting $all_reads & $reference_genome temporary docker files
rm -r $output_folder/fastq_tmp
rm -r $output_folder/reference_tmp

# End of file
