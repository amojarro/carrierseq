#!/bin/bash

# carrierseq.sh
# Angel Mojarro <mojarro at mit dot edu>
# angelmojarro.com

# Usage info -h
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

# Getops error checking -i, -r, and -o are required
if ( ! getopts "i:r:o:" opt); then
	echo "Usage: `basename $0` options [-i INPUT] [-r REFERENCE] [-o OUTPUT] are required. Use -h for help";
	exit $E_OPTERROR;
fi

OPTIND=1   

while getopts "h?i:r:t:q:p:o:" opt; do
    case "$opt" in
    h)
        show_help
        exit 0
        ;;
    i)  all_reads=$OPTARG
        ;;
    \?)
        exit 1
        ;;    
    r)  reference_genome=$OPTARG
        ;;
    \?)
        exit 1
        ;;    
    t)  bwa_threads=$OPTARG
        ;;
    q)  q_score=$OPTARG
        ;;
    p)  p_value=$OPTARG
        ;;
    o)  output_folder=$OPTARG
        ;;
    \?)
        exit 1
        ;;           
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

echo "all_reads='$all_reads', reference_genome='$reference_genome', bwa_threads='$bwa_threads', q_score='$q_score', p_value='$p_value', output_folder='$output_folder'"

# Make output directories
echo Creating output directories...
mkdir -p $output_folder/fastq_tmp # Temporarily copies reads into docker working folder
mkdir -p $output_folder/reference_tmp # Temporarily copies reference genome into docker working folder
mkdir -p $output_folder/00_bwa # map all reads to carrier reference genome
mkdir -p $output_folder/01_samtools # extract unmapped sam file
mkdir -p $output_folder/02_seqtk # extract unmapped reads
mkdir -p $output_folder/03_fastqc # discard low-quality reads < qx
mkdir -p $output_folder/03_01_low_quality_reads # save low-quality reads
mkdir -p $output_folder/04_fqtrim_dusted # discard low complexity reads
mkdir -p $output_folder/04_01_low_complexity_reads # save low-complexity reads
mkdir -p $output_folder/05_reads_of_interest # filtered reads to poisson calculation
mkdir -p $output_folder/06_poisson_calculation # calculations for sorting "real" reads versus "possible noise"
mkdir -p $output_folder/07_hqnrs # "high-quality noise reads" 
mkdir -p $output_folder/08_target_reads # final output reads to be analyzed if target is unknown

##################
#                #
#    Mapping     #
#                #
##################

# Docker settings
CarrierSeq="mojarro/carrierseq:latest" 
DockerPath="/carrierseq"
DockerOptions="-v $output_folder:$DockerPath $CarrierSeq"

# Copy reads and reference genome into temporary folders for Docker
echo Copying reads and reference genome into temporary docker file...
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
echo Unmapped reads saved to 02_seqtk!

# 02.2 grep - count reads from fasta file
grep -c ">" $output_folder/02_seqtk/unmapped_reads.fasta > $output_folder/02_seqtk/unmapped_reads.txt
echo Total unmapped reads:
cat $output_folder/02_seqtk/unmapped_reads.txt

# 03 quality_score_filter.py - discard low-quality reads
echo Applying quality filter...
python python/quality_score_filter.py $output_folder/02_seqtk/unmapped_reads.fastq $output_folder/03_fastqc/unmapped_reads_qc $q_score
echo Reads saved to 03_fastqc!

# 03.1 grep - count "high-quality" reads
grep -c ">" $output_folder/03_fastqc/unmapped_reads_qc.fa > $output_folder/03_fastqc/unmapped_reads_qc.txt
echo Reads ≥ $q_score quality score:
cat $output_folder/03_fastqc/unmapped_reads_qc.txt

# 03.01 grep - identify discarded low-quality reads
echo Saving low-quality reads...
grep -e ">" $output_folder/03_fastqc/unmapped_reads_qc.fa | awk '{print $1}' | sed 's/^.//' > $output_folder/03_fastqc/unmapped_reads_qc.lst
grep -Fxvf $output_folder/03_fastqc/unmapped_reads_qc.lst $output_folder/01_samtools/bwa_unmapped_reads.lst > $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.lst
echo Reads saved to 03_01_low_quality_reads!

# 03.01.1 seqtk - save discarded low-quality read
Cmd="$DockerOptions seqtk subseq"
docker run $Cmd $DockerPath/02_seqtk/unmapped_reads.fastq $DockerPath/03_01_low_quality_reads/low_quality_unmapped_reads.lst > $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.fastq

# 03.01.2 seqtk - make fasta file
Cmd="$DockerOptions seqtk seq -a"
docker run $Cmd $DockerPath/03_01_low_quality_reads/low_quality_unmapped_reads.fastq > $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.fasta

# 03.01.3 grep - count low-quality reads from fasta file
grep -c ">" $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.fasta > $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.txt
echo Total discarded low-quality reads:
cat $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.txt

# 04 fqtrim - discard low complexity reads
echo Applying DUST filter...
Cmd="$DockerOptions fqtrim -D"
docker run $Cmd $DockerPath/03_fastqc/unmapped_reads_qc.fq > $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.fastq

# 04.1 seqtk - make fasta file
Cmd="$DockerOptions seqtk seq -a"
docker run $Cmd $DockerPath/04_fqtrim_dusted/unmapped_reads_qc_dusted.fastq > $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.fasta
echo Reads saved to 04_fqtrim_dusted!

# 04.2 grep - count dusted reads from fasta file
grep -c ">" $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.fasta > $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.txt
echo Reads after DUST filter:
cat $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.txt

# 04.01 - identify discarded low-complexity reads
echo Saving low-complexity reads...
grep -e ">" $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.fasta | awk '{print $1}' | sed 's/^.//' > $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.lst
grep -Fxvf $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.lst $output_folder/03_fastqc/unmapped_reads_qc.lst > $output_folder/04_01_low_complexity_reads/low_complexity_reads_qc.lst

# 04.01.1 seqtk - save discarded low-complexity reads
Cmd="$DockerOptions seqtk subseq"
docker run $Cmd $DockerPath/03_fastqc/unmapped_reads_qc.fq $DockerPath/04_01_low_complexity_reads/low_complexity_reads_qc.lst > $output_folder/04_01_low_complexity_reads/low_complexity_reads_qc.fastq

# 04.01.2 seqtk - make fasta file
Cmd="$DockerOptions seqtk seq -a"
docker run $Cmd $DockerPath/04_01_low_complexity_reads/low_complexity_reads_qc.fastq > $output_folder/04_01_low_complexity_reads/low_complexity_reads_qc.fasta
echo Reads saved to 04_01_low_complexity_reads!

# 04.01.3 grep - count low-complexity reads from fasta file
grep -c ">" $output_folder/04_01_low_complexity_reads/low_complexity_reads_qc.fasta > $output_folder/04_01_low_complexity_reads/low_complexity_reads_qc.txt
echo Total discarded low-complexity reads:
cat $output_folder/04_01_low_complexity_reads/low_complexity_reads_qc.txt

# 05 copy "high-quality reads" qc and higher + complex reads to 05_target_reads for further analysis
echo Organizing reads of interest...
cp $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.fastq $output_folder/05_reads_of_interest/carrierseq_roi.fastq
cp $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.fasta $output_folder/05_reads_of_interest/carrierseq_roi.fasta 
echo Reads saved to 05_reads_of_interest!

# 05.1 grep - count target reads from fasta file
grep -c ">" $output_folder/05_reads_of_interest/carrierseq_roi.fasta > $output_folder/05_reads_of_interest/carrierseq_roi.txt
echo Total reads of interest:
cat $output_folder/05_reads_of_interest/carrierseq_roi.txt

##################
#                #
#  Poisson Calc  #
#                #
##################

echo Starting Poisson calculation...

ChannelsInUse="$output_folder/06_poisson_calculation/03_channels_in_use.txt"
TotalROIs="$output_folder/05_reads_of_interest/carrierseq_roi.txt"
LambdaValue="$output_folder/06_poisson_calculation/04_lambda_value.txt"
ROIChannels="$output_folder/06_poisson_calculation/08_roi_channels_clean.lst"
XCrit="$output_folder/06_poisson_calculation/06_xcrit_threshold_for_dictionary_search.txt"

# 06 grep - extract all channels used, delete duplicates to count unique (n/512) channels used
echo 'Counting total channels in use...'
grep -Eo '_ch[0-9]+_|ch=[0-9]+' $all_reads > $output_folder/06_poisson_calculation/01_reads_channels.lst
awk '!seen[$0]++' $output_folder/06_poisson_calculation/01_reads_channels.lst > $output_folder/06_poisson_calculation/02_channels_used.lst

# 06.01 - count unique channels (n/512)
grep -c "ch" $output_folder/06_poisson_calculation/02_channels_used.lst > $output_folder/06_poisson_calculation/03_channels_in_use.txt
echo 'Channels in use:'
cat $output_folder/06_poisson_calculation/03_channels_in_use.txt

# 06.02 python - calculate lambda for poisson calculation
echo Calculating lambda and x_crit values...
python python/calculate_lambda.py $TotalROIs $ChannelsInUse > $output_folder/06_poisson_calculation/04_lambda_value.txt

# 06.02.1 python - calculate x_critical
python python/xcrit.py $LambdaValue $p_value > $output_folder/06_poisson_calculation/05_read_channel_threshold.txt
sed -n 6p $output_folder/06_poisson_calculation/05_read_channel_threshold.txt > $output_folder/06_poisson_calculation/06_xcrit_threshold_for_dictionary_search.txt
cat $output_folder/06_poisson_calculation/05_read_channel_threshold.txt

# 06.03 grep - get channel list from carrierseq_roi.fasta (now compatible with poretools and albacore fastqs)
echo 'Extracting read IDs and channels from reads of interest....'
grep -Eo '_ch[0-9]+_' $output_folder/05_reads_of_interest/carrierseq_roi.fasta | sed 's/_//g' | sed 's/ch//g' > $output_folder/06_poisson_calculation/07_poretools_roi_channels.lst # Get Channel List from poretools output
awk 'NR % 2 == 0' $output_folder/06_poisson_calculation/07_poretools_roi_channels.lst | sed 's/_//g' | sed 's/ch//g' > $output_folder/06_poisson_calculation/08_roi_channels_clean.lst # Remove duplicate channels from poretools fastq
grep -Eo 'ch=[0-9]+' $output_folder/05_reads_of_interest/carrierseq_roi.fasta | sed 's/ch=//g' >> $output_folder/06_poisson_calculation/08_roi_channels_clean.lst # Get channel list from albacore fastq

# 06.03 pyton - Calculate Frequency and create channel dictionary
echo Creating channel frequency dictionaries...
python python/frequency_calc.py $ROIChannels $XCrit $output_folder/06_poisson_calculation/xx_roi_channel_dictionary.txt $output_folder/06_poisson_calculation/xx_hqnr_channel_dictionary.txt $output_folder/06_poisson_calculation/xx_target_channel_dictionary.txt $output_folder/06_poisson_calculation/09_target_channels.lst
echo Poisson caculation complete! Files saved to 06_poisson_calculation.

##################
#                #
#    Sorting     #
#                #
##################

echo 'Preparing files to begin Poisson sorting...'

PSorter="$output_folder/06_poisson_calculation/09_target_channels.lst" # target reads "good" channel to grep -f with
ROIids="$output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.lst" # reads of interest read id only

# Make grep search compatible for albacore and poretools header format
sed 's/^/ch=/' $PSorter > $output_folder/06_poisson_calculation/10_albacore_target_channels.lst
sed 's/^/_ch/' $PSorter | sed 's/$/_/' > $output_folder/06_poisson_calculation/10_poretools_target_channels.lst

Poretools="$output_folder/06_poisson_calculation/10_poretools_target_channels.lst"
Albacore="$output_folder/06_poisson_calculation/10_albacore_target_channels.lst"

# Dump reads of interest header read id, channel, etc
grep -e '>' $output_folder/05_reads_of_interest/carrierseq_roi.fasta | sed 's/>//g' > $output_folder/05_reads_of_interest/carrierseq_roi_header.lst # dump all rois with channels

ROIHeader="$output_folder/05_reads_of_interest/carrierseq_roi_header.lst"

echo Identifying target reads...
# 07 grep - identify target read ids from channel grep poretools and albacore format
grep -f $Poretools $ROIHeader | awk '{print $1}' > $output_folder/08_target_reads/carrierseq_target_reads.lst
grep -f $Albacore $ROIHeader | awk '{print $1}' >> $output_folder/08_target_reads/carrierseq_target_reads.lst

echo Identifying HQNRs...
# 07.01 grep use grep diff to get hqnr channels
grep -Fxvf $output_folder/08_target_reads/carrierseq_target_reads.lst $ROIids > $output_folder/07_hqnrs/carrierseq_hqnrs.lst

# 07.02 seqtk - extract target reads and make fasta file
echo 'Extracting target reads from carrierseq_roi.fastq ...'
seqtk subseq $output_folder/05_reads_of_interest/carrierseq_roi.fastq $output_folder/08_target_reads/carrierseq_target_reads.lst > $output_folder/08_target_reads/carrierseq_target_reads.fastq # use seqtk to extract target reads
seqtk seq -a $output_folder/08_target_reads/carrierseq_target_reads.fastq > $output_folder/08_target_reads/carrierseq_target_reads.fasta # make target reads fasta
echo 'Target reads saved to 08_target_reads!'
echo Target target reads: 

# 07.02.1 grep - count target reads
grep -c ">" $output_folder/08_target_reads/carrierseq_target_reads.fasta > $output_folder/08_target_reads/carrierseq_target_reads.txt # count target reads
cat $output_folder/08_target_reads/carrierseq_target_reads.txt

# 07.03 seqtk - extract hqnr reads and make fasta file
echo 'Extracting HQNRs from carrierseq_roi.fastq ...'
seqtk subseq $output_folder/05_reads_of_interest/carrierseq_roi.fastq $output_folder/07_hqnrs/carrierseq_hqnrs.lst > $output_folder/07_hqnrs/carrierseq_hqnrs.fastq # use seqtk to extract target reads
seqtk seq -a $output_folder/07_hqnrs/carrierseq_hqnrs.fastq > $output_folder/07_hqnrs/carrierseq_hqnrs.fasta # make hqnrs reads fasta
echo 'HQNRs saved to 07_hqnrs!'
echo Total HQNRs: 

# 07.03.1 grep - count hqnr reads
grep -c ">" $output_folder/07_hqnrs/carrierseq_hqnrs.fasta > $output_folder/07_hqnrs/carrierseq_hqnrs.txt # count hqnrs reads
cat $output_folder/07_hqnrs/carrierseq_hqnrs.txt

# Cleaning up Docker
echo Deleting temporary reads and reference genome docker files...
rm -r $output_folder/fastq_tmp
rm -r $output_folder/reference_tmp

echo 'Thank you for using CarrierSeq!'

# End of file!
