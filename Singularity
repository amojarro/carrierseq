# To build the container (Singularity version 2.4):
# singularity build carrierseq.img Singularity
#
# This is an example of a Singularity container delivering a specific analysis
# pipeline. A set of different steps are installed to each app in the container 
# based on performing a componenent of the analysis. This strategy allows the
# creator to better represent his or her specific pipeline via the container, 
# and makes it easier for a user to run (or even swap out) different steps w/o
# needing to know the subtle details about the software relevant to each step.
# Notably, software that is shared between steps is installed globally to the
# container.

Bootstrap: docker
From: ubuntu:14.04


%help

    CarrierSeq is a sequence analysis workflow for low-input nanopore 
    sequencing which employs a genomic carrier.

    Github Contributors: Angel Mojarro (@amojarro), 
                         Srinivasa Aditya Bhattaru (@sbhattaru),   
                         Christopher E. Carr (@CarrCE),
                         Vanessa Sochat (@vsoch).

    fastq-filter from: https://github.com/nanoporetech/fastq-filter

    To run a typical pipeline, you might do:
    #TODO: add here

%labels
bioRxiv-doi https://doi.org/10.1101/175281


%post
apt-get update

# Install Shared Dependencies
apt-get install -y git wget build-essential 
apt-get install -y gcc-multilib apt-utils
apt-get install -y zlib1g-dev samtools 
apt-get install -y pkg-config python-setuptools 
apt-get install -y python-dev python-tables python-pip python-tk
apt-get install -y python-numpy python-scipy vim tree

# Install seqtk
git clone https://github.com/lh3/seqtk.git build
cd build && git checkout v1.2 && make
mv seqtk /usr/local/bin

# Install biopython
pip install biopython

# The client will map to the respective data folders, always


#%appinstall download
    #wget https://www.dropbox.com/sh/vyor82ulzh7n9ke/AAC4W8rMe4z5hdb7j4QhF_IYa?dl=1
    #tar -xzvf sratoolkit.2.8.2-1-ubuntu64.tar.gz
    #mv sratoolkit.2.8.2-1-ubuntu64/bin/* bin/

%apphelp download
  This module includes the entire sra-toolkit for ubuntu. For this container,
  it is (hard coded) to get the example data, and download to /scif/data. Note
  that you must run with a local folder mapped to /scif/data for this to work,
  as the image runs in read only mode.

      singularity run --app sra-toolkit --bind data:/scif/data carrierseq.img

  This will obtain SRR5935058, use fastq-dump to extract it to /scif/data. You can
  also access any of the commands (eg, prefetch) by exec-ing to the container

      singularity exec --app sra-toolkit --bind data:/scif/data carrierseq.img prefetch

  Or you can obtain your own data with a different method, and map to the container.
  When you have data, look at the next step:

      singularity run --app mapping --bind data:/scif/data carrierseq.img



%apprun download
   wget https://sra-download.ncbi.nlm.nih.gov/traces/sra51/SRR/005795/SRR5935058 -O /scif/data/SRR5935058
   fastq-dump -I  --split-files SRR5935058 -v --outdir /scif/data

%apphelp readme
Print the repository's README.md to the console

%appfiles readme
README.md

%apprun readme
    cat ${SINGULARITY_APPROOT}/README.md


%appfiles reference
reference/lambda_ecoli.fa

%apphelp reference
This module only exists to define the location of the reference, at
/scif/apps/reference/lambda_ecoli.fa

# To get externally:
singularity exec --app reference env | grep REFERENCE_FILE

# To get internally
REFERENCE=("/scif/apps/reference/"`ls /scif/apps/reference`)

%appenv reference
REFERENCE_FILE=/scif/apps/reference/lambda_ecoli.fa
export REFERENCE_FILE


%appsetup mapping
    DATAROOT="${SINGULARITY_ROOTFS}/scif/data/mapping"
    mkdir -p $DATAROOT/00_bwa # map all reads to carrier reference genome
    mkdir -p $DATAROOT/01_samtools # extract unmapped sam file
    mkdir -p $DATAROOT/02_seqtk # extract unmapped reads
    mkdir -p $DATAROOT/03_fastqc # discard low-quality reads < qx
    mkdir -p $DATAROOT/03_01_low_quality_reads # save low-quality reads
    mkdir -p $DATAROOT/04_fqtrim_dusted # discard low complexity reads
    mkdir -p $DATAROOT/04_01_low_complexity_reads # save low-complexity reads

    # Input folder for poisson
    mkdir -p $DATAROOT/05_reads_of_interest # filtered reads to poisson calculation
    

%appinstall mapping

    # Install bwa
    git clone https://github.com/lh3/bwa.git build
    cd build && git checkout v0.7.15
    make
    mv -t ../bin bwa bwakit 
    
    # Install fqtrim
    cd .. && wget http://ccb.jhu.edu/software/fqtrim/dl/fqtrim-0.9.5.tar.gz
    tar xvfz fqtrim-0.9.5.tar.gz && cd fqtrim-0.9.5 && make release
    mv fqtrim ../bin

%appenv mapping
    all_reads=${CSEQ_ALLREADS:-/scif/data/SRR5935058_1.fastq}
    reference_genome="${CSEQ_REF:-}"
    bwa_threads="${CSEQ_BWATHREADS:-1}"
    q_score="${CSEQ_QSCORE:-9}"
    output_folder="${SINGULARITY_APPDATA}"
    export all_reads reference_genome bwa_threads q_score output_folder

%appfiles mapping
python/quality_score_filter.py bin/quality_score_filter.py

%apprun mapping

    echo Indexing reference genome...
    # -01 bwa - Index carrier reference genome
    bwa index $reference_genome

    # 00 bwa - map $all_reads to the $reference_genome
    echo Mapping all reads to reference genome...
    bwa mem -x ont2d -t $bwa_threads $reference_genome $all_reads > $output_folder/00_bwa/bwa_mapped.sam

    # 01 samtools - extract unmapped reads as sam file
    echo Extracting unmapped reads...
    samtools view -S -f4 $output_folder/00_bwa/bwa_mapped.sam > $output_folder/01_samtools/bwa_unmapped.sam

    # 01.1 samtools - identify unmapped reads
    cut -f1 $output_folder/01_samtools/bwa_unmapped.sam | sort | uniq > $output_folder/01_samtools/bwa_unmapped_reads.lst

    # 02 seqtk - extract unmapped reads as fastq file
    seqtk subseq $all_reads $output_folder/01_samtools/bwa_unmapped_reads.lst > $output_folder/02_seqtk/unmapped_reads.fastq

    # 02.1 seqtk - make fasta file
    seqtk seq -a $output_folder/02_seqtk/unmapped_reads.fastq > $output_folder/02_seqtk/unmapped_reads.fasta
    echo Unmapped reads saved to 02_seqtk!

    # 02.2 grep - count reads from fasta file
    grep -c ">" $output_folder/02_seqtk/unmapped_reads.fasta > $output_folder/02_seqtk/unmapped_reads.txt
    echo Total unmapped reads:
    cat $output_folder/02_seqtk/unmapped_reads.txt

    # 03 quality_score_filter.py - discard low-quality reads
    echo Applying quality filter...
    python ${SINGULARITY_APPROOT}/bin/quality_score_filter.py $output_folder/02_seqtk/unmapped_reads.fastq $output_folder/03_fastqc/unmapped_reads_qc $q_score
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
    seqtk subseq $output_folder/02_seqtk/unmapped_reads.fastq $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.lst > $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.fastq

    # 03.01.2 seqtk - make fasta file
    seqtk seq -a $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.fastq > $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.fasta

    # 03.01.3 grep - count low-quality reads from fasta file
    grep -c ">" $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.fasta > $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.txt
    echo Total discarded low-quality reads:
    cat $output_folder/03_01_low_quality_reads/low_quality_unmapped_reads.txt

    # 04 fqtrim - discard low complexity reads
    echo Applying DUST filter...
    fqtrim -D $output_folder/03_fastqc/unmapped_reads_qc.fq > $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.fastq

    # 04.1 seqtk - make fasta file
    seqtk seq -a $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.fastq > $output_folder/04_fqtrim_dusted/unmapped_reads_qc_dusted.fasta
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
    seqtk subseq $output_folder/03_fastqc/unmapped_reads_qc.fq $output_folder/04_01_low_complexity_reads/low_complexity_reads_qc.lst > $output_folder/04_01_low_complexity_reads/low_complexity_reads_qc.fastq

    # 04.01.2 seqtk - make fasta file
    seqtk seq -a $output_folder/04_01_low_complexity_reads/low_complexity_reads_qc.fastq > $output_folder/04_01_low_complexity_reads/low_complexity_reads_qc.fasta
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

%applabels mapping
    bwa-version v0.7.15
    seqtk-version v1.2
    fqtrim-version v0.9.5

%appenv poisson
    MAPPING=${MAPPING:-/scif/data/mapping}
    all_reads=${CARRIERSEQ_ALLREADS:-}
    p_value="${CARRIERSEQ_PVALUE:-0.0001}"
    output_folder="${SINGULARITY_APPDATA}"
    ChannelsInUse="$output_folder/06_poisson_calculation/03_channels_in_use.txt"
    TotalROIs="$MAPPING/05_reads_of_interest/carrierseq_roi.txt"
    LambdaValue="$output_folder/06_poisson_calculation/04_lambda_value.txt"
    ROIChannels="$output_folder/06_poisson_calculation/08_roi_channels_clean.lst"
    XCrit="$output_folder/06_poisson_calculation/06_xcrit_threshold_for_dictionary_search.txt"
    export XCrit ROIChannels LambdaValue TotalROIs ChannelsInUse output_folder p_value MAPPING
 
%appsetup poisson
    DATAROOT="${SINGULARITY_ROOTFS}/scif/data/poisson"
    mkdir -p $DATAROOT/06_poisson_calculation # calculations for sorting "real" reads versus "possible noise"
    mkdir -p $DATAROOT/07_hqnrs # "high-quality noise reads" 
    mkdir -p $DATAROOT/08_target_reads # final output reads to be analyzed if target is unknown

%appfiles poisson
python/calculate_lambda.py bin/calculate_lambda.py
python/xcrit.py bin/xcrit.py
python/frequency_calc.py bin/frequency_calc.py

%apprun poisson
    echo Starting Poisson calculation...
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
    python ${SINGULARITY_APPROOT}/bin/calculate_lambda.py $TotalROIs $ChannelsInUse > $output_folder/06_poisson_calculation/04_lambda_value.txt

    # 06.02.1 python - calculate x_critical
    python xcrit.py $LambdaValue $p_value > $output_folder/06_poisson_calculation/05_read_channel_threshold.txt
    sed -n 6p $output_folder/06_poisson_calculation/05_read_channel_threshold.txt > $output_folder/06_poisson_calculation/06_xcrit_threshold_for_dictionary_search.txt
    cat $output_folder/06_poisson_calculation/05_read_channel_threshold.txt

    # 06.03 grep - get channel list from carrierseq_roi.fasta (now compatible with poretools and albacore fastqs)
    echo 'Extracting read IDs and channels from reads of interest....'
    grep -Eo '_ch[0-9]+_' $MAPPING/05_reads_of_interest/carrierseq_roi.fasta | sed 's/_//g' | sed 's/ch//g' >  $output_folder/06_poisson_calculation/07_poretools_roi_channels.lst # Get Channel List from poretools output
    awk 'NR % 2 == 0' $output_folder/06_poisson_calculation/07_poretools_roi_channels.lst | sed 's/_//g' | sed 's/ch//g' > $output_folder/06_poisson_calculation/08_roi_channels_clean.lst # Remove duplicate channels from poretools fastq
    grep -Eo 'ch=[0-9]+' $MAPPING/05_reads_of_interest/carrierseq_roi.fasta | sed 's/ch=//g' >> $output_folder/06_poisson_calculation/08_roi_channels_clean.lst # Get channel list from albacore fastq

    # 06.03 pyton - Calculate Frequency and create channel dictionary
    echo Creating channel frequency dictionaries...
    python frequency_calc.py $ROIChannels $XCrit $output_folder/06_poisson_calculation/xx_roi_channel_dictionary.txt $output_folder/06_poisson_calculation/xx_hqnr_channel_dictionary.txt $output_folder/06_poisson_calculation/xx_target_channel_dictionary.txt $output_folder/06_poisson_calculation/09_target_channels.lst
    echo Poisson caculation complete! Files saved to 06_poisson_calculation.


%appenv sorting
    MAPPING=${MAPPING:-/scif/data/mapping}
    POISSON=${POISSON:-/scif/data/poisson}
    output_folder="${SINGULARITY_APPDATA}"
    PSorter="$POISSON/06_poisson_calculation/09_target_channels.lst"
    ROIids="$MAPPING/04_fqtrim_dusted/unmapped_reads_qc_dusted.lst"
    Poretools="$POISSON/06_poisson_calculation/10_poretools_target_channels.lst"
    Albacore="$POISSON/06_poisson_calculation/10_albacore_target_channels.lst"
    ROIHeader="$MAPPING/05_reads_of_interest/carrierseq_roi_header.lst"
    export output_folder PSorter ROIids Pretools Albacore MAPPING POISSON ROIHeader

%appsetup sorting
    DATAROOT= "${SINGULARITY_ROOTFS}/scif/data/sorting"
    mkdir -p $DATAROOT/07_hqnrs
    mkdir -p $DATAROOT/08_target_reads

%apprun sorting
    echo 'Preparing files to begin Poisson sorting...'

    # Make grep search compatible for albacore and poretools header format
    sed 's/^/ch=/' $PSorter > $POISSON/06_poisson_calculation/10_albacore_target_channels.lst
    sed 's/^/_ch/' $PSorter | sed 's/$/_/' > $POISSON/06_poisson_calculation/10_poretools_target_channels.lst

    # Dump reads of interest header read id, channel, etc
    grep -e '>' $MAPPING/05_reads_of_interest/carrierseq_roi.fasta | sed 's/>//g' > $MAPPING/05_reads_of_interest/carrierseq_roi_header.lst # dump all rois with channels

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
    echo Total target reads: 

    # 07.02.1 grep - count target reads
    grep -c ">" $output_folder/08_target_reads/carrierseq_target_reads.fasta > $output_folder/08_target_reads/carrierseq_target_reads.txt # count target reads
    cat $output_folder/08_target_reads/carrierseq_target_reads.txt

    # 07.03 seqtk - extract hqnr reads and make fasta file
    echo 'Extracting HQNRs from carrierseq_roi.fastq ...'
    seqtk subseq $output_folder/05_reads_of_interest/carrierseq_roi.fastq $output_folder/07_hqnrs/carrierseq_hqnrs.lst > $output_folder/07_hqnrs/carrierseq_hqnrs.fastq
    seqtk seq -a $output_folder/07_hqnrs/carrierseq_hqnrs.fastq > $output_folder/07_hqnrs/carrierseq_hqnrs.fasta # make hqnrs reads fasta
    echo 'HQNRs saved to 07_hqnrs!'
    echo Total HQNRs: 

    # 07.03.1 grep - count hqnr reads
    grep -c ">" $output_folder/07_hqnrs/carrierseq_hqnrs.fasta > $output_folder/07_hqnrs/carrierseq_hqnrs.txt # count hqnrs reads
    cat $output_folder/07_hqnrs/carrierseq_hqnrs.txt

    echo 'Thank you for using CarrierSeq Singularity!'
