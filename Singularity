# To build the container:
# singularity create --size 2000 carrierseq.img
# sudo singularity bootstrap carrierseq.img Singularity

Bootstrap: docker
From: ubuntu:14.04

%help
    echo "This container is installed with the following apps:"
    tree /scif/apps
    echo "Type singularity --app <appname> <container>.img help to see help


%post
apt-get update &&

# Install Shared Dependencies
apt-get install -y git wget build-essential gcc-multilib apt-utils zlib1g-dev samtools pkg-config python-setuptools python-dev python-tables python-pip python-tk python-numpy python-scipy vim tree

# Install biopython
pip install biopython

%apprun setup
#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo "Please supply an output folder."
    exit 0
fi

output_folder="$1"

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
    tree $output_folder



%appinstall bwa
    git clone https://github.com/lh3/bwa.git build
    cd build && git checkout v0.7.15 && make
    mkdir ../bin && mv -t ../bin bwa bwakit

%apprun bwa
    exec bwa "$@"
    
%applabels bwa
    version v0.7.15

%apphelp bwa
    exec bwa --help


%appinstall seqtk
    git clone https://github.com/lh3/seqtk.git build
    cd build && git checkout v1.2 && make
    mkdir ../bin && mv seqtk ../bin

%applabels seqtk
    version v1.2

%apprun seqtk
    exec seqtk "$@"

%appinstall python
    mkdir bin

%apprun python
    exec /usr/bin/python "$@"

%appfiles python
    python/quality_score_filter.py bin/quality_score_filter.py


%appinstall fqtrim
    wget http://ccb.jhu.edu/software/fqtrim/dl/fqtrim-0.9.5.tar.gz
    tar xvfz fqtrim-0.9.5.tar.gz && cd fqtrim-0.9.5 && make release
    mkdir ../bin && mv fqtrim ../bin

%apphelp fqtrim
    exec fqtrim --help

%applabels fqtrim
    version v0.9.5

%apprun fqtrim
    exec fqtrim "$@"
