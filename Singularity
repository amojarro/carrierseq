Bootstrap: docker
From: ubuntu:14.04

%post
apt-get update &&

# Install Shared Dependencies
apt-get install -y git wget build-essential gcc-multilib apt-utils zlib1g-dev samtools pkg-config python-setuptools python-dev python-tables python-pip python-tk python-numpy python-scipy

# Install biopython
RUN pip install biopython

%labels
MAINTAINER  Angel Mojarro <mojarro at mit dot edu>
REPRODUCER  Vanessa Sochat <vsochat@stanford.edu>

%appinstall bwa
    git clone https://github.com/lh3/bwa.git .
    git checkout v0.7.15
    make
    mkdir bin
    cp -p bwa /bin

%appinstall seqtk
    git clone https://github.com/lh3/seqtk.git .
    git checkout v1.2
    make
    mkdir bin
    cp -p seqtk bin

%appinstall fqtrim
    wget http://ccb.jhu.edu/software/fqtrim/dl/fqtrim-0.9.5.tar.gz
    tar xvfz fqtrim-0.9.5.tar.gz -C .
    make release
    mkdir bin
    cp -p fqtrim bin
