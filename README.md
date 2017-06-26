# CarrierSeq: 

CarrierSeq and CarrierSeqXL are sequence analysis 

## Requirements

The CarrierSeq scripts require the following dependencies to be installed on your local machine.

bwa - https://github.com/lh3/bwa</br>
seqtk - https://github.com/lh3/seqtk</br>
samtools - https://github.com/samtools/samtools</br>
fqtrim - https://ccb.jhu.edu/software/fqtrim/</br>

Alternatively, use Docker and the Docker scripts.

## Using Docker
### Building Your Own Docker Image

1. Download & install Docker - https://www.docker.com/
2. Start Docker and increase threads/memory if possible, the carrierseq(XL).sh scripts are set to 6 CPUs and 18 GB ram.
3. Save Dockerfile to your directory 
4. ```cd to/your/directory```
5. ```docker build -t <name-your-image> .```

### Or Pull from DockerHub

1. Start docker
2. run ```docker pull mojarro/carrierseq:latest```

## Using CarrierSeq and CarrierSeqXL
## Local Machine

## Docker
