# CarrierSeq: 

CarrierSeq and CarrierSeqXL are sequence analysis 

## Requirements 

>bwa - https://github.com/lh3/bwa
>seqtk - https://github.com/lh3/seqtk
>samtools - https://github.com/samtools/samtools
>fqtrim - https://ccb.jhu.edu/software/fqtrim/

### Building Your Own Docker Image

1. Download & install Docker - https://www.docker.com/
2. Start Docker and increase threads/memory if possible, the carrierseq.sh script is set to 6 CPUs and 18 GB ram.
3. Save Dockerfile to your directory 
4. ```$ cd to/your/directory```
5. ```$ docker build -t <name-your-image> .```

### Or Pull from Docker Hub

1. Start docker
2. run ```$ docker pull mojarro/carrierseq:latest```
