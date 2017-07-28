# CarrierSeq

## About

bioRxiv doi: xxx.xxx.xxx

CarrierSeq is a sequence analysis workflow for low-input nanopore sequencing. For many environmental samples, the total extractable DNA is far below the current input requirements of nanopore sequencing, preventing “sample to sequence” metagenomics from low-biomass or recalcitrant samples. One approach is to employ carrier sequencing, a method to sequence low-input DNA by preparing the target DNA with a genomic carrier to achieve ideal library preparation and sequencing stoichiometry without amplification. We can then apply CarrierSeq to identify the low-input target reads from the genomic carrier

## Requirements

The CarrierSeq scripts requires the following dependencies to be installed on your local machine.

bwa - https://github.com/lh3/bwa</br>
seqtk - https://github.com/lh3/seqtk</br>
samtools - https://github.com/samtools/samtools</br>
fqtrim - https://ccb.jhu.edu/software/fqtrim/</br>

Alternatively, use Docker and the Docker script.

## Using Docker
### Building Your Own Docker Image

1. Download & install Docker - https://www.docker.com/
2. Start Docker and increase threads/memory if possible, the carrierseq(XL).sh scripts are set to 6 CPUs and 18 GB RAM.
3. Save Dockerfile to your directory.
4. ```cd to/your/directory```
5. ```docker build -t <name-your-image> .```

### Or Pull from DockerHub

1. Start docker.
2. run ```docker pull mojarro/carrierseq:latest```

## Using CarrierSeq 

Reads to be analyzed must be compiled into a single fastq file and the carrier reference genome must be in fasta format.

### Locally

Run CarrierSeq with:

```./carrierseq.sh <all_reads.fastq> <reference_genome.fasta> <q-score> <p-value>```

### Docker

```./carrierseq_docker.sh <all_reads.fastq> <reference_genome.fasta> <q-score> <p-value>```

