# CarrierSeq

## About

bioRxiv doi: xxx.xxx.xxx

CarrierSeq is a sequence analysis workflow for low-input nanopore sequencing which utlizes a genomic carrier.

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

Run CarrierSeq with:

```./carrierseq.sh``` or with docker ```./carrierseq_docker.sh```

Then select your working directoty, fastq file, reference genome, custom q-score and p-value.

CarrierSeq will generate the following folders and files within your working directory:
```
00_bwa/bwa_mapped.sam #

01_samtools/bwa_unmapped_reads.lst #
           /bwa_unmapped.sam       #

02_seqtk/unmapped_reads.fasta #  
        /unmapped_reads.fastq #
        /unmapped_reads.txt   #

03_01_low_quality_reads/low_quality_unmapped_reads.fasta #
                       /low_quality_unmapped_reads.fastq #
                       /low_quality_unmapped_reads.lst   #
                       /low_quality_unmapped_reads.txt   #

03_fastq9/unmapped_reads_q9.fa  #
         /unmapped_reads_q9.fq  #
         /unmapped_reads_q9.lst #
         /unmapped_reads_q9.txt #

04_01_low_complexity_reads/low_complexity_reads_q9.fasta #
                          /low_complexity_reads_q9.fastq #
                          /low_complexity_reads_q9.lst   #
                          /low_complexity_reads_q9.txt   #

04_fqtrim_dusted/unmapped_reads_q9_dusted.fasta #
                /unmapped_reads_q9_dusted.fastq #
                /unmapped_reads_q9_dusted.lst   #
                /unmapped_reads_q9_dusted.txt   #

05_reads_of_interest/x #
                    /x #
                    /x #

06_poisson_caculation/x #
                     /x #
                     /x #

07_hqnrs/x #
        /x #
        /x #

08_target_reads/carrierseq_out.fasta #
               /carrierseq_out.fastq #
               /carrierseq_out.txt #
```
