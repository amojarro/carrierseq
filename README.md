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

## Using Docker and Dockerhub

1. Download & install Docker - https://www.docker.com/
1. Start docker
2. run ```docker pull mojarro/carrierseq:latest```

That's it, no installing dependencies!

## Using CarrierSeq 

Reads to be analyzed must be compiled into a single fastq file and the carrier reference genome must be in fasta format.

Run CarrierSeq with:

```./carrierseq.sh -i <input.fastq> -r <reference.fasta> -q <q-score> -p <p-value> -o <output_directory>```

or with Docker...

```./carrierseq_docker.sh -i <input.fastq> -r <reference.fasta> -q <q-score> -p <p-value> -o <output_directory>```

CarrierSeq will use the default q-score and p-value if -q and -p are not defined:

```
q_score = 9
p_value = 0.0001 or 0.05/512 active channels
```

Also, you may need to make the script executable with:

```chmod +x path/to/carrierseq.sh```

## CarrierSeq Output 

CarrierSeq will generate the following folders and files within your working directory:

```
00_bwa/bwa_mapped.sam 

01_samtools/bwa_unmapped_reads.lst 
           /bwa_unmapped.sam       

02_seqtk/unmapped_reads.fasta 
        /unmapped_reads.fastq 
        /unmapped_reads.txt   

03_01_low_quality_reads/low_quality_unmapped_reads.fasta 
                       /low_quality_unmapped_reads.fastq 
                       /low_quality_unmapped_reads.lst   
                       /low_quality_unmapped_reads.txt   

03_fastqc/unmapped_reads_qc.fa  
         /unmapped_reads_qc.fq  
         /unmapped_reads_qc.lst 
         /unmapped_reads_qc.txt 

04_01_low_complexity_reads/low_complexity_reads_qc.fasta 
                          /low_complexity_reads_qc.fastq 
                          /low_complexity_reads_qc.lst   
                          /low_complexity_reads_qc.txt   

04_fqtrim_dusted/unmapped_reads_qc_dusted.fasta 
                /unmapped_reads_qc_dusted.fastq
                /unmapped_reads_qc_dusted.lst 
                /unmapped_reads_qc_dusted.txt

05_reads_of_interest/carrierseq_roi.fasta
                    /carrierseq_roi.fastq
                    /carrierseq_roi.txt

06_poisson_caculation/channels_used.lst
                     /channels_in_use.txt
                     /lambda_value.txt
                     /read_channel_threshold.txt

07_hqnrs/carrierseq_hqnrs.fasta
        /carrierseq_hqnrs.fastq
        /carrierseq_hqnrs.txt

08_target_reads/carrierseq_out.fasta
               /carrierseq_out.fastq
               /carrierseq_out.txt
```
