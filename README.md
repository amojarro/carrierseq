# CarrierSeq

## About

bioRxiv doi: xxx.xxx.xxx

CarrierSeq is a sequence analysis workflow for low-input nanopore sequencing which utlizes a genomic carrier.

Contributors: Angel Mojarro (@amojarro) and Christopher E. Carr (@CarrCE).

### Motivation
Long-read nanopore sequencing technology is of particular significance for taxonomic identification at or below the species level. For many environmental samples, the total extractable DNA is far below the current input requirements of nanopore sequencing, preventing “sample to sequence” metagenomics from low-biomass or recalcitrant samples.

### Results
Here we address this problem by employing carrier sequencing, a method to sequence low-input DNA by preparing the target DNA with a genomic carrier to achieve ideal library preparation and sequencing stoichiometry without amplification. We then use CarrierSeq, a sequence analysis workflow to identify the low-input target reads from the genomic carrier

### Methods
CarrierSeq implements bwa-mem (Li, 2013) to first map all reads to the genomic carrier then extracts unmapped reads by using samtools (Li et al., 2009) and seqtk (Li, 2012). Thereafter, the user can define a quality score threshold and CarrierSeq proceeds to discard low-complexity reads with fqtrim (Pertea, 2015). This set of unmapped and filtered reads are labeled “reads of interest” and should theoretically comprise target reads and likely contamination. However, reads of interest may also include “high-quality noise reads” (HQNRs), defined as reads that satisfy quality score and complexity filters yet do not match to any database and disproportionately originate from specific channels. By treating reads as a Poisson arrival process, CarrierSeq models the expected reads of interest channel distribution and rejects data from channels exceeding a reads/channels threshold (xcrit). Reads of interest are then sorted in ```08_target_reads``` (reads/channel ≤ xcrit) or ```07_hqnrs``` (reads/channel > xcrit).

## Requirements

The CarrierSeq scripts requires the following packages to be installed on your local machine.

Biopython - http://biopython.org/</br>
SciPy - https://www.scipy.org/</br>
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

Note: You may need to first make the script executable with:

```chmod +x path/to/carrierseq.sh```

Reads to be analyzed must be compiled into a single fastq file and the carrier reference genome must be in fasta format.

Run CarrierSeq with:

```./carrierseq.sh -i <input.fastq> -r <reference.fasta> -q <q_score> -p <p_value> -o <output_directory> -t <bwa_threads>```

or with Docker...

```./carrierseq_docker.sh -i <input.fastq> -r <reference.fasta> -q <q_score> -p <p_value> -o <output_directory> -t <bwa_threads>```

-i -r and -o are mandatory flags, CarrierSeq will use the default values if -q -p or -t are not defined:

```
bwa_threads = 1 
q_score = 9
p_value = 0.0001 or 0.05/512 active channels
```

## CarrierSeq Output 

CarrierSeq will generate the following folders and files within your working directory:

```
# All reads mapped to the carrier reference genome.
00_bwa/bwa_mapped.sam 

# Unmapped reads to the carrier.
01_samtools/bwa_unmapped_reads.lst 
           /bwa_unmapped.sam       

02_seqtk/unmapped_reads.fasta 
        /unmapped_reads.fastq 
        /unmapped_reads.txt   
        
# Reads equal to or greater than a given q-score threshold (default = 9).
03_fastqc/unmapped_reads_qc.fa  
         /unmapped_reads_qc.fq  
         /unmapped_reads_qc.lst 
         /unmapped_reads_qc.txt 

# Discarded reads below the given q-score threshold
03_01_low_quality_reads/low_quality_unmapped_reads.fasta 
                       /low_quality_unmapped_reads.fastq 
                       /low_quality_unmapped_reads.lst   
                       /low_quality_unmapped_reads.txt   
                       
# Reads with less than 50% of its length detected as low complexity
04_fqtrim_dusted/unmapped_reads_qc_dusted.fasta 
                /unmapped_reads_qc_dusted.fastq
                /unmapped_reads_qc_dusted.lst 
                /unmapped_reads_qc_dusted.txt
                
# Discarded reads with over than 50% of its length detected as low complexity               
04_01_low_complexity_reads/low_complexity_reads_qc.fasta 
                          /low_complexity_reads_qc.fastq 
                          /low_complexity_reads_qc.lst   
                          /low_complexity_reads_qc.txt   

# Reads of Interest - should theoretically consist of target reads and contamination,
# but may also include "high-quality noise reads" HQNRs which originate from specific channels.
05_reads_of_interest/carrierseq_roi.fasta
                    /carrierseq_roi.fastq
                    /carrierseq_roi.txt

# By treating reads as a Poisson arrival process, CarrierSeq models the expected reads-of-interest 
# channel distribution and rejects data from channels exceeding a reads/channels threshold (xcrit)
06_poisson_caculation/channels_used.lst
                     /channels_in_use.txt
                     /lambda_value.txt
                     /read_channel_threshold.txt
                     
# Likely HQNRs (reads/channel > xcrit) 
07_hqnrs/carrierseq_hqnrs.fasta
        /carrierseq_hqnrs.fastq
        /carrierseq_hqnrs.txt
        
# Likely Target Reads (reads/channel ≤ xcrit)
08_target_reads/carrierseq_out.fasta
               /carrierseq_out.fastq
               /carrierseq_out.txt
```
## Known Issues

CarrierSeq is currently only compatible with fastq files generated using the new header implemented after Albacore 1.2 or the latest MinKNOW live basecalling option. 

Example fastq header:
```
{read_id} runid={run_id} read={read_number} ch={channel_id} start_time={start_time_utc}
```
or 
```
@006b70e0-2dd9-4f1a-a18a-521b3d4668d7 runid=ebfdf272ef6e469d83f5d045e4351889384fef6f read=179 ch=106 start_time=2017-05-07T04:28:48Z
```

