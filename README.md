# CarrierSeq

## About

bioRxiv doi: https://doi.org/10.1101/175281

CarrierSeq is a sequence analysis workflow for low-input nanopore sequencing which employs a genomic carrier.

Github Contributors: Angel Mojarro (@amojarro), Srinivasa Aditya Bhattaru (@sbhattaru), and Christopher E. Carr (@CarrCE).</br> 
fastq-filter from: https://github.com/nanoporetech/fastq-filter

### Motivation
Long-read nanopore sequencing technology is of particular significance for taxonomic identification at or below the species level. For many environmental samples, the total extractable DNA is far below the current input requirements of nanopore sequencing, preventing “sample to sequence” metagenomics from low-biomass or recalcitrant samples.

### Results
Here we address this problem by employing carrier sequencing, a method to sequence low-input DNA by preparing the target DNA with a genomic carrier to achieve ideal library preparation and sequencing stoichiometry without amplification. We then use CarrierSeq, a sequence analysis workflow to identify the low-input target reads from the genomic carrier.

### Methods
CarrierSeq implements ```bwa-mem``` (Li, 2013) to first map all reads to the genomic carrier then extracts unmapped reads by using ```samtools``` (Li et al., 2009) and ```seqtk``` (Li, 2012). Thereafter, the user can define a quality score threshold and CarrierSeq proceeds to discard low-complexity reads with ```fqtrim``` (Pertea, 2015). This set of unmapped and filtered reads are labeled “reads of interest” and should theoretically comprise target reads and likely contamination. However, reads of interest may also include “high-quality noise reads” (HQNRs), defined as reads that satisfy quality score and complexity filters yet do not match to any database and disproportionately originate from specific channels. By treating reads as a Poisson arrival process, CarrierSeq models the expected reads of interest channel distribution and rejects data from channels exceeding a reads/channels threshold (xcrit). Reads of interest are then sorted into ```08_target_reads``` (reads/channel ≤ xcrit) or ```07_hqnrs``` (reads/channel > xcrit).

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

That's it!

## Using CarrierSeq 

Note: You may first need to make the script executable with:

```chmod +x path/to/carrierseq.sh```
or
```chmod +x path/to/carrierseq_docker.sh```

Reads to be analyzed must be compiled into a single fastq file and the carrier reference genome must be in fasta format.

```cd``` into your CarrierSeq folder containing the bash and python scripts and run CarrierSeq with:

```./carrierseq.sh -i <input.fastq> -r <reference.fasta> -q <q_score> -p <p_value> -o <output_directory> -t <bwa_threads>```

or with Docker...

```./carrierseq_docker.sh -i <input.fastq> -r <reference.fasta> -q <q_score> -p <p_value> -o <output_directory> -t <bwa_threads>```

-i, -r, and -o are mandatory flags, CarrierSeq will use the default values if -q, -p, or -t are not defined:

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

# Discarded reads below the given q-score threshold.
03_01_low_quality_reads/low_quality_unmapped_reads.fasta 
                       /low_quality_unmapped_reads.fastq 
                       /low_quality_unmapped_reads.lst   
                       /low_quality_unmapped_reads.txt   
                       
# Reads with less than 50% of its length detected as low complexity.
04_fqtrim_dusted/unmapped_reads_qc_dusted.fasta 
                /unmapped_reads_qc_dusted.fastq
                /unmapped_reads_qc_dusted.lst 
                /unmapped_reads_qc_dusted.txt
                
# Discarded reads with over than 50% of its length detected as low complexity.               
04_01_low_complexity_reads/low_complexity_reads_qc.fasta 
                          /low_complexity_reads_qc.fastq 
                          /low_complexity_reads_qc.lst   
                          /low_complexity_reads_qc.txt   

# Reads of Interest - should theoretically consist of target reads and contamination,
# but may also include "high-quality noise reads" HQNRs which originate from specific channels.
05_reads_of_interest/carrierseq_roi_header.lst
                    /carrierseq_roi.fasta
                    /carrierseq_roi.fastq
                    /carrierseq_roi.txt

# By treating reads as a Poisson arrival process, CarrierSeq models the expected reads-of-interest 
# channel distribution and rejects data from channels exceeding a reads/channels threshold (xcrit).
06_poisson_caculation/01_reads_channels.lst # all channels used during sequencing.
                     /02_channels_used.lst # Unique channels used during sequencing.
                     /03_channels_in_use.txt # Number of unique channels.
                     /04_lambda_value.txt # Lambda = Unkown Reads / Used Channels.
                     /05_read_channel_threshold.txt # Critical read/channel (xcrit) threshold calculation summary.
                     /06_xcrit_threshold_for_dictionary_search.txt # xcrit value.
                     /07_poretools_roi_channels.lst # Channels used in reads of interest from fastq generated using poretools.
                     /08_roi_channels_clean.lst # Channels used in reads of interest from fastq generated using albacore or minknow or formatted channels from 07_poretools_roi_channels.lst.
                     /09_target_channels.lst # "Good" channels used to sort target reads.
                     /10_albacore_target_channels.lst # "Good" channels list formatted for poretools fastq files.
                     /10_poretools_target_channels.lst # "Good" channel list formatted for albacore/minknow fastq files.
                     /xx_hqnr_channel_dictionary.txt # HQNRs read/channel frequency dictionary for python.
                     /xx_roi_channel_dictionary.txt # Reads of interest read/channel frequency dictionary for python.
                     /xx_target_channel_dictionary.txt # Target reads read/channel frequency dictionary for python.
                     
# Likely HQNRs (reads/channel > xcrit). 
07_hqnrs/carrierseq_hqnrs.fasta
        /carrierseq_hqnrs.fastq
        /carrierseq_hqnrs.lst
        /carrierseq_hqnrs.txt
        
# Likely Target Reads (reads/channel ≤ xcrit).
08_target_reads/carrierseq_target_reads.fasta
               /carrierseq_target_readst.fastq
               /carrierseq_target_reads.lst
               /carrierseq_target_reads.txt
```

## CarrierSeq Example
Supplementary sequencing data available from NCBI </br>
BioProject: https://www.ncbi.nlm.nih.gov/bioproject/398368 </br>
BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN07509071 </br>
SRA Download: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run=SRR5935058 </br>

### Library preparation and sequencing
0.2 ng of B. subtilis DNA was prepared with 1 µg of Lambda DNA using the Oxford Nanopore Technologies (ONT) ligation sequencing kit (LSK-SQK108). The library was then sequenced on a MinION Mark-1B sequencer and R9.4 flowcell for 48 hours and basecalled using ONT’s Albacore (v1.10) offline basecaller.

### CarrierSeq Parameters
q-score = 9 (default) and p-value = 0.05. 

### Sequencing and CarrierSeq Summary
At Q9, the expected B. subtilis abundance is 590 reads for this sequencing data. The xcrit value was calculated to be 7 reads/channel.

```
All Reads (Lambda + B. subtilis + Contamination + Noise)
Total Reads: 547,478 reads
Total Bases: 4,914,693,436 bases
###
Reads of Interest (B. subtilis + Contamination + HQNRs) [05_reads_of_interest]
Total Reads: 1,811 reads
Total Bases: 8,132,374 bases
###
HQNRS [07_hqnrs]
Total Reads: 1,179 reads (including 17 false negative B. subtilis reads)
Total Bases: 7,282,767 bases
###
Target Reads [08_target_reads]
Total Reads: 632 reads (including 574 true positive B. subtilis reads, 4 true positive contamination reads, and 54 false positive HQNRs)
Total Bases: 849,607 bases
```

### ROI Pore Occupancy
The matrix illustrates the reads/channel distribution of B. subtilis, contamination, and HQNRs across all 512 nanopore channels. Here we are able to visually identify overly productive channels (e.g., 191 reads/channel, etc) producing likely HQNRs.
![alt text](https://github.com/amojarro/carrierseq/blob/master/example/carrierseq_roi_q9_p005.png)

### HQNR Pore Occupancy
“Bad” channels identified by CarrierSeq as HQNR-associated (reads/channel > 7).
![alt text](https://github.com/amojarro/carrierseq/blob/master/example/carrierseq_hqnrs_q9_p005.png)

### Target Reads Pore Occupancy
“Good” channels identified by CarrierSeq as non-HQNR-associated (reads/channel ≤ 7). Channels producing 6 or more reads yield HQNRs that have satisfied our CarrierSeq parameters. By imposing a stricter p-value, CarrierSeq may be able to reject more HQNRs (e.g., xcrit = 5).
![alt text](https://github.com/amojarro/carrierseq/blob/master/example/carrierseq_target_reads_q9_p005.png)
