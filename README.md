# CarrierSeq: 

CarrierSeq and CarrierSeqXL.

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
2. Start Docker and increase threads/memory if possible, the carrierseq(XL).sh scripts are set to 6 CPUs and 18 GB RAM.
3. Save Dockerfile to your directory.
4. ```cd to/your/directory```
5. ```docker build -t <name-your-image> .```

### Or Pull from DockerHub

1. Start docker.
2. run ```docker pull mojarro/carrierseq:latest```

## Using CarrierSeq and CarrierSeqXL
### Local Machine

scripts/

Edit the appropriate *.sh file and configure your working directory as such:

```
DataFolder="<your/working/directory>" 

$DataFolder/fastq # fastq reads file to be analyzed - "all_reads.fastq"
$DataFolder/reference # reference genome(s). example - "lambda/lambda.fa"
$DataFolder/python # fastq quality filter python script by Michael Micorescu <Michael dot Micorescu at nanoporetech dot com> and edited by Angel Mojarro
```
Find and replace <reference_1> and <reference_2> (if using CarrierSeqXL) with your reference genome(s).

Then ```./carrierseq.sh``` or ```./carrierseqxl.sh```

python/

(local machine only)

Edit the appropriate *.py file and configure your working directory as such:

calculate_lamnda.py
```
reads_txt = open('<your/working/directory>/05_target_reads/carrierseq_out.txt', 'r')
channels_txt = open('<your/working/directory>/06_poisson_calculation/channels_in_use.txt', 'r')
```



### Docker

Follow the same steps as above except edit ```CarrierSeq="mojarro/carrierseq:latest"``` if youd decide to build your own Docker Image.

Make sure Docker is running then simply ```./carrierseq_docker.sh``` or ```./carrierseqxl_docker.sh```

