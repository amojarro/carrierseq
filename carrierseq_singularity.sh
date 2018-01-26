#!/bin/sh
# For all steps, we will be mapping a data folder (on host) into the container
# The folder should have all_reads.fastq 
# Download: https://www.dropbox.com/sh/vyor82ulzh7n9ke/AAC4W8rMe4z5hdb7j4QhF_IYa?dl=0

# CarrierSeq Singularity:
# 1. Download data from  link above
# 2. mapping: Perform mapping step of pipeline, mapping the same folder.
# 3. poisson: perform poisson regression on filtered reads
# 4. sorting: Finally, sort results


if [ $# -eq 0 ]
  then
    echo "Please provide a local data folder to map to the container."
    exit 0
fi

DATA=$1

# Check that singularity installed
if [ ! -f "cseq" ]; then

    if type singularity 2>/dev/null; then
        sudo singularity build cseq Singularity
    else
        echo "Please install singularity 2.4+ and run again"
        exit 0
    fi
else
    echo "Found carrierseq.img"
fi

if [ ! -f "$DATA/all_reads.fastq" ]; then
    echo "Please download data to $DATA from https://www.dropbox.com/sh/vyor82ulzh7n9ke/AAC4W8rMe4z5hdb7j4QhF_IYa?dl=0"
    exit 0
fi


#singularity run --bind $DATA:/scif/data cseq run download
singularity run --bind $DATA:/scif/data cseq run mapping
singularity run --bind $DATA:/scif/data cseq run poisson
singularity run --bind $DATA:/scif/data cseq run sorting
