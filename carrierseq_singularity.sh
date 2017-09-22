#!/bin/sh

# CarrierSeq Singularity:
# 1. download: Download data, map the data base to an empty folder on our local machine
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
if [ ! -f "carrierseq.img" ]; then

    if type singularity 2>/dev/null; then
        sudo singularity build carrierseq.img Singularity
    else
        echo "Please install singularity 2.4+ and run again"
        exit 0
    fi
else
    echo "Found carrierseq.img"
fi


singularity run --app download --bind $DATA:/scif/data carrierseq.img
singularity run --app mapping --bind $DATA:/scif/data carrierseq.img
singularity run --app poisson --bind $DATA:/scif/data carrierseq.img
singularity run --app sorting --bind $DATA:/scif/data carrierseq.img
