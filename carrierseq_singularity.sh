#!/bin/sh

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

# 1. Download data, map the data base to an empty folder on our local machine
singularity run --app sra-toolkit --bind $DATA:/scif/data carrierseq.img

# 2. Perform mapping step of pipeline, mapping the same folder.
singularity run --app mapping --bind $DATA:/scif/data carrierseq.img

# 3. perform poisson regression on filtered reads
singularity run --app poisson --bind $DATA:/scif/data carrierseq.img

# 4. Finally, sort results
singularity run --app sorting --bind $DATA:/scif/data carrierseq.img
