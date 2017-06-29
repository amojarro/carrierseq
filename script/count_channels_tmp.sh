#!/bin/bash

# Extract all channels used, delete duplicates to count unique (n/512) channels used
grep -Eio "_ch[0-9]+_" $DataFolder/fastq/all_reads.fastq | awk '!seen[$0]++' > $DataFolder/06_poisson_calculation/channels_used.lst

# Count unique channels (n/512)
grep -c "ch" $DataFolder/06_poisson_calculation/channels_used.lst > $DataFolder/06_poisson_calculation/channels_in_use.txt
 
# Calculate lambda
python $DataFolder/python/calculate_lambda.py > $DataFolder/06_poisson_calculation/lambda_value.txt

# Calculate x_critical
python $DataFolder/python/xcrit.py > DataFolder/06_poisson_calculation/read_channel_threshold.txt