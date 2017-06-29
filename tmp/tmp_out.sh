# get channel list from carrierseq_out.fastq
grep -Eio "_ch[0-9]+_" carrierseq_out.fastq | sed 's/_//g' > carrierseq_channels.lst # Get Channel List
sed 's/ch//g' carrierseq_channels.lst > carrierseq_channels_clean.lst # Remove Characters


