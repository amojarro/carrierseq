import math

print 'Lambda = Unkown Reads / Used Channels'

reads_txt = open('/DataFolder/05_target_reads/carrierseq_out.txt', 'r')
reads = reads_txt.read()
print 'Unknown Reads (05_target_reads/carrierseq_out.txt)'
print reads

channels_txt = open('/DataFolder/06_poisson_calculation/channels_in_use.txt', 'r')
channels = channels_txt.read()
print 'Used Channels (channels_in_use.txt)'
print channels

lambda_value = (float(reads) / float(channels))
print 'Lambda:'
print lambda_value

