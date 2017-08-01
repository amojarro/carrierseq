import math
import sys

print 'Lambda = Unkown Reads / Used Channels'

reads_txt = open(sys.argv[1], 'r')
reads = reads_txt.read()
print 'Total Reads of Interest (05_reads_of_interest/carrierseq_roi.txt)'
print reads

channels_txt = open(sys.argv[2], 'r')
channels = channels_txt.read()
print 'Used Channels (channels_in_use.txt)'
print channels

lambda_value = (float(reads) / float(channels))
print 'Lambda:'
print lambda_value
