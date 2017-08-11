# Primary author : Michael Micorescu (Michael.Micorescu(at)nanoporetech.com)
# Modified by AM for CarrierSeq application

from Bio import SeqIO
import math
from Tkinter import Tk
import sys

name = sys.argv[1]
qs = float(sys.argv[3])
output = sys.argv[2]

count = 0
for rec in SeqIO.parse(name, "fastq"):
    count += 1

qual_sequences = []

cnt = 0
for rec in SeqIO.parse(name, "fastq"):
    rec.letter_annotations["phred_quality"]
    probs = []
    for q in rec.letter_annotations["phred_quality"]:
        e = float(math.pow(10.0,-1*(float(q)/10.0)))
#        print q, e
        probs.append(e)
    av_prob = float(sum(probs))/float(len((rec.letter_annotations["phred_quality"])))
#    print av_prob
    av_q = float(-10.0*(math.log10(float(av_prob))))
 #   print av_prob, av_q
    if av_q >= qs:
	cnt += 1
        qual_sequences.append(rec)

output_handle = open(output +'.fa', "w")
SeqIO.write(qual_sequences, output_handle, "fasta")
output_handle.close()

output_handle = open(output +'.fq', "w")
SeqIO.write(qual_sequences, output_handle, "fastq")
output_handle.close()
