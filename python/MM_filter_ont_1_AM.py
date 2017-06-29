from Bio import SeqIO
import math
from Tkinter import Tk


name = '/$DataFolder/02_seqtk/unmapped_reads.fastq'
qs = 9
output = '/$DataFolder/03_fastq9/unmapped_reads_q9'

count = 0
for rec in SeqIO.parse(name, "fastq"):
    count += 1
print("%i reads in fastq file" % count)

qual_sequences = [] # Setup an empty list

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

print cnt,'quality reads saved'

output_handle = open(output +'.fa', "w")
SeqIO.write(qual_sequences, output_handle, "fasta")
output_handle.close()

output_handle = open(output +'.fq', "w")
SeqIO.write(qual_sequences, output_handle, "fastq")
output_handle.close()

