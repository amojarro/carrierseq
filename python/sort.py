import sys
# Arguments: in, out1, out2
    
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def writefq(fp,header,seq,qual):
	fp.write("@" + header + "\n")
	fp.write(seq + "\n")
	fp.write("+" + "\n")
	fp.write(qual + "\n")
	return

in_h = open(sys.argv[1],"rU")
out1_h = open(sys.argv[2],"w")
out2_h = open(sys.argv[3],"w")

for record in readfq(in_h):
	# record[0] is read name
	# record[1] is sequence
	# record[2] is quality
	print(record[0])
	# Determine channel for this sequence
	# If channel is on the keep list, write to out1
	writefq(out1_h,record[0],record[1],record[2])
	# If channel is on the reject list, write to out2    
	writefq(out2_h,record[0],record[1],record[2])

#if __name__ == "__main__":
#    import sys
#    n, slen, qlen = 0, 0, 0
#    for name, seq, qual in readfq(sys.stdin):
#        n += 1
#        slen += len(seq)
#        qlen += qual and len(qual) or 0
#    print n, '\t', slen, '\t', qlen
    

