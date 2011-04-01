#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#03/31/2011

import os, sys, re
from optparse import OptionParser

def writeNewSeqs(seq, f, name):
    seqs = re.split('N{10,}', seq)
    seqStart = 0
    for s in seqs:
	if s == '':
	    continue
	size = len(s)
	seqname = '.'.join([name, str(seqStart), str(size)])
	seqStart += size
	f.write('>' + seqname + '\n')
	for i in range(0, size,100):
	    f.write(s[i:i+100] + '\n')
    return

def removeNs(files, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    for file in files:
        outfile = "%s/%s" %(outdir, file)
        f = open(outfile, "w")
        
        seq = ""
        name = ""
	for line in f.readlines():
	    if re.match(">", line):
	        if seq != "":
		    writeNewSeqs(seq, f, name)
	        name = line.strip().strip('>')
		seq = ""
            elif re.search("\w", line):
	        seq += line.strip()
	    
	#last sequence in current file:
	if seq != "":
	    writeNewSeqs(seq, f, name)
	f.close()

def getList(file):
    f = open(file, "r")
    list = []
    for line in f.readlines():
        list.append(line[:len(line) -1])
    f.close()
    return list

def main():
    usage = "usage: %prog [options] <list of sequences>"
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--outdir", dest="outdir", help="Output directory", default="outputs")
    parser.add_option("-s", "--sequences", dest="seqList", help="File that contains list of interest species. Specify this option if do not want to list the species in the arguments", default='')
    (options, args) = parser.parse_args()
    if options.speciesList == '':
        removeNs(args, options.outdir)
    else:
        files = getList(options.seqList)
        removeNs(files, options.outdir)

if __name__ == "__main__":
    main()
