#!/usr/bin/env python

import os, re, sys, subprocess
from optparse import OptionParser

class Sample():
    def __init__( self, line ):
        self.desc = line.strip()
        items = line.strip().split('\t')
	if len(items) < 7:
	    sys.stderr.write( "Each sample line is expect to have 7 fields. Only see %d fields\n" % len(items) )
	    sys.exit( 1 )
        self.name = items[0]
	self.contigNum = int(items[1])
	self.n50 = int(items[2])
	self.maxContigLen = int(items[3])
	self.totalLen = int(items[4])
	self.usedReads = int(items[5])
	self.totalReads = int(items[6])
	if self.totalReads > 0:
	    self.percentReadUsed = self.usedReads/float(self.totalReads)*100
	else:
	    self.percentReadUsed = 0

def prettyInt( number ):
    numstr = str( number )
    prettyStr = ''
    for i in range( len(numstr) ):
        if i > 0 and (len(numstr) - i) %3 == 0:
	    prettyStr = prettyStr + ','
	prettyStr = prettyStr + numstr[i]

    return prettyStr

def prettyFloat( num ):
    return "%.2f" %num

def writeSample( f, items, cellColor ):
    for item in items:
        if cellColor == 1:
	    f.write("& \\cellcolor[gray]{0.9} ")
        else:
	    f.write("& ")
        f.write( item )
    f.write("\\\\\n")
    return
    
def expTab( f, expname, samples, cellColor ):    
    if len(samples) < 1:
        return

    for i in range( len(samples) ):
        s = samples[i]
        items = [s.name, prettyInt(s.contigNum), prettyInt(s.n50), \
	         prettyInt(s.maxContigLen), prettyInt(s.totalLen), \
		 prettyInt(s.usedReads), prettyInt(s.totalReads), \
		 prettyFloat(s.percentReadUsed) ]

	if i == 0: #first line
	    f.write( "\\multirow{%d}{*}{%s} " % (len(samples) + 1, expname) )
	writeSample( f, items, cellColor )
        cellColor = 1 - cellColor

    #Calculate the Average:
    items = [ 'Average', 0, 0, 0, 0, 0, 0, 0 ]
    for s in samples:
        items[1] += s.contigNum
	items[2] += s.n50
	items[3] += s.maxContigLen
	items[4] += s.totalLen
	items[5] += s.usedReads
	items[6] += s.totalReads
    for i in range(1, 7):
        items[i] /= len(samples)
    items[7] = prettyFloat( items[5]/float(items[6])*100 )
    for i in range(1, 7):
        items[i] = prettyInt( items[i] )
    writeSample( f, items, cellColor )
    cellColor = 1 - cellColor
    
    f.write("\\hline\n")
    f.write("\n")
    
    return cellColor

def readTab( file ):
    if not os.path.exists( file ):
        sys.stderr.write( "File %s does not exist\n" % file)
        sys.exit(1)
    f = open( file, 'r' )
    f.readline()
    samples = []
    for line in f.readlines():
        samples.append( Sample(line) )
    
    samples = sorted( samples, key=lambda s:s.n50 )
    f.close()
    return samples

def tab( args, outfh ):
    cellColor = 1
    for file in args:
        expname = file.split('.')[0]
        samples = readTab( file )
	cellColor = expTab( outfh, expname, samples, cellColor )


def tableCloser(f, captionStr, label):
    f.write("\\end{tabular}\n")
    f.write("}\n")
    f.write("\\caption{%s}\n" %captionStr)
    f.write("\\end{center}\n")
    f.write("\\label{%s}" %label)
    f.write("\\end{table}\n\n")
    
def tableHeader(f, title):
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n") 
    f.write("\\scalebox{0.7}{%\n")
    #f.write("\\begin{tabular}{|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|}\n")
    f.write("\\begin{tabular}{ccccccccc}\n")
    f.write("\\multicolumn{9}{c}{%s} \\\\\n" %title)
    f.write("\\hline\n")
    f.write("\\hline\n")
    f.write("Experiment & Sample & NumContigs & N50 & MaxContigLength & TotalLength & UsedReads & TotalReads & PercentReadsUsed\\\\\n")
    #f.write("\\cline{3-8}\n")
    f.write("\\hline\n")

def writeDocumentStart(f):
    f.write("\\documentclass[11pt]{article}\n")
    f.write("\\usepackage{epsfig}\n")
    f.write("\\usepackage{multirow}\n")
    f.write("\\usepackage{graphicx}\n")
    f.write("\\usepackage{array}\n")
    f.write("\\usepackage{slashbox}\n")
    f.write("\\usepackage{color}\n")
    f.write("\\usepackage[table]{xcolor}\n")
    f.write("\\usepackage{rotating}\n")
    f.write("\n")

    f.write("\\newcommand{\\figref}[1]{Figure~\\ref{fig:#1}}\n")
    f.write("\\newcommand{\\tabref}[1]{Table~\\ref{tab:#1}}\n")
    f.write("\n")

    f.write("\\textwidth=6.5in\n")
    f.write("\\textheight=9in\n")
    f.write("\\oddsidemargin=0in\n")
    f.write("\\evensidemargin=0in\n")
    f.write("\\topmargin=0in\n")
    f.write("\\topskip=0in\n")
    f.write("\\headheight=0in\n")
    f.write("\\headsep=0in\n")
    f.write("\n")

    f.write("\\begin{document}\n")
    return

def writeDocumentEnd(f):
    f.write("\\end{document}\n")

def main():
    usg = "usage: %prog [options] stats1 [stats2 stats3 ...]"
    parser = OptionParser(usage=usg)
    parser.add_option("-n", "--title", dest = "title", help="Title of the table\n", default="Assembly Stats")
    parser.add_option("-o", "--output", dest = "out", default="assemblyStats.tex", help="Name of output file")
    (options, args) = parser.parse_args()

    f = open( options.out, 'w' )
    writeDocumentStart(f)
    tableHeader(f, options.title)

    tab( args, f )
    
    captionStr = ""
    label = "%s" %options.title
    tableCloser(f, captionStr, label)
    
    writeDocumentEnd(f)
    f.close()

if __name__ == "__main__":
    main()


