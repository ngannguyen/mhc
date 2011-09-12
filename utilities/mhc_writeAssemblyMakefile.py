#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#06/22/2011

import os, sys, re

def writeMakefile(f, options):
    f.write("dataDir = %s\n" %options.datadir)
    f.write("sample = %s\n" %options.sample)
    f.write("outdir = %s\n" %options.outdir)
    f.write("ref = %s\n" %options.ref)
    f.write("\n")
    #f.write("all: sample.lst reads ${sample}-mhc.sam runVelveth runVelvetg\n")
    f.write("all: sample.lst reads runVelveth runVelvetg\n")
    f.write("\n")
    f.write("sample.lst:\n")
    f.write("\tls -1 ${dataDir}/${sample}*.bam |gawk 'match($$0, /.+\/([^\/]+)\.bam/, ary) {print ary[1]}' > sample.lst\n")
    f.write("\n")
    f.write("samples = $(shell cat sample.lst)\n")
    f.write("targets = ${samples:%=%-mhcSorted.bam}\n")
    #f.write("readfiles = ${samples:%=%-mhcSorted.bam}\n")
    f.write("range = %s\n" %options.coor)
    f.write("kmer = %s\n" %options.kmer)
    f.write("\n")
    f.write("reads: ${targets}\n")
    f.write("\n")
    f.write("%-mhcSorted.bam:\n")
    f.write("\tsamtools view -b ${dataDir}/$*.bam ${range} > $*-mhc.bam\n")
    f.write("\tsamtools sort -n $*-mhc.bam $*-mhcSorted\n")
    f.write("\trm -f $*-mhc.bam\n")
    f.write("\n")
    #f.write("${sample}-mhc.sam:\n")
    #f.write("\tcat *-mhcSorted.bam > ${sample}-mhc.bam\n")
    #f.write("\tsamtools view ${sample}-mhc.bam > ${sample}-mhc.sam\n")
    f.write("\n")
    f.write("runVelveth:\n")
    #f.write("\tvelveth ${outdir} ${kmer} -reference ${ref} -shortPaired -sam ${sample}-mhc.sam\n")
    f.write("\tvelveth ${outdir} ${kmer} -reference ${ref} -shortPaired -bam ${targets}\n")
    f.write("\n")
    f.write("runVelvetg:\n")
    f.write("\tvelvetg ${outdir} -read_trkg yes -exp_cov auto -cov_cutoff auto\n")
    f.write("\n")
    f.write("clean:\n")
    f.write("\trm -f sample.lst\n")
    f.write("\trm -f *.bam\n")
    f.write("\trm -f *.sam\n")
    f.write("\trm -Rf ${outdir}\n")
    f.write("\n")

def main():
    from optparse import OptionParser
    #Adding options
    usage = "usage: %prog [options]\n"
    parser = OptionParser(usage=usage)
    parser.add_option("-d", "--dataDirectory", dest="datadir", help="Location of the raw reads") 
    parser.add_option("-s", "--sample", dest="sample", help="Sample name")
    parser.add_option("-o", "--outdir", dest="outdir", help="Output directory")
    parser.add_option("-r", "--reference", dest="ref", help="Reference file")
    parser.add_option("-c", "--coordinate", dest="coor", help="Coordinate range of the reads to be extracted. E.g: 6:28,477,797-33,428,773")
    parser.add_option("-k", "--kmer", dest="kmer", help="size of kmer for velvet to use")
    
    (options, args) = parser.parse_args()
    writeMakefile(sys.stdout, options)


if __name__ == "__main__":
    main()
