#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#04/14/2010
#
#Convert cactus MAF blocks to genome's coordinate
#1/ If names of (cactus) sequences in the input MAF have the format:
#species.chr.chrSize.start.length.strand
#Where: 
#chr: is the chromosome (or contig) the cactus sequence is from
#chrSize: is total Length of 'chr'
#start: is start position of the input sequence on chromosome 'chr'
#length: length of input sequence
#strand: = 1 if the sequence is on the forward strand of 'chr', otherwise = 0
#
#Example:  hg19.chr7.159138663.55083851.195218.1
#
#(a side note: strand is in 0 and 1 instead of -+ because for some other pipeline 
#which use repeatMasker, '-' or '+' in the sequence name would yield errors
#
#2/ If sequences names are not in the format above, their location on the corresponding
#genomes need to be inputed from a file with the following format:
#SequenceName\tspcecies\tchr\tchrSize\tstart\tlength\tstrand
#with SequenceName has to be the same with the input MAF file
#
#Example:
#MySequence	hg19	chr7	159138663	55083851	195218	1
#
#3/ The input coordinates information have priority over the info in the sequences names

import os
import sys
import re
from optparse import OptionParser

#def mapCoorPos(chrsize, start, totalLen, strand, coor):
#mapCoorPos is equivalent with mapCoorNeg --> just mapCoor
#def mapCoor(chrsize, start, totalLen, strand, coor, len):
def mapCoor(chrsize, start, strand, coor, len, mafstrand):
   #(cactus) input sequence is on the positive strand of the chromosome (reference sequence)
   #Map 'coor' to the coordinates of the reference sequence
   #refsize: size of reference sequence, start: location of current seq on the ref-seq
   #strand: strand of the maf fragment relatively to the (cactus) input sequence
   #strand = 1 if current seq is on the forward strand of the ref-seq, 0 otherwise
   if strand == '0':#NEED TO take into account mafStrand HERE
      sys.stderr.write("mapCoor, refStrand == 0\n")
      return chrsize -1 -(start + coor + len)
   else:
      if mafstrand == '+':
         return coor + start 
      else:
         return chrsize -(start + len) + coor

#def mapCoorNeg(chrsize, start, totalLen, strand, coor):
   #(cactus) input sequence is on the negative strand of the chromosome (reference sequence)
   #
   ## posStart = chrsize -1 -start
   #if strand == '-': #-- = +
      ##posCoor = totalLen -1 -coor
      #return chrsize - totalLen - start + coor #posStart - posCoor
   #else:
      #return start + coor #'start' on chr. neg strand + 'coor' = mapped coor on chr. neg strand

def mapStrand(refstrand, strand):
   if refstrand == '1':
      if strand == '+':
         return '+'
      else:
         return '-'
   else:
      #sys.stderr.write("mapStrand, refStrand == 0\n")
      if strand == '+':
         return '-'
      else:
         return '+'

def getModifiedName(file):
   f = open(file, 'r')
   dict = {} #key = orignal seq name, val = species.chr.order
   chrDict = {} #key = chr.spc, val = list of visited threads
   for line in f.readlines():
      if re.match('s', line):
         list = line.split('\t')
	 nameOrder = list[1].split('_')
	 name = nameOrder[0]
	 nameli = name.split('.')
	 namechrom = nameli[0] + '.' + nameli[1]
	 if len(nameOrder) < 2 and namechrom not in dict:
	    dict[list[1]] = namechrom 
	 if len(nameOrder) == 2:
	    if namechrom not in chrDict:
	       chrDict[namechrom] = [list[1]]
	    else:
	       threadList = chrDict[namechrom]
	       if list[1] not in threadList:
	          threadList.append(list[1])
   f.close()

   for seq in chrDict:
      threadList = chrDict[seq]
      threadList.sort()
      for i in range(0, len(threadList)):
         dict[threadList[i]] = seq + '.' + str(i)

   for seq in sorted(dict):
      sys.stderr.write("%s: %s\n" % (seq, dict[seq]))
   return dict



def maf_mapCoor(inputfile, outputfile, infoLi):
   """Map coordinate of input MAF records to the genome coordinate"""
   f = open(outputfile, 'w')
   #dict = getModifiedName(inputfile)
   inf = open(inputfile, 'r')

   #go through input file (stdin) and get the sequences
   for line in inf.readlines():
      ps = re.compile('s')
      pi = re.compile('i')
      pe = re.compile('e')

      if pi.match(line) or pe.match(line) or ps.match(line):
         #0    1     2      3       4          5          6
         #s, name, start, length, strand, totalLength, sequence
         list = line.split('\t')
	 #if len(list) != 7:
	 #   sys.stderr.write("wrong input format")
	 #   continue
         
         #li = genome location info 
	 #(spc = li[0], chr = li[1], chrSize = li[2], genomeStart = li[3], 
	 # totalLength = li[4], strand = li[5])
         if list[1] not in infoLi:
	    li = list[1].split('.') #split the name to get species.chr
         else:
            li = infoLi[list[1]].split('.')

	 if len(li) < 6:
            #convert the underscore:
            li1 = li[1].split('_')
            if len(li1) < 2:
               f.write(line)
            else:
               list[1] = li[0] + li1[1] + '.' + li1[0]
	       seq = '\t'.join(list)
	       f.write(seq)
	    #sys.stderr.write(li[0] + " don't have coor info\n" )
	    continue

	 #list[1] = dict[list[1]]
	 #if li[5] == '1':#refStrand
	    #list[2] = mapCoorPos(int(li[2]), int(li[3]), int(li[4]), list[4], int(list[2]))
	 #else:
	    #list[2] = mapCoorNeg(int(li[2]), int(li[3]), int(li[4]), list[4], int(list[2]))
	 #list[2] = mapCoor(int(li[2]), int(li[3]), int(li[4]), list[4], int(list[2]), int(list[3]))

         #list[1] = li[0] + '.' + li[1] #name 
         li5 = li[5].split('_')
         if len(li5) >= 2:
            #list[1] = list[1] + '.' + li5[len(li5) -1]
            list[1] = li[0] + li5[len(li5) -1] + '.' + li[1] #NAME
         else:
            list[1] = li[0] + '.' + li[1] #name 
         if ps.match(line) or pe.match(line):
            li[5] = li5[0]
	    list[2] = mapCoor(int(li[2]), int(li[3]), li[5], int(list[2]), int(list[5]), list[4])
            list[2] = str(list[2])
            list[4] = mapStrand(li[5],list[4])
            list[5] = li[2]

	 seq = '\t'.join(list)
	 f.write(seq)
      else:
         f.write(line)

   inf.close()
   f.close()

def getCoorInfo(file):
   f = open(file,'r')
   li = {} #key=SequenceName, val = spc.chr.chrSize.start.len.strand
   i = 0
   for line in f.readlines():
      list = line.rstrip().split('\t')
      if len(list) < 7:
        sys.stderr.write("infoFile must be in the following format:\n"); 
        sys.stderr.write("SequenceName\tspcecies\tchr\tchrSize\tstart\tlength\tstrand\n");
      if list[0] not in li:
        li[list[0]] = '.'.join(list[1:7])
   return li

if __name__ == "__main__" :
   import os
   import sys
   import re
   from optparse import OptionParser
   #Adding Options
   usage = "usage: %prog [options] < input maf file"
   parser = OptionParser(usage=usage)
   parser.add_option("-i", "--input", dest="inputFile", help="intput maf file")
   parser.add_option("-o", "--output", dest="outputFile", help="MAFs with new coordinates will be printed to this file", default="./out.maf")
   parser.add_option("-a", "--info", dest="infoFile", help="", default='')
   (options, args) = parser.parse_args()
   infoLi = {}
   if options.infoFile != '':
      infoLi = getCoorInfo(options.infoFile)
   maf_mapCoor(options.inputFile, options.outputFile, infoLi)

