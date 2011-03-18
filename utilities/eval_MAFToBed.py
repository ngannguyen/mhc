#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#04/14/2010
#

import os
import sys
import re
from optparse import OptionParser

def isOverlap(start1, end1, start2, end2):
   if end1 <= start2 or end2 <= start1:
      return 'false'
   else:
      return 'true'

def mafToBed(species):
   #numOverlap = 0
   #blocks = []
   #go through input file (stdin) and get the sequences
   for line in sys.stdin.readlines():
      ps = re.compile('s')
      pspc = re.compile(species)
      chr = "chr"
      if ps.match(line) and pspc.search(line):
	 #0,  1  ,   2  ,   3   ,   4   ,      5     ,    6
         #s, name, start, length, strand, totalLength, sequence
         list = line.split('\t')
	 if len(list) != 7:
	    sys.stderr.write("wrong input format")
	    continue
	 li = list[1].split('.') #split the name to get species.chr
	 if len(li) >= 2:
            chr = li[1]
	    #sys.stderr.write("sequences name must be species.chr\n")
	    #continue
         chrlen = int(list[5])
	 start = int(list[2]) 
         if list[4] == '-':
            start = chrlen - start - int(list[3])
            #start = start - int(list[3])
         end = start + int(list[3])
	 seq = '\t'.join([chr,str(start),str(end)])
	 sys.stdout.write(seq + "\n")
         #testing..
         #for block in blocks:
         #   se = block.split('-')
         #   if isOverlap(int(se[0]), int(se[1]), start, end) == 'true':
         #      numOverlap += 1
         #blck = str(start) + '-' + str(end)
         #blocks.append(blck)
  # sys.stderr.write("Number of overlaps: " + str(numOverlap) + "\n")
         

if __name__ == "__main__" :
   import os
   import sys
   import re
   from optparse import OptionParser
   #Adding Options
   usage = "usage: %prog -s speciesName < input maf file > outputFile"
   parser = OptionParser(usage=usage)
   parser.add_option("-s", "--species", dest="spc", help="species to extract BED records from", default="hg19")
   (options, args) = parser.parse_args()
   mafToBed(options.spc)

