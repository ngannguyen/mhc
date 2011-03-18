#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#05/12/2010
#

import os
import sys
import re
from optparse import OptionParser

def mapCoor(start, convertFlag):
   #go through input file (stdin) and get the sequences
   for line in sys.stdin.readlines():
      list = line.split('\t')
      if len(list) < 3:
         sys.stderr.write("BED record must have format: chr\tstart\tend\n")
         sys.stderr.write("current line: " + line)
      else:
         if not convertFlag:
            list[1] = str(int(list[1]) + start)  
            list[2] = str(int(list[2]) + start) 
	 else:
            list[1] = str(int(list[1]) - start + 2)  
            list[2] = str(int(list[2]) - start + 2)
         if len(list) >= 8:
            list[6] = list[1]
            list[7] = list[2] 
         seq = '\t'.join(list)
	 sys.stdout.write(seq )
         

if __name__ == "__main__" :
   #Adding Options
   usage = "usage: %prog start < input Bed file > outputFile"
   parser = OptionParser(usage=usage)
   parser.add_option("-c", dest="convertFlag", action="store_true", help="If specified, convert coordinates to CactusCoordinates.\n")
   (options, args) = parser.parse_args()
   if len(args) >= 1:
      #sys.stderr.write("start = *" + args[0] + "*\n")
      mapCoor(int(args[0]), options.convertFlag)

