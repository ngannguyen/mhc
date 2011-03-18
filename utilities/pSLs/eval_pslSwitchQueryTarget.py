#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#04/30/2010
#
#This scripts switches query to target, and target to query for input psls

def switch():
   #go through input file (stdin) and get the psl
   for line in sys.stdin.readlines():
      line = (line.split('\n'))[0]
      list = line.split('\t')
      if len(list) != 21:
	 #sys.stderr.write("wrong input format.Should be [Chr\\tStart\\tEnd]\n")
         continue
      #switch strand (list[8]):
      if len(list[8]) == 2:
         list[8] = (list[8])[1] + (list[8])[0]
      elif len(list[8]) == 1:
         list[8] = '+' + (list[8])[0]
      else:
         sys.stderr.write("strand with wrong format: *" + list[8] + "*\n");
      seq = '\t'.join(list[0:4])
      seq = seq + '\t' + '\t'.join(list[6:8])
      seq = seq + '\t' + '\t'.join(list[4:6])
      seq = seq + '\t' + list[8]
      seq = seq + '\t' + '\t'.join(list[13:17])
      seq = seq + '\t' + '\t'.join(list[9:13])
      seq = seq + '\t' + '\t'.join(list[17:19])
      seq = seq + '\t' + list[20]
      seq = seq + '\t' + list[19]
      sys.stdout.write(seq + "\n")

if __name__ == "__main__" :
   import os
   import sys
   import re
   from optparse import OptionParser
   #Adding Options
   usage = "usage: %prog < inputPSLFile > outputFile"
   parser = OptionParser(usage=usage)
   (options, args) = parser.parse_args()
   switch()

