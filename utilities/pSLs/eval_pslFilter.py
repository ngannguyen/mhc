#!/usr/bin/env python

#nknguyen@soe.ucsc.edu

#Read file from STDIN and output to STDOUT
def filter():
   lines = []
   for line in sys.stdin.readlines():
      check = 0
      for l in lines:
         if line == l:
            #sys.stderr.write(line)
            check = 1
            break
      if check == 0:
         sys.stdout.write(line)
         lines.append(line)      

if __name__ == "__main__":
   import os
   import sys
   import re
   from optparse import OptionParser
   #Adding Options
   usage = "usage: %prog < input psl file > outputFile"
   parser = OptionParser(usage=usage)
   (options, args) = parser.parse_args()
   filter()


