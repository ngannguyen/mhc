#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#01/07/2010
#Extract portions (fragments) of a maf file that are within the inputted range
#E.g: extract the alignments of human chr6 from position 1 to position 2
#Input: species.chr start end


#import sys
#import re

#Read file from STDIN and output to STDOUT
def extract(name, start, end):
   check = 0 #check if the alignment is in range or not, default=0(false)
   curralign = []#current alignment block 
   for line in sys.stdin.readlines():
      pcomment = re.compile("#")
      if pcomment.match(line): #comment line
         sys.stdout.write(line)
         continue
      pa = re.compile("a")
      ps = re.compile("s")
      if pa.match(line): #start a new alignment
         if check == 1: #if previous alignment in range, print to stdout
            for a in curralign:
	       sys.stdout.write(a)
	 curralign = [line]
         check = 0
         continue
      curralign.append(line)
      if ps.match(line): #an alignment line
	 list = line.split()
	 if list[1] == name:
            currstrand = list[4]
            if currstrand == '-':
               currend = int(list[5]) - int(list[2])
               currstart = currend - int(list[3])
            else:
               currstart = int(list[2])
               currend = currstart + int(list[3])
	    #if currstart >= start and  currend <= end:
	    if currend > start and  currstart < end:
	       check = 1
   if check == 1:
      for a in curralign:
         sys.stdout.write(a)
	

if __name__ == "__main__":
   import sys
   import re
   #print "input arguments: " + sys.argv[1] + ";" + sys.argv[2] + ";" + sys.argv[3] + "\n"
   if len(sys.argv) < 4:
      sys.stderr.write("Usage: maf-extract.py <name> <start> <end> < <original_file> > <output_file>\n")
   else:
      extract(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))



