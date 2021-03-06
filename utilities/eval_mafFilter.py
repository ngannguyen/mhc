#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#08/10/2010
#Read in maf file from STDIN, 
#printed out to STDOUT only species of interest

#Read file from STDIN and output to STDOUT
def filter(speciesList):
   check = 0 #equal to 1 if the MAF has at least one spcies in speciesList
   curralign = []#current alignment block 
   currSpeciesLi = []#list of species in current MAF
   for line in sys.stdin.readlines():
      pcomment = re.compile("#")
      if pcomment.match(line): #comment line
         #sys.stdout.write(line)
         continue
      pa = re.compile("a")
      ps = re.compile("s")
      if pa.match(line): #start a new MAF
         check = 0
         for spc in currSpeciesLi:
            if spc in speciesList:
               check = 1
               break
         if check == 1: #if previous MAF has at least one species in specieslist, print to stdout
            sys.stdout.write('\n')
            for a in curralign:
               sys.stdout.write(a)
	 curralign = [line]
         currSpeciesLi = []
         continue
      if ps.match(line): #an alignment line
	 list = line.split()
         li = list[1].split('.')
         if li[0] in speciesList:
            curralign.append(line)
            if li[0] not in currSpeciesLi:
               currSpeciesLi.append(li[0])
   check = 0
   for spc in currSpeciesLi:
      if spc in speciesList:
         check = 1
         break
   if check == 1:
      for a in curralign:
         sys.stdout.write(a)

	
def getSpeciesList(file):
   f = open(file,'r')
   li = []
   for line in f.readlines():
      li.append(line[:len(line)-1])
   return li

if __name__ == "__main__":
   import os
   import sys
   import re
   from optparse import OptionParser
   #Adding Options
   usage = "usage: %prog [options] [list of species, e.g: hg19 panTor2] < input maf file > output maf file"
   parser = OptionParser(usage=usage)
   parser.add_option("-s", "--species", dest="speciesList", help="File that contains list of interest species. File format: each species per line. Specify this option if do not want to list the species in the arguments", default='')
   (options, args) = parser.parse_args()
   if options.speciesList == '':
      filter(args)
   else:
      spcli = getSpeciesList(options.speciesList)
      filter(spcli)

