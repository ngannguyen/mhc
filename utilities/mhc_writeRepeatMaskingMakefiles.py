#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#03/13/2010

def modify_dirname(dir):
   """Add slash / at the end of the directory name if it doesnt have yet"""
   if (not re.search('/$', dir)): #not end with /
      dir = dir + '/'
   return dir

def check_dir(filepath):
   """Check if directories on the path of a file exist, and create them if not."""
   d = os.path.dirname(filepath)
   if not os.path.exists(d):
      os.makedirs(d)

def writeMakefiles(outdir, file):
   f = open(file, 'r')
   for line in f.readlines():
      li = line.rstrip().split('\t')
      if len(li) < 2:
         sys.stderr.write(file + " possibly has wrong format. Should be: <species\\tScientific_name> for each line\n")
      else:
         writeMakefile( outdir, li[0], li[1])
   

def writeMakefile(outdir, species, speciesName):
   filename = outdir + species + "/Makefile"
   check_dir(filename)
   f = open(filename, 'w')
   f.write("outputDir = ../..\n")
   f.write("species = " + species + "\n")
   f.write("speciesName = " + speciesName + "\n")
   f.write("files = $(wildcard *.fa)\n\n")
   f.write(".PHONY: all clean\n\n")
   f.write("all: ${files:%=${outputDir}/${species}/%} dirs\n\n")
   f.write("dirs:\n")
   f.write("\tmkdir repeatMasking\n")    
   f.write("\tmkdir ${outputDir}/${species}\n\n")
   f.write("${outputDir}/${species}/%: % dirs\n")    
   f.write("\tmkdir repeatMasking/$(basename $<)\n")
   f.write("\trepeatMasking_doCluster.py --genome=$< --workDir=repeatMasking/$(basename $<) --species=\"${speciesName}\"\n")    
   f.write("\ttwoBitToFa repeatMasking/$(basename $<)/*.rmsk.2bit ${outputDir}/${species}/$<\n\n")
   f.write("clean:\n")
   f.write("\trm -R ./repeatMasking\n")
   f.close()
   return

if __name__ == "__main__" :
   import os
   import sys
   import re
   from optparse import OptionParser
   #Adding Options
   usage = "usage: %prog [options] <speciesList>\n"
   parser = OptionParser(usage=usage)
   parser.add_option("-o", "--outputs", dest="outdir", help="Outputs directory", default="./outputs/")
   (options, args) = parser.parse_args()
   if len(args) < 1:
      sys.stderr.write(usage)
   else:
      writeMakefiles(modify_dirname(options.outdir), args[0])

