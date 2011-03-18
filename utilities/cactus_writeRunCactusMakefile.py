#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#01/09/2011
import os, sys, re

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
   

def writeMakefile(options):
   filename = options.outdir + "/Makefile"
   check_dir(filename)
   f = open(filename, 'w')
   f.write("include %s\n" %(options.include))
   f.write("binPath=%s\n" %(options.binPath))
   f.write("libPath=%s\n" %(options.libPath))
   f.write("\n")
   f.write("experimentXml = experiment.xml\n")
   f.write("databaseString = '%s'\n" %(options.databaseString))

   list = options.refseq.split(".")
   if len(list) < 6:
      sys.stderr.write("Wrong refseq format, must have 6 fields, instead of %s\n" %(options.refseq))
   refseq = list[0]
   chrom = list[1]
   start = int(list[3])
   end = start + int(list[4])

   f.write("refseq = %s\n" %(refseq))
   f.write("species = \"%s\"\n" %(options.species))
   f.write("\n")
   f.write("all: runCactus\n")
   f.write("\n")

   f.write("runCactus: ${experimentXml} ${binPath}/cactus_workflow.py\n")
   f.write("\tcactus_workflow.py --batchSystem %s --experiment ${experimentXml} --buildReference --setupAndBuildAlignments --logDebug --jobTree ./jobTree\n" %(options.batchSystem))
   f.write("\n")

   f.write("runCactusMaf: ${binPath}/cactus_MAFGenerator\n")
   f.write("\tcactus_MAFGenerator -c ${databaseString} -d 0 -e \"cactus.maf\" -f\n")
   f.write("\teval_mafMapCoor.py -o cactusM.maf  -i cactus.maf\n")
   f.write("\n")

   f.write("runCactusBed: ${binPath}/cactus_bedGenerator\n")
   f.write("\tcactus_bedGenerator -b \"${refseq}\" -c ${databaseString} -d 0 -e \"${refseq}.bed\"\n")
   f.write("\n")

   f.write("runCactusAugMaf: ${binPath}/cactus_augmentedMaf\n")
   f.write("\tcactus_augmentedMaf -b ${species} -c ${databaseString} -d 0 -e \"${refseq}-A.maf\"\n")
   f.write("\teval_mafMapCoor.py -i \"${refseq}-A.maf\" -o \"${refseq}-AM.maf\"\n")
   f.write("\n")

   f.write("runCactusRef: ${binPath}/cactus_MAFGenerator ${binPath}/cactus_bedGenerator ${binPath}/cactus_getReferenceSeq\n")
   f.write("\tcactus_MAFGenerator -c ${databaseString} -d 0 -e \"ref.maf\" -f -g \"reference\"\n")
   f.write("\tcactus_bedGenerator -b \"reference\" -c ${databaseString} -d 0 -e \"ref.bed\"\n")
   f.write("\tcactus_getReferenceSeq -b \"reference\" -c ${databaseString} -d 0 -e \"ref.fa\"\n")
   f.write("\n")

   f.write("refGene.bed:\n")
   f.write("\thgsql -e \"SELECT * FROM refGene WHERE chrom = '%s' AND txStart >= %d AND txEnd <= %d AND name LIKE 'NM%%'\" ${refseq} | sed '1d' | cut -f 2- | sort -u > refGene.gp\n" %(chrom, start, end))
   f.write("\tgenePredToCdsBed.py -c -l < refGene.gp > refGene.bed\n")
   f.write("\n")

   f.write("refGene1.bed: refGene.bed\n")
   f.write("\teval_bedMapCoor.py -c %d < refGene.bed > refGene1.bed\n" %(start))
   f.write("\n")
   
   f.write("refGene2.bed: refGene.bed\n")
   f.write("\tcat refGene.bed | sort -n --key=2 > refGene2.bed\n")
   f.write("\n")
   
   
   f.write("runCactusGenemap: ${binPath}/cactus_geneMap refGene1.bed\n")
   f.write("\tcactus_geneMap -c ${databaseString} -o \"cdsGeneMap.xml\" -s \"%s\" -g \"refGene1.bed\"\n" %(options.refseq))
   f.write("\n")
   #f.write("\tgeneMap.py -o genemap -b stats -c cat -d ntb -e cns -f breakLen < cdsGenemap.xml\n")

   f.write("runCactusGenemap2: ${binPath}/cactus_geneMap2 refGene2.bed\n")
   f.write("\tcactus_geneMap2 -c ${databaseString} -o \"cdsGeneMap2\" -s \"%s\" -g \"refGene2.bed\" > cdsGeneMap2.verbose\n" %(refseq))
   f.write("\n")
   
   f.write("clean:\n")
   f.write("\trm -Rf jobTree\n")
   f.close()
   return

if __name__ == "__main__" :
   from optparse import OptionParser
   #Adding Options
   usage = "usage: %prog [options]\n"
   parser = OptionParser(usage=usage)
   parser.add_option("-o", "--outdir", dest="outdir", help="Directory where the Makefile will be")
   parser.add_option("-i", "--include", dest="include", help="Path to include.mk")
   parser.add_option("-b", "--binPath", dest="binPath", help="Bin path" )
   parser.add_option("-l", "--libPath", dest="libPath", help="Library path" )
   parser.add_option("-n", "--speciesNames", dest="species", help="List of species. E.g: \"hg18 panTro2 ponAbe2 rheMac2\"" )
   parser.add_option("-d", "--databaseString", dest="databaseString", help="Cactus database string")
   parser.add_option("-r", "--refseq", dest="refseq", help="Reference sequence name, e.g hg18")
   parser.add_option("-s", "--batchSystem", dest="batchSystem", help="The type of batch system to run the job(s) with, currently can be 'singleMachine'/'parasol'/'acidTest'/'gridEngine'")
   

   (options, args) = parser.parse_args()
   writeMakefile(options)

