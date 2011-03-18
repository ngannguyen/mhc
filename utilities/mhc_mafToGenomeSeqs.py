#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#01/13/2010
#
#Extract sequences from a maf file of desired (input) species into fasta format
#Read in maf file from STDIN, get the coordinates (start pos., end pos.) of aligned fragments for each species
#If two fragments are within and prox (proximity) basepairs (default = 1000), sort them together in a "chunk"
#For each chunk, find min and max coordinates, look up the orignal genome of that species, and pull out the sequence.
#
#Note:1/ Name of sequences in the maf file should be in the format: "species.chromosome"
#     2/ Genome sequences of the species are in a directory called "input" that is specified by the option -i
#	Within this directory, each species has its own directory
#	Each species-directory contains fasta files of that species', with filenames in the format: chromosome.fa
#	E.g: if -i ./genomes/ then:
#	genomes/
#	   hg19/
#	      chr1.fa
#	      chr2.fa	
#	   panTro2/
#             chr8.fa
#	   ponAbe2/
#	      contig1.fa
#The outputs are written out to the directory specified in -o, with the same structure as the input dir
#

import os
import sys
import re
from optparse import OptionParser

def reversecomplement(s):
   """Return the reverse complement of the dna string."""   
   complement = {'n':'n','N':'N','A':'T', 'a':'t', 'C':'G', 'c':'g', 'g':'c', 'G':'C', 'T':'A','t':'a'}
   s = s[::-1]
   list = [complement[base] for base in s]
   return ''.join(list)

#Classes
class Fragment:
   """Represent a sequence fragment, with source, start, end, strand..."""
   #def setvals(self, src, s, e, str, totalLen, seq):
   def setvals(self, src, s, l, str, totalLen, seq):
      items = src.split('.')
      if len(items) > 2:
          self.source = items[0] + "." + "_".join(items[1:])
      else:
          self.source = src #usually: species.chromosome
      self.start = s
      self.size = l
      self.strand = str
      self.seq = seq
      self.srcsize = totalLen
   def mergePos (self, frag):
      """Merge (take min start & max end between the two frags) this fragment with frag"""
      """Two frags have to be on same chr & same strand"""
      selfend = self.start + self.size
      fragend = frag.start + frag.size
      if self.start > frag.start:
         self.start = frag.start
      if selfend < fragend:
         selfend = fragend
      self.size = selfend - self.start
   def reverse (self): #get the reverse complement of the fragment
      temp = self.start
      self.start = self.srcsize - (self.start + self.size)
      #self.end = self.srcsize - temp
      if self.strand == '-':
         self.strand = '+'
      else:
         self.strand = '-'
      self.seq = reversecomplement(self.seq)
   def printfrag(self):
      #print "%s\t%d\t%d\t%s\t%d\t%s\n" % (self.source, self.start, self.size, self.strand, self.srcsize, self.seq) 
      print "%s\t%d\t%d\t%d\t%s\n" % (self.source, self.srcsize, self.start, self.size, self.strand) 

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

	    
def inProx (frag1, frag2, prox):
   """Return True if frag1 and frag2 are in proximity of each other."""
   """Otherwise return False."""
   end1 = frag1.start + frag1.size
   end2 = frag2.start + frag2.size
   if end2 >= frag1.start - prox and frag2.start <= end1 + prox:
      return True
   else:
      return False

def exhaustMerge(list, prox):
   """Merge list elements until no two fragments in the list"""
   """ are within each other proximity"""
   for i in range(0,len(list)-1):
      for j in range(i+1, len(list)):
         if inProx(list[i], list[j], prox):
	    list[i].mergePos(list[j])
	    list.pop(j)
	    exhaustMerge(list,prox)
	    return
   return

def getGenomeFilename(indir, name):
   '''Get the fasta-file's name that contains sequence <name>'''
   spc_chr = name.split('.')
   spc = spc_chr[0]
   if len(spc_chr) > 2:
       chr = "_".join(spc_chr[1:])
   else:
       chr = spc_chr[1]
   filename = indir + spc + "/"+ chr + ".fa"
   try:
      open(filename, 'r')
   except IOError:
      print "Cannot open %s\n" % filename
      filename = indir + spc + "/" + spc +".fa"
   return filename
      
def extractSeq(dir, source, start, size):
   """Extract nucleotide (or a.a) sequence of a source sequence"""
   spc_chr = source.split('.')
   #filename = dir + spc_chr[0] + "/" + spc_chr[1] + ".fa"
   filename = getGenomeFilename(dir, source)
   print "Current input file: %s\n" %filename
   print "Searching for sequence: %s\n" %spc_chr[1]
   try:
      f = open(filename, 'r')
   except IOError:
      print "File %s does not exist\n" %filename
      return ''
   #pn = re.compile('>')
   #pmyseq = re.compile('>' + spc_chr[1])
   if start < 0:
      start = 0
   if size < 0:
      size = 0
   currpos = 0 #sum of length of seqs read in so far
   #currpos = -1 #sum of length of seqs read in so far
   seq = ''
   foundSeq = False #True when found the source's sequence in the input file
   for line in f.readlines():
      if re.match('>',line):
         #if re.search(spc_chr[1]+'\n',line):
         if re.search(spc_chr[1],line) or re.search(line[1:].strip(),spc_chr[1]):
	    print "Line:Found matched sequence!\n"
            foundSeq = True
            continue
         else:
            if foundSeq:
	       print "End matched sequence, break\n"
               break
      if not foundSeq:
         continue
      if len(seq) == size:
         break
      line = line[:len(line)-1] #get rid of '\n'
      currpos += len(line)
      if currpos < start:
         continue
      currstart = len(line) - (currpos - start)
      if currstart < 0: 
         currstart = 0
      currend = currstart + (size - len(seq))
      if currend > len(line):
         currend = len(line)
      seq += line[currstart: currend]
      #seq += line
   return seq

#def writeFrag(name, outdir, files, frag, indir, removeN, lenN):
def getChrom(name):
    items = name.split('.')
    chr = "_".join(items[1:])
    return chr

def writeFrag(name, outdir, files, frag, indir, removeN):
   li = name.split('.')
   #filename = outdir + li[0] + "/" + name + ".fa"
   #filename = outdir + li[0] + "/" + li[1] + ".fa"
   chr = getChrom(name) 
   filename = outdir + li[0] + "/" + chr + ".fa"
   print filename
   wa = 'w'
   if name not in files:
      check_dir(filename)
      files.append(name)
   else:
      wa = 'a'
   f = open(filename, wa)
   
   if frag.strand == '-':
      frag.reverse()
   strand = "1" #always extract sequence from the + strand
   if frag.seq == '':
      frag.seq = extractSeq(indir, name, frag.start, frag.size)
   frag.seq = frag.seq.replace('-','')

   #Remove Ns (unavailable data) in the soft-masked inputed sequence
   #and separate the sequences on two sides of Ns. E.g: ...acgNNNNNNNNNTGCA...
   #will be split into ...acg and TGCA...
   if removeN:
      #p = re.compile('[^N]+')
      #seqs = p.finditer(frag.seq)
      seqs = re.split('(N{10,})', frag.seq)
      #for seq in seqs:
      seqStart = 0
      for i in range(0, len(seqs), 2):
         if seqs[i] == '':
             continue
         if i >= 2:
             seqStart += len(seqs[i-2]) + len(seqs[i-1])
         start = seqStart + frag.start
         size = len(seqs[i])
         seqname = '.'.join([frag.source, str(frag.srcsize), str(start), str(size), strand])
         f.write('>' + seqname + '\n')
         #f.write('>' + frag.source + '\n')
         #print the sequence
         s = seqs[i]
         for i in range(0, size,100):
            f.write(s[i:i+100] + '\n')
   else:
      #seqname = '.'.join([frag.source, str(frag.srcsize), str(frag.start), str(len(frag.seq)), strand])
      seqname = '.'.join([frag.source, str(frag.srcsize), str(frag.start), str(frag.size), strand])
      f.write('>' + seqname + '\n')
      for i in range(0,len(frag.seq),100):
         f.write(frag.seq[i:i+100] + '\n')

   f.close()
   return

#ref = reference source
#def getseqs(indir, outdir, spcList, prox, removeN):
def getseqs(indir, outdir, spcList, prox, removeN, summary):
   #dir = "./outputs/" #may want to change this
   files = [] 
   frags = {} #key = spc.chr, val = list of fragments Objs
   currFrag = Fragment()
   #go through input file (stdin) and get the sequences
   for line in sys.stdin.readlines():
      ps = re.compile('s')
      #pi = re.compile('i')
      if ps.match(line):
         #s, src (spc.chr), start, length, strand, totalLength, sequence
         list = line.split()
	 li = list[1].split('.') #split the name to get species.chr
         species = li[0]
	 if species in spcList:
	    currFrag = Fragment()
            #Set empty sequence 
            currFrag.setvals(list[1],int(list[2]), int(list[3]), list[4], int(list[5]), '')
            #Reverse to positive strand if negative:
            if currFrag.strand == '-':
               currFrag.reverse()  
	    if list[1] not in frags: #if the source chromosome (or contig) is new 
               frags[list[1]] = [currFrag]
            else:
	       flist = frags[list[1]]
               #Merge current fragment to an existing one if it is within 
               #proximity of that fragment. If cannot append it to any 
               #existing frag, add it to the fragment list of this chr
	       checkMerge = False
	       for f in flist:
		  if inProx(currFrag, f, prox):
		     f.mergePos(currFrag)
		     checkMerge = True
		     break
	       if not checkMerge:
	          flist.append(currFrag)

   #Check again if any of the fragments in frags is mergable. If yes merge it
   for chr in frags:
      exhaustMerge(frags[chr], prox)
      #write to output files:
      for frag in frags[chr]:
         #frag.printfrag()
         if summary:
             items = frag.source.split('.')
             sys.stdout.write("%s\t%s\t%d\t%d\n" %(items[0], items[1], frag.start, frag.start + frag.size))
         else:
             writeFrag(frag.source, outdir, files, frag, indir, removeN)

def getSpeciesList(file):
   f = open(file,'r')
   li = []
   for line in f.readlines():
      li.append(line[:len(line)-1])
   return li

def test():
   #Test exhaustMerge:
   list = []
   frag1 = Fragment()
   frag1.setvals('A', 1, 10, '+', 1000, '')
   list.append(frag1)
   frag2 = Fragment()
   frag2.setvals('B', 89, 100, '+', 1000, '')
   list.append(frag2)
   frag3 = Fragment()
   frag3.setvals('C', 55, 15, '+', 1000, '')
   list.append(frag3)
   frag4 = Fragment()
   frag4.setvals('D', 25, 15, '+', 1000, '')
   list.append(frag4)
   print "Before:\n"
   for l in list:
      l.printfrag()
   exhaustMerge(list, 20)
   print "After:\n"
   for l in list:
      l.printfrag()

if __name__ == "__main__" :
   import os
   import sys
   import re
   from optparse import OptionParser
   #Adding Options
   usage = "usage: %prog [options] [list of species, e.g: HUMAN CHIMP GORILLA] < input maf file"
   parser = OptionParser(usage=usage)
   parser.add_option("-o", "--outputs", dest="outdir", help="Outputs directory", default="./outputs/")
   parser.add_option("-n", "--removeN", dest="rmN", action="store_true", help="remove Ns and split the sequence into separate sequences")
   #parser.add_option("-l", "--lenN", dest="lenN", help="minimum size of stretches of Ns to be removed\n", default=10)
   parser.add_option("-i", "--inputs", dest="indir", help="Inputs directory, where whole-chromosome-fasta sequences are.", default="./")
   parser.add_option("-s", "--species", dest="speciesList", help="File that contains list of interest species. Specify this option if do not want to list the species in the arguments", default='')
   parser.add_option("-p", "--proximity",type="int", dest="prox",help="Maximum distance for two fragments to be put in the same region of the genome. Default = 100000 base pairs", default=1000000)
   #parser.add_option("-m", "--merge", dest="merge", action="store_true", default=False, help="Merge consecutive fragments together")
   parser.add_option("-a", "--summary", dest="summary", action="store_true", default=False, help="If specified, only print to stdout the summary sequence ranges")
   (options, args) = parser.parse_args()
   if options.speciesList == '':
      getseqs(modify_dirname(options.indir), modify_dirname(options.outdir), args, options.prox, options.rmN, options.summary)
      #getseqs(modify_dirname(options.indir), modify_dirname(options.outdir), args, options.prox, options.rmN, options.lenN)
   else:
      spcli = getSpeciesList(options.speciesList)
      getseqs(modify_dirname(options.indir), modify_dirname(options.outdir), spcli, options.prox, options.rmN, options.summary)
      #getseqs(modify_dirname(options.indir), modify_dirname(options.outdir), spcli, options.prox, options.rmN, options.lenN)
   #test()
   #seq = extractSeq("./", "my.seq1", 21, 40)
   #print seq + "\n"

