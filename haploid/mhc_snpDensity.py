#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#04/01/2011
#

import sys, os, re
from optparse import OptionParser

class Coord:
    def __init__(self, chr, coord):
        self.chr = chr
	self.coord = coord

    def __cmp__(self, other):
        if self.chr < other.chr:
            return -1
        elif self.chr > other.chr:
            return 1
        else:
            return cmp(self.coord, other.coord)
    
    def getStr(self):
        str = "%s:%d" %(self.chr, self.coord)
        return str


class Region:
    def __init__(self, region):
        (chr, start, end) = re.split("[:-]", region)
        self.start = Coord(chr, int(start))
        self.end = Coord(chr, int(end))
    
class Het:
    def __init__(self, line):
        self.desc = line.strip()
	items = self.desc.split('\t')
	self.chr = items[0]
	self.start = int(items[1])
	self.end = int(items[2])
	self.ref = items[3]
	self.norm1 = items[4]
	self.norm2 = items[5]
	self.tumor1 = items[6]
	self.tumor2 = items[7]
	self.id = items[8]
	self.tissue = items[9]
	self.qual = items[10]
	self.ploidy_status = items[11]
	self.variant_type = items[12]

    def inRange(self, region):
        start = Coord(self.chr, self.start)
        end = Coord(self.chr, self.end) 
        if cmp(region.start, start) < 1 and cmp(end, region.end) < 1:
            return True
        else:
	    return False
        
def maxVal(dict):
    max = -1
    for k in dict:
        if max < dict[k]:
            max = dict[k]
    return max

def getTotal(dict):
    total = 0
    for k in dict:
        total += dict[k]
    return total

def printDist(fout, hetDist):
    width = 100
    for type in hetDist:
        dist = hetDist[type]
	fout.write("\n%s\n"% type)
        #max = maxVal(dist)
        total = getTotal(dist)
	for k in sorted(dist.keys()):
            #N = int(width*dist[k]/max)
            N = int(width*dist[k]/total)
            str = '*'*N
	    fout.write("%6d | %6d | %s\n" %(k, dist[k], str))
    return

def getDensity(fin, fout, region):
    if region:
        reg = Region(region)
    else:
        reg = None
    prevCoor = {'LOH':-1, 'SNP': -1, 'Deletion': -1, 'Insertion': -1} 
    hetDist = {'LOH':{}, 'SNP':{}, 'Deletion':{}, 'Insertion':{}}

    for line in fin.readlines():
        het = Het(line)
	if reg and not het.inRange(reg):
	    continue
	
	if het.variant_type not in prevCoor:
	    sys.stderr.write("Unrecognized variant_type: %s\n" % (het.variant_type))
	    continue

	prev = prevCoor[het.variant_type]
	if prev > 0:
	    dist = het.start - prev
	    if dist in hetDist[het.variant_type]:
	        hetDist[het.variant_type][dist] += 1
	    else:
	        hetDist[het.variant_type][dist] = 1
	prevCoor[het.variant_type] = het.end
        
    printDist(fout, hetDist)	        
    return

def main():
    #usage = ""
    parser = OptionParser()
    parser.add_option("-r", "--region", dest="region", help="Location of the region of interest. E.g: chr6:28585776-33536751" )
    #parser.add_option("")

    (options, args) = parser.parse_args()
    getDensity(sys.stdin, sys.stdout, options.region)



if __name__ == "__main__":
    main()
