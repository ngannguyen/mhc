#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#03/29/2011

"""
Wrapper to run velvet for multiple samples, and refine the parameters of velvet until get nearly-optimal assembly

The script can do assemblies for inputed list of 'experiments', each corresponding to one dataset. E.g: aml, gbm, colon, ov...
    Each experiment contains multiple samples, which are required to be listed in [experiment]/"sample.lst". E.g: aml/sample.lst
    Each sample (usually) contains one germline genome and one tumor genome, whose names are required to be listed in [experiment]/[sample]/sample.lst. Eg: aml/TCGA-AB-2907/sample.lst  (TCGA-AB-2907-03A-01D-0739-09_whole)
    

"""

import os, sys, re, time
from optparse import OptionParser
import xml.etree.ElementTree as ET

#from jobTree.src.jobTree import runJobTree
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import nameValue
from sonLib.bioio import getTempDirectory
from sonLib.bioio import setLogLevel

from mhc.ref.lib import *

class Setup(Target):
    """Sets up velvet runs for all samples
    """
    def __init__(self, options):
        Target.__init__(self, time=0.00025)
        self.options = options

    def run(self):
        setLogLevel("DEBUG")
	logger.info("Adding experiments to jobTree\n")
        experiments = getList(self.options.experimentList)
        for exp in experiments:
            self.addChildTarget(RunExperiment(exp, self.options))
        
	#self.setFollowOnTarget(Cleanup())


class RunExperiment(Target):
    """Add children jobs (each child = one sample) for the experiment
    """
    def __init__(self, experiment, options):
        Target.__init__(self, time=0.00025)
        self.exp = experiment
        self.options = options

    def run(self):
        sampleListName = "%s/%s/%s" %(self.options.dataDir, self.exp, "sample.lst")
	logger.info("sampleListName: %s\n" %(sampleListName))
	assert os.path.exists(sampleListName)
        samples = getList(sampleListName)
        for sample in samples:
	    sampleDir = "%s/%s" %(self.exp, sample)
	    self.addChildTarget(RunSample(self.exp, sample, sampleDir, self.options))
	
	#self.setFollowOnTarget(CleanupExperiment(self.exp))

class RunSample(Target):
    """
    """
    def __init__(self, experiment, sample, dir, options):
        Target.__init__(self, time=500)
        self.exp = experiment
        self.sample = sample
        self.dir = dir
        self.options = options

    def run(self):
        readFileList = "%s/%s/%s" %(self.options.dataDir, self.dir, "sample.lst")
	assert os.path.exists(readFileList)
	readFiles = getList(readFileList)
	suffix = "-mhcSorted.bam"
	velvetInputFiles = ""
	for file in readFiles:
	    velvetInputFiles += self.options.dataDir + "/" + self.dir + "/" + file + suffix + " "
	
	#velvetOut = "%s/%s" %(self.dir, "assembly")
	velvetOut = "%s/%s/%s" %(self.options.outdir, self.dir, "assembly")
        if not os.path.exists(self.options.outdir + "/" + self.dir):
            system("mkdir -p %s/%s" %(self.options.outdir, self.dir))

        if not self.options.cont:
	    logger.info("Running velveth...\n")
            velvethCmd = "velveth %s %s -reference %s -shortPaired -bam %s" %(velvetOut, self.options.kmer, self.options.ref, velvetInputFiles)
            #sys.stderr.write("%s\n" %velvethCmd)
	    system(velvethCmd)
            
	    #system("velveth %s %s -reference %s -shortPaired -bam %s" %(velvetOut, self.options.kmer, self.options.ref, velvetInputFiles))
	    logger.info("Finished velveth\n")

	    #Initial parameters of velvetg:
	    #params = {"-read_trkg":"yes", "-exp_cov":"auto", "-cov_cutoff":"auto"}
	    paramStr = "-read_trkg yes -exp_cov auto -cov_cutoff auto"
        else:
	    paramStr = ""

	self.addChildTarget(RunVelvetg(velvetOut, paramStr, self.options.maxIter, 0))
	self.setFollowOnTarget(CleanupSample(self.exp, self.sample, velvetOut, self.options.outdir, self.options.mask))

class CleanupSample(Target):
    """
    """
    def __init__(self, experiment, sample, dir, outdir, mask):
        Target.__init__(self, time=0.00025)
        self.exp = experiment
        self.sample = sample
        self.dir = dir
	self.outdir = outdir
        self.mask = mask

    def run(self):
        outfile = "%s/%s.stats" %(self.outdir, self.exp)
	if os.path.exists(outfile):
	    f = open(outfile, "a")
	else:
	    f = open(outfile, "w")
	    f.write("Sample\tNodes\tN50\tMaxContigLength\tTotal\treadUsed\treadTotal\n")
        logFile = "%s/Log" %(self.dir)
	stats = parseLogLine(getLastLine(logFile))
	f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(self.sample, stats["nodes"], stats["n50"], stats["max"], stats["total"], stats["readUsed"], stats["readTotal"]))
	f.close()

	system("rm -f %s/expCov.txt" %(self.dir))
	system("rm -f %s/insLen.txt" %(self.dir))
        
        if self.mask:
            logger.info("repeatMasking the assembled contigs...\n")
            self.addChildTarget(RepeatMask(self.dir))


class RunVelvetg(Target):
    """
    """
    def __init__(self, dir, paramStr, maxIter, iterNum):
        Target.__init__(self, time=35000)
        self.dir = dir
	self.paramStr = paramStr
        self.maxIter = maxIter
        self.iterNum = iterNum #This velvetg run is the iteration number 'iterNum'
	
	#getting specific velvet parameters:
	items = paramStr.split()
	self.params = {}
	for i in range(0, len(items) -1, 2):
	    self.params[items[i]] = items[i+1]

    def run(self):
        if self.paramStr != "":
	    logger.info("Running velvetg with paramStr: '%s', velvetDir: %s\n" %(self.paramStr, self.dir))
	    velvetgCmd = "velvetg %s %s" %(self.dir, self.paramStr)
	    #system("velvetg %s %s" %(self.dir, self.paramStr))
	    system(velvetgCmd)
	    sys.stderr.write("%s\n" %velvetgCmd)
	    logger.info("Finished velvetg\n")
        
        iterNum = self.iterNum + 1
        if self.maxIter == iterNum:
            return

	assert os.path.exists(self.dir)

	logger.info("Getting new expected Coverage\n")
	expCovFile = "%s/expCov.txt" %(self.dir)
        #assert os.path.exists(self.dir + "/stats.txt")
	expCovCmd = "velvet-estimate-exp_cov.pl %s/stats.txt > %s" %(self.dir, expCovFile)
	system(expCovCmd)
	#sys.stderr.write("%s\n" %(expCovCmd))
	#system("velvet-estimate-exp_cov.pl %s/stats.txt > %s" %(self.dir, expCovFile))
        expCovParams = parseExpCovLine(getLastLine(expCovFile))

        logger.info("Getting new insert-length\n")
	insLenFile = "%s/insLen.txt" %(self.dir)
	insLenCmd = "observed-insert-length.pl %s > %s" %(self.dir, insLenFile)
	system(insLenCmd)
	#sys.stderr.write("%s\n" % insLenCmd)
	#system("observed-insert-length.pl %s > %s" %(self.dir, insLenFile))
	insLenParams = parseInsLenLine(getLastLine(insLenFile))

	#compare current with new parameters to see if need to rerun velvet:
	newExpCov = int(expCovParams["-exp_cov"])
	newInsLen = int(insLenParams["-ins_length"])
	#paramStr = "-read_trkg yes -exp_cov %s -cov_cutoff %s -ins_length %s -ins_length_sd %s" \
	#           %(expCovParams["-exp_cov"], expCovParams["-cov_cutoff"], insLenParams["-ins_length"], insLenParams["-ins_length_sd"])
	paramStr = "-read_trkg yes -exp_cov %s -cov_cutoff auto -ins_length %s -ins_length_sd %s" \
	           %(expCovParams["-exp_cov"], insLenParams["-ins_length"], insLenParams["-ins_length_sd"])
	
	if self.paramStr == "" :
	    logger.info("Run velvet with new parameters: %s\n" %(paramStr))
	    self.addChildTarget(RunVelvetg(self.dir, paramStr, self.maxIter, 0))
	else:    
	    currExpCov = self.params["-exp_cov"]
	    if "-ins_length" in self.params:
		currInsLen = int(self.params["-ins_length"])
            
	    if currExpCov == "auto" or abs(newExpCov - int(currExpCov)) > 5 or \
	       "-ins_length" not in self.params or abs(newInsLen - currInsLen) >  5:
		logger.info("Rerun velvet with new parameters: %s\n" %(paramStr))
		self.addChildTarget(RunVelvetg(self.dir, paramStr, self.maxIter, iterNum))
	    else:
		logger.info("Done with velvetg runs. Final parameters: %s\n" %(paramStr))

class RepeatMask(Target):
    """repeatMask the assembled contigs
    """
    def __init__(self, dir):
        Target.__init__(self, time=3600)
        self.dir = dir

    def run(self):
        contigFile = os.path.join(self.dir, "contigs.fa")
        repeatMaskDir = os.path.join(self.dir, "repeatMask")
        system("mkdir -p %s" %(repeatMaskDir))
        system("RepeatMasker -species human --xsmal -dir %s %s" %(repeatMaskDir, contigFile))
        #system("")

def main():
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)

    parser.add_option("-e", "--experimentList",dest= "experimentList", help= "File containing list of experiments (sample sets). E.g: \naml\ngbm\n")
    parser.add_option("-d", "--dataDir",dest= "dataDir", help= "Location of the reads. E.g: /inside/depot/aml/wg")
    parser.add_option("-r", "--ref",dest= "ref", help= "Fasta file of the reference sequence")
    parser.add_option("-k", "--kmer",dest= "kmer", help= "kmer for velvet", default="25")
    parser.add_option("-o", "--output", dest="outdir", help="Output directory")
    parser.add_option("-c", "--continue", dest="cont", action="store_true", help="If specified, skip the velveth step, continue running velvetg", default=False)
    parser.add_option("-m", "--mask", dest="mask", action="store_true", help="If specified, soft-repeatMask the assembled contigs", default=False)
    parser.add_option("--maxIter", dest="maxIter", help="The number of maximun velvetg iterations", default=2)

    options, args = parser.parse_args()

    if options.experimentList == None:
        raise RuntimeError("No experiment list (sample sets) given")

    #i = Stack( Setup(options.experimentList, options.outdir, options.kmer, options.ref, options.cont) ).startJobTree(options)
    i = Stack( Setup(options) ).startJobTree(options)
    if i:
        raise RuntimeError("The jobTree contains %d failed jobs\n" %i)

if __name__ == "__main__":
    from mhc.ref.mhc_runVelvet import *    
    main()


