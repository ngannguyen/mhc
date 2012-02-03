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
from jobTree.src.bioio import getTempFile

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
	if self.options.inputInfo:
	    self.addChildTarget( PreProcess(self.options) )
	else:
	    self.addChildTarget( Start(self.options) )
        
	#self.setFollowOnTarget(Cleanup())

class PreProcess(Target):
    """Extract the MHC (or any other range that specified in the ref file) from original files
    Sort files into dataDir/experiment/sample (or patient)
    Sort bam files by name
    Also print out 'dataDir/exp.lst', dataDir/experiment/sample.lst, and dataDir/experiment/sample/file.lst
    """
    def __init__(self, options):
        Target.__init__(self, time=0.00025)
	self.options = options

    def run(self):
	localTempDir = self.getLocalTempDir()
        options = self.options
        e2p2s = options.exp2patient2samples
	explst = os.path.join(localTempDir, 'exp.lst')
	f = open(explst, 'w')
	for exp in e2p2s:
	    expdir = os.path.join(options.dataDir, exp)
	    system("mkdir -p %s" %( expdir ))
	    self.addChildTarget( PreProcessExp(expdir, e2p2s[exp], options.ref2info) )
            f.write("%s\n" %exp)
	f.close()
	system("mv %s %s" %(explst, options.dataDir))

        self.setFollowOnTarget( Start(self.options) )

class PreProcessExp(Target):
    def __init__(self, expdir, sample2files, ref2info):
        Target.__init__(self, time=0.00025)
	self.expdir = expdir
	self.sample2files = sample2files
	self.ref2info = ref2info

    def run(self): 
	localTempDir = self.getLocalTempDir()
        samplelst = os.path.join(localTempDir, 'sample.lst')
	f = open(samplelst, 'w')
	for s in self.sample2files:
	    f.write('%s\n' %s)
	    sampledir = os.path.join(self.expdir, s)
	    system('mkdir -p %s' %(sampledir))
	    self.addChildTarget( PreProcessSample(sampledir, self.sample2files[s], self.ref2info) )
	f.close()
	system("mv %s %s" %(samplelst, self.expdir))

class PreProcessSample(Target):
    def __init__(self, sampledir, files, ref2info):
        Target.__init__(self)
        self.sampledir = sampledir
	self.files = files
	self.ref2info = ref2info

    def run(self):
	localTempDir = self.getLocalTempDir()
        filelst = os.path.join(localTempDir, 'file.lst')
	f = open(filelst, 'w')
	for file in self.files:
	    #Copy necessary file to local tempdir first:
	    localbam = os.path.join( localTempDir, os.path.basename(file.path) )
	    #localbambai = os.path.join( localTempDir, "%s.bai" % os.path.basename(file.path) )
	    system("ln -s %s %s" %(file.path, localbam))

	    range = self.ref2info[file.ref][1]

	    filename = os.path.basename(file.path).rstrip('.bam')
	    localout = os.path.join( localTempDir, "%s-sorted" %(filename) )

	    logger.info("Pre-processing sample %s\n" %(filename))
	    f.write( "%s\t%s\n" %(filename, self.ref2info[file.ref][0]) )
	    #Extract range and sort by name:
            if os.path.exists( "%s.bai" %file.path ):
                system("ln -s %s.bai %s.bai" %(file.path, localbam))
            else:
                system("samtools index %s" %(localbam))
	    cmd = "samtools view -b %s %s | samtools sort -n - %s" %(localbam, range, localout)
	    #cmd = "samtools view -b %s %s | samtools sort -n - %s" %(file.path, range, localout)
	    system(cmd)
            system("scp -C %s.bam %s" %(localout, self.sampledir))
	    
	    #Clean up right away:
	    system("rm -f %s.bam" %localout)
	    #system("rm -f %s" %localbam)
        system("mv %s %s" %(filelst, self.sampledir))	
	f.close()

class Start(Target):
    def __init__(self, options):
        Target.__init__(self, time=0.00025)
	self.options = options
    
    def run(self):
        experiments = getList( os.path.join(self.options.dataDir, 'exp.lst') )
        for exp in experiments:
	    self.addChildTarget(RunExperiment(exp, self.options))

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
	globalTempDir = self.getGlobalTempDir()
        readFileList = "%s/%s/%s" %(self.options.dataDir, self.dir, "file.lst")
	assert os.path.exists(readFileList)
	readFiles = getList(readFileList)
	suffix = "-sorted.bam"
	velvetInputFiles = ""
	reffile = ''
	for file in readFiles:
	    fileitems = file.split('\t')
	    #if len(fileitems) < 2:
	    #    raise RuntimeError("Wrong file.lst format. Require: <filename>\t<ref>")
	    if len(fileitems) == 2:
	        if reffile == '':
		   reffile = fileitems[1]
		elif reffile != fileitems[1]:
		    raise RuntimeError("Files of the same sample must be mapped to the same reference (%s)." %(file))
	    velvetInputFiles += self.options.dataDir + "/" + self.dir + "/" + fileitems[0] + suffix + " "
        
	#Copy reference to local dir
	localref = ''
	if reffile != '':
	    localref = os.path.join( globalTempDir, os.path.basename(reffile) )
	    system("scp -C %s %s" %(reffile, localref))
	
	velvetOut = os.path.join(self.options.outdir, self.dir, "assembly")
	if not os.path.exists( velvetOut ):
            system("mkdir -p %s" %( velvetOut ))

	localOut = os.path.join(globalTempDir, 'assembly')
	system("mkdir -p %s" %(localOut))

        if not self.options.cont:
	    logger.info("Running velveth...\n")
            velvethCmd = "velveth %s %s -reference %s -shortPaired -bam %s" %(localOut, self.options.kmer, localref, velvetInputFiles)
	    if localref == '':
	        logger.info("No reference file - de novo assembly instead of ref-assisted assembly\n")
                velvethCmd = "velveth %s %s -shortPaired -bam %s" %(localOut, self.options.kmer, velvetInputFiles)
	    system(velvethCmd)
	    #system("velveth %s %s -reference %s -shortPaired -bam %s" %(velvetOut, self.options.kmer, self.options.ref, velvetInputFiles))
	    logger.info("Finished velveth\n")
            #raise RuntimeError("DEBUGGING -- DID IT RUN VELVETH?")
	    #Initial parameters of velvetg:
	    paramStr = "-read_trkg yes -exp_cov auto -cov_cutoff auto"
        else:
	    paramStr = ""

	#self.addChildTarget(RunVelvetg(velvetOut, paramStr, self.options.maxIter, 0))
	#self.addChildTarget(RunVelvetg(globalTempDir, velvetOut, paramStr, self.options.maxIter, 0))
	self.addChildTarget(RunVelvetg(localOut, velvetOut, paramStr, self.options.maxIter, 0))
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
    def __init__(self, globalTempDir, dir, paramStr, maxIter, iterNum):
        Target.__init__(self, time=35000)
	self.tempdir = globalTempDir
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
	#localTempDir = self.getLocalTempDir()
        if self.paramStr != "":
	    logger.info("Running velvetg with paramStr: '%s', velvetDir: %s\n" %(self.paramStr, self.dir))
	    velvetgCmd = "velvetg %s %s" %(self.tempdir, self.paramStr)
	    #system("velvetg %s %s" %(self.dir, self.paramStr))
	    system(velvetgCmd)
	    sys.stderr.write("%s\n" %velvetgCmd)
	    logger.info("Finished velvetg\n")
        
        iterNum = self.iterNum + 1
        if self.maxIter == iterNum:
	    system("mv %s/* %s" %(self.tempdir, self.dir))
            return

	assert os.path.exists(self.tempdir)

	logger.info("Getting new expected Coverage\n")
	expCovFile = "%s/expCov.txt" %(self.tempdir)
        #assert os.path.exists(self.dir + "/stats.txt")
	expCovCmd = "velvet-estimate-exp_cov.pl %s/stats.txt > %s" %(self.tempdir, expCovFile)
	system(expCovCmd)
	#sys.stderr.write("%s\n" %(expCovCmd))
	#system("velvet-estimate-exp_cov.pl %s/stats.txt > %s" %(self.dir, expCovFile))
        expCovParams = parseExpCovLine(getLastLine(expCovFile))

        logger.info("Getting new insert-length\n")
	insLenFile = "%s/insLen.txt" %(self.tempdir)
	insLenCmd = "observed-insert-length.pl %s > %s" %(self.tempdir, insLenFile)
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
	    self.addChildTarget(RunVelvetg(self.tempdir, paramStr, self.maxIter, 0))
	else:    
	    currExpCov = self.params["-exp_cov"]
	    if "-ins_length" in self.params:
		currInsLen = int(self.params["-ins_length"])
            
	    if currExpCov == "auto" or abs(newExpCov - int(currExpCov)) > 5 or \
	       "-ins_length" not in self.params or abs(newInsLen - currInsLen) >  5:
		logger.info("Rerun velvet with new parameters: %s\n" %(paramStr))
		self.addChildTarget(RunVelvetg(self.tempdir, self.dir, paramStr, self.maxIter, iterNum))
	    else:
	        system("mv %s/* %s" %(self.tempdir, self.dir))
		logger.info("Done with velvetg runs. Final parameters: %s\n" %(paramStr))

class RepeatMask(Target):
    """repeatMask the assembled contigs
    """
    def __init__(self, dir):
        Target.__init__(self, time=3600)
        self.dir = dir

    def run(self):
        localTempDir = self.getLocalTempDir()
        contigFile = os.path.join(self.dir, "contigs.fa")
        localContigFile = os.path.join(localTempDir, "contigs.fa")
        repeatMaskDir = os.path.join(localTempDir, "repeatMask")
        system("mkdir -p %s" %(repeatMaskDir))

	#Copy necessary files to local temp dir:
	system("scp -C %s %s" %(contigFile, localContigFile))
        system("RepeatMasker -species human --xsmal -dir %s %s" %(repeatMaskDir, localContigFile))
	#Copy to output dir
	outdir = os.path.join(self.dir, "repeatMask")
	system("mkdir -p %s" %outdir)
	system("mv %s/* %s" %(repeatMaskDir, outdir))
        #system("")

class Sample():
    def __init__(self, line):
        items = line.strip().split(',')
        if len(items) < 20:
            raise RuntimeError("inputInfo file does not have enough field\n")
	    
        self.barcode = items[0] #name: TCGA-AA-3518-11A-01D-1518-05
	nameItems = self.barcode.split('-')
        self.patient = items[5]	
	self.experiment = items[8] #diseaseAbreviation, eg: COAD
	self.ref = items[11] #NCBI36/hg18 or NCBI37/hg19
	self.exptype = items[13] #DNA
	self.type = items[14] #normal/tumor
        self.path = items[19] #file path: /inside/depot3/cwilks/new_order/depot_redos/staging/TCGA-AA-3518-11A-01D.mdups.bam

def readInputInfo(file):
    exp2patient2samples = {} #key = exp, val = {patient: [samples]}
    
    f = open(file, 'r')
    for line in f:
        if re.search('barcode', line):
	    continue
        sample = Sample(line)
        if sample.experiment not in exp2patient2samples:
            exp2patient2samples[sample.experiment] = {sample.patient: [sample]}
	elif sample.patient not in exp2patient2samples[ sample.experiment ]:
            exp2patient2samples[sample.experiment][sample.patient] = [sample]
	else:
            exp2patient2samples[sample.experiment][sample.patient].append( sample )
    f.close()
    return exp2patient2samples

def readRef(file):
    ref2info = {} #key = ref, val = [path, range]
    f = open(file, 'r')
    for line in f:
        items = line.strip().split('\t')
	if len(items) < 3:
	    raise RuntimeError("ref file does not have enough fields. Expected 3, only see %d\n" %(len(items)))
	ref2info[items[0]] = [items[1], items[2]]
    f.close()
    return ref2info

def initOptions(parser):
    parser.add_option("-p", "--preprocess", dest="inputInfo", default=None, help="The cvs file that contains sample information (e.g: filename, experiment type, path to file etc). If specified, the program will sort and extract portion of files listed in this csv into dataDir/experiment/sample/file-sorted.bam(s). NOTE: if specify this option, experimentList can be skipped.")
    parser.add_option("-r", "--ref",dest= "ref", help= "File that maps referenceGenome to path_to_ref_fasta_file, with format:<name>\t<path>\trange. Eg:<HG19>\t</inside/depot/users/nknguyen/mhc/data/hg19/chr6.fa>\t<6:28,477,797-33,428,773>")
    parser.add_option("-e", "--experimentList",dest= "experimentList", help= "File containing list of experiments (sample sets). E.g: \naml\ngbm\n")
    parser.add_option("-d", "--dataDir",dest= "dataDir", help= "Location of the reads. E.g: /inside/grotto/nknguyen/mhcRef/data")
    parser.add_option("-k", "--kmer",dest= "kmer", help= "kmer for velvet, default=25", default="25")
    parser.add_option("-o", "--output", dest="outdir", help="Output directory")
    parser.add_option("-c", "--continue", dest="cont", action="store_true", help="If specified, skip the velveth step, continue running velvetg", default=False)
    parser.add_option("-m", "--mask", dest="mask", action="store_true", help="If specified, soft-repeatMask the assembled contigs", default=False)
    parser.add_option("--maxIter", dest="maxIter", help="The number of maximun velvetg iterations. Default=2", default=2)

def checkOptions(parser, options, args):
    if not os.path.exists(options.dataDir):
        system('mkdir -p %s' %options.dataDir)

    if options.inputInfo:
        if not os.path.exists(options.inputInfo):
	    parser.error("File %s does not exist\n" %options.inputInfo)
        options.exp2patient2samples = readInputInfo(options.inputInfo)

    if options.experimentList == None and options.inputInfo == None:
        parser.error("Neither experiment list nor preprocess was given\n")
    
    if not options.ref or not os.path.exists(options.ref):
        parser.error("Reference-Info File was not specified or does not exist (%s)\n" %(options.ref))
    else:
        options.ref2info = readRef(options.ref)

def main():
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    initOptions(parser)
    options, args = parser.parse_args()
    checkOptions(parser, options, args)

    #i = Stack( Setup(options.experimentList, options.outdir, options.kmer, options.ref, options.cont) ).startJobTree(options)
    i = Stack( Setup(options) ).startJobTree(options)
    if i:
        raise RuntimeError("The jobTree contains %d failed jobs\n" %i)

if __name__ == "__main__":
    from mhc.ref.runVelvet import *    
    main()


