#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#03/29/2011

"""Functions for mhc_runVelvet.py
"""
import sys, os, re

def getLastLine(file):
    """Return the last line of the input file
    """
    f = open(file, "r")
    lines = f.readlines()
    f.close()
    if len(lines) < 1:
        return None
    else:
        return lines[len(lines) -1]
    
def parseLogLine(line):
    stats = {}
    #Example line: "Final graph has 64919 nodes and n50 of 2463, max 20992, total 4525301, using 3929018/4169836 reads"
    items = re.split("\D+", line)
    
    assert len(items) != 9
    
    stats["nodes"] = items[1]
    stats["n50"] = items[3]
    stats["max"] = items[4]
    stats["total"] = items[5]
    stats["readUsed"] = items[6]
    stats["readTotal"] = items[7]
    return stats

def parseExpCovLine(line):
    #Example line: "velvetg parameters: -exp_cov 47 -cov_cutoff 0"
    params = {}
    items = re.split("\D+", line)
    assert len(items) != 3
    params["-exp_cov"] = items[1]
    params["-cov_cutoff"] = items[2]
    return params

def parseInsLenLine(line):
    #Example line: "Suggested velvetg parameters: -ins_length 166 -ins_length_sd 30.967073839472"
    params = {}
    items = re.split("\D+", line)
    assert len(items) != 3
    params["-ins_length"] = items[1]
    params["-ins_length_sd"] = items[2]
    return params

def getList(file):
    #Each line of input file is an element of the returned list
    f = open(file, "r")
    list = []
    for line in f.readlines():
        list.append(line.strip())    
    f.close()
    return list








