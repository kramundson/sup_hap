#! /usr/bin/env python

import os, sys, math
from optparse import OptionParser

#Comai Lab, Ucdavis Genome Center
#Meric Lieberman, 2011
# This work is the property of UC Davis Genome Center - Comai Lab

# Use at your own risk. 
# We cannot provide support.
# All information obtained/inferred with this script is without any 
# implied warranty of fitness for any purpose or use whatsoever. 
#------------------------------------------------------------------------------

#Part 1: run-mpileup.py

#This program is meant to be run on a directory of sorted.bam files. It will generate a mpileup file with columns for each library.

#INPUT:
#This program is run in a folder full of .sorted.bam files as input

#OUTPUT:
#This program outputs a mpileup file.

#NOTE:
#If the program samtools is not in /usr/bin, then the path to samtools must be specified using the command line parameters

#PARAMETERS, default value in []:
#1. REQUIRED:
#-r or reference_file, The alignment reference (fasta format) [required]
#-o or--output_file, The output mpileup.txt filename [required]
#2. OPTIONAL:
#-q or --mapqual, Minimum mapping quality for an alignment to be used [20]
#-Q or --basequal, Minimum base quality for a base to be considered [20]
#-d or --maxdepth, Max per-BAM depth coverage to avoid excessive memory usage [8000]
#-s or --samtools, File path to Samtools [/usr/bin/samtools]

usage = "\n\n%prog -r reference.fa -o output.txt [-q y] [-Q x] [-s path to Samtools}"
usage += "\nRun in a directory only full of sorted.bam files, will generate a mplieup from all bams."
parser = OptionParser(usage=usage)
parser.add_option("-r", "--reference_file", dest="ref", help="Alignment reference file.")
parser.add_option("-o", "--output_file", dest="dest", help="Output file name.")
parser.add_option("-q", "--mapqual", dest="mapqual", default="20", help="(OPTIONAL, default = 20) Minimum mapping quality for an alignment to be used")
parser.add_option("-Q", "--basequal", dest="basequal", default="20", help="(OPTIONAL, default = 20) Minimum base quality for a base to be considered")
parser.add_option("-d", "--maxdepth", dest="maxdepth", default="8000", help="(OPTIONAL, default = 8000) Max per-BAM depth to avoid excessive memory usage")
parser.add_option("--samtools", "-s", dest="pathSAM",  type = "str", default='/share/apps/samtools-github-1.18/samtools', help="File path to Samtools")

(opt, args) = parser.parse_args()
mapqual = opt.mapqual
basequal = opt.basequal
maxdepth = opt.maxdepth

try:
   file = opt.ref
   o = open(opt.dest, 'w')
except:
   parser.error("Please check your command line paramters with -h or --help")

li = os.listdir(os.getcwd())
ind = filter(lambda x: ".sorted.bam" in x, li)
ind.sort()
a = map(lambda x: ["Cov-"+'-'.join(x.split('_')[:-1]).replace('lib','').replace('Lib',''),'Call-'+'-'.join(x.split('_')[:-1]).replace('lib','').replace('Lib',''),'Qual-'+'-'.join(x.split('_')[:-1]).replace('lib','').replace('Lib','')], ind)
b = [item for sublist in a for item in sublist]
header = '\t'.join(['Chrom', 'Pos', 'Ref']+b)+'\n'

runline = ' '.join(ind)
o.write(header)
o.close()

line = opt.pathSAM + " mpileup -d "+ maxdepth +" -Q "+basequal+" -q "+mapqual+" -f "+file+" "+runline+" >> "+opt.dest
print line
os.system(line)

