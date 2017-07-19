#! /usr/bin/env python

import os, sys, math
import re, subprocess

if not os.environ.has_key('MODULE_VERSION'):
	os.environ['MODULE_VERSION_STACK'] = '3.2.10'
	os.environ['MODULE_VERSION'] = '3.2.10'
else:
	os.environ['MODULE_VERSION_STACK'] = os.environ['MODULE_VERSION']
os.environ['MODULESHOME'] = '/software/modules/3.2.10/x86_64-linux-ubuntu14.04/Modules/3.2.10'

if not os.environ.has_key('MODULEPATH'):
	f = open(os.environ['MODULESHOME'] + "/init/.modulespath", "r")
	path = []
	for line in f.readlines():
		line = re.sub("#.*$", '', line)
		if line is not '':
			path.append(line)
	os.environ['MODULEPATH'] = ':'.join(path)

if not os.environ.has_key('LOADEDMODULES'):
	os.environ['LOADEDMODULES'] = ''
	
def module(*args):
	if type(args[0]) == type([]):
		args = args[0]
	else:
		args = list(args)
	(output, error) = subprocess.Popen(['/software/modules/3.2.10/x86_64-linux-ubuntu14.04/Modules/%s/bin/modulecmd' % os.environ['MODULE_VERSION'], 'python'] + 
			args, stdout=subprocess.PIPE).communicate()
	exec output
	
module("load samtools/0.1.18")
module("load bwa/0.6.2")

#Comai Lab, Ucdavis Genome Center
#Meric Lieberman, 2014
# This work is the property of UC Davis Genome Center - Comai Lab

# Use at your own risk. 
# We cannot provide support.
# All information obtained/inferred with this script is without any 
# implied warranty of fitness for any purpose or use whatsoever. 
#------------------------------------------------------------------------------

#This program processes library FASTQ files through .sorted.bam files.
#This process includes the following:
#1. If option is selected to remove chimeric reads, this is done first. New library files are generated
#   with nc appended for no chimeric, the orignial fastq files are still kept, and a statistic file
#   is generated for the removal called rescan-cut-log.txt
#2a. If a paired end overamp is specified, the INTERLEAVED pair ended Fastq file is then split to forward/reverse
#2b. Specified overamp (single or paired) is performed after bwa alignment is done SINGLE ENDED.
#3. Overamp is run, this gives unique and aligned numbers in the file master-OverAmp.txt.
#   If the -o/--overamp option is used, the unique reads sam files will be used to continue instead of the original
#4. The .bam files are generated, .sorted.bam are made from these. The orignial unsorted bams are removed
#
#INPUT:
#This program must be run in a folder of .fq fastq files. It will look at the contents of the folder and run on all
#files that end with ".fq".
#
#OUTPUT:
#Folders for sam, non-overamp sams, sai, sorted bams, bai, original fq files, uncut fq files (if remove chimeric is used) 
#are generated and relevent files moved to the correct result directory.
#
#PARAMETERS:
#1. REQUIRED:
#   i. -d or --database, the reference database for mapping
#2. OPTIONAL:
#   i. -c or --chimeric, remove chimeric reads, generated new .fq files
#   ii. -o or --overamp, remove duplicate reads, use unique .sam files
#   iii. -p or -- paired, align single ended and overamp pair ended, requires input .fq files be interleaved
#   iv. -t or --thread, default == 1, use additioanl threads for alignment
#   v. -q or --trimqual, default == 20, Default mapping quality, for use as the bwa aln -q X during alignment
#   
#NOTE:
#The default parameters assume overamp-4.py is in the same directory as this program, and that bwa and samtools have been installed to /usr/bin/
#If this is not true, the path of all three (overamp, bwa, samtools) must be specified with the following command line parameters:
#   i. --bwa or -b, path for bwa
#   ii. --samtools or -a, path for samtools
#   iii. --scritps or -s, path to overamp-2
#
#The sorted.bam files that are generated here can then be used to create a mpileup file (using our mplieup script package), 
#which is then used for mutation detection / genotyping (using our MAPS package).





from optparse import OptionParser
usage = "USAGE: bwa-samtools.py -d database_file.fa [-c] [-p] [-o] [-t #threads] <Note overamp can not be down with bwa meme or bwa sw>"
parser = OptionParser(usage=usage)
parser.add_option("-d", "--database-file", dest="database", help="Input database file for mapping.")
parser.add_option("-c", "--chimeric", dest="chimeric",  action="store_true", default = False, help="Remove Chimeric Reads")
parser.add_option("-o", "--overamp", dest="overamp",  action="store_true", default = False, help="Remove Overamplified Reads")
parser.add_option("-O", "--supressoveramp", dest="nooveramp",  action="store_true", default = False, help="Suppress all overamp actions")
parser.add_option("-m", "--mode", dest="mode", type = "str", default='s', help="Read type, s = SE (default), p = PE no salvage, ps = PE with SE salvage, pb = PE with PE salvage")
parser.add_option("-D", "--disableSD", dest="disableSD", action="store_false", default = True, help="If running in ps or pb mode, this turns off collecting pairs mapping in the same direction (0,0) or (16,16)")
parser.add_option("-t", "--thread", dest="threads", default="8", help="How many threads to use during alignment.")
parser.add_option("-q", "--trimqual", dest="trimqual", default="20", help="Default mapping quality, for use as the bwa aln -q X during alignment.")
parser.add_option("--bwa", "-b", dest="pathBWA",  type = "str", default='', help="File path to BWA")
parser.add_option("--samtools", "-a", dest="pathSAM",  type = "str", default='', help="File path to Samtools")
parser.add_option("--scripts", "-s", dest="pathScript",  type = "str", default='/isner/share/scripts/', help="File path to all package scripts.")
parser.add_option("-M", "--mem", dest="memmode", action="store_true", default = False, help="Use the bwa mem algorithm instead of aln defualt (bwa-backtack) ")
parser.add_option("-S", "--sw", dest="swmode", action="store_true", default = False, help="Use the bwa sw algorithm instead of aln defualt (bwa-backtack) ")

(opt, args) = parser.parse_args()




if opt.memmode == True or opt.swmode == True:
   nooveramp = True


opt.mode = str.lower(opt.mode)
if opt.mode not in ['s', 'p', 'ps', 'pb']:
   parser.error("Please check your command line mode paramter, must be s, p, ps, or pb")
   

if opt.disableSD == True:
   DISaddon = ''
else:
   DISaddon = ' -D '


database = opt.database
remchim = opt.chimeric
threads = opt.threads
trimQual = opt.trimqual
pathBWA = opt.pathBWA
pathSam = opt.pathSAM
pathElse = opt.pathScript

#from Isabelle Henry, this looks for chimeric reads to be removed
def rescanCutter(fname):
   print "Removing chimeric"
   all = os.listdir(os.getcwd())
   fqs = filter(lambda x: ".fq" in x or ".fastq" in x, all)
   fqs.sort()
   out = open(fname,'w')
   out.write('File\tIntact\tCut\tRemoved\tTotal\t%intact\t%cut\t%removed\n')
   os.system("mkdir orig-fq")
   for fq in fqs:
      print str(fq)
      f = open(fq)
      newname = str(fq).split('.')[0]+'nc.fq'
      o = open(str(fq).split('.')[0]+'nc.fq','w')
      countyes = 0
      countrem = 0
      countno = 0 
      while True:
         name = f.readline()
         if name == "":
            break
         seq = f.readline()
         plus = f.readline()
         qual = f.readline()
         if 'CATG' in seq[4:]:
            seq = seq[:4]+seq[4:].split('CATG')[0]+'CATG'+'\n'
            if len(seq) > 35:
               o.write(name+seq+plus+qual[:len(seq)-1]+'\n')
               countyes +=1
            else:
               countrem +=1
         else:
            o.write(name+seq+plus+qual)
            countno +=1
      f.close()
      o.close()
   
      total = countno + countyes + countrem
      data = [fq, countno, countyes, countrem, total, round(100.00*countno/total,3), round(100.00*countyes/total,3), round(100.00*countrem/total,3)]
      data = map(lambda x: str(x), data)
      out.write('\t'.join(data)+'\n')   
      os.system("mv "+fq+" orig-fq/")
   out.close()



#/////////// do index for bwa if have not
if '/' in database:
   dbfiles = os.listdir(database[:database.rindex('/')+1])
   dname = database[database.rindex('/')+1:]
else:
   dbfiles = os.listdir(os.getcwd())
   dname = database
ind = filter(lambda x: dname in x, dbfiles)

index = ""
dsize = os.path.getsize(database)
if dsize < 1048576000:
   index = " is "
else:
   index = " bwtsw "   
#makes sure ref has been indexed, will skip if already done
ref = ['', '.pac', '.ann', '.amb', '.bwt', '.sa']
check = list((x in ref for x in map(lambda x: x[len(database):],ind)))
if check.count('False') > 1 or len(check) < len(ref):
   print(pathBWA+"bwa index -a"+index+" "+database) 
   os.system(pathBWA+"bwa index -a"+index+" "+database)


#/////////////////////////

if remchim == True:
   rescanCutter("rescan-cut-log.txt")

#setup result directories
li = os.listdir(os.getcwd())
if opt.memmode == True or opt.swmode == True:
   dirs = ['sam', 'bam', 'bai', 'fq', 'usam']
   if False in list((x in li for x in dirs)):
      os.system("mkdir sam bam bai fq usam")
else:
   dirs = ['sam', 'sai', 'bam', 'bai', 'fq', 'usam']
   if False in list((x in li for x in dirs)):
      os.system("mkdir sam sai bam bai fq usam")



if opt.nooveramp == True:
   os.system("rm -r usam")



#make header for overamp run
if opt.nooveramp == False:
   temp = open('master-OverAmp.txt', 'w')
   temp.write('\t'.join(["File", "Type", "#Unique", "#Aligned", "Percent", "All"])+'\n')
   temp.close()

todo = filter(lambda x: x[-3:] == '.fq', li)
todo.sort()
print "Aligning"
for file in todo:
   print file
   name = file.split('.')[0]
   if opt.memmode == False and opt.swmode == False:
      os.system(pathBWA+"bwa aln -t "+threads+" -q "+trimQual+" "+database+" "+file+" > "+name+"_aln.sai")
      os.system(pathBWA+"bwa samse "+database+" "+name+"_aln.sai "+file+" > "+name+"_aln.sam")
   elif opt.memmode == True:
      os.system(pathBWA+"bwa mem -t "+threads+" "+database+" "+file+" > "+name+"_aln.sam")
   else:
      os.system(pathBWA+"bwa bwasw -t "+threads+" "+database+" "+file+" > "+name+"_aln.sam") 
   if opt.nooveramp == False:
      os.system(pathElse+"overamp-5.py -f "+name+"_aln.sam -o u"+name+"_aln.sam "+ " -m "+opt.mode+DISaddon+" >> master-OverAmp.txt")
   if opt.overamp == True:
      name = 'u'+name


   os.system(pathSam+"samtools view -bS "+name+"_aln.sam"+" > "+name+"_aln.bam")
   os.system(pathSam+"samtools sort "+name+"_aln.bam"+" "+name+"_aln.sorted")
   os.system("rm -f "+name+"_aln.bam")
   os.system(pathSam+"samtools index "+name+"_aln.sorted.bam")
   if opt.nooveramp == False:
      os.system('mv u*_aln.sam usam/') 
   os.system('mv *_aln.sam sam/')     
   os.system("mv "+file+" fq/") 
   if opt.memmode == False and opt.swmode == False:
      os.system("mv *.sai sai/") 
   os.system("mv *.bam bam/")
   os.system("mv *.bai bai/")







