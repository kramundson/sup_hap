#! /usr/bin/env python
import sys, math, os, time
from optparse import OptionParser
from collections import defaultdict

#Comai Lab, Ucdavis Genome Center
#Meric Lieberman, 2011
# This work is the property of UC Davis Genome Center - Comai Lab

# Use at your own risk. 
# We cannot provide support.
# All information obtained/inferred with this script is without any 
# implied warranty of fitness for any purpose or use whatsoever. 
#------------------------------------------------------------------------------

#This program looks at a sam file and picks out the unique reads.
#It can be used on a single or paired ended alignment.
#
#INPUT:
#THis program takes a samfile as input mapped SE
#
#OUTPUT:
#This program outputs a unique read sma file, and prints alignment statistics to the screen.
#The output columns are [Filename, #unique, #aligned, #unique, total reads].
#
#NOTE:
#This program is designed to be used in conjuction wit bwa-samtools-do-all.py, 
#with its printed output being redirected to a combined all lib files, 
#that is why there are no headers for the printed output.
#Also, overamp generates the unique files while doing the counts, not independently.
#
#PARAMETERS:
#1. REQUIRED:
#   i. -f or --samfile, The input.sam file
#   ii.  -o or--outfile, The output unique .sam file name
#2. OPTIONAL:
#   i. -s or --sort, sort sam by chrom and start position
#   ii. -p or --paired, Switches to pair-end sam file. Default is single ended

usage = "USAGE: overAmp.py -f samfile.sam -o outputfile [-p] [-s]"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--samfile", dest="sam", help="Input .sam file.")
parser.add_option("-o", "--outfile", dest="out", help="Output file.")
parser.add_option("-m", "--mode", dest="mode", type = "str", default='NA', help="Read type, s = SE, p = PE no salvage, ps = PE with SE salvage, pb = PE with PE salvage")
parser.add_option("-D", "--disableSD", dest="disableSD", action="store_false", default = True, help="If running in ps or pb mode, this turns off collecting pairs mapping in the same direction (0,0) or (16,16)")

(opt, args) = parser.parse_args()

opt.mode = str.lower(opt.mode)
if opt.mode not in ['s', 'p', 'ps', 'pb']:
   parser.error("Please check your command line mode paramter, must be s, p, ps, or pb")



#print time.ctime()
try: 
   f = open(opt.sam)
   o = open(opt.out,'w')
except:
   parser.error("Please check your command line paramters with -h or --help")

def form(flo):
   return str(flo).split('.')[0]+'.'+str(flo).split('.')[1][:3]

#Read in header names and lengths
refseq = []
maxlength = 0
while True:
   x = f.readline()
   if x[0] != '@':
      break
   if x[:3] == "@SQ":
      o.write(x)
      refseq.append(x.split('\t')[1].split(':')[1])
      maxlength = max(maxlength, int(x[:-1].split('\t')[2].split(':')[1]))
      continue 

#Take header names, find out all possible lengths and pick the longest one
#then going through all possible cuts of the longest name, figure out what the optimal cut is
#the optimal cut will be the one wil the smallest differnece in the possible front and bac keys
lens = list(set(map(lambda x: len(x), refseq)))
lens.sort(reverse = True)
poss = []
for i in range(lens[0]):
   front = set(map(lambda x: x[:i], refseq))
   back = set(map(lambda x: x[i:], refseq))
   poss.append([i, len(front), len(back), abs(len(front) - len(back))])
poss.sort(key = lambda x: x[3])
cut = poss[0][0]

#for the length cut, take the square root of the longest one, should always be less than ~7k
numcut = int(math.sqrt(maxlength))+1

reads = 0
pealign = 0
peunique = 0
sealign = 0
seunique = 0
oealign = 0
oeunique = 0
f.seek(0)



#Treatign as PE mapped
if opt.mode in ['p', 'ps', 'pb']:
   peres = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))))
   #            Type                key1                key2               pos1mul           pos1mod              pos2m   pos2mod
   seres = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
   #              F/R               key1                 key2              posmul   posmod
   while True:
      x = f.readline()
      if x == "":
         break
      if x[0] == '@':
         continue
      y = f.readline()
      l = x.split('\t')
      m = y.split('\t')
      reads +=2
      if l[2] == '*' and m[2] == '*':
         continue       
      
      flagset = [l[1],m[1]]
      # IF CORRECT PE FLAGS
      if (flagset == ['0','16'] or flagset == ['16','0']) and l[2] == m[2]:
         typer = 'reg'
         #good pe
         if flagset == ['0','16']:
            pos1 = int(l[3])
            pos2 = int(m[3])+len(m[9])-1
         elif flagset == ['16','0']:
            pos1 = int(m[3])
            pos2 = int(l[3])+len(l[9])-1         
         name = l[2]
         key0 = typer[:]        
         key1 = name[:cut]
         key2 = name[cut:]
         key3 = pos1 / numcut
         key4 = pos1 % numcut
         key5 = pos2 / numcut
         key6 = pos2 % numcut    
         pealign +=2     
         try: 
            if key6 in peres[key0][key1][key2][key3][key4][key5]:
               continue
            else:
               peres[key0][key1][key2][key3][key4][key5].append(key6)
         except:
            x[45678987654]         
         
         o.write(x+y)
         peunique+=2
      elif (flagset == ['0','0'] or flagset == ['16','16']) and l[2] == m[2] and opt.mode not in ['s', 'p'] and opt.disableSD:
         #good pe
         if flagset == ['0','0']:
            pos1 = int(l[3])
            pos2 = int(m[3])
            typer = 'bfor'
         elif flagset == ['16','16']:
            pos1 = int(l[3])+len(l[9])-1 
            pos2 = int(m[3])+len(m[9])-1 
            typer = 'brev'       
         name = l[2]      
         key0 = typer[:]  
         key1 = name[:cut]
         key2 = name[cut:]
         key3 = pos1 / numcut
         key4 = pos1 % numcut
         key5 = pos2 / numcut
         key6 = pos2 % numcut    
         pealign +=2     
         try: 
            if key6 in peres[key0][key1][key2][key3][key4][key5]:
               continue
            else:
               peres[key0][key1][key2][key3][key4][key5].append(key6)
         except:
            x[45678987654]         
         
         o.write(x+y)
         peunique+=2         
      # OTERWISE TREAT AS SE   
      elif opt.mode == 'ps':
         for temp in [l, m]:      
            if temp[2] == '*':
               continue   
                           
            if temp[1] == '0':
               key1 = "F"
               pos = temp[3]   
            elif temp[1] == '16':   
               key1 = "R"
               pos = int(temp[3])+len(temp[9])-1
            else:
               continue
            sealign +=1                  
            key2 = temp[2][:cut]
            key3 = temp[2][cut:]
            key4 = int(pos) / numcut
            key5 = int(pos) % numcut
                
            if key5 in seres[key1][key2][key3][key4]:
               continue
            else:
               seres[key1][key2][key3][key4].append(key5)
            o.write(x)
            seunique+=1
            
      elif opt.mode == 'pb':
         good = 0
         for temp in [l, m]:      
            if temp[2] == '*':
               continue   
                           
            if temp[1] == '0':
               key1 = "F"
               pos = temp[3]   
            elif temp[1] == '16':   
               key1 = "R"
               pos = int(temp[3])+len(temp[9])-1
            else:
               continue
            oealign +=1                  
            key2 = temp[2][:cut]
            key3 = temp[2][cut:]
            key4 = int(pos) / numcut
            key5 = int(pos) % numcut
                
            if key5 in seres[key1][key2][key3][key4]:
               continue
            else:
               seres[key1][key2][key3][key4].append(key5)
               good +=1
         if good != 0:   
            o.write(x+y)
            oeunique+=good

 
else:
   seres = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
   #              F/R               key1                 key2              posmul   posmod
   while True:
      l = f.readline()
      if l == "":
         break
      if l[0] == '@':
         continue
      x = l.split('\t')
      reads +=1
    
      if x[2] == '*':
         continue   
                     
      if x[1] == '0':
         key1 = "F"
         pos = x[3]   
      elif x[1] == '16':   
         key1 = "R"
         pos = int(x[3])+len(x[9])-1
      else:
         continue
      sealign +=1                  
      key2 = x[2][:cut]
      key3 = x[2][cut:]
      key4 = int(pos) / numcut
      key5 = int(pos) % numcut
            
      if key5 in seres[key1][key2][key3][key4]:
         continue
      else:
         seres[key1][key2][key3][key4].append(key5)
      o.write(l)
      seunique+=1      
        
o.close()
#total = sum([len(refs[x]) for x in refs])
if opt.mode not in ['p', 'pb']:
   sepercent = 100*seunique/float(sealign)
   print '\t'.join([opt.sam+"\tSE", str(seunique), str(sealign), str(round(sepercent,3)), str(reads)])
if opt.mode == 'pb':
   oepercent = 100*oeunique/float(oealign)
   print '\t'.join([opt.sam+"\tPSE", str(oeunique), str(oealign), str(round(oepercent,3)), str(reads)])   
if opt.mode in ['p', 'pb', 'ps']:
   pepercent = 100*peunique/float(pealign)
   print '\t'.join([opt.sam+"\tPE", str(peunique), str(pealign), str(round(pepercent,3)), str(reads)])
allpercent = 100*(seunique+peunique+oeunique)/float(pealign+sealign+oealign)
if opt.mode not in ['s','p']:
   print '\t'.join([opt.sam+"\tALL", str(peunique+seunique+oeunique), str(pealign+sealign+oealign), str(round(allpercent,3)), str(reads)])
 
#print "For file : "+opt.sam
#print "The total number of aligned reads is " + str(align)
#print "The percent unique reads is " +str(round(percent,3))
#print time.ctime()

