# readcut.py
# Kirk Amundson
# Hard cut fastq-formatted reads to specified number of nucleotides
# Note: requires uninterleaved file. I use this to sample randomly from the forward reads
# of a 150PE .fq and cut readlength to 50 for feeding into alignment + bin-by-sam

# USAGE: python script.py reads.fq length

import sys
# import random

# def write_random_records(fq, N=100000):
#     """get N random headers from a fastq file without reading whole thing into memory"""
#     records = sum(1 for _ in open(fq)) / 4
#     print(records)
#     rand_records = sorted([random.randint(0,records-1) for _ in xrange(N)])
#     fh = open(fq)
#     pref = fq.split('.')[0]
#     sub = open(pref+"_subset"+N+"_bp.fq", 'w')
#     rec_no = - 1
#     
#     for rr in rand_records:
#         while rec_no < rr:
#             rec_no += 1
#             name = fh.readline()
#             seq = fh.readline()
#             plus = fh.readline()
#             qual = fh.readline()
#             scut = seq[:keep]
#             qcut = qual[:keep]
#             out.write(name+scut+"\n"+plus+qcut+"\n")

# write_random_records(sys.argv[1])

pref = sys.argv[1].split('.')[0]
print pref
keep = int(sys.argv[2])
fh = open(sys.argv[1],'r')
out = open(pref+"_"+str(keep)+".fq", 'w')

while True:
    # read four lines at a time
    name = fh.readline()
    if name == "":
        break
    seq = fh.readline().rstrip()
    plus = fh.readline()
    qual = fh.readline().rstrip()
    
    # scut = seq[:keep-1] # bug keeps one shorter than asked for.
    scut = seq[:keep] # all better
    qcut = qual[:keep] # same for qual
    #qcut = qual[:keep-1] # comment out and fix 4/7/16
    out.write(name+scut+"\n"+plus+qcut+"\n")

fh.close()
out.close()
