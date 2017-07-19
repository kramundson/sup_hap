#!/bin/bash
# Kirk Amundson
# potato instance setup
# Install software packages: sratoolkit, aspera-connect, samtools, bwa, bin-by-sam
# Also get wrapper scripts from comai lab website

# This should be used once on an AMI and the AMI stored thereafter to prevent excessive 
# accumulation of gray hairs

# Keep all software in spud-software and symlink to /usr/local/bin maybe?
cd
mkdir ~/spud-software
cd ~/spud-software

# install ascp
wget http://download.asperasoft.com/download/sw/ascp-client/3.5.4/ascp-install-3.5.4.102989-linux-64.sh
chmod +x ascp-install-3.5.4.102989-linux-64.sh
sudo ./ascp-install-3.5.4.102989-linux-64.sh

# install bwa
cd ~/spud-software
curl -L https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.15.tar.bz2/download > bwa-0.7.15.tar.bz2
tar xjvf bwa-0.7.15.tar.bz2
cd bwa-0.7.15
make
sudo cp bwa /usr/local/bin
echo 'export PATH=$PATH:/usr/local/bin' >> $HOME/.profile
source ~/.profile

# install samtools
cd
sudo apt-get -y install samtools

# install sratoolkit
cd
wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -zxvf sratoolkit.tar.gz
echo 'export PATH=$PATH:$HOME/sratoolkit.2.8.2-1-ubuntu64/bin' >> $HOME/.profile
source ~/.profile
which fastq-dump # should return a path where fastq dump would likely live
fastq-dump -X 5 -Z SRR390728

# install bwa-doall
mkdir ~/bin
cd ~/bin
wget http://comailab.genomecenter.ucdavis.edu/download.php?f=f/ff/BWA-DoAll-1.0.tgz
tar xzvf download.php?f=f%2Fff%2FBWA-DoAll-1.0.tgz

# install bin by sam
cd ~/bin
wget http://comailab.genomecenter.ucdavis.edu/images/6/66/Bin-by-Sam-tool.tar.gz
tar xzvf Bin-by-Sam-tool.tar.gz

# install run-mpileup.py and mpileup-parser-v2.py
cd ~/bin
wget http://comailab.genomecenter.ucdavis.edu/images/6/6b/Mpileup-tools-2.0.tgz
tar xzvf Mpileup-tools-2.0.tgz

# should produce the lines below and nothing else, but uncommented:

# Read 5 spots for SRR390728
# Written 5 spots for SRR390728
# @SRR390728.1 1 length=72
# CATTCTTCACGTAGTTCTCGAGCCTTGGTTTTCAGCGATGGAGAATGACTTTGACAAGCTGAGAGAAGNTNC
# +SRR390728.1 1 length=72
# ;;;;;;;;;;;;;;;;;;;;;;;;;;;9;;665142;;;;;;;;;;;;;;;;;;;;;;;;;;;;;96&&&&(
# @SRR390728.2 2 length=72
# AAGTAGGTCTCGTCTGTGTTTTCTACGAGCTTGTGTTCCAGCTGACCCACTCCCTGGGTGGGGGGACTGGGT
# +SRR390728.2 2 length=72
# ;;;;;;;;;;;;;;;;;4;;;;3;393.1+4&&5&&;;;;;;;;;;;;;;;;;;;;;<9;<;;;;;464262
# @SRR390728.3 3 length=72
# CCAGCCTGGCCAACAGAGTGTTACCCCGTTTTTACTTATTTATTATTATTATTTTGAGACAGAGCATTGGTC
# +SRR390728.3 3 length=72
# -;;;8;;;;;;;,*;;';-4,44;,:&,1,4'./&19;;;;;;669;;99;;;;;-;3;2;0;+;7442&2/
# @SRR390728.4 4 length=72
# ATAAAATCAGGGGTGTTGGAGATGGGATGCCTATTTCTGCACACCTTGGCCTCCCAAATTGCTGGGATTACA
# +SRR390728.4 4 length=72
# 1;;;;;;,;;4;3;38;8%&,,;)*;1;;,)/%4+,;1;;);;;;;;;4;(;1;;;;24;;;;41-444//0
# @SRR390728.5 5 length=72
# TTAAGAAATTTTTGCTCAAACCATGCCCTAAAGGGTTCTGTAATAAATAGGGCTGGGAAAACTGGCAAGCCA
# +SRR390728.5 5 length=72
# ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;9445552;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;446662
