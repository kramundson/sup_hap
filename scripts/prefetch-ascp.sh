#!/bin/bash
# Kirk Amundson
# prefetch-ascp.sh
# 19 Juy 2017

# wrapper for calling SRA prefetch without entering in ascp path every time
# Also runs vdb-validate on each downloaded .sra and writes out to file

# Usage: prefetch-ascp.sh <SRR12345678>

prefetch --ascp-path '/home/ubuntu/.aspera/connect/bin/ascp|/home/ubuntu/.aspera/connect/etc/asperawebid_dsa.openssh' $1 
vdb-validate $1 > $1'-vdbvalidate.out' 2> $1'-vdbvalidate.err'
