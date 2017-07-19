#!/bin/bash
# fastq-dump-edwards.sh
# Kirk Amundson
# Wrapper script for calling fastq-dump with recommended settings from Edwards lab
# 19 July 2017

# Usage: fasterq-dump.sh </path/to/.sra>

fastq-dump --gzip --skip-technical --readids --dumpbase --clip --split-files $1
