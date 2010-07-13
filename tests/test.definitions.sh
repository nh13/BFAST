#!/bin/sh

error()
{
  echo $1
  exit 1
}

# Stores all the definitions required for the tests
CMD_PREFIX="../bfast/";
#CMD_PREFIX="valgrind --leak-check=yes --log-file=pipe.valgrind ../bfast/";
REF_ID_CC="corner.cases";
OUTPUT_ID_NT="phi.x.174";
OUTPUT_ID_CC_NT="corner.cases.nt";
OUTPUT_ID_CS="DH10B";
OUTPUT_ID_CC_CS="corner.cases.cs";
NUM_THREADS="1";
DATA_DIR="data/";
INPUT_DIR="data/";
OUTPUT_DIR="output/";
TMP_DIR="tmp/";
S2F="../scripts/solid2fastq"

# Try to find any version of md5sum
MD5BIN="md5sum"
MLIST="Please, contact: bfast-help@lists.sourceforge.net"
[ `which $MD5BIN` ] || MD5BIN="gmd5sum"
[ `which $MD5BIN` ] || error "I can't find md5sum in your system. $MLIST"
