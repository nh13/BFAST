bin_PROGRAMS = solid2fastq

solid2fastq_SOURCES = \
					  ../bfast/BError.c	../bfast/BError.h \
					  ../bfast/RGIndex.c  ../bfast/RGIndex.h \
					  ../bfast/BLib.c ../bfast/BLib.h \
					  ../bfast/RGBinary.c ../bfast/RGBinary.h \
					  ../bfast/RGRanges.c ../bfast/RGRanges.h \
					  ../bfast/RGMatch.c ../bfast/RGMatch.h \
					  solid2fastq.c

solid2fastq_LDADD = -lm -lz

dist_bin_SCRIPTS = bfast.submit.pl solid2fastq.pl qseq2fastq.pl
