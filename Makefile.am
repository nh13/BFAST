bin_PROGRAMS = bfast 

bfast_SOURCES = \
				AlignedEnd.c AlignedEnd.h \
				AlignedEntry.c AlignedEntry.h \
				AlignedRead.c AlignedRead.h \
				AlignedReadConvert.c AlignedReadConvert.h \
				BError.c BError.h \
				BLib.c BLib.h \
				BLibDefinitions.h \
				RGBinary.c RGBinary.h \
				RGIndex.c RGIndex.h \
				RGIndexAccuracy.c RGIndexAccuracy.h \
				RGIndexExons.c RGIndexExons.h \
				RGIndexLayout.c RGIndexLayout.h \
				RGMatch.c RGMatch.h \
				RGMatches.c RGMatches.h \
				RGRanges.c RGRanges.h \
				RGReads.c RGReads.h \
				ScoringMatrix.c ScoringMatrix.h \
				FindMatches.c FindMatches.h \
				MatchesReadInputFiles.c MatchesReadInputFiles.h \
				RunAligner.c RunAligner.h \
				Align.c Align.h \
				AlignNTSpace.c AlignNTSpace.h \
				AlignColorSpace.c AlignColorSpace.h \
				InputOutputToFiles.c InputOutputToFiles.h \
				FilterAlignments.c FilterAlignments.h \
				BfastFasta2BRG.c BfastFasta2BRG.h \
				BfastIndex.c BfastIndex.h \
				BfastMatch.c BfastMatch.h \
				BfastLocalAlign.c BfastLocalAlign.h \
				BfastPostProcess.c BfastPostProcess.h \
				BfastBAFConvert.c \
				BfastHeader.c \
				BfastBMFConvert.c \
				BfastBRG2Fasta.c \
				Main.c Main.h

bfast_LDADD = -lm -lz

EXTRA_DIST = LICENSE \
			 autogen.sh \
			 manual/bfast-book.pdf \
			 scripts/solid2fastq.pl \
			 scripts/qseq2fastq.pl \
			 scripts/bfast.submit.pl \
			 tests/test.definitions.sh \
			 tests/data/data.tar.gz \
			 tests/test.initialize.sh \
			 tests/test.fasta2brg.sh \
			 tests/test.index.sh \
			 tests/test.match.sh \
			 tests/test.localalign.sh \
			 tests/test.postprocess.sh \
			 tests/test.diff.sh \
			 tests/test.cleanup.sh

SUBDIRS = butil scripts tests
