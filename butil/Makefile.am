bin_PROGRAMS = balignmentscoredistribution \
			   balignsim \
			   bevalsim \
			   bexonify \
			   bgeneratereads \
			   berrordistribution \
			   bindexdist \
			   bindexhist \
			   bindexhistchr \
			   bmergesorted \
			   bpairedenddist \
			   brepeat \
			   bsort \
			   btestindexes \
			   btranslocations


balignmentscoredistribution_SOURCES = \
									  ../bfast/BError.c	../bfast/BError.h \
									  ../bfast/RGIndex.c	../bfast/RGIndex.h \
									  ../bfast/BLib.c	../bfast/BLib.h \
									  ../bfast/RGBinary.c ../bfast/RGBinary.h \
									  ../bfast/RGRanges.c ../bfast/RGRanges.h \
									  ../bfast/AlignedRead.c	../bfast/AlignedRead.h \
									  ../bfast/AlignedReadConvert.c	../bfast/AlignedReadConvert.h \
									  ../bfast/AlignedEnd.c	../bfast/AlignedEnd.h \
									  ../bfast/AlignedEntry.c	../bfast/AlignedEntry.h \
									  balignmentscoredistribution.c balignmentscoredistribution.h

balignmentscoredistribution_LDADD = -lz

balignsim_SOURCES = \
					SimRead.c	SimRead.h \
					../bfast/BError.c	../bfast/BError.h \
					../bfast/RGIndexExons.c  ../bfast/RGIndexExons.h \
					../bfast/RGIndex.c	../bfast/RGIndex.h \
					../bfast/BLib.c	../bfast/BLib.h \
					../bfast/RGBinary.c ../bfast/RGBinary.h \
					../bfast/RGMatch.c ../bfast/RGMatch.h \
					../bfast/RGMatches.c	../bfast/RGMatches.h \
					../bfast/RGRanges.c ../bfast/RGRanges.h \
					../bfast/AlignedRead.c	../bfast/AlignedRead.h \
					../bfast/AlignedEnd.c	../bfast/AlignedEnd.h \
					../bfast/RGIndexAccuracy.c	../bfast/RGIndexAccuracy.h \
					../bfast/AlignedEntry.c	../bfast/AlignedEntry.h \
					../bfast/RunAligner.c	../bfast/RunAligner.h \
					../bfast/ScoringMatrix.c	../bfast/ScoringMatrix.h \
					../bfast/Align.c	../bfast/Align.h \
					../bfast/AlignColorSpace.c	../bfast/AlignColorSpace.h \
					../bfast/AlignNTSpace.c	../bfast/AlignNTSpace.h \
					balignsim.c	balignsim.h

balignsim_LDADD = -lz

bevalsim_SOURCES = \
				   ../bfast/BError.c	../bfast/BError.h \
				   ../bfast/RGIndexExons.c  ../bfast/RGIndexExons.h \
				   ../bfast/RGIndex.c	../bfast/RGIndex.h \
				   ../bfast/BLib.c	../bfast/BLib.h \
				   ../bfast/RGBinary.c ../bfast/RGBinary.h \
				   ../bfast/RGRanges.c ../bfast/RGRanges.h \
				   ../bfast/RGMatches.c ../bfast/RGMatches.h \
				   ../bfast/RGMatch.c ../bfast/RGMatch.h \
				   ../bfast/AlignedRead.c	../bfast/AlignedRead.h \
				   ../bfast/AlignedEnd.c	../bfast/AlignedEnd.h \
				   ../bfast/AlignedEntry.c	../bfast/AlignedEntry.h \
				   bevalsim.c	bevalsim.h

bevalsim_LDADD = -lz 

bexonify_SOURCES = \
				   ../bfast/BError.c	../bfast/BError.h \
				   ../bfast/RGIndexExons.c  ../bfast/RGIndexExons.h \
				   ../bfast/RGIndex.c	../bfast/RGIndex.h \
				   ../bfast/BLib.c	../bfast/BLib.h \
				   ../bfast/RGBinary.c ../bfast/RGBinary.h \
				   ../bfast/RGMatch.c ../bfast/RGMatch.h \
				   ../bfast/RGRanges.c ../bfast/RGRanges.h \
				   bexonify.c	bexonify.h

bexonify_LDADD = -lz

bgeneratereads_SOURCES = \
						 ../bfast/BError.c	../bfast/BError.h \
						 ../bfast/RGIndexExons.c  ../bfast/RGIndexExons.h \
						 ../bfast/RGIndex.c	../bfast/RGIndex.h \
						 ../bfast/BLib.c	../bfast/BLib.h \
						 ../bfast/RGBinary.c ../bfast/RGBinary.h \
						 ../bfast/RGMatch.c ../bfast/RGMatch.h \
						 ../bfast/RGRanges.c ../bfast/RGRanges.h \
						 SimRead.c	SimRead.h \
						 bgeneratereads.c	bgeneratereads.h 

bgeneratereads_LDADD = -lz  

berrordistribution_SOURCES = \
							 ../bfast/AlignedEnd.c	../bfast/AlignedEnd.h \
							 ../bfast/AlignedEntry.c	../bfast/AlignedEntry.h \
							 ../bfast/AlignedRead.c	../bfast/AlignedRead.h \
							 ../bfast/BError.c	../bfast/BError.h \
							 ../bfast/RGIndexExons.c  ../bfast/RGIndexExons.h \
							 ../bfast/RGIndex.c	../bfast/RGIndex.h \
							 ../bfast/BLib.c	../bfast/BLib.h \
							 ../bfast/RGBinary.c ../bfast/RGBinary.h \
							 ../bfast/RGMatch.c ../bfast/RGMatch.h \
							 ../bfast/RGRanges.c ../bfast/RGRanges.h \
							 berrordistribution.c	berrordistribution.h 

berrordistribution_LDADD = -lz  

bindexdist_SOURCES = \
					 ../bfast/BError.c	../bfast/BError.h \
					 ../bfast/RGIndexExons.c  ../bfast/RGIndexExons.h \
					 ../bfast/RGIndex.c	../bfast/RGIndex.h \
					 ../bfast/BLib.c	../bfast/BLib.h \
					 ../bfast/RGBinary.c ../bfast/RGBinary.h \
					 ../bfast/RGMatch.c ../bfast/RGMatch.h \
					 ../bfast/RGRanges.c ../bfast/RGRanges.h \
					 ../bfast/RGReads.c	../bfast/RGReads.h \
					 bindexdist.c	bindexdist.h

bindexdist_LDADD = -lz

bindexhistchr_SOURCES = \
					 ../bfast/BError.c	../bfast/BError.h \
					 ../bfast/RGIndexExons.c  ../bfast/RGIndexExons.h \
					 ../bfast/RGIndex.c	../bfast/RGIndex.h \
					 ../bfast/BLib.c	../bfast/BLib.h \
					 ../bfast/RGBinary.c ../bfast/RGBinary.h \
					 ../bfast/RGMatch.c ../bfast/RGMatch.h \
					 ../bfast/RGRanges.c ../bfast/RGRanges.h \
					 ../bfast/RGReads.c	../bfast/RGReads.h \
					 bindexhistchr.c	bindexhistchr.h

bindexhistchr_LDADD = -lz


bindexhist_SOURCES = \
					 ../bfast/BError.c	../bfast/BError.h \
					 ../bfast/RGIndexExons.c  ../bfast/RGIndexExons.h \
					 ../bfast/RGIndex.c	../bfast/RGIndex.h \
					 ../bfast/BLib.c	../bfast/BLib.h \
					 ../bfast/RGBinary.c ../bfast/RGBinary.h \
					 ../bfast/RGMatch.c ../bfast/RGMatch.h \
					 ../bfast/RGRanges.c ../bfast/RGRanges.h \
					 ../bfast/RGReads.c	../bfast/RGReads.h \
					 bindexhist.c	bindexhist.h

bindexhist_LDADD = -lz

bmergesorted_SOURCES = \
					   ../bfast/BError.c	../bfast/BError.h \
					   ../bfast/RGIndexExons.c  ../bfast/RGIndexExons.h \
					   ../bfast/RGIndex.c	../bfast/RGIndex.h \
					   ../bfast/BLib.c	../bfast/BLib.h \
					   ../bfast/RGBinary.c ../bfast/RGBinary.h \
					   ../bfast/RGMatch.c ../bfast/RGMatch.h \
					   ../bfast/RGMatches.c ../bfast/RGMatches.h \
					   ../bfast/RGRanges.c ../bfast/RGRanges.h \
					   ../bfast/RGReads.c	../bfast/RGReads.h \
					   ../bfast/AlignedRead.c	../bfast/AlignedRead.h \
					   ../bfast/AlignedEnd.c	../bfast/AlignedEnd.h \
					   ../bfast/AlignedEntry.c	../bfast/AlignedEntry.h \
					   bmergesorted.c	bmergesorted.h

bmergesorted_LDADD = -lz

bpairedenddist_SOURCES = \
						 ../bfast/BError.c	../bfast/BError.h \
						 ../bfast/RGIndexExons.c  ../bfast/RGIndexExons.h \
						 ../bfast/RGIndex.c	../bfast/RGIndex.h \
						 ../bfast/BLib.c	../bfast/BLib.h \
						 ../bfast/RGBinary.c ../bfast/RGBinary.h \
						 ../bfast/RGMatch.c ../bfast/RGMatch.h \
						 ../bfast/RGMatches.c ../bfast/RGMatches.h \
						 ../bfast/RGRanges.c ../bfast/RGRanges.h \
						 ../bfast/RGReads.c	../bfast/RGReads.h \
						 ../bfast/AlignedRead.c	../bfast/AlignedRead.h \
						 ../bfast/AlignedEnd.c	../bfast/AlignedEnd.h \
						 ../bfast/AlignedEntry.c	../bfast/AlignedEntry.h \
						 bpairedenddist.c	bpairedenddist.h

bpairedenddist_LDADD = -lz

brepeat_SOURCES = \
				  ../bfast/BError.c	../bfast/BError.h \
				  ../bfast/RGIndexExons.c  ../bfast/RGIndexExons.h \
				  ../bfast/RGIndex.c	../bfast/RGIndex.h \
				  ../bfast/BLib.c	../bfast/BLib.h \
				  ../bfast/RGBinary.c ../bfast/RGBinary.h \
				  ../bfast/RGMatch.c ../bfast/RGMatch.h \
				  ../bfast/RGRanges.c ../bfast/RGRanges.h \
				  brepeat.c	brepeat.h

brepeat_LDADD = -lz 

bsort_SOURCES = \
				../bfast/AlignedEnd.c	../bfast/AlignedEnd.h \
				../bfast/AlignedEntry.c	../bfast/AlignedEntry.h \
				../bfast/AlignedRead.c	../bfast/AlignedRead.h \
				../bfast/BError.c	../bfast/BError.h \
				../bfast/RGIndexExons.c  ../bfast/RGIndexExons.h \
				../bfast/RGIndex.c	../bfast/RGIndex.h \
				../bfast/BLib.c	../bfast/BLib.h \
				../bfast/RGBinary.c ../bfast/RGBinary.h \
				../bfast/RGMatch.c ../bfast/RGMatch.h \
				../bfast/RGMatches.c ../bfast/RGMatches.h \
				../bfast/RGRanges.c ../bfast/RGRanges.h \
				../bfast/RGReads.c	../bfast/RGReads.h \
				bsort.c	bsort.h

bsort_LDADD = -lz

btestindexes_SOURCES = \
					   ../bfast/RGRanges.c ../bfast/RGRanges.h \
					   ../bfast/RGBinary.c ../bfast/RGBinary.h \
					   ../bfast/RGIndex.c	../bfast/RGIndex.h \
					   ../bfast/RGIndexAccuracy.c	../bfast/RGIndexAccuracy.h \
					   ../bfast/BError.c	../bfast/BError.h \
					   ../bfast/BLib.c	../bfast/BLib.h \
					   btestindexes.c	btestindexes.h

btestindexes_LDADD = -lz

btranslocations_SOURCES = \
						  ../bfast/AlignedEnd.c	../bfast/AlignedEnd.h \
						  ../bfast/AlignedEntry.c	../bfast/AlignedEntry.h \
						  ../bfast/AlignedRead.c	../bfast/AlignedRead.h \
						  ../bfast/BError.c	../bfast/BError.h \
						  ../bfast/RGIndexExons.c  ../bfast/RGIndexExons.h \
						  ../bfast/RGIndex.c	../bfast/RGIndex.h \
						  ../bfast/BLib.c	../bfast/BLib.h \
						  ../bfast/RGBinary.c ../bfast/RGBinary.h \
						  ../bfast/RGMatch.c ../bfast/RGMatch.h \
						  ../bfast/RGMatches.c ../bfast/RGMatches.h \
						  ../bfast/RGRanges.c ../bfast/RGRanges.h \
						  ../bfast/RGReads.c	../bfast/RGReads.h \
						  btranslocations.c btranslocations.h

btranslocations_LDADD = -lz 
