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
									  ../BError.c	../BError.h \
									  ../RGIndex.c	../RGIndex.h \
									  ../BLib.c	../BLib.h \
									  ../RGBinary.c ../RGBinary.h \
									  ../RGRanges.c ../RGRanges.h \
									  ../AlignedRead.c	../AlignedRead.h \
									  ../AlignedReadConvert.c	../AlignedReadConvert.h \
									  ../AlignedEnd.c	../AlignedEnd.h \
									  ../AlignedEntry.c	../AlignedEntry.h \
									  balignmentscoredistribution.c balignmentscoredistribution.h

balignmentscoredistribution_LDADD = -lz

balignsim_SOURCES = \
					SimRead.c	SimRead.h \
					../BError.c	../BError.h \
					../RGIndexExons.c  ../RGIndexExons.h \
					../RGIndex.c	../RGIndex.h \
					../BLib.c	../BLib.h \
					../RGBinary.c ../RGBinary.h \
					../RGMatch.c ../RGMatch.h \
					../RGMatches.c	../RGMatches.h \
					../RGRanges.c ../RGRanges.h \
					../AlignedRead.c	../AlignedRead.h \
					../AlignedEnd.c	../AlignedEnd.h \
					../RGIndexAccuracy.c	../RGIndexAccuracy.h \
					../AlignedEntry.c	../AlignedEntry.h \
					../RunAligner.c	../RunAligner.h \
					../ScoringMatrix.c	../ScoringMatrix.h \
					../Align.c	../Align.h \
					../AlignColorSpace.c	../AlignColorSpace.h \
					../AlignNTSpace.c	../AlignNTSpace.h \
					balignsim.c	balignsim.h

balignsim_LDADD = -lz

bevalsim_SOURCES = \
				   ../BError.c	../BError.h \
				   ../RGIndexExons.c  ../RGIndexExons.h \
				   ../RGIndex.c	../RGIndex.h \
				   ../BLib.c	../BLib.h \
				   ../RGBinary.c ../RGBinary.h \
				   ../RGRanges.c ../RGRanges.h \
				   ../RGMatches.c ../RGMatches.h \
				   ../RGMatch.c ../RGMatch.h \
				   ../AlignedRead.c	../AlignedRead.h \
				   ../AlignedEnd.c	../AlignedEnd.h \
				   ../AlignedEntry.c	../AlignedEntry.h \
				   bevalsim.c	bevalsim.h

bevalsim_LDADD = -lz 

bexonify_SOURCES = \
				   ../BError.c	../BError.h \
				   ../RGIndexExons.c  ../RGIndexExons.h \
				   ../RGIndex.c	../RGIndex.h \
				   ../BLib.c	../BLib.h \
				   ../RGBinary.c ../RGBinary.h \
				   ../RGMatch.c ../RGMatch.h \
				   ../RGRanges.c ../RGRanges.h \
				   bexonify.c	bexonify.h

bexonify_LDADD = -lz

bgeneratereads_SOURCES = \
						 ../BError.c	../BError.h \
						 ../RGIndexExons.c  ../RGIndexExons.h \
						 ../RGIndex.c	../RGIndex.h \
						 ../BLib.c	../BLib.h \
						 ../RGBinary.c ../RGBinary.h \
						 ../RGMatch.c ../RGMatch.h \
						 ../RGRanges.c ../RGRanges.h \
						 SimRead.c	SimRead.h \
						 bgeneratereads.c	bgeneratereads.h 

bgeneratereads_LDADD = -lz  

berrordistribution_SOURCES = \
							 ../AlignedEnd.c	../AlignedEnd.h \
							 ../AlignedEntry.c	../AlignedEntry.h \
							 ../AlignedRead.c	../AlignedRead.h \
							 ../BError.c	../BError.h \
							 ../RGIndexExons.c  ../RGIndexExons.h \
							 ../RGIndex.c	../RGIndex.h \
							 ../BLib.c	../BLib.h \
							 ../RGBinary.c ../RGBinary.h \
							 ../RGMatch.c ../RGMatch.h \
							 ../RGRanges.c ../RGRanges.h \
							 berrordistribution.c	berrordistribution.h 

berrordistribution_LDADD = -lz  

bindexdist_SOURCES = \
					 ../BError.c	../BError.h \
					 ../RGIndexExons.c  ../RGIndexExons.h \
					 ../RGIndex.c	../RGIndex.h \
					 ../BLib.c	../BLib.h \
					 ../RGBinary.c ../RGBinary.h \
					 ../RGMatch.c ../RGMatch.h \
					 ../RGRanges.c ../RGRanges.h \
					 ../RGReads.c	../RGReads.h \
					 bindexdist.c	bindexdist.h

bindexdist_LDADD = -lz

bindexhistchr_SOURCES = \
					 ../BError.c	../BError.h \
					 ../RGIndexExons.c  ../RGIndexExons.h \
					 ../RGIndex.c	../RGIndex.h \
					 ../BLib.c	../BLib.h \
					 ../RGBinary.c ../RGBinary.h \
					 ../RGMatch.c ../RGMatch.h \
					 ../RGRanges.c ../RGRanges.h \
					 ../RGReads.c	../RGReads.h \
					 bindexhistchr.c	bindexhistchr.h

bindexhistchr_LDADD = -lz


bindexhist_SOURCES = \
					 ../BError.c	../BError.h \
					 ../RGIndexExons.c  ../RGIndexExons.h \
					 ../RGIndex.c	../RGIndex.h \
					 ../BLib.c	../BLib.h \
					 ../RGBinary.c ../RGBinary.h \
					 ../RGMatch.c ../RGMatch.h \
					 ../RGRanges.c ../RGRanges.h \
					 ../RGReads.c	../RGReads.h \
					 bindexhist.c	bindexhist.h

bindexhist_LDADD = -lz

bmergesorted_SOURCES = \
					   ../BError.c	../BError.h \
					   ../RGIndexExons.c  ../RGIndexExons.h \
					   ../RGIndex.c	../RGIndex.h \
					   ../BLib.c	../BLib.h \
					   ../RGBinary.c ../RGBinary.h \
					   ../RGMatch.c ../RGMatch.h \
					   ../RGMatches.c ../RGMatches.h \
					   ../RGRanges.c ../RGRanges.h \
					   ../RGReads.c	../RGReads.h \
					   ../AlignedRead.c	../AlignedRead.h \
					   ../AlignedEnd.c	../AlignedEnd.h \
					   ../AlignedEntry.c	../AlignedEntry.h \
					   bmergesorted.c	bmergesorted.h

bmergesorted_LDADD = -lz

bpairedenddist_SOURCES = \
						 ../BError.c	../BError.h \
						 ../RGIndexExons.c  ../RGIndexExons.h \
						 ../RGIndex.c	../RGIndex.h \
						 ../BLib.c	../BLib.h \
						 ../RGBinary.c ../RGBinary.h \
						 ../RGMatch.c ../RGMatch.h \
						 ../RGMatches.c ../RGMatches.h \
						 ../RGRanges.c ../RGRanges.h \
						 ../RGReads.c	../RGReads.h \
						 ../AlignedRead.c	../AlignedRead.h \
						 ../AlignedEnd.c	../AlignedEnd.h \
						 ../AlignedEntry.c	../AlignedEntry.h \
						 bpairedenddist.c	bpairedenddist.h

bpairedenddist_LDADD = -lz

brepeat_SOURCES = \
				  ../BError.c	../BError.h \
				  ../RGIndexExons.c  ../RGIndexExons.h \
				  ../RGIndex.c	../RGIndex.h \
				  ../BLib.c	../BLib.h \
				  ../RGBinary.c ../RGBinary.h \
				  ../RGMatch.c ../RGMatch.h \
				  ../RGRanges.c ../RGRanges.h \
				  brepeat.c	brepeat.h

brepeat_LDADD = -lz 

bsort_SOURCES = \
				../AlignedEnd.c	../AlignedEnd.h \
				../AlignedEntry.c	../AlignedEntry.h \
				../AlignedRead.c	../AlignedRead.h \
				../BError.c	../BError.h \
				../RGIndexExons.c  ../RGIndexExons.h \
				../RGIndex.c	../RGIndex.h \
				../BLib.c	../BLib.h \
				../RGBinary.c ../RGBinary.h \
				../RGMatch.c ../RGMatch.h \
				../RGMatches.c ../RGMatches.h \
				../RGRanges.c ../RGRanges.h \
				../RGReads.c	../RGReads.h \
				bsort.c	bsort.h

bsort_LDADD = -lz

btestindexes_SOURCES = \
					   ../RGRanges.c ../RGRanges.h \
					   ../RGBinary.c ../RGBinary.h \
					   ../RGIndex.c	../RGIndex.h \
					   ../RGIndexAccuracy.c	../RGIndexAccuracy.h \
					   ../BError.c	../BError.h \
					   ../BLib.c	../BLib.h \
					   btestindexes.c	btestindexes.h

btestindexes_LDADD = -lz

btranslocations_SOURCES = \
						  ../AlignedEnd.c	../AlignedEnd.h \
						  ../AlignedEntry.c	../AlignedEntry.h \
						  ../AlignedRead.c	../AlignedRead.h \
						  ../BError.c	../BError.h \
						  ../RGIndexExons.c  ../RGIndexExons.h \
						  ../RGIndex.c	../RGIndex.h \
						  ../BLib.c	../BLib.h \
						  ../RGBinary.c ../RGBinary.h \
						  ../RGMatch.c ../RGMatch.h \
						  ../RGMatches.c ../RGMatches.h \
						  ../RGRanges.c ../RGRanges.h \
						  ../RGReads.c	../RGReads.h \
						  btranslocations.c btranslocations.h

btranslocations_LDADD = -lz 
