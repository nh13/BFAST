#!/usr/bin/perl
use strict;
use Getopt::Long;
# TODO 
# Clean up this code
# Comment the code

my $iniFileName;

#Command-line options
GetOptions("iniFileName=s"=>\$iniFileName)
	or die("Error parsing command-line options.  Terminating!\n");
if(!defined($iniFileName)) {
	die("Option iniFileName was not given");
}

Main($iniFileName);

# TODO 
sub Main {
	my $iniFileName = shift;
	local *FH;

	# Create data object to hold all options
	my %data = ();

	# open the file for reading
	open(FH, "$iniFileName") || die("Error.  Could not open $iniFileName for reading.  Terminating!\n"); 

	####################################################################
	# Parse the file 
	# Read in by line
	####################################################################
	while(defined(my $line = <FH>)) {
		chomp($line);
		if($line =~ m/^#/) {
			# Ignore comment line 
		}
		elsif($line =~ m/^([A-Z]*)="(.*)"/) {
			# Read in the name and a data pair
			my $name = $1;
			my $data = $2;
			#print("name:$name\tdata:$data\n");
			$data{$name} = $data;
		}
		else {
			die("Error.  Could not recognize line in $iniFileName:\n$line\nTerminating!\n");
		}
	}

	# close the file
	close(FH);

	# Check that all options are specified
	my @options = (
		"GENERATETREESSH",
		"GENERATEINDEXESSH",
		"FINDMATCHESSH",
		"ALIGNSH",
		"COMMANDDIR",
		"OUTPUTDIR",
		"OUTPUTID",
		"RGLISTANDLENGTHFILENAME",
		"REGIONS",
		"INDEXMATCHLENGTH",
		"TREEMATCHLENGTH",
		"GAPS",
		"OFFSETS",
		"READSFILENAME",
		"NUMREADS",
		"NUMPROCESSES",
		"NUMMISMATCHES",
		"NUMINSERTIONS",
		"NUMDELETIONS",
		"PAIREDEND",
		"SCORINGMATRIXFILENAME",
		"ALGORITHM",
		"ALIGNOFFSET");
	for(my $i=0;$i<scalar(@options);$i++) {
		if(!defined($data{$options[$i]})){ 
			die("Error.  ".$options[$i]." was not found in $iniFileName.  Terminating!\n");
		}
	}

	# TODO - check that all paths have an '/' at the end


	# Read in chromsomes with list and lengths
	# open the file for reading
	open(FH, $data{"OUTPUTDIR"}."".$data{"RGLISTANDLENGTHFILENAME"}) || die("Error.  Could not open ".$data{"OUTPUTDIR"}."".$data{"RGLISTANDLENGTHFILENAME"}." for reading.  Terminating!\n"); 
	my @chrLengths = ();
	my @chrFileNames = ();
	while(defined(my $line = <FH>)) {
		chomp($line);
		if($line =~ m/^(.*)\t(.*)$/) {
			push(@chrFileNames, $1);
			push(@chrLengths, $2);
		}
		else {
			die("Error.  Could not understand line in ".$data{"RGLISTANDLENGTHFILENAME"}.":\n$line\nTerminating!\n");
		}
	}
	# close the file
	close(FH);

	# Create rgListFileName
	$data{"RGLISTFILENAME"} = $data{"OUTPUTDIR"}."rgListFileName.txt"; # insert the file name into our data
	open(FH, ">".$data{"RGLISTFILENAME"}) || die("Error.  Could not open ".$data{"RGLISTFILENAME"}." for writing.  Terminating!\n");
	for(my $i=0;$i<scalar(@chrFileNames);$i++) {
		print FH $chrFileNames[$i]."\n";
	}
	# close the file
	close(FH);

	# TODO 
	# Create regions for the indexes and trees
	my @startChr = ();
	my @startPos = ();
	my @endChr = ();
	my @endPos = ();
	my $curPos = $data{"REGIONS"};
	my $cont = 1;
	my $i=0;
	my $length = $chrLengths[0];
	push(@startChr, 1);
	push(@startPos, 1);
	while(1==$cont) {
		while($i<scalar(@chrLengths) && $curPos > $length) {
			$i++;
			if($i<scalar(@chrLengths)) {
				$curPos = $curPos - $length;
				$chrLengths[$i-1];
				$length = $chrLengths[$i];
			}
		}
		if($i < scalar(@chrLengths) && $curPos <= $length) {
			push(@endChr, $i+1);
			push(@endPos, $curPos);
			if($curPos == $length) {
				die("Error.  Could not create regions.  This is a bug.  Please report.  Terminating!\n");
			}
			else {
				push(@startChr, $i+1);
				push(@startPos, $curPos+1);
			}
		}
		else {
			$cont = 0;
		}
		$curPos += $data{"REGIONS"};
	}
	push(@endChr, $i);
	push(@endPos, $length);
	#for(my $i=0;$i<scalar(@startChr);$i++) {
	#	print "chr".$startChr[$i].":".$startPos[$i]." to chr".$endChr[$i].":".$endPos[$i]."\n";
	#}

	# TODO
	# Create a gap file
	my @gaps = split(/,/, $data{"GAPS"});
	$data{"GAPSFILENAME"} = $data{"OUTPUTDIR"}."gaps.txt";
	open(FH, ">".$data{"GAPSFILENAME"}) || die("Error.  Could not open ".$data{"GAPSFILENAME"}." for writing.  Terminating!\n");
	for(my $i=0;$i<scalar(@gaps);$i++) {
		print FH $gaps[$i]."\t";
	}
	print FH "\n";
	# close the file
	close(FH);

	
	local *FHindex;
	local *FHtree;
	open(FHindex, ">".$data{"OUTPUTDIR"}.$data{"GENERATEINDEXESSH"}) || die("Error.  Could not open ".$data{"OUTPUTDIR"}.$data{"GENERATEINDEXESSH"}." for writing.  Terminating!\n");
	open(FHtree, ">".$data{"OUTPUTDIR"}.$data{"GENERATETREESSH"}) || die("Error.  Could not open ".$data{"OUTPUTDIR"}.$data{"GENERATETREESSH"}." for writing.  Terminating!\n");
	my $command = "";
	my @indexFileNames = ();
	my @treeFileNames = ();
	# Generate indexes and trees for each starting/ending chr/pos
	for(my $i=0;$i<scalar(@startChr);$i++) {
		####################################################################
		# BPREPROCESS - Create Indexes.
		####################################################################
		$command = "";
		$command = $data{"COMMANDDIR"}; # Command directory 
		$command .= "bpreprocess/bpreprocess"; # Command 
		$command .= " -i"; # Create Index 
		$command .= " -r ".$data{"RGLISTFILENAME"}; # The reference genome file list 
		$command .= " -l ".$data{"INDEXMATCHLENGTH"}; # The length of reads in the index 
		$command .= " -s ".$startChr[$i];
		$command .= " -S ".$startPos[$i];
		$command .= " -e ".$endChr[$i]; 
		$command .= " -E ".$endPos[$i];
		$command .= " -o ".$data{"OUTPUTID"};
		$command .= " -d ".$data{"OUTPUTDIR"};
		print FHindex $command."\n";
		push(@indexFileNames, sprintf("%sblatter.index.file.%s.%d.%d.%d.%d.%d.%s",
				$data{"OUTPUTDIR"},
				$data{"OUTPUTID"},
				$startChr[$i],
				$startPos[$i],
				$endChr[$i],
				$endPos[$i],
				$data{"INDEXMATCHLENGTH"},
				"bif"));

		####################################################################
		# BPREPROCESS - Create Trees.
		####################################################################
		$command = $data{"COMMANDDIR"}; # Command directory 
		$command .= "bpreprocess/bpreprocess"; # Command 
		$command .= " -r ".$data{"RGLISTFILENAME"}; # The reference genome file list 
		$command .= " -l ".$data{"TREEMATCHLENGTH"}; # The length of reads in the index 
		$command .= " -s ".$startChr[$i];
		$command .= " -S ".$startPos[$i];
		$command .= " -e ".$endChr[$i]; 
		$command .= " -E ".$endPos[$i];
		$command .= " -g ".$data{"GAPSFILENAME"};
		$command .= " -o ".$data{"OUTPUTID"};
		$command .= " -d ".$data{"OUTPUTDIR"};
		print FHtree $command."\n";

		# One file per gap
		for(my $j=0;$j<scalar(@gaps);$j++) {
			push(@treeFileNames, sprintf("%sblatter.tree.file.%s.%d.%d.%d.%d.%d.%d.%s",
					$data{"OUTPUTDIR"},
					$data{"OUTPUTID"},
					$startChr[$i],
					$startPos[$i],
					$endChr[$i],
					$endPos[$i],
					$gaps[$j],
					$data{"TREEMATCHLENGTH"},
					"btf"));

		}
	}
	close(FHindex);
	close(FHtree);

	# Print the indexes to a file 
	$data{"INDEXESFILENAME"} = $data{"OUTPUTDIR"}."indexes.txt";
	open(FH, ">".$data{"INDEXESFILENAME"}) || die("Error.  Could not open ".$data{"INDEXESFILENAME"}." for writing.  Terminating!\n");
	for(my $i=0;$i<scalar(@indexFileNames);$i++) {
		print FH $indexFileNames[$i]."\n";
	}
	print FH "\n";
	# close the file
	close(FH);

	# Print the trees to a file
	$data{"TREESFILENAME"} = $data{"OUTPUTDIR"}."trees.txt";
	open(FH, ">".$data{"TREESFILENAME"}) || die("Error.  Could not open ".$data{"TREESFILENAME"}." for writing.  Terminating!\n");
	for(my $i=0;$i<scalar(@treeFileNames);$i++) {
		print FH $treeFileNames[$i]."\n";
	}
	print FH "\n";
	# close the file
	close(FH);

	# Create offsets file
	$data{"OFFSETSFILENAME"} = $data{"OUTPUTDIR"}."offsets.txt";
	my @offsets = split(/,/, $data{"OFFSETS"});
	open(FH, ">".$data{"OFFSETSFILENAME"}) || die("Error.  Could not open ".$data{"OFFSETSFILENAME"}." for writing.  Terminating!\n");
	for(my $i=0;$i<scalar(@offsets);$i++) {
		print FH $offsets[$i]."\t";
	}
	print FH "\n";
	# close the file
	close(FH);

	# Run once for each set of reads
	local *FHmatch;
	open(FHmatch, ">".$data{"OUTPUTDIR"}.$data{"FINDMATCHESSH"}) || die("Error.  Could not open ".$data{"OUTPUTDIR"}.$data{"FINDMATCHESSH"}." for writing.  Terminating!\n");
	my @matchesFileNames = ();
	my $incReads = $data{"NUMREADS"}/$data{"NUMPROCESSES"};
	for(my $i=$incReads;$i<=$data{"NUMREADS"};$i+=$incReads) {
		####################################################################
		# BMATCHES - find candidate matches
		####################################################################
		$command = "";
		$command = $data{"COMMANDDIR"}; # Command directory 
		$command .= "bmatches/bmatches"; # Command
		$command .= " -I ".$data{"INDEXESFILENAME"};
		$command .= " -T ".$data{"TREESFILENAME"};
		$command .= " -R ".$data{"OUTPUTDIR"}.$data{"READSFILENAME"}; 
		$command .= " -O ".$data{"OFFSETSFILENAME"};
		$command .= " -s ".($i-$incReads+1);
		$command .= " -e ".$i;
		$command .= " -m ".$data{"NUMMISMATCHES"};
		$command .= " -i ".$data{"NUMINSERTIONS"};
		$command .= " -d ".$data{"NUMDELETIONS"};
		if($data{"PAIREDEND"} == 1) {
			$command .= " -2";
		}
		$command .= " -o ".$data{"OUTPUTID"};
		$command .= " -d ".$data{"OUTPUTDIR"};
		print FHmatch $command."\n";
		push(@matchesFileNames, sprintf("%sblatter.matches.file.%s.%d.%d.%d.%d.%d.%d.%s",
				$data{"OUTPUTDIR"},
				$data{"OUTPUTID"},
				($i-$incReads+1),
				$i,
				$data{"NUMMISMATCHES"},
				$data{"NUMINSERTIONS"},
				$data{"NUMDELETIONS"},
				$data{"PAIREDEND"},
				"bmf"));

	}
	close(FHmatch);

	$data{"OFFSETSFILENAME"} = $data{"OUTPUTDIR"}."offsets.txt";
	my @offsets = split(/,/, $data{"OFFSETS"});
	open(FH, ">".$data{"OFFSETSFILENAME"}) || die("Error.  Could not open ".$data{"OFFSETSFILENAME"}." for writing.  Terminating!\n");
	for(my $i=0;$i<scalar(@offsets);$i++) {
		print FH $offsets[$i]."\t";
	}
	print FH "\n";
	# close the file
	close(FH);

	# Run once for each set of reads
	local *FHalign;
	open(FHalign, ">".$data{"OUTPUTDIR"}.$data{"ALIGNSH"}) || die("Error.  Could not open ".$data{"OUTPUTDIR"}.$data{"ALIGNSH"}." for writing.  Terminating!\n");
	my $j=0;
	for(my $i=$incReads;$i<=$data{"NUMREADS"};$i+=$incReads) {
		# Create a matches file name
		my $matchesFileName = $data{"OUTPUTDIR"}."matches.".$data{"OUTPUTID"}.".".($i-$incReads+1).".$i.txt";
		open(FH, ">$matchesFileName") || die("Error.  Could not open $matchesFileName for reading.  Terminating!\n");
		print FH $matchesFileNames[$j]."\n";
		# close the file
		close(FH);

		####################################################################
		# BALIGN - align candidate matches
		####################################################################
		$command = "";
		$command = $data{"COMMANDDIR"}; # Command directory 
		$command .= "balign/balign"; # Command
		$command .= " -r ".$data{"RGLISTFILENAME"}; # The reference genome file list 
		$command .= " -m $matchesFileName";
		$command .= " -s ".$data{"OUTPUTDIR"}.$data{"SCORINGMATRIXFILENAME"};
		$command .= " -a ".$data{"ALGORITHM"};
		$command .= " -O ".$data{"ALIGNOFFSET"};
		if($data{"PAIREDEND"} == 1) {
		$command .= " -2";
		}
		$command .= " -o ".$data{"OUTPUTID"};
		$command .= " -d ".$data{"OUTPUTDIR"};
		print FHalign $command."\n";

	}
	close(FHalign);
}
