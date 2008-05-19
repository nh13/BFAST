#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

# TODO 
############################################################
# bcommand.pl
#
# Usage:
# 	perl bcommand.pl -i <ini file name>
#	perl bcommand.pl -o <output file name>
#
# For each of the options below, specify their values 
# by including a line:
# OPTION="VALUE"
# where the OPTION specifies the option, and in quotes
# is the OPTION's value (VALUE).
#
# The -o option will create a blank bfast.ini file for use
#
# Please see accompanying bfast.ini example file and man
# page.
#
############################################################

# The options allowed by the .ini file. Each 
my @options = (
	# General
	"NUMPROCESSES",
	"NUMTHREADS",
	"PAIREDEND",
	"TIMING",

	# Output
	"COMMANDDIR",
	"OUTPUTDIR",
	"OUTPUTID",

	# Input
	"RGLISTANDLENGTHFILENAME",

	# bpreprocess
	"BPREPROCESS",
	"LAYOUTFILENAME",
	"NUMLAYOUTS",
	"BPREPROCESSINPUTBINARY",
	"REGIONS",
	"STARTCHR",
	"STARTPOS",
	"ENDCHR",
	"ENDPOS",
	"BPREPROCESSOUTPUTBINARY",


	# bmatches
	"BMATCHES",
	"READSFILENAME",
	"OFFSETS",
	"BMATCHESINPUTBINARY",
	"NUMREADS",
	"NUMMISMATCHES",
	"NUMINSERTIONS",
	"NUMDELETIONS",
	"NUMGAPINSERTIONS",
	"NUMGAPDELETIONS",
	"BMATCHESOUTPUTBINARY",

	# balign
	"BALIGN",
	"SCORINGMATRIXFILENAME",
	"ALGORITHM",
	"ALIGNOFFSET",
	"MAXNUMMATCHES");

# Global variables
my $BFAST_INDEX_FILE_EXTENSION = "bif";
my $BFAST_TREE_FILE_EXTENSION = "btf";
my $BFAST_MATCHES_FILE_EXTENSION = "bmf";

# Command line variables
my $iniFileName;
my $outputFileName;

#Command-line options
GetOptions("iniFileName=s"=>\$iniFileName,
	"outputFileName=s"=>\$outputFileName)
	or PrintUsage();
if(!defined($iniFileName) && !defined($outputFileName)) {
	print "Option iniFileName was not given.\n";
	print "Option outputFileName was not given.\n";
	print "One of the two options must be given.\n";
	die("***** Exiting due to errors *****\n");
}
if(defined($iniFileName)) {
	Main($iniFileName);
}
elsif(defined($outputFileName)) {
	PrintIniFile($outputFileName);
}

sub PrintUsage {
	print "\nUsage: perl bcommand.pl [options]\n\n";
	print "-i <ini file name>\n";
	exit(1);
}


############################################################
# Main function
# TODO 
#
# params:
#	data - a reference to an empty hash to store the comands
#		from the .ini file.
#	iniFileName - the name of the .ini file.
#
# return:
# 	null.
############################################################
sub Main {
	my $iniFileName = shift;
	my $indexCommandsFileName;
	my $matchesCommandsFileName;
	my $alignCommandsFileName;

	############################################################
	# Create data object to hold all options
	my %data = ();

	############################################################
	# Create arrays for holding the reference genome 
	# chromosomes and their respective lengths
	my @chrLengths = ();
	my @chrFileNames = ();

	############################################################
	# We will break up the genome into regions and align
	# them seperately.  Create arrays to store the start
	# and end locations.
	my @startChr = ();
	my @startPos = ();
	my @endChr = ();
	my @endPos = ();

	############################################################
	# Create an array to hold the user specified offsets.
	my @offsets = ();

	############################################################
	# Parse the file 
	ParseIniFile($iniFileName, \%data);

	############################################################
	# Read in chromsomes with list and lengths
	ReadChrAndLengths($data{"RGLISTANDLENGTHFILENAME"}, \@chrFileNames, \@chrLengths); 

	############################################################
	# Add a new command to our data 
	$data{"RGLISTFILENAME"} = $data{"OUTPUTDIR"}."rgListFileName.txt"; 

	############################################################
	# Create rgListFileName
	CreateReferenceGenomeFileList($data{"RGLISTFILENAME"}, \@chrFileNames);

	############################################################
	# Create regions for the indexes 
	CreateRegions(\@startChr, 
		\@startPos, 
		\@endChr, 
		\@endPos, 
		\@chrLengths, 
		$data{"REGIONS"},
		$data{"STARTCHR"},
		$data{"STARTPOS"},
		$data{"ENDCHR"},
		$data{"ENDPOS"});

	############################################################
	# Create file names to store the indexes created by
	# the bpreprocess commands.
	$indexCommandsFileName = $data{"OUTPUTDIR"}.$data{"BPREPROCESS"};
	$data{"MAININDEXLISTFILENAME"} = $data{"OUTPUTDIR"}."main.indexes.txt";
	$data{"INDEXLISTFILENAME"} = $data{"OUTPUTDIR"}."indexes.txt";
	# Create and output to file the commans, as well as 
	# outputting the main index and index file names created by the
	# commands.
	CreateBProcessCommandsAndFiles(
		$indexCommandsFileName,
		\@startChr,
		\@startPos,
		\@endChr,
		\@endPos,
		\%data);


	############################################################
	# Create offsets file name
	$data{"OFFSETSFILENAME"} = $data{"OUTPUTDIR"}."offsets.txt";
	# Get the offsets
	@offsets = split(/,/, $data{"OFFSETS"});
	# Create the offsets file
	CreateOffsetsFile($data{"OFFSETSFILENAME"}, \@offsets);

	############################################################
	# BMATCHES
	$matchesCommandsFileName = $data{"OUTPUTDIR"}.$data{"BMATCHES"};
	$data{"MATCHESLISTFILENAME"} = $data{"OUTPUTDIR"}."matches.txt";
	CreateBMatchesCommandsAndFiles(
		$matchesCommandsFileName,
		\%data);

	############################################################
	# BALIGN
	$alignCommandsFileName = $data{"OUTPUTDIR"}.$data{"BALIGN"};
	CreateBAlignCommandsAndFiles(
		$alignCommandsFileName,
		\@startChr,
		\@startPos,
		\@endChr,
		\@endPos,
		\%data);

}

############################################################
# ParseIniFile
#
# Parses the specified .ini file and stores the options in
# the given hash.
#
# params:
#	data - a reference to a hash to store the comands from
#		the .ini file.
#	iniFileName - the name of the .ini file.
#
# return:
#	null.
############################################################
sub ParseIniFile {
	my $iniFileName = shift;
	my $data = shift;
	local *FH;

	# Initialize data 
	%$data = ();

	# Open the file for reading
	open(FH, "$iniFileName") || die("Error.  Could not open $iniFileName for reading.  Terminating!\n"); 

	# Read through every line of the .ini file
	while(defined(my $line = <FH>)) {
		chomp($line);
		if($line =~ m/^#/) {
			# Ignore comment line 
		}
		elsif($line =~ m/^([A-Z]*)="(.*)"/) {
			# Read in the name and a data pair
			my $name = $1;
			my $cmd = $2;
			$$data{$name} = $cmd;
		}
		else {
			# Ignore 
		}
	}

	# close the file
	close(FH);

	# Check that all options are specified
	for(my $i=0;$i<scalar(@options);$i++) {
		if(!defined($$data{$options[$i]})){ 
			die("Error.  ".$options[$i]." was not found in $iniFileName.  Terminating!\n");
		}
	}

	# Check that all paths have an '/' at the end
	if(!($$data{"COMMANDDIR"} =~ m/\/$/)) {
		die("Error.  Within $iniFileName, COMMANDDIR must end with a '/'.  Terminating!\n");
	}
	if(!($$data{"OUTPUTDIR"} =~ m/\/$/)) {
		die("Error.  Within $iniFileName, OUTPUTDIR must end with a '/'.  Terminating!\n");
	}
}

############################################################
# ReadChrAndLengths
#
# Reads the file names of files that contain the reference
# genome for each chromosome respectively.  Additionally, we
# read the length of each chromosome.   The format for a 
# line in the file must be:
# NAME	LENGTH
#
# params:
#	fileName - the name of the input file to be read.
#	chrFileNames - a reference to an array to store the 
#		the chromosome file name.s
#	chrLength - a reference to an array to store the length
#		of each chromosome.
#
# return:
#	null.
############################################################
sub ReadChrAndLengths {
	my $fileName = shift;
	my $chrFileNames = shift;
	my $chrLengths = shift;
	local *FH;

	# Initialize arrays
	@$chrFileNames = ();
	@$chrLengths = ();

	# open the file for reading
	open(FH, "$fileName") || die("Error.  Could not open $fileName for reading.  Terminating!\n");
	while(defined(my $line = <FH>)) {
		chomp($line);
		if($line =~ m/^(.*)\t(.*)$/) {
			push(@$chrFileNames, $1);
			push(@$chrLengths, $2);
		}
		else {
			die("Error.  Could not understand line in $fileName:\n$line\nTerminating!\n");
		}
	}
	# close the file
	close(FH);
}

############################################################
# CreateReferenceGenomeFileList
#
# Writes the chromosome file names to the specified file.
#
# params:
#	fileName - the name of the output file to be written.
#	chrFileNames - a reference to an array to store the 
#		the chromosome file name.s
#
# return:
#	null.
############################################################
sub CreateReferenceGenomeFileList {
	my $fileName = shift;
	my $chrFileNames = shift;
	local *FH;

	# open the output file
	open(FH, ">$fileName") || die("Error.  Could not open $fileName for writing.  Terminating!\n");

	# Write each file name to file
	for(my $i=0;$i<scalar(@$chrFileNames);$i++) {
		print FH $$chrFileNames[$i]."\n";
	}

	# close the file
	close(FH);
}

############################################################
# CreatRegions
#
# Splits the specified region int subregions of a given 
# length.  This is further complicated by the fact that
# chromosomes have different length, and therefore we may
# span one or more chromosomes.
#
# params:
#	startChrArr - a reference to an array in which
#		to store the start chromosome for reach region.
#	startPosArr - a reference to an array in which
#		to store the start positionfor reach region.
#	endChrArr - a reference to an array in which
#		to store the end chromosome for reach region.
#	endPosArr - a reference to an array in which
#		to store the start position for reach region.
#	chrLengths - a reference to an array that has the
#		lengths of every chromosome.
#	regionLength - the length of the subregions.
#	startChr - the chromosome at which to start.
#	startPos - the position at which to start.
#	endChr - the position at which to end.
#	endPos - the position at which to end.
#
# return:
#	null.
############################################################
sub CreateRegions {
	my $startChrArr = shift;
	my $startPosArr = shift;
	my $endChrArr = shift;
	my $endPosArr = shift;
	my $chrLengths = shift;
	my $regionLength = shift;
	my $startChr = shift;
	my $startPos = shift;
	my $endChr = shift;
	my $endPos = shift;

	if($regionLength == 0) {
		$regionLength = 4000000000;
	}

	my $cont = 1;
	my $curChr=0;
	my $curPos=0;
	my $prevChr=0;
	my $prevPos=0;

	# Check start and end bounds 
	if($startChr <= 0) {
		$startChr = 1;
		$startPos = 1;
	}
	elsif($startPos == 0) {
		$startPos = 1;
	}
	if($endChr <= 0) {
		$endChr = scalar(@$chrLengths);
	}
	if($endPos <= 0) {
		$endPos = $$chrLengths[$endChr-1];
	}

	# Initialize arrays
	@$startChrArr = ();
	@$startPosArr = ();
	@$endChrArr = ();
	@$endPosArr = ();

	# Initialize start positions and add them to the array
	$prevPos = $startPos;
	$prevChr = $startChr;
	$curChr = $startChr;
	$curPos = $startPos;
	$startPos = $startPos;

	# While it is reasonable to continue
	while($prevChr < $endChr || ($prevChr == $endChr && $prevPos <= $endPos)) {
		# Update 
		$curPos += $regionLength;
		# Find the next end position.  Move to the next chromosome if we the 
		# region to cover is larger than the remaining positions on the current
		# chromosome.
		while($curChr<$endChr && $curPos > $$chrLengths[$curChr-1]) {
			# Update the current chromosome
			$curChr++;
			# Check that we have chromosomes left
			if($curChr<$endChr) {
				# Update the current position by using the final positions
				# on the previous chromosome and the positions on the new
				# chromosome.
				$curPos = $curPos - $$chrLengths[$curChr-1];
			}
		}
		# Now, if we have found a legal end position, update 
		if($curChr < $endChr || ($curChr == $endChr && $curPos <= $endPos)) {
			if(($curChr != $prevChr) || ($curPos != $prevPos)) {
				# add to arrays 
				push(@$startChrArr, $prevChr);
				push(@$startPosArr, $prevPos);
				push(@$endChrArr, $curChr);
				push(@$endPosArr, $curPos);
			}
		}
		elsif(($endChr != $prevChr) || ($endPos != $prevPos)) {
			push(@$startChrArr, $prevChr);
			push(@$startPosArr, $prevPos);
			push(@$endChrArr, $endChr);
			push(@$endPosArr, $endPos);
		}
		# Update prev
		$prevChr = $curChr;
		$prevPos = $curPos;
	}

	# Print regions
	#for(my $i=0;$i<scalar(@$startChrArr);$i++) {
	#		print sprintf("From chr%d:%d to chr:%d:%d\n",
	#		$$startChrArr[$i],
	#		$$startPosArr[$i],
	#		$$endChrArr[$i],
	#		$$endPosArr[$i]);
	#}
}

############################################################
# CreateBProcessCommandsAndFiles 
# 
# Creates all the bpreprocess commands storing the commands
# in two separate files for the indexes
# respectively.  Also creates two files that list 
# the output files created by the two groups of bpreprocess
# commands.
#
# params:
#	indexCommandsFileName - the name of the file in which 
#		to store all bpreprocess commands to create indexes.
#	startChrArr - a reference to an array in which
#		to store the start chromosome for reach region.
#	startPosArr - a reference to an array in which
#		to store the start positionfor reach region.
#	endChrArr - a reference to an array in which
#		to store the end chromosome for reach region.
#	endPosArr - a reference to an array in which
#		to store the start position for reach region.
#	data - a reference to a hash to store the comands from
#		the .ini file.
#
# return:
#	null.
############################################################
sub CreateBProcessCommandsAndFiles {
	my $indexCommandsFileName = shift;
	my $startChrArr = shift;
	my $startPosArr = shift;
	my $endChrArr = shift;
	my $endPosArr = shift;
	my $data = shift;

	local *FHMainIndexes;
	local *FHIndexes; 
	local *FHIndexCommands;
	my $command = "";
	my $mainIndexListFileName = $$data{"MAININDEXLISTFILENAME"};
	my $indexListFileName = $$data{"INDEXLISTFILENAME"};

	# Open output files
	open(FHMainIndexes, ">$mainIndexListFileName") || die("Error.  Could not open $mainIndexListFileName for writing.  Terminating!\n");
	open(FHIndexes, ">$indexListFileName") || die("Error.  Could not open $indexListFileName for writing.  Terminating!\n");
	open(FHIndexCommands, ">$indexCommandsFileName") || die("Error.  Could not open $indexCommandsFileName for writing.  Terminating!\n");

	# Generate indexes for each starting/ending chr/pos
	for(my $i=0;$i<scalar(@$startChrArr);$i++) {

		############################################################
		# BPREPROCESS - Create Indexes.
		$command = $$data{"COMMANDDIR"}; # Command directory 
		$command .= "bpreprocess/bpreprocess"; # Command 
		$command .= " -r ".$$data{"RGLISTFILENAME"}; # The reference genome file list 
		$command .= " -i ".$$data{"OUTPUTDIR"}.$$data{"LAYOUTFILENAME"}; # The layout file name
		$command .= " -a 0"; # Create Index 
		if($$data{"BPREPROCESSINPUTBINARY"}==1) {
			$command .= " -b";
		}
		$command .= " -s ".$$startChrArr[$i];
		$command .= " -S ".$$startPosArr[$i];
		$command .= " -e ".$$endChrArr[$i]; 
		$command .= " -E ".$$endPosArr[$i];
		$command .= " -n ".$$data{"NUMTHREADS"};
		$command .= " -o ".$$data{"OUTPUTID"};
		$command .= " -d ".$$data{"OUTPUTDIR"};
		if($$data{"BPREPROCESSOUTPUTBINARY"}==1) {
			$command .= " -B";
		}
		if($$data{"TIMING"}==1) {
			$command .= " -t";
		}
		# Print the command
		print FHIndexCommands $command."\n";

		# Print the index file names that this command will create. 
		for(my $j=0;$j<$$data{"NUMLAYOUTS"};$j++) {

			if($j==0) {
				print FHMainIndexes sprintf("%sbfast.index.file.%s.%d.%d.%d.%d.%d.%s\n",
					$$data{"OUTPUTDIR"},
					$$data{"OUTPUTID"},
					$$startChrArr[$i],
					$$startPosArr[$i],
					$$endChrArr[$i],
					$$endPosArr[$i],
					$j+1,
					$BFAST_INDEX_FILE_EXTENSION);
			}
			else {
				print FHIndexes sprintf("%sbfast.index.file.%s.%d.%d.%d.%d.%d.%s\n",
					$$data{"OUTPUTDIR"},
					$$data{"OUTPUTID"},
					$$startChrArr[$i],
					$$startPosArr[$i],
					$$endChrArr[$i],
					$$endPosArr[$i],
					$j+1,
					$BFAST_INDEX_FILE_EXTENSION);
			}
		}
	}
	# Close the output files 
	close(FHMainIndexes);
	close(FHIndexes);
	close(FHIndexCommands);
}

############################################################
# CreateOffsetsFile
#
# Writes the array of offsets to the offsets file.
#
# params:
#	fileName - the name of the output file to be written.
#	offsets - a reference to an array of integers representing
#		the offsets.
#
# return:
#	null.
############################################################
sub CreateOffsetsFile {
	my $fileName = shift;
	my $offsets = shift;
	local *FH;

# Open the output file
	open(FH, ">$fileName") || die("Error.  Could not open $fileName for writing.  Terminating!\n");

# Print each offset
	for(my $i=0;$i<scalar(@$offsets);$i++) {
		print FH $$offsets[$i]."\t";
	}
	print FH "\n";

# close the file
	close(FH);
}

############################################################
# CreateBMatchesCommandsAndFiles 
# 
# Creates all the bmatches commands storing the commands
# in a given file. Also creates a list of the output files
# generated by the commands.
#
# params:
#	matchesCommandsFileName - the name of the file in which 
#		to store all bmatches commands.
#	matchesListFileName - the name of the file in which to 
#		store all matches files generated.
#	data - a reference to a hash to store the comands from
#		the .ini file.
#
# return:
#	null.
############################################################
sub CreateBMatchesCommandsAndFiles {
	my $matchesCommandsFileName = shift;
	my $data = shift;

	local *FHMatchesCommands;
	local *FHMatchesList;

	my $incReads;
	my $matchesListFileName = $$data{"MATCHESLISTFILENAME"};
	my $numReads = $$data{"NUMREADS"};
	my $numProcesses = $$data{"NUMPROCESSES"};
	my $command;

	# Open the file in which to store the commands
	open(FHMatchesCommands, ">$matchesCommandsFileName") || die("Error.  Could not open $matchesCommandsFileName for writing.  Terminating!\n");
	open(FHMatchesList, ">$matchesListFileName") || die("Error.  Could not open $matchesListFileName for writing.  Terminating!\n");

	# Initialize the number of reads per bmatches command to process.
	$incReads = $numReads/$numProcesses; 

	# Create each bmatches command
	for(my $i=$incReads;$i<=$numReads;$i+=$incReads) {
		############################################################
		# BMATCHES - find candidate matches
		############################################################
		$command = $$data{"COMMANDDIR"}; # Command directory 
		$command .= "bmatches/bmatches"; # Command
		$command .= " -r ".$$data{"RGLISTFILENAME"}; # The reference genome file list 
		$command .= " -i ".$$data{"MAININDEXLISTFILENAME"};
		$command .= " -I ".$$data{"INDEXLISTFILENAME"};
		$command .= " -R ".$$data{"OUTPUTDIR"}.$$data{"READSFILENAME"}; 
		$command .= " -O ".$$data{"OFFSETSFILENAME"};
		if($$data{"BMATCHESINPUTBINARY"} == 1) {
			$command .= " -b";
		}
		$command .= " -s ".($i-$incReads+1);
		$command .= " -e ".$i;
		$command .= " -x ".$$data{"NUMMISMATCHES"};
		$command .= " -y ".$$data{"NUMINSERTIONS"};
		$command .= " -z ".$$data{"NUMDELETIONS"};
		$command .= " -Y ".$$data{"NUMGAPINSERTIONS"};
		$command .= " -Z ".$$data{"NUMGAPDELETIONS"};
		if($$data{"PAIREDEND"} == 1) {
			$command .= " -2";
		}
		$command .= " -n ".$$data{"NUMTHREADS"};
		$command .= " -o ".$$data{"OUTPUTID"};
		$command .= " -d ".$$data{"OUTPUTDIR"};
		if($$data{"BMATCHESOUTPUTBINARY"}==1) {
			$command .= " -B";
		}
		if($$data{"TIMING"}==1) {
			$command .= " -t";
		}
		# Print the command
		print FHMatchesCommands $command."\n";
		# Print the file generated
		print FHMatchesList sprintf("%sbfast.matches.file.%s.%d.%d.%d.%d.%d.%d.%d.%d.%s\n",
			$$data{"OUTPUTDIR"},
			$$data{"OUTPUTID"},
			($i-$incReads+1),
			$i,
			$$data{"NUMMISMATCHES"},
			$$data{"NUMINSERTIONS"},
			$$data{"NUMDELETIONS"},
			$$data{"NUMGAPINSERTIONS"},
			$$data{"NUMGAPDELETIONS"},
			$$data{"PAIREDEND"},
			"bmf");

	}
	close(FHMatchesCommands);
	close(FHMatchesList);
}

############################################################
# CreateBAlignCommandsAndFiles 
# 
# Creates all the balign commands storing the commands
# in a given file. 
#
# params:
#	alignCommandsFileName - the name of the file in which 
#		to store all balign commands.
#	startChrArr - a reference to an array in which
#		to store the start chromosome for reach region.
#	startPosArr - a reference to an array in which
#		to store the start positionfor reach region.
#	endChrArr - a reference to an array in which
#		to store the end chromosome for reach region.
#	endPosArr - a reference to an array in which
#		to store the start position for reach region.
#	data - a reference to a hash to store the comands from
#		the .ini file.
#
# return:
#	null.
############################################################
sub CreateBAlignCommandsAndFiles {
	my $alignCommandsFileName = shift;
	my $startChrArr = shift;
	my $startPosArr = shift;
	my $endChrArr = shift;
	my $endPosArr = shift;
	my $data = shift;

	my $numReads = $$data{"NUMREADS"};
	my $numProcesses = $$data{"NUMPROCESSES"};
	my $incReads = $numReads/$numProcesses;
	my $j=0;
	my $command;
	my $alignFileName;

	local *FHAlignCommands;

# Open the file in which to store the commands
	open(FHAlignCommands, ">$alignCommandsFileName") || die("Error.  Could not open $alignCommandsFileName for writing.  Terminating!\n");

# Create balign command

############################################################
# BALIGN - align candidate align
####################################################################
	$command = $$data{"COMMANDDIR"}; # Command directory 
	$command .= "balign/balign"; # Command
	$command .= " -r ".$$data{"RGLISTFILENAME"}; # The reference genome file list 
	$command .= " -m ".$$data{"MATCHESLISTFILENAME"};
	$command .= " -x ".$$data{"OUTPUTDIR"}.$$data{"SCORINGMATRIXFILENAME"};
	$command .= " -a ".$$data{"ALGORITHM"};
	$command .= " -s ".$$startChrArr[0];
	$command .= " -S ".$$startPosArr[0];
	$command .= " -e ".$$endChrArr[scalar(@$startPosArr)-1]; 
	$command .= " -E ".$$endPosArr[scalar(@$startPosArr)-1];
	$command .= " -O ".$$data{"ALIGNOFFSET"};
	$command .= " -M ".$$data{"MAXNUMMATCHES"};
	if($$data{"PAIREDEND"} == 1) {
		$command .= " -2";
	}
	$command .= " -n ".$$data{"NUMTHREADS"};
	$command .= " -o ".$$data{"OUTPUTID"};
	$command .= " -d ".$$data{"OUTPUTDIR"};
	if($$data{"TIMING"}==1) {
		$command .= " -t";
	}
	print FHAlignCommands $command."\n";

# Close the output file
	close(FHAlignCommands);
}

############################################################
# PrintIniFile
#
# Prints a sample .ini file to the output file
#
# params:
#	outputFileName - the name of the output file
#
# return:
#	null.
############################################################
sub PrintIniFile {
	my $outputFileName = shift;
	local *FH;

	# Open the output file
	open(FH, ">$outputFileName") || die("Error.  Could not open $outputFileName for writing.  Terminating!\n");

	# Print the options
	for(my $i=0;$i<scalar(@options);$i++) {
		print FH $options[$i]."=\"\"\n";
	}

	# Close the output file
	close(FH);
}
