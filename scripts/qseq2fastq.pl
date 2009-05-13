#!/usr/bin/perl -w

# Author: Nils Homer

use strict;
use warnings;
use Getopt::Std;

select STDERR; $| = 1;  # make unbuffered

my %opts;
my $version = '0.1.1';
my $usage = qq{
Usage: qseq2fastq.pl <input .qseq prefix>

	This script will convert Illumina output files (*qseq files)
	to the BFAST fastq multi-end format.  For single-end reads 
	that do not have more than one end (neither mate-end nor
	paired end), the output format is strictly fastq format.
	For multi-end data (mate-end, paired end ect.) each end is
	outputted consecutively in 5'->3' ordering in fastq blocks.
	Therefore, the output reamins strictly conforming to the 
	fastq format, but keeps intact the relationships between
	sequences that all originate from the same peice of DNA.
	We assume for paired end data that all data is able to be
	paired.

	All .qseq files will be inferred from the specified directory 
	and accompanying prefix.  For paired end data, reads in
		<prefix>_1_XXXX_qseq.txt and <prefix>_2_XXXX_qseq.txt will be
	paired.  For example, if we wish to create a paired end fastq
	file the first lane, the command should be:

		qseq2fastq.pl s_1_
};
my $ROTATE_NUM = 100000;

getopts('', \%opts);
die($usage) if (@ARGV < 1);

my $input_prefix = shift @ARGV;

# Get input files
my @files_one = ();
my @files_two = ();
GetDirContents($input_prefix."_1", \@files_one);
GetDirContents($input_prefix."_2", \@files_two);

if(0 == scalar(@files_one)) {
	die("Error.  Could not find any files.  Terminating!\n");
}
elsif(0 < scalar(@files_two) && scalar(@files_one) != scalar(@files_two)) {
	die("Error.  Did not find an equal number of paired end files.  Terminating!\n");
}

my $FH_index = 0;

if(0 == scalar(@files_two)) { # Single end
	while($FH_index < scalar(@files_one)) {
		open(FH_one, "$files_one[$FH_index]") || die;
		while(defined(my $line = <FH_one>)) {
			my ($name, $seq, $qual) = parse_line($line, 1);
			print STDOUT "$name\n$seq\n+\n$qual\n";
		}
		close(FH_one);
		$FH_index++;
	}
}
else { # Paired end
	while($FH_index < scalar(@files_one)) {
		my $min_read_name = "";
		print STDOUT "Opening $files_one[$FH_index]\n";
		open(FH_one, "$files_one[$FH_index]") || die;
		open(FH_two, "$files_two[$FH_index]") || die;
		my %read_one = ();
		my %read_two = ();
		while(1 == get_read(*FH_one, \%read_one, 1) &&
			1 == get_read(*FH_two, \%read_two, 2)) {
			if(!($read_one{"NAME"} eq $read_two{"NAME"})) {
				die;
			}
			print STDOUT "".$read_one{"NAME"}."\n".$read_one{"SEQ"}."\n+\n".$read_one{"QUAL"}."\n";
			print STDOUT "".$read_two{"NAME"}."\n".$read_two{"SEQ"}."\n+\n".$read_two{"QUAL"}."\n";
		}
		close(FH_one);
		close(FH_two);
		$FH_index++;
	}
}

sub cmp_read_names {
	my $a = shift;
	my $b = shift;

	# characters, the numbers, then characters, then numbers, ...
	# recursion is for sissies
	while(0 != length($a) &&
		0 != length($b)) {
		my $a_char = "";
		my $b_char = "";
		my $a_num = 0;
		my $b_num = 0;
		if($a =~ m/^(\d+)(.*?)$/) {
			$a_num = $1;
			$a = $2;
		}
		elsif($a =~ m/^(\D+)(.*?)$/) {
			$a_char = $1;
			$a = $2;
		}
		if($b =~ m/^(\d+)(.*?)$/) {
			$b_num = $1;
			$b = $2;
		}
		elsif($b =~ m/^(\D+)(.*?)$/) {
			$b_char = $1;
			$b = $2;
		}
		# Compare numbers then letters
		if(!($a_char eq $b_char)) {
			return ($a_char cmp $b_char);
		}
		elsif($a_num != $b_num) {
			return ($a_num <=> $b_num);
		}
	}

	return (length($a) <=> length($b));
}

sub GetDirContents {
	my $prefix = shift;
	my $dirs = shift;

	my $dir = "";
	if($prefix =~ m/^(.+\/)([^\/]*)/) {
		$dir = $1;
		$prefix = $2;
	}
	else {
		$dir = "./";
	}

	print STDOUT "dir=$dir\n";
	local *DIR;
	opendir(DIR, "$dir") || die("Error.  Could not open $dir.  Terminating!\n");
	@$dirs = grep !/^\.\.?$/, readdir DIR;
	@$dirs = grep /^$prefix.*qseq\.txt/, @$dirs;
	for(my $i=0;$i<scalar(@$dirs);$i++) {
		@$dirs[$i] = $dir."".@$dirs[$i];
	}
	close(DIR);
}

sub parse_line {
	my $line = shift;
	my $end = shift;

	my @arr = split(/\s+/, $line);

	# Convert Illumina PHRED to Sanger PHRED
	my $qual = "";
	for(my $i=0;$i<length($arr[9]);$i++) {
		my $Q = ord(substr($arr[9], $i, 1)) - 64;
		$qual .= chr(($Q<=93? $Q : 93) + 33);
	}

	my $name = "@".$arr[0]."_".$arr[1]."_".$arr[2]."_".$arr[3]."_".$arr[4]."_".$arr[5]."_".$arr[6]."";
	my $seq = $arr[8];

	if(1 != $end) {
		$seq = reverse($seq);
		$qual = reverse($qual);
	}

	return ($name, $seq, $qual);
}
		
sub get_read {
	my $FH = shift;
	my $read = shift;
	my $end = shift;

	if(defined(my $line = <$FH>)) {
		my ($name, $seq, $qual) = parse_line($line, $end);

		$read->{"NAME"} = $name;
		$read->{"SEQ"} = $seq;
		$read->{"QUAL"} = $qual;

		return 1;
	}
	else {
		return 0;
	}
}
