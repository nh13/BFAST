#!/usr/bin/perl
#
# Get file offsets when splitting csfasta files into roughly even sizes
#
# Kevin Squire (23 Sept 2010), ksquire@mednet.ucla.edu

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use File::Path qw(mkpath);
use POSIX qw/ceil/;
use Fcntl ":seek";
use Cwd 'abs_path';
use Time::HiRes qw( gettimeofday tv_interval );

my $path = dirname(abs_path($0));
#my $extract_region = "$path/extract_region";
my $solid2fastq = "$path/solid2fastq";

my $help;
my $man;

my $dryrun = 0;
my $split_size;
my $chunks;
my $bwa_output;
my $single_output;
my $gzip_output;
my $output_prefix;
my $readname_prefix;
my $basename;
my $wait = 0;

GetOptions('help' => \$help,
           man => \$man,
           'dryrun' => \$dryrun,
           'split-size=i' => \$split_size,
           'chunks=i' => \$chunks,
	   'output-prefix=s' => \$output_prefix,
	   'readname-prefix=s' => \$readname_prefix,
           'wait' => \$wait,
	   'b' => \$bwa_output,
	   'w' => \$single_output,
	   'Z' => \$gzip_output) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose =>2) if $man;
pod2usage("Please specify exactly one of --split_size or --chunks.") if ((!defined($split_size) && !defined($chunks)) || ($split_size && $chunks));
pod2usage("Wrong number of arguments.") if ($#ARGV != 1 && $#ARGV != 3);
pod2usage("Please specify at most one of -b or -w.") if ($bwa_output && $single_output);

my $FAKEQSUBID = 0;

################################################

sub SubmitJob {

   my ($run_file, $dryrun, $command, $outputID, $threads, $dependent_job_ids, $sync) = @_;

   print STDERR "[psolid2fastq] RUNFILE=$run_file\n";

   my $output = <<END_OUTPUT;
run ()
{
        echo "running: \$*" 2>&1;
        eval \$*;
        if test \$? != 0 ; then
        echo "error: while running '\$*'";
        exit 100;
        fi
}
END_OUTPUT

   $output .= "\nrun \"hostname\";\n";
   $output .= "run \"$command\";\n";
   $output .= "exit 0;\n";

   open(FH, ">$run_file") or die("Error.  Could not open $run_file for writing!\n");
   print FH "$output";
   close(FH);

   # Create qsub command
   my $qsub = "qsub -terse";

   if (1 < $threads) {
      $qsub .= " -pe serial ".$threads;
   }

   if(0 < scalar(@$dependent_job_ids)) {
       $qsub .= " -hold_jid ".join(",", @$dependent_job_ids);
   }

   if (0 < $sync) {
      $qsub .= " -sync y";
   }

   $qsub .= " -N $outputID -o $run_file.out -e $run_file.err $run_file";

   if(1 == $dryrun) {
      $FAKEQSUBID++;
      print $qsub."\n";
      print STDERR "[psolid2fastq] NAME=$outputID QSUBID=$FAKEQSUBID\n";
      return $FAKEQSUBID;
   }

   my $qsub_id = `$qsub`;
   $qsub_id = "$qsub_id";
   chomp($qsub_id);

   die("Error submitting QSUB_COMMAND=$qsub\nQSUB_ID=$qsub_id\n") unless (0 < length($qsub_id));
   if($qsub_id !~ m/^\S+$/) {
      die("Error submitting QSUB_COMMAND=$qsub\nQSUB_ID=$qsub_id\n") unless (0 < length($qsub_id));
   }

   print STDERR "[psolid2fastq] NAME=$outputID QSUBID=$qsub_id\n";

   return $qsub_id;
}


sub estimate_avg_read_size {

   my $FILE = shift;
   my $total_size = 0;

   # skip header
   while (<$FILE>) {
      last if (!m/^#/);
   }

   # First read line (read above)
   $total_size = length($_);

   for (my $i = 1; $i < 64; $i++) {
      $total_size += length(<$FILE>);
   }

   seek ($FILE, -2*$total_size, SEEK_END);

   my $read_name;
   get_next_read($FILE, \$read_name);

   for (my $i = 0; $i < 64; $i++) {
      $total_size += length(<$FILE>);
   }

   my $avg_size = $total_size / 64;

   seek ($FILE, 0, SEEK_SET);

   return $avg_size;

}

sub decode {

   my $read_name = shift;

   my ($f1, $f2, $f3, $junk) = split (/_/, substr($read_name, 1));

   return int(sprintf("%04d%04d%04d", $f1, $f2, $f3));
}

sub get_next_read {

   my ($FILE, $readname) = @_;
   my @lines = ();
   my $i;
   my $line;
   my $pos;

   for ($i = 0; $i < 3; $i++) {
      $pos = tell($FILE);
      $line = <$FILE>;
      if ($line =~ m/^>/) {
         $$readname = decode($line);
         return $pos;
      }
   }

   die "Can't find next line!";
}

sub get_matching_read {
   my ($FILE, $readname, $qual_seek_pos, $step_size) = @_;
   my $pos;
   my $qual_readname;

   my $last = $#$qual_seek_pos;

   my $step_back = - $step_size/250;

   if (@$qual_seek_pos < 2) {
       $pos = $qual_seek_pos->[$last] + $step_size;
   }
   elsif (@$qual_seek_pos == 2) {
       $pos = $qual_seek_pos->[1] * 2 - $qual_seek_pos->[0];
   }
   else {
       my ($a, $b, $c);
       $a = $qual_seek_pos->[$last-1];
       $c = ($qual_seek_pos->[$last-2] + $qual_seek_pos->[$last])/2 - $a;
       $b = $a + $c - $qual_seek_pos->[$last-2];

       $pos = $a + 2*$b + 4*$c;
   }

   seek($FILE, $pos, 0);

   #print STDERR "Readname: $readname\r";

   do {
      # find the next readname
      while(<$FILE>) {
	 if (m/^>/) {
	    last;
	 }
      }

      $qual_readname = decode($_);

      if($qual_readname > $readname) {
	 seek($FILE, $step_back, 1);
      }

   } until ($qual_readname <= $readname);

   #print STDERR "Readname: $readname   Qual readname: $qual_readname\r";

   do {
      if (m/^>/) {
	 $qual_readname = decode($_);
	 #print STDERR "Readname: $readname   Qual readname: $qual_readname\r";
      }

      if (m/^>/ && decode($_) == $readname) {
         $pos = tell($FILE) - length;
	 #print STDERR "Qual file pos: $pos                                                           \n";

         return $pos;
      }
   } while (<$FILE>);

   return -1;

}

################################################

my $paired = 0;
$paired = 1 if ($#ARGV == 3);

my @csfasta_files = ();
my @csfasta_file_sizes = ();
my @csfasta_chunk_sizes = ();
my @csfasta_fhs = ();
my @csfasta_pos = (0,0);
my @csfasta_seek_pos = ([0],[0]);

my @qual_files = ();
my @qual_file_sizes = ();
my @qual_chunk_sizes = ();
my @qual_fhs = ();
my @qual_pos = (0,0);
my @qual_seek_pos = ([0],[0]);
my $i;

for ($i = 0; $i <= $paired; $i++) {
   $csfasta_files[$i] = rel2abs(shift @ARGV);
}

for ($i = 0; $i <= $paired; $i++) {
   $qual_files[$i] = rel2abs(shift @ARGV);
}

# Get file sizes, and open the files
my $file;
my $filesize;

for $file (@csfasta_files) {
   my $filesize = -s $file;
   push(@csfasta_file_sizes, $filesize);
   open(my $fh, "<", $file) || die "Could not open $file!";
   push(@csfasta_fhs, $fh);

   if (!defined($chunks)) {
      my $size = estimate_avg_read_size($fh);
      $chunks = ceil($filesize / $size / $split_size);
      print STDERR "Splitting file into $chunks chunks.\n";
   }

   push(@csfasta_chunk_sizes, int($filesize / $chunks));
}

for $file (@qual_files) {
   $filesize = -s $file;
   push(@qual_file_sizes, $filesize);
   push(@qual_chunk_sizes, int($filesize / $chunks));

   open(my $fh, "<", $file) || die "Could not open $file!";
   push(@qual_fhs, $fh);
}

my $input_dir = rel2abs(dirname($csfasta_files[0]));

if (!defined($output_prefix)) {
   $basename = fileparse($csfasta_files[0], "_F[35].*\.csfasta");
   $output_prefix = $input_dir."/".$basename;
}
else {
   $output_prefix = rel2abs($output_prefix);
   $basename = fileparse($output_prefix);
   my $output_dir = dirname($output_prefix);
   if (!-d$output_dir) {
      mkpath($output_dir);
   }
}

# set run dir
my $RUNDIR = "$input_dir/run.files";
if (!-d$RUNDIR) {
   mkpath($RUNDIR);
}

# set up start and end strings
my @start_strings = ();
my @end_strings = ();

my $start_string = "0";
for ($i = 0; $i <= $paired*2; $i++) {
   $start_string = "$start_string,0";
}

push (@start_strings, $start_string);

my $chunk;

# For each chunk > 1
for ($chunk = 1; $chunk < $chunks; $chunk++) {

   print STDERR "Finding boundaries of chunk $chunk...\r";

   # get next read_name_1 at pos_csfasta_1

   my $read_name_1;
   my $read_name_2;

   seek ($csfasta_fhs[0], $chunk * $csfasta_chunk_sizes[0], 0);   
   $csfasta_pos[0] = get_next_read($csfasta_fhs[0], \$read_name_1);

   if ($paired) {
      seek ($csfasta_fhs[1], $chunk * $csfasta_chunk_sizes[1], 0);
      $csfasta_pos[1] = get_next_read($csfasta_fhs[1], \$read_name_2);

      do {
         # find read_name_2 >= read_name_1 in csfasta_2
	 #print "Read name 1: $read_name_1\n";
         while($read_name_2 < $read_name_1) {
            $csfasta_pos[1] = get_next_read($csfasta_fhs[1], \$read_name_2);
	    #print "Read name 2: $read_name_2\r";
         }
         # find read_name_1 >= read_name_2 in csfasta_1
	 #print "Read name 2: $read_name_2\n";
         while($read_name_1 < $read_name_2) {
            $csfasta_pos[0] = get_next_read($csfasta_fhs[0], \$read_name_1);
	    #print "Read name 1: $read_name_1\r";
         }
      } while ($read_name_1 != $read_name_2);
      # while read_name_1 != read_name_2
   }

   push(@{$csfasta_seek_pos[0]}, $csfasta_pos[0]);
   if ($paired) {
       push(@{$csfasta_seek_pos[1]}, $csfasta_pos[1]);
   }

   # find read_name_1 in qual_1

   $qual_pos[0] = $chunk * $qual_chunk_sizes[0];   
   $qual_pos[0] = get_matching_read($qual_fhs[0], $read_name_1, $qual_seek_pos[0], $qual_chunk_sizes[0]);
   push(@{$qual_seek_pos[0]}, $qual_pos[0]);

   if ($paired) {
      # find read_name_2 in qual_2
      $qual_pos[1] = $chunk * $qual_chunk_sizes[1];
      $qual_pos[1] = get_matching_read($qual_fhs[1], $read_name_2, $qual_seek_pos[1], $qual_chunk_sizes[1]);
      push(@{$qual_seek_pos[1]}, $qual_pos[1]);
   }

   # save offsets
   if (!$paired) {
      $start_string = $csfasta_pos[0].",".$qual_pos[0];
   }
   else {
      $start_string = $csfasta_pos[0].",".$csfasta_pos[1].",".$qual_pos[0].",".$qual_pos[1];
   }

   push (@start_strings, $start_string);

}

# set up end strings

@end_strings = @start_strings;
shift(@end_strings);

my $end_string;

if (!$paired) {
   $end_string = $csfasta_file_sizes[0].",".$qual_file_sizes[0];
}
else {
   $end_string =$csfasta_file_sizes[0].",".$csfasta_file_sizes[1].",".$qual_file_sizes[0].",".$qual_file_sizes[1];
}

push (@end_strings, $end_string);

################
# create jobs...
################

my $command;
my $run_file;
my $id;
my $solid2fastq_args;
my $compress_arg;

my $job_id;
my @job_ids;

$chunk = 0;

for my $start (@start_strings) {
   my $end = shift(@end_strings);

   $chunk++;

   $solid2fastq_args = "-N $chunk";

   if ($bwa_output) {
      $solid2fastq_args .= " -b";
   }
   elsif ($single_output) {
      $solid2fastq_args .= " -w";
   }

   if ($gzip_output) {
      $solid2fastq_args .= " -Z";
   }

   if (defined($readname_prefix)) {
      $solid2fastq_args .= " -p $readname_prefix";
   }

   $command = "$solid2fastq -o $output_prefix $solid2fastq_args -s $start -e $end ";
   if (!$paired) {
      $command = $command.$csfasta_files[0]." ".$qual_files[0];
   }
   else {
      $command = $command.$csfasta_files[0]." ".$csfasta_files[1]." ".$qual_files[0]." ".$qual_files[1];
   }
   $id = "solid2fastq.$basename.$chunk";
   $run_file = "$RUNDIR/$id.sh";
   my @a_sub = ();
   $job_id = SubmitJob($run_file, $dryrun, $command, $id, 1, \@a_sub, 0);

   push (@job_ids, $job_id);
}

if (0 < $wait) {
   $command = "echo 'Splitting Complete.'";

   $id = "psolid2fastq.$basename";
   $run_file = "$RUNDIR/$id.sh";

   $job_id = SubmitJob($run_file, $dryrun, $command, $id, 1, \@job_ids, 1);
}

__END__

=head1 NAME

psolid2fastq.pl - determine file offsets of split starts when splitting SOLiD csfasta, qual files

=head1 SYNOPSIS

psolid2fastq.pl [options] [csfasta_files] [qual_files]

 Options:
   --help                              brief help message
   --man                               full documentation
   --dryrun                            do a dry run only
   --split-size SPLIT_SIZE             approximate size of splits
   --chunks CHUNKS                     number of chunks to split the file into
   --output-prefix PREFIX              output prefix
   --readname-prefix READNAME_PREFIX   read name prefix
   --wait                              wait for last job to finish before exiting
   -b                                  Enable bwa output (for 'bwa aln', not for 'bfast bwaaln').
   -w                                  Create a single file to dump reads with only one end.
   -Z                                  Output files are gzip comressed

Only specify one of --split-size or --chunks

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=item B<--dryrun>

Does a dry run.  Run scripts are created, but not executed.

=item B<--split-size SPLIT_SIZE>

Approximate size of splits.  For speed purposes, much of the file is not
read, so each split will only be approximately this size.  If specified,
the number of chunks is automatically determined.

=item B<--chunks CHUNKS>

Number of chunks to split the file into.

=item B<--output-prefix>

Output prefix

=item B<--readname-prefix>

Read name prefix.  This will be prepended to all read names in the file.

=item B<--wait>

Wait for last job to finish running before exiting

=item B<-b>

Enable bwa output (for 'bwa aln', not for 'bfast bwaaln') (solid2fastq option)

=item B<-w>

Create a single file to dump reads with only one end (solid2fastq option)

=item B<-Z>

Output files are gzip compressed (solid2fastq option)

=back

=head1 DESCRIPTION

This script takes in approximate split sizes or a chunk count, and determines
seek positions in csfasta and qual files for future splitting.

=cut

