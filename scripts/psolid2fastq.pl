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
use List::Util qw( min );

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
my $merge_single = 0;

my $output_dir;

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
	   'Z' => \$gzip_output,
	   'merge-single' => \$merge_single) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose =>2) if $man;
pod2usage("Please specify exactly one of --split_size or --chunks.") if ((!defined($split_size) && !defined($chunks)) || ($split_size && $chunks));
pod2usage("Wrong number of arguments.") if ($#ARGV != 1 && $#ARGV != 3);
pod2usage("Please specify at most one of -b or -w.") if ($bwa_output && $single_output);
pod2usage("--merge-single is only useful for -b or -w modes.") if ($merge_single && !$bwa_output && !$single_output);

my $FAKEQSUBID = 0;

################################################

sub SubmitArrayJob {
   my ($run_file, $job_file, $dryrun, $commands, $outputID, $threads, $dependent_job_ids, $sync) = @_;

   print STDERR "[psolid2fastq] RUNFILE=$run_file  JOBFILE=$job_file\n";

   open(FH, ">$job_file") or die("Error.  Could not open $job_file for writing!\n");

   # create jobs file

   my $job_num = 0;
   for my $command (@$commands) {
      $job_num++;
      print FH "############# Job $job_num start #############\n";
      print FH "$command\n";
      print FH "############# Job $job_num end #############\n";
   }
   close(FH);

   # Create script file

   my $output = <<'END_OUTPUT';
run ()
{
    echo \"running: $*\" 2>&1;
    eval $*;
    if test $? != 0 ; then
        echo "error: while running '$*'";
        exit 100;
    fi
}

echo "==========="
echo "Job $SGE_TASK_ID"
echo "==========="
run "date";
run "hostname";
lines=`awk -v job_start="Job $SGE_TASK_ID start" -v job_end="Job $SGE_TASK_ID end" \
'$0 ~ job_start {
    while (getline > 0 && $0 !~ job_end) print;
    exit;
}' $job_file`

while read -r line
do
    run "$line";
done <<< "$lines"

exit 0;
END_OUTPUT

   $output =~ s/\$job_file/$job_file/g;

   open(FH, ">$run_file") or die("Error.  Could not open $run_file for writing!\n");
   print FH "$output";
   close(FH);

   # Create qsub command
   my $qsub = "qsub -terse";

   $qsub .= " -t 1-$job_num";

   if (1 < $threads) {
      $qsub .= " -pe serial ".$threads;
   }

   if(0 < scalar(@$dependent_job_ids)) {
       $qsub .= " -hold_jid ".join(",", @$dependent_job_ids);
   }

   if (0 < $sync) {
      $qsub .= " -sync y";
   }

   $qsub .= " -N $outputID -o $run_file.\$TASK_ID.out -e $run_file.\$TASK_ID.err $run_file";

   if(1 == $dryrun) {
      $FAKEQSUBID++;
      print $qsub."\n";
      print STDERR "[psolid2fastq] NAME=$outputID QSUBID=$FAKEQSUBID\n";
      return $FAKEQSUBID;
   }

   my $qsub_id = `$qsub`;
   $qsub_id = "$qsub_id";
   chomp($qsub_id);
   $qsub_id =~ s/(\d+)\..*/$1/;

   die("Error submitting QSUB_COMMAND=$qsub\nQSUB_ID=$qsub_id\n") unless (0 < length($qsub_id));
   if($qsub_id !~ m/^\S+$/) {
      die("Error submitting QSUB_COMMAND=$qsub\nQSUB_ID=$qsub_id\n") unless (0 < length($qsub_id));
   }

   print STDERR "[psolid2fastq] NAME=$outputID QSUBID=$qsub_id\n";

   return $qsub_id;

}

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

   if (eof($FILE)) {
      return -1;
   }

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
   my ($FILE, $readname, $new_readname, $seek_pos, $step_size, $file_size, $exact) = @_;
   my $pos;

   my $last = $#$seek_pos;

   # Determined emperically.  If the estimated position is overshot,
   # most reads require only one step back.
   my $step_back = - min($step_size/256, 4096);

   if (@$seek_pos < 3) {
       $pos = $seek_pos->[$last] + $step_size;
   }
   else {
       # Estimate next position based on last two positions and
       # assuming a quadratic curve.  Empirically, it works better 
       # than not doing it, but it's not perfect.
       my ($a, $b, $c);
       $a = $seek_pos->[$last-1];
       $c = ($seek_pos->[$last-2] + $seek_pos->[$last])/2 - $a;
       $b = $a + $c - $seek_pos->[$last-2];

       $pos = $a + 2*$b + 4*$c;
   }

   # To protect against craziness near the end of a file...
   if ($pos > $file_size + $step_back) {
       $pos = $file_size + $step_back;
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

      $$new_readname = decode($_);
      $pos = tell($FILE) - length;

      if($$new_readname > $readname) {
	 seek($FILE, $step_back, 1);
      }
      else {
	 $pos = tell($FILE) - length;
      }
   } until ($$new_readname <= $readname);

   if (0 == $exact) {
      return $pos;
   }

   do {
      if (m/^>/) {
	 $$new_readname = decode($_);
	 if ($$new_readname == $readname) {
	    $pos = tell($FILE) - length;

	    return $pos;
	 }
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
   $output_dir = rel2abs(dirname($output_prefix));
}
else {
   $output_prefix = rel2abs($output_prefix);
   $basename = fileparse($output_prefix);
   $output_dir = rel2abs(dirname($output_prefix));
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

   my $qual_read_name_1;
   my $qual_read_name_2;

   seek ($csfasta_fhs[0], $chunk * $csfasta_chunk_sizes[0], 0);
   $csfasta_pos[0] = get_next_read($csfasta_fhs[0], \$read_name_1);

   if ($paired) {
      $csfasta_pos[1] = get_matching_read($csfasta_fhs[1], $read_name_1, \$read_name_2,
					  $csfasta_seek_pos[1], $csfasta_chunk_sizes[1],
					  $csfasta_file_sizes[1], 0);

      do {
         # find read_name_2 >= read_name_1 in csfasta_2
         while($read_name_2 < $read_name_1) {
            $csfasta_pos[1] = get_next_read($csfasta_fhs[1], \$read_name_2);
	    if ($csfasta_pos[1] == -1) {
	       print STDERR "Warning! Reached end of $csfasta_files[1] while searching for $read_name_1\n";
	       last;
	    }
         }
         # find read_name_1 >= read_name_2 in csfasta_1
         while($read_name_1 < $read_name_2) {
            $csfasta_pos[0] = get_next_read($csfasta_fhs[0], \$read_name_1);
	    if ($csfasta_pos[0] == -1) {
	       print STDERR "Warning! Reached end of $csfasta_files[0] while searching for $read_name_2\n";
	       last;
	    }
         }
      } while ($read_name_1 != $read_name_2);
   }

   push(@{$csfasta_seek_pos[0]}, $csfasta_pos[0]);
   if ($paired) {
       push(@{$csfasta_seek_pos[1]}, $csfasta_pos[1]);
   }

   # find read_name_1 in qual_1

   $qual_pos[0] = $chunk * $qual_chunk_sizes[0];
   $qual_pos[0] = get_matching_read($qual_fhs[0], $read_name_1, \$qual_read_name_1,
				    $qual_seek_pos[0], $qual_chunk_sizes[0], $qual_file_sizes[0], 1);

   push(@{$qual_seek_pos[0]}, $qual_pos[0]);
   if ($paired) {
      # find read_name_2 in qual_2
      $qual_pos[1] = $chunk * $qual_chunk_sizes[1];
      $qual_pos[1] = get_matching_read($qual_fhs[1], $read_name_1, \$qual_read_name_2,
				       $qual_seek_pos[1], $qual_chunk_sizes[1], $qual_file_sizes[1], 1);
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

   print STDERR "$start_string\n";

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

my @commands = ();
my $command;
my $run_file;
my $job_file;
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
#   $id = "solid2fastq.$basename.$chunk";
#   $run_file = "$RUNDIR/$id.sh";
#   my @a_sub = ();
#   $job_id = SubmitJob($run_file, $dryrun, $command, $id, 1, \@a_sub, 0);
#   push (@job_ids, $job_id);

   push (@commands, $command);
}

$id = "solid2fastq.$basename";
$run_file = "$RUNDIR/$id.sh";
$job_file = "$RUNDIR/$id.jobs";
my @a_sub = ();
$job_id = SubmitArrayJob($run_file, $job_file, $dryrun, \@commands, $id, 1, \@a_sub, $wait && !$merge_single);
push (@a_sub, $job_id);

# if (0 < $wait) {
#    $command = "echo 'Splitting Complete.'";

#    $id = "psolid2fastq.$basename";
#    $run_file = "$RUNDIR/$id.sh";

#    $job_id = SubmitJob($run_file, $dryrun, $command, $id, 1, \@job_ids, 1);
# }

if ($merge_single) {
   my $ext = "fastq";
   if ($gzip_output) {
      $ext .= ".gz";
   }
   $command = "(cat $output_dir/*.single.* > $output_prefix.SINGLE.$ext && rm $output_dir/*.single.* && mv $output_prefix.SINGLE.$ext $output_prefix.single.$ext)";
   $id = "concat.$basename";
   $run_file = "$RUNDIR/$id.sh";
   $job_id = SubmitJob($run_file, $dryrun, $command, $id, 1, \@a_sub, $wait);
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
   -Z                                  Output files are gzip compressed
   --merge-single                      merge single-end read files into one file (-b or -w only)

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

=item B<--merge-single>

bwa (-b) and single (-w) output modes produce many (potentially small) files with single 
end reads.  This option allows them to be merged together after splitting.

=back

=head1 DESCRIPTION

This script takes in approximate split sizes or a chunk count, and determines
seek positions in csfasta and qual files for future splitting.

=cut

