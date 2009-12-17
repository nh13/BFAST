#!/usr/bin/perl

# Please see the LICENSE accompanying this distribution for 
# # details.  Report all bugs to nhomer@cs.ucla.edu or 
# # bfast-help@lists.sourceforge.net.  For documentation, use the
# # -man option.

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# Resubmits and Eqw job on the cluster, then updates all
# depdent jobs to point to this new job.

my $config;
my ($man, $username, $delete, $help, $jids) = (0, "", 0, 0, 0);
my $version = "0.1.1";
my $CORE = "[bfast resubmit]";

GetOptions('help|?' => \$help,
	'man' => \$man,
	'username=s' => \$username,
	'jids' => \$jids,
	'delete' => \$delete)
	or pod2usage(1);
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(1) if ($help);

my @jids_to_do = ();
my %jids_completed = ();

if($jids && 0 < length($username)) {
	print STDERR "Error: both jids and username options cannot be used.\n";
	exit(1);
}
elsif($jids) {
	foreach my $jid (@ARGV) {
		push(@jids_to_do, $jid);
	}
}
elsif(0 < length($username)) {
	if(0 < scalar(@ARGV)) {
		print STDERR "Error: jids should not be supplied with the username option.\n";
		exit(1);
	}
	@jids_to_do = get_jids_from_username($username);
}
else {
	pod2usage(1);
}

while(0 < scalar(@jids_to_do)) {
	resub(\@jids_to_do, \%jids_completed, $delete);
}


sub get_jids_from_username {
	my $username = shift;

	my $qstat=`qstat`;
	my @lines = split(/\n/, $qstat);
	my @jids = ();
	foreach my $line (@lines) {
		$line =~ s/^\s+//;
		if($line =~ m/Eqw/ && $line =~ m/$username/) {
			if($line =~ m/^(\d+)/) {
				push(@jids, $1);
			}
		}
	}
	return @jids;
}

sub resub {
	my ($jids, $jids_completed, $delete) = @_;
	my $job_name = "";
	my $jid_successors = ();
	my @lines = ();
	my $new_jid = "";
	my $out = "";
	my %params = ();
	my $old_jid = shift(@$jids);

	# ignore previously processed ids
	if(defined($jids_completed{$old_jid})) {
		return;
	}

	printf(STDOUT "%s processing %s\n", $CORE, $old_jid); 

	get_params($old_jid, \%params);

	if(0 == $delete) {
		# Resub with user hold
		$new_jid=`qresub -h u $old_jid`;
		$new_jid =~ s/^.*\n//;
		$new_jid =~ s/^.*job\s+(\S+).*?$/$1/;
		chomp($new_jid);
		printf(STDOUT "%s resubbed with hold on %s\n", $CORE, $new_jid); 

		if(defined($params{'jid_successor_list'})) {
			# Alter the successors
			my @jid_successors = split(/,/, $params{'jid_successor_list'});
			foreach my $jid_successor (@jid_successors) {
				my $new_hold_list = get_new_hold_list($jid_successor, $old_jid, $new_jid);
				$out=`qalter $jid_successor -hold_jid $new_hold_list`;
				printf(STDOUT "%s altered successor %s\n", $CORE, $jid_successor);
			}
		}

		# Delete the old job
		$out=`qdel $old_jid`;
		printf(STDOUT "%s deleted %s\n", $CORE, $old_jid); 

		# Remove user hold
		$out=`qalter $new_jid -h U`;
		printf(STDOUT "%s removed hold on %s\n", $CORE, $new_jid); 
	}
	else {

		if(defined($params{'jid_successor_list'})) {
			# Alter the successors
			my @jid_successors = split(/,/, $params{'jid_successor_list'});
			foreach my $jid_successor (@jid_successors) {
				if(!defined($jids_completed->{$jid_successor})) {
					push(@$jids, $jid_successor);
					printf(STDOUT "%s added successor for deletion %s\n", $CORE, $jid_successor);
				}
				else {
					printf(STDOUT "%s successor is queued for deletion %s\n", $CORE, $jid_successor);
				}
			}
		}

		# Delete the old job
		$out=`qdel $old_jid`;
		printf(STDOUT "%s deleted %s\n", $CORE, $old_jid); 
	}
	$jids_completed->{$old_jid} = 1;
}

sub get_params {
	my ($jid, $params) = @_;
	my $qstat=`qstat -j $jid`;
	my @lines = split(/\n/, $qstat);
	foreach my $line (@lines) {
		if($line =~ m/^(\S+):\s+(\S+)/) {
			$params->{$1} = $2;
		}
	}
}

sub get_new_hold_list {
	my ($jid_successor, $old_jid, $new_jid) = @_;

	my %params = ();
	get_params($jid_successor, \%params);

	# Check that the job parameters were ok
	die unless defined($params{'jid_predecessor_list'});

	my $new_jid_predecessor_list = $params{'jid_predecessor_list'};
	$new_jid_predecessor_list =~ s/^$old_jid$/$new_jid/; # one predecessor
	$new_jid_predecessor_list =~ s/^$old_jid,/$new_jid,/; # first
	$new_jid_predecessor_list =~ s/,$old_jid,/,$new_jid,/; # middle
	$new_jid_predecessor_list =~ s/,$old_jid$/,$new_jid/; # last
	die unless ($new_jid_predecessor_list ne $params{'jid_predecessor_list'});

	return $new_jid_predecessor_list;
}

__END__
=head1 SYNOPSIS

bfast.resubmit.pl [options] <jobids> 

=head1 OPTIONS

=over 8

=item B<-help>
Print a brief help message and exits.

=item B<-man>
Prints the manual page and exits.

=item B<-username>
Process all jobs in the error state from the given username. 

=item B<-jids>
Process all given job ids.

=item B<-delete>
Delete all given jobs as well as their successors.

=back

=head1 DESCRIPTION

B<bfast.resubmit.pl> will resubmit each of the given jobs and modify any 
successor job(s) (i.e. jobs dependencies) to hold on the new resubmitted
job.  This script can process all jobs submitted by a given user, or the
jobs given by the specified job identifiers.  If desired, all jobs and
their successors can be deleted using the B<-delete> option.

Please report all bugs to nhomer@cs.ucla.edu or bfast-help@lists.sourceforge.net.

=head1 REQUIREMENTS

This has only been tested on SGE clusters.  PBS clusters have not been 
tested.

=cut
