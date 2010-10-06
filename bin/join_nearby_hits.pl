#!/usr/bin/perl -w

my $cutoff = 4000;

my $lastcontig = "no";
my $laststart; 
my $lastend;
my $lastmotif;

my $joincount = 0;
my $totals = 0;

$outfile = $ARGV[0];
open OUT, ">".$outfile;

while(my $r = <STDIN>) {

    chomp $r; 

    if( ($contig, $start, $end,$motif) = ($r =~ /^(\S+)\t(\d+)\t(\d+)\t(\S+)/) ) {

	$totals++;
	if (($contig eq $lastcontig) and ($motif ne $lastmotif) ) {
	    #check distance
	    if ( $start - $lastend < $cutoff or $end - $laststart < $cutoff) {
		print "$contig $start $end $motif JOINS $lastcontig $laststart $lastend $lastmotif\n";
		print OUT $lastcontig."_".$laststart."-".$lastend."\t".$lastmotif."\n";
		print OUT $contig."_".$start."-".$end."\t".$motif."\n";
		$joincount++;
	    }
	}

	$lastcontig =$contig;
	$laststart=$start;
	$lastend = $end;
	$lastmotif = $motif;
    }
    
}
close OUT;
print "---> joins $joincount of $totals lines\n";
