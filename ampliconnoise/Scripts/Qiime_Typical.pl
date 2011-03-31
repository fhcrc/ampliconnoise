#!/usr/bin/perl

use strict;

my $d         = $ARGV[0];
my $fastaFile = $ARGV[1];
my $listFile  = $ARGV[2];

my @ids     = ();
my @seqs    = ();
my @weights = ();
my @cWeights = ();
my %hashID = {};

open(FILE, $fastaFile) or die "Can't open $fastaFile\n";

my $count = 0;
while(my $line = <FILE>){
    if($line =~ />(.*)/){
	my $name = $1;
	$ids[$count] = $name;
	$ids[$count] =~ /^.*_\d+_(\d+)/;
	#print "$1\n";
	$weights[$count] = $1;
	$hashID{$name} = $count;
	$count++;
    }
    else{
	chomp($line);
	if($seqs[$count - 1] ne undef){
	    $seqs[$count - 1] = $seqs[$count - 1] . $line;
	}
	else{
	    $seqs[$count - 1] = $line;
	}
    }
}

close(FILE);

#for(my $i = 0; $i < $count; $i++){
 #   print "$i $ids[$i] $weights[$i]\n";
#}

my @consensus = ();

open(FILE, $listFile) or die "Can't open $listFile\n";


while(my $line = <FILE>){
    my @tokens = split(/ /,$line);

    my $dLocal = shift(@tokens);
    
    my $N = shift(@tokens);

    if($d == $dLocal){
	for(my $i = 0; $i < $N; $i++){
	    my $clust = shift(@tokens);

	    my @cluster = split(/,/,$clust);
	    my $cSize = scalar(@cluster);

	    my $otuWeight = $weights[$hashID{$cluster[0]}];

	    my $bestWeight = $weights[$hashID{$cluster[0]}];
	    my $bestJ      = $hashID{$cluster[0]};

	    for(my $j = 1; $j < $cSize; $j++){
		if($weights[$hashID{$cluster[$j]}] > $bestWeight){
		    $bestJ = $hashID{$cluster[$j]};
		    $bestWeight = $weights[$hashID{$cluster[$j]}];
		}

		$otuWeight += $weights[$hashID{$cluster[$j]}];
	    } 
	    
	    $consensus[$i] = $seqs[$bestJ];
	    $cWeights[$i] = $otuWeight;
	    

	    #print "$i $consensus[$i]\n";
	}

	for(my $i = 0; $i < $N; $i++){
	    print ">${i}\n";

	    print "$consensus[$i]\n";
	}


    }

}

close(FILE);

