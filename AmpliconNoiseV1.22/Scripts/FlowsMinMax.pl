#!/usr/bin/perl
use POSIX;
my $nTotal = 0;
my $nClean = 0;
my $id;
my $bGood   = 0;
my $flowSeq = "TACG";

my $primer    = $ARGV[0];

my $out       = $ARGV[1];

my $minFlows = 360;

my $maxFlows = 720;

my $ffile = "${out}.fa";
my $dfile = "${out}.dat";
open(FFILE, ">${ffile}") or die "Can't open ${out}.fa\n";
open(DFILE, ">${dfile}") or die "Can't open ${out}.dat\n";

while(my $line = <STDIN>){
    chomp($line);
 
    if($line =~ /\# of Reads:\s+(\d+)/){
	$nReads = $1;
    }

    if($line =~ /\# of Flows:\s+(\d+)/){
	$nFlows = $1;
    }

    if($line =~ />(.*)/){
	$id = $1;
    }

    if($line =~ /Clip Qual Right:\s+(\d+)/){
	my $clipright = $1;
	my $flowRef  = -1;
	my $length   = -1;
	my $sequence = -1;

	$bGood = 1;

	$line = <STDIN>; $line = <STDIN>; $line = <STDIN>; $line = <STDIN>; 
	chomp($line);
	if($line =~ /Flowgram:(.*)/){
	    
	    my @tokens = split(/\s+/,$1);
	
	    shift(@tokens);

	    my $test = scalar(@tokens);
	
	    if($test != $nFlows){
		$bGood = 0;
#print "Problem... $test $nFlows\n";
	    }
	    else{
		$flowRef = \@tokens;
	    }
	}
	else{
	    $bGood = 0;
	    #printf "Format error $line\n";
	}
	$line = <STDIN>; chomp($line);
	if($line =~ /Flow Indexes:(.*)/){
	    my @tokens = split(/\s+/,$1);
	
	    shift(@tokens);
	    
	    $length = $tokens[$clipright - 1];
	    #print "$clipright $length\n";
	    
	}
	else{
	    $bGood = 0;
	    #printf "Format error $line\n";
	}
	$line = <STDIN>; chomp($line);
	if($line =~ /Bases:\s+(\w+)$/){
#	    print "$1\n";
	    $sequence = $1;
#
	}
	else{
	    $bGood = 0;
	    #printf "Format error $line\n";
	}

	if($bGood == 1){
	    my @flows = @{$flowRef};
	    my $read;
	    my $newLength = 0;
	    #print "$length @flows\n";
	    while($newLength*4 < $length){
		my $signal = 0;
		my $noise  = 0;

		for($j = 0; $j < 4; $j++){
		    my $f = $flows[$j + 4*$newLength];
		    if($f > 0.5){
			$signal++;
			if($f < 0.7 || $f > 9.49){
			    $noise++;
			}
		    }
	    
		}
    
		if($noise > 0 || $signal == 0){
		    last;
		}

		$newLength++;
	    }

	    $length = $newLength*4;

	    my $read = flowToSeq($length,@flows);
	#    print "$length $read\n";
	    if($length >= $minFlows && $read=~/^TCAG(${primer}.*)/){
		
		if($length > $maxFlows){
		    $length = $maxFlows;
		}
		    
		printf DFILE "$id $length ";
		for($j = 0; $j < $maxFlows; $j++){
		    printf DFILE "$flows[$j] ";
		}
		
		printf DFILE "\n";
		    

		
		printf FFILE ">$id\n";
		printf FFILE "$1\n";
		
		$nClean++;
		goto found;
	    }

	    found : $nTotal++;
	}

    }
}

close(DFILE);
close(FFILE);

open my $in,  '<',  $dfile      or die "Can't read old file: $!";
open my $out, '>', "$dfile.new" or die "Can't write new file: $!";

print $out "$nClean $maxFlows\n";

while( <$in> )
{     
    print $out $_;
}

close $out;

rename("$dfile.new", $dfile);

sub flowToSeq()
{
    my ($length, @flowgram) = @_;
    my $retSeq = "";

    for(my $i = 0; $i < $length; $i++){
	my $signal = floor($flowgram[$i] + 0.5);
	my $base   = substr($flowSeq, $i % 4, 1);

	for(my $j = 0; $j < $signal; $j++){
	    $retSeq .= $base;
	}
    }

    return $retSeq;
}
