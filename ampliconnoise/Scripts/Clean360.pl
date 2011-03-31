#!/usr/bin/perl

use POSIX;

my $nTotal = 0;
my $nClean = 0;
my $flowSeq = "TACG";

my $primer    = $ARGV[0];
my $out       = $ARGV[1];

my $ffile = "${out}.fa";
my $dfile = "${out}.dat";
open(FFILE, ">${ffile}") or die "Can't open ${out}.fa\n";
open(DFILE, ">${dfile}") or die "Can't open ${out}.dat\n";

while(my $line = <STDIN>){
    chomp($line);
 
    if($line =~ />(.*)/){
	my $id = $1;
	#print "$id\n";
	$line = <STDIN>;
	chomp($line);
	#print "$id $line\n";
	my @flows = split(/ /,$line);
	
	#print "@flows\n";
	my $length = shift(@flows);
	my $read;
	my $newLength = 0;

	#print "$length\n";
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
	   
	if($length >= 360 && $read=~/^TCAG.*(${primer}.*)/){
	    $length = 360;
		    
	    printf DFILE "$id $length ";
	    for($j = 0; $j < 360; $j++){
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

close(DFILE);
close(FFILE);

open my $in,  '<',  $dfile      or die "Can't read old file: $!";
open my $out, '>', "$dfile.new" or die "Can't write new file: $!";

print $out "$nClean 360\n";

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
