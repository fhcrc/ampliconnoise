#!/usr/bin/perl

my $fastaFile = shift;
my $distFile = shift;

my $cutoff = shift;

my @Seq = ();
my $total = 0;
my @id = ();

open(FILE, $fastaFile) or die "Can't open $fastaFile\n";

my $seq = "";

while($line = <FILE>){ 
    chomp($line);
    
    if($line =~ />(.*)/){
	
	$id[$count] = $1;
	
	if($seq ne ""){
	    $Seq[$count - 1] = $seq;

	    $seq = "";
	}

	$count++;
    }
    else{
	$seq .= $line;
    }
}

$Seq[$count - 1] = $seq;

$total = $count;

open(FILE, $distFile) or die "Can't open $distFile\n";

for($i = 0; $i < $total; $i++){
    $line = <FILE>;
    
    chomp($line);

    @tokens = split(/ /,$line);

    $id = $tokens[0];

    $id =~ /.*_(.*)/;

    $freq = $1;
    
    $x = $tokens[1]; 
    $y = $tokens[2];

    $p = $tokens[4];

    if($p < $cutoff){
	printf STDERR ">$id[$i]\n$Seq[$i]\n";
    }
    else{
	printf STDOUT ">$id[$i]\n$Seq[$i]\n";
    }
}

close(FILE);
