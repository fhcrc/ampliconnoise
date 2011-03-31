#!/usr/bin/perl

use POSIX;

my $fastaFile = shift;
my $split     = shift;
my $offset    = shift;

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

$size = floor(($total - $offset)/$split);
$sizeZero = $size + (($total - $offset) % $split);

for($i = 0; $i < $split; $i++){
	my $sizeJ = 0;
	$dirName = "Split$i";
	print "$dirName\n";
	mkdir $dirName;

	if($i == 0){
		$sizeJ = $sizeZero;
	}
	else{
		$sizeJ = $size;
	}
	my $file = ">${dirName}\/${dirName}.fa";

	open(FILE, $file) or die "Can't open $file\n";
	for($j = 0; $j < $sizeJ; $j++){
		print FILE ">$id[$offset + $j]\n$Seq[$offset+$j]\n"; 
	}
	close(FILE);
	$offset += $sizeJ;
}


