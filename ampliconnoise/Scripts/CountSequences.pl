#!/usr/bin/perl

$total = 0;

while($line = <STDIN>){
    chomp($line);

    if($line=~/>(.*)/){
	$total ++;
    }
}
if($ARGV[0] eq "c"){
	print "$total,";
}
else{
	print "$total\n";
}
