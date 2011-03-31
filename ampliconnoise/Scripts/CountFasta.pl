#!/usr/bin/perl

$total = 0;

while($line = <STDIN>){
    chomp($line);

    if($line=~/>(.*)/){
	$id = $1;
	@tokens = split(/_/,$id);
	$size = scalar(@tokens);
	$total += $tokens[$size - 1];
    }
}
if($ARGV[0] eq "c"){
	print "$total,";
}
else{
	print "$total\n";
}
