#!/usr/bin/perl

my $newLength = shift(@ARGV);

while($line = <STDIN>){
    chomp($line);

    if($line =~ />(.*)/){
	print "$line\n";
    }
    else{
	$truncate = substr($line, 0, $newLength);

	print "$truncate\n";
    }
}
