#!/usr/bin/perl
my $tag = shift(@ARGV);
my $newLength = shift(@ARGV);

$tag =~ /G*(.*)/;
$tag2 = $1;

#print "$tag2\n";

while($line = <STDIN>){
    chomp($line);

    if($line =~ />(.*)/){
	print "$line\n";
    }
    else{
	$line =~ /$tag2(.*)/; 	
	$parse = $1;
	$truncate = substr($parse, 0, $newLength);

	print "$truncate\n";
    }
}
