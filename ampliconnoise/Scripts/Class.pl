#!/usr/bin/perl

my $distFile = shift;

my $alpha = shift;
my $beta =  shift;

my $total = 0;

open(FILE, $distFile) or die "Can't open $distFile\n";

$line = <FILE>;
chomp($line);
@tokens = split(/ /,$line);
$total = $tokens[0];

for($i = 0; $i < $total; $i++){
    $line = <FILE>;
    
    chomp($line);

    @tokens = split(/ /,$line);

    $id = $tokens[1];

    $id =~ /.*_(.*)/;

    $freq = $1;
    
    $x = $tokens[11]; 
    $y = $tokens[12];
    $z = $tokens[13];

    if($x >= 0.15 || $y > 0.0){
	$p = 0.0;
    }
    else{
	$r = $alpha + $beta*$z;
	$p = 1.0/(1.0 + exp(-$r));
    }
    
    print "$id $x $y $z $p\n";
}

close(FILE);
