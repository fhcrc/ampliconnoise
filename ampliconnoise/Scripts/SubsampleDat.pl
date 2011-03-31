#!/usr/bin/perl

my $datFile = shift;

my $N = shift;

my @lines = ();
my $total = 0;
my @id = ();

srand(78126723);

open(FILE, $datFile) or die "Can't open $datFile\n";

$line = <FILE>;
chomp($line);

($total, $M) = split(/ /,$line);

for($i = 0; $i < $total; $i++){
	$line = <FILE>;
	chomp($line);
	push(@lines, $line);
}

close(FILE);
my @select = ();
for($i = 0; $i < $total; $i++){
	$select[$i] = 0;
}
$i = 0;
while($i < $N){
	my $random_number = int(rand($total));

	while($select[$random_number] == 1){		
		$random_number = int(rand($total));
	}
	$select[$random_number] = 1;

	$i++;
}

print "$N $M\n";
for($i = 0; $i < $total; $i++){
	if($select[$i] == 1){
		printf("$lines[$i]\n");
	}
}

