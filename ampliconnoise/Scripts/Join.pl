#!/usr/bin/perl
my $nK = 0;
my $count = 0;
my @lines = ();
foreach $file(@ARGV){
	open (FILE, $file) or die "Can't open $file\n";
	$line = <FILE>;
	if($count == 0){
		chomp($line);
		@tokens = split(/ /,$line);
		$nK = $tokens[1];
	}
	
	while($line = <FILE>){
		chomp($line);

		@tokens = split(/ /,$line);

		shift(@tokens);

		my $temp = "$count @tokens\n";
		push(@lines, $temp);
		$count++;
	}

	close(FILE);
}

print "$count $nK\n";

for($i = 0; $i < $count; $i++)
{
	print "$lines[$i]";
}
