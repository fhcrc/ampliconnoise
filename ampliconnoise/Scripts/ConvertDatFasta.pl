#!/usr/bin/perl
use POSIX;
my $flowSeq = "TACG";

@datFiles = @ARGV;

foreach $datFile (@datFiles){
	open (DATFILE, $datFile) or die "Can't open $datFile\n";
	$line = <DATFILE>;
	while($line = <DATFILE>){
		@tokens = split(/ /,$line);

		$id = shift(@tokens);
#		print "@tokens\n";
		my $length = shift(@tokens);
		my $seq = &flowToSeq($length, @tokens);
		$seq =~ s/^TCAG//;
		print ">$id\n$seq\n";
	}
	close(DATFILE);
}

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

