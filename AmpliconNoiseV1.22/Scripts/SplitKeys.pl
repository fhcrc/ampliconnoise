#!/usr/bin/perl
use POSIX; 
use IO::File;

my $nReads = 0;
my $nFlows = 400;

my %counts = {};

my $bGood   = 0;

my %keys = {};

my $count = 0;
my $fileName = 0;

my $primer  = $ARGV[0]; 
my $keyFile = $ARGV[1];

open(FILE, "$keyFile") or die;

while($line = <FILE>){
    chomp($line);

    @tokens = split(/,/,$line);

    my $tag = $tokens[1];

    $fileName = $tokens[0].".raw";
    #my $fileName = "Seq${count}.fa";

    my $fh = IO::File->new();

    open($fh, ">$fileName") or die "Can't open $fileName\n";

    print $fh "$count $tag $fileName\n";

    $keys{$tag} = $fh;

    #close($fh);

    $count++;
}
close(FILE);

my $flowSeq = "TACG";


while(my $line = <STDIN>){
    chomp($line);
 
    if($line =~ /\# of Reads:\s+(\d+)/){
	$nReads = $1;
    }

    if($line =~ /\# of Flows:\s+(\d+)/){
	$nFlows = $1;
    }

    if($line =~ />(.*)/){
	$id = $1;
    }

    if($line =~ /Clip Qual Right:\s+(\d+)/){
	my $clipright = $1;
	my $flowRef  = -1;
	my $length   = -1;
	my $sequence = -1;

	$bGood = 1;

	$line = <STDIN>; $line = <STDIN>; $line = <STDIN>; $line = <STDIN>; 
	chomp($line);
	if($line =~ /Flowgram:(.*)/){
	    
	    my @tokens = split(/\s+/,$1);
	
	    shift(@tokens);

	    my $test = scalar(@tokens);
	
	    if($test != $nFlows){
		$bGood = 0;
#print "Problem... $test $nFlows\n";
	    }
	    else{
		$flowRef = \@tokens;
	    }
	}
	else{
	    $bGood = 0;
	    #printf "Format error $line\n";
	}
	$line = <STDIN>; chomp($line);
	if($line =~ /Flow Indexes:(.*)/){
	    my @tokens = split(/\s+/,$1);
	
	    shift(@tokens);
	    
	    $length = $tokens[$clipright - 1];
	    #print "$clipright $length\n";
	    
	}
	else{
	    $bGood = 0;
	    #printf "Format error $line\n";
	}
	$line = <STDIN>; chomp($line);
	if($line =~ /Bases:\s+(\w+)$/){
#	    print "$1\n";
	    $sequence = $1;
#
	}
	else{
	    $bGood = 0;
	    #printf "Format error $line\n";
	}

	my $read = flowToSeq($length,@{$flowRef});
	#print "$bGood $read\n";
	if($bGood == 1 && $read =~ /^TCAG(\w*)$primer.*$/){
	    my $tag = $1;
	    if($keys{$tag} ne undef){
		my $fh = $keys{$tag};
		
		print $fh ">$id\n$length @{$flowRef}\n";

		$counts{$tag}++;
	    }
	    else{
		print STDERR ">$id\n$sequence\n";
	    }
	}
#	else{
	#    printf STDERR "Format error\n";
	#}

    }
}

my $total = 0;
foreach $key (keys %counts){
    my $count = $counts{$key};
    if($count > 0){
	print "$key $count\n";
	$total += $count;
    }
}
print "$nReads $total\n";

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
