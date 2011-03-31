#!/usr/bin/perl

use strict;

my $cutoff = $ARGV[0];
my $taxonomyFile = $ARGV[1];
my $suffix = $ARGV[2];

my $minFreq = 0;

my $nOTUs  = 0;
my $minval = 0;
my %taxHash = {};
open(TAXFILE, $taxonomyFile) or die "Can't open $taxonomyFile\n";

while(my $line = <TAXFILE>){
	chomp($line);
	my @tokens = split(/\t/,$line);

	$taxHash{$tokens[0]} = $tokens[1]; 
}

close(TAXFILE);


my %OTUVectors = {};

while(my $line = <STDIN>){
    my @tokens = split(/ /,$line);
    my $d = shift(@tokens);

    if($d == $cutoff){
	$nOTUs = shift(@tokens);

	for(my $i = 0; $i < $nOTUs; $i++){
	    my %hashCluster = {};

	    my @Cluster = split(/,/,$tokens[$i]);
	    foreach my $entry(@Cluster){
		$entry =~ /^(.*)_\d+_(\d+)$/;
		if($1 ne ""){
		    $hashCluster{$1} += $2;
		}
	    }

	    foreach my $sample(keys %hashCluster){
		my $ref = $OTUVectors{$sample};
		if($hashCluster{$sample} > 0){
		    #print "good $sample\n";
		    if($ref ne undef){
			${$OTUVectors{$sample}}[$i] = $hashCluster{$sample}; 
		    }
		    else{
			my @vector = ();
			$vector[$i] = $hashCluster{$sample};
			$OTUVectors{$sample} = \@vector; 
		    }
		}
		else{
		   # print "bad $sample\n";
		}
	    }
	}
    }
}

my @Vectors = ();
my @names = ();
my @totals = ();
my @goodKeys = ();
foreach my $testKey(keys %OTUVectors){
       
	if($OTUVectors{$testKey} ne undef){
		push(@goodKeys,$testKey);
	}
}

@goodKeys = sort {keySort()} (@goodKeys);

foreach my $sample(@goodKeys){
    if($OTUVectors{$sample} ne undef){
	my @vector = @{$OTUVectors{$sample}};
		
	for(my $i = 0; $i < $nOTUs; $i++){
	    if($vector[$i] == undef){
		$vector[$i] = 0;
	    }
	}
	
	push(@Vectors, \@vector);
	push(@names, $sample);
    }
}

my $nSamples = scalar(@Vectors);

my @print = ();

for(my $j = 0; $j < $nOTUs; $j++){

    $print[$j] = 0;

    for(my $i = 0; $i < $nSamples; $i++){
	if($Vectors[$i][$j] > $minFreq){
	    $print[$j] = 1;
	}
    }

}

for(my $i = 0; $i < $nSamples; $i++){
	$totals[$i] = 0.0;

	for(my $j = 0; $j < $nOTUs; $j++){
		if($print[$j] == 1){
			$totals[$i] += $Vectors[$i][$j];
		}	
	}
}

print "#Full OTU Counts\n";
print "#OTU ID\t";

my $stringSamples = join("\t",@goodKeys);

print "$stringSamples";

print "\tConsensus Lineage\n";

my @printOTUs = ();

for(my $j = 0; $j < $nOTUs; $j++){
	if($print[$j] == 1){
	  printf("%d\t",$j);
	  my @temp = ();

	  for(my $i = 0; $i < $nSamples; $i++){
	    push(@temp, $Vectors[$i][$j]);
	  }

	  my $tempString = join("\t", @temp);
	  printf("$tempString\t$taxHash{$j}\n");
	}	
}


sub keySort
{
    $a=~/${suffix}(\d+)/; 
    my $a1 = $1;
     
    $b=~/${suffix}(\d+)/; 
    my $b1 = $1;
	#print "$a $a1 $b $b1\n"; 
    return $a1 <=> $b1;   
}

sub min
{
    my ($x, $y) = @_;

    if($x < $y){
	return $x;
    }
    else{
	return $y;
    }
}
