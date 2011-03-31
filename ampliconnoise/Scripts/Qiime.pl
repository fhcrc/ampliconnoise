#!/usr/bin/perl

$listFile = shift(@ARGV);

$fastaFile = shift(@ARGV);

$d = shift(@ARGV);

my %fastaMap = {};
my $count = 1;

open(FASTAFILE, $fastaFile) or die "Can't open $fastaFile\n";

while($line = <FASTAFILE>){
	chomp($line);
	if($line=~/>(.*)/){
		$fastaMap{$1} = $count;
		$count++;
	}
}

close(FASTAFILE);

@mapFiles = @ARGV;
my %readMap = {};
foreach $mapFile(@mapFiles)
{
	my @map = ();
	open (MAPFILE, $mapFile) or die "Can't open $mapFile\n";
	while($line = <MAPFILE>){
		chomp($line);
		my @stokens = split(/ /,$line);
		my @ctokens = split(/,/,$stokens[1]);
		my @maptokens = ();
		#$stokens[1]=~tr/,/\t/;
		foreach $ctok(@ctokens){
#			print "$ctok $fastaMap{$ctok}\n";
			push(@maptokens,$fastaMap{$ctok});
		}
		push(@map,\@maptokens);
	}
	close(MAPFILE);
	$mapFile=~/${mapDir}\/(.*)_S1000_s60_c01_T220_s30_c08_cd.mapping/;
	$sample=$1;
	$readMap{$sample} = \@map;
#	print "$mapFile\n";
}

#foreach $sample (keys %readMap){
#	print "$sample\n";

#	my @array = @{$readMap{$sample}};
#	my $tsize = scalar(@array);

#	for($i = 0; $i < $tsize; $i++){
#		my @finalidx = @{$array[$i]};
#		print "$i @finalidx\n";
#	}
#}
open(FILE, $listFile) or die "Can't open $file\n";

while($line = <FILE>){
        chomp($line);

        @tokens = split(/ /,$line);

        $locald = shift(@tokens);
        
        if($d == $locald){
		$N = shift(@tokens);

		for($i = 0; $i < $N; $i++){
			print "$i\t";		
			my @clust = split(/,/,$tokens[$i]);
			my @temp = ();
			foreach $read (@clust){
				$read=~/(TS\d+).*_(\d+)_\d+/;
				my $sample = $1;
				my $idx = $2;
				my @array = @{$readMap{$sample}};

				my @finalidx = @{$array[$idx]};
#				print "$sample $idx @finalidx\n";
				foreach my $c(@finalidx){
#					print "$sample\n";
					
					push(@temp,"${sample}_${c}");
				}
			}
			my $tempString = join("\t",@temp);
			print "$tempString\n";
		}
                
        }
}

close(FILE);

