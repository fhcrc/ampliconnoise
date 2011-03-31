#!/bin/bash

defaultPrimer="ATTAGATACCC[ACTG]GGTAG" #default primer
nodes=8                              #no. of cluster nodes to use

sfffile=$1; #first argument name of sff file (necessary)
primer=${2:-$defaultPrimer} #second argument primer as a Perl regular expression

stub=${sfffile%.sff};

echo $stub $sfffile $primer

# first generate sff text file if necessary
if [ ! -f ${sfffile}.txt ]; then
    echo "Generating .sff.txt file"
    sffinfo $sfffile > ${sfffile}.txt
fi

#generate flowgram and fasta files
if [ ! -f ${stub}.dat ]; then
    echo "Generating .dat file"
    FlowsFA360.pl $primer $stub < ${sfffile}.txt
fi

#get unique sequences
if [ ! -f ${stub}_U.fa ]; then
    echo "Getting unique sequences"
    FastaUnique -in ${stub}.fa > ${stub}_U.fa
fi

#use NDist to get sequence distances 
if [ ! -f ${stub}_U_I.ndist ]; then
    echo "Calculating sequence distances"
    mpirun -np $nodes NDist -i -in ${stub}_U.fa > ${stub}_U_I.ndist
fi

#use NDist to get sequence distances 
if [ ! -f ${stub}_U_I.list ]; then
    echo "Cluster sequences..";
#cluster sequences using average linkage and sequence weights
    FCluster -a -w -in ${stub}_U_I.ndist -out ${stub}_U_I > ${stub}_U_I.fcout
fi

rm ${stub}_U_I.ndist

SplitClusterEven -din ${stub}.dat -min ${stub}.map -tin ${stub}_U_I.tree -s 5000 -m 1000 > ${stub}_split.stats

echo "Calculating .fdist files"
for c in C*
do
        if [ -d $c ] ; then
                mpirun -np $nodes PyroDist -in ${c}/${c}.dat -out ${c}/${c} > ${c}/${c}.fout
        fi
done

echo "Clustering .fdist files"

for c in C*
do
        if [ -d $c ] ; then
                FCluster -in ${c}/${c}.fdist -out ${c}/${c}_X > ${c}/${c}.fout
		rm ${c}/${c}.fdist
        fi
done

echo "Running PyroNoise"
for dir in C*
do
        if [ -d $dir ] ; then
                mpirun -np $nodes PyroNoiseM -din ${dir}/${dir}.dat -out ${dir}/${dir}_s60_c01 -lin ${dir}/${dir}_X.list -s 60.0 -c 0.01 > ${dir}/${dir}_s60_c01.pout
        fi
done

cat C*/C*_s60_c01_cd.fa > All_s60_c01_cd.fa
cat C*/C*_s60_c01.mapping > All_s60_c01.mapping

Truncate.pl 220 < All_s60_c01_cd.fa > All_s60_c01_T220.fa

echo "Running SeqDist"
mpirun -np $nodes SeqDist -in All_s60_c01_T220.fa > All_s60_c01_T220.seqdist

FCluster -in All_s60_c01_T220.seqdist -out All_s60_c01_T220_S > All_s60_c01_T220.fcout

echo "Running SeqNoise"
mpirun -np $nodes SeqNoise -in All_s60_c01_T220.fa -din All_s60_c01_T220.seqdist -lin All_s60_c01_T220_S.list -out All_s60_c01_T220_s30_c08 -s 30.0 -c 0.08 -min All_s60_c01.mapping > All_s60_c01_T220_s30_c08.snout

rm All_s60_c01_T220.seqdist
