#!/bin/bash


defaultBarCode="ACACACGTCG" #default primer

fastafile=$1; 
barcode=${2:-$defaultBarCode} 
#second argument primer as a Perl regular expression

stub=${fastafile%.fa};
parseFile=${stub}_P.fa
echo $stub $fastafile $barcode

sed "s/^${barcode}//" $fastafile > $parseFile

Perseus -sin $parseFile > ${stub}_P.per

Class.pl ${stub}_P.per -7.5 0.5  > ${stub}_P.class

FilterGoodClass.pl ${stub}_P.fa ${stub}_P.class 0.5 2> ${stub}_P_Good.fa > ${stub}_P_Chi.fa
