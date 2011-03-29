#!/bin/bash

for file in *.raw
do
  stub=${file%.raw}
  ./CleanMinMax.pl "TGCTGCCTCCCGTAGGAGT" $stub < $file
  echo "$stub $i";
done
