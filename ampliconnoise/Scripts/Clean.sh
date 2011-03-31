#!/bin/bash

for file in *.raw
do
  stub=${file%.raw}
  ./Clean360.pl "TGCTGCCTCCCGTAGGAGT" $stub < $file
  echo "$stub $i";
done
