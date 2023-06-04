#!/bin/bash

echo "chr1	2	100" > Test.bed

for file in *.fa
	
do	
	bedtools getfasta -fi $file -bed Test.bed -fo Test.whatawaste
	

done

rm Test.whatawaste
rm Test.bed