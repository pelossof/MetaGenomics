#!/bin/bash
filename=$1
awk -F "\t" 'NR !=2 {{out=$1}{for (i = 2; i <= NF; i+=3) out = out "\t" $i} {print out}}' ${filename} > ${filename}.counts.txt
