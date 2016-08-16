#!/bin/bash
filename=$1
awk -F "\t" 'NR !=2 {print $1"\t"$3"\t"$4}' ${filename} > ${filename}.genes.txt &
awk -F "\t" 'NR !=2 {{out=$1}{for (i = 2; i <= NF; i+=4) out = out "\t" $i} {print out}}' ${filename} > ${filename}.meth.txt
