#/bin/bash

if [ "$#" -le 1 ]; then
    echo "usage: mergemaf.sh inputpath outputfile"
    exit 1
fi

> $2

firstfile=`ls $1/*.maf.txt | head -1`
head -1 $firstfile >> $2

for f in $1/*.maf.txt
do 
	tail -n +2 "$f" >> $2
done

