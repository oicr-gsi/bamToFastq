#!/bin/bash
cd $1

files=$( ls | grep -v samstats.txt)
for file in $files
do 
    md5sums+=$(md5sum $file | awk '{print $1" "}')
done

printf '%s\n' ${md5sums[*]} | sort
