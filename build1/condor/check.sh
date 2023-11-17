#!/bin/bash

dir="./out/"
str="Good bye"
n=0
for file in $(ls  ${dir}*)
	do 

	if [ `grep -c "${str}" $file` -eq 0 ];then
	echo "Error in $file: Does NOT Finished!!!"
	((n++))
	fi
done

echo "File number of NOT Finished: ${n}"
