#!/bin/bash

temp_arr=(3213 34 3 534 5 46 54 6 54 6 546 54)
len=${#temp_arr[@]}

echo "${len}"


for (( i=0; i<$len; i++  )); do
	echo "${i}  ${temp_arr[$i]}"
done



val=$( echo "3*5" | bc -l)
echo "${val}"
