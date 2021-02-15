mkdir -p Data
ref_site=9

file_out="denden_refsite${ref_site}.txt"
rm Data/${file_out}

file=
line=$(awk '/denden/ {print NR}' ${file})
line=$(echo "${line}+2+${ref_site}" | bc -l)

#echo "line=${line}"

for i in {1..16}
do
val=$(awk -v row=${line} -v col=${i} 'NR==row {print $col}' ${file} )

#echo "${val}"
val=$(echo "${val}" | cut -f1 -d "," | tr -d "(");val=$(printf "%1.10f" ${val})

site=$(echo "${i}-1" | bc -l)
echo "0 0 ${site} ${val} 0" >> Data/${file_out}
done


