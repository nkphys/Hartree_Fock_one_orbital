mkdir -p Data
theta=4.0
ref_site=9

for eps in 2.0 #4.5 6.5 7.0 7.5 11.0 13.0 20.0 #2.5 3.5 4.5 5.0 5.5 6.5 7.0 7.5 8.5 9.0 9.5 10.5 11.0 12.0 13.0 14.0 #1.0 2.0 4.0 6.0 8.0 10.0 15.0 18.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 70.0 80.0 90.0 100.0
do

file_out="SS_refsite${ref_site}_theta${theta}_eps${eps}.txt"
rm Data/${file_out}

#file=theta${theta}/eps${eps}/out.txt
file=out.txt
lineSzSz=$(awk '/SzSz/ {print NR}' ${file})
lineSzSz=$(echo "${lineSzSz}+2+${ref_site}" | bc -l)
lineSpSm=$(awk '/SpSm/ {print NR}' ${file})
lineSpSm=$(echo "${lineSpSm}+2+${ref_site}" | bc -l)
lineSmSp=$(awk '/SmSp/ {print NR}' ${file})
lineSmSp=$(echo "${lineSmSp}+2+${ref_site}" | bc -l)
for i in {1..16}
do
val_SzSz=$(awk -v row=${lineSzSz} -v col=${i} 'NR==row {print $col}' ${file} )
val_SzSz=$(echo "${val_SzSz}" | cut -f1 -d "," | tr -d "(");val_SzSz=$(printf "%1.10f" ${val_SzSz})

val_SpSm=$(awk -v row=${lineSpSm} -v col=${i} 'NR==row {print $col}' ${file} )
val_SpSm=$(echo "${val_SpSm}" | cut -f1 -d "," | tr -d "(");val_SpSm=$(printf "%1.10f" ${val_SpSm})

val_SmSp=$(awk -v row=${lineSmSp} -v col=${i} 'NR==row {print $col}' ${file} )
val_SmSp=$(echo "${val_SmSp}" | cut -f1 -d "," | tr -d "(");val_SmSp=$(printf "%1.10f" ${val_SmSp})

val=$(echo "${val_SzSz} + 0.5*(${val_SpSm} + ${val_SmSp})" | bc -l)

site=$(echo "${i}-1" | bc -l)
echo "0 0 ${site} ${val} 0" >> Data/${file_out}
done


done #eps
