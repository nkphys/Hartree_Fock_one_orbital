
for band in {1..19} {21..28} # {1..30}
do
i=$(echo "2 + (${band}-1)*7" | bc -l)
#echo "${i}"

col=$(echo "${i}+1" | bc -l)
string="\$${col}"
for j in 3 5
do
col=$(echo "${i}+${j}" | bc -l)
string="${string}+\$${col}"
done
for j in 2 4 6
do
col=$(echo "${i}+${j}" | bc -l)
string="${string}-\$${col}"
done

#echo "${i}  ${string}"

printf "\"Bands0.0001000000.txt\" u 1:${i}:(${string}) w lp palette pt 4 lw 4 ti \"${band}\", "
#printf "\"Bands0.0001000000.txt\" u 1:${i}:(${string}) w lp palette pt 4 lw 4 ti \"${band}\", "
done

echo ""
