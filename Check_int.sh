for site in {0..35}
do
grep "^${site} " dat.txt > Site${site}_DS.txt
E=0.0;
for i in {1..35}
do 
val=$(awk -v row=${i} 'NR==row {print $3}' Site${site}_DS.txt)
E=$(echo "${val} + ${E}" | bc -l)
done
E=$(echo "0.0*960 + ${E}" | bc -l)
echo "${site}  ${E}"
done
