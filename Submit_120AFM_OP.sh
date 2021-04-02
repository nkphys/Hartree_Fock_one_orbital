
m=0.0
file_energy_out="energy_vs_m_U8_A1_120AFM_XZ.txt"
echo "#m energies .." > ${file_energy_out}
for i in {0..100}
do
m=$(echo "(0.5/100)*${i}" | bc -l)
#m=0.44294040

#sed -i -e "s/N_TOTAL_VALUE/${N}/g" $input


OP_input="OUTPUT_OP_ZX120_run.txt"
rm $OP_input
cp OUTPUT_OP_ZX120_template.txt $OP_input

nup0=$(echo "(1.0 + (2*${m}))/2.0" | bc -l)
ndn0=$(echo "(1.0 - (2*${m}))/2.0" | bc -l)
cupdcdn0=$(echo "0" | bc -l)

nup1=$(echo "(1.0 - (${m}))/2.0" | bc -l)
ndn1=$(echo "(1.0 + (${m}))/2.0" | bc -l)
cupdcdn1=$(echo "-1.0*(sqrt(3.0)/2.0)*${m}" | bc -l)

nup2=$(echo "(1.0 - (${m}))/2.0" | bc -l)
ndn2=$(echo "(1.0 + (${m}))/2.0" | bc -l)
cupdcdn2=$(echo "1.0*(sqrt(3.0)/2.0)*${m}" | bc -l)

##############
sed -i -e "s/n_up_0/(${nup0},0)/g" $OP_input
sed -i -e "s/n_dn_0/(${ndn0},0)/g" $OP_input
sed -i -e "s/n_up_1/(${nup1},0)/g" $OP_input
sed -i -e "s/n_dn_1/(${ndn1},0)/g" $OP_input
sed -i -e "s/n_up_2/(${nup2},0)/g" $OP_input
sed -i -e "s/n_dn_2/(${ndn2},0)/g" $OP_input
sed -i -e "s/cupD_cdn_0/(${cupdcdn0},0)/g" $OP_input
sed -i -e "s/cupD_cdn_1/(${cupdcdn1},0)/g" $OP_input
sed -i -e "s/cupD_cdn_2/(${cupdcdn2},0)/g" $OP_input
###############


time ./k_space_SelfConsistency TriangularLattice modelinput_TL_k_space_SC.inp > out_${m}.txt

val=$(awk "NR==55" out_${m}.txt)
echo "${m}   ${val}" >> ${file_energy_out}
echo "m=${m} done"
done

