
txy=-1.0
txztoyz_by_txy=0.0

txz_x_by_txy=1.3
txz_y_by_txz_x=0.0



txy_x=$(printf "%1.2f" ${txy})
txy_y=$(printf "%1.2f" ${txy})

txz_x_val=$(echo "${txz_x_by_txy}*${txy}" | bc -l)
txz_y_val=$(echo "${txz_y_by_txz_x}*${txz_x_val}" | bc -l)
txz_x=$(printf "%1.2f" ${txz_x_val})
txz_y=$(printf "%1.2f" ${txz_y_val})
tyz_x=${txz_y}
tyz_y=${txz_x}

txztoyz_val=$(echo "${txztoyz_by_txy}*${txy}" | bc -l)
txztoyz=$(printf "%1.2f" ${txztoyz_val})


cp t0_mat.txt t0_mat_3orb.txt
cp t1_minus_a2_mat.txt t1_minus_a2_mat_3orb.txt
cp t1_plus_a1_mat.txt t1_plus_a1_mat_3orb.txt
cp t1_plus_a1_minus_a2_mat.txt t1_plus_a1_minus_a2_mat_3orb.txt
cp t1_plus_a1_plus_a2_mat.txt t1_plus_a1_plus_a2_mat_3orb.txt


sed -i "s/txz_x/${txz_x}/" t1_plus_a1_mat_3orb.txt
sed -i "s/tyz_x/${tyz_x}/" t1_plus_a1_mat_3orb.txt
sed -i "s/txy_x/${txy_x}/" t1_plus_a1_mat_3orb.txt

sed -i "s/txz_y/${txz_y}/" t1_minus_a2_mat_3orb.txt
sed -i "s/tyz_y/${tyz_y}/" t1_minus_a2_mat_3orb.txt
sed -i "s/txy_y/${txy_y}/" t1_minus_a2_mat_3orb.txt

sed -i "s/txztoyz/${txztoyz}/" t1_plus_a1_mat_3orb.txt
sed -i "s/txztoyz/${txztoyz}/" t1_minus_a2_mat_3orb.txt
