i=69;j=70;ip=$(echo "${i}+1" | bc -l);jp=$(echo "${j}+1" | bc -l);awk -v row=${ip} -v col=${jp} 'NR==row {print $col}' YCn_ZZ_REALMatrixUPPERTRIANGLE_form_SPINUPTL_Hoppings_12x12_OBCXOBC.txt
