set terminal epslatex size 5,3 standalone color colortext 8
set output 'Fig_bands.tex'


set xlabel '\fontsize{16}{60}\selectfont$$'
set xlabel offset 17,1.5
set ylabel '\fontsize{16}{60}\selectfont$E(meV)$'
set ylabel offset 5.5,2.5

set key at 140,5.6
set key maxrows 13
set key spacing 2
set key width -1.7
set key samplen 1.0

#set label 2 '\fontsize{12}{60}\selectfont$$'
#set label 2 at 0.01,5

#set arrow 1 from 3.3,-1 to 3.3,15 nohead dt 4 lw 2 lc "black"
#set arrow 2 from 2.9,-1 to 2.9,20 nohead dt 4 lw 2 lc "black"


#set arrow 1 from 48,0.6 to 58,0.6 nohead lw 8 lc "black" front

#set arrow 1 from 0.1,11.5 to 0.1,9.5
#set label 6 at 0.2,0.46
#set label 6 '\fontsize{12}{60}\selectfont $\textcolor{black}{(a)\; \textrm{3-orb model, L=16}} $' front


#set label 7 at 2.8,8.0
#set label 7 '\fontsize{16}{60}\selectfont{\color{red} $ ? $}'


#set label 8 at 6,18
#set label 8 '\fontsize{8}{60}\selectfont ZigZag Phase'


#set label 7 at 0.2,0.75
#set label 7 '\fontsize{7}{60}\selectfont$$'

#set xlabel offset 1.35,1.1
#set xlabel '\fontsize{8}{60}\selectfont$i$'
#set format x ""
#set format y "%1.0f"
#set xtics ('0.1' 0.1, '' 0.2, ''0.3, '' 0.4, '0.5' 0.5, '' 0.6, '' 0.7, '' 0.8, '' 0.9, '1.0' 1.0, '2.0' 2.0 , '' 3.0, '' 4.0, '5.0' 5.0, '' 6.0, '' 7.0, '' 8.0, '' 9.0, '10.0' 10.0)

#set xr [0.5:10]
#set yr [*:0.08]

#set logscale x

#set xr [0.2:10]


#set ytics 1.0
set xtics offset 0,-0.5
set ytics offset 0.4,0
#set ytics ( '\fontsize{12}{60}\selectfont$0.0$' 0.0 , '\fontsize{12}{60}\selectfont$0.1$' 0.1 ,'\fontsize{12}{60}\selectfont$0.2$' 0.2  ,'\fontsize{12}{60}\selectfont$0.3$' 0.3, '\fontsize{12}{60}\selectfont$0.4$' 0.4  ,'\fontsize{12}{60}\selectfont$0.5$' 0.5)
#set xtics ( '\fontsize{12}{60}\selectfont$0$' 0 , '\fontsize{12}{60}\selectfont$\pi/2$' 1.5707 ,'\fontsize{12}{60}\selectfont$\pi$' 3.14159  ,'\fontsize{12}{60}\selectfont$3\pi/2$' 4.7123850, '\fontsize{12}{60}\selectfont$2\pi$' 6.28318)
#set cbtics ( '\fontsize{12}{60}\selectfont$0$' 0 ,'\fontsize{12}{60}\selectfont$5$' 5  ,'\fontsize{12}{60}\selectfont$10$' 10, '\fontsize{12}{60}\selectfont$15$' 15  ,'\fontsize{12}{60}\selectfont$20$' 20  )

set ytics ( '\fontsize{16}{60}\selectfont-$6$' -6,  '\fontsize{16}{60}\selectfont-$3$' -3,'\fontsize{16}{60}\selectfont$0$' 0,  '\fontsize{16}{60}\selectfont$3$' 3, '\fontsize{16}{60}\selectfont$6$' 6)

set xtics ( '\fontsize{14}{60}\selectfont$\Gamma$' 0.0, '\fontsize{14}{60}\selectfont$K$' 48, '\fontsize{14}{60}\selectfont$M$' 72, '\fontsize{14}{60}\selectfont$\Gamma$' 144 )
#set ytics offset 0.7

#set pm3d map
#set pm3d corners2color c1
set xr [0.0:144]
set yr [-6:6]
#set cbr [0:20]
#set palette define (0 "white", 0.05 "white", 0.15 "#C5E1EF", 0.25 "#9EC9E2", 0.38 "#3C93C2", 0.6 "#0D4A70", 1.0 "black")

#LzLzkw_OBC_lambda0.0_L36_eta0.02.txt  MzMzkw_OBC_lambda0.0_L36_eta0.02.txt
#sp "MzMzkw_OBC_lambda0.45_L16_eta0.05_3orb_Parzen.txt" u 1:2:3 w pm3d ti ''

#p "Bands_energy_MoSe2_theta_1.5.txt" u 1:6 w l lw 16 lc "red" ti "", "Bands_energy_MoSe2_theta_1.5.txt" u 1:7 w l lw 16 lc "red" ti "",  "Bands_energy_MoSe2_theta_1.5.txt" u 1:(($4-3)*0.17) w l dt 3 lw 8 lc "black" ti "", "Bands_energy_MoSe2_theta_1.5.txt" u 1:(($5-3)*0.17) w l dt 3 lw 8 lc "black" ti ""

p "Bands_Path5_0.0001000000.txt" u ($1*6):4 w p pt 4 ps 2 lc "black" lw 2 ti '\fontsize{14}{60}\selectfont tight-binding model', "Bands_Path5_0.0001000000.txt" u ($1*6):6 w p pt 4 ps 2 lc "black" lw 2  ti "", "Bands_Path5_0.0001000000.txt" u ($1*6):8 w p pt 4 ps 2 lc "black" lw 2  ti "", "Bands_Path5_0.0001000000.txt" u ($1*6):10 w p pt 4 ps 2 lc "black" lw 2  ti "", "Bands_energy_144_kpoints.txt" u 1:($8+36.2) w l lw 6 lc "red" ti "", "Bands_energy_144_kpoints.txt" u 1:($9+36.2) w l lw 6 lc "red" ti "", "Bands_energy_144_kpoints.txt" u 1:($10+36.2) w l lw 6 lc "red" ti "", "Bands_energy_144_kpoints.txt" u 1:($11+36.2) w l lw 6 lc "red" ti '\fontsize{14}{60}\selectfont continuum model'


unset label 7
unset label 8
unset label 6
unset arrow 1
unset arrow 2


