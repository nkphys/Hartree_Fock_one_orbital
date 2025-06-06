set term epslatex standalone color 
unset key
# unset xtics
# unset ytics

#set cbrange [-1.0:1.0] 
set output 'plot_TL.tex'

#A1_Len=60
#A2_Len=200

#set border 4095 lw 3 

#unset border 
#unset xtics
#unset ytics

#n1y=sqrt(3.0)*0.5*A2_Len
#n1x=0.5*A2_Len
#n2x=n1x+A1_Len
#n3x=n2x-n1x

#x_off=4
#y_off=4

#xrange_val_min=0-x_off
#xrange_val_max=n2x+x_off

#yrange_val_min=0-y_off
#yrange_val_max=n1y+y_off

#set xrange [xrange_val_min:xrange_val_max]
#set yrange [yrange_val_min:yrange_val_max]
#set yrange [xrange_val_min:xrange_val_max]

#ratio=(yrange_val_max-yrange_val_min+4.7)/(xrange_val_max-xrange_val_min)

ratio=1.0
#ratio=2.0

set size 1,ratio
set style arrow 1 head filled size screen 0.0,20 lw 10 lc pal 
set style arrow 3 head filled size screen 0.1,1,1  lw 1 lc pal 


#set arrow 4 nohead from 0,0 to n1x,n1y lw 2 lc "black"
#set arrow 5 nohead from n1x,n1y to n2x,n1y lw 2 lc "black"
#set arrow 6 nohead from n2x,n1y to n3x,0 lw 2 lc "black"
#set arrow 7 nohead from n3x,0 to 0,0 lw 2 lc "black"

scale = 0.7 
#p 'MySkyrmion.txt' u (( $1 + (0.5*$2)  )-0.4*$4):(( sqrt(3.0)*0.5*($2) )-0.4*$5):(($4)*scale):(($5)*scale):(($3)) w vec arrowstyle 3 notitle

#set yr [-0.7:0.7]
#set xr [0:2]
unset cbr
set cbr [-0.1:0.1]
#set palette define (0 "red", 0.3 "red", 1.0 "blue")
set xr [-2:16]
set yr [-2:16]
set palette define (0 "red", 0.3 "red", 1.0 "blue")
#set palette 2 define (0 "grey", 0.3 "grey", 1.0 "black")

set style fill transparent solid 0.6 noborder
#p "layer_1_RealSpace_OParams.txt" u ($3):($4):(($5)*50.0):($5) with circles lc palette notitle
#p "Local_spins0.0001000000.txt" u ($5):($6):((($11))*0.1):($11) with circles lc palette notitle

p "Local_spins0.0001000000.txt" u ($5):($6):($8*3):($9*3):(($9)*4) w vec arrowstyle 3 notitle, "Local_spins0.0001000000.txt" u ($5):($6):((($11))*0.14):($11) with circles lc palette notitle
#sp "layer_0_Temp_0.0001000000RealSpace_OParams.txt" u ($5):($6):((($11))*50.0) w pm3d map lc palette notitle


#p "RealSpace_OParams.txt" u ($3):($4):($6*8000*1):($7*8000*1):($5) w vec arrowstyle 3 notitle 
