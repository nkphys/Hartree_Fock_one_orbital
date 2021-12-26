#! /usr/bin/env wolframscript

data = Cases[Import["/home/nitin/Documents/Codes/Hartree_Fock_one_orbital/Local_spins0.0001000000.txt","Table"], {_?NumberQ, ___}];
#Export["check.dat", data]

g1 = Graphics3D[{RGBColor[4.0*#[[6]], 4.0*#[[7]],4.0*#[[8]]], Opacity[1],Arrowheads[0.1*(Sqrt[#[[7]]*#[[7]]+#[[8]]*#[[8]]+#[[6]]*#[[6]]])],Arrow[Tube[{{#[[4]],#[[5]],0}, {#[[4]] + 4.0*#[[7]],#[[5]]+ 4.0*#[[8]], 4.0*#[[6]]}},0.04]], Text [Style[#[[1]], 5, Black], {#[[4]] - 0.1,#[[5]] + 0.2, 0.1}], Text [Style[Cell == 12 ", "Seed == 10 ", " Subscript[U,0] ==  6.0 ", " Subscript[U,1] \[Backslash] Subscript[U,0] == 0.0 ", " Energy == "VALUE_ETOT" , 8, Blue], {18,15,0.1}],
PointSize[0.01*#[[6]]],
RGBColor[0.0*#[[6]], 0.0*#[[7]], 0.0*#[[8]]],
Point[{#[[4]],#[[5]], 0}]} & /@ data,
Axes -> {True, True, False}, Boxed -> False,
AspectRatio -> {1.0, 1.0, 1.0}, ViewPoint -> {0,0, 100} ]
Export["out1.pdf", Graphics[Inset[g1, Automatic, Automatic, Scaled[1]]], ImageResolution->400]
