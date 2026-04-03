
reset
# Plot Settings {{{ 
    set terminal postscript enh col
    set output "plot_length_scale_Re600_Pe60_B100.eps"
    set xlabel font "Roman,25"
    set ylabel font "Roman,25"
    set label font "Roman,25"
    set title font "Roman,30"
    set tics font "Roman,22"
    set key font "Roman,25"
# }}}
# {{{ Length scale settings

Fr = (100.**-0.5)
uh = 2.205
EU = 2.448
Re = 600.
mdisp = 145.820
PI = 3.14159

KF = 0.5
KForce = 2.**(0.5)
LF = 2*3.14159/KF
LZ = LF*Fr
KZ = 2*PI/LZ
LO = LF*Fr*(Fr*mdisp/Re)**0.5
KO = 2*PI/LO
LK = LF*(1./(Re**2*mdisp))**0.25
KK = 2*PI/LK
LT = LZ*uh*(15./mdisp)**0.5
KT = 2*PI/LT
# LH
# }}}

# ---------------------------------------------------------------

perform_block_1 = 1

# {{{ Block 1: h_spec 

if (perform_block_1){

#set format y "10^{%T}"
set xlabel "t"
set log y
#set yrange[1e-14:1]

#set arrow from KForce, 1e-14 to KForce, 1 nohead lw 3 dt 2 lc rgb "blue"
#set label "k_F" at KForce,1e-12 offset -4,0 font "Roman,25"
#set arrow from KZ, 1e-14 to KZ, 1e-3 nohead lw 3 dt 2 lc rgb "forest-green"
#set label "k_Z" at KZ,1e-12 offset -4,0 font "Roman,25"
#set arrow from KO, 1e-14 to KO, 1e-3 nohead lw 3 dt 2 lc rgb "goldenrod"
#set label "k_O" at KO,1e-12 offset -4,0 font "Roman,25"
#set arrow from KT, 1e-14 to KT, 1e-3 nohead lw 3 dt 2 lc rgb "dark-violet"
#set label "k_T" at KT,1e-12 offset -4,0 font "Roman,25"
#set arrow from KK, 1e-14 to KK, 1e-3 nohead lw 3 dt 2 lc rgb "red"
#set label "k_K" at KK,1e-4 offset -4,0 font "Roman,25"

L = 4*3.14159
#LO = L*Fr*((Fr/Re)*$47)**0.5
#LK = L*(1/(Re**2*$47))**0.25
#LT = L*(15.*($36**2+$37**2)/$47)**0.5
#LT = LK*(15.**0.5)*(Re*($36**2+$37**2)**0.5)**0.25
#LZ = L*Fr
#($36**2+$37**2)**0.5, uh

Fr = 0.1
Re = 600.
Pe = 60.


plot "<cat OUT*" u 2:(L*Fr*($36**2+$37**2)**0.5) w lp pt 7 ps 1 lc rgb "black" title "l_z",\
"" u 2:(L*Fr*(($36**2+$37**2)**0.5)*(((Fr/Re)*$47)**0.5)) w lp pt 7 ps 1 lc rgb "blue" title "l_O",\
"" u 2:(L*(1/(Re**2*($36**2+$37**2)*$47))**0.25)*(15.**0.5)*(Re*($36**2+$37**2)**0.5)**0.25 \
w lp pt 7 ps 1 lc rgb "red" title "l_T",\
"" u 2:(L*(1/(Re**2*($36**2+$37**2)*$47))**0.25) w lp pt 7 ps 1 lc rgb "forest-green" title "l_K",\

#LO,\
#LZ,\
#LT,\
#LK


} # End of Block 1 }}}

# ---------------------------------------------------------------


