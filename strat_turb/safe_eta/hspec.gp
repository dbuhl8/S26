
reset
# Plot Settings {{{ 
    set terminal postscript enh col
    set output "plot_hspec.eps"
    set xlabel font "Roman,25"
    set ylabel font "Roman,25"
    set label font "Roman,25"
    set title font "Roman,30"
    set tics font "Roman,22"
    set key font "Roman,25"
# }}}
# {{{ Length scale settings

Fr = (1.**-0.5)
uh = 2.040
EU = 2.289
Re = 600.
mdisp = 211.810 
PI = 3.14159

KF = 2.**0.5
LF = 2*3.14149/KF
LO = LF*Fr*(Fr*mdisp/Re)**0.5
KO = 2*PI/LO
LK = LF*(1./(Re**2*mdisp))**0.25
KK = 2*PI/LK
LT = LO*(10*EU/mdisp)**0.5
KT = 2*PI/LT
KT = 10
# LH
# }}}

# ---------------------------------------------------------------

perform_block_1 = 1

# {{{ Block 1: h_spec 

if (perform_block_1){

set format y "10^{%T}"
set xlabel "k_h"
set ylabel "E" rotate by 0 
set log xy
set yrange[1e-14:1]
set arrow from KF, 1e-14 to KF, 1 nohead lw 3 dt 2 lc rgb "blue"
set label "k_F" at KF,1e-12 offset -4,0 font "Roman,25"
set arrow from KO, 1e-14 to KO, 1 nohead lw 3 dt 2 lc rgb "forest-green"
set label "k_O" at KO,1e-12 offset -4,0 font "Roman,25"
set arrow from KT, 1e-14 to KT, 1 nohead lw 3 dt 2 lc rgb "goldenrod"
set label "k_T" at KT,1e-12 offset -4,0 font "Roman,25"
set arrow from KK, 1e-14 to KK, 1 nohead lw 3 dt 2 lc rgb "red"
set label "k_K" at KK,1e-12 offset -4,0 font "Roman,25"

plot "XY_SPEC9" i 25 u (($1**2 + $2**2)**0.5):(($3 + $4)/(uh**2)) pt 7 ps 1 lc rgb "black" title "E(U_h)",\
"" i 25 u (($1**2 + $2**2)**0.5):(($5)/(uh**2)) pt 7 ps 1 lc rgb "dark-turquoise" title "E(w)",\
"" i 25 u (($1**2 + $2**2)**0.5):($11) pt 7 ps 1 lc rgb "sienna1" title "E(b)",\
[KO:KT] 0.1*x**(-3) lw 2 lc rgb "blue" notitle,\
[KF:KO] 0.2*x**(-3) lw 2 lc rgb "red" notitle

} # End of Block 1 }}}

# ---------------------------------------------------------------


