reset

idx = 0

# {{{ Plot Settings 
png_output = 0
eps_output = 1
multiplot_mode = 0

if (png_output){
    set terminal png size 800,600
    set output "plot_mdisp.png"
} else if (eps_output) {
    set terminal postscript enh col
    set output "ReG_VTurb.eps"
}
set tics font "Roman,22"
set title font "Roman,35"
set key font "Roman,20"
set xlabel font "Roman,25"
set ylabel font "Roman,25"

if (multiplot_mode) {
set multiplot layout 3,1 columnsfirst 
}

# }}}
# {{{file key: 
#
# idx 0: Re300_Pe30, steady
# idx 1: Re600_Pe30, steady
# idx 2: Re600_Pe60, steady
# idx 3: Re1000_Pe10, steady
# idx 4: Re1000_Pe100, steady
# idx 5: Re600_Pe600, steady
# idx 6: Re600_Pe60, stoch
# idx 7: Re1000_Pe100, stoch

# }}}

# -------------------------------------------------------------

perform_block_1 = 0

# {{{ First Plot

if (perform_block_1) {


set xlabel "Re_{B}^*"
set key bottom left
set ylabel "VLam"
#set title "EffecRe_B v VTurb"
set log xy
set yrange [0.1:1]

plot "tavg_bflux.dat" \
   i 0 u (($30**3)*$1/$3):(1-$20):($31*3*($30**2)*$1/$3):21 w xyerrorbars pt 7 ps 2 lc rgb "blue" title 'Steady (300,30)',\
"" i 1 u (($30**3)*$1/$3):(1-$20):($31*3*($30**2)*$1/$3):21 w xyerrorbars pt 7 ps 2 lc rgb "red" title 'Steady (600,30)',\
"" i 2 u (($30**3)*$1/$3):(1-$20):($31*3*($30**2)*$1/$3):21 w xyerrorbars pt 7 ps 2 lc rgb "black" title 'Steady (600,60)',\
"" i 3 u (($30**3)*$1/$3):(1-$20):($31*3*($30**2)*$1/$3):21 w xyerrorbars pt 7 ps 2 lc rgb "purple" title 'Steady (1000,10)',\
"" i 4 u (($30**3)*$1/$3):(1-$20):($31*3*($30**2)*$1/$3):21 w xyerrorbars pt 7 ps 2 lc rgb "forest-green" title 'Steady (1000,100)',\
"" i 5 u (($30**3)*$1/$3):(1-$20):($31*3*($30**2)*$1/$3):21 w xyerrorbars pt 7 ps 2 lc rgb "#CC79A7" title 'Steady (600,600)',\
"" i 6 u (($30**3)*$1/$3):(1-$20):($31*3*($30**2)*$1/$3):21 w xyerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch (600,60)',\
"" i 7 u (($30**3)*$1/$3):(1-$20):($31*3*($30**2)*$1/$3):21 w xyerrorbars pt 9 ps 2 lc rgb "forest-green" title 'Stoch (1000,100)',\
0.8*(x/100.)**(-0.25) w l lw 2 dt 2 title "Re_G^(-0.25)"
#[0.01:1] 0.05*(x/100)**4 w l lw 2 dt 4 title "Re_G^4"


} # end of block 1 }}}

# -------------------------------------------------------------

perform_block_2 = 1

# {{{ Second Plot

if (perform_block_2) {

set xlabel "Re_{G}"
set key bottom right
set ylabel "VTurb"
set log xy
set yrange [0.001:1]

plot "tavg_bflux.dat" \
   i 0 u ($18/$3):20:($19/$3):21 w xyerrorbars pt 7 ps 2 lc rgb "blue" title 'Steady (300,30)',\
"" i 1 u ($18/$3):20:($19/$3):21 w xyerrorbars pt 7 ps 2 lc rgb "red" title 'Steady (600,30)',\
"" i 2 u ($18/$3):20:($19/$3):21 w xyerrorbars pt 7 ps 2 lc rgb "black" title 'Steady (600,60)',\
"" i 3 u ($18/$3):20:($19/$3):21 w xyerrorbars pt 7 ps 2 lc rgb "purple" title 'Steady (1000,10)',\
"" i 4 u ($18/$3):20:($19/$3):21 w xyerrorbars pt 7 ps 2 lc rgb "forest-green" title 'Steady (1000,100)',\
"" i 5 u ($18/$3):20:($19/$3):21 w xyerrorbars pt 7 ps 2 lc rgb "#CC79A7" title 'Steady (600,600)',\
"" i 6 u ($18/$3):20:($19/$3):21 w xyerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch (600,60)',\
"" i 7 u ($18/$3):20:($19/$3):21 w xyerrorbars pt 9 ps 2 lc rgb "forest-green" title 'Stoch (1000,100)',\
(x**4.)/(100 + x**4 + x**(15./4)) w l lw 2 dt 2 lc rgb "black"
#1 - x**(-0.25) w l lw 2 dt 2 title "1 - Re_G^(-0.25)",\
#[0.01:1] 0.05*x**4 w l lw 2 dt 4 title "Re_G^4"



} # end of block 2 }}}

# -------------------------------------------------------------

perform_block_3 = 0

# {{{ Third Plot

if (perform_block_3) {


set xlabel "Re_{B, eff}"
set key top left
set ylabel "VTurb"
set title "Effective Re_B v VTurb"
set log x
set yrange [0:1]

plot "tavg_bflux_corrected.dat" \
   i 0 u ((($6**2+$8**2)**(1.5))*$1/$3):20:21 w yerrorbars pt 7 ps 2 lc rgb "blue" title '(300,30)',\
"" i 1 u ((($6**2+$8**2)**(1.5))*$1/$3):20:21 w yerrorbars pt 7 ps 2 lc rgb "red" title '(600,30)',\
"" i 2 u ((($6**2+$8**2)**(1.5))*$1/$3):20:21 w yerrorbars pt 7 ps 2 lc rgb "black" title '(600,60)',\
"" i 3 u ((($6**2+$8**2)**(1.5))*$1/$3):20:21 w yerrorbars pt 7 ps 2 lc rgb "purple" title '(1000,10)',\
"" i 4 u ((($6**2+$8**2)**(1.5))*$1/$3):20:21 w yerrorbars pt 7 ps 2 lc rgb "forest-green" title '(1000,100)',\
"" i 5 u ((($6**2+$8**2)**(1.5))*$1/$3):20:21 w yerrorbars pt 7 ps 2 lc rgb "black" title '(600,60)',\
"" i 6 u ((($6**2+$8**2)**(1.5))*$1/$3):20:21 w yerrorbars pt 7 ps 2 lc rgb "forest-green" title '(1000,100)',\


} # end of block 3 }}}

# -------------------------------------------------------------



