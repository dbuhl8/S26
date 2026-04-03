reset


# {{{ Plot Settings 
png_output = 0
eps_output = 1
multiplot_mode = 0

if (png_output){
    set terminal png size 800,600
    set output "plot_mdisp.png"
} else if (eps_output) {
    set terminal postscript enh col
    set output "length_scale_ratios.eps"
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
# steady_tavg_eta.dat: 
# idx 0: Re600_Pe60
# idx 1: Re600_Pe30
# idx 2: Re1000_Pe100
# idx 3: Re300_Pe30
# idx 4: Re1000_Pe10

# stoch_tavg_eta.dat:
# idx 0: Re600_Pe60
# idx 1: Re1000_Pe100
# }}}
# {{{ Columns 
# in each file:
# Re(1)
# B(2) 
# Pr(3)
# Pe(4)
# BPe(5)
# lb(6)
# ub(7)
# uh_rms(8/9)
# vortz_rms(10/11)   
# wrms(12/13)
# tdisp(14/15)
# mdisp(16/17)
# eta (local)(18/19)
# eta (global)(20/21)   
# lam_wrms(22/23)
# lam tdisp(24/25)
# lam mdisp(26/27)
# lam eta (local)(28/29)
# lam eta (global)(30/31)    
# turb_wrms (32/33)
# turb tdisp (34/35)
# turb mdisp (36/37)
# turb eta (local)(38/39)
# turb eta (global)(40/41)
# lam_Fr_wrms(42/43) 
# lam_Fr tdisp(44/45)
# lam_Fr mdisp (46/47)
# lam_Fr eta (local) (48/49)
# lam_Fr eta (global)(50/51)
# turb_Fr_wrms(52/53) 
# turb_Fr tdisp (54/55)
# turb_Fr mdisp (56/57)
# turb_Fr eta (local) (58/59)
# turb_Fr eta (global)(60/61)
# lam_Fr_wmrs(62/63)
# lam_Fr_vortz tdisp(64/65)
# lam_Fr_vortz mdisp(66/67)
# lam_Fr_vortz eta (local)(68/69)
# lam_Fr_vortz eta (global)(70/71)
# turb_Fr_wrms(72/73)
# turb_Fr_vortz tdisp (74/75)
# turb_Fr_vortz mdisp(76/77)
# turb_Fr_vortz eta(local)(78/79)
# turb_Fr_vortz eta(glbal)(80/81) 
# vlam(82/83)
# vturb(84/85)
# vlam_Fr(86/87)
# vturb_Fr(88/89)
# vlam_Fr_vortz (90/91)
# vturb_Fr_vortz (92/93)
# wlam_wght (94/95) 
# wturb_wght (96/97)
# wlam_eff_wght(98/99)
# wturb_eff_wght (100/101) 
# }}}

# -------------------------------------------------------------

perform_block_1 = 1

# {{{ First Plot

if (perform_block_1) {

set xlabel "l_z/l_T"
set key top right
set ylabel "l_O/l_K"
set title "Correlation Between Length Scales"
set log xy

LO = 4*3.14159

plot "steady_tavg_eta.dat" \
   i 0 u ((15**0.5)*$8/($16**0.5)):(($2**(-0.75))*($16**.75)) w yerrorbars pt 5 ps 2 lc rgb "black" title 'Steady (Re,Pr) = (600,0.1)',\
"" i 1 u ((15**0.5)*$8/($16**0.5)):(($2**(-0.75))*($16**.75)) w yerrorbars pt 5 ps 2 lc rgb "blue" title 'Steady (Re,Pr) = (600,0.05)',\
"" i 2 u ((15**0.5)*$8/($16**0.5)):(($2**(-0.75))*($16**.75)) w yerrorbars pt 5 ps 2 lc rgb "red" title 'Steady (Re,Pr) = (1000,0.1)',\
"" i 3 u ((15**0.5)*$8/($16**0.5)):(($2**(-0.75))*($16**.75)) w yerrorbars pt 5 ps 2 lc rgb "dark-violet" title 'Steady (Re,Pr) = (300,0.1)',\
"" i 4 u ((15**0.5)*$8/($16**0.5)):(($2**(-0.75))*($16**.75)) w yerrorbars pt 5 ps 2 lc rgb "green" title 'Steady (Re,Pr) = (1000,0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u ((15**0.5)*$8/($16**0.5)):(($2**(-0.75))*($16**.75)) w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. (Re,Pr) = (600,0.1)',\
"" i 1 u ((15**0.5)*$8/($16**0.5)):(($2**(-0.75))*($16**.75)) w yerrorbars pt 9 ps 2 lc rgb "red" title 'Stoch. (Re,Pr) = (1000,0.1)',\
"" i 2 u ((15**0.5)*$8/($16**0.5)):(($2**(-0.75))*($16**.75)) w yerrorbars pt 9 ps 2 lc rgb "green" title 'Stoch. (Re,Pr) = (300,0.1)',\

} # end of block 1 }}}

# -------------------------------------------------------------

perform_block_2 = 0

# {{{ Second Plot

if (perform_block_2) {

set xlabel "Re_{B, eff}"
set key top right
set ylabel "MDisp_{Lam}"
set title "Effective MDisp v Re_{B}"
set log x


plot "steady_tavg_eta.dat" \
   i idx u (($8**3)*$1/$2):($66/($1*$8**pwr)):($67/($1*$8**pwr)) w yerrorbars pt 4 ps 2 lc rgb "red" title 'Steady MDisp_{Lam}',\
"stoch_tavg_eta.dat" \
   i idx u (($8**3)*$1/$2):($66/($1*$8**pwr)):($67/($1*$8**pwr)) w yerrorbars pt 9 ps 2 lc rgb "red" title 'Stoch. MDisp_{Lam}',\

} # end of block 2 }}}

# -------------------------------------------------------------

perform_block_3 = 0

# {{{ Third Plot

if (perform_block_3) {


set xlabel "Re_{B, eff}"
set key top right
set ylabel "MDisp_{Turb}"
set title "Effective MDisp_{Turb} v Re_{B}"
set log x


plot "steady_tavg_eta.dat" \
   i idx u (($8**3)*$1/$2):($76/($1*$8**pwr)):($77/($1*$8**pwr)) w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady MDisp_{Turb}',\
"stoch_tavg_eta.dat" \
   i idx u (($8**3)*$1/$2):($76/($1*$8**pwr)):($77/($1*$8**pwr)) w yerrorbars pt 9 ps 2 lc rgb "blue" title 'Stoch. MDisp_{Turb}',\

} # end of block 3 }}}

# -------------------------------------------------------------

perform_block_4 = 0

# {{{ Third Plot

if (perform_block_4) {


set xlabel "Re_{B, eff}"
set key top right
set ylabel "MDisp_{Turb}"
set title "Effective MDisp_{Turb} v Re_{B}"
set log xy

plot "steady_tavg_eta.dat" \
   i idx u (($8**3)*$1/$2):($76/($1*$8**pwr)):($77/($1*$8**pwr)) w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady MDisp_{Turb}',\
"" i idx u (($8**3)*$1/$2):($66/($1*$8**pwr)):($67/($1*$8**pwr)) w yerrorbars pt 4 ps 2 lc rgb "red" title 'Steady MDisp_{Lam}',\
"" i idx u (($8**3)*$1/$2):($16/($1*$8**pwr)):($17/($1*$8**pwr)) w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady MDisp',\
"stoch_tavg_eta.dat" \
   i idx u (($8**3)*$1/$2):($76/($1*$8**pwr)):($77/($1*$8**pwr)) w yerrorbars pt 9 ps 2 lc rgb "blue" title 'Stoch. MDisp_{Turb}',\
"" i idx u (($8**3)*$1/$2):($66/($1*$8**pwr)):($67/($1*$8**pwr)) w yerrorbars pt 9 ps 2 lc rgb "red" title 'Steady MDisp_{Lam}',\
"" i idx u (($8**3)*$1/$2):($16/($1*$8**pwr)):($17/($1*$8**pwr)) w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Steady MDisp',\

} # end of block 4 }}}

# -------------------------------------------------------------

perform_block_5 = 0

# {{{ Fifth Plot

if (perform_block_5) {


set xlabel "Re_{B, eff}"
set key bottom right
set ylabel "Re_{G}"
#set title "Effective Re_{B} v Re_{G}"
set format x "10^{%T}"
set format y "10^{%T}"
set log xy

plot "steady_tavg_eta.dat" \
   i 0 u (($8**3)*$1/$2):($16/$2):($17/$2) w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady (Re,Pr) = (600,0.1)',\
"" i 1 u (($8**3)*$1/$2):($16/$2):($17/$2) w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady (Re,Pr) = (600,0.05)',\
"" i 2 u (($8**3)*$1/$2):($16/$2):($17/$2) w yerrorbars pt 4 ps 2 lc rgb "red" title 'Steady (Re,Pr) = (1000,0.1)',\
"" i 3 u (($8**3)*$1/$2):($16/$2):($17/$2) w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady (Re,Pr) = (300,0.1)',\
"" i 4 u (($8**3)*$1/$2):($16/$2):($17/$2) w yerrorbars pt 4 ps 2 lc rgb "green" title 'Steady (Re,Pr) = (1000,0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u (($8**3)*$1/$2):($16/$2):($17/$2) w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. (Re,Pr) = (600,0.1)',\
"" i 1 u (($8**3)*$1/$2):($16/$2):($17/$2) w yerrorbars pt 9 ps 2 lc rgb "red" title 'Stoch. (Re,Pr) = (1000,0.1)',\
0.03*x title "0.03x"


} # end of block 5 }}}

# -------------------------------------------------------------

perform_block_6 = 0

# {{{ Sixth Plot

if (perform_block_6) {


set xlabel "B"
set key bottom right
set ylabel "Mdisp"
#set title "Effective Re_{B} v Re_{G}"
#set format x "10^{%T}"
#set format y "10^{%T}"
set log x
set yrange [61.7:61.8]

plot "steady_tavg_eta.dat" \
   i 0 u 2:16:17 w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady (Re,Pr) = (600,0.1)',\
"stoch_tavg_eta.dat" \
   i 0 u 2:16:17 w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. (Re,Pr) = (600,0.1)',\
61.7521


} # end of block 5 }}}

# -------------------------------------------------------------



