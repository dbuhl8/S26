unset multiplot
reset

# {{{ Plot Settings

png_output = 0
eps_output = 1
multiplot_mode = 0

if (png_output){
    if (multiplot_mode){
        set terminal png size 900,1600
    } else {
        set terminal png size 800,600
    }
    set output "plot_vturb.png"
} else if (eps_output) {
    set terminal postscript enh col
    set output "plot_ReG_v_Vturb.eps"
}
set tics font "Roman,22"
set title font "Roman,35"
set key font "Roman,22"
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

perform_block_1 = 0

# {{{ First Plot

if (perform_block_1) {

set xlabel "Re_{b, eff}"
set key top right
set ylabel "Eta" rotate by 0 
set title "Normal"
set log x
set yrange[0:1]

plot "steady_tavg_eta.dat" \
 i idx u (($8**3)*$1/$2):38 pt 5 ps 2 lc rgb "red" title 'Steady VTurb',\
"stoch_tavg_eta.dat" \
 i idx u (($8**3)*$1/$2):38 pt 7 ps 2 lc rgb "red" title 'Stoch. VTurb',\
 i idx u (($8**3)*$1/$2):37 pt 7 ps 2 lc rgb "blue" title 'Stoch. VLam',\
 i idx u (($8**3)*$1/$2):37 pt 5 ps 2 lc rgb "blue" title 'Steady VLam',\

}

# }}} End of block 1

# -------------------------------------------------------------

perform_block_2 = 0

# {{{ Second Plot

if (perform_block_2) {

set xlabel "Re_{b, eff}"
set key top right
set ylabel "Eta" rotate by 0 
set title "Using Effective Fr"
set log x
set yrange[0:1]

plot "steady_tavg_eta.dat" \
 i idx u (($8**3)*$1/$2):40 pt 5 ps 2 lc rgb "red" title 'Steady VTurb',\
"stoch_tavg_eta.dat" \
 i idx u (($8**3)*$1/$2):40 pt 7 ps 2 lc rgb "red" title 'Stoch. VTurb',\
 i idx u (($8**3)*$1/$2):39 pt 7 ps 2 lc rgb "blue" title 'Stoch. VLam',\
 i idx u (($8**3)*$1/$2):39 pt 5 ps 2 lc rgb "blue" title 'Steady VLam',\

} 

# }}} End of block 2

# -------------------------------------------------------------

perform_block_3 = 0

# {{{ Third Plot

if (perform_block_3) {

set xlabel "Re_{B, eff}"
set key top left
set ylabel "V_{Turb, eff}"
set title "Effective Re_{B, eff} v V_{Turb, eff}"
set log x
set yrange[0:1]

plot "steady_tavg_eta.dat" \
   i 0 u (($8**3)*$1/$2):92 pt 4 ps 2 lc rgb "red" title 'Steady Re=600,Pr=0.1',\
"" i 1 u (($8**3)*$1/$2):92 pt 4 ps 2 lc rgb "forest-green" title 'Steady Re=600,Pr=0.05',\
"" i 2 u (($8**3)*$1/$2):92 pt 4 ps 2 lc rgb "blue" title 'Steady Re=1000,Pr=0.1',\
"" i 3 u (($8**3)*$1/$2):92 pt 4 ps 2 lc rgb "dark-violet" title 'Steady Re=300,Pr=0.1',\
"" i 4 u (($8**3)*$1/$2):92 pt 4 ps 2 lc rgb "black" title 'Steady Re=1000,Pr=0.01',\
"stoch_tavg_eta.dat" \
   i 0 u (($8**3)*$1/$2):92 pt 9 ps 2 lc rgb "red" title 'Stoch. Re=600,Pr=0.1',\
"" i 1 u (($8**3)*$1/$2):92 pt 9 ps 2 lc rgb "blue" title 'Stoch. Re=1000,Pr=0.1'

}

# }}}

# -------------------------------------------------------------

perform_block_4 = 0 

# {{{ Fourth Plot

if (perform_block_4) {

set xlabel "Vturb_{old}"
set key top right
set ylabel "VTurb_{new,eff}" rotate by 0 
set yrange[0:1]

plot "steady_tavg_eta.dat" \
   i 0:2:2 u 38:42 pt 4 ps 2 lc rgb "red" title 'Steady VTurb',\
"stoch_tavg_eta.dat" \
   i 0:1 u 38:42 pt 9 ps 2 lc rgb "red" title 'Stoch. VTurb',\
   x

}

# }}} End of block 4

# -------------------------------------------------------------

perform_block_5 = 0

# {{{ Fifth Plot

if (perform_block_5) {

set xlabel "Re_{G}"
set key bottom right offset 3,0
set ylabel "V_{Turb}"
#set format x "10^{%T}"
#set format y "10^{%T}"
set log x
set yrange[0:1]

set arrow from 1,0 to 1,1 nohead dt 2 lw 3 lc rgb "black"
set label "Re_G = 1" at 1,0.8 offset 1,0 font "Roman,25"

plot "steady_tavg_eta.dat" \
   i 0 u ($16/$2):92 pt 4 ps 2 lc rgb "black" title 'Steady (Re, Pr) = (600, 0.1)',\
"" i 1 u ($16/$2):92 pt 4 ps 2 lc rgb "blue" title 'Steady (Re, Pr) = (600, 0.05)',\
"" i 2 u ($16/$2):92 pt 4 ps 2 lc rgb "red" title 'Steady (Re, Pr) = (1000, 0.1)',\
"" i 3 u ($16/$2):92 pt 4 ps 2 lc rgb "dark-violet" title 'Steady (Re, Pr) = (300, 0.1)',\
"" i 4 u ($16/$2):92 pt 4 ps 2 lc rgb "green" title 'Steady (Re, Pr) = (1000, 0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u ($16/$2):92 pt 9 ps 2 lc rgb "black" title 'Stoch. (Re, Pr) = (600, 0.1)',\
"" i 1 u ($16/$2):92 pt 9 ps 2 lc rgb "red" title 'Stoch. (Re, Pr) = (1000, 0.1)',\
"" i 2 u ($16/$2):92 pt 9 ps 2 lc rgb "dark-violet" title 'Stoch. (Re, Pr) = (300, 0.1)'

}

# }}}

# -------------------------------------------------------------

perform_block_6 = 1

# {{{ Sixth Plot

if (perform_block_6) {

set xlabel "Re_{B, eff}"
set key bottom right offset 3,0
set ylabel "V_{Turb}"
set log x
set yrange[0:1]

unset arrow
unset label

plot "steady_tavg_eta.dat" \
   i 0 u (($8**3)*$1/$2):92 pt 4 ps 2 lc rgb "black" title 'Steady (Re, Pr) = (600, 0.1)',\
"" i 1 u (($8**3)*$1/$2):92 pt 4 ps 2 lc rgb "blue" title 'Steady (Re, Pr) = (600, 0.05)',\
"" i 2 u (($8**3)*$1/$2):92 pt 4 ps 2 lc rgb "red" title 'Steady (Re, Pr) = (1000, 0.1)',\
"" i 3 u (($8**3)*$1/$2):92 pt 4 ps 2 lc rgb "dark-violet" title 'Steady (Re, Pr) = (300, 0.1)',\
"" i 4 u (($8**3)*$1/$2):92 pt 4 ps 2 lc rgb "green" title 'Steady (Re, Pr) = (1000, 0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u (($8**3)*$1/$2):92 pt 9 ps 2 lc rgb "black" title 'Stoch. (Re, Pr) = (600, 0.1)',\
"" i 1 u (($8**3)*$1/$2):92 pt 9 ps 2 lc rgb "red" title 'Stoch. (Re, Pr) = (1000, 0.1)',\
"" i 2 u (($8**3)*$1/$2):92 pt 9 ps 2 lc rgb "dark-violet" title 'Stoch. (Re, Pr) = (300, 0.1)'

}

# }}}

# -------------------------------------------------------------

