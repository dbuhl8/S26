unset multiplot
reset

ub = 0.7
idx = 0
pwr = 1

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
    set output "plot_eta.png"
} else if (eps_output) {
    set terminal postscript enh col
    set output "plot_ReB_eta.eps"
}
set tics font "Roman,22"
set title font "Roman,30"
set key font "Roman,20"
set xlabel font "Roman,25"
set ylabel font "Roman,25"

if (multiplot_mode) {
set multiplot layout 2,1 columnsfirst 
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
# idx 2: Re300_Pe30
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
# {{{ Length scale settings
LF_Steady = 1
LF_Stoch = 2**0.5
#}}}

# -------------------------------------------------------------

perform_block_1 = 1

# {{{ First Plot

if (perform_block_1) {

set xlabel "Re_{G}" font "Roman,20"
set key top right
set ylabel "n"
#set title "Effective  Re_G v n"
set format x "10^{%T}"
set log x
set yrange [0:1]

set arrow from 1,0 to 1,1 nohead lw 3 dt 2 lc rgb "black"
set label "Re_G = 1" at 1,0.8 offset 1,0 font "Roman,25"

plot "steady_tavg_eta.dat" \
   i 0 u ($16/$2):20:21 w yerrorbars pt 5 ps 2 lc rgb "dark-violet" title 'Steady Eta',\
"" i 0 u ($16/$2):70:71 w yerrorbars pt 5 ps 1.6 lc rgb "red" title 'Steady Eta_{Lam}',\
"" i 0 u ($16/$2):80:81 w yerrorbars pt 5 ps 1.6 lc rgb "blue" title 'Steady Eta_{Turb}',\
"" i 2 u ($16/$2):20:21 w yerrorbars pt 4 ps 2 lc rgb "dark-violet" notitle,\
"" i 2 u ($16/$2):70:71 w yerrorbars pt 4 ps 1.6 lc rgb "red" notitle,\
"" i 2 u ($16/$2):80:81 w yerrorbars pt 4 ps 1.6 lc rgb "blue" notitle,\
"stoch_tavg_eta.dat" \
   i 0 u ($16/$2):20:21 w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Steady Eta',\
"" i 0 u ($16/$2):70:71 w yerrorbars pt 9 ps 1.6 lc rgb "red" title 'Steady Eta_{Lam}',\
"" i 0 u ($16/$2):80:81 w yerrorbars pt 9 ps 1.6 lc rgb "blue" title 'Steady Eta_{Turb}',\
"" i 1 u ($16/$2):20:21 w yerrorbars pt 8 ps 2 lc rgb "dark-violet" notitle ,\
"" i 1 u ($16/$2):70:71 w yerrorbars pt 8 ps 1.6 lc rgb "red" notitle,\
"" i 1 u ($16/$2):80:81 w yerrorbars pt 8 ps 1.6 lc rgb "blue" notitle,\


} # end of block 1 }}}

# -------------------------------------------------------------

perform_block_2 = 0

# {{{ Second Plot

if (perform_block_2) {

set xlabel "B_{eff}"
set key bottom right
set ylabel "{/Symbol h}"
set log xy

plot "steady_tavg_eta.dat" \
   i 0 u ($2/($8**2)):20:21 w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady (Re, Pr) = (600, 0.1)',\
"" i 1 u ($2/($8**2)):20:21 w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady (Re, Pr) = (600, 0.05)',\
"" i 2 u ($2/($8**2)):20:21 w yerrorbars pt 4 ps 2 lc rgb "red" title 'Steady (Re, Pr) = (1000, 0.1)',\
"" i 3 u ($2/($8**2)):20:21 w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady (Re, Pr) = (300, 0.1)',\
"" i 4 u ($2/($8**2)):20:21 w yerrorbars pt 4 ps 2 lc rgb "green" title 'Steady  (Re, Pr) = (1000, 0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u ($2/($8**2)):20:21 w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. (Re, Pr) = (600, 0.1)',\
"" i 1 u ($2/($8**2)):20:21 w yerrorbars pt 9 ps 2 lc rgb "red" title 'Stoch. (Re, Pr) = (1000, 0.1)',\
"" i 2 u ($2/($8**2)):20:21 w yerrorbars pt 9 ps 2 lc rgb "green" title 'Stoch. (Re, Pr) = (300, 0.1)',\

} # end of block 2 }}} 

# -------------------------------------------------------------

perform_block_3 = 0

# {{{ Third Plot

if (perform_block_3) {

set xlabel "B_{eff}"
set key top right
set ylabel "n"
#set title "Effective Re_{B} v n_{Lam}"
set log xy
##set yrange [0:ub]

plot "steady_tavg_eta.dat" \
   i 0 u ($2/($8**2)):70:71 w yerrorbars pt 4 ps 2 lc rgb "red" title 'Steady Eta_{Lam}',\
"" i 2 u ($2/($8**2)):70:71 w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady Eta_{Lam}',\
"stoch_tavg_eta.dat" \
   i 0 u ($2/($8**2)):70:71 w yerrorbars pt 9 ps 2 lc rgb "red" title 'Stoch. Eta_{Lam}',\
"" i 1 u ($2/($8**2)):70:71 w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. Eta_{Lam}',\

} # end of block 3 }}} 

# -------------------------------------------------------------

perform_block_4 = 0

# {{{ B_eff v eta_turb

if (perform_block_4) {

set xlabel "B_{eff}"
set key top right
set ylabel "{/Symbol h}"
set log xy
#set yrange [0:ub]

plot "steady_tavg_eta.dat" \
   i 0 u ($2/($8**2)):80:81 w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady Eta_{Turb}',\
"" i 2 u ($2/($8**2)):80:81 w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady Eta_{Turb}',\
"stoch_tavg_eta.dat" \
   i 0 u ($2/($8**2)):80:81 w yerrorbars pt 9 ps 2 lc rgb "blue" title 'Stoch. Eta_{Turb}',\
"" i 1 u ($2/($8**2)):80:81 w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. Eta_{Turb}',\

} # end of block 3 }}} 

# -------------------------------------------------------------

perform_block_5 = 0

# {{{ VTurb v Eta

if (perform_block_5) {

set xlabel "V_{Turb, eff}"
set key bottom center
set ylabel "n"
set title "V_{Turb, eff} v n"
set xrange [0.5*1e-2:1]
set log x

plot "steady_tavg_eta.dat" \
   i 0 u 92:20:21 w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady Re=600, Pr=0.1',\
"" i 2 u 92:20:21 w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady Re=1000, Pr=0.1',\
"" i 1 u 92:20:21 w yerrorbars pt 4 ps 2 lc rgb "forest-green" title 'Steady Re=600, Pr=0.05',\
"" i 3 u 92:20:21 w yerrorbars pt 4 ps 2 lc rgb "navy-blue" title 'Steady Re=300, Pr=0.1',\
"" i 4 u 92:20:21 w yerrorbars pt 4 ps 2 lc rgb "dark-red" title 'Steady Re=1000, Pr=0.01',\
"stoch_tavg_eta.dat" \
   i 0 u 92:20:21 w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Stoch. Re=600',\
"" i 1 u 92:20:21 w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. Re=1000',\
[1e-2:0.7] 0.06*log(x) + 0.46 w l dt 2 lc rgb "red" title "0.06 log(V_{Turb}) + 0.46",\
[1e-2:0.7] 0.06*log(x) + 0.58 w l dt 2 lc rgb "blue" title "0.06 log(V_{Turb}) + 0.58"

} # end of block 5 }}} 

# -------------------------------------------------------------

perform_block_6 = 0

# {{{ Gibson v eta

if (perform_block_6) {

set xlabel "Re_{G}"
set ylabel "{/Symbol h}"
set log xy
set key bottom left offset -15,0

set arrow from 1,0.03 to 1,1 nohead dt 2 lw 3 lc rgb "black"
set label "Re_G = 1" at 1,0.05 offset 1,0 font "Roman, 25"

plot "steady_tavg_eta.dat" \
   i 0 u ($16/$2):20:21 w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady  (Re, Pr) = (600, 0.1)',\
"" i 1 u ($16/$2):20:21 w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady  (Re, Pr) = (600, 0.05)',\
"" i 2 u ($16/$2):20:21 w yerrorbars pt 4 ps 2 lc rgb "red" title 'Steady  (Re, Pr) = (1000, 0.1)',\
"" i 3 u ($16/$2):20:21 w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady  (Re, Pr) = (300, 0.1)',\
"" i 4 u ($16/$2):20:21 w yerrorbars pt 4 ps 2 lc rgb "green" title 'Steady  (Re, Pr) = (1000, 0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u ($16/$2):20:21 w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. (Re, Pr) = (600, 0.1)',\
"" i 1 u ($16/$2):20:21 w yerrorbars pt 9 ps 2 lc rgb "red" title 'Stoch. (Re, Pr) = (1000, 0.1)',\
"" i 2 u ($16/$2):20:21 w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Stoch (Re, Pr) = (300, 0.1)'


} # end of block 6 }}} 

# -------------------------------------------------------------

perform_block_7 = 0

# {{{ Gibson v Weird Eta

if (perform_block_7) {

set xlabel "Re_{B, emergent}"
set key top right
set ylabel "n_{emergent}"
set title "Re_{B, emergent} v n_{emergent}"
set log x

plot "steady_tavg_eta.dat" \
   i 0 u ($16/$2):20:21 w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady Re=600, Pr=0.1',\
"" i 2 u ($16/$2):20:21 w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady Re=1000, Pr=0.1',\
"" i 1 u ($16/$2):20:21 w yerrorbars pt 4 ps 2 lc rgb "forest-green" title 'Steady Re=600, Pr=0.05',\
"" i 3 u ($16/$2):20:21 w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady Re=300, Pr=0.1',\
"" i 4 u ($16/$2):20:21 w yerrorbars pt 4 ps 2 lc rgb "dark-red" title 'Steady Re=1000, Pr=0.01',\
"stoch_tavg_eta.dat" \
   i 0 u ($16/$2):20:21 w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Stoch. Re=600',\
"" i 1 u ($16/$2):20:21 w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. Re=1000',\

} # end of block 6 }}} 

# -------------------------------------------------------------

perform_block_8 = 0

# {{{ 8: Inertial Range v Eta

if (perform_block_8) {

set xlabel "Inertial Range"
set key bottom right
set ylabel "{/Symbol h}"
set log xy

plot "steady_tavg_eta.dat" \
   i 0 u ($2/($8**2)):20:21 w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady (Re, Pr) = (600, 0.1)',\
"" i 1 u ($2/($8**2)):20:21 w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady (Re, Pr) = (600, 0.05)',\
"" i 2 u ($2/($8**2)):20:21 w yerrorbars pt 4 ps 2 lc rgb "red" title 'Steady (Re, Pr) = (1000, 0.1)',\
"" i 3 u ($2/($8**2)):20:21 w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady (Re, Pr) = (300, 0.1)',\
"" i 4 u ($2/($8**2)):20:21 w yerrorbars pt 4 ps 2 lc rgb "green" title 'Steady  (Re, Pr) = (1000, 0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u ($2/($8**2)):20:21 w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. (Re, Pr) = (600, 0.1)',\
"" i 1 u ($2/($8**2)):20:21 w yerrorbars pt 9 ps 2 lc rgb "red" title 'Stoch. (Re, Pr) = (1000, 0.1)',\
"" i 2 u ($2/($8**2)):20:21 w yerrorbars pt 9 ps 2 lc rgb "green" title 'Stoch. (Re, Pr) = (300, 0.1)',\

}

# }}} End of block 8

# -------------------------------------------------------------

perform_block_9 = 0

# {{{ l_T/l_K v eta

if (perform_block_9) {

set xlabel "Re_{G}"
set ylabel "(l_T/l_K)^{3/4}"
set log xy
set key top right

#set arrow from 1,0.03 to 1,1 nohead dt 2 lw 3 lc rgb "black"
#set label "l_T/L_K = 1" at 1,0.05 offset 1,0 font "Roman, 25"

plot "steady_tavg_eta.dat" \
   i 0 u ($16/$2):((5*($8**2 + $12**2)/$16)**0.5/(1./($1**2*$16))**0.25)**(3./4) w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady  (Re, Pr) = (600, 0.1)',\
"" i 1 u ($16/$2):((5*($8**2 + $12**2)/$16)**0.5/(1./($1**2*$16))**0.25)**(3./4) w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady  (Re, Pr) = (600, 0.05)',\
"" i 2 u ($16/$2):((5*($8**2 + $12**2)/$16)**0.5/(1./($1**2*$16))**0.25)**(3./4) w yerrorbars pt 4 ps 2 lc rgb "red" title 'Steady  (Re, Pr) = (1000, 0.1)',\
"" i 3 u ($16/$2):((5*($8**2 + $12**2)/$16)**0.5/(1./($1**2*$16))**0.25)**(3./4) w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady  (Re, Pr) = (300, 0.1)',\
"" i 4 u ($16/$2):((5*($8**2 + $12**2)/$16)**0.5/(1./($1**2*$16))**0.25)**(3./4) w yerrorbars pt 4 ps 2 lc rgb "green" title 'Steady  (Re, Pr) = (1000, 0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u ($16/$2):((5*($8**2 + $12**2)/$16)**0.5/(1./($1**2*$16))**0.25)**(3./4) w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. (Re, Pr) = (600, 0.1)',\
"" i 1 u ($16/$2):((5*($8**2 + $12**2)/$16)**0.5/(1./($1**2*$16))**0.25)**(3./4) w yerrorbars pt 9 ps 2 lc rgb "red" title 'Stoch. (Re, Pr) = (1000, 0.1)',\
"" i 2 u ($16/$2):((5*($8**2 + $12**2)/$16)**0.5/(1./($1**2*$16))**0.25)**(3./4) w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Stoch (Re, Pr) = (300, 0.1)'


} # end of block 6 }}} 

# -------------------------------------------------------------

