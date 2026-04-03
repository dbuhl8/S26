unset multiplot
reset

# {{{ Plot Settings

png_output = 1
eps_output = 0
multiplot_mode = 0

if (png_output){
    set terminal png size 800,600
    set output "wrms_wz_sep_all.png"
} else if (eps_output) {
    set terminal postscript enh col
    set output "wrms_wz_sep_all.eps"
}
set tics font "Roman,22"
set title font "Roman,35"
set key font "Roman,20"
set xlabel font "Roman,30"
set ylabel font "Roman,30"

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
# lam_Fr_vortz wmrs(62/63)
# lam_Fr_vortz tdisp(64/65)
# lam_Fr_vortz mdisp(66/67)
# lam_Fr_vortz eta (local)(68/69)
# lam_Fr_vortz eta (global)(70/71)
# turb_Fr_vortz wrms(72/73)
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

# Plot variables
# 1: true, 2: false
plot_Re1000 = 1

# -------------------------------------------------------------

perform_block_1 = 0

# {{{ First Plot - set perform_block_1 to 1 to run this block

if (perform_block_1) {

set xlabel "Fr^{-1}" font "Roman,20"
set key bottom left
set ylabel "w_{rms}"
set title "Fr^{-1} v w_{rms}"
set log xy
#set yrange [1e-3:2]

plot "steady_tavg_eta.dat" \
   i 0:2:2 u ($2**0.5):12:13   w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady: Wrms',\
"" i 0:2:2 u ($2**0.5):94:95   w yerrorbars pt 4 ps 2 lc rgb "red" title '       Wlam',\
"" i 0:2:2 u ($2**0.5):96:97   w yerrorbars pt 4 ps 2 lc rgb "blue" title '       Wturb',\
"stoch_tavg_eta.dat" \
   i 0:1 u ($2**0.5):12:13   w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Stoch.: Wrms',\
"" i 0:1 u ($2**0.5):94:95   w yerrorbars pt 9 ps 2 lc rgb "red" title '       Wlam',\
"" i 0:1 u ($2**0.5):96:97   w yerrorbars pt 9 ps 2 lc rgb "blue" title '       Wturb',\

} # end block 1 }}}

# -------------------------------------------------------------

perform_block_2 = 1

# {{{ Second Plot - set perform_block_2 to 1 to run this block


if (perform_block_2){

set xlabel "Fr_{eff}^{-1}"
set key bottom left
set ylabel "w_{rms}/u_{h,rms}"
set title "Vortz Weighting: Effective Fr_{eff}^{-1} v w_{rms, eff}"
set log xy
#set xrange [0.1:10]

if (plot_Re1000) {

plot "steady_tavg_eta.dat" \
   i 0 u ($2**0.5/$8):($12/$8):($13/$8)   w yerrorbars pt 5 ps 2 lc rgb "dark-violet" title 'Steady: Wrms',\
"" i 0 u ($2**0.5/$8):($94/$8):($95/$8)   w yerrorbars pt 5 ps 2 lc rgb "red" title '       Wlam',\
"" i 0 u ($2**0.5/$8):($96/$8):($97/$8)   w yerrorbars pt 5 ps 2 lc rgb "blue" title '       Wturb',\
"" i 2 u ($2**0.5/$8):($12/$8):($13/$8)   w yerrorbars pt 4 ps 2 lc rgb "dark-violet" notitle,\
"" i 2 u ($2**0.5/$8):($94/$8):($95/$8)   w yerrorbars pt 4 ps 2 lc rgb "red" notitle,\
"" i 2 u ($2**0.5/$8):($96/$8):($97/$8)   w yerrorbars pt 4 ps 2 lc rgb "blue" notitle,\
"stoch_tavg_eta.dat" \
   i 0 u ($2**0.5/$8):($12/$8):($13/$8)   w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Stoch.: Wrms',\
"" i 0 u ($2**0.5/$8):($94/$8):($95/$8)   w yerrorbars pt 9 ps 2 lc rgb "red" title '       Wlam',\
"" i 0 u ($2**0.5/$8):($96/$8):($97/$8)   w yerrorbars pt 9 ps 2 lc rgb "blue" title '       Wturb',\
"" i 1 u ($2**0.5/$8):($12/$8):($13/$8)   w yerrorbars pt 8 ps 2 lc rgb "dark-violet" notitle,\
"" i 1 u ($2**0.5/$8):($94/$8):($95/$8)   w yerrorbars pt 8 ps 2 lc rgb "red" notitle,\
"" i 1 u ($2**0.5/$8):($96/$8):($97/$8)   w yerrorbars pt 8 ps 2 lc rgb "black" notitle,\
[1:8] 0.4*x**(-0.5) w l dt 2 lw 2 lc rgb "blue" title "Fr^{-1/2}",\
[1:8] 0.4*x**(-1) w l dt 2 lw 2 lc rgb "red" title "Fr^{-1}"



} else {

plot "steady_tavg_eta.dat" \
   i 0 u ($2**0.5/$8):($12/$8):($13/$8)   w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady: Wrms',\
"" i 0 u ($2**0.5/$8):($94/$8):($95/$8)   w yerrorbars pt 4 ps 2 lc rgb "red" title '       Wlam',\
"" i 0 u ($2**0.5/$8):($96/$8):($97/$8)   w yerrorbars pt 4 ps 2 lc rgb "blue" title '       Wturb',\
"stoch_tavg_eta.dat" \
   i 0 u ($2**0.5/$8):($12/$8):($13/$8)   w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Stoch.: Wrms',\
"" i 0 u ($2**0.5/$8):($94/$8):($95/$8)   w yerrorbars pt 9 ps 2 lc rgb "red" title '       Wlam',\
"" i 0 u ($2**0.5/$8):($96/$8):($97/$8)   w yerrorbars pt 9 ps 2 lc rgb "blue" title '       Wturb',\

}

} # end block 2 }}}

# -------------------------------------------------------------

perform_block_3 = 0

# {{{ Third Plot - set perform_block_3 to 1 to run this block

if (perform_block_3) {

set xlabel "Fr_{eff}^{-1}"
set key bottom left
set ylabel "w_{rms}/u_{h,rms}"
set title "Targeted Avg: Effective Fr^{-1} v w_{rms}"
set log xy
#set yrange [1e-3:10]
set xrange [.1:100]

if (plot_Re1000) {

plot "steady_tavg_eta.dat" \
   i 0 u (($2**0.5)/$8):($12/$8):($13/$8) w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady: Wrms',\
"" i 0 u (($2**0.5)/$8):($62/$8):($63/$8) w yerrorbars pt 4 ps 2 lc rgb "red"         title '        Wlam',\
"" i 0 u (($2**0.5)/$8):($72/$8):($73/$8) w yerrorbars pt 4 ps 2 lc rgb "blue"        title '        Wturb',\
"" i 2 u (($2**0.5)/$8):($12/$8):($13/$8) w yerrorbars pt 4 ps 2 lc rgb "black"       title 'Re1000: Wrms',\
"" i 2 u (($2**0.5)/$8):($62/$8):($63/$8) w yerrorbars pt 4 ps 2 lc rgb "black"       title '        Wlam',\
"" i 2 u (($2**0.5)/$8):($72/$8):($73/$8) w yerrorbars pt 4 ps 2 lc rgb "black"       title '        Wturb',\
"stoch_tavg_eta.dat" \
   i 0 u (($2**0.5)/$8):($12/$8):($13/$8) w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Stoch.: Wrms',\
"" i 0 u (($2**0.5)/$8):($62/$8):($63/$8) w yerrorbars pt 9 ps 2 lc rgb "red"         title '        Wlam',\
"" i 0 u (($2**0.5)/$8):($72/$8):($73/$8) w yerrorbars pt 9 ps 2 lc rgb "blue"        title '        Wturb',\
"" i 1 u (($2**0.5)/$8):($12/$8):($13/$8) w yerrorbars pt 9 ps 2 lc rgb "black"       title 'Re1000: Wrms',\
"" i 1 u (($2**0.5)/$8):($62/$8):($63/$8) w yerrorbars pt 9 ps 2 lc rgb "black"       title '        Wlam',\
"" i 1 u (($2**0.5)/$8):($72/$8):($73/$8) w yerrorbars pt 9 ps 2 lc rgb "black"       title '        Wturb',\

} else {

plot "steady_tavg_eta.dat" \
   i 0 u (($2**0.5)/$8):($12/$8):($13/$8) w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady: Wrms',\
"" i 0 u (($2**0.5)/$8):($62/$8):($63/$8) w yerrorbars pt 4 ps 2 lc rgb "red"         title '        Wlam',\
"" i 0 u (($2**0.5)/$8):($72/$8):($73/$8) w yerrorbars pt 4 ps 2 lc rgb "blue"        title '        Wturb',\
"stoch_tavg_eta.dat" \
   i 0 u (($2**0.5)/$8):($12/$8):($13/$8) w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Stoch.: Wrms',\
"" i 0 u (($2**0.5)/$8):($62/$8):($63/$8) w yerrorbars pt 9 ps 2 lc rgb "red"         title '        Wlam',\
"" i 0 u (($2**0.5)/$8):($72/$8):($73/$8) w yerrorbars pt 9 ps 2 lc rgb "blue"        title '        Wturb',\

}


} # end of block 3 }}}

# -------------------------------------------------------------

perform_block_4 = 0

# {{{ Fourth Plot - set perform_block_4 to 1 to run this block


if (perform_block_4) {

set xlabel "Fr_{eff}^{-1}"
set key top right
set ylabel "Delta w_{rms}"
set title "Extraction Routine Comparsion"
set log xy

plot "steady_tavg_eta.dat" \
   i 0 u (($2**0.5)/$8):(abs($62-$94)/$8):(($63-$95)/$8) w yerrorbars pt 4 ps 2 lc rgb "red"         title '        Wlam',\
"" i 0 u (($2**0.5)/$8):(abs($72-$96)/$8):(($73-$97)/$8) w yerrorbars pt 4 ps 2 lc rgb "blue"        title '        Wturb',\
"stoch_tavg_eta.dat" \
   i 0 u (($2**0.5)/$8):(abs($62-$94)/$8):(($63-$95)/$8) w yerrorbars pt 9 ps 2 lc rgb "red"         title '        Wlam',\
"" i 0 u (($2**0.5)/$8):(abs($72-$96)/$8):(($73-$97)/$8) w yerrorbars pt 9 ps 2 lc rgb "blue"        title '        Wturb',\

#plot "steady_tavg_eta.dat" \
   #i 0 u (($2**0.5)/$8):($43/$8) pt 5 ps 2 lc rgb "dark-violet" title "Steady: w_{rms}",\
#"" i 0 u (($2**0.5)):(0.25*(($2**0.5)**(-1))*$41 + 0.35*(($2**0.5)**(-0.5))*$42) pt 9 ps 2 lc rgb "green" title "Steady: w_{recon}",\
#"stoch_tavg_eta.dat" \
   #i 0 u (($2**0.5)/$8):($43/$8) pt 9 ps 2 lc rgb "dark-violet" title "Stoch: w_{rms}",\
#"" i 0 u (($2**0.5)):(0.25*(($2**0.5)**(-1))*$41 + 0.35*(($2**0.5)**(-0.5))*$42) pt 9 ps 2 lc rgb "green" title "Stoch: w_{recon}",\

} # end of block 4 }}}

# -------------------------------------------------------------

perform_block_5 = 0

# {{{ Fifth Plot - set perform_block_5 to 1 to run this block

if (perform_block_5) {

set xlabel "Fr_{eff} V_{Turb, eff}"
set key bottom right
set ylabel "w_{rms, eff}"
set title "Targeted Avg: Effective V_{Turb} v w_{rms}"
set log xy

if (plot_Re1000) {

plot "steady_tavg_eta.dat" \
   i 0 u ($92/($2**0.5/$8)):($12/$8):($13/$8) w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady: Wrms',\
"" i 0 u ($92/($2**0.5/$8)):($62/$8):($63/$8) w yerrorbars pt 4 ps 2 lc rgb "red"         title '        Wlam',\
"" i 0 u ($92/($2**0.5/$8)):($72/$8):($73/$8) w yerrorbars pt 4 ps 2 lc rgb "blue"        title '        Wturb',\
"" i 2 u ($92/($2**0.5/$8)):($12/$8):($13/$8) w yerrorbars pt 4 ps 2 lc rgb "black"       title 'Re1000: Wrms',\
"" i 2 u ($92/($2**0.5/$8)):($62/$8):($63/$8) w yerrorbars pt 4 ps 2 lc rgb "black"       title '        Wlam',\
"" i 2 u ($92/($2**0.5/$8)):($72/$8):($73/$8) w yerrorbars pt 4 ps 2 lc rgb "black"       title '        Wturb',\
"stoch_tavg_eta.dat" \
   i 0 u ($92/($2**0.5/$8)):($12/$8):($13/$8) w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Stoch.: Wrms',\
"" i 0 u ($92/($2**0.5/$8)):($62/$8):($63/$8) w yerrorbars pt 9 ps 2 lc rgb "red"         title '        Wlam',\
"" i 0 u ($92/($2**0.5/$8)):($72/$8):($73/$8) w yerrorbars pt 9 ps 2 lc rgb "blue"        title '        Wturb',\
"" i 1 u ($92/($2**0.5/$8)):($12/$8):($13/$8) w yerrorbars pt 9 ps 2 lc rgb "black"       title 'Re1000: Wrms',\
"" i 1 u ($92/($2**0.5/$8)):($62/$8):($63/$8) w yerrorbars pt 9 ps 2 lc rgb "black"       title '        Wlam',\
"" i 1 u ($92/($2**0.5/$8)):($72/$8):($73/$8) w yerrorbars pt 9 ps 2 lc rgb "black"       title '        Wturb',\
0.37*x**0.5 w l dt 2 lc rgb "red" title '0.37 V_{Turb}^{0.5}'

} else {

plot "steady_tavg_eta.dat" \
   i 0 u ($92/($2**0.5/$8)):($12/$8):($13/$8) w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady: Wrms',\
"" i 0 u ($92/($2**0.5/$8)):($62/$8):($63/$8) w yerrorbars pt 4 ps 2 lc rgb "red"         title '        Wlam',\
"" i 0 u ($92/($2**0.5/$8)):($72/$8):($73/$8) w yerrorbars pt 4 ps 2 lc rgb "blue"        title '        Wturb',\
"stoch_tavg_eta.dat" \
   i 0 u ($92/($2**0.5/$8)):($12/$8):($13/$8) w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Stoch.: Wrms',\
"" i 0 u ($92/($2**0.5/$8)):($62/$8):($63/$8) w yerrorbars pt 9 ps 2 lc rgb "red"         title '        Wlam',\
"" i 0 u ($92/($2**0.5/$8)):($72/$8):($73/$8) w yerrorbars pt 9 ps 2 lc rgb "blue"        title '        Wturb',\
[1e-4:1e-1] 0.8*x**0.5 w l dt 2 lc rgb "red" title '0.37 V_{Turb}^{0.5}'

}


} # end of block 5 }}}

# -------------------------------------------------------------

perform_block_6 = 0

# {{{ Sixth Plot - set perform_block_6 to 1 to run this block

if (perform_block_6) {

set xlabel "V_{Turb, eff}"
set key bottom right
set ylabel "w_{rms, eff}"
set title "Targeted Avg: Effective V_{Turb} v w_{rms}"
set log xy

if (plot_Re1000) {

plot "steady_tavg_eta.dat" \
   i 0 u ($92):($12/$8):($13/$8) w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady: Wrms',\
"" i 0 u ($92):($62/$8):($63/$8) w yerrorbars pt 4 ps 2 lc rgb "red"         title '        Wlam',\
"" i 0 u ($92):($72/$8):($73/$8) w yerrorbars pt 4 ps 2 lc rgb "blue"        title '        Wturb',\
"" i 2 u ($92):($12/$8):($13/$8) w yerrorbars pt 4 ps 2 lc rgb "black"       title 'Re1000: Wrms',\
"" i 2 u ($92):($62/$8):($63/$8) w yerrorbars pt 4 ps 2 lc rgb "black"       title '        Wlam',\
"" i 2 u ($92):($72/$8):($73/$8) w yerrorbars pt 4 ps 2 lc rgb "black"       title '        Wturb',\
"stoch_tavg_eta.dat" \
   i 0 u ($92):($12/$8):($13/$8) w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Stoch.: Wrms',\
"" i 0 u ($92):($62/$8):($63/$8) w yerrorbars pt 9 ps 2 lc rgb "red"         title '        Wlam',\
"" i 0 u ($92):($72/$8):($73/$8) w yerrorbars pt 9 ps 2 lc rgb "blue"        title '        Wturb',\
"" i 1 u ($92):($12/$8):($13/$8) w yerrorbars pt 9 ps 2 lc rgb "black"       title 'Re1000: Wrms',\
"" i 1 u ($92):($62/$8):($63/$8) w yerrorbars pt 9 ps 2 lc rgb "black"       title '        Wlam',\
"" i 1 u ($92):($72/$8):($73/$8) w yerrorbars pt 9 ps 2 lc rgb "black"       title '        Wturb',\
0.37*x**0.5 w l dt 2 lc rgb "red" title '0.37 V_{Turb}^{0.5}'

} else {

plot "steady_tavg_eta.dat" \
   i 0 u ($92):($12/$8):($13/$8) w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady: Wrms',\
"" i 0 u ($92):($62/$8):($63/$8) w yerrorbars pt 4 ps 2 lc rgb "red"         title '        Wlam',\
"" i 0 u ($92):($72/$8):($73/$8) w yerrorbars pt 4 ps 2 lc rgb "blue"        title '        Wturb',\
"stoch_tavg_eta.dat" \
   i 0 u ($92):($12/$8):($13/$8) w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Stoch.: Wrms',\
"" i 0 u ($92):($62/$8):($63/$8) w yerrorbars pt 9 ps 2 lc rgb "red"         title '        Wlam',\
"" i 0 u ($92):($72/$8):($73/$8) w yerrorbars pt 9 ps 2 lc rgb "blue"        title '        Wturb',\
[1e-4:1] 0.37*x**0.5 w l dt 2 lc rgb "red" title '0.37 V_{Turb}^{0.5}'

}


} # end of block 6 }}}

# -------------------------------------------------------------

perform_block_7 = 0

# {{{ Seventh Plot - set perform_block_7 to 1 to run this block


if (perform_block_7) {

set xlabel "w_{Lam, Targeted Avg}"
set ylabel "w_{No-Turb, Vortz Weighted}"
set yrange [1e-2:1]
set xrange [1e-2:1]
set key bottom right
#set title "Extraction Routine Comparsion"
set log xy

plot "steady_tavg_eta.dat" \
   i 0 u 62:94 w yerrorbars pt 4 ps 2 lc rgb "black"            title 'Steady (Re, Pr) = (600, 0.1)',\
"" i 1 u 62:94 w yerrorbars pt 4 ps 2 lc rgb "blue"             title 'Steady (Re, Pr) = (600, 0.05)',\
"" i 2 u 62:94 w yerrorbars pt 4 ps 2 lc rgb "red"              title 'Steady (Re, Pr) = (1000, 0.1)',\
"" i 3 u 62:94 w yerrorbars pt 4 ps 2 lc rgb "dark-violet"      title 'Steady (Re, Pr) = (300 ,0.1)',\
"" i 4 u 62:94 w yerrorbars pt 4 ps 2 lc rgb "green"            title 'Steady (Re, Pr) = (1000, 0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u 62:94 w yerrorbars pt 9 ps 2 lc rgb "black"            title 'Stoch. (Re, Pr) = (600, 0.1)',\
"" i 1 u 62:94 w yerrorbars pt 9 ps 2 lc rgb "red"              title 'Stoch. (Re, Pr) = (1000, 0.1)',\
x dt 2 lc rgb "black"

#plot "steady_tavg_eta.dat" \
   #i 0 u (($2**0.5)/$8):($43/$8) pt 5 ps 2 lc rgb "dark-violet" title "Steady: w_{rms}",\
#"" i 0 u (($2**0.5)):(0.25*(($2**0.5)**(-1))*$41 + 0.35*(($2**0.5)**(-0.5))*$42) pt 9 ps 2 lc rgb "green" title "Steady: w_{recon}",\
#"stoch_tavg_eta.dat" \
   #i 0 u (($2**0.5)/$8):($43/$8) pt 9 ps 2 lc rgb "dark-violet" title "Stoch: w_{rms}",\
#"" i 0 u (($2**0.5)):(0.25*(($2**0.5)**(-1))*$41 + 0.35*(($2**0.5)**(-0.5))*$42) pt 9 ps 2 lc rgb "green" title "Stoch: w_{recon}",\

} # end of block 7 }}}

# -------------------------------------------------------------

perform_block_8 = 0

# {{{ Eighth Plot - set perform_block_8 to 1 to run this block


if (perform_block_8) {

set xlabel "Fr_{eff}"
set key bottom right
set ylabel "Fr_{h}"
#set title "Fr_{eff} v.s. Fr_{h}"
set format x "10^{%T}"
set format y "10^{%T}"
set log xy

plot "steady_tavg_eta.dat" \
   i 0 u (($2**-0.5)*$8):(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2))) w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady (Re,Pr) = (600,0.1)',\
"" i 1 u (($2**-0.5)*$8):(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2))) w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady (Re,Pr) = (600,0.05)',\
"" i 2 u (($2**-0.5)*$8):(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2))) w yerrorbars pt 4 ps 2 lc rgb "red" title 'Steady (Re,Pr) = (1000,0.1)',\
"" i 3 u (($2**-0.5)*$8):(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2))) w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady (Re,Pr) = (300,0.1)',\
"" i 4 u (($2**-0.5)*$8):(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2))) w yerrorbars pt 4 ps 2 lc rgb "green" title 'Steady (Re,Pr) = (1000,0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u (($2**-0.5)*$8):(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2))) w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stochastic (Re,Pr) = (600,0.1)',\
"" i 1 u (($2**-0.5)*$8):(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2))) w yerrorbars pt 9 ps 2 lc rgb "red" title 'Stochastic (Re,Pr) = (1000,0.1)',\
0.06*x title '0.06x'

} # end of block 8 }}}

# -------------------------------------------------------------

perform_block_9 = 0

# {{{ Ninth Plot - set perform_block_9 to 1 to run this block


if (perform_block_9) {

set xlabel "Fr_{h, emergent}^{-1}"
set ylabel "w_{rms}"
set title "Fr_{h, emergent}^{-1} v.s. w_{rms}"
set key bottom left
set log xy

plot "steady_tavg_eta.dat" \
   i 0 u (1./(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2)))):($12/$8):($13/$8)\
   w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady Re=600 w_{rms}',\
"" i 0 u (1./(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2)))):($94/$8):($95/$8)\
   w yerrorbars pt 4 ps 2 lc rgb "red" title 'Steady Re=600 w_{Lam}',\
"" i 0 u (1./(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2)))):($96/$8):($97/$8)\
   w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady Re=600 w_{Turb}',\
"" i 2 u (1./(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2)))):($12/$8):($13/$8)\
   w yerrorbars pt 6 ps 2 lc rgb "dark-violet" title 'Steady Re=1000 w_{rms}',\
"" i 2 u (1./(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2)))):($94/$8):($95/$8)\
   w yerrorbars pt 6 ps 2 lc rgb "red" title 'Steady Re=1000 w_{Lam}',\
"" i 2 u (1./(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2)))):($96/$8):($97/$8)\
   w yerrorbars pt 6 ps 2 lc rgb "blue" title 'Steady Re=1000 w_{Turb}',\
"stoch_tavg_eta.dat" \
   i 0 u (1./(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2)))):($12/$8):($13/$8)\
   w yerrorbars pt 9 ps 2 lc rgb "dark-violet" title 'Stoch. Re=600 w_{rms}',\
"" i 0 u (1./(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2)))):($94/$8):($95/$8)\
   w yerrorbars pt 9 ps 2 lc rgb "red" title 'Stoch. Re=600 w_{Lam}',\
"" i 0 u (1./(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2)))):($96/$8):($97/$8)\
   w yerrorbars pt 9 ps 2 lc rgb "blue" title 'Stoch. Re=600 w_{Turb}',\
"" i 1 u (1./(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2)))):($12/$8):($13/$8)\
   w yerrorbars pt 13 ps 2 lc rgb "dark-violet" title 'Stoch. Re=1000 w_{rms}',\
"" i 1 u (1./(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2)))):($94/$8):($95/$8)\
   w yerrorbars pt 13 ps 2 lc rgb "red" title 'Stoch. Re=1000 w_{Lam}',\
"" i 1 u (1./(($2**-0.5)*$16/($1*0.5*($8**2 + $12**2)))):($96/$8):($97/$8)\
   w yerrorbars pt 13 ps 2 lc rgb "blue" title 'Stoch. Re=1000 w_{Turb}',\

} # end of block 9 }}}

# -------------------------------------------------------------

perform_block_10 = 0

# {{{ tenth Plot - set perform_block_10 to 1 to run this block


if (perform_block_10) {

set xlabel "w_{Turb, Targetted Avg}"
set key bottom right
set ylabel "w_{Turb, Vortz Weighted}"
set yrange [1e-2:1]
set xrange [1e-2:1]
#set title "Extraction Routine Comparsion"
set log xy

plot "steady_tavg_eta.dat" \
   i 0 u 72:96 w yerrorbars pt 4 ps 2 lc rgb "black"            title 'Steady (Re, Pr) = (600, 0.1)',\
"" i 1 u 72:96 w yerrorbars pt 4 ps 2 lc rgb "blue"             title 'Steady (Re, Pr) = (600, 0.05)',\
"" i 2 u 72:96 w yerrorbars pt 4 ps 2 lc rgb "red"              title 'Steady (Re, Pr) = (1000, 0.1)',\
"" i 3 u 72:96 w yerrorbars pt 4 ps 2 lc rgb "dark-violet"      title 'Steady (Re, Pr) = (300 ,0.1)',\
"" i 4 u 72:96 w yerrorbars pt 4 ps 2 lc rgb "green"            title 'Steady (Re, Pr) = (1000, 0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u 72:96 w yerrorbars pt 9 ps 2 lc rgb "black"            title 'Stoch. (Re, Pr) = (600, 0.1)',\
"" i 1 u 72:96 w yerrorbars pt 9 ps 2 lc rgb "red"              title 'Stoch. (Re, Pr) = (1000, 0.1)',\
x dt 2 lc rgb "black"

#plot "steady_tavg_eta.dat" \
   #i 0 u (($2**0.5)/$8):($43/$8) pt 5 ps 2 lc rgb "dark-violet" title "Steady: w_{rms}",\
#"" i 0 u (($2**0.5)):(0.25*(($2**0.5)**(-1))*$41 + 0.35*(($2**0.5)**(-0.5))*$42) pt 9 ps 2 lc rgb "green" title "Steady: w_{recon}",\
#"stoch_tavg_eta.dat" \
   #i 0 u (($2**0.5)/$8):($43/$8) pt 9 ps 2 lc rgb "dark-violet" title "Stoch: w_{rms}",\
#"" i 0 u (($2**0.5)):(0.25*(($2**0.5)**(-1))*$41 + 0.35*(($2**0.5)**(-0.5))*$42) pt 9 ps 2 lc rgb "green" title "Stoch: w_{recon}",\

} # end of block 7 }}}

# -------------------------------------------------------------


