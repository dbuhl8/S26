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
    set output "bflux_model_rms.eps"
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
# idx 2: Re300_Pe30
# }}}

# -------------------------------------------------------------

perform_block_1 = 0
# <wb>/(wb)

# {{{ First Plot

if (perform_block_1) {

set xlabel "Re_G"
set key top right
set ylabel "<wb>/wb" rotate by 0
set format x "10^{%T}"
#set yrange [0:1]
set log x

plot "steady_tavg_eta.dat" \
   i 0 u ($16/$2):($68/($12*$54)) w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady (Re,Pr) = (600,0.1)',\
"" i 1 u ($16/$2):($68/($12*$54)) w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady (Re,Pr) = (600,0.05)',\
"" i 2 u ($16/$2):($68/($12*$54)) w yerrorbars pt 4 ps 2 lc rgb "red" title 'Steady (Re,Pr) = (1000,0.1)',\
"" i 3 u ($16/$2):($68/($12*$54)) w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady (Re,Pr) = (300,0.1)',\
"" i 4 u ($16/$2):($68/($12*$54)) w yerrorbars pt 4 ps 2 lc rgb "forest-green" title 'Steady (Re,Pr) = (1000,0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u ($16/$2):($68/($12*$54)) w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. (Re,Pr) = (600,0.1)',\
"" i 1 u ($16/$2):($68/($12*$54)) w yerrorbars pt 9 ps 2 lc rgb "red" title 'Stoch. (Re,Pr) = (1000,0.1)',\
#"" i 2 u 2:($68/($12*$54)) w yerrorbars pt 9 ps 2 lc rgb "forest-green" title 'Steady MDisp',\

f(x) = a
fit [20:1e6] f(x) "stoch_tavg_eta.dat" u ($16/$2):($68/($12*$54)) via a
stoch_c = a
fit [20:1e6] f(x) "steady_tavg_eta.dat" u ($16/$2):($68/($12*$54)) via a
steady_c = a

print "Stochastic Coefficient: ", stoch_c
print "Steady Coefficient: ", steady_c


} # end of block 1 }}}

# -------------------------------------------------------------

perform_block_2 = 0

# <wb>_turb/(w_turb*b_turb)

# {{{ Second Plot 

if (perform_block_2) {

set xlabel "Re_G"
set key top right
set ylabel "<wb>_{Turb}/w_{Turb}b_{Turb}" rotate by 0
set format x "10^{%T}"
#set yrange [0:1]
set log x

plot "steady_tavg_eta.dat" \
   i 0 u ($16/$2):($72/($32*$58)) w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady (Re,Pr) = (600,0.1)',\
"" i 1 u ($16/$2):($72/($32*$58)) w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady (Re,Pr) = (600,0.05)',\
"" i 2 u ($16/$2):($72/($32*$58)) w yerrorbars pt 4 ps 2 lc rgb "red" title 'Steady (Re,Pr) = (1000,0.1)',\
"" i 3 u ($16/$2):($72/($32*$58)) w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady (Re,Pr) = (300,0.1)',\
"" i 4 u ($16/$2):($72/($32*$58)) w yerrorbars pt 4 ps 2 lc rgb "forest-green" title 'Steady (Re,Pr) = (1000,0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u ($16/$2):($72/($32*$58)) w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. (Re,Pr) = (600,0.1)',\
"" i 1 u ($16/$2):($72/($32*$58)) w yerrorbars pt 9 ps 2 lc rgb "red" title 'Stoch. (Re,Pr) = (1000,0.1)',\
#"" i 2 u 2:($68/($12*$54)) w yerrorbars pt 9 ps 2 lc rgb "forest-green" title 'Steady MDisp',\

f(x) = a
fit [20:1e6] f(x) "stoch_tavg_eta.dat" u ($16/$2):($72/($32*$58)) via a
stoch_c = a
fit [20:1e6] f(x) "steady_tavg_eta.dat" u ($16/$2):($72/($32*$58)) via a
steady_c = a

print "Stochastic Coefficient: ", stoch_c
print "Steady Coefficient: ", steady_c


} # end of block 2 }}}

# -------------------------------------------------------------

perform_block_3 = 0

# <wb>_lam/(w_lam*b_lam)

# {{{ Third Plot 

if (perform_block_3) {

set xlabel "Re_G"
set key top right
set ylabel "<wb>_{Lam}/w_{Lam}b_{Lam}" rotate by 0
set format x "10^{%T}"
#set yrange [0:1]
set log x

plot "steady_tavg_eta.dat" \
   i 0 u ($16/$2):($70/($22*$56)) w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady (Re,Pr) = (600,0.1)',\
"" i 1 u ($16/$2):($70/($22*$56)) w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady (Re,Pr) = (600,0.05)',\
"" i 2 u ($16/$2):($70/($22*$56)) w yerrorbars pt 4 ps 2 lc rgb "red" title 'Steady (Re,Pr) = (1000,0.1)',\
"" i 3 u ($16/$2):($70/($22*$56)) w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady (Re,Pr) = (300,0.1)',\
"" i 4 u ($16/$2):($70/($22*$56)) w yerrorbars pt 4 ps 2 lc rgb "forest-green" title 'Steady (Re,Pr) = (1000,0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u ($16/$2):($70/($22*$56)) w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. (Re,Pr) = (600,0.1)',\
"" i 1 u ($16/$2):($70/($22*$56)) w yerrorbars pt 9 ps 2 lc rgb "red" title 'Stoch. (Re,Pr) = (1000,0.1)',\
#"" i 2 u 2:($68/($12*$54)) w yerrorbars pt 9 ps 2 lc rgb "forest-green" title 'Steady MDisp',\

f(x) = a
fit [20:1e6] f(x) "stoch_tavg_eta.dat" u ($16/$2):($70/($22*$56)) via a
stoch_c = a
fit [20:1e6] f(x) "steady_tavg_eta.dat" u ($16/$2):($70/($22*$56)) via a
steady_c = a

print "Stochastic Coefficient: ", stoch_c
print "Steady Coefficient: ", steady_c


} # end of block 3 }}}

# -------------------------------------------------------------

perform_block_4 = 0

# Reconstruction v.s. BFlux (turb)

# {{{ Fourth Plot 

if (perform_block_4) {

set xlabel "Reconstruction"
set key bottom right
set ylabel "<wb>" rotate by 0
set format x "10^{%T}"
set format y "10^{%T}"
#set yrange [0:1]
set log xy

# model = sqrt(VLam(c_lam*w_lam*b_lam)^2 + VTurb(c_turb*w_turb*b_turb))
# = (($42*(1.2*$56*$22)**2 + $44*(1.2*$58*$32)**2)**0.5)

plot "steady_tavg_eta.dat" \
   i 0 u (($42*(1.2*$56*$22)**2 + $44*(1.2*$58*$32)**2)**0.5):70:71 w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady (Re,Pr) = (600,0.1)',\
"" i 1 u (($42*(1.2*$56*$22)**2 + $44*(1.2*$58*$32)**2)**0.5):70 w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady (Re,Pr) = (600,0.05)',\
"" i 2 u (($42*(1.2*$56*$22)**2 + $44*(1.2*$58*$32)**2)**0.5):70 w yerrorbars pt 4 ps 2 lc rgb "red" title 'Steady (Re,Pr) = (1000,0.1)',\
"" i 3 u (($42*(1.2*$56*$22)**2 + $44*(1.2*$58*$32)**2)**0.5):70 w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady (Re,Pr) = (300,0.1)',\
"" i 4 u (($42*(1.2*$56*$22)**2 + $44*(1.2*$58*$32)**2)**0.5):70 w yerrorbars pt 4 ps 2 lc rgb "forest-green" title 'Steady (Re,Pr) = (1000,0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u (($42*(1.2*$56*$22)**2 + $44*(1.2*$58*$32)**2)**0.5):70 w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. (Re,Pr) = (600,0.1)',\
"" i 1 u (($42*(1.2*$56*$22)**2 + $44*(1.2*$58*$32)**2)**0.5):70 w yerrorbars pt 9 ps 2 lc rgb "red" title 'Stoch. (Re,Pr) = (1000,0.1)',\
x w l lw 2 dt 2 lc rgb "black"
#"" i 2 u 2:($68/($12*$54)) w yerrorbars pt 9 ps 2 lc rgb "forest-green" title 'Steady MDisp',\

f(x) = a
fit [20:1e6] f(x) "stoch_tavg_eta.dat" u ($16/$2):($68/($12*$54)) via a
stoch_c = a
fit [20:1e6] f(x) "steady_tavg_eta.dat" u ($16/$2):($68/($12*$54)) via a
steady_c = a

print "Stochastic Coefficient: ", stoch_c
print "Steady Coefficient: ", steady_c


} # end of block 4 }}}

# -------------------------------------------------------------

perform_block_5 = 1

# Model v.s. BFlux (turb)

# {{{ Fifth Plot 

if (perform_block_5) {

set xlabel "Model"
set key bottom right
set ylabel "<wb>" rotate by 0
set format x "10^{%T}"
set format y "10^{%T}"
#set yrange [0:1]
set log xy

# model = sqrt(VLam(c_lam*w_lam*b_lam)^2 + VTurb(c_turb*w_turb*b_turb))
# model = sqrt(VLam(c_lam*Fr**3)^2 + VTurb(c_turb*Fr**1.5)^2)
#       = (($42*(1.2*(($2**-0.5)*$8)**3)**2 + $44*(1.2*(($2**-0.5)*$8)**2)**2)**0.5)

plot "steady_tavg_eta.dat" \
   i 0 u (($42*(1.2*(($2**-0.5)*$8)**3)**2 + $44*(1.2*(($2**-0.5)*$8)**2)**2)**0.5):70:71 w yerrorbars pt 4 ps 2 lc rgb "black" title 'Steady (Re,Pr) = (600,0.1)',\
"" i 1 u (($42*(1.2*(($2**-0.5)*$8)**3)**2 + $44*(1.2*(($2**-0.5)*$8)**2)**2)**0.5):70 w yerrorbars pt 4 ps 2 lc rgb "blue" title 'Steady (Re,Pr) = (600,0.05)',\
"" i 2 u (($42*(1.2*(($2**-0.5)*$8)**3)**2 + $44*(1.2*(($2**-0.5)*$8)**2)**2)**0.5):70 w yerrorbars pt 4 ps 2 lc rgb "red" title 'Steady (Re,Pr) = (1000,0.1)',\
"" i 3 u (($42*(1.2*(($2**-0.5)*$8)**3)**2 + $44*(1.2*(($2**-0.5)*$8)**2)**2)**0.5):70 w yerrorbars pt 4 ps 2 lc rgb "dark-violet" title 'Steady (Re,Pr) = (300,0.1)',\
"" i 4 u (($42*(1.2*(($2**-0.5)*$8)**3)**2 + $44*(1.2*(($2**-0.5)*$8)**2)**2)**0.5):70 w yerrorbars pt 4 ps 2 lc rgb "forest-green" title 'Steady (Re,Pr) = (1000,0.01)',\
"stoch_tavg_eta.dat" \
   i 0 u (($42*(1.2*(($2**-0.5)*$8)**3)**2 + $44*(1.2*(($2**-0.5)*$8)**2)**2)**0.5):70 w yerrorbars pt 9 ps 2 lc rgb "black" title 'Stoch. (Re,Pr) = (600,0.1)',\
"" i 1 u (($42*(1.2*(($2**-0.5)*$8)**3)**2 + $44*(1.2*(($2**-0.5)*$8)**2)**2)**0.5):70 w yerrorbars pt 9 ps 2 lc rgb "red" title 'Stoch. (Re,Pr) = (1000,0.1)',\
x w l lw 2 dt 2 lc rgb "black"
#"" i 2 u 2:($68/($12*$54)) w yerrorbars pt 9 ps 2 lc rgb "forest-green" title 'Steady MDisp',\

f(x) = a
fit [20:1e6] f(x) "stoch_tavg_eta.dat" u ($16/$2):($68/($12*$54)) via a
stoch_c = a
fit [20:1e6] f(x) "steady_tavg_eta.dat" u ($16/$2):($68/($12*$54)) via a
steady_c = a

print "Stochastic Coefficient: ", stoch_c
print "Steady Coefficient: ", steady_c


} # end of block 5 }}}

# -------------------------------------------------------------





