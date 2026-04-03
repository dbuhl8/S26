reset
set multiplot layout 3,1 columnsfirst 

# file key: 
# steady_tavg_eta.dat: 
# idx 0: Re600_Pe60
# idx 1: Re600_Pe30
# idx 2: Re1000_Pe100
# idx 3: Re300_Pe30
# idx 4: Re1000_Pe10

# stoch_tavg_eta.dat:
# idx 0: Re600_Pe60
# idx 1: Re1000_Pe100

# in each file:
# 43 - wrms
# 44 - wlam
# 45 - wturb
# 46 - wlam_Fr
# 47 - wturb_Fr
# 48 - wlam_Fr_vortz
# 49 - wturb_Fr_vortz
# 50 - wlam_weighted
# 51 - wturb_weighted
# 52 - wlam_Fr_weighted
# 53 - wturb_Fr_weighed

idx = 0


# -------------------------------------------------------------

# First Plot

set xlabel "Fr^{-1}"
set key top right
#set ylabel "Eta" rotate by 0 
set title "Normal"
set log xy

plot "steady_tavg_eta.dat" \
   i idx u ($2**0.5):43 pt 5 ps 2 lc rgb "magenta" title 'Steady: Wrms',\
"" i idx u ($2**0.5):44 pt 5 ps 2 lc rgb "red" title '       Wlam',\
"" i idx u ($2**0.5):45 pt 5 ps 2 lc rgb "blue" title '       Wturb',\
"" i idx u ($2**0.5):50 pt 4 ps 2 lc rgb "red" title '       Wlam wght',\
"" i idx u ($2**0.5):51 pt 4 ps 2 lc rgb "blue" title '       Wlam wght',\
"stoch_tavg_eta.dat" \
   i idx u ($2**0.5):43 pt 7 ps 2 lc rgb "magenta" title 'Stoch: Wrms',\
"" i idx u ($2**0.5):44 pt 7 ps 2 lc rgb "red" title '       Wlam',\
"" i idx u ($2**0.5):45 pt 7 ps 2 lc rgb "blue" title '       Wturb',\
"" i idx u ($2**0.5):50 pt 6 ps 2 lc rgb "red" title '       Wlam wght',\
"" i idx u ($2**0.5):51 pt 6 ps 2 lc rgb "blue" title '       Wlam wght',\



# -------------------------------------------------------------

# Second Plot

set xlabel "Fr_{eff}^{-1}"
set key top right
#set ylabel "Eta" rotate by 0 
set title "Using Effective Fr"
set log x


plot "steady_tavg_eta.dat" \
   i idx u (($2**0.5)/$8):43 pt 5 ps 2 lc rgb "magenta" title 'Steady: Wrms',\
"" i idx u (($2**0.5)/$8):46 pt 5 ps 2 lc rgb "red" title '       Wlam',\
"" i idx u (($2**0.5)/$8):47 pt 5 ps 2 lc rgb "blue" title '       Wturb',\
"" i idx u (($2**0.5)/$8):52 pt 4 ps 2 lc rgb "red" title '       Wlam wght',\
"" i idx u (($2**0.5)/$8):53 pt 4 ps 2 lc rgb "blue" title '       Wlam wght',\
"stoch_tavg_eta.dat" \
   i idx u (($2**0.5)/$8):43 pt 7 ps 2 lc rgb "magenta" title 'Stoch: Wrms',\
"" i idx u (($2**0.5)/$8):46 pt 7 ps 2 lc rgb "red" title '       Wlam',\
"" i idx u (($2**0.5)/$8):47 pt 7 ps 2 lc rgb "blue" title '       Wturb',\
"" i idx u (($2**0.5)/$8):52 pt 6 ps 2 lc rgb "red" title '       Wlam wght',\
"" i idx u (($2**0.5)/$8):53 pt 6 ps 2 lc rgb "blue" title '       Wlam wght',\



# -------------------------------------------------------------

# Third Plot

set xlabel "Fr_{eff}^{-1}"
set key top right
#set ylabel "Eta" rotate by 0 
set title "Using Effective Fr & w"
set log x

plot "steady_tavg_eta.dat" \
   i idx u (($2**0.5)/$8):43 pt 5 ps 2 lc rgb "magenta" title 'Steady: Wrms',\
"" i idx u (($2**0.5)/$8):48 pt 5 ps 2 lc rgb "red" title '       Wlam',\
"" i idx u (($2**0.5)/$8):49 pt 5 ps 2 lc rgb "blue" title '       Wturb',\
"" i idx u (($2**0.5)/$8):52 pt 4 ps 2 lc rgb "red" title '       Wlam_wght',\
"" i idx u (($2**0.5)/$8):53 pt 4 ps 2 lc rgb "blue" title '       Wlam_wght',\
"stoch_tavg_eta.dat" \
   i idx u (($2**0.5)/$8):43 pt 7 ps 2 lc rgb "magenta" title 'Stoch: Wrms',\
"" i idx u (($2**0.5)/$8):48 pt 7 ps 2 lc rgb "red" title '       Wlam',\
"" i idx u (($2**0.5)/$8):49 pt 7 ps 2 lc rgb "blue" title '       Wturb',\
"" i idx u (($2**0.5)/$8):52 pt 6 ps 2 lc rgb "red" title '       Wlam_wght',\
"" i idx u (($2**0.5)/$8):53 pt 6 ps 2 lc rgb "blue" title '       Wlam_wght',\



