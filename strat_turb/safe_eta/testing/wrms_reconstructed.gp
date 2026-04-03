reset
#set multiplot layout 3,1 columnsfirst 

idx = 0

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


# -------------------------------------------------------------

# First Plot

set xlabel "Fr^{-1}"
set key top right
#set ylabel "Eta" rotate by 0 
set title "Normal"
set log xy

plot "steady_tavg_eta.dat" \


#   i idx u (($2**0.5)/$8):($42*(($2**0.5)/$8)**-0.5 + 1.9*$41*(($2**0.5)/$8)**-1) pt 5 ps 2 lc rgb "red" title 'Steady Wrms (reconstructed)',\
#"stoch_tavg_eta.dat" \
#   i idx u (($2**0.5)/$8):($42*(($2**0.5)/$8)**-0.5 + 1.9*$41*(($2**0.5)/$8)**-1) pt 7 ps 2 lc rgb "blue" title 'Stoch. Wrms (reconstructed)'

#"" i idx u (($2**0.5)/$8):43 pt 5 ps 2 lc rgb "black" title 'Steady Wrms',\
#"" i idx u (($2**0.5)/$8):43 pt 5 ps 2 lc rgb "black" title 'Steady Wrms',\


# -------------------------------------------------------------

# Second Plot

#set xlabel "Fr_{eff}^{-1}"
#set key top right
#set ylabel "Eta" rotate by 0 
#set title "Using Effective Fr"
#set log x

#plot "steady_tavg_eta.dat" \
#   i idx u (($2**0.5)/$8):($38*x**-0.5 + 1.9*$37*x**-1) pt 5 ps 2 lc rgb "green" title 'Steady Wrms (reconstructed)',\
#"" i idx u (($2**0.5)/$8):43 pt 5 ps 2 lc rgb "black" title 'Steady Wrms',\
#"stoch_tavg_eta.dat" \
#   i idx u (($2**0.5)/$8):($38*x**-0.5 + 1.9*$37*x**-1) pt 5 ps 2 lc rgb "green" title 'Steady Wrms (reconstructed)',\
#"" i idx u (($2**0.5)/$8):43 pt 5 ps 2 lc rgb "black" title 'Steady Wrms',\


# -------------------------------------------------------------

# Third Plot

#set xlabel "Fr_{eff}^{-1}"
#set key top right
#set ylabel "Eta" rotate by 0 
#set title "Using Effective Fr & w"
#set log x



