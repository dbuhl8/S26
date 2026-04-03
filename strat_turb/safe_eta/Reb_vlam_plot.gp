

reset
set multiplot layout 1,3 columnsfirst 

idx = 0


# -------------------------------------------------------------

# First Plot

set xlabel "Re_{b, eff}"
set key top right
#set ylabel "Eta" rotate by 0 
set title "Normal"

plot "steady_tavg_eta.dat" \
 i idx u (($8**3)*$1/$2):37 pt 5 ps 2 lc rgb "blue" title 'Steady VLam',\
"stoch_tavg_eta.dat" \
 i idx u (($8**3)*$1/$2):37 pt 7 ps 2 lc rgb "blue" title 'Stoch. VLam',\
#"" i idx u (($8**3)*$1/$2):38 pt 5 ps 2 lc rgb "red" title 'Steady VTurb',\
#"" i idx u (($8**3)*$1/$2):38 pt 7 ps 2 lc rgb "red" title 'Stoch. VTurb',\


# -------------------------------------------------------------

# Second Plot

set xlabel "Re_{b, eff}"
set key top right
#set ylabel "Eta" rotate by 0 
set title "Using Effective Fr"

plot "steady_tavg_eta.dat" \
 i idx u (($8**3)*$1/$2):39 pt 5 ps 2 lc rgb "blue" title 'Steady VLam',\
"stoch_tavg_eta.dat" \
 i idx u (($8**3)*$1/$2):39 pt 7 ps 2 lc rgb "blue" title 'Stoch. VLam',\
#"" i idx u (($8**3)*$1/$2):40 pt 5 ps 2 lc rgb "red" title 'Steady VTurb',\
#"" i idx u (($8**3)*$1/$2):40 pt 7 ps 2 lc rgb "red" title 'Stoch. VTurb',\


# -------------------------------------------------------------

# Third Plot

set xlabel "Re_{b, eff}"
set key top right
#set ylabel "Eta" rotate by 0 
set title "Using Effective Fr & {/Symbol w}_z"

plot "steady_tavg_eta.dat" \
 i idx u (($8**3)*$1/$2):41 pt 5 ps 2 lc rgb "blue" title 'Steady VLam',\
"stoch_tavg_eta.dat" \
 i idx u (($8**3)*$1/$2):41 pt 7 ps 2 lc rgb "blue" title 'Stoch. VLam',\
#"" i idx u (($8**3)*$1/$2):42 pt 5 ps 2 lc rgb "red" title 'Steady VTurb',\
#"" i idx u (($8**3)*$1/$2):42 pt 7 ps 2 lc rgb "red" title 'Stoch. VTurb',\
