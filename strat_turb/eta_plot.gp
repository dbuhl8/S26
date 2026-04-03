




# plotting routine to plot Ri frac / ReB
set xlabel "Re_b"
set log x
set ylabel '{/Symbol h}' rotate by 0
set key top left
plot "nonrotating_data_new.dat" \
i 0 u ($1/$2):26:27 w yerrorbars pt 2 ps 2 lc rgb "red" title "Stoch. Re = 600", \
"" i 1 u ($1/$2):26:27 w yerrorbars pt 2 ps 2 lc rgb "blue" title "Stoch. Re = 1000", \
"Final_data.dat" i 5 u ($1/$2):42:43 w yerrorbars pt 5 ps 2 lc rgb "red" title "Steady Re = 600", \
"" i 9 u ($1/$2):42:43 w yerrorbars pt 5 ps 2 lc rgb "blue" title "Steady Re = 1000",\
"" i 8 u ($1/$2):42:43 w yerrorbars pt 5 ps 2 lc rgb "green" title "Steady Pr = 0.01"

#i 0 u ((($10**2+$12**2)**1.5)*$1/$2):26:27 w yerrorbars pt 2 ps 2 lc rgb "red" title "Stoch. Re = 600", \
#"" i 1 u ((($10**2+$12**2)**1.5)*$1/$2):26:27 w yerrorbars pt 2 ps 2 lc rgb "blue" title "Stoch. Re = 1000", \
#"Final_data.dat" i 5 u ((($10**2+$12**2)**1.5)*$1/$2):42:43 w yerrorbars pt 5 ps 2 lc rgb "red" title "Steady Re = 600", \
#"" i 9 u ((($10**2+$12**2)**1.5)*$1/$2):42:43 w yerrorbars pt 5 ps 2 lc rgb "blue" title "Steady Re = 1000",\
#"" i 8 u ((($10**2+$12**2)**1.5)*$1/$2):42:43 w yerrorbars pt 5 ps 2 lc rgb "green" title "Steady Pr = 0.01"


#"Final_data.dat" i 5 u ((($10**2+$12**2)**1.5)*$1/$2):(($2*$18/$4)/($2*$18/$4 + $24/$1)):(($2*$19/$4)/($2*$19/$4 + $25/$1)) w yerrorbars pt 5 ps 2 lc rgb "red" title "Steady Re = 600", \
#"" i 9 u ((($10**2+$12**2)**1.5)*$1/$2):(($2*$18/$4)/($2*$18/$4 + $24/$1)):(($2*$19/$4)/($2*$19/$4 + $25/$1)) w yerrorbars pt 5 ps 2 lc rgb "blue" title "Steady Re = 1000", \

