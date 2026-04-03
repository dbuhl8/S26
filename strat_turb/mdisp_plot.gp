




# plotting routine to plot Ri frac / ReB
set xlabel "Re_b"
set ylabel "M Disp"
set log xy
set key top left
plot "nonrotating_data_new.dat" \
i 0 u ($1/$2):($22/$1):($23/$1) w yerrorbars pt 2 ps 2 lc rgb "red" title "Stoch. Re = 600", \
"" i 1 u ($1/$2):($22/$1):($23/$1) w yerrorbars pt 2 ps 2 lc rgb "blue" title "Stoch. Re = 1000",\
"Final_data.dat" i 5 u ($1/$2):($24/$1):($25/$1) w yerrorbars pt 5 ps 2 lc rgb "red" title "Steady Re = 600", \
"" i 9 u ($1/$2):($24/$1):($25/$1) w yerrorbars pt 5 ps 2 lc rgb "blue" title "Steady Re = 1000"

#i 0 u ((($10**2+$12**2)**1.5)*$1/$2):($22/$1):($23/$1) w yerrorbars pt 2 ps 2 lc rgb "red" title "Stoch. Re = 600", \
#"" i 1 u ((($10**2+$12**2)**1.5)*$1/$2):($22/$1):($23/$1) w yerrorbars pt 2 ps 2 lc rgb "blue" title "Stoch. Re = 1000",\
#"Final_data.dat" i 5 u ((($10**2+$12**2)**1.5)*$1/$2):($24/$1):($25/$1) w yerrorbars pt 5 ps 2 lc rgb "red" title "Steady Re = 600", \
#"" i 9 u ((($10**2+$12**2)**1.5)*$1/$2):($24/$1):($25/$1) w yerrorbars pt 5 ps 2 lc rgb "blue" title "Steady Re = 1000"


