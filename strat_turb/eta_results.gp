




# plotting routine to plot Ri frac / ReB
set xlabel "Re_b"
set ylabel "{/Symbol w}_z" rotate by 0
set key top left
plot "B30Re600Pe60_stoch_eta.dat" \
 i 0 u 6:12 w yerrorbars pt 7 ps 2 lc rgb "red" title "Stoch Lam Eta (global)", \
"" i 0 u 6:16 w yerrorbars pt 7 ps 2 lc rgb "blue" title "Stoch Turb Eta (global)", \
"<cat B30Re600Pe60_stoch/OUT*" u 2:($45) w l lc rgb "green" title "RMS ETA",\
"" u 2:9 w l lc "black" title "RMS WT",\
"B30Re600Pe60_steady_eta.dat" \
 i 0 u 6:12 w yerrorbars pt 7 ps 2 lc rgb "red" title "Steady Lam Eta (global)", \
"" i 0 u 6:16 w yerrorbars pt 7 ps 2 lc rgb "blue" title "Steady Turb Eta (global)", \
"<cat B30Re600Pe60_steady/OUT*" u 2:(($34/2)/($34/2 + $36/600)) w l lc rgb "green" title "RMS ETA"


# weighted adding togetht
#plot "B30Re600Pe60_stoch_eta.dat" \
#i 0 u 6:($12*$19 + $16*$20) w yerrorbars pt 5 ps 2 lc rgb "red" title "Stoch. Lam Eta (global)", \
#"<cat B30Re600Pe60_stoch/OUT*" u 2:(($45/2)/($45/2 + $47/600)) w l lc rgb "green" title "RMS ETA",\
#"B30Re600Pe60_steady_eta.dat" \
#i 0 u 6:($12*$19 + $16*$20) w yerrorbars pt 7 ps 2 lc rgb "red" title "Steady Lam Eta (global)", \
#"<cat B30Re600Pe60_steady/OUT*" u 2:(($34/2)/($34/2 + $36/600)) w l lc rgb "green" title "RMS ETA"


#"" i 0 u 6:16 w yerrorbars pt 5 ps 2 lc rgb "blue" title "Stoch. Turb Eta (global)", \
#i 0 u 6:11 w yerrorbars pt 5 ps 2 lc rgb "red" title "Stoch. Lam Eta (Local)", \
#"" i 0 u 6:15 w yerrorbars pt 5 ps 2 lc rgb "blue" title "Stoch. Turb Eta (local)", \
#"" i 0 u 6:15 w yerrorbars pt 5 ps 2 lc rgb "black" title "Stoch. Turb Eta (local)", \
#i 0 u 6:11 w yerrorbars pt 7 ps 2 lc rgb "red" title "Steady Lam Eta (Local)", \
#"" i 0 u 6:15 w yerrorbars pt 7 ps 2 lc rgb "blue" title "Steady Turb Eta (local)", \


#"Final_data.dat" i 5 u 2:($2*$18/$4):($2*$19/$4) w yerrorbars pt 5 ps 2 lc rgb "red" title "Steady Re = 600", \
#"" i 9 u 2:($2*$18/$4):($2*$19/$4) w yerrorbars pt 5 ps 2 lc rgb "blue" title "Steady Re = 1000"


#i 0 u ((($10**2+$12**2)**1.5)*$1/$2):($2*$16/$4):($2*$17/$4) w yerrorbars pt 2 ps 2 lc rgb "red" title "Stoch. Re = 600", \
#"" i 1 u ((($10**2+$12**2)**1.5)*$1/$2):($2*$16/$4):($2*$17/$4) w yerrorbars pt 2 ps 2 lc rgb "blue" title "Stoch. Re = 1000",\
#"Final_data.dat" i 5 u ((($10**2+$12**2)**1.5)*$1/$2):($2*$18/$4):($2*$19/$4) w yerrorbars pt 5 ps 2 lc rgb "red" title "Steady Re = 600", \
#"" i 9 u ((($10**2+$12**2)**1.5)*$1/$2):($2*$18/$4):($2*$19/$4) w yerrorbars pt 5 ps 2 lc rgb "blue" title "Steady Re = 1000"


