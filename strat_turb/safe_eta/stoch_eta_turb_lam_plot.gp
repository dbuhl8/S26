reset
set xlabel "t"
set key top right
#set yrange [0:1]
lb = 30
ub = 80

idx = 5

f(x) = a
fit [lb:ub] f(x) "safe_nonrotating_stoch_eta.dat"\
   i idx u 6:11 via a
eta = a
fit [lb:ub] f(x) "safe_nonrotating_stoch_eta.dat"\
   i idx u 6:15 via a
eta_lam = a
fit [lb:ub] f(x) "safe_nonrotating_stoch_eta.dat"\
   i idx u 6:19 via a
eta_turb = a
fit [lb:ub] f(x) "safe_nonrotating_stoch_eta.dat"\
   i idx u 6:23 via a
eta_lam_Fr = a
fit [lb:ub] f(x) "safe_nonrotating_stoch_eta.dat"\
   i idx u 6:27 via a
eta_turb_Fr = a
print "Mean of VTurb_Fr: ", a
fit [lb:ub] f(x) "safe_nonrotating_stoch_eta.dat"\
   i idx u 6:31 via a
eta_lam_Fr_vortz = a
fit [lb:ub] f(x) "safe_nonrotating_stoch_eta.dat"\
   i idx u 6:35 via a
eta_turb_Fr_vortz = a
fit [lb:ub] f(x) "safe_nonrotating_steady_eta.dat"\
i idx u 6:7 via a
uh = a


plot "safe_nonrotating_stoch_eta.dat" \
   i idx u 6:7 pt 5 ps 2 lc rgb "green" title 'Eta',\
"" i idx u 6:15 pt 5 ps 2 lc rgb "blue" title 'Eta Lam',\
"" i idx u 6:19 pt 5 ps 2 lc rgb "red" title 'Eta Turb',\
"" i idx u 6:23 pt 9 ps 2 lc rgb "blue" title 'Eta Lam_{Fr}',\
"" i idx u 6:27 pt 9 ps 2 lc rgb "red" title 'Eta Turb_{Fr}',\
"" i idx u 6:31 pt 7 ps 2 lc rgb "blue" title 'Eta Lam_{Fr,vortz}',\
"" i idx u 6:35 pt 7 ps 2 lc rgb "red" title 'Eta Turb_{Fr_vortz}',\


print "Mean of Eta: ", eta
print "Mean of Eta Lam: ", eta_lam
print "Mean of Eta Lam_Fr: ", eta_lam_Fr
print "Mean of Eta Lam_Fr_vortz: ", eta_lam_Fr_vortz
print "Mean of Eta Turb: ", eta_turb
print "Mean of Eta Turb_Fr: ", eta_turb_Fr
print "Mean of Eta Turb_Fr_vortz: ", eta_turb_Fr_vortz
print "Mean of Uh: ", uh
