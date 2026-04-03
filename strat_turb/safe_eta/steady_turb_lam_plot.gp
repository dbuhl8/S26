reset
set xlabel "t"
set key center right

lb = 0
#ub = 300

idx = 7

f(x) = a
fit [lb:ub] f(x) "safe_nonrotating_steady_eta.dat"\
   i idx u 6:36 via a
vlam = a
fit [lb:ub] f(x) "safe_nonrotating_steady_eta.dat"\
   i idx u 6:37 via a
vturb = a
fit [lb:ub] f(x) "safe_nonrotating_steady_eta.dat"\
   i idx u 6:38 via a
vlam_Fr = a
fit [lb:ub] f(x) "safe_nonrotating_steady_eta.dat"\
   i idx u 6:39 via a
vturb_Fr = a
print "Mean of VTurb_Fr: ", a
fit [lb:ub] f(x) "safe_nonrotating_steady_eta.dat"\
   i idx u 6:40 via a
vlam_Fr_vortz = a
fit [lb:ub] f(x) "safe_nonrotating_steady_eta.dat"\
   i idx u 6:41 via a
vturb_Fr_vortz = a
fit [lb:ub] f(x) "safe_nonrotating_steady_eta.dat"\
   i idx u 6:7 via a
uh = a

plot "safe_nonrotating_steady_eta.dat" \
   i idx u 6:36 pt 5 ps 2 lc rgb "blue" title 'VLam',\
"" i idx u 6:37 pt 5 ps 2 lc rgb "red" title 'VTurb',\
"" i idx u 6:38 pt 9 ps 2 lc rgb "blue" title 'VLam_{Fr}',\
"" i idx u 6:39 pt 9 ps 2 lc rgb "red" title 'VTurb_{Fr}',\
"" i idx u 6:40 pt 7 ps 2 lc rgb "blue" title 'VLam_{Fr,vortz}',\
"" i idx u 6:41 pt 7 ps 2 lc rgb "red" title 'VTurb_{Fr_vortz}',\
"" i idx u 6:7 pt 5 ps 2 lc rgb "green" title 'U_h',\


print "Mean of VLam: ", vlam
print "Mean of VLam_Fr: ", vlam_Fr
print "Mean of VLam_Fr_vortz: ", vlam_Fr_vortz
print "Mean of VTurb: ", vturb
print "Mean of VTurb_Fr: ", vturb_Fr
print "Mean of VTurb_Fr_vortz: ", vturb_Fr_vortz
print "Mean of Uh: ", uh
