import numpy as np
import matplotlib.pyplot as plt


# extracted data files
steady_tavg_file = "steady_tavg_eta.dat"
stoch_tavg_file = "stoch_tavg_eta.dat"

print("Begin gathering time averaged data")
Forcing = []
Re = []
Fr = []
Pr = []
Uh = []
W = []
EU = []
mdisp = []
bdisp = []
# getting steady data
steady_data = open(steady_tavg_file, 'r')
for line in steady_data.readlines():
    if "#" in line:
        continue
    elif len(line.split()) == 0:
        continue
    data = [float(x) for x in line.split()]
    Forcing.append("Steady")
    Re.append(data[0])
    Fr.append(1./np.sqrt(data[1]))
    Pr.append(data[2])
    Uh.append(data[7])
    W.append(data[11])
    EU.append(0.5*(data[7]**2 + data[11]**2))
    mdisp.append(data[15])
    bdisp.append(data[13])
steady_data.close()

stoch_data = open(stoch_tavg_file, 'r')
for line in stoch_data.readlines():
    if "#" in line:
        continue
    elif len(line.split()) == 0:
        continue
    data = [float(x) for x in line.split()]
    Forcing.append("Stoch")
    Re.append(data[0])
    Fr.append(1./np.sqrt(data[1]))
    Pr.append(data[2])
    Uh.append(data[7])
    W.append(data[11])
    EU.append(0.5*(data[7]**2 + data[11]**2))
    mdisp.append(data[15])
    bdisp.append(data[13])
stoch_data.close()

print("Finished gathering time averaged data")
Forcing = np.array(Forcing)
Re = np.array(Re)
Fr = np.array(Fr)
Pr = np.array(Pr)
Uh = np.array(Uh)
W = np.array(W)
EU = np.array(EU)
mdisp = np.array(mdisp)
bdisp = np.array(bdisp)
idx = np.where(abs(Re - 1000) <= 2)[0]
print(idx)
Forcing = Forcing[idx]
Re = Re[idx]
Fr = Fr[idx]
Pr = Pr[idx]
Uh = Uh[idx]
W = W[idx]
EU = EU[idx]
mdisp = mdisp[idx]
bdisp = bdisp[idx]


headers = ["Forcing","$Re$","$Fr$","$Pr$","$u_{h, rms}$",\
    "$w_{rms}$","$\mathcal{K}$","$|\\nabla u|^2$",\
    "$|\\nabla b|^2$"]

textabular = "||"+'c '*(len(headers)-1)+"c||"
texheader = "\\hline " + " & ".join(headers)
texdata = "\n\\\\ \\hline\n"
newline = texdata
for i in range(len(Forcing)):
   texdata += f"{Forcing[i]} & {int(Re[i])} & {Fr[i]:.3f} & {Pr[i]:.3f} &"+\
   f"{Uh[i]:.3f} & {W[i]:.3f} & {EU[i]:.3f} & {mdisp[i]:.3f} & {bdisp[i]:.3f}"\
   + newline
print("\\begin{tabular}{"+textabular+"}")
print(texheader)
print(texdata,end="")
print("\\end{tabular}")


