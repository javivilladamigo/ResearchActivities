
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import mplhep as hep
plt.style.use([hep.style.CMS])
plt.clf()
plt.close("all")
QCD_color, Z_color, H_color = ["#af272f", "#608fbe", "darkorange"]



def quotient(a, b): # quotient uncertainty
    y = a / b
    uy = y * (1. / a + 1. / b) ** 0.5
    return y, uy


## DATA EXTRACTED FROM ANALYSIS IN efficiencies.ipynb ## 


Nevents     =   {
                "Ntot" :
                    {"QCD" : 21335910, "Z" : 459000, "H" : 448000},
                "kin" :
                    {"Mu18 Ph32"                        : {"QCD" : 34,  "Z" : 98630,    "H" : 160587},

                    "Mu15 Mu10 Ph15"                    : {"QCD" : 94,  "Z" : 94330,    "H" : 128356},
                    "Mu18 Ph24 dR04"                    : {"QCD" : 29,  "Z" : 113505,   "H" : 168665},
                    "Mu15 Ph20 2mumuM4"                 : {"QCD" : 28,  "Z" : 125836,   "H" : 173615},
                    "Mu10 Mu05 Ph15 2mumuM4"            : {"QCD" : 97,  "Z" : 136929,   "H" : 178371},
                    "Mu10 Mu05 Ph22.7 2mumuM4"          : {"QCD" : 35,  "Z" : 129513,   "H" : 174286}}

                } 

efficiency  =   {
                "selection" :
                    {"Mu18 Ph32"                                : {"QCD" : None, "Z" : None, "H" : None}, # must have the same names as Nevents["kin"]
                    

                    "Mu15 Mu10 Ph15"                            : {"QCD" : None, "Z" : None, "H" : None},
                    "Mu18 Ph24 dR04"                            : {"QCD" : None, "Z" : None, "H" : None},
                    "Mu15 Ph20 2mumuM4"                         : {"QCD" : None, "Z" : None, "H" : None},
                    "Mu10 Mu05 Ph15 2mumuM4"                    : {"QCD" : None, "Z" : None, "H" : None},
                    "Mu10 Mu05 Ph22.7 2mumuM4"                  : {"QCD" : None, "Z" : None, "H" : None}}
                }



for kin in efficiency["selection"]:
    for nature in ("QCD", "Z", "H"):
        efficiency["selection"][kin][nature] = quotient(Nevents["kin"][kin][nature], Nevents["Ntot"][nature])




fig, (ax, ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize = (10, 10)); fig.patch.set_facecolor("xkcd:white")


ax.tick_params(which = 'both', axis = 'y', direction='in', right = True, top = True, length = 6)
ax.tick_params(which = 'major', axis = 'y', labelsize = 30)
ax.minorticks_on()
ax.tick_params(which = 'minor', axis = 'y', direction='in', right = True, top = True, length = 3)
ax.tick_params(which = "minor", axis = "x", direction = "in", length = 0)
ax.tick_params(axis = 'x', labelbottom = False)

ax1.tick_params(which = 'both', axis = 'y', direction='in', right = True, top = True, length = 6)
ax1.tick_params(which = 'major', axis = 'y', labelsize = 25)
ax1.minorticks_on()
ax1.tick_params(which = 'minor', axis = 'y', direction='in', right = True, top = True, length = 3)
ax1.tick_params(which = "minor", axis = "x", direction = "in", length = 0)

ax.tick_params(which = 'major', axis = 'x', direction='in', bottom = True, top = True, length = 6)
ax1.tick_params(which = 'major', axis = 'x', direction='in', bottom = True, top = True, length = 6)


_labels = ["Mu18 Ph32", "Mu15 Mu10 Ph15", "Mu18 Ph24 dR04", "Mu15 Ph20 2mumuM4", "Mu10 Mu05 Ph15 2mumuM4", "Mu10 Mu05 Ph22.7 2mumuM4"]
_x = np.arange(len(_labels))

for i, kin in enumerate(_labels):
    # QCD
    ax1.plot(i, 100 * efficiency["selection"][kin]["QCD"][0], marker = 'o', markersize = 8., c = QCD_color, lw = 0)
    ax1.errorbar(i, 100 * efficiency["selection"][kin]["QCD"][0], yerr = 100 * efficiency["selection"][kin]["QCD"][1] if 100 * efficiency["selection"][kin]["QCD"][1] + 100 * efficiency["selection"][kin]["QCD"][0] <= 100. else 100 - 100 * efficiency["selection"][kin]["QCD"][0], fmt = "o", capsize = 3., elinewidth = 1., lw = 0, color = QCD_color, label = "QCD" if i == 0 else None)

    # Z
    ax.plot(i, 100 * efficiency["selection"][kin]["Z"][0], marker = 'o', markersize = 8., c = Z_color, lw = 0)
    ax.errorbar(i, 100 * efficiency["selection"][kin]["Z"][0], yerr = 10 * 100 * efficiency["selection"][kin]["Z"][1] if 100 * efficiency["selection"][kin]["Z"][1] + 100 * efficiency["selection"][kin]["Z"][0] <= 100. else 100 - 100 * efficiency["selection"][kin]["Z"][0], fmt = "o", capsize = 3., elinewidth = 1., lw = 0, color = Z_color, label = "Z" if i == 0 else None)

    # H 
    ax.plot(i, 100 * efficiency["selection"][kin]["H"][0], marker = 'o', markersize = 8., c = H_color, lw = 0)
    ax.errorbar(i, 100 * efficiency["selection"][kin]["H"][0], yerr = 10 * 100 * efficiency["selection"][kin]["H"][1] if 100 * efficiency["selection"][kin]["H"][1] + 100 * efficiency["selection"][kin]["H"][0] <= 100. else 100 - 100 * efficiency["selection"][kin]["H"][0], fmt = "o", capsize = 3., elinewidth = 1., lw = 0, color = H_color, label = "H" if i == 0 else None)


ax.axvline(x = 0.5, linestyle = "dashed", color = "darkgrey", lw = 5.)
ax1.axvline(x = 0.5, linestyle = "dashed", color = "darkgrey", lw = 5.)

fig.text(x = 0.096, y = 0.575, s = "Simulated", fontsize = 30)
fig.text(x = 0.50, y = 0.575, s = "Calculated", fontsize = 30)

ax.set_ylabel("Efficiency (%)", labelpad = 30, fontsize = 30)
ax1.set_xlabel("Selection criteria", labelpad = 30, fontsize = 30)

ax.set_xlim(-0.5, _x[-1]+0.5)
ax1.set_xlim(-0.5, _x[-1]+0.5)

ax.set_ylim(18, 42)
ax1.set_ylim(0.00005, 0.00055)



rect_ax = patches.Rectangle((-0.5, 18), width = 0.5 - (-0.5), height = 42-18, linewidth = 0.5, color = None, facecolor = 'darkgrey', alpha = 0.35)
rect_ax1 = patches.Rectangle((-0.5, 0.00005), width = 0.5 - (-0.5), height = 0.00055-0.00005, linewidth = 0.5, color = None, facecolor = 'darkgrey', alpha = 0.35)

ax.add_patch(rect_ax)
ax1.add_patch(rect_ax1)

#fig.legend(loc = (0.84,0.73), fontsize = 30)
fig.legend(loc = "upper right", fontsize = 30)
ax1.set_xticks(_x, _labels, rotation = 45, fontsize = 25)
ax1.set_yticks(np.arange(0.0001, 0.0006, 0.0002))

ax.grid(which = "minor", axis = "y")
ax1.grid(which = "minor", axis = "y")

fig.subplots_adjust(left = 0.076, bottom = 0.449, right = 0.888, top = 0.97, wspace = 0.2, hspace = 0.2)

fig.show()