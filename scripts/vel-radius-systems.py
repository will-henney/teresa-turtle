TAB=[["System", "Type", "PA", "Axis", "V", "eV", "R", "eR", "t", "et"], ["Inner shell A NNW", 1, 340, "A^+", 49, 7, 0.07, 0.01, 1.4, 0.2], ["Outer lobe I+/VII+", 3, 340, "A^+", 89, 27, 0.17, 0.01, 1.9, 0.7], ["North knot", 2, 340, "A^+", 67, 4, 0.19, 0.01, 2.7, 0.2], ["Inner shell A SSE", 1, 160, "A^-", 55, 7, 0.08, 0.01, 1.2, 0.2], ["Outer lobe I-", 3, 160, "A^-", 81, 26, 0.17, 0.03, 2.0, 0.8], ["Outer lobe VIII-", 3, 160, "A^-", 68, 9, 0.23, 0.01, 3.2, 0.4], ["Outer lobe II+", 3, 170, "C^+", 65, 11, 0.27, 0.02, 3.9, 0.9], ["Outer lobe II-", 3, 355, "C^-", 25, 25, 0.25, 0.08, 0.0, 0.0], ["Inner shell B ENE", 1, 65, "B^-", 43, 4, 0.06, 0.01, 1.2, 0.4], ["NE Blue complex", 2, 60, "", 65, 2, 0.14, 0.04, 2.1, 0.6], ["IS ENE", 1, 50, "B^-", 35, 12, 0.12, 0.04, 3.2, 2.0], ["Inner shell B WSW", 1, 240, "B^+", 36, 4, 0.06, 0.01, 1.5, 0.4], ["IS WSW", 1, 260, "B^+", 38, 20, 0.16, 0.06, 3.8, 3.5], ["SW Red complex", 2, 250, "", 49, 2, 0.2, 0.06, 3.7, 1.2], ["N(E) Red complex", 2, 30, "", 30, 2, 0.26, 0.21, 8.0, 6.7], ["SW Faint Blue complex", 2, 225, "", 41, 3, 0.11, 0.02, 2.6, 0.5], ["Outer lobe X+", 3, 285, "D+", 51, 43, 0.19, 0.02, 3.7, 3.3], ["Outer lobe IV+", 3, 315, "D+", 37, 6, 0.11, 0.01, 2.9, 0.6], ["NW knot ([N II])", 2, 310, "", 24, 4, 0.18, 0.1, 6.7, 4.5], ["N(W) Red complex", 2, 300, "", 27, 3, 0.16, 0.08, 6.0, 3.1], ["Outer lobe X-/IV-/V+", 3, 120, "D-", 26, 13, 0.14, 0.04, 5.2, 3.3], ["SE Blue complex", 2, 135, "", 54, 8, 0.07, 0.02, 1.1, 0.4]]
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_color_codes()

figfile = "vel-radius-systems.pdf"
t = Table(rows=TAB[1:], names=TAB[0])
t["PA"] = t["PA"] % 180

mshell = t["Type"] == 1
mknot = t["Type"] == 2
mlobe = t["Type"] == 3


fig, ax = plt.subplots()

ax.errorbar("R", "V", xerr="eR", yerr="eV",
            fmt="none", data=t, lw=0.3, zorder=-100)
for mask, marker in [
        [mshell, "o"],
        [mknot, "P"],
        [mlobe, "d"],
]:
    scat = ax.scatter("R", "V",
                      marker=marker,
                      c="PA", cmap="twilight_shifted", vmin=0.0, vmax=180.0,
                      edgecolors="k", linewidths=0.5,
                      s=100,
                      data=t[mask])

R0 = 1.0
for time in 1000, 2000, 4000, 8000, 16000:
    # 1 pc in 1000 yr means speed of 977.79 km/s
    V0 = 977.79*1000/time
    ax.plot([0, R0], [0, V0], color="k", lw=0.5, zorder=-100)

fig.colorbar(scat, ax=ax).set_label(r"Position angle modulo $180^{\circ}$")
ax.set(
    xlim=[0, 0.35],
    ylim=[0, 120],
    ylabel="Velocity, km/s",
    xlabel="Distance, pc",
)
sns.despine()
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end="")
