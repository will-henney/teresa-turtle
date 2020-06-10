# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import numpy as np
from astropy.table import Table
from pathlib import Path

datadir = Path("../cspn-tables")

byte_by_byte_description = """
Byte-by-byte Description of file:
--------------------------------------------------------------------------------
   Bytes Format Units     Label       Explanations
--------------------------------------------------------------------------------
   1-  5  I5    ---       N           Track point number
   7- 15  F9.6  [Lsun]    logL        logarithm of the stellar luminosity
  17- 25  F9.6  [K]       logTeff     logarithm of the effective temperature
  27- 35  F9.6  [cm/s2]   logg        logarithm of the surface gravity
  40- 51  F12.4 yr        t           Age since the point at LogTeff=3.85
  53- 61  F9.6  ---       Menv        Fractional mass of the envelope
  63- 71  F9.6  Msun      Mstar       Total mass of the star
  73- 82  F10.6 [Msun/yr] log(-dM/dt)  Logarithm of the Mass Loss Rate,
                                       log(-dMstar/dt)
--------------------------------------------------------------------------------
"""


def read_tracks(datafile):
    """Read each Miller–Bertolami track into a separate astropy.table
    
    Input argument `datafile` is a CDS file containing all tracks 
    for a given metallicity, e.g., "0100_t03.dat"
    
    Returns list of tables. Each table has a metadata "comments" field 
    that contains additional info (mass and surface composition). 
    """
    with open(datafile) as f:
        # Each track is separated by two blank lines
        tracks = f.read().split("\n\n\n")[:-1]
        tables = []
        for track in tracks:
            lines = track.split("\n")
            metadata = lines[:6]
            data = lines[7:]
            datastring = "\n".join(
                [byte_by_byte_description] + data
            )
            table = Table.read(datastring, format="ascii.cds")
            table.meta["comments"] = metadata
            tables.append(table)
    return tables


tabs = read_tracks(datadir / "0100_t03.dat")
[_.meta for _ in tabs]

10**3.589592

tabs[2].show_in_notebook()

from matplotlib import pyplot as plt
import seaborn as sns
# %matplotlib inline
sns.set_context("talk")
sns.set_color_codes()


# Plot of effective temperature versus gravity, which we compare with the turtle observed values.

def extract_masses(data):
    _, Mi, Mf, _ = data.meta["comments"][2].split()
    return round(float(Mi), 2), round(float(Mf), 3)


fig, ax = plt.subplots(figsize=(8, 8))
ax.axhline(4.875, color="k", lw=0.5)
ax.axhline(4.792, color="k", lw=0.5)
ax.axvspan(4.6, 5.0, color="k", alpha=0.1)
ax.axvline(4.8, color="k", lw=0.5)
for data in tabs:
    try:
        Mi, Mf = extract_masses(data)
        label = f"({Mi}, {Mf})"
    except:
        continue
    ax.plot(
        "logg", "logTeff",
        data=data, label=label,
    )
ax.legend()
ax.set(
    xlabel="$\log_{10}\, g$",
    ylabel="$\log_{10}\, T_{\mathrm{eff}}$",
    xlim=[4.0, 6.0],
    ylim=[4.7, 5.0],
)
sns.despine()
None

# From this we conclude that the $2\,M_\odot$ model is the best fit, but we cannot rule out $1.25\,M_\odot$ to $3\,M_\odot$.

# Next, look at the timescales.

fig, [axL, axMd, ax] = plt.subplots(3, 1, sharex=True, figsize=(12, 8))
ax.axhline(62000, color="k", lw=0.5)
ax.axhline(75000, color="k", lw=0.5)
ax.axhline(20000, color="r", ls="--", lw=0.8)
for axx in axL, axMd, ax:
    axx.axvline(0.0, color="k", lw=0.5)
    axx.axvspan(-100, 100, color="m", alpha=0.05)
for data in tabs:
    try:
        Mi, Mf = extract_masses(data)
        label = f"({Mi}, {Mf})"
    except:
        continue
    data["Teff"] = 10**data["logTeff"]
    data["L"] = 10**data["logL"]
    data["Mdot"] = 10**data["log(-dM/dt)"]
    ax.plot(
        "t", "Teff",
        data=data, label=label,
    )
    axMd.plot(
        "t", "Mdot",
        data=data, label=label,
    )
    axL.plot(
        "t", "L",
        data=data,
    )
ax.legend()
ax.set(
    xlabel="time, years",
    ylabel="$T_{\mathrm{eff}}$, K",
    xlim=[-100000, 100000],
    ylim=[3000, 3e5],
)
ax.set_xscale("symlog", linthreshx=100)
ax.set_yscale("log")
axL.set(
    ylim=[0.0, 1.5e4],
    ylabel="$L / L_\odot$",
)
axMd.set(
    ylim=[1e-8, 1e-5],
    ylabel="$d M / dt$, $M_\odot$/yr",
    yscale="log",
)
sns.despine()


# So, conclusion from this is that we need the mass to be 1.25 to 1.5.  If it is 1.25, then the MB models have it not going through a C-star phase, which is needed to explain the low C/O ratio in nebula.
#
# On the other hand, if it was a triple interaction that ejected the envelope, it might have happened before the AGB got to end of its evolution. 

def make_table_of_times(tabs, Teff):
    logTeff = np.log10(Teff)
    tTkey = f"t({Teff})"
    rslts = {
        "Mi": [],
        "Mf": [],
        tTkey: [],
        "t_cross": [],
        "t_tr": [],
    }
    for data in tabs:
        Mi, Mf = extract_masses(data)
        rslts["Mi"].append(Mi)
        rslts["Mf"].append(Mf)
        # First time to reach given Teff
        mask = data["logTeff"] >= logTeff
        tT = data[mask]["t"].min()
        rslts[tTkey].append(tT)
        # Time to cross to maximum Teff
        icross = data["logTeff"].argmax()
        rslts["t_cross"].append(data[icross]["t"])
        # Transition time before t = 0
        rslts["t_tr"].append(-data["t"].min())
    return Table(rslts)


times = make_table_of_times(tabs, 75000)

times

make_table_of_times(tabs_tb2, 75000)

make_table_of_times(tabs_Z0010, 75000)

# ## HR diagram with cousin nebulae

# +
coustab = Table.read("../doc/cousins.ecsv")
coustab = coustab[:-1]
coustab["log T"] = np.round(3.0 + np.log10(0.5*(coustab["T_eff"] + coustab["T_Z(He II)"])), 2)
coustab["log L"] = np.round(np.log10(coustab["L"]), 2)

def fcol(label):
    peimb, symm = label.split('-')
    color = {"III": "b", "IIb": "c", "IIa": "r"}[peimb]
    return color

def fmark(label):
    peimb, symm = label.split('-')
    marker = {"S": "o", "C": "d", "M": "P", "A": "*"}[symm]
    return marker



coustab["c"] = [fcol(_) for _ in coustab["WJH"]]
coustab["marker"] = [fmark(_) for _ in coustab["WJH"]]


m_use = coustab["Use"] == "Y"
coustab[m_use]
# -

for _x, _y, _c, _m in coustab[m_use][["log L", "log T", "c", "marker"]]:
    print(_x, _y, _c, _m)

beartab = Table.read("../doc/bear-triples.ecsv")
beartab

# +
fig, ax = plt.subplots(figsize=(8, 8))
#ax.axvspan(4.7, 5.0, 0.6, 0.9, color="k", alpha=0.1)
lw = 0.5
for data in tabs:
    try:
        Mi, Mf = extract_masses(data)
        label = f"({Mi}, {Mf})"
    except:
        continue
    ax.plot(
        "logTeff", "logL",
        data=data, label=label,
        zorder=-100, c="k", lw=lw,
    )
    lw += 0.2
    


m = beartab["Bear"] == "Triple"
ax.scatter(
    "log T", "log L", data=beartab[m], 
    color="g", marker="v", s=100, label="Triple",
    facecolors="none",
)
ax.scatter(
    "log T", "log L", data=beartab[~m], 
    color="g", marker="v", s=50, label="Likely",
    facecolors="none",
)
for _n, _x, _y, _c, _m in coustab[m_use][["Name", "log T", "log L", "c", "marker"]]:
    size = 800 if "6210" in _n else 100
    ax.scatter(_x, _y, c=_c, marker=_m, s=size)
ax.legend()
ax.set(
    ylabel="$\log_{10}\, L/L_\odot$",
    xlabel="$\log_{10}\, T_{\mathrm{eff}}$",
    xlim=[5.5, 3.8],
    ylim=[2.0, None],
)
sns.despine()
fig.savefig("hr-planetaries.pdf")
None
# -

m = beartab["Bear"] == "Triple"
~m

# ls



# ### Less important stuff

# Now, we try the appendix B tracks (how are they different?)

tabs_tb2 = read_tracks(datadir / "0100_tb2.dat")

tabs_tb2[0].meta

fig, ax = plt.subplots(figsize=(8, 8))
ax.axhline(4.875, color="k", lw=0.5)
ax.axhline(4.792, color="k", lw=0.5)
ax.axvspan(4.6, 5.0, color="k", alpha=0.1)
ax.axvline(4.8, color="k", lw=0.5)
for data in tabs_tb2:
    try:
        _, Mi, Mf, _ = data.meta["comments"][2].split()
        Mi = round(float(Mi), 2)
        Mf = round(float(Mf), 3)
        label = f"({Mi}, {Mf})"
    except:
        continue
    ax.plot(
        "logg", "logTeff",
        data=data, label=label,
    )
ax.legend()
ax.set(
    xlabel="$\log_{10}\, g$",
    ylabel="$\log_{10}\, T_{\mathrm{eff}}$",
    xlim=[4.0, 6.0],
    ylim=[4.7, 5.0],
)
sns.despine()
None

fig, ax = plt.subplots(figsize=(8, 8))
ax.axhline(62000, color="k", lw=0.5)
ax.axhline(75000, color="k", lw=0.5)
ax.axvline(0.0, color="k", lw=0.5)
for data in tabs_tb2:
    try:
        _, Mi, Mf, _ = data.meta["comments"][2].split()
        Mi = round(float(Mi), 2)
        Mf = round(float(Mf), 3)
        label = f"({Mi}, {Mf})"
    except:
        continue
    data["Teff"] = 10**data["logTeff"]
    ax.plot(
        "t", "Teff",
        data=data, label=label,
    )
ax.legend()
ax.set(
    xlabel="time, years",
    ylabel="$T_{\mathrm{eff}}$, K",
    xlim=[-10000, 10000],
    ylim=[3000, 3e5],
)
ax.set_xscale("symlog", linthreshx=100)
ax.set_yscale("log")
sns.despine()

# Now, the low metallicity tracks

tabs_Z0010 = read_tracks(datadir / "0010_t03.dat")

fig, ax = plt.subplots(figsize=(8, 8))
ax.axhline(4.875, color="k", lw=0.5)
ax.axhline(4.792, color="k", lw=0.5)
ax.axvspan(4.6, 5.0, color="k", alpha=0.1)
ax.axvline(4.8, color="k", lw=0.5)
for data in tabs_Z0010:
    try:
        _, Mi, Mf, _ = data.meta["comments"][2].split()
        Mi = round(float(Mi), 2)
        Mf = round(float(Mf), 3)
        label = f"({Mi}, {Mf})"
    except:
        continue
    ax.plot(
        "logg", "logTeff",
        data=data, label=label,
    )
ax.legend()
ax.set(
    xlabel="$\log_{10}\, g$",
    ylabel="$\log_{10}\, T_{\mathrm{eff}}$",
    xlim=[4.0, 6.0],
    ylim=[4.7, 5.0],
)
sns.despine()
None

fig, ax = plt.subplots(figsize=(8, 8))
ax.axhline(62000, color="k", lw=0.5)
ax.axhline(75000, color="k", lw=0.5)
ax.axvline(0.0, color="k", lw=0.5)
for data in tabs_Z0010:
    try:
        _, Mi, Mf, _ = data.meta["comments"][2].split()
        Mi = round(float(Mi), 2)
        Mf = round(float(Mf), 3)
        label = f"({Mi}, {Mf})"
    except:
        continue
    data["Teff"] = 10**data["logTeff"]
    ax.plot(
        "t", "Teff",
        data=data, label=label,
    )
ax.legend()
ax.set(
    xlabel="time, years",
    ylabel="$T_{\mathrm{eff}}$, K",
    xlim=[-10000, 10000],
    ylim=[3000, 3e5],
)
ax.set_xscale("symlog", linthreshx=100)
ax.set_yscale("log")
sns.despine()


