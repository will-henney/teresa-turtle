TAB=[
    ["Component", "V", "Ratio", "E(R)", "R'", "E(R')"],
    ["N knot", "+24", 0.89, 0.03, 0.7, 0.02],
    ["N Lobe VI-/IX+", "+20", 0.89, 0.05, 0.7, 0.04],
    ["Inner shell", -16, 0.72, 0.06, 0.57, 0.05],
    ["N Red complex", "+23", 0.67, 0.08, 0.53, 0.06],
    ["NE Blue complex", -43, 0.82, 0.02, 0.65, 0.02],
    ["SW Red complex", "+34", 0.79, 0.13, 0.62, 0.1],
    ["SE Blue complex", -28, 0.88, 0.04, 0.69, 0.03],
    ["Total", 0, 0.76, 0.0, 0.6, 0.0]
]
import numpy as np
import pyneb as pn
from astropy.table import Table

data = Table(rows=TAB[1:], names=TAB[0])

T = 1.07e4

S2 = pn.Atom(atom='S2')

dens = S2.getTemDen(
    data["R'"],
    tem=T,
    to_eval="L(6716) / L(6731)",
)

dens_hi = S2.getTemDen(
    data["R'"] - data["E(R')"],
    tem=T,
    to_eval="L(6716) / L(6731)",
)
dens_lo = S2.getTemDen(
    data["R'"] + data["E(R')"],
    tem=T,
    to_eval="L(6716) / L(6731)",
    )

edens = 0.5*(dens_hi - dens_lo)

data["Dens"] = np.round(dens)
data["E+(Dens)"] = np.round(dens_hi - dens)
data["E-(Dens)"] = np.round(dens - dens_lo)

data.remove_columns(['Ratio', 'E(R)'])

outtab = [
    data.colnames,
    None
] + [list(row) for row in data]
