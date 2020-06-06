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
from pathlib import Path
from astropy.table import Table


datadir = Path("../quireza-tables")

tab1 = Table.read(
    str(datadir / "table1.dat"),
    format="ascii.cds",
    readme=str(datadir / "ReadMe"),
)
tab1

tab6 = Table.read(
    str(datadir / "table6.dat"),
    format="ascii.cds",
    readme=str(datadir / "ReadMe"),
)
tab6

gaiadir = Path("../stang-2020-tables")
tab_g = Table.read(
    str(gaiadir / "apjab59e4t1_mrt.dat"),
    format="ascii.cds",
)
tab_g

import pandas as pd

df = pd.merge(tab1.to_pandas(), tab6.to_pandas(), "outer", "Name")
df

df["R"] = np.round(4.84813681e-3*df["Theta"]*df["dhel"], 3)
df["L"] = df["S5GHz"]*df["dhel"]**2


def fix_colnames(s):
    return s.replace("(", "_").replace(")", "").replace("-", "_")


df.rename(columns=fix_colnames, inplace=True)

# Add some new columns, proportional to physical size and luminosity.

df



qstring = "p_typeIIb + p_typeIII > 0.9 & p_typeIII >  2*p_typeIIb"
qIII = "`Post_class` == 'III'"
qIIb = "`Post_class` == 'IIb'"

df.query(qstring).describe()

df.query(qIII).describe()

df.query(qIIb).describe()

# So this is a sample selected to be mainly Type III, with a soft border with the Type IIb.
#
#

pd.options.display.max_rows = 999
pd.options.display.max_columns = 999

df.query(qstring)

df.corr()

df.query('Name == "NGC 6210"')

df.query('Name == "NGC 6826"')

df.query('Name == "J320"')

df.query('Name == "Fg1"')

df.query('Name == "NGC 2392"')

df.query('Name == "NGC 6543"')

df.query('Name == "NGC 7009"')

df.query('Name == "IC4634"')

# A new selection that has R = 0.05 to 0.2 and L = 300 to 1000 and D < 3 kpc

qstring_physical = "0.08 < R < 0.2 & 500 < L < 3000 & dhel < 4.2"

df.query(qIII).query(qstring_physical)

df.query(qIIb).query(qstring_physical)

# Cross-reference with the Gaia data

dfIII = df.query(qIII).query(qstring_physical)

dfG = tab_g.to_pandas().rename(columns={"ID": "PNG"})

dfGIII = dfG.query("PNG in @dfIII['PNG']")

df_all = pd.merge(dfIII, dfGIII, "outer", "PNG")

df_all.corr()

df_all["D"] = 1.0/df_all["plx"]

df_all["Rpc"] = 10**df_all["logRpc"]

df_all[["Name", "dhel", "D", "theta", "Theta", "R", "Rpc"]]

df_all

dfG.query("PNG in ['197.8+17.3', '037.7-34.5', '096.4+29.9']")

# ## Testing stuff

d = pd.DataFrame({"a": [1, 2], "x_y+12": [3, 4]})
d

d.query("`x_y+12` < 4")

pd.__version__


