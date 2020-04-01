import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from astropy.io import fits
from astropy.wcs import WCS

DATADIR = "data/maps20"
HDU = "SCALED"
LINEID = "heii"
STRETCH = 0.4
Q = 5

def load_and_scale_channel(vel, v1, v2, ihdu=HDU, nancolor=0.5):
    """Load channel given central velocity wrt systemic

    Then scale so that [0 -> 1] maps to [v1 -> v2]
    """
    velstring = f"{int(-40.0 + vel - 10.0):+04d}{int(-40.0 + vel + 10.0):+04d}"
    hdu = fits.open(f"{DATADIR}/turtle-slits-{LINEID}{velstring}.fits")[ihdu]
    im = (hdu.data - v1) / (v2 - v1)
    im[~np.isfinite(im)] = nancolor
    return im, WCS(hdu)


image_r, w = load_and_scale_channel(20, 0.0, 0.004, nancolor=0.1)
image_g, _ = load_and_scale_channel(0, 0.0, 0.005, nancolor=0.1)
image_b, _ = load_and_scale_channel(-20, 0.0, 0.002, nancolor=0.02)
im = make_lupton_rgb(image_r, image_g, image_b, stretch=STRETCH, Q=Q)

hdu_hst = fits.open("propermotions/oiii/oiii-ago1997-win.fits")[0]
whst = WCS(hdu_hst)

figfile = "turtle-rgb-map-heii-shell.pdf"

fig, ax = plt.subplots(figsize=(5, 5), subplot_kw=dict(projection=w))
ax.imshow(im)
levs = [2.5, 5, 10, 20, 40, 80] 
lws = np.array([0.3, 0.5, 0.7, 1.0, 1.3, 1.6])
ax.contour(
    hdu_hst.data,
    levels=levs,
    linewidths=lws + 0.8,
    colors="k",
    transform=ax.get_transform(whst),
)
ax.contour(
    hdu_hst.data,
    levels=levs,
    linewidths=lws,
    colors="y",
    transform=ax.get_transform(whst),
)
ax.set(
    xlim=[190, 320],
    ylim=[200, 330],
    xlabel="RA (J2000)",
    ylabel="Dec (J2000)",
)


fig.tight_layout(
    rect=[0.18, 0.04, 1.0, 1.0]
)
fig.savefig(figfile, dpi=200)
print(figfile, end="")
