import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from astropy.io import fits
from astropy.wcs import WCS

DATADIR = "data/maps20"
HDU = "SCALED"
LINEID = "nii"
STRETCH = 0.4
Q = 5

def load_and_scale_channel(vel, v1, v2, ihdu=HDU):
    """Load channel given central velocity wrt systemic

    Then scale so that [0 -> 1] maps to [v1 -> v2]
    """
    velstring = f"{int(-40.0 + vel - 10.0):+04d}{int(-40.0 + vel + 10.0):+04d}"
    hdu = fits.open(f"{DATADIR}/turtle-slits-{LINEID}{velstring}.fits")[ihdu]
    return (hdu.data - v1) / (v2 - v1), WCS(hdu)


image_r, w = load_and_scale_channel(-10, 0.0, 0.03)
image_g, _ = load_and_scale_channel(-30, 0.0, 0.06)
image_b, _ = load_and_scale_channel(-50, 0.0, 0.08)
image = make_lupton_rgb(image_r, image_g, image_b, stretch=STRETCH, Q=Q)

figfile = "turtle-rgb-map-nii-knots.png"

fig, ax = plt.subplots(subplot_kw=dict(projection=w))
ax.imshow(image)

ax.set(
    xlim=[190, 320],
    ylim=[210, 300],
    xlabel="RA (J2000)",
    ylabel="Dec (J2000)",
)
fig.savefig(figfile, dpi=200)
print(figfile, end="")
