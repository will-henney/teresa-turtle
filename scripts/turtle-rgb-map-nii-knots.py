import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from astropy.io import fits
from astropy.wcs import WCS

DATADIR = "data/maps20"
HDU = "SCALED"
LINEID = "nii"
STRETCH = 0.15
Q = 8

def load_and_scale_channel(vel, v1, v2, ihdu=HDU, nancolor=0.5):
    """Load channel given central velocity wrt systemic

    Then scale so that [0 -> 1] maps to [v1 -> v2]
    """
    velstring = f"{int(-40.0 + vel - 10.0):+04d}{int(-40.0 + vel + 10.0):+04d}"
    hdu = fits.open(f"{DATADIR}/turtle-slits-{LINEID}{velstring}.fits")[ihdu]
    im = (hdu.data - v1) / (v2 - v1)
    im[~np.isfinite(im)] = nancolor
    return im, WCS(hdu)


image_r, w = load_and_scale_channel(-10, 0.0, 0.04, nancolor=0.02)
image_g, _ = load_and_scale_channel(-30, 0.0, 0.05, nancolor=0.02)
image_b, _ = load_and_scale_channel(-50, 0.0, 0.06, nancolor=0.02)
# im_blue = make_lupton_rgb(image_r, image_g, image_b, stretch=STRETCH, Q=Q)
im_blue = make_lupton_rgb(
    0.0*image_r,
    (image_g + 2*image_r)/3,
    (2*image_b + image_g)/3,
    stretch=STRETCH, Q=Q)

image_r, w = load_and_scale_channel(50, 0.0, 0.08, nancolor=0.02)
image_g, _ = load_and_scale_channel(30, 0.0, 0.2, nancolor=0.02)
image_b, _ = load_and_scale_channel(10, 0.0, 0.1, nancolor=0.02)
#im_red = make_lupton_rgb(image_r, image_g, image_b, stretch=STRETCH, Q=Q)
im_red = make_lupton_rgb(
    (2*image_r + image_g)/3,
    (image_g + 2*image_b)/3,
    0.0*image_b,
    stretch=STRETCH, Q=Q)

hdu_hst = fits.open("propermotions/nii/n6210n-rescaled-1998.fits")[0]
whst = WCS(hdu_hst)

figfile = "turtle-rgb-map-nii-knots.pdf"

fig, [ax1, ax2] = plt.subplots(
    1, 2,
    figsize=(10, 5),
    sharex=True, sharey=True, 
    subplot_kw=dict(projection=w)
)
ax1.imshow(im_blue)
ax1.contour(
    hdu_hst.data,
    levels=[0.05, 0.1, 0.2, 0.5, 1.0, 2.0],
    linewidths=[0.3, 0.5, 0.7, 1.0, 1.3, 1.6],
    colors="pink",
    transform=ax1.get_transform(whst),
)
ax1.set(
    xlim=[190, 320],
    ylim=[200, 330],
    xlabel="RA (J2000)",
    ylabel="Dec (J2000)",
)

ax2.imshow(im_red)
ax2.contour(
    hdu_hst.data,
    levels=[0.05, 0.1, 0.2, 0.5, 1.0, 2.0],
    linewidths=[0.3, 0.5, 0.7, 1.0, 1.3, 1.6],
    colors="cyan",
    transform=ax2.get_transform(whst),
)
lon, lat = ax2.coords
lat.set_ticklabel_visible(False)
lon.set_axislabel(" ")


fig.tight_layout(rect=[0.1, 0.08, 1.0, 1.0])
fig.savefig(figfile, dpi=200)
print(figfile, end="")
