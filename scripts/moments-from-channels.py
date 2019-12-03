import sys
import numpy as np
from astropy.io import fits

try:
    lineid = sys.argv[1]
except:
    sys.exit(f"Usage {sys.argv[0]} LINEID")


prefix = f"../data/maps/turtle-slits-{lineid}"
vels = np.arange(-75.0, +15.0, 10.0)
nv = len(vels)

# make a cube from the stack of isovel images
imstack = []
for vel in vels:
    vstring = f"{int(vel-5.0):+04d}{int(vel+5.0):+04d}"
    filename = f"{prefix}{vstring}.fits"
    hdu = fits.open(filename)["SCALED"]
    imstack.append(hdu.data)

imcube = np.stack(imstack, axis=0)
vcube = vels[:, None, None]*np.ones_like(imcube)
vmean = np.average(vcube, weights=imcube, axis=0)
sigsq = np.average((vcube - vmean)**2, weights=imcube, axis=0)
sigma = np.sqrt(sigsq)

fits.PrimaryHDU(header=hdu.header, data=vmean).writeto(
    f"{prefix}-vmean.fits", overwrite=True,
)
fits.PrimaryHDU(header=hdu.header, data=sigma).writeto(
    f"{prefix}-sigma.fits", overwrite=True,
)
