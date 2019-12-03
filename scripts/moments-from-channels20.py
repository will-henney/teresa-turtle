import sys
import numpy as np
from astropy.io import fits

try:
    lineid = sys.argv[1]
    try:
        suffix = sys.argv[2]
    except:
        suffix = ""
except:
    sys.exit(f"Usage {sys.argv[0]} LINEID [SUFFIX]")

mapdir = "maps-oiii-20" if lineid == "oiii" else "maps20"

prefix = f"../data/{mapdir}/turtle-slits-{lineid}"
vels = np.arange(-90.0, +50.0, 20.0)
nv = len(vels)

# make a cube from the stack of isovel images
imstack = []
for vel in vels:
    vstring = f"{int(vel-10.0):+04d}{int(vel+10.0):+04d}"
    filename = f"{prefix}{vstring}{suffix}.fits"
    hdu = fits.open(filename)["SCALED"]
    imstack.append(hdu.data)

imcube = np.stack(imstack, axis=0)
vcube = vels[:, None, None]*np.ones_like(imcube)
vmean = np.sum(vcube*imcube, axis=0)/np.sum(imcube, axis=0)
# vmean = np.average(vcube, weights=imcube, axis=0)
sigsq = np.sum((vcube - vmean)**2*imcube, axis=0)/np.sum(imcube, axis=0)
# sigsq = np.average((vcube - vmean)**2, weights=imcube, axis=0)
sigma = np.sqrt(sigsq)

fits.PrimaryHDU(header=hdu.header, data=vmean).writeto(
    f"{prefix}-vmean{suffix}.fits", overwrite=True,
)
fits.PrimaryHDU(header=hdu.header, data=sigma).writeto(
    f"{prefix}-sigma{suffix}.fits", overwrite=True,
)
