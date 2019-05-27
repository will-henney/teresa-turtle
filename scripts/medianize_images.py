import sys
import os
import glob
import numpy as np
from astropy.io import fits

try:
    datapath = sys.argv[1]
    fnlist = glob.glob(f"{datapath}/imslit-??.fits")
except:
    sys.exit(f"Usage: {sys.argv[0]} DATAPATH")

imlist = []
for fitsname in fnlist:
    hdu, = fits.open(fitsname)
    imlist.append(hdu.data)
imstack = np.dstack(imlist)
median = np.median(imstack, axis=-1)
fits.PrimaryHDU(header=hdu.header,
                data=median).writeto(f'{datapath}/imslit-median.fits', overwrite=True)

ratcombo = np.zeros_like(median)
combo = np.zeros_like(median)
for im, fn in zip(imlist, fnlist):
    combo = combo + im
    head, tail = os.path.split(fn)
    outname = os.path.join(head, tail.replace('imslit', 'imslit-ratio'))
    ratio = im/median
    ratcombo = ratcombo + ratio
    fits.PrimaryHDU(header=hdu.header,
                    data=ratio).writeto(outname, overwrite=True)
fits.PrimaryHDU(header=hdu.header,
                data=ratcombo).writeto(f'{datapath}/imslit-ratcombo.fits',
                                       overwrite=True)
fits.PrimaryHDU(header=hdu.header,
                data=combo).writeto(f'{datapath}/imslit-combo.fits',
                                    overwrite=True)
