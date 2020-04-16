import sys
import numpy as np
from astropy.io import fits
from pathlib import Path

datadir = Path("~/Dropbox/Papers/LL-Objects/NGC6210/HST").expanduser()
fitsfiles = tuple(datadir.glob("*_drz.fits"))

try:
    FILTER, YEAR = sys.argv[1], sys.argv[2]
except:
    sys.exit(f"Usage: {sys.argv[0]} FILTER YEAR")

# Find all observations that are for the requested filter and year
obslist = []
for fitsfile in fitsfiles:
    hdulist = fits.open(fitsfile)
    try:
        h = hdulist[0].header
        hi = hdulist["SCI"].header
        if FILTER in hi["PHOTMODE"] and h["DATE-OBS"].startswith(YEAR):
           obslist.append(hdulist) 
    except:
        pass

print(f"Found {len(obslist)} files with FILTER={FILTER} YEAR={YEAR}")

outfile = f"turtle-{FILTER}-{YEAR}.fits"
outdir = Path("../data/hst-cr-combine")


data = None
THRESH = 1.5
# supposedly already in counts/s
image_stack = np.stack(
    # [_["SCI"].data/float(_[0].header["EXPTIME"]) for _ in obslist]
    [_["SCI"].data for _ in obslist]
)
times = np.array([float(_[0].header["EXPTIME"]) for _ in obslist])
times = times[:, None, None] * np.ones_like(image_stack)
if len(obslist) >= 3:
    # Easy case - median masking for CRs
    median_image = np.median(image_stack, axis=0, keepdims=True)
    cr_mask = (image_stack > 0.0) & (image_stack > THRESH*median_image)
    times[cr_mask] = 0.0
    image = np.average(image_stack, axis=0, weights=times)
else:
    # Only two images - take the minimum and hope!
    image = np.nanmin(image_stack, axis=0)

obslist[0]["SCI"].data = image
obslist[0].writeto(outdir / outfile, overwrite=True)
