tab=[["spec", "imslit", "run", "t", "wa", "islit", "ij"], ["spm244", "spm242_bcr", "2004-06", 1800, 1, 275, 1], ["spm253", "spm251_bcr", "2004-06", 1800, 1, 280, 1], ["spm258", "spm256_bcr", "2004-06", 1800, 1, 285, 1], ["spm294", "spm293_bcr", "2004-06", 1800, 1, 262, 1], ["spm297", "spm296_bcr", "2004-06", 1800, 1, 265, 1], ["spm114", "spm115_bcr", "2011-05", 600, 2, 263.5, 1], ["spm0033", "spm0035o_b", "2015-08", 900, 2, 485, 2], ["spm0039", "spm0041o_b", "2015-08", 900, 2, 482, 2], ["spm0045", "spm0047o_b", "2015-08", 900, 2, 479, 2], ["spm0051", "spm0053o_b", "2015-08", 900, 2, 476, 2], ["spm0059", "spm0061o_b", "2015-08", 900, 2, 475, 2], ["spm0112", "spm0114o_b", "2015-08", 1200, 2, 486, 2], ["spm0118", "spm0120o_b", "2015-08", 1800, 2, 479.5, 2], ["spm0126", "spm0128o_b", "2015-08", 1200, 2, 476.5, 2], ["spm0175", "spm0177o_b", "2015-08", 900, 2, 489.5, 2], ["spm0181", "spm0183o_b", "2015-08", 1800, 2, 482, 2], ["spm0188", "spm0190o_b", "2015-08", 1800, 2, 512, 2]]
import numpy as np
from scipy.interpolate import griddata
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table


#
# First set up WCS for the output image
# We use capital letters for the output variables
#

NX, NY = 512, 512
# 0.5 arcsec pixels
dRA, dDec = -0.3/3600., 0.3/3600.
# Center on central star of NGC 6210
RA0, Dec0 = 251.122998321, 23.7998586853
W = WCS(naxis=2)
W.wcs.cdelt = [dRA, dDec]
W.wcs.crpix = [0.5*(1 + NX), 0.5*(1 + NY)]
W.wcs.crval = [RA0, Dec0]
W.wcs.ctype = ['RA---TAN', 'DEC--TAN']

outimage = np.zeros((NY, NX))
# Create world coord arrays for output image
II, JJ = np.meshgrid(np.arange(NX), np.arange(NY))
RA, Dec = W.all_pix2world(II, JJ, 0)

#
# Read in the list of slits
#
table = Table(rows=tab[1:], names=tab[0])

for i, row in enumerate(table):
    hdu, = fits.open(f"data/imslit/{row['imslit']}-wcs.fits")
    # image = (hdu.data - row['bias']) / (row['core'] - row['bias'])
    image = hdu.data
    outfilename = f'data/imslit-oiii/imslit-{i:02d}.fits'
    ny, nx = image.shape
    #hdu.header.remove('@EPOCH')
    w = WCS(hdu.header)
    # Create world coord arrays for input image
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny))
    ra, dec = w.all_pix2world(ii, jj, 0)
    # Do the interpolation
    points = np.array(list(zip(ra.ravel(), dec.ravel())))
    xi = np.array(list(zip(RA.ravel(), Dec.ravel())))
    outimage = griddata(points, image.ravel(), xi, method='nearest').reshape((NY, NX))
    bg = np.nanmedian(outimage[420:, 420:])
    core = np.nanmedian(outimage[230:290, 230:290])
    print(core, bg)
    outimage = (outimage - bg)/(core - bg)
    # Save the output image
    fits.PrimaryHDU(header=W.to_header(), data=outimage).writeto(outfilename, overwrite=True)
