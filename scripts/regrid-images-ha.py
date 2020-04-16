tab=[["spec", "imslit", "run", "t", "wa", "islit", "ij", "s", "id"], ["obj1003_bcrx", "obj1002_bcr", "1998-06", 60, 1, 390, 1, 1, "o"], ["obj1007_bcrx", "obj1006_bcr", "1998-06", 300, 1, 388, 1, 1, "n"], ["obj1010_bcrx", "obj1009_bcr", "1998-06", 1200, 1, 388, 1, 1, "m"], ["obj1015_bcrx", "obj1012_bcr", "1998-06", 1200, 1, 390.5, 1, 1, "p"], ["obj1018_bcrx", "obj1017_bcr", "1998-06", 1200, 1, 390.5, 1, 1, "q"], ["spm112_bcrx", "spm110_u", "2003-06", 1800, 1, 257, 1, 1, "f"], ["spm116_bcrx", "spm115_bcr", "2003-06", 1800, 1, 263, 1, 1, "e"], ["spm119_bcrx", "spm118_bcr", "2003-06", 1800, 1, 265, 1, 1, "d"], ["spm122_bcrx", "spm121_bcr", "2003-06", 1800, 1, 267, 1, 1, "h"], ["spm125_bcrx", "spm124_bcr", "2003-06", 1800, 1, 269.5, 1, 1, "i"], ["spm130_bcrx", "spm127_bcr", "2003-06", 1800, 1, 272.5, 1, 1, "k"], ["spm135_bcrx", "spm134_bcr", "2003-06", 1800, 1, 277.5, 1, 1, "b"], ["spm224_bcrx", "spm223_bcr", "2003-10", 1800, 1, 294.5, 1, 1, "t"], ["spm230_bcrx", "spm229_bcr", "2003-10", 1800, 1, 294.5, 1, 1, "v"], ["spm329_bcrx", "spm328_bcr", "2003-10", 1800, 1, 273, 1, 1, "w"], ["spm112-2_bcrx", "spm111_cr", "2011-05", 1800, 2, 531, 2, 1, "g"], ["spm018-bcrx", "spm017_b", "2013-07", 1800, 1, 505, 2, -1, "c"], ["spm020-bcrx", "spm019_b", "2013-07", 1800, 1, 483, 2, -1, "i"], ["spm023-bcrx", "spm021_b", "2013-07", 1800, 1, 496, 2, -1, "j"], ["spm0031o_bcrx", "spm0030o_b", "2015-08", 900, 2, 485, 2, 1, "c"], ["spm0037o_bcrx", "spm0036o_b", "2015-08", 900, 2, 484, 2, 1, "d"], ["spm0043o_bcrx", "spm0042o_b", "2015-08", 900, 2, 481, 2, 1, "e"], ["spm0049o_bcrx", "spm0048o_b", "2015-08", 900, 2, 480, 2, 1, "f"], ["spm0057o_bcrx", "spm0056o_b", "2015-08", 900, 2, 476, 2, 1, "g"], ["spm0110o_bcrx", "spm0109o_b", "2015-08", 1200, 2, 490.5, 2, 1, "b"], ["spm0116o_bcrx", "spm0115o_b", "2015-08", 1800, 2, 485.5, 2, 1, "a"], ["spm0124o_bcrx", "spm0123o_b", "2015-08", 1200, 2, 478, 2, 1, "i"], ["spm0173o_bcrx", "spm0172o_b", "2015-08", 900, 2, 493, 2, 1, "h"], ["spm0179o_bcrx", "spm0178o_b", "2015-08", 1800, 2, 489.5, 2, 1, "j"], ["spm0186o_bcrx", "spm0185o_b", "2015-08", 1800, 2, 515.5, 2, 1, "k"]]
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
    outfilename = f'data/imslit-ha/imslit-{i:02d}.fits'
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
