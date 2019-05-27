tab=[["spec", "imslit", "run", "t", "wa"], ["obj1003", "obj1002_bcr", "1998-06", 60, 1], ["obj1007", "obj1006_bcr", "1998-06", 300, 1], ["obj1010", "obj1009_bcr", "1998-06", 1200, 1], ["obj1015", "obj1012_bcr", "1998-06", 1200, 1], ["obj1018", "obj1017_bcr", "1998-06", 1200, 1], ["spm112", "spm110_u", "2003-06", 1800, 1], ["spm116", "spm115_bcr", "2003-06", 1800, 1], ["spm119", "spm118_bcr", "2003-06", 1800, 1], ["spm122", "spm121_bcr", "2003-06", 1800, 1], ["spm125", "spm124_bcr", "2003-06", 1800, 1], ["spm130", "spm127_bcr", "2003-06", 1800, 1], ["spm135", "spm134_bcr", "2003-06", 1800, 1], ["spm224", "spm223_bcr", "2003-10", 1800, 1], ["spm230", "spm229_bcr", "2003-10", 1800, 1], ["spm329", "spm328_bcr", "2003-10", 1800, 1], ["spm112-2", "spm111_cr", "2011-05", 1800, 2], ["spm018", "spm017_b", "2013-07", 1800, 1], ["spm020", "spm019_b", "2013-07", 1800, 1], ["spm023", "spm021_b", "2013-07", 1800, 1], ["spm0031", "spm0030o_b", "2015-08", 900, 2], ["spm0037", "spm0036o_b", "2015-08", 900, 2], ["spm0043", "spm0042o_b", "2015-08", "?", 2], ["spm0049", "spm0048o_b", "2015-08", 900, 2], ["spm0057", "spm0056o_b", "2015-08", 900, 2], ["spm0110", "spm0109o_b", "2015-08", 1200, 2], ["spm0116", "spm0115o_b", "2015-08", 1800, 2], ["spm0124", "spm0123o_b", "2015-08", 1200, 2], ["spm0173", "spm0172o_b", "2015-08", 900, 2], ["spm0179", "spm0178o_b", "2015-08", 1800, 2], ["spm0186", "spm0185o_b", "2015-08", 1800, 2]]
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
