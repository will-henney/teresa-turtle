import glob
import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
import astropy.units as u
from astropy.coordinates import SkyCoord
sys.path.append('/Users/will/Dropbox/OrionWest')
from helio_utils import waves2vels

# Center on central star of NGC 6210
RA0, Dec0 = 251.122998321, 23.7998586853
c0 = SkyCoord(RA0, Dec0, unit="deg")

wav0_dict = {
    'oiii': 5006.84,
    'ha': 6562.79,
    'nii': 6583.45,
    'nii_s': 6548.05,
    'heii': 6560.10,
    'cii': 6578.15,
    'hei': 5015.68,
}

# Which spectrum file do we *really* need to look in?
parent = {
    'oiii': 'oiii',
    'ha': 'ha',
    'nii': 'nii',
    'nii_s': 'nii_s',
    'heii': 'ha',
    'cii': 'nii',
    'hei': 'oiii',
}

try:
    lineid = sys.argv[1]
except:
    sys.exit(f"Usage: {sys.argv[0]} {{{'|'.join(wav0_dict)}}}")

datadir = "pvextract-oiii" if lineid in ["oiii", "hei"] else "pvextract"

speclist = glob.glob(f'data/{datadir}/*-{parent[lineid]}.fits')
# Rest wavelength in meters
wav0 = 1e-10*wav0_dict[lineid]

for fn in speclist:
    print('Processing', fn)
    spechdu, = fits.open(fn)
    wspec = WCS(spechdu.header, key='A')

    # Eliminate degenerate 3rd dimension from data array and trim off bad bits
    spec2d = spechdu.data[0]

    # Convert to heliocentric velocity
    [[wav1, _, _], [wav2, _, _]] = wspec.all_pix2world([[0, 0, 0], [1, 0, 0]], 0)
    [v1, v2] = waves2vels(np.array([wav1, wav2]), wav0, spechdu.header, usewcs="A")


    # sequence of pixels along the slit spatial axis
    ipixels = np.arange(wspec.array_shape[1])
    # sky coordinate that corresponds to each pixel
    coords = pixel_to_skycoord(ipixels, 0, wcs=wspec)
    # separation of each pixel from star
    radii = c0.separation(coords).to(u.arcsec).value
    # reference pixel has minimum separation from star
    iref = radii.argmin()
    cref = coords[iref]
    # slit offset is minimum of radii
    offset = radii.min()
    # but sign depends on PA
    pa_offset = c0.position_angle(cref)
    if np.sin(pa_offset) < 0.0:
        offset *= -1.0

    # separation from reference pixel gives coordinate along slit
    s = cref.separation(coords).to(u.arcsec).value
    # inter-pixel separation
    ds = np.abs(np.diff(s)).mean()
    # PA of slit
    pa = cref.position_angle(coords[-1])
    if np.cos(pa) < 0.0:
        # Make sure positive offset is to N
        ds *= -1.0
        # And flip slit PA to compensate
        pa += np.pi*u.rad

    wnew = WCS(naxis=2)
    wnew.wcs.ctype = ['VHEL', 'LINEAR']
    wnew.wcs.crpix = [1, iref+1]
    wnew.wcs.cunit = ['km/s', 'arcsec']
    wnew.wcs.crval = [v1.to('km/s').value, 0.0]
    wnew.wcs.cdelt = [(v2 - v1).to('km/s').value, ds]

    newsuffix = f"-PA{int(pa.to(u.deg).value)%360:03d}-sep{int(offset):+04d}"    
    newhdr = wnew.to_header()
    newhdr["PA"] = pa.to(u.deg).value, "Position angle of slit, degrees"
    newhdr["OFFSET"] = offset, "Perpendicular offset of slit from star, arcsec"
    newhdr["WEIGHT"] = spechdu.header["WEIGHT"]
    new_fn = fn.replace(f'data/{datadir}/', 'data/PVoffset/')
    new_fn = new_fn.replace(parent[lineid], lineid)
    new_fn = new_fn.replace('.fits', newsuffix + '.fits')
    print('Writing', new_fn)
    fits.PrimaryHDU(data=spec2d, header=newhdr).writeto(new_fn, overwrite=True)
