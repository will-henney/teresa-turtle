import os
import sys
import numpy as np
import astropy
from astropy.table import Table, Row, hstack
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from astropy import units as u
from astropy.coordinates import SkyCoord
import seaborn as sns
import turtle_utils
from turtle_utils import (
    slit_profile,
    extract_full_profile_from_pv,
    extract_slit_profile_from_imslit,
    get_orig_folder,
    find_slit_coords,
    subtract_sky_and_trim,
    make_three_plots,
)

try:
    choice = int(sys.argv[1])
except:
    choice = None

restwavs = {'ha': 6562.79, 'nii': 6583.45, 'nii_s': 6548.05}

# Position of star
RA0, Dec0 = 251.122998321, 23.7998586853

saturation = 6e4


sns.set_palette('RdPu_d', 3)

table1 = Table.read('data/ha-slits.tab', format="ascii.tab")
table2 =  Table.read('data/align-ha.tab', format="ascii.tab")
# The align-ha table takes precedence if islit has been modified
table1.remove_column("islit")
# We already have spec in the ha-slits table
table2.remove_column("spec") 
table = hstack([table1, table2], join_type="exact")



# Photometric reference image
photom, = fits.open('data/imslit-ha/imslit-median.fits')
wphot = WCS(photom.header)
turtle_utils.VERBOSE = 1
neighbors = [-2, -1, 1, 2]
for row in table:
    if choice is not None and row["id"] != choice:
        # If we asked for a single spectrum, then skip all others
        continue
    spec_hdu, = fits.open(get_orig_folder(row["run"]) + "/" + row["spec"] + ".fits")
    im_hdu, = fits.open("data/imslit/" + row["imslit"] + "-wcs.fits")
    # Mask out saturated pixels with NaN
    spec_hdu.data[spec_hdu.data > saturation] = np.nan
    # trim the edge or arrays since sometimes the outer pixels contain garbage
    spec_hdu.data = subtract_sky_and_trim(spec_hdu.data, row)
    spec_profile = extract_full_profile_from_pv(
        spec_hdu,
        wavaxis=row["wa"],
        bandwidth=90.0,
        linewavs=restwavs.values())
    imslit_profile = extract_slit_profile_from_imslit(im_hdu.data, row)
    print(row)
    jslit = np.arange(len(spec_profile))
    # jslit0_spec = np.average(jslit, weights=spec_profile)
    # jslit0_imslit = np.average(jslit, weights=imslit_profile)
    # jslit0_spec = np.nanargmax(spec_profile)
    # jslit0_imslit = np.nanargmax(imslit_profile)
    jslit0_spec = row["j0_s"]
    jslit0_imslit = row["j0_i"]
    print(jslit0_spec, jslit0_imslit, 'shift =', row["shift"])
    slit_coords = find_slit_coords(row, im_hdu.header, spec_hdu.header)
    calib_profile = slit_profile(slit_coords['RA'], slit_coords['Dec'],
                                 photom.data, wphot)


    # Look at neighboring slit positions
    nb_calib_profiles = {}
    for nb in neighbors:
        nbrow = Table(row)[0]   # This is the trick to get a copy of the row
        nbrow["islit"] += nb
        nb_slit_coords = find_slit_coords(nbrow, im_hdu.header, spec_hdu.header)
        nb_calib_profiles[nb] = slit_profile(
            nb_slit_coords['RA'], nb_slit_coords['Dec'], photom.data, wphot)


    # Offset in arcsec along the slit
    slit_points = (np.arange(len(spec_profile)) - jslit0_spec)*slit_coords["ds"]
    # Extra correction for optical halos that show up at +/- 40 arcsec
    halo_mask = np.abs(np.abs(slit_points) - 40.0) < 5.0
    spec_profile -= np.median(spec_profile[halo_mask])

    # Take a window about profile peak to normalize spec_profile
    jslice0 = slice(jslit0_spec-20, jslit0_spec+20)
    # propagate saturated pixels to the calibration profile
    calib_profile_nan = calib_profile.copy()
    calib_profile_nan[~np.isfinite(spec_profile)] = np.nan
    rat0 = np.nansum(spec_profile[jslice0])/np.nansum(calib_profile_nan[jslice0])
    print('Coarse calibration: ratio =', rat0)
    spec_profile /= rat0


    # Make a figure comparing the profiles
    plt_prefix = f"figs/{row.index:03d}-calib"
    ratio = make_three_plots(spec_profile, calib_profile, plt_prefix,
                             slit_points=slit_points,
                             neighbors=nb_calib_profiles, db=row, sdb=slit_coords)
