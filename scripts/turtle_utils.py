import numpy as np
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import pixel_to_skycoord
from astropy import units as u
from astropy.modeling import models, fitting
from matplotlib import pyplot as plt
import seaborn as sns

VERBOSE = 0

ORIG_DATA_ROOT = "../Papers/LL-Objects/NGC6210"

def get_orig_folder(temporada):
    if "2015" in temporada:
        return ORIG_DATA_ROOT + "/Temporada2015"
    else:
        return  ORIG_DATA_ROOT + "/Varias-temporadas"

def slit_profile(ra, dec, image, wcs):
    """
    Find the image intensity for a list of positions (ra and dec)
    """
    xi, yj = wcs.all_world2pix(ra, dec, 0)
    # Find nearest integer pixel
    ii, jj = np.floor(xi + 0.5).astype(int), np.floor(yj + 0.5).astype(int)
    if VERBOSE > 0:
        print(ra[::100], dec[::100])
        print(ii[::100], jj[::100])
    ny, nx = image.shape
    return np.array([image[j, i]
                     if (0 < i < nx and 0 < j < ny)
                     else np.nan
                     for i, j in list(zip(ii, jj))])

def find_slit_coords(db, hdr, shdr):
    """Find the coordinates of all the pixels along a spectrograph slit

    Input arguments are a dict-like 'db' of hand-measured values (must
    contain 'wa', 'ij', islit' and 'shift') and a FITS headers 'hdr' from
    the image+slit exposure and 'shdr' from a spectrum exposure

    Returns a dict of 'ds' (slit pixel scale), 'PA' (slit position
    angle), 'RA' (array of RA values in degrees along slit), 'Dec'
    (array of Dec values in degrees along slit)

    """

    # Decide on axis order for both spectrum and image. Note that
    # values of 'wa' and 'ij' give the axis that is perpendicular to
    # the slit length (wavelength or position, respectively). Hence we
    # subtract from 3 to get the slit length axis
    jstring_i = str(3 - db['ij'])  # which image (I+S) axis lies along slit
    jstring_s = str(3 - db['wa'])  # which spec (pv) axis lies along slit

    dRA_arcsec = hdr['CD1_'+jstring_i]*3600*np.cos(np.radians(hdr['CRVAL2']))
    dDEC_arcsec = hdr['CD2_'+jstring_i]*3600
    ds = np.hypot(dRA_arcsec, dDEC_arcsec)
    PA = np.degrees(np.arctan2(dRA_arcsec, dDEC_arcsec)) % 360.0

    # Deal with parameters that depend on orientation of the PV image
    if jstring_s == '1':
        # PV slit has spatial axis horizontal in IMAGE coords
        ns = shdr['NAXIS1']
        # Mezcal has used two different ways of specifying on-chip binning
        try:
            # Older way
            spec_binning = shdr['CBIN']
        except KeyError:
            try:
                # Newer way
                spec_binning = shdr['CCDXBIN']
            except KeyError:
                # And the very old data don't have it at all
                spec_binning = 1
    elif jstring_s == '2':
        # PV slit has spatial axis vertical in IMAGE coords
        ns = shdr['NAXIS2']
        try:
            spec_binning = shdr['RBIN']
        except KeyError:
            try:
                spec_binning = shdr['CCDYBIN']
            except KeyError:
                spec_binning = 1
    else:
        raise ValueError('PV slit axis (3 - wa) must be 1 or 2')

    # Pixel coords of each slit pixel on image (in 0-based convention)
    # Deal with parameters that depend on orientation of the I+S image
    if jstring_i == '1':
        # Slit is horizontal in IMAGE coords
        iarr = np.arange(ns) - float(db['shift'])
        jarr = np.ones(ns)*float(db['islit'])
        # Mezcal has used two different ways of specifying on-chip binning
        try:
            # Older way
            image_binning = hdr['CBIN']
        except KeyError:
            try:
                # Newer way
                image_binning = hdr['CCDXBIN']
            except KeyError:
                # And the very old data don't have it at all
                image_binning = 1
        # correct for difference in binning between the image+slit and the spectrum
        iarr *= spec_binning/image_binning
    elif jstring_i == '2':
        # Slit is vertical in IMAGE coords
        iarr = np.ones(ns)*float(db['islit'])
        jarr = np.arange(ns) - float(db['shift'])
        try:
            image_binning = hdr['RBIN']
        except KeyError:
            try:
                image_binning = hdr['CCDYBIN']
            except KeyError:
                image_binning = 1
        jarr *= spec_binning/image_binning
    else:
        raise ValueError('I+S slit axis (3 - ij) must be 1 or 2')

    if db['s'] < 0:
        # Slit pixel axis has opposite sense in I+S and PV
        iarr = iarr[::-1]
        jarr = jarr[::-1]

    print('iarr =', iarr[::100], 'jarr =', jarr[::100])
    # Also correct the nominal slit plate scale
    ds *= spec_binning/image_binning

    # Convert to world coords, using the native frame
    w = WCS(hdr)
    observed_frame = w.wcs.radesys.lower()
    # Note it is vital to ensure the pix2world transformation returns
    # values in the order (RA, Dec), even if the image+slit may have
    # Dec first
    coords = SkyCoord(*w.all_pix2world(iarr, jarr, 0, ra_dec_order=True),
                      unit=(u.deg, u.deg), frame=observed_frame)
    print('coords =', coords[::100])
    print('Binning along slit: image =', image_binning, 'spectrum =', spec_binning)
    # Make sure to return the coords in the ICRS frame
    return {'ds': ds, 'PA': PA,
            'RA': coords.icrs.ra.value,
            'Dec': coords.icrs.dec.value}

def subtract_sky_and_trim(data, db, trim=3, margin=10):
    """Assume that pixels within `trim` of edge might be bad, and use
    average sky within margin of edge in spatial direction to define
    the bg
    """
    # convert axis notation from FITS to python convention
    wav_axis = 2 - db["wa"]
    if wav_axis == 0:
        bg = 0.5*(data[:, trim:margin] +
                  data[:, -margin:-trim]).mean(axis=1, keepdims=True)
    else:
        bg = 0.5*(data[trim:margin, :] +
                  data[-margin:-trim, :]).mean(axis=0, keepdims=True)
    # Remove sky background
    newdata = data - bg
    # And only then can we clean up the trim zone
    newdata[:trim, :] = 0.0
    newdata[-trim:, :] = 0.0
    newdata[:, :trim] = 0.0
    newdata[:, -trim:] = 0.0
    return newdata

def extract_full_profile_from_pv(spec_hdu, wavaxis, bandwidth, linewavs):
    assert(wavaxis in [1, 2]) # wavaxis is in FITS convention
    w = WCS(spec_hdu.header)
    if wavaxis == 1:
        nwav = spec_hdu.header['NAXIS1']
        im = spec_hdu.data[:, :]
        wavs, _ = w.all_pix2world(np.arange(nwav), [0], 0)
    else:
        nwav = spec_hdu.header['NAXIS2']
        im = spec_hdu.data[:, :].T
        _, wavs = w.all_pix2world([0], np.arange(nwav), 0)

    # im should have wavelength as last axis (python convention)
    assert(nwav == im.shape[-1])

    wavmask = np.ones((nwav,)).astype(bool)
    # remove from continuum mask +/- 150 km/s around each line
    for wav0 in linewavs:
        vels = 3e5*(wavs - wav0)/wav0
        wavmask = wavmask & (np.abs(vels) > 150.0)
    # broadcast to 2 dimensions
    imwts = np.ones_like(im)*wavmask[None, :]
    av_cont_profile = np.average(im, weights=imwts, axis=-1)
    # find how much extra continuum to add
    dwav = abs(wavs[1] - wavs[0])
    pv_bw = abs(wavs[-1] - wavs[0])
    missing_cont_profile = av_cont_profile*(bandwidth - pv_bw)/dwav
    # Add to the profile summed over the PV bandwidth
    full_profile = im.sum(axis=-1) + missing_cont_profile
    return full_profile


def extract_slit_profile_from_imslit(data, db, slit_width=1):
    print(db["islit"])
    i1, i2 = int(db["islit"]) - slit_width, int(db["islit"]) + slit_width
    if db["ij"] == 1:
        return data[:, i1:i2].sum(axis=1)
    elif db["ij"] == 2:
        return data[i1:i2, :].sum(axis=0 )
    else:
        raise ValueError("ij must be 1 or 2")

def fit_cheb(x, y, npoly=3, mask=None):
    """Fits a Chebyshev poly to y(x) and returns fitted y-values"""
    fitter = fitting.LinearLSQFitter()
    p_init = models.Chebyshev1D(npoly, domain=[x.min(), x.max()])
    if mask is None:
        mask = np.ones_like(x).astype(bool)
    p = fitter(p_init, x[mask], y[mask])
    if VERBOSE > 0:
        print(p)
    return p(x)

def make_three_plots(spec, calib, prefix,
                     slit_points=None, niirat=None, neighbors=None, db=None, sdb=None):
    assert spec.shape == calib.shape
    fig, axes = plt.subplots(3, 1)

    if slit_points is None:
        ypix = np.arange(len(calib))
        xlabel = "Slit pixel"
        xlim = None
    else:
        ypix = slit_points
        xlabel = "Slit position, arcsec"
        xlim = -80, 80

    xlim = xlim or (ypix.min(), ypix.max())

    # vmax = np.percentile(calib, 95) + 2*calib.std()
    vmax = 20.0
    vmin = -0.01
    ratio = spec/calib

    alpha = 0.8

    # First, plot two profiles against each other to check for zero-point offsets
    # axes[0].plot(calib, spec/ratio_fit, '.', alpha=alpha)
    axes[0].plot(calib, spec, '.', alpha=alpha)
    axes[0].plot([vmin, vmax], [vmin, vmax], '-', alpha=alpha)
    axes[0].set_xlim(vmin, vmax)
    axes[0].set_ylim(vmin, vmax)
    axes[0].set_xlabel('Calibration Image')
    axes[0].set_ylabel('Uncorrected Integrated Spectrum')
    axes[0].set_xscale('symlog', linthreshx=0.01)
    axes[0].set_yscale('symlog', linthreshy=0.01)

    # Second, plot each against slit pixel to check spatial offset
    axes[1].plot(ypix, spec, alpha=alpha, lw=1,
                 label='Integrated Spectrum')
    axes[2].plot(ypix, spec/np.nanmax(calib), alpha=alpha, lw=1,
                 label='Integrated Spectrum')
    axes[1].plot(ypix, calib, alpha=alpha, label='Calibration Image')
    axes[2].plot(ypix, calib/np.nanmax(calib), alpha=alpha, label='Calibration Image')
    if neighbors is not None:
        for nb, calib_nb in neighbors.items():
            lw = 0.5 + 0.1*nb
            label = f"Slit $\Delta x = {nb:+1d}$"
            axes[1].plot(ypix, calib_nb,
                         alpha=0.3*alpha, lw=lw, color="k", label=label)
            axes[2].plot(ypix, calib_nb/np.nanmax(calib_nb),
                         alpha=0.3*alpha, lw=lw, color="k", label=label)
    # axes[1].plot(ypix, spec/ratio_fit, alpha=alpha, lw=1.0,
    #              label='Corrected Integrated Spectrum')
    axes[1].set_xlim(*xlim)
    axes[1].set_ylim(vmin, vmax)
    axes[1].legend(fontsize='xx-small', loc='upper right')
    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel('Profile (absolute log scale)')
    axes[1].set_yscale('symlog', linthreshy=0.01)

    # # Third, plot ratio to look for spatial trends
    # axes[2].plot(ypix, ratio, alpha=alpha)
    # axes[2].plot(ypix, ratio_fit, alpha=alpha)
    # if niirat is not None:
    #     axes[2].plot(ypix, niirat, 'b', lw=0.5, alpha=0.5)
    axes[2].set_xlim(-40, 40)
    axes[2].set_ylim(-0.05, 1.05)
    axes[2].set_xlabel(xlabel)
    # axes[2].set_ylabel('Ratio: Spec / Calib')
    axes[2].set_ylabel('Profile (relative linear scale)')

    info = ""
    if db is not None:
        # Add some info to the graphs
        info += fr"H$\alpha$ slit {db['id']:02d}" + "\n"
        info += f"PV: {db['spec']}" +"\n"
        info += f"I+S: {db['imslit']}" + "\n"
        info += f"Date: {db['run']}, t = {db['t']} s" + "\n"
    if sdb is not None:
        info += fr"Slit PA = ${sdb['PA']:.1f}^\circ$, ds = {sdb['ds']:.2f} arcsec" + "\n"
    if info:
        axes[0].text(0.95, 0.05, info,
                     fontsize="small",
                     ha="right", va="bottom", transform=axes[0].transAxes)
    fig.set_size_inches(5, 8)
    fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    fig.savefig(prefix+'.png', dpi=300)

    return None
