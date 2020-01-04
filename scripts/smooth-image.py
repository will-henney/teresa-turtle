import numpy as np
import argh
from astropy.io import fits
from astropy.convolution import convolve_fft, Gaussian2DKernel

def main(infile, width=16, twopass=False, threshold=1.5):
    hdu = fits.open(infile)[0]
    im = hdu.data
    gauss_kernel = Gaussian2DKernel(width)
    smoothim = convolve_fft(im, gauss_kernel)
    sharpim = im/smoothim

    if twopass:
        mask = (sharpim > threshold) | (im == 0.0)
        im[mask] = np.nan
        print('Starting second pass: N(masked) =', mask.sum())
        smoothim = convolve_fft(im, gauss_kernel, normalize_kernel=True)
        sharpim = im/smoothim

    outhdu = fits.PrimaryHDU(data=smoothim, header=hdu.header)
    outfile = infile.replace(".fits", f"_smooth_{width}.fits")
    outhdu.writeto(outfile, overwrite=True)

    outhdu.data = sharpim
    outhdu.writeto(outfile.replace("_smooth_", "_sharp_"), overwrite=True)

if __name__ == "__main__":
    argh.dispatch_command(main)
