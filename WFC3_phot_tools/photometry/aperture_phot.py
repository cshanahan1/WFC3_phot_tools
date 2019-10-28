from astropy.stats import mad_std
from ginga.util import zscale
from photutils import DAOStarFinder
from photutils import aperture_photometry, CircularAperture, CircularAnnulus
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import DAOStarFinder, aperture_photometry,  deblend_sources, detect_sources, detect_threshold, source_properties

import warnings
warnings.filterwarnings("ignore")

def show_source_detection_plot(data, coo_tab):
    """
    Displays a plot of `data` (z-scaled), with sources in `coo_tab` overplotted.

    Parameters
    ----------
    data : array
        Image array.
    coo_tab : `astropy.table.Table`
        Table with source positions as columns 'xcentroid' and 'ycentroid'.

    Returns
    -------
    Window with plot.

    """

    z1, z2 = zscale.zscale(data)
    
    plt.imshow(data, origin='lower', cmap='Greys_r', vmin=z1, vmax=z2)
    if coo_tab:
        plt.scatter(coo_tab['xcentroid'], coo_tab['ycentroid'], c='r')
        plt.title('{} Sources Detected'.format(len(coo_tab)))
    else:
        plt.title('No Sources Detected')
    plt.show()

def detect_sources_segmap(data, threshold, npixels, kernel_fwhm = 1.8, 
                           bkgrnd_threshold=True, show=False):
    """
    Runs image segmentation to detect sources in `data`. 

    Parameters
    ----------
    data : array
        Image array.
    threshold : float or array
        Detection threshold value, or pixel-wise threshold image (must be same
        shape as `data`.)
    npixels : int
        Positive integer number of connected pixels, each greater that 
        `threshold` that an object must have to be detected.
    kernel_fwhm : float
        FWHM of gaussian kernel used to smooth image before segmentation. 

    Returns
    -------
    coo_tab : `astropy.table.Table` or int
        Table with detected source(s). Returns '0' if no sources are detected.

    """

    sigma = kernel_fwhm * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(data, threshold=threshold, npixels=npixels, 
                          filter_kernel=kernel)
    if segm:
        coo_tab = source_properties(data, segm).to_table()
        if show:
            show_source_detection_plot(data, coo_tab)
        return coo_tab

    if not segm:
        if show:
            show_source_detection_plot(data, None)
        return 0

def calc_sky_annulus(data, x, y, r_in, r_out, sky_method='median', 
                     sigma_clip=True):
    """ 
        Calculates the background level in a circular annulus centered
        at `x, y` with an inner radius of `r_in` pixels and an outer radius of
        `r_out` pixels. 

        Sky level and standard deviation in annulus are returned. 
        Options for `sky_method` are `mean`, `median`, or `mode`. 
        Uses a sigma clip with sigma=`n_sigma_clip` to compute stats. 
        If n_sigma_clip is set to 0, sigma clipping won't be done. 

        Parameters
        ----------
            data : array
                Science data array.
            x, y : float
                Pixel position for annulus center.
            r_in, r_out : int
                Inner, outer circular annulus radii
            sky_method : str
                'mean', 'median', or 'mode' 
            n_sigma_clip : int
                Used for sigma clipping. 0 to turn off sigma clipping

        Returns
        -------
            (back, backstd) : tuple
                Sky level in annulus, std of sky.

    """
    from .background_median import aperture_stats_tbl

    sky_apertures = CircularAnnulus((x, y), r_in, r_out)

    sky_tbl = aperture_stats_tbl(data, sky_apertures, sigma_clip=sigma_clip)

    return(sky_tbl['aperture_'+sky_method].item(), sky_tbl['aperture_std'].item())

def compute_phot_err_daophot(flux, back, back_rms, phot_ap_area,
                             sky_ap_area, gain=1.0):

    """Calculates flux errors in the same manner as IRAF/DAOphot. 

    The error terms in this model represent Poisson noise from the source, 
    Poisson noise in the sky and readout noise, and error in the sky measurement
    err = sqrt((flux - back*phot_ap_area)/epadu + phot_ap_area * backrms**2 +
    phot_ap_area**2 * backrms**2 / sky_ap_area)
    Where flux is the aperture sum, back is the per-pixel background level,
    epadu is the conversion factor between e- and adu (gain), phot_ap_area is
    the area of the photometric aperture, backrms is the uncertainty in the sky,
    and sky_ap_area is the sky aperture area. Note that the flux/background
    in the above equation are in ADU, but the input to this function is in e-.
    The conversion is done internally.
    Parameters
    ----------
    flux : float
        Non-sky subtracted flux of star, in electrons.
    back : float
        Per-pixel background level.
    back_rms : float
        Background RMS.
    phot_ap_area : int or float
        Area of photometric aperture.
    sky_ap_area : int or float
        Area of sky aperture.
    gain : float
        CCD gain, default 1.0.
     Returns
     --------
     errs : tuple
        (error in insturmental magnitudes, error in flux)
    """

    #convert values input in e- to ADU, as equation expects
    #sky subtract flux to isolate poisson noise from source
    flux = (flux - back*phot_ap_area) / gain
    back = back / gain
    back_rms = back_rms / gain

    err1 = (flux/gain)
    err2 = (phot_ap_area*back_rms**2)
    err3 = (phot_ap_area**2 * back_rms**2) / sky_ap_area

    flux_err_adu = np.sqrt(np.abs(err1 + err2 + err3)) # in ADU
    flux_err = flux_err_adu * gain

    return flux_err


def circular_aperture_photometry(data, xc, yc, aperture_radii):
    """
    Aperture photometry with circular apertures on single source in `data` 
    located at (xc, yc).

    Returns a table of aperture sums for every aperture size in 
    `aperture_radii`.

    Parameters
    ----------
    data : array 
        Image array.
    xc, yc : float
        x, y location of source in `data`.
    aperture_radii : list of floats
        Desired aperture radii, in pixels.

    Returns
    -------
    phot_table : `astropy.table.Table`
        Table with columns for x&y source position, and sum in every circular
        aperture in `aperture_radii`. 

    """
    aps = [CircularAperture((xc, yc), r=rad) for rad in aperture_radii]
    phot_table = Table(aperture_photometry(data, aps))
    phot_table.remove_column('id')

    table_order = []
    for i, rad in enumerate(aperture_radii):
        rad = str(rad)
        phot_table.rename_column('aperture_sum_{}'.format(str(i)), \
                                 'countrate_{}'.format(rad))
        table_order.append('countrate_{}'.format(rad))

    table_order = ['xcenter','ycenter'] + table_order

    phot_table = phot_table[table_order]

    return(phot_table)
