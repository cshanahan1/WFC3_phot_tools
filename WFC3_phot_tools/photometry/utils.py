""" Contains utility functions for photometry, including PAM correction. 

	Authors
	-------
	Clare Shanahan, Oct 2019

"""
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.io import fits
import numpy as np

def query_simbad_by_targname(targname):

	"""Querys SIMBAD for `targname`. Constructs a `astropy.coordinates.SkyCoord`
	object containing the RA, Dec, proper motion, radial velocity, and distance
	to the object.
	
	Paramters
	---------
		targname : str
			Target name. Simbad is queried by target name.
	Returns
	-------
		c : `astropy.coordinates.SkyCoord`
			SkyCoord object with query results.
		"""

	from astroquery.simbad import Simbad
	
	customSimbad = Simbad()
	customSimbad.add_votable_fields('pmra', 'pmdec', 'distance', 'rv_value')
	result_table = customSimbad.query_object(targname)

	RA = str(result_table['RA'].item())
	RA = RA[0:2]+'h'+RA[3:5]+'m'+RA[6:]+'s'
	
	DEC = str(result_table['DEC'].item())
	DEC = DEC[0:3]+'d'+DEC[4:6]+'m'+DEC[7:]+'s'
	
	PM_RA = result_table['PMRA'].item()*u.mas/u.yr
	PM_DEC = result_table['PMDEC'].item()*u.mas/u.yr

	c = SkyCoord(RA, DEC, frame='icrs', obstime=Time('J2000'),
				 distance=result_table['Distance_distance'].item()*u.pc,
				 pm_ra_cosdec=PM_RA,
				 pm_dec=PM_DEC)
	
	return(c)


def apply_proper_motion_targ(targname, mjd):
	
	"""Querys Simbad for ra, dec, and proper motions for 'targname'.
	   Returns (RA, and Dec) at date `mjd` considering the proper motion.

		Parameters
		---------
		targname : str
			Target name. Simbad is queried by target name. 
		mjd : float
			MJD for which to calculate RA, and Dec. 

		Returns
		-------
		(ra_new, dec_new) : tuple of floats
			RA and Dec for `targname` after applying proper motion. In degrees.
			
			"""

	c = query_simbad_by_targname(targname)

	delta_t_yr = (mjd - 51544)/365. * u.yr
	ra_new = c.ra.deg*u.deg+(c.pm_ra_cosdec * delta_t_yr)
	dec_new = c.dec.deg*u.deg+(c.pm_dec * delta_t_yr)

	return(ra_new.value, dec_new.value)


def make_PAMcorr_image_UVIS(data, prihdr, scihdr, pamdir):
    """Creates the Pixel Area Map (PAM) image.
    Parameters:
        data : array
            Name of FITS file.
        pri : header
            Primary header of file for data.
        scihdr : header
            Header from science extension of data.
        pamdir : str
            Path to where pixel area maps for UVIS1 and/or UVIS2 are located.
    Returns:
        pamcorr_data : array
            PAM corrected data
    """

    x0 = int(np.abs(scihdr['LTV1']))
    y0 = int(np.abs(scihdr['LTV2']))
    x1 = int(x0 + scihdr['NAXIS1'])
    y1 = int(y0 + scihdr['NAXIS2'])

    if scihdr['CCDCHIP'] == 1:
        pam=fits.getdata(pamdir+'UVIS1wfc3_map.fits')
        pamcorr_data = data * pam[y0:y1,x0:x1]

    elif scihdr['CCDCHIP'] == 2:
        pam=fits.getdata(pamdir+'UVIS2wfc3_map.fits')
        pamcorr_data = data * pam[y0:y1,x0:x1]
    else:
        raise Exception('Chip case not handled.')

    return pamcorr_data

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
