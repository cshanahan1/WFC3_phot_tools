from astropy.io import fits
import glob
import os

from lacosmic_source import *
from lacosmic_best_params import best_params

def run_lacosmic(input_file, params='best', save_mask=True):

	""" Runs LACosmic algorithm on `input_file` to produce a CR rejected fits
		file. Optionally, saves the mask identifying CRs in image. CR corrected
		image and mask are output in same directory as `input_file`

		Parameters
		----------

		input_file : str
			Full path to input fits file.

		params : 'best' or list
			Parameters for LACosmic algorithm. If 'best', a set of predtermined 
			parameters for LACosmic will be used. Otherwise, a tuple can be 
			passed in with (gain, readnoise, sigclip, sigfrac, objlim, niter)
			in that order. 
		save_mask : bool
			If mask image with detected CRs should be saved. 

		Outputs
		-------
		input_file.clean.fits, and optionally input_file.mask.fits, in same
		directory as input.
	"""
	hdu = fits.open(f)
	if params == 'best':
		filt = hdu[0].header['filter']
		lacosmic_param_dict = best_params()
		sigclip, sigfrac, objlim, niter, sigclip_pf = lacosmic_param_dict[filt]

		post_flash = hdu[0].header['flshcorr']
		if post_flash == 'COMPLETE':
			sigclip = sigclip_pf
		params = lacosmic_best_params_dict[filt]
	
		gain, readnoise = 1.5, 3.0
	else:
		sigclip, sigfrac, objlim, niter = params

	print(input_file, sigclip)
	# # Now run LACosmic 

 #    c = cosmics.cosmicsimage(sci_array, gain=gain, readnoise=readnoise, 
 #    						 sigclip = sigclip, sigfrac = sigfrac, 
 #    						 objlim = objlim, verbose=False)
 #    c.run(maxiter = niter, verbose=False)
  


