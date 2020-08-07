from astroquery.mast import Observations
from astropy.table import Table, vstack
import glob
import os
import shutil

""" 
	Functions to search for and download WFC3 data from MAST using Astroquery. Intended to be 
	imported by other scripts. 

	Authors
	-------
	Clare Shanahan, Oct 2019
""" 


def query_by_propid_targ_filter(prop_ids, targnames='any', filters='any',
								file_type='any'):

	""" Astroquery query for data from `instrument` by target name, filter, 
	and proposal ID. Returns table of data products of file type(s) in 
	`file_type`.

	Parameters
	----------
	prop_ids : str or list of str
		Proposal ID(s)
	targnames : str list of str
		Exact target names that should be returned in query, others
		will be excluded. All spelling/name variations that might appear in MAST
		should be provided. If 'any', all available targets will be returned.
	filters: str of list of str.
		Filters that should be returned in query, other will be excluded. If
		'any', all available filters will be returned. 
	file_type : str or list or 
		File extention type(s) desired (i.e flt, flc, drz...), as a string for 
		a single type or a list for many. If 'any', all available file types 
		will be returned.

	Returns
	--------
	query_products : `astropy.table.Table`
		Table of products returned from query. 

	"""

	if type(prop_ids) != list:
		prop_ids = [prop_ids]
	if targnames == 'any':
		targnames = '*'
	if filters == 'any':
		filters = '*'

	query_products = Table()
	j = 0
	for i, prop_id in enumerate(prop_ids): #iterate to avoid timeout
		print('Querying for data from {}.'.format(prop_id))

		obsTable = Observations.query_criteria(obs_collection='HST', 
											   proposal_id=prop_id,
											   target_name=targnames,
											   filters=filters)
		if file_type == 'any':
			query_products = Observations.get_product_list(obsTable)
		else:
			if type(file_type) == str:
				file_type=file_type.upper()
			if type(file_type) == list:
				file_type=[x.upper() for x in file_type]
			query_products = Observations.get_product_list(obsTable)
			query_products = Observations.filter_products(query_products, 
							 productSubGroupDescription=file_type)

		if len(query_products) == 0:
			print('No records found in query.')
			j = 0
			if len(prop_ids) == 1:
				return
			continue

		if (i == 0) & (j == 0):
			query_products_total = query_products
		else:
			query_products_total = vstack([query_products, query_products_total])
		print('{} records found'.format(len(query_products)))
		j = 1

	return query_products_total


def query_by_data_id(dataset_ids, file_type):

	""" Astroquery query by file rootname(s) or ASN ID. Query will return 
		all records found for IDs in `dataset_ids` list (or string if 
		single ID) of type `file_type`. `file_type` can be set to a string 
		('FLT'), a list of strings (['FLT', DRZ']), or 'any'. 

		Parameters
		----------
		dataset_ids : str or list of str
			9-digit dataset id(s). Note, if these end in 'j' or 's', which
			occasionally happens due to failed observations, the id may not
			appear in the database. Change these files to end with 'q'
			instead
		file_type : str
			Fits file type (e.g 'flt', 'flc').

		Returns
		-------
		query_products : `astropy.table.Table`
			Results of query.
	"""

	obsTable = Observations.query_criteria(obstype='all', 
										   obs_collection='HST', 
										   obs_id=dataset_ids)
	if file_type == 'any':
		query_products = Observations.get_product_list(obsTable)
	else:
		if type(file_type) == str:
			file_type=file_type.upper()
		if type(file_type) == list:
			file_type=[x.upper() for x in file_type]
		query_products = Observations.get_product_list(obsTable)
		query_products = Observations.filter_products(query_products, 
						 productSubGroupDescription=file_type)

	# Initially all visit files were returned. Now select only specified IDs
	remove_rows = []
	for i, obs_id in enumerate(query_products['obs_id']):
		if obs_id not in dataset_ids:
			remove_rows.append(i)
	query_products.remove_rows(remove_rows)		

	return query_products

def download_products(query_products, output_dir=''):

	""" Downloads all products in `query_products` to `output_dir`. 

	Parameters
	----------
	query_products : `astropy.table.Table`
		Table of data products to download. 


	Notes
	-----
	Files are initially downladed temporary directory within `output_dir` 
	called 'temp', so if a subdirectory `temp` already exists within 
	`output_dir` an error is raised. ."""

	# format path for output dir
	if output_dir == '':
		output_dir = os.getcwd()
	output_dir = os.path.join(output_dir, '')

	# make temp dir in `output_dir`. error if it exists.
	assert(os.path.isdir(output_dir + 'temp') is False)
	os.makedirs(output_dir + 'temp')

	print('Downloading {} files.'.format(len(query_products)))
	Observations.download_products(query_products, 
								   download_dir=output_dir+'temp', 
								   mrp_only=False)

	# move files from mast download directories in `temp` to `output_dir`
	files = glob.glob(output_dir + 'temp/mastDownload/HST/*/*.fits')
	if len(files) > 0:
		print('Cleaning up temp directory.')
		for f in files:
			shutil.copy(f, output_dir+os.path.basename(f))

		# remove temp directory
		shutil.rmtree(output_dir + 'temp') # remove mast download dir 
