import numpy as np



def convert(jj):
	assert len(jj)!= 0	
	entries = len(jj['rad'])
	dt = np.dtype([('rad',float), ('exbv_solar',float), ('teff',float), ('mag2',float), ('mag1',float), ('mag0',float), ('satid',int), ('vx',float), ('vy',float), ('vz',float), ('mtip',float), ('pz',float), ('px',float), ('py',float), ('feh',float), ('exbv_schlegel',float), ('lum',float), ('exbv_schlegel_inf',float), ('mact',float), ('glon',float), ('popid',int), ('glat',float), ('alpha',float), ('smass',float), ('ubv_r',float), ('ubv_u',float), ('partid',int), ('ubv_v',float), ('age',float), ('grav',float), ('ubv_b',float), ('ubv_i',float), ('ubv_h',float), ('ubv_k',float), ('ubv_j',float), ('fieldid',int), ('sdss_u',float), ('sdss_g',float), ('sdss_r',float), ('sdss_i',float), ('sdss_z',float), ('gaia_g',float), ('gaia_g_bp',float), ('gaia_g_rp',float)])

	x = np.empty((entries,), dtype = dt)
	for i,item in enumerate(jj.keys()):
		if not (item == 'log' or item == 'center'):
			x[item] = jj[item]
	return(x)