import numpy as np

def apply_bovy_extinction(data):
	import mwdust
	from healpy import ang2pix	
	## This updates the Galaxia extinction to the bovy2015 extinction.
	## the resolution is given in the next line (healpix in factors of 2) 128 takes approximately 3hours 256 takes approximately 13hours
	## The extinction is averaged over healpix pixel to make use of arrays, otherwise each star needs to be calculated individually.
	## Still the distance information of the extinction map is used.
		
	healpix_resolution = 128
	print(len(data))	
	dust = mwdust.Combined15()
	phi = data['glon'] * np.pi/180.
	theta = (90.- data['glat']) * np.pi/180.
	number = ang2pix(healpix_resolution,theta,phi)
	index = np.arange(len(number))
	sorted_indices = np.argsort(number)
	number = number[sorted_indices]
	data = data[sorted_indices]
	index_list = [0]
	temp = min(number)
	for i in range(len(number)):
		if number[i] != temp:
			index_list.append(i)	
		temp = number[i]	
	index_list.append(len(number))
	for i in range(len(index_list)-1):
		if i%1000 == 0:		
			print i, len(index_list)
		data['exbv_schlegel'][index_list[i]:index_list[i+1]] = dust(np.median(data['glon'][index_list[i]:index_list[i+1]]),np.median(data['glat'][index_list[i]:index_list[i+1]]),data['rad'][index_list[i]:index_list[i+1]])

	print(np.median(data['exbv_schlegel']))
		
	return(data)



