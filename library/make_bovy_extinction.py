import numpy as np
import mwdust
from healpy import ang2pix



def apply_bovy_extinction(data):
	## This makes the bovy2015 map extinction.
	## the resolution is given in the next line (healpix in factors of 2) 128 takes approximately 3hours 256 takes approximately 13hours
		
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

	print('median extinction: ',np.median(data['exbv_schlegel']))
		
	return(data)


