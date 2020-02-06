import numpy as np

def apply_gdr3mock_extinction(data, nside):
    from numpy.lib.recfunctions import rename_fields
    from healpy import ang2pix
    def distmod2dist(distmod):
        """distance modulus to distance in kpc"""
        return 10.**(distmod/5.-2.)
    distmod= np.linspace(4.,18.875,120)
    dists= np.hstack((0,distmod2dist(distmod)))

    # Load the right extinction map
    def get_hpx_number(nside):
        """
        returns hpx level for given nside
        """
        return(int(np.log2(nside)))

    hpx_level = get_hpx_number(nside)
    extmap = np.load("../input/ext_map/combined1719_hpxlevel_%d.npy" %(hpx_level)) 
    # prepending 0 extinction for the 0 distance bin
    zero = np.zeros(shape = len(extmap),dtype = extmap.dtype)
    extmap = np.concatenate((zero[:,None],extmap), axis = 1)
    
    nHpx = ang2pix(nside,data["glon"],data["glat"],nest = True, lonlat = True)
    extmap = extmap[nHpx]
    for i in range(len(extmap)):
        if i%1000000==0:
            print(i,len(extmap))
        data["exbv_schlegel"][i] = np.interp(data["rad"][i], dists, extmap[i])
    a0 = data['exbv_schlegel'] * 0.86 * 3.1
    a0[a0<0.] = 0.
    data['exbv_schlegel'] = a0
    data = rename_fields(data, {'exbv_schlegel':'a0'})	
    return(data)


def apply_bovy_extinction(data,hpx_res = 128):
    import mwdust
    from healpy import ang2pix
    from numpy.lib.recfunctions import rename_fields
    ## This updates the Galaxia extinction to the bovy2015 extinction.
    ## the resolution is given in the next line (healpix in factors of 2) 128 takes approximately 3hours 256 takes approximately 13hours
    ## The extinction is averaged over healpix pixel to make use of arrays, otherwise each star needs to be calculated individually.
    ## Still the distance information of the extinction map is used.
        
    healpix_resolution = hpx_res
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
            print(i, len(index_list))
        data['exbv_schlegel'][index_list[i]:index_list[i+1]] = dust(np.median(data['glon'][index_list[i]:index_list[i+1]]),np.median(data['glat'][index_list[i]:index_list[i+1]]),data['rad'][index_list[i]:index_list[i+1]])
    # 14% less than the SFD value and Rv of 3.1        
    a0 = data['exbv_schlegel'] * 0.86 * 3.1
    a0[a0<0.] = 0.
    data['exbv_schlegel'] = a0
    data = rename_fields(data, {'exbv_schlegel':'a0'})	
    return(data)

