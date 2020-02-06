import numpy as np
from healpy import ang2pix
from numpy.lib.recfunctions import append_fields, stack_arrays, drop_fields

def get_nside(healpix_number = 1):
    """
    returns the nside for specific healpix level
    """
    return(np.power(2,healpix_number))

def gaia_hpx_factor(healpix_number = 1):
    """
    returns the number by which to divide the source_id in order to get a hpx number of a specific hpx level
    INPUT:
       healpix_number: the healpix level, ranging from 0 to 12, an integer
    OUTPUT:
       the gaia source id factor to get a specific hpx dicretization
    """
    return(np.power(2,35)*np.power(4,12-healpix_number))

def number_of_healpixels(healpix_number = 1):
    """
    returns the number of pixels for a specific level
    """
    return(np.power(4,healpix_number)*12)

def cumhpx_prepend(healpix_number = 1):
    """
    returns the number of cumulative healpixels below the given healpix level
    """
    return(int((number_of_healpixels(healpix_number)-12)/3))

def cumhpx(healpix_level, hpx):
    """
    returns the unique cumulative healpix number
    """
    assert (hpx < number_of_healpixels(healpix_level)), "healpixnumber %d is too great for healpix level %d" %(hpx,healpix_level)
    return(cumhpx_prepend(healpix_level)+hpx)

def sources_per_hpx(hpx, hpx_level):
    """
    this function returns the number of sources per healpix, hpx need to be ordered
    """
    sort = np.searchsorted(hpx,np.arange(number_of_healpixels(hpx_level)), side = "right")
    return(np.hstack((sort[0],np.diff(sort))))

def hpx_indices(hpx,hpx_level):
    """
    returns the index borders of the sorted healpixel numbers adding a 0 upfront, all hpx numbers need to be there and sorted by hpx
    """
    return(np.hstack((0,np.searchsorted(hpx,np.arange(number_of_healpixels(hpx_level)), side = "right"))))

def cumhpx_indices(cumhpx):
    """
    returns the index borders of the sorted cumhealpixel numbers adding a 0 upfront, array need to be sorted by cumhpx
    """
    return(np.hstack((0,np.searchsorted(cumhpx,np.unique(cumhpx), side = "right"))))

def assign_cumhpx(data,threshold = 3000, verbose = False):
    """
    this function adds cumhpx to a catalogue with ra and dec entries
    the threshold defines how much sources a hpx should contain at most. otherwise it goes a level deeper
    """
    data = append_fields(data, ["hpx","cumhpx","mask"], (np.zeros(len(data)),np.zeros(len(data)),np.zeros(len(data))),dtypes = int,usemask = False, asrecarray = True)
    hpxl_list = []
    # Loop through healpix levels in increasing order, stop when all healpixels contain less sources than the threshold
    for healpix_level in np.arange(12):
        if (len(data)==0):
            break
        if verbose:
            print(healpix_level)
        # assign healpix of that oder
        nside = get_nside(healpix_level)
        data.hpx = ang2pix(nside,data['ra'],data['dec'],nest = True, lonlat = True)
        data = np.sort(data, kind = 'mergesort', order = 'hpx')
        src_per_hpx = sources_per_hpx(data.hpx, healpix_level)
        src_ind = hpx_indices(data.hpx,healpix_level)
        # if Nsources is below threshold assign the cumhpx of that healpix and delete these data from the next iteration
        for i,item in enumerate(src_per_hpx):
            if (item < threshold):
                if verbose:
                    print(healpix_level,i,item,src_ind[i],src_ind[i+1])
                chpx = cumhpx(healpix_level=healpix_level,hpx=i)
                data.cumhpx[src_ind[i]:src_ind[i+1]] = chpx
                data.mask[src_ind[i]:src_ind[i+1]] = 1
        cut = (data.mask == 1)
        hpxl_list.append(data[cut])
        data = data[~cut]
    data = stack_arrays(hpxl_list,usemask=False, asrecarray=True,autoconvert=True)
    assert (np.unique(data.mask)[0]==1),"not all sources have been assigned"
    # Cleaning of useless fields and assigning good hpx
    data = drop_fields(data,"mask", usemask=False, asrecarray=True)
    max_level = highest_hpx_level_from_largest_cumhpx(max(data.cumhpx))
    nside = get_nside(max_level)
    data.hpx = ang2pix(nside,data['ra'],data['dec'],nest = True, lonlat = True)
    data = np.sort(data, kind = 'mergesort', order = 'hpx')                                                  
    return(data)

def highest_hpx_level_from_largest_cumhpx(max_cumhpx):
    """
    given the largest cumhpx value it returns the needed hpx level in order to map the grid
    """
    for i in range(12):
        if (max_cumhpx < cumhpx_prepend(i)):
            break
    return(i-1)

def hpx2cumhpx(cat):
    """
    pointing from hpx to cumhpx, this is how the cumhpx structure of one catalogue can be transferred to another, array needs to be sorted by hpx
    """
    # Slicing by hpx level which should be as fine as the finest cumhpx level
    max_hpx_level = highest_hpx_level_from_largest_cumhpx(max(cat.cumhpx))
    scr_idx = hpx_indices(cat.hpx,max_hpx_level)
    result = []
    for i in range(len(scr_idx)-1):
        temp = cat[scr_idx[i]:scr_idx[i+1]]
        # taking the unique cumhpx per hpx pixel
        cumhpx_pointer = np.unique(temp.cumhpx)
        if (len(cumhpx_pointer)==0):
            print('hpx %d is empty' %(i))
            result.append(cumhpx(max_hpx_level, i))
            continue
        assert(len(cumhpx_pointer)==1),"more cumhpx for single hpx, cumhpxpointer_len = %d" %(len(cumhpx_pointer))
        result.append(cumhpx_pointer[0])
    result = np.hstack(result)
    return(result)

def highest_hpx_level_from_largest_hpx(max_hpx):
    """
    given the largest hpx value it returns the corresponding hpx level
    """
    for i in range(12):
        if (max_hpx < number_of_healpixels(i)):
            break
    return(i)

def transfer_cumhpx_structure(data_with_cumhpx,new_cat):
    """
    This function gives a new catalogue the same cumhpx structure as an existing one
    """
    # Extract the structure from old cat
    cumhpx_array = hpx2cumhpx(data_with_cumhpx) ## corresponding cumhpx
    max_hpx_level = highest_hpx_level_from_largest_cumhpx(max(data_with_cumhpx.cumhpx))
    hpx_array = np.arange(number_of_healpixels(max_hpx_level)) ## hpx from 0 to hpx_max
    healpix_level = highest_hpx_level_from_largest_hpx(max(data_with_cumhpx.hpx))
    assert(max_hpx_level==healpix_level), "healpixlevel from cumhpx and hpx are not the same"    
    # Prepare the new catalogue
    new_cat = append_fields(new_cat, ["hpx","cumhpx"], (np.zeros(len(new_cat)),np.zeros(len(new_cat))),dtypes = int,usemask = False, asrecarray = True)
    nside = get_nside(healpix_level)
    new_cat.hpx = ang2pix(nside,new_cat['ra'],new_cat['dec'],nest = True, lonlat = True)
    new_cat = np.sort(new_cat, kind = 'mergesort', order = 'hpx')
    
    # use the cumhpx with hpx_array identification to assign the same cumhpx structure to the second catalogue
    src_ind = hpx_indices(new_cat.hpx,healpix_level)
    for i in range(len(src_ind)-1):
        new_cat.cumhpx[src_ind[i]:src_ind[i+1]] = cumhpx_array[i]
    return(new_cat)

def likelihood(ndata,nmodel):
    """
    equation 1 of https://arxiv.org/abs/1904.04350
    !what happens for nmodel or dmodel == 0?
    !same proportion deviations in higher density bins are valued stronger into the likelihood
    """
    if (ndata!=0 and nmodel !=0):
        return(-1* (nmodel - ndata + ndata * np.log(np.divide(ndata,nmodel))))
    elif (ndata==0 and nmodel !=0):
        return(-1* (nmodel))
    elif (ndata == 0 and nmodel == 0):
        return(0.)
    elif (ndata !=0 and nmodel == 0):
        return(-1*(ndata))

def hess_likelihood(cumhpx_cell_data,cumhpx_cell_model, bins = 10, avoid_bp = True):
    """
    calculates the hess diagram likelihood for a cumhpx cell, g bp rp needed. 
    binlimits according to data ranges, might miss model data! bin resolution affects likelihood significantly
    """
    # same edges for both data sets
    if avoid_bp:
        nd,xedges,yedges = np.histogram2d(cumhpx_cell_data.phot_g_mean_mag - cumhpx_cell_data.phot_rp_mean_mag,cumhpx_cell_data.phot_g_mean_mag, bins = bins)
        nm,_,_ = np.histogram2d(cumhpx_cell_model.phot_g_mean_mag - cumhpx_cell_model.phot_rp_mean_mag,cumhpx_cell_model.phot_g_mean_mag, bins = [xedges,yedges])
    else:
        nd,xedges,yedges = np.histogram2d(cumhpx_cell_data.phot_bp_mean_mag - cumhpx_cell_data.phot_rp_mean_mag,cumhpx_cell_data.phot_g_mean_mag, bins = bins)
        nm,_,_ = np.histogram2d(cumhpx_cell_model.phot_bp_mean_mag - cumhpx_cell_model.phot_rp_mean_mag,cumhpx_cell_model.phot_g_mean_mag, bins = [xedges,yedges])
    nd = np.hstack(nd)
    nm = np.hstack(nm)
    assert(len(nm)==len(nd)),"not the same binning"
    # likelihood per cmd bin, might need matrix form in the future
    ln_list = []
    for i in range(len(nd)):
        ln_list.append(likelihood(nd[i],nm[i]))
    result = sum(ln_list)
    diff = len(cumhpx_cell_model) - len(cumhpx_cell_data)
    return(result, diff)

def hess_likelihood_per_cumhpx(data,model,bins = 10):
    """
    splits the catalogue into cumhpx cells and calculates hess likelihood per cell. 
    returns array of likelihoods corresponding to cumhpx ids
    Beware usually in this package catalogs are sorted by hpx, here it is cumhpx!
    """
    # cumhpx structure in indices
    data = np.sort(data, kind = 'mergesort', order = 'cumhpx')
    model = np.sort(model, kind = 'mergesort', order = 'cumhpx')
    data_cumhpx = np.unique(data.cumhpx)
    model_cumhpx = np.unique(model.cumhpx)
    # This throws out all cumhpx in the model that are not filled in the data
    if not (np.all(data_cumhpx==model_cumhpx)):
        throw_out = []
        for item in model_cumhpx:
            if item not in data_cumhpx:
                print('cumhpx %d is deleted in the model because of no representation in the data, (might affect likelihood)' %(item))
                cut = (model.cumhpx==item)
                model = model[~cut]
    data_idx = cumhpx_indices(data.cumhpx)
    model_idx = cumhpx_indices(model.cumhpx)
    # In order to fix this behaviour we will need to see that data and model idx are pointing to the same cumhpx
    assert(len(data_idx)==len(model_idx)),"data and model have different cumhpx structure"
    overall_like = []
    overall_diff = []
    # cut data and model separately in cumhpx bins and get the hess likelihood per bin
    for i in range(len(data_idx)-1):
        temp_data = data[data_idx[i]:data_idx[i+1]]
        temp_model = model[model_idx[i]:model_idx[i+1]]
        # Also we need to think about the likelihood if data or model are empty
        hl, diff = hess_likelihood(temp_data,temp_model,bins = bins)
        overall_like.append(hl)
        overall_diff.append(diff)
    return(np.hstack(overall_like), np.hstack(overall_diff))

def plot_mollweide_linear(data,label = 'lnlike', savefile=False, name = 'likelihood.png', title = 'likelihood'):
    import healpy as hp
    import matplotlib.pylab as plt
    
    vrange = np.nanmax(np.abs(data))
    map_mollweide = data
    hp.mollview(map_mollweide, cbar = True, nest = True, coord= "CG", unit = label,notext =True, min = -vrange, max = vrange, cmap=plt.get_cmap("bwr"))
    plt.title(title)
    if savefile:
        plt.savefig(name)
        plt.close()
    else:
        plt.show()

def plot_mollweide_log(data, label = 'lnlike', savefile=False, name = 'likelihood.png', title = 'likelihood'):
    import healpy as hp
    from matplotlib.colors import LogNorm
    import matplotlib.pylab as plt
    norm = LogNorm()
    data[data==0]=np.nan
    map_mollweide = data        
    total = np.nansum(map_mollweide)
    hp.mollview(map_mollweide, cbar = True, min=None, max=None, nest = True,norm = norm, coord= "CG", unit = label,notext =True)
    plt.title(title)
    if savefile:
        plt.savefig(name)
        plt.close()
    else:
        plt.show()
    
def assess_likelihood(data,model,threshold=10000,cmd_bins=4):
    """
    This function returns the loglikelihood of one catalogue compared to the other based on the CMD likelihood per cumhpx
    INPUT
       data = one catalogue needs to contain ra, dec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag
       model = another catalogue with the same entries
       threshold = the highest number of stars which shpould occur per cumhpx
       cmd_bins = number of bins in the CMD, the same along both axis
    OUTPUT
       cumhpx_like = the log likelihood per cumhpx
    To get the corresponding cumhpx ID and starcount per cumhpx one can do:
    cumhealpix, ct = np.unique(data.cumhpx, return_counts=True)
    """
    ## Assign a cumhpx structure
    if 'cumhpx' in data.dtype.names:
        print('data already has cumhpx structure')
    else:
        data = assign_cumhpx(data, threshold)
    ## Transfer cumhpx structure
    if 'cumhpx' in model.dtype.names:
        print('model already has cumhpx structure, might be incompatible with data')
    else:
        model = transfer_cumhpx_structure(data,model)
    cumhpx_like, cumhpx_diff = hess_likelihood_per_cumhpx(data,model,bins = cmd_bins)
    return(data,cumhpx_like, cumhpx_diff)

def get_hpx_structure_for_plotting(cumhpx_like, data):
    '''
    for plotting purposes the cumhpx structure needs to be written back into hpx compatible structure.
    This is done here for the likelihood and the starcounts per cumhpx which is not equal area any more.
    INPUT
       cumhpx_like = the likelihood per cumhpx, i.e. output from assess_likelihood
       data = a catalogue that contains hpx
    OUTPUT
       hpx_like = likelihood per hpx
       hpx_dens = starcount per hpx (not exactly true because the starcount is per cumhpx but for visualisation the number is repeated)
    '''
    data = np.sort(data, kind = 'mergesort', order = 'hpx')
    # pointing of cumhpx to stellar density and to loglike
    cumhealpix, ct = np.unique(data.cumhpx, return_counts=True)
    # pointing of cumhpx to hpx 
    cum2hpx = hpx2cumhpx(data)
    healpix_level = highest_hpx_level_from_largest_hpx(max(data.hpx))
    hpx = np.arange(number_of_healpixels(healpix_level)) ## hpx from 0 to hpx_max
    hpx_like = np.ones_like(hpx)*np.nan
    hpx_dens = np.ones_like(hpx)*np.nan
    # plot stellar density and loglike in mollweide projection
    for i,item in enumerate(cumhealpix):
        cut = (cum2hpx == item)
        hpx_like[cut] = cumhpx_like[i]
        hpx_dens[cut] = ct[i]
    return(hpx_like,hpx_dens)    

def likelihood_assessment(x,fSample, nside = 128,threshold = 10000, cmd_bins = 5, comp_file = "../output/GDR2_207/GDR2_207_cleaned_0.0025sampling.fits"):
    """
    This function returns the likelihood compared to the GDR2 catalogue
    INPUT
       x = catalogue like Galaxia
       fSample = the sampling of that catalogue. Up to 0.0025 is represented by the GDR2data
    OUTPUT
       cumhpx_like_new = the likelihood in cumhpx structure (assess likelihood from the sum of this)
       hpx_like = likelihood per hpx for plottingpurposes (simply summing would give a wrong number)
       hpx_dens = the number of stars per hpx for plotting purposes (simply summing would give a wrong number)
    """
    from astropy.io import fits
    from cumhpx import assess_likelihood, get_hpx_structure_for_plotting, delete_sources_below_maglim, delete_magellanic_clouds
    x = x.view(np.recarray)
    # get GDR2data and apply downsampling and magcut per hpx
    data = fits.getdata(comp_file)
    downsampling = fSample / 0.0025
    if downsampling > 1:
        print('could not assess likelihood, because GDR2data only has sampling of 0.0025')
    choice = np.random.choice(np.arange(len(data)), size = int(len(data)*downsampling), replace = False)
    data = data[choice]
    data = delete_sources_below_maglim(data,nside,'g',True)
    # taking the maglims of rp does throw away many sources
    data = delete_sources_below_maglim(data,nside,'rp',True)
    # throw away magellanic clouds
    data = delete_magellanic_clouds(data)
    # assess likelihood and prepare arrays for plotting
    data, cumhpx_like_new, cumhpx_diff = assess_likelihood(data,x, threshold = threshold, cmd_bins = cmd_bins)
    print('likelihood: ', sum(cumhpx_like_new))
    print('total diff: ', sum(np.abs(cumhpx_diff)))
    print('too many: %d, too little: %d' %(sum(cumhpx_diff[cumhpx_diff>0]) ,sum(cumhpx_diff[cumhpx_diff<0]) ))
    hpx_diff, hpx_dens = get_hpx_structure_for_plotting(cumhpx_diff,data)
    hpx_like, hpx_dens = get_hpx_structure_for_plotting(cumhpx_like_new,data)
    
    return(cumhpx_like_new,hpx_like, hpx_dens, len(data),cumhpx_diff, hpx_diff)

def delete_sources_below_maglim(data,nside,band = 'g', verbose = False):
    """
    This function deletes all sources that are below a specific magnitude limit
    INPUT
       data = a catalog with ra and dec and phot_g_mean_mag
       nside = nside level of the maglim file
       verbose = print number of stars lost
    OUTPUT
       data = only the stars within the maglim
    """
    import healpy as hp
    maglim = np.load("../input/gdr2_maglim/maglim_%s_in_nside%d.npy" %(band,nside))
    hpx = hp.ang2pix(nside,data['ra'],data['dec'],nest = True, lonlat = True)
    maglim_star = maglim[hpx]
    visible = (data['phot_%s_mean_mag' %(band)] <= maglim_star) 
    lendata = len(data)
    data = data[visible]
    if verbose:
        print("nstar before cleaning: %d" %(lendata))
        print("nstar after %s-band cleaning: %d" %(band,len(data)))
    return(data)

def do_likelihood(outputDir, outputFile = "GDR2mock_20.7Gmag",fSample = 0.001, comp_file = "../output/GDR2_207/GDR2_207_cleaned_0.0025sampling.fits", savefile = True):
    """
    applies the likelihood including maglim
    OUTPUT
       starcounts = the number of stars, after applying observational cuts
       likelihood = the likelihood from hess diagrams
       diffplus = how many stars model overpredicts over the sky
       diffminus = how many stars model underpredicts over the sky
    """
    from astropy.io import fits
    from cumhpx import plot_mollweide_log, delete_sources_below_maglim, likelihood_assessment, delete_magellanic_clouds

    x = fits.getdata(outputDir + "/" + outputFile + '.fits')
    nside_maglim = 128
    print('GDR2mock')
    x = delete_sources_below_maglim(x,nside_maglim,'g',True)
    x = delete_sources_below_maglim(x,nside_maglim,'rp',True)
    x = delete_magellanic_clouds(x)
    print('GDR2')
    cumhpx_like,hpx_like,hpx_dens,lengdr2data,cumhpx_diff, hpx_diff = likelihood_assessment(x,fSample,comp_file=comp_file)
    print(len(x),lengdr2data)
    # add the numbers into the plot.
    plot_mollweide_linear(hpx_diff,'diff (m-d)',savefile=savefile, name = outputDir + "/" + outputFile + "_starcount_diff.png", title = "total deviation %d" %(sum(np.abs(cumhpx_diff))))
    plot_mollweide_log(-1*hpx_like,'-lnlike %.0f' %(sum(cumhpx_like)),savefile=savefile, name = outputDir + "/" + outputFile + "_likelihood.png", title = "gdr2: %d, mock: %d" %(lengdr2data,len(x)))
    starcounts = len(x)
    likelihood = sum(cumhpx_like)
    diffplus = sum(cumhpx_diff[cumhpx_diff>0])
    diffminus = sum(-cumhpx_diff[cumhpx_diff<0])
    return(starcounts,likelihood,diffplus,diffminus)

def delete_magellanic_clouds(data, verbose = True):
    """
    this function deletes the magellanic clouds from the sky and returns the cleaned catalogue
    """
    smc_l = 302.8
    smc_b = -44.3
    smc_radius = 7
    lmc_l = 280.5
    lmc_b = -32.9
    lmc_radius = 12
    if verbose:
        print(len(data))
    cut = (np.sqrt((smc_l-data['l'])**2+(smc_b-data['b'])**2)<smc_radius)
    data = data[~cut]
    if verbose:
        print("after smc cut: %d" %(len(data)))
    cut = (np.sqrt((lmc_l-data['l'])**2+(lmc_b-data['b'])**2)<lmc_radius)
    data = data[~cut]
    if verbose:
        print("after lmc cut: %d" %(len(data)))
    return(data)

def local_normalisation(catName, verbose = False):
    '''
    calculates the local densities from the 100pc sample mock catalogue and returns the numbers
    INPUT
       catName = name of the catalogue
    OUTPUT
       stellardiskmass = the local mass density of thick and thin disk stars in Msun/pc^3*10^-3
       wddiskmass = the local mass density of thick and thin disk white dwarfs in Msun/pc^3*10^-3
       totalcount = total starcount within 100pc, simply all sources in the mock, no cuts applied
    '''
    from astropy.io import fits
    if verbose:
        print('## Local normalisation, compare to Czekaj+14 table 7 ##')
    outputDir = '../output/%s_100pc' %(catName)
    #outputDir = '../output/mockcat_new_SFR_new_IMF_100pc_more_thick'
    x = fits.getdata(outputDir + '/' + 'GDR2mock_20.7Gmag.fits')
    totalcount = len(x)
    smaller_sphere = 50
    sphere = (4/3)*np.pi*smaller_sphere**3
    local_sum = 0
    for popid in np.arange(9):
        selection = (x.popid==popid) & (x.logg < 7) & (np.divide(1,x.parallax)<smaller_sphere/1000)
        if (popid!=8):
            if verbose:
                print('popid=',popid,' %.1f Msun/pc^3*10^-3' %(sum(x.smass[selection])/sphere*1000))
            local_sum+=sum(x.smass[selection])/sphere*1000
        else:
            if verbose:
                print('popid=',popid,' %.1f Msun/pc^3*10^-6' %(sum(x.smass[selection])/sphere*1e6))
    if verbose:
        print("## total sum of stars: %.1f" %(local_sum))
    stellardiskmass = local_sum
    selection = (x.popid<=6) & (x.logg > 7)& (np.divide(1,x.parallax)<smaller_sphere/1000)
    wddiskmass = 0
    if verbose:
        print('thindisc WD',' %.1f Msun/pc^3*10^-3 %.1f mact' %(sum(x.smass[selection])/sphere*1000,sum(x.mact[selection])/sphere*1000))
    wddiskmass+=sum(x.mact[selection])/sphere*1000
    selection = (x.popid==7) & (x.logg > 7)& (np.divide(1,x.parallax)<smaller_sphere/1000)
    if verbose:
        print('thick disc WD',' %.1f Msun/pc^3*10^-3 %.1f mact' %(sum(x.smass[selection])/sphere*1000,sum(x.mact[selection])/sphere*1000))
        print('####################################')
    wddiskmass+=sum(x.mact[selection])/sphere*1000
    return(stellardiskmass,wddiskmass,totalcount)

    
