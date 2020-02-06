from sklearn.ensemble import ExtraTreesRegressor
import pickle
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from numpy.lib.recfunctions import append_fields

# define a function that returns the ext in a band from model, av, and linear or quadratic fitting,
def av2ext(band_name,log_teff,log_grav,meh_ini,log_lum,av, av_axis = np.array([0,1,2,3,5,10,20])):
    """
    get extinctions for a specific band for an array of Av,and stellar parameters
    INPUT
       band_name = specific band (this will load a pre-trained extinction model)
       log_teff = teff
       log_grav = logg
       meh_ini = metallicity
       log_lum = luminosity
           the above stellar parameters have been used in the trainig
       av = av values of the stars
       av_axis = model grid of Av values (on these values the extinction model has been trained)
    OUTPUT
       ext_array are the extinctions in the band for which the model was trained
    """
    from scipy.interpolate import interp1d
    # Load the extinction model
    filename = "../input/ext_models/ext_model_%s_band" %(band_name)
    ext_model = pickle.load(open(filename, 'rb'))
    
    # prepare the features
    X = np.vstack((log_teff,log_grav,meh_ini,log_lum)).T
    # predict the extinction grid
    prediction = ext_model.predict(X)
    #interpolate cubic to a finer grid
    x_axis = np.linspace(0,20,41)
    f = interp1d(av_axis,prediction, kind = 'cubic', copy = False, bounds_error=False, fill_value='extrapolate', assume_sorted = True)
    result = f(x_axis)
    # prepare ext array
    ext_array = np.zeros_like(av)
    assert(len(ext_array)==len(result)),"Problem with dimensionality of ext_model input cat and av array"
    # for each star linearly interpolate the av value from the finer grid
    for i,item in enumerate(result):
        ext_array[i] = np.interp(av[i],x_axis,result[i], left = 0.)
    # if av > 20 then extrapolate linearly
    cut = (av>20)
    ext_array[cut] = ext_array[cut] * (av[cut]/20.)
    return(ext_array)

def apply_extinction_curves(data,indexing,ext_curve,band):
    """
    this function converts A0 to extinction in specific bands.
    INPUT
       data = Galaxia catalog data
       indexing = a mapping from the parsec_index of each source in the catalogue onto the isochrone/extinction grid
       ext_curve = an array of the form [parsec_index,photometric_band,extinction_grid] that has been prepared from isochrones beforehand
       band = name of the photometric band
    OUTPUT
       ext_array = the absorption in the specified band according to the a0 provided in data
    """
    from scipy.interpolate import interp1d    
    # These bands correpond to the ones trained on and to the bands in the galaxia catalog
    bands = ['g','bpft','bpbr','rp','rvs']
    for i, item in enumerate(bands):
        if item == band:
            break
    # ext_curve is projected onto that band only
    ext_curve = ext_curve[:,i,:]
    avs = [0,1,2,3,5,10,20]
    # prepend zeros and bring onto the data shape
    ext_curve = np.concatenate((np.zeros(shape=(len(ext_curve[indexing])),dtype = np.float32)[:,None],ext_curve[indexing][:,:]), axis=1)
    # interpolate the ext_curve cubically onto a finer grid
    x_axis = np.linspace(0,20,41)
    f = interp1d(avs,ext_curve, kind = 'cubic', copy = False, bounds_error=False, fill_value='extrapolate', assume_sorted = True)
    result = f(x_axis)
    # for each star linearly interpolate the av value from the finer grid
    ext_array = np.zeros(shape = (len(data)),dtype = np.float32)
    for i,item in enumerate(result):
        ext_array[i] = np.interp(data[i]['a0'],x_axis,result[i], left = 0.)
    # if av > 20 then extrapolate linearly
    cut = (data['a0']>20)
    ext_array[cut] = ext_array[cut] * (data['a0'][cut]/20.)
    return(ext_array)


def parsec_index_mapper(data):
    """
    creates an index array that can be applied on the isochrones or extinction_curve grid
    """
    iso = np.load("../input/isochrones/parsec.npy")
    return()

def train_ext(feature,labels,band_name):
    """
    Training an extinction model from stellar parameters for a specific band
    """
    # creating a subsample on which to train
    test = np.random.choice(np.arange(len(feature)),size = int(len(feature)/100), replace = False)
    training_feature = feature[test]
    training_label = g_ext[test]
    # the rest is stored for validation
    test_feature = np.delete(feature, test)
    test_label = np.delete(g_ext,test, axis = 0)
    # prepare input for regression
    X = np.vstack((training_feature["log_teff"],training_feature["log_grav"],training_feature["meh_ini"],training_feature["log_lum"])).T
    y = training_label#[:,-1]
    X_test = np.vstack((test_feature["log_teff"],test_feature["log_grav"],test_feature["meh_ini"],test_feature["log_lum"])).T
    y_test = test_label#[:,-1]
    # train the model
    model = ExtraTreesRegressor(verbose = 1, n_jobs=-1, max_depth=20, criterion = 'mse',n_estimators=10, max_features = 4)
    model.fit(X,y)
    # validate
    y_pred = model.predict(X_test)
    print("90% of the validation set from Av = 1 has less than: ",np.percentile(np.abs(y_pred[:,1] - y_test[:,1]),90), "error in Av")
    print("model feature importance in the order: teff, logg, metallicity, lum: ",model.feature_importances_)
    # Save and reload model
    filename = "ext_model_" + band_name
    pickle.dump(model, open(filename, 'wb'))



def add_parsec_index(data,iso, verbose = False):
    """
    This function assigns an index to each star in Galaxia. Beware that not all indexes will be represented by the parsec isochrone table
    """
    from numpy.lib.recfunctions import append_fields
    def return_index_feh(feh):
        dfeh = 0.05
        offset = 1.5
        return(int((feh+offset)/dfeh))
    min_value = min(iso['meh_ini'])
    max_value = max(iso['meh_ini'])
    min_index_feh = return_index_feh(min_value)
    max_index_feh = return_index_feh(max_value)
    stretch = max_index_feh - min_index_feh
    if verbose:
        print('teff values in parsec')
        print('teff min max value: ', min_value, max_value )
        print('teff min max index: ', min_index_feh, max_index_feh)
        print('dimension cut into %d pieces' %(stretch))

    def return_index_teff(teff):
        dteff = 0.02
        return(int(teff/dteff))
    min_value = min(iso['log_teff'])
    max_value = max(iso['log_teff'])
    min_index_teff = return_index_teff(min_value)
    max_index_teff = return_index_teff(max_value)
    stretch = max_index_teff - min_index_teff
    if verbose:
        print('teff values in parsec')
        print('teff min max value: ', min_value, max_value )
        print('teff min max index: ', min_index_teff, max_index_teff)
        print('dimension cut into %d pieces' %(stretch))
    def return_index_lum(lum):
        dlum = 0.05
        offset = 5.
        return(int((lum+offset)/dlum))
    min_value = min(iso['log_lum'])
    max_value = max(iso['log_lum'])
    min_index_lum = return_index_lum(min_value)
    max_index_lum = return_index_lum(max_value)
    stretch = max_index_lum - min_index_lum
    if verbose:
        print('lum values in parsec')
        print('lum min max value: ', min_value, max_value )
        print('lum min max index: ', min_index_lum, max_index_lum )
        print('dimension cut into %d pieces' %(stretch))    
    
    def return_index(feh,teff,lum):
        """
        feh given in dex
        teff given in log teff
        lum given in log lum
        """
        index_feh = return_index_feh(feh)
        if index_feh < min_index_feh:
            index_feh = min_index_feh
        elif index_feh > max_index_feh:
            index_feh = max_index_feh

        index_teff = return_index_teff(teff)
        if index_teff < min_index_teff:
            index_teff = min_index_teff
        elif index_teff > max_index_teff:
            index_teff = max_index_teff
        index_teff *= 1000

        index_lum = return_index_lum(lum)
        if index_lum > max_index_lum:
            index_lum = max_index_lum
        elif index_lum < min_index_lum:
            index_lum = min_index_lum
        index_lum *= 1000 * 1000
        assert(index_feh >= 0)
        assert(index_teff >= 0)
        assert(index_lum >= 0)
        return (index_feh + index_teff + index_lum)
    
    indizes_g = np.zeros(len(data),dtype = np.int32)
    for i in range(len(data)):
        indizes_g[i] = return_index(data['feh'][i],data['teff'][i],data['lum'][i])
    data = append_fields(data,"parsec_index",indizes_g,usemask=False)
    return(data)

def map_galaxia_indexes_onto_parsec(data, iso, verbose = False):
    #check problems with Galaxia index cells to parsec index cells
    problems = 0
    total_fail=0
    distances = []
    indizes = iso['parsec_index']
    unique = np.unique(indizes)
    indizes_g = data['parsec_index']
    for i,item in enumerate(np.unique(indizes_g)):
        if i % 1000 == 0:
            if verbose:
                print(i, len(np.unique(indizes_g)))
        try:
            assert(item in unique)
        except:
            problems += 1
            temp = data[np.where(indizes_g == item)]
            #We are staying with the same metallicity
            temp_cut = np.where(indizes%1000==item%1000)
            temp_x = iso[temp_cut]
            temp_indizes = indizes[temp_cut]
            total_fail += len(temp)
            distance = (np.median(temp['teff'])-temp_x['log_teff'])**2 + (np.median(temp['lum']) - temp_x['log_lum'])**2
            min_distance = min(distance)
            if verbose:
                print(i,np.sqrt(min_distance))
            distances.append(np.sqrt(min_distance))
            indizes_g[np.where(indizes_g == item)] = temp_indizes[np.where(distance == min_distance)][0]  
    if verbose:
        print("Problems with %d cells in galaxia with a total of %d stars" %(problems,total_fail))
        print('distances of shifted cells in loglum logteff space, median and 95 percentile, and 95 and max')
        print(np.median(distances), np.percentile(distances,95), np.percentile(distances,99), np.percentile(distances,100))
    for i,item in enumerate(np.unique(indizes_g)):
        assert(item in unique)
    data['parsec_index'] = indizes_g
    return(data)

def plot_sky_map(x,nside,fSample,appMagLimits1,outputDir,outputFile):
    """
    Plots a skymap for this catalogue
    """
    import healpy as hp
    import matplotlib.cm as cm
    from matplotlib.colors import LogNorm
    import matplotlib
    import matplotlib.pylab as plt
    from defaults import getLogTickMarks
    NSIDE = nside
    oversampling = 1./fSample
    total = len(x)
    x['l'] = x['l'] * (np.pi/180.)
    x['b'] = (90.-x['b']) * (np.pi/180.)
    print('total number of stars = %.d' %(total * oversampling))
    count = hp.ang2pix(NSIDE,x['b'],x['l'])
    sqdegs = 41253
    pixels = NSIDE * NSIDE * 12
    pixel_per_sqdeg = pixels / float(sqdegs)
    min_density = oversampling * pixel_per_sqdeg 
    m = np.arange(hp.nside2npix(NSIDE))
    density = np.zeros(hp.nside2npix(NSIDE))
    for item in count:
        density[item] += oversampling * pixel_per_sqdeg
    print(np.min(density))
    print(sum(density)/pixel_per_sqdeg)
    print(sum(density))
    cmap = cm.jet
    cmap.set_under(cmap(0.0))
    cmap.set_over(cmap(1.0))
    norm=LogNorm()
    minVal = np.nanmin(density[density>0])
    maxVal = np.nanmax(density[density<+np.inf])
    density[density<minVal] = minVal
    x['l'] = x['l'] * (180./np.pi)
    cbLabel=r'$n_{\rm stars}$ [sq.deg$^{-1}$]'
    hp.mollview(density, unit=cbLabel, min=minVal, max=maxVal, nest=False, title='', norm=norm, cmap=cmap, cbar=None)
    fig = plt.gcf()
    ax = plt.gca()
    pos1 = ax.get_position() # get the original position 
    pos2 = [pos1.x0, pos1.y0 + 0.06,  pos1.width, pos1.height]
    ax.set_position(pos2) # set a new position
    im = ax.get_images()[0]
    cbAx = fig.add_axes([0.1, 0.12, 0.8, 0.03])
    cb = plt.colorbar(im, cax=cbAx, orientation='horizontal', )
    cb.ax.minorticks_on()
    tickMarks = getLogTickMarks(minVal, maxVal)
    minorticks = im.norm(tickMarks)
    cb.ax.xaxis.set_ticks(minorticks, minor=True)
    cb.solids.set_edgecolor("face")
    cb.set_label(cbLabel)
    plt.title(r"%.1fm stars with m$_\mathrm{G}$<%fmag" %(total*oversampling/1e6,appMagLimits1), fontsize = 12)
    plt.savefig(outputDir + "/" + outputFile + ".png")
    plt.clf()
    plt.close()




def adding_velocities(x):
	from astropy.coordinates import ICRS, CartesianRepresentation, CartesianDifferential, Galactic, Galactocentric
	
	#old routine 
	#icrs = ICRS(x= (x['px'])*u.kpc, y=x['py']*u.kpc, z=x['pz']*u.kpc, v_x=x['vx']*u.km/u.s, v_y=x['vy']*u.km/u.s, v_z=x['vz']*u.km/u.s, representation=CartesianRepresentation, differential_cls=CartesianDifferential)
	#galactic = icrs.transform_to(Galactic)
	#icrs = galactic.transform_to(ICRS)

	galactocentric = Galactocentric(x = (x['px']-8.0)*u.kpc,
                          y = x['py']*u.kpc,
                          z = (x['pz']+0.015)*u.kpc, 
                          v_x=(x['vx']+11.1)*u.km/u.s,
                          v_y=(x['vy']+239.08)*u.km/u.s, 
                          v_z=(x['vz']+7.25)*u.km/u.s,
                          galcen_distance = 8.0*u.kpc,
                          z_sun = 15*u.pc,
                          galcen_v_sun=CartesianDifferential(d_x=11.1*u.km/u.s, d_y=239.08*u.km/u.s, d_z=7.25*u.km/u.s))#, unit=None,)
	icrs = galactocentric.transform_to(ICRS)

	radial_velocity = icrs.radial_velocity.value # km/s
	pm_ra_cosdec = icrs.pm_ra_cosdec.value # mas/yr
	pm_dec = icrs.pm_dec.value # mas/yr
	#total_pm = np.sqrt(pm_ra_cosdec * pm_ra_cosdec + pm_dec * pm_dec) # mas/yr

	pc_to_km = 3.086e13
	ARCSECONDS2RADIANS = 4.84813681109536e-6
	year_to_second = 3.154e+7
	conversion = ARCSECONDS2RADIANS * (1./year_to_second) * pc_to_km
	#transverse_velocity = total_pm * x['rad'] * conversion # km/s 

	#error_pm = x['parallax'] * x['fractional_parallax_error'] * 0.526 # pm error in mas/yr

	x = append_fields(x, 'radial_velocity', (radial_velocity), usemask = False)
	#x = append_fields(x, 'transverse_velocity', (transverse_velocity), usemask = False)
	#x = append_fields(x, 'pm_error', (error_pm), usemask = False)
	#x = append_fields(x, 'total_pm', (total_pm), usemask = False)
	x = append_fields(x, 'pm_ra', (pm_ra_cosdec), usemask = False)	
	x = append_fields(x, 'pm_dec', (pm_dec), usemask = False)		
	return(x)

