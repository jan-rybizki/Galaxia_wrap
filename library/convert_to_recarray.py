import numpy as np

def create_gdr2mock_mag_limited_survey(nside = 16, 
                                       outputFile = "GDR2mock_20.7Gmag",
                                       modelFile = "Model/population_parameters_BGM_update.ebf",
                                       codeDataDir = "/home/rybizki/Programme/GalaxiaData",
                                       outputDir = '../output/mockcat_old',
                                       photoSys = "parsec1/GAIADR3",
                                       bandName = "gaia_g",
                                       colorName = "gaia_bpft-gaia_rp",
                                       appMagLimits0 = -1000, appMagLimits1 = 20.7,
                                       absMagLimits0 = -1000, absMagLimits1 = 1000,
                                       colorLimits0 = -1000, colorLimits1 = 1000,
                                       fSample = 0.0001, seed = 1, verbose = True,
                                       use_previous = False, delete_ebf = True,
                                       make_likelihood_asessment = False,
                                       warpFlareOn = 1, r_max = 1000):
    """
    This function creates the all sky mock survey catalogue via galaxia
    main input parameters are the Galaxia parameter file parameters
    additionally we can specify the nside of the extinction map and the all sky plot
    and also if additional output should be given via 'verbose'
    
    modelFile = "Model/population_parameters_mrtd5.ebf" # sharma+19 update
    modelFile = "Model/population_parameters_mptd.ebf" # old besancon
    modelFile = "Model/population_parameters_BGM_update.ebf" # my new
    """
    import subprocess, os, shutil, ebf, sys,time
    import numpy as np
    from healpy import ang2pix
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from numpy.lib.recfunctions import append_fields, rename_fields, drop_fields
    from astropy.io import fits
    path = os.path.abspath('../library/')
    if path not in sys.path:
        sys.path.append(path)
    from convert_to_recarray import convert
    from util import adding_velocities,av2ext, plot_sky_map, add_parsec_index, map_galaxia_indexes_onto_parsec, apply_extinction_curves
    
    start = time.time()
    magcolorNames = bandName + ',' + colorName
    # usually unchanged parameters
    geometryOption = 0 #int
    longitude = 0
    latitude = 90
    surveyArea = 1000
    popID = -1 #int
    starType = 0 #int
    photoError = 0 #int
    

    if not use_previous:
        # Creating the folder (will be overwritten if existed before)
        folder_create = outputDir + '/'
        if os.path.exists(folder_create):
            shutil.rmtree(folder_create)
            os.mkdir(folder_create)
            print(folder_create, "existed and was recreated")
        else:
            os.mkdir(folder_create)
        ###################
        # Creating the parameterfile
        filedata = 'outputFile                          %s\nmodelFile                          %s\ncodeDataDir                          %s\noutputDir                           %s\nphotoSys                            %s\nmagcolorNames                       %s\nappMagLimits[0]                     %f\nappMagLimits[1]                     %f\nabsMagLimits[0]                     %f\nabsMagLimits[1]                     %f\ncolorLimits[0]                      %f\ncolorLimits[1]                      %f\ngeometryOption                      %d\nlongitude                           %f\nlatitude                            %f\nsurveyArea                          %f\nfSample                             %f\npopID                               %d\nwarpFlareOn                         %d\nseed                                %d\nr_max                               %f\nstarType                            %d\nphotoError                          %d\n' %(outputFile,modelFile,codeDataDir,outputDir,photoSys,magcolorNames,appMagLimits0,appMagLimits1,absMagLimits0,absMagLimits1,colorLimits0,colorLimits1,geometryOption,longitude,latitude,surveyArea,fSample,popID,warpFlareOn,seed,r_max,starType,photoError)
        myparameterfile = outputDir + '/' + outputFile + '.log'
        file = open(myparameterfile, "w")
        file.write(filedata)
        file.close()

        # Creating mock catalogue
        args = ['galaxia', '-r', myparameterfile]
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)	
        print("Galaxia spawns catalogue")
        (output, err) = p.communicate()
        #p_status = p.wait()
        if verbose:
            print(output)
    
    print("########################################################################################")
    print("############################# GALAXIA OUTPUT END ##################")
    print("########################################################################################")
    # Reading in the Catalogue and converting it to npy file. (If you added photometric bands then you will have to edit the convert routine)
    start1 = time.time()
    bes = ebf.read(outputDir + "/" + outputFile + ".ebf",'/')
    x = convert(bes)
    print(len(x))
    print(x.dtype.names)
    # We add ra and dec
    x = append_fields(x,('ra','dec'),(np.zeros(len(x)),np.zeros(len(x))),usemask=False)
    c = SkyCoord(b=x['glat'], l=x['glon'], frame = 'galactic', unit=(u.deg, u.deg))
    x['ra'] = c.icrs.ra.deg
    x['dec'] = c.icrs.dec.deg
    
    conv = time.time()
    print("converting to npy and appending ra and dec took %.1f sec" %(conv-start1))
    ## Turning absolute into apparent magnitudes
    #filternames = ['gaia_g','gaia_gbp','gaia_grp']	
    filternames = ['gaia_g','gaia_bpft','gaia_bpbr','gaia_rp','gaia_rvs']
    for band in filternames:
        x[band] += 5 * np.log10(x['rad']*1000.) - 5

    use_bovy = False
    ## Optional you can use extinctions from Bovy 2016 python package mwdust (using Green2015, Marshall2006 and Schlegel)
    if use_bovy:
        from add_extinction import apply_bovy_extinction
        x = apply_bovy_extinction(x,nside)
    else:
        from add_extinction import apply_gdr3mock_extinction
        x = apply_gdr3mock_extinction(x,nside)
    exptime = time.time()
    print("converting time and applying extinction map for %d sources in nside = %d took %.1f sec" %(len(x),nside, exptime-conv))
       
    ## adding parsec index
    iso = np.load("../input/isochrones/parsec.npy")
    x = add_parsec_index(x,iso)
    x = map_galaxia_indexes_onto_parsec(x, iso)
    indextime = time.time()
    print("indexing and remapping to isochrones took %.1f sec" %(indextime-exptime))
    
    ## Adding extinction to apparent magnitudes
    x = append_fields(x,('a_g','a_bpft','a_bpbr','a_rp','a_rvs'),(np.zeros(len(x)),np.zeros(len(x)),np.zeros(len(x)),np.zeros(len(x)),np.zeros(len(x))),usemask=False)
    extratree_based = False
    if extratree_based:
        for band in filternames:
            ext_array = av2ext(band[5:],x["teff"],x["grav"],x["feh"],x["lum"],x["a0"])
            x[band] += ext_array
    else:
        ext_curve = np.load("../input/isochrones/ext.npy")
        indexing = np.searchsorted(iso['parsec_index'],x['parsec_index'])
        for i,band in enumerate(filternames):
            ext_array = apply_extinction_curves(x,indexing,ext_curve,band[5:])
            x[band] += ext_array
            ext_band = "a_%s" %(band[5:])
            x[ext_band] = ext_array
    extcurvetime = time.time()
    print("calculating extinction curve for all bands took %.1f sec" %(extcurvetime- indextime))
            
    # Now create a mag limited sample roughly three quarter of the sample will be lost due to extinction
    print(len(x))
    cut = (x[bandName]<appMagLimits1)
    x = x[cut]
    print(len(x))

    # clean fields
    x['teff'] = np.power(10,x['teff'])
    x = rename_fields(x, {'grav':'logg'})
    x = rename_fields(x, {'glon':'l'})
    x = rename_fields(x, {'glat':'b'})
    x = rename_fields(x, {'rad':'parallax'})
    x['parallax'] = np.divide(1.,x['parallax'])
    x['age'] = np.power(10,x['age'])
    # merge bpfaint and bpbright at G = 10.87
    bp = np.zeros_like(x['gaia_g'])
    a_bp = np.zeros_like(x['gaia_g'])
    bright_end = (x['gaia_g']<10.87)
    faint_end = (x['gaia_g']>=10.87)
    bp[bright_end] = x['gaia_bpbr'][bright_end]
    bp[faint_end] = x['gaia_bpft'][faint_end]
    a_bp[bright_end] = x['a_bpbr'][bright_end]
    a_bp[faint_end] = x['a_bpft'][faint_end]

    x['gaia_bpft'] = bp
    x['a_bpft'] = a_bp
    x = drop_fields(x, 'gaia_bpbr')
    x = drop_fields(x, 'a_bpbr')
    x = rename_fields(x, {'gaia_bpft':'phot_bp_mean_mag'})
    x = rename_fields(x, {'a_bpft':'a_bp'})
    x = rename_fields(x, {'gaia_g':'phot_g_mean_mag'})
    x = rename_fields(x, {'gaia_rp':'phot_rp_mean_mag'})
    x = rename_fields(x, {'gaia_rvs':'phot_rvs_mean_mag'})

    ## Nested healpix resolution of GAIA
    healpix_resolution = 4096
    number = ang2pix(healpix_resolution,x['ra'],x['dec'], nest = True, lonlat = True)
    number = number.view(np.int64)
    number *= 34359738368
    x = append_fields(x,'source_id',number,usemask = False)
    print('calculated healpix')
    x = adding_velocities(x)
    print('calculated pmdec pmra and rv')
    
    cleaningtime = time.time()
    print("cleaning of data took %.1f sec" %(cleaningtime - extcurvetime))
    
    
    
    # Save output
    if os.path.isfile(outputDir + "/" + outputFile + ".fits"):
        os.remove(outputDir + "/" + outputFile + ".fits")
    fits.writeto(outputDir + "/" + outputFile + ".fits",x)
    # Delete ebf file
    if delete_ebf:
        os.remove(outputDir + "/" + outputFile + ".ebf")
    
    ### Plot the stellar density across the sky
    if verbose:
        ## Using plotting routines from Tri Astraatmadja
        plot_sky_map(x,nside,fSample,appMagLimits1,outputDir,outputFile)
    plottingtime = time.time()
    print("plotting time took %.1f sec" %(plottingtime - cleaningtime))
    
    if make_likelihood_asessment:
        from cumhpx import plot_mollweide_log, delete_sources_below_maglim, likelihood_assessment
        nside_maglim = 128
        x = delete_sources_below_maglim(x,nside_maglim,'g',True)
        x = delete_sources_below_maglim(x,nside_maglim,'rp',True)
        cumhpx_like,hpx_like,hpx_dens,lengdr2data = likelihood_assessment(x,fSample)
        # add the numbers into the plot.
        plot_mollweide_log(-1*hpx_like,'-lnlike %.0f' %(sum(cumhpx_like)),savefile=True, name = outputDir + "/" + outputFile + "_likelihood.png", title = "gdr2: %d, mock: %d" %(lengdr2data,len(x)))
        likelihoodtime = time.time()
        print('assessing likelihood took %.1f sec' %(likelihoodtime-plottingtime))
        print("Total time in minutes: %.1f" %((likelihoodtime-start)/60))
    else:
        print("Total time in minutes: %.1f" %((plottingtime-start)/60))

def create_gdr2mock_mag_limited_survey_from_nbody(nbody_filename = 'nbody_file',nside = 16, 
                                       outputFile = "nbody",
                                       modelFile = "Model/population_parameters_BGM_update.ebf",
                                       codeDataDir = "/home/rybizki/Programme/GalaxiaData",
                                       outputDir = '../output/nbody',
                                       photoSys = "parsec1/GAIADR3",
                                       bandName = "gaia_g",
                                       colorName = "gaia_bpft-gaia_rp",
                                       appMagLimits0 = -1000, appMagLimits1 = 20.7,
                                       absMagLimits0 = -1000, absMagLimits1 = 1000,
                                       colorLimits0 = -1000, colorLimits1 = 1000,
                                       fSample = 1.0, seed = 1, verbose = True,
                                       use_previous = False, delete_ebf = True,
                                       make_likelihood_asessment = False,
                                       warpFlareOn = 1, r_max = 1000, popid=10):
    """
    This function creates the all sky mock survey catalogue via galaxia
    main input parameters are the Galaxia parameter file parameters
    additionally we can specify the nside of the extinction map and the all sky plot
    and also if additional output should be given via 'verbose'
    
    modelFile = "Model/population_parameters_mrtd5.ebf" # sharma+19 update
    modelFile = "Model/population_parameters_mptd.ebf" # old besancon
    modelFile = "Model/population_parameters_BGM_update.ebf" # my new
    """
    import subprocess, os, shutil, ebf, sys,time
    import numpy as np
    from healpy import ang2pix
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from numpy.lib.recfunctions import append_fields, rename_fields, drop_fields
    from astropy.io import fits
    path = os.path.abspath('../library/')
    if path not in sys.path:
        sys.path.append(path)
    from convert_to_recarray import convert
    from util import adding_velocities,av2ext, plot_sky_map, add_parsec_index, map_galaxia_indexes_onto_parsec, apply_extinction_curves
    
    start = time.time()
    magcolorNames = bandName + ',' + colorName
    # usually unchanged parameters
    geometryOption = 0 #int
    longitude = 0
    latitude = 90
    surveyArea = 1000
    popID = -1 #int
    starType = 0 #int
    photoError = 0 #int
    

    if not use_previous:
        # Creating the folder (will be overwritten if existed before)
        folder_create = outputDir + '/'
        if os.path.exists(folder_create):
            shutil.rmtree(folder_create)
            os.mkdir(folder_create)
            print(folder_create, "existed and was recreated")
        else:
            os.mkdir(folder_create)
        ###################
        # Creating the parameterfile
        filedata = 'outputFile                          %s\nmodelFile                          %s\ncodeDataDir                          %s\noutputDir                           %s\nphotoSys                            %s\nmagcolorNames                       %s\nappMagLimits[0]                     %f\nappMagLimits[1]                     %f\nabsMagLimits[0]                     %f\nabsMagLimits[1]                     %f\ncolorLimits[0]                      %f\ncolorLimits[1]                      %f\ngeometryOption                      %d\nlongitude                           %f\nlatitude                            %f\nsurveyArea                          %f\nfSample                             %f\npopID                               %d\nwarpFlareOn                         %d\nseed                                %d\nr_max                               %f\nstarType                            %d\nphotoError                          %d\n' %(outputFile,modelFile,codeDataDir,outputDir,photoSys,magcolorNames,appMagLimits0,appMagLimits1,absMagLimits0,absMagLimits1,colorLimits0,colorLimits1,geometryOption,longitude,latitude,surveyArea,fSample,popID,warpFlareOn,seed,r_max,starType,photoError)
        myparameterfile = outputDir + '/' + outputFile + '.log'
        file = open(myparameterfile, "w")
        file.write(filedata)
        file.close()

        # Creating mock catalogue
        args = ['galaxia', '-r', '--nfile=%s' %(nbody_filename), '--hdim=6', myparameterfile]
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)	
        print("Galaxia spawns catalogue")
        (output, err) = p.communicate()
        #p_status = p.wait()
        if verbose:
            print('output: ',output)
            print('error: ',err)
    
    print("########################################################################################")
    print("############################# GALAXIA OUTPUT END ##################")
    print("########################################################################################")
    # Reading in the Catalogue and converting it to npy file. (If you added photometric bands then you will have to edit the convert routine)
    start1 = time.time()
    bes = ebf.read(outputDir + "/" + outputFile + ".ebf",'/')
    x = convert(bes)
    print(len(x))
    print(x.dtype.names)
    # We add ra and dec
    x = append_fields(x,('ra','dec'),(np.zeros(len(x)),np.zeros(len(x))),usemask=False)
    c = SkyCoord(b=x['glat'], l=x['glon'], frame = 'galactic', unit=(u.deg, u.deg))
    x['ra'] = c.icrs.ra.deg
    x['dec'] = c.icrs.dec.deg
    
    conv = time.time()
    print("converting to npy and appending ra and dec took %.1f sec" %(conv-start1))
    ## Turning absolute into apparent magnitudes
    #filternames = ['gaia_g','gaia_gbp','gaia_grp']	
    filternames = ['gaia_g','gaia_bpft','gaia_bpbr','gaia_rp','gaia_rvs']
    for band in filternames:
        x[band] += 5 * np.log10(x['rad']*1000.) - 5

    use_bovy = False
    ## Optional you can use extinctions from Bovy 2016 python package mwdust (using Green2015, Marshall2006 and Schlegel)
    if use_bovy:
        from add_extinction import apply_bovy_extinction
        x = apply_bovy_extinction(x,nside)
    else:
        from add_extinction import apply_gdr3mock_extinction
        x = apply_gdr3mock_extinction(x,nside)
    exptime = time.time()
    print("converting time and applying extinction map for %d sources in nside = %d took %.1f sec" %(len(x),nside, exptime-conv))
       
    ## adding parsec index
    iso = np.load("../input/isochrones/parsec.npy")
    x = add_parsec_index(x,iso)
    x = map_galaxia_indexes_onto_parsec(x, iso)
    indextime = time.time()
    print("indexing and remapping to isochrones took %.1f sec" %(indextime-exptime))
    
    ## Adding extinction to apparent magnitudes
    x = append_fields(x,('a_g','a_bpft','a_bpbr','a_rp','a_rvs'),(np.zeros(len(x)),np.zeros(len(x)),np.zeros(len(x)),np.zeros(len(x)),np.zeros(len(x))),usemask=False)
    extratree_based = False
    if extratree_based:
        for band in filternames:
            ext_array = av2ext(band[5:],x["teff"],x["grav"],x["feh"],x["lum"],x["a0"])
            x[band] += ext_array
    else:
        ext_curve = np.load("../input/isochrones/ext.npy")
        indexing = np.searchsorted(iso['parsec_index'],x['parsec_index'])
        for i,band in enumerate(filternames):
            ext_array = apply_extinction_curves(x,indexing,ext_curve,band[5:])
            x[band] += ext_array
            ext_band = "a_%s" %(band[5:])
            x[ext_band] = ext_array
    extcurvetime = time.time()
    print("calculating extinction curve for all bands took %.1f sec" %(extcurvetime- indextime))
            
    # Now create a mag limited sample roughly three quarter of the sample will be lost due to extinction
    print(len(x))
    cut = (x[bandName]<appMagLimits1)
    x = x[cut]
    print(len(x))

    # clean fields
    x['teff'] = np.power(10,x['teff'])
    x = rename_fields(x, {'grav':'logg'})
    x = rename_fields(x, {'glon':'l'})
    x = rename_fields(x, {'glat':'b'})
    x = rename_fields(x, {'rad':'parallax'})
    x['parallax'] = np.divide(1.,x['parallax'])
    x['age'] = np.power(10,x['age'])
    # merge bpfaint and bpbright at G = 10.87
    bp = np.zeros_like(x['gaia_g'])
    a_bp = np.zeros_like(x['gaia_g'])
    bright_end = (x['gaia_g']<10.87)
    faint_end = (x['gaia_g']>=10.87)
    bp[bright_end] = x['gaia_bpbr'][bright_end]
    bp[faint_end] = x['gaia_bpft'][faint_end]
    a_bp[bright_end] = x['a_bpbr'][bright_end]
    a_bp[faint_end] = x['a_bpft'][faint_end]

    x['gaia_bpft'] = bp
    x['a_bpft'] = a_bp
    x = drop_fields(x, 'gaia_bpbr')
    x = drop_fields(x, 'a_bpbr')
    x = rename_fields(x, {'gaia_bpft':'phot_bp_mean_mag'})
    x = rename_fields(x, {'a_bpft':'a_bp'})
    x = rename_fields(x, {'gaia_g':'phot_g_mean_mag'})
    x = rename_fields(x, {'gaia_rp':'phot_rp_mean_mag'})
    x = rename_fields(x, {'gaia_rvs':'phot_rvs_mean_mag'})
    x['popid'] = popid

    ## Nested healpix resolution of GAIA
    healpix_resolution = 4096
    number = ang2pix(healpix_resolution,x['ra'],x['dec'], nest = True, lonlat = True)
    number = number.view(np.int64)
    number *= 34359738368
    x = append_fields(x,'source_id',number,usemask = False)
    print('calculated healpix')
    x = adding_velocities(x)
    print('calculated pmdec pmra and rv')
    
    cleaningtime = time.time()
    print("cleaning of data took %.1f sec" %(cleaningtime - extcurvetime))
    
    
    
    # Save output
    if os.path.isfile(outputDir + "/" + outputFile + ".fits"):
        os.remove(outputDir + "/" + outputFile + ".fits")
    fits.writeto(outputDir + "/" + outputFile + ".fits",x)
    # Delete ebf file
    if delete_ebf:
        os.remove(outputDir + "/" + outputFile + ".ebf")
    
    ### Plot the stellar density across the sky
    if verbose:
        ## Using plotting routines from Tri Astraatmadja
        plot_sky_map(x,nside,fSample,appMagLimits1,outputDir,outputFile)
    plottingtime = time.time()
    print("plotting time took %.1f sec" %(plottingtime - cleaningtime))
    
    if make_likelihood_asessment:
        from cumhpx import plot_mollweide_log, delete_sources_below_maglim, likelihood_assessment
        nside_maglim = 128
        x = delete_sources_below_maglim(x,nside_maglim,'g',True)
        x = delete_sources_below_maglim(x,nside_maglim,'rp',True)
        cumhpx_like,hpx_like,hpx_dens,lengdr2data = likelihood_assessment(x,fSample)
        # add the numbers into the plot.
        plot_mollweide_log(-1*hpx_like,'-lnlike %.0f' %(sum(cumhpx_like)),savefile=True, name = outputDir + "/" + outputFile + "_likelihood.png", title = "gdr2: %d, mock: %d" %(lengdr2data,len(x)))
        likelihoodtime = time.time()
        print('assessing likelihood took %.1f sec' %(likelihoodtime-plottingtime))
        print("Total time in minutes: %.1f" %((likelihoodtime-start)/60))
    else:
        print("Total time in minutes: %.1f" %((plottingtime-start)/60))


def convert(jj):
	assert len(jj)!= 0	
	entries = len(jj['rad'])
	dt = np.dtype([('rad',float), ('teff',float), ('vx',float), ('vy',float), ('vz',float), ('pz',float), ('px',float), ('py',float), ('feh',float), ('exbv_schlegel',float), ('lum',float), ('glon',float), ('glat',float), ('smass',float), ('age',float), ('grav',float), ('gaia_g',float), ('gaia_bpft',float), ('gaia_bpbr',float),
('gaia_rp',float), ('gaia_rvs',float), ('popid',np.int32), ('mact',float)])

	x = np.empty((entries,), dtype = dt)
	for i,item in enumerate(jj.keys()):
		if item not in ['log','center','label','mag1','mag2', 'mag0','fieldid','satid','mtip','partid','exbv_solar', 'exbv_schlegel_inf', 'alpha']:
			x[item] = jj[item]
	return(x)
