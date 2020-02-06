import numpy as np
from astropy.io import fits
from sklearn.ensemble import ExtraTreesRegressor
import pickle
from numpy.lib.recfunctions import append_fields
import os
from time import sleep

nobs_model = pickle.load(open('errors/lb2vpunobs_model_bigger','rb'))
pege_model = pickle.load(open('errors/gbprpvpunobs2pege_model_bigger','rb'))
rvs_model = pickle.load(open('errors/gbprpteff2rvse_model','rb'))
for i in range(400):
    print(i,400)
    x = fits.getdata("../data/GDR3mock207_%d.fits" %(i))
    # train nobs
    l = x.l
    b = x.b
    X = np.vstack((l,b)).T
    y = nobs_model.predict(X)
    x = append_fields(x,["visibility_periods_used","phot_g_n_obs"],(np.round(y[:,0]*(34/22)),np.round(y[:,1]*(34/22))), usemask = False, asrecarray = True)
    # train parallax and gmag error
    g = x.phot_g_mean_mag
    bprp = x.phot_bp_mean_mag - x.phot_rp_mean_mag
    vpu = x.visibility_periods_used
    nobs = x.phot_g_n_obs
    X = np.vstack((g,bprp,vpu,nobs)).T
    y = pege_model.predict(X)
    x = append_fields(x,["parallax_error","phot_g_mean_mag_error"],(y[:,0]*np.sqrt(22/34),y[:,1]*np.sqrt(22/34)), usemask = False, asrecarray = True)
    # train rvs error
    teff = x.teff
    X = np.vstack((g,bprp,teff)).T
    y = rvs_model.predict(X)
    x = append_fields(x,"radial_velocity_error",np.hstack(y)*np.sqrt(22/34), usemask = False, asrecarray = True)
    fits.writeto("../data/GDR3mock_207_%d.fits" %(i),x)
    os.remove("../data/GDR3mock207_%d.fits" %(i))

