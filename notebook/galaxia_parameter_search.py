
# coding: utf-8

# In[1]:


import numpy as np
import os, sys
path = os.path.abspath('../library/')
if path not in sys.path:
    sys.path.append(path)
from convert_to_recarray import create_gdr2mock_mag_limited_survey
import ebf, subprocess
from cumhpx import do_likelihood, local_normalisation


# In[2]:


### change files
'''
default start = 0
default end = 10
0-6 == thin disk
7 = thick disk
8 = halo
9 = bulge
what can be changed:
/norm_thin_disk
/norm_thick_disk
/norm_halo
/norm_bulge
'''
modelFile = '/home/rybizki/Programme/GalaxiaData/Model/population_parameters_BGM_update.ebf'   
thin_norms = np.linspace(0.6,1.0,41)

thick_norms = [1.0,1.15,1.3,1.45,1.6,1.75,1.9,2.05]
result = np.zeros(shape = (len(thin_norms),len(thick_norms),7))
for i,thin_norm in enumerate(thin_norms):
    print(thin_norm)
    print('#######################################')
    ebf.update_ind(modelFile, '/popidstart', 0, ind=0)
    ebf.update_ind(modelFile, '/popidend', 7, ind=0)
    ebf.update_ind(modelFile, '/norm_thin_disk', thin_norm, ind=0)
    ### make new bhtree
    args = ['galaxia', '-s', 'warp']
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)	
    print("Galaxia spawns new thin disk BHtree")
    (output, err) = p.communicate()
    if True:
        print(output)

    for j,thick_norm in enumerate(thick_norms):
        print(thick_norm)
        print('#######################################')
        ebf.update_ind(modelFile, '/popidstart', 7, ind=0)
        ebf.update_ind(modelFile, '/popidend', 8, ind=0)
        ebf.update_ind(modelFile, '/norm_thick_disk', thick_norm, ind=0)
        ### make new bhtree
        args = ['galaxia', '-s', 'warp']
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)	
        print("Galaxia spawns new thick disk BHtree")
        (output, err) = p.communicate()
        if True:
            print(output)

        ### calculate mock catalogue
        name = 'raster_%.1f_%.1f' %(thin_norm,thick_norm)
        create_gdr2mock_mag_limited_survey(nside = 512, outputDir = '../output/%s_100pc' %(name),
                                           use_previous = False, delete_ebf = True,
                                          fSample = 1, make_likelihood_asessment=False, r_max = 0.1,
                                          verbose = False)
        create_gdr2mock_mag_limited_survey(nside = 512, outputDir = '../output/%s_0.001' %(name),
                                           use_previous = False, delete_ebf = True,
                                          fSample = 0.001, make_likelihood_asessment=False)
        ### get likelihoods
        stellardiskmass,wddiskmass,totalcount = local_normalisation(name, verbose = False)
        outputDir = '../output/%s_0.001' %name
        fSample = 0.001
        starcounts,likelihood,diffplus,diffminus = do_likelihood(outputDir,fSample= fSample)
        result[i,j,0] = stellardiskmass
        result[i,j,1] = wddiskmass
        result[i,j,2] = totalcount
        result[i,j,3] = starcounts
        result[i,j,4] = likelihood
        result[i,j,5] = diffplus
        result[i,j,6] = diffminus

np.save('result.npy',result)
print(result)
