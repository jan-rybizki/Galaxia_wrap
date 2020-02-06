import numpy as np
from astropy.io import fits
import os
from numpy.lib.recfunctions import stack_arrays, drop_fields
from time import sleep

fields_to_drop = ['px','py','pz','vx','vy','vz']
edges = np.load('edges.npy')

#MW
for i in np.arange(1000):
    x = fits.getdata('../output/GDR3mock/%d/GDR3mock207.fits' %(i) )
    x = drop_fields(x,fields_to_drop,usemask= False, asrecarray = True)
    print(i,len(x))
    x=np.sort(x,order='source_id')
    indexes = np.searchsorted(x.source_id,edges)
    for j, jtem in enumerate(indexes[:-1]):
        np.save('../data/%d_%d.npy' %(i,j) ,x[jtem:indexes[j+1]])
    os.remove('../output/GDR3mock/%d/GDR3mock207.fits' %(i) )
	
#MC
for i in np.arange(10):
    x = fits.getdata('../output/GDR3mock_extra/MCs_%d/nbody.fits' %(i) )
    x = drop_fields(x,fields_to_drop,usemask= False, asrecarray = True)
    print(i,len(x))
    x=np.sort(x,order='source_id')
    indexes = np.searchsorted(x.source_id,edges)
    for j, jtem in enumerate(indexes[:-1]):
        np.save('../data/%d_%d.npy' %(i+1000,j) ,x[jtem:indexes[j+1]])

#CL
for i in np.arange(1):
    x = fits.getdata('../output/GDR3mock_extra/Clusters_1/nbody.fits')
    x = drop_fields(x,fields_to_drop,usemask= False, asrecarray = True)
    print(i,len(x))
    x=np.sort(x,order='source_id')
    indexes = np.searchsorted(x.source_id,edges)
    for j, jtem in enumerate(indexes[:-1]):
        np.save('../data/%d_%d.npy' %(i+1010,j) ,x[jtem:indexes[j+1]])


count = 0
for i in np.arange(400):
    array_stack = []
    for j in np.arange(1011):
        array_stack.append(np.load('../data/%d_%d.npy' %(j,i)))
    x = stack_arrays(array_stack, usemask = False, asrecarray=True,autoconvert = True)
    for j in np.arange(1011):
        os.remove('../data/%d_%d.npy' %(j,i))    
    count+=len(x)
    print(i,len(x))
    x=np.sort(x,order = 'source_id', kind = 'mergesort')
    fits.writeto('../data/%d.fits' %(i),x)
np.save('total_starcount.npy', count)
    

