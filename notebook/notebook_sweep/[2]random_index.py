import numpy as np
from astropy.io import fits
import os
from numpy.lib.recfunctions import stack_arrays, append_fields

count = np.load('total_starcount.npy')
# generated on think-rybizki

"""
print(count)
indexes = np.arange(1,count+1,dtype = np.int32)
ri = np.random.choice(indexes,size = count, replace=False)
"""
ri = np.load('random_index.npy')
assert count == len(ri)

for i in np.arange(400):
    x = fits.getdata('../data/%d.fits' %(i))
    length = len(x)
    print(i,length, len(ri))
    x = append_fields(x,'random_index',ri[:length],usemask = False)
    ri = ri[length:]

    fits.writeto('../data/GDR3mock_207_%d.fits' %(i),x)
    os.remove("../data/%d.fits" %(i))
print(ri)
print("finished!")
