import numpy as np
from astropy.io import fits
import os,sys
i = int(sys.argv[1])
x = fits.getdata('../data/GDR3mock_207_%d.fits' %(i))
length = len(x)
print(i,length)
temp = x["source_id"][0]
counter = 0
for j in range(len(x)):
    if x["source_id"][j] == temp:
        x['source_id'][j] += counter
        counter += 1
    else:
        temp = x['source_id'][j]
        counter = 1    
fits.writeto('../data/GDR3mock207_%d.fits' %(i),x)
os.remove("../data/GDR3mock_207_%d.fits" %(i))
print("finished!")
