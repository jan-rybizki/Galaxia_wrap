
# coding: utf-8

# In[1]:


import numpy as np
import os, sys
path = os.path.abspath('../library/')
if path not in sys.path:
    sys.path.append(path)
from convert_to_recarray import create_gdr2mock_mag_limited_survey


# In[ ]:

seed = int(sys.argv[1])
create_gdr2mock_mag_limited_survey(nside = 512, outputDir = '/home/rybizki/Desktop/Galaxia_wrap-master/output/GDR3mock/%d' %(seed),outputFile = "GDR3mock207",
                                   use_previous = False, delete_ebf = True,
                                  fSample = 0.001, make_likelihood_asessment=False, seed = seed)

