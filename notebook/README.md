## Galaxia wrap notebooks roadmap

You can load the notebooks using ``jupyter notebook`` and then change the content interactively, once you downloaded this repository.
You can also just have a look at them here on github.

**Legacy GDR2mock notebooks**

------------------------------------

In **[notebook 1](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/%5B1%5DCreating%20magnitude%20limited%20synthetic%20catalogues%20with%20Galaxia.ipynb)** we will use Galaxia to produce a magnitude limited survey of the Milky Way Galaxy.

**[Notebook 2](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/%5B2%5DProviding%20n-body%20particle%20input%20for%20Galaxia.ipynb)** is about the n-body input to Galaxia. We will not use real n-body simulations but just create a CMD by producing SSPs with a constant SFR. The output can then be used to check e.g. the age distribution of specific stellar tracer population.

------------------------------------

**New notebooks in order to reproduce GeDR3mock results**

**[Notebook 3](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/%5B3%5DComparing%20Catalogs%20CMDs%20per%20healpix.ipynb)** Can be used to create a likelihood function based on CMDs per HEALpix.

**[Notebook 4](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/%5B4%5Dmag_limited_survey_function.ipynb)** Shows the different functionalities of the wrapper function that creates magnitude limited all-sky surveys.

**[Notebook 5](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/%5B5%5D%20parameter%20search.ipynb)([a](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/%5B5a%5Dinspect_parameter_search.ipynb))** Galactic parameter search using Galaxia as a forward model (notebook 4) together with different likelihoods (notebook 3)

**[Notebook 6](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/%5B6%5DProviding%20n-body%20particle%20input%20of%20SMC%20LMC%20for%20GDR3mock.ipynb)** How to create n-body input for the SMC LMC.

**[Notebook 7](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/%5B7%5Dcluster%20mockup.ipynb)([a](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/%5B7a%5Dprepare_cluster_catalog.ipynb))** Mocking up cluster data (based on real observations) and interfacing it with the galaxia n-body functionality in order to produce mock catalogs of the clusters.

**[Notebook 8](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/%5B8%5DTrain%20parallax%20and%20proper%20motion%20error_alternative_projection.ipynb)([a](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/%5B8a%5DCalculate%20the%20error%20relations.ipynb))** Empirical training of the parallax uncertainty and relations for proper motion error. 

**[isochrone_generation](https://github.com/jan-rybizki/Galaxia_wrap/tree/master/notebook/isochrone_generation)** Here the conversion and filtering of the raw stellar evolutionary tracks to the final isochrone tables can be inspected.

**[notebook_sweep](https://github.com/jan-rybizki/Galaxia_wrap/tree/master/notebook/notebook_sweep)** shows the final compilation of GeDR3mock for upload to the Virtual Observatory

**[notebook_aida](https://github.com/jan-rybizki/Galaxia_wrap/tree/master/notebook/notebook_aida)** is similar to the notebooks here but also includes a notebook which computes the new 3d extinction map used when generating GeDR3mock.
