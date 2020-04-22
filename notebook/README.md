## Galaxia wrap notebooks roadmap

You can load the notebooks using ``jupyter notebook`` and then change the content interactively, once you downloaded this repository.
You can also just have a look at them here on github.

**Legacy GDR2mock notebooks**

------------------------------------

In **notebook 1** we will use Galaxia to produce a magnitude limited survey of the Milky Way Galaxy.

**Notebook 2** is about the n-body input to Galaxia. We will not use real n-body simulations but just create a CMD by producing SSPs with a constant SFR. The output can then be used to check e.g. the age distribution of specific stellar tracer population.

------------------------------------

**New notebooks in order to reproduce GeDR3mock results**

**Notebook 3** Can be used to create a likelihood function based on CMDs per HEALpix.

**[Notebook 4](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/%5B4%5Dmag_limited_survey_function.ipynb)** Shows the different functionalities of the wrapper function that creates magnitude limited all-sky surveys.

**Notebook 5(a)** Galactic parameter search using Galaxia as a forward model (notebook 4) together with different likelihoods (notebook 3)

**Notebook 6** How to create n-body input for the SMC LMC.

**Notebook 7(a)** Mocking up cluster data (based on real observations) and interfacing it with the galaxia n-body functionality in order to produce mock catalogs of the clusters.

**Notebook 8(a)** Empirical training of the parallax uncertainty and relations for proper motion error. 
