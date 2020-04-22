## Galaxia wrap notebooks roadmap

You can load the notebooks using ``jupyter notebook`` and then change the content interactively, once you downloaded this repository.
You can also just have a look at them here on github.

------------------------------------

**New notebooks in order to generate the isochrone set for GeDR3mock mock**

**[Notebook 0](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/isochrone_generation/notebook/%5B0%5Dextract%20isochrones.ipynb)** shows how to convert the raw stellar evolutionary tracks into npy tables.

**[Notebook 1](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/isochrone_generation/notebook/%5B1%5Ddelete_bad_extinctions_from_isochrones.ipynb)** shows how we deleted stellar models with too high extinction values which compromised the grid interpolation.

**[Notebook 2new](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/isochrone_generation/notebook/%5B2old%5Dcreate_parsec_index_isochrones_fast.ipynb)** shows how to index the isochrone tables and generates the isochrone grid (for different photometric systems).

**[Notebook 3new](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/isochrone_generation/notebook/%5B3new%5Dinspect_parsec_index.ipynb)** evaluates and checks the resulting grid.

**[Notebook 4](https://github.com/jan-rybizki/Galaxia_wrap/blob/master/notebook/isochrone_generation/notebook/%5B4%5Dcreate%20parsec_and_ext_for_upload.ipynb)** prepares the isochrone grid for upload to Virtual Observatory.
