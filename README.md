# Galaxia_wrap

A few python scripts and a tutorial notebook to create and manipulate Galaxia mock catalogues


**What is Galaxia**

Galaxia is a software to create mock stellar catalogues [http://adsabs.harvard.edu/abs/2011ApJ...730....3S]. You can plug in n-body simulations (which we will do in another notebook) or use a Galaxy model. Galaxia uses the old 2003 Besancon Model [http://adsabs.harvard.edu/abs/2003A%26A...409..523R] and populates stars accordingly using Padova Isochrones [http://adsabs.harvard.edu/abs/2008A%26A...482..883M]. You can add extinctions that are coming from a 3D extinction map [http://adsabs.harvard.edu/abs/2003A%26A...409..205D].

**Prerequisits**

Download Galaxia [https://sourceforge.net/projects/galaxia/files/] and follow the instructions [http://galaxia.sourceforge.net/Galaxia3pub.html]. You will also have to install ebf into your python distribution via

```
pip install ebfpy
```

Before creating your first mock catalogue you will need to run

```
galaxia -s warp
```
which will run for some time.

**Start with the tutorial**

Now you can dive into the tutorials which will show you (a) how to create a magnitude limited catalogue and (b) how to input your n-body data to create mock observations.

This tutorial was setup for the Gaia Sprint 2017 [http://gaia.lol/].
