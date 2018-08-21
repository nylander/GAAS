
[<img align="center" src="NBIS.png" width="200" height="100" />](https://nbis.se) 
<h2 ><em>G</em>enome <em>A</em>ssembly <em>A</em>nnotation <i>S</i>ervice (GAAS)</h2>

---------------------------

Contains development done in the GAA (Genome Assembly Annotation) Service.

[__annotation__](annotation) directory contains development related to annotation side of the service.</br>

To use the tools in this repository, clone the directory and run `make install` to update some paths:
```
git clone https://github.com/NBISweden/GAAS.git
cd GAAS
make install   # Updates paths in the environment profiles to point to the correct GAAS repository location
```

Two profiles are available to setup the necessary environment variables to use the scripts:

If you are on Rackham:
```
source GAAS/profiles/activate_rackham_env
```

If you are on the NBIS production server:
```
source GAAS/profiles/activate_nbis_env
```
