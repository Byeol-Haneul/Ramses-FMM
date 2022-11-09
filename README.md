[1]: https://bitbucket.org/rteyssie/ramses/
[2]: https://bitbucket.org/rteyssie/mini-ramses/

## mini-ramses ##

The mini-ramses repository is a fork of the [main RAMSES repository][1]. It was created as a stripped-down version of the main code base created in order to facilitate the development of major updates of RAMSES' core routines. This stripped-down version is still available in the `master` branch.

A rapidly-evolving development of mini-ramses based on a completely changed data-structure is contained in the `develop` branch. This branch is work in progress, use it only if you know what you're doing ;)


You can download the code by cloning the git repository using 
```
$ git clone https://bitbucket.org/rteyssie/mini-ramses.git
```

If you want to contribute to mini-ramses, you can either ask me (romain.teyssier@gmail.com) for your personal new branch in this repository which I will give you write access to, or you can fork this repository. To bring changes back into the `develop` branch of mini-ramses, simply issue a pull request.

To compile and execute the standard test cases, please follow these steps.

1- Molecular core test in 3D:

```
$ cd bin
$ make clean
$ make
$ cd ..
$ bin/ramses3d namelist/coeur.nml
```

2- Sedov explosion test in 3D:

```
$ cd bin
$ make clean
$ make UNITS=default CONDINIT=default
$ cd ..
$ bin/ramses3d namelist/sedov3d.nml
```

3- Sedov explosion in 2D:

```
$ cd bin
$ make clean
$ make NDIM=2 UNITS=default CONDINIT=default
$ cd ..
$ bin/ramses2d namelist/sedov2d.nml
```

4- Shock tube test in 1D:

```
$ cd bin
$ make clean
$ make NDIM=1 UNITS=default CONDINIT=default
$ cd ..
$ bin/ramses1d namelist/tube1d.nml
```

5- Cosmological N body simulation in 3D

```
$ cd bin
$ make clean
$ make UNITS=cosmo CONDINIT=default
$ cd ..
$ utils/script/load_cosmo_ic.sh
$ bin/ramses3d namelist/dmo.nml
```

You get the picture now ;-)

To visualize the 2D (resp. 3D) results, compile the map making executable in the bin directory right after you compiled ramses using NDIM=2 (resp. 3).

```
$ cd bin
$ make NDIM=2 amr2map
$ cd ..
$ bin/amr2map2d -inp output_00002 -out dens,map
$ utils/py/map2img.py dens.map --log
```

In the molecular cloud collapse case, you can also explore the movie1 directory and use the python function directly on any of the maps in there.

