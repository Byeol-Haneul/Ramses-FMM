[1]: https://bitbucket.org/rteyssie/ramses/
[2]: https://bitbucket.org/rteyssie/mini-ramses/

# mini-ramses-AGN
A repository for developing sink particle / black hole / AGN models for mini-ramses.
The original mini-ramses trunk is available on bitbucket at: https://bitbucket.org/rteyssie/mini-ramses/src/develop/

## Here is the to-do list for sink particle development:
### Sink Refine
- [ ] Add sink-refine
- [ ] Add ability to modify the size of the sink particle region based on the mass of the sink particle

### Dyanmics
- [ ]  Under construction

### Accretion
- [ ] Add simple bondi accretion
- [ ] Add a weighting scheme (mass, volume, gaussian kernel)
- [ ] Add flux accretion method (Bleuler+2014)
- [ ] Add a time-step mechanism, such that the accretion is limited
- [ ] Design a new angular momentum-based weighting scheme, (e.g. Hopkins,Quataert 2011, Angles-Alcazar+2013)
- [ ] Investigate the need for a sector based approach, (see Talbot+21, Fiacconi+18)

### Feedback
- [ ] Implement a simple version of two-mode AGN feedback (Teyssier+2011, Dubois+2012)
- [ ] Implement a different weighting scheme for feedback and accretion
- [ ] Play around with having a different (larger) feedback region than accretion region

### Internal Black Hole evolution
- [ ] Finalise the accretion disk (AD) and black hole (BH) internal model, (Talbot+21)
- [ ] Design an efficient workflow for all of these processes
- [ ] Implement a first version
- [ ] Include a switch between the 'simple' BH model and this one
- [ ] Leave room for other physics (e.g. disk alpha, R/H, phi_BH to be dynamically updated, see Chris Bambic's paper suggestions)

### Couple the Internal Black hole model to a Blandford-Znajek Jet
- [ ] Understand exactly what's happening in the Talbot+21 model
- [ ] Implement into the code, adding a switch between this feedback mode and the original two-mode model
- [ ] find an optimum injection geometry

### Couple the Black Hole to RT
- [ ] Settle on a full model to be used (e.g. Trebitsch+19,21,23)
- [ ] Implement this model
- [ ] Couple to the rest of feedback, ensure energy is conserved


## Intended tests
- [ ] Sink refine: Orbiting sink particle with refining mesh
- [ ] Accretion: Typical Bondi-Hoyle-Lyttleton tests
- [ ] Accretion: Idealised disk accretion
- [ ] Feedback: Compare jet/quasar mode properties to RAMSES OG
- [ ] Advanced Feedback: Same as above, redo some of the experiments from Talbot+21
- [ ] RT Feedback: Idealised RT feedback with no jet, test if this produces a good Narrow Line Region (requires RTZ)
- [ ] RT Feedback: Same jet tests as Talbot+21, compare the loading factors with the non-RT version above

# Original README:
## mini-ramses ##

The mini-ramses repository is a fork of the [main RAMSES repository][1]. It was created as a stripped-down version of the main code base created in order to facilitate the development of major updates of RAMSES' core routines. This stripped-down version is still available in the `master` branch.

A rapidly-evolving development of mini-ramses based on a completely changed data-structure is contained in the `develop` branch. This branch is work in progress, use it only if you know what you're doing ;)

You can download the code by cloning the git repository using 
```
$ git clone https://bitbucket.org/rteyssie/mini-ramses.git
```

If you want to contribute to mini-ramses, you can either ask me (romain.teyssier@gmail.com) for your personal new branch in this repository which I will give you write access to, or you can fork this repository. To bring changes back into the `develop` branch of mini-ramses, simply issue a pull request.

To compile and execute the standard test cases, please follow these steps.

1- Shock tube test in 1D:

```
$ cd bin
$ make clean
$ make NDIM=1 HYDRO=1
$ cd ..
$ bin/ramses1d namelist/tube1d.nml
```

2- Sedov explosion in 2D:

```
$ cd bin
$ make clean
$ make NDIM=2 HYDRO=1
$ cd ..
$ bin/ramses2d namelist/sedov2d.nml
```

3- Magnetic loop advection in 2D:

```
$ cd bin
$ make clean
$ make NDIM=2 HYDRO=1 MHD=1 INIT=LOOP
$ cd ..
$ bin/ramses2d namelist/loop.nml
```

4- Sedov explosion test in 3D:

```
$ cd bin
$ make clean
$ make NDIM=3 HYDRO=1
$ cd ..
$ bin/ramses3d namelist/sedov3d.nml
```

5- Cosmological N body simulation in 3D

```
$ cd bin
$ make clean
$ make NDIM=3 HYDRO=0 GRAV=1 UNITS=COSMO
$ cd ..
$ utils/script/load_cosmo_ic.sh
$ bin/ramses3d namelist/dmo.nml
```

6- Molecular core test in 3D:

```
$ cd bin
$ make clean
$ make NDIM=3 HYDRO=1 GRAV=1 UNITS=COEUR INIT=COEUR
$ cd ..
$ bin/ramses3d namelist/coeur.nml
```

You get the picture now ;-)

To visualize the 2D and 3D results, compile the map making executable in the utils/f90 directory.

```
$ cd utils/f90
$ gfortran amr2map.f90 -o amr2map
$ cd ../..
$ utils/f90/amr2map -inp output_00002 -out dens.map -typ 1
$ utils/py/map2img.py dens.map --log
```

In the molecular cloud collapse case, you can also explore the movie1 directory and use the python function directly on any of the maps in there. 
If you have the MPI library properly installed on your system, you can repeat all the tests above using the MPI=1 compilation option.

