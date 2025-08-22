This folder contains grafic IC generation scripts, courtesy of Eric Moseley (emoseley@stanford.edu). 
(Yes, I know that "grafic IC" is like saying ATM machine.)

Included in this directory are several examples of IC generation scripts, which can
then be modified for your own purposes. Read the code to learn more, but put simply, if I want to generate e.g. a box with decaying turbulence in it with a resolution of e.g. 2^6 in 2D, with modes ranging from kmin=1 to kmax=2, with a 50/50 mix of solenoidal and compressive modes, an initial r.m.s. velocity of 1.0, and a magnetic field strength of 1.0 in the y direction, then I would write:

python turb.py 6 --ndim 2 --kmin 1 --kmax 2 --alpha 0.5 --vrms 1.0 --by 1.0

You would then see a folder created where you ran the script (unless you specify another output directory):

ic_turb/ic_turb_6_2d

The other two scripts, khi.py and uniform.py work similarly. khi.py will generate ICs for a Kelvin-Helmholtz intsability, while uniform.py generates uniform initial conditions. Both uniform.py and turb.py include the option to add particles to the ICs.

You can test that the code works properly by using pip to install pytest (if you haven't already):
pip install pytest

Then you test the code within the repository folder with:
pytest -q

Hopefully, it passes the 11 tests!

Happy IC generation!
