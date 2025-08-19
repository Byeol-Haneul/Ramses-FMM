This is the github repository for my new grafic IC generation scripts.
(Yes, I know that "grafic IC" is like saying ATM machine.)
I grew quite tired of using MUSIC for something that didn't need it, 
so I rewrote it all as just python scripts. Now, there's no need to deal with
f2py compilation or anything like that; everything's just in python. 

It should be fairly clear how to use the code. By default, if you simply write (e.g.):
./generateIC 8

you'll generate initial conditions for a uniform box with a resolution of (2^8)^3 = 256^3. 
By default (though this can be changed), it will be in a folder in the directory you ran the script in:
./ic_dturb/ic_dturb_8

Within generateIC, you can also specify the magnetic field strength. As well, if you want to modify the densities, velocities, etc.,
you can go into generate.py and make your own custom version of it. As it stands, the code just generates uniform ICs. 

You can test that the code works properly by using pip to install pytest (if you haven't already):
pip install pytest

Then you test the code within the repository folder with:
pytest -q

Hopefully, it passes the 11 tests!

Happy IC generation!
