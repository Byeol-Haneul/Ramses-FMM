#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import argparse
import miniramses as ram

parser = argparse.ArgumentParser()
parser.add_argument("nout", help="enter output number")
parser.add_argument("--path", help="specify a path")
parser.add_argument("--log", help="plot log variable",action="store_true")
parser.add_argument("--out", help="output a png image")
parser.add_argument("--pref", help="specify a file prefix")
parser.add_argument("--min", help="specify a minimum variable value")
parser.add_argument("--var", help="specify a variable number")
parser.add_argument("--xcen", help="specify the image center x-coordinate")
parser.add_argument("--ycen", help="specify the image center y-coordinate")
parser.add_argument("--zcen", help="specify the image center z-coordinate")
parser.add_argument("--rad", help="specify the image radius")
parser.add_argument("--clump", help="specify the image radius")
args = parser.parse_args()
# path the the file
path = args.path
prefix = args.pref
ivar = args.var
vmin = args.min
radius = args.rad
xcenter = args.xcen
ycenter = args.ycen
zcenter = args.zcen
clump = args.clump
if clump==None:
    clump=False
if xcenter==None:
    xcenter=0.5
else:
    xcenter=float(xcenter)
if ycenter==None:
    ycenter=0.5
else:
    ycenter=float(ycenter)
if zcenter==None:
    zcenter=0.5
else:
    zcenter=float(zcenter)
if radius==None:
    radius=0.1
else:
    radius=float(radius)
center=np.array([xcenter,ycenter,zcenter])

if ivar==None:
    ivar=0
else:
    ivar=int(ivar)-1
if prefix==None:
    prefix="hydro"
if path==None:
    path="./"
else:
    path=path+"/"
if prefix=="hydro":
    isort=0
if prefix=="peak":
    isort=1
if prefix=="grav":
    isort=0

nout = args.nout
print("Reading output number ",nout)
print(path)

print(center)
print(radius)

c=ram.rd_cell(nout,path=path,prefix=prefix,center=center,radius=radius)
ram.visu(c.x[0],c.x[1],c.dx,c.u[ivar],sort=c.u[isort],log=True,vmin=vmin)

if clump:
    h=ram.rd_clump(nout)
    r = np.sqrt((h.x-center[0])**2+(h.y-center[1])**2+(h.z-center[2])**2)
    nn = np.count_nonzero(r < radius)
    xx = h.x[r < radius]
    yy = h.y[r < radius]
    zz = h.z[r < radius]
    mm = h.m[r < radius]
    plt.plot(xx,yy,'r.')

if args.out:
    plt.savefig(args.out)

plt.show()

