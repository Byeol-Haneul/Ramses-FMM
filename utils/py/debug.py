#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import argparse
import miniramses as ram

parser = argparse.ArgumentParser()
parser.add_argument("nout", help="enter output number")
parser.add_argument("--path", help="specify a path")
parser.add_argument("--bkp", help="plot log variable",action="store_true")
parser.add_argument("--prefix", help="specify a file prefix")
parser.add_argument("--xcen", help="specify the image center x-coordinate")
parser.add_argument("--ycen", help="specify the image center y-coordinate")
parser.add_argument("--zcen", help="specify the image center z-coordinate")
parser.add_argument("--rad", help="specify the image radius")
args = parser.parse_args()
# path the the file
path = args.path
prefix = args.prefix
radius = args.rad
xcenter = args.xcen
ycenter = args.ycen
zcenter = args.zcen
backup = args.bkp

if xcenter==None:
    xcenter=None
else:
    xcenter=float(xcenter)
if ycenter==None:
    ycenter=None
else:
    ycenter=float(ycenter)
if zcenter==None:
    zcenter=None
else:
    zcenter=float(zcenter)
if radius==None:
    radius=None
else:
    radius=float(radius)
center=np.array([xcenter,ycenter,zcenter])
if path==None:
    path="./"
else:
    path=path+"/"

nout = args.nout
print("Reading output number ",nout)

i=ram.rd_info(nout,backup=backup)
c=ram.rd_cell(nout,backup=backup)

d=c.u[0]
if(backup):
    ux=c.u[1]/c.u[0]
    uy=c.u[2]/c.u[0]
    uz=c.u[3]/c.u[0]
else:
    ux=c.u[1]
    uy=c.u[2]
    uz=c.u[3]
if(backup):
    ekin=0.5*d*(ux*ux+uy*uy+uz*uz)
    p=2/3*(c.u[4]-ekin)
else:
    p=c.u[4]
cs=np.sqrt(5/3*p/d)
dt=c.dx/(3*cs+np.abs(ux)+np.abs(uy)+np.abs(uz))

print("min d=",np.min(d)," max d=",np.max(d))
print("max ux=",np.max(np.abs(ux)))
print("max uy=",np.max(np.abs(uy)))
print("max uz=",np.max(np.abs(uz)))
print("max cs=",np.max(np.abs(cs)))
print("min dt=",np.min(np.abs(dt)))




