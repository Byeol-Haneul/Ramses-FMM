#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import argparse
import miniramses as ram

parser = argparse.ArgumentParser()
parser.add_argument("nout", help="enter output number")
parser.add_argument("--path", help="specify a path")
parser.add_argument("--log", help="plot log SFR",action="store_true")
parser.add_argument("--out", help="output a png image")
parser.add_argument("--xcen", help="specify the image center x-coordinate")
parser.add_argument("--ycen", help="specify the image center y-coordinate")
parser.add_argument("--zcen", help="specify the image center z-coordinate")
parser.add_argument("--rad", help="specify the image radius")
args = parser.parse_args()
# path the the file
path = args.path
radius = args.rad
xcenter = args.xcen
ycenter = args.ycen
zcenter = args.zcen
log = args.log

if path==None:
    path="./"
else:
    path=path+"/"
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

nout = args.nout
print("Reading output number ",nout)

s=ram.rd_part(nout,path=path,prefix='star',center=center,radius=radius)
i=ram.rd_info(nout,path=path)
bins=np.linspace(0,1,100)
unit_m=i.unit_d*i.unit_l**3/2e33/(bins[1]-bins[0])/1e9
plt.hist(s.tp*i.unit_t/3e7/1e9,weights=s.mp*unit_m,bins=bins)
if log:
    plt.yscale("log")
plt.xlabel('t [Gyr]')
plt.ylabel('SFR [Msol/yr]')

if args.out:
    plt.savefig(args.out)

plt.show()

