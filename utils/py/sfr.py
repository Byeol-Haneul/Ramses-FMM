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
parser.add_argument("--xcen", help="specify the region center x-coordinate")
parser.add_argument("--ycen", help="specify the region center y-coordinate")
parser.add_argument("--zcen", help="specify the region center z-coordinate")
parser.add_argument("--rad", help="specify the region radius")
parser.add_argument("--bin", help="specify the bin size in Myr")
args = parser.parse_args()
# path the the file
path = args.path
radius = args.rad
xcenter = args.xcen
ycenter = args.ycen
zcenter = args.zcen
log = args.log
dt = args.bin

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
if dt==None:
    bin_size=0.1  #in Gyr
else:
    bin_size=float(dt)/1000

nout = args.nout
print("Reading output number ",nout)

s=ram.rd_part(nout,path=path,prefix='star',center=center,radius=radius)
i=ram.rd_info(nout,path=path)
time=abs(s.tp*i.unit_t/i.aexp**2/(365*24*3600*1e9))
n_bin=int(np.max(time)/bin_size)
bins=np.linspace(0,np.max(time),n_bin)
unit_m=i.unit_d*i.unit_l**3/2e33/(bins[1]-bins[0])/1e9
plt.hist(time,weights=s.mp*unit_m,bins=bins)
if log:
    plt.yscale("log")
plt.xlabel('t [Gyr]')
plt.ylabel('SFR [Msol/yr]')

if args.out:
    plt.savefig(args.out)

plt.show()

