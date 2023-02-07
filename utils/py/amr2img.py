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
args = parser.parse_args()
# path the the file
path = args.path
if path==None:
    path="./"
else:
    path=path+"/"

nout = args.nout
print("Reading output number ",nout)
print(path)


c=ram.rd_cell(nout,path=path)

#vmin = kwargs.get("vmin",None)
#vmax = kwargs.get("vmax",None)
#sort = kwargs.get("sort",None)

x=c.x[0]
y=c.x[1]
v=c.u[0]
dx=c.dx

xmin=np.min(x-dx/2)
xmax=np.max(x+dx/2)
ymin=np.min(y-dx/2)
ymax=np.max(y+dx/2)

if args.log:
    v = np.log10(abs(v))

vmin=np.min(v)
vmax=np.max(v)

print("min=",vmin," max=",vmax)

plt.rcParams['figure.dpi'] = 58
plt.rcParams.update({'font.size': 22})

px = 1/plt.rcParams['figure.dpi']

fig, ax = plt.subplots(figsize=(1000*px,1000*px))
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
plt.scatter(x,y,s=0.0001)
ax.set_aspect("equal")
rescale=np.maximum(xmax-xmin,ymax-ymin)
plt.scatter(x,y,c=v,s=(dx*800/rescale)**2,cmap="viridis",marker="s",vmin=vmin,vmax=vmax)

if args.log:
    plt.colorbar(shrink=0.8,label="log density")
else:
    plt.colorbar(shrink=0.8,label="density")

plt.xlabel("x")
plt.ylabel("y")

if args.out:
    plt.savefig(args.out)

plt.show()

