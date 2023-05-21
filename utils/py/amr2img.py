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
ram.visu(c.x[0],c.x[1],c.dx,c.u[0],sort=c.u[0],log=True)

if args.out:
    plt.savefig(args.out)

plt.show()

