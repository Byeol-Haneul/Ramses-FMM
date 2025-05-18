#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from astropy.io import ascii
import argparse
import os
parser = argparse.ArgumentParser()
parser.add_argument("file", help="enter filename run.log")
parser.add_argument("--log", help="plot log variable",action="store_true")
parser.add_argument("--sym", help="use a circle for each cell",action="store_true")
parser.add_argument("--out", help="output a png image")
args = parser.parse_args()
print("Reading "+args.file)

# path the the file
path_to_output = args.file

# read log file
cmd="grep Main "+path_to_output+" | cut -b 95-109 > /tmp/emag.txt"
os.system(cmd)
cmd="grep Main -A2 "+path_to_output+" | grep dt | cut -b 19-33 > /tmp/time.txt"
os.system(cmd)
cmd="grep Main -A2 "+path_to_output+" | grep dt | cut -b 34-47 > /tmp/dt.txt"
os.system(cmd)
cmd="grep Main -A2 "+path_to_output+" | grep dt | cut -b 48-60 > /tmp/aexp.txt"
os.system(cmd)
do_mag=True
try:
    lines = ascii.read("/tmp/emag.txt")
    emag = lines["col2"]
except Exception as e:
    print("no magnetic energy")
    do_mag=False
do_time=True
try:
    lines = ascii.read("/tmp/time.txt")
    time = lines["col2"]
except Exception as e:
    print("no proper time")
    do_time=False
lines = ascii.read("/tmp/dt.txt")
dt = lines["col2"]
lines = ascii.read("/tmp/aexp.txt")
aexp = lines["col2"]

if do_mag:
    if args.sym:
        plt.plot(time,emag,"o")
    else:
        plt.plot(time,emag)
    if args.log:
        plt.yscale("log")
    plt.xlabel('time')
    plt.ylabel('magnetic energy')

    if args.out:
        plt.savefig(args.out)
    plt.show()

if do_time:
    if args.sym:
        plt.plot(time,dt,"o")
    else:
        plt.plot(time,dt)
else:
    if args.sym:
        plt.plot(aexp,dt,"o")
    else:
        plt.plot(aexp,dt)
    
if args.log:
    plt.yscale("log")
plt.xlabel('time')
plt.ylabel('time step')

if args.out:
    plt.savefig(args.out)
plt.show()


