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
cmd="grep -n Output "+path_to_output+" > /tmp/out.txt"
os.system(cmd)
lines = ascii.read("/tmp/out.txt")
i = int(lines["col1"][-1][:-1])
n = int(lines["col3"][-1])

data = ascii.read(path_to_output,header_start=i-3,data_start=i-2,data_end=i+n-2)

if args.sym:
    plt.plot(data["x"],data["d"],"o")
else:
    plt.plot(data["x"],data["d"])

if args.out:
    plt.savefig(args.out)
plt.show()


