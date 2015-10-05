#!/usr/bin/env python
"""Usage:
   extract_energies.py <dir> [--print | --save <fname>]

Extract SCF energy from multiple Gaussian output files in a given dir

Arguments:
   <dir>            Directory to look for Gaussian files

Options:
   --save <fname>   Save to file
   --print          Print table to screen

pv278@cam.ac.uk, 31/07/15
"""
import numpy as np
import pandas as pd
import sys, os, re
from docopt import docopt

args = docopt(__doc__, version=0.0)
#print args
outdir = args["<dir>"]
outfiles = [file for file in os.listdir(outdir) if file.split(".")[-1] == "out"]
positions, energies = [], []

for file in outfiles:
    filepath = outdir + "/" + file
    if not os.path.exists(filepath):
        continue
    E_line = [line for line in open(filepath) if "SCF Done" in line]
    if E_line:       # if energy line exists
        E_line = E_line[0]
        energies.append(float(E_line.split()[4]))
        pos = re.search(r"(?<=_d)[^}]*(?=.out)", file)   # WHAT IS THIS DOING?
        pos = pos.group()
        positions.append(float(pos))

if not energies:     # if the array is empty
    sys.exit()

tbl = pd.DataFrame(energies, index=positions).sort_index()
tbl.columns = ["E"]
print tbl

if args["--save"]:
    datapath = args["--save"]
    tbl.to_csv(datapath, sep="\t", header=False)
    print "Table saved in ", datapath
