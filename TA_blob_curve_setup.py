#!/usr/bin/env python
"""Usage:
    TA_blob_curve_setup.py <position> <dmin> <dmax> <N> [--blobnum <bn>]

[AD HOC] Generate gjf files of one water blob on triflic acid

Arguments:
    <position>   Where on TA is the blob placed TA (number 1..6, see README)
    <dmin>       Minimum distance of two blobs
    <dmax>       Maximum distance of two blobs
    <N>          Number of points in between

Options:
    --blobnum <bn>   Which blob to use [default: 0]

pv278@cam.ac.uk, 05/10/15
"""
import numpy as np
import os, sys
from math import pi
from docopt import docopt
from xyzlib import Atoms


def save_configs(TA, blob, Drange, theta, phi, outdir, blobnum="0"):
    for d in Drange:
        xyzname = "run_blob_" + blobnum + "_d" + str(d) + ".xyz"
        xyzpath = os.path.join(outdir, xyzname)
        dist = np.array([[d, 0, 0]])
        dist = Atoms(["A"], dist)
        dist.rotate(theta, phi)
        dist = dist.coords[0].round(2)
        blob.shift(dist)
        TA_blob = triflic + blob
        blob.shift(-dist)   # STUPID TRICK
        TA_blob.save(xyzpath)


args = docopt(__doc__)
#print args
try:
    pos = args["<position>"]
except int(pos) not in range(1, 7):
    print "Position numbers are 1...6"
print "pos:", pos, type(pos)
maindir = os.path.expanduser("~/DPDcoeffs/TA_WaterBlob")
blobdir = os.path.expanduser("~/DPDcoeffs/Files/Waterblobs")
blobnum = args["--blobnum"]   # blob number used, default is 0
blobfile = os.path.join(blobdir, "waterblob_" + blobnum + ".xyz")
blob = Atoms().read(blobfile)
blob.shift_com()

triflic = Atoms().read(os.path.expanduser("~/DPDcoeffs/Files/triflic.xyz"))
Ccoord = triflic.coords[[i for i, x in enumerate(triflic.names) if x == "C"][0]]

Dmin = float(args["<dmin>"])
Dmax = float(args["<dmax>"])
N = int(args["<N>"])
Drange = np.linspace(Dmin, Dmax, N).round(2)
print "Range of distances: ", Drange

outdir = os.path.join(maindir, "Pos_" + pos)
if not os.path.exists(outdir):
    os.mkdir(outdir)

if pos == "1":       # above C
    Drange += Ccoord[2]  # shift by z-coord
    Drange = Drange.round(2)
    theta, phi = pi/2.0, 0.0
    save_configs(triflic, blob, Drange, theta, phi, outdir, blobnum)
elif pos == "2":     # below S
    theta, phi = -pi/2.0, 0.0
    save_configs(triflic, blob, Drange, theta, phi, outdir, blobnum)
elif pos == "3":     # towards H
    theta, phi = 0.0, 0.0
    save_configs(triflic, blob, Drange, theta, phi, outdir, blobnum)
elif pos == "4":     # opposite H
    theta, phi = 0.0, pi
    save_configs(triflic, blob, Drange, theta, phi, outdir, blobnum)
elif pos == "5":     # towards O w/o H
    theta, phi = 0.0, pi/3.0
    save_configs(triflic, blob, Drange, theta, phi, outdir, blobnum)
elif pos == "6":     # between O w/ H and O w/o H
    theta, phi = 0.0, 2*pi/3.0
    save_configs(triflic, blob, Drange, theta, phi, outdir, blobnum)

     

