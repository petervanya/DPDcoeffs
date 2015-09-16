#!/usr/bin/env python
"""Usage:
    rdf.py <dir> <fname> <binsize>

Read all LAMMPS data files and compute the pair correlation function.

Arguments:
    <dir>          Directory with output files
    <fname>        Core of the filename, extension is 'xyz'
    <binsize>      Size of the bins for the histogram [default: 0.05]

pv278@cam.ac.uk, 02/09/15
"""
from docopt import docopt
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from math import pi
import os
import glob

def read_outfile(outfile):
    """Read one xyz outfile into a np matrix"""
    A = open(outfile, "r").readlines()[2:]
    A = [line.split() for line in A]
    A = np.array(A).astype(float)
    return A

def get_one_rdf(mat, dr=0.05):
    """Compute pair correl'n fcn from matrix of positions"""
    N = len(mat)
    Npairs = N*(N-1)/2
    d = []
    for i in range(N):
        for j in range(i+1, N):
            d.append(norm(mat[i] - mat[j]))
    bins = np.arange(0, 10+dr, dr)         # assume box size is 10
    (rdf, r, c) = plt.hist(d, bins=bins)
    r = r[:-1] + np.diff(r)/2
    return r, rdf/(4*pi*r**2*dr)

def get_master_rdf(outfiles):
    Npairs = Nparts*(Nparts-1)/2
    PCFmat = np.zeros((Npairs, Nfiles))
    for i in range(Nfiles):
        A = read_outfile(outfiles[i])
        r, PCFmat[:, i] = get_one_rdf(A[:, 1:])
        print outfiles[i], "done"
    rdf = np.sum(PCFmat, 1)/len(outfiles)
    return r, rdf

def save_data(vec1, vec2, outfile):
    with open(outfile, "w") as f:
        for i in range(len(vec1)):
            f.write(str(vec1[i]) + "\t" + str(vec2[i]) + "\n")


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args
    dr = float(args["<binsize>"])
    outfiles = glob.glob(args["<dir>"] + "/" + args["<fname>"] + "*.xyz")
    Nfiles = len(outfiles)
    Nparts = int(open(outfiles[0], "r").readline())
    print "Particles:", Nparts

    A = read_outfile(outfiles[0])
    r, vals = get_one_rdf(A[:, 1:], dr=dr)
    print outfiles[0]
#    rdf = get_master_rdf(outfiles)

    plt.cla()
    plt.plot(r, vals)
    plt.savefig("hist.png")
    save_data(r, vals, "rdf.out")



