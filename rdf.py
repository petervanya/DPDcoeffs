#!/usr/bin/env python
"""Usage:
    rdf.py <fnames> [--binsize <bs>] [--picname <pn>]

Read all LAMMPS data files in the directory and compute the pair correlation function.
Produce a histogram plot and save data into rdf.out.
Uses Fortran routines from mat_ops.f95 produced with f2py.

Arguments:
    <fnames>          Core of the filename, extension is 'xyz'

Options:
    --binsize <bs>   Size of the bins for the histogram [default: 0.05]
    --picname <pn>   Name of the histogram plot [default: hist.png]

pv278@cam.ac.uk, 02/09/15
"""
from docopt import docopt
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from math import pi
import os, glob, time
import mat_ops

def read_outfile(outfile):
    """Read one xyz outfile into a numpy matrix"""
    A = open(outfile, "r").readlines()[2:]
    A = [line.split() for line in A]
    A = np.array(A, order="F").astype(float)
    return A

def get_one_rdf(mat, dr=0.05):
    """Compute pair correl'n fcn from the xyz matrix"""
    N = len(mat)
    Npairs = N*(N-1)/2
#    d = np.zeros(Npairs)
#    here = 0
#    for i in range(N):
#        for j in range(i+1, N):
#            d[here] = norm(mat[i] - mat[j])
#            here += 1
    d = mat_ops.get_pair_dist(mat)

    bins = np.arange(0, 10+dr, dr)         # assume box size is 10
    (rdf_raw, r, c) = plt.hist(d, bins=bins)  # c is useless
    r = r[:-1] + np.diff(r)/2
    rdf = rdf_raw/(4*pi*r**2*dr)
    return r, rdf


def get_master_rdf(outfiles, N):
    """Construct an rdf from all the available xyz files"""
    Nfiles = len(outfiles)
    Npairs = N*(N-1)/2
#    rdf_mat = np.zeros((Npairs, Nfiles))
    rdf_mat = []
    for i in range(Nfiles):
        A = read_outfile(outfiles[i])
        r, rdf_i = get_one_rdf(A[:, 1:])
        rdf_mat.append(rdf_i)
        print outfiles[i], "done"
    rdf = np.sum(np.array(rdf_mat), 0)/len(outfiles)
    return r, rdf


def save_data(vec1, vec2, outfile):
    with open(outfile, "w") as f:
        for i in range(len(vec1)):
            f.write(str(vec1[i]) + "\t" + str(vec2[i]) + "\n")


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args
    dr = float(args["--binsize"])
    outfiles = glob.glob(args["<fnames>"])
    print outfiles
    Nfiles = len(outfiles)
    N = int(open(outfiles[0], "r").readline())
    print "Particles:", N
    
    r, vals = get_master_rdf(outfiles, N)
#    A = read_outfile(outfiles[0])
#    r, vals = get_one_rdf(A[:, 1:], dr=dr)

    plt.cla()
    plt.plot(r, vals)
    plt.savefig(args["--picname"])
    save_data(r, vals, "rdf.out")



