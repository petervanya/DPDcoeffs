#!/usr/bin/env python
"""Usage:
    rdf.py gen_data <fnames> [--beads <b>] [--binsize <bs>] [--picname <pn>]
    rdf.py plot <fnames>

Read all LAMMPS data files in the directory and compute the pair correlation function.
Uses Fortran routines from mat_ops.f95 produced with f2py.
* gen_data  generate the histogram of distances and save it in rdf.out
* plot      plot rdf.out

Arguments:
    <fnames>         Regex for all the required xyz files

Options:
    --binsize <bs>   Size of the bins for the histogram [default: 0.05]
    --picname <pn>   Name of the histogram plot [default: hist]
    --beads <b>      Bead type, 'all' or 1,...,n [default: all]

pv278@cam.ac.uk, 02/09/15
"""
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from math import pi
import os, glob
from docopt import docopt
import mat_ops           # Fortran .so file

def read_outfile(outfile):
    """Read one xyz outfile into a numpy matrix"""
    A = open(outfile, "r").readlines()[2:]
    A = [line.split() for line in A]
    A = np.array(A, order="F").astype(float)
    return A


def save_data(outfile, *args): #vec1, vec2):
    """Save two vectors into file"""
    m, n = len(args), len(args[0])   # cols, rows
    args = zip(*args)
    with open(outfile, "w") as f:
        for i in range(n):
            line = ""
            for j in range(m):
                line += str(args[i][j]) + "\t"
            line += "\n"
            f.write(line)


def plot_hist(r, vals, picname="hist.png", title="title"):
    """Plot a rdf function from given values"""
    plt.clf()
    plt.plot(r, vals)
    ymax = np.max(vals[1:])*1.2/200*200     # round to 200
    plt.ylim([0, ymax])
    plt.xlabel("$r$ (DPD units)")
    plt.ylabel("$g(r)$")
    plt.title(title)
    plt.savefig(picname)
    print "rdf saved into", picname


def get_one_rdf_old(outfile, dr=0.05, bead_type="all"):
    """Compute radial dist'n fcn from the xyz matrix using Python"""
    A = read_outfile(outfile)
    if bead_type == "all":
        xyz_mat = A[:, 1:]
    else:
        xyz_mat = A[A[:, 0] == int(bead_type)][:, 1:]
    N = len(xyz_mat)
    Npairs = N*(N-1)/2
    d = np.zeros(Npairs)
    here = 0
    for i in range(N):
        for j in range(i+1, N):
            d[here] = norm(xyz_mat[i] - xyz_mat[j])
            here += 1
    bins = np.arange(0, 10+dr, dr)            # assume box size is 10
    (rdf_raw, r, c) = plt.hist(d, bins=bins)  # c is useless here
    r = r[:-1] + np.diff(r)/2.0
    rdf = rdf_raw/(4*pi*r**2*dr)
    return r, rdf


def get_one_rdf(outfile, dr=0.05, bead_type="all"):
    """Compute radial dist'n fcn from the xyz matrix using Fortran routine"""
    A = read_outfile(outfile)
    if bead_type == "all":
        xyz_mat = A[:, 1:]
    else:
        xyz_mat = A[A[:, 0] == int(bead_type)][:, 1:]
    N = len(xyz_mat)
    d = mat_ops.get_pair_dist(xyz_mat)        # call Fortran

    bins = np.arange(0, 10+dr, dr)            # assume box size is 10
    (rdf_raw, r, c) = plt.hist(d, bins=bins)  # c is useless here
    r = r[:-1] + np.diff(r)/2.0
    rdf = rdf_raw/(4*pi*r**2*dr)
    return r, rdf


def get_master_rdf(outfiles, dr=0.05, bead_type="all"):
    """Construct an rdf from all the available xyz files"""
    Nfiles = len(outfiles)
    rdf_mat = []
    for outfile in outfiles:
        r, rdf_i = get_one_rdf(outfile, dr, bead_type)
        rdf_mat.append(rdf_i)
        print outfile, "done"
    rdf_mat = np.array(rdf_mat).T
    np.savetxt("rdf_mat.out", rdf_mat, fmt="%.0f", delimiter="\t")
    rdf = np.array(np.sum(rdf_mat, 1)/len(outfiles))
    return r, rdf


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args
    if args["gen_data"]:
        dr = float(args["--binsize"])
        bead_type = args["--beads"]
        outfiles = glob.glob(args["<fnames>"])
        print outfiles
        Nfiles = len(outfiles)
        N = int(open(outfiles[0], "r").readline())
        print "Particles:", N
        r, vals = get_master_rdf(outfiles, dr, bead_type)
        save_data("rdf_" + bead_type + ".out", r, vals)

    elif args["plot"]:   # reads rdf.out and creates histograms
        outfiles = glob.glob(args["<fnames>"])
        for outfile in outfiles:
            A = np.loadtxt(outfile)
            r, vals = A[:, 0], A[:, 1]
            outfile = outfile[:-4]
            bead_type = outfile.split("_")[1]    # "all", "1", etc
            dirname = "_".join(outfile.split("_")[2:])
            picname = args["--picname"] + "_" + bead_type + "_" + dirname + ".png"
            title = "beads " + bead_type
            plot_hist(r, vals, picname, title=title)
        

