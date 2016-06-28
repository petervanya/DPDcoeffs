#!/usr/bin/env python2
"""Usage:
    rdf2.py <files> [--rc <rc> --boxsize <L> --types <n> --bins <nbins>] 

Read xyz files and compute radial distribution function for a given atom type.
Specify cutoff rc to consider only pairs up to this distance.
Correct implementation of PBCs.
f2py does not work with Python3!

Arguments:
    <files>             Regex match xyz files

Options:
    --types <n>         Atom name in number [default: 1]
    --rc <rc>           Only pairs up to this <rc> in in DPD units
    --bins <nbins>      Number of bins [default: 1000]
    --boxsize <L>       Box size in arb units [default: 20]

pv278@cam.ac.uk, 28/06/16
"""
import numpy as np
import glob, sys
from docopt import docopt
from f_rdf import f_rdf      # Fortran module
import lmp_lib as ll


def rdf_1type(dumpfiles, atom_type, rc, cell, bins, bc=1.3):
    """construct an rdf for water beads from given xyz frames.
    Compute radial dist'n fcn from the xyz frame using 
    Fortran routine for pair distances and numpy binning
    """
    # FIX IF rc = L
    rdf = np.zeros(len(bins)-1)
    Nd = len(dumpfiles)
    for dump in dumpfiles:
        A = ll.read_xyzfile(dump)
        xyz = A[A[:, 0] == atom_type][:, 1:]
        N = len(xyz)
        print("Calculating rdf...")

        r = bins[:-1] + np.diff(bins)/2.0
        dr = bins[1] - bins[0]
        L = cell[0, 0]

        nn = int(2*N**2 * (rc/L)**3 * bc)  # Bulgarian const
        d = f_rdf.dist_vec_cut(xyz, rc, L, cell, nn)
        d = d[d != 0.0]
        rdf_raw, _ = np.histogram(d, bins)
        rdf += rdf_raw / (4*np.pi*r**2 * dr) * L**3 / (N*(N-1)/2)
        print("Done: %s" % dump)
    rdf /= Nd
    return rdf


def rdf_2types(dumpfiles, atom_types, rc, cell, bins, bc=1.3):
    """construct an rdf for water beads from given xyz frames"""
    # FIX IF rc = L
    rdf = np.zeros(len(bins)-1)
    Nd = len(dumpfiles)
    for dump in dumpfiles:
        A = ll.read_xyzfile(dump)
        xyz1 = A[A[:, 0] == atom_types[0]][:, 1:]
        xyz2 = A[A[:, 0] == atom_types[1]][:, 1:]
        N1 = len(xyz1)
        N2 = len(xyz2)
        print("Atoms: %i %i | Calculating rdf..." % (N1, N2))

        r = bins[:-1] + np.diff(bins)/2.0
        dr = bins[1] - bins[0]
        L = cell[0, 0]

        nn = int(2*N1**2 * (rc/L)**3 * bc)  # Bulgarian const
        d = f_rdf.dist_vec_cut_2mat(xyz1, xyz2, rc, L, cell, nn)
        d = d[d != 0.0]
        rdf_raw, _ = np.histogram(d, bins)
        rdf += rdf_raw /(4*np.pi*r**2 * dr) * L**3 / (N1*N2)
        print("Done: %s" % dump)
    rdf /= Nd
    return rdf


if __name__ == "__main__":
    args = docopt(__doc__)
#    print(args)
    dumpfiles = glob.glob(args["<files>"])
    L = float(args["--boxsize"])
    Nbins = int(args["--bins"])
    Nfiles = len(dumpfiles)
    N = int(open(dumpfiles[0], "r").readline())
 
    if args["--rc"]:
        rc = float(args["--rc"])
        bins = np.linspace(0, rc, Nbins+1)
        r = bins[:-1] + np.diff(bins)/2.0
    else:
        rc = L
    bins = np.linspace(0, L, Nbins+1)
    r = bins[:-1] + np.diff(bins)/2.0
    
    atom_types = args["--types"]
    atom_types = [int(i) for i in atom_types.split()]
    if len(atom_types) > 2:
        print("Only two atom types allowed for rdf calculation.")
        sys.exit()
    N_xyz_types = set(ll.read_xyzfile(dumpfiles[0])[:, 0])
    if not set(atom_types).issubset(N_xyz_types):
        print("Requested atom types not present in the xyz files.")
        sys.exit()

    cell = L*np.eye(3)
    if len(dumpfiles) == 0:
        raise ValueError("No xyz files captured, aborting.")
    
    print("===== Calculating rdf =====")
    print("Atom types: %s | Bins: %i | Cutoff: %.2f | xyz files: %i" % \
         (repr(atom_types), Nbins, rc, len(dumpfiles)))

    if len(atom_types) == 1:
        vals = rdf_1type(dumpfiles, atom_types[0], rc, cell, bins, bc=1.3)
    else:
        vals = rdf_2types(dumpfiles, atom_types, rc, cell, bins, bc=1.3)

    fname = "rdf_water.out"
    np.savetxt(fname, np.vstack((r, vals)).T)
    print("rdf saved in", fname)


