#!/usr/bin/env python
"""Usage:
    rdf_general.py <files> [--rc <rc> --boxsize <L> --types <n> --bins <nbins>] 

Read xyz files and compute radial distribution function for a given atom type.
Specify cutoff rc to consider only pairs up to this distance.
Correct implementation of PBCs.

Arguments:
    <files>             Regex match xyz files

Options:
    --types <n>         Atom name in number [default: 1]
    --rc <rc>           Only pairs up to this <rc> in in DPD units
    --bins <nbins>      Number of bins [default: 500]
    --boxsize <L>       Box size in arb units [default: 20]

pv278@cam.ac.uk, 28/06/16
"""
import numpy as np
from numpy import pi
import glob, sys
from docopt import docopt
from f_rdf import f_rdf      # Fortran module
import lmp_lib as ll


class mydict(dict):
    """A container for all the system constants"""
    def __getattr__(self, key):
        return self[key]


def rdf_1type(dumpfiles, sp):
    """Construct an rdf for water beads from given xyz frames.
    Compute radial dist'n fcn from the xyz frame using 
    Fortran routine for pair distances and numpy binning
    for one particle type.
    Input:
    * dumpflies: list of xyz frame names
    * sp: system params, containing:
        * atom_type: integer
        * rc: cutoff for atom neighbourhood
        * cell: (3, 3) matrix
        * bins: vector
        * bc: Bulgarian const to guess the distance vector size
    """
    rdf = np.zeros(len(sp.bins)-1)
    L = sp.cell[0, 0]
    r = sp.bins[:-1] + np.diff(sp.bins)/2.0
    dr = r[1] - r[0]
    norm = sp.N*(sp.N-1)/2

    if sp.rc == L/2:
        nn = norm
    else:
        nn = int(2*sp.N**2 * (sp.rc/L)**3 * sp.bc)  # estimate at dist vec length

    for dump in dumpfiles:
        A = ll.read_xyzfile(dump)
        xyz = A[A[:, 0] == sp.atom_types[0]][:, 1:]
        print("Calculating rdf...")
        d = f_rdf.dist_vec_cut(xyz, sp.rc, L, sp.cell, nn)
        d = d[d != 0.0]
        rdf_raw, _ = np.histogram(d, sp.bins)
        rdf += rdf_raw / (4*pi*r**2 * dr) * L**3 / norm
        print("Done: %s" % dump)
    return rdf / sp.Nd     # normalise


def rdf_2types(dumpfiles, sp):
    """Construct an rdf for water beads from given xyz frames
    Compute radial dist'n fcn from the xyz frame using 
    Fortran routine for pair distances and numpy binning
    for one particle type.
    Similar to rdf_1type."""
    rdf = np.zeros(len(bins)-1)
    L = cell[0, 0]
    r = sp.bins[:-1] + np.diff(sp.bins)/2.0
    dr = r[1] - r[0]
    A = ll.read_xyzfile(dumpfiles[0])
    N1 = len(A[A[:, 0] == sp.atom_types[0]][:, 1:])
    N2 = len(A[A[:, 0] == sp.atom_types[1]][:, 1:])
    norm = N1*N2

    if rc == L/2:
        nn = norm
    else:
        nn = int(2*max(N1, N2)**2 * (sp.rc/L)**3 * sp.bc)

    for dump in dumpfiles:
        A = ll.read_xyzfile(dump)
        xyz1 = A[A[:, 0] == sp.atom_types[0]][:, 1:]
        xyz2 = A[A[:, 0] == sp.atom_types[1]][:, 1:]
        print("Atoms: %i %i | Calculating rdf..." % (N1, N2))

        d = f_rdf.dist_vec_cut_2mat(xyz1, xyz2, sp.rc, L, sp.cell, nn)
        d = d[d != 0.0]
        rdf_raw, _ = np.histogram(d, sp.bins)
        rdf += rdf_raw /(4*pi*r**2 * dr) * L**3 / norm/2
        print("Done: %s" % dump)
    return rdf / sp.Nd


if __name__ == "__main__":
    args = docopt(__doc__)
#    print(args)
    dumpfiles = glob.glob(args["<files>"])
    L = float(args["--boxsize"])
    Nbins = int(args["--bins"])
    Nd = len(dumpfiles)
    N = int(open(dumpfiles[0], "r").readline())
 
    if args["--rc"]:
        rc = float(args["--rc"])
        bins = np.linspace(0, rc, Nbins+1)
        r = bins[:-1] + np.diff(bins)/2.0
    else:
        rc = L/2        # max distance between atoms
    bins = np.linspace(0, L/2, Nbins+1)
    r = bins[:-1] + np.diff(bins)/2.0
    
    atom_types = [int(i) for i in args["--types"].split()]
    if len(atom_types) > 2:
        print("Only two atom types allowed for rdf calculation.")
        sys.exit()
    xyz_types = set(ll.read_xyzfile(dumpfiles[0])[:, 0])
    if not set(atom_types).issubset(xyz_types):
        print("Requested atom types not present in the xyz files.")
        sys.exit()

    cell = L*np.eye(3)
    if len(dumpfiles) == 0:
        raise ValueError("No xyz files captured, aborting.")

    sp = mydict(N=N, Nd=Nd, rc=rc, cell=cell, bins=bins,\
                bc=1.3, atom_types=atom_types)
    
    print("===== Calculating rdf =====")
    print("Atom types: %s | Bins: %i | Cutoff: %.2f | xyz files: %i" % \
         (repr(atom_types), Nbins, rc, Nd))

    if len(atom_types) == 1:
        vals = rdf_1type(dumpfiles, sp)
        fname = "rdf_%i.out" % atom_types[0]
    else:
        vals = rdf_2types(dumpfiles, sp)
        fname = "rdf_%i_%i.out" % (atom_types[0], atom_types[1])

    np.savetxt(fname, np.vstack((r, vals)).T)
    print("rdf saved in", fname)


