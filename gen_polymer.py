#!/usr/bin/env python
"""Usage:
    gen_polymer.py <N> <f> [--dist <a>] [--FAA <AA>] [--FAB <AB>] 
                           [--seed <s>] [--boxsize <L>] [--save <fname>] 

Generate a binary mixture from A and B monomers: AA...ABB...B
UNFINISHED

Arguments:
    <N>              Polymerisation (# of monomers)
    <f>              fraction of A monomers

Options:
    --FAA <AA>       Force between AA (and BB) beads [default: 1.0]
    --FAB <AB>       Force between AB beads [default: 1.0]
    --dist <a>       Initial distance between successive beads
    --save <fname>   Save into file [default: poly.txt]
    --seed <s>       Random seed [default: 123]
    --boxsize <L>    Size of the simulation box L**3 [default: 10]

pv278@cam.ac.uk, 28/09/15
"""
import numpy as np
from docopt import docopt
from xyzlib import Atoms


def get_copoly_xyz(N, f, a, start=np.zeros(3), angles=[0, 0], sigma=0.1):
    """generate the xyz matrix of the polymer with coordinates 
    wiggled by gaussian noise.
    Atom_style molecular
    * sigma -- standard deviation for the noise"""
    NA = int(f*N)
    types = [1]*NA + [2]*(N-NA)
    coords = np.hstack(( np.zeros((N, 2)), np.matrix(np.linspace(0, (N-1)*a, N)).T ))  # z-coord
    xyz = Atoms(names=[], coords=coords)
    xyz.rotate(angles[0], angles[1])           # rotate the polymer
    xyz.coords += np.random.randn(N, 3)*0.1    # gaussian noise
    xyz.coords += start                        # shift
    xyz_mat = np.hstack((np.arange(1, N+1, dtype=int).reshape((N, 1)), \
                         np.ones((N, 1), dtype=int), np.ones((N, 1), dtype=int), \
                         xyz.coords))
    return xyz_mat


def get_copoly_bonds(N):
    """Generate the bond matrix: [num, type, part1, part2]
    Stiffness constants (type) are the same"""
    bond_mat = np.vstack((np.arange(1, N), np.ones(N-1, dtype=int), \
                          np.arange(1, N), np.arange(2, N+1) )).T
    return bond_mat


def get_header(xyz_mat, bond_mat, ang_mat=np.array([]), \
               dih_mat=np.array([]), imp_mat=np.array([])):
    """Generate LAMMPS header"""
    s = "#blabla\n"
    s += str(len(xyz_mat)) + " atoms\n"
    s += str(len(bond_mat)) + " bonds\n"
#    s += str(len(ang_mat)) + " angles\n"
    s += "\n"
    s += str(len(set(xyz_mat[:, 1]))) + " atom types\n"   # select unique
    s += str(len(set(bond_mat[:, 1]))) + " bond types\n"
    s += "\n"
    s += "0.0 " + str(L) + " xlo xhi\n"
    s += "0.0 " + str(L) + " ylo yhi\n"
    s += "0.0 " + str(L) + " zlo zhi\n"
    return s


def get_masses(id_mass):
    """Save list of masses in a string"""
    s = "Masses\n\n"
    for k, v in id_mass.iteritems():
        s += str(k) + "\t" + str(v) + "\n"
    return s


def get_bond_coeffs(id_bondcoeffs):
    """Save bond coefficients into a string"""
    s = "Bond Coeffs\n\n"
    for k, v in id_bondcoeffs.iteritems():
        str_v = ""
        for i in len(v):
            str_v += str(i) + "\t"
        s += str(k) + "\t" + str_v + "\n"
    return s


def mat_to_str(mat):
    """Convert matrix to string"""
    s = ""
    M, N = mat.shape
    for i in range(M):
        for j in range(N):
            s += str(mat[i, j]) + "\t"
        s += "\n"
    return s


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args
    N = int(args["<N>"])
    f = int(args["<f>"])
    a = float(args["--dist"])
    L = float(args["--boxsize"])
    FAA = float(args["--FAA"])
    FAB = float(args["--FAB"])
    fname = args["--save"]
    seed = int(args["--seed"])
    np.random.seed(seed)
    
    xyz_mat = get_copoly_xyz(N, f, a) 
    bond_mat = get_copoly_bonds(N)

    final_string = get_header(xyz_mat, bond_mat) + \
                   "\nAtoms\n\n" + mat_to_str(xyz_mat) + \
                   "\nBonds\n\n" + mat_to_str(bond_mat)
    
    print final_string
