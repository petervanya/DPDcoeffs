#!/usr/bin/env python
"""Usage:
    gen_polymer.py <N> <f> [--L <L> --rho <rho> --save <fname> --xyz <xyz>]

Generate a binary mixture from A and B monomers: AA...ABB...B

Arguments:
    <N>              Polymerisation (number of monomers in chain)
    <f>              fraction of A monomers in chain

Options:
    --L <L>          Size of the simulation box L**3 [default: 10]
    --rho <rho>      Density of the system [default: 5.0]
    --save <fname>   Output save file [default: diblock.data]
    --xyz <xyz>      Create xyz file

pv278@cam.ac.uk, 15/03/16
"""
import numpy as np
from math import *
from docopt import docopt
import lmp_lib as ll


def grow_polymer(L, f, n, Nc, mu=1.0, sigma=0.1):
    """Generate coordinates of matrix polymer chains (taken from gen_pmma.py)
    return (n*Nc, 5) matrix, columns: [mol_ids, bead_type, xyz]
    Input:
    * L: cubic box size
    * n: polymerisation
    * Nc: number of chains
    * mu: mean of bead distance
    * sigma: deviation of bead distance"""
    xyz = np.zeros((Nc*n, 5))
    atom_ids = np.matrix( [1]*int(n*f) + [2]*int(n*(1-f)) ).T
    for i in range(Nc):
        xyz[i*n : (i+1)*n] = np.hstack(( (i+1)*np.matrix(np.ones(n)).T,\
                                         atom_ids,\
                                         grow_one_chain(L, n, Nc, mu, sigma) ))
    return xyz


def grow_one_chain(L, n, Nc, mu, sigma):
    """Return (n, 3) xyz matrix of one chain (taken from gen_pmma.py)"""
    xyz = np.zeros((n, 3))
    xyz[0] = np.random.rand(3)*L
    for i in range(1, n):
        theta = np.random.rand()*pi
        phi = np.random.rand()*2*pi
        r = mu + np.random.randn()*L*sigma
        new_bead_pos = [r*cos(theta), r*sin(theta)*cos(phi), r*sin(theta)*sin(phi)]
        xyz[i] = xyz[i-1] + new_bead_pos
#        xyz[i] = np.where(xyz[i] > L, L, xyz[i])       # set coord to L or 0 on the boundary
#        xyz[i] = np.where(xyz[i] < 0.0, 0.0, xyz[i])
    return xyz


def gen_bonds(n, Nc):
    """Create bond matrix (taken from gen_pmma.py)
    return (n*Nc, 3) matrix, columns: [bond_type, atom1, atom2]
    Input:
    * n: polymerisation
    * Nc: number of chains"""
    mat = np.zeros(((n-1)*Nc, 3), dtype=int)
    for i in range(Nc):
        one_chain_bonds = np.hstack(( np.matrix([1]*(n-1)).T,\
                                      np.matrix( np.arange(n*i+1, n*(i+1)) ).T,\
                                      np.matrix( np.arange(n*i+2, n*(i+1)+1) ).T ))
        mat[i*(n-1) : (i+1)*(n-1)] = one_chain_bonds
    return mat


if __name__ == "__main__":
    args = docopt(__doc__)
    N = int(args["<N>"])
    f = float(args["<f>"])
    L = float(args["--L"])
    rho = float(args["--rho"])
    np.random.seed(1234)

    Nb = int(rho * L**3)
    Nc = int(Nb/N)
    rc = 1.0
    mu, sigma = rc/2.0, rc/10.0
    
    print("=== Creating LAMMPS data file for diblock copolymer melt ===")
    print("Set interaction params in the input file")
    print("Box size:", L, "| Density:", rho, "| Polymerisation:", N, "| A beads fraction:", f)
    
    poly_xyz = grow_polymer(L, f, N, Nc, mu, sigma)  # SET mu CAREFULLY
    xyz_str = ll.atoms2str(poly_xyz)
    print(len(poly_xyz), "beads created, density:", len(poly_xyz)/L**3)

    bonds = gen_bonds(N, Nc)
    bonds_str = ll.bonds2str(bonds)
    print(len(bonds), "bonds created")

    final_string = ll.header2str(len(poly_xyz), len(bonds), 2, 1, L) + \
                   ll.mass2str({1: 1.0, 2: 1.0}) + \
                   "\nAtoms\n\n" + xyz_str +  \
                   "\nBonds\n\n" + bonds_str
    
    fname = args["--save"]
    open(fname, "w").write(final_string)
    print("Data file written in", fname)

    if args["--xyz"]:
        fname = args["--xyz"]
        ll.save_xyzfile(fname, poly_xyz[:, 1:])
        print("xyz file saved in", fname)



