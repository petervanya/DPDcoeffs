#!/usr/bin/env python
"""Usage:
    gen_binmixt2.py [--f <f> --L <L> --rho <rho --save <fname> --xyz <xyz>]

Create a LAMMPS data file for A/B binary mixture with fraction f of A beads.
atom_style atomic

Options:
    --f <f>          Fraction of A beads [default: 0.5]
    --rho <rho>      Density of the system [default: 3.0]
    --L <L>          Box size [default: 10]
    --save <fname>   Save into file [default: binmixt.data]
    --xyz <xyz>      Create xyz file

pv278@cam.ac.uk, 13/06/16
"""
import numpy as np
import lmp_lib as ll
from docopt import docopt


def get_pos(N, f, L):
    """Return position matrix and type vector of particles"""
    xyz = np.random.rand(N, 3)*L
    names = [1]*int(f*N) + [2]*int((1-f)*N)
    return names, xyz


def header2str(N, L):
    """Generate LAMMPS header"""
    s = "#binmixt\n"
    s += str(N) + " atoms\n"
    s += "2 atom types\n"
    s += "\n"
    s += "0.0 " + str(L) + " xlo xhi\n"
    s += "0.0 " + str(L) + " ylo yhi\n"
    s += "0.0 " + str(L) + " zlo zhi\n\n"
    return s


def mass2str(m1, m2):
    s = "Masses\n\n"
    s += "1 " + str(m1) + "\n"
    s += "2 " + str(m2) + "\n\n"
    return s


def atoms2str(names, xyz):
    """Convert atomic matrix to string"""
    M = len(xyz)
    s = "Atoms\n\n"
    for i in range(M):
        s += "%i\t%i\t%f\t%f\t%f\n" % (i+1, names[i], xyz[i, 0], xyz[i, 1], xyz[i, 2])
    s += "\n"
    return s


if __name__ == "__main__":
    args = docopt(__doc__)
    np.random.seed(1234)
    f = float(args["--f"])
    L = float(args["--L"])
    rho = float(args["--rho"])
    N = int(rho * L**3)
    if f < 0.0 or f > 1.0:
        print("Fraction f of A beads must be between 0 and 1.")
        sys.exit()

    print("=== LAMMPS data file for binary mixture ====")
    print("L: %.1f | rho: %.1f | f: %.2f" % (L, rho, f))

    names, xyz = get_pos(N, f, L)
    header = header2str(N, L)
    final_string = header + \
                   mass2str(1.0, 1.0) + \
                   atoms2str(names, xyz)

    fname = args["--save"]
    open(fname, "w").write(final_string)
    print("Data file written in", fname)

    if args["--xyz"]:
        fname = args["--xyz"]
        ll.save_xyzfile(fname, np.hstack((np.matrix(names).T, xyz)) )
        print("xyz file written in", fname)


