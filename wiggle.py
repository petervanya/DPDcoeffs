#!/usr/bin/env python
"""Usage:
    wiggle.py <infile> (--all | --bymol) [-s <outfile>] [--sigma <s>] [--seed <seed>]

Wiggle atomic coordinates by a Gaussian noise with a given sigma

Arguments:
    <infile>         Input xyz file

Options:
    -s <outfile>     Save into xyz file
    --sigma <s>      Gaussian sigma [default: 0.1]
    --seed <seed>    Random seed [default: 123]
    --all            Wiggle all atoms randomly
    --bymol          Wiggle each water molecule separately

pv278@cam.ac.uk, 14/08/15
"""
from docopt import docopt
import numpy as np
from math import degrees
from xyzlib import Atoms


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args
    np.random.seed(int(args["--seed"]))
    sigma = float(args["--sigma"])
    A = Atoms().read(args["<infile>"])

    if args["--all"]:
        noise = np.random.randn(len(A), 3)*sigma
        A.coords += noise
    if args["--bymol"]:
        Nmols = len(A)/3
        for i in range(Nmols):
            mol = Atoms(coords=A.coords[i:i+3, :])
            theta = degrees(np.random.randn()*sigma)
            phi = degrees(np.random.randn()*sigma)
            mol.rotate(theta, phi)
            mol.shift(np.random.randn(3, 3)*sigma)
            A.coords[i:i+3, :] = mol.coords

    if args["-s"]:
        A.save(args["-s"])
    else:
        print A

