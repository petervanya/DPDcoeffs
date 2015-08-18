#!/usr/bin/env python
"""Usage:
    gen_gjf.py xyz <file> <outfile>
    gen_gjf.py g09out <file> <outfile>

Read coords from either xyz file or g09 outfile and create a new gjf file

Arguments:
    <file>         The input xyz file or g09 outfile to extract coords from
    <outfile>      File to output

pv278@cam.ac.uk, 10/08/15
"""
from docopt import docopt
import numpy as np
from g09lib import gen_g09_header, gen_g09_script
from xyzlib import Atoms
from parse_coords import parse_last_coords

SCFcycles = 1000
OPTcycles = 200
NP = 16
g09params = {"np" : NP, "method" : "B3LYP", "basis" : "6-31G*",\
             "text" : "Some silly text",\
             "opt" : "Opt=(maxcycles=" + str(OPTcycles) + ")",\
             "scf" : "Scf=(direct, maxcycle=" + str(SCFcycles) + ")"}


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args
    outfile = args["<outfile>"]
    header = gen_g09_header(params=g09params)

    if args["xyz"]:
        infile = args["<file>"]
        assert(infile[-3:] == "xyz")
        A = Atoms().read(infile)
        gen_g09_script(header, str(A), outfile)

    if args["g09out"]:
        infile = args["<file>"]
        assert(infile[-3:] == "out")
        A = parse_last_coords(infile)
        gen_g09_script(header, str(A), outfile)

