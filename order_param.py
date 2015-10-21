#!/usr/bin/env python
"""Usage:
    order_param.py calc <systems> <rc> [--save <fname>] [--gnuplot]
    order_param.py plot <file>

Calculate the order parameter for a binary mixture phase transition
(number of AB contacts divited by number of AA plus BB contacts
in a given cutoff distance rc). 
Look for fA_*_AB_* dirs and collect all the dump files of the form "Dump/dump*.xyz".

Arguments:
    <systems>       The bin. mixture systems stored in dirs e.g. fA_0.1_AB_4.0, use regex
    <rc>            Cutoff distance in which to consider pairing

Options:
    --save <fname>  Save the final op matrix into file [default: temp.out]
    --gnuplot       Print in Gnuplot data format

pv278@cam.ac.uk, 12/10/15
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os, sys, glob
from docopt import docopt
import mat_ops


def read_dumpfile(dumpfile):
    """Read a LAMMPS xyz dumpfile into a numpy matrix"""
    A = open(dumpfile, "r").readlines()[2:]
    A = [line.split() for line in A]
    A = np.array(A, order="F").astype(float)
    return A


def save_matrix(mat, fname):
    m, n = mat.shape
    with open(fname, "w") as f:
        for i in range(m):
            for j in range(n):
                f.write(str(mat[i, j]) + "\t")
            f.write("\n")
    print "Matrix written in file", fname


def get_op(xyz_mat, rc):
    """Get number of 'contacts' between beads, DEPRECATED"""
    AB, AA, BB = 0, 0, 0
    N = len(xyz_mat)
    arr = mat_ops.get_pair_dist2(xyz_mat)        # call Fortran
    arr = arr[arr[:, 2] <= rc]
    mask = (arr[:, 0] == 1) & (arr[:, 1] == 1)
    AA = sum(arr[mask, 2])
    mask = (arr[:, 0] == 2) & (arr[:, 1] == 2)
    BB = sum(arr[mask, 2])
    AB = sum(arr[:, 2][arr[:, 0] != arr[:, 1]])
    return float(AB)/(AA+BB)


def get_average_op(dumpfiles, rc, fA):
    """From given xyz files extract average order param"""
    order_params = []
    for dumpfile in dumpfiles:
        A = read_dumpfile(dumpfile)
#        order_params.append(get_op(A, rc))                # common sense, NOT GOOD
#        order_params.append(mat_ops.get_local_op(A, rc))   # Goyal thesis, calling Fortran
        n1, n2 = mat_ops.get_local_op2(A, rc)               # Goyal thesis, alternative call
        n1, n2 = n1/fA, n2/(1-fA)                           # rescale to make sense
        op = float(np.dot(n1-n2, n1-n2))/np.dot(n1+n2, n1+n2)
        order_params.append(op)
    return np.average(order_params)


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args

    if args["calc"]:
        rc = float(args["<rc>"])
        systems = glob.glob(args["<systems>"])
        fAs = list(set([float(system.split("/")[0].split("_")[1]) for system in systems]))
        ABs = list(set([float(system.split("/")[0].split("_")[3]) for system in systems]))
        fAs.sort()
        ABs.sort()
        print fAs, "\n", ABs
 
        op_mat = np.zeros((len(ABs), len(fAs)))
        for i in range(len(ABs)):
            for j in range(len(fAs)):
                system = "fA_" + str(fAs[j]) + "_AB_" + str(ABs[i])
                path = os.path.join(os.path.expanduser("~/DPDcoeffs/BinMixt2"),\
                       system, "Dump/dump9?00.xyz")  # MAKE THIS REGEX CMD LN OPTION
                xyzfiles = glob.glob(path)
                op_mat[i, j] = get_average_op(xyzfiles, rc, fAs[j])
                if args["--gnuplot"]:
                    print ABs[i], fAs[j], op_mat[i, j]
        save_matrix(op_mat, args["--save"])

    elif args["plot"]:   # NOT TESTED YET
        data = np.loadtxt(args["<file>"])
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        Axes3D.plot_surface(data)
        plt.show()


