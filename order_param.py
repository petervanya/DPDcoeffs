#!/usr/bin/env python
"""Usage:
    order_param.py calc <systems> <rc> [--save <fname>]
    order_param.py plot <file>

Calculate the order parameter for a binary mixture phase transition
(number of AB contacts divited by number of AA plus BB contacts
in a given cutoff distance rc).

Arguments:
    <systems>       The bin. mixture systems stored in dirs e.g. fA_0.1_AB_4.0, use regex
    <rc>            Cutoff distance in which to consider pairing

Options:
    --save <fname>  Save the final op matrix into file [default: temp.out]

pv278@cam.ac.uk, 12/10/15
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os, sys, glob
from docopt import docopt
import mat_ops


def read_outfile(outfile):
    """Read a LAMMPS xyz outfile into a numpy matrix"""
    A = open(outfile, "r").readlines()[2:]
    A = [line.split() for line in A]
    A = np.array(A, order="F").astype(float)
    return A


def get_op(xyz_mat, rc):
    """Get number of 'contacts' between beads"""
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


def get_average_op(outfiles, rc):
    """From given xyz files extract average order param"""
    order_params = []
    for outfile in outfiles:
        A = read_outfile(outfile)
        order_params.append(get_op(A, rc))
    return np.average(order_params)


def savemat(mat, fname):
    m, n = mat.shape
    with open(fname, "w") as f:
        for i in range(m):
            for j in range(n):
                f.write(str(mat[i, j]) + "\t")
            f.write("\n")
    print "Matrix written in file", fname


def pipeline_old(outfiles, rc):
    """(Old working of this script)
    From the outfiles extract the order param and print it"""
    print get_average_op(outfiles, rc)


if __name__ == "__main__":
    args = docopt(__doc__)
    print args

    if args["calc"]:
        rc = float(args["<rc>"])
##       outfiles = glob.glob(args["<files>"])
        systems = glob.glob(args["<systems>"])
##       print systems
        fAs = list(set([float(system.split("_")[1]) for system in systems]))
        ABs = list(set([float(system.split("_")[3]) for system in systems]))
        fAs.sort()
        ABs.sort()
        print fAs, "\n", np.array([ABs]).T
 
        op_mat = np.zeros((len(ABs), len(fAs)))
        for i in range(len(ABs)):
            for j in range(len(fAs)):
##               print fAs[j], ABs[i]
                system = "fA_" + str(fAs[j]) + "_AB_" + str(ABs[i])
                path = os.path.join(os.path.expanduser("~/DPDcoeffs/BinMixt"), system, "Dump/dump*00.xyz")
                xyzfiles = glob.glob(path)
                op_mat[i, j] = get_average_op(xyzfiles, rc)
 
        if args["--save"]:
            savemat(op_mat, args["--save"])
        else:
            print op_mat

    elif args["plot"]:
        data = np.loadtxt(args["<file>"])
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        Axed3D.plot_surface(data)
        plt.show()



