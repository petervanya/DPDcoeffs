#!/usr/bin/env python
"""Usage:
    analyse_curves.py <infiles> <Efile> [--T <T>] [--plot] [--fit [--cutoff <rc>]]

1. Read in all the energy curves for all blob combinations
and produce a temperature-weighted master curve
2. Fit the master curve by a parabola

Options:
    --T <T>             Temperature [default: 300]
    --plot              Plot master curve
    --fit               Fit the curve within a given distance cutoff
    --cutoff <rc>       Distance cutoff for fitting [default: 6e-10]

pv278@cam.ac.uk
"""
import numpy as np
import matplotlib.pyplot as plt
import glob, sys
from docopt import docopt


kB = 1.38e-23
J_H = 1.602e-19 * 27.211     # Hartrees to Joules
m_AA = 1e-10


def get_master_curve(infiles, E, T):
    """Produce master curve from list of curves, energies of separate clusters
    and temperature
    Arguments:
    * infiles: list of files storing (N, 2) matrix of energy (eV) vs. distance (AA)
    * E: (N, 2) array of blob number and energy
    * T: temperature
    """
    master_curve = np.asarray(np.loadtxt(infiles[0]))
    master_curve[:, 1] = 0.0
    Z = 0.0                  # normalisation constant
    E *= J_H

    for i in range(len(infiles)):
        blob1, blob2 = [int(n) for n in infiles[i].rstrip(".out").split("_")[1:]]
        c = np.loadtxt(infiles[i])
        scaled_E = (E[blob1, 1] + E[blob2, 1]) - 2*np.min(E[:, 1])  # rescale wrt energy minimum
        boltzmann_fact = np.exp(-(scaled_E / (kB*T)))
        master_curve[:, 1] += (c[:, 1] - (E[blob1, 1] + E[blob2, 1])) * boltzmann_fact
        Z += boltzmann_fact
    
    master_curve[:, 1] /= Z
    return master_curve


def fit_curve(mc, cutoff):
    """
    Fit master curve by a quadratic function
    Arguments:
    * cutoff: discard all values beyond this one
    TO DO: recreate the right side of the parabola for correct fit
    """
    new_c = mc[mc[:, 0] <= cutoff]
    coeffs = np.polyfit(new_c[:, 0], new_c[:, 1], 2)
    return coeffs


if __name__ == "__main__":
    args = docopt(__doc__)
#    print args
    infiles = glob.glob(args["<infiles>"])
    energies = np.loadtxt(args["<Efile>"])
    T = float(args["--T"])
    
    if not infiles:
        print "No files captured, aborting."
        sys.exit()
    print len(infiles), "energy curves will be analysed."

    master_curve = get_master_curve(infiles, energies, T)
    master_curve[:, 0] *= m_AA
    master_curve[:, 1] *= J_H
    
    print master_curve
    outname = "master_curve.out"
    np.savetxt(outname, master_curve)
    print "Master curve saved in", outname
    
    if args["--fit"]:
        cutoff = float(args["--cutoff"])
        print "Fitting curve, cutoff:", cutoff
        coeffs = fit_curve(master_curve, cutoff)
        print "Fit:", coeffs
        print "The a_ij coeff: %.4f" % (coeffs[0]/2)


    if args["--plot"]:
        plt.plot(master_curve[:, 0], master_curve[:, 1])
        plt.plot(master_curve[:, 0], np.polyval(coeffs, master_curve[:, 0]), "r-")
        plt.xlim([4e-10, 8e-10])
#        plt.ylim([0.0, 1e-19])
        plotname = "plot_master_curve.png"
        plt.savefig(plotname)
        print "Plot of master curve saved in", plotname
 



