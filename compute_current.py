import argparse
import cv2
import numpy as np
import sys, os
import biot_savart_v4_3 as bs
import cvxopt
from matplotlib import pyplot as plt
from matplotlib.colors import TwoSlopeNorm as diverge

# given the basis functions of the coils at 1 amp, and the offresonance artifact
# both should be in the same FOV and RES
# compute the current values for each of the coils and output
# save them to a
# output the estimated corrected basis map

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("coilsfolder", default="default", help="")
    parser.add_argument("offres", default="default", help="")
    parser.add_argument("volume", default="default", help="")
    parser.add_argument("--numcoils", type=int, default=15, help="")
    parser.add_argument("--R", type=float, default=18., help="")

    args = parser.pargse_args()


    offres = np.load(args.offres)

    coilsfolder = os.path.join("headcaps", args.coilsfolder)
    volume = np.loadtxt(os.path.join(coilsfolder, 'volume_def.txt'))

    bases = []
    for coil in os.listdir(coilsfolder):
        else:
            basis = np.loadtxt(os.path.join(coilsfolder, coil, "field.csv"), delimiter=',')[:,2]
            basis.reshape(volume[1])
            bases.append(basis)
    basis = np.array(bases)


    if basis[0].shape != offres.shape:
        upsamples = offress.shape / basis.shape[1:]
        re = offress.shape % basis.shape[1:]

        if re != 0:
            raise Exception("need integer upscale factor")
        for i in range(len(basis)):
            for j in range(len(upsamples)):
                basis[i] = np.repeat(basis[i], upsamples[j], axis=j)

        assert basis[0].shape == offres.shape, "still not upsample correct"

    # masking 
    start = volume[0]
    end = volume[0] + volume[1]
    x = np.linspace(start[0], end[0], volume[0]+1)[:-1]
    y = np.linspace(start[1], end[1], volume[1]+1)[:-1]
    z = np.linspace(start[2], end[2], volume[2]+1)[:-1]

    x = x + (x[1] - x[0])/2
    y = y + (y[1] - y[0])/2
    z = z + (z[1] - z[0])/2
    


    



