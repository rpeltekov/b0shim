from pydicom import dcmread 
import dicom2nifti      
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="2 folders of phase compute"
                                                 +" phasemap")
    parser.add_argument("folder1", type=str,default=0, help="")
    parser.add_argument("folder2", type=str,default=0, help="")
    #parser.add_argument("o", type=str,default=0, help="")
    #parser.add_argument("thresh", type=int,default=0, help="% value 0 to 100")
    parser.add_argument("--debug", action="store_true", help="")

    args = parser.parse_args()

    folder1 = os.path.join("processed", args.folder1)
    folder2 = os.path.join("processed", args.folder2)
    b0_fold = "b0maps"
    b0_dir = os.path.join(os.getcwd(), b0_fold)
    if not os.path.exists(b0_dir):
        os.makedirs(b0_dir)

    #o_filename = args.o

    num_slice = len(os.listdir(folder1)) // 2
    phase_pref = "ph.MRDC."
    mag_pref = "mag.MRDC."

    for i in range(num_slice):

        mag1file = os.path.join(folder1, mag_pref+str(i))
        phase1file = os.path.join(folder1, phase_pref+str(i))
        phase2file = os.path.join(folder2, phase_pref+str(i))

        print(mag1file)
        mag1 = dcmread(mag1file)
        phase1 = dcmread(phase1file)
        phase2 = dcmread(phase2file)
        
        mag1a = mag1.pixel_array
        phase1a = phase1.pixel_array
        phase2a = phase2.pixel_array

        te1 = phase1[0x0018, 0x0081].value
        te2 = phase2[0x0018, 0x0081].value
        print(f"[INFO] TE1 = {te1} TE2 = {te2}")

        phase1a = phase1a * np.pi / 32767
        phase2a = phase2a * np.pi / 32767

        delt_phase = phase2a - phase1a

        m, n = mag1a.shape
        mask = (mag1a > np.max(np.abs(mag1a)) / 10).astype(int)
        
        deltafreq = delt_phase / (2*np.pi*(te2-te1)*10e-3)

        gradient = np.linspace(0, 1, 256)
        colors = np.zeros((256, 3))
        colors[:, 1] = gradient  # Green channel from 0 to 1
        colors[:, 0] = 1 - gradient  # Red channel from 1 to 0

        plt.imshow(deltafreq, cmap='gray')
        plt.colorbar()
        plt.show()

        
        

        
        






