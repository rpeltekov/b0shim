from pydicom import dcmread 
import dicom2nifti      
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def MRPhaseUnwrap(theta, voxelsize=[1.5, 1.5, 1.6], padsize=[12, 12, 12]):
    # Pad the input data with zeros
    theta = np.pad(theta, [(padsize[0], padsize[0]), (padsize[1], padsize[1]),
                           (padsize[2], padsize[2]), (0,0)], mode='constant')

    # Get the size of the padded data
    NP = theta.shape
    if len(theta.shape) == 3:
        print(theta.shape)
        print("we should not be here")
        theta = theta[:,:,:,np.newaxis]
        NP = theta.shape

    # Make the size of the third dimension even
    if theta.shape[2] % 2 == 1:
        print("we should not be here")
        theta = np.concatenate((theta, np.zeros_like(theta[:, :, -1, np.newaxis])), axis=2)

    # Create coordinates and wave numbers
    yy, xx, zz = np.meshgrid(np.arange(1, NP[1] + 1), 
                             np.arange(1, NP[0] + 1), 
                             np.arange(1, NP[2] + 1))
    if len(voxelsize) > 1:
        voxelsize = np.array(voxelsize)

    FOV = voxelsize * NP[0:3]
    xx = (xx - NP[0] / 2 - 1) / FOV[0]
    yy = (yy - NP[1] / 2 - 1) / FOV[1]
    zz = (zz - NP[2] / 2 - 1) / FOV[2]
    k2 = xx**2 + yy**2 + zz**2

    # Initialize arrays for phi and Laplacian
    phi = np.zeros(NP)
    Laplacian = np.zeros(NP).astype(np.cdouble)

    def ifftnc(d):
        scale = np.sqrt(len(d.flatten()))
        return np.fft.fftshift(np.fft.ifftn(np.fft.fftshift(d)))*scale

    def fftnc(d):
        scale = np.sqrt(len(d.flatten()))
        return np.fft.fftshift(np.fft.fftn(np.fft.fftshift(d)))/scale

    for iii in range(NP[3]):
        Laplacian0 = np.cos(theta[:,:,:,iii]) * ifftnc(k2 * fftnc(np.sin(theta[:,:,:,iii])))
        Laplacian0 = Laplacian0 - np.sin(theta[:,:,:,iii]) * ifftnc(k2*fftnc(np.cos(theta[:,:,:,iii])))
        Laplacian[:,:,:,iii] = Laplacian0
        phi0 = fftnc(Laplacian0) / k2
        phi0[NP[0] // 2, NP[1] // 2, NP[2] // 2] = 0
        phi0 = np.real(ifftnc(phi0))
        phi[:,:,:,iii] = phi0

    # Remove padding to get the final results
    phi = phi[padsize[0]:-padsize[0], padsize[1]:-padsize[1], padsize[2]:-padsize[2], :]
    Laplacian = Laplacian[padsize[0]:-padsize[0], padsize[1]:-padsize[1], padsize[2]:-padsize[2], :]

    return phi, Laplacian


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="3d basis")
    parser.add_argument("folder1", type=str,default=0, help="")
    parser.add_argument("folder2", type=str,default=0, help="")
    parser.add_argument("--voxX", type=int,default=1.5, help="")
    parser.add_argument("--voxY", type=int,default=1.5, help="")
    parser.add_argument("--voxZ", type=int,default=1.6, help="")
    parser.add_argument("--debug", action="store_true", help="")

    args = parser.parse_args()

    folder1 = os.path.join("processed", args.folder1)
    folder2 = os.path.join("processed", args.folder2)

    b0_fold = "b0maps"
    b0_dir = os.path.join(os.getcwd(), b0_fold)

    comb_fold = f"{args.folder1}_{args.folder2}"
    comb_dir = os.path.join(b0_dir, comb_fold)

    if not os.path.exists(comb_dir):
        os.makedirs(comb_dir)

    tes = np.zeros(2)
    for file_name in os.listdir(folder1):
        if file_name.endswith("_3dPhase.npy"):
            # Extract TE value from the file name
            tes[0] = float(file_name.split("_")[1])
            tes[1] = float(file_name.split("_")[2])
            # Load and store the phase array
            phases1 = np.load(os.path.join(folder1, file_name))
         
            phases2 = np.load(os.path.join(folder2, file_name))

    # unwrap
    unwrapf1, lap1 = MRPhaseUnwrap(phases1, [args.voxX, args.voxY, args.voxZ])
    unwrapf2, lap2 = MRPhaseUnwrap(phases2, [args.voxX, args.voxY, args.voxZ])

    delt_phase1 = unwrapf1[:,:,:,1] - unwrapf1[:,:,:,0]
    delt_phase2 = unwrapf2[:,:,:,1] - unwrapf2[:,:,:,0]

    deltfreq1 = delt_phase1 / (2*np.pi*(tes[1]-tes[0])*10e-3)
    deltfreq2 = delt_phase2 / (2*np.pi*(tes[1]-tes[0])*10e-3)

    tms_basis = deltfreq2-deltfreq1
    np.save(os.path.join(comb_dir,"basis.npy"), tms_basis)

    for slice in range(deltfreq1.shape[2]):
        fig, ax = plt.subplots(1,3, figsize=(15, 5))
        z = (slice-deltfreq1.shape[2]/2)* 1.6/10
        plt.title(f"Off resonance, at z={z}cm")
        im0 = ax[0].imshow(deltfreq1[:,:,slice], cmap='RdYlGn')
        ax[0].set_title(f"Hz of {args.folder1}")
        fig.colorbar(im0, ax=ax[0])
        im1 = ax[1].imshow(deltfreq2[:,:,slice], cmap='RdYlGn')
        ax[1].set_title(f"Hz of {args.folder2}")
        fig.colorbar(im1, ax=ax[1])
        im2 = ax[2].imshow(tms_basis[:,:,slice], cmap='RdYlGn')
        fig.colorbar(im2, ax=ax[2])
        ax[2].set_title(f"Off Resonance Due to TMS (Hz)") 
        #plt.show()
        plt.savefig(os.path.join(comb_dir, f"offres_{z}z.png"))
        plt.close()
