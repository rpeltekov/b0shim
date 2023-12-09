from pydicom import dcmread 
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from apng import APNG

def get_ticks(volume):
    start = volume[0]
    size = volume[1]
    end = start + size
    res = volume[2]
    x = np.linspace(start[0], end[0], int(res[0])+1)[:-1]
    y = np.linspace(start[1], end[1], int(res[1])+1)[:-1]
    z = np.linspace(start[2], end[2], int(res[2])+1)[:-1]

    x = x + (x[1] - x[0])/2
    y = y + (y[1] - y[0])/2
    z = z + (z[1] - z[0])/2
    return x, y, z

def MRPhaseUnwrap(theta, voxelsize=[1.5, 1.5, 1.6], padsize=[12, 12, 12]):
    # Make the size of the third dimension even
    print("before unwrap size", theta.shape)
    made_even = False
    if theta.shape[2] % 2 == 1:
        theta = np.concatenate((theta, np.zeros_like(theta[:, :, -1, np.newaxis])), axis=2)
        made_even = True

    # Pad the input data with zeros
    theta = np.pad(theta, [(padsize[0], padsize[0]), (padsize[1], padsize[1]),
                           (padsize[2], padsize[2]), (0,0)], mode='constant')

    # Get the size of the padded data
    NP = theta.shape
    if len(theta.shape) == 3:
        print(theta.shape)
        theta = theta[:,:,:,np.newaxis]
        NP = theta.shape

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
        # centered ifft
        scale = np.sqrt(len(d.flatten()))
        return np.fft.fftshift(np.fft.ifftn(np.fft.fftshift(d)))*scale

    def fftnc(d):
        # centered fft
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
    print("after unwrap size", phi.shape)

    return phi, Laplacian


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="3d basis")
    parser.add_argument("folder1", type=str,default=0, help="")
    parser.add_argument("folder2", type=str,default=0, help="")
    parser.add_argument("--voxX", type=int,default=1.5625, help="")
    parser.add_argument("--voxY", type=int,default=1.5625, help="")
    parser.add_argument("--voxZ", type=int,default=2, help="")
    parser.add_argument("--debug", action="store_true", help="")
    parser.add_argument("--channel", action="store_true", help="")
    parser.add_argument("--nocrop", action="store_true", help="")
    parser.add_argument("-R", type=float,default=1.6, help="")

    args = parser.parse_args()

    folder1 = os.path.join("processed", args.folder1)
    folder2 = os.path.join("processed", args.folder2)

    b0_fold = "b0maps"
    b0_dir = os.path.join(os.getcwd(), b0_fold)

    comb_fold = f"{args.folder1}_{args.folder2}"
    comb_dir = os.path.join(b0_dir, comb_fold)

    if not os.path.exists(comb_dir):
        os.makedirs(comb_dir)

    if not os.path.exists(os.path.join(comb_dir,"basis.npy")):
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
        #unwrapf1, lap1 = MRPhaseUnwrap(phases1, [args.voxX, args.voxY, args.voxZ])
        #unwrapf2, lap2 = MRPhaseUnwrap(phases2, [args.voxX, args.voxY, args.voxZ])

        unwrapf1 = phases1
        unwrapf2 = phases2

        delt_phase1 = unwrapf1[:,:,:,1] - unwrapf1[:,:,:,0]
        delt_phase2 = unwrapf2[:,:,:,1] - unwrapf2[:,:,:,0]

        deltfreq1 = delt_phase1 / (2*np.pi*(tes[1]-tes[0])*10e-3)
        deltfreq2 = delt_phase2 / (2*np.pi*(tes[1]-tes[0])*10e-3)

        tms_basis = deltfreq2-deltfreq1
        np.save(os.path.join(comb_dir,"basis.npy"), tms_basis)
    else:
        tms_basis = np.load(os.path.join(comb_dir,"basis.npy"))
    
    b0offres = os.path.join(comb_dir, "differencemap")
    if not os.path.exists(b0offres):
        os.makedirs(b0offres)

    #TODO: change this ifelse shit and fix naming from tms_basis...
    
    if (not args.channel):
        for slice in range(deltfreq1.shape[2]):

            fig, ax = plt.subplots()
            #fig, ax = plt.subplots(1,3, figsize=(15, 5))
            #z = (slice-deltfreq1.shape[2]/2)* 1.6/10
            #plt.title(f"Off resonance, at z={z}cm")
            #im0 = ax[0].imshow(deltfreq1[:,:,slice], cmap='jet')
            #ax[0].set_title(f"Hz of {args.folder1}")
            #fig.colorbar(im0, ax=ax[0])
            #im1 = ax[1].imshow(deltfreq2[:,:,slice], cmap='jet')
            #ax[1].set_title(f"Hz of {args.folder2}")
            #fig.colorbar(im1, ax=ax[1])
            im2 = ax.imshow(tms_basis[:,:,slice], cmap='jet')
            fig.colorbar(im2, ax=ax)
            ax.set_title(f"Off Resonance Due to TMS (Hz)") 
            #plt.show()
            plt.savefig(os.path.join(b0offres, f"offres_{z}z.png"))
            plt.close()

            
            plt.imshow(tms_basis[:,:,slice], cmap='jet')
            plt.title(f"offres due to tms z = {z}")
            plt.clim(-100, 100)
            plt.colorbar()
            plt.savefig(os.path.join(b0offres, f"{slice}.png"))
            plt.close()
    else:

        # This is to simply output the b0maps for every slice in the difference
        # image from the two inputs given.
            
        volume = np.array([[-12.0, -12.0, -9.919999999999999929e+00],
                  [24.0, 24.0, 19.839],
                  [256, 256, 124]])
        x, y, z = get_ticks(volume)

        axial = os.path.join(b0offres, "axial")
        if not os.path.exists(axial):
            os.makedirs(axial)

        apng = APNG()
        print(f'[INFO] Max field value seen {np.max(np.abs(tms_basis))}')
        for sl in range(tms_basis.shape[2]):
            xx, yy = np.meshgrid(x, y)
            z_lev = sl #z[sl]

            radius = np.sqrt(xx**2 + yy**2 + z_lev**2*np.ones(yy.shape))

            in_phantom = np.argwhere(radius.flatten() <= args.R).flatten()

            if len(in_phantom) == 0 and not args.nocrop:
                continue

            num_pos = int(tms_basis.shape[0]*tms_basis.shape[1])
            refslice = tms_basis[:,:,sl]
            zsl = refslice.reshape((num_pos,))
            alphas = np.zeros((num_pos,))
            alphas[in_phantom] = 1.0
            zsl = zsl.reshape(refslice.shape)
            alphas = alphas.reshape(refslice.shape)
            if args.nocrop:
                alphas = None

            fig, ax = plt.subplots()
            im = ax.imshow(zsl, alpha=alphas, cmap='jet')
            im.set_clim(-300, 300)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            plt.axis('off')
            plt.title(f'Axial Z ~ {round(z_lev)}')
            imgname = os.path.join(axial, f"{sl}.png")#_z{z_lev.round(2)}.png")
            plt.savefig(imgname, transparent=True)

            apng.append_file(imgname)
            plt.close()

        apng.save(os.path.join(axial, f'anim.apng'))
