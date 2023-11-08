import argparse
import numpy as np
import sys, os
from cvxopt import solvers, matrix
from matplotlib import pyplot as plt

# given the basis functions of the coils at 1 amp, and the offresonance artifact
# both should be in the same FOV and RES
# compute the current values for each of the coils and output
# save them to a
# output the estimated corrected basis map

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("coilsfolder", default="default", help="")
    parser.add_argument("offres", default="default", help="")
    parser.add_argument("--gradient", action='store_true', help="")
    parser.add_argument("--perslice", action='store_true', help="")
    parser.add_argument("--roi", type=float, default=8.0, help="")
    parser.add_argument("-v", action='store_true', help="")

    args = parser.parse_args()

    headcap = args.coilsfolder
    obj = args.offres.split("/")
    obj = "_".join(obj[0:-1])

    R = args.roi

    infodir = f"shim_results/"+obj+f"-roi{R}-"+headcap
    outputdir = f"shim_results/"+obj+f"-roi{R}-"+headcap+("/perslice/" if args.perslice else "")
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    b0offres = os.path.join(outputdir, "streamline")
    if not os.path.exists(b0offres):
        os.makedirs(b0offres)

    offres = np.load(args.offres)
    if args.v:
        print(f"[INFO] offres: {offres.shape}")

    coilsfolder = os.path.join("headcaps", args.coilsfolder)
    volume = np.loadtxt(os.path.join(coilsfolder, 'volume_def.txt'))
    numcoils = int(volume[3,2])

    print(f"SHIMMING ROI RADIUS {R} cm with {numcoils} coils")

    if args.v:
        print("[INFO] collecting the basis functions")
    bases = []
    for i in range(numcoils):
        coiltxt = os.path.join(coilsfolder, f"coil_{i}", "field.csv")
        coilnpy = os.path.join(coilsfolder, f"coil_{i}", "field.npy")
        if not os.path.exists(coilnpy):
            basis = np.loadtxt(coiltxt, delimiter=',')[:,2]
            np.save(coilnpy, basis)
        else:
            basis = np.load(coilnpy)
        basis = basis.reshape(volume[2].astype(int))
        bases.append(basis)
    basis = np.array(bases)

    def upsample(this, tothis):
        if args.v:
            print(f"[INFO] tothis: {tothis.shape}, this[0]: {this[0].shape}")
            print("[INFO] upsampling this functions")
        upsamples = np.array(tothis.shape) / np.array(this.shape[1:])
        re = np.array(tothis.shape) % np.array(this.shape[1:])
        upsamples = upsamples.astype(int)
        if args.v:
            print(upsamples)

        if (re != np.zeros(re.shape)).any():
            raise Exception("need integer upscale factor")
        for j in range(len(upsamples)):
            newthis = np.zeros((this.shape[0],
                                 this.shape[1]*(upsamples[j] if j == 0 else 1),
                                 this.shape[2]*(upsamples[j] if j == 1 else 1),
                                 this.shape[3]*(upsamples[j] if j == 2 else 1)))
            for i in range(numcoils):
                newthis[i] = np.repeat(this[i], upsamples[j], axis=j)
            this = newthis.copy()

        assert this[0].shape == tothis.shape, "still not upsample correct"
        return this


    if basis[0].shape != offres.shape:
        basis = upsample(basis, offres)

    # masking 
    if not args.perslice:
        if args.v:
            print("[INFO] forming QP problem")
        num_pos = int(offres.shape[0]*offres.shape[1]*offres.shape[2])
        A = basis.reshape((numcoils, num_pos)).T
        Y = offres.reshape((num_pos,))

        if args.v:
            print(f"[INFO] A: {A.shape}, y: {Y.shape}")

        start = volume[0].astype(int)
        end = volume[0].astype(int) + volume[1].astype(int)
        xx = np.linspace(start[0], end[0], int(offres.shape[0])+1)[:-1]
        yy = np.linspace(start[1], end[1], int(offres.shape[1])+1)[:-1]
        zz = np.linspace(start[2], end[2], int(offres.shape[2])+1)[:-1]

        xx = xx + (xx[1] - xx[0])/2
        yy = yy + (yy[1] - yy[0])/2
        zz = zz + (zz[1] - zz[0])/2

        xx, yy, zz = np.meshgrid(xx, yy, zz)
        xx = xx.flatten()
        yy = yy.flatten()
        zz = zz.flatten()

        radius = np.sqrt(xx**2 + yy**2 + zz**2)
        in_phantom = np.argwhere(radius < R).flatten()
        if args.v:
            print(in_phantom)

        A_masked = A[in_phantom]
        Y_masked = Y[in_phantom]
        if args.v:
            print(f"[info] a_masked: {A_masked.shape}, y_masked: {Y_masked.shape}")
            print(f"[info] a_masked: {A_masked}, y_masked: {Y_masked}")

        p = 2 * A_masked.T @ A_masked
        q = 2 * Y_masked.T @ A_masked #todo: add lasso term
        g = np.vstack((np.eye(numcoils), -np.eye(numcoils)))
        if args.gradient:
            h = 2 * np.ones(numcoils-3)
            h = np.concatenate((np.ones(3)*3.2934, h))
            h = np.concatenate((h, h))
        else:
            h = 2 * np.ones(2*numcoils)

        if args.v:
            print(f"[info] p: {p.shape}\n       q: {q.shape}\n"+
              f"       g: {g.shape}\n       h: {h.shape}")
        res = solvers.qp(matrix(p), matrix(q), matrix(g), matrix(h))

        if args.v:
            print("[info] solved currents: \n", res['x'])

        x = np.array(res['x']).flatten()
        np.savetxt(os.path.join(outputdir, "currents.txt"), x)

        if args.v:
            print("[info] applying solved currents")

        corrections = A_masked.dot(x).flatten()
        if args.v:
            print(f"[info] corrections: {corrections.shape}")
            print(f"[info] y_masked: {Y_masked.shape}")
        
        outputmap = np.zeros(Y.shape)
        reference = np.zeros(Y.shape)
        if args.v:
            print(f"[info] inputs: {outputmap[in_phantom].shape}")
        outputmap[in_phantom] = corrections + Y_masked
        reference[in_phantom] = Y_masked
        outputmap = outputmap.reshape(offres.shape)

        zz = zz.reshape(offres.shape)
        reference = reference.reshape(offres.shape)
        ma = max(np.max(outputmap), np.max(reference))
        mi = min(np.min(outputmap), np.min(reference))
        for z in range(offres.shape[2]):
            fig, ax = plt.subplots(1,2, figsize=(10, 5))

            im0 = ax[0].imshow(reference[:,:,z], cmap='jet')
            ax[0].set_title("original")
            fig.colorbar(im0, ax=ax[0], fraction=0.046, pad=0.04)
            im0.set_clim(mi,ma)
            
            im1 = ax[1].imshow(outputmap[:,:,z], cmap='jet')
            ax[1].set_title("corrections")
            fig.colorbar(im1, ax=ax[1], fraction=0.046, pad=0.04)
            im1.set_clim(mi,ma)

            plt.savefig(os.path.join(outputdir, f"{z}corr_z_{zz[0,0,z].round(2)}.png"))
            plt.close()

            plt.imshow(outputmap[:,:,z], cmap='jet')
            plt.title(f"corrected fieldmap at z = {zz[0,0,z]}")
            plt.colorbar(fraction=0.046, pad=0.04)
            plt.clim(mi,ma)
            plt.savefig(os.path.join(b0offres, f"{z}.png"))
            plt.close()
    else:
        start = volume[0].astype(float)
        end = volume[0].astype(float) + volume[1].astype(float)

        xx = np.linspace(start[0], end[0], int(offres.shape[0])+1)[:-1]
        yy = np.linspace(start[1], end[1], int(offres.shape[1])+1)[:-1]
        zz = np.linspace(start[2], end[2], int(offres.shape[2])+1)[:-1]

        xx = xx + (xx[1] - xx[0])/2
        yy = yy + (yy[1] - yy[0])/2
        zz = zz + (zz[1] - zz[0])/2

        xx, yy = np.meshgrid(xx, yy)
        xx = xx.flatten()
        yy = yy.flatten()

        ma = 0
        mi = 0
        outputmaps = []
        references = []
        zis = []
        currents = []

        solvers.options['show_progress'] = False

        print("[INFO] forming QP problem PER slice")
        for zi in range(offres.shape[2]):
            z = zz[zi]
            radius = np.sqrt(xx**2 + yy**2 + z**2*np.ones(yy.shape))

            basis_slice = basis[:,:,:,zi]
            ref_slice = offres[:,:,zi]

            in_phantom = np.argwhere(radius < R).flatten()
            if len(in_phantom) == 0:
                continue
            if args.v:
                print(in_phantom)

            num_pos = int(offres.shape[0]*offres.shape[1])
            A = basis_slice.reshape((numcoils, num_pos)).T
            y = ref_slice.reshape((num_pos,))

            if args.v:
                print(f"[INFO] A: {A.shape}, y: {y.shape}")

            A_masked = A[in_phantom]
            y_masked = y[in_phantom]
            if args.v:
                print(f"[INFO] A_masked: {A_masked.shape}, y_masked: {y_masked.shape}")

            P = 2 * A_masked.T @ A_masked
            q = 2 * y_masked.T @ A_masked #TODO: add lasso term
            G = np.vstack((np.eye(numcoils), -np.eye(numcoils)))
            if args.gradient:
                h = 2 * np.ones(numcoils-3)
                h = np.concatenate((np.ones(3)*3.2934, h))
                h = np.concatenate((h, h))
            else:
                h = 2 * np.ones(2*numcoils)

            if args.v:
                print(f"[INFO] P: {P.shape}\n       q: {q.shape}\n"+
                  f"       G: {G.shape}\n       h: {h.shape}")
            res = solvers.qp(matrix(P), matrix(q), matrix(G), matrix(h))

            if args.v:
                print("[INFO] Solved Currents: \n", res['x'])

            x = np.array(res['x']).flatten()
            np.savetxt(os.path.join(outputdir, "currents.txt"), x)

            if args.v:
                print("[INFO] Applying solved currents")

            corrections = A_masked.dot(x).flatten()
            if args.v:
                print(f"[INFO] corrections: {corrections.shape}")
            if args.v:
                print(f"[INFO] y_masked: {y_masked.shape}")
            
            outputmap = np.zeros(y.shape)
            reference = np.zeros(y.shape)
            if args.v:
                print(f"[INFO] inputs: {outputmap[in_phantom].shape}")
            outputmap[in_phantom] = corrections + y_masked
            reference[in_phantom] = y_masked

            outputmaps.append(outputmap.reshape(ref_slice.shape))
            references.append(reference.reshape(ref_slice.shape))
            zis.append(zi)
            currents.append(x)

        print('[INFO] solved all slices')
            
        outputmaps = np.array(outputmaps)
        references = np.array(references)
        currents = np.array(currents)
        
        xx = np.linspace(start[0], end[0], int(offres.shape[0])+1)[:-1]
        yy = np.linspace(start[1], end[1], int(offres.shape[1])+1)[:-1]
        XX, YY, ZZ = np.meshgrid(xx,yy,zz)
        XX = XX.flatten()
        YY = YY.flatten()
        ZZ = ZZ.flatten()
        RR = np.argwhere(np.sqrt(XX**2 + YY**2 + ZZ**2) < R)

        outputmaps_f = outputmaps.flatten()
        references_f = references.flatten()

        rmse_outputmaps = np.sqrt(1/len(RR) * np.sum((outputmaps_f)**2))
        rmse_references = np.sqrt(1/len(RR) * np.sum((references_f)**2))
        rmse_difference = rmse_references - rmse_outputmaps
        percent = (rmse_difference/rmse_references*100).round(2)
        print(f'[INFO] RMSE reference ROI = {rmse_references}')
        print(f'[INFO] RMSE shimmed ROI   = {rmse_outputmaps}')
        print(f'[INFO] increase in RMSE   = {rmse_difference}')
        print(f'[INFO] %increase in RMSE  = {percent}%')
        
        np.savetxt(os.path.join(infodir, f"performance_{percent}.txt"),
                   [R, rmse_references, rmse_outputmaps, rmse_difference, percent])

        ma = max(np.max(outputmaps), np.max(references))
        mi = min(np.min(outputmaps), np.min(references))

        np.save(os.path.join(infodir, "corrected_cropped_offres.npy"), outputmaps)
        np.save(os.path.join(infodir, "original_cropped_offres.npy"), references)

        print('[INFO] Saving figures')

        for i, zi in enumerate(zis):
            fig, ax = plt.subplots(1,2, figsize=(14, 5))

            im0 = ax[0].imshow(references[i], cmap='jet')
            ax[0].set_title("original off resonance (Hz)")
            ax[0].set_xlabel("X (cm)")
            ax[0].set_xticks(np.linspace(-.5,int(offres.shape[0])-.5,9), 
                         np.linspace(start[0], end[0], 9))
            ax[0].set_ylabel("Y (cm)")
            ax[0].set_yticks(np.linspace(-.5,int(offres.shape[1])-.5,9), 
                       np.linspace(start[1], end[1], 9))
            im0.set_clim(mi, ma)
            fig.colorbar(im0, ax=ax[0], fraction=0.046, pad=0.04)
            
            im1 = ax[1].imshow(outputmaps[i], cmap='jet')
            ax[1].set_title(f"corrected off resonance (Hz) z={zz[zi].round(2)}")
            ax[1].set_xlabel("X (cm)")
            ax[1].set_xticks(np.linspace(-.5,int(offres.shape[0])-.5,9),
                             np.linspace(start[0], end[0], 9))
            ax[1].set_ylabel("Y (cm)")
            ax[1].set_yticks(np.linspace(-.5,int(offres.shape[1])-.5,9), 
                             np.linspace(start[1], end[1], 9))
            im1.set_clim(mi, ma)
            fig.colorbar(im1, ax=ax[1], fraction=0.046, pad=0.04)

            plt.savefig(os.path.join(outputdir, f"{zi}corr_z_{zz[zi].round(2)}.png"))
            plt.close()
 
            #plt.imshow(outputmaps[i], cmap='jet')
            #plt.title(f"corrected fieldmap at z = {zz[zi]}")
            #plt.colorbar(fraction=0.046, pad=0.04)
            #plt.clim(mi, ma)
            #plt.savefig(os.path.join(b0offres, f"{zi}.png"))
            #plt.close()
        

