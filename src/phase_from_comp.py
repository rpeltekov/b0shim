from pydicom import dcmread # to read the dicom images
import dicom2nifti      # to write to nifti images
import os               # to iterate through a folder
import re               # to extract the number / order from a dcm folder
import numpy as np      # to modify the images
import matplotlib.pyplot as plt # to visualize the images
import argparse
import sys

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Converting dicoms and "
                                                 +"getting their new phase")
    parser.add_argument("folder", type=str,default=0, help="")
    parser.add_argument("echos", type=int,default=1, help="")
    parser.add_argument("--pic", action="store_true", help="")
    parser.add_argument("--save", action="store_true", help="")
    parser.add_argument("--nifti", action="store_true", help="")
    parser.add_argument("--debug", action="store_true", help="")

    args = parser.parse_args()

    folder = args.folder
    proc_fold = "processed"
    proc_dir = os.path.join(folder, "processed")
    if args.save and not os.path.exists(proc_dir):
        os.makedirs(proc_dir)

    pattern = r'\.(\d+)'
    
    for scan_folder in os.listdir(folder):
        if scan_folder in ("processed", "b0maps"):
            continue
        print(f"[INFO] Iterating on scan {scan_folder}")

        scan_dir = os.path.join(folder, scan_folder)
        save_dir = os.path.join(proc_dir, scan_folder)

        if not os.path.exists(save_dir):
            print("[INFO] saving to:", save_dir)
            os.makedirs(save_dir)

        if len(os.listdir(scan_dir)) == 0:
            print(scan_dir, "exiting")

        # separate the real and imaginary images from the dcm folder that they
        # interleaved in

        tes = []
        images = {}
        names = {}
        
        for file_name in os.listdir(scan_dir):
            file_path = os.path.join(scan_dir, file_name)
            ds = dcmread(file_path)
            file_index = int(re.findall(pattern, file_name)[0])
            names[file_index] = file_name
            images[file_index] = ds
            te = float(ds[0x0018, 0x0081].value)
            if te not in tes:
                tes.append(te)

        sorted_items = sorted(images.items())
        images_list = np.array([img.pixel_array for idx, img in sorted_items])

        tes.sort()

        reals = np.array(images_list[2::4])
        imags = np.array(images_list[3::4])
        mags = np.array(images_list[::4])

        if (args.echos == 2):
            real = np.moveaxis([reals[::2], reals[1::2]], 0, -1)
            imag = np.moveaxis([imags[::2], imags[1::2]], 0, -1)
            mags = np.moveaxis([mags[::2], mags[1::2]], 0, -1)

            phases = np.arctan2(imag, real)

            if (args.debug):
                print(f"[DEBUG] Max {np.max(phase)} min {np.min(phase)}")

            testring = "".join(["te_"]+[f"{te}_" for te in tes])
            name = testring + "3dPhase.npy"
            phases = np.moveaxis(phases, 0, 2)
            np.save(os.path.join(save_dir, name), phases)

            namemag = testring + "3dMag.npy"
            mags = np.moveaxis(mags, 0, 2)
            np.save(os.path.join(save_dir, namemag), mags)

        else:
            phases = np.arctan2(imags, reals)
            name = f"te_{tes[0]}"+"_3dPhase.npy"
            phases = np.moveaxis(phases, 0, 2)
            np.save(os.path.join(save_dir, name), phases)


print("[INFO] Done.")
