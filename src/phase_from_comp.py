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
    parser.add_argument("--pic", action="store_true", help="")
    parser.add_argument("--save", action="store_true", help="")
    parser.add_argument("--nifti", action="store_true", help="")

    args = parser.parse_args()

    folder = args.folder
    raw_dir = os.path.join(folder, "raw")
    proc_fold = "processed"
    proc_dir = os.path.join(folder, "processed")
    if args.save and not os.path.exists(proc_dir):
        os.makedirs(proc_dir)

    nifti_fold = "niftis"
    nifti_dir = os.path.join(folder, "niftis")
    if args.nifti and not os.path.exists(nifti_dir):
        os.makedirs(nifti_dir)

    pattern = r'\.(\d+)'
    
    for scan_folder in os.listdir(raw_dir):
        print("iterating on scan {scan_folder}")

        scan_dir = os.path.join(raw_dir, scan_folder)
        save_dir = os.path.join(proc_dir, scan_folder)
        nifti_save_dir = os.path.join(nifti_dir, scan_folder)

        # separate the real and imaginary images from the dcm folder that they
        # interleaved in

        images = {}
        names = {}
        
        for file_name in os.listdir(scan_dir):
            file_path = os.path.join(scan_dir, file_name)
            ds = dcmread(file_path)
            file_index = int(re.findall(pattern, file_name)[0])
            names[file_index] = file_name
            images[file_index] = ds
            #print(ds.pixel_array.dtype)


        sorted_items = sorted(images.items())
        images_list = [img.pixel_array for idx, img in sorted_items]

        # from looking at all of them, it looks like real and imaginary images are
        # every 3rd and 4th image in order. 
        #
        # so if there are 16 slices, then the order from all of them will be 
        # magnitude slice 1, phase slice 1, real slice 1, imag slice 1, mag slice 2
        # ... etc

        # i think it goes real = 3rd and imag = 4th and not the other way around

        real = np.array(images_list[2::4])
        imag = np.array(images_list[3::4])
        og_phases = np.array(images_list[1::4])

        og_phases = np.clip(og_phases, -4096, 4095)
        og_phases = og_phases / np.max(np.abs(og_phases))

        phase = np.arctan2(imag, real)
        phase = phase / np.pi
        phase = (phase * 4095).astype(np.int16)
        print(f"max {np.max(phase)} min {np.min(phase)}")
                #arr = arr / np.max(np.abs(arr))


        if (args.save):
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            for i in range(0,16):
                #save the new phases
                images[i*4+2].PixelData = phase[i].tobytes()
                images[i*4+2].save_as(os.path.join(save_dir, names[i*4+1]))
                # save the magnitudes
                images[i*4+1].save_as(os.path.join(save_dir, names[i*4]))

        if (args.nifti):
            if not os.path.exists(nifti_save_dir):
                os.makedirs(nifti_save_dir)
            dicom2nifti.convert_directory(save_dir, nifti_save_dir, reorient=True)



        if (args.pic):
            for i in range(0,16):
                fig, ax = plt.subplots(2, 2, figsize=(10, 10))  # 1 row, 2 columns

                im1 = ax[0,0].imshow(phase[i], cmap='gray')
                ax[0,0].set_title('phase from arctan(imag/real)')
                fig.colorbar(im1, ax=ax[0,0])  # Add a color bar for the first image

                arr = og_phases[i]
                im2 = ax[0,1].imshow(arr, cmap='gray')
                ax[0,1].set_title('original phase image')
                fig.colorbar(im2, ax=ax[0,1])  # Add a color bar for the first image

                arr = images_list[i*4+2]
                im4 = ax[1,0].imshow(arr, cmap='gray')
                ax[1,0].set_title('original "real" image')

                arr = images_list[i*4+3]
                im5 = ax[1,1].imshow(arr, cmap='gray')
                ax[1,1].set_title('original "imag" image')

                plt.show()



