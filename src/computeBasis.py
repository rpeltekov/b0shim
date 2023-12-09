import argparse
import threading
import numpy as np
import sys, os
sys.path.append('../tools/biot-savart-master')
import biot_savart_v4_3 as bs
from matplotlib import pyplot as plt
from matplotlib.colors import TwoSlopeNorm as diverge
from apng import APNG


gyro = 4258

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


def coil_from_ancor(d, r, R, n=100):
    """ 
    this function will create a coil geometry to be placed on a sphere of radius
    r with provided normal ancor from the center of the coil. 
    it will spit out n wire segments that can be pasted into a file and used on
    the biot savart library.

    inputs:
    d:       3x1 [float] the normal ancor of the coil plane from sphere center
    r:       [float] the radius of the coil
    R:       [float] the radius of the sphere
    n:       [int] the number of segments to make for the coil
    """

    d = np.array(d)
    d = d / np.linalg.norm(d, ord=2) # normalize the dection ancor

    # parametric equation of 2d circle in R^3
    # x0 +r(cost)v1+r(sint)v2
    #   x0 is the center of the circle, v1 and v2 are the basis ancors

    # get the center of the circle, that sits in the sphere, because we have to
    # drop the coil onto the surface of the sphere and the sphere is convex
    x0 = np.sqrt(R**2 - r**2) * d

    # get a single ancor that is orthogonal to x0 and normalize

    # need x0 dot v1 = 0
    # x0[0]*v1[0] + x0[1]*v1[1] + x0[2]*v1[2] = 0
    # set two components of v1 to 1, and then solve for the 4th

    # i.e. 
    # v1[0] = 1
    # v1[1] = 1
    # v1[2] = (-x0[0]-x0[1])/x0[2]

    # if x0[2] is 0 then use another one

    if (x0[2] != 0):
        v1 = np.array([1, 1, (-x0[0]-x0[1])/x0[2]])
    elif (x0[1] != 0):
        v1 = np.array([1, (-x0[0]-x0[2])/x0[1], 1])
    else:
        v1 = np.array([(-x0[1]-x0[2])/x0[0], 1, 1])
    v1 = v1 / np.linalg.norm(v1, ord=2)

    # get a third ancor that is orthogonal to x0 and v1 using crossproduct 
    # and normalize
    v2 = np.cross(x0, v1)
    v2 = v2 / np.linalg.norm(v2, ord=2)

    t = np.linspace(0, 2*np.pi, n)

    x0s = np.tile(x0, (n,1))
    points = x0s + np.cos(t).reshape((n,1))@(r*v1).reshape((1,3)) \
               + np.sin(t).reshape((n,1))@(r*v2).reshape((1,3))
    return points

def gen_loop_ancors(n=0):
    """
      generate the positions for the coils. probably use some fibonacci 
      spiral sequence or something
      
      inputs:
       n: [int] the number of points to pattern onto the surface
       n=0 -> means that we generate the basic octohedron
    """
    ancors = []
    if n==0:
        ancors = [[ 1, 0, 0], [-1, 0, 0], [ 0, 1, 0], [ 0,-1, 0], [ 0, 0, 1],
                  [ 0, 0,-1], [ 1, 1, 0], [ 1,-1, 0], [-1, 1, 0], [-1,-1, 0],
                  [ 1, 0, 1], [ 1, 0,-1], [-1, 0, 1], [-1, 0,-1], [ 0, 1, 1],
                  [ 0, 1,-1], [ 0,-1, 1], [ 0,-1,-1], [ 1, 1, 1], [-1, 1, 1],
                  [ 1,-1, 1], [-1,-1, 1], [ 1, 1,-1], [-1, 1,-1], [ 1,-1,-1],
                  [-1,-1,-1]]

    # 7 aroung the middle, a row of 7 above that, and then 1 at the top
    pistep = 2*np.pi/7
    pis = np.arange(0, 2*np.pi, pistep)

    if n>=15:

        bottom = np.vstack((np.cos(pis), np.sin(pis), np.zeros(7))).T
        # for now place the middle row at 45 degrees on top of the bottom row
        middle = np.vstack((np.cos(pis+pistep/2), np.sin(pis+pistep/2), np.ones(7)*.9)).T
        top = [0, 0, 1]
        
        ancors = np.vstack((bottom, middle, top))

    if n==1:
        ancors = [[0,0,1]]

    return ancors

def gen_coils(ancors, r, R, folder="", gradient=False, n=501):
    """
    create geometries of the coils with 1 amp set to a folder with the given
    folder
    """
    coils = []
    if not os.path.exists(folder):
        os.makedirs(folder)

    if gradient:
        for i in range(3):
            coil_folder_filepath = os.path.join(folder, f"coil_{i}")
            if not os.path.exists(coil_folder_filepath):
                os.makedirs(coil_folder_filepath)

    for i, ancor in enumerate(ancors):
        if gradient:
            coil_folder_filepath = os.path.join(folder, f"coil_{i+3}")
        else:
            coil_folder_filepath = os.path.join(folder, f"coil_{i}")

        coil = coil_from_ancor(ancor, r, R, n);
        # add the current values of 1 to each of the segments
        coil = np.hstack((coil, np.ones((coil.shape[0], 1))))
        coils.append(coil)

        if not os.path.exists(coil_folder_filepath):
            os.makedirs(coil_folder_filepath)

        coil_filename = os.path.join(coil_folder_filepath, f'coil_{i}.txt')
        np.savetxt(coil_filename, coil, delimiter=", ")
    
    return coils

def compute_bzfields(coils, volume, folder, vis=False, gradient=False, debug=False):
    """
        generate the delt bz magnetization that the coil generates at pos z

        if z is None, it will return the magnetization for the whole volume.
    """
    x, y, z = get_ticks(volume)

    X, Y, Z = np.meshgrid(x, y, z)
    Z = Z.flatten()
    Y = Y.flatten()
    X = X.flatten()

    fields = []
    if gradient:
        g = np.array([[1, 0, 0],[0, 1, 0], [0, 0, 1]]) # gaus per cm
        together = np.vstack((X, Y, Z))
        for i in range(3):
            coilname = f"coil_{i}"
            field_folder_filepath = os.path.join(folder, coilname)
            field = g[i].dot(together) * gyro
            plot_save_field_slices(np.vstack((np.zeros(field.shape),
                                             np.zeros(field.shape), 
                                             field)).T, 
                                   x, y, z, Z, field_folder_filepath,
                                   coilname, vis)
            np.save(os.path.join(field_folder_filepath, "field.npy"), field)

    def compute_coil_field(i, coil):
        if gradient:
            i = i+3
        coilname = f"coil_{i}"
        print(f"[INFO]  Computing field for {coilname}")
        field_folder_filepath = os.path.join(folder, coilname)

        centerpoint = bs.calculate_field(coil.T, 0, 0, 5)
        print("[INFO]: the center of the coil sees this magnetic field gaus:", centerpoint)
        mag = np.sqrt(np.sum(centerpoint**2))
        print("[INFO]: magnitude", mag)

        # if i am iterating the plots, and data is available, reuse it
        if os.path.exists(os.path.join(field_folder_filepath, "field.npy")):
            field = np.load(os.path.join(field_folder_filepath, "field.npy"))
        else:
            field = bs.calculate_field(coil.T, X,Y,Z)
            field = field * gyro

            field = field[:,2]
            np.save(os.path.join(field_folder_filepath, "field.npy"), field)

        plot_save_field_slices(field, x, y, z, X, Y, Z, 
                               field_folder_filepath, coilname, vis)

    if len(coils) < 3:
        for i, coil in enumerate(coils):
            compute_coil_field(i, coil)
    else:
        threads = []
        for i, coil in enumerate(coils):
            threads.append(threading.Thread(target=compute_coil_field, 
                                            args=(i,coil)))
            threads[i].start()
        for thread in threads:
            thread.join()
        

def plot_save_field_slices(field, x, y, z, X,Y,Z, folder, coilname, vis=False):
    zfolder = os.path.join(folder, "z")
    if not os.path.exists(zfolder):
        os.makedirs(zfolder)
    radius = np.sqrt(X**2 + Y**2 + Z**2)
    r_outside_i = np.argwhere(radius >= 8)
    field[r_outside_i] = 0
    alpha = np.ones(field.shape)
    alpha[r_outside_i] = 0

    print('[INFO] Saving axial slices of coil to png and apng')
    print(f'[INFO] Max value: {np.max(np.abs(field))}')
    apng = APNG()
    for i, zlevel in enumerate(z):
        level = np.where(Z == zlevel)
        bslice = field[level].reshape(len(x), len(y))
        aslice = alpha[level].reshape(len(x), len(y))
        bslice = -np.flip(bslice, axis=0)
        if not aslice.any():
            continue
        x_label,y_label = "x","y"
        x_array,y_array = x,y
        plt.xlabel('x (cm)')
        plt.ylabel('y (cm)')
        plt.xticks(np.linspace(0,len(x)-1,5), np.linspace(x[-1], x[0], 5))
        plt.yticks(np.linspace(0,len(y)-1,11), np.linspace(y[-1], y[0], 11))
        plt.title(f'Bz Field (Hz) of {coilname} at Z={zlevel} cm')
        plt.imshow(bslice, alpha=aslice, cmap='jet')
        plt.clim(-300,300)
        plt.colorbar()
        filename = os.path.join(zfolder,f"{i}_{round(zlevel, 2)}.png")
        if vis:
            plt.show()
        plt.savefig(filename, transparent=True)
        apng.append_file(filename)
        plt.close()
    apng.save(os.path.join(zfolder,f'anim.apng'))


def plot_coils(coils, volume=None, field=None, pos=None, folder=None,
               name=None, vis=False):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    if volume is not None:
        maxdim = np.max(volume[1])*np.ones(3)
        start = -maxdim/2
        end = start + maxdim
        ax.set_xlim(start[0], end[0])
        ax.set_ylim(start[1], end[1])
        ax.set_zlim(start[2], end[2])

    if field is not None:
        print("pos shape,", pos.shape)
        print("field shape,", field.shape)
        ax.plot(coils[:,0], coils[:,1], coils[:,2], label='parametric curve')
        ax.quiver(pos[0], pos[1], pos[2], field[:,0], field[:,1],
                  field[:,2], length=1, normalize=True, color='b')
        plt.show()
        plt.close()
        return
    
    if name is not None:
        ax.plot(coils[:,0], coils[:,1], coils[:,2], label='parametric curve')
        plt.savefig(os.path.join(folder, name+"_isolated.png"))
        plt.close()
        return

    for coil in coils:
        ax.plot(coil[:,0], coil[:,1], coil[:,2], label='parametric curve')


    if vis:
        plt.show()
    plt.savefig(os.path.join(folder, "headcap.png"))
    plt.close()

def plot_ancors(ancors, r, R, volume, folder, vis):

    # add the sphere to the drawing
    # u = np.linspace(0, 2 * np.pi, 100)
    # v = np.linspace(0, np.pi, 100)
    # x = 1.9*R * np.outer(np.cos(u), np.sin(v))
    # y = 1.9*R * np.outer(np.sin(u), np.sin(v))
    # z = 1.9*R * np.outer(np.ones(np.size(u)), np.cos(v))
    #ax.plot_surface(x, y, z)

    coils = []
    for ancor in ancors:
        coils.append(coil_from_ancor(ancor,r,R,100))
        
    plot_coils(coils, volume=volume, folder=folder, vis=vis)

def injest_coils(filename):
  pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--folder", default="default", help="")
    parser.add_argument("--resX", type=int, default=256, 
                        help="NUMBER OF VOXELS NOT VOXEL SIZE")
    parser.add_argument("--resY", type=int, default=256,
                        help="NUMBER OF VOXELS NOT VOXEL SIZE")
    parser.add_argument("--resZ", type=int, default=16,
                        help="NUMBER OF VOXELS NOT VOXEL SIZE")
    parser.add_argument("--fovX", type=float, default=20, help="")
    parser.add_argument("--fovY", type=float, default=20, help="")
    parser.add_argument("--fovZ", type=float, default=20, help="")
    parser.add_argument("--numcoils", type=int, default=0, help="")
    parser.add_argument("--R", type=float, default=10., help="")
    parser.add_argument("--r", type=float, default=5., help="")
    parser.add_argument("--gradient", action="store_true", help="")
    parser.add_argument("--vis", action="store_true", help="")
    parser.add_argument("--nosave", action="store_false", help="")

    args = parser.parse_args()

    folder = os.path.join("../data/headcaps", args.folder)
    if not os.path.exists(folder):
        os.makedirs(folder)

    r = args.r
    R = args.R

    res = np.array([args.resX, args.resY, args.resZ])
    fov = np.array([args.fovX, args.fovY, args.fovZ])
    start = -fov/2
    volume_def = np.vstack((start, fov, res))

    print(f"[INFO] Headcap Setup:\n   VOLUME={volume_def}\n   R={R}, r={r}") 
    np.savetxt(os.path.join(folder, "volume_def.txt"),
               np.vstack((volume_def, [R,r, args.numcoils + 3*args.gradient])))

    if folder == "default":
        print("[INFO] doing default coil calculations: octohedron")
    else:
        print(f"[INFO] calculating basis maps in {folder}")

    ancors = gen_loop_ancors(args.numcoils)
    plot_ancors(ancors, r, R, volume_def, folder, args.vis)   
    coils = gen_coils(ancors, r, R, folder, args.gradient, 201)
    compute_bzfields(coils, volume_def, folder, args.vis, args.gradient)


