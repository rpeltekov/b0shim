import argparse
import numpy as np
import sys, os
sys.path.append('../../tools/biot-savart-master')
import biot_savart_v4_3 as bs
from matplotlib import pyplot as plt

gyro = 4258

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
    ancors_basic = [[ 1, 0, 0], [-1, 0, 0], [ 0, 1, 0], [ 0,-1, 0], [ 0, 0, 1],
                    [ 0, 0,-1], [ 1, 1, 0], [ 1,-1, 0], [-1, 1, 0], [-1,-1, 0],
                    [ 1, 0, 1], [ 1, 0,-1], [-1, 0, 1], [-1, 0,-1], [ 0, 1, 1],
                    [ 0, 1,-1], [ 0,-1, 1], [ 0,-1,-1], [ 1, 1, 1], [-1, 1, 1],
                    [ 1,-1, 1], [-1,-1, 1], [ 1, 1,-1], [-1, 1,-1], [ 1,-1,-1],
                    [-1,-1,-1]]
    return ancors_basic

def gen_coils(ancors, r, R, folder="", n=200):
    """
    create geometries of the coils with 1 amp set to a folder with the given
    folder
    """
    coils = []
    if not os.path.exists(folder):
        os.makedirs(folder)
    for i, ancor in enumerate(ancors):
        coil_folder_filepath = os.path.join(folder, f"coil_{i}")
        coil = coil_from_ancor(ancor, r, R, n);
        # add the current values of 1 to each of the segments
        coil = np.hstack((coil, np.ones((coil.shape[0], 1))))
        coils.append(coil)

        if not os.path.exists(coil_folder_filepath):
            os.makedirs(coil_folder_filepath)

        coil_filename = os.path.join(coil_folder_filepath, 'coil.txt')
        np.savetxt(coil_filename, coil, delimiter=", ")
    
    return coils

def compute_bzfields(coils, volume, folder, vis=False, debug=False):
    """
        generate the delt bz magnetization that the coil generates at pos z

        if z is None, it will return the magnetization for the whole volume.
    """
    np.savetxt(os.path.join(folder, "volume_def.txt"), volume)

    start = volume[0]
    size = volume[1]
    end = start + size
    res = volume[2]
    x = np.linspace(start[0], end[0], int(res[0]+1))
    y = np.linspace(start[1], end[1], int(res[1]+1))
    z = np.linspace(start[2], end[2], int(res[2]+1))

    # now we need to place the sample points in the middle of each voxel
    x = x[:-1]
    y = y[:-1]
    z = z[:-1]
    x = x + (x[1] - x[0])/2
    y = y + (y[1] - y[0])/2
    z = z + (z[1] - z[0])/2

    Z, Y, X = np.meshgrid(z, y, x, indexing='ij')

    fields = []
    for i, coil in enumerate(coils):
        coilname = f"coil_{i}"
        field_folder_filepath = os.path.join(folder, coilname)

        if vis:
            plot_coils([coil[:,:-1]]) # cut off the current value for vis
            
        if i >4:
            field = bs.calculate_field(coil, X,Y,Z) * gyro
            fields.append(field)

            plot_save_field_slices(field, x, y, z, field_folder_filepath, coilname, True)
        

        #np.savetxt(folder+"/"+volume_name+".txt", field, delimeter=", ")

def plot_save_field_slices(field, x, y, z, folder, coilname, viz=False):
    for zlevel in z:
        level = np.where(z == zlevel)
        bslice = field[:,:,level[0][0],2]
        x_label,y_label = "x","y"
        x_array,y_array = x,y
        plt.xlabel('x (cm)')
        plt.xlabel('y (cm)')
        plt.xticks(x)
        plt.yticks(y)
        plt.title(f'Bz Field (Hz) of {coilname} at Z={zlevel} cm')
        plt.imshow(bslice, cmap='RdYlGn')
        plt.colorbar()
        filename = os.path.join(folder, coilname+f"_{zlevel}.png")
        if viz:
            plt.show()
        #plt.savefig(filename)

def plot_coils(coils, volume=None):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    if volume is not None:
        start = volume[0]
        end = volume[0] + volume[1]
        ax.set_xlim(start[0], end[0])
        ax.set_ylim(start[1], end[1])
        ax.set_zlim(start[2], end[2])

    for coil in coils:
        ax.plot(coil[:,0], coil[:,1], coil[:,2], label='parametric curve')
    plt.show()

def plot_ancors_on_circle(ancors, r, R, volume=None):

    # add the sphere to the drawing
    # TODO: currently covers the actual coils not nicely
    # u = np.linspace(0, 2 * np.pi, 100)
    # v = np.linspace(0, np.pi, 100)
    # x = 1.9*R * np.outer(np.cos(u), np.sin(v))
    # y = 1.9*R * np.outer(np.sin(u), np.sin(v))
    # z = 1.9*R * np.outer(np.ones(np.size(u)), np.cos(v))
    #ax.plot_surface(x, y, z)

    coils = []
    for ancor in ancors:
        coils.append(coil_from_ancor(ancor,r,R,100))
        
    plot_coils(coils, volume)


def injest_coils(filename):
  pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--folder", default="default", help="")
    parser.add_argument("--resX", default=256, 
                        help="NUMBER OF VOXELS NOT VOXEL SIZE")
    parser.add_argument("--resY", default=256,
                        help="NUMBER OF VOXELS NOT VOXEL SIZE")
    parser.add_argument("--resZ", default=16,
                        help="NUMBER OF VOXELS NOT VOXEL SIZE")
    parser.add_argument("--fovX", default=16, help="")
    parser.add_argument("--fovY", default=16, help="")
    parser.add_argument("--fovZ", default=16, help="")
    parser.add_argument("--vis", action="store_true", help="")
    parser.add_argument("--debug", action="store_true", help="")
    parser.add_argument("--nosave", action="store_false", help="")

    args = parser.parse_args()

    res = np.array([args.resX, args.resY, args.resZ])
    fov = np.array([args.fovX, args.fovY, args.fovZ])
    start = -fov/2
    volume_def = np.vstack((start, fov, res))
    
    folder = args.folder

    if folder == "default":
        print("[INFO] doing default coil calculations: octogon")
        #TODO generate the mask volume, which will be used to mask effects?
        # is this necessary tho, because if we create these basis maps, then
        # adding them to some existing inhomogeneity will simply mask it there
        # when we compute the ridge regression problem... i guess we do not need
        # to do this...
        ancors = gen_loop_ancors()
        r = 3.5
        R = 10

    coil_folder = os.path.join(folder, "coils")

    if args.vis:
        plot_ancors_on_circle(ancors, r, R, volume_def)   

    coils = gen_coils(ancors, r, R, coil_folder)
    
    compute_bzfields(coils, volume_def, folder, True, args.vis)


