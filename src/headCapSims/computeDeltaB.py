import numpy as np
import sys, os
sys.path.append('biot-savart-master')
import biot_savart_v4_3 as bs
from matplotlib import pyplot as plt

def coil_from_vector(d, r, R, n=100):
  """ 
  this function will create a coil geometry to be placed on a sphere of radius
  r with provided normal vector from the center of the coil. 
  it will spit out n wire segments that can be pasted into a file and used on
  the biot savart library.
  
  inputs:
   d:       3x1 [float] the normal vector of the coil plane from the circle center
   r:       [float] the radius of the coil
   R:       [float] the radius of the sphere
   n:       [int] the number of segments to make for the coil
  """
  
  d = np.array(d)
  d = d / np.linalg.norm(d, ord=2) # normalize the dection vector

  # parametric equation of 2d circle in R^3
  # x0 +r(cost)v1+r(sint)v2
  #   x0 is the center of the circle, v1 and v2 are the basis vectors

  # get the center of the circle, that sits in the sphere, because we have to
  # drop the coil onto the surface of the sphere and the sphere is convex
  x0 = np.sqrt(R**2 - r**2) * d
  
  # get a single vector that is orthogonal to x0 and normalize
  v1 = np.array([1, 1, (-x0[0]-x0[1])/x0[2]])
  v1 = v1 / np.linalg.norm(v1, ord=2)

  # get a third vector that is orthogonal to x0 and v1 using crossproduct and normalize
  v2 = np.cross(x0, v1)
  v2 = v2 / np.linalg.norm(v2, ord=2)

  t = np.linspace(0, 2*np.pi, n)
  
  x0s = np.tile(x0, (n,1))
  points = x0s + np.cos(t).reshape((n,1))@(r*v1).reshape((1,3)) \
               + np.sin(t).reshape((n,1))@(r*v2).reshape((1,3))
  return points

def gen_loop_ancors(n, R):
  """
    generate the positions for the coils. probably use some fibonacci 
    spiral sequence or something
    
    inputs:
     n: [int] the number of points to pattern onto the surface
     R: [float] the size of the sphere to pattern them onto.
  """
  pass

def save_generated_coils(ancors, r, R, folder):
  """
    save the geometries of the coils with 1 amp set to a folder with the given
    folder
  """
  if not os.path.exists(folder):
    os.makedirs(folder)

  for i, ancor in enumerate(ancors):
    coil = coil_from_vector(ancor, r, R, 200);
    # add the current values of 1 to each of the segments
    coil = np.hstack((coil, np.ones((coil.shape[0], 1))))

    #append(folder+'/coil_'+str(i)+'.txt')
    np.savetxt(folder+'/coil_'+str(i)+'.txt', coil, delimiter=", ")
    

def compute_and_save_bzfields(folder, R, z: None):
  """
    generate the delt bz magnetization that the coil generates at pos z

    if z is None, it will return the magnetization for the whole volume.
  """
  fields = []
  for file_path in os.listdir(folder):
    coil = file_path[:-4] # remove the .txt
    coil_filepath = folder+"/"+file_path
    bs.plot_coil(coil_filepath)
    volume_name = coil+"vol"
    volume = (2.1*R, 2.1*R, 2.1*R)
    start = (-1.05*R, -1.05*R, -1.05*R)
    bs.write_target_volume(coil_filepath, volume_name, volume, start, .5, .5)
    field = bs.read_target_volume(volume_name)
    fields.append(field)
    bs.plot_fields(field, volume, start, .5)
    #np.savetxt(folder+"/"+volume_name+".txt", field, delimeter=", ")

  final_field = fields[0]
  for field in fields[1:]:
    final_field += field
  bs.plot_fields(final_field, volume, start, .5)

def test_plot_circle():
  r = 3.5
  R = 10
  fig = plt.figure()
  ax = fig.add_subplot(projection='3d')
  
  # # add the sphere to the drawing
  # u = np.linspace(0, 2 * np.pi, 100)
  # v = np.linspace(0, np.pi, 100)
  # x = 2*R * np.outer(np.cos(u), np.sin(v))
  # y = 2*R * np.outer(np.sin(u), np.sin(v))
  # z = 2*R * np.outer(np.ones(np.size(u)), np.cos(v))
  # ax.plot_surface(x, y, z)

  # add the coil
  vectors = np.array(
             [[0,110,-20],
             [-30,86,-80]])
  #             [35,85.5,-80],
  #             [0.9,-105,40.3],
  #             [-58.3,-88.8,16.9],
  #             [-32.2,-70.7,81.7],
  #             [33.9,-71.4,81.9],
  #             [57,-91.3,18.9],
  #             [-86.7,-36.1,-7],
  #             [-78,-43.4,51.7],
  #             [-53.7,-4.6,91.3],
  #             [0.4,-27.3,104.4],
  #             [53.2,-5.3,89.7],
  #             [73.4,-44,50.4],
  #             [86.8,-37,-6.3],
  #             [-87.7,-19.9,-70.3],
  #             [-88.7,27.8,-23.6],
  #             [-87.8,10.1,40.6],
  #             [-52.5,60.2,72.2],
  #             [0.9,31.9,101.2],
  #             [55.6,58.8,69.9],
  #             [86.7,12.1,36],
  #             [89.8,26.3,-25.5],
  #             [91.7,-23.7,-70.9],
  #             [-72,39.8,-86.4],
  #             [-60,80.5,-32.8],
  #             [-31.6,102.4,23.7],
  #             [0.9,84.6,71],
  #             [34.4,100,22.8],
  #             [76.4,63.6,17.6],
  #             [58.1,80.3,-35.7],
  #             [75.9,37.1,-86.4]])

  for i in range(len(vectors)):
    vector = vectors[i]
    coil = coil_from_vector(vector,r,R,100)
    ax.plot(coil[:,0], coil[:,1], coil[:,2], label='parametric curve')

  
  save_generated_coils(list(vectors), r, R, "try1")
  compute_and_save_bzfields("try1", R, 5)
  #plt.show()



def injest_coils(filename):
  pass

test_plot_circle()

