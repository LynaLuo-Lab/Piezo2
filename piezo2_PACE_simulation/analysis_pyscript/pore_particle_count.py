
import pandas as pd
import numpy as np

import numpy as np
from numpy.linalg import norm
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import sys


def point_in_triangle(x, y, a1, b1, a2, b2, a3, b3):
    def area(x1, y1, x2, y2, x3, y3):
        return abs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2))
    area_ABC = area(a1, b1, a2, b2, a3, b3)
    area_PAB = area(x, y, a1, b1, a2, b2)
    area_PBC = area(x, y, a2, b2, a3, b3)
    area_PCA = area(x, y, a3, b3, a1, b1)
    return area_ABC == area_PAB + area_PBC + area_PCA

def point_in_circle(x, y, a1, b1, r):
    distance_squared = (x - a1) ** 2 + (y - b1) ** 2
    if distance_squared <= r ** 2:
        return True
    else:
        return False

def point_in_cylinder(point_pos, center_pos, r, h):
    z_min = center_pos[2] - h
    z_max = center_pos[2] + h
    dist_xy = np.linalg.norm(point_pos[:2] - center_pos[:2])
    if z_min <= point_pos[2] <= z_max and dist_xy <= r:
        return True
    else:
        return False


def point_in_UL(z, c1):
    z_posit = ""
    if z > c1:
        z_posit = "upper"
    elif z < c1:
        z_posit = "lower"
    elif z == c1:
        z_posit = "middle"
    return z_posit


u = mda.Universe('./md.pdb',sys.argv[1])

waters = u.select_atoms('resname W') 
channel = u.select_atoms('protein and segid J T d and (resid 515) and name CA')
lipid_heads = u.select_atoms('resname POPC and name PO4 NC3  GL1 GL2') 
lipid_tails = u.select_atoms('resname POPC and name C1A C2A C3A C4A C1B C2B C3B C4B')
threshold = 30.0  


water_index_start = 191182   #  POPC 11355   protein residue number 5859
water_index_resid_start = 17214  # POPC 11355   protein residue number 5859

waters_pass_last_frame = np.zeros(896114, dtype=int)

for ts in u.trajectory[::2]:
    channel_center = channel.center_of_geometry()
    channel_coordinate = channel.positions
    vx1, vy1, vz1 = channel_coordinate[0][0],channel_coordinate[0][1],channel_coordinate[0][2]
    vx2, vy2, vz2 = channel_coordinate[1][0],channel_coordinate[1][1],channel_coordinate[1][2]
    vx3, vy3, vz3 = channel_coordinate[2][0],channel_coordinate[2][1],channel_coordinate[2][2]
    vcx, vcy, vcz = channel_center[0], channel_center[1], channel_center[2]


    M = u.select_atoms("segid J and (resid 515) and name CB").center_of_geometry()
    N = u.select_atoms("segid T and (resid 515) and name CB").center_of_geometry()
    Z = u.select_atoms("segid d and (resid 515) and name CB").center_of_geometry()
    ZN = np.linalg.norm(Z - N)
    NM = np.linalg.norm(N - M)
    MZ = np.linalg.norm(M - Z)
    pore_dist = (ZN + NM + MZ)/3

    val_distance_to_center = distances.distance_array(channel_coordinate, channel_center[None, :])
    distances_to_center = distances.distance_array(waters.positions, channel_center[None, :])
    waters_in_channel = waters[distances_to_center[:, 0] < threshold]
    waters_in_channel_index = waters_in_channel.indices - water_index_start
    coordinates = waters_in_channel.positions
   
    # cylinder parameters
    r = np.mean(val_distance_to_center)           # radius
    h = 15                                        # height

    heads_in_cylinder = []
    tails_in_cylinder = []
    for atom in lipid_heads:
        pos = atom.position
        lipid_point = point_in_cylinder(pos, channel_center, r, h)
        if lipid_point:
            heads_in_cylinder.append(atom.id)
            # if atom.resid not in heads_in_cylinder:
                # heads_in_cylinder.append(atom.resid)

    for atom in lipid_tails:
        pos = atom.position
        lipid_point = point_in_cylinder(pos, channel_center, r, h)
        if lipid_point:
            tails_in_cylinder.append(atom.id) 
            # if atom.resid not in tails_in_cylinder:
                # tails_in_cylinder.append(atom.resid) 
        
    num = 0
    cylinder_index = []
    upper_index = []
    lower_index = []
    middle_index = []
    for line in coordinates:
        x, y, z = line[0], line[1], line[2]
  
        water_in_cylinder = point_in_cylinder(line, channel_center, r, h )
        water_in_UL = point_in_UL(z, vcz)
        if water_in_cylinder:
           cylinder_index.append(waters_in_channel_index[num])

    print("%5s %4s %4s %4s %4s %6.2f\n"%(ts.frame/2, len(heads_in_cylinder), len(tails_in_cylinder), len(cylinder_index), pore_dist))
    #print "%5s %4s %4s %4s %6.2f"%(ts.frame/2, len(heads_in_cylinder), len(tails_in_cylinder), len(cylinder_index), pore_dist)
