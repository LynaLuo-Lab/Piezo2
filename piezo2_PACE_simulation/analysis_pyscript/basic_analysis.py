import os

import pandas as pd
import numpy as np
import math
from numpy.linalg import norm
import MDAnalysis as mda
from scipy.stats import pearsonr


def angle(x1, y1, x2, y2, x3, y3):
    # Calculate vectors BA and BC
    BA_x = x1 - x2
    BA_y = y1 - y2
    BC_x = x3 - x2
    BC_y = y3 - y2
    
    # Calculate dot product of BA and BC
    dot_product = BA_x * BC_x + BA_y * BC_y
    
    # Calculate magnitudes of vectors BA and BC
    magnitude_BA = math.sqrt(BA_x**2 + BA_y**2)
    magnitude_BC = math.sqrt(BC_x**2 + BC_y**2)
    
    # Calculate cosine of the angle
    cos_theta = dot_product / (magnitude_BA * magnitude_BC)
    
    # Calculate angle in radians
    theta_radians = math.acos(cos_theta)
    
    # Convert radians to degrees
    theta_degrees = math.degrees(theta_radians)
    
    return theta_degrees

def compute_dihedral(center1, center2, center3, center4):
    """
    Calculate the dihedral angle (in radians) defined by four centers of mass.

    Args:
    - center1, center2, center3, center4: Numpy arrays representing the positions
      of the four centers of mass in 3D space.

    Returns:
    - dihedral_angle: The computed dihedral angle in radians.
    """

    # Calculate the vectors between the centers of mass
    vector1 = center2 - center1
    vector2 = center3 - center2
    vector3 = center4 - center3

    # Calculate the normal vectors to the planes defined by the vectors
    normal1 = np.cross(vector1, vector2)
    normal2 = np.cross(vector2, vector3)

    # Calculate the angle between the normals
    dot_product = np.dot(normal1, normal2)
    cross_product = np.cross(normal1, normal2)
    dihedral_angle = np.arctan2(np.linalg.norm(cross_product), dot_product)

    return dihedral_angle


# use domz1 and distI to compute radius of curvature in nm
def domerad(rad, h):
    R = (rad*rad + h*h)/(2*h)
    return R

def projA(rad):
    area=np.pi*rad*rad
    return area
        
# inradius of a triangle
def rad(a,b,c):
    return ((a*b*c)/np.sqrt((a+b+c)*(a+b-c)*(a-b+c)*(b+c-a))) 

def centroid(a, b, c):
    # Given triangle side lengths a, b, c
    # Calculate semi-perimeter
    s = (a + b + c) / 2

    # Calculate centroid coordinates
    centroid_x = (a**2 + b**2 - c**2) / (4 * math.sqrt(a**2 - (b - c)**2))
    centroid_y = math.sqrt(a**2 - centroid_x**2)

    # Calculate distance from centroid to any vertex (let's use the first vertex as an example)
    distance_to_centroid = math.sqrt(centroid_x**2 + centroid_y**2)

    return distance_to_centroid



######################
def flatten_angle(u):
    """beam and Pore helix angle flatting"""
    A1 = u.select_atoms("segid A B and name CA").center_of_geometry() # arm1
    A2 = u.select_atoms("segid K L and name CA").center_of_geometry() # arm1
    A3 = u.select_atoms("segid U V and name CA").center_of_geometry() # arm1
    B = u.select_atoms("resid 505-526 and segid J T d and name CA").center_of_geometry() # CTD
    C = u.select_atoms("resid 395-405 and segid J T d and name CA").center_of_geometry() # cap
    BA1 = B - A1
    BA2 = B - A2
    BA3 = B - A3
    CB =  B - C
    theta1 = np.rad2deg(np.arccos(np.dot(BA1, CB)/(norm(BA1)*norm(CB))))
    theta2 = np.rad2deg(np.arccos(np.dot(BA2, CB)/(norm(BA2)*norm(CB))))
    theta3 = np.rad2deg(np.arccos(np.dot(BA3, CB)/(norm(BA3)*norm(CB))))
    return (theta1 + theta2 + theta3)/3

def flatten_angle_repeatF(u):
    """beam and Pore helix angle flatting"""
    A1 = u.select_atoms("segid F and name CA").center_of_geometry() # arm1
    A2 = u.select_atoms("segid P and name CA").center_of_geometry() # arm1
    A3 = u.select_atoms("segid Z and name CA").center_of_geometry() # arm1
    B = u.select_atoms("resid 505-526 and segid J T d and name CA").center_of_geometry() # CTD
    C = u.select_atoms("resid 395-405 and segid J T d and name CA").center_of_geometry() # cap
    BA1 = B - A1
    BA2 = B - A2
    BA3 = B - A3
    CB =  B - C
    theta1 = np.rad2deg(np.arccos(np.dot(BA1, CB)/(norm(BA1)*norm(CB))))
    theta2 = np.rad2deg(np.arccos(np.dot(BA2, CB)/(norm(BA2)*norm(CB))))
    theta3 = np.rad2deg(np.arccos(np.dot(BA3, CB)/(norm(BA3)*norm(CB))))
    return (theta1 + theta2 + theta3)/3

def hrot1(u):
    # Calculate the rotation of beam
    beam = u.select_atoms("segid J and (resid 235 to 260) and name CA").center_of_geometry()[:2]
    IH = u.select_atoms("segid J and resid 505 to 526 and name CA").center_of_geometry()[:2]
    com= u.select_atoms("segid J T d and resid 505 to 526 and name CA").center_of_geometry()[:2]
 
    combined= np.concatenate((beam, com, IH))
    return angle(*combined)

def hrot2(u):
    beam = u.select_atoms("segid d and (resid 235 to 260) and name CA").center_of_geometry()[:2]
    IH = u.select_atoms("segid d and resid 505 to 526 and name CA").center_of_geometry()[:2]
    com= u.select_atoms("segid J T d and resid 505 to 526 and name CA").center_of_geometry()[:2]
 
    combined= np.concatenate((beam, com, IH))
    return angle(*combined)

def hrot3(u):
    beam = u.select_atoms("segid T and (resid 235 to 260) and name CA").center_of_geometry()[:2]
    IH = u.select_atoms("segid T and resid 505 to 526 and name CA").center_of_geometry()[:2]
    com= u.select_atoms("segid J T d and resid 505 to 526 and name CA").center_of_geometry()[:2]
 
    combined= np.concatenate((beam, com, IH))
    return angle(*combined)

def domz1(u):    
    """Repeat I height"""
    M = u.select_atoms("segid A  K  U and name CA").center_of_geometry()
    N = u.select_atoms("segid J T d and (resid 505 to 526) and name CA").center_of_geometry()
    NM = M[2] - N[2]
    return NM

def domz2(u):    
    """Repeat F height"""
    M = u.select_atoms("segid F P Z and name CA").center_of_geometry()
    N = u.select_atoms("segid J T d and (resid 505 to 526) and name CA").center_of_geometry()
    NM = M[2] - N[2]
    return NM

def distValPore(u):    
    """V2750 distance"""
    M = u.select_atoms("segid J and (resid 515) and name CB").center_of_geometry()
    N = u.select_atoms("segid T and (resid 515) and name CB").center_of_geometry()
    Z = u.select_atoms("segid d and (resid 515) and name CB").center_of_geometry()
    ZN = np.linalg.norm(Z - N)
    NM = np.linalg.norm(N - M)
    MZ = np.linalg.norm(M - Z)
    return (ZN + NM + MZ)/3

def distF1(u):    
    """repeat F distance"""
    M = u.select_atoms("segid F and (resid 101) and name CA").center_of_geometry()
    N = u.select_atoms("segid P and (resid 101) and name CA").center_of_geometry()
    Z = u.select_atoms("segid Z and (resid 101) and name CA").center_of_geometry()
    ZN = Z - N
    NM = N - M
    MZ = M - Z
    return np.linalg.norm(ZN)

def distF2(u):    
    M = u.select_atoms("segid F and (resid 101) and name CA").center_of_geometry()
    N = u.select_atoms("segid P and (resid 101) and name CA").center_of_geometry()
    Z = u.select_atoms("segid Z and (resid 101) and name CA").center_of_geometry()
    ZN = Z - N
    NM = N - M
    MZ = M - Z
    return np.linalg.norm(NM)

def distF3(u):    
    M = u.select_atoms("segid F and (resid 101) and name CA").center_of_geometry()
    N = u.select_atoms("segid P and (resid 101) and name CA").center_of_geometry()
    Z = u.select_atoms("segid Z and (resid 101) and name CA").center_of_geometry()
    ZN = Z - N
    NM = N - M
    MZ = M - Z
    return np.linalg.norm(MZ)


def dist67a(u):    
    """last repeat 67S distance"""
    M = u.select_atoms("segid A and resid 62 and name CA").center_of_geometry()
    N = u.select_atoms("segid K and resid 62 and name CA").center_of_geometry()
    Z = u.select_atoms("segid U and resid 62 and name CA").center_of_geometry()
    ZN = Z - N
    NM = N - M
    MZ = M - Z
    return np.linalg.norm(ZN)

def dist67b(u):    
    M = u.select_atoms("segid A and resid 62 and name CA").center_of_geometry()
    N = u.select_atoms("segid K and resid 62 and name CA").center_of_geometry()
    Z = u.select_atoms("segid U and resid 62 and name CA").center_of_geometry()
    ZN = Z - N
    NM = N - M
    MZ = M - Z
    return np.linalg.norm(NM)

def dist67c(u):    
    M = u.select_atoms("segid A and resid 62 and name CA").center_of_geometry()
    N = u.select_atoms("segid K and resid 62 and name CA").center_of_geometry()
    Z = u.select_atoms("segid U and resid 62 and name CA").center_of_geometry()
    ZN = Z - N
    NM = N - M
    MZ = M - Z
    return np.linalg.norm(MZ)


def distI1(u):    
    """last repeat I distance"""
    M = u.select_atoms("segid A B and name CA").center_of_geometry()
    N = u.select_atoms("segid K L and name CA").center_of_geometry()
    Z = u.select_atoms("segid U V and name CA").center_of_geometry()
    ZN = Z - N
    NM = N - M
    MZ = M - Z
    return np.linalg.norm(ZN)

def distI2(u):    
    M = u.select_atoms("segid A B and name CA").center_of_geometry()
    N = u.select_atoms("segid K L and name CA").center_of_geometry()
    Z = u.select_atoms("segid U V and name CA").center_of_geometry()
    ZN = Z - N
    NM = N - M
    MZ = M - Z
    return np.linalg.norm(NM)

def distI3(u):    
    M = u.select_atoms("segid A B and name CA").center_of_geometry()
    N = u.select_atoms("segid K L and name CA").center_of_geometry()
    Z = u.select_atoms("segid U V and name CA").center_of_geometry()
    ZN = Z - N
    NM = N - M
    MZ = M - Z
    return np.linalg.norm(MZ)


import sys

xtc_file = []
path = "/home/shuli/work/piezo/martini-pace/total_analysis/1_xtc/"
for rl in open(sys.argv[1],'r'):
    xtc_file.append(rl.split(" ")[0][:-1])


analysisR = []

for xtc in xtc_file:
  print(xtc)
  ru  = mda.Universe('./md_protein.pdb', path + xtc + ".xtc")
  num = 0
  data = []
  fw = open( xtc + ".txt", 'w')
  fw.write("#time, _distVP, _dist67, _distF,  _domz1, _flattenA, _projA, _domeA, _hrot, _projA_RF, _rad\n")
  '''
       #time      : simulation time, unit: ns \n 
       #_distVP   : pore size\n 
       #_dist67   : arm-arm distance S67\n 
       #_distF    : arm-arm distance K776\n 
       #_domz1    : dome height\n 
       #_flattenA : flatten angle\n 
       #_projA    : projected area\n 
       #_domeA    : dome area\n 
       #_hrot     : rotation_angle\n 
       #_projA_RF : projected area repeatF\n
       #_rad      : dome_radius\n
  '''
  
  for fr in ru.trajectory:
    # flatten angle
    _flattenA = flatten_angle(ru)
    
    _flattenA_rF = flatten_angle_repeatF(ru)
    
    # pore distance 
    _distVP = distValPore(ru)/10
    
    # rotation angle between tm37, center of tm38 and tm38 
    _hrot1 = hrot1(ru)
    _hrot2 = hrot2(ru)
    _hrot3 = hrot3(ru)
    _hrot  = (_hrot1+_hrot2+_hrot3) / 3
    
    # repeatF distance
    _distF1 = distF1(ru)
    _distF2 = distF2(ru)
    _distF3 = distF3(ru)
    _distF  = (_distF1 + _distF2 + _distF3) / 30
    
    _radF   = rad(_distF1, _distF2, _distF3) / 10
    _projA_rF = projA(_radF) 
    
    # repeatI distance
    _distI1 = distI1(ru)
    _distI2 = distI2(ru)
    _distI3 = distI3(ru)
    _distI  = (_distI1 + _distI2 + _distI3) / 30
    
    _rad   = rad(_distI1, _distI2, _distI3) / 10

    # project area
    _projA = projA(_rad)

    # dome Z distance between repeatI and tm38
    _domz1 = domz1(ru) / 10
    
    # domerad
    _domeR = domerad(_rad, _domz1)
    
    # dome area
    _domeA = 2*3.14*_domz1*_domeR

    # repeatI SER67 distance
    _dist67a = dist67a(ru)
    _dist67b = dist67b(ru)
    _dist67c = dist67c(ru)
    _dist67  = (_dist67a + _dist67b + _dist67c)/30
    
    fw.write(f"{num:6d} {_distVP:7.2f} {_dist67:7.2f} {_distF:7.2f} {_domz1:7.2f} {_flattenA:7.2f} {_projA:7.2f} {_domeA:7.2f} {_hrot:7.2f} {_projA_rF:7.2f} {_domeR:7.2f}\n")
    num += 1
  fw.close()

