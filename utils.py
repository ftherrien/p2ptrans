import numpy as np
import numpy.linalg as la
from config import *

def lcm(x, y):
   """This function takes two
   integers and returns the L.C.M."""

   lcm = (x*y)//gcd(x,y)
   return lcm

def gcd(x, y):
    """This function implements the Euclidian algorithm
    to find G.C.D. of two numbers"""
    while(y):
        x, y = y, x % y
    return x

def normal(A):
    return A/np.ones((3,1)).dot(la.norm(A, axis=0).reshape((1,np.shape(A)[1])))

def find_uvw(stretch_dir, basis = np.eye(3)):
    min_dist = np.ones(3)*1000
    min_uvw = np.zeros((3,np.shape(stretch_dir)[1]))
    stretch_dir = normal(stretch_dir)
    for l,vec in enumerate(stretch_dir.T):
        for i in np.arange(-10,10):
            for j in np.arange(-10,10):
                for k in np.arange(-10,10):
                    if [i,j,k] != [0,0,0]:
                        unit_vec = basis.dot(np.array([i,j,k]))
                        unit_vec = unit_vec/la.norm(unit_vec)
                        cur_dist = la.norm(unit_vec-vec)
                        if cur_dist < min_dist[l]:
                            min_dist[l] = cur_dist
                            min_uvw[:,l] = np.array([i,j,k], np.int).T

        gcd3 = gcd(min_uvw[0,l],gcd(min_uvw[1,l], min_uvw[2,l]))
        min_uvw[:,l] = min_uvw[:,l]/gcd3

    return min_uvw


