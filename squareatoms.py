 # -*- coding: utf-8 -*-
"""
ATOM SIMULATOR squareatoms
Created on Mon Nov 15 14:45:06 2021
@author: Asari
"""

# Use notebook for zoom/pan
# %matplotlib notebook
#Use inline for interactive feats
# %matplotlib inline
# from matplotlib.gridspec import GridSpec
# from ipywidgets import interactive
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# import pylab as pl
import numpy as np

from numpy import fft as npf
import math as m
from numpy import cos as cos
from numpy import sin as sin
from numpy import tan as tan
from numpy import arcsin as arcsin
from numpy import arccos as arccos
from numpy import arctan as arctan

from numpy import sqrt as sqrt
from numpy import log as log



def squareatoms(pix, L, a, theta, e11, e12, e22, honeycomb, origin_x = 0, origin_y = 0):
    """
    Simulates a square atomic lattice and takes the FFT of the lattice

    Creates a (pix x pix) image array of a square atomic lattice with periodicity 'a' (nm),
    rotated from the x-direction by an angle theta. The image has side length L (nm).
    If honeycomb == 1, it creates the same image but with honeycomb pattern. 
    For either choice, the image height is normalized to [0,1].

    Args:

    pix [int]: # of pixels, like 256
    L [float]: length of simulated image in nanometers
    a [float]: atomic lattice constant in nanometers
    theta [degrees]: angle to rotate the image
    e11, e12, e22 [float]: strain tensor elements
    honeycomb [0 or 1]: is asking whether you want the atoms as dots (0) or holes (1)
    origin is optional, you can give a lattice vector [o1,o2]
        where o1,o2 are pixel positions of your origin
    """
  
    # Reciprocal lattice vectors for square crystal WITH STRAIN
    k1 = (2*np.pi/a)*np.array([1 + e11, e12])
    k2 = (2*np.pi/a)*np.array([e12, 1 + e22]) 

    # Create rotation matrix to the lattice by theta (assuming theta is in degrees)
    theta_rad = np.deg2rad(theta) # convert theta to radians
    rotmat = np.array([[cos(theta_rad), -sin(theta_rad)], [sin(theta_rad), cos(theta_rad)]])
    
    # Rotate the reciprocal lattice vector by the rotation matrix
    k1 = np.matmul(k1,rotmat)
    k2 = np.matmul(k2,rotmat)
    
    # Create a meshgrid of 'pix' number of points, values from 0 --> L
    xx = np.linspace(0, L, pix)
    [X, Y] = np.meshgrid(xx,xx)

    # This just sets the center pixel of the image
    origin = [origin_x, origin_y]
    if origin == [0,0]:
        ctrX = xx[int(np.floor(pix/2))]
        ctrY = ctrX
    else:
        ctrX = xx[origin[0]-1]
        ctrY = xx[origin[1]-1]

    # Calculate the distance from a point on the meshgrid to the center of image
    # Distance from point on meshgrid to center of x/y axes
    x = X - ctrX  # these are also meshgrids
    y = Y - ctrY

    # Set amplitudes to plot lattice as honeycomb or dots ... This basically just changes the phase of the atoms 
    if honeycomb == 1:
        A = 1/2
        B = -1/4
    else:
        A = 1/2
        B = 1/4

    # Create square lattice
    Z = A + B * (cos((k1[0]*x) + (k1[1]*y)) + cos((k2[0]*x) + (k2[1]*y)))

    # FFT of lattice:
    fftZ = np.abs(npf.fftshift(npf.fft2(Z- np.mean(np.mean(Z))))) # subtract by the mean to get rid of bright constant background peak at k=0
  

    return Z, fftZ