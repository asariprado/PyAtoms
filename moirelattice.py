 # -*- coding: utf-8 -*-
"""
ATOM SIMULATOR Squareatoms
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

from numpy import minimum as minn
from numpy import maximum as maxx 
from scipy.ndimage import gaussian_filter

from hexatoms import hexatoms
from squareatoms import squareatoms


def moirelattice(pix, L, a1, a2, a3, moireBtn, lattice1, lattice2, lattice3, theta_offset, theta_tw, theta_tw2, e11, e12, e22, d11, d12, d22, f11, f12, f22,  alpha1, beta1, alpha2, beta2, alpha3, beta3, origin1, origin2, origin3, filter_bool, sigma):
   
    ### Define the rotation angles

    # Rotate the first lattice (offset rotation of the whole image)
    theta_im = theta_offset 
    
    # Rotate the second lattice wrt to the first one (add theta_offset + theta_twist)
    theta_tw12 = theta_tw + theta_im

    # Rotate the third lattice wrt to the second one 
    theta_tw23 = (theta_tw2 + theta_tw + theta_im)


    ## CREATE FIRST LATTICE ##
    if lattice1 == 'Hexagonal':
        Z1, fftZ1 = hexatoms(pix, L, a1, theta_im, e11, e12, e22, alpha1, beta1, origin1)

    elif lattice1 == 'Square':
        Z1, fftZ1 = squareatoms(pix, L, a1, theta_im, d11, d12, d22)


    ## CREATE SECOND LATTICE ##
    if lattice2 == 'Hexagonal':
        Z2, fftZ2 = hexatoms(pix, L, a2, theta_tw12, d11, d12, d22, alpha2, beta2, origin2)
   
    elif lattice2 == 'Square':
        Z2, fftZ2 = squareatoms(pix, L, a2, theta_tw12, d11, d12, d22)
  

    # If only doing bilayer, create the moire lattice: Z = Z1*Z2 and filter it if filter btn is checked
    if moireBtn == 'Bilayer': 
        if filter_bool == True: 
            # Filter the stacked lattices: multiply the two and then add each one individually. divide by 3 to normalize bc each one has max of 1, so 1*1 + 1 + 1 = 3 
            Z = gaussian_filter((Z1 * Z2 + Z1 + Z2)/3, sigma,mode='mirror')
            fftZ = np.abs(npf.fftshift(npf.fft2(Z - np.mean(np.mean(Z)))))
            # Normalize the FFT:
            # fftZ_norm = (fftZ - np.min(np.min((fftZ))))/(np.max(np.max(fftZ)) - np.min(np.min(fftZ)))


        else:
            # Just stack the 2 lattices to create moire superlattice
            Z = (Z1 * Z2 + Z1 + Z2)/3
            # Take the fft     
            fftZ = np.abs(npf.fftshift(npf.fft2(Z - np.mean(np.mean(Z)))))    
            # Normalize the FFT to be betweeen 0-1. this is so the plots dont move around much 
            # fftZ_norm = (fftZ - np.min(np.min((fftZ))))/(np.max(np.max(fftZ)) - np.min(np.min(fftZ)))

    
    ## CREATE THIRD LATTICE IF TRILAYER MOIRE BUTTON IS CHOSEN ## 
    elif moireBtn == 'Trilayer': 
        if lattice3 == 'Hexagonal':
            Z3, fftZ3 = hexatoms(pix, L, a3, theta_tw23, f11, f12, f22, alpha3, beta3, origin3)

            
        elif lattice3 == 'Square':
            Z3, fftZ3 = squareatoms(pix, L, a3, theta_tw23, f11, f12, f22)
      
      # Low pass filter the image if the button is checked:
        if filter_bool == True:
            Z = gaussian_filter((Z1 * Z2 * Z3 + Z1 + Z2 + Z3)/4, sigma,mode='mirror') # filter the stacked 3 lattices
            fftZ = np.abs(npf.fftshift(npf.fft2(Z - np.mean(np.mean(Z)))))

        else:
            # Stack the 3 lattices to create moire superlattice
            Z = (Z1 * Z2 * Z3 + Z1 + Z2 + Z3)/4 #Normalize after this line, but need if statement for alpha AND beta ==0
            

            # Take the fft     
            fftZ = np.abs(npf.fftshift(npf.fft2(Z - np.mean(np.mean(Z)))))   
            # fftZ_norm = (fftZ - np.min(np.min((fftZ))))/(np.max(np.max(fftZ)) - np.min(np.min(fftZ)))
 


    return Z, np.abs(fftZ)


# explicit function to normalize array. from https://www.geeksforgeeks.org/how-to-normalize-an-array-in-numpy-in-python/
def normalize_2d(matrix):
    norm = np.linalg.norm(matrix)
    matrix = matrix/norm  # normalized matrix
    return matrix




