 # -*- coding: utf-8 -*-
"""
ATOM SIMULATOR Squareatoms
Created on Mon Nov 15 14:45:06 2021
@author: Asari
"""

import numpy as np

from numpy import fft as npf
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


def moirelattice(pix, L, a1, a2, a3, moireBtn, lattice1, lattice2, lattice3, theta_offset, theta_tw, theta_tw2, e11, e12, e22, d11, d12, d22, f11, f12, f22,  alpha1, beta1, alpha2, beta2, alpha3, beta3, eta, origin1, origin2, origin3, filter_bool, sigma):
   
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

        Z = (eta * Z1 * Z2) + (1-eta)*(Z1 + Z2)

        if filter_bool == True: 

            Z = gaussian_filter(Z, sigma,mode='mirror')
       
        fftZ = np.abs(npf.fftshift(npf.fft2(Z - np.mean(np.mean(Z)))))

        fftZ_norm = (fftZ - np.min(np.min((fftZ))))/(np.max(np.max(fftZ)) - np.min(np.min(fftZ)))
        fftZ = fftZ_norm

    
    
    ## CREATE THIRD LATTICE IF TRILAYER MOIRE BUTTON IS CHOSEN ## 
    elif moireBtn == 'Trilayer': 
        if lattice3 == 'Hexagonal':
            Z3, fftZ3 = hexatoms(pix, L, a3, theta_tw23, f11, f12, f22, alpha3, beta3, origin3)

            
        elif lattice3 == 'Square':
            Z3, fftZ3 = squareatoms(pix, L, a3, theta_tw23, f11, f12, f22)
      

        Z = (eta * Z1 * Z2 * Z3) + (1-eta)*(Z1 + Z2 + Z3)

      # Low pass filter the image if the button is checked:
        if filter_bool == True: 

            Z = gaussian_filter(Z, sigma,mode='mirror') # filter the stacked 3 lattices

        fftZ = np.abs(npf.fftshift(npf.fft2(Z - np.mean(np.mean(Z)))))


        fftZ_norm = (fftZ - np.min(np.min((fftZ))))/(np.max(np.max(fftZ)) - np.min(np.min(fftZ)))
        fftZ = fftZ_norm
 
    # Normalize the real space image
    Z = (Z - np.min(np.min(Z)))/(np.max(np.max(Z)) - np.min(np.min(Z))) 


    return Z, np.abs(fftZ)


# explicit function to normalize array. from https://www.geeksforgeeks.org/how-to-normalize-an-array-in-numpy-in-python/
def normalize_2d(matrix):
    norm = np.linalg.norm(matrix)
    matrix = matrix/norm  # normalized matrix
    return matrix




