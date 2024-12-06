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


def moirelattice(pix, L, a1, a2, a3, moireBtn, modeBtn, lattice1, lattice2, lattice3, theta_offset, theta_tw, theta_tw2, e11, e12, e22, d11, d12, d22, f11, f12, f22,  alpha1, beta1, alpha2, beta2, alpha3, beta3, eta, xi, origin1, origin2, origin3, filter_bool, sigma,center):
   
    ### Define the rotation angles

    # Rotate the first lattice (offset rotation of the whole image)
    theta_im = theta_offset 
    
    # Rotate the second lattice wrt to the first one (add theta_offset + theta_twist)
    theta_tw12 = theta_tw + theta_im

    # Rotate the third lattice wrt to the second one 
    theta_tw23 = (theta_tw2 + theta_tw + theta_im)


    ## CREATE FIRST LATTICE ##
    if lattice1 == 'Hexagonal':
        Z1, fftZ1 = hexatoms(pix, L, a1, theta_im, e11, e12, e22, alpha1, beta1, origin1,center)

    elif lattice1 == 'Square':
        Z1, fftZ1 = squareatoms(pix, L, a1, theta_im, d11, d12, d22,center)


    ## CREATE SECOND LATTICE ##
    if lattice2 == 'Hexagonal':
        Z2, fftZ2 = hexatoms(pix, L, a2, theta_tw12, d11, d12, d22, alpha2, beta2, origin2,center)
   
    elif lattice2 == 'Square':
        Z2, fftZ2 = squareatoms(pix, L, a2, theta_tw12, d11, d12, d22,center)
  

    # If only doing bilayer, create the moire lattice: Z = Z1*Z2 and filter it if filter btn is checked
    if moireBtn == 'Bilayer':
        if modeBtn == 'Simple':
            if eta<0.0:
                eta = 0.0
            elif eta>1.0:
                eta = 1.0

            Z = (eta * Z1 * Z2) + (1-eta)*(Z1 + Z2)

        if modeBtn == 'Log':
            if xi<0.0:
                xi = 0.0
            elif xi>10.0:
                xi = 10.0

            # Add small non-zero constant to avoid -inf from log 0
            Z = np.log(1e-6 + Z1 + Z2*np.exp(-xi))


        if filter_bool == True: 

            Z = gaussian_filter(Z, sigma,mode='mirror')
       
        fftZ = np.abs(npf.fftshift(npf.fft2(Z - np.mean(np.mean(Z)))))

        fftZ = mat2gray(fftZ)

    
    
    ## CREATE THIRD LATTICE IF TRILAYER MOIRE BUTTON IS CHOSEN ## 
    elif moireBtn == 'Trilayer':
        if lattice3 == 'Hexagonal':
            Z3, fftZ3 = hexatoms(pix, L, a3, theta_tw23, f11, f12, f22, alpha3, beta3, origin3,center)
            
        elif lattice3 == 'Square':
            Z3, fftZ3 = squareatoms(pix, L, a3, theta_tw23, f11, f12, f22,center)
        
        if modeBtn == 'Simple':
            if eta<0.0:
                eta = 0.0
            elif eta>1.0:
                eta = 1.0  
      
            Z = (eta * Z1 * Z2 * Z3) + (1-eta)*(Z1 + Z2 + Z3)

        if modeBtn == 'Log':
            if xi<0.0:
                xi = 0.0
            elif xi>10.0:
                xi = 10.0

            Z = np.log(1e-6 + Z1 + Z2*np.exp(-xi) + Z3*np.exp(-2*xi))

      # Low pass filter the image if the button is checked:
        if filter_bool == True: 

            Z = gaussian_filter(Z, sigma,mode='mirror') # filter the stacked 3 lattices

        fftZ = np.abs(npf.fftshift(npf.fft2(Z - np.mean(np.mean(Z)))))


        fftZ = mat2gray(fftZ)
 
    # Normalize the real space image
    Z = mat2gray(Z)


    return Z, np.abs(fftZ)


# explicit function to normalize the 2D matrix.
def mat2gray(Z_un):
    Z = (Z_un - np.min(np.min(Z_un)))/(np.max(np.max(Z_un)) - np.min(np.min(Z_un)))
    return Z




