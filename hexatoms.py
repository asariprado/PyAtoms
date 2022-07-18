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


def hexatoms(pix, L, a, theta, e11, e12, e22, honeycomb, origin_x = 0, origin_y = 0): 
    """
    
    Returns np.arrays of  hexagonal atomic lattice and its FFT

    Creates a (pix x pix) image array of a hexagonal atomic lattice with periodicity 'a' (nm),
    rotated from the x-direction by an angle theta. The image has side length L (nm).
    If honeycomb == 1, it creates the same image but with honeycomb pattern. 
    For either choice, the image height is normalized to [0,1].

    
    Args:
        pix [int]: number of pixels to make the simulated lattice (will be a NxN pixel image)
        L [float]: side length of the sample (in real space/ real units) [nm]
        a [float]: periodicity/spacing between atoms/lattice parameter (in real space/ real units) [nm]
        theta [degrees]: angle to rotate lattice by theta degrees from the x-direction. Default is 0, no rotation
        e11, e22, e12 [float]: STRAIN tensor components. the amount of strain/distortion to add to each component (default is 0 for each - no distortion)
                               The strain tensor is:   e = (e11  e12) 
                                                           (-e12 e22)
        honeycomb [0 or 1]: Creates the lattice with (1) or without (0) a honeycomb pattern
                            0: no honeycomb pattern. 1: honeycomb pattern
        origin_x, origin_y [int]: sets the origin of the lattice. Default is set to 0,0 ((not really that important?))
        
    """
    
    
    # Create 2D meshgrid of points
    x = np.linspace(0, L, pix) 
    [X,Y] = np.meshgrid(x,x) 
    ## ^^ IMPORTANT: make sure its typed [X,Y] instead of [Y,X], 
    ##    otherwise you might need to transpose the output
    
    
    # Set the origin/center of the image to be at (the pixel corresponding to?) L/2
    ctrX, ctrY = x[int(np.floor(pix/2))], x[int(np.floor(pix/2))]


    # Calculate the distance from a point on the meshgrid to the center of image
    # Distance from point on meshgrid to center of x/y axes
    x_i = X - ctrX    # these are also meshgrids
    y_i = Y - ctrY

    ## Reciprocal lattice vectors for hexagonal crystal WITH STRAIN
    k1 = (2*np.pi/a)* np.array([1 - e11 - (e12/m.sqrt(3)), (1/m.sqrt(3)) - e12 - (e22/m.sqrt(3))])
    k2 = (2*np.pi/a)* np.array([-1 + e11 - (e12/m.sqrt(3)), (1/m.sqrt(3)) + e12 - (e22/m.sqrt(3))])
    k3 = -(k1 + k2) # Get the 3rd wavevector by taking a linear combo of k1,k2 
    ## ADD STRAIN:  to see how the strain tensor is defined and how this affects the recip lattice vectors, 
    ## look at paper: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.80.045401 
 

    # Rotate the lattice by theta (assuming theta is in degrees)
    theta_rad = np.deg2rad(theta) # convert theta to radians

    # Create rotation matrix to multiply reciprocal lattice vectors to rotate the image 
    # Rotation matrix =  [[cos -sin 
    #                    [sin, cos]]
    rotmat = np.array([[cos(theta_rad), -sin(theta_rad)], 
                       [sin(theta_rad), cos(theta_rad)]])

    # Rotate wavevectors by angle theta, using rotation matrix 
    k1 = np.matmul(k1,rotmat)
    k2 = np.matmul(k2,rotmat)
    k3 = np.matmul(k3, rotmat)


    # Create variables for each wavevector component - For CONVENIENCE only 
    # (just to make the code for Z below look less cluttered)
    k1x, k1y = k1[0], k1[1]
    k2x, k2y = k2[0], k2[1]
    k3x, k3y = k3[0], k3[1]



    # Create the image with or without a honeycomb pattern:
    if honeycomb == 1: # Normalized amplitudes
        A1 = 1
        B1 = -2/9
    else:
        A1 = 0
        B1 = 2/9
    # ## This is un normalized -- comment out
    # if honeycomb == 1: # Creates the image with a honeycomb pattern
    #     A = -1
    # else:              # No honeycomb pattern (equivalent to plotting dots as atoms)
    #     A = 1


    # Create hexagonal lattice
    # Remember that electron wavefunctions are complex exponentials
    # Taking the real part only we get a cosine bc e^ikx = cos(kx) + isin(kx)
    # Multiply by 2 to account for both +k and -k frequencies]
    # Z = 2 * A * (cos((k1x * x_i) + (k1y * y_i)) 
    #            + cos((k2x * x_i) + (k2y * y_i)) 
    #            + cos((k3x * x_i) + (k3y * y_i))) 

    ## should it be this one or the commented out one above????
    Z = A1 + B1 * (3/2 + (cos((k1[0]*x_i) + (k1[1]*y_i)) + (cos((k2[0]*x_i) + (k2[1]*y_i))) + (cos((k3[0]*x_i) + (k3[1]*y_i)))))




    #######################################   
    ### Take the 2D FFT of the lattice ####
    ####################################### 
    fftZ = np.abs((npf.fftshift(npf.fft2(Z - np.mean(np.mean(Z))))))
    # subtract by mean(mean(Z)) to remove the strong peak at k=0 (DC/constant background)

    # Normalize the FFT to be in between 0-1 values
    # Subtract by the minimum value and then divide by the maximum value:
    # fftZ_norm = (fftZ - np.min(np.min((fftZ))))/(np.max(np.max(fftZ)) - np.min(np.min(fftZ)))


    return Z, np.abs(fftZ)# fftZ_norm

