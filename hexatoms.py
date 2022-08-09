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
from numpy import pi as pi
from numpy import exp as exp


def hexatoms(pix, L, a, theta, e11, e12, e22, alpha, beta, origin): 
    """
    
    Returns np.arrays of  hexagonal atomic lattice and its FFT (absolute value)

    Creates a (pix x pix) image array of a hexagonal atomic lattice with periodicity 'a' (nm),
    rotated from the x-direction by an angle theta. The image has side length L (nm).
    The image height is normalized to [0,1].

    
    Args:
        pix [int]: number of pixels to make the simulated lattice (will be a NxN pixel image)
        L [float]: side length of the sample (in real space/ real units) [nm]
        a [float]: periodicity/spacing between atoms/lattice parameter (in real space/ real units) [nm]
        theta [degrees]: angle to rotate lattice by theta degrees from the x-direction. Default is 0, no rotation
        e11, e22, e12 [float]: STRAIN tensor components. the amount of strain/distortion to add to each component (default is 0 for each - no distortion)
                               The strain tensor is:   e = (e11  e12) 
                                                           (-e12 e22)
        alpha, beta [float]: amplitudes of A and B sublattices
        origin [str]: "Hollow", "A-site" or "B-site". to choose whether a hollow spot or an A/B atom is at the center of the image 
        
    """
    
    
    ## Create 2D meshgrid of points. Make sure there is *always* a dedicated pixel for the origin (0,0) for pix even or odd
    x = np.arange(-(pix//2),(pix-1)//2 + 1)*L/(pix-1) 
    [X,Y] = np.meshgrid(x,x) 
    # ^^ IMPORTANT: make sure its typed [X,Y] instead of [Y,X], otherwise you might need to transpose the output
    
    ## Reciprocal lattice vectors for hexagonal crystal WITH STRAIN defined in local coordinates
    # to see how the strain tensor is defined and how this affects the recip lattice vectors, 
    # look at paper: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.80.045401 
    k1 = (2*np.pi/a)* np.array([1 - e11 - (e12/m.sqrt(3)), (1/m.sqrt(3)) - e12 - (e22/m.sqrt(3))])
    k2 = (2*np.pi/a)* np.array([-1 + e11 - (e12/m.sqrt(3)), (1/m.sqrt(3)) + e12 - (e22/m.sqrt(3))])
    k3 = -(k1 + k2) # Get the 3rd wavevector by taking a linear combo of k1,k2 
    
 
    ## Rotate the lattice by theta (in degrees)
    # Convert theta to radians
    theta_rad = np.deg2rad(theta) 

    # Create rotation matrix to multiply reciprocal lattice vectors to rotate the image 
    rotmat = np.array([[cos(theta_rad), -sin(theta_rad)],        # Rotation matrix =  [[cos -sin] 
                       [sin(theta_rad), cos(theta_rad)]])		 #                     [sin, cos]]

    # Rotate wavevectors by angle theta, using rotation matrix 
    k1 = np.matmul(k1,rotmat)
    k2 = np.matmul(k2,rotmat)
    k3 = np.matmul(k3, rotmat)


    # Create variables for each wavevector component - For CONVENIENCE only 
    # (just to make the code for T below look less cluttered)
    k1x, k1y = k1[0], k1[1]
    k2x, k2y = k2[0], k2[1]
    k3x, k3y = k3[0], k3[1]


    ## Create plain lattice
    T = exp(1j*(k1x*X + k1y*Y)) + exp(1j*(k2x*X + k2y*Y)) + exp(1j*(k3x*X + k3y*Y))

    ## Amplitude (alpha, beta) and phase-offsets on the A and B sublattices, respectively.
    phase = alpha + beta*exp(-1j*2*pi/3)


    ## To shift origin to hollowsite or B-sublattice
    if origin == "A-site":
        phase = phase # Keep as is
        
    elif origin == "Hollow":
        phase *= exp(-1j*2*pi/3) # Muliply by phase shift to bring hollow site down to origin

    elif origin == "B-site":
        phase *= exp(-1j*2*pi/3)**2 # Multiply by phase shift again to bring B-site to origin


    ## Multiply the lattice by the phase , add the conjugate to get rid of imaginary.
    # (The real part is taken since the result has a 0j complex component.)
    Z_un = np.real(T*phase + np.conjugate(T*phase))

    # relative strength of the sublattices, put a number btw  0 & 1

    ## Normalize image 
    # First, avoid division by zero
    if alpha == 0 and beta == 0:
        Z = np.ones((pix,pix))

    # Then normalize so that image values are between 0<Z<1.
    else:
        Z = (Z_un - np.min(np.min(Z_un)))/(np.max(np.max(Z_un)) - np.min(np.min(Z_un))) 



    ## Take 2d FFT of the lattice 
    fftZ = np.abs((npf.fftshift(npf.fft2(Z - np.mean(np.mean(Z))))))
    # subtract by mean(mean(Z)) to remove the strong peak at k=0 (DC/constant background)

    # Normalize the FFT to be in between 0-1 values
    # Subtract by the minimum value and then divide by the maximum value:
    # fftZ_norm = (fftZ - np.min(np.min((fftZ))))/(np.max(np.max(fftZ)) - np.min(np.min(fftZ)))
    ## ^^ commented out bc the normalization is done in moirelattice.py for bilayer/trilayer, and in _Widgets.py for single layer

    return Z, fftZ
