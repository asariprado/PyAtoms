# PyAtoms


Simulates scanning probe microscopy (SPM) images

(Formerly named: SPM Simulator, Atom Simulator)

<img width="900" alt="image" src="https://user-images.githubusercontent.com/62832051/186005147-f2597f19-aa75-4967-9a7f-ddda06bf7613.png">


### Dependencies:
- Python 3
- NumPy
- Matplotlib (optimized for versions < 3.5.1)
- PyQt5 

## Installation instructions - Mac OS and Windows. (Not tested on Linux systems)
1. Download AtomSimulator folder (via Code-> Download .zip)
2. Open terminal (Mac OSX) or command line (Windows).
3. Navigate to the directory where the AtomSimulator folder is located. Example:
    ```
    ~/Downloads/AtomSimulator
    ```
5. To run program, type: 
    ```
    python AtomSimulator_GUI.py
    ```
## Known issues

The most recent version of Matplotlib (3.5.1) on Mac OS systems creates figures with overly large fonts. To correct this, AtomSimulator detects the user's operating system and version of Matplotlib and corrects the fonts and figure properties "by hand." We are monitoring this for a better permanent fix.

### For windows users:
- Make sure python is installed and that its path is set in your environment
- To check if it is, open the command line and type
    ```
    python -V
    ```
- Alternatively, if you installed python, NumPy, SciPy, etc. through Anaconda for Windows, you can run the above code through the Anaconda prompt.

For any other issues or crash reports, contact asariprado@physics.ucla.edu

##
## How to use 

Note that all fields accept typical mathematical operations in python and NumPy such as `+` `-` `*` `/` `sqrt` `log` and all valid NumPy functions `func` can be called via `np.func()`. 

1. Moiré lattice
    - Choose to simulate a 1, 2 or 3 layer lattice
       - Lattice 1 parameters change the single/first layer.
       - Lattice 2 only works if bilayer/trilayer are selected. These change the second lattice.
       - Lattice 3 only works if trilayer is selected. These change the third lattice.
   - `eta`, $\eta$ : A phenomenological parameter we use to weigh the relative strength of the sum of lattices, $Z_1 + Z_2$, versus the product of lattices, $Z_1 * Z_2$. $\eta$ is a real number between 0 and 1: The moiré image for $\eta=1$ is purely the product and $\eta=0$ purely the sum.
   - If $Z_1$ and $Z_2$ are two periodic lattices, the superimposed moiré/superlattice image is approximated as $Z_{moire} = (1-\eta)(Z_1 + Z_2) + \eta Z_1Z_2$ This toy model provides a good match to both experimental STM images and their Fourier transforms.
   - Similarly, for three periodic lattices, the moiré/superlattice image is approximated as $Z_{moire} = (1-\eta)(Z_1 + Z_2 + Z_3) + \eta Z_1Z_2Z_3$.
   - For an alternative method for simulating SPM images that may be implemented in future releases, see F. Joucken, et al. <a href="https://doi.org/10.1016/j.carbon.2014.11.030" target="_blank" rel="noopener noreferrer"> *Carbon* **83**, 48 (2015). </a>

2. Image parameters
    - `pixels`: number of pixels.
    - `L`: length of the image window in nanometers.
    - `theta`: Rotation (counter-clockwise) of the atomic lattice.

3. Colormap
    -  `Real space image`: colormap of the real space image / simulated atomic lattice.
    -  `FFT`: colormap of the 2D fast Fourier transform.
    
4. Low pass filtering
    - `sigma`: radius of Gaussian mask in real space (pixels) units. The half-width at half-max of the gaussian is shown as a white circle in the bottom left corner of the image.
    - The radius of the gaussian mask in nanometers is shown in the text box.

5. Save files
    - Click to save the files to a specific directory.
      - Clicking will open the file explorer.
      - Navigate to the directory you want to save the files in
      - Input a `filename` to save as
    - It will save a folder called `filename` with 4 files:
      - .png images of the real space and FFT, .txt file of the real space image, .txt file of the parameter values
      - `filename.png`, `filename_FFT.png`, `filename.txt`, `filename_params.txt`
      - Example of folder with saved files: 
  
        <img width="700" alt="Screen Shot 2022-08-04 at 1 30 29 PM" src="https://user-images.githubusercontent.com/62832051/182946811-ba2a1e4d-04d7-4658-b013-38dac1c8ef42.png">
      - Example of params .txt file: 
  
        <img width="446" alt="Screen Shot 2022-08-03 at 3 40 46 PM" src="https://user-images.githubusercontent.com/62832051/182724722-b820f3b3-a2e2-413c-8c9d-cb03da7b78ce.png">

      
6. SPM image time estimator
    - `Tip velocity`: Velocity of the tip, $v_t$ across 1 line in nanometers/second. The total image time is calculated by considering the total pixels, scanning left/right and in a single slow-scan direction (upwards, for example). This time is estimated via $T_{im} = (2 * N_{pix} * L) / v_t$, where $N_{pix}$ is the number of pixels in the image, $L$ is the length of the image in nanometers, and the factor 2 takes into account the left/right fast-scan direction. 
    
7. Spectroscopy map Time estimator
    - `Time per spectra`: Total time in seconds (including overhead) to record a single spectrum, $T_{spec}$ i.e. one dI/dV(V) sweep. We assume that spectra are recorded in one direction along both the fast- and slow-scan direction. The time is estimated via $T_{map} = N_{pix}^2 * T_{spec} + (2 * N_{pix} * L) / v_t$.

8. Lattices
    - Parameters tab
       - `symmetry`: Choose to simulate either a hexagonal (triangular/honeycomb) or square lattice.
       - `Lattice constant`: periodicity/spacing between atoms in nanometers.
       - `Twist angle`: twists the second lattice with respect to the first lattice (in Lattice 2 params) // twists the third lattice with respect to the second lattice (in Lattice 3 params)

    - Sublattices tab -- only affects hexagonal (triangular/honeycomb) lattices
       - Lattice site at the image `origin`: choose whether the `origin` should be a hollow site, an A-site atom or a B-site atom
          - To test this, set `L = 1` and click the different options for the origin
          <img width="400" alt="Screen Shot 2022-10-03 at 4 26 18 PM" src="https://user-images.githubusercontent.com/62832051/193703355-855b46de-f020-428f-af81-0ae4fee0bf57.png">

          
      - Weight of sublattices:
          - `alpha1`, $\alpha_1$ : weight of A sublattice
          - `beta1`, $\beta_1$: weight of B sublattice 
      - For a triangular lattice: `alpha = 1`, `beta = 0`
      
        ![triangular](https://user-images.githubusercontent.com/62832051/183219252-90edd400-bd36-46c0-9e39-e5d0c4b0e4c0.png)

      - For a honeycomb lattice: `alpha = 1`, `beta = 1`
      
        ![honeycomb](https://user-images.githubusercontent.com/62832051/183219271-329a51b9-b0b9-4e44-a7b2-04bf34960c7c.png)
     
    - Strain tab
        - Apply the 2D strain tensor, $e_{xy}$, where the $x$-axis is defined as the *local* direction of that lattice, i.e. the strain tensor rotates with the local axes set by `theta` or `twist angle`. See https://doi.org/10.1103/PhysRevB.80.045401.


##
## Examples
1. Twisted bilayer graphene with $1.1^\circ$ twist angle.
    - `Moire lattice`: bilayer
    - `L = 35`
    - `Pixels = 1024`
    - Lattice 1:
      - `Hexagonal`, `a = 0.3`, `A-site`, `alpha1 = 1`, `beta1 = 1`
    - Lattice 2:
       - `Hexagonal`, `b = 0.3`, `Twist angle = 1.1`, `A-site`, `alpha2 = 1`, `beta2 = 1` 
      
    <img width="260" alt="Screen Shot 2022-08-05 at 1 41 35 PM" src="https://user-images.githubusercontent.com/62832051/183159495-dc4b4c38-5e67-4cbc-bbe9-55ff9b696ec6.png">

    <img width="280" alt="Screen Shot 2022-08-03 at 3 43 14 PM" src="https://user-images.githubusercontent.com/62832051/182725026-1a462df7-9372-4b7e-bd02-6041134966b7.png">



2. 1T-TaS2 with $(\sqrt{13}\times\sqrt{13})R13.9^\circ$ charge density wave superlattice.
    - `Moire lattice`: bilayer
    - `L = 7`
    - `Pixels = 256`
    - Lattice 1:
      - `Hexagonal`, `a = 0.3`, `A-site`, `alpha1 = 1`, `beta1 = 0`
    - Lattice 2:
       - `Hexagonal`, `b = 0.3*np.sqrt(13)`, `Twist angle = 13.9`, `A-site`, `alpha2 = 1`, `beta2 = 0` 
       
    ![1T-TaS2](https://user-images.githubusercontent.com/62832051/182723975-b59e6b83-545a-47fe-8a68-59146fc1879b.png)
    ![1T-TaS2_FFT](https://user-images.githubusercontent.com/62832051/182723980-90b7689d-55c0-466d-8fd7-f993574f8955.png)


3. 2H-NbSe2 with $(3\times 3)R0^\circ$ charge density wave superlattice.
    - `Moire lattice`: bilayer
    - `L = 7`
    - `Pixels = 256`
    - Lattice 1:
      - `Hexagonal`, `a = 0.3`, `A-site`, `alpha1 = 1`, `beta1 = 0`
    - Lattice 2:
       - `Hexagonal`, `b = 0.3*3`, `Twist angle = 0`, `A-site`, `alpha2 = 1`, `beta2 = 0` 
      
   ![2H-NbSe2](https://user-images.githubusercontent.com/62832051/182723639-dc7b7277-1328-4ecd-8913-8428cc38331f.png)
   ![2H-NbSe2_FFT](https://user-images.githubusercontent.com/62832051/182723651-0ced0fed-5f33-4a78-bbc4-d9bf321e9811.png)


4. Kekule-O (trivial) distorted graphene with $(\sqrt{3}\times\sqrt{3})R30^\circ$ superlattice.
    - `Moire lattice`: bilayer
    - `L = 7`
    - `Pixels = 256`
    - Lattice 1:
      - `Hexagonal`, `a = 0.3`,  `Hollow`, `alpha1 = 1`, `beta1 = 1`
    - Lattice 2:
       - `Hexagonal`, `b = 0.3*sqrt(3)`, `Twist angle = 30`, `Hollow`, `alpha2 = 1`, `beta2 = 0` 
       
    ![Kekule-O trivial](https://user-images.githubusercontent.com/62832051/182722880-113f3ada-2199-4fe0-926f-73e645fda904.png)
    ![Kekule-O trivial_FFT](https://user-images.githubusercontent.com/62832051/182722894-8a1e5cc1-afaa-41c3-8d5a-ee098aed254d.png)


 5. Kekule-O (topological) distorted graphene with $(\sqrt{3}\times\sqrt{3})R30^\circ$ superlattice.
    - `Moire lattice`: bilayer
    - `L = 7`
    - `Pixels = 256`
    - Lattice 1:
      - `Hexagonal`, `a = 0.3`,  `Hollow`, `alpha1 = 1`, `beta1 = 1`
    - Lattice 2:
       - `Hexagonal`, `b = 0.3*sqrt(3)`, `Twist angle = 30`, `Hollow`, `alpha2 = 1`, `beta2 = 1`
       
    ![Topological](https://user-images.githubusercontent.com/62832051/182723054-ede96db1-1f19-4eb9-ab59-3f8a60c52b32.png)
    ![Topological_FFT](https://user-images.githubusercontent.com/62832051/182723064-db5399fb-3dcb-4a23-890f-95ff5b65e9d8.png)

