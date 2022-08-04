# AtomSimulator

Simulates SPM images ,,,,,,,,,,,,,,,,,,,,,,,,,,
<img width="900" alt="Screen Shot 2022-08-03 at 5 10 41 PM" src="https://user-images.githubusercontent.com/62832051/182736096-5b4b863d-10f3-4cc2-91c6-e5fd9825a752.png">

### Dependencies:
- python
- numpy
- matplotlib
- PyQt5


## Installation instructions
1. Download folder
2. Open terminal (Mac OSX) or command line (Windows).
3. Navigate to the directory where the AtomSimulator folder is located. Example:
    ```
    ~/Downloads/AtomSimulator
    ```
5. To run program, type: 
    ```
    python AtomSimulator_GUI.py
    ```


### For windows users:
- Make sure python is installed and that its path is set in your environment
- To check if it is, open the command line and type
    ```
    python -V
    ```
- Alternatively, if you installed python, NumPy, SciPy, etc. through Anaconda for Windows, you can run the above code through the Anaconda prompt.

##
## How to use ?

Note that all fields accept typical mathematical operations in python and NumPy such as `+` `-` `*` `/` `sqrt` `log` and all valid NumPy functions `func` can be called via `np.func()`. 

1. Moiré lattice
    - Choose to simulate a 1, 2 or 3 layer lattice
       - Lattice 1 parameters change the single/first layer
       - Lattice 2 only works if bilayer/trilayer are selected. These change the second lattice
       - Lattice 3 only works if trilayer is selected. These change the third lattice


2. Image parameters
    - `pixels`: number of pixels
    - `L`: length of image in nanometers
    - `theta`: Rotation of image

3. Colormap
    -  `Real space image`: colormap of the real space image / simulated atomic lattice
    -  `FFT`: colormap of the FFT
    
4. Low pass filtering
    - `sigma`: radius of Gaussian mask in real space (pixels)
    - The number on top tells you the radius of the mask in nanometers

5. Save files
    - Click to save the files to a specific directory. Will save a folder with 4 files:
      - .png images of the real space and FFT, .txt file of the real space image, .txt file of the parameter values
      - Example of params .txt file: <img width="446" alt="Screen Shot 2022-08-03 at 3 40 46 PM" src="https://user-images.githubusercontent.com/62832051/182724722-b820f3b3-a2e2-413c-8c9d-cb03da7b78ce.png">

      

6. SPM Time estimator
    - `Tip velocity`: Velocity of the tip across 1 line in nanometers/second.
    
7. dI/dV Time estimator
    - `Time per spectra`: Total time in seconds (including overhead) to record a single spectrum, i.e. one dI/dV(V) sweep. 

8. Lattices
    - Parameters tab
       - `symmetry`: Choose to simulate either a triangular/honeycomb or square lattice.
       - `Lattice constant`: periodicity/spacing between atoms in nanometers.
       - `Twist angle`: twists the second lattice with respect to the first lattice (in Lattice 2 params) // twists the third lattice with respect to the second lattice (in Lattice 3 params)

    - Sublattices tab (only affects hexagonal lattices)
       - Lattice site at the image origin: choose whether the origin should be a hollow site, an A-site atom or a B-site atom
          - To test this, set `L = 1 nm` and click the different options for the origin
      - Weight of sublattices:
          - `alpha1`: weight of A sublattice
          - `beta1`: weight of B sublattice 
      - If `alpha = 1`, `beta = 0`, the lattice is triangular
      - If `alpha = beta = 1`, the lattice is honeycomb
    - Strain tab
        - Apply the 2D strain tensor, $e_{xy}$, where the $x$-axis is defined as the *local* direction of that lattice, i.e. the strain tensor rotates with the local axes set by `theta` or `twist angle`. See https://doi.org/10.1103/PhysRevB.80.045401.


##
## Examples
1. Twisted bilayer graphene with $1.1^\circ$ twist angle.
    - `Moire lattice`: bilayer
    - `L = 35`
    - `Pixels = 1024`
    - Lattice 1:
      - `Hexagonal`, `a = 0.3`, `A-site`, `alpha1 = 1`, `beta1 = 0`
    - Lattice 2:
       - `Hexagonal`, `b = 0.3`, `Twist angle = 1.1`, `A-site`, `alpha2 = 1`, `beta2 = 0` 
       
    <img width="256" alt="Screen Shot 2022-08-03 at 3 42 15 PM" src="https://user-images.githubusercontent.com/62832051/182725015-be33b834-783e-45d2-b6f8-b9d12843235b.png"> <img width="280" alt="Screen Shot 2022-08-03 at 3 43 14 PM" src="https://user-images.githubusercontent.com/62832051/182725026-1a462df7-9372-4b7e-bd02-6041134966b7.png">


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

