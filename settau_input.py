import numpy as np
import sys
import os
import subprocess
import time

subprocess.call(["f2py", "kmc_compute_settau.f90", "-h", "kmc_rates.pyf", "-m", "kmc_rates"])  # signature file generation

subprocess.call(["f2py", "-c", "-m", "kmc_rates", "--f90exec=ifort", "--opt=-O3", "--f90flags='-mcmodel=medium'", "kmc_compute_settau.f90"])
 
time.sleep(5)

import kmc_rates
from kmc_rates import kmc_compute_rates

coords_npy = np.load("input.npy") # sample input system corresponding to a PEO length
box_len = int((len(coords_npy))**(float(1.0)/float(3.0))) + 1 # Taking the cube root returns a float ending in .99999...
Lx = int(box_len) # Number of lattice points in the x, y, or z directions
Ly = int(box_len)
Lz = int(box_len) 
delta = np.ndarray(shape = (3), dtype = float) # values for lattice spacing in the x, y, and z directions 
delta[0] = 1.0 
delta[1] = 1.0
delta[2] = 1.0
Domain = np.empty([Lx,Ly,Lz])  # PS block = 0, Solvent/PEO block = 1 
density_input = np.empty([Lx,Ly,Lz])
dsite_input = np.empty([Lx,Ly,Lz])

# coords_npy is given in the shape(299^3, 6) and must be mapped onto the multidimensional array formatt
x = 0
for i in range(Lx):
	for j in range(Ly):
		for k in range(Lz):
		
			Domain[i,j,k] = coords_npy[x,3]
                        density_input[i,j,k] = coords_npy[x,4]
			dsite_input[i,j,k] = coords_npy[x,5]
			x += 1

 
def compute_tau_rates(Domain, dsite_input, density_input, Lx, Ly, Lz, delta):
	
	kmc_rates.kmc_compute_rates.delta = delta # delta is inputed directly outside of the FORTRAN subroutine due to problems trying to pass it into the subroutine
	rates = kmc_rates.kmc_compute_rates.settau_calc(Domain, dsite_input, density_input, Lx, Ly, Lz) # passing values and calling the subroutine 
        return rates 
	

settau_rates = compute_tau_rates(Domain, dsite_input, density_input, Lx, Ly, Lz, delta) 
np.save("settau_rates", settau_rates) # saves calculated rates in .npy formatt



