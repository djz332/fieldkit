import numpy as np
import sys
import os
import subprocess
import time
import shlex
 
subprocess.call(["f2py", "-c", "-m", "kmc_test", "--f90exec=ifort", "--opt=-O3", "--f90flags='-fopenmp -mcmodel=large'", "kmc_main_RW.f90"])    

time.sleep(5)

import kmc_test
from kmc_test import kmc_rw

coords_npy = np.load("coords.npy")
box_len = len(coords_npy[0])
totalmoves_input = 80000000                   #Total number of moves

coords_npy = np.reshape(coords_npy,(box_len**3,6)) # Contains hopping rates     

kmc_test.kmc_rw.l = box_len
kmc_test.kmc_rw.l_unscale = box_len
kmc_test.kmc_rw.settau_input = coords_npy
kmc_test.kmc_rw.totalmoves = totalmoves_input
kmc_test.kmc_rw.case_parameter = 0           #A parameter used to differentiate different cases E.g. changing this to 1 will create files likes trajectory_1_1.csv, trajectory_1_2.csv, etc.
kmc_test.kmc_rw.numtracer = 8192             #Total number of walkers
kmc_test.kmc_rw.num_status = 2               #Total number of states for walkers. In out system, this is to represent different blocks.
kmc_test.kmc_rw.msdflag = [1,1,1]            #Flags used to denote whether to consider transport in particular direction to calculate MSDs. E.g. 0 1 1 for lamellar, 0 0 1 for cylinder  
kmc_test.kmc_rw.msd_freq = 100000            #Interval at which trajectories are updated
kmc_test.kmc_rw.origin_freq = 6
kmc_test.kmc_rw.prob_freq = 1000000
kmc_test.kmc_rw.print_freq = 500000000
kmc_test.kmc_rw.dfit_points = 2500
kmc_test.kmc_rw.num_runs = 15
kmc_test.kmc_rw.msdfile_write_freq = 100
kmc_test.kmc_rw.normalizeflag = 0


print '             '

kmc_test.kmc_rw.create_files()

print 'Creating files'

kmc_test.kmc_rw.grid_generator()

print 'Lattice sites are assigned'

kmc_test.kmc_rw.status_allocation()

print 'Allocating site status'

kmc_test.kmc_rw.set_tau()

print 'Setting up hopping rates equal to 1 (for random walk)'

kmc_test.kmc_rw.randomly_insert_tracers()

print 'Insert tracer done'

kmc_test.kmc_rw.algorithm_variable_setting()

print 'Algorithm set'

kmc_test.kmc_rw.dynamics()

print 'Successful dynamics!'

print 'RUN SUCESSFUL!!!'

sys.exit()
