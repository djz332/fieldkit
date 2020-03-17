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
totalmoves_input = 80000000
#print box_len
#print np.shape(coords_npy)

coords_npy = np.reshape(coords_npy,(box_len**3,6))
#print np.shape(coords_npy)

kmc_test.kmc_rw.l = box_len
kmc_test.kmc_rw.l_unscale = box_len
kmc_test.kmc_rw.settau_input = coords_npy
kmc_test.kmc_rw.totalmoves = totalmoves_input
kmc_test.kmc_rw.numtracer = 8192
kmc_test.kmc_rw.num_status = 2
kmc_test.kmc_rw.msdflag = [1,1,1]
kmc_test.kmc_rw.msd_freq = 100000
kmc_test.kmc_rw.origin_freq = 6
kmc_test.kmc_rw.prob_freq = 1000000
kmc_test.kmc_rw.print_freq = 500000000
kmc_test.kmc_rw.dfit_points = 2500
kmc_test.kmc_rw.num_runs = 15
kmc_test.kmc_rw.msdfile_write_freq = 100
kmc_test.kmc_rw.normalizeflag = 0


kmc_test.kmc_rw.read_input()

print 'Read input'

kmc_test.kmc_rw.grid_generator()

print 'Lattice sites are assigned'

kmc_test.kmc_rw.status_allocation()

print 'Allocating site status'
#print kmc_test.algorithm_rw.totalsites

kmc_test.kmc_rw.set_tau()

print 'Differentials are caclulated/settau = 1'
#print kmc_test.inputs.settau_input
#print kmc_test.kmc_subs_rw.siteid_settau 

kmc_test.kmc_rw.randomly_insert_tracers()

print 'Intert tracer done'

kmc_test.kmc_rw.algorithm_variable_setting()

print 'Algorithm set'
#sys.exit()

kmc_test.kmc_rw.dynamics()

print 'Successful dynamics!'

kmc_test.kmc_rw.output_generator()

print 'RUN SUCESSFUL!!!'

sys.exit()
