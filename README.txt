This version of the code and coords.npy generated are sepecifically for the random walk only. 
Certain calculations and sections of the generic code are deleted for slimplicity.

To run the random walk

1. run input_gen.py to generate the coords.npy input (an array with hopping rates values of 1)
2. run set_input.py to compile the kmc_main_RW.f90 code into a wrapper and run it
