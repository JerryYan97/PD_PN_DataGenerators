# PD_PN_DataGenerators

### Overview

This project is intended for generating training data for Projected Newton method FEM simulation
and Projective Dynamics simulation in 3D. 

### How to build it

It is a CMAKE project with several 3rd party libraries attached. 
Besides, you also need to download and install the Intel's OneAPI to get MKL and TBB supports. 
In addition, you also need to download and build the libigl lib.
Then, you need to specify the folder of the MKL, TBB and libigl lib folder in your system. In order to link this 
project to these packages, you will need to specify three environment variables:

1. ```MKL_ROOT=/opt/intel/oneapi/mkl/YOUR VERSION NUMBER```

2. ```TBB_ROOT="/opt/intel/oneapi/tbb/YOUR VERSION NUMBER"```
   
3. ```LIBIGL_DIR = YOUR Libigl PROJECT ROOT PATH```

### Reference, Third Parties and Acknowledgments

1. LIBIGL
2. Intel MKL, TBB
3. JSON parser
4. https://github.com/ipc-sim/IPC
5. UPenn CG Lab