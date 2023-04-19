# PhaseField
To load and compile the code the script $ . Paramesh.scr

To execute: create a run directory, say, "RUNS/C01"

Copy files from directory, Campfire:

campfire]$ cp IM_AlNi3D.intel IM_AlNi3D-SRC/CASE.inp IM_AlNi3D-SRC/amr_runtime_parameters ../RUNS/C01/.

We need mpi ...

C01]$ module load mpi/openmpi-x86_64

and an intel compiler

C01]$ module load intel-compiler

Then to run:

C01]$ mpirun -n 10 ./IM_AlNi3D.intel 

or

C01]$ mpirun -n 10 ./IM_AlNi3D.intel > out.txt&
