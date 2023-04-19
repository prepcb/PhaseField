# PhaseField
To load and compile the code the script $ . Paramesh.scr

To execute: create a run directory, say, "RUNS/C01"

Copy files from directory, Campfire:

campfire]$ cp IM_AlNi3D.intel IM_AlNi3D-SRC/CASE.inp IM_AlNi3D-SRC/amr_runtime_parameters ../RUNS/C01/.

We need mpi ...

C01]$ module load mpi/openmpi-x86_64

and an intel compiler

C01]$ module load intel-compiler

Then to run a 10 core job:

C01]$ mpirun -n 10 ./IM_AlNi3D.intel 

or to send output to out.txt

C01]$ mpirun -n 10 ./IM_AlNi3D.intel > out.txt&

This will create files: out.txt,tip_loc.txt, paramesh_chk_000001, CHK.out and

iso_06_XXX000_YYY.vtk

where XXX is in the range 001 to 999 (max for paraview .vtk files), and YYY is from 000 to 009 (for the 10 cores)

The result seen in the .vtk files for the given CASE.inp and amr_runtime_parameters input files is a cubic "hopper" crystal
