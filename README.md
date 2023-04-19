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

the "06" in this case relates to the mesh size determined from the parameter g_max_lvl=6 in CASE.inp. g_max_lvl=7 being finer.

The result seen in the .vtk files for the given CASE.inp and amr_runtime_parameters input files is (eventually) a cubic "hopper" crystal

The CHK.out file contains a single number, e.g. 97000. This number can be use to restart a run from time step 97000 by first copying:

cp paramesh_chk_000001 paramesh_chk_097000

and then running...

C01]$ mpirun -n 10 ./IM_AlNi3D.intel 2 > out.txt2& 

Note the "2" in the command line indicating a restart from the value held in the CHK.out file. 

The reason for this (copying of the restart file paramesh_chk_000001) procedure is to avoid build up of many restart file tyaking up hard disk space.


