!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  input_parameter_file_reading 
!!!!  * Reads in parameters from the file 'CASE.inp' that may be set at runtime
!!!!  * To add new ones just include any module variable in the set defined by 
!!!!    "namelist/input_variables/"
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Latest version: Chris Goodyer, September 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine input_parameter_file_reading()
  ! include modules for governing solution parameters
  use multigrid_parameters
  use time_dep_parameters
  use refinement_parameters
  use generic_parameters

  ! include file to define physical quantities of the model
  use solution_parameters
  ! include paramesh data for nvar
  use paramesh_dimensions

  implicit none

  include 'mpif.h'

  integer :: pid, ierr
  logical :: there

  ! Make list of parameters available
     ! multigrid_parameters
  namelist/input_variables/ weight, smooth_count_max, solve_count_max, mg_min_lvl, mg_max_level, &
                            max_v_cycle_count, defect_tol_min, defect_too_big, defect_tol_max, &
                            low_vcycle_count, high_vcycle_count, C2F_start_level, C2F_desired_level, &
                            g_lrefine_min,g_grid_max,gg_Nbx
     ! time_dep_parameters
  namelist/input_variables/ dt, max_dt, min_dt, simulation_time, max_time, max_time_iter, shampine_test, &
                               shampine_tol, stepsize_increase_frequency, lee_tol_failing, lee_tol_halving, &
                               lee_tol_decrease, lee_tol_increase
     ! refinement_parameters
  namelist/input_variables/ allsolid_test, ctore, ctode
     ! generic_parameters
  namelist/input_variables/ local_adaptation, numglobalref, grow, verbose, output_rate, multigrid_on, C2F_switch, paraview_output
     ! solution_parameters
  namelist/input_variables/ lambda, delta, epsilon, ke, le, mcinf, nucleate_x, nucleate_y, nucleate_z, &
                            nuc_radius, anti_trapping_mod, solute

  
  inquire(file="CASE.inp", exist=there)
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)

  if (there) then
     if (pid.eq.0.and.verbose.ge.1) print *,'File CASE.inp exists'
     open(file="CASE.inp", unit=30,status="unknown")
     read(30,nml=input_variables)
     if ((pid.eq.0.and.verbose.ge.2).or.verbose.ge.4) then
        write (6,nml=input_variables)
     endif
     close(30)
  else
     if (pid.eq.0.and.verbose.ge.1) print *,'No input file (CASE.inp) exists'
  endif

end subroutine input_parameter_file_reading
