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
  namelist/input_variables/g_lambda,g_beta,g_c0,g_t0,heat_c_term,Grad_DW_term,heat_phi_term,cross_term,g_Le,g_linear,g_weight,g_diagonal,multigrid_on,g_triangular,g_npoints,new_potential,g_Isotropic,g_LUT,g_Quad,g_right_limit,g_left_limit,g_gamma,g_cinfty,g_axis_symmetric,g_slider
  namelist/input_variables/  nuc_radius,  delta,  epsilon,thermal,g_max_time,smooth_count_max,g_max_dt,g_mu,g_1D,g_vel,g_ratio,g_mu_kinetic,g_obstacle,g_Dch,g_Halley
  namelist/input_variables/g_sym,defect_tol_min,g_max_lvl,g_output_rate,g_grid,g_nbx,g_nby,g_nbz,AlNi,hex_facet,g_alpha,g_kurv,g_curve,g_approx,PCB_order,g_needle,g_x0,g_y0
  namelist/input_variables/ g_kinetic,g_k_mult,g_c_max,g_c_min,g_q0,g_skew_potential,g_power_skew,g_power,g_at_factor,g_pure,uniform_mesh,g_s_skew,g_s,g_AMM,g_plapp,g_ellipse,g_equal_delta,g_r_factor,uniformstart,solve_count_max,verbose,smooth_count_max,low_vcycle_count,high_vcycle_count,max_v_cycle_count,g_drive
  
  inquire(file="CASE.inp", exist=there)
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)

  if (there) then
     if (pid.eq.0.and.verbose.ge.1) print *,'File CASE.inp exists'
     open(file="CASE.inp", unit=30,status="unknown")
     write(*,*)"opened"
     read(30,nml=input_variables)
     write(*,*)"read"
     close(30)
  else
     if (pid.eq.0.and.verbose.ge.1) print *,'No input file (CASE.inp) exists'
  endif

end subroutine input_parameter_file_reading
