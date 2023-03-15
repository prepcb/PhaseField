!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  input_parameter_file_reading 
!!!!  * Reads in parameters from the file 'CASE.inp' that may be set at runtime
!!!!  * To add new ones just include any module variable in the set defined by 
!!!!    "namelist/input_variables/"
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Latest version: Chris Goodyer, September 2012
!!!!  Latest version: Peter Bollada, July 2020
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
  namelist/input_variables/g_lambda,g_beta,g_c0,g_t0,heat_c_term,Grad_DW_term,heat_phi_term,cross_term,g_Le,g_linear,g_weight,g_diagonal,multigrid_on,g_triangular,g_npoints,g_Wheeler,g_Dch, g_infinity, g_mu_infinity,g_d_facets,g_amax,g_pot_h,g_shape
  namelist/input_variables/  nuc_radius,  delta,  epsilon,thermal,g_max_time,smooth_count_max,g_max_dt,g_mu,g_1D,g_vel,g_ratio,paraview,g_obstacle, g_gamma,g_drive,g_Radius_min
  namelist/input_variables/g_sym,defect_tol_min,g_max_lvl,g_output_rate,g_grid,g_nbx,g_nby,g_nbz,AlNi,hex_facet,g_alpha,g_xx1,g_xx0,g_barrier
  namelist/input_variables/g_kurv,g_curve,g_approx,PCB_order,g_needle,g_x0,g_y0,g_z0,aspect_ratio,kinetic_force,g_abs,g_abs_linear
  namelist/input_variables/ g_kinetic,g_k_mult,g_c_max,g_c_min,g_q0,g_skew_potential,g_power_skew,g_power,g_at_factor,g_pure,g_antitrapping,g_mesh_dependence,g_hex
  namelist/input_variables/uniform_mesh,g_s_skew,g_s,g_AMM,g_plapp,g_ellipse,g_equal_delta,g_r_factor,solve_count_max
  namelist/input_variables/g_accel_factor,g_accel_factor_liq,g_slope_factor,g_no_log_term,convection_predictor, n_sharpen,g_multiple,g_grid_factor,g_time_end
  namelist/input_variables/verbose,smooth_count_max,low_vcycle_count,high_vcycle_count,max_v_cycle_count,g_toy_model,New_potential, g_velocity,g_parallel_scaling
  namelist/input_variables/g_mol_n,g_mol_d, g_mol_epsilon
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
