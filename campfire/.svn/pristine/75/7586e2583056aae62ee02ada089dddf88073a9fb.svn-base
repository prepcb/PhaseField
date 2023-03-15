!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_parameter_set                                            REQUIRED
!!!!  * Called at start of program execution.
!!!!  * Paramesh setup needed first to get nvar set
!!!!  * Sets case specific parameters
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!  
!!!!  app_domain_setup                                             REQUIRED
!!!!  * Sets up the bounding boxes and boundary conditions for the initial 
!!!!    coarsest block distribution
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Latest version: Chris Goodyer, May 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!
!  Control file for paramesh time dependant non-linear multigrid solver
!  Written by J. Green
!  Based upon: http://www.physics.drexel.edu/~olson/paramesh-doc/Users_manual/amr_users_guide.html#main_p
!

subroutine app_parameter_set
  ! include modules for governing solution parameters
  use multigrid_parameters
  use generic_parameters
  ! include file to define physical quantities of the model
  use solution_parameters
  ! include paramesh data for nvar
  use paramesh_dimensions
  use time_dep_parameters

  implicit none


  integer :: i,j
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set parameters that affect global generic parameters
  tstep_increase_factor = 1.1
  tstep_decrease_factor = 0.8
  multigrid_on=1
  weight =1.0
  smooth_count_max = 3! if explicit mp_smooth set to 1
  solve_count_max = 3

  defect_tol_min = 1e-7 !if explicit set to ~1e-1
!  defect_tol_max = 1e-1
  defect_too_big = 1e1
  min_dt = 1.0e-16
!  output_rate = 100
  max_time_iter = 400000
  verbose=1            ! How much terminal output
  numglobalref=1       ! How many levels of global refinement are needed initially

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  dt = 1.0e-6

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Physical phase field case - variables declared in pf_modules.f90

!  namelist/input_variables/g_lambda,g_beta,g_c0,g_t0,heat_c_term,Grad_DW_term,heat_phi_term,cross_term,g_Le
!
! Lead-Tin   Pb-Sn   Pb rich solid  has the lower values of c and Sn the higher.
! phi(1) liquid
! phi(2) Pb
! phi(3) Sn 

! from mp_modules

  g_lambda=4.0
  g_beta = 0.0
  g_c0 = 0.5
  g_T0 = (g_L(1)*(1-g_c0)+g_L(2)*g_c0)/(g_L(1)*(1-g_c0)/g_T(1)+g_L(2)*g_c0/g_T(2))
  heat_c_term = .true.
  Grad_DW_term=.true. 
  heat_phi_term = .true.
  cross_term = .false.
  g_Le=100.0
  nuc_radius = 20.0!2.2d-8/g_delta(1)
  delta=(175.+273.-g_T0)/g_T0
  epsilon = 0.02!(1d-6)/4.
  g_mu(1)=0.04
  g_mu(2)=0.04

  call input_parameter_file_reading() ! replace the 8 above if CASE.inp present
  low_vcycle_count =    high_vcycle_count - 2
  max_v_cycle_count =   high_vcycle_count
  output_rate = g_output_rate
  g_delta(1)=g_lambda*g_d0
  g_sigma(1)=0.033 !Pb
  g_sigma(2)=0.059 !Sn
  g_T(1) = 600.
  g_T(2) = 505.
  g_L(1) = 11340*23000 !density*Latent heat per kg
  g_L(2) = 7280*59000
  g_MolPerVol(1)=54730.
  g_MolPerVol(2)=61505.
  g_delta(2)=g_delta(1)*(g_T(2)/g_T(1))*(g_sigma(1)/g_sigma(2))
  g_M(1)=g_T(1)*g_mu(1)/(sqrt(72.)*g_L(1)*g_delta(1))
  g_M(2)=g_T(2)**2/g_T(1)*g_mu(2)/(sqrt(72.)*g_L(2)*g_delta(2))
  g_W(1)=3.*g_sigma(1)/(sqrt(2.)*g_delta(1))
  g_W(2)=3.*g_sigma(2)*g_T(1)/(sqrt(2.)*g_delta(2)*g_T(2))
  g_epsilon(1) = sqrt(72.)*g_sigma(1)*g_delta(1)
  g_epsilon(2) = sqrt(72.)*g_sigma(2)*g_delta(2)*g_T(1)/g_T(2)
  g_tau=1./(g_M(1)*g_W(1)*g_lambda**2)
  g_lengthSQ=g_epsilon(1)/(g_W(1)*g_lambda**2)
  g_Dch=g_lengthSQ/g_tau
  if(zero_start)then
     g_D(2)=g_D(1)
     g_D(3)=g_D(1)
     g_eps22= 1.0
  else


     g_D(2)=g_D(1)*1d-4
     g_D(3)=g_D(1)*1d-4
     g_eps22(1,2)= 1.0  !sigma_12 = 33.3 ergs per sq cm for Pb-Liq
     g_eps22(3,1)= g_epsilon(2)/g_epsilon(1)  !sigma_12 = 59.0 ergs per sq cm for Sn-Liq
     g_eps22(2,3)= 1.0!0.5*(g_eps22(1,2)+g_eps22(3,1))
  endif
  g_D_tilde=g_Le*g_D(1)/g_Dch







  hsa(1, 1) = -2977.961 ! J/mol 
  hsa(1, 2) = 93.949561 
  hsa(1, 3) = -24.5242231 
  hsa(1, 4) = -0.365895e-2 
  hsa(1, 5) = -0.24395e-6 
  hsa(1, 6) = 0. 
  hsa(1, 7) = -0.6019e-18 
  hsa(1, 8) = 0. 
  hsb(1, 1) = 1247.957 +000.
  hsb(1, 2) = 51.355548 
  hsb(1, 3) = -15.961 
  hsb(1, 4) = -0.188702e-1 
  hsb(1, 5) = 0.3121167e-5 
  hsb(1, 6) = -61960.0 
  hsb(1, 7) = 0.147031e-17 
  hsb(1, 8) = 0. 
  hsrk(1, 1) = 6204.5 
  hsrk(1, 2) = -.67981 
  hsrk(1, 3) = 791.7 
  hsrk(1, 4) = -1.5219 
  hsrk(1, 5) = 0. 
  hsrk(1, 6) = 0. 
  hsrk(1, 7) = 0. 
  hsrk(1, 8) = 0. 
  hsa(2, 1) = -7650.085 
  hsa(2, 2) = 101.700244 
  hsa(2, 3) = -24.5242231 
  hsa(2, 4) = -0.365895e-2 
  hsa(2, 5) = -0.24395e-6 
  hsa(2, 6) = 0. 
  hsa(2, 7) = 0. 
  hsa(2, 8) = 0. 
  hsb(2, 1) = -345.135 
  hsb(2, 2) = 56.983315 
  hsb(2, 3) = -15.961 
  hsb(2, 4) = -0.188702e-1 
  hsb(2, 5) = 0.3121167e-5 
  hsb(2, 6) = -61960.0 
  hsb(2, 7) = 0. 
  hsb(2, 8) = 0. 
  hsrk(2, 1) = 7145.3 
  hsrk(2, 2) = -2.30237 
  hsrk(2, 3) = 0. 
  hsrk(2, 4) = 0. 
  hsrk(2, 5) = 0. 
  hsrk(2, 6) = 0. 
  hsrk(2, 7) = 0. 
  hsrk(2, 8) = 0. 
  hsa(3, 1) = -7161.085 
  hsa(3, 2) = 105.220244 
  hsa(3, 3) = -24.5242231 
  hsa(3, 4) = -0.365895e-2 
  hsa(3, 5) = -0.24395e-6 
  hsa(3, 6) = 0. 
  hsa(3, 7) = 0. 
  hsa(3, 8) = 0. 
  hsb(3, 1) = -5855.135 
  hsb(3, 2) = 65.443315 
  hsb(3, 3) = -15.961 
  hsb(3, 4) = -0.188702e-1 
  hsb(3, 5) = 0.3121167e-5 
  hsb(3, 6) = -61960.0 
  hsb(3, 7) = 0. 
  hsb(3, 8) = 0. 
  hsrk(3, 1) = 19700.0 
  hsrk(3, 2) = -15.89 
  hsrk(3, 3) = 0. 
  hsrk(3, 4) = 0. 
  hsrk(3, 5) = 0. 
  hsrk(3, 6) = 0. 
  hsrk(3, 7) = 0. 
  hsrk(3, 8) = 0. 









  ! Do we really need all these values?
  N_unk=nvar/nunkvbles -1 !some spare capacity for output

  N_phases = N_unk-2

  M_unk = N_unk*6

  total_vars=N_unk
  NN_unk= N_unk*4  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Derived assignments - do not change
  epsilon_tilde = 4.0*epsilon/(1.0-3.0*epsilon)  
  ! See page 26 of Jan's thesis
  A_0 = 1.0-3.0*epsilon                           ! See page 26 of Jan's thesis


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   anything you do to the following will almost certainly break it
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rrblock =reshape((/0,N_unk,2*N_unk,3*N_unk, 0,N_unk, 0,2*N_unk, N_unk,3*N_unk, 2*N_unk,3*N_unk,&
                   0,N_unk,2*N_unk,3*N_unk, N_unk,0, 2*N_unk,0, 3*N_unk,N_unk, 3*N_unk,2*N_unk,&
                   1,2,3,4,1,2,1,3,2,4,3,4, &
                   2,2,2,2,3,1,2,2,2,2,3,1,&
                   2,2,2,2,2,2,3,1,3,1,2,2/),(/12,5/))



rrarray = reshape((/0,1,0,1,&
                    0,0,1,1/),(/4,2/))
! read as 
!        x  | 0    1    0    1|
!        y  | 0    0    1    1|





end subroutine app_parameter_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_domain_setup(pid, noprocs)
  ! include file to define physical qualities of the model and mesh
  use paramesh_dimensions
  use physicaldata
  ! include workspace so interp_mask_work can be set
  use workspace
  ! include file defining the tree
  use tree
  use solution_parameters
  use multigrid_parameters
  use time_dep_parameters
  use generic_parameters

  ! IMPORTANT: Expose the PARAMESH interface blocks to your code
  use paramesh_interfaces
  use paramesh_mpi_interfaces
  use solution_parameters
  
  ! Set implicit to none
  implicit none
  
  ! include file required for mpi library.
  include 'mpif.h'
  
  integer, intent(in) :: pid, noprocs
  ! VARIABLES NEEDED LOCALLY GO HERE
  integer :: ierr
  integer nBx, nBy, nBz, i, j, k
  double precision dBx, dBy, dBz
  double precision :: grid_min= 0.,grid_max=12.5

  


  eighthdmn = 0

  if (eighthdmn.eq.1) grid_min = 0.0


  nBx = 20
  nBy = 20
  nBz = 1

  grid_max=g_grid
  nbx=g_nbx
  nby=g_nby
  nbz=g_nbz

  if (ndim.eq.2) nBz = 1

  dBx = grid_max-grid_min
  dBy = dBx
  dBz = dBx

  max_x = dBx*nBx
  ! set a limit on the refinement level
  lrefine_min = 1


  lrefine_max = g_max_lvl


   numglobalref=1






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!          Leave untouched from here down for phase field code                                       !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


! CEG Initialise my startlist for the blocklist for that grid level
  allocate(block_starts(lrefine_max+1))
  
  ! set up a single block covering the whole cube domain
  lnblocks = 0
  block_starts(1) = 1
  if (pid.eq.0) then
     do k = 1, nBz
        do j = 1, nBy
           do i = 1, nBx
              lnblocks = lnblocks+1
              coord(1,lnblocks) = nBx*grid_min + (i-0.5)*dBx
              coord(2,lnblocks) = nBy*grid_min + (j-0.5)*dBy

              bnd_box(1,1,lnblocks) = nBx*grid_min + (i-1.0)*dBx
              bnd_box(2,1,lnblocks) = nBx*grid_min + i*dBx
!              write(*,*)bnd_box(1,1,lnblocks),bnd_box(2,1,lnblocks),nBx*dBx
              bnd_box(1,2,lnblocks) = nBy*grid_min + (j-1.0)*dBy
              bnd_box(2,2,lnblocks) = nBy*grid_min + j*dBy
               
!               if (ndim.eq.3) then
!                  coord(3,lnblocks) = nBz*grid_min+(k-0.5)*dBz! Fix z coord at middle of layer
!                  bnd_box(1,3,lnblocks) = nBz*grid_min+(k-1.0)*dBy
!                  bnd_box(2,3,lnblocks) = nBz*grid_min+k*dBy
!               endif

              bsize(:,lnblocks) = bnd_box(2,:,lnblocks) - bnd_box(1,:,lnblocks)
              nodetype(lnblocks) = 1
              lrefine(lnblocks) = 1
     
              ! boundary conditions need to represent whether far field or symmetry
              
            neigh(1,1,lnblocks) = lnblocks-1
            neigh(1,2,lnblocks) = lnblocks+1
            neigh(1,3,lnblocks) = lnblocks-nBx
            neigh(1,4,lnblocks) = lnblocks+nBx
            if (i.eq.1)  neigh(:,1,lnblocks)=bc_cond_sym
            if (i.eq.nBx)neigh(:,2,lnblocks)=bc_cond_sym
            if (j.eq.1)neigh(:,3,lnblocks)=bc_cond_sym
            if (j.eq.nBy)neigh(:,4,lnblocks)=bc_cond_sym
            neigh(2,:,lnblocks) = 0
              refine(lnblocks)=.true.
           end do
        end do
     end do
     min_dx = (bnd_box(2,1,1)-bnd_box(1,1,1))/(2**(lrefine_max-1)*nxb)
     if (pid.eq.0.and.verbose.ge.2) print *,"Min available dx=", min_dx
      max_dt = g_max_dt

   
      
  endif
  boundary_index(1) = bc_cond_sym 
  boundary_index(2) = bc_cond_sym 
  boundary_index(3) = bc_cond_sym  
  boundary_index(4) = bc_cond_sym  

  ! ensure all processors have min_dx and max_dt set
  call MPI_Bcast(min_dx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(max_dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
end subroutine app_domain_setup

