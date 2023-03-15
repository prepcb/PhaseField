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
  smooth_count_max = 1! if explicit mp_smooth set to 1
  solve_count_max = 1
  low_vcycle_count = 3 
  high_vcycle_count =5 
  max_v_cycle_count = 5
  defect_tol_min = 1e-7 !if explicit set to ~1e-1
!  defect_tol_max = 1e-1
  defect_too_big = 1e1
  min_dt = 1.0e-16
  output_rate = 100
  max_time_iter = 400000
  verbose=1            ! How much terminal output
  numglobalref=1       ! How many levels of global refinement are needed initially

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  dt = 1.0e-6

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Physical phase field case - variables declared in pf_modules.f90


! from mp_modules

  g_delta(1)=5d-10
  g_delta(2)=g_delta(1)*(g_T(2)/g_T(1))*(g_sigma(1)/g_sigma(2))
  g_M(1)=g_T(1)*g_beta(1)/(sqrt(72.)*g_L(1)*g_delta(1))
  g_M(2)=g_T(2)**2/g_T(1)*g_beta(2)/(sqrt(72.)*g_L(2)*g_delta(2))
  g_W(1)=3.*g_sigma(1)/(sqrt(2.)*g_delta(1))
  g_W(2)=3.*g_sigma(2)*g_T(1)/(sqrt(2.)*g_delta(2)*g_T(2))
  g_epsilon(1) = sqrt(72.)*g_sigma(1)*g_delta(1)
  g_epsilon(2) = sqrt(72.)*g_sigma(2)*g_delta(2)*g_T(1)/g_T(2)
  g_tau=1./(g_M(1)*g_W(1))
  g_lengthSQ=g_epsilon(1)/g_W(1)
  g_Dch=g_lengthSQ/g_tau
  g_D_tilde=1d4*g_D(1)/g_Dch
  g_c0 = 0.6
  g_T0 = (g_L(1)*(1-g_c0)+g_L(2)*g_c0)/(g_L(1)*(1-g_c0)/g_T(1)+g_L(2)*g_c0/g_T(2))
  

  hsa(2, 1) = -g_C(1)*g_T(1) 
  hsa(2, 2) = g_C(1)*log(g_T(1))
  hsa(2, 3) = -g_C(1)
  
  hsb(2, 1) = -g_C(2)*g_T(2) 
  hsb(2, 2) = g_C(2)*log(g_T(2))
  hsb(2, 3) = -g_C(2)
  

  hsa(1, 1) = g_L(1) + hsa(2, 1)
  hsa(1, 2) = -g_L(1)/g_T(1) + hsa(2, 2)
  hsa(1, 3) = hsa(2, 3)
  
  hsb(1, 1) = g_L(2) + hsb(2, 1)
  hsb(1, 2) = -g_L(2)/g_T(2) + hsb(2, 2)
  hsb(1, 3) = hsb(2, 3)



  nuc_radius = 5.
  delta=-80.0/g_T0
  epsilon = 0d0!4d-2!(1d-6)/4.
  heat_c_term = .true.
  heat_phi_term = .true.
  Grad_DW_term=.true.
  ! Do we really need all these values?
   total_vars=nvar/nunkvbles

  N_phases = total_vars-2
  M_unk = total_vars*6
  N_unk = total_vars
  NN_unk= N_unk*4
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Derived assignments - do not change
  epsilon_tilde = 4.0*epsilon/(1.0-3.0*epsilon)  
  ! See page 26 of Jan's thesis
  A_0 = 1.0-3.0*epsilon                           ! See page 26 of Jan's thesis


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
  double precision :: grid_min= 0.,grid_max=200.


  eighthdmn = 0

  if (eighthdmn.eq.1) grid_min = 0.0
  mg_min_lvl = 1

  nBx = 2
  nBy = 2
  nBz = 1

  if (ndim.eq.2) nBz = 1

  dBx = grid_max-grid_min
  dBy = dBx
  dBz = dBx

  max_x = dBx*nBx
  ! set a limit on the refinement level
  lrefine_min = mg_min_lvl


  lrefine_max = 8


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
!      max_dt = 0.2*min_dx
      max_dt = 1d0 
   
      
  endif
  boundary_index(1) = bc_cond_sym 
  boundary_index(2) = bc_cond_sym 
  boundary_index(3) = bc_cond_sym  
  boundary_index(4) = bc_cond_sym  

  ! ensure all processors have min_dx and max_dt set
  call MPI_Bcast(min_dx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(max_dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
end subroutine app_domain_setup

