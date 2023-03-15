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
  use solution_parameters
  use paramesh_dimensions
  use time_dep_parameters

  implicit none
  double precision :: qq = 0.0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set parameters that affect global generic parameters
  tstep_increase_factor = 1.1
  tstep_decrease_factor = 0.8
  multigrid_on=1
  weight = 1.0
  smooth_count_max = 3! if explicit mp_smooth set to 1
  solve_count_max = 3
  low_vcycle_count = 3 
  high_vcycle_count =5 
  max_v_cycle_count = 12
  defect_tol_min = 1e-10 !if explicit set to ~1e-1
  defect_tol_min_f = 1.0
  defect_tol_max = 1e1
  defect_too_big = 1e1
  min_dt = 1.0e-16
  output_rate = 1000
  max_time_iter =5000000
  multiphase = .false.

  v_mesh = 0.0
  verbose=1            ! How much terminal output
  local_adaptation = 1 ! Are we using local_adaptation?  0=no 1=yes 2=user supplied
  numglobalref=1       ! How many levels of global refinement are needed initially
  grow = 0             ! Are we letting the domain grow as adaptation spreads?  1=yes
  C2F_switch = 0       ! Are we performing coarse to fine opration? 1=yes

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  dt = 1.0e-6

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Physical phase field case - variables declared in pf_modules.f90
  eps2(1) = 1.*(1.-qq)
  eps2(2) = 1.*qq
  eps2(3) = 1.*(0.77+qq)
  eps22(1,2) = 1.0 
  eps22(1,3) = 1.77
  eps22(3,2) = 0.77
  eps22(2,1) = 1.0 
  eps22(3,1) = 1.77
  eps22(2,3) = 0.77
  scale_factor=1.0
  nuc_radius = 7d0
  lambda=4.0
  beta=0.0045  !beta=0.006
  alpha=0.0
  T_scale = 1.0
  delta = -0.02192309379
  delta = 0.12
  delta = delta*T_scale
  epsilon = 0.02
  ke = 0.30
  le = 40.0                  ! Only if thermal-solutal coupled
  mcinf = 0.08	             ! Only if thermal-solutal coupled
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Derived assignments - do not change
  epsilon_tilde = 4.0*epsilon/(1.0-3.0*epsilon)  
  ! See page 26 of Jan's thesis
  A_0 = 1.0-3.0*epsilon                           ! See page 26 of Jan's thesis


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  
  ! Set implicit to none
  implicit none
  
  ! include file required for mpi library.
  include 'mpif.h'
  
  integer, intent(in) :: pid, noprocs
  ! VARIABLES NEEDED LOCALLY GO HERE
  integer :: ierr
  integer nBx, nBy, nBz, i, j, k
  double precision dBx, dBy, dBz
!  double precision :: grid_min=-3200.0,grid_max=3200.0      ! grid 11
!  double precision :: grid_min=-1600.0,grid_max=1600.0      ! grid 10
!  double precision :: grid_min=-1200.0,grid_max=1200.0      ! grid 9
!  double precision :: grid_min=-800.0,grid_max=800.0      ! grid 9
!  double precision :: grid_min=-400.0,grid_max=400.0      ! grid 8
!  double precision :: grid_min=-200.0,grid_max=200.0     ! grid 7
!  double precision :: grid_min= 0.,grid_max= 37.05242907  ! grid 6   !Recommended coarse block resolution
  double precision :: grid_min= 0.,grid_max= 1.25  ! grid 6   !Recommended coarse block resolution
!  double precision :: grid_min=-50.0,grid_max=40.0       ! grid 5
!  double precision :: grid_min=-25.0,grid_max=25.0       ! grid 4
!  double precision :: grid_min=-12.50,grid_max=12.50       ! grid 3

  write(*,*)total_vars, nvar, nunkvbles

  eighthdmn = 1
  grid_max=grid_max/scale_factor
  if (eighthdmn.eq.1) grid_min = 0.0
  mg_min_lvl = 1

  nBx = 320
  nBy = 1
  nBz = 1


  if (ndim.eq.2) nBz = 1

  dBx = grid_max-grid_min
  dBy = dBx
  dBz = dBx

  ! set a limit on the refinement level
  lrefine_min = mg_min_lvl

!  lrefine_max = 12
!  lrefine_max = 11
!  lrefine_max = 10
!  lrefine_max = 9
!  lrefine_max = 8
  lrefine_max = 2
!  lrefine_max = 6                       ! For G6 this is dx=0.39
!  lrefine_max = 5                       ! For G6 this is dx=0.78
!  lrefine_max = 4
!  lrefine_max = 3
!  lrefine_max = 2

   numglobalref=1

   !!--tim


   if(C2F_switch.eq.1)then
     !!--tim
     !! user edit bit
     C2F_start_level = lrefine_max-1
     !!--tim
     !! leave is along
     C2F_desired_level = lrefine_max
     if(C2F_desired_level.lt.C2F_start_level)then
       if(pid.eq.0)write(*,*) "ERROR C2F is on and desired level is less then current level"
       stop
     endif
   endif

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
      max_dt = 1d-2 /scale_factor**2

  endif
  boundary_index(1) = bc_cond_sym 
  boundary_index(2) = bc_cond_sym 
  boundary_index(3) = bc_cond_sym  
  boundary_index(4) = bc_cond_sym  

  ! ensure all processors have min_dx and max_dt set
  call MPI_Bcast(min_dx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(max_dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
end subroutine app_domain_setup

