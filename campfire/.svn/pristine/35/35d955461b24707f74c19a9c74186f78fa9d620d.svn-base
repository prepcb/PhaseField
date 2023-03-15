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
  double precision :: qq = 0.0
  integer :: i,j
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set parameters that affect global generic parameters
  tstep_increase_factor = 1.1
  tstep_decrease_factor = 0.8
  multigrid_on=1
  weight = 1.0
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
  multiphase = .true.
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

  eps22  = 0d0
  eps22(1,2) = 1.00 
  eps22(1,3) = 1.77
  eps22(2,3) = 1.0 !0.2
  eps22(2,1) = eps22(1,2) 
  eps22(3,1) = eps22(1,3)
  eps22(3,2) = eps22(2,3)
  

  Mobility(1,3)=1d0
  Mobility(3,1)=1d0
  Mobility(1,2)=1d0
  Mobility(2,1)=1d0
  Mobility(2,3)=1d0
  Mobility(3,2)=1d0





  scale_factor=1.0
  nuc_radius = 5/scale_factor
  lambda=4.0
  T_scale = 1.0
!  delta = -0.02192309379
  delta=-0.022
  delta = delta*T_scale
  epsilon = 0d0!4d-2!(1d-6)/4.
  width=20d0
  
  ! Do we really need all these values?
   total_vars=nvar/nunkvbles

  N_phases = total_vars-2


  ! How big and where should the initial seed be?
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Derived assignments - do not change
  epsilon_tilde = 4.0*epsilon/(1.0-3.0*epsilon)  
  ! See page 26 of Jan's thesis
  A_0 = 1.0-3.0*epsilon                           ! See page 26 of Jan's thesis


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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

rrblock =reshape((/0,5,10,15,0,5,0,10,5,15,10,15,&
                   0,5,10,15,5,0,10,0,15,5,15,10,&
                   1,2,3,4,1,2,1,3,2,4,3,4, &
                   2,2,2,2,3,1,2,2,2,2,3,1,&
                   2,2,2,2,2,2,3,1,3,1,2,2/),(/12,5/))

rrarray = reshape((/0,1,0,1,&
                    0,0,1,1/),(/4,2/))
! read as 
!        x  | 0    1    0    1|
!        y  | 0    0    1    1|

rrrblock = reshape((/&
     0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75,&
     0, 5, 0, 10, 5, 15, 10, 15, 20, 25, 20, 30, 25, 35, 30, 35,&
     40, 45, 40, 50, 45, 55, 50, 55, 60, 65, 60, 70, 65, 75, 70, 75,&
     5, 20, 15, 30, 10, 40, 15, 45, 45, 60, 55, 70, 30, 60, 35, 65,&
     0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75,&
     5, 0, 10, 0, 15, 5, 15, 10, 25, 20, 30, 20, 35, 25, 35, 30,&
     45, 40, 50, 40, 55, 45, 55, 50, 65, 60, 70, 60, 75, 65, 75, 70,&
     20, 5, 30, 15, 40, 10, 45, 15, 60, 45, 70, 55, 60, 30, 65, 35,&
     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,&
     1, 2, 1, 3, 2, 4, 3, 4, 5, 6, 5, 7, 6, 8, 7, 8,&
     9, 10, 9, 11, 10, 12, 11, 12, 13, 14, 13, 15, 14, 16, 15, 16,&
     2, 5, 4, 7, 3, 9, 4, 10, 10, 13, 12, 15, 7, 13, 8, 14,&
     2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,&
     3, 1, 2, 2, 2, 2, 3, 1, 3, 1, 2, 2, 2, 2, 3, 1,&
     3, 1, 2, 2, 2, 2, 3, 1, 3, 1, 2, 2, 2, 2, 3, 1,&
     3, 1, 3, 1, 2, 2, 2, 2, 3, 1, 3, 1, 2, 2, 2, 2,&
     2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,&
     2, 2, 3, 1, 3, 1, 2, 2, 2, 2, 3, 1, 3, 1, 2, 2,&
     2, 2, 3, 1, 3, 1, 2, 2, 2, 2, 3, 1, 3, 1, 2, 2,&
     2, 2, 2, 2, 3, 1, 3, 1, 2, 2, 2, 2, 3, 1, 3, 1/),(/64,5/))

rrrarray= reshape((/&
     0, 1, 0, 1, 2, 3, 2, 3,&
     0, 1, 0, 1, 2, 3, 2, 3,&
     0, 0, 1, 1, 0, 0, 1, 1,&
     2, 2, 3, 3, 2, 2, 3, 3/),(/16,2/))



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

!  double precision :: grid_min= 0.,grid_max=37.0524290
  double precision :: grid_min= 0.,grid_max

  grid_max=width*4.1d0
  eighthdmn = 0
  grid_max=grid_max/scale_factor
  if (eighthdmn.eq.1) grid_min = 0.0
  mg_min_lvl = 1

  nBx = 5
  nBy = 1
  nBz = 1

  if (ndim.eq.2) nBz = 1

  dBx = grid_max-grid_min
  dBy = dBx
  dBz = dBx

  max_x = dBx*nBx
  ! set a limit on the refinement level
  lrefine_min = mg_min_lvl


  lrefine_max = 6


   numglobalref=1

   !!--tim
   C2F_switch = 0
   
   if(C2F_switch.eq.1)then
     !!--tim
     !! user edit bit
     C2F_start_level = 5
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
      max_dt = 1d0 /scale_factor**2
   
      
  endif
  boundary_index(1) = bc_cond_sym 
  boundary_index(2) = bc_cond_sym 
  boundary_index(3) = bc_cond_sym  
  boundary_index(4) = bc_cond_sym  

  ! ensure all processors have min_dx and max_dt set
  call MPI_Bcast(min_dx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(max_dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
end subroutine app_domain_setup

