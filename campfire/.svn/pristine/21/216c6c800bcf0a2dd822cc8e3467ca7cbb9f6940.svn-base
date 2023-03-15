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
  low_vcycle_count = 3 
  high_vcycle_count =5 
  max_v_cycle_count = 5
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

  output_rate = g_output_rate

  N_unk=nvar/nunkvbles -1 !some spare capacity for output

  N_phases = 1
  M_unk = N_unk*6

  total_vars=N_unk
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


 d_10_10(1,1) = -0.1273035714D2
 d_10_10(1,2) = 0.4049999961D2
 d_10_10(1,3) = -0.8100000161D2
 d_10_10(1,4) = 0.1259999836D3
 d_10_10(1,5) = -0.1417499919D3
 d_10_10(1,6) = 0.1134000069D3
 d_10_10(1,7) = -0.6299999614D2
 d_10_10(1,8) = 0.2314285451D2
 d_10_10(1,9) = -0.5062499749D1
 d_10_10(1,10) = 0.5000000571D0
 d_10_10(2,1) = -0.499999996D0
 d_10_10(2,2) = -0.7730357237D1
 d_10_10(2,3) = 0.1799999982D2
 d_10_10(2,4) = -0.2100000443D2
 d_10_10(2,5) = 0.2100000238D2
 d_10_10(2,6) = -0.1574999903D2
 d_10_10(2,7) = 0.8400001523D1
 d_10_10(2,8) = -0.3000000548D1
 d_10_10(2,9) = 0.642857171D0
 d_10_10(2,10) = -0.6249998964D-1
 d_10_10(3,1) = 0.6250000164D-1
 d_10_10(3,2) = -0.1125000024D1
 d_10_10(3,3) = -0.4917857063D1
 d_10_10(3,4) = 0.1049999904D2
 d_10_10(3,5) = -0.7874999362D1
 d_10_10(3,6) = 0.524999993D1
 d_10_10(3,7) = -0.2624999467D1
 d_10_10(3,8) = 0.8999998774D0
 d_10_10(3,9) = -0.1874999996D0
 d_10_10(3,10) = 0.1785714413D-1
 d_10_10(4,1) = -0.1785714261D-1
 d_10_10(4,2) = 0.2410714255D0
 d_10_10(4,3) = -0.192857141D1
 d_10_10(4,4) = -0.2775000107D1
 d_10_10(4,5) = 0.6750000098D1
 d_10_10(4,6) = -0.3375000067D1
 d_10_10(4,7) = 0.1500000079D1
 d_10_10(4,8) = -0.4821428676D0
 d_10_10(4,9) = 0.9642857002D-1
 d_10_10(4,10) = -0.8928571316D-2
 d_10_10(5,1) = 0.8928571426D-2
 d_10_10(5,2) = -0.1071428572D0
 d_10_10(5,3) = 0.6428571413D0
 d_10_10(5,4) = -0.3000000004D1
 d_10_10(5,5) = -0.899999995D0
 d_10_10(5,6) = 0.4499999998D1
 d_10_10(5,7) = -0.1499999994D1
 d_10_10(5,8) = 0.4285714244D0
 d_10_10(5,9) = -0.8035714292D-1
 d_10_10(5,10) = 0.7142857165D-2
 d_10_10(6,1) = -0.7142857139D-2
 d_10_10(6,2) = 0.803571427D-1
 d_10_10(6,3) = -0.4285714275D0
 d_10_10(6,4) = 0.1499999998D1
 d_10_10(6,5) = -0.4499999999D1
 d_10_10(6,6) = 0.900000004D0
 d_10_10(6,7) = 0.2999999993D1
 d_10_10(6,8) = -0.6428571381D0
 d_10_10(6,9) = 0.1071428569D0
 d_10_10(6,10) = -0.892857144D-2
 d_10_10(7,1) = 0.8928571397D-2
 d_10_10(7,2) = -0.964285718D-1
 d_10_10(7,3) = 0.4821428573D0
 d_10_10(7,4) = -0.1500000008D1
 d_10_10(7,5) = 0.3375000024D1
 d_10_10(7,6) = -0.675000002D1
 d_10_10(7,7) = 0.2775000019D1
 d_10_10(7,8) = 0.1928571414D1
 d_10_10(7,9) = -0.2410714282D0
 d_10_10(7,10) = 0.1785714292D-1
 d_10_10(8,1) = -0.178571434D-1
 d_10_10(8,2) = 0.1874999976D0
 d_10_10(8,3) = -0.9000000167D0
 d_10_10(8,4) = 0.2625000019D1
 d_10_10(8,5) = -0.5249999942D1
 d_10_10(8,6) = 0.787499994D1
 d_10_10(8,7) = -0.1050000002D2
 d_10_10(8,8) = 0.4917857131D1
 d_10_10(8,9) = 0.1125000005D1
 d_10_10(8,10) = -0.6250000005D-1
 d_10_10(9,1) = 0.6249999837D-1
 d_10_10(9,2) = -0.6428571438D0
 d_10_10(9,3) = 0.2999999953D1
 d_10_10(9,4) = -0.8399999911D1
 d_10_10(9,5) = 0.1575000033D2
 d_10_10(9,6) = -0.210000004D2
 d_10_10(9,7) = 0.2099999997D2
 d_10_10(9,8) = -0.1799999987D2
 d_10_10(9,9) = 0.7730357135D1
 d_10_10(9,10) = 0.5000000024D0
 d_10_10(10,1) = -0.5000000099D0
 d_10_10(10,2) = 0.5062499968D1
 d_10_10(10,3) = -0.2314285703D2
 d_10_10(10,4) = 0.6300000018D2
 d_10_10(10,5) = -0.1133999993D3
 d_10_10(10,6) = 0.1417499991D3
 d_10_10(10,7) = -0.1259999992D3
 d_10_10(10,8) = 0.8099999995D2
 d_10_10(10,9) = -0.4049999997D2
 d_10_10(10,10) = 0.1273035716D2



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

