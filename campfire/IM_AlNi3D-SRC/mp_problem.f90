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
  use workspace
  use physicaldata
  implicit none


  integer :: i,j,i1,j1
  double precision :: FreeEnergy,FE,pi = 3.1415926535897932384626433
  LOGICAL :: file_exists
  double precision :: dh,h,a1,a2,dc0,dc1,c0,c1,S,x0,x1,a0,SS,xx0,xx1
  integer :: ii,n,k,np,l
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision ::  child1(10,10),child2(10,10),child3(10,10),child4(10,10),parent(10,10)
  double precision :: pcb_LagrangePoly,pcb_dLagrangePoly,pcb_d2LagrangePoly
  doubleprecision :: y10(10)=(/-1.,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1./)
  doubleprecision :: z10(10)=(/-1.080466092,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1.080466092/)
  doubleprecision :: y6(10)=(/-1.,-0.765055323929464692851,-0.2852315164806450963142,0.2852315164806450963142,0.765055323929464692851,1.,0.,0.,0.,0./)
  doubleprecision ::  y4(10)=(/-1.,-0.447213595499957939282,0.447213595499957939282,1.,0.,0.,0.,0.,0.,0./)
  double precision :: Trans(10,10),diff(10,10),psi=1.618033988749895
  double precision :: tmp_g_T0,RootFinderBiSection2,T_melt,cS,cL,RootFinderBiSection,cSmin,cLmin,df,radius
!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set parameters that affect global generic parameters



  tstep_increase_factor = 1.1
  tstep_decrease_factor = 0.8
  multigrid_on=1
  weight =1.0
  smooth_count_max = 1! if explicit mp_smooth set to 1
  solve_count_max = 1
  low_vcycle_count = 2 
  high_vcycle_count =3
  max_v_cycle_count = 3
  defect_tol_min = 1e-7 !if explicit set to ~1e-1
!  defect_tol_max = 1e-1
  defect_too_big = 1e1
  min_dt = 1.0e-16
!  output_rate = 100
  max_time_iter = 999000
  verbose=2           ! How much terminal output
  numglobalref=1       ! How many levels of global refinement are needed initially
  g_linear=.true.     !you also need to change inter_mask in paramesh/source/amr_initialize.F90 for .false.



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  dt = 1.0e-6
  g_D(1) = 1d0
  g_D(2) = 1d-4
  call input_parameter_file_reading() ! replace the 8 above if CASE.inp present


  interp_mask_unk(:) = 1
  interp_mask_work(:) = 1
  interp_mask_unk_res(:)   = 1
  interp_mask_work_res(:)  = 1



  g_gpure=g_pure
  if(AlNi.eq.7)then
     g_pure=.true.
     g_c0=0d0
  endif
  g_vel_max = g_vel
  INQUIRE(FILE="Antitrapping.txt", EXIST=file_exists) 
  if(file_exists)then
     open(unit=199,file="Antitrapping.txt")
     read(199,'(3F14.7)')g_alpha,g_c_min,g_c_max
     close(199)
  endif

  INQUIRE(FILE="Common_tan_params.txt", EXIST=file_exists) 

  if(file_exists.and.g_plapp)then

     open(unit=199,file="Common_tan_params.txt")
     read(199,'(8F14.5)')g_c_min,g_f_c_min,g_f1_c_min,g_f2_c_min,g_c_max,g_f_c_max,g_f1_c_max,g_f2_c_max
     close(199)
     write(*,*) "read Common_tan_params.txt"
     write(*,'(8F14.5)')g_c_min,g_f_c_min,g_f1_c_min,g_f2_c_min,g_c_max,g_f_c_max,g_f1_c_max,g_f2_c_max
  endif
  if(g_alpha.gt.0d0)g_alpha0 = g_alpha
  output_rate = g_output_rate
  g_delta(1)=g_lambda*g_d0
  !  if(g_s_skew)g_delta(1)=g_delta(1)*2d0 
  
  if(AlNi.eq.8)then !AlFe
     g_delta(1)=g_lambda*g_d0
     g_delta(2) = g_delta(1)
     g_delta(2)=g_delta(1)
     g_W(1)=6.*g_sigma(1)/g_delta(1)
     g_W(2)=6.*g_sigma(2)/g_delta(2)
     
!     g_Dch=g_D(1)*g_Dch
     g_tau = g_delta(1)**2/g_Dch
     g_lengthSQ=g_delta(1)**2
     g_D_tilde=g_Le*g_D(1)/g_Dch
     
     
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
     
     g_sigma(1)=1.0
     g_sigma(2)=1.0
     g_T(1) = g_T0
     g_T(2) = g_T0
     g_MolPerVol(1)=1.0
     g_MolPerVol(2)=1.0


     tmp_g_T0 = g_T0
     if(.not.g_toy_model)write(*,*)"melting point is",RootFinderBiSection2(.false.)
     T_melt = g_T0
     g_T0 = tmp_g_T0

     if(g_T0.ge.T_melt)then
        write(*,*)"Melting. T melt =",T_melt

     endif

!!!!!!!!!!!!!

! Even if I say so myself , this is a nifty routine
! to find the common tangent values c_Liq, and c_Sol.
! the method first finds the minimums and then quickly iterates (6 times) 
     if(.not.g_toy_model)then
        cS=RootFinderBiSection(.false.,0d0,0.)
        cSmin = cS
        cL=RootFinderBiSection(.false.,1d0,0.)
        !     cL = cS-0.1 !common tangent is to the left!
        cLmin = cL
        
        df = (FreeEnergy(cS,0.,0d0,0,.false.)-FreeEnergy(cL,0.,1d0,0,.false.))/(cS-cL)
        write(*,*)cL,cS,df
        
        do k = 1, 10
           !        df=(f3(cS)-f4(cL))/(cS-cL)
           df = (FreeEnergy(cS,0.,0d0,0,.false.)-FreeEnergy(cL,0.,1d0,0,.false.))/(cS-cL)
           cS=RootFinderBiSection(.false.,0d0,df)
           cL=RootFinderBiSection(.false.,1d0,df)
           write(*,'(6E17.9)')cL,cS,df
        enddo
        
        xx0 = cS
        xx1 = cL
        
        df = df*g_slope_factor
        g_df = df
        
        
        ! quadratic
        !y=1/2*g_AlFe_c(1,4)(x-g_AlFe_c(1,1))^2+g_AlFe_c(1,2)(x-g_AlFe_c(1,1))+g_AlFe_c(1,3)
        if(g_xx0.gt.0d0)then
           
           g_AlFe_c(1,1) = g_xx0!xx0*(1d0-g_xx1) !solid  
           g_AlFe_c(2,1) = g_xx1!+xx1*(1d0-g_xx1) !liquid
           
           
           g_AlFe_c(1,2) = FreeEnergy(xx0,0d0,0d0,2,.false.)
           g_AlFe_c(2,2) = FreeEnergy(xx1,0d0,1d0,2,.false.)
           g_AlFe_c(1,3) = FreeEnergy(xx0,0d0,0d0,0,.false.) + df*(g_xx0 - xx0)
           g_AlFe_c(2,3) = FreeEnergy(xx1,0d0,1d0,0,.false.) + df*(g_xx1 - xx1)
           g_AlFe_c(1,4) = g_accel_factor*(FreeEnergy(xx0+1d-6,0d0,0d0,2,.false.) - FreeEnergy(xx0,0d0,0d0,2,.false.))/1d-6
           g_AlFe_c(2,4) = g_accel_factor_liq*(FreeEnergy(xx1+1d-6,0d0,1d0,2,.false.) - FreeEnergy(xx1,0d0,1d0,2,.false.))/1d-6
!!$
           
           
        else
           g_AlFe_c(1,1) = xx0
           g_AlFe_c(2,1) = xx1
           g_AlFe_c(1,2) = FreeEnergy(xx0,0d0,0d0,2,.false.)
           g_AlFe_c(2,2) = FreeEnergy(xx1,0d0,1d0,2,.false.)
           g_AlFe_c(1,3) = FreeEnergy(xx0,0d0,0d0,0,.false.)
           g_AlFe_c(2,3) = FreeEnergy(xx1,0d0,1d0,0,.false.)
           g_AlFe_c(1,4) = g_accel_factor*(FreeEnergy(xx0+1d-6,0d0,0d0,2,.false.) - FreeEnergy(xx0,0d0,0d0,2,.false.))/1d-6
           g_AlFe_c(2,4) = (FreeEnergy(xx1+1d-6,0d0,1d0,2,.false.) - FreeEnergy(xx1,0d0,1d0,2,.false.))/1d-6
        endif
     endif
     if(g_toy_model)then
        xx0 = g_xx0
        xx1 = g_xx1

        g_AlFe_c(1,1) = g_xx0 !solid                                                                                                                                                
        g_AlFe_c(2,1) = g_xx1 !liquid                                                                                                                                               
        g_AlFe_c(1,2) = g_slope_factor
        g_AlFe_c(2,2) = g_slope_factor
        g_AlFe_c(1,3) = 0.0
        g_AlFe_c(2,3) = g_slope_factor*(g_xx1-g_xx0)
        g_AlFe_c(1,4) = g_accel_factor
        g_AlFe_c(2,4) = g_accel_factor_liq
     endif

!!$                  write(*,'(2F17.5)') g_AlFe_c
!!$                  stop
     
     g_approx = .true.
     
! original 3D code
!!$     g_L(1) = abs(FreeEnergy(xx0,0d0,0d0,0,.true.)-FreeEnergy(xx1,0d0,1d0,0,.true.))
!!$     g_L(2)=g_L(1)
!!$
!!$     g_L = g_R*g_T0
!!$     g_L = T_melt/ (T_melt-g_t0)*g_L
!!!!!!!!!!!!!!!

! from 2D code
!        g_L = abs((g_AlFe_c(1,1)-g_AlFe_c(2,1))*df*T_melt/(T_melt-g_T0))/g_drive 
! using the above at T=1192.8 which is near T_melt 1192.9 we obtain...
!     g_L = 2.19634466191911

!!$     g_L(1) =  abs((g_AlFe_c(1,3)- g_AlFe_c(2,3)))
!!$     g_L(2) =  abs((g_AlFe_c(1,3)- g_AlFe_c(2,3)))
!!%     g_L = 6d0*abs((g_AlFe_c(1,3)- g_AlFe_c(2,3)))/(g_lambda*g_velocity) 
!!$

!!$     g_L = 3d0*abs(g_AlFe_c(1,3)- g_AlFe_c(2,3))*g_Radius_min/g_lambda 
!!$     g_L = g_L *g_accel_factor/4d0*(g_xx1-g_xx0)/g_slope_factor

     g_L = 3d0*g_Radius_min/g_lambda
     g_L = g_L*(g_AlFe_c(2,1) - g_AlFe_c(1,1))**2 * g_accel_factor/4d0
!     g_L = g_L*(g_AlFe_c(2,1) - g_AlFe_c(1,1))**2 * g_accel_factor/2d0


!     g_L = g_L 
!     g_L = g_L*g_Dch

!!     g_Dch = g_Dch *g_velocity/(g_lambda*6d0)
!!!!!!!!!!
     write(*,*)"Latent heat",g_L(1)
     
     write(*,*)xx0,xx1

!     write(*,*)g_L(1),abs(FreeEnergy(xx0,0d0,0d0,0,.true.)-FreeEnergy(xx1,0d0,1d0,0,.true.))
 !    stop
!!$    open (unit=199,file="energy_out.txt")
!!$            dh = 0.27246-0.2350
!!$     do ii=1,999
!!$        h= ii*0.001
!!$        a1 = FreeEnergy(h,0d0,1d0,0,.false.)/g_L(1)
!!$        a2 = FreeEnergy(0.235+h*dh,0d0,0d0,0,.false.)/g_L(2)
!!$        
!!$        write(199,'(4F17.5)')h,a1,0.235+h*dh,a2
!!$     enddo
!!$     close(199)
!!$     
!!$     open (unit=199,file="energy_out_quad.txt")
!!$     !           dh = 0.27246-0.2350
!!$     dh=0.1
!!$     do ii=1,999
!!$        h= ii*0.001
!!$        a1 = FreeEnergy(h,0d0,1d0,0,.false.)/g_L(1)
!!$        !              a2 = FreeEnergy(0.235+h*dh,0d0,0d0,0,0,0,.true.)/g_L(2)
!!$!        a2 = FreeEnergy(0.2+h*dh,0d0,0d0,0,.true.)/g_L(2)
!!$        a2 = FreeEnergy(h,0d0,0d0,0,.false.)/g_L(2)
!!$        
!!$!        write(199,'(4F17.5)')h,a1,0.2+h*dh,a2
!!$        write(199,'(4F17.5)')h,a1,a2   
!!$     enddo
!!$     close(199)
!!$     open (unit=199,file="energy_params.txt")
!!$     write(199,'(5F17.5)')xx0,FreeEnergy(xx0,0d0,0d0,0,.false.)/g_L(1),xx1,FreeEnergy(xx1,0d0,1d0,0,.false.)/g_L(1),g_L(1)
!!$     close(199)




   l=0
   do i=-1,1,2
      do j=-1,1,2
         do k=-1,1,2
            l=l+1
            g_Dodec(1,l)=i*1d0
            g_Dodec(2,l)=j*1d0
            g_Dodec(3,l)=k*1d0
         enddo
      enddo
   enddo
   
   do i=-1,1,2
      do j=-1,1,2
         l=l+1
         g_Dodec(1,l)=i*psi
         g_Dodec(2,l)=j*(psi-1d0)
         g_Dodec(3,l)=0d0
      enddo
   enddo

   
   do i=-1,1,2
      do j=-1,1,2
         l=l+1
         g_Dodec(2,l)=i*psi
         g_Dodec(3,l)=j*(psi-1d0)
         g_Dodec(1,l)=0d0
      enddo
   enddo

   
   do i=-1,1,2
      do j=-1,1,2
         l=l+1
         g_Dodec(3,l)=i*psi
         g_Dodec(1,l)=j*(psi-1d0)
         g_Dodec(2,l)=0d0
      enddo
   enddo
   
   g_Dodec = g_Dodec/sqrt(3.)
   
!!!!!!!!!!!!! Truncated Octahedron !!!!!!!!!!
   l = 0
   do i=-2,2
      do j=-2,2
         do k=-2,2
            if(abs(i).ne.abs(j).and.abs(i).ne.abs(k).and.abs(k).ne.abs(j))then
               l = l +1
               g_TrunOct(1,l)=i
               g_TrunOct(2,l)=j
               g_TrunOct(3,l)=k
            endif
         enddo
      enddo
   enddo
   
   g_TrunOct = g_TrunOct/sqrt(5.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   do k=0,1
      do i=1,6
         g_hexcol(1,i+k*6)=sin(pi*(i-1d0)/3d0)
         g_hexcol(2,i+k*6)=cos(pi*(i-1d0)/3d0)
         g_hexcol(3,i+k*6)=aspect_ratio*(k-0.5)
      enddo
   enddo
   radius = sqrt(g_HexCol(1,1)**2+g_HexCol(2,1)**2+g_HexCol(3,1)**2)
   g_HexCol = G_HexCol/radius

   call set_mollifier

   g_cube(1,1)=1d0
   g_cube(2,1)=1d0
   g_cube(3,1)=1d0

   g_cube(1,2)=-1d0
   g_cube(2,2)=1d0
   g_cube(3,2)=1d0

   g_cube(1,3)=1d0
   g_cube(2,3)=-1d0
   g_cube(3,3)=1d0

   g_cube(1,4)=-1d0
   g_cube(2,4)=-1d0
   g_cube(3,4)=1d0

   g_cube(1,5)=1d0
   g_cube(2,5)=1d0
   g_cube(3,5)=-1d0

   g_cube(1,6)=-1d0
   g_cube(2,6)=1d0
   g_cube(3,6)=-1d0

   g_cube(1,7)=1d0
   g_cube(2,7)=-1d0
   g_cube(3,7)=-1d0

   g_cube(1,8)=-1d0
   g_cube(2,8)=-1d0
   g_cube(3,8)=-1d0

!g_cube = g_cube/sqrt(3d0) ! put vertices on unit sphere


else
   write(*,*)"only AlNi=8 now coded"
   write(*,*)" Alni=1 is Al3Ni2;AlNi=2 is AlNi, AlNi=3 is SiGe AlNi=4 is AlSi AlNi=5 is PbSn, AlNi=6 is NiCu, AlNi=7 is pure Al"
   stop
endif

  min_dx = g_grid/real(2**(g_max_lvl-1)*nxb)
!!$  do i=1,6
!!$     g_phex(1,i)=cos(pi*(i-1d0)/3d0+pi/30.)
!!$     g_phex(2,i)=sin(pi*(i-1d0)/3d0+pi/30.)
!!$  enddo
!!$  g_phex(1,3)=cos(0.8*3.141592654+pi/30.)
!!$  g_phex(2,3)=sin(0.8*3.141592654+pi/30.)










  ! Do we really need all these values?
  N_unk=nvar/nunkvbles  !some spare capacity for output

  N_phases = 1
  M_unk = N_unk*nunkvbles

  total_vars=N_unk



  NN_unk= N_unk*nxb*nxb*nxb  



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Derived assignments - do not change
  epsilon_tilde = 4.0*epsilon/(1.0-3.0*epsilon)  
  ! See page 26 of Jan's thesis
  A_0 = 1.0-3.0*epsilon                           ! See page 26 of Jan's thesis



end subroutine app_parameter_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_domain_setup14May2020(pid, noprocs)
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

  grid_max=g_grid
  eighthdmn = 0

  if (eighthdmn.eq.1) grid_min = 0.0
  mg_min_lvl = 1

  nBx = g_nbx
  nBy = g_nby
  nBz = g_nbz

  if (ndim.eq.2) nBz = 1

  dBx = grid_max-grid_min
  dBy = dBx
  dBz = dBx

  max_x = dBx*nBx
  ! set a limit on the refinement level
  lrefine_min = mg_min_lvl


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
               
               if (ndim.eq.3) then
                  coord(3,lnblocks) = nBz*grid_min+(k-0.5)*dBz! Fix z coord at middle of layer
                  bnd_box(1,3,lnblocks) = nBz*grid_min+(k-1.0)*dBy
                  bnd_box(2,3,lnblocks) = nBz*grid_min+k*dBy
               endif

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

!!!! ------- MCH initialize block_list ---------
  if(morton_flag==.true.) then
    allocate(block_list(lnblocks))
    write(*,*) 'add block_list'
  endif
!!!! -------------------------------------------

end subroutine app_domain_setup14may2020

!!!!!!!!
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

  grid_max=g_grid
  eighthdmn = 0

  if (eighthdmn.eq.1) grid_min = 0.0
  mg_min_lvl = 1

  nBx = g_nbx
  nBy = g_nby
  nBz = g_nbz

  if (ndim.eq.2) nBz = 1

  dBx = grid_max-grid_min
  dBy = dBx
  dBz = dBx

  max_x = dBx*nBx
  ! set a limit on the refinement level
  lrefine_min = mg_min_lvl


  lrefine_max = g_max_lvl


   numglobalref=1






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!          Leave untouched from here down for phase field code                                       !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


! CEG Initialise my startlist for the blocklist for that grid level
  allocate(block_starts(lrefine_max+1))


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

              if (ndim.eq.3) then
                 coord(3,lnblocks) = nBz*grid_min+(k-0.5)*dBz        ! Fix z coord at middle of layer
                 bnd_box(1,3,lnblocks) = nBz*grid_min+(k-1.0)*dBy
                 bnd_box(2,3,lnblocks) = nBz*grid_min+k*dBy
              endif

              bsize(:,lnblocks) = bnd_box(2,:,lnblocks) - bnd_box(1,:,lnblocks)
              nodetype(lnblocks) = 1
              lrefine(lnblocks) = 1
     
              ! boundary conditions need to represent whether far field or symmetry
!              neigh(2,:,lnblocks) = 0
              if (i.eq.1.and.eighthdmn.eq.1) then
                 neigh(:,1,lnblocks) = bc_cond_sym
              else if (i.eq.1) then
                 neigh(:,1,lnblocks) = bc_cond_far
              else
                 neigh(1,1,lnblocks) = lnblocks-1
                 neigh(2,2,lnblocks) = 0
              endif
              if (i.eq.nBx) then
                 neigh(:,2,lnblocks) = bc_cond_far
              else
                 neigh(1,2,lnblocks) = lnblocks+1
                 neigh(2,2,lnblocks) = 0
              endif

              if (j.eq.1.and.eighthdmn.eq.1) then
                 neigh(:,3,lnblocks) = bc_cond_sym
              else if (j.eq.1) then
                 neigh(:,3,lnblocks) = bc_cond_far
              else
                 neigh(1,3,lnblocks) = lnblocks-nBx
                 neigh(2,3,lnblocks) = 0
              endif
              if (j.eq.nBy) then
                 neigh(:,4,lnblocks) = bc_cond_far
              else
                 neigh(1,4,lnblocks) = lnblocks+nBx
                 neigh(2,4,lnblocks) = 0
              endif

              if (ndim.eq.3) then
                 if (k.eq.1.and.eighthdmn.eq.1) then
                    neigh(:,5,lnblocks) = bc_cond_sym
                 else if (k.eq.1) then
                    neigh(:,5,lnblocks) = bc_cond_far
                 else
                    neigh(1,5,lnblocks) = lnblocks-nBx*nBy
                    neigh(2,5,lnblocks) = 0
                 endif
                 if (k.eq.nBz) then
                    neigh(:,6,lnblocks) = bc_cond_far
                 else
                    neigh(1,6,lnblocks) = lnblocks+nBx*nBy
                    neigh(2,6,lnblocks) = 0
                 endif
              endif
              neigh(2,:,lnblocks) = 0
              refine(lnblocks)=.true.
           end do
        end do
     end do
     
     min_dx = (bnd_box(2,1,1)-bnd_box(1,1,1))/(2**(lrefine_max-1)*nxb)
     if (pid.eq.0.and.verbose.ge.2) print *,"Min available dx=", min_dx
     max_dt = min_dx
     max_dt = g_max_dt
  endif

  ! Set boundary conditions
  if (eighthdmn.eq.1) then
     boundary_index(1)       = bc_cond_sym
  else  
     boundary_index(1)       = bc_cond_far
  endif
  boundary_index(2)       = bc_cond_far
  
  ! y boundaries (boundary boxes 3 and 4)
  if (ndim .ge. 2) then
     if (eighthdmn.eq.1) then
        boundary_index(3)       = bc_cond_sym
     else  
        boundary_index(3)       = bc_cond_far
     endif
     boundary_index(4)     = bc_cond_far
  end if
  
  ! z boundaries
  if (ndim .eq. 3) then
     if (eighthdmn.eq.1) then
        boundary_index(5)       = bc_cond_sym
     else  
        boundary_index(5)       = bc_cond_far
     endif
     boundary_index(6)          = bc_cond_far
  end if

  ! ensure all processors have min_dx and max_dt set
  call MPI_Bcast(min_dx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(max_dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

!!!! ------- MCH initialize block_list ---------                                                                                                                                                                                             
  if(morton_flag==.true.) then
    allocate(block_list(lnblocks))
    write(*,*) 'add block_list'
  endif
!!!! -------------------------------------------   



end subroutine app_domain_setup



!!!!!!

