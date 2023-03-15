!mp_problem.f90 sets up the problem parameters

!!$        g_AlFe_c(1,1) = xx0 !solid                                                                                                                                               !!$        g_AlFe_c(2,1) = xx1 !liquid                                                                                                                                              !!$        g_AlFe_c(1,2) = cSmin !solid                                                                                                                                             !!$        g_AlFe_c(2,2) = cLmin !liquid                                                                                                                                            !!$        g_AlFe_c(1,3) = FreeEnergy(cSmin,0d0,0d0,0,.false.)
!!$        g_AlFe_c(2,3) = FreeEnergy(cLmin,0d0,1d0,0,.false.)
!!$        g_AlFe_c(1,4) = FreeEnergy(xx0,0d0,0d0,0,.false.)
!!$        g_AlFe_c(2,4) = FreeEnergy(xx1,0d0,1d0,0,.false.)
!!$        g_AlFe_c(1,5) = (FreeEnergy(cSmin+1d-7,0d0,0d0,2,.false.)-FreeEnergy(cSmin-1d-7,0d0,0d0,2,.false.))/2d-7
!!$        g_AlFe_c(2,5) = (FreeEnergy(cLmin+1d-7,0d0,1d0,2,.false.)-FreeEnergy(cLmin-1d-7,0d0,1d0,2,.false.))/2d-7
!!$        g_AlFe_c(1,6) = df !solid common tangents                                                                                                                                !!$        g_AlFe_c(2,6) = df !liquid        


double precision function U_CurveL(c,d)
 use solution_parameters
  implicit none
  integer, intent(in) :: d
  double precision, intent(in) :: c
  double precision :: t,dc,f0,ddf0,f1,df1


  t = (c-g_AlFe_c(2,2))/(g_AlFe_c(2,1)-g_AlFe_c(2,2))
  dc = g_AlFe_c(2,1)-g_AlFe_c(2,2)

  f0 = g_AlFe_c(2,3) !energy at min
  f1 = g_AlFe_c(2,4) !energy at common tangent
  ddf0 = g_AlFe_c(2,5)*dc**2 !accel at min
  df1 = g_AlFe_c(2,6)*dc !tan at common tangent

  if(d.eq.0)then
     U_CurveL = f0 + ddf0*t**2/2 + (-4*f0 - ddf0 + 4*f1 - df1)*t**3 + (3*f0 + ddf0/2 - 3*f1 + df1)*t**4
  elseif(d.eq.1)then
     U_CurveL = (    ddf0*t    +3* (-4*f0 - ddf0 + 4*f1 - df1)*t**2 + 4*(3*f0 + ddf0/2 - 3*f1 + df1)*t**3)/dc
  elseif(d.eq.2)then
     U_CurveL = (    ddf0      +6* (-4*f0 - ddf0 + 4*f1 - df1)*t + 12*(3*f0 + ddf0/2 - 3*f1 + df1)*t**2)/dc**2
  else
     write(*,*)"d> 2 not coded, mp_problem.f90"
     stop
  endif
end function U_CurveL

double precision function U_CurveS(c,d)
 use solution_parameters
  implicit none
  integer, intent(in) :: d
  double precision, intent(in) :: c
  double precision :: t,dc,f0,ddf0,f1,df1


  t = (c-g_AlFe_c(1,2))/(g_AlFe_c(1,1)-g_AlFe_c(1,2))
  dc = g_AlFe_c(1,1)-g_AlFe_c(1,2)

  f0 = g_AlFe_c(1,3)
  f1 = g_AlFe_c(1,4)
  ddf0 = g_AlFe_c(1,5)*dc**2
  df1 = g_AlFe_c(1,6)*dc

  if(d.eq.0)then
     U_CurveS = f0 + 0.5*ddf0*t**2 + (-4*f0 - ddf0 + 4*f1 - df1)*t**3 + (3*f0 + 0.5*ddf0 - 3*f1 + df1)*t**4
  elseif(d.eq.1)then
     U_CurveS = (    ddf0*t    +3* (-4*f0 - ddf0 + 4*f1 - df1)*t**2 + 4*(3*f0 + 0.5*ddf0 - 3*f1 + df1)*t**3)/dc
  elseif(d.eq.2)then
     U_CurveS = (    ddf0      +6* (-4*f0 - ddf0 + 4*f1 - df1)*t + 12*(3*f0 + 0.5*ddf0 - 3*f1 + df1)*t**2)/dc**2
  else
     write(*,*)"d> 2 not coded, mp_problem.f90"
     stop
  endif
end function U_CurveS


double precision function Hermite(c,i,d)
  use solution_parameters
  implicit none
  integer, intent(in) :: i,d
  double precision, intent(in) :: c
  double precision :: t,dt
  
  t = (c-g_AlFe_c(1,2))/(g_AlFe_c(1,1)-g_AlFe_c(1,2))
  dt = 1./(g_AlFe_c(1,1)-g_AlFe_c(1,2))  

  if(d.eq.0)then
     if(i.eq.1)Hermite = -1.25*t**4 + 9./2.*t**3 - 17./4.*t**2 + 1d0
     if(i.eq.2)Hermite =t**4 - 4*t**3 + 4*t**2
     if(i.eq.3)Hermite =0.25*t**4 - 0.5*t**3 + 0.25*t**2
     if(i.eq.4)Hermite =(-t**4 + 3*t**3 - 2*t**2)/dt
  elseif(d.eq.1)then

     if(i.eq.1)Hermite = -5*t**3 + 27./2.*t**2 - 17./2.*t
     if(i.eq.2)Hermite =4*t**3 - 12*t**2 + 8*t
     if(i.eq.3)Hermite =t**3 - 1.5*t**2 + 0.5*t
     if(i.eq.4)Hermite =(-4*t**3 + 9*t**2 - 4*t)/dt

     Hermite = Hermite*dt
  else
     write(*,*)"Hermite"
     stop
  end if
end function Hermite


double precision function HermiteL(c,i,d)
  use solution_parameters
  implicit none
  integer, intent(in) :: i,d
  double precision, intent(in) :: c
  double precision :: t,dt

  t = (c-g_AlFe_c(2,2))/(g_AlFe_c(2,1)-g_AlFe_c(2,2))
  dt = 1./(g_AlFe_c(2,1)-g_AlFe_c(2,2))  


  if(d.eq.0)then
     if(i.eq.1)HermiteL =t**4 - 2.*t**2 + 1.
     if(i.eq.2)HermiteL =-1.25*t**4 + 0.5*t**3 + 1.75*t**2
     if(i.eq.3)HermiteL =0.25*t**4 - 0.5*t**3 + 0.25*t**2
     if(i.eq.4)HermiteL =(0.5*t**4 - 0.5*t**2)/dt
  elseif(d.eq.1)then
     if(i.eq.1)HermiteL =4*t**3 - 4*t
     if(i.eq.2)HermiteL =-5*t**3 + 1.5*t**2 + 3.5*t
     if(i.eq.3)HermiteL =t**3 - 1.5*t**2 + 0.5*t
     if(i.eq.4)HermiteL =(2*t**3 - t)/dt

     HermiteL = HermiteL*dt

  else
     write(*,*)"HermiteL"
     stop
  end if
end function HermiteL

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
!  use ExampleFuncs
  implicit none


  integer :: i,j,i1,j1
  double precision :: FreeEnergy,FE,pi = 3.1415926535897932384626433
  LOGICAL :: file_exists
  double precision :: dh,h,a1,a2,a3,a4,a5,dc0,dc1,c0,c1,S,x0,x1,a0,SS,xx0,xx1
  integer :: ii,n,k
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision ::  child1(10,10),child2(10,10),child3(10,10),child4(10,10),parent(10,10)
!  double precision :: pcb_LagrangePoly,pcb_dLagrangePoly,pcb_d2LagrangePoly
  doubleprecision :: y10(10)=(/-1.,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1./)
  doubleprecision :: z10(10)=(/-1.080466092,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1.080466092/)
  doubleprecision :: y6(10)=(/-1.,-0.765055323929464692851,-0.2852315164806450963142,0.2852315164806450963142,0.765055323929464692851,1.,0.,0.,0.,0./)
  doubleprecision ::  y4(10)=(/-1.,-0.447213595499957939282,0.447213595499957939282,1.,0.,0.,0.,0.,0.,0./)
  double precision :: Trans(10,10),diff(10,10)
  double precision :: RootFinderBiSection,RootFinderBiSection2,Matrix33(3,3),Vector3(3), Ans3(3),g_phi,Y_func,Cphi_f
  double precision :: df,cL,cS,cLmin,cSmin,potential,tmp_g_T0,T_melt,AlFe_c(2,4),c_i, Y_phi
  double precision :: Gibbs_FE_e2_AlFe!(c,T,d_c,d_T)

!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set parameters that affect global generic parameters
  tstep_increase_factor = 1.1
  tstep_decrease_factor = 0.5
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
  max_time_iter = 9999000
  verbose=2           ! How much terminal output
  numglobalref=1       ! How many levels of global refinement are needed initially
  g_linear=.false.     !you also need to change inter_mask in paramesh/source/amr_initialize.F90 for .false.
  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








  dt = 1.0e-6

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Physical phase field case - variables declared in pf_modules.f90

!  namelist/input_variables/g_lambda,g_beta,g_c0,g_t0,heat_c_term,Grad_DW_term,heat_phi_term,cross_term,g_Le
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! from mp_modules

!  g_lambda=80.0!4.0 giving g_delta(1)=4e-8
!  g_beta = 0.0
!  g_c0 = 0.5
!  g_T0 = (g_L(1)*(1-g_c0)+g_L(2)*g_c0)/(g_L(1)*(1-g_c0)/g_T(1)+g_L(2)*g_c0/g_T(2))
!  heat_c_term = .false.
!  Grad_DW_term=.false. 
!  heat_phi_term = .false.
!  cross_term = .false.
!  g_Le=100.0
!  nuc_radius = 20.0!2.2d-8/g_delta(1)
!  delta=(175.+273.-g_T0)/g_T0 !over written in CASE.inp
!  epsilon = 0.02!(1d-6)/4.
!  g_mu(1)=0.002640!0.04
!  g_mu(2)=0.003120!0.04
  g_D(1) = 1d-9
  g_D(2) = 1d-12
  call input_parameter_file_reading() ! replace the 8 above if CASE.inp present

  if(g_linear)then
     interp_mask_unk(:) = 1
     interp_mask_work(:) = 1
     interp_mask_unk_res(:)   = 1
     interp_mask_work_res(:)  = 1
  else
!!$     interp_mask_unk(:)   = nxb+22
!!$     interp_mask_work(:)  = nxb+22
!!$     interp_mask_unk_res(:)   = nxb+22
!!$     interp_mask_work_res(:)  = nxb+22
!!$     g_order=nxb+2
!!$
!!$     Trans=0d0
!!$     do i=1,g_order
!!$        Trans(i,i)=1d0
!!$     enddo
!!$     Trans(1,1)=0.5
!!$     Trans(2,1)=0.5
!!$     Trans(g_order-1,g_order)=0.5
!!$     Trans(g_order,g_order)=0.5
!!$     if(g_linear)then
!!$     elseif(g_order.eq.10)then
!!$        z10=y10
!!$        g_y10=y10
!!$     elseif(g_order.eq.6)then
!!$        z10=y6
!!$        g_y10=y6
!!$     elseif(g_order.eq.4)then
!!$        z10=y4
!!$        g_y10=y4
!!$     else
!!$        write(*,*)"mp_problem.f90. order=4,6,10 only"
!!$        stop
!!$     endif
!!$     do i=1,g_order
!!$        do j=1,g_order
!!$           diff(j,i)=pcb_LagrangePoly(g_order,z10(i),j)
!!$        enddo
!!$     enddo
!!$     g_diff0=0d0
!!$     do i=1,g_order
!!$        do j=1,g_order
!!$           do k=1,g_order
!!$              g_diff0(i,j)  = g_diff0(i,j)+Trans(i,k)*diff(k,j)
!!$           enddo
!!$        enddo
!!$     enddo
!!$     
!!$     
!!$     do i=1,g_order
!!$        do j=1,g_order
!!$           diff(j,i)=pcb_dLagrangePoly(g_order,z10(i),j)
!!$        enddo
!!$     enddo
!!$     g_diff1=0d0
!!$     do i=1,g_order
!!$        do j=1,g_order
!!$           do k=1,g_order
!!$              g_diff1(i,j)  = g_diff1(i,j)+Trans(i,k)*diff(k,j)
!!$           enddo
!!$        enddo
!!$     enddo
!!$     
!!$     do i=1,g_order
!!$        do j=1,g_order
!!$           diff(j,i)=pcb_d2LagrangePoly(g_order,z10(i),j)
!!$        enddo
!!$     enddo
!!$     g_diff2=0d0
!!$     do i=1,g_order
!!$        do j=1,g_order
!!$           do k=1,g_order
!!$              g_diff2(i,j)  = g_diff2(i,j)+Trans(i,k)*diff(k,j)
!!$           enddo
!!$        enddo
!!$     enddo
  endif



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
     
     

     g_Dch=g_D(1)*g_Dch !g_Dch from input file
     g_tau = g_delta(1)**2/g_Dch
     g_lengthSQ=g_delta(1)**2
     g_D_tilde=g_Le*g_D(1)/g_Dch
     
     
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
     
     g_sigma(1)=1.0!0.16895d0!0.033
     g_sigma(2)=1.0!0.35241d0!
     g_T(1) = g_T0!1.0!933.14d0
     g_T(2) = g_T0!1.0!1687d0
     g_MolPerVol(1)=1.0!2.7d3/26.98d-3
     g_MolPerVol(2)=1.0!2.57d3/28.085d-3
     g_left_tan = 0d0
     g_right_tan = 0d0
     g_left_acc =0d0
     g_right_acc = 0d0     
     g_right_3 = 0d0
     tmp_g_T0 = g_T0
     write(*,*)"melting point is",RootFinderBiSection2(.false.)
     T_melt = g_T0
     g_T0 = tmp_g_T0

!!$     if(g_T0.ge.T_melt)then
!!$        write(*,*)"Melting. T melt =",T_melt
!!$        stop
!!$     endif

!!!!!!!!!!!!!

! Even if I say so myself , this is a nifty routine
! to find the common tangent values c_Liq, and c_Sol.
! the method first finds the minimums and then quickly iterates (6 times) 

     cS=RootFinderBiSection(.false.,0d0,0.)
     cSmin = cS
     cL=RootFinderBiSection(.false.,1d0,0.)
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

!!!!!!!!!!!!!!!!!!
! Having found cL and cS now construct quadratic approximation to real free energy functions     

     g_L = g_R*g_T0

     xx0 = cS
     xx1 = cL
     if(g_Quad)then
        g_AlFe_c(1,1) = xx0 !solid  
        g_AlFe_c(2,1) = xx1 !liquid
        g_AlFe_c(1,1) = (1d0-g_slider)*g_AlFe_c(1,1)+g_slider*g_AlFe_c(2,1)
        g_AlFe_c(1,2) = FreeEnergy(xx0,0d0,0d0,2,.false.)!The tangents
        g_AlFe_c(2,2) = FreeEnergy(xx1,0d0,1d0,2,.false.)
        g_AlFe_c(1,3) = FreeEnergy(xx0,0d0,0d0,0,.false.)
        g_AlFe_c(2,3) = FreeEnergy(xx1,0d0,1d0,0,.false.)
        g_AlFe_c(1,3) = (1d0-g_slider)*g_AlFe_c(1,3)+g_slider*g_AlFe_c(2,3)
        g_AlFe_c(1,4) = (FreeEnergy(xx0+1d-7,0d0,0d0,2,.false.) - FreeEnergy(xx0,0d0,0d0,2,.false.))/1d-7!the 2nd derivatives
        g_AlFe_c(2,4) = (FreeEnergy(xx1+1d-7,0d0,1d0,2,.false.) - FreeEnergy(xx1,0d0,1d0,2,.false.))/1d-7

!        g_L(1) = g_AlFe_c(2,3) -g_AlFe_c(1,3)! abs(FreeEnergy(xx1,0d0,1d0,0,.false.)-FreeEnergy(xx0,0d0,0d0,0,.false.))
        
     else !quatic approximation using min values, tangent values and one other
        g_AlFe_c(1,1) = xx0 !solid 
        g_AlFe_c(2,1) = xx1 !liquid    
        g_AlFe_c(1,2) = cSmin !solid 
        g_AlFe_c(2,2) = cLmin !liquid
        g_AlFe_c(1,3) = FreeEnergy(cSmin,0d0,0d0,0,.false.)
        g_AlFe_c(2,3) = FreeEnergy(cLmin,0d0,1d0,0,.false.)
        g_AlFe_c(1,4) = FreeEnergy(xx0,0d0,0d0,0,.false.)
        g_AlFe_c(2,4) = FreeEnergy(xx1,0d0,1d0,0,.false.)
        g_AlFe_c(1,5) = (FreeEnergy(cSmin+1d-7,0d0,0d0,2,.false.)-FreeEnergy(cSmin-1d-7,0d0,0d0,2,.false.))/2d-7
        g_AlFe_c(2,5) = (FreeEnergy(cLmin+1d-7,0d0,1d0,2,.false.)-FreeEnergy(cLmin-1d-7,0d0,1d0,2,.false.))/2d-7
        g_AlFe_c(1,6) = df !solid common tangents
        g_AlFe_c(2,6) = df !liquid



!        g_L(1) = (g_AlFe_c(2,4) -g_AlFe_c(1,4))! abs(FreeEnergy(xx1,0d0,1d0,0,.false.)-FreeEnergy(xx0,0d0,0d0,0,.false.))
     endif



!!$                  write(*,'(2F17.5)') g_AlFe_c
!!$                  stop
!!$     
! compute c(phi) and c'(phi)
!!$     do i=-2,g_Ncphi+2
!!$        g_Cphi(i,1) = real(i*1.0/g_Ncphi)
!!$        g_Cphi(i,2) = RootFinderBiSection(.true.,g_Cphi(i,1),df)
!!$     enddo
!!$     do i=-1,g_Ncphi+1
!!$        g_Cphi(i,3) = g_Ncphi*0.5*(g_Cphi(i+1,2)-g_CPhi(i-1,2))
!!$     enddo



     g_df = df !make df global

     g_approx = .true. ! possibly redundant

     open (unit=199,file="potential.txt")
     do i=0,101
        write(199,'(1E12.5)')potential(i*0.01,0)
     enddo
     close(199)




     open (unit=199,file="energy_out.txt")
     do ii=1,999
        h= ii*0.001
        a1 = potential(0d0,0) +FreeEnergy(h,0d0,0d0,0,.true.)/g_L(1)
        a2 = potential(0.25,0)+FreeEnergy(h,0d0,Y_phi(0.25d0,0),0,.true.)/g_L(1)
        a3 = potential(0.5,0) +FreeEnergy(h,0d0,Y_phi(0.5,0),0,.true.)/g_L(1)
        a4 = potential(0.75,0)+FreeEnergy(h,0d0,Y_phi(0.75,0),0,.true.)/g_L(1)
        a5 = potential(1d0,0) +FreeEnergy(h,0d0,1d0,0,.true.)/g_L(1)
        
        write(199,'(6E17.5)')h,a1,a2,a3,a4,a5

!        write(199,'(3E17.5)')h,Y_phi(h,0),Y_phi(h,1)
     enddo
     close(199)


     open (unit=199,file="energy_out_quad.txt")
        write(199,'(2F17.7)')g_AlFe_c(:,1:2)
        write(199,'(2F17.7)')g_AlFe_c(:,3:6)/g_L(1)
     close(199)

     open (unit=199,file="energy_params.txt")
     write(199,'(5F17.5)')xx0,FreeEnergy(xx0,0d0,0d0,0,.false.)/g_L(1),xx1,FreeEnergy(xx1,0d0,1d0,0,.false.)/g_L(1),g_L(1)
     close(199)



     write(*,*)FreeEnergy(g_left_limit,0d0,0d0,2,.false.),FreeEnergy(g_right_limit,0d0,0d0,2,.false.)
     write(*,*)1d3*(FreeEnergy(g_left_limit+1d-3,0d0,0d0,2,.false.)-FreeEnergy(g_left_limit,0d0,0d0,2,.false.)),1d3*FreeEnergy(g_right_limit,0d0,0d0,2,.false.)-FreeEnergy(g_right_limit-1d-3,0d0,0d0,2,.false.)







     if(g_LUT)then !use Look Up Table LUT. Option here to change .true. to .false.. The latter uses actual FE.
        g_left_tan = FreeEnergy(g_left_limit,0d0,0d0,2,.false.)
        g_right_tan = FreeEnergy(g_right_limit,0d0,0d0,2,.false.)
        g_left_acc = 1d6*(FreeEnergy(g_left_limit+1d-6,0d0,0d0,2,.false.)-FreeEnergy(g_left_limit,0d0,0d0,2,.false.))
        g_right_acc = 1d6*(FreeEnergy(g_right_limit,0d0,0d0,2,.false.)-FreeEnergy(g_right_limit-1d-6,0d0,0d0,2,.false.))
        g_right_3 =1d12*(FreeEnergy(g_right_limit,0d0,0d0,2,.false.)-2d0*FreeEnergy(g_right_limit-1d-6,0d0,0d0,2,.false.)+FreeEnergy(g_right_limit-2d-6,0d0,0d0,2,.false.))
        write(*,*) g_right_tan, g_right_acc, g_right_3

        g_LUT = .false. ! suspend LUT to establish the table
        g_AlFe_LUT = 0d0
        do i = 1,g_Ncphi-1
           c_i = real(i*1.0/g_Ncphi)
              g_AlFe_LUT(i,1) = FreeEnergy(c_i,0d0,0d0,0,.false.)!Solid
              g_AlFe_LUT(i,2) = FreeEnergy(c_i,0d0,0d0,2,.false.)!Solid 1st derivative
              g_AlFe_LUT(i,3) = 0.5*(FreeEnergy(c_i+1d-7,0d0,0d0,2,.false.) - FreeEnergy(c_i-1d-7,0d0,0d0,2,.false.))/1d-7
              g_AlFe_LUT(i,4) = FreeEnergy(c_i,0d0,1d0,0,.false.)!Liquid
              g_AlFe_LUT(i,5) = FreeEnergy(c_i,0d0,1d0,2,.false.)!Liquid 1st derivatve
              g_AlFe_LUT(i,6) = 0.5*(FreeEnergy(c_i+1d-7,0d0,1d0,2,.false.) - FreeEnergy(c_i-1d-7,0d0,1d0,2,.false.))/1d-7
        enddo

        g_LUT=.true. !now using LU table.
     endif
     g_L(1) = abs(FreeEnergy(xx1,0d0,1d0,0,.false.) - FreeEnergy(xx0,0d0,0d0,0,.false.) )
     g_L(2)=g_L(1)
     g_L = T_melt/ (T_melt-g_t0)*g_L/g_drive
     write(*,*)"Latent heat",g_L(1)
     write(*,*)xx0,xx1

!!$


     open (unit=199,file="AlFeIM.txt")
     do ii=1,9999
        write(199,'(3F17.8)')0.0001*ii,FreeEnergy(0.0001*ii,0d0,0d0,0,.false.)/g_L(1),FreeEnergy(0.0001*ii,0d0,0d0,2,.false.)/g_L(1)
     enddo

     open (unit=199,file="AlFeIM_approx.txt")
     do ii=1,9999
        write(199,'(3F17.8)')0.0001*ii,FreeEnergy(0.0001*ii,0d0,0d0,0,.true.)/g_L(1),FreeEnergy(0.0001*ii,0d0,0d0,2,.true.)/g_L(1)
     enddo
     close(199)
!     stop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
  else
     write(*,*)"only AlNi=8 now coded"
     write(*,*)" Alni=1 is Al3Ni2;AlNi=2 is AlNi, AlNi=3 is SiGe AlNi=4 is AlSi AlNi=5 is PbSn, AlNi=6 is NiCu, AlNi=7 is pure Al"
     stop
  endif
!  g_Dch=1d-9
  min_dx = g_grid/real(2**(g_max_lvl-1)*nxb)
  do i=1,6
     g_phex(1,i)=cos(pi*(i-1d0)/3d0+pi/30.)
     g_phex(2,i)=sin(pi*(i-1d0)/3d0+pi/30.)
  enddo
  g_phex(1,3)=cos(0.8*3.141592654+pi/30.)
  g_phex(2,3)=sin(0.8*3.141592654+pi/30.)










  ! Do we really need all these values?
  N_unk=nvar/nunkvbles  !some spare capacity for output

  N_phases = 1
  M_unk = N_unk*nunkvbles

  total_vars=N_unk



  NN_unk= N_unk*nxb*nxb  



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Derived assignments - do not change
  epsilon_tilde = 4.0*epsilon/(1.0-3.0*epsilon)  
  ! See page 26 of Jan's thesis
  A_0 = 1.0-3.0*epsilon                           ! See page 26 of Jan's thesis



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

!!!! ------- MCH initialize block_list ---------
  if(morton_flag==.true.) then
    allocate(block_list(lnblocks))
    write(*,*) 'add block_list'
  endif
!!!! -------------------------------------------

end subroutine app_domain_setup







double precision function pcb_dLagrangePoly(order,x,j)
  use solution_parameters
implicit none
integer, intent (in) :: order,j
double precision, intent (in) :: x
integer :: i
!double precision :: y(10),pcb_LagrangePoly0
doubleprecision :: y4(10)=(/-1.,-0.447213595499957939282,0.447213595499957939282,1.,0.,0.,0.,0.,0.,0./)
doubleprecision :: y6(10)=(/-1.,-0.765055323929464692851,-0.2852315164806450963142,0.2852315164806450963142,0.765055323929464692851,1.,0.,0.,0.,0./)
doubleprecision :: y8(10)=(/-1.,-0.8717401485096066153375,-0.5917001814331423021445,-0.2092992179024788687687,0.2092992179024788687687,0.5917001814331423021445,0.8717401485096066153375  ,1.,0.,0./)
doubleprecision :: y10(10)=(/-1.,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1./)
!!$
!!$if(order.eq.4)then
!!$   y=y4
!!$elseif(order.eq.6)then
!!$   y=y6
!!$elseif(order.eq.8)then
!!$   y=y8
!!$elseif(order.eq.10)then
!!$   y=g_y10
!!$else
!!$   write(*,*)"order out of range in pcb_LagrangePoly: 4,6,8,10 only"
!!$   stop
!!$endif
!!$pcb_dLagrangePoly=0d0
!!$
!!$do i = 1,order
!!$   if(i.ne.j)pcb_dLagrangePoly = pcb_dLagrangePoly + pcb_LagrangePoly0(order,x,j,i,0)
!!$!   if(i.ne.j.and.i.ne.k)pcb_dLagrangePoly = pcb_dLagrangePoly + pcb_LagrangePoly(order,x,j,i)/(y(j)-y(i))
!!$enddo
!!$pcb_dLagrangePoly = pcb_dLagrangePoly/pcb_LagrangePoly0(order,y(j),j,0,0)
pcb_dLagrangePoly = 0d0
end function pcb_dLagrangePoly

double precision function pcb_d2LagrangePoly(order,x,j)
  use solution_parameters
implicit none
integer, intent (in) :: order,j
double precision, intent (in) :: x
integer :: i,k
!double precision :: y(10),pcb_LagrangePoly0
doubleprecision :: y4(10)=(/-1.,-0.447213595499957939282,0.447213595499957939282,1.,0.,0.,0.,0.,0.,0./)
doubleprecision :: y6(10)=(/-1.,-0.765055323929464692851,-0.2852315164806450963142,0.2852315164806450963142,0.765055323929464692851,1.,0.,0.,0.,0./)
doubleprecision :: y8(10)=(/-1.,-0.8717401485096066153375,-0.5917001814331423021445,-0.2092992179024788687687,0.2092992179024788687687,0.5917001814331423021445,0.8717401485096066153375  ,1.,0.,0./)
doubleprecision :: y10(10)=(/-1.,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1./)
!!$
!!$if(order.eq.4)then
!!$   y=y4
!!$elseif(order.eq.6)then
!!$   y=y6
!!$elseif(order.eq.8)then
!!$   y=y8
!!$elseif(order.eq.10)then
!!$   y=g_y10
!!$else
!!$   write(*,*)"order out of range in pcb_LagrangePoly: 4,6,8,10 only"
!!$   stop
!!$endif
!!$pcb_d2LagrangePoly=0d0
!!$
!!$do i = 1,order
!!$   if(i.ne.j)then
!!$   do k = 1,order
!!$      if(k.ne.j .and. k.ne.i)pcb_d2LagrangePoly = pcb_d2LagrangePoly + pcb_LagrangePoly0(order,x,j,i,k)
!!$   enddo
!!$   endif
!!$enddo
pcb_d2LagrangePoly = 0d0!pcb_d2LagrangePoly/pcb_LagrangePoly0(order,y(j),j,0,0)
end function pcb_d2LagrangePoly
