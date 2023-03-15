!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  Utilities for phase field stencils in 2-d and 3-d
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  History: James Green 2008-9, Chris Goodyer, 2009-2012
!!!!  Latest version: Peter Bollada, April 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function A_grad(phi,grad_phi)
use solution_parameters
implicit none
double precision, intent (in)::phi(3),grad_phi(3,2)
double precision :: r3(3,3,2),A3(3,3,2),r2(3,3)
integer :: i,j,a




do i=1,3
do j=1,3
do a=1,2
   r3(i,j,a) = phi(i)*grad_phi(j,a)-phi(j)*grad_phi(i,a)
enddo
enddo
enddo

do i=1,3
do j=1,3
   r2(i,j)=r3(i,j,1)**2+r3(i,j,2)**2
enddo
enddo



do i=1,3
do j=1,3
do a=1,2
   if(r2(i,j).gt.1d-20)then
      A3(i,j,a) = (1d0+epsilon_tilde*(r3(i,j,1)**4/r2(i,j)**2 + r3(i,j,2)**4/r2(i,j)**2))*r3(i,j,a) 
   else
      A3(i,j,a) = r3(i, j, a)
   end if
enddo
enddo
enddo

A_grad = 0d0

do i=2,3
do j=1,i-1
do a=1,2
   A_grad = A_grad + eps22(i,j)*0.5*A3(i,j,a)**2 
enddo
enddo
enddo
end function A_grad


DOUBLE PRECISION FUNCTION A_GRAD_phi(phi, phid, grad_phi, a_grad)
  USE SOLUTION_PARAMETERS
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: phi(3), grad_phi(3, 2)
  DOUBLE PRECISION, INTENT(IN) :: phid(3)
  DOUBLE PRECISION :: r3(3, 3, 2), a3(3, 3, 2), r2(3, 3)
  DOUBLE PRECISION :: r3d(3, 3, 2), a3d(3, 3, 2), r2d(3, 3)
  INTEGER :: i, j, a
  REAL :: result1
  DOUBLE PRECISION :: a_grad

  r3d = 0.D0
  DO i=1,3
    DO j=1,3
      DO a=1,2
        r3d(i, j, a) = grad_phi(j, a)*phid(i) - grad_phi(i, a)*phid(j)
        r3(i, j, a) = phi(i)*grad_phi(j, a) - phi(j)*grad_phi(i, a)
      END DO
    END DO
  END DO
  r2d = 0.D0
  DO i=1,3
    DO j=1,3
      r2d(i, j) = 2*r3(i, j, 1)*r3d(i, j, 1) + 2*r3(i, j, 2)*r3d(i, j, 2&
&       )
      r2(i, j) = r3(i, j, 1)**2 + r3(i, j, 2)**2
    END DO
  END DO
  a3d = 0.D0
  DO i=1,3
    DO j=1,3
      DO a=1,2
        IF (r2(i, j) .GT. 1d-20) THEN
          a3d(i, j, a) = epsilon_tilde*((4*r3(i, j, 1)**3*r3d(i, j, 1)*&
&           r2(i, j)**2-r3(i, j, 1)**4*2*r2(i, j)*r2d(i, j))/(r2(i, j)**&
&           2)**2+(4*r3(i, j, 2)**3*r3d(i, j, 2)*r2(i, j)**2-r3(i, j, 2)&
&           **4*2*r2(i, j)*r2d(i, j))/(r2(i, j)**2)**2)*r3(i, j, a) + (&
&           1d0+epsilon_tilde*(r3(i, j, 1)**4/r2(i, j)**2+r3(i, j, 2)**4&
&           /r2(i, j)**2))*r3d(i, j, a)
          a3(i, j, a) = (1d0+epsilon_tilde*(r3(i, j, 1)**4/r2(i, j)**2+&
&           r3(i, j, 2)**4/r2(i, j)**2))*r3(i, j, a)
        ELSE
          a3d(i, j, a) = r3d(i, j, a)
          a3(i, j, a) = r3(i, j, a)
        END IF
      END DO
    END DO
  END DO
  a_grad = 0d0
  a_grad_phi = 0.D0
  DO i=2,3
    DO j=1,i-1
      DO a=1,2
        result1 = EPS22(i, j)
        a_grad_phi = a_grad_phi + result1*0.5*2*a3(i, j, a)*a3d(i, j, a)
        a_grad = a_grad + result1*0.5*a3(i, j, a)**2
      END DO
    END DO
  END DO
END FUNCTION A_GRAD_phi

DOUBLE PRECISION FUNCTION A_GRAD_D_phi(phi, grad_phi, grad_phid, a_grad)
  USE SOLUTION_PARAMETERS
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: phi(3), grad_phi(3, 2)
  DOUBLE PRECISION, INTENT(IN) :: grad_phid(3, 2)
  DOUBLE PRECISION :: r3(3, 3, 2), a3(3, 3, 2), r2(3, 3)
  DOUBLE PRECISION :: r3d(3, 3, 2), a3d(3, 3, 2), r2d(3, 3)
  INTEGER :: i, j, a
  REAL :: result1
  DOUBLE PRECISION :: a_grad

  r3d = 0.D0
  DO i=1,3
    DO j=1,3
      DO a=1,2
        r3d(i, j, a) = phi(i)*grad_phid(j, a) - phi(j)*grad_phid(i, a)
        r3(i, j, a) = phi(i)*grad_phi(j, a) - phi(j)*grad_phi(i, a)
      END DO
    END DO
  END DO
  r2d = 0.D0
  DO i=1,3
    DO j=1,3
      r2d(i, j) = 2*r3(i, j, 1)*r3d(i, j, 1) + 2*r3(i, j, 2)*r3d(i, j, 2&
&       )
      r2(i, j) = r3(i, j, 1)**2 + r3(i, j, 2)**2
    END DO
  END DO
  a3d = 0.D0
  DO i=1,3
    DO j=1,3
      DO a=1,2
        IF (r2(i, j) .GT. 1d-20) THEN
          a3d(i, j, a) = epsilon_tilde*((4*r3(i, j, 1)**3*r3d(i, j, 1)*&
&           r2(i, j)**2-r3(i, j, 1)**4*2*r2(i, j)*r2d(i, j))/(r2(i, j)**&
&           2)**2+(4*r3(i, j, 2)**3*r3d(i, j, 2)*r2(i, j)**2-r3(i, j, 2)&
&           **4*2*r2(i, j)*r2d(i, j))/(r2(i, j)**2)**2)*r3(i, j, a) + (&
&           1d0+epsilon_tilde*(r3(i, j, 1)**4/r2(i, j)**2+r3(i, j, 2)**4&
&           /r2(i, j)**2))*r3d(i, j, a)
          a3(i, j, a) = (1d0+epsilon_tilde*(r3(i, j, 1)**4/r2(i, j)**2+&
&           r3(i, j, 2)**4/r2(i, j)**2))*r3(i, j, a)
        ELSE
          a3d(i, j, a) =r3d(i,j,a)
          a3(i, j, a) = r3(i,j,a)
        END IF
      END DO
    END DO
  END DO
  a_grad = 0d0
  a_grad_d_phi = 0.D0
  DO i=2,3
    DO j=1,i-1
      DO a=1,2
        result1 = EPS22(i, j)
        a_grad_d_phi = a_grad_d_phi + result1*0.5*2*a3(i, j, a)*a3d(i, j, a)
        a_grad = a_grad + result1*0.5*a3(i, j, a)**2
      END DO
    END DO
  END DO
END FUNCTION A_GRAD_D_phi


double precision function GradientEnergy_FV(LapPhi,vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::LapPhi(3),dx,vbles(5,3,3)
integer, intent (in):: lp
double precision ::phi_x(3),phi_y(3),phi(3),phi_xx(3),phi_yy(3),grad_phi(3,2)
double precision :: A_grad_phi,A_grad_d_phi,phid(3),A_grad,grad_phid(3,2)
integer :: ii,kk


!A_GRAD_phi(phi, phid, grad_phi, a_grad)
!A_GRAD_D_phi(phi, grad_phi, grad_phid, a_grad)
GradientEnergy_FV = 0.






!if(epsilon_tilde.lt.-1d-6)then
!isotropic

   phi(1) = vbles(1,2,2)
   phi(2) = vbles(2,2,2)
   phi(3) = vbles(3,2,2)

   phi_x(1) =  (vbles(1,3,2)-vbles(1,1,2))*0.5/dx
   phi_x(2) =  (vbles(2,3,2)-vbles(2,1,2))*0.5/dx
   phi_x(3) =  (vbles(3,3,2)-vbles(3,1,2))*0.5/dx
   
   phi_y(1) =  (vbles(1,2,3)-vbles(1,2,1))*0.5/dx
   phi_y(2) =  (vbles(2,2,3)-vbles(2,2,1))*0.5/dx
   phi_y(3) =  (vbles(3,2,3)-vbles(3,2,1))*0.5/dx
   
   
   do ii=1,N_phases
      if(ii.ne.lp)then
         GradientEnergy_FV = GradientEnergy_FV + eps22(lp,ii)*(2.*((phi(lp)*phi_x(ii)-phi(ii)*phi_x(lp))*phi_x(ii)&
              + (phi(lp)*phi_y(ii)-phi(ii)*phi_y(lp))*phi_y(ii)) -&
              (phi(ii)*LapPhi(lp)-phi(lp)*LapPhi(ii))*phi(ii))
      endif
   enddo

   GradientEnergy_FV = -A_0*GradientEnergy_FV
   if(vbles(lp,2,2).gt.0.4.and.vbles(lp,2,2).lt.0.6)then
      write(*,*)GradientEnergy_FV
   endif
!else
!anisotropic

   GradientEnergy_FV = 0.

   phid = 0d0
   phid(lp) = 1d0

   
   do ii=1,3
      phi(ii)=vbles(ii,2,2)
      grad_phi(ii,1) =  (vbles(ii,3,2)-vbles(ii,1,2))*0.5/dx
      grad_phi(ii,2) =  (vbles(ii,2,3)-vbles(ii,2,1))*0.5/dx
   enddo

   GradientEnergy_FV = GradientEnergy_FV + A_GRAD_phi(phi, phid, grad_phi, a_grad)
   
   
   grad_phid = 0d0
   grad_phid(lp,1) = 1d0 !x component of phi(lp)

   ! right point((1/2,0)   
   do ii=1,3
      phi(ii) = (vbles(ii,2,2)+vbles(ii,3,2))*0.5
      grad_phi(ii,1) =  (vbles(ii,3,2)-vbles(ii,2,2))/dx
      grad_phi(ii,2) =  (vbles(ii,2,3)-vbles(ii,2,1)+vbles(ii,3,3)-vbles(ii,3,1))*0.25/dx
   enddo

   GradientEnergy_FV = GradientEnergy_FV - A_GRAD_D_phi(phi, grad_phi, grad_phid, a_grad)/dx



   !left point(-1/2,0)   
   do ii=1,3
      phi(ii) = (vbles(ii,2,2)+vbles(ii,1,2))*0.5
      grad_phi(ii,1) =  (vbles(ii,2,2)-vbles(ii,1,2))/dx
      grad_phi(ii,2) =  (vbles(ii,2,3)-vbles(ii,2,1)+vbles(ii,1,3)-vbles(ii,1,1))*0.25/dx
   enddo

   GradientEnergy_FV = GradientEnergy_FV + A_GRAD_D_phi(phi, grad_phi, grad_phid, a_grad)/dx
   
   grad_phid = 0d0
   grad_phid(lp,2) = 1d0 !y component of phi(lp)

   !up point(0,1/2)   
   do ii=1,3
      phi(ii) = (vbles(ii,2,2)+vbles(ii,2,3))*0.5
      grad_phi(ii,1) =  (vbles(ii,2,3)-vbles(ii,2,2))/dx
      grad_phi(ii,2) =  (vbles(ii,3,2)-vbles(ii,1,2)+vbles(ii,3,3)-vbles(ii,1,3))*0.25/dx
   enddo

   GradientEnergy_FV = GradientEnergy_FV - A_GRAD_D_phi(phi, grad_phi, grad_phid, a_grad)/dx
   
   !down point (0,-1/2)   
   do ii=1,3
      phi(ii) = (vbles(ii,2,2)+vbles(ii,2,1))*0.5
      grad_phi(ii,1) =  (vbles(ii,2,2)-vbles(ii,2,1))/dx
      grad_phi(ii,2) =  (vbles(ii,3,2)-vbles(ii,1,2)+vbles(ii,3,1)-vbles(ii,1,1))*0.25/dx
   enddo

   GradientEnergy_FV = GradientEnergy_FV + A_GRAD_D_phi(phi, grad_phi, grad_phid, a_grad)/dx
!end if
GradientEnergy_FV = -A_0*GradientEnergy_FV
   if(vbles(lp,2,2).gt.0.4.and.vbles(lp,2,2).lt.0.6)then
      write(*,*)GradientEnergy_FV
      write(*,*)"---------------------------"
   endif





end function GradientEnergy_FV




double precision function A_func(X,Y)
double precision, intent(in) :: X,Y
double precision :: epsilon_tilde=0d0
A_func = 1d0+epsilon_tilde
if (X+Y.lt.1d-11)return
A_func = 1d0 + epsilon_tilde*(X**2 + Y**2)/(X+Y)**2
end function A_func




!T_K converts non-dimensional parameter T to dimensioned T_K(T)
double precision function T_K(T)
implicit none
double precision, intent (in)::T
double precision :: T_eut= 456.14
!Alternate non-dimensionalisation
! double precision T_b,C_b,L_b, 
!  T_b = 505.07
!  C_b = 210.
!  L_b = 7029.
!  T_K = L_b/C_b*T+T_b
  
 T_eut = 456.14
 
 T_K = T_eut*(1.+T)

 end function T_K
 double precision function TemperatureRHS(vbles,dx,c,T,phi,c_dot,phi_dot)
!This is not the most efficient routine. It takes Div (D*c*(1-c)Grad(h)) by
!examining its neighboring values. Quicker is to evaluate all points in
!one sweep then use them for this routine. But that approach would have to be !done outside this routine.
implicit none
double precision, intent (in):: dx,c,T,phi(3),c_dot,phi_dot(3),vbles(5,3,3)
!integer, intent (in)::df
double precision, dimension(3,3)::stencil =reshape((/0.,1.,0.,1.,-4.,1.,0.,1.,0./),(/3,3/))
double precision, dimension(3)::stencil_x =(/-.5,.0,.5/)
double precision:: FreeEnergy,x,a,xa,dFdc(3,3) !x coordinate
double precision:: Cp= 28.5!29.9758394     !J/mol/K
! double precision:: D=2.110795e-5      !m^2/s
double precision:: T_eut=456.14  !K
! double precision ::D_c = 2.976894361E-9 !m^2/s characteristic diffusivity
! double precision ::D_ch = 1.126392460619016E-011 !hunt Jackson
double precision ::D_ch = 2.83e-10,R = 8.31 !J/(K mol) !hunt Jackson
double precision, dimension(3):: D1=(/6.7e-10,6.7e-14,6.7e-14/)!Jackson Hunt
! double precision :: D_tilde = 7090.594237 !dimensionless diffusivity
double precision :: D_tilde = 2630.6602369, D_heat_c !dimensionless debug Le=10.
double precision :: T_,c_,phi_(3)
logical :: heat_c_term=.true.
integer :: ii,jj,kk, N_phases=3
!  rhoa(1)=10660. !liquid Lead kg/m^3
!  rhoa(2)=11340. !Solid Lead
!  rhoa(3)=11340. !The same (for want of a better option)
!  rhob(1)=6999.  !liquid Tin
!  rhob(2)=7365.  !Solid Tin (White)
!  rhob(3)=7365.  !The same

!  rho=0.25*(11340.+ 7365.+6999.+10660.)!kg/m^3 keep it simple
!  kappa =35.3 !J/(s m)
  TemperatureRHS=0.
  return
  D_tilde = 40.
!    return !delete for isothermal

 TemperatureRHS = FreeEnergy(c,T,phi,c_dot,phi_dot,N_phases+1) !J K/mol
!  TemperatureRHS = TemperatureRHS/Cp !means LHS is T_dot only
! heat capacity is in FreeEnergy function now Cp = -T(d/dT)dF/dT
 TemperatureRHS = TemperatureRHS/T_eut !non dimensional T
 TemperatureRHS = TemperatureRHS
 
 do ii=1,3
 do jj=1,3
!  T_=unk(1+N_phases*nunkvbles,i-2+ii,j-2+jj,k,lb)*stencil(ii,jj)
  T_=vbles(4,ii,jj)*stencil(ii,jj)
 TemperatureRHS = TemperatureRHS + D_tilde*T_/(dx*dx)
 enddo
 enddo

! replaces c_dot_term an heat_flux_term which does not appear to operate well
 if(heat_c_term)then
  do ii=1,3
 do jj=1,3
!   T_=unk(1+N_phases*nunkvbles,i-2+ii,j-2+jj,k,lb)
    T_=vbles(4,ii,jj)
!   c_=unk(1+(N_phases+1)*nunkvbles,i-2+ii,j-2+jj,k,lb)
    c_=vbles(5,ii,jj)
   do kk=1,N_phases
!    phi_(kk)=unk(1+(kk-1)*nunkvbles,i-2+ii,j-2+jj,k,lb)
      phi_(kk)=vbles(kk,ii,jj)
   enddo
   dFdc(ii,jj)=((1.+T)*(log(c)-log(1.-c))+FreeEnergy(c_,T_,Phi_,c_dot,phi_dot,N_phases+2))
 enddo
 enddo
 D_heat_c = 0.
 do ii=1,N_phases
  D_heat_c = D_heat_c + D1(ii)*phi(ii)
 enddo
 D_heat_c = (c*(1.-c)*D_heat_c/D_ch)/(Cp*T_eut)
 
!   grad(dfdc) dot grad(dFdc)
 TemperatureRHS = TemperatureRHS + D_heat_c*((dFdc(3,2)-dFdc(1,2))**2 + (dFdc(2,3)-dFdc(2,1))**2)/(4.*dx*dx)

 endif
!  TemperatureRHS = TemperatureRHS + v_mesh*(unk(1+3*nunkvbles,i+1,j,k,lb)-unk(1+3*nunkvbles,i,j,k,lb))/(dx)

 

 end function TemperatureRHS

!

double precision function SoluteRHS(vbles,dx,c_dot,phi_dot,c_c)
! this uses volume element method. This guarantees conservation of solute. The main source of solute is the variation of free energy with phase: d (d F/dc) = ...+ (d^2 F/dc d phi) d phi + ..., which is revealed by a plot of c_dot.

implicit none
integer :: N_phases=3
double precision, intent (in):: dx,c_dot,c_c,vbles(5,3,3)
!integer, intent (in)::i,j,k,lb,df
integer ii,jj,ip,ih,jv
double precision c,T,Phi(3),FreeEnergy,D,phi_dot(3),phi_,phi_tmp(3)
double precision, dimension(3,3)::stencil11=reshape((/0.,6.,0.,6.,-24.,6.,0.,6.,0./),(/3,3/))
double precision, dimension(3)::stencil1=(/-.5,0.,.5/)
double precision, dimension(3):: D1=(/6.7e-10,6.7e-12,6.7e-12/)!Jackson Hunt
! double precision, dimension(3):: D1=(/1e-9,1e-11,1e-11/)
double precision :: DFDc,potential,Lap_c,D_, D_c,D_cc,D_c1,D_cc1,Grad_c(2),y0,g_phi,f_c,f_cc
double precision ::D_ch =  2.83e-10! (m^2/s)    characteristic diffusivity delta2/tau=M*eps2(1)
double precision :: phaseRHS_c,beta=1d1
! 
beta=0.5
SoluteRHS=0.
D_ch =  2.83d-10
!
!if(df.eq.0)then
!x-coordinate transform values
!! Right
 ih=3
 jv=2
 include 'c_stencil.h'
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Left
 ih=1
 jv=2
 include 'c_stencil.h'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!up
 ih=2
 jv=3
 include 'c_stencil.h'
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !down
 ih=2
 jv=1
 include 'c_stencil.h'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SoluteRHS = SoluteRHS/D_ch/(dx*dx) 
! SoluteRHS = SoluteRHS+v_mesh*(unk(1+4*nunkvbles,i+1,j,k,lb)-unk(1+4*nunkvbles,i,j,k,lb))/(dx)
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

return

end function SoluteRHS

double precision function potential(phi,c,lp)

implicit none

double precision, intent (in)::phi(3),c
double precision, dimension(3,3)::sigma_a,sigma_b
double precision, dimension(3)::Wa,Wb
double precision, dimension(3,3)::beta_a,beta_b,delta_a,delta_b
! double precision ::R = 8.31,T_eut=456.14,MolPerVol= 54185.3329242099
double precision :: VolPerJoule = 4.87E-009 !=1/(R*T_eut*MolPerVol)
integer, intent(in):: lp
double precision t_a,t_b,DX,g_phi,bump,qq,eps2(3)
integer :: ip,jp, N_phases=3
bump=0.
potential =0.
qq=0.!-0.335
eps2(1) = 1.*(1.-qq)
eps2(2) = 1.*qq
eps2(3) = 1.*(0.77+qq)

Wa = eps2



!Folch-Plapp potential
 potential=0.
 if(lp.le.N_phases)then 
   potential =  Wa(lp)*((2.+(-6.+4.*phi(lp))*phi(lp))*phi(lp))
 endif

end function potential


double precision function PhaseRHS(vbles,dx,lp,Phi,c,T,c_dot,phi_dot)
use solution_parameters
implicit none
!integer, intent(in)::i,j,k,lb
double precision, intent(in)::dx,c,T,c_dot,phi_dot(3),Phi(3),vbles(5,3,3)
integer, intent(in):: lp
double precision :: Phi11(3,2),Phi12(3,2,2),LapPhi(3)
double precision GradientEnergy,potential,FreeEnergy,GE,dd_fi
double precision rf4,rf5,rf6,rfactor
double precision ::M_a= 0.948190 ,M_b= 0.5680384!m^3/(J s) mobility
double precision ::t_a = 600.61,t_b = 505.07
! double precision, dimension (3)::eps2=(/0d0,1d0,1.771771772d0/)
! double precision, dimension (3)::eps2=(/.5d0,.5d0,1d0/)
! double precision, dimension (3)::eps2=(/1.187939392393404E-011,2.706804758382105E-010,4.887522071561416E-010/)

integer ii,jj
! double precision ::W_a(3,3)
! double precision ::W_b = 5.7337864e7,xtmp,MolPerVol
! double precision ::eps2=3.139554108e-8! (J/m)

double precision ::W_a(3,3)
! double precision ::W_b(3,3)
! double precision ::eps2(3,3)
double precision ::M_tilde !s characteristic time
!double precision ::delta2 = 4.0000000E-016 !m^2 characteristic area (original)
double precision :: delta2 = 4E-18!decrease to sharpen. This i.c. must increase by similar factor and 
!D_solute <- D_solute/scale_factor
double precision Cross_term

integer ip,n

!!!!! Calculate characteristic time and space (squared)

! tau_a=1./(M_a*W_a)
 PhaseRHS=0.

 w_a(1,1)=7e7 ! see notes = 0.07 / delta = 0.07/1e-9 

 call Calc_GradPhi(LapPhi,vbles,dx)
!assemble and non-dimensionlise phaseRHS 

 M_tilde = 1.0-0.159*c ! see notes Anisotropic2d3d.tex !debug
!
  PhaseRHS=M_tilde*(GradientEnergy(LapPhi,vbles,lp,dx)-potential(Phi,c,lp)*scale_factor**2/lambda**2&
       -scale_factor**2/lambda*FreeEnergy(c,T,Phi,c_dot,phi_dot,lp)/W_a(1,1))
  if(lp.eq.1)PhaseRHS=PhaseRHS - Cross_term(dx,c_dot,phi_dot,vbles)
end function PhaseRHS

double precision function Cross_term(dx,c_dot,phi_dot,vbles)
use solution_parameters
implicit none
double precision, intent(in)::dx,c_dot,phi_dot(3),vbles(5,3,3)
double precision FreeEnergy,A,D,D_c,D_cc,SoluteRHS,f_c,f_cc,c_c,c,T,phi(3),potential
double precision, dimension(3):: D1=(/6.7e-10,6.7e-12,6.7e-12/)!Jackson Hunt
double precision ::beta=0d0,phaseRHS_c,beta_=0.5
integer ih,jv,ip
!x-coordinate transform values
!! Right
 ih=3
 jv=2
 include 'c_stencil.h'
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Left
 ih=1
 jv=2
 include 'c_stencil.h'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!up
 ih=2
 jv=3
 include 'c_stencil.h'
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !down
 ih=2
 jv=1
 include 'c_stencil.h'




Cross_term = soluteRHS*beta_/2.83d-10/(dx*dx)
return
end function Cross_term

double precision function PhaseRHS_c(vbles,dx,lp,Phi,c,T,c_dot,phi_dot)
use solution_parameters
implicit none
!integer, intent(in)::i,j,k,lb
double precision, intent(in)::dx,c,T,c_dot,phi_dot(3),Phi(3),vbles(5,3,3)
integer, intent(in):: lp
double precision :: Phi11(3,2),Phi12(3,2,2),LapPhi(3)
double precision GradientEnergy,potential,FreeEnergy,GE,dd_fi
double precision rf4,rf5,rf6,rfactor
double precision ::M_a= 0.948190 ,M_b= 0.5680384!m^3/(J s) mobility
double precision ::t_a = 600.61,t_b = 505.07
! double precision, dimension (3)::eps2=(/0d0,1d0,1.771771772d0/)
! double precision, dimension (3)::eps2=(/.5d0,.5d0,1d0/)
! double precision, dimension (3)::eps2=(/1.187939392393404E-011,2.706804758382105E-010,4.887522071561416E-010/)

integer ii,jj
! double precision ::W_a(3,3)
! double precision ::W_b = 5.7337864e7,xtmp,MolPerVol
! double precision ::eps2=3.139554108e-8! (J/m)

double precision ::W_a(3,3)
! double precision ::W_b(3,3)
! double precision ::eps2(3,3)
double precision ::M_tilde !s characteristic time
!double precision ::delta2 = 4.0000000E-016 !m^2 characteristic area (original)
double precision :: delta2 = 4E-18!decrease to sharpen. This i.c. must increase by similar factor and 
!D_solute <- D_solute/scale_factor


integer ip,n

!!!!! Calculate characteristic time and space (squared)

! tau_a=1./(M_a*W_a)
 PhaseRHS_c=0.

 w_a(1,1)=7e7 ! see notes = 0.07 / delta = 0.07/1e-9 

 call Calc_GradPhi(LapPhi,vbles,dx)
!assemble and non-dimensionlise phaseRHS 


!
  PhaseRHS_c=-(lambda**2*GradientEnergy(LapPhi,vbles,lp,dx)-potential(Phi,c,lp)*scale_factor**2&
       -scale_factor**2*lambda*FreeEnergy(c,T,Phi,c_dot,phi_dot,lp)/W_a(1,1))

end function PhaseRHS_c


subroutine Project_Phase(vbles,dx,phi_rhs,Phi,c,T,c_dot,phi_dot)
implicit none
!integer, intent(in)::df
double precision, intent(in)::dx,Phi(3),c,T,c_dot,phi_dot(3),vbles(5,3,3)
double precision, intent(out) :: phi_rhs(3)
integer :: N_phases=3
double precision PhaseRHS,Wk(3), dWk(5,5),H,e(5,5),tip,Wk_T(3),Wk_S(3),pphi(3),tr
integer ip,jp,kp,lp,ap,bp,izero,iflag
double precision, dimension(3,3)::I2 = reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))
double precision, dimension(3,3)::Projector = 0.
double precision, dimension(3,3,3)::D_Projector=0. 
double precision, dimension(3)::q_factor=(/-2e-12,1e-12,1e-12/)
double precision, dimension(3,3)::L2 = reshape((/2.,-1.,-1.,-1.,2.,-1.,-1.,-1.,2./),(/3,3/))
!  perturbs growing tip
!   tip = 16.*unk(1, i, j, k, lb)**2*(1.-unk(1, i, j, k, lb))**2
!   tip = tip*(dmod(time*tip,exp(1.)*1.3e-10)/exp(1.)/1.3e-10-5d-1)
  phi_rhs=0.
!assemble Projector
  izero=0
 do ip=1,N_phases
  if(1d0-Phi(ip).lt.1d-16)izero=ip
 enddo
 D_Projector=0.
 
 
  if(izero.eq.0)then

! maple generated A
Projector(1,1)  =  phi(1) + phi(2) * phi(1) ** 2 / (1 - phi(2))** 2 &
     + phi(3) * phi(1) ** 2 / (1 - phi(3)) ** 2
Projector(1,2)  =  -phi(1) / (1 - phi(1)) * phi(2) - phi(2) * phi(1)&
     / (1 - phi(2)) + phi(3) * phi(1) * phi(2) / (1 - phi(3)) ** 2
Projector(1,3)  =  -phi(1) / (1 - phi(1)) * phi(3) + phi(2) * phi(1) &
     * phi(3) / (1 - phi(2)) ** 2 - phi(3) * phi(1) / (1 - phi(3))
Projector(2,1)  =  -phi(1) / (1 - phi(1)) * phi(2) - phi(2) * phi(1) &
     / (1 - phi(2)) + phi(3) * phi(1) * phi(2) / (1 - phi(3)) ** 2
Projector(2,2)  =  phi(1) * phi(2) ** 2 / (1 - phi(1)) ** 2 + phi(2)&
     + phi(3) * phi(2) ** 2 / (1 - phi(3)) ** 2
Projector(2,3)  =  phi(1) * phi(2) * phi(3) / (1 - phi(1)) ** 2 - phi(2) &
     / (1 - phi(2)) * phi(3) - phi(3) * phi(2) / (1 - phi(3))
Projector(3,1)  =  -phi(1) / (1 - phi(1)) * phi(3) + phi(2) * phi(1) *&
     phi(3) / (1 - phi(2)) ** 2 - phi(3) * phi(1) / (1 - phi(3))
Projector(3,2)  =  phi(1) * phi(2) * phi(3) / (1 - phi(1)) ** 2 - phi(2)&
     / (1 - phi(2)) * phi(3) - phi(3) * phi(2) / (1 - phi(3))
Projector(3,3)  =  phi(1) * phi(3) ** 2 / (1 - phi(1)) ** 2 + phi(2) *&
     phi(3) ** 2 / (1 - phi(2)) ** 2 + phi(3)
! model B maple generated

  else
!     Projector=0.
   do kp=1,N_phases
   do lp=1,N_phases
    Projector(kp,lp)=2.*I2(kp,lp) - 1.
   enddo
   enddo
  endif

    
  
  do ip=1,N_phases
   Wk(ip) =  PhaseRHS(vbles,dx,ip,Phi,c,T,c_dot,phi_dot) 
  enddo 

  
!Project phi_vector
  do ip=1,N_phases
   phi_rhs(ip)=0.!q_factor(ip)*tip*dt
   do jp=1,N_phases
    phi_rhs(ip)=phi_rhs(ip)+Projector(ip,jp)*Wk(jp)
   enddo
  enddo
!   do ip=1,N_phases
!     phi_rhs(ip)=phi_rhs(ip)+v_mesh*(unk(1+(ip-1)*nunkvbles,i+1,j,k,lb)-unk(1+(ip-1)*nunkvbles,i,j,k,lb))/(dx) 
!   enddo
 
end subroutine Project_Phase





subroutine get_RHS_Fi_2d(i,j,k,lb,dx,Fi)
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none

integer, intent(in)::i,j,k,lb
double precision, intent(in)::dx
double precision, intent(out):: Fi(5)
double precision :: c,T,phi(3),c_dot,phi_dot(3)
double precision phi_rhs(3),dphi_rhs(3),rfactor,rf4,rf5,rf6,S
double precision TemperatureRHS,SoluteRHS,T_k,phi_,phi_tmp(3)
double precision vbles(5,3,3),unk_(30)
logical Not_number
integer ip,n,ii,jj,kk
  Fi=0.

  do ii=1,3
     do jj=1,3
        do ip=1,N_phases+2
           vbles(ip,ii,jj)=unk(1+(ip-1)*nunkvbles,i+ii-2,j+jj-2,k,lb)
        enddo
     enddo
  enddo
do ii=1,30
 unk_(ii) = unk(ii,i,j,k,lb)
 enddo
 call get_RHS_Fi(vbles,dx,Fi,mode_1,N_phases,nunkvbles,dt,dtold,lnblocks,unk_)
 do ip=1,N_phases+2
    unk(1+(ip-1)*nunkvbles,i,j,k,lb)=vbles(ip,2,2)
 enddo



 
end subroutine get_RHS_Fi_2d

subroutine get_RHS_Fi(vbles,dx,Fi,mode_1,N_phases,nunkvbles,dt,dtold,lnblocks,unk_)
  implicit none

integer, intent(in)::N_phases,mode_1,nunkvbles,lnblocks
double precision, intent(in)::dx,dt,dtold,unk_(30)
double precision, intent(out):: Fi(5)
double precision :: c,T,phi(3),c_dot,phi_dot(3)
double precision phi_rhs(3),dphi_rhs(3),rfactor,rf4,rf5,rf6,S
double precision TemperatureRHS,SoluteRHS,T_k,phi_,phi_tmp(3)
double precision vbles(5,3,3)
logical Not_number
integer ip,n,ii,jj,kk
  if (mode_1.eq.1) then

   do ip=1,N_phases
    n=(ip-1)*nunkvbles
!    phi_dot(ip) = (unk (1+n, i, j, k, lb)-unk (n+2, i, j, k, lb))/dt
    phi_dot(ip) = (vbles(ip,2,2)-unk_(n+2))/dt

!     phi_dot(ip) = phi_dot(ip) - v_mesh*(unk(1+n,i+1,j,k,lb)-unk(1+n,i,j,k,lb))/(dx)
   enddo
   n=(N_phases+1)*nunkvbles
!   c_dot = (unk (1+n, i, j, k, lb)-unk (2+n, i, j, k, lb))/dt
   c_dot = (vbles(5,2,2)-unk_(2+n))/dt
!    c_dot = c_dot - v_mesh*(unk(1+4*nunkvbles,i+1,j,k,lb)-unk(1+4*nunkvbles,i,j,k,lb))/(dx)
  elseif(mode_1.eq.2)then
   rfactor = dt/dtold
   rf4 = rfactor*rfactor/(rfactor+1.0)
   rf5 = (rfactor+1.0)
   rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
 
   do ip=1,N_phases
   n=(ip-1)*nunkvbles
!   phi_dot(ip) = (rf4*unk(5+n,i,j,k,lb)-rf5*unk(4+n,i,j,k,lb)+rf6*unk (1+n, i, j, k, lb))/dt
   phi_dot(ip) = (rf4*unk_(5+n)-rf5*unk_(4+n)+rf6*vbles(ip,2,2))/dt
!     phi_dot(ip) = phi_dot(ip) - v_mesh*(unk(1+n,i+1,j,k,lb)-unk(1+n,i,j,k,lb))/(dx)
   enddo
   n=(N_phases+1)*nunkvbles
!   c_dot = (rf4*unk(5+n,i,j,k,lb)-rf5*unk(4+n,i,j,k,lb)+rf6*unk (1+n, i, j, k, lb))/dt
   c_dot = (rf4*unk_(5+n)-rf5*unk_(4+n)+rf6*vbles(5,2,2))/dt
!    c_dot = c_dot - v_mesh*(unk(1+4*nunkvbles,i+1,j,k,lb)-unk(1+4*nunkvbles,i,j,k,lb))/(dx)
  else
   write(*,*)"mode_1 not = 1 or 2"
   stop
  endif 
!   S=unk (1, i, j, k, lb)+unk (1+nunkvbles, i, j, k, lb)+unk (1+2*nunkvbles, i, j, k, lb)
  do ip=1,N_phases
   n=(ip-1)*nunkvbles
!   phi(ip)=max(0d0,unk (1+n, i, j, k, lb))
   phi(ip)=max(0d0,vbles(ip,2,2))
  enddo
  do ip=1,N_phases
     n=(ip-1)*nunkvbles
     phi(ip)=min(1d0,phi(ip)/(phi(1)+phi(2)+phi(3)))
!     unk (1+n, i, j, k, lb)=phi(ip)
     vbles(ip,2,2)=phi(ip)
  enddo
 
!   if(phi(1)+phi(2)+phi(3)-1..gt.1d0)then
!   write(*,*) phi(1)+phi(2)+phi(3)-1d0
!   write(*,*)"in get_RHS_Fi_2d"
!   stop
!   endif
  T = vbles(4,2,2)!debug pcb
  c = min(1d0-1d-10,vbles(5,2,2))
  
!   if(phi(2)+phi(3).lt.1d-14)then
!   if(abs(c-0.74).le.1d-6)return
!   endif
  
!  unk(1+(N_phases+1)*nunkvbles, i, j, k, lb)=c
  vbles(5,2,2)=c
!  call Project_Phase( i,j,k,lb,dx,phi_rhs,dphi_rhs,phi,c,T,c_dot,phi_dot,0)
  call Project_Phase( vbles,dx,phi_rhs,phi,c,T,c_dot,phi_dot)
  do ip=1,N_phases
   Fi(ip)=phi_rhs(ip)
  enddo
!  Fi(N_phases+1) = TemperatureRHS(i,j,k,lb,dx,c,T,phi,c_dot,phi_dot,0)
  Fi(N_phases+1) = TemperatureRHS(vbles,dx,c,T,phi,c_dot,phi_dot)
!  Fi(N_phases+2) = SoluteRHS(i,j,k,lb,dx,c_dot,phi_dot,c,0)
  Fi(N_phases+2) = SoluteRHS(vbles,dx,c_dot,phi_dot,c)
 
!  if(Not_number(Fi(5)))then
!   write(*,*)Fi,i,j,lb
!   write(*,*)"phi=", phi
!   write(*,*)"dx,c,T,etc",dx,c,T,c_dot,phi_dot
!   write(*,*)unk(1,i,j,1,lb),unk(7,i,j,1,lb),unk(13,i,j,1,lb)
!   write(*,*)"T and c =",unk(19,i,j,1,lb),unk(25,i,j,1,lb)
!   write(*,*)"Temperature=",T_k(T),T
!   write(*,*)"===================================="
!   write(*,*)"start",i,j,k,lb,dx,c_dot,phi_dot,c,"end"
!   write(*,*)"SoluteRHS=",SoluteRHS(i,j,k,lb,dx,c_dot,phi_dot,c,0)
!   write(*,*)"===================================="
!  do ii=1,10
!  do jj=1,10
!!  do kk=1,lnblocks
!   c=vbles(5,2,2)
!   if(Not_number(c))then 
!   write(*,*)ii,jj,kk,c
!   write(*,*)"not number"
!   endif
!  enddo
!  enddo
!  enddo
!  write(*,*)"in  get_RHS_Fi_2, trap",nunkvbles*(N_phases+2)
!  write(*,'(5e12.5)') Fi
!  do ii=1,5
!     write(*,'(3e12.5)') vbles(ii,:,:)
!     write(*,*)"=================================="
!  enddo
!  stop
!  end if
end subroutine get_RHS_Fi

double precision function GradientEnergy(LapPhi,vbles,lp,dx)

implicit none
double precision, intent (in)::LapPhi(3),dx,vbles(5,3,3)
integer, intent (in):: lp
logical :: folch=.true.
double precision ::phi_x(3),phi_y(3),phi(3),phi_xx(3),phi_yy(3),eps22(3,3)
integer :: ii,kk, N_phases=3

  eps22(1,2) = 1.0 
  eps22(1,3) = 1.77
  eps22(3,2) = 0.77
  eps22(2,1) = 1.0
  eps22(3,1) = 1.77
  eps22(2,3) = 0.77



  GradientEnergy = 0.
!if(folch)then
!   GradientEnergy = LapPhi(lp)
!  return
! endif
!if(df.eq.0)then

phi(1) = vbles(1,2,2)
phi(2) = vbles(2,2,2)
phi(3) = vbles(3,2,2)

phi_x(1) =  (vbles(1,3,2)-vbles(1,1,2))*0.5/dx
phi_x(2) =  (vbles(2,3,2)-vbles(2,1,2))*0.5/dx
phi_x(3) =  (vbles(3,3,2)-vbles(3,1,2))*0.5/dx

phi_y(1) =  (vbles(1,2,3)-vbles(1,2,1))*0.5/dx
phi_y(2) =  (vbles(2,2,3)-vbles(2,2,1))*0.5/dx
phi_y(3) =  (vbles(3,2,3)-vbles(3,2,1))*0.5/dx


do ii=1,N_phases
 if(ii.ne.lp)then
  GradientEnergy = GradientEnergy + eps22(lp,ii)*(2.*((phi(lp)*phi_x(ii)-phi(ii)*phi_x(lp))*phi_x(ii)&
       + (phi(lp)*phi_y(ii)-phi(ii)*phi_y(lp))*phi_y(ii)) -&
       (phi(ii)*LapPhi(lp)-phi(lp)*LapPhi(ii))*phi(ii))
 endif
enddo

GradientEnergy = -GradientEnergy
end function GradientEnergy

! for now only calculates LapPhi
subroutine Calc_GradPhi(LapPhi,vbles,dx)

implicit none
!integer, intent(in)::i,j,k,lb
double precision, intent(in)::dx,vbles(5,3,3)
double precision, intent(out)::LapPhi(3)
double precision, dimension(3,3)::stencil=reshape((/1.,4.,1.,4.,-20.,4.,1.,4.,1./),(/3,3/)),&
     stencil5=reshape((/0.,1.,0.,1.,-4.,1.,0.,1.,0./),(/3,3/))

double precision, dimension(3,3)::stencilyy =reshape((/0.,1.,0.,0.,-2.,0.,0.,1.,0./),(/3,3/))
double precision, dimension(3,3)::stencilxx =reshape((/0.,0.,0.,1.,-2.,1.,0.,0.,0./),(/3,3/))
double precision, dimension(3,3)::stencilx =reshape((/0.,0.,0.,-.5,0.,.5,0.,0.,0./),(/3,3/))
integer ip,n,ii,jj
double precision dphi_r,dphi_l,x,a,xa
double precision XX,YY, phi_x,phi_y,A_func,A_
logical Not_number
logical :: nine_point=.true.
integer :: N_phases=3
double precision :: epsilon_tilde=0d0 
 
 LapPhi=0.
 do ip=1,N_phases
  if(epsilon_tilde.gt.1d-10)then
   phi_x = 0.5*(vbles(ip,3,2) - vbles(ip,1,2))/dx
   phi_y = 0.5*(vbles(ip,2,3) - vbles(ip,2,1))/dx
   XX = phi_x**2
   YY = phi_y**2
   A_=A_func(XX,YY)
  else
   XX=0.
   YY=0.
   A_=1.0
  endif
  if(XX+YY.lt.1d-8)then
   do ii=1,3
   do jj=1,3
    if(nine_point)then
      LapPhi(ip)=LapPhi(ip)+stencil(ii,jj)*vbles(ip,ii,jj)
    else
      LapPhi(ip)=LapPhi(ip)+stencil5(ii,jj)*vbles(ip,ii,jj)
    endif
   enddo
   enddo
   if(nine_point)then
    LapPhi(ip)=A_*A_*LapPhi(ip)/(6.*dx*dx)
   else
    LapPhi(ip)=A_*A_*LapPhi(ip)/(dx*dx)
   endif
  else
!  write(*,*) "in Calc_GradPhi. No anisotropy"
!  stop
!   finite volume method

! x right
   phi_x = (vbles(ip,3,2) - vbles(ip,2,2))/dx
   phi_y = 0.25*(vbles(ip,2,3) - vbles(ip,2,1) + vbles(ip,3,3) - vbles(ip,3,1))/dx
   XX = phi_x**2
   YY = phi_y**2
   A_=A_func(XX,YY)
   LapPhi(ip) = LapPhi(ip) +A_*(phi_x*(epsilon_tilde*2.*YY*(XX-YY)/(XX+YY)**2+A_))/dx
! x left
   phi_x = (vbles(ip,1,2) -vbles(ip,2,2))/dx
   phi_y = 0.25*(vbles(ip,2,3) - vbles(ip,2,1) + vbles(ip,1,3) - vbles(ip,1,1))/dx
   XX = phi_x**2
   YY = phi_y**2
   A_=A_func(XX,YY)
   LapPhi(ip) = LapPhi(ip) +A_*(phi_x*(epsilon_tilde*2.*YY*(XX-YY)/(XX+YY)**2+A_))/dx
! y up
   phi_y = (vbles(ip,2,3) -vbles(ip,2,2))/dx
   phi_x = 0.25*(vbles(ip,3,2) - vbles(ip,1,2) + vbles(ip,3,3) - vbles(ip,1,3))/dx
   XX = phi_x**2
   YY = phi_y**2
   A_=A_func(XX,YY)
   LapPhi(ip) = LapPhi(ip) +A_*(phi_y*(epsilon_tilde*2.*XX*(YY-XX)/(XX+YY)**2+A_))/dx
! y down
   phi_y = (vbles(ip,2,1) -vbles(ip,2,2))/dx
   phi_x = 0.25*(vbles(ip,3,2) - vbles(ip,1,2) + vbles(ip,3,1) - vbles(ip,1,1))/dx
   XX = phi_x**2
   YY = phi_y**2
   A_=A_func(XX,YY)
   LapPhi(ip) = LapPhi(ip) +A_*(phi_y*(epsilon_tilde*2.*XX*(YY-XX)/(XX+YY)**2+A_))/dx
!    these terms to reduce grid anisotropy
   if(nine_point)then
!    top right
   phi_x = 0.5*(vbles(ip,3,3) - vbles(ip,2,3) + vbles(ip,3,2) -vbles(ip,2,2))/dx
   phi_y = 0.5*(vbles(ip,3,3) - vbles(ip,3,2) + vbles(ip,2,3) -vbles(ip,2,2))/dx
   XX = phi_x**2
   YY = phi_y**2
   A_=A_func(XX,YY)
   LapPhi(ip) = LapPhi(ip) +0.5*A_*(phi_x*(epsilon_tilde*2.*YY*(XX-YY)/(XX+YY)**2+A_))/dx
   LapPhi(ip) = LapPhi(ip) +0.5*A_*(phi_y*(epsilon_tilde*2.*XX*(YY-XX)/(XX+YY)**2+A_))/dx
!    top left
   phi_x = 0.5*(vbles(ip,1,3) - vbles(ip,2,3) + vbles(ip,1,2) -vbles(ip,2,2))/dx
   phi_y = 0.5*(vbles(ip,1,3) - vbles(ip,1,2) + vbles(ip,2,3) -vbles(ip,2,2))/dx
   XX = phi_x**2
   YY = phi_y**2
   A_=A_func(XX,YY)
   LapPhi(ip) = LapPhi(ip) +0.5*A_*(phi_x*(epsilon_tilde*2.*YY*(XX-YY)/(XX+YY)**2+A_))/dx
   LapPhi(ip) = LapPhi(ip) +0.5*A_*(phi_y*(epsilon_tilde*2.*XX*(YY-XX)/(XX+YY)**2+A_))/dx
!    bottom left
   phi_x = 0.5*(vbles(ip,1,1) - vbles(ip,2,1) + vbles(ip,1,2) -vbles(ip,2,2))/dx
   phi_y = 0.5*(vbles(ip,1,1) - vbles(ip,1,2) + vbles(ip,2,1) -vbles(ip,2,2))/dx
   XX = phi_x**2
   YY = phi_y**2
   A_=A_func(XX,YY)
   LapPhi(ip) = LapPhi(ip) +0.5*A_*(phi_x*(epsilon_tilde*2.*YY*(XX-YY)/(XX+YY)**2+A_))/dx
   LapPhi(ip) = LapPhi(ip) +0.5*A_*(phi_y*(epsilon_tilde*2.*XX*(YY-XX)/(XX+YY)**2+A_))/dx
!    bottom right
   phi_x = 0.5*(vbles(ip,3,1) - vbles(ip,2,1) + vbles(ip,3,2) -vbles(ip,2,2))/dx
   phi_y = 0.5*(vbles(ip,3,1) - vbles(ip,3,2) + vbles(ip,2,1) -vbles(ip,2,2))/dx
   XX = phi_x**2
   YY = phi_y**2
   A_=A_func(XX,YY)
   LapPhi(ip) = LapPhi(ip) +0.5*A_*(phi_x*(epsilon_tilde*2.*YY*(XX-YY)/(XX+YY)**2+A_))/dx
   LapPhi(ip) = LapPhi(ip) +0.5*A_*(phi_y*(epsilon_tilde*2.*XX*(YY-XX)/(XX+YY)**2+A_))/dx
   LapPhi(ip) = LapPhi(ip)/3.
   endif !9-point stencil
  endif
 enddo
! if(Not_number(LapPhi(ip)))then
!  write(*,*)i,j,lb, "LapPhi"
!  stop
! endif
 return
end subroutine Calc_GradPhi



double precision function FreeEnergy(c,T,phi,c_dot,phi_dot,lp)
implicit none
double precision ta,tb,L_a,L_b,R,mola,molb,rhoa(3),rhob(3)
double precision,dimension (3,8)::hsa,hsb,hsrk
double precision, intent (in):: T,c,phi(3),c_dot,phi_dot(3)
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha,hb,hr,h(3),f_phi,EntropyMixing,t_a,E_const
double precision :: Tk,T_K,T_eut=456.14,Cp,RE
double precision :: MolPerVol=54185.!=(rhoa(1)+rhob(1))/(mola+molb)
! double precision :: MolPerVol=1.347708895e5 !debug
integer :: i,j,TT, N_phases=3
logical :: c_dot_term=.true. ! for extra term for PRL paper
!non dimension T to dimensioned Tk for gibbs expansion etc.

Tk=T_K(T)
R = 8.31 !J/(K mol)
! NU = 7.42e-6 !m^3 !taken from phase_8.08/init_eut_Pb_Sn_8.h
! Latent heat
L_a = 4773.925095 !J/mol 
L_b = 7029.016991
! Setting Lead diffusivity =1 then implies characteristic length
!delta_a= 1.489764340*10^(-7)*m
! and characteristic time
!tau_a = 9.809908287 10^(-10) s,

! not used: reference only
!!!!!!!!!!!!!!!!!!!!!!!
! rhoa = 10660. !kg/m^3 liquid Lead at melting point
! rhob = 6999.  !              Tin
! ta = 600.61 !Lead
! tb = 505.07 !Tin
!!!!!!!!!!!!!!!!!!!!!!!!
!  mola=.2072     !kg/mol
!  molb= 0.1187  
! rhoa(1)=10660. !liquid Lead kg/m^3
! rhob(1)=6999.  !liquid Tin
! rhoa(2)=11340. !Solid Lead
! rhoa(3)=11340. !The same (for want of a better option)
! rhob(2)=7365.  !Solid Tin (White)
! rhob(3)=7365.  !The same
! t_a=600.61 
! E_const=(R*rhoa(1)*t_a)/mola !~ 2.5 E8 J/m^3
! write(*,*)E_const
! stop
!note: but perhaps Delta T, undercooling is better than t_a
! E_const= 2.56779906e8 !J/m^3 for non-dimensionlising D F/D c

!convert Latent heat yo J/m^3 used in temperature equation.
!This approximation will be superceeded
!Not needed as can work with J/mol in temperature equation
! do i=1,3
! L_a(i)=L_a(i)*rhoa(i)/mola !J/m^3
! L_b(i)=L_b(i)*rhob(i)/molb
! enddo



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
if(lp.eq.0)then !free energy proper
 call GibbsVector(V,Tk,0)
 call MixingVector(RK,c,0)
 FreeEnergy=0.
 do i=1,3
 ha=0.
 hb=0.
 hr=0.
 do j = 1,8
  ha = ha + hsa(i,j)*V(j)
  hb = hb + hsb(i,j)*V(j)
 enddo
 do j = 1,4
   hr = hr + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
 enddo
!Convert to J/m^3 
!  ha=ha*MolPerVol
!  hb=hb*MolPerVol 
!  hr=hr*MolPerVol
 
 FreeEnergy = FreeEnergy + (c*hb+(1.-c)*ha+hr)*f_phi(phi,i,0,0)
 enddo
  FreeEnergy = FreeEnergy + R*Tk*EntropyMixing(c,0)

elseif(lp.le.N_phases)then !phase = current_var
 FreeEnergy = 0.
 do i=1,N_phases
 call GibbsVector(V,Tk,0)
 call MixingVector(RK,c,0)

 ha=0.
 hb=0.
 hr=0.
 do j = 1,8
  ha = ha + hsa(i,j)*V(j)
  hb = hb + hsb(i,j)*V(j)
 enddo
 do j = 1,4
   hr = hr + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
 enddo
!Convert to J/m^3 
 ha=ha*MolPerVol
 hb=hb*MolPerVol 
 hr=hr*MolPerVol
 FreeEnergy = FreeEnergy +(c*hb+(1.-c)*ha+hr)*f_phi(phi,i,lp,0)
 enddo
 return
elseif(lp.eq.N_phases+2)then !Solute
! the below term does not discretise well. So this term is combined
!with the diffusivity (dependent on c and phi_i and T)
!  FreeEnergy = R*Tk*EntropyMixing(c,1) !dimension (J/mol)
 FreeEnergy = 0.!debug
 call GibbsVector(V,Tk,0)
 call MixingVector(RK,c,1)

 do i= 1,N_phases
  ha=0.
  hb=0.
  hr=0.
  do j = 1,8
   ha = ha + hsa(i,j)*V(j)
   hb = hb + hsb(i,j)*V(j)
  enddo
  do j = 1,4
   hr = hr + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
  enddo
! add in other factor and non dimensionalise
! Folch Plapp
  FreeEnergy = FreeEnergy+(hb-ha+hr)*f_phi(phi,i,0,0)
 enddo
! add in other factor and non dimensionalise

 
 
 ! non-dimensionalise Free energy
  FreeEnergy = FreeEnergy/(R*T_eut)
elseif(lp.eq.N_phases+1)then !Temperature Keep in J/mol
!  T (d/dT)dF/dphi
  call GibbsVector(V,Tk,1) !first derivative of free energy
  call MixingVector(RK,c,0)
  FreeEnergy = 0.
  do i= 1,N_phases
   ha=0.
   hb=0.
   hr=0.
   do j = 1,8
    ha = ha + hsa(i,j)*V(j)
    hb = hb + hsb(i,j)*V(j)
   enddo
   do j = 1,4
    hr = hr + hsrk(i, 2*j)*RK(j)
   enddo
   do j=1,N_phases
    FreeEnergy = FreeEnergy + Tk*(c*hb+(1.-c)*ha+hr)*f_phi(phi,i,j,0)*phi_dot(j)
   enddo
  enddo
!    pure -dF/dphi term  
  call GibbsVector(V,Tk,0)
  call MixingVector(RK,c,0)
  do i= 1,N_phases
   ha=0.
   hb=0.
   hr=0.
   do j = 1,8
    ha = ha + hsa(i,j)*V(j)
    hb = hb + hsb(i,j)*V(j)
   enddo
   do j = 1,4
    hr = hr + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
   enddo
   do j=1,N_phases
    FreeEnergy = FreeEnergy - (c*hb+(1.-c)*ha+hr)*f_phi(phi,i,j,0)*phi_dot(j)
   enddo
  enddo
 if(c_dot_term)then !there is a c_dot term see PRL small for large Le
!    T(d/dT)dF/dc =-T (ds/dc) 
  call GibbsVector(V,Tk,1)
  call MixingVector(RK,c,1)
  RE = R*EntropyMixing(c,1)
  do i= 1,N_phases
   ha=0.
   hb=0.
   hr=0.
   do j = 1,8
    ha = ha + hsa(i,j)*V(j)
    hb = hb + hsb(i,j)*V(j)
   enddo
   do j = 1,4
    hr = hr +hsrk(i, 2*j)*RK(j)
   enddo
   FreeEnergy = FreeEnergy + Tk*(RE+(hb-ha+hr)*f_phi(phi,i,0,0))*c_dot
  enddo
 endif !c_dot_term
!   heat capacity = -T (d/dT) dF/dT
  call GibbsVector(V,Tk,2)
  Cp=0.
  do i= 1,N_phases
   ha=0.
   hb=0.
   do j = 1,8
    ha = ha + hsa(i,j)*V(j)
    hb = hb + hsb(i,j)*V(j)
   enddo
   Cp = Cp - (c*hb+(1.-c)*ha)*f_phi(phi,i,0,0)
  enddo
   Cp = Cp*Tk
  FreeEnergy = FreeEnergy/Cp
else
 write(*,*)"bad current_var"
 stop
endif
end function FreeEnergy

double precision function EntropyMixing(c,d)
double precision, intent (in)::c
integer, intent (in)::d

! EntropyMixing =0.
! if((1.-c).le.1e-1)then
! if(d.eq.1) EntropyMixing = 2.197224577
! return
! endif
 
 if(d.eq.1)then
  EntropyMixing = log(c/(1.-c))
 else
  EntropyMixing = c*log(c)+(1.-c)*log(1.-c)
 endif
end function EntropyMixing


subroutine GibbsVector(V,T,d)
double precision, intent (out)::V(8)
double precision, intent (in)::T
integer, intent (in)::d
if(d.eq.0)then
 V(1)=1.0
 V(2)=T
 V(3)=T*log(T)
 V(4)=T**2
 V(5)=T**3
 V(6)=1./T
 V(7)=T**7
 V(8)=1.0/T**9
else if(d.eq.1)then
 V(1)=0.
 V(2)=1.
 V(3)=1.+log(T)
 V(4)=2.*T
 V(5)=3.*T**2
 V(6)=-1./T**2
 V(7)=7.*T**6
 V(8)=-9.0/T**10
elseif(d.eq.2)then
 V(1)=0.
 V(2)=0.
 V(3)=1./T
 V(4)=2.
 V(5)=6.*T
 V(6)=2./T**3
 V(7)=42.*T**5
 V(8)=90./T**11
else
 write(*,*)"bad diff order in GibbsVector"
stop
endif
end subroutine GibbsVector
subroutine MixingVector(RK,c,d)
double precision, intent (out)::RK(5)
double precision, intent (in)::c
integer, intent (in)::d
integer i
if(d.eq.0)then
 RK(1)=c*(1.-c)
 RK(2)=RK(1)*(1.-2.*c)
 RK(3)=RK(2)*(1.-2.*c)
 RK(4)=RK(3)*(1.-2.*c)
 RK(5)=RK(4)*(1.-2.*c)
elseif(d.eq.1)then
 RK(1)=1.-2.*c
 RK(2)=1.+(-6.+6.*c)*c
 RK(3)=1.+(-10.+(24.-16.*c)*c)*c
 RK(4)=1.+(-14.+(54.+(-80.+40.*c)*c)*c)*c
 RK(5)=1.+(-18.+(96.+(-224.+(240.-96.*c)*c)*c)*c)*c
elseif(d.eq.2)then
 RK(1)=-2.
 RK(2)=12.*c-6. 
 RK(3)=-10.+(48.-48.*c)*c
 RK(4)=-14.+(108.+(-240.+160.*c)*c)*c
 RK(5)=-18.+(192.+(-672.+(960.-480.*c)*c)*c)*c
elseif(d.eq.3)then
 RK(1)=0.
 RK(2)=12.
 RK(3)=-96.*c+48.
 RK(4)=108.+(-480.+480.*c)*c
 RK(5)=192.+(-1344.+(2880.-1920.*c)*c)*c
else
 write(*,*)"bad diff order in MixingVector"
 stop
endif
end subroutine MixingVector

! double precision function g_phi(phi,d)
!   double precision, intent (in):: phi
!   integer, intent (in):: d
!   if(d.eq.1)then
!    g_phi = 30.*phi**2*(1.0-phi)**2
!   elseif(d.eq.2)then
!    g_phi = 60.0*phi*(2.0*phi-1.0)*(phi-1.0)
!   else !assume 0
!    g_phi = 6.0*phi**5 - 15.0*phi**4 + 10.0*phi**3
!   endif
! end function g_phi
double precision function g_phi(phi,d)
  double precision, intent (in):: phi
  integer, intent (in):: d
  if(d.eq.1)then
   g_phi = 6.*phi*(1.-phi)
!    g_phi = 30.*phi**4-60.*phi**3+30.*phi**2
  elseif(d.eq.2)then
   g_phi = 6.-12.*phi
!    g_phi = 120.*phi**3-180.*phi**2+60.*phi
  else !assume 0
   g_phi = 3.*phi**2-2.*phi**3
!    g_phi = 6.*phi**5-15.*phi**4+10.*phi**3
  endif
end function g_phi
double precision function gg_phi(phi,i,j,k)
  double precision, intent (in):: phi(3)
  integer, intent (in):: i,j,k
  double precision g_phi 
   
   if(k.eq.0)then
    if(j.eq.0)then
     gg_phi =  g_phi(phi(i),0)
    elseif(j.eq.i)then
     gg_phi =  g_phi(phi(i),1)
    else
     gg_phi = 0.
    endif
   elseif(k.eq.i.and.i.eq.j)then
     gg_phi =  g_phi(phi(i),2)
   else
    gg_phi = 0.
   endif
      
end function gg_phi
double precision function f_phi(phi,i,j,k)
  double precision, intent (in):: phi(3)
  double precision g_phi,s,gg_phi
  integer, intent (in):: i,j,k
  s=g_phi(phi(1),0)+g_phi(phi(2),0)+g_phi(phi(3),0)
  if(j.eq.0)then
   f_phi = g_phi(phi(i),0)
  elseif(k.eq.0)then
   f_phi = gg_phi(phi,i,j,0)/s-g_phi(phi(i),0)*gg_phi(phi,j,j,0)/s**2
  else
   f_phi = gg_phi(phi,i,j,k)/s - gg_phi(phi,i,j,0)*gg_phi(phi,k,k,0)/s**2&
        - gg_phi(phi,i,k,0)*gg_phi(phi,j,j,0)/s**2 - g_phi(phi(i),0)*gg_phi(phi,j,j,k)/s**2 &
        + 2.*gg_phi(phi,j,j,0)*g_phi(phi(i),0)*gg_phi(phi,k,k,0)/s**3
  endif
end function f_phi

logical function Not_number(x)
double precision, intent(in):: x
if(x.gt.1e20)then
Not_number=.true.
elseif(x.lt.-1e20)then
Not_number=.true.
elseif(x.le.1e20.or.x.ge.-1e20)then
Not_number=.false.
else
Not_number=.true.
endif
end function Not_number

!*********************************************************
!*   SOLVE A LINEAR SYSTEM BY TRIANGULARIZATION METHOD   *
!* ----------------------------------------------------- *
!* SAMPLE RUN:                                           *
!*                                                       *
!* Linear system AX = B                                  *
!*                                                       *
!* Matrix A:                                             *
!* 2.00    1.00   -2.00                                  *
!* 3.00   -1.00    1.00                                  *
!* 7.00    5.00   -3.00                                  *
!*                                                       *
!* Right side:                                           *
!* B(1) =    1.00                                        *
!* B(2) =    0.00                                        *
!* B(3) =    0.00                                        *
!*                                                       *
!* Solution of AX=B:                                     *
!* X(1) =    0.0625                                      *
!* X(2) =   -0.0500                                      *
!* X(3) =   -0.6875                                      *
!*                                                       *
!* ----------------------------------------------------- *
!* Reference: "Méthodes de calcul numérique - Tome 1 By  *
!*             Claude Nowakowski, PS1 1981" [BIBLI 07].  *
!*                                                       *
!*                    F90 Release By J-P Moreau, Paris.  *
!*********************************************************
subroutine Tlinear(AA,BB,X)
double precision,intent(inout) :: AA(5,5),BB(5)
double precision,intent(out) :: X(5)
double precision ::  A(5,5),B(5)
INTEGER, PARAMETER :: N = 5
do i=1,N
B(i)=BB(i)
do j=1,N
A(i,j)=AA(i,j)
enddo
enddo

! A is given by columns
! DATA A /2.d0,3.d0,7.d0,1.d0,-1.d0,5.d0,-2.d0,1.d0,-3.d0/
! DATA B /1.d0,0.d0,0.d0/

! PRINT *,' '
! PRINT *,' Linear system AX = B'
! PRINT *,' '
! PRINT *,' Matrix A:'
! DO I = 1, N
!   WRITE(*,20) (A(I, J), J=1,N)
! END DO
! PRINT *,' '
! PRINT *,' Right Side:'
! DO I=1, N
!   WRITE(*,10)  I, B(I)
! END DO

! Transform A into triangular matrix
DO K = 1, N-1
  DO I = K + 1, N
    B(I) = B(I) - A(I, K) / A(K, K) * B(K)
    DO J = K + 1, N
      A(I, J) = A(I, J) - A(I, K) / A(K, K) * A(K, J)
    END DO
  END DO
END DO

! Solve triangular system
X(N) = B(N) / A(N, N)
DO I = N - 1, 1, -1
  S = 0.d0
  DO K = I + 1, N
    S = S + A(I, K) * X(K)
  END DO
  X(I) = (B(I) - S) / A(I, I)
END DO
! X(4)=0.!BB(4)/AA(4,4)
! X(5)=BB(5)/AA(5,5)
! Print results:
! PRINT *,' '
! PRINT *,' Solution of AX=B'
! DO I=1, N
!   WRITE(*,11)  I, X(I)
! END DO
! 
! 10 format('  B(',I2,') = ',F6.2)
! 11 format('  X(',I2,') = ',F8.4)
! 20 format(5(F6.2,'  '))
! 
! print *,' '
!
END
double precision function flatten(x)
double precision, intent (in)::x
double precision:: h=0.01 !decrease to flatten more
! flatten=h*tanh(x/h)
flatten=max(-h,min(h,x))
end function flatten
double precision function f_c(x,df)
double precision, intent (in)::x
integer, intent(in)::df

f_c=0.
if(df.eq.0)then
if(x*(1-x).gt.0.)f_c=x*(1.-x)
elseif(df.eq.1)then
 if(x*(1-x).gt.0.)f_c=1.-2.*x
endif
end function f_c
double precision function f_cc(x,df)
double precision, intent (in)::x
integer, intent(in)::df
if(x*(1.-x).le.0.)then
f_cc=0.
else
 if(df.eq.0)then
  f_cc=x*(1.-x)*(log(x)-log(1.-x))
 elseif(df.eq.1)then
  f_cc=log(x)-log(1.-x)-2.*x*log(x)+2.*x*log(1.-x)+1.
 endif
endif
end function f_cc


double precision function phi_(phi,ip)
double precision, intent (in)::phi(3)
integer,intent (in):: ip
! phi_=(3.*phi(ip)-phi(1)-phi(2)-phi(3)+1.)/3.
phi_=phi(ip)/(phi(1)+phi(2)+phi(3))
end function phi_
