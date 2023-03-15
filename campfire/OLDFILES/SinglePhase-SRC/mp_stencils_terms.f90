!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  Utilities for phase field stencils in 2-d and 3-d
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  History: James Green 2008-9, Chris Goodyer, 2009-2012
!!!!  Latest version: Peter Bollada, Jun 2014 single phase wheeler type model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function Ramp(x)
implicit none
double precision, intent (in):: x
if(x.lt.1d2)then
Ramp = x*0.002
else
Ramp = 0.2
end if
end function Ramp



double precision function GradientEnergy(LapPhi,vbles,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(3,3,3),LapPhi
double precision ::phix,phiy,X,Y,u,v,A
double precision :: G11,G22,G12,M11,M12,g,GG


phix = 0.5*(vbles(1,3,2) - vbles(1,1,2))/dx 
phiy = 0.5*(vbles(1,2,3) - vbles(1,2,1))/dx 
u = phix*phix
v = phiy*phiy
g  = u + v
if(g.gt.1d-20)then

   X = u/g
   Y= v/g
   A = A_0*(1d0+epsilon_tilde*(X**2+Y**2))

   GG = 16*A_0*epsilon_tilde*X*Y*(A_0*epsilon_tilde*(X-Y)**2+A)
   G11 = Y*(GG +4*epsilon_tilde*A*A_0*(X-Y))+ A**2
   G22 = X*(GG +4*epsilon_tilde*A*A_0*(Y-X))+ A**2
   G12 = -phix*phiy*GG/g



   M11 = 0.5*(vbles(1,3,2)+vbles(1,1,2)-vbles(1,2,3)-vbles(1,2,1))/dx**2  !(phi_xx - phi_yy)/2
   M12 = 0.25*(vbles(1,3,3)+vbles(1,1,1)-vbles(1,1,3)-vbles(1,3,1))/dx**2 !phi_xy

   GradientEnergy = 0.5*LapPhi*(G11+G22)+M11*(G11-G22)+2*M12*G12


else
   GradientEnergy = (A_0*(1d0+epsilon_tilde))**2*LapPhi
endif


end function GradientEnergy

double precision function GradientEnergy_d(vbles,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(3,3,3)
double precision ::phix,phiy,X,Y,u,v,A
double precision :: G11,G22,G12,M11,M12,g,GG


phix = 0.5*(vbles(1,3,2) - vbles(1,1,2))/dx 
phiy = 0.5*(vbles(1,2,3) - vbles(1,2,1))/dx 
u = phix*phix
v = phiy*phiy
g  = u + v
if(g.gt.1d-20)then

   X = u/g
   Y= v/g
   A = A_0*(1d0+epsilon_tilde*(X**2+Y**2))

   GG = 16d0*A_0*epsilon_tilde*X*Y*(A_0*epsilon_tilde*(X-Y)**2+A)
   G11 = Y*(GG +4*epsilon_tilde*A*A_0*(X-Y))+ A**2
   G22 = X*(GG +4*epsilon_tilde*A*A_0*(Y-X))+ A**2
   G12 = -phix*phiy*GG/g

   
   GradientEnergy_d = -5d0/3d0*(G11+G22)/dx**2
else
   GradientEnergy_d = -(A_0*(1d0+epsilon_tilde))**2*1d1/3d0/dx**2
endif


end function GradientEnergy_d










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


!

double precision function SoluteRHS(vbles,dx,c_dot,phi_dot,T)
! this uses volume element method. This guarantees conservation of solute. The main source of solute is the variation of free energy with phase: d (d F/dc) = ...+ (d^2 F/dc d phi) d phi + ..., which is revealed by a plot of c_dot.
use solution_parameters

implicit none

double precision, intent (in):: dx,c_dot,vbles(3,3,3),phi_dot,T


integer ii,jj,ip,ih,jv
double precision c,Phi,psi,FreeEnergy,D,phi_,phi_tmp
double precision, dimension(3,3)::stencil11=reshape((/0.,6.,0.,6.,-24.,6.,0.,6.,0./),(/3,3/))
double precision, dimension(3)::stencil1=(/-.5,0.,.5/)
double precision, dimension(3):: D1=(/6.70e-10,6.70e-12,6.70e-12/)
double precision :: DFDc,potential,Lap_c,D_, D_c,D_cc,D_c1,D_cc1,Grad_c(2),y0,g_phi,f_c,f_cc
double precision ::D_ch =  2.83e-10! (m^2/s)    characteristic diffusivity delta2/tau=M*eps2(1)
double precision :: beta1,beta2,M_tilde,phaserhs_c,ramp


SoluteRHS=0.
D_ch =  2.83d-10




!if(df.eq.0)then
!x-coordinate transform values
!! Right
 ih=3
 jv=2
 include 'c_stencil.h'
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Left !
 ih=1
 jv=2
 include 'c_stencil.h'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!up
 ih=2
 jv=3
 include 'c_stencil.h'
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !down
 ih=2
 jv=1
 include 'c_stencil.h'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SoluteRHS = SoluteRHS/D_ch/(dx*dx) 

end function SoluteRHS



double precision function cross_term(vbles,dx,c_dot,phi_dot,phi,c_c,psi,T)
! this uses volume element method. This guarantees conservation of solute. The main source of solute is the variation of free energy with phase: d (d F/dc) = ...+ (d^2 F/dc d phi) d phi + ..., which is revealed by a plot of c_dot.
use solution_parameters
implicit none

double precision, intent (in):: dx,c_dot,c_c,psi,vbles(3,3,3),T
integer ii,jj,ip,ih,jv
double precision c,Phi,FreeEnergy,D,phi_dot,phi_,phi_tmp
double precision, dimension(3,3)::stencil11=reshape((/0.,6.,0.,6.,-24.,6.,0.,6.,0./),(/3,3/))
double precision, dimension(3)::stencil1=(/-.5,0.,.5/)
double precision, dimension(3):: D1=(/6.70e-10,6.70e-12,6.70e-12/)
double precision :: DFDc,potential,Lap_c,D_, D_c,D_cc,D_c1,D_cc1,Grad_c(2),y0,g_phi,f_c,f_cc
double precision ::D_ch =  2.83e-10! (m^2/s)    characteristic diffusivity delta2/tau=M*eps2(1)
double precision :: beta1,beta2,M_tilde,phaserhs_c



cross_term=0.
D_ch =  2.83d-10

!if(df.eq.0)then
!x-coordinate transform values
!! Right
 ih=3
 jv=2
 include 'c_stencil2.h'
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Left !
 ih=1
 jv=2
 include 'c_stencil2.h'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!up
 ih=2
 jv=3
 include 'c_stencil2.h'
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !down
 ih=2
 jv=1
 include 'c_stencil2.h'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


cross_term = cross_term/D_ch/(dx*dx) 

end function cross_term

double precision function potential(phi,c,lp)

implicit none

double precision, intent (in)::phi,c
integer, intent (in):: lp
double precision, dimension(3,3)::sigma_a,sigma_b
double precision, dimension(3)::Wa,Wb
double precision, dimension(3,3)::beta_a,beta_b,delta_a,delta_b
! double precision ::R = 8.31,T_eut=456.14,MolPerVol= 54185.3329242099
double precision :: VolPerJoule = 4.87E-009 !=1/(R*T_eut*MolPerVol)
potential=0d0
if(lp.eq.1) potential = 4*phi**3-6*phi**2+2*phi
end function potential

double precision function PsiRHS(Phi,c,T,c_dot,phi_dot,vbles,dx,unk_,mode_1,dt,dtold,nunkvbles)
use solution_parameters
implicit none
double precision, intent(in)::dx,vbles(3,3,3),phi,c,T,c_dot,phi_dot,unk_(18),dt,dtold
integer, intent(in):: mode_1,nunkvbles
double precision :: dfdphi,PhaseRHS_c,rfactor,rf4,rf5,rf6
integer ii,jj
 
 PsiRHS=0.

 dfdphi = PhaseRHS_c(Phi,c,T,c_dot,phi_dot,vbles,dx)


 if (mode_1.eq.1) then
    
    psiRHS = (dfdphi-unk_(2+2*nunkvbles))/dt
    
 elseif(mode_1.eq.2)then
    rfactor = dt/dtold
    rf4 = rfactor*rfactor/(rfactor+1.0)
    rf5 = (rfactor+1.0)
    rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
    



    psiRHS = (rf6*dfdphi+rf4*unk_(5+2*nunkvbles)-rf5*unk_(4+2*nunkvbles))/dt
    
 else
    write(*,*)"mode_1 not = 1 or 2",mode_1
    stop
 endif
 
 psiRHS = psiRHS*1d-3



end function PsiRHS






double precision function PhaseRHS(vbles,dx,Phi,c,psi,T,c_dot,phi_dot)
use solution_parameters
implicit none
double precision, intent(in)::dx,c,T,c_dot,phi_dot,Phi,psi,vbles(3,3,3)
double precision :: LapPhi
double precision GradientEnergy,potential,FreeEnergy,GE,dd_fi
double precision rf4,rf5,rf6,rfactor
double precision ::M_a= 0.948190 ,M_b= 0.5680384!m^3/(J s) mobility
double precision ::t_a = 600.61,t_b = 505.07
integer ii,jj
double precision ::W_a(3,3),beta_,cross_term
double precision ::M_tilde !s characteristic time
double precision :: delta2 = 4E-18!decrease to sharpen. This i.c. must increase by similar factor and 
integer ip,n

!!!!! Calculate characteristic time and space (squared)

! tau_a=1./(M_a*W_a)
 PhaseRHS=0.

 w_a(1,1)=7e7 ! see notes = 0.07 / delta = 0.07/1e-9 

 call Calc_GradPhi(LapPhi,vbles,dx)
!assemble and non-dimensionlise phaseRHS 

 M_tilde = 1.0-0.159*c ! see notes Anisotropic2d3d.tex !!

!  M_tilde = M_tilde*(1+beta*alpha*lambda)

!
  PhaseRHS=M_tilde*(GradientEnergy(LapPhi,vbles,dx)-potential(Phi,c,1)/lambda**2&
       -FreeEnergy(c,T,Phi,c_dot,phi_dot,1)/lambda/W_a(1,1))
!
  if(beta.gt.1d-6)then
     phaserhs = phaserhs - cross_term(vbles,dx,c_dot,phi_dot,phi,c,psi,T)
  endif




end function PhaseRHS




double precision function PhaseRHS_c(Phi,c,T,c_dot,phi_dot,vbles,dx)

use solution_parameters
implicit none
double precision, intent(in)::c,T,c_dot,phi_dot,Phi,vbles(3,3,3),dx
double precision potential,FreeEnergy,GE,dd_fi,gradientEnergy,LapPhi
double precision ::W_a(3,3)

 w_a(1,1)=7e7 ! see notes = 0.07 / delta = 0.07/1e-9 

!
!assemble and non-dimensionlise phaseRHS 

!phaseRHS_c is (D F/D phi)/W_a(1,1)

 call Calc_GradPhi(LapPhi,vbles,dx)
 
 PhaseRHS_c = -lambda*GradientEnergy(LapPhi,vbles,dx)+potential(Phi,c,1)/lambda +FreeEnergy(c,T,Phi,c_dot,phi_dot,1)/W_a(1,1)


end function PhaseRHS_c


   

!
subroutine Project_Phase(vbles,dx,phi_rhs,Phi,c,psi,T,c_dot,phi_dot)
implicit none

double precision, intent(in)::dx,Phi,c,psi,T,c_dot,phi_dot,vbles(3,3,3)
double precision, intent(out) :: phi_rhs
double precision PhaseRHS



   phi_rhs =  PhaseRHS(vbles,dx,Phi,c,psi,T,c_dot,phi_dot) 

 
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
double precision, intent(out):: Fi(3)
double precision vbles(3,3,3),unk_(18),phi,c,psi,dd_phi,s
double precision rf4,rf5,rf6,rfactor,x
integer ip,ii,jj,iflag
  Fi=0.



  do ii=1,3
     do jj=1,3
        do ip=1,total_vars
           vbles(ip,ii,jj)=unk(1+(ip-1)*nunkvbles,i+ii-2,j+jj-2,k,lb)
        enddo
     enddo
  enddo


  do ii=1,nvar
     unk_(ii) = unk(ii,i,j,k,lb)
  enddo
  phi = vbles(1,2,2)
  c=vbles(2,2,2)
  psi=vbles(3,2,2)
  x = bnd_box(1,1,lb)+(i-1.5)*dx
  call get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_,phi,c,psi,time,x)

  
  
  do ip=1,total_vars
     unk(1+(ip-1)*nunkvbles,i,j,k,lb)=vbles(ip,2,2)
  enddo
  
  

 
end subroutine get_RHS_Fi_2d

subroutine get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_,phi,c,psi,time,x)
use solution_parameters


  implicit none

integer, intent(in)::mode_1,nunkvbles,lnblocks
double precision, intent(in)::dx,dt,dtold,unk_(18),phi,c,psi,time,x
double precision, intent(out):: Fi(3)
double precision :: T,c_dot,phi_dot,psi_dot
double precision phi_rhs,dphi_rhs,rfactor,rf4,rf5,rf6,S
double precision SoluteRHS,PsiRHS,T_k,phi_,vbles(3,3,3),ff,gg,u1,u2,ramp
logical Not_number
integer ip,n,ii,jj,kk ,i,j,lb


!these three lines are to help tapenade differentiation
vbles(1,2,2)=phi
vbles(2,2,2)=c
vbles(3,2,2)=psi


if (mode_1.eq.1) then
   
   
   
   phi_dot = (phi-unk_(2))/dt
   
   c_dot = (c-unk_(2+nunkvbles))/dt
   


elseif(mode_1.eq.2)then
   rfactor = dt/dtold
   rf4 = rfactor*rfactor/(rfactor+1.0)
   rf5 = (rfactor+1.0)
   rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
   

   phi_dot = (rf4*unk_(5)-rf5*unk_(4)+rf6*phi)/dt
   c_dot = (rf4*unk_(5+nunkvbles)-rf5*unk_(4+nunkvbles)+rf6*c)/dt


else
   write(*,*)"mode_1 not = 1 or 2",mode_1
   stop
endif



vbles(1,2,2)=max(0d0,phi)

vbles(1,2,2)=min(1d0,phi)
!T = min(delta,delta - (time - 3d2)*3d-2/3d3)
T = delta
!T = ramp(x+4d1-0.04*time)



call Project_Phase(vbles,dx,phi_rhs,Phi,c,psi,T,c_dot,phi_dot)


Fi(1)=phi_rhs


Fi(2) = SoluteRHS(vbles,dx,c_dot,phi_dot,T)

Fi(3)=PsiRHS(Phi,c,T,c_dot,phi_dot,vbles,dx,unk_,mode_1,dt,dtold,nunkvbles)
end subroutine get_RHS_Fi



! for now only calculates LapPhi
subroutine Calc_GradPhi(LapPhi,vbles,dx)

implicit none

double precision, intent(in)::dx,vbles(3,3,3)
double precision, intent(out)::LapPhi
double precision, dimension(3,3)::stencil=reshape((/1.,4.,1.,4.,-20.,4.,1.,4.,1./),(/3,3/))
integer ii,jj
 
LapPhi=0.0
do ii=1,3
   do jj=1,3
      LapPhi=LapPhi+stencil(ii,jj)*vbles(1,ii,jj)
   enddo
enddo
LapPhi = LapPhi/6./dx**2

end subroutine Calc_GradPhi



double precision function FreeEnergy(c,T,phi,c_dot,phi_dot,lp)
implicit none
double precision ta,tb,L_a,L_b,R,mola,molb,rhoa(3),rhob(3)
double precision,dimension (3,8)::hsa,hsb,hsrk
double precision, intent (in):: T,c,phi,c_dot,phi_dot
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha,hb,hr,h(3),f_phi,EntropyMixing,t_a,E_const
double precision :: Tk,T_K,T_eut=456.14,Cp,RE
double precision :: MolPerVol=54185.!=(rhoa(1)+rhob(1))/(mola+molb)
! double precision :: MolPerVol=1.347708895e5 !debug
integer :: i,j,TT

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
 do i=1,2
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


 
 FreeEnergy = FreeEnergy + (c*hb+(1.-c)*ha+hr)*f_phi(phi,i,0)
 enddo
  FreeEnergy = FreeEnergy + R*Tk*EntropyMixing(c,0)

elseif(lp.eq.1)then !phase = current_var
 FreeEnergy = 0.

 call MixingVector(RK,c,0)
 i=1
 call GibbsVector(V,Tk,0)


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

 FreeEnergy = FreeEnergy +(c*hb+(1.-c)*ha+hr)*f_phi(phi,i,1)*molpervol

i=2
 call GibbsVector(V,Tk,0)


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

 FreeEnergy = FreeEnergy -(c*hb+(1.-c)*ha+hr)*f_phi(phi,i,1)*molpervol




 return
elseif(lp.eq.2)then !Solute
! the below term does not discretise well. So this term is combined
!with the diffusivity (dependent on c and phi_i and T)
!  FreeEnergy = R*Tk*EntropyMixing(c,1) !dimension (J/mol)
 FreeEnergy = 0.!debug
 call GibbsVector(V,Tk,0)
 call MixingVector(RK,c,1)

 do i= 1,2
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
  FreeEnergy = FreeEnergy+(hb-ha+hr)*f_phi(phi,i,0)
 enddo
! add in other factor and non dimensionalise

 
 
 ! non-dimensionalise Free energy
  FreeEnergy = FreeEnergy/(R*T_eut)
else
 write(*,*)"bad current_variable",lp
 stop
endif
end function FreeEnergy

double precision function EntropyMixing(c,d)
double precision, intent (in)::c
integer, intent (in)::d

 
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


double precision function g_phi(phi,d)
  double precision, intent (in):: phi
  integer, intent (in):: d
  if(d.eq.1)then
   g_phi = 6.*phi*(1.-phi)
!    g_phi = 30.*phi**4-60.*phi**3+30.*phi**2
  else !assume 0
   g_phi = 3.*phi**2-2.*phi**3
!    g_phi = 6.*phi**5-15.*phi**4+10.*phi**3
  endif
end function g_phi

double precision function f_phi(phi,i,j)
  double precision, intent (in):: phi
  double precision g_phi,s,gg_phi
  integer, intent (in):: i,j
  if(j.eq.0)then
     if(i.eq.1)then
        f_phi = g_phi(phi,0) !solute or free energy
     else
        f_phi = g_phi(1-phi,0) !solute or free energy
     endif
  else
     if(i.eq.1)then
        f_phi = g_phi(phi,1) ! phase so: d F/(d phi)
     else
        f_phi = g_phi(1-phi,1) ! phase so: d F/(d phi)
     endif
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
double precision,intent(inout) :: AA(3,3),BB(3)
double precision,intent(out) :: X(3)
double precision ::  A(3,3),B(3)
INTEGER, PARAMETER :: N = 3
do i=1,N
B(i)=BB(i)
do j=1,N
A(i,j)=AA(i,j)
enddo
enddo

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

!
END

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



