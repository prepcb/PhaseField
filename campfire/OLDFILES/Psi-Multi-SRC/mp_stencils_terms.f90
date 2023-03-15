!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  Utilities for phase field stencils in 2-d and 3-d
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function phi_psi(psi,i)
use solution_parameters
double precision, intent(in)::psi(2)
integer, intent(in):: i
double precision :: a,x
a=g_alpha
x=4*a**2*psi(1)*psi(2)-2*a**2*psi(1)-2*a**2*psi(2)+a**2-2*a*psi(1)+2*a*psi(2)+3
if(i.eq.1)then



   phi_psi =4.D0*a*psi(1)*psi(2)-2.D0*a*psi(1)-2.D0*psi(1)-4.D0*psi(2)+4.D0
elseif(i.eq.2)then
   phi_psi=4.D0*a*psi(1)*psi(2)-4.D0*a*psi(1)-2.D0*a*psi(2)+2.D0*a+4.D0*psi(1)+2.D0*psi(2)-2.D0
elseif(i.eq.3)then
   phi_psi=4.D0*a**2*psi(1)*psi(2)-2.D0*a**2*psi(1)-2.D0*a**2*psi(2)-8.D0*a*psi(1)*psi(2)+a**2&
        +4.D0*a*psi(1)+4.D0*a*psi(2)-2.D0*a-2.D0*psi(1)+2.D0*psi(2)+1.D0
else
   write(*,*)"i>3 in phi_psi"
   stop
endif
phi_psi=phi_psi/x
end function phi_psi


double precision function d_phi_psi(psi,j,i)
double precision, intent(in)::psi(2)
integer, intent(in):: i,j
double precision :: x,a !a=0 is linear model, a=1 is the max
a=g_alpha
x = 4*a**2*psi(1)*psi(2)-2*a**2*psi(1)-2*a**2*psi(2)+a**2-2*a*psi(1)+2*a*psi(2)+3
x=x*x
if(x.lt.1d-10)then
   write(*,*)"function not defined in d_phi_psi",psi
   stop
endif

if(i.eq.1)then
   if(j.eq.1)then
     d_phi_psi=  -2*(2*a*psi(2)-a-1)*(2*a**2*psi(2)-a**2-6*a*psi(2)+4*a-3)/x
  elseif(j.eq.2)then
     d_phi_psi= -8*a**2*psi(1)+4*a**2+24*a*psi(1)-8*a-12/x
  else
     write(*,*)"j>2",i
     stop
  endif
elseif(i.eq.2)then
   if(j.eq.1)then  
      d_phi_psi=-8*a**2*psi(2)+4*a**2+24*a*psi(2)-16*a+12/x
   elseif(j.eq.2)then
      d_phi_psi=2*(2*a*psi(1)-a+1)*(2*a**2*psi(1)-a**2-6*a*psi(1)+2*a+3)/x
  else
     write(*,*)"j>2",i
     stop
  endif
elseif(i.eq.3)then
   if(j.eq.1)then
      d_phi_psi=8*a**3*psi(2)**2-8*a**3*psi(2)-24*a**2*psi(2)**2+2*a**3+32*a**2*psi(2)-10*a**2-24*a*psi(2)+14*a-6/x
   elseif(j.eq.2)then
      d_phi_psi=-8*a**3*psi(1)**2+8*a**3*psi(1)+24*a**2*psi(1)**2-2*a**3-16*a**2*psi(1)+2*a**2-24*a*psi(1)+10*a+6/x
  else
     write(*,*)"j>2",i
     stop
  endif
else
   write(*,*)"i>3"
   stop
endif



end function d_phi_psi




double precision function GradientEnergy(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision :: psix,psiy,LapPsi,eps,d_eps,gradientEnergy0,psi(2),s,psi_(2),Psi3,psi_3(3,3)
integer :: i,ii,jj


if(lp.le.0)then
   call Calc_GradPhi(LapPsi,vbles(lp,:,:),dx)
   do i=1,2
      psix = 0.5*(vbles(i,3,2)-vbles(i,1,2))/dx
      psiy = 0.5*(vbles(i,2,3)-vbles(i,2,1))/dx
      GradientEnergy0 = GradientEnergy0 +  0.5*(psix**2 +psiy**2)
   enddo
elseif(lp.le.2)then
      GradientEnergy= -LapPsi 
      do ii=1,3
         do jj=1,3
            psi(1)=vbles(1,ii,jj)
            psi(2)=vbles(2,ii,jj)
            psi_3(ii,jj)=psi3(psi,0)
         enddo
      enddo
      call Calc_GradPhi(LapPsi,psi_3,dx)
      psi(1)=vbles(1,2,2)
      psi(2)=vbles(2,2,2)
      GradientEnergy = GradientEnergy - LapPsi*psi3(psi,lp) 
      
else
   write(*,*)"in GradientEnergy: lp > 2"
endif



end function GradientEnergy








!T_K converts non-dimensional parameter T to dimensioned T_K(T)
double precision function T_K(T)
use solution_parameters
implicit none
double precision, intent (in)::T


 
 T_K = g_T0*(1.+T)

 end function T_K
 double precision function TemperatureRHS(vbles,dx,c,T,phi,c_dot,phi_dot)
use solution_parameters
implicit none
double precision, intent (in):: dx,c,T,phi(n_phases),c_dot,phi_dot(n_phases),vbles(N_unk+1,3,3)

double precision, dimension(3,3)::stencil =reshape((/0.,1.,0.,1.,-4.,1.,0.,1.,0./),(/3,3/))
double precision, dimension(3)::stencil_x =(/-.5,.0,.5/)
double precision:: FreeEnergy,GradientEnergy,potential,x,a,xa,dFdc(3,3),d_H(3,3) !x coordinate
double precision:: Cp= 28.5!29.9758394     !J/mol/K
double precision :: D_heat_c 
double precision :: T_,c_,phi_(n_phases),MolPerVol,X_energy,Y_energy,div,Grad_DW
integer :: ii,jj,kk,ip,jp
stop
  MolPerVol = g_MolPerVol(1)*(1.-c)+g_MolPerVol(2)*c  !mol/m^3



! to help with floating point arithmetic
  X_energy=g_R*g_T0  !J/mol



  TemperatureRHS=0.
  g_C=0d0 !set heat capacity to zero. All binary interface from FE_ij contribute
  TemperatureRHS = FreeEnergy(c,T,phi(1),c_dot,phi_dot(1),N_unk-1,1,2) !J K/mol
  if(n_phases.eq.3)then
     TemperatureRHS = TemperatureRHS + FreeEnergy(c,T,phi(1),c_dot,phi_dot(1),N_unk-1,1,3) !J K/mol
     TemperatureRHS = TemperatureRHS + FreeEnergy(c,T,phi(2),c_dot,phi_dot(2),N_unk-1,2,3) !J K/mol
  endif
  Cp = (1.-c)*g_C(1)+c*g_C(2)  !J/K/mol

  Y_energy=Cp*g_T0  !J/mol  
  TemperatureRHS = TemperatureRHS/Y_energy  !no dimension

  
  do ii=1,3
     do jj=1,3
        T_=vbles(N_unk-1,ii,jj)*stencil(ii,jj)
        TemperatureRHS = TemperatureRHS + g_D_tilde*T_/(dx*dx)
     enddo
  enddo
  
! div(D_heat_c grad(dFdc))
  if(heat_c_term)then
     ip=1
     jp=2
     do ii=1,3
        do jj=1,3
           T_=vbles(N_unk-1,ii,jj)
           c_=vbles(N_unk,ii,jj)
           do kk=1,n_phases
              phi_(kk)=vbles(kk,ii,jj)
           enddo
           dFdc(ii,jj)=(FreeEnergy(c_,T_,Phi_(ip),c_dot,phi_dot(ip),N_unk,ip,jp)/X_energy)**2
           D_heat_c =  g_D(1)*phi(ip)+(1d0-phi(ip))*g_D(2) !Change
           d_H(ii,jj) = 0.5*c_*(1.-c_)*(D_heat_c/g_Dch)*(X_energy/Y_energy)   
        enddo
     enddo
     div=0d0
     div = div  + (dFdc(3,2)-dFdc(2,2))*(d_H(3,2)+d_H(2,2)) 
     div = div  + (dFdc(1,2)-dFdc(2,2))*(d_H(1,2)+d_H(2,2)) 
     div = div  + (dFdc(2,3)-dFdc(2,2))*(d_H(2,3)+d_H(2,2)) 
     div = div  + (dFdc(2,1)-dFdc(2,2))*(d_H(2,1)+d_H(2,2)) 
     
     div =0.5*div/(dx*dx)

     heat_size(3)=max(heat_size(3),abs(div*Y_energy)) !how big is this term
     Q_heat(3)=div*Y_energy
     TemperatureRHS = TemperatureRHS + div
! for 3 phases there are two more interface terms (13) and (23)
     if(n_phases.eq.3)then
        ip=1
        jp=3
        do ii=1,3
           do jj=1,3
              T_=vbles(N_unk-1,ii,jj)
              c_=vbles(N_unk,ii,jj)
              do kk=1,n_phases
                 phi_(kk)=vbles(kk,ii,jj)
              enddo
              dFdc(ii,jj)=(FreeEnergy(c_,T_,Phi_(ip),c_dot,phi_dot(ip),N_unk,ip,jp)/X_energy)**2
              D_heat_c =  g_D(1)*phi(ip)+(1d0-phi(ip))*g_D(2) !Change
              d_H(ii,jj) = 0.5*c_*(1.-c_)*(D_heat_c/g_Dch)*(X_energy/Y_energy)   
           enddo
        enddo
        div=0d0
        div = div  + (dFdc(3,2)-dFdc(2,2))*(d_H(3,2)+d_H(2,2)) 
        div = div  + (dFdc(1,2)-dFdc(2,2))*(d_H(1,2)+d_H(2,2)) 
        div = div  + (dFdc(2,3)-dFdc(2,2))*(d_H(2,3)+d_H(2,2)) 
        div = div  + (dFdc(2,1)-dFdc(2,2))*(d_H(2,1)+d_H(2,2)) 
        
        div =0.5*div/(dx*dx)
        
        heat_size(3)=max(heat_size(3),abs(div*Y_energy)) !how big is this term
        Q_heat(3)=div*Y_energy
        TemperatureRHS = TemperatureRHS + div
        
        ip=2
        jp=3
        do ii=1,3
           do jj=1,3
              T_=vbles(N_unk-1,ii,jj)
              c_=vbles(N_unk,ii,jj)
              do kk=1,n_phases
                 phi_(kk)=vbles(kk,ii,jj)
              enddo
              dFdc(ii,jj)=(FreeEnergy(c_,T_,Phi_(ip),c_dot,phi_dot(ip),N_unk,ip,jp)/X_energy)**2
              D_heat_c =  g_D(1)*phi(ip)+(1d0-phi(ip))*g_D(2) !Change
              d_H(ii,jj) = 0.5*c_*(1.-c_)*(D_heat_c/g_Dch)*(X_energy/Y_energy)   
           enddo
        enddo
        div=0d0
        div = div  + (dFdc(3,2)-dFdc(2,2))*(d_H(3,2)+d_H(2,2)) 
        div = div  + (dFdc(1,2)-dFdc(2,2))*(d_H(1,2)+d_H(2,2)) 
        div = div  + (dFdc(2,3)-dFdc(2,2))*(d_H(2,3)+d_H(2,2)) 
        div = div  + (dFdc(2,1)-dFdc(2,2))*(d_H(2,1)+d_H(2,2)) 
        
        div =0.5*div/(dx*dx)
        
        heat_size(3)=max(heat_size(3),abs(div*Y_energy)) !how big is this term
        Q_heat(3)=div*Y_energy
        TemperatureRHS = TemperatureRHS + div
     endif
  endif




end function TemperatureRHS

double precision function SoluteRHS(vbles,dx,c_c)!
! this uses volume element method. This guarantees conservation of solute. The main source of solute is the variation of free energy with phase: d (d F/dc) = ...+ (d^2 F/dc d phi) d phi + ..., which is revealed by a plot of c_dot.
use solution_parameters
implicit none
double precision, intent (in):: dx,c_c,vbles(N_unk+1,3,3)
integer ii,jj,ip,ih,jv
double precision c,T,Psi(2),FreeEnergyPhi,D
double precision, dimension(3,3)::stencil11=reshape((/0.,6.,0.,6.,-24.,6.,0.,6.,0./),(/3,3/))
double precision, dimension(3)::stencil1=(/-.5,0.,.5/)
double precision :: potential,D_c,f_c,MolPerVol,phi_psi
double precision :: c_dot=0d0,phi_dot=0d0,beta,phi(3),s1,s2,vbles_(3,3,3)
!

do ii=1,3
do jj=1,3
   psi(1)=vbles(1,ii,jj)
   psi(2)=vbles(2,ii,jj)
   vbles_(1,ii,jj)= phi_psi(psi,1)
   vbles_(2,ii,jj)= phi_psi(psi,2)
   vbles_(3,ii,jj)= phi_psi(psi,3)
enddo
enddo
SoluteRHS=0.


!x-coordinate transform values
!! Right
 ih=3
 jv=2
 include 'c_stencil.f90'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Left
 ih=1
 jv=2
 include 'c_stencil.f90'


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!up
 ih=2
 jv=3
 include 'c_stencil.f90'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !down
 ih=2
 jv=1
 include 'c_stencil.f90'

 g_c_force = SoluteRHS/(g_Dch*dx*dx*g_R*g_T0)
 SoluteRHS = SoluteRHS/(g_Dch*dx*dx*g_R*g_T0) 

end function SoluteRHS

double precision function potential(psi,c,lp)
use solution_parameters
implicit none

double precision, intent (in)::psi(2),c
integer, intent(in):: lp
integer :: i,j
double precision :: S,potential0,eps,d_eps,psi_(2),Psi3,psi_3
potential=0d0
if(lp.eq.0)then
   potential0 = 0d0
   do i=1,2
      potential0 = potential0 + psi(i)**2*(1d0-psi(i))**2
   enddo
   potential=potential0
   return
elseif(lp.le.2)then
   potential = 4*psi(lp)**3-6*psi(lp)**2+2*psi(lp)
   psi_3 = Psi3(psi,0)
   potential = potential + (4*psi_3**3-6*psi_3**2+2*psi_3)*psi3(psi,lp) 
else
   continue
endif


end function potential



double precision function PhaseRHS(vbles,dx,lp,psi,c,T)
use solution_parameters
implicit none
double precision, intent(in)::dx,c,T,Psi(2),vbles(N_unk+1,3,3)
integer, intent(in):: lp
double precision GradientEnergy,potential,FreeEnergyPhi,PsiSoluteRHS
double precision ::M_tilde 
double precision :: c_dot=0,psi_dot(2),phi(3)
double precision :: MolPerVol,beta,d_phi_psi,phi_psi

psi_dot=0d0
MolPerVol=(1.-c)*g_MolPerVol(1)+c*g_MolPerVol(2)

 !assemble and non-dimensionlise phaseRHS 

phi(1)=phi_psi(psi,1)
phi(2)=phi_psi(psi,2)
phi(3)=phi_psi(psi,3)

M_tilde = (1.-c)+c*g_M(2)/g_M(1)
phaseRHS=- M_tilde*(GradientEnergy(vbles,lp,dx)+potential(Psi,c,lp)/g_lambda**2&
     +(d_phi_psi(psi,lp,1)*FreeEnergyPhi(c,T,Phi,c_dot,psi_dot,1)&
     + d_phi_psi(psi,lp,2)*FreeEnergyPhi(c,T,Phi,c_dot,psi_dot,2)&
     + d_phi_psi(psi,lp,3)*FreeEnergyPhi(c,T,Phi,c_dot,psi_dot,3))*MolPerVol/(g_W(1)*g_lambda**2))


!phaseRHS=- M_tilde*(GradientEnergy(vbles,lp,dx)+potential(Psi,c,lp)/g_lambda**2&
!     +FreeEnergy(c,T,Phi,c_dot,psi_dot,lp)*MolPerVol/(g_W(1)*g_lambda**2))


end function PhaseRHS






!


subroutine get_RHS_Fi_2d(i,j,k,lb,dx,Fi)
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use physicaldata
  use tree
  implicit none

integer, intent(in)::i,j,k,lb
double precision, intent(in)::dx
double precision, intent(out):: Fi(N_unk)
double precision :: c,T,phi(2),c_dot,phi_dot(2)
double precision phi_rhs(2),rfactor,rf4,rf5,rf6,S
double precision TemperatureRHS,SoluteRHS,T_k,phi_
double precision vbles(N_unk+1,3,3),unk_(M_unk),MolPerVol
logical Not_number
integer ip,n,ii,jj,kk
  Fi=0.

  do ii=1,3
     do jj=1,3
        do ip=1,N_unk
           vbles(ip,ii,jj)=unk(1+(ip-1)*nunkvbles,i+ii-2,j+jj-2,1,lb)
        enddo
     enddo
  enddo

!(1.-c)+c*g_M(2)/g_M(1)
!  if(neigh(1,2,lb).eq.bc_cond_sym.and. i.eq.9)then
!     ii = 3
!     do jj=1,3
!        vbles(2,ii,jj)=delta
!     enddo
!  endif



  do ii=1,M_unk
     unk_(ii) = unk(ii,i,j,k,lb)
  enddo


  call get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_)


!
  if(thermal.and.heat_c_term)then
     unk(1+N_unk*nunkvbles,i,j,1,lb)=Q_heat(1)
     unk(2+N_unk*nunkvbles,i,j,1,lb)=Q_heat(2)
     unk(3+N_unk*nunkvbles,i,j,1,lb)=Q_heat(3)
     unk(4+N_unk*nunkvbles,i,j,1,lb)=Q_heat(4)
  endif

 
end subroutine get_RHS_Fi_2d

subroutine get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_)

  use solution_parameters

  implicit none

!
integer, intent(in)::mode_1,nunkvbles,lnblocks
double precision, intent(in)::dx,dt,dtold,unk_(M_unk)
double precision, intent(out):: Fi(N_unk)
double precision, intent(inout):: vbles(N_unk+1,3,3)
double precision :: c,T,psi(2),c_dot,psi_dot(2),psi_rhs(2)
double precision rfactor,rf4,rf5,rf6,S,GradientEnergy,potential,Projector(2,2)
double precision TemperatureRHS,SoluteRHS,PhaseRHS,T_k,psi_,FreeEnergy,MolPerVol,Sum
logical Not_number
integer ip,n

Fi=0d0

  if (mode_1.eq.1) then
     do ip=1,2
        psi_dot(ip) = (vbles(ip,2,2)-unk_(ip+1))/dt
     enddo
     n=(N_unk-1)*nunkvbles !puts pointer at beginning of solute part
     c_dot = (vbles(N_unk,2,2)-unk_(2+n))/dt
  elseif(mode_1.eq.2)then
     rfactor = dt/dtold

     rf4 = rfactor*rfactor/(rfactor+1.0)
     rf5 = (rfactor+1.0)
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
     do ip=1,2
        psi_dot(ip) = (rf4*unk_(ip+4)-rf5*unk_(ip+3)+rf6*vbles(ip,2,2))/dt
     enddo
     n=(N_unk-1)*nunkvbles
     c_dot = (rf4*unk_(5+n)-rf5*unk_(4+n)+rf6*vbles(N_unk,2,2))/dt
  else
     write(*,*)"mode_1 not = 1 or 2"
     stop
  endif
  do ip=1,2
     psi(ip)=vbles(ip,2,2)
  enddo
  T = vbles(N_unk-1,2,2)
  c = vbles(N_unk,2,2)

  do ip=1,2
     Fi(ip) =  PhaseRHS(vbles,dx,ip,psi,c,T)
  enddo

  if(thermal)then
     Fi(N_unk-1) = TemperatureRHS(vbles,dx,c,T,psi,c_dot,psi_dot)
  else
     Fi(N_unk-1)=0d0
  endif
  Fi(N_unk) = SoluteRHS(vbles,dx,c)
!
  do ip=1,N_unk
     n=(ip-1)*nunkvbles
     if(mode_1.le.1)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)-unk_(n+2)
     else if(mode_1.eq.2)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)/rf6-unk_(n+2)
     end if
  enddo





end subroutine get_RHS_Fi

subroutine Calc_GradPhi(LapPhi,dphi,dx)
use solution_parameters
implicit none
double precision, intent(in)::dx,dphi(3,3)
double precision, intent(out)::LapPhi
double precision, dimension(3,3)::stencil=reshape((/1.,4.,1.,4.,-20.,4.,1.,4.,1./),(/3,3/))
double precision, dimension(3,3):: stencil5=reshape((/0.,1.,0.,1.,-4.,1.,0.,1.,0./),(/3,3/))
integer ip,ii,jj

 LapPhi=0.
 do ii=1,3
    do jj=1,3
       LapPhi=LapPhi+stencil(ii,jj)*dphi(ii,jj)/(6.*dx*dx)
    end do
 end do
  
end subroutine Calc_GradPhi!



double precision function EntropyMixing(c,d)
double precision, intent (in)::c
integer, intent (in)::d

if(c.eq.1d0.or.c.eq.0d0)then
   write(*,*)"c=",c, "in EntropyMixing"
   stop
endif
 
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
  elseif(d.eq.2)then
   g_phi = 6.-12.*phi
!    g_phi = 120.*phi**3-180.*phi**2+60.*phi
  else !assume 0
   g_phi = 3.*phi**2-2.*phi**3
!    g_phi = 6.*phi**5-15.*phi**4+10.*phi**3
  endif
end function g_phi


logical function Not_number(x)
double precision, intent(in):: x
if(x.eq.x)then
Not_number=.false.
else
Not_number=.true.
endif
end function Not_number

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


double precision function Psi3(psi,i)
use solution_parameters
implicit none
double precision, intent (in) :: psi(2)
integer, intent(in):: i
double precision :: a
a=g_alpha
if(i.eq.0)then
   Psi3 = (2.D0*a**2*psi(1)*psi(2)-a**2*psi(1)-a**2*psi(2)+0.5D0*a**2-psi(1)-psi(2)+0.15D1)&
        /(4.D0*a**2*psi(1)*psi(2)-2.D0*a**2*psi(1)-2.D0*a**2*psi(2)+a**2+1.D0)
elseif(i.eq.1)then
   psi3 = (4.D0*a**2*psi(2)**2-4.D0*a**2*psi(2)+a**2-1.D0)&
        /(4.D0*a**2*psi(1)*psi(2)-2.D0*a**2*psi(1)-2.D0*a**2*psi(2)+a**2+1.D0)**2
elseif(i.eq.2)then
   psi3 = ((4.D0*a**2*psi(1)**2-4.D0*a**2*psi(1)+a**2-1.D0)&
        /(4.D0*a**2*psi(1)*psi(2)-2.D0*a**2*psi(1)-2.D0*a**2*psi(2)+a**2+1.D0)**2)
else
   write(*,*)"Psi12"
   stop
endif
end function Psi3


double precision function FreeEnergy(c,T,psi,c_dot,phi_dot,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,psi(2),c_dot,phi_dot(n_phases)
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(n_phases),hb(n_phases),hr(n_phases),f_phi,g_phi,EntropyMixing
double precision :: Tk,T_K,Cp,R,a
double precision :: FE,S,dS,nn(4)=(/1,2,3,1/),phi(3),ip,jp,FE12,psi3
integer :: i,j
a=g_alpha
!phi is not the bulk variable. It is psi extended to include psi(3)
phi(1)=psi(1)
phi(2)=psi(2)
phi(3)=Psi3(psi,0)


FreeEnergy=0d0



R=g_R
Tk=T_K(T)

if(lp.eq.0)then !free energy proper J/mol

  call GibbsVector(V,Tk,0)
  call MixingVector(RK,c,0)
   FreeEnergy=0.
   ha=0.
   hb=0.
   hr=0.
   do i=1,n_phases
      do j = 1,8
         ha(i) = ha(i) + hsa(i,j)*V(j)
         hb(i) = hb(i) + hsb(i,j)*V(j)
      enddo
   enddo
   do i=1,3
      do j = 1,4
         hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
      enddo
   enddo


   FreeEnergy =            ((c*hb(1)+(1.-c)*ha(1)+hr(1)) -(c*hb(2)+(1.-c)*ha(2)+hr(2)))*g_phi(phi(3),0)
   FreeEnergy = FreeEnergy+((c*hb(2)+(1.-c)*ha(2)+hr(2)) -(c*hb(3)+(1.-c)*ha(3)+hr(3)))*g_phi(phi(1),0)
   FreeEnergy = FreeEnergy+((c*hb(3)+(1.-c)*ha(3)+hr(3)) -(c*hb(1)+(1.-c)*ha(1)+hr(1)))*g_phi(phi(2),0)


   FreeEnergy = FreeEnergy + R*Tk*EntropyMixing(c,0)

   return
elseif(lp.le.2)then !phase = current_var
   if(lp.eq.1)then
      ip=2
      jp=3
   elseif(lp.eq.2)then
      ip=3
      jp=1
   else
      write(*,*)"FreeEnergy: lp=",lp
      stop
   endif
      
   if(zero_start)return
   call GibbsVector(V,Tk,0)
   call MixingVector(RK,c,0)
   FreeEnergy=0.
   ha=0.
   hb=0.
   hr=0.
   do i=1,n_phases
      do j = 1,8
         ha(i) = ha(i) + hsa(i,j)*V(j)
         hb(i) = hb(i) + hsb(i,j)*V(j)
      enddo
   enddo
   do i=1,n_phases
      do j = 1,4
         hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
      enddo
   enddo


   FreeEnergy = FreeEnergy +( (c*hb(ip)+(1.-c)*ha(ip)+hr(ip))-(c*hb(jp)+(1.-c)*ha(jp)+hr(jp)))*g_phi(phi(lp),1)

   FE12 = ( (c*hb(1)+(1.-c)*ha(1)+hr(1))-(c*hb(2)+(1.-c)*ha(2)+hr(2)))*g_phi(phi(3),1)



   FreeEnergy = FreeEnergy +Psi3(psi,lp) * FE12



  
 return
elseif(lp.eq.N_unk)then !Solute

   FreeEnergy=R*Tk*EntropyMixing(c,1)
   if(zero_start)return
   call GibbsVector(V,Tk,0)
   call MixingVector(RK,c,1)

   ha=0.
   hb=0.
   hr=0.
   do i=1,n_phases
      do j = 1,8
         ha(i) = ha(i) + hsa(i,j)*V(j)
         hb(i) = hb(i) + hsb(i,j)*V(j)
      enddo
   enddo
   do i=1,n_phases
      do j = 1,4
         hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
      enddo
   enddo


   FreeEnergy =            ((hb(1)-ha(1)+hr(1)) -(hb(2)-ha(2)+hr(2)))*g_phi(phi(3),0)
   FreeEnergy = FreeEnergy+((hb(2)-ha(2)+hr(2)) -(hb(3)-ha(3)+hr(3)))*g_phi(phi(1),0)
   FreeEnergy = FreeEnergy+((hb(3)-ha(3)+hr(3)) -(hb(1)-ha(1)+hr(1)))*g_phi(phi(2),0)


   return

elseif(lp.eq.N_unk-1)then !Temperature Keep in J/mol
   write(*,*)"in Temp FE"
   stop
   FreeEnergy = 0d0
   FE = 0d0
   if(heat_phi_term)then
!  T (d/dT)dF/dphi
      call GibbsVector(V,Tk,1) !first derivative of free energy
      call MixingVector(RK,c,0)      
      ha=0.
      hb=0.
      hr=0.
      do i=1,2
         do j = 1,8
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      do i=1,2
         do j = 1,4
            hr(i) = hr(i) + hsrk(i, 2*j)*RK(j)
         enddo!
      enddo
!      FE = FE + (c*hb(i)+(1.-c)*ha(i)+hr(i))*g_phi(phi,1)*phi_dot*Tk
!      FE = FE - (c*hb(j)+(1.-c)*ha(j)+hr(j))*g_phi(phi,1)*phi_dot*Tk


      
      
      !    pure -dF/dphi term  
      call GibbsVector(V,Tk,0)
      
      
      ha=0.
      hb=0.
      hr=0.
      do i=1,2
         do j = 1,8
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      do i=1,2
         do j = 1,4
            hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
         enddo
      enddo

!      FE = FE - (c*hb(i)+(1.-c)*ha(i)+hr(i))*g_phi(phi,1)*phi_dot
!      FE = FE + (c*hb(j)+(1.-c)*ha(j)+hr(j))*g_phi(phi,1)*phi_dot
      Q_heat(1)=FE
      heat_size(1)=max(heat_size(1),abs(FE))
!      FreeEnergy=FE
   endif !end phi_dot_term
   if(heat_c_term)then 
      !    T(d/dT)dF/dc =-T (ds/dc) ! see AT14.pdf for details

      call GibbsVector(V,Tk,1)
      call MixingVector(RK,c,1)      
      FE=0d0
      ha=0.
      hb=0.
      hr=0.
      do i=1,2
         do j = 1,8
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      do i=1,2
         do j = 1,4
            hr(i) = hr(i) + hsrk(i, 2*j)*RK(j)
         enddo
      enddo
      FE = FE + (hb(i)-ha(i)+hr(i))*g_phi(phi,0)*Tk*c_dot 
      FE = FE + (hb(j)-ha(j)+hr(j))*g_phi(1d0-phi,0)*Tk*c_dot 


      call GibbsVector(V,Tk,0)
      ha=0.
      hb=0.
      hr=0.
      do i=1,2
         do j = 1,8
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      do i=1,2
         do j = 1,4
            hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
         enddo
      enddo
      FE = FE - (hb(i)-ha(i)+hr(i))*g_phi(phi,0)*c_dot
      FE = FE - (hb(j)-ha(j)+hr(j))*g_phi(1d0-phi,0)*c_dot
      
      heat_size(2)=max(heat_size(2),abs(FE))
      Q_heat(2)=FE
      FreeEnergy=FreeEnergy + FE
   endif !c_dot_term
!   heat capacity = -T (d/dT) dF/dT
   call GibbsVector(V,Tk,2)

   g_C=0d0
   ha=0.
   hb=0.
   do j = 1,8
      ha = ha + hsa(i,j)*V(j)
      hb = hb + hsb(i,j)*V(j)
   enddo
   g_C(1) = g_C(1) - ha(i)*g_phi(phi,0) - ha(j)*(1.-g_phi(phi,0))
   g_C(2) = g_C(2) - hb(i)*g_phi(phi,0) - hb(j)*(1.-g_phi(phi,0))
   g_C=g_C*Tk
   return
else
   write(*,*)"bad current_var"
   stop
endif
end function FreeEnergy




double precision function FreeEnergyPhi(c,T,phi,c_dot,phi_dot,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi(n_phases),c_dot,phi_dot(n_phases)
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(n_phases),hb(n_phases),hr(n_phases),f_phi,g_phi,EntropyMixing
double precision :: Tk,T_K,Cp,R
double precision :: FE,S,dS
integer :: i,j

FreeEnergyPhi=0d0



R=g_R
Tk=T_K(T)

if(lp.eq.0)then !free energy proper J/mol

  call GibbsVector(V,Tk,0)
  call MixingVector(RK,c,0)
   FreeEnergyPhi=0.
   ha=0.
   hb=0.
   hr=0.
   do i=1,n_phases
      do j = 1,8
         ha(i) = ha(i) + hsa(i,j)*V(j)
         hb(i) = hb(i) + hsb(i,j)*V(j)
      enddo
   enddo
   do i=1,n_phases
      do j = 1,4
         hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
      enddo
   enddo
   S=0d0
   do i=1,n_phases
      S=S+g_phi(phi(i),0)
   enddo
   do i=1,n_phases
      FreeEnergyPhi = FreeEnergyPhi + (c*hb(i)+(1.-c)*ha(i)+hr(i))*g_phi(phi(i),0)/S
   enddo

   FreeEnergyPhi = FreeEnergyPhi + R*Tk*EntropyMixing(c,0)

   return
elseif(lp.le.n_phases)then !phase = current_var
   if(zero_start)return
   call GibbsVector(V,Tk,0)
   call MixingVector(RK,c,0)
   FreeEnergyPhi=0.
   ha=0.
   hb=0.
   hr=0.
   do i=1,n_phases
      do j = 1,8
         ha(i) = ha(i) + hsa(i,j)*V(j)
         hb(i) = hb(i) + hsb(i,j)*V(j)
      enddo
   enddo
   do i=1,n_phases
      do j = 1,4
         hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
      enddo
   enddo

   S=0d0
   do i=1,n_phases
      S=S+g_phi(phi(i),0)
   enddo
   FreeEnergyPhi = (c*hb(lp)+(1.-c)*ha(lp)+hr(lp))*g_phi(phi(lp),1)/S
   do i = 1, n_phases
      dS = -1d0/S**2*g_phi(phi(lp),1)
      FreeEnergyPhi = FreeEnergyPhi+ (c*hb(i)+(1.-c)*ha(i)+hr(i))*g_phi(phi(i),0)*dS
   enddo
  
 return
elseif(lp.eq.N_unk)then !Solute

   FreeEnergyPhi=R*Tk*EntropyMixing(c,1)
   if(zero_start)return
   call GibbsVector(V,Tk,0)
   call MixingVector(RK,c,1)

   ha=0.
   hb=0.
   hr=0.
   do i=1,n_phases
      do j = 1,8
         ha(i) = ha(i) + hsa(i,j)*V(j)
         hb(i) = hb(i) + hsb(i,j)*V(j)
      enddo
   enddo
   do i=1,n_phases
      do j = 1,4
         hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
      enddo
   enddo
   S=0d0
   do i=1,n_phases
      S=S+g_phi(phi(i),0)
   enddo
   do i=1,n_phases
      FreeEnergyPhi = FreeEnergyPhi + (hb(i) - ha(i)+hr(i))*g_phi(phi(i),0)/S
   enddo

   return

elseif(lp.eq.N_unk-1)then !Temperature Keep in J/mol
   write(*,*)"in Temp FE"
   stop
   FreeEnergyPhi = 0d0
   FE = 0d0
   if(heat_phi_term)then
!  T (d/dT)dF/dphi
      call GibbsVector(V,Tk,1) !first derivative of free energy
      call MixingVector(RK,c,0)      
      ha=0.
      hb=0.
      hr=0.
      do i=1,2
         do j = 1,8
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      do i=1,2
         do j = 1,4
            hr(i) = hr(i) + hsrk(i, 2*j)*RK(j)
         enddo!
      enddo
!      FE = FE + (c*hb(i)+(1.-c)*ha(i)+hr(i))*g_phi(phi,1)*phi_dot*Tk
!      FE = FE - (c*hb(j)+(1.-c)*ha(j)+hr(j))*g_phi(phi,1)*phi_dot*Tk


      
      
      !    pure -dF/dphi term  
      call GibbsVector(V,Tk,0)
      
      
      ha=0.
      hb=0.
      hr=0.
      do i=1,2
         do j = 1,8
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      do i=1,2
         do j = 1,4
            hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
         enddo
      enddo

!      FE = FE - (c*hb(i)+(1.-c)*ha(i)+hr(i))*g_phi(phi,1)*phi_dot
!      FE = FE + (c*hb(j)+(1.-c)*ha(j)+hr(j))*g_phi(phi,1)*phi_dot
      Q_heat(1)=FE
      heat_size(1)=max(heat_size(1),abs(FE))
!      FreeEnergyPhi=FE
   endif !end phi_dot_term
   if(heat_c_term)then 
      !    T(d/dT)dF/dc =-T (ds/dc) ! see AT14.pdf for details

      call GibbsVector(V,Tk,1)
      call MixingVector(RK,c,1)      
      FE=0d0
      ha=0.
      hb=0.
      hr=0.
      do i=1,2
         do j = 1,8
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      do i=1,2
         do j = 1,4
            hr(i) = hr(i) + hsrk(i, 2*j)*RK(j)
         enddo
      enddo
      FE = FE + (hb(i)-ha(i)+hr(i))*g_phi(phi,0)*Tk*c_dot 
      FE = FE + (hb(j)-ha(j)+hr(j))*g_phi(1d0-phi,0)*Tk*c_dot 


      call GibbsVector(V,Tk,0)
      ha=0.
      hb=0.
      hr=0.
      do i=1,2
         do j = 1,8
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      do i=1,2
         do j = 1,4
            hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
         enddo
      enddo
      FE = FE - (hb(i)-ha(i)+hr(i))*g_phi(phi,0)*c_dot
      FE = FE - (hb(j)-ha(j)+hr(j))*g_phi(1d0-phi,0)*c_dot
      
      heat_size(2)=max(heat_size(2),abs(FE))
      Q_heat(2)=FE
      FreeEnergyPhi=FreeEnergyPhi + FE
   endif !c_dot_term
!   heat capacity = -T (d/dT) dF/dT
   call GibbsVector(V,Tk,2)

   g_C=0d0
   ha=0.
   hb=0.
   do j = 1,8
      ha = ha + hsa(i,j)*V(j)
      hb = hb + hsb(i,j)*V(j)
   enddo
   g_C(1) = g_C(1) - ha(i)*g_phi(phi,0) - ha(j)*(1.-g_phi(phi,0))
   g_C(2) = g_C(2) - hb(i)*g_phi(phi,0) - hb(j)*(1.-g_phi(phi,0))
   g_C=g_C*Tk
   return
else
   write(*,*)"bad current_var"
   stop
endif
end function FreeEnergyPhi
