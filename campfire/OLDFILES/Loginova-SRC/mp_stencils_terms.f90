!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  Utilities for phase field stencils in 2-d and 3-d
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




double precision function GradientEnergy(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk,3,3)
integer, intent (in):: lp
double precision :: LapPhi(1)

!isotropic   
call Calc_GradPhi(LapPhi,vbles,dx)
                                                                                             
GradientEnergy = - LapPhi(lp)
  
  
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
double precision, intent (in):: dx,c,T,phi(N_phases),c_dot,phi_dot(N_phases),vbles(N_unk,3,3)

double precision, dimension(3,3)::stencil =reshape((/0.,1.,0.,1.,-4.,1.,0.,1.,0./),(/3,3/))
double precision, dimension(3)::stencil_x =(/-.5,.0,.5/)
double precision:: FreeEnergy,GradientEnergy,potential,x,a,xa,dFdc(3,3),d_H(3,3) !x coordinate
double precision:: Cp= 28.5!29.9758394     !J/mol/K
double precision :: D_heat_c 
double precision :: T_,c_,phi_(N_phases),MolPerVol,X_energy,Y_energy,div,Grad_DW
integer :: ii,jj,kk

  MolPerVol = g_MolPerVol(1)*0.8+g_MolPerVol(2)*0.2
  Cp = (1.-c)*g_C(1)+c*g_C(2)


! to help with floating point arithmetic
  X_energy=g_R*g_T0*MolPerVol
  Y_energy=Cp*g_T0


  TemperatureRHS=0.




  
  TemperatureRHS = FreeEnergy(c,T,phi,c_dot,phi_dot,N_phases+1) !J K/m^3

  TemperatureRHS = TemperatureRHS/(Cp*g_T0)

  
  do ii=1,3
     do jj=1,3
        T_=vbles(N_phases+1,ii,jj)*stencil(ii,jj)
        TemperatureRHS = TemperatureRHS + g_D_tilde*T_/(dx*dx)
     enddo
  enddo
  
! div(D_heat_c grad(dFdc))
  if(heat_c_term)then
     do ii=1,3
        do jj=1,3
           T_=vbles(2,ii,jj)
           c_=vbles(3,ii,jj)
           phi_(1)=vbles(1,ii,jj)
           dFdc(ii,jj)=(FreeEnergy(c_,T_,Phi_,c_dot,phi_dot,3)/X_energy)**2
           D_heat_c =  g_D(1)*phi(1)+(1d0-phi(1))*g_D(2)
           d_H(ii,jj) = 0.5*c_*(1.-c_)*(D_heat_c/g_Dch)*(X_energy/Y_energy)   
        enddo
     enddo
     div=0d0
     div = div  + (dFdc(3,2)-dFdc(2,2))*(d_H(3,2)+d_H(2,2)) 
     div = div  + (dFdc(1,2)-dFdc(2,2))*(d_H(1,2)+d_H(2,2)) 
     div = div  + (dFdc(2,3)-dFdc(2,2))*(d_H(2,3)+d_H(2,2)) 
     div = div  + (dFdc(2,1)-dFdc(2,2))*(d_H(2,1)+d_H(2,2)) 
     
     div =0.5*div/(dx*dx)

     heat_size(3)=max(heat_size(3),abs(div*Cp*g_T0)) !how big is this term

     TemperatureRHS = TemperatureRHS + div
  endif
! grad double well term in phase equation
  if(grad_DW_term)then
     Grad_DW = (GradientEnergy(vbles,1,dx)+potential(Phi,c,1))*g_W(1)
     heat_size(4)=max(heat_size(4),abs(Grad_DW))
    TemperatureRHS = TemperatureRHS + Grad_DW/(Cp*g_T0)
  endif


!
 end function TemperatureRHS

double precision function SoluteRHS(vbles,dx,c_c)
! this uses volume element method. This guarantees conservation of solute. The main source of solute is the variation of free energy with phase: d (d F/dc) = ...+ (d^2 F/dc d phi) d phi + ..., which is revealed by a plot of c_dot.
use solution_parameters
implicit none
double precision, intent (in):: dx,c_c,vbles(N_unk,3,3)
integer ii,jj,ip,ih,jv
double precision c,T,Phi(N_phases),FreeEnergy,D
double precision, dimension(3,3)::stencil11=reshape((/0.,6.,0.,6.,-24.,6.,0.,6.,0./),(/3,3/))
double precision, dimension(3)::stencil1=(/-.5,0.,.5/)
double precision :: potential,D_c,f_c
double precision::  MolPerVol
double precision :: c_dot=0d0,phi_dot(1)=0d0



SoluteRHS=0.


!!
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


 MolPerVol = g_MolPerVol(1)*0.8+g_MolPerVol(2)*0.2



 SoluteRHS = SoluteRHS/(g_Dch*dx*dx)/(g_R*g_T0*MolPerVol) 

return

end function SoluteRHS



double precision function potential(phi,c,lp)
use solution_parameters
implicit none

double precision, intent (in)::phi(N_phases),c
integer, intent(in):: lp
integer :: ip,jp

potential=0d0
if(lp.eq.1)then
   potential = ((1.-c)+c*g_W(2)/g_W(1))* 2d0*phi(1)*(1d0-phi(1))*(1d0-2d0*phi(1))
elseif(lp.eq.3)then
   potential = (g_W(2)/g_W(1)-1d0)* phi(1)**2*(1d0-phi(1))**2
else
   stop
endif

return


end function potential




double precision function PhaseRHS(vbles,dx,lp,Phi,c,T)
use solution_parameters
implicit none
double precision, intent(in)::dx,c,T,Phi(N_phases),vbles(N_unk,3,3)
integer, intent(in):: lp
double precision :: LapPhi(N_phases)
double precision GradientEnergy,potential,FreeEnergy
double precision ::M_tilde 
double precision :: c_dot=0,phi_dot=0




 !assemble and non-dimensionlise phaseRHS 

 M_tilde = (1.-c)+c*g_M(2)/g_M(1)
 phaseRHS=- M_tilde*(GradientEnergy(vbles,lp,dx)+potential(Phi,c,lp)&
      +FreeEnergy(c,T,Phi,c_dot,phi_dot,1)/g_W(1))


end function PhaseRHS









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
double precision :: c,T,phi(N_phases),c_dot,phi_dot(N_phases)
double precision phi_rhs(N_phases),rfactor,rf4,rf5,rf6,S
double precision TemperatureRHS,SoluteRHS,T_k,phi_
double precision vbles(N_unk,3,3),unk_(M_unk)
logical Not_number
integer ip,n,ii,jj,kk
  Fi=0.

  do ii=1,3
     do jj=1,3
        do ip=1,N_unk
           vbles(ip,ii,jj)=unk(1+(ip-1)*nunkvbles,i+ii-2,j+jj-2,k,lb)
        enddo
     enddo
  enddo




  if(neigh(1,2,lb).eq.bc_cond_sym.and. i.eq.9)then
     ii = 3
     do jj=1,3
        vbles(2,ii,jj)=delta
     enddo
  endif


  do ii=1,M_unk
     unk_(ii) = unk(ii,i,j,k,lb)
  enddo

  
  call get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_)




 
end subroutine get_RHS_Fi_2d

subroutine get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_)

  use solution_parameters

  implicit none


integer, intent(in)::mode_1,nunkvbles,lnblocks
double precision, intent(in)::dx,dt,dtold,unk_(M_unk)
double precision, intent(out):: Fi(N_unk)
double precision, intent(in):: vbles(N_unk,3,3)
double precision :: c,T,phi(N_phases),c_dot,phi_dot(N_phases)
double precision rfactor,rf4,rf5,rf6,S
double precision TemperatureRHS,SoluteRHS,PhaseRHS,T_k,phi_
logical Not_number
integer ip,n

Fi=0d0

  if (mode_1.eq.1) then
     phi_dot(1) = (vbles(1,2,2)-unk_(2))/dt
     n=2*nunkvbles
     c_dot = (vbles(3,2,2)-unk_(2+n))/dt
  elseif(mode_1.eq.2)then
     rfactor = dt/dtold

!     rf1 = (rfactor+1.0)*dt/(2.0*rfactor+1.0)
!     rf2 = (rfactor+1.0)*(rfactor+1.0)/(2.0*rfactor+1.0)
!     rf3 = rfactor*rfactor/(2.0*rfactor+1.0)

     rf4 = rfactor*rfactor/(rfactor+1.0)
     rf5 = (rfactor+1.0)
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
     
     phi_dot(1) = (rf4*unk_(5)-rf5*unk_(4)+rf6*vbles(1,2,2))/dt

     n=2*nunkvbles
     c_dot = (rf4*unk_(5+n)-rf5*unk_(4+n)+rf6*vbles(3,2,2))/dt
  else
     write(*,*)"mode_1 not = 1 or 2"
     stop
  endif


  phi(1)=vbles(1,2,2)
  T = vbles(2,2,2)!
  c = vbles(3,2,2)

  


  Fi(1) = PhaseRHS(vbles,dx,1,Phi,c,T)
  Fi(2) = TemperatureRHS(vbles,dx,c,T,phi,c_dot,phi_dot)
  Fi(3) = SoluteRHS(vbles,dx,c)

  do ip=1,N_unk
     n=(ip-1)*nunkvbles
     if(mode_1.le.1)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)-unk_(n+2)
     else if(mode_1.eq.2)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)/rf6-unk_(n+2)
     end if
  enddo


end subroutine get_RHS_Fi

subroutine Calc_GradPhi(LapPhi,vbles,dx)
use solution_parameters
implicit none
double precision, intent(in)::dx,vbles(N_unk,3,3)
double precision, intent(out)::LapPhi(N_phases)
double precision, dimension(3,3)::stencil=reshape((/1.,4.,1.,4.,-20.,4.,1.,4.,1./),(/3,3/))
double precision, dimension(3,3):: stencil5=reshape((/0.,1.,0.,1.,-4.,1.,0.,1.,0./),(/3,3/))
integer ip,ii,jj

 
 LapPhi=0.
 do ip=1,N_phases
    do ii=1,3
       do jj=1,3
          LapPhi(ip)=LapPhi(ip)+stencil(ii,jj)*vbles(ip,ii,jj)/(6.*dx*dx)
       end do
    end do
 end do
 return
 
end subroutine Calc_GradPhi


double precision function FreeEnergy(c,T,phi,c_dot,phi_dot,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi(N_phases),c_dot,phi_dot(N_phases)
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(2),hb(2),hr,f_phi,g_phi,EntropyMixing
double precision :: Tk,T_K,Cp,R
double precision :: MolPerVol,FE
integer :: i,j


R=g_R
MolPerVol= g_MolPerVol(1)*0.8+g_MolPerVol(2)*0.2
Tk=T_K(T)






if(lp.eq.0)then !free energy proper J/m^3

  call GibbsVector(V,Tk,0)
   FreeEnergy=0.
   ha=0.
   hb=0.
   
   do i=1,2
      do j = 1,3
         ha(i) = ha(i) + hsa(i,j)*V(j)
         hb(i) = hb(i) + hsb(i,j)*V(j)
      enddo
   enddo
   FreeEnergy = FreeEnergy + (c*hb(1)+(1.-c)*ha(1))*g_phi(phi(1),0)
   FreeEnergy = FreeEnergy + (c*hb(2)+(1.-c)*ha(2))*g_phi(1d0-phi(1),0)
   FreeEnergy = FreeEnergy + R*Tk*(c*g_MolPerVol(2)+(1-c)*g_MolPerVol(1))*EntropyMixing(c,0)
   return
elseif(lp.eq.1)then !phase = current_var
   call GibbsVector(V,Tk,0)
   FreeEnergy=0.
   ha=0.
   hb=0.
   
   do i=1,2
      do j = 1,3
         ha(i) = ha(i) + hsa(i,j)*V(j)
         hb(i) = hb(i) + hsb(i,j)*V(j)
      enddo
   enddo
   FreeEnergy = FreeEnergy + (c*hb(1)+(1.-c)*ha(1))*g_phi(phi(1),1)
   FreeEnergy = FreeEnergy - (c*hb(2)+(1.-c)*ha(2))*g_phi(phi(1),1)

 return
elseif(lp.eq.3)then !Solute



   FreeEnergy =0d0

   call GibbsVector(V,Tk,0)

   ha=0.
   hb=0.
   
   do i=1,2
      do j = 1,3
         ha(i) = ha(i) + hsa(i,j)*V(j)
         hb(i) = hb(i) + hsb(i,j)*V(j)
      enddo
   enddo


   FreeEnergy = FreeEnergy + (hb(1) - ha(1))*g_phi(phi(1),0)
   FreeEnergy = FreeEnergy + (hb(2) - ha(2))*(1.-g_phi(phi(1),0))
   FreeEnergy = FreeEnergy + R*Tk*((g_MolPerVol(2) - g_MolPerVol(1))*EntropyMixing(c,0)&
        +(c*g_MolPerVol(2)+(1d0-c)*g_MolPerVol(1))*EntropyMixing(c,1))

  return
elseif(lp.eq.2)then !Temperature Keep in J/mol
   FreeEnergy = 0d0
   FE = 0d0
   if(heat_phi_term)then
!  T (d/dT)dF/dphi
      call GibbsVector(V,Tk,1) !first derivative of free energy
      ha=0.
      hb=0.
      
      do i=1,2
         do j = 1,3 !,2  since terms linear in T do not contribute
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      FE = FE + (c*hb(1)+(1.-c)*ha(1))*g_phi(phi(1),1)*phi_dot(1)*Tk
      FE = FE - (c*hb(2)+(1.-c)*ha(2))*g_phi(phi(1),1)*phi_dot(1)*Tk


      
      
      !    pure -dF/dphi term  
      call GibbsVector(V,Tk,0)
      
      
      ha=0.
      hb=0.
      
      do i=1,2
         do j = 1,3!,2  since terms linear in T do not contribute
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      FE = FE - (c*hb(1)+(1.-c)*ha(1))*g_phi(phi(1),1)*phi_dot(1)
      FE = FE + (c*hb(2)+(1.-c)*ha(2))*g_phi(phi(1),1)*phi_dot(1)

      heat_size(1)=max(heat_size(1),abs(FE))
      FreeEnergy=FE
   endif !end phi_dot_term
   if(heat_c_term)then 
      !    T(d/dT)dF/dc =-T (ds/dc) ! see AT14.pdf for details

      call GibbsVector(V,Tk,1)
      
      FE=0d0
      ha=0.
      hb=0.
      
      do i=1,2
         do j = 1,3!,2  since terms linear in T do not contribute
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      FE = FE + (hb(1)+ha(1))*g_phi(phi(1),0)*Tk*c_dot 
      FE = FE + (hb(2)+ha(2))*g_phi(1d0-phi(1),0)*Tk*c_dot 


      call GibbsVector(V,Tk,0)
      ha=0.
      hb=0.
      
      do i=1,2
         do j = 1,3!,2  since terms linear in T do not contribute
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      FE = FE - (hb(1)+ha(1))*g_phi(phi(1),0)*c_dot
      FE = FE - (hb(2)+ha(2))*g_phi(1d0-phi(1),0)*c_dot
      
      heat_size(2)=max(heat_size(2),abs(FE))
      FreeEnergy=FreeEnergy + FE
   endif !c_dot_term

   return
else
   write(*,*)"bad current_var"
   stop
endif
end function FreeEnergy


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
