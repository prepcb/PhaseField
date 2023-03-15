!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  Utilities for phase field stencils in 2-d and 3-d
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





double precision function GradientEnergy(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision ::phix,phiy
double precision :: M11,M12,M22,g,q
double precision :: g0,h0,gx,hx,gxx,hxx,gy,hy,gyy,hyy,hxy,phixx,phiyy,phixy

q=.2d0/g_lambda

phix = 0.5*(vbles(1,3,2) - vbles(1,1,2))/dx 
phiy = 0.5*(vbles(1,2,3) - vbles(1,2,1))/dx 
phixx = (vbles(1,3,2) - 2d0*vbles(1,2,2) + vbles(1,1,2))/dx**2
phiyy = (vbles(1,2,3) - 2d0*vbles(1,2,2) + vbles(1,2,1))/dx**2
phixy = 0.25*(vbles(1,3,3)+vbles(1,1,1)-vbles(1,3,1)-vbles(1,1,3))/dx**2



g0  = phix*phix +phiy*phiy +q**2
gx = 2d0*phix
gy = 2d0*phiy
gxx=2d0
gyy=2d0


 h0 = (phix**3-3*phix*phiy**2)**2 + q**6
 hx = 2*(phix**3-3*phix*phiy**2)*(3*phix**2-3*phiy**2)
 hxx = 2*(3*phix**2-3*phiy**2)**2+12*(phix**3-3*phix*phiy**2)*phix
 hy = -12*(phix**3-3*phix*phiy**2)*phix*phiy
 hyy = 72*phix**2*phiy**2-12*(phix**3-3*phix*phiy**2)*phix
 hxy = -12*phix*phiy*(3*phix**2-3*phiy**2)-12*(phix**3-3*phix*phiy**2)*phiy







m11 = -(-2*h0*g0**2*hxx+g0**2*hx**2+2*h0*g0*hx*gx+2*h0**2*g0*gxx-3*h0**2*gx**2)/h0/g0**3/sqrt(h0/g0)*epsilon/3+gxx/2

     
m22 = -(-2*h0*g0**2*hyy+g0**2*hy**2+2*h0*g0*gy*hy+2*h0**2*g0*gyy-3*h0**2*gy**2)/h0/g0**3/sqrt(h0/g0)*epsilon/3+gyy/2


m12 = -(-2*h0*g0**2*hxy+g0**2*hx*hy+h0*g0*gx*hy+h0*g0*gy*hx-3*h0**2*gx*gy)/h0/g0**3/sqrt(h0/g0)*epsilon/3





GradientEnergy = phixx*m11+2d0*phixy*m12+phiyy*m22


GradientEnergy = - GradientEnergy




end function GradientEnergy


double precision function GradientEnergy4(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision ::phix,phiy,X,Y,u,v,A
double precision :: G11,G22,G12,M11,M12,g,GG
double precision :: LapPhi(1)

call Calc_GradPhi(LapPhi(1),vbles,dx)

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

   GradientEnergy4 = 0.5*LapPhi(1)*(G11+G22)+M11*(G11-G22)+2*M12*G12


else
   GradientEnergy4 = (A_0*(1d0+epsilon_tilde))**2*LapPhi(1)
endif

GradientEnergy4 = - GradientEnergy4

!
end function GradientEnergy4








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
double precision, intent (in):: dx,c,T,phi(N_phases),c_dot,phi_dot(N_phases),vbles(N_unk+1,3,3)

double precision, dimension(3,3)::stencil =reshape((/0.,1.,0.,1.,-4.,1.,0.,1.,0./),(/3,3/))
double precision, dimension(3)::stencil_x =(/-.5,.0,.5/)
double precision:: FreeEnergy,GradientEnergy,potential,x,a,xa,dFdc(3,3),d_H(3,3) !x coordinate
double precision:: Cp= 28.5!29.9758394     !J/mol/K
double precision :: D_heat_c 
double precision :: T_,c_,phi_(N_phases),MolPerVol,X_energy,Y_energy,div,Grad_DW
integer :: ii,jj,kk

  MolPerVol = g_MolPerVol(1)*(1.-c)+g_MolPerVol(2)*c  !mol/m^3



! to help with floating point arithmetic
  X_energy=g_R*g_T0  !J/mol



  TemperatureRHS=0.

  TemperatureRHS = FreeEnergy(c,T,phi,c_dot,phi_dot,N_phases+1) !J K/mol
  Cp = g_C  !J/K/mol

  Y_energy=Cp*g_T0  !J/mol  
  TemperatureRHS = TemperatureRHS/Y_energy  !no dimension

  
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

     heat_size(3)=max(heat_size(3),abs(div*Y_energy)) !how big is this term
     Q_heat(3)=div*Y_energy
     TemperatureRHS = TemperatureRHS + div
  endif
! grad double well term in phase equation
  if(grad_DW_term)then
     Grad_DW = (g_lambda*GradientEnergy(vbles,1,dx)+potential(Phi,c,1)/g_lambda)*g_W(1)/MolPerVol*phi_dot(1) !J/mol
     heat_size(4)=max(heat_size(4),abs(Grad_DW))
     Q_heat(4)=-Grad_DW
    TemperatureRHS = TemperatureRHS - Grad_DW/(Cp*g_T0)   !no dimension  
  endif


end function TemperatureRHS

double precision function SoluteRHS(vbles,dx,c_c)
! this uses volume element method. This guarantees conservation of solute. The main source of solute is the variation of free energy with phase: d (d F/dc) = ...+ (d^2 F/dc d phi) d phi + ..., which is revealed by a plot of c_dot.
use solution_parameters
implicit none
double precision, intent (in):: dx,c_c,vbles(N_unk+1,3,3)
integer ii,jj,ip,ih,jv
double precision c,T,Phi(N_phases),FreeEnergy,D
double precision, dimension(3,3)::stencil11=reshape((/0.,6.,0.,6.,-24.,6.,0.,6.,0./),(/3,3/))
double precision, dimension(3)::stencil1=(/-.5,0.,.5/)
double precision :: potential,D_c,f_c,MolPerVol
double precision :: c_dot=0d0,phi_dot(1)=0d0,beta,phi_,PhiSoluteRHS=0d0



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

 PhiSoluteRHS=0d0
 if(cross_term)then
    !!
    !x-coordinate transform values
    !! Right
    ih=3
    jv=2
    
    include 'phi_stencil.f90'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!Left
    ih=1
    jv=2
    
    include 'phi_stencil.f90'
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !up
    ih=2
    jv=3
    
    include 'phi_stencil.f90'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !down
    ih=2
    jv=1
    
    include 'phi_stencil.f90'
 endif
 
 g_anti_trap = PhiSoluteRHS/(g_Dch*dx*dx*g_R*g_T0)
 g_c_force = SoluteRHS/(g_Dch*dx*dx*g_R*g_T0)
 SoluteRHS = (SoluteRHS+PhiSoluteRHS)/(g_Dch*dx*dx*g_R*g_T0) 

return

end function SoluteRHS

double precision function PhiSoluteRHS(vbles,dx,c_c)
! this uses volume element method. This guarantees conservation of solute. The main source of solute is the variation of free energy with phase: d (d F/dc) = ...+ (d^2 F/dc d phi) d phi + ..., which is revealed by a plot of c_dot.
use solution_parameters
implicit none
double precision, intent (in):: dx,c_c,vbles(N_unk+1,3,3)
integer ii,jj,ip,ih,jv
double precision c,T,Phi(N_phases),FreeEnergy,D
double precision, dimension(3,3)::stencil11=reshape((/0.,6.,0.,6.,-24.,6.,0.,6.,0./),(/3,3/))
double precision, dimension(3)::stencil1=(/-.5,0.,.5/)
double precision :: potential,D_c,f_c,SoluteRHS
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

 PhiSoluteRHS = SoluteRHS/(g_Dch*dx*dx)

return

end function PhiSoluteRHS

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
double precision, intent(in)::dx,c,T,Phi(N_phases),vbles(N_unk+1,3,3)
integer, intent(in):: lp
double precision :: LapPhi(N_phases)
double precision GradientEnergy,potential,FreeEnergy,PhiSoluteRHS
double precision ::M_tilde 
double precision :: c_dot=0,phi_dot=0
double precision :: MolPerVol,beta

MolPerVol=(1.-c)*g_MolPerVol(1)+c*g_MolPerVol(2)

 !assemble and non-dimensionlise phaseRHS 

 M_tilde = (1.-c)+c*g_M(2)/g_M(1)
 phaseRHS=- M_tilde*(GradientEnergy(vbles,lp,dx)+potential(Phi,c,lp)/g_lambda**2&
      +FreeEnergy(c,T,Phi,c_dot,phi_dot,1)*MolPerVol/(g_W(1)*g_lambda**2))
! if(cross_term)then
!    beta = g_beta/g_lambda !
!    phaseRHS = phaseRHS + g_sym* beta*PhiSoluteRHS(vbles,dx,c)/(g_R*g_T0) 
! endif
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
double precision :: c,T,phi(N_phases),c_dot,phi_dot(N_phases)
double precision phi_rhs(N_phases),rfactor,rf4,rf5,rf6,S
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

  ip=N_unk+1
  do ii=1,3
     do jj=1,3
        if(mode_1.eq.1) then
           vbles(ip,ii,jj)=(unk(1,i+ii-2,j+jj-2,1,lb)-&
                unk(2,i+ii-2,j+jj-2,1,lb))/dt
        else
           rfactor = dt/dtold
           rf4 = rfactor*rfactor/(rfactor+1.0)
           rf5 = (rfactor+1.0)
           rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
           vbles(ip,ii,jj)=(rf4*unk(5,i+ii-2,j+jj-2,1,lb)-&
                rf5*unk(4,i+ii-2,j+jj-2,1,lb)+&
                rf6*unk(1,i+ii-2,j+jj-2,1,lb))/dt
        endif
           c=vbles(N_unk,ii,jj)
           vbles(ip,ii,jj)=-vbles(ip,ii,jj)*g_W(1)/((1.-c)+c*g_M(2)/g_M(1))
           MolPerVol=g_MolPerVol(1)*(1.-c)+g_MolPerVol(2)*c
           vbles(ip,ii,jj)=vbles(ip,ii,jj)/MolPerVol*g_lambda**2
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

!  unk(2+3*nunkvbles,i,j,1,lb)=g_anti_trap
!  unk(3+3*nunkvbles,i,j,1,lb)=g_c_force
!
  unk(1+N_unk*6,i,j,1,lb)=Q_heat(1)
  unk(2+N_unk*6,i,j,1,lb)=Q_heat(2)
  unk(3+N_unk*6,i,j,1,lb)=Q_heat(3)
  unk(4+N_unk*6,i,j,1,lb)=Q_heat(4)


 
end subroutine get_RHS_Fi_2d

subroutine get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_)

  use solution_parameters

  implicit none

!
integer, intent(in)::mode_1,nunkvbles,lnblocks
double precision, intent(in)::dx,dt,dtold,unk_(M_unk)
double precision, intent(out):: Fi(N_unk)
double precision, intent(in):: vbles(N_unk+1,3,3)
double precision :: c,T,phi(N_phases),c_dot,phi_dot(N_phases)
double precision rfactor,rf4,rf5,rf6,S,GradientEnergy,potential
double precision TemperatureRHS,SoluteRHS,PhaseRHS,T_k,phi_,FreeEnergy,MolPerVol
logical Not_number
integer ip,n

Fi=0d0

  if (mode_1.eq.1) then
     phi_dot(1) = (vbles(1,2,2)-unk_(2))/dt
     n=2*nunkvbles
     c_dot = (vbles(3,2,2)-unk_(2+n))/dt
  elseif(mode_1.eq.2)then
     rfactor = dt/dtold

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
  if(thermal)then
     Fi(2) = TemperatureRHS(vbles,dx,c,T,phi,c_dot,phi_dot)
  else
     Fi(2)=0d0
  endif
  Fi(3) = SoluteRHS(vbles,dx,c)
!  if(grad_DW_term) then
!     MolPerVol=(1.-c)*g_MolPerVol(1)+c*g_MolPerVol(2)
!     Fi(4) = (g_lambda**2*GradientEnergy(vbles,1,dx)+potential(phi,c,1))/MolPerVol*g_W(1) +FreeEnergy(c,T,Phi,c_dot,phi_dot,1)
!  else
!     Fi(4)=0d0
!  endif

  do ip=1,N_unk
     n=(ip-1)*nunkvbles
     if(mode_1.le.1)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)-unk_(n+2)
     else if(mode_1.eq.2)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)/rf6-unk_(n+2)
     end if
  enddo

!  Fi(4)=vbles(4,2,2)-Fi(4)



end subroutine get_RHS_Fi

subroutine Calc_GradPhi(LapPhi,vbles,dx)
use solution_parameters
implicit none
double precision, intent(in)::dx,vbles(N_unk+1,3,3)
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


double precision function FreeEnergy(c0,T,phi,c_dot,phi_dot,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c0,phi(N_phases),c_dot,phi_dot(N_phases)
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(2),hb(2),hr(2)
double precision :: Tk,Cp,R
double precision :: FE,c
double precision Gibbs_FE_liq,Gibbs_FE_sol,Gibbs_FE_l,Gibbs_FE_e,T_K,g_phi,EntropyMixing !functions
integer :: i,j
double precision :: a = 1d1


!if(c0.lt.1d0-1d0/a)then
!   c = c0
!else
!   c =1d0-1d0/a*exp(-1d0+a*(1d0-c0))
!endif
c=c0



if(AlNi.eq.1)then
   write(*,*)"AlNi flag must be 2 for AlSi"
   stop
   if(lp.eq.0)then !free energy proper J/mol
      
      FreeEnergy =              Gibbs_FE_liq(c,T,0,0)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_FE_sol(c,T,0,0)*g_phi(1d0-phi(1),0)
      
      return
   elseif(lp.eq.1)then !phase
      
      FreeEnergy =              Gibbs_FE_liq(c,T,0,0)*g_phi(phi(1),1)
      FreeEnergy = FreeEnergy - Gibbs_FE_sol(c,T,0,0)*g_phi(phi(1),1)
      
      return
   elseif(lp.eq.3)then !Solute
      
      FreeEnergy = FreeEnergy + Gibbs_FE_liq(c,T,1,0)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_FE_sol(c,T,1,0)*(1.-g_phi(phi(1),0))
      
      return
   elseif(lp.eq.2)then !Temperature Keep in J/mol
      Tk=T_k(T)
      
      !-(1-T*d/dT)dF/dphi
      if(heat_phi_term)then
         !  T (d/dT)dF/dphi
         FreeEnergy =  Tk*(Gibbs_FE_liq(c,T,0,1) - Gibbs_FE_sol(c,T,0,1))
         ! -Df/Dphi
         FreeEnergy = FreeEnergy - Gibbs_FE_liq(c,T,0,0) + Gibbs_FE_sol(c,T,0,0)
         
         FreeEnergy = FreeEnergy*phi_dot(1)*g_phi(phi(1),1)
         
         Q_heat(1) = FreeEnergy
      endif !end phi_dot_term
      if(heat_c_term)then 
         !   -(1- T(d/dT))dF/dc 
         FE=(Tk*Gibbs_FE_liq(c,T,1,1)*g_phi(phi(1),0)- Gibbs_FE_liq(c,T,1,0))*g_phi(phi(1),0)
         FE=FE+(Tk*Gibbs_FE_sol(c,T,1,1)- Gibbs_FE_sol(c,T,1,0))*(1.-g_phi(phi(1),0))
         FE=FE*c_dot
         FreeEnergy = FreeEnergy + FE
         Q_heat(2)= FE
         !      FreeEnergy = FreeEnergy + Tk*Gibbs_FE_liq(c,T,1,1)*g_phi(phi(1),0)*c_dot
         !      FreeEnergy = FreeEnergy + Tk*Gibbs_FE_sol(c,T,1,1)*(1.-g_phi(phi(1),0))*c_dot
         !      FreeEnergy = FreeEnergy - Gibbs_FE_liq(c,T,1,0)*g_phi(phi(1),0)*c_dot
         !      FreeEnergy = FreeEnergy - Gibbs_FE_sol(c,T,1,0)*(1.-g_phi(phi(1),0))*c_dot
      endif !c_dot_term
      !   heat capacity = -T (d/dT) dF/dT
      g_C =  -Tk*(g_phi(phi(1),0)*Gibbs_FE_liq(c,T,0,2)+g_phi(1d0-phi(1),0)*Gibbs_FE_sol(c,T,0,2))
      return
   else
      write(*,*)"bad current_var"
      stop
   endif
elseif(AlNi.eq.2)then
   if(lp.eq.0)then !free energy proper J/mol
      
      FreeEnergy =              Gibbs_FE_l(c,T,0,0)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_FE_e(c,T,0,0)*g_phi(1d0-phi(1),0)
      
      return
   elseif(lp.eq.1)then !phase
      
      FreeEnergy =              Gibbs_FE_l(c,T,0,0)*g_phi(phi(1),1)
      FreeEnergy = FreeEnergy - Gibbs_FE_e(c,T,0,0)*g_phi(phi(1),1)
      
      return
   elseif(lp.eq.3)then !Solute
      
      FreeEnergy = FreeEnergy + Gibbs_FE_l(c,T,1,0)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_FE_e(c,T,1,0)*(1.-g_phi(phi(1),0))
      
      return
   elseif(lp.eq.2)then !Temperature Keep in J/mol
      Tk=T_k(T)
      
      !-(1-T*d/dT)dF/dphi
      if(heat_phi_term)then
         !  T (d/dT)dF/dphi
         FreeEnergy =  Tk*(Gibbs_FE_l(c,T,0,1) - Gibbs_FE_e(c,T,0,1))
         ! -Df/Dphi
         FreeEnergy = FreeEnergy - Gibbs_FE_l(c,T,0,0) + Gibbs_FE_e(c,T,0,0)
         
         FreeEnergy = FreeEnergy*phi_dot(1)*g_phi(phi(1),1)
         
         Q_heat(1) = FreeEnergy
      endif !end phi_dot_term
      if(heat_c_term)then 
         !   -(1- T(d/dT))dF/dc 
         FE=(Tk*Gibbs_FE_l(c,T,1,1)*g_phi(phi(1),0)- Gibbs_FE_l(c,T,1,0))*g_phi(phi(1),0)
         FE=FE+(Tk*Gibbs_FE_e(c,T,1,1)- Gibbs_FE_e(c,T,1,0))*(1.-g_phi(phi(1),0))
         FE=FE*c_dot
         FreeEnergy = FreeEnergy + FE
         Q_heat(2)= FE
         !      FreeEnergy = FreeEnergy + Tk*Gibbs_FE_l(c,T,1,1)*g_phi(phi(1),0)*c_dot
         !      FreeEnergy = FreeEnergy + Tk*Gibbs_FE_e(c,T,1,1)*(1.-g_phi(phi(1),0))*c_dot
         !      FreeEnergy = FreeEnergy - Gibbs_FE_l(c,T,1,0)*g_phi(phi(1),0)*c_dot
         !      FreeEnergy = FreeEnergy - Gibbs_FE_e(c,T,1,0)*(1.-g_phi(phi(1),0))*c_dot
      endif !c_dot_term
      !   heat capacity = -T (d/dT) dF/dT
      g_C =  -Tk*(g_phi(phi(1),0)*Gibbs_FE_l(c,T,0,2)+g_phi(1d0-phi(1),0)*Gibbs_FE_e(c,T,0,2))
      return
   else
      write(*,*)"bad current_var"
      stop
   endif
else
   write(*,*)"models: 1,2 "
   stop
endif
end function FreeEnergy



double precision function EntropyMixing(c,d)
double precision, intent (in)::c
integer, intent (in)::d


!if(c.eq.1d0.or.c.eq.0d0)then
!   write(*,*)"c=",c, "in EntropyMixing"
!   stop
!endif
 



 if(d.eq.1)then
  EntropyMixing =  log(c/(1.-c+g_alpha))
 else
  EntropyMixing = c*log(c)+(1.-c+g_alpha)*log(1.-c+g_alpha)
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
 
subroutine get_y1_y2(y,d,y1_y2)
integer, intent(in)::d
double precision, intent(in) :: y(2)
double precision, intent(out) :: y1_y2(6)



y1_y2=0d0
if(d.eq.0)then
   y1_y2(1)=y(1)*y(2)
   y1_y2(2)=(1d0-y(1))*y(2)
   y1_y2(3)=(1d0-y(2))*y(1)
   y1_y2(4)=(1d0-y(1))*(1d0-y(2))
   y1_y2(5)=y(1)*(1d0-y(1))
   y1_y2(6)=y(2)*(1d0-y(2))
elseif(d.eq.1)then
   y1_y2(1)=y(2)
   y1_y2(2)=-y(2)
   y1_y2(3)=(1d0-y(2))
   y1_y2(4)=-(1d0-y(2))
   y1_y2(5)=1d0-2d0*y(1)
   y1_y2(6)=0d0
elseif(d.eq.2)then
   y1_y2(1)=y(1)
   y1_y2(2)=(1d0-y(1))
   y1_y2(3)=-y(1)
   y1_y2(4)=-(1d0-y(1))
   y1_y2(5)=0d0
   y1_y2(6)=1d0-2d0*y(2)
else
   write(*,*)"d=0 to 2"
   stop
endif
end subroutine get_y1_y2




double precision function y1(c,d)
integer, intent(in)::d
double precision, intent(in) :: c
if(c.gt.1d0-1d-9.or.c.lt.1d-9)then
   write(*,*)"c out of bound in y1()"
   stop
endif
if(d.eq.0)then
   if(c.ge.0.4d0)then
      y1=1d0-1d-5
      return
   else
      y1=2.5*c
      return
   endif
elseif(d.eq.1)then
   if(c.ge.0.4d0)then
      y1=1d-5
      return
   else
      y1=2.5
      return
!   else
!      write(*,*)"d y1/dc not determined at c=",c
!      stop
   endif
else
   write(*,*)"d =0 or 1"
endif
end function y1

double precision function y2(c,d)
integer, intent(in)::d
double precision, intent(in) :: c
if(c.gt.1d0-1d-9.or.c.lt.1d-9)then
   write(*,*)"c out of bound in y2()"
   stop
endif

if(d.eq.0)then
   if(c.ge.0.5d0)then
      y2=1d-5
      return
   elseif(c.gt.0.4d0)then
      y2=1d0-(2d0-5d0*c)/(c-1d0)
      return
   else
      y2=1d0-1d-5
      return
   endif
elseif(d.eq.1)then
!   if(c.eq.0.5.or.c.eq.0.4)then
!       write(*,*)"d y2/dc not determined at c =",c
!       stop
   if(c.gt.0.5d0)then
      y2=1d-5
      return
   elseif(c.gt.0.4d0)then
      y2=-3d0/(1d0-c)**2
      return
   else
      y2=1d-5
      return
   endif
else
   write(*,*)"d =0 or 1"
endif
end function y2





double precision function y3(c,d)
integer, intent(in)::d
double precision, intent(in) :: c
if(c.gt.1d0-1d-9.or.c.lt.1d-9)then
   write(*,*)"c out of bound in y3()"
   stop
endif
if(d.eq.0)then
   if(c.gt.0.5d0)then
      y3=2d0*(1d0-c)
      return
   else
      y3=1d0-1d-5
      return
   endif
elseif(d.eq.1)then
   if(c.gt.0.5d0)then
      y3=-2d0
      return
   else
      y3=0d0
      return
!   else
!      write(*,*)"d y3/dc not determined at c=",c
!      stop
   endif
else
   write(*,*)"d =0 or 1"
endif
end function y3


double precision function y4(c,d)
integer, intent(in)::d
double precision, intent(in) :: c
if(c.gt.1d0-1d-9.or.c.lt.1d-9)then
   write(*,*)"c out of bound in y4()"
   stop
endif
if(d.eq.0)then
   if(c.ge.0.5d0)then
      y4=1d0-1d-5
      return
   else
      y4=c/(1d0-c)
      return
   endif
elseif(d.eq.1)then
   if(c.ge.0.5d0)then
      y4=0d0
      return
   else
      y4=1d0/(1d0-c)**2
      return
!   else
!      write(*,*)"d y4/dc not determined at c=",c
!      stop
   endif
else
   write(*,*)"d =0 or 1"
endif
end function y4



double precision function gibbs_T(h,T,d)
double precision, intent(in) :: T,h(8)
integer, intent(in) :: d
double precision V(8)
integer i
gibbs_T=0d0
call GibbsVector(V,T,d)
do i= 1,8
   gibbs_T=gibbs_T+h(i)*V(i)
enddo
end function gibbs_T


double precision function Gibbs_FE_l(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5), gs1A,gs1B, gs2A,gs2B

  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gl_,Tk,Gibbs_T
  double precision :: hlA_(8),hlB_(8),hs1A_(8),hs1B_(8),hs2A_(8),hs2B_(8)


! AlSi below




Tk=T_k(T)
if (Tk .gt. 1687.0) then
  hlA_=hlA3 
  hlB_=hlB2 
elseif(Tk .gt. 933.47) then
  hlA_=hlA3 
  hlB_=hlB1 
elseif (Tk .gt. 700.0) then
  hlA_=hlA2 
  hlB_=hlB1 
else
  hlA_=hlA1 
  hlB_=hlB1 
endif

if(d_c.eq.0.and.d_T.eq.0)then
   glA=gibbs_T(hlA_,Tk,0)
   glB=gibbs_T(hlB_,Tk,0)
   call MixingVector(RK,c,0)
   RKl_=         (hrkl0(1)+hrkl0(2)*Tk)*RK(1)
   RKl_=RKl_ + (hrkl1(1)+hrkl1(2)*Tk)*RK(2)
   RKl_=RKl_ + (hrkl2(1)+hrkl2(2)*Tk)*RK(3)
   
   ent_=g_R*Tk*EntropyMixing(c,0)
   gl_=glA*(1d0-c)+glB*c+ent_+RKl_
   Gibbs_FE_l=gl_
elseif(d_c.eq.0.and.d_T.eq.1)then
   glA=gibbs_T(hlA_,Tk,1)
   glB=gibbs_T(hlB_,Tk,1)
   call MixingVector(RK,c,0)
   RKl_=         hrkl0(2)*RK(1)
   RKl_=RKl_ + hrkl1(2)*RK(2)
   RKl_=RKl_ + hrkl2(2)*RK(3)
   
   ent_=g_R*EntropyMixing(c,0)
   gl_=glA*(1d0-c)+glB*c+ent_+RKl_
   Gibbs_FE_l=gl_
elseif(d_c.eq.0.and.d_T.eq.2)then
   glA=gibbs_T(hlA_,Tk,2)
   glB=gibbs_T(hlB_,Tk,2)
   gl_=glA*(1d0-c)+glB*c
   Gibbs_FE_l=gl_
elseif(d_c.eq.1.and.d_T.eq.0)then
   glA=gibbs_T(hlA_,Tk,0)
   glB=gibbs_T(hlB_,Tk,0)
   call MixingVector(RK,c,1)
   RKl_=         (hrkl0(1)+hrkl0(2)*Tk)*RK(1)
   RKl_=RKl_ + (hrkl1(1)+hrkl1(2)*Tk)*RK(2)
   RKl_=RKl_ + (hrkl2(1)+hrkl2(2)*Tk)*RK(3)
   
   ent_=g_R*Tk*EntropyMixing(c,1)
   gl_=-glA+glB+ent_+RKl_
   Gibbs_FE_l=gl_
elseif(d_c.eq.1.and.d_T.eq.1)then
   glA=gibbs_T(hlA_,Tk,1)
   glB=gibbs_T(hlB_,Tk,1)
   call MixingVector(RK,c,1)
   RKl_=       hrkl0(2)*RK(1)
   RKl_=RKl_ + hrkl1(2)*RK(2)
   RKl_=RKl_ + hrkl2(2)*RK(3)
   
   ent_=g_R*EntropyMixing(c,1)
   gl_=-glA+glB+ent_+RKl_
   Gibbs_FE_l=gl_
elseif(d_c.eq.1.and.d_T.eq.2)then
   glA=gibbs_T(hlA_,Tk,2)
   glB=gibbs_T(hlB_,Tk,2)
   gl_=-glA+glB
   Gibbs_FE_l=gl_
else
   write(*,*)"unexpected d_c,d_T combination",d_c,d_T
   stop
endif

end function Gibbs_FE_l


double precision function Gibbs_FE_e(c,T,d_c,d_T)
  use solution_parameters
  use time_dep_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5), gs1A,gs1B, gs2A,gs2B

  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gs2_,Tk,Gibbs_T
  double precision :: hlA_(8),hlB_(8),hs1A_(8),hs1B_(8),hs2A_(8),hs2B_(8)
  double precision :: t_factor



! AlSi below




Tk=T_k(T)
if (Tk .gt. 1687.0) then
  hs2A_=hs2A3 
  hs2B_=hs2B2 
elseif(Tk .gt. 933.47) then
  hs2A_=hs2A3 
  hs2B_=hs2B1 
elseif (Tk .gt. 700.0) then
  hs2A_=hs2A2 
  hs2B_=hs2B1 
else
  hs2A_=hs2A1 
  hs2B_=hs2B1 
endif

if(d_c.eq.0.and.d_T.eq.0)then
   gs2A=gibbs_T(hs2A_,Tk,0)
   gs2B=gibbs_T(hs2B_,Tk,0)
   call MixingVector(RK,c,0)
   RKl_=       (hrks20(1)+hrks20(2)*Tk)*RK(1)
   RKl_=RKl_ + (hrks21(1)+hrks21(2)*Tk)*RK(2)
   RKl_=RKl_ + (hrks22(1)+hrks22(2)*Tk)*RK(3)
   ent_=g_R*Tk*EntropyMixing(c,0)
   gs2_=gs2A*(1d0-c)+gs2B*c+ent_+RKl_
   Gibbs_FE_e=gs2_
elseif(d_c.eq.0.and.d_T.eq.1)then
   gs2A=gibbs_T(hs2A_,Tk,1)
   gs2B=gibbs_T(hs2B_,Tk,1)
   call MixingVector(RK,c,0)
   RKl_=       hrks20(2)*RK(1)
   RKl_=RKl_ + hrks21(2)*RK(2)
   RKl_=RKl_ + hrks22(2)*RK(3)
   ent_=g_R*EntropyMixing(c,0)
   gs2_=gs2A*(1d0-c)+gs2B*c+ent_+RKl_
   Gibbs_FE_e=gs2_
elseif(d_c.eq.0.and.d_T.eq.2)then
   gs2A=gibbs_T(hs2A_,Tk,2)
   gs2B=gibbs_T(hs2B_,Tk,2)
   gs2_=gs2A*(1d0-c)+gs2B*c
   Gibbs_FE_e=gs2_
elseif(d_c.eq.1.and.d_T.eq.0)then
   gs2A=gibbs_T(hs2A_,Tk,0)
   gs2B=gibbs_T(hs2B_,Tk,0)
   call MixingVector(RK,c,1)
   RKl_=       (hrks20(1)+hrks20(2)*Tk)*RK(1)
   RKl_=RKl_ + (hrks21(1)+hrks21(2)*Tk)*RK(2)
   RKl_=RKl_ + (hrks22(1)+hrks22(2)*Tk)*RK(3)
   ent_=g_R*Tk*EntropyMixing(c,1)
   gs2_=-gs2A+gs2B+ent_+RKl_
   Gibbs_FE_e=gs2_
elseif(d_c.eq.1.and.d_T.eq.1)then
   gs2A=gibbs_T(hs2A_,Tk,1)
   gs2B=gibbs_T(hs2B_,Tk,1)
   call MixingVector(RK,c,1)
   RKl_=       hrks20(2)*RK(1)
   RKl_=RKl_ + hrks21(2)*RK(2)
   RKl_=RKl_ + hrks22(2)*RK(3)
   ent_=g_R*EntropyMixing(c,1)
   gs2_=-gs2A+gs2B+ent_+RKl_
   Gibbs_FE_e=gs2_
elseif(d_c.eq.1.and.d_T.eq.2)then
   gs2A=gibbs_T(hs2A_,Tk,2)
   gs2B=gibbs_T(hs2B_,Tk,2)
   gs2_=-gs2A+gs2B
   Gibbs_FE_e=gs2_
else
   write(*,*)"unexpected d_c,d_T combination",d_c,d_T
   stop
endif

end function Gibbs_FE_e






double precision function Gibbs_fe_eX(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in)::d_c,d_T
  double precision :: SLe(6),yy(2),yy_d(2),y1_y2(6),y1_y2_d1(6),y1_y2_d2(6)
  double precision :: y3,y4,EntropyMixing,T_k
  double precision ::entSL_,entSL_d,ge1_,ge1_d1,ge1_d2,ge2_,ge2_d1,ge2_d2,ge_,ge_d,Tk,Gibbs_T
  double precision :: hlA_(8),hlB_(8),haA_(8),haB_(8),hslb1_(8),hslb2_(8),hslb3_(8),hslb4_(8),hc_(8),hsle1_(8),hsle2_(8),hsle3_(8),hsle4_(8)
  
 Tk=T_k(T)
  if (Tk > 1728.0) then
  hlA_=hlA3 
  hlB_=hlB2 
  haA_=haA3 
  haB_=haB2 
  hslb1_=hslb1_4 
  hslb2_=hslb2_3 
  hslb3_=hslb3_4 
  hslb4_=hslb4_4 
  hc_=hc4 
  hsle1_=hsle1_4 
  hsle2_=hsle2_3 
  hsle3_=hsle3_2 
  hsle4_=hsle4_2 
elseif (Tk > 933.47)then
  hlA_=hlA3 
  hlB_=hlB1 
  haA_=haA3 
  haB_=haB1 
  hslb1_=hslb1_3 
  hslb2_=hslb2_3 
  hslb3_=hslb3_3 
  hslb4_=hslb4_3 
  hc_=hc3 
  hsle1_=hsle1_3 
  hsle2_=hsle2_3 
  hsle3_=hsle3_1 
  hsle4_=hsle4_1 
elseif (Tk > 700.0)then
  hlA_=hlA2 
  hlB_=hlB1 
  haA_=haA2 
  haB_=haB1 
  hslb1_=hslb1_2 
  hslb2_=hslb2_2 
  hslb3_=hslb3_2 
  hslb4_=hslb4_2 
  hc_=hc2 
  hsle1_=hsle1_2 
  hsle2_=hsle2_2 
  hsle3_=hsle3_1 
  hsle4_=hsle4_1 
else
  hlA_=hlA1 
  hlB_=hlB1 
  haA_=haA1 
  haB_=haB1 
  hslb1_=hslb1_1 
  hslb2_=hslb2_1 
  hslb3_=hslb3_1 
  hslb4_=hslb4_1 
  hc_=hc1 
  hsle1_=hsle1_1 
  hsle2_=hsle2_1 
  hsle3_=hsle3_1 
  hsle4_=hsle4_1 
endif

  if(d_c.eq.0.and.d_T.eq.0)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(hsle1_,Tk,0)
     SLe(2)=gibbs_T(hsle2_,Tk,0)
     SLe(3)=gibbs_T(hsle3_,Tk,0)
     SLe(4)=gibbs_T(hsle4_,Tk,0)
     SLe(5)=hsle5(1)+hsle5(2)*Tk
     SLe(6)=hsle6(1)+hsle6(2)*Tk

     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     entSL_=g_R*Tk*(EntropyMixing(yy(1),0)+EntropyMixing(yy(2),0))
     call get_y1_y2(yy,0,y1_y2)
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge2_=SLe(5)*y1_y2(5)+SLe(6)*y1_y2(6)
     ge_=(ge1_+ge2_+entSL_)/(1d0+yy(1))
     Gibbs_fe_eX=ge_
  elseif(d_c.eq.1.and.d_T.eq.0)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(hsle1_,Tk,0)
     SLe(2)=gibbs_T(hsle2_,Tk,0)
     SLe(3)=gibbs_T(hsle3_,Tk,0)
     SLe(4)=gibbs_T(hsle4_,Tk,0)
     SLe(5)=hsle5(1)+hsle5(2)*Tk
     SLe(6)=hsle6(1)+hsle6(2)*Tk

     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     yy_d(1)=y4(c,1)
     yy_d(2)=y3(c,1)
     entSL_=g_R*Tk*(EntropyMixing(yy(1),0)+EntropyMixing(yy(2),0))
     entSL_d=g_R*Tk*(EntropyMixing(yy(1),1)*yy_d(1)+EntropyMixing(yy(2),1)*yy_d(2))


     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1)
     call get_y1_y2(yy,2,y1_y2_d2)
     
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge2_=SLe(5)*y1_y2(5)+SLe(6)*y1_y2(6)
     ge1_d1=(SLe(1)*y1_y2_d1(1)+SLe(2)*y1_y2_d1(2)+SLe(3)*y1_y2_d1(3)+SLe(4)*y1_y2_d1(4))*yy_d(1)
     ge1_d2=(SLe(1)*y1_y2_d2(1)+SLe(2)*y1_y2_d2(2)+SLe(3)*y1_y2_d2(3)+SLe(4)*y1_y2_d2(4))*yy_d(2)
     ge2_d1=(SLe(5)*y1_y2_d1(5)+SLe(6)*y1_y2_d1(6))*yy_d(1)
     ge2_d2=(SLe(5)*y1_y2_d2(5)+SLe(6)*y1_y2_d2(6))*yy_d(2)
     ge_=(ge1_+ge2_+entSL_)/(1d0+yy(1))
     ge_d=-ge_/(1d0+yy(1))*yy_d(1)+(ge1_d1+ge1_d2+ge2_d1+ge2_d2+entSL_d)/(1d0+yy(1))
     Gibbs_fe_eX=ge_d
  elseif(d_c.eq.1.and.d_T.eq.1)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(hsle1_,Tk,1)
     SLe(2)=gibbs_T(hsle2_,Tk,1)
     SLe(3)=gibbs_T(hsle3_,Tk,1)
     SLe(4)=gibbs_T(hsle4_,Tk,1)
     SLe(5)=hsle5(2)
     SLe(6)=hsle6(2)

     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     yy_d(1)=y4(c,1)
     yy_d(2)=y3(c,1)
     entSL_=g_R*(EntropyMixing(yy(1),0)+EntropyMixing(yy(2),0))
     entSL_d=g_R*(EntropyMixing(yy(1),1)*yy_d(1)+EntropyMixing(yy(2),1)*yy_d(2))


     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1)
     call get_y1_y2(yy,2,y1_y2_d2)
     
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge2_=SLe(5)*y1_y2(5)+SLe(6)*y1_y2(6)
     ge1_d1=(SLe(1)*y1_y2_d1(1)+SLe(2)*y1_y2_d1(2)+SLe(3)*y1_y2_d1(3)+SLe(4)*y1_y2_d1(4))*yy_d(1)
     ge1_d2=(SLe(1)*y1_y2_d2(1)+SLe(2)*y1_y2_d2(2)+SLe(3)*y1_y2_d2(3)+SLe(4)*y1_y2_d2(4))*yy_d(2)
     ge2_d1=(SLe(5)*y1_y2_d1(5)+SLe(6)*y1_y2_d1(6))*yy_d(1)
     ge2_d2=(SLe(5)*y1_y2_d2(5)+SLe(6)*y1_y2_d2(6))*yy_d(2)
     ge_=(ge1_+ge2_+entSL_)/(1d0+yy(1))
     ge_d=-ge_/(1d0+yy(1))*yy_d(1)+(ge1_d1+ge1_d2+ge2_d1+ge2_d2+entSL_d)/(1d0+yy(1))
     Gibbs_fe_eX=ge_d
  elseif(d_c.eq.1.and.d_T.eq.2)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(hsle1_,Tk,2)
     SLe(2)=gibbs_T(hsle2_,Tk,2)
     SLe(3)=gibbs_T(hsle3_,Tk,2)
     SLe(4)=gibbs_T(hsle4_,Tk,2)


     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     yy_d(1)=y4(c,1)
     yy_d(2)=y3(c,1)



     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1)
     call get_y1_y2(yy,2,y1_y2_d2)
     
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge1_d1=(SLe(1)*y1_y2_d1(1)+SLe(2)*y1_y2_d1(2)+SLe(3)*y1_y2_d1(3)+SLe(4)*y1_y2_d1(4))*yy_d(1)
     ge1_d2=(SLe(1)*y1_y2_d2(1)+SLe(2)*y1_y2_d2(2)+SLe(3)*y1_y2_d2(3)+SLe(4)*y1_y2_d2(4))*yy_d(2)
     ge_=ge1_/(1d0+yy(1))
     ge_d=-ge_/(1d0+yy(1))*yy_d(1)+(ge1_d1+ge1_d2)/(1d0+yy(1))
     Gibbs_fe_eX=ge_d
  elseif(d_c.eq.0.and.d_T.eq.1)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(hsle1_,Tk,1)
     SLe(2)=gibbs_T(hsle2_,Tk,1)
     SLe(3)=gibbs_T(hsle3_,Tk,1)
     SLe(4)=gibbs_T(hsle4_,Tk,1)
     SLe(5)=hsle5(2)
     SLe(6)=hsle6(2)

     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     entSL_=g_R*(EntropyMixing(yy(1),0)+EntropyMixing(yy(2),0))
     call get_y1_y2(yy,0,y1_y2)
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge2_=SLe(5)*y1_y2(5)+SLe(6)*y1_y2(6)
     ge_=(ge1_+ge2_+entSL_)/(1d0+yy(1))
     Gibbs_fe_eX=ge_
  elseif(d_c.eq.0.and.d_T.eq.2)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(hsle1_,Tk,2)
     SLe(2)=gibbs_T(hsle2_,Tk,2)
     SLe(3)=gibbs_T(hsle3_,Tk,2)
     SLe(4)=gibbs_T(hsle4_,Tk,2)
     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     call get_y1_y2(yy,0,y1_y2)
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge_=ge1_/(1d0+yy(1))
     Gibbs_fe_eX=ge_
  endif
end function Gibbs_fe_eX

double precision function Gibbs_FE_liq(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5)
  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gl_,Tk,Gibbs_T

  
  if(d_c.eq.0.and.d_T.eq.0)then
     Tk=T_k(T)
     glA=gibbs_T(hlA,Tk,0)
     glB=gibbs_T(hlB,Tk,0)
     call MixingVector(RK,c,0)
     RKl_=         (hrkl0(1)+hrkl0(2)*Tk)*RK(1)
     RKl_=RKl_ + (hrkl1(1)+hrkl1(2)*Tk)*RK(2)
     RKl_=RKl_ + (hrkl2(1)+hrkl2(2)*Tk)*RK(3)
     RKl_=RKl_ + (hrkl3(1)+hrkl3(2)*Tk)*RK(4)
     RKl_=RKl_ + (hrkl4(1)+hrkl4(2)*Tk)*RK(5)
     ent_=g_R*Tk*EntropyMixing(c,0)
     gl_=glA*(1d0-c)+glB*c+ent_+RKl_
     Gibbs_FE_liq=gl_
  elseif(d_c.eq.0.and.d_T.eq.1)then
     Tk=T_k(T)
     glA=gibbs_T(hlA,Tk,1)
     glB=gibbs_T(hlB,Tk,1)
     call MixingVector(RK,c,0)
     RKl_=         hrkl0(2)*RK(1)
     RKl_=RKl_ +   hrkl1(2)*RK(2)
     RKl_=RKl_ +   hrkl2(2)*RK(3)
     RKl_=RKl_ +   hrkl3(2)*RK(4)
     RKl_=RKl_ +   hrkl4(2)*RK(5)
     ent_=g_R*EntropyMixing(c,0)
     gl_=glA*(1d0-c)+glB*c+ent_+RKl_
     Gibbs_FE_liq=gl_
  elseif(d_c.eq.0.and.d_T.eq.2)then
     Tk=T_k(T)
     glA=gibbs_T(hlA,Tk,2)
     glB=gibbs_T(hlB,Tk,2)
     call MixingVector(RK,c,0)
     gl_=glA*(1d0-c)+glB*c
     Gibbs_FE_liq=gl_
  elseif(d_c.eq.1.and.d_T.eq.0)then
     Tk=T_k(T)
     glA=gibbs_T(hlA,Tk,0)
     glB=gibbs_T(hlB,Tk,0)
     call MixingVector(RK,c,1)
     RKl_=         (hrkl0(1)+hrkl0(2)*Tk)*RK(1)
     RKl_=RKl_ + (hrkl1(1)+hrkl1(2)*Tk)*RK(2)
     RKl_=RKl_ + (hrkl2(1)+hrkl2(2)*Tk)*RK(3)
     RKl_=RKl_ + (hrkl3(1)+hrkl3(2)*Tk)*RK(4)
     RKl_=RKl_ + (hrkl4(1)+hrkl4(2)*Tk)*RK(5)
     ent_=g_R*Tk*EntropyMixing(c,1)
     gl_=-glA+glB+ent_+RKl_
     Gibbs_FE_liq=gl_
  elseif(d_c.eq.1.and.d_T.eq.1)then
     Tk=T_k(T)
     glA=gibbs_T(hlA,Tk,1)
     glB=gibbs_T(hlB,Tk,1)
     call MixingVector(RK,c,0)
     RKl_=         hrkl0(2)*RK(1)
     RKl_=RKl_ +   hrkl1(2)*RK(2)
     RKl_=RKl_ +   hrkl2(2)*RK(3)
     RKl_=RKl_ +   hrkl3(2)*RK(4)
     RKl_=RKl_ +   hrkl4(2)*RK(5)
     ent_=g_R*EntropyMixing(c,1)
     gl_=-glA+glB+ent_+RKl_
     Gibbs_FE_liq=gl_
  elseif(d_c.eq.1.and.d_T.eq.2)then
     Tk=T_k(T)
     glA=gibbs_T(hlA,Tk,2)
     glB=gibbs_T(hlB,Tk,2)
!     call MixingVector(RK,c,0)
     gl_=-glA+glB
     Gibbs_FE_liq=gl_
  else
     write(*,*)"unexpected d_c,d_T combination",d_c,d_T
     stop
  endif

  end function Gibbs_FE_liq




double precision function Gibbs_FE_sol(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in)::d_c,d_T
  double precision :: SLb(6),yy(2),yy_d(2),y1_y2(6),y1_y2_d1(6),y1_y2_d2(6)
  double precision :: y1,y2,EntropyMixing,T_k
  double precision ::entSL_,entSL_d,gb1_,gb1_d1,gb1_d2,gb2_,gb2_d1,gb2_d2,gb_,gb_d,Tk,Gibbs_T
  
  if(d_c.eq.0.and.d_T.eq.0)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(hslb1,Tk,0)
     SLb(2)=gibbs_T(hslb2,Tk,0)
     SLb(3)=gibbs_T(hslb3,Tk,0)
     SLb(4)=gibbs_T(hslb4,Tk,0)
     SLb(5)=hslb5(1)+hslb5(2)*Tk
     SLb(6)=hslb6(1)+hslb6(2)*Tk
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     entSL_=g_R*Tk*(EntropyMixing(yy(1),0)/3.+EntropyMixing(yy(2),0)/6.)
     call get_y1_y2(yy,0,y1_y2)
     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)
     gb2_=SLb(5)*y1_y2(5)+SLb(6)*y1_y2(6)
     gb_=(gb1_+gb2_+entSL_)/(1-yy(2)/6d0)
     Gibbs_FE_sol=gb_
  elseif(d_c.eq.0.and.d_T.eq.1)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(hslb1,Tk,1)
     SLb(2)=gibbs_T(hslb2,Tk,1)
     SLb(3)=gibbs_T(hslb3,Tk,1)
     SLb(4)=gibbs_T(hslb4,Tk,1)
     SLb(5)=hslb5(2)
     SLb(6)=hslb6(2)
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     entSL_=g_R*(EntropyMixing(yy(1),0)/3.+EntropyMixing(yy(2),0)/6.)
     call get_y1_y2(yy,0,y1_y2)
     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)
     gb2_=SLb(5)*y1_y2(5)+SLb(6)*y1_y2(6)
     gb_=(gb1_+gb2_+entSL_)/(1-yy(2)/6d0)
     Gibbs_FE_sol=gb_
  elseif(d_c.eq.0.and.d_T.eq.2)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(hslb1,Tk,2)
     SLb(2)=gibbs_T(hslb2,Tk,2)
     SLb(3)=gibbs_T(hslb3,Tk,2)
     SLb(4)=gibbs_T(hslb4,Tk,2)
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     call get_y1_y2(yy,0,y1_y2)
     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)
     gb_=gb1_/(1-yy(2)/6d0)
     Gibbs_FE_sol=gb_
  elseif(d_c.eq.1.and.d_T.eq.0)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(hslb1,Tk,0)
     SLb(2)=gibbs_T(hslb2,Tk,0)
     SLb(3)=gibbs_T(hslb3,Tk,0)
     SLb(4)=gibbs_T(hslb4,Tk,0)
     SLb(5)=hslb5(1)+hslb5(2)*Tk
     SLb(6)=hslb6(1)+hslb6(2)*Tk
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     yy_d(1)=y1(c,1)
     yy_d(2)=y2(c,1)

     entSL_=g_R*Tk*(EntropyMixing(yy(1),0)/3.+EntropyMixing(yy(2),0)/6.)
     entSL_d=g_R*Tk*(yy_d(1)*EntropyMixing(yy(1),1)/3.+yy_d(2)*EntropyMixing(yy(2),1)/6.)


     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1) !d(y1-y2)/d(y1)
     call get_y1_y2(yy,2,y1_y2_d2) !d(y1-y2)/d(y2)

     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)

     gb1_d1=SLb(1)*y1_y2_d1(1)+SLb(2)*y1_y2_d1(2)+SLb(3)*y1_y2_d1(3)+SLb(4)*y1_y2_d1(4)
     gb1_d1=gb1_d1*yy_d(1) !d(y1)/d(c)*d(gb1)/d(y1)

     gb1_d2=SLb(1)*y1_y2_d2(1)+SLb(2)*y1_y2_d2(2)+SLb(3)*y1_y2_d2(3)+SLb(4)*y1_y2_d2(4)
     gb1_d2=gb1_d2*yy_d(2)!d(y2)/d(c)*d(gb1)/d(y2)


     gb2_=SLb(5)*y1_y2(5)+SLb(6)*y1_y2(6)

     gb2_d1=SLb(5)*y1_y2_d1(5)+SLb(6)*y1_y2_d1(6)
     gb2_d1=gb2_d1*yy_d(1) !d(y1)/d(c)*d(gb2)/d(y1) 

     gb2_d2=SLb(5)*y1_y2_d2(5)+SLb(6)*y1_y2_d2(6)
     gb2_d2=gb2_d2*yy_d(2)!d(y2)/d(c)*d(gb2)/d(y2) 


     gb_=(gb1_+gb2_+entSL_)/(1-yy(2)/6d0)
     gb_d=(gb1_d1+gb1_d2+gb2_d1+gb2_d2+entSL_d)/(1-yy(2)/6d0)-(gb_)/(1-yy(2)/6d0)*(-yy_d(2)/6d0)
     Gibbs_FE_sol=gb_d
  elseif(d_c.eq.1.and.d_T.eq.1)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(hslb1,Tk,1)
     SLb(2)=gibbs_T(hslb2,Tk,1)
     SLb(3)=gibbs_T(hslb3,Tk,1)
     SLb(4)=gibbs_T(hslb4,Tk,1)
     SLb(5)=hslb5(2)
     SLb(6)=hslb6(2)
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     yy_d(1)=y1(c,1)
     yy_d(2)=y2(c,1)

     entSL_=g_R*(EntropyMixing(yy(1),0)/3.+EntropyMixing(yy(2),0)/6.)
     entSL_d=g_R*(yy_d(1)*EntropyMixing(yy(1),1)/3.+yy_d(2)*EntropyMixing(yy(2),1)/6.)


     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1) !d(y1-y2)/d(y1)
     call get_y1_y2(yy,2,y1_y2_d2) !d(y1-y2)/d(y2)

     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)

     gb1_d1=SLb(1)*y1_y2_d1(1)+SLb(2)*y1_y2_d1(2)+SLb(3)*y1_y2_d1(3)+SLb(4)*y1_y2_d1(4)
     gb1_d1=gb1_d1*yy_d(1) !d(y1)/d(c)*d(gb1)/d(y1)

     gb1_d2=SLb(1)*y1_y2_d2(1)+SLb(2)*y1_y2_d2(2)+SLb(3)*y1_y2_d2(3)+SLb(4)*y1_y2_d2(4)
     gb1_d2=gb1_d2*yy_d(2)!d(y2)/d(c)*d(gb1)/d(y2)


     gb2_=SLb(5)*y1_y2(5)+SLb(6)*y1_y2(6)

     gb2_d1=SLb(5)*y1_y2_d1(5)+SLb(6)*y1_y2_d1(6)
     gb2_d1=gb2_d1*yy_d(1) !d(y1)/d(c)*d(gb2)/d(y1) 

     gb2_d2=SLb(5)*y1_y2_d2(5)+SLb(6)*y1_y2_d2(6)
     gb2_d2=gb2_d2*yy_d(2)!d(y2)/d(c)*d(gb2)/d(y2) 


     gb_=(gb1_+gb2_+entSL_)/(1-yy(2)/6d0)
     gb_d=(gb1_d1+gb1_d2+gb2_d1+gb2_d2+entSL_d)/(1-yy(2)/6d0)-(gb_)/(1-yy(2)/6d0)*(-yy_d(2)/6d0)
     Gibbs_FE_sol=gb_d
  elseif(d_c.eq.1.and.d_T.eq.2)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(hslb1,Tk,2)
     SLb(2)=gibbs_T(hslb2,Tk,2)
     SLb(3)=gibbs_T(hslb3,Tk,2)
     SLb(4)=gibbs_T(hslb4,Tk,2)
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     yy_d(1)=y1(c,1)
     yy_d(2)=y2(c,1)

     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1) !d(y1-y2)/d(y1)
     call get_y1_y2(yy,2,y1_y2_d2) !d(y1-y2)/d(y2)

     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)

     gb1_d1=SLb(1)*y1_y2_d1(1)+SLb(2)*y1_y2_d1(2)+SLb(3)*y1_y2_d1(3)+SLb(4)*y1_y2_d1(4)
     gb1_d1=gb1_d1*yy_d(1) !d(y1)/d(c)*d(gb1)/d(y1)

     gb1_d2=SLb(1)*y1_y2_d2(1)+SLb(2)*y1_y2_d2(2)+SLb(3)*y1_y2_d2(3)+SLb(4)*y1_y2_d2(4)
     gb1_d2=gb1_d2*yy_d(2)!d(y2)/d(c)*d(gb1)/d(y2)

     gb_=gb1_/(1-yy(2)/6d0)
     gb_d=gb1_d1/(1-yy(2)/6d0)-(gb_)/(1-yy(2)/6d0)*(-yy_d(2)/6d0)
     Gibbs_FE_sol=gb_d

  else
     write(*,*)"unexpected d_c d_T combination",d_c,d_T
  endif
  
end function Gibbs_FE_sol
