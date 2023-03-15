double precision function SoluteRHS(vbles,dx,c_dot,phi_dot,c_c)
! this uses volume element method. This guarantees conservation of solute. The main source of solute is the variation of free energy with phase: d (d F/dc) = ...+ (d^2 F/dc d phi) d phi + ..., which is revealed by a plot of c_dot.
use solution_parameters
implicit none

double precision, intent (in):: dx,c_dot,c_c,vbles(2,3,3)
!integer, intent (in)::i,j,k,lb,df
integer ii,jj,ip,ih,jv
double precision c,T,Phi,FreeEnergy,D,phi_dot,phi_,phi_tmp
double precision, dimension(3,3)::stencil11=reshape((/0.,6.,0.,6.,-24.,6.,0.,6.,0./),(/3,3/))
double precision, dimension(3)::stencil1=(/-.5,0.,.5/)
!double precision, dimension(3):: D1=(/6.7e-10,6.7e-12,6.7e-12/)!Jackson Hunt
double precision, dimension(3):: D1=(/6.70e-10,6.70e-12,6.70e-12/)
double precision :: DFDc,potential,Lap_c,D_, D_c,D_cc,D_c1,D_cc1,Grad_c(2),y0,g_phi,f_c,f_cc
double precision ::D_ch =  2.83e-10! (m^2/s)    characteristic diffusivity delta2/tau=M*eps2(1)
double precision :: phaseRHS_c,beta=1d1
! 
beta=0.0 !changePCB
!

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
!!Left !
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
