!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  Utilities for phase field stencils in 2-d and 3-d
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





double precision function GradientEnergy(vbles,lp,dx,i,j)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,10,10)
integer, intent (in):: lp,i,j
double precision ::phix,phiy,X,Y,u,v,A
double precision :: G11,G22,G12,M11,M12,g,GG
double precision :: LapPhi(1)

call Calc_GradPhi(LapPhi(1),vbles,dx,i,j)
!!$
!!$phix = 0.5*(vbles(1,3,2) - vbles(1,1,2))/dx 
!!$phiy = 0.5*(vbles(1,2,3) - vbles(1,2,1))/dx 
!!$u = phix*phix
!!$v = phiy*phiy
!!$g  = u + v
!!$if(g.gt.1d-20)then
!!$
!!$   X = u/g
!!$   Y= v/g
!!$   A = A_0*(1d0+epsilon_tilde*(X**2+Y**2))
!!$
!!$   GG = 16*A_0*epsilon_tilde*X*Y*(A_0*epsilon_tilde*(X-Y)**2+A)
!!$   G11 = Y*(GG +4*epsilon_tilde*A*A_0*(X-Y))+ A**2
!!$   G22 = X*(GG +4*epsilon_tilde*A*A_0*(Y-X))+ A**2
!!$   G12 = -phix*phiy*GG/g
!!$
!!$
!!$
!!$   M11 = 0.5*(vbles(1,3,2)+vbles(1,1,2)-vbles(1,2,3)-vbles(1,2,1))/dx**2  !(phi_xx - phi_yy)/2
!!$   M12 = 0.25*(vbles(1,3,3)+vbles(1,1,1)-vbles(1,1,3)-vbles(1,3,1))/dx**2 !phi_xy
!!$
!!$   GradientEnergy = 0.5*LapPhi(1)*(G11+G22)+M11*(G11-G22)+2*M12*G12
!!$
!!$
!!$else
!!$   GradientEnergy = (A_0*(1d0+epsilon_tilde))**2*LapPhi(1)
!!$endif

GradientEnergy = LapPhi(1)


GradientEnergy = - GradientEnergy

!
end function GradientEnergy



double precision function potential(phi,lp)
use solution_parameters
implicit none

double precision, intent (in)::phi(N_phases)
integer, intent(in):: lp
integer :: ip,jp

potential=0d0
!potential =  2d0*phi(1)*(1d0-phi(1))*(1d0-2d0*phi(1))
potential = 2*phi(1)-6*phi(1)**2+4*phi(1)**3



end function potential




double precision function PhaseRHS(vbles,dx,lp,Phi,i,j)
use solution_parameters
implicit none
double precision, intent(in)::dx,Phi(N_phases),vbles(N_unk+1,10,10)
integer, intent(in):: lp,i,j
double precision :: LapPhi(N_phases)
double precision GradientEnergy,potential
double precision :: phi_dot=0





 phaseRHS=- GradientEnergy(vbles,lp,dx,i,j)-potential(Phi,lp)/g_lambda**2 -7d-3

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
double precision vbles(N_unk+1,10,10),unk_(M_unk),MolPerVol
logical Not_number
integer ip,n,ii,jj,kk
  Fi=0.

!  write(*,*)"in get_RHS_Fi_2d"


  do ii=1,10
     do jj=1,10
        do ip=1,N_unk
!           vbles(ip,ii,jj)=unk(1+(ip-1)*nunkvbles,i+ii-2,j+jj-2,1,lb)
           vbles(ip,ii,jj)=unk(1+(ip-1)*nunkvbles,ii,jj,1,lb)
        enddo
     enddo
  enddo




  do ii=1,M_unk
     unk_(ii) = unk(ii,i,j,k,lb)
  enddo


  call get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_,i,j)



!  write(*,*)"exit get_RHS_Fi_2d"
 
end subroutine get_RHS_Fi_2d

subroutine get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_,i,j)

  use solution_parameters

  implicit none

!
integer, intent(in)::mode_1,nunkvbles,lnblocks,i,j
double precision, intent(in)::dx,dt,dtold,unk_(M_unk)
double precision, intent(out):: Fi(N_unk)
double precision, intent(in):: vbles(N_unk+1,10,10)
double precision :: c,T,phi(N_phases),c_dot,phi_dot(N_phases)
double precision rfactor,rf4,rf5,rf6,S,GradientEnergy,potential
double precision TemperatureRHS,SoluteRHS,PhaseRHS,T_k,phi_,FreeEnergy,MolPerVol
logical Not_number
integer ip,n

!  write(*,*)"in get_RHS_Fi"


Fi=0d0

  if (mode_1.eq.1) then
     phi_dot(1) = (vbles(1,i,j)-unk_(2))/dt

  elseif(mode_1.eq.2)then
     rfactor = dt/dtold

     rf4 = rfactor*rfactor/(rfactor+1.0)
     rf5 = (rfactor+1.0)
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
     
     phi_dot(1) = (rf4*unk_(5)-rf5*unk_(4)+rf6*vbles(1,i,j))/dt
  else
     write(*,*)"mode_1 not = 1 or 2"
     stop
  endif


  phi(1)=vbles(1,i,j)
  


  Fi(1) = PhaseRHS(vbles,dx,1,Phi,i,j)

  do ip=1,N_unk
     n=(ip-1)*nunkvbles
     if(mode_1.le.1)then
        Fi(ip)=vbles(ip,i,j)-dt*Fi(ip)-unk_(n+2)
     else if(mode_1.eq.2)then
        Fi(ip)=vbles(ip,i,j)-dt*Fi(ip)/rf6-unk_(n+2)
     end if
  enddo




!  write(*,*)"exit get_RHS_Fi"
end subroutine get_RHS_Fi

subroutine Calc_GradPhi(LapPhi,vbles,dx,i,j)
use solution_parameters
implicit none
integer, intent(in) :: i,j
double precision, intent(in)::dx,vbles(N_unk+1,10,10)
double precision, intent(out)::LapPhi(N_phases)
double precision, dimension(3,3)::stencil=reshape((/1.,4.,1.,4.,-20.,4.,1.,4.,1./),(/3,3/))
double precision, dimension(3,3):: stencil5=reshape((/0.,1.,0.,1.,-4.,1.,0.,1.,0./),(/3,3/))
integer ip,ii,jj

 write(*,*)"write this code"
stop
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


