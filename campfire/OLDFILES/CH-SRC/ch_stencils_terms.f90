!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  Utilities for phase field stencils in 2-d and 3-d
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  History: James Green 2008-9, Chris Goodyer, 2009-2012
!!!!  Latest version: Peter Bollada, April 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_w_aniso_2(phix,phiy,phiz,grad_phi_sq,A_func)
  use solution_parameters
  implicit none
  double precision, intent(out) :: A_func
  double precision, intent(in) :: phix,phiy,phiz,grad_phi_sq
  if(epsilon_tilde.eq.0.0)then
     A_func=A_0
  else
     if(grad_phi_sq.gt.0)then
        A_func=A_0*(1+epsilon_tilde*(phix*phix*phix*phix+phiy*phiy*phiy*phiy+phiz*phiz*phiz*phiz)/(grad_phi_sq*grad_phi_sq))
     else
        A_func=A_0
     end if
  end if
  if (abs(A_func).lt.1.0e-14) A_func=0.0
end subroutine get_w_aniso_2

subroutine get_dw_dphix(phix,phiy,phiz,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phiz,grad_phi_sq

  val=4.0*A_0*epsilon_tilde*phix*(phix*phix*grad_phi_sq-phix*phix*phix*phix-phiy*phiy*phiy*phiy-phiz*phiz*phiz*phiz)/(grad_phi_sq*grad_phi_sq*grad_phi_sq)

end subroutine get_dw_dphix

subroutine get_dw_dphiy(phix,phiy,phiz,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phiz,grad_phi_sq

  val=4.0*A_0*epsilon_tilde*phiy*(phiy*phiy*grad_phi_sq-phix*phix*phix*phix-phiy*phiy*phiy*phiy-phiz*phiz*phiz*phiz)/(grad_phi_sq*grad_phi_sq*grad_phi_sq)
end subroutine get_dw_dphiy

subroutine get_dw_dphiz(phix,phiy,phiz,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phiz,grad_phi_sq
        val=4.0*A_0*epsilon_tilde*phiz*(phiz*phiz*grad_phi_sq-phix*phix*phix*phix-phiy*phiy*phiy*phiy-phiz*phiz*phiz*phiz)/(grad_phi_sq*grad_phi_sq*grad_phi_sq)
end subroutine get_dw_dphiz

subroutine get_dw_dx(phix,phiy,phiz,phixx,phixy,phixz,phiyy,phiyz,phizz,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phiz,phixx,phixy,phixz,phiyy,phiyz,phizz,grad_phi_sq
        val=4.0*A_0*epsilon_tilde*(grad_phi_sq*(phix*phix*phix*phixx+phiy*phiy*phiy*phixy+phiz*phiz*phiz*phixz)&
             -(phix*phix*phix*phix+phiy*phiy*phiy*phiy+phiz*phiz*phiz*phiz)*(phix*phixx+phiy*phixy+phiz*phixz))/(grad_phi_sq*grad_phi_sq*grad_phi_sq)
end subroutine get_dw_dx

subroutine get_dw_dy(phix,phiy,phiz,phixx,phixy,phixz,phiyy,phiyz,phizz,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phiz,phixx,phixy,phixz,phiyy,phiyz,phizz,grad_phi_sq
        val=4.0*A_0*epsilon_tilde*(grad_phi_sq*(phix*phix*phix*phixy+phiy*phiy*phiy*phiyy+phiz*phiz*phiz*phiyz)&
             -(phix*phix*phix*phix+phiy*phiy*phiy*phiy+phiz*phiz*phiz*phiz)*(phix*phixy+phiy*phiyy+phiz*phiyz))/(grad_phi_sq*grad_phi_sq*grad_phi_sq)
end subroutine get_dw_dy

subroutine get_dw_dz(phix,phiy,phiz,phixx,phixy,phixz,phiyy,phiyz,phizz,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phiz,phixx,phixy,phixz,phiyy,phiyz,phizz,grad_phi_sq
        val=4.0*A_0*epsilon_tilde*(grad_phi_sq*(phix*phix*phix*phixz+phiy*phiy*phiy*phiyz+phiz*phiz*phiz*phizz)&
             -(phix*phix*phix*phix+phiy*phiy*phiy*phiy+phiz*phiz*phiz*phiz)*(phix*phixz+phiy*phiyz+phiz*phizz))/(grad_phi_sq*grad_phi_sq*grad_phi_sq)
end subroutine get_dw_dz

subroutine get_dw_dpsi(sin_psi,cos_psi,sin_theta,cos_theta,dw_dpsi_val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: dw_dpsi_val
  double precision, intent(in) :: sin_psi,cos_psi,sin_theta,cos_theta
  if(epsilon_tilde.eq.0.0)then
     dw_dpsi_val=0.0
  else
     dw_dpsi_val=4.0*A_0*epsilon_tilde*cos_psi*sin_psi*(sin_psi*sin_psi*(1.0-2.0*sin_theta*sin_theta*cos_theta*cos_theta)-cos_psi*cos_psi)
  end if
end subroutine get_dw_dpsi

subroutine get_dunk_n_dx(n,i,j,k,dx,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dx
  logical, intent(in) :: workf
  double precision, intent(out) :: val
  if(workf.eqv..true.)then
     val = (work(i+1,j,k,lb,n)-work(i-1,j,k,lb,n))/(2.0*dx)
  else
     val = (unk(n,i+1,j,k,lb)-unk(n,i-1,j,k,lb))/(2.0*dx)
  end if
  if (abs(val).lt.1.0e-14) val=0.0
end subroutine get_dunk_n_dx

subroutine get_dunk_n_dy(n,i,j,k,dy,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dy
  logical, intent(in) :: workf
  double precision, intent(out) :: val
  if(workf)then
     val = (work(i,j+1,k,lb,n)-work(i,j-1,k,lb,n))/(2.0*dy) 
  else
     val = (unk(n,i,j+1,k,lb)-unk(n,i,j-1,k,lb))/(2.0*dy)
  end if
  if (abs(val).lt.1.0e-14) val=0.0
end subroutine get_dunk_n_dy

subroutine get_dunk_n_dz(n,i,j,k,dz,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dz
  logical, intent(in) :: workf
  double precision, intent(out) :: val
  if(workf)then
     val = (work(i,j,k+1,lb,n)-work(i,j,k-1,lb,n))/(2.0*dz) 
  else
     val = (unk(n,i,j,k+1,lb)-unk(n,i,j,k-1,lb))/(2.0*dz)
  end if
end subroutine get_dunk_n_dz

subroutine get_d2unk_n_dx2(n,i,j,k,dx,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dx
  logical, intent(in) :: workf
  double precision, intent(out) :: val
  if(workf)then
     val = (work(i+1,j,k,lb,n)-2.0*work(i,j,k,lb,n)+work(i-1,j,k,lb,n))/(dx*dx) 
  else
     val = (unk(n,i+1,j,k,lb)-2.0*unk(n,i,j,k,lb)+unk(n,i-1,j,k,lb))/(dx*dx)
  end if
  if (abs(val).lt.1.0e-14) val=0.0
end subroutine get_d2unk_n_dx2

subroutine get_d2unk_n_dy2(n,i,j,k,dy,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dy
  logical, intent(in) :: workf
  double precision, intent(out) :: val
  if(workf)then
     val = (work(i,j+1,k,lb,n)-2.0*work(i,j,k,lb,n)+work(i,j-1,k,lb,n))/(dy*dy)
  else
     val = (unk(n,i,j+1,k,lb)-2.0*unk(n,i,j,k,lb)+unk(n,i,j-1,k,lb))/(dy*dy)
  end if
  if (abs(val).lt.1.0e-14) val=0.0
end subroutine get_d2unk_n_dy2

subroutine get_d2unk_n_dz2(n,i,j,k,dz,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dz
  logical, intent(in) :: workf
  double precision, intent(out) :: val
  if(workf)then
     val = (work(i,j,k+1,lb,n)-2.0*work(i,j,k,lb,n)+work(i,j,k-1,lb,n))/(dz*dz)
  else
     val = (unk(n,i,j,k+1,lb)-2.0*unk(n,i,j,k,lb)+unk(n,i,j,k-1,lb))/(dz*dz)
  end if
end subroutine get_d2unk_n_dz2

subroutine get_d2unk_n_dxdy(n,i,j,k,dx,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dx
  double precision, intent(out) :: val
  logical, intent(in) :: workf
  double precision :: temp1,temp2
  call get_dunk_n_dx(n,i,j+1,k,dx,lb,temp1,workf)
  call get_dunk_n_dx(n,i,j-1,k,dx,lb,temp2,workf)
  val = (temp1-temp2)/(2.0*dx)
  if (abs(val).lt.1.0e-14) val=0.0
end subroutine get_d2unk_n_dxdy

subroutine get_d2unk_n_dxdz(n,i,j,k,dx,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dx
  double precision, intent(out) :: val
  logical, intent(in) :: workf
  double precision :: temp1,temp2
  call get_dunk_n_dx(n,i,j,k+1,dx,lb,temp1,workf)
  call get_dunk_n_dx(n,i,j,k-1,dx,lb,temp2,workf)
  val = (temp1-temp2)/(2.0*dx)
end subroutine get_d2unk_n_dxdz

subroutine get_d2unk_n_dydz(n,i,j,k,dy,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dy
  double precision, intent(out) :: val
  logical, intent(in) :: workf
  double precision :: temp1,temp2
  call get_dunk_n_dy(n,i,j,k+1,dy,lb,temp1,workf)
  call get_dunk_n_dy(n,i,j,k-1,dy,lb,temp2,workf)
  val = (temp1-temp2)/(2.0*dy)
end subroutine get_d2unk_n_dydz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Alternative versions for 2-d case                                  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_w_aniso_2_2d(phix,phiy,grad_phi_sq,A_func)
  use solution_parameters
  implicit none
  double precision, intent(out) :: A_func
  double precision, intent(in) :: phix,phiy,grad_phi_sq
  if(epsilon_tilde.eq.0.0)then
     A_func=A_0
  else
     if(grad_phi_sq.gt.0)then
        A_func=A_0*(1+epsilon_tilde*(phix*phix*phix*phix+phiy*phiy*phiy*phiy)/(grad_phi_sq*grad_phi_sq))
     else
        A_func=A_0
     end if
  end if
end subroutine get_w_aniso_2_2d

subroutine get_dw_dphix_2d(phix,phiy,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,grad_phi_sq
  val=4.0*A_0*epsilon_tilde*phix*(phix*phix*grad_phi_sq-phix*phix*phix*phix-phiy*phiy*phiy*phiy)/(grad_phi_sq*grad_phi_sq*grad_phi_sq)
end subroutine get_dw_dphix_2d

subroutine get_dw_dphiy_2d(phix,phiy,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,grad_phi_sq
  val=4.0*A_0*epsilon_tilde*phiy*(phiy*phiy*grad_phi_sq-phix*phix*phix*phix-phiy*phiy*phiy*phiy)/(grad_phi_sq*grad_phi_sq*grad_phi_sq)
end subroutine get_dw_dphiy_2d

subroutine get_dw_dx_2d(phix,phiy,phixx,phixy,phiyy,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phixx,phixy,phiyy, grad_phi_sq
  val=4.0*A_0*epsilon_tilde*(grad_phi_sq*(phix*phix*phix*phixx + phiy*phiy*phiy*phixy)&
        -(phix*phix*phix*phix+phiy*phiy*phiy*phiy)*(phix*phixx+phiy*phixy))/(grad_phi_sq*grad_phi_sq*grad_phi_sq)
end subroutine get_dw_dx_2d

subroutine get_dw_dy_2d(phix,phiy,phixx,phixy,phiyy,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phixx,phixy,phiyy,grad_phi_sq
  val=4.0*A_0*epsilon_tilde*(grad_phi_sq*(phix*phix*phix*phixy+phiy*phiy*phiy*phiyy)&
             -(phix*phix*phix*phix+phiy*phiy*phiy*phiy)*(phix*phixy+phiy*phiyy))/(grad_phi_sq*grad_phi_sq*grad_phi_sq)
end subroutine get_dw_dy_2d

double precision function Lap2d(dx,d_flag,i,j,k,lb,n)
  use solution_parameters
  use time_dep_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  include 'mpif.h'
  integer, intent (in) :: i,j,k,lb,n  
  double precision, intent (in):: dx
  logical, intent(in) :: d_flag
  double precision Lap_normal,Lap_diag,dLap,alpha 
  !alpha = 0 5_point,alpha=2/3 = 9_point, alpha=1 diagonal Lap  
  alpha = 1./3.
  if(d_flag) then  
     Lap2d = (2.0*alpha-4.0)/(dx*dx) 
  else  
     Lap_normal =   unk(1+n,i+1,j,  k,lb)&
               +unk(1+n,i  ,j+1,k,lb)&
               +unk(1+n,i-1,j  ,k,lb)&
               +unk(1+n,i  ,j-1,k,lb)&
               -4.0* unk(1+n,i ,j,k,lb)
     Lap_diag =     unk(1+n,i+1,j+1,k,lb)&
               +unk(1+n,i-1,j+1,k,lb)&
               +unk(1+n,i-1,j-1,k,lb)&
               +unk(1+n,i+1,j-1,k,lb)&
               -4.0* unk(1+n,i ,j,k,lb)
     Lap_diag = Lap_diag*0.5

     Lap2d = ((1.0-alpha)*Lap_normal+alpha*Lap_diag)/(dx*dx)
  end if
end function Lap2d

double precision function A_func2d(X,Y)
  use solution_parameters
  implicit none
  double precision, intent(in) :: X,Y
  double precision :: p,q

  if(epsilon_tilde.eq.0.0)then
     A_func2d=A_0
  else
     if(X+Y.gt.0)then
       p=X/(X+Y)
       q=Y/(X+Y)
       A_func2d=A_0*(1.0+epsilon_tilde*(1.0-2.*p*q))
     else
        A_func2d=A_0
     end if
  end if
  if (abs(A_func2d).lt.1.0e-14) then
     write(*,*)"In A_func2d"
     stop
  endif
end function A_func2d

double precision function A_func3d(X,Y,Z)
  use solution_parameters
  implicit none
  double precision, intent(in) :: X,Y,Z

  if(epsilon_tilde.eq.0.0)then
     A_func3d=A_0
  else
     A_func3d=A_0*(1.0+epsilon_tilde*(X*X+Y*Y+Z*Z))
  end if
end function A_func3d

double precision function Ai_func3d(phix,X,Y,Z)
  use solution_parameters
  implicit none
  double precision, intent(in) :: phix,X,Y,Z

  if(epsilon_tilde.eq.0.0)then
     Ai_func3d=0.0
  else
     Ai_func3d=4.0*A_0*epsilon_tilde*phix*(X*(Y+Z)-Y*Y-Z*Z)
  end if
end function Ai_func3d

double precision function Aij_func3d(phix,phiy,X,Y,Z)
  use solution_parameters
  implicit none
  double precision, intent(in) :: phix,phiy,X,Y,Z

  if(epsilon_tilde.eq.0.0)then
     Aij_func3d=0.0
  else
     Aij_func3d=8.0*A_0*epsilon_tilde*phix*phiy*(X*X+Y*Y-4.0*X*Y-2.0*Z*(X+Y)+3.0*Z*Z)
  end if
end function Aij_func3d

double precision function Aii_func3d(X,Y,Z)
  use solution_parameters
  implicit none
  double precision, intent(in) :: X,Y,Z

  if(epsilon_tilde.eq.0.0)then
     Aii_func3d=0.0
  else
     Aii_func3d=-4.0*A_0*epsilon_tilde*(3.0*X*X*(Y+Z)-8.0*X*(Y*Y+Z*Z)-6.0*X*Y*Z+Y*Y*Y+Y*Y*Z+Z*Z*Y+Z*Z*Z)
  end if
end function Aii_func3d
