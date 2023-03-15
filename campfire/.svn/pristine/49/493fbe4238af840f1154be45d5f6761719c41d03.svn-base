!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  Routines for phase field stencils in 2-d
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  History: Chris Goodyer, 2011-2012
!!!!  Latest version: Peter Bollada, April 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_stencil_diff_2d(i,j,dx,lb,m,dFi_dvi,LapMult,tau_val)
  use multigrid_parameters
  use solution_parameters
  use time_dep_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  include 'mpif.h'
  integer, intent(in) :: i,j,lb,m
  double precision, intent(in) :: dx,LapMult,tau_val
  double precision, intent(out) :: dFi_dvi
  double precision :: U,phi,c, divisor,  alewis, phix, phiy, gradphisq
  double precision :: rfactor, rf4, rf5, rf6,Lap2d
  integer k
  k=1

  if (mode_1.eq.2) then 
     rfactor = dt/dtold
     rf4 = rfactor*rfactor/(rfactor+1.0)/dt
     rf5 = (rfactor+1.0)/dt
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)/dt
  else
     rf4=0.0
     rf5=1.0/dt
     rf6=1.0/dt
  endif

  ! dFi_dvi=d(N(v))/dvi
  ! N(v) is "whatever equals f" for node i,j,k
  ! Ie solving N(v)=f. N() can be non-linear
  ! d(N(v))/dvi is the result of differentiating N(v) wrt the node
  !    currently being solved (not the same as d(N(v))/dv)
  ! Eg if N(v)=v(x+1)-2v(x)+v(x-1) (1D d2v/dx2)
  ! d(N(v))/dvi=-2
  ! n is the variable number
  ! m is the variable type (cell/face/edge/vertex centered)
  ! m=1 -> cell centered
  !   2 -> vertex
  !   3 -> face
  !   4 -> edge
  ! so get_dFi_dvi(i,j,k,dx,dy,lb,2,1,dFi_dvi)
  !    corresponds to unk(2,i,j,k,lb)
  ! USER EDIT CODE
  U=unk(1+nunkvbles,i,j,k,lb)
  phi=unk(1,i,j,k,lb)
  if(total_vars.eq.3)then
     c=unk(1+2*nunkvbles,i,j,k,lb)
     alewis = 1.0/le    
  ! Handler code for when running phase/thermal or phase/solute only (ie not fully coupled)
  else
     if(solute)then
        U=delta
        c=unk(1+nunkvbles,i,j,k,lb)
        if(current_var.eq.2)then
           current_var=3
        end if
     else
        c=0.0
     end if
  end if
  if(m.eq.1)then
     if(current_var.eq.1)then
       dFi_dvi=LapMult*Lap2d(dx,.true.)-(3.0*phi*phi-1.0+lambda*(U+mcinf*c)*4.0*phi*(phi*phi-1))/tau_val
     else if(current_var.eq.2)then
        dFi_dvi=D_therm*Lap2d(dx,.true.)
     else if(current_var.eq.3)then
        if((1+ke)/2-(1-ke)/2*unk(1,i,j,k,lb).eq.0)then 
           dFi_dvi=0.
        else if (mode_1.eq.1) then
           dFi_dvi=D_solute*(1-unk(1,i,j,k,lb))/2.0*(-2.0/(dx*dx)-2.0/(dx*dx))/((1+ke)/2-(1-ke)/2*unk(1,i,j,k,lb))+&
                0.5*(1-ke)*(phi-unk(2,i,j,k,lb))/dt/((1+ke)/2-(1-ke)/2*phi)
        else           
           dFi_dvi=D_solute*(1-unk(1,i,j,k,lb))/2.0*(-2.0/(dx*dx)-2.0/(dx*dx))/((1+ke)/2-(1-ke)/2*unk(1,i,j,k,lb))+&
                0.5*(1-ke)*(rf4*unk(5,i,j,k,lb)-rf5*unk(4,i,j,k,lb)+rf6*phi)/((1+ke)/2-(1-ke)/2*phi)
        end if
     else
        print *,"Unhandled m/cvar combination in subroutine get_dFi_dvi m=",m,"n=",current_var
     end if

  else
     print *,"Unhandled m/n combination in subroutine get_dFi_dvi m=",m,"n=",current_var
  end if
  if((total_vars.lt.3).and.(solute).and.(current_var.eq.3))then
     current_var=2
  end if
  ! END USER EDIT
end subroutine get_stencil_diff_2d

subroutine get_stencil_2d(i,j,dx,lb,m,Fi, dx2inv,LapMult,tau_val)
  use solution_parameters
  use time_dep_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  include 'mpif.h'
  integer, intent(in) :: i,j,lb,m
  double precision, intent(in) :: dx, dx2inv
  double precision, intent(out) :: Fi,LapMult,tau_val
  double precision :: phix,phiy
  double precision :: phixx,phixy,phiyy
  double precision :: A_func,t1,t2,t3
  double precision :: dw_dphix,dw_dphiy,dw_dx,dw_dy,gradphisq,cx,cy,cxx,cyy,phi_dot_x,phi_dot_y
  double precision :: dh,phi,U,c, divisor, laplacian
  double precision :: rfactor, rf4, rf5, rf6
  double precision :: alewis
  double precision :: A_func2d
  double precision :: M11,M22,M12,A1,A2,A11,A22,A12, xx,yy,Lap2d,Lap_val
  integer k
  k=1

  if (mode_1.eq.2) then 
     rfactor = dt/dtold
     rf4 = rfactor*rfactor/(rfactor+1.0)
     rf5 = (rfactor+1.0)
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
  else
     rf4=0.0
     rf5=1.0
     rf6=1.0
  endif

  ! Fi=N(v)
  ! f is the right hand side for node i,j,k
  ! N(v) is "whatever equals f" for node i,j,k
  ! Ie solving N(v)=f. N() can be non-linear
  ! n is the variable number
  ! m is the variable type (cell/face/edge/vertex centered)
  ! m=1 -> cell centered
  !   2 -> vertex
  !   3 -> face
  !   4 -> edge
  ! so get_Fi(i,j,k,dx,dy,lb,2,1,Fi)
  !    corresponds to unk(2,i,j,k,lb)
  ! Ignore above, had such trouble getting mg working for non cell centered node I've not
  ! coded any of the others.
  ! USER EDIT CODE
  U=unk(1+nunkvbles,i,j,k,lb)
  phi=unk(1,i,j,k,lb)
  if(total_vars.eq.3)then
     c=unk(1+2*nunkvbles,i,j,k,lb)
     alewis = 1.0/le    
  else
  ! Handler code for when running phase/thermal or phase/solute only (ie not fully coupled)
     if(solute)then
        alewis = 0.0
        U = delta
        if(current_var.eq.2)then
           current_var=3
        end if
        c=unk(1+nunkvbles,i,j,k,lb)
     else
        alewis = 1.0/le    
        c=0.0
     end if
  end if
  if(m.eq.1)then
    if(current_var.eq.1)then
      call get_dunk_n_dx(1,i,j,k,dx,lb,phix,.false.)
      call get_dunk_n_dy(1,i,j,k,dx,lb,phiy,.false.)
      gradphisq=phix*phix+phiy*phiy       
      call get_w_aniso_2_2d(phix,phiy,gradphisq,A_func)      
      if((total_vars.lt.3))then
       divisor=1.0
      else
       divisor=(alewis+mcinf*(1.0+(1.0-ke)*c))
      end if
      tau_val=A_func*A_func*divisor
      LapMult=1.0/divisor
      Lap_val=Lap2d(dx,.false.,i,j,k,lb,0)
      Fi=-(phi*(phi*phi-1.0)+lambda*(U+mcinf*c)*(1.0-phi*phi)**2)
! see Docs/AnisotropyPCB.tex for explanation of following      
      if(gradphisq.gt.0.0.and.epsilon_tilde.gt.0.0)then 
	call get_d2unk_n_dx2(1,i,j,k,dx,lb,phixx,.false.)
	call get_d2unk_n_dy2(1,i,j,k,dx,lb,phiyy,.false.)
	call get_d2unk_n_dxdy(1,i,j,k,dx,lb,phixy,.false.)
	xx=phix*phix/gradphisq
	yy=phiy*phiy/gradphisq
	A1=4.*A_0*epsilon_tilde*phix*yy*(xx-yy)       
	A2=4.*A_0*epsilon_tilde*phiy*xx*(yy-xx)
	A11=-4.*A_0*epsilon_tilde*yy*(3.*xx*xx-8.*xx*yy+yy*yy)
	A22=-4.*A_0*epsilon_tilde*xx*(xx*xx-8.*xx*yy+3.*yy*yy)
	A12=8.*A_0*epsilon_tilde*phix*phiy*(xx*xx-4.*xx*yy+yy*yy)/gradphisq
	M11=(A1*A1+4.*A_func*A1*phix)/gradphisq+A_func*A11
	M22=(A2*A2+4.*A_func*A2*phiy)/gradphisq+A_func*A22
	M12=(A1*A2+2.*A_func*(A1*phiy+A2*phix))/gradphisq+A_func*A12
 	LapMult=LapMult+(M11+M22)*0.5/tau_val	
        Fi=Fi+(M11+M22)*0.5*Lap_val
	Fi=Fi+(M11-M22)*0.5*(phixx-phiyy)+2.*M12*phixy
      endif
      Fi=Fi/tau_val+Lap_val/divisor
    else if(current_var.eq.2)then
      call get_dh_dpsi(phi,dh) 
      if (mode_1.eq.1) then
       Fi=D_therm*Lap2d(dx,.false.,i,j,k,lb,nunkvbles)  +0.5*dh*(phi-unk(2,i,j,k,lb))/dt
      else
       Fi=D_therm*Lap2d(dx,.false.,i,j,k,lb,nunkvbles) +0.5*dh*(rf4*unk(5,i,j,k,lb)-rf5*unk(4,i,j,k,lb)+rf6*phi)/dt
      endif
    else if(current_var.eq.3)then
        call get_dunk_n_dx(1,i,j,k,dx,lb,phix,.false.)
        call get_dunk_n_dy(1,i,j,k,dx,lb,phiy,.false.)
        call get_d2unk_n_dx2(1,i,j,k,dx,lb,phixx,.false.)
        call get_d2unk_n_dy2(1,i,j,k,dx,lb,phiyy,.false.)
        call get_d2unk_n_dxdy(1,i,j,k,dx,lb,phixy,.false.)
        call get_dunk_n_dx(current_unk,i,j,k,dx,lb,cx,.false.)
        call get_dunk_n_dy(current_unk,i,j,k,dx,lb,cy,.false.)
        call get_d2unk_n_dx2(current_unk,i,j,k,dx,lb,cxx,.false.)
        call get_d2unk_n_dy2(current_unk,i,j,k,dx,lb,cyy,.false.)
        if (mode_1.ge.1) then
           phi_dot_x=1.0/(2.0*dt*dx)*(unk(1,i+1,j,k,lb)-unk(1,i-1,j,k,lb)-(unk(4,i+1,j,k,lb)-unk(4,i-1,j,k,lb)))
           phi_dot_y=1.0/(2.0*dt*dx)*(unk(1,i,j+1,k,lb)-unk(1,i,j-1,k,lb)-(unk(4,i,j+1,k,lb)-unk(4,i,j-1,k,lb)))
        else 
           phi_dot_x=1.0/(2.0*dt*dx)*(  rf6*(unk(1,i+1,j,k,lb)-unk(1,i-1,j,k,lb)) &
                                      - rf5*(unk(4,i+1,j,k,lb)+unk(4,i-1,j,k,lb))&
                                      + rf4*(unk(5,i+1,j,k,lb)-unk(5,i-1,j,k,lb)))
           phi_dot_y=1.0/(2.0*dt*dx)*(  rf6*(unk(1,i,j+1,k,lb)-unk(1,i,j-1,k,lb)) &
                                      - rf5*(unk(4,i,j+1,k,lb)+unk(4,i,j-1,k,lb))&
                                      + rf4*(unk(5,i,j+1,k,lb)-unk(5,i,j-1,k,lb)))
        endif
        Fi=0
! 1st antitrapping term
        t1=phix*cx+phiy*cy
        if(phix*phix+phiy*phiy.lt.1.0e-15)then
           t2=0
        else
           if (mode_1.eq.1) then
!           if (mode_1.ge.1) then
              t2=-D_solute/2.0+1.0*anti_trapping_mod/(2.0*sqrt(2.0))*(1-ke)*&
                   (unk(1,i,j,k,lb)-unk(2,i,j,k,lb))/dt/sqrt(phix*phix+phiy*phiy)
           else
              t2=-D_solute/2.0+1.0*anti_trapping_mod/(2.0*sqrt(2.0))*(1-ke)*&
                   (rf6*unk(1,i,j,k,lb)- rf5*unk(4,i,j,k,lb)+ &
                    rf4*unk(5,i,j,k,lb))/dt/sqrt(phix*phix+phiy*phiy)
           endif
        end if
        Fi=t1*t2

! 2nd antitrapping term
        t1=cxx+cyy
        if (ABS(t1).lt.1.0e-15) t1=0.0
        t2=D_solute/2*(1-unk(1,i,j,k,lb))
        Fi=Fi+t1*t2
!        print *,i,j,lb,Fi

        if(phix*phix+phiy*phiy.lt.1.0e-15)then
           t1=0
        else
           t1=(phixx+phiyy)/sqrt(phix*phix+phiy*phiy)
           t1=t1-(phix*(phix*phixx+phiy*phixy)&
                +phiy*(phix*phixy+phiy*phiyy)&
                )/sqrt(phix*phix+phiy*phiy)**3
        end if
        if (mode_1.eq.1) then
!        if (mode_1.ge.1) then
           t2=(1.0+(1.0-ke)*unk(current_unk,i,j,k,lb))*(unk(1,i,j,k,lb)-unk(2,i,j,k,lb))/(dt*2.0*sqrt(2.0))
        else
           t2=(1.0+(1.0-ke)*unk(current_unk,i,j,k,lb))*(rf6*unk(1,i,j,k,lb) -        &
                    rf5*unk(4,i,j,k,lb)  +  rf4*unk(5,i,j,k,lb))/(dt*2.0*sqrt(2.0))
!           t2=(1.0+(1.0-ke)*unk(current_unk,i,j,k,lb))*(unk(1,i,j,k,lb)-unk(2,i,j,k,lb))/(dt*2.0*sqrt(2.0))
!           if (rf6*unk(1,i,j,k,lb) -rf5*unk(4,i,j,k,lb)  +  rf4*unk(5,i,j,k,lb)- (unk(1,i,j,k,lb)-unk(2,i,j,k,lb)).gt.1.0e-12)&
!           print *, mype, rf6*unk(1,i,j,k,lb) -rf5*unk(4,i,j,k,lb)  +  rf4*unk(5,i,j,k,lb), unk(1,i,j,k,lb)-unk(2,i,j,k,lb),&
!                 rf6*unk(1,i,j,k,lb) -rf5*unk(4,i,j,k,lb)  +  rf4*unk(5,i,j,k,lb)- (unk(1,i,j,k,lb)-unk(2,i,j,k,lb))
        endif
        Fi=Fi+anti_trapping_mod*t1*t2


! 3rd antitrapping term
        t1=phix*phi_dot_x+phiy*phi_dot_y
        if(phix*phix+phiy*phiy.lt.1.0e-15)then
           t2=0
        else
           t2=1.0/(2.0*sqrt(2.0))*(1.0+(1.0-ke)*unk(current_unk,i,j,k,lb))/sqrt(phix*phix+phiy*phiy)
        end if
        Fi=Fi+anti_trapping_mod*t1*t2
        t3=(1+ke)/2-(1-ke)/2*phi
        if (mode_1.eq.1) then
!        if (mode_1.ge.1) then
           Fi=Fi+0.5*(1+(1-ke)*c)*(phi-unk(2,i,j,k,lb))/dt
        else
           Fi=Fi+0.5*(1+(1-ke)*c)*(rf4*unk(5,i,j,k,lb)-rf5*unk(4,i,j,k,lb)+rf6*phi)/dt
        endif
        if(t3.eq.0.0)then
           Fi=0.0
        else
           Fi=Fi/t3
        end if
     else
        print *,"Unhandled m/c_var combination in subroutine get_Fi m=",m,"cvar=",current_var
     end if   
   
   

  else
     print *,"Unhandled m/n combination in subroutine get_Fi m=",m,"n=",current_var
  end if
  if((total_vars.lt.3).and.(solute).and.(current_var.eq.3))then
     current_var=2
  end if
  ! END USER EDIT
end subroutine get_stencil_2d



