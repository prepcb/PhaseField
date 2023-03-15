!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  Routines for phase field stencils in 3-d
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  History: James Green 2008-9, Chris Goodyer, 2009-2012
!!!!  Latest version: Peter Bollada, April 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!  Note that get_stencil_diff and get_stencil are now 3-d only

subroutine get_stencil_diff(i,j,k,dx,lb,m,dFi_dvi,LapMult)
  use multigrid_parameters
  use solution_parameters
  use time_dep_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use multigrid_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: i,j,k,lb,m
  double precision, intent(in) :: dx,LapMult
  double precision, intent(out) :: dFi_dvi
  double precision :: U,phi,c, divisor, A_func,alewis, phix, phiy, phiz, gradphisq
  
  double precision :: rfactor, rf4, rf5, rf6
  double precision :: Lap_mid_27=-128./30. !27 point
  double precision :: Lap_mid_19=-4. !19 point
  double precision :: Lap_mid_7=-1. !7 point
  double precision :: f_,g_,df_,dg_
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
  ! so get_dFi_dvi(i,j,k,dx,lb,2,1,dFi_dvi)
  !    corresponds to unk(2,i,j,k,lb)
  ! USER EDIT CODE
  U=unk(1+nunkvbles,i,j,k,lb)
  phi=unk(1,i,j,k,lb)
  if(total_vars.eq.3)then
     c=unk(1+2*nunkvbles,i,j,k,lb)
     alewis=1.0/le
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
        if((total_vars.lt.3))then
           ! For thermal only model mcinf=0, alewis=1;
           ! For solute only model alewis=0 and c is actually U_0 which boils term down to 1
           divisor=1.0
        else
           divisor=alewis+mcinf*(1.0+(1.0-ke)*c)
        end if
        call get_dunk_n_dx(1,i,j,k,dx,lb,phix,.false.)
        call get_dunk_n_dy(1,i,j,k,dx,lb,phiy,.false.)
        call get_dunk_n_dz(1,i,j,k,dx,lb,phiz,.false.)
        gradphisq = (phix*phix+phiy*phiy+phiz*phiz)
        call get_w_aniso_2(phix,phiy,phiz,gradphisq,A_func)
        dFi_dvi=LapMult*Lap_mid_27/(dx*dx*divisor)-(3.0*phi*phi-1.0+lambda*(U+mcinf*c)*4.0*phi*(phi*phi-1))/(A_func*A_func*divisor)
     else if(current_var.eq.2)then
        dFi_dvi=D_therm*Lap_mid_27/(dx*dx)
     else if(current_var.eq.3)then
        if((1+ke)/2-(1-ke)/2*unk(1,i,j,k,lb).eq.0)then 
           dFi_dvi=0.
        else if (mode_1.eq.1) then
           dFi_dvi=D_solute*(1.-unk(1,i,j,k,lb))/2.0*(Lap_mid_27/(dx*dx))/((1+ke)/2-(1-ke)/2*unk(1,i,j,k,lb))+&
                0.5*(1-ke)*(phi-unk(2,i,j,k,lb))/dt/((1.+ke)/2.-(1.-ke)/2.*phi)
        else           
           dFi_dvi=D_solute*(1-unk(1,i,j,k,lb))/2.0*(Lap_mid_27/(dx*dx))/((1+ke)/2-(1-ke)/2*unk(1,i,j,k,lb))+&
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
end subroutine get_stencil_diff

subroutine get_stencil(i,j,k,dx,lb,m,Fi,LapMult)
  use solution_parameters
  use time_dep_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  include 'mpif.h'
  integer, intent(in) :: i,j,k,lb,m
  double precision, intent(in) :: dx
  double precision, intent(out) :: Fi,LapMult
!  double precision :: sin_theta,cos_theta,sin_phi,cos_phi,phix,phiy,phiz
  double precision :: phix,phiy,phiz
!  double precision :: thetax,thetay,thetaz,phix,phiy,phiz,phixx,phixy,phixz,phiyy,phiyz,phizz
  double precision :: phixx,phixy,phixz,phiyy,phiyz,phizz
!  double precision :: A_func,dw_dw_dphi,t1,t2,t3,t4,t5,t6
  double precision :: A,A_func,t1,t2,t3, laplacian
  double precision :: dw_dphix,dw_dphiy,dw_dphiz,dw_dx,dw_dy,dw_dz,gradphisq,cx,cy,cz,cxx,cyy,czz,phi_dot_x,phi_dot_y,phi_dot_z
  double precision :: dh, phi, U, c, alewis, divisor
  double precision :: rfactor, rf4, rf5, rf6
  double precision :: A_func3d,Ai_func3d,Aij_func3d,Aii_func3d,xx,yy,zz
  double precision :: Ai(3),Aij(3,3),Gij(3,3),trace_Gij
  double precision, dimension(3,3,3)::stencil27 =reshape((/1., 3., 1., 3., 14., 3., 1., 3., 1., 3., 14., 3., 14., -128., 14., 3., 14., 3., 1., 3., 1., 3., 14., 3., 1., 3., 1./),(/3,3,3/))
  double precision, dimension(3,3,3)::stencil19 =reshape((/0., 1., 0., 1., 2., 1., 0., 1., 0., 1., 2., 1., 2., -24., 2., 1., 2., 1., 0., 1., 0., 1., 2., 1., 0., 1., 0./),(/3,3,3/))
  double precision, dimension(3,3,3)::stencil7 =reshape((/0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 1., -6., 1., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0./),(/3,3,3/))
  double precision :: Lap_27=30.,Lap_19=6.,Lap_7=1.
  integer ii,jj,kk
  
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
  ! so get_Fi(i,j,k,dx,lb,2,1,Fi)
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
        call get_dunk_n_dz(1,i,j,k,dx,lb,phiz,.false.)
        gradphisq=(phix*phix+phiy*phiy+phiz*phiz)
        call get_w_aniso_2(phix,phiy,phiz,gradphisq,A_func)
        if((total_vars.lt.3))then
           ! For thermal only model mcinf=0, alewis=1;
           ! For solute only model alewis=0 and c is actually U_0 which boils term down to 1
!           divisor=A_func*A_func
           divisor=1.0
        else
!           divisor=A_func*A_func*(alewis+mcinf*(1.0+(1.0-ke)*c))
           divisor=alewis+mcinf*(1.0+(1.0-ke)*c)
        end if
         
!   
!        27-point stencil - from page 64 of Jan's thesis  
!             CEG - still a bug as of 24/5/10, probably in amr_1blk_bcset.f90
        laplacian = 0.
        do ii=1,3
        do jj=1,3
        do kk=1,3
          laplacian = laplacian + stencil27(ii,jj,kk)*unk(1,i+ii-2,j+jj-2,k+kk-2,lb)
        enddo
        enddo
        enddo
        laplacian = laplacian/(Lap_27*dx*dx)
        Fi=-(phi*(phi*phi-1.0)+lambda*(U+mcinf*c)*(1-phi*phi)**2) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! new computation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       LapMult=1.0
        if(gradphisq.gt.0.0.and.epsilon_tilde.gt.0.0) then
           call get_d2unk_n_dx2(1,i,j,k,dx,lb,phixx,.false.)
           call get_d2unk_n_dxdy(1,i,j,k,dx,lb,phixy,.false.)
           call get_d2unk_n_dxdz(1,i,j,k,dx,lb,phixz,.false.)
           call get_d2unk_n_dy2(1,i,j,k,dx,lb,phiyy,.false.)
           call get_d2unk_n_dydz(1,i,j,k,dx,lb,phiyz,.false.)
           call get_d2unk_n_dz2(1,i,j,k,dx,lb,phizz,.false.)
 

           xx=phix*phix/gradphisq
           yy=phiy*phiy/gradphisq
           zz=phiz*phiz/gradphisq
           A = A_0*(1d0+epsilon_tilde*(xx**2+yy**2+zz**2))

           Gij(1,1) = A**2+4*(A-A_0)*(xx*(6*A-4*A_0)-A)-4*epsilon_tilde*A_0*xx*(xx*(12*A-8*A_0)-3*A)+16*epsilon_tilde**2*A_0**2*xx**3
           Gij(2,2) = A**2+4*(A-A_0)*(yy*(6*A-4*A_0)-A)-4*epsilon_tilde*A_0*yy*(yy*(12*A-8*A_0)-3*A)+16*epsilon_tilde**2*A_0**2*yy**3
           Gij(3,3) = A**2+4*(A-A_0)*(zz*(6*A-4*A_0)-A)-4*epsilon_tilde*A_0*zz*(zz*(12*A-8*A_0)-3*A)+16*epsilon_tilde**2*A_0**2*zz**3
           Gij(1,2) = 8*phix*phiy/gradphisq*(8*(3*A-2*A_0)*(A-A_0)-epsilon_tilde*(xx+yy)*(3*A-2*A_0)+2*epsilon_tilde**2*A_0**2)
           Gij(1,3) = 8*phix*phiz/gradphisq*(8*(3*A-2*A_0)*(A-A_0)-epsilon_tilde*(xx+zz)*(3*A-2*A_0)+2*epsilon_tilde**2*A_0**2)
           Gij(3,2) = 8*phiz*phiy/gradphisq*(8*(3*A-2*A_0)*(A-A_0)-epsilon_tilde*(zz+yy)*(3*A-2*A_0)+2*epsilon_tilde**2*A_0**2)

           trace_Gij=(Gij(1,1)+Gij(2,2)+Gij(3,3))/3.0 

           Gij(1,1)=Gij(1,1)-trace_Gij                
           Gij(2,2)=Gij(2,2)-trace_Gij                
           Gij(3,3)=Gij(3,3)-trace_Gij                
                                                      
           LapMult = LapMult+trace_Gij/(A_func*A_func)
                                                      
           Fi = Fi+trace_Gij*laplacian                
           Fi = Fi +2.0*(Gij(1,2)*phixy+Gij(3,2)*phiyz+Gij(1,3)*phixz)
           Fi = Fi +Gij(1,1)*phixx+Gij(2,2)*phiyy+Gij(3,3)*phizz


           
        endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! end new computation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 






!        LapMult=1.0
!        if(gradphisq.gt.0.0.and.epsilon_tilde.gt.0.0) then
!           call get_d2unk_n_dx2(1,i,j,k,dx,lb,phixx,.false.)
!           call get_d2unk_n_dxdy(1,i,j,k,dx,lb,phixy,.false.)
!           call get_d2unk_n_dxdz(1,i,j,k,dx,lb,phixz,.false.)
!           call get_d2unk_n_dy2(1,i,j,k,dx,lb,phiyy,.false.)
!           call get_d2unk_n_dydz(1,i,j,k,dx,lb,phiyz,.false.)
!           call get_d2unk_n_dz2(1,i,j,k,dx,lb,phizz,.false.)
           
!           xx=phix*phix/gradphisq
!           yy=phiy*phiy/gradphisq
!           zz=phiz*phiz/gradphisq
           
!           Ai(1)=Ai_func3d(phix,xx,yy,zz)
!           Ai(2)=Ai_func3d(phiy,yy,zz,xx)
!           Ai(3)=Ai_func3d(phiz,zz,xx,yy)
           
!           Aij(1,1)=Aii_func3d(xx,yy,zz)
!           Aij(2,2)=Aii_func3d(yy,zz,xx)
!           Aij(3,3)=Aii_func3d(zz,xx,yy)
           
!           Aij(1,2)=Aij_func3d(phix,phiy,xx,yy,zz)/gradphisq
!           Aij(2,3)=Aij_func3d(phiy,phiz,yy,zz,xx)/gradphisq
!           Aij(3,1)=Aij_func3d(phiz,phix,zz,xx,yy)/gradphisq
           
!           Gij(1,1)=(Ai(1)*Ai(1)+2.0*A_func*(Ai(1)*phix+Ai(1)*phix))/gradphisq+A_func*Aij(1,1)
!           Gij(2,2)=(Ai(2)*Ai(2)+2.0*A_func*(Ai(2)*phiy+Ai(2)*phiy))/gradphisq+A_func*Aij(2,2)
!           Gij(3,3)=(Ai(3)*Ai(3)+2.0*A_func*(Ai(3)*phiz+Ai(3)*phiz))/gradphisq+A_func*Aij(3,3)
!           Gij(1,2)=(Ai(1)*Ai(2)+2.0*A_func*(Ai(1)*phiy+Ai(2)*phix))/gradphisq+A_func*Aij(1,2)
!           Gij(2,3)=(Ai(2)*Ai(3)+2.0*A_func*(Ai(2)*phiz+Ai(3)*phiy))/gradphisq+A_func*Aij(2,3)
!           Gij(3,1)=(Ai(3)*Ai(1)+2.0*A_func*(Ai(3)*phix+Ai(1)*phiz))/gradphisq+A_func*Aij(3,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reproduce old version by removing these 6 lines        !
!           trace_Gij=(Gij(1,1)+Gij(2,2)+Gij(3,3))/3.0    !
                                                         !
!           Gij(1,1)=Gij(1,1)-trace_Gij                   !
!           Gij(2,2)=Gij(2,2)-trace_Gij                   !
!           Gij(3,3)=Gij(3,3)-trace_Gij                   !
                                                         !
!           LapMult = LapMult+trace_Gij/(A_func*A_func)   !
                                                         !
!           Fi = Fi+trace_Gij*laplacian                   !
                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!           Fi = Fi +2.0*(Gij(1,2)*phixy+Gij(2,3)*phiyz+Gij(3,1)*phixz)
!           Fi = Fi +Gij(1,1)*phixx+Gij(2,2)*phiyy+Gij(3,3)*phizz
           
!        end if
        Fi=Fi/(divisor*A_func*A_func)
        Fi = Fi + laplacian/divisor
     else if(current_var.eq.2)then
        laplacian = 0.
        do ii=1,3
        do jj=1,3
        do kk=1,3
          laplacian = laplacian + stencil27(ii,jj,kk) *unk(1+nunkvbles,i+ii-2,j+jj-2,k+kk-2,lb)
        enddo
        enddo
        enddo
        laplacian = laplacian/(Lap_27*dx*dx)
        
        call get_dh_dpsi(phi,dh)
        if (mode_1.eq.1) then
          Fi=D_therm*laplacian + 0.5*dh*(phi-unk(2,i,j,k,lb))
        else
          Fi=D_therm*laplacian + 0.5*dh*(rf4*unk(5,i,j,k,lb)-rf5*unk(4,i,j,k,lb)+rf6*phi)/dt
        endif
     else if(current_var.eq.3)then
        call get_dunk_n_dx(1,i,j,k,dx,lb,phix,.false.)
        call get_dunk_n_dy(1,i,j,k,dx,lb,phiy,.false.)
        call get_dunk_n_dz(1,i,j,k,dx,lb,phiz,.false.)
        call get_d2unk_n_dx2(1,i,j,k,dx,lb,phixx,.false.)
        call get_d2unk_n_dy2(1,i,j,k,dx,lb,phiyy,.false.)
        call get_d2unk_n_dz2(1,i,j,k,dx,lb,phizz,.false.)
        call get_d2unk_n_dxdy(1,i,j,k,dx,lb,phixy,.false.)
        call get_d2unk_n_dxdz(1,i,j,k,dx,lb,phixz,.false.)
        call get_d2unk_n_dydz(1,i,j,k,dx,lb,phiyz,.false.)
        call get_dunk_n_dx(current_unk,i,j,k,dx,lb,cx,.false.)
        call get_dunk_n_dy(current_unk,i,j,k,dx,lb,cy,.false.)
        call get_dunk_n_dz(current_unk,i,j,k,dx,lb,cz,.false.)
! INCORRECTLY USING BDF1 here to make 3rd antitrapping work
!        if (mode_1.eq.1) then
        if (mode_1.ge.1) then
           phi_dot_x=1.0/(2.0*dt*dx)*(unk(1,i+1,j,k,lb)-unk(1,i-1,j,k,lb)-(unk(4,i+1,j,k,lb)-unk(4,i-1,j,k,lb)))
           phi_dot_y=1.0/(2.0*dt*dx)*(unk(1,i,j+1,k,lb)-unk(1,i,j-1,k,lb)-(unk(4,i,j+1,k,lb)-unk(4,i,j-1,k,lb)))
           phi_dot_z=1.0/(2.0*dt*dx)*(unk(1,i,j,k+1,lb)-unk(1,i,j,k-1,lb)-(unk(4,i,j,k+1,lb)-unk(4,i,j,k-1,lb)))
        else 
           phi_dot_x=1.0/(2.0*dt*dx)*(  rf6*(unk(1,i+1,j,k,lb)-unk(1,i-1,j,k,lb)) &
                                      - rf5*(unk(4,i+1,j,k,lb)+unk(4,i-1,j,k,lb))&
                                      + rf4*(unk(5,i+1,j,k,lb)-unk(5,i-1,j,k,lb)))
           phi_dot_y=1.0/(2.0*dt*dx)*(  rf6*(unk(1,i,j+1,k,lb)-unk(1,i,j-1,k,lb)) &
                                      - rf5*(unk(4,i,j+1,k,lb)+unk(4,i,j-1,k,lb))&
                                      + rf4*(unk(5,i,j+1,k,lb)-unk(5,i,j-1,k,lb)))
           phi_dot_z=1.0/(2.0*dt*dx)*(  rf6*(unk(1,i,j,k+1,lb)-unk(1,i,j,k-1,lb)) &
                                      - rf5*(unk(4,i,j,k+1,lb)+unk(4,i,j,k-1,lb))&
                                      + rf4*(unk(5,i,j,k+1,lb)-unk(5,i,j,k-1,lb)))
        endif
        Fi=0
! 1st antitrapping term
        t1=phix*cx+phiy*cy+phiz*cz
        if(phix*phix+phiy*phiy+phiz*phiz.eq.0)then
           t2=0
        else
           if (mode_1.eq.1) then
!           if (mode_1.ge.1) then
              t2=-D_solute/2.0+1.0*anti_trapping_mod/(2.0*sqrt(2.0))*(1-ke)*&
                   (unk(1,i,j,k,lb)-unk(2,i,j,k,lb))/dt/sqrt(phix*phix+phiy*phiy+phiz*phiz)
           else
              t2=-D_solute/2.0+1.0*anti_trapping_mod/(2.0*sqrt(2.0))*(1-ke)*&
                   (rf6*unk(1,i,j,k,lb)- rf5*unk(4,i,j,k,lb)+ &
                    rf4*unk(5,i,j,k,lb))/dt/sqrt(phix*phix+phiy*phiy+phiz*phiz)
           endif
        end if
        Fi=t1*t2

! 2nd antitrapping term
!         t1=cxx+cyy+czz
        laplacian=0.
        do ii=1,3
        do jj=1,3
        do kk=1,3
          laplacian = laplacian + stencil27(ii,jj,kk)*unk(1+2*nunkvbles,i+ii-2,j+jj-2,k+kk-2,lb)
        enddo
        enddo
        enddo
        laplacian = laplacian/(Lap_27*dx*dx)
        
        t2=0.5*D_solute*(1.-unk(1,i,j,k,lb))
        Fi=Fi+laplacian*t2
        if(phix*phix+phiy*phiy+phiz*phiz.eq.0)then
           t1=0
        else
           t1=(phixx+phiyy+phizz)/sqrt(phix*phix+phiy*phiy+phiz*phiz)
           t1=t1-(phix*(phix*phixx+phiy*phixy+phiz*phixz)&
                +phiy*(phix*phixy+phiy*phiyy+phiz*phiyz)&
                +phiz*(phix*phixz+phiy*phiyz+phiz*phizz))/sqrt(phix*phix+phiy*phiy+phiz*phiz)**3
        end if
        if (mode_1.eq.1) then
!        if (mode_1.ge.1) then
           t2=(1.0+(1.0-ke)*unk(current_unk,i,j,k,lb))*(unk(1,i,j,k,lb)-unk(2,i,j,k,lb))/(dt*2.0*sqrt(2.0))
        else
           t2=(1.0+(1.0-ke)*unk(current_unk,i,j,k,lb))*(rf6*unk(1,i,j,k,lb) -        &
                    rf5*unk(4,i,j,k,lb)  +  rf4*unk(5,i,j,k,lb))/(dt*2.0*sqrt(2.0))
        endif
        Fi=Fi+anti_trapping_mod*t1*t2


! 3rd antitrapping term
        t1=phix*phi_dot_x+phiy*phi_dot_y+phiz*phi_dot_z
        if(phix*phix+phiy*phiy+phiz*phiz.eq.0)then
           t2=0
        else
           t2=1.0/(2.0*sqrt(2.0))*(1.0+(1.0-ke)*unk(current_unk,i,j,k,lb))/sqrt(phix*phix+phiy*phiy+phiz*phiz)
        end if
        Fi=Fi+anti_trapping_mod*t1*t2
        t3=(1+ke)/2-(1-ke)/2*phi
        if (mode_1.eq.1) then
!        if (mode_1.ge.1) then
           Fi=Fi+0.5*(1+(1-ke)*c)*(phi-unk(4,i,j,k,lb))/dt
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
end subroutine get_stencil

