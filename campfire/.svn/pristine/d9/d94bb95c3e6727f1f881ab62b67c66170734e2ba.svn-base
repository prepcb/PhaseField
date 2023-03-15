  
!!$  D_c=0.5* ( (vbles(1,ih,jv)*g_D(1)+(1.-vbles(1,ih,jv))*g_D(2))*f_c(vbles(2,ih,jv),0)&
!!$       +(vbles(1,2,2)*g_D(1)+(1.-vbles(1,2,2))*g_D(2))*f_c(vbles(2,2,2),0))

  D_c=0.5* ( vbles(1,ih,jv)*g_D(1)*f_c(vbles(2,ih,jv),0)+(1.-vbles(1,ih,jv))*g_D(2)*f_c(AlFey1(vbles(2,ih,jv),0),0)&
       +vbles(1,2,2)*g_D(1)*f_c(vbles(2,2,2),0)+(1.-vbles(1,2,2))*g_D(2)*f_c(AlFey1(vbles(2,2,2),0),0))
  
  FE=0d0
  if(g_alpha.ne.0d0)then
     if((.not.g_quad).or.New_potential)then
        write(*,*)"anti trapping only works for g_Quad=.true.c_stencil.f90",g_quad,new_potential
        stop
     end if
     
     abs_grad_phi=sqrt((vbles(1,ih,jv)-vbles(1,2,2))**2&
          +0.0625*(   vbles(1,4-jv,ih)+vbles(1,2+ih-jv,ih+jv-2)&
          -vbles(1,jv,4-ih)-vbles(1,ih+jv-2,2+jv-ih)    )**2)/dx+1d-9 
     
     u =  g_phi_dot/abs_grad_phi
     
     FE = 0.5*(FreeEnergy(vbles(2,ih,jv),0.,vbles(1,ih,jv),-1,.true.)+FreeEnergy(vbles(2,2,2),0.,vbles(1,2,2),-1,.true.))
     
  endif

  if(g_xy_position(2).le.0d0)then
     write(*,*)"y position < 0"
     write(*,*)"c_stencil.f90"
     stop
  endif
  radial_factor=0d0
  if(g_axis_symmetric)radial_factor = (jv-2)*0.5*dx/g_xy_position(2)
  !  centre
  c=vbles(2,2,2)
  T=0d0
  Phi=vbles(1,2,2)
  SoluteRHS = SoluteRHS - D_c*FreeEnergy(c,T,Phi,2,.true.)*(1d0+radial_factor)
  if(g_quad.and.(.not. new_potential))then
     if(g_alpha.ge.0d0)then
        SoluteRHS = SoluteRHS - D_c*a_trap*u*FE*phi
     endif
  endif
  ! neighbours
  c=vbles(2,ih,jv)
  T=0d0
  
  phi = vbles(1,ih,jv)
  SoluteRHS = SoluteRHS + D_c*FreeEnergy(c,T,Phi,2,.true.)*(1d0+radial_factor)
  
  if(g_quad.and.(.not. new_potential))then
     if(g_alpha.ge.0d0)then
        SoluteRHS = SoluteRHS + D_c*a_trap*u*FE*phi
     endif
  endif
  
