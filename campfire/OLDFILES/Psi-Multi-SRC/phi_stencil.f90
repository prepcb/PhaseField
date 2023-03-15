


!  phi_=0.5*(vbles(1,ih,jv)+vbles(1,2,2))
!  beta = g_beta*phi_**2*(1.0-phi_)**2*16.0*g_lambda
  beta = g_beta*g_lambda


  D_c=0.5* ( (vbles(1,ih,jv)*g_D(1)+(1.-vbles(1,ih,jv))*g_D(2))*f_c(vbles(N_unk,ih,jv),0)&
       +(vbles(1,2,2)*g_D(1)+(1.-vbles(1,2,2))*g_D(2))*f_c(c_c,0))

!  D_c=0.5* ( g_D(1)*f_c(vbles(3,ih,jv),0)&
!       +g_D(1)*f_c(vbles(3,2,2),0))


  ! neighbours
  c=vbles(N_unk,ih,jv)
  T=vbles(N_unk-1,ih,jv)
    Phi(1)=vbles(1,ih,jv)

  PhiSoluteRHS = PhiSoluteRHS + D_c*beta*vbles(4,ih,jv)
  
  
  
  
  !  centre
  c=vbles(N_unk,2,2)
  T=vbles(N_unk-1,2,2)
    Phi(1)=vbles(1,2,2)

  PhiSoluteRHS = PhiSoluteRHS - D_c*beta*vbles(4,2,2)
  





