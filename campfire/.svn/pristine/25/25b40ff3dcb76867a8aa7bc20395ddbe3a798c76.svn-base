
  D_c=0.5* ( (vbles(1,ih,jv,kd)*g_D(1)+(1.-vbles(1,ih,jv,kd))*g_D(2))*f_c(vbles(N_phases+1,ih,jv,kd),0)&
       +(vbles(1,2,2,2)*g_D(1)+(1.-vbles(1,2,2,2))*g_D(2))*f_c(c_c,0))


  ! neighbours
  c=vbles(N_phases+1,ih,jv,kd)
  T=vbles(N_phases+2,ih,jv,kd)
  
  Phi=vbles(1,ih,jv,kd)

  SoluteRHS = SoluteRHS+D_c*FreeEnergy(c,T,Phi,c_dot,phi_dot,2)


  
  
  
  !  centre
  c=vbles(N_phases+1,2,2,2)
  T=vbles(N_phases+2,2,2,2)
  
  Phi=vbles(1,2,2,2)

  SoluteRHS = SoluteRHS - D_c*FreeEnergy(c,T,Phi,c_dot,phi_dot,2)





