  phi_=(vbles(1,2,2)+vbles(1,ih,jv))*0.5
  D_c=0.5* ( (vbles(1,ih,jv)*g_D(1)+(1.-vbles(1,ih,jv))*g_D(2))*f_c(vbles(N_phases+2,ih,jv),0)&
       +(vbles(1,2,2)*g_D(1)+(1.-vbles(1,2,2))*g_D(2))*f_c(c_c,0))*3d1*phi_**2*(1d0-phi_)**2


  ! neighbours
  c=vbles(N_phases+2,ih,jv)
  T=vbles(N_phases+1,ih,jv)
  
  Phi(1)=vbles(1,ih,jv)

  SoluteRHS = SoluteRHS+D_c*FreeEnergy(c,T,Phi,c_dot,phi_dot,3)


  
  
  
  !  centre
  c=vbles(N_phases+2,2,2)
  T=vbles(N_phases+1,2,2)
  
  Phi(1)=vbles(1,2,2)

  SoluteRHS = SoluteRHS - D_c*FreeEnergy(c,T,Phi,c_dot,phi_dot,3)





