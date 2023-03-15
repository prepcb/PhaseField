  write(*,*)
  stop
  phi_=(vbles(1,2,2)+vbles(1,ih,jv))*0.5
  D_c=0.5* ( (vbles(1,ih,jv)*g_D(1)+(1.-vbles(1,ih,jv))*g_D(2))*f_c(vbles(N_unk,ih,jv),0)&
       +(vbles(1,2,2)*g_D(1)+(1.-vbles(1,2,2))*g_D(2))*f_c(c_c,0))*3d1*phi_**2*(1d0-phi_)**2


  ! neighbours
  c=vbles(N_unk,ih,jv)
  T=vbles(N_unk-1,ih,jv)
  
  Phi(1)=vbles(1,ih,jv)

  SoluteRHS = SoluteRHS+D_c*FreeEnergy_ij(c,T,Phi,c_dot,phi_dot,3,1,2)


  
  
  
  !  centre
  c=vbles(N_unk,2,2)
  T=vbles(N_unk-1,2,2)
  
  Phi(1)=vbles(1,2,2)

  SoluteRHS = SoluteRHS - D_c*FreeEnergy_ij(c,T,Phi,c_dot,phi_dot,3,1,2)





