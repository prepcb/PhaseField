

  s1 = 0d0
  s2 = 0d0
  do ip=1,n_phases
     s1 = s1 + vbles(ip,ih,jv)*g_D(ip)
     s2 = s2 + vbles(ip,2,2)*g_D(ip)
  enddo
  D_c = 0.5*(s1*f_c(vbles(N_unk,ih,jv),0)+s2*f_c(c_c,0))
  ! neighbours
  c=vbles(N_unk,ih,jv)
  T=vbles(N_unk-1,ih,jv)
  do ip=1,n_phases
     Phi(ip)=vbles(ip,ih,jv)
  enddo
     SoluteRHS = SoluteRHS+D_c*FreeEnergy(c,T,Phi,c_dot,phi_dot,N_unk)
 
  
  !  centre
  c=vbles(N_unk,2,2)
  T=vbles(N_unk-1,2,2)
  do ip=1,n_phases
     Phi(ip)=vbles(ip,2,2)
  enddo
     SoluteRHS = SoluteRHS-D_c*FreeEnergy(c,T,Phi,c_dot,phi_dot,N_unk)

