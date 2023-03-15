

!I need to address this properly PCB
     s1 = vbles_(1,ih,jv)*g_D(2)+ vbles_(2,ih,jv)*g_D(1)
     s2 = vbles_(1,2,2)*g_D(2)+ vbles_(2,2,2)*g_D(1)


  D_c = 0.5*(s1*f_c(vbles(N_unk,ih,jv),0)+s2*f_c(c_c,0))
  ! neighbours
  c=vbles(N_unk,ih,jv)
  T=vbles(N_unk-1,ih,jv)
  do ip=1,3
     Phi(ip)=vbles_(ip,ih,jv)
  enddo
     SoluteRHS = SoluteRHS+D_c*FreeEnergyPhi(c,T,Phi,c_dot,phi_dot,N_unk)
  
  
  !  centre
  c=vbles(N_unk,2,2)
  T=vbles(N_unk-1,2,2)
  do ip=1,3
     Phi(ip)=vbles_(ip,2,2)
  enddo
     SoluteRHS = SoluteRHS-D_c*FreeEnergyPhi(c,T,Phi,c_dot,phi_dot,N_unk)
