! The 0.34 = W/(R*T_melt*molpervol)=7e7/(8.31*456*54185)

 D_c=0.

 do ip=1,N_phases
  D_c=D_c + 0.5* D1(ip)*( vbles(ip,ih,jv)*f_c(vbles(5,ih,jv),0)+vbles(ip,2,2)*f_c(c_c,0))
 enddo
 c=vbles(5,ih,jv)
 T=vbles(4,ih,jv)
 do ip=1,N_phases
  Phi(ip)=vbles(ip,ih,jv)
 enddo
      SoluteRHS = SoluteRHS + D_c*(FreeEnergy(c,T,Phi,c_dot,phi_dot,N_phases+2) + potential(phi,c,N_phases+2)) 
      

!  centre
 c=vbles(5,2,2)
 T=vbles(4,2,2)
 do ip=1,N_phases
  Phi(ip)=vbles(ip,2,2)
 enddo
      SoluteRHS = SoluteRHS - D_c*(FreeEnergy(c,T,Phi,c_dot,phi_dot,N_phases+2) + potential(phi,c,N_phases+2)) 
				   
