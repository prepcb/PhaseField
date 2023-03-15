
 D=0.
 D_c=0.
 D_cc=0.

 do ip=1,N_phases
  D= D + 0.5 * D1(ip)*( (1.+vbles(N_phases+1,ih,jv))*vbles(ip,ih,jv)+(1.+vbles(N_phases+1,2,2))*vbles(ip,2,2))
  D_c=D_c + 0.5* D1(ip)*( vbles(ip,ih,jv)*f_c(vbles(N_phases+2,ih,jv),0)+vbles(ip,2,2)*f_c(c_c,0))
  D_cc=D_cc + 0.5* D1(ip)*( vbles(ip,ih,jv)*f_cc(vbles(N_phases+2,ih,jv),0)+vbles(ip,2,2)*f_cc(c_c,0))
 enddo
 c=vbles(N_phases+2,ih,jv)
 T=vbles(N_phases+1,ih,jv)

  Phi(1)=vbles(1,ih,jv)

      SoluteRHS = SoluteRHS+D_c*(FreeEnergy(c,T,Phi,0d0,0d0,N_phases+2) + potential(phi,c,N_phases+2)&
				 ) +D*(c)+D_cc*T



!  centre
 c=vbles(N_phases+2,2,2)
 T=vbles(N_phases+1,2,2)

  Phi(1)=vbles(1,2,2)

      SoluteRHS = SoluteRHS - D_c*(FreeEnergy(c,T,Phi,0d0,0d0,N_phases+2) + potential(phi,c,N_phases+2)&
				   ) - D*(c)-D_cc*T




