! The 0.34 = W/(R*T_melt*molpervol)=7e7/(8.31*456*54185)
 D=0.
 D_c=0.
 D_cc=0.

 do ip=1,N_phases
  D= D + 0.5 * D1(ip)*( (1.+vbles(4,ih,jv))*vbles(ip,ih,jv)+(1.+vbles(4,2,2))*vbles(ip,2,2))
      D_c=D_c + 0.5* D1(ip)*( vbles(ip,ih,jv)*f_c(vbles(5,ih,jv),0)+vbles(ip,2,2)*f_c(vbles(5,2,2),0))
      D_cc=D_cc + 0.5* D1(ip)*( vbles(ip,ih,jv)*f_cc(vbles(5,ih,jv),0)+vbles(ip,2,2)*f_cc(vbles(5,2,2),0))
 enddo
 c=vbles(5,ih,jv)
 T=vbles(4,ih,jv)
 do ip=1,N_phases
  Phi(ip)=vbles(ip,ih,jv)
 enddo
      SoluteRHS = SoluteRHS+D_c*(FreeEnergy(c,T,Phi,0d0,0d0,N_phases+2) + potential(phi,c,N_phases+2) ) +D*(c)+D_cc*T
      lamda=vbles(9,ih,jv)
       SoluteRHS = SoluteRHS+D_c*lamda

!  centre
 c=vbles(5,2,2)
 T=vbles(4,2,2)
 do ip=1,N_phases
  Phi(ip)=vbles(ip,2,2)
 enddo
      SoluteRHS = SoluteRHS - D_c*(FreeEnergy(c,T,Phi,0d0,0d0,N_phases+2) + potential(phi,c,N_phases+2)) - D*(c)-D_cc*T
      lamda=vbles(9,2,2)
       SoluteRHS = SoluteRHS - D_c*lamda



       

				   
