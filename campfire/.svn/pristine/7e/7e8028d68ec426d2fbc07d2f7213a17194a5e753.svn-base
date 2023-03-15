! The 0.34 = W/(R*T_melt*molpervol)=7e7/(8.31*456*54185)

  D=0.
  D_c=0.
  D_cc=0.
  
  D= 0.5 * D1(1)*( (1.+delta)*vbles(1,ih,jv)+(1.+delta)*vbles(1,2,2))
  D_c=0.5* D1(1)*( vbles(1,ih,jv)*f_c(vbles(2,ih,jv),0) + vbles(1,2,2)*f_c(vbles(2,2,2),0))
  D_cc= 0.5* D1(1)*( vbles(1,ih,jv)*f_cc(vbles(2,ih,jv),0) + vbles(1,2,2)*f_cc(vbles(2,2,2),0))

  if(beta.gt.1d-6) beta1 = beta*lambda*8d0*(vbles(1,2,2)**2*(1d0-vbles(1,2,2))**2+vbles(1,ih,jv)**2*(1d0-vbles(1,ih,jv))**2)
  
  

  
  
  c=vbles(2,ih,jv)
!  T= delta
  Phi=vbles(1,ih,jv)
  psi=vbles(3,ih,jv)
		     
  


  SoluteRHS = SoluteRHS+D_c*(FreeEnergy(c,T,Phi,c_dot,phi_dot,2)) +D*c+D_cc*T

 if(beta.gt.1d-6)SoluteRHS = SoluteRHS + D_c*0.34*beta1*psi*1d3


!  centre
  c=vbles(2,2,2)
  Phi=vbles(1,2,2)
  psi=vbles(3,2,2)

  SoluteRHS = SoluteRHS - D_c*(FreeEnergy(c,T,Phi,c_dot,phi_dot,2)) - D*c-D_cc*T

  if(beta.gt.1d-6)SoluteRHS = SoluteRHS - D_c*0.34*beta1*psi*1d3


				   

