
!!$  D_c=0.5* ( (vbles(1,ih,jv,kd)*g_D(1)+(1.-vbles(1,ih,jv,kd))*g_D(2))*f_c(vbles(2,ih,jv,kd),0)&
!!$       +(vbles(1,2,2,2)*g_D(1)+(1.-vbles(1,2,2,2))*g_D(2))*f_c(vbles(2,2,2,2),0))

!!$  D_c=0.5* ( vbles(1,ih,jv,kd)*g_D(1)*f_c(vbles(2,ih,jv,kd),0)+(1.-vbles(1,ih,jv,kd))*g_D(2)*f_c(vbles(2,ih,jv,kd),0)&
!!$       +vbles(1,2,2,2)*g_D(1)*f_c(vbles(2,2,2,2),0)+(1.-vbles(1,2,2,2))*g_D(2)*f_c(vbles(2,2,2,2),0))
!!$
!!$  D_c=0.5* ( vbles(1,ih,jv,kd)*g_D(1)*f_c(vbles(2,ih,jv,kd),0)+(1.-vbles(1,ih,jv,kd))*g_D(2)*f_c(AlFey1(vbles(2,ih,jv,kd),0),0)&
!!$       +vbles(1,2,2,2)*g_D(1)*f_c(vbles(2,2,2,2),0)+(1.-vbles(1,2,2,2))*g_D(2)*f_c(AlFey1(vbles(2,2,2,2),0),0))!  6th July


! since we use quadratic free energy it is not necessary to have c(1-c) as a factor in the diffusion
if(g_plapp)then
   c_neigh = vbles(2,ih,jv,kd)
   c_centr = vbles(2,2,2,2)
      
!!$   y_neigh = max(0d0,min(1d0,(c_neigh - g_AlFe_c(1,1))/(1d0 - g_AlFe_c(1,1))))
!!$   y_centr = max(0d0,min(1d0,(c_centr - g_AlFe_c(1,1))/(1d0 - g_AlFe_c(1,1))))
   phi_centr = vbles(1,2,2,2)
   phi_neigh = vbles(1,ih,jv,kd)
   D_c=0.5* ( phi_neigh*g_D(1)+(1.-phi_neigh)*g_D(2) + phi_centr*g_D(1) + (1.-phi_centr)*g_D(2)) !no c(1-c) term

   SoluteRHS = SoluteRHS + D_c*(c_neigh - c_centr)

else
   if(g_no_log_term)then
      D_c=0.5* ( vbles(1,ih,jv,kd)*g_D(1)+(1.-vbles(1,ih,jv,kd))*g_D(2)&
           +vbles(1,2,2,2)*g_D(1) + (1.-vbles(1,2,2,2))*g_D(2)) ! additional condition added 24th Sep
   else
      c_neigh = vbles(2,ih,jv,kd)
      c_centr = vbles(2,2,2,kd)
      
      y_neigh = max(0d0,min(1d0,(c_neigh - g_AlFe_c(1,1))/(1d0 - g_AlFe_c(1,1))))
      y_centr = max(0d0,min(1d0,(c_centr - g_AlFe_c(1,1))/(1d0 - g_AlFe_c(1,1))))
      
      if(g_antitrapping)then
         
         phi_centr=sharpen(vbles(1,2,2,2), n_sharpen)
         phi_neigh=sharpen(vbles(1,ih,jv,kd), n_sharpen)
         
         
!!$     if(vbles(1,2,2,2).gt.0.6)then
!!$        phi_centr = 1d0
!!$     elseif(vbles(1,2,2,2).lt.0.4)then
!!$        phi_centr=0d0
!!$     else
!!$        phi_centr = (vbles(1,2,2,2)-0.5)/0.2+0.5
!!$     endif
!!$     if(vbles(1,ih,jv,kd).gt.0.6)then
!!$        phi_neigh = 1d0
!!$     elseif(vbles(1,ih,jv,kd).lt.0.4)then
!!$        phi_neigh=0d0
!!$     else
!!$        phi_neigh = (vbles(1,ih,jv,kd)-0.5)/0.2+0.5
!!$     endif
!!$  else  
         
         phi_centr = vbles(1,2,2,2)
         phi_neigh = vbles(1,ih,jv,kd)
      endif
      
      
!!$      D_c = 0.5* ( (phi_neigh*g_D(1)+(1.-phi_neigh)*g_D(2))*f_c(y_neigh,0)&                                                                             
!!$           +(phi_centr*g_D(1)+(1.-phi_centr)*g_D(2))*f_c(y_centr,0)) 

      D_c = 0.5* ( (phi_neigh*g_D(1)+(1.-phi_neigh)*g_D(2))&                                                                             
           +(phi_centr*g_D(1)+(1.-phi_centr)*g_D(2))) 
   endif
   
   
   
   !  centre
   c=vbles(2,2,2,2)
   T=0d0
   Phi=vbles(1,2,2,2)
   
   if(g_antitrapping)then
      phi = sharpen(vbles(1,2,2,2), n_sharpen)
!!$
!!$
!!$     if(vbles(1,2,2,2).gt.0.6)then
!!$        phi = 1d0
!!$     elseif(vbles(1,2,2,2).lt.0.4)then
!!$        phi=0d0
!!$     else
!!$        phi = (vbles(1,2,2,2)-0.5)/0.2+0.5
!!$     endif
   endif
   
   SoluteRHS = SoluteRHS - D_c*FreeEnergy(c,T,Phi,2,g_approx)
   
!!$  if(g_antitrapping.and..false.)then
!!$
!!$     norm =   (vbles(1,2,2,2)-vbles(1,ih,2,2))**2
!!$     norm = norm + (vbles(1,2,2,2)-vbles(1,2,jv,2))**2
!!$     norm = norm + (vbles(1,2,2,2)-vbles(1,2,2,kd))**2
!!$     norm = (sqrt(norm)+1d-3)/dx
!!$
!!$     if (mode_1.eq.1) then
!!$        phi_t = (unk(1,i,j,k,lb)-unk(2,i,j,k,lb))/dt
!!$        phi_t = phi_t+(unk(1,i-2+ih,j-2+jv,k-2+kd,lb)-unk(2,i-2+ih,j-2+jv,k-2+kd,lb))/dt
!!$        phi_t = phi_t*0.5
!!$        
!!$     elseif(mode_1.eq.2)then
!!$        phi_t = (rf4*unk(5,i,j,k,lb)-rf5*unk(4,i,j,k,lb)+rf6*unk(1,i,j,k,lb))/dt
!!$        phi_t = phi_t+(rf4*unk(5,i-2+ih,j-2+jv,k-2+kd,lb)-rf5*unk(4,i-2+ih,j-2+jv,k-2+kd,lb)+rf6*unk(1,i-2+ih,j-2+jv,k-2+kd,lb))/dt
!!$        phi_t = phi_t*0.5
!!$     endif
!!$g_width = max(1d-6,g_width)
!!$     SoluteRHS = SoluteRHS + D_c/g_D(1)*g_Dch*a_trap*g_velocity*phi
!!$
!!$  endif

  ! neighbours
   c=vbles(2,ih,jv,kd)
   T=0d0
   Phi=vbles(1,ih,jv,kd)
   if(g_antitrapping)then
      
      phi = sharpen(vbles(1,ih,jv,kd), n_sharpen)
!!$     if(vbles(1,ih,jv,kd).gt.0.6)then
!!$        phi = 1d0
!!$     elseif(vbles(1,ih,jv,kd).lt.0.4)then
!!$        phi=0d0
!!$     else
!!$        phi = (vbles(1,ih,jv,kd)-0.5)/0.2+0.5
!!$     endif
   endif
   
   
   
   SoluteRHS = SoluteRHS + D_c*FreeEnergy(c,T,Phi,2,g_approx)
   
   if(g_antitrapping.and..false.)then
      SoluteRHS = SoluteRHS - D_c/g_D(1)*g_Dch*a_trap*g_velocity*phi
   endif
   
!!$  if(g_alpha.ge.0d0)then
!!$     SoluteRHS = SoluteRHS + D_c*a_trap*u*FE*phi
!!$  endif
endif
