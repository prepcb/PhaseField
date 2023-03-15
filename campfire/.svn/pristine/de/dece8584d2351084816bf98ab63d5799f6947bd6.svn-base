!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_initial_soln_blk                                         REQUIRED
!!!!  * Sets the initial solution on a block
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_test_refinement                                          REQUIRED
!!!!  * Application specific refinement tests - unused in phase field code
!!!!  * Only called if local_adaptation=2
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Latest version: Chris Goodyer, May 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



double precision function phi_f(x,d,c)
implicit none
double precision, intent(in) :: x,d
integer, intent(in)::c
double precision :: phi
if(c.eq.0) then
!   phi_f = 0.5*(1d0+tanh(2d0*x/d))
   phi_f = 1./(1.+exp(-4.*x/d))
else
!   phi_f = (1d0-tanh(2d0*x/d)**2)/d
   phi = 1./(1.+exp(-4.*x/d))
   phi_f = 4.*phi*(1.-phi)/d

   
endif
end function phi_f
subroutine app_initial_soln_blk(lb)
  
!--------------------------------------------------------------
! include files 
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use solution_parameters
  implicit none
  include 'mpif.h'
  
  integer, intent(in) :: lb
  integer :: nguard0, i, j, k,ip,ii,jj
  double precision :: radius, pi, dx, dy, dz, xi, yi, zi,u1
  double precision c_ave,ci,f11,f12,f21,f22,c_mid,c0,c1,dc0,dc1,c2,a1,a2,a0,x0,x1,xx0,xx1
  double precision :: phi_,d_phi,c_,FreeEnergy,S,h,dh,SS,hh, f_s
  double precision :: d=4.27,x_L,x_R,alpha,phi_f,D_sol,f_,f_cc,f_c_phi,c_f(80000),ff(80000),u=0.01,a_L,a_R
  integer :: n=999,nn=9999
  double precision :: Gibbs_FE_l_AlSi,Gibbs_FE_e_AlSi,c_l,c_s,c_a,c_b,T_,gh_phi
  double precision :: Gibbs_FE_e,Gibbs_FE_l,c_min_,c_max_,phi0,Gibbs_FE_e_PbSn,Gibbs_FE_l_PbSn,fe(2,1000)
  logical :: AdamsBashforth = .false.

  c_ave = g_c0

! g_vel is the estimated final tip speed
! g_ratio is (actual width/ equilibrium width)  usually < 1. modified periodically.
  u = g_vel

  if(g_obstacle.le.0d0)then
     d = sqrt(8d0) * g_lambda * g_ratio
  else
     d = 1.5*2d0/sqrt(g_obstacle) * g_lambda * g_ratio
  endif

!FreeEnergy(c,T,phi,c_dot,phi_dot,lp,approx)
  if(AlNi.ne.8)then
     if(g_c0.gt.0d0)then
        do i=1,999
           c0=i*0.001
           phi_=1d0
           fe(1,i)=FreeEnergy(c0,0d0,phi_,0,.false.)
           phi_=0d0
           fe(2,i)=FreeEnergy(c0,0d0,phi_,0,.false.)
           
        enddo
        
        open (unit=199,file="energy_out.txt")
        do ii=1,999
           write(199,'(3F14.5)')(fe(jj,ii),jj=1,2)
        enddo
        close(199)
     endif
  endif


  if(init_c)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     T_=0d0!i*0.01 - 0.2
     ! C_Min c_max finder
     
     dh=1d0/(real(n)+1d0)
     h=1d9
     hh=1d9
     
     

     
     
!!!!!!!!!!
     if(AlNi.ne.8)then
        write(*,*)"unknown AlNi value mp_initial_soln"
        stop
     endif
     
     
     init_c=.false.
  endif





!!!!!!!!!!!!
  g_alpha0 = g_alpha

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nguard0 = nguard*npgs
  unk(:,:,:,:,lb) =0.
  ! set values for unk

  dx = bsize(1,lb)/real(nxb)



  k=1
  if(g_linear)then
     do j=jl_bnd+nguard0*k2d,ju_bnd-nguard0*k2d
        yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5)
        do i=il_bnd+nguard0,iu_bnd-nguard0
           xi = bnd_box(1,1,lb)+ dx*(real(i-nguard0)-.5)
           do ii=1,nunkvbles-1
              unk(ii,i,j,k,lb)=1./(1+exp(-4*(sqrt((xi-g_x0)**2+(yi-g_y0)**2)-nuc_radius)/d) &
!!$                    exp(-4*(sqrt((xi-2*g_x0)**2+(yi-g_y0)**2)-nuc_radius)/d) + &
!!$                   exp(-4*(sqrt((xi-3*g_x0)**2+(yi-g_y0)**2)-1.5*nuc_radius)/d)+&
!!$                   exp(-4*(sqrt((xi-4*g_x0)**2+(yi-g_y0)**2)-.5*nuc_radius)/d)+&
!!$                   exp(-4*(sqrt((xi-5*g_x0)**2+(yi-g_y0)**2)-nuc_radius)/d) + &
!!$                   exp(-4*(sqrt((xi-6*g_x0)**2+(yi-g_y0)**2)-1.6*nuc_radius)/d)+&
!!$                   exp(-4*(sqrt((xi-7*g_x0)**2+(yi-g_y0)**2)-.6*nuc_radius)/d)+&
!!$                   exp(-4*(sqrt((xi-8*g_x0)**2+(yi-g_y0)**2)-nuc_radius)/d) + &
!!$                   exp(-4*(sqrt((xi-9*g_x0)**2+(yi-g_y0)**2)-1.5*nuc_radius)/d)+&
!!$                   exp(-4*(sqrt((xi-10*g_x0)**2+(yi-g_y0)**2)-2.5*nuc_radius)/d)&
)



!!$              if(xi.gt.nuc_radius.or. yi.gt.8*nuc_radius)then
!!$                 unk(ii,i,j,k,lb)=1d0
!!$              else
!!$                 unk(ii,i,j,k,lb)=0d0
!!$              endif
              c_ave = g_cinfty*g_AlFe_c(1,1)+(1.0-g_cinfty)*g_AlFe_c(2,1)
!              unk(ii+nunkvbles,i,j,k,lb) = (1d0-unk(ii,i,j,k,lb))*g_AlFe_c(1,1)+unk(ii,i,j,k,lb)*c_ave
!              unk(ii+nunkvbles,i,j,k,lb) = (1.0-gh_phi(0.5*unk(ii,i,j,k,lb),0))*g_AlFe_c(1,1)+gh_phi(0.5*unk(ii,i,j,k,lb),0)*g_AlFe_c(2,1)
              unk(ii+nunkvbles,i,j,k,lb) = (1.0-gh_phi(g_cinfty*unk(ii,i,j,k,lb),0))*g_AlFe_c(1,1)+gh_phi(g_cinfty*unk(ii,i,j,k,lb),0)*g_AlFe_c(2,1)
!              unk(ii+nunkvbles,i,j,k,lb) = 0.5*(g_AlFe_c(1,1)+g_AlFe_c(2,1))
           enddo
        enddo
     enddo
  else
     do j=2,g_order-1
        yi = bnd_box(1,2,lb) +(0.5+g_y10(j)/(1.+g_y10(g_order-1)))*nxb*dx
        do i=2,g_order-1
           xi = bnd_box(1,1,lb)+(0.5+g_y10(i)/(1.+g_y10(g_order-1)))*nxb*dx
           do ii=1,nunkvbles-1
!              unk(ii,i,j,k,lb)=1./(1+exp(-4*(sqrt((xi-g_x0)**2+(yi-g_y0)**2)-nuc_radius)/d))

              unk(ii,i,j,k,lb)=1./(1+exp(-4*(sqrt((xi-g_x0)**2+(yi-g_y0)**2)-nuc_radius)/d))

              unk(ii,i,j,k,lb)=1./(1+exp(-4*(sqrt((xi-2*g_x0)**2+(yi-2*g_y0)**2)-nuc_radius)/d))
              c_ave = 0.5*(g_AlFe_c(1,1)+g_AlFe_c(2,1))
              unk(ii+nunkvbles,i,j,k,lb) = (1d0-unk(ii,i,j,k,lb))*g_AlFe_c(1,1)+unk(ii,i,j,k,lb)*c_ave
           enddo
        enddo
     enddo
  endif


  return

end subroutine app_initial_soln_blk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_initial_soln_end                                   
  return
end subroutine app_initial_soln_end


subroutine app_test_refinement
  return
end subroutine app_test_refinement


!!!!!!!!!!!!!!!!
subroutine find_antitrapping_AlFe
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use solution_parameters
  implicit none
  include 'mpif.h'


  integer :: nguard0, i, j, k,ip,ii,jj
  double precision :: radius, pi, dx, dy, dz, xi, yi, zi,u1
  double precision c_ave,ci,f11,f12,f21,f22,c_mid,c0,c1,c2,xx
  double precision :: phi_,d_phi,c_,FreeEnergy,S,h,dh,SS,hh, f_s
  double precision :: x_L,x_R,alpha,phi_f,D_sol,f_,f_cc,f_c_phi,c_f(80000),ff(80000),a_L,a_R,d,u
  integer :: n=999,nn=9999
  double precision :: Gibbs_FE_l_AlSi,Gibbs_FE_e_AlSi,c_l,c_s,c_a,c_b,T_,alpha_tmp,T_i
  double precision :: Gibbs_FE_e,Gibbs_FE_l,c_min_,c_max_,phi0,Gibbs_FE_e_PbSn,Gibbs_FE_l_PbSn
  logical :: AdamsBashforth = .false.,phi_method=.false.



  d = sqrt(8d0) * g_lambda * g_ratio

  u = g_vel
  write(*,*)"velocity =",u



  if(u.lt.0d0)then
     g_alpha=1d-6
     return
  endif
  
  if(g_alpha.ne.0d0)then !else use input value 
     
     g_alpha=0
     alpha_tmp = g_alpha
     
     
     !ODE solver for 1D
     if(phi_method)then
        x_L = 1d-9
        x_R =  1d0 - 1d-9
     else
        x_R =  4d0*d
        x_L = -4d0*d
     endif
     
     
     dh = (x_R-x_L)/real(nn)
     !        dh = g_grid*g_nbx*2d0**(-1-g_max_lvl)/nxb
     a_L = -1d1
     a_R = 1d1
     do jj=1,50
        g_alpha =0.5*(a_L+a_R)
        
        !ODE solver main loop  
        g_c_max = g_AlFe_c(2,1)
        g_c_min = g_AlFe_c(1,1)
        
        c_ = g_c_min
        c_f(1) = c_
        
        do ii =2,nn
           h = x_L + ii*dh !x coord
           
           if(phi_method)then
              phi_=h !phi
              d_phi = 4.*h*(1.-h)/d
           else
              phi_=phi_f(h,d,0) !phi
              d_phi = phi_f(h,d,1) !d_phi/d_x
           endif
           
           
           
           D_sol = ((1d0-phi_)*g_D(2) + g_D(1)*phi_)*c_f(ii-1)*(1d0-c_f(ii-1))/g_Dch
           f_cc =FreeEnergy(c_f(ii-1),0d0,phi_,-2,.true.)/(g_R*g_T0)
           f_c_phi = FreeEnergy(c_f(ii-1),0d0,phi_,-1,.true.)/(g_R*g_T0)
           
           if(phi_method)then
              ff(ii) = (dh/f_cc)*(f_c_phi*(1d0+g_alpha*u*d/(g_D(1))*g_Dch)&
                   +u/d_phi*(c_f(ii-1)-g_c_max)/D_sol)
           else
              ff(ii) = (dh/f_cc)*(f_c_phi*(1d0+g_alpha*u*d/(g_D(1))*g_Dch)*d_phi&
                   +u*(c_f(ii-1)-g_c_min)/D_sol)
           endif
           
           
           c_f(ii) = c_f(ii-1) - ff(ii)
           
           if(c_f(ii).gt.c_)then
              c_ = c_f(ii)
              if(phi_method)then
                 xx = h!phi_
              else
                 xx = phi_
              endif
           endif
        enddo
        
        if(c_.ge.g_c_max)then
           a_R = g_alpha
        else
           a_L = g_alpha
        endif

     enddo
     write(*,'(4F17.5)')c_,g_alpha,xx,u

     if(g_alpha.gt.9d0)g_alpha=alpha_tmp
  endif
  open (unit=199,file="energy_cf.txt")
  do ii=2,nn
     write(199,'(3F17.7)')x_L+ii*dh,c_f(ii),phi_f(x_L + ii*dh,d,0)
  enddo
  close(199)


end subroutine

