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


subroutine app_initial_soln_blk(lb)
integer, intent(in) :: lb
logical :: single_phase=.false.,multi_phase=.true.
if(single_phase)then
   call mp_initial_soln_single(lb)
elseif(multi_phase)then
   call mp_initial_soln_multi(lb)
else
   write(*,*)"initial condition"
endif
end subroutine app_initial_soln_blk

subroutine mp_initial_soln_multi(lb)
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use solution_parameters
  implicit none
  include 'mpif.h'
  
  integer, intent(in) :: lb
  integer :: nguard0, i, j, k,ip,ii,jj,i1,i2,ismooth
  double precision :: radius, pi, dx, dy, dz, xi, yi, zi,u1,u2,u3
  double precision :: steepness=1.0,snap21,eut_ic
   double precision,dimension(32) ::p =0
!   double precision :: C_min=0.252,C_ave=0.73,C_max=.98, C_ratio
   double precision :: C_ave=0.70, C_ratio
  double precision :: phi_(3),c_,FreeEnergyPhi,FF(3,10000),df(10000),dg(10000),S,h,n=9999,dh
  double precision :: psi(2),x,a=0.5,phi_psi
  a=g_alpha
  c_ave=g_c0
  if(init_c)then
  dh=1d0/(n+1d0)
  do ii=1,n
     do jj=1,3
        phi_= 0d0
        phi_(jj) = 1d0
        c_ = ii*dh
        FF(jj,ii)=FreeEnergyPhi(c_,delta*0.5,phi_,0d0,0d0,0)
     enddo
!     write(*,'(3F12.2X)')FF
  enddo
  do ii=2,n-2
     df(ii)=(FF(3,ii+1)-FF(3,ii-1))/dh/2
     dg(ii)=(FF(2,ii+1)-FF(2,ii-1))/dh/2
  enddo
  h=1d9
  do ii=2,n/3
     do jj=ii+1,n-2
!        (df[a]-dg[b])^2+(df[a]-(g[b]-f[a])/(b*dx-a*dx))^2

        S = (df(jj)-dg(ii))**2+(df(jj)-(FF(2,ii)-FF(3,jj))/(ii*dh-jj*dh))**2
        if(S.lt.h)then
           h = S
           c_min = ii*dh
           c_max = jj*dh
        endif
     enddo
  enddo
  write(*,*) c_min,c_max
  init_c=.false.
  endif



  C_ratio=1.-(c_max-c_ave)/(C_max-C_min)



  
  pi=4.0*datan(1.0)
  snap21=0.5+0.5*tanh(steepness*nuc_radius)
  
  nguard0 = nguard*npgs
  unk(:,:,:,:,lb) =0.
! set values for unk
  dx = bsize(1,lb)/real(nxb)
  do k=kl_bnd+nguard0*k3d,ku_bnd-nguard0*k3d
     zi = bnd_box(1,3,lb) + dx*(real(k-nguard0)-.5) 
     if (ndim.eq.2) zi = 0.0
     do j=jl_bnd+nguard0*k2d,ju_bnd-nguard0*k2d
        yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5)
        do i=il_bnd+nguard0,iu_bnd-nguard0
           xi = bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5) 
           do ii=1,nunkvbles-1
           

!!!Eutectic

              u1 =1.0-(.5 -0.5*tanh(steepness*(xi-nuc_radius)))
              u2 =(1.0  -u1) *eut_ic(yi,c_ratio)
              u3=1.0 - u1 - u2

!              u1 =1.0-(.5 -0.5*tanh(steepness*(sqrt(xi**2+yi**2)-nuc_radius)))
!              u2 = 1d0-u1
!              u3=0d0
      
              unk(ii+0*nunkvbles,i,j,k,lb)=(u2-u3)/(a*(u2+u3)+1d0-a)/2d0+0.5d0
              unk(ii+1*nunkvbles,i,j,k,lb)=(u3-u1)/(a*(u3+u1)+1d0-a)/2d0+0.5d0
              psi(1)=unk(ii+0*nunkvbles,i,j,k,lb)
              psi(2)=unk(ii+1*nunkvbles,i,j,k,lb)




              x=4*a**2*psi(1)*psi(2)-2*a**2*psi(1)-2*a**2*psi(2)+a**2-2*a*psi(1)+2*a*psi(2)+3

              if(x.lt.1d-10)then
                 write(*,*)"function not defined in initial condition",u1,u2,u3
                 stop
              endif
!              unk(ii+2*nunkvbles,i,j,k,lb)=u3
              ! Temperature            
              !            unk(ii+(N_phases)*nunkvbles,i,j,k,lb) = delta
              
              unk(ii+(N_unk-2)*nunkvbles,i,j,k,lb) =delta
              
              ! solute
              if(int(g_grid/width)*width.ge. yi)then
                 unk(ii+(N_unk-1)*nunkvbles,i,j,k,lb) =C_ave+u3*(c_max-c_ave)+u2*(c_min-c_ave)
              else
                 unk(ii+(N_unk-1)*nunkvbles,i,j,k,lb) =C_ave
              endif
           enddo
        enddo
     enddo
  enddo
  return

end subroutine mp_initial_soln_multi




subroutine mp_initial_soln_single(lb)  
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
  integer :: nguard0, i, j, k,ip,ii
  double precision :: radius, pi, dx, dy, dz, xi, yi, zi,u1
  double precision c_ave
  c_ave=g_c0




  
  nguard0 = nguard*npgs
  unk(:,:,:,:,lb) =0.
! set values for unk
  dx = bsize(1,lb)/real(nxb)
  do k=kl_bnd+nguard0*k3d,ku_bnd-nguard0*k3d
     zi = bnd_box(1,3,lb) + dx*(real(k-nguard0)-.5) 
     if (ndim.eq.2) zi = 0.0
     do j=jl_bnd+nguard0*k2d,ju_bnd-nguard0*k2d
        yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5)
        do i=il_bnd+nguard0,iu_bnd-nguard0
           xi = bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5)
           do ii=1,nunkvbles-1
!phase (liquid)
              if(g_1D)then
                 unk(ii,i,j,k,lb) = 0.5 +0.5*tanh(xi-nuc_radius)
              else
                 unk(ii,i,j,k,lb) = 0.5 +0.5*tanh((xi**2+yi**2)**0.5-nuc_radius)
              endif
              unk(ii+nunkvbles,i,j,k,lb) = 1d0 - unk(ii,i,j,k,lb)
              unk(ii+2*nunkvbles,i,j,k,lb) = 0d0

              
! Temperature            
              unk(ii+(N_unk-2)*nunkvbles,i,j,k,lb) =delta !overwrites if N_phases=2
! solute
              unk(ii+(N_unk-1)*nunkvbles,i,j,k,lb) =c_ave


           enddo
        enddo
     enddo
  enddo
 
  return

end subroutine mp_initial_soln_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_initial_soln_end                                   
  return
end subroutine app_initial_soln_end


subroutine app_test_refinement
  return
end subroutine app_test_refinement

double precision function eut_ic(x,phi2)
use solution_parameters
double precision,intent (in):: x,phi2
double precision ::phi1,xx,q(12)=0
!double precision ::width=37.05242907,phi1,xx,q(12)=0



phi1=1.-phi2

eut_ic =1.
xx=x-int(x/width)*width
if(xx.lt.width*phi2*0.5)eut_ic = 0. 
if(xx.gt.width*(phi2*0.5+phi1))eut_ic = 0.
return


end function eut_ic

