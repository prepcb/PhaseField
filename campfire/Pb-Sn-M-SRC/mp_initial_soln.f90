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
  double precision :: radius, pi, dx, dy, dz, xi, yi, zi,u4,u5,u1,u2,u3,v1,v2,v3
  double precision :: steepness=2.0,snap21,eut_ic,angle_ic
  double precision,dimension(32) ::p =0
  double precision :: C_ave=0.70, C_ratio
  double precision :: phi_(3),c_,c1_,FreeEnergy,FF(3,100000),df(100000),dg(100000),S,h,n=9999.0,dh
  logical :: angle =.true.
  double precision :: phi1(3)=(/1,0,0/),phi2(3)=(/0,1,0/),phi3(3)=(/0,0,1/)
  double precision :: n23(2)=(/0.0,1.0/),n32(2)=(/0.0,-1.0/)
  double precision :: n31(2)=(/ .8660254035,-0.5/),n13(2)=(/ -.8660254035,0.5/)
  double precision :: n12(2)=(/ -.8660254035,-0.5/),n21(2)=(/ .8660254035,0.5/)
  double precision :: xx,yy,theta, T_k,Tk,Latent(3),c1=0.35,dc1=0.05,c2=0.45,dc2=0.05

theta = 180 -180./3.14159*2*asin(1/sqrt(2+2*g_d1**2))

write(*,*)"theta=",theta

theta = theta/180.*3.14159*0.5

  n31(1)=sin(theta)
  n31(2)=-cos(theta)
  n13=-n31
  n12(1)=-sin(theta)
  n12(2)=-cos(theta)
  n21 = -n12
  n23 = n21 + n13
  n32 = -n23


  steepness = 8.! 1.0/g_lambda

!  steepness = sqrt(2.)/g_lambda
  k=int(n)
  write(*,*)"n=",k
  c_ave=g_c0
  write(*,*)"model",g_model

!  if(g_model.eq.0)then
  if(.false.)then
     if(init_c)then
        dh=1d0/(n+1d0)
        do ii=1,k
!!$           do jj=1,3
!!$              phi_= 0d0
!!$              phi_(jj) = 1d0
!!$              c_ = c1 + ii*dh
!!$              FF(jj,ii)=FreeEnergy(c_,delta,phi_,0d0,0d0,0)
!!$           enddo
           FF(1,ii)=FreeEnergy(c1+ii*dh*dc1,delta,phi1,0d0,0d0,0)
           FF(2,ii)=FreeEnergy(c2+ii*dh*dc2,delta,phi2,0d0,0d0,0)
        enddo
        if(.false.)then
           open (unit=199,file="energy_out.txt")
           do ii=1,k
              write(199,'(3F14.5)')(FF(jj,ii),jj=1,3)
           enddo
           close(199)
           stop
        endif
        do ii=2,n-1
!           df(ii)=0.5*(FF(1,ii+1)-FF(1,ii-1))/dh
!           dg(ii)=0.5*(FF(2,ii+1)-FF(2,ii-1))/dh
           df(ii)=FreeEnergy(c1+ii*dh*dc1,delta,phi1,0d0,0d0,5)
           dg(ii)=FreeEnergy(c2+ii*dh*dc2,delta,phi2,0d0,0d0,5)
        enddo
        h=1d9
        !  do ii=2,n/3
        !     do jj=ii+1,n-2
        do ii=2,k-1
           do jj=2,k-1
!              if(ii.ne.jj)then
                 !        (df[a]-dg[b])^2+(df[a]-(g[b]-f[a])/(b*dx-a*dx))^2
                 
!                 S = (df(jj)-dg(ii))**2+(df(jj)-(FF(1,ii)-FF(2,jj))/(ii*dh-jj*dh))**2
                 S = (df(jj)-dg(ii))**2+(dg(ii)-(FF(1,ii)-FF(2,jj))/(c1+ii*dh*dc1-c2-jj*dh*dc2))**2
                 if(S.lt.h)then
                    h = S
!                    write(*,*)h
                    c_max = c2+ii*dh*dc2
                    c_min = c1+jj*dh*dc1
 !                endif
              endif
           enddo
        enddo
        write(*,*) c_min,c_max

        init_c=.false.
     endif



  C_ratio=1.-(c_max-c_ave)/(c_max-c_min)


  endif


!FreeEnergy(c_,delta*0.5,phi_,0d0,0d0,0)

!    (c,T,phi,c_dot,phi_dot,lp)
if(.true.)then
open(unit=99,file='c_out.txt')
!write(99,'(4E14.7)')c_min,FreeEnergy(c_min,0d0,phi1,0d0,0d0,0),c_max,FreeEnergy(c_max,0d0,phi2,0d0,0d0,0)
do ii = 1,999
   c_=ii*0.001!*(0.27246-0.2350)+0.2350
   c1_=ii*0.001*(0.27246-0.2350)+0.2350

!!$   write(*,*)FreeEnergy(c_,0d0,phi1,0d0,0d0,0),FreeEnergy(c_,0d0,phi2,0d0,0d0,0),FreeEnergy(c_,0d0,phi3,0d0,0d0,0)
!!$   g_Tflag = .true.
!!$   Tk=T_k(0d0)
!!$   Latent(3) = Tk*(FreeEnergy(c_,0d0,phi1,0d0,0d0,1)-FreeEnergy(c_,0d0,phi2,0d0,0d0,1))
!!$   Latent(1) = Tk*(FreeEnergy(c_,0d0,phi2,0d0,0d0,1)-FreeEnergy(c_,0d0,phi3,0d0,0d0,1))
!!$   Latent(2) = Tk*(FreeEnergy(c_,0d0,phi3,0d0,0d0,1)-FreeEnergy(c_,0d0,phi1,0d0,0d0,1))
!!$   g_Tflag = .false.
!!$   Latent(3) = Latent(3) - (FreeEnergy(c_,0d0,phi1,0d0,0d0,0)-FreeEnergy(c_,0d0,phi2,0d0,0d0,0))
!!$   Latent(1) = Latent(1) - (FreeEnergy(c_,0d0,phi2,0d0,0d0,0)-FreeEnergy(c_,0d0,phi3,0d0,0d0,0))
!!$   Latent(2) = Latent(2) - (FreeEnergy(c_,0d0,phi3,0d0,0d0,0)-FreeEnergy(c_,0d0,phi1,0d0,0d0,0))
 write(99,'(4E14.7)')c_,FreeEnergy(c_,0d0,phi1,0d0,0d0,0),FreeEnergy(c_,0d0,phi2,0d0,0d0,0),FreeEnergy(c1_,0d0,phi3,0d0,0d0,0)
!!$ write(99,'(4E14.7)')c_,abs(Latent)
enddo
close(99)
write(*,*)"done c_out.txt"

endif






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
!     do j=jl_bnd,ju_bnd
        yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5)
       do i=il_bnd+nguard0,iu_bnd-nguard0
!       do i=il_bnd,iu_bnd
           xi = bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5) 
           do ii=1,nunkvbles-1

              if(g_model.eq.3)then
                 c_max=0.75
                 c_min=0.25
                 c_ave=0.5
!!$                    u1 =1.0-(.5 -0.5*tanh(steepness*(xi-nuc_radius)))
!!$                    u2 =(1.0  -u1) *eut_ic(yi,0.5)
!!$                    u3 = 1-u1-u2

                    u1 = 1.0/(1+exp(1.0/(g_dlta(1)*g_lambda)*((xi-nuc_radius)*n12(1)+(yi-g_grid*0.5)*n12(2)))+exp(1.0/(g_dlta(1)*g_lambda)*((xi-nuc_radius)*n13(1)+(yi-g_grid*0.5)*n13(2))))
                    u2 = 1.0/(1+exp(1.0/(g_dlta(2)*g_lambda)*((xi-nuc_radius)*n23(1)+(yi-g_grid*0.5)*n23(2)))+exp(1.0/(g_dlta(2)*g_lambda)*((xi-nuc_radius)*n21(1)+(yi-g_grid*0.5)*n21(2))))
                    u3 = 1.0/(1+exp(1.0/(g_dlta(3)*g_lambda)*((xi-nuc_radius)*n31(1)+(yi-g_grid*0.5)*n31(2)))+exp(1.0/(g_dlta(3)*g_lambda)*((xi-nuc_radius)*n32(1)+(yi-g_grid*0.5)*n32(2))))



                 S= u1+u2+u3

                 unk(ii+0*nunkvbles,i,j,k,lb) = u1/S
                 unk(ii+1*nunkvbles,i,j,k,lb) = u3/S
                 unk(ii+2*nunkvbles,i,j,k,lb) = u2/S
                 unk(ii+(N_unk-2)*nunkvbles,i,j,k,lb) =delta
                 unk(ii+(N_unk-1)*nunkvbles,i,j,k,lb) =c_ave+u2*(c_max-c_ave)+u3*(c_min-c_ave)
              elseif(g_model.eq.4)then
!                 c_min=0.381
!                 c_max=0.52
!                 c_ave=(c_min+3*c_max)/4.
!                 g_c0 = c_ave
                 u1 = 1d0/(1d0+exp(-steepness*(sqrt(xi**2+yi**2)-nuc_radius)))
!                 u1 = 0.5+0.5*tanh(1*(sqrt(xi**2+yi**2)-nuc_radius))
                 u2 = 0d0
                 u3 = 1d0-u1
                 unk(ii+0*nunkvbles,i,j,k,lb) = u1
                 unk(ii+1*nunkvbles,i,j,k,lb) = u2
                 unk(ii+2*nunkvbles,i,j,k,lb) = u3
                 unk(ii+(N_unk-2)*nunkvbles,i,j,k,lb) = 0d0
!                 unk(ii+(N_unk-1)*nunkvbles,i,j,k,lb) = u1*0.25 + u3*0.2351
                 unk(ii+(N_unk-1)*nunkvbles,i,j,k,lb)  =  g_c0
              elseif(g_model.eq.2.or.g_model.eq.1)then
                 if(angle)then
                    theta= pi/6d0
                    u1 = 0.5+0.5*tanh(steepness*(sqrt(xi**2+yi**2)-nuc_radius))
                    u2 = 1d0/(1+exp(steepness*(xi*cos(theta)-yi*sin(theta))))
                    u3 = 1d0 - u2
                    u2=u2*(1-u1)
                    u3=u3*(1-u1)

                 else
!                    u1 =1.0-(.5 -0.5*tanh(steepness*(xi-nuc_radius)))
!                    u2 =(1.0  -u1) *eut_ic(yi,0.5)
!                    u3 = 1-u1-u2


                    u1 = 1.0/(1+exp(1.0/(g_dlta(1)*g_lambda)*((xi-2.5)*n12(1)+(yi-2.5)*n12(2)))+exp(1.0/(g_dlta(1)*g_lambda)*((xi-2.5)*n13(1)+(yi-2.5)*n13(2))))
                    u2 = 1.0/(1+exp(1.0/(g_dlta(2)*g_lambda)*((xi-2.5)*n23(1)+(yi-2.5)*n23(2)))+exp(1.0/(g_dlta(2)*g_lambda)*((xi-2.5)*n21(1)+(yi-2.5)*n21(2))))
                    u3 = 1.0/(1+exp(1.0/(g_dlta(3)*g_lambda)*((xi-2.5)*n31(1)+(yi-2.5)*n31(2)))+exp(1.0/(g_dlta(3)*g_lambda)*((xi-2.5)*n32(1)+(yi-2.5)*n32(2))))


                 endif

                 S= u1+u2+u3

                 unk(ii+0*nunkvbles,i,j,k,lb) = u1/S
                 unk(ii+1*nunkvbles,i,j,k,lb) = u2/S
                 unk(ii+2*nunkvbles,i,j,k,lb) = u3/S
                 unk(ii+(N_unk-2)*nunkvbles,i,j,k,lb) =delta
!                 if(g_model.eq.0)unk(ii+(N_unk-1)*nunkvbles,i,j,k,lb) =C_ave+u3*(c_max-c_ave)+u2*(c_min-c_ave)
!                 if(g_model.eq.2) unk(ii+(N_unk-1)*nunkvbles,i,j,k,lb) =g_c0!u2*c_max+u1*c_ave
                 if(g_model.eq.2) unk(ii+(N_unk-1)*nunkvbles,i,j,k,lb) =g_c0+u2*(0.998-g_c0)+u3*(0.1-g_c0)
                 if(g_model.eq.1) unk(ii+(N_unk-1)*nunkvbles,i,j,k,lb) =g_c0 + u2*(c_max-g_c0)+u3*(c_min-g_c0)

                 !             unk(ii+2*nunkvbles,i,j,k,lb) =(1d0-unk(ii,i,j,k,lb))*g_c_max+unk(ii,i,j,k,lb)*c_ave
              elseif(g_model.eq.0.or.g_model.eq.5.or.g_model.eq.6)then
!!!Eutectic
                 write(*,*)"model KKS for Pb-Sn",g_model

!!$                    u1 =1.0-(.5 -0.5*tanh(steepness*(xi-nuc_radius*(1d0 +slant*yi/g_grid/real(g_nby)))))
!!$
!!$                 u1 = 1. - 1/(1+exp(4*steepness*(sqrt((xi)**2+(yi)**2)-nuc_radius)))
!!$!                 u1 = 0.5+0.5*tanh(steepness*(sqrt(xi**2+yi**2)-nuc_radius))
!!$                 u2 = 1d0/(1+exp(4*steepness*((xi)*cos(0.41*pi/2)-(yi)*sin(0.41*pi/2))))
!!$                 u3 = 1d0 - u2
!!$                 u2=u2*(1-u1)
!!$                 u3=u3*(1-u1)

!!$
!                    u1 =1.0-(.5 -0.5*tanh(steepness*(xi-nuc_radius)))
                    
                    u1 = 1-1/(1+exp(2d0*(xi-nuc_radius))+exp(2d0*(yi-4*width)))
                    u2 =(1.0  -u1) *eut_ic(yi,g_c0)
                    u3 = 1-u1-u2

                 S= u1+u2+u3


                 unk(ii+0*nunkvbles,i,j,k,lb) = u1/S
                 unk(ii+1*nunkvbles,i,j,k,lb) = u2/S
                 unk(ii+2*nunkvbles,i,j,k,lb) = u3/S

                 unk(ii+(N_unk-2)*nunkvbles,i,j,k,lb) =delta  !Temperature = 0
                 ! solute
                 c_min = g_aa(2)
                 c_max = g_aa(3)
                 if(g_c0.gt.0)then
                    c_ave =g_c0
                 else
                    c_ave = 0.5*c_min + 0.5*c_max
                    g_c0 = c_ave
                 endif

                 unk(ii+(N_unk-1)*nunkvbles,i,j,k,lb) =C_ave+u3/S*(c_max-c_ave)+u2/S*(c_min-c_ave)
!                 unk(ii+(N_unk-1)*nunkvbles,i,j,k,lb) =C_ave+u3*(c_max-c_ave)+u2*(c_min-c_ave)
              elseif(g_model.eq.1)then
                  u3 = 0.5-0.5*tanh(steepness*(sqrt((xi)**2+(yi-0*nuc_radius*2d0)**2)-nuc_radius)) !AlNi
                  u3=min(1d0,max(0d0,u3))
                  u2=0d0
                  u1 = min(1d0,max(0d0,1d0 - u2 -u3))                                  !AlNi
                  S=u1+u2+u3
                  u1=u1/S
                  u2=u2/S
                  u3=u3/S

                  unk(ii+0*nunkvbles,i,j,k,lb)=u1
                  unk(ii+1*nunkvbles,i,j,k,lb)=u2
                  unk(ii+2*nunkvbles,i,j,k,lb)=u3

                  unk(ii+(N_unk-2)*nunkvbles,i,j,k,lb) =delta
                  unk(ii+(N_unk-1)*nunkvbles,i,j,k,lb) =g_c0 + (0.38-g_c0)*(0.5-0.5*tanh(steepness*(sqrt(xi**2+yi**2)-2d2)))
!                  unk(ii+(N_unk-1)*nunkvbles,i,j,k,lb) =g_c0
                  
              else
                 write(*,*)"i.c. for g_model not coded"
                 stop
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


double precision function angle_ic(x,phi2,angle)
use solution_parameters
double precision,intent (in):: x,phi2,angle
double precision ::phi1,xx




phi1=1.-phi2

angle_ic =1.
xx=x-int(x/angle)*angle
if(xx.lt.angle*phi2*0.5)angle_ic = 0. 
if(xx.gt.angle*(phi2*0.5+phi1))angle_ic = 0.
return


end function angle_ic




