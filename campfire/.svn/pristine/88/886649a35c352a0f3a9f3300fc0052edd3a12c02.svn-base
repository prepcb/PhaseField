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
  integer :: nguard0, i, j, k,ip,ii,jj,i1,i2,ismooth
  double precision :: radius, pi, dx, dy, dz, xi, yi, zi,u1,u2,u3
  double precision :: steepness=1.0,snap21,eut_ic
   double precision,dimension(32) ::p =0
!   double precision :: C_min=0.252,C_ave=0.73,C_max=.98, C_ratio
   double precision :: C_ave=0.72, C_ratio
  double precision :: phi_(3),c_,FreeEnergy,FF(3,10000),df(10000),dg(10000),S,h,n=9999,dh


  if(init_c)then
  dh=1d0/(n+1d0)
  do ii=1,n
     do jj=1,3
        phi_= 0d0
        phi_(jj) = 1d0
        c_ = ii*dh
        FF(jj,ii)=FreeEnergy(c_,delta*0.5,phi_,0d0,0d0,0)
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
!  write(*,*)C_ratio,"hello"
  
  pi=4.0*datan(1.0)
  snap21=0.5+0.5*tanh(steepness*nuc_radius*scale_factor)
  
  nguard0 = nguard*npgs
  unk(:,:,:,:,lb) =0.
! set values for unk
  dx = bsize(1,lb)/real(nxb)
  do k=kl_bnd+nguard0*k3d,ku_bnd-nguard0*k3d
     zi = bnd_box(1,3,lb) + dx*(real(k-nguard0)-.5) - nucleate_z
     if (ndim.eq.2) zi = 0.0
     do j=jl_bnd+nguard0*k2d,ju_bnd-nguard0*k2d
        yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5) - nucleate_y
        do i=il_bnd+nguard0,iu_bnd-nguard0
           xi = bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5) - nucleate_x
           do ii=1,nunkvbles-1
           

!!!Eutectic

        u1 =1.0-(.5 -0.5*tanh(steepness*(xi-nuc_radius)*scale_factor))/snap21
        u3 =(1.0  -u1) *eut_ic(yi,c_ratio)
        u2=1.0 - u1 - u3
        unk(ii+0*nunkvbles,i,j,k,lb)=u1
        unk(ii+1*nunkvbles,i,j,k,lb)=u2
        unk(ii+2*nunkvbles,i,j,k,lb)=u3
! Temperature            
!            unk(ii+(N_phases)*nunkvbles,i,j,k,lb) = delta

           unk(ii+(N_phases)*nunkvbles,i,j,k,lb) =delta

! solute
            unk(ii+(N_phases+1)*nunkvbles,i,j,k,lb) =C_ave+min(1.0,exp(10.*(nuc_radius+0.0*sin(yi*.1)-xi)*scale_factor))*((c_max-c_min)*eut_ic(yi,c_ratio)+c_min-c_ave)
           enddo
        enddo
     enddo
  enddo
 
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



double precision function eut_ic(x,phi2)
use solution_parameters
double precision,intent (in):: x,phi2
double precision ::phi1,xx,q(12)=0
!double precision ::width=37.05242907,phi1,xx,q(12)=0



phi1=1.-phi2
eut_ic =1.
xx=x
do i=1,12
 if(xx.gt.width+q(i))then
   xx=xx-(width+q(i))
 else
  if(xx.lt.(width+q(i))*phi1*0.5)eut_ic = 0. 
  if(xx.gt.(width+q(i))*(phi1*0.5+phi2))eut_ic = 0.
  return
 end if
enddo
end function eut_ic



