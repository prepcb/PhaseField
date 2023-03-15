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
   double precision :: C_ave=0.72, C_ratio,lamda,c1,c2,c3
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
           !the ci s  

           ! the idea is to keep the solutions from the initial condition to use as subsequent guesses
           c1 = c_ave
           c2 = c_min
           c3 = c_max
           lamda = 0d0
           phi_(1)=u1
           phi_(2)=u2
           phi_(3)=u3
           call get_ci(c1,c2,c3,lamda,phi_,c_ave,delta)

           unk(total_vars*nunkvbles+1,i,j,k,lb)=c1
           unk(total_vars*nunkvbles+2,i,j,k,lb)=c2
           unk(total_vars*nunkvbles+3,i,j,k,lb)=c3
           unk(total_vars*nunkvbles+4,i,j,k,lb)=lamda
           
        enddo
     enddo
  enddo




  return

end subroutine app_initial_soln_blk


!iterative solver to find solute compositions, c1,c2,c3 and df/dc = lamda
subroutine get_ci(c1,c2,c3,lamda,phi,c,T)
implicit none
double precision, intent(in)::phi(3),c,T
double precision,intent(inout)::c1,c2,c3,lamda
double precision :: H(4,4),g1,g2,g3,f_phi,r(4),FreeEnergy,d(4),v
double precision :: q1,q2,q3,cc
integer :: i,j,k


 if(abs(1d0-phi(1)).lt.1d-20)then
    lamda = FreeEnergy(c,T,phi,0d0,0d0,6)
    c1 = c
    return
 elseif(abs(1d0-phi(2)).lt.1d-20)then
    lamda = FreeEnergy(c,T,phi,0d0,0d0,7)
    c2 = c
    return
 elseif(abs(1d0-phi(3)).lt.1d-20)then
    lamda = FreeEnergy(c,T,phi,0d0,0d0,8)
    c3 = c
    return
 endif

! interpolating functions
 g1 =  f_phi(phi,1,0,0)
 g2 =  f_phi(phi,2,0,0)
 g3 =  f_phi(phi,3,0,0)

!iterative loop - usually about 4
 do i=1,100
    v =g1*c1+g2*c2+g3*c3  != c

    ! system to solve
    r(1) = FreeEnergy(c1,T,phi,0d0,0d0,6) - lamda
    r(2) = FreeEnergy(c2,T,phi,0d0,0d0,7) - lamda
    r(3) = FreeEnergy(c3,T,phi,0d0,0d0,8) - lamda
    r(4) = c-v


       
    ! for jacobi entries - see H(i,j)
    q1 = FreeEnergy(c1,T,phi,0d0,0d0,9)
    q2 = FreeEnergy(c2,T,phi,0d0,0d0,10)
    q3 = FreeEnergy(c3,T,phi,0d0,0d0,11)


    !the inverse Jacobi - from maple

    H(1,1) = (g2 * q3 + g3 * q2) / (g1 * q2 * q3 + g2 * q1 * q3 + g3 * q1 * q2)
    H(1,2) = -g2 / (g1 * q2 * q3 + g2 * q1 * q3 + g3 * q1 * q2) * q3
    H(1,3) = -g3 / (g1 * q2 * q3 + g2 * q1 * q3 + g3 * q1 * q2) * q2
    H(1,4) = -q2 * q3 / (g1 * q2 * q3 + g2 * q1 * q3 + g3 * q1 *   q2)
    H(2,1) = -g1 / (g1 * q2 * q3 + g2 * q1 * q3 + g3 * q1 * q2) * q3
    H(2,2) = (g1 * q3 + g3 * q1) / (g1 * q2 * q3 + g2 * q1 *  q3 + g3 * q1 * q2)
    H(2,3) = -g3 / (g1 * q2 * q3 + g2 * q1 * q3 + g3 * q1 * q2) * q1
    H(2,4) = -q1 * q3 / (g1 * q2 * q3 + g2 * q1 * q3 + g3 * q1 *   q2) 
    H(3,1) = -g1 / (g1 * q2 * q3 + g2 * q1 * q3 + g3 * q1 * q2) * q2
    H(3,2) = -g2 / (g1 * q2 * q3 + g2 * q1 * q3 + g3 * q1 * q2) * q1
    H(3,3) = (g1 * q2 + g2 * q1) / (g1 * q2 * q3 + g2 * q1 * q3 + g3 * q1 * q2)
    H(3,4) = -q1 * q2 / (g1 * q2 * q3 + g2 * q1 * q3 + g3 * q1 *  q2)  
    H(4,1) = -g1 / (g1 * q2 * q3 + g2 * q1 * q3 + g3 * q1 * q2) * q2 * q3
    H(4,2) = -g2 / (g1 * q2 * q3 + g2 * q1 * q3 + g3 * q1 * q2) * q1 * q3
    H(4,3) = -g3 / (g1 * q2 * q3 + g2 * q1 * q3 + g3 * q1 * q2) * q1 * q2
    H(4,4) = -q1 * q2 * q3 / (g1 * q2 * q3 + g2 * q1 * q3 + g3 * q1 * q2)     


    !update defect d(4)
    d=0d0
    do j=1,4
       do k=1,4
          d(j) = d(j) + H(j,k)*r(k)
       enddo
    enddo

    ! tolerance - note if say phi(3)=0 then c3 is arbitrary
    if(abs(phi(1)*d(1))+abs(phi(2)*d(2))+abs(phi(3)*d(3)).lt.1d-8)exit ! tolerance reached


    ! improve the guesses
    c1 = max(1d-9,min(1d0-1d-9,c1 - d(1))) 
    c2 = max(1d-9,min(1d0-1d-9,c2 - d(2)))
    c3 = max(1d-9,min(1d0-1d-9,c3 - d(3)))
    lamda = max(-1d9,min(1d9,lamda - d(4)))

 enddo

 if(i.ge.200)then
    write(*,*)"maximum iterations reached"
    stop
 endif


return
end subroutine get_ci

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



