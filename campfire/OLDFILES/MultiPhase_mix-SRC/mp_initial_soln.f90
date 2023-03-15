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
  integer :: nguard0, i, j, k,ip,ii,jj,i1,i2,ismooth,iter
  double precision :: radius, pi, dx, dy, dz, xi, yi, zi,u1,u2,u3
  double precision :: steepness=1.0,snap21,eut_ic
   double precision,dimension(32) ::p =0
!   double precision :: C_min=0.252,C_ave=0.73,C_max=.98, C_ratio
   double precision :: C_ave=0.72, C_ratio
  double precision :: phi_(3),c_,FreeEnergy,FF(3,10000),df(10000),dg(10000),S,h,n=9999,dh
  double precision :: w(4)
!!!!!!!!!!!!!!!




!!!!!!!!!!!!!

  if(init_c)then
     w(1)=-1d0
  call get_mix(delta,2,3,w,iter)


  c_min=w(1)
  c_max=w(2)

  c23(1)=w(1)
  c23(2)=w(2)

  c23(3)=w(4) !w(4) is free energy at c23(1)
  lamda(2,3)=w(3)
  lamda(3,2)=w(3)
  cc23(1)=w(1)
  cc23(2)=w(2)
  cc23(3)=w(3)
  call get_mix(delta,1,3,w,iter)
  c31(1)=w(1)
  c31(2)=w(2)
  c31(3)=w(4)
  cc31(1)=w(1)
  cc31(2)=w(2)
  cc31(3)=w(3)
  lamda(1,3)=w(3)
  lamda(3,1)=w(3)
  call get_mix(delta,2,1,w,iter)
  c12(1)=w(1)
  c12(2)=w(2)
  c12(3)=w(4)
  cc12(1)=w(1)
  cc12(2)=w(2)
  cc12(3)=w(3)
  lamda(1,2)=w(3)
  lamda(2,1)=w(3)
  init_c=.false.
  endif



  C_ratio=1.-(c_max-c_ave)/(C_max-C_min)



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

subroutine get_mix(T,p1,p2,w,iter)
use solution_parameters
implicit none
double precision, intent (in)::T
double precision, intent (inout):: w(4)
integer, intent (in):: p1,p2
integer, intent (out)::iter
double precision :: FreeEnergy,phi(3),e1,e2,EntropyMixing,JJ(2,2),a,b,c,d,e,f,g,detJ, invJ(2,2),tol,cc(2)
integer :: i,j,max_iter=40
double precision :: ha,hb,hr,V(8),RK(5),T_k,Tk,R=8.31,T_eut=456.14


if(w(1).lt.0d0)then
   w(1)=0.01
   w(2)=0.99
endif

do i = 1,max_iter

   !e1 e2 are the equations to solve
   !1st guess w(1) and w(2)
   phi=0d0
   phi(p1)=1d0
   a = FreeEnergy(w(1),T,phi,0d0,0d0,5)!+(1d0+T)*EntropyMixing(w(1),1)
   c = FreeEnergy(w(1),T,phi,0d0,0d0,0)+(1d0+T)*EntropyMixing(w(1),0)
   e = FreeEnergy(w(1),T,phi,0d0,0d0,6)+(1d0+T)*EntropyMixing(w(1),2)
   phi=0d0
   phi(p2)=1d0
   b =  FreeEnergy(w(2),T,phi,0d0,0d0,5)!+ (1d0+T)*EntropyMixing(w(2),1)
   d = FreeEnergy(w(2),T,phi,0d0,0d0,0) + (1d0+T)*EntropyMixing(w(2),0)
   f = FreeEnergy(w(2),T,phi,0d0,0d0,6) + (1d0+T)*EntropyMixing(w(2),2)
   g = (w(2)- w(1))
   
   !The Eqs to solve with initial guesses
   e1 = a - b
   e2 = (d - c)/g - a 
   



   !jacobi matrix
   JJ(1,1) = e                          !d e1/d w(1)
   JJ(1,2) = -f                         !d e1/d w(2)
   JJ(2,1) = -a/g + (d-c)/g**2 - e      !d e2/ d w(1)
   JJ(2,2) =  b/g - (d-c)/g**2          !d e2/ d w(2)
   

   detJ = JJ(1,1)*JJ(2,2) - JJ(1,2)*JJ(2,1)
   
   invJ(1,1) = JJ(2,2)
   invJ(2,2) = JJ(1,1)
   invJ(1,2) =-JJ(1,2)
   invJ(2,1) =-JJ(2,1)
   invJ = invJ/detJ



   
   cc(1) = invJ(1,1)*e1+invJ(1,2)*e2
   cc(2) = invJ(2,1)*e1+invJ(2,2)*e2
   w(1) = min(0.999,max(0.001,w(1)-cc(1)))
   w(2) = min(0.999,max(0.001,w(2)-cc(2)))
   w(3) = a
   w(4) = c

   tol = abs(cc(1))+abs(cc(2))
   if(tol.lt.1d-5)exit
   if(i.eq.max_iter)then
      write(*,*)"iterations for tangents equals max = ",i
      stop
   endif
      
enddo
iter=i

return
end subroutine get_mix

