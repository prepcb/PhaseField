!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  Utilities for phase field stencils in 2-d and 3-d
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! 168 
integer function imaximum(p,n)
implicit none
double precision, intent (in):: p(100)
integer, intent (in) :: n
integer i,j
double precision maximum
maximum = p(1)
j=1
do i = 2, n
if(p(i).ge.p(j))then
j=i
endif
enddo
imaximum=j

end function imaximum

double precision function maximum(p,n)
implicit none
double precision, intent (in):: p(100)
integer, intent (in) :: n
integer i,j,k
double precision:: q(100),amaximum

maximum=amaximum(p,n)
return
maximum = p(1)
j=1
do i = 2, n
   if(p(i).ge.p(j))then
      j=i
   endif
enddo
maximum=p(j)

end function maximum

double precision function amaximum(p,n)
implicit none
double precision, intent (in):: p(100)
integer, intent (in) :: n
integer i,j
amaximum = 0d0

do i = 1, n
   if(p(i).gt.0d0)amaximum = amaximum + p(i)**64
enddo
amaximum=amaximum**(0.01562500000)

end function amaximum

double precision function SurfaceEnergy(p,x,n)
use solution_parameters
implicit none
double precision, intent (in):: p(3,100),x(3)
integer, intent (in) :: n
integer i,j
double precision :: maximum,q(100),qmax

q=0d0
do j=1,n
   do i=1,3
      q(j) = q(j) + p(i,j)*x(i)
   enddo
enddo

qmax=0.04+maximum(q,n)

SurfaceEnergy =0.5*qmax**2
end function SurfaceEnergy


double precision function iSurfaceEnergy(p,x,n)
implicit none
double precision, intent (in):: p(3,100),x(3)
integer, intent (in) :: n
integer i,j
double precision :: imaximum,q(12),qmax

q=0d0
do j=1,n
   do i=1,3
      q(j) = q(j) + p(i,j)*x(i)
   enddo
enddo

iSurfaceEnergy =imaximum(q,n)
end function iSurfaceEnergy



double precision function GradientEnergyHexPrismHess(vbles,lp,dx)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3,3)
integer, intent (in):: lp
double precision :: phi_(3),phi(3,3),phi_tmp(3,3),A_(3),A(3,3),AA(3,3),fgf
double precision :: yy(3),y,xx(3),x,px,LapS,H,Lap,d=4d-2,c(6)=(/1,2,3,1,2,3/),rr(3)=(/1,1,2/)
double precision :: p(3,6 ),pi=3.141592654,q(12),dh(3),pi6=.5235987758
integer :: i,j,np=6,k,l
double precision :: SurfaceEnergy
double precision :: ID(3,3)=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))

do i=-1,1
p(1,i+2)=cos(pi*i/3d0+pi6)
p(2,i+2)=sin(pi*i/3d0+pi6)
p(3,i+2)=1d0
p(1,i+5)=cos(pi*i/3d0+pi6)
p(2,i+5)=sin(pi*i/3d0+pi6)
p(3,i+5)=-1d0
enddo



GradientEnergyHexPrismHess=0
   phi_(1)=0.5*(vbles(1,3,2,2)-vbles(1,1,2,2))/dx
   phi_(2)=0.5*(vbles(1,2,3,2)-vbles(1,2,1,2))/dx
   phi_(3)=0.5*(vbles(1,2,2,3)-vbles(1,2,2,1))/dx
!!$   phi(1,1)=(vbles(1,3,2,2)-2d0*vbles(1,2,2,2)+vbles(1,1,2,2))/dx**2
!!$   phi(2,2)=(vbles(1,2,3,2)-2d0*vbles(1,2,2,2)+vbles(1,2,1,2))/dx**2
!!$   phi(3,3)=(vbles(1,2,2,3)-2d0*vbles(1,2,2,2)+vbles(1,2,2,1))/dx**2
!!$
!!$   phi(1,2) = 0.25*(vbles(1,3,3,2)+vbles(1,1,1,2)-vbles(1,3,1,2)-vbles(1,1,3,2))/dx**2
!!$   phi(1,3) = 0.25*(vbles(1,3,2,3)+vbles(1,1,2,1)-vbles(1,3,2,1)-vbles(1,1,2,3))/dx**2
!!$   phi(2,3) = 0.25*(vbles(1,2,3,3)+vbles(1,2,1,1)-vbles(1,2,3,1)-vbles(1,2,1,3))/dx**2
!!$   phi(2,1) = phi(1,2)
!!$   phi(3,1) = phi(1,3)
!!$   phi(3,2) = phi(2,3)


   phi(1,1)=0.5*(vbles(1,3,2,2)-2d0*vbles(1,2,2,2)+vbles(1,1,2,2))
   phi(1,1) = phi(1,1) +0.25*(vbles(1,3,1,1)-2d0*vbles(1,2,1,1)+vbles(1,1,1,1))
   phi(1,1) = phi(1,1) +0.25*(vbles(1,3,3,3)-2d0*vbles(1,2,3,3)+vbles(1,1,3,3))

   phi(2,2)=0.5*(vbles(1,2,3,2)-2d0*vbles(1,2,2,2)+vbles(1,2,1,2))
   phi(2,2) = phi(2,2) +0.25*(vbles(1,1,3,1)-2d0*vbles(1,1,2,1)+vbles(1,1,1,1))
   phi(2,2) = phi(2,2) +0.25*(vbles(1,3,3,3)-2d0*vbles(1,3,2,3)+vbles(1,3,1,3))

   phi(3,3)=0.5*(vbles(1,2,2,3)-2d0*vbles(1,2,2,2)+vbles(1,2,2,1))
   phi(3,3) = phi(3,3) +0.25*(vbles(1,1,1,3)-2d0*vbles(1,1,1,2)+vbles(1,1,1,1))
   phi(3,3) = phi(3,3) +0.25*(vbles(1,3,3,3)-2d0*vbles(1,3,3,2)+vbles(1,3,3,1))
   

   phi(1,2) = 0.5*0.25*(vbles(1,3,3,2)+vbles(1,1,1,2)-vbles(1,3,1,2)-vbles(1,1,3,2))
   phi(1,2) = phi(1,2) + 0.25*0.25*(vbles(1,3,3,1)+vbles(1,1,1,1)-vbles(1,3,1,1)-vbles(1,1,3,1))
   phi(1,2) = phi(1,2) + 0.25*0.25*(vbles(1,3,3,3)+vbles(1,1,1,3)-vbles(1,3,1,3)-vbles(1,1,3,3))



   phi(1,3) = 0.5*0.25*(vbles(1,3,2,3)+vbles(1,1,2,1)-vbles(1,3,2,1)-vbles(1,1,2,3))
   phi(1,3) = phi(1,3) + 0.25*0.25*(vbles(1,3,1,3)+vbles(1,1,1,1)-vbles(1,3,1,1)-vbles(1,1,1,3))
   phi(1,3) = phi(1,3) + 0.25*0.25*(vbles(1,3,3,3)+vbles(1,1,3,1)-vbles(1,3,3,1)-vbles(1,1,3,3))
   phi(2,3) = 0.5*0.25*(vbles(1,2,3,3)+vbles(1,2,1,1)-vbles(1,2,3,1)-vbles(1,2,1,3))
   phi(2,3) = phi(2,3) + 0.25*0.25*(vbles(1,1,3,3)+vbles(1,1,1,1)-vbles(1,1,3,1)-vbles(1,1,1,3))
   phi(2,3) = phi(2,3) + 0.25*0.25*(vbles(1,3,3,3)+vbles(1,3,1,1)-vbles(1,3,3,1)-vbles(1,3,1,3))
   phi(2,1) = phi(1,2)
   phi(3,1) = phi(1,3)
   phi(3,2) = phi(2,3)

   phi = phi/dx**2




   do i=1,3
      xx(i) = phi_(i)
   enddo


   x = sqrt(xx(1)**2+xx(2)**2+xx(3)**2+1d-9)



!   d=0.1/min_dx*4d-2
   d=max(x*4d-2,1d-5)


   do i=1,3
      dh=0d0
      dh(i)=d
      yy=xx-dh
      A(i,i) = SurfaceEnergy(p,yy,np)
      A(i,i) = A(i,i) -2d0*SurfaceEnergy(p,xx,np)
      yy=xx+dh
      A(i,i) = A(i,i) + SurfaceEnergy(p,yy,np)
      A(i,i) = A(i,i)/d**2
   enddo

   do i=1,2
      do j=i+1,3

         dh=0d0
         dh(i)=d
         dh(j)=d
         yy=xx + dh
         A(i,j) = SurfaceEnergy(p,yy,np)

         yy=xx - dh
         A(i,j) = A(i,j) + SurfaceEnergy(p,yy,np)

         dh=0d0
         dh(i)=d
         dh(j)=-d
         yy=xx+dh
         A(i,j) = A(i,j) - SurfaceEnergy(p,yy,np)

         yy=xx-dh
         A(i,j) = A(i,j) - SurfaceEnergy(p,yy,np)

         A(i,j)=A(i,j)*0.25/d**2
         A(j,i)=A(i,j)
      enddo
   enddo

if(.false.)then
   xx =xx/x

   AA = A
   phi_tmp = phi
   do i=1,3
      do j=1,3
         do k=1,3
            phi(i,j) = phi(i,j)-phi_tmp(i,k)*xx(k)*xx(j)-phi_tmp(j,k)*xx(k)*xx(i)
            A(i,j) = A(i,j)-AA(i,k)*xx(k)*xx(j)-AA(j,k)*xx(k)*xx(i)
            do l=1,3
               A(i,j) = A(i,j)+2*xx(i)*xx(k)*AA(k,l)*xx(l)*xx(j)
               phi(i,j) = phi(i,j)+2*xx(i)*xx(k)*phi_tmp(k,l)*xx(l)*xx(j)
            enddo
         enddo
      enddo
   enddo
endif



   H = 0d0
   do i=1,3
      do j=1,3
         H = H + phi(i,j)*A(i,j)
      enddo
   enddo

   GradientEnergyHexPrismHess = - H





end function GradientEnergyHexPrismHess

double precision function GradientEnergyHexPrism(vbles,lp,dx)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3,3)
integer, intent (in):: lp
double precision :: phi_(3),phi(3,3),phi_tmp(3,3),A_(3),A(3,3),AA(3,3),fgf
double precision :: yy(3),y,xx(3),x,px,LapS,H,Lap,d=4d-2,c(6)=(/1,2,3,1,2,3/),rr(3)=(/1,1,2/)
double precision :: p(3,12),pi=3.141592654,q(12),dh(3),pi6=.5235987758,t1(3),t2(3),t3(3)
integer :: i,j,np=12,k,l,iSurfaceEnergy
double precision :: SurfaceEnergy,maxA
double precision :: ID(3,3)=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))

do i=1,6
p(2,i)=cos(pi*i/3d0+pi6)
p(3,i)=sin(pi*i/3d0+pi6)
p(1,i)=1d0
p(2,i+6)=cos(pi*i/3d0+pi6)
p(3,i+6)=sin(pi*i/3d0+pi6)
p(1,i+6)=-1d0
enddo

p=p/sqrt(2.)


GradientEnergyHexPrism=0
   phi_(1)=0.5*(vbles(1,3,2,2)-vbles(1,1,2,2))/dx
   phi_(2)=0.5*(vbles(1,2,3,2)-vbles(1,2,1,2))/dx
   phi_(3)=0.5*(vbles(1,2,2,3)-vbles(1,2,2,1))/dx
   phi(1,1)=(vbles(1,3,2,2)-2d0*vbles(1,2,2,2)+vbles(1,1,2,2))/dx**2
   phi(2,2)=(vbles(1,2,3,2)-2d0*vbles(1,2,2,2)+vbles(1,2,1,2))/dx**2
   phi(3,3)=(vbles(1,2,2,3)-2d0*vbles(1,2,2,2)+vbles(1,2,2,1))/dx**2

   phi(1,2) = 0.25*(vbles(1,3,3,2)+vbles(1,1,1,2)-vbles(1,3,1,2)-vbles(1,1,3,2))/dx**2
   phi(1,3) = 0.25*(vbles(1,3,2,3)+vbles(1,1,2,1)-vbles(1,3,2,1)-vbles(1,1,2,3))/dx**2
   phi(2,3) = 0.25*(vbles(1,2,3,3)+vbles(1,2,1,1)-vbles(1,2,3,1)-vbles(1,2,1,3))/dx**2
   phi(2,1) = phi(1,2)
   phi(3,1) = phi(1,3)
   phi(3,2) = phi(2,3)
   do i=1,3
      xx(i) = phi_(i)
   enddo


   x = sqrt(xx(1)**2+xx(2)**2+xx(3)**2+1d-9)

   xx=xx/x


  d=1d-4


   do i=1,3
      dh=0d0
      dh(i)=d
      yy=xx-dh
      A(i,i) = SurfaceEnergy(p,yy,np)
      A(i,i) = A(i,i) -2d0*SurfaceEnergy(p,xx,np)
      yy=xx+dh
      A(i,i) = A(i,i) + SurfaceEnergy(p,yy,np)
      A(i,i) = A(i,i)/d**2
   enddo

   do i=1,2
      do j=i+1,3

         dh=0d0
         dh(i)=d
         dh(j)=d
         yy=xx + dh
         A(i,j) = SurfaceEnergy(p,yy,np)

         yy=xx - dh
         A(i,j) = A(i,j) + SurfaceEnergy(p,yy,np)

         dh=0d0
         dh(i)=d
         dh(j)=-d
         yy=xx+dh
         A(i,j) = A(i,j) - SurfaceEnergy(p,yy,np)

         yy=xx-dh
         A(i,j) = A(i,j) - SurfaceEnergy(p,yy,np)

         A(i,j)=A(i,j)*0.25/d**2
         A(j,i)=A(i,j)
      enddo
   enddo


if(.false.)then
   AA = A
   phi_tmp = phi
   do i=1,3
      do j=1,3
         do k=1,3
            phi(i,j) = phi(i,j)-phi_tmp(i,k)*xx(k)*xx(j)-phi_tmp(j,k)*xx(k)*xx(i)
            A(i,j) = A(i,j)-AA(i,k)*xx(k)*xx(j)-AA(j,k)*xx(k)*xx(i)
            do l=1,3
               A(i,j) = A(i,j)+2*xx(i)*xx(k)*AA(k,l)*xx(l)*xx(j)
               phi(i,j) = phi(i,j)+2*xx(i)*xx(k)*phi_tmp(k,l)*xx(l)*xx(j)
            enddo
         enddo
      enddo
   enddo
endif



   

   H = 0d0
   do i=1,3
      do j=1,3
         H = H + phi(i,j)*A(i,j) 
      enddo
   enddo


   GradientEnergyHexPrism = - H





end function GradientEnergyHexPrism

double precision function GradientEnergyPCB3(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3,3)
integer, intent (in):: lp
double precision :: phi_(3),phi(3,3),A_(3),A(3,3),AA,fgf
double precision :: xx(3),x,px,LapS,H,Lap,d=4d-2,c(6)=(/1,2,3,1,2,3/),rr(3)=(/1,1,2/)
integer :: i,j

GradientEnergyPCB3=0
   phi_(1)=0.5*(vbles(1,3,2,2)-vbles(1,1,2,2))/dx
   phi_(2)=0.5*(vbles(1,2,3,2)-vbles(1,2,1,2))/dx
   phi_(3)=0.5*(vbles(1,2,2,3)-vbles(1,2,2,1))/dx
   phi(1,1)=(vbles(1,3,2,2)-2d0*vbles(1,2,2,2)+vbles(1,1,2,2))/dx**2
   phi(2,2)=(vbles(1,2,3,2)-2d0*vbles(1,2,2,2)+vbles(1,2,1,2))/dx**2
   phi(3,3)=(vbles(1,2,2,3)-2d0*vbles(1,2,2,2)+vbles(1,2,2,1))/dx**2

   phi(1,2) = 0.25*(vbles(1,3,3,2)+vbles(1,1,1,2)-vbles(1,3,1,2)-vbles(1,1,3,2))/dx**2
   phi(1,3) = 0.25*(vbles(1,3,2,3)+vbles(1,1,2,1)-vbles(1,3,2,1)-vbles(1,1,2,3))/dx**2
   phi(2,3) = 0.25*(vbles(1,2,3,3)+vbles(1,2,1,1)-vbles(1,2,3,1)-vbles(1,2,1,3))/dx**2
   phi(2,1) = phi(1,2)
   phi(3,1) = phi(1,3)
   phi(3,2) = phi(2,3)
   x = sqrt(phi_(1)**2 +phi_(2)**2 + phi_(3)**2 + 1d-9)
   Lap = 0d0

   do i=1,3
      do j=1,3
         A(i,j)=rr(i)*rr(j)
      enddo
   enddo
   do i=1,3
      xx(i) = phi_(i)/x
   enddo
   fgf=0d0
   do i=1,3
      fgf = fgf + xx(i)*rr(i)
   enddo
   do i=1,3
      px = abs(xx(i))
      if(px.lt.d)then
         A(i,i)=A(i,i)+fgf*rr(i)/d
!         A(i,i)=sqrt(2d0)/d*(2d0*(d-px)/d)

!quatic approx
!         A(i,i)=15./16.*(1-px/d)**2*(1+px/d)**2/d
!         A(i,i) = A(i,i)*sqrt(2d0)

         !         A(i,i)= A(i,i) +rr(i)*(rr(c(i+1))*xx(c(i+1))+rr(c(i+2))*xx(c(i+2)))/d
         !         A(i,i)= A(i,i) +rr(i)*(rr(c(i+1))*xx(c(i+1))+rr(c(i+2))*xx(c(i+2)))/d&
         !              *2d0*(1d0-px/d)
         !         A(i,i)=2d0*(d-px)/d**2                                                                                                                        

!         A(i,c(i+1))=0d0
!         A(i,c(i+2))=0d0
!         A(c(i+1),i)=0d0
!         A(c(i+2),i)=0d0

      endif
   enddo


   H = 0d0
   do i=1,3
      do j=1,3
         H = H + phi(i,j)*A(i,j)
      enddo
   enddo

   GradientEnergyPCB3 = - H




end function GradientEnergyPCB3




double precision function GradientEnergyPCB2(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3,3)
integer, intent (in):: lp
double precision :: phi_(3),phi(3,3),A_(3),A(3,3),AA,XA,XXAA,phiRR
double precision :: xx(3),x,LapS,H,Lap
integer :: i,j

GradientEnergyPCB2=0
   phi_(1)=0.5*(vbles(1,3,2,2)-vbles(1,1,2,2))/dx
   phi_(2)=0.5*(vbles(1,2,3,2)-vbles(1,2,1,2))/dx
   phi_(3)=0.5*(vbles(1,2,2,3)-vbles(1,2,2,1))/dx
   phi(1,1)=(vbles(1,3,2,2)-2d0*vbles(1,2,2,2)+vbles(1,1,2,2))/dx**2
   phi(2,2)=(vbles(1,2,3,2)-2d0*vbles(1,2,2,2)+vbles(1,2,1,2))/dx**2
   phi(3,3)=(vbles(1,2,2,3)-2d0*vbles(1,2,2,2)+vbles(1,2,2,1))/dx**2

   phi(1,2) = 0.25*(vbles(1,3,3,2)+vbles(1,1,1,2)-vbles(1,3,1,2)-vbles(1,1,3,2))/dx**2
   phi(1,3) = 0.25*(vbles(1,3,2,3)+vbles(1,1,2,1)-vbles(1,3,2,1)-vbles(1,1,2,3))/dx**2
   phi(2,3) = 0.25*(vbles(1,2,3,3)+vbles(1,2,1,1)-vbles(1,2,3,1)-vbles(1,2,1,3))/dx**2
   phi(2,1) = phi(1,2)
   phi(3,1) = phi(1,3)
   phi(3,2) = phi(2,3)
   x = sqrt(phi_(1)**2 +phi_(2)**2 + phi_(3)**2 + 1d-9)
   Lap = 0d0
   do i=1,3
      xx(i) = phi_(i)/x
      Lap = Lap + phi(i,i)
   enddo
   AA = 1.-3.*epsilon + 4.*epsilon*(xx(1)**4 + xx(2)**4 + xx(3)**4)
   A = 0d0
   XA=0d0
   do i = 1,3
      A_(i) = 16.*epsilon*xx(i)**3
      A(i,i)= 48.*epsilon*xx(i)**2
      XA=XA+2d0*AA*xx(i)*A_(i)
   enddo
   H = 0d0
   XXAA=0d0
   phiRR=0d0
   do i=1,3
      do j=1,3
         H = H + 2d0*phi(i,j)*(A_(i)*A_(j)+AA*A(i,j))
         XXAA=XXAA + 2d0*(A_(i)*A_(j)+AA*A(i,j))*xx(i)*XX(j)
         phiRR=phiRR + phi(i,j)*xx(i)*xx(j)
      enddo
   enddo
   
   GradientEnergyPCB2 =  (AA**2-0.5*XA)*Lap + 0.5*H + 0.5*(XA - XXAA)*phiRR
   GradientEnergyPCB2 =  -   GradientEnergyPCB2 




      
end function GradientEnergyPCB2

double precision function GradientEnergyPCB(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3,3)
integer, intent (in):: lp
double precision :: phi_(3),phi(3,3),A_(3),A(3,3),AA
double precision :: xx(3),x,LapS,H,Lap
integer :: i,j

GradientEnergyPCB=0
   phi_(1)=0.5*(vbles(1,3,2,2)-vbles(1,1,2,2))/dx
   phi_(2)=0.5*(vbles(1,2,3,2)-vbles(1,2,1,2))/dx
   phi_(3)=0.5*(vbles(1,2,2,3)-vbles(1,2,2,1))/dx
   phi(1,1)=(vbles(1,3,2,2)-2d0*vbles(1,2,2,2)+vbles(1,1,2,2))/dx**2
   phi(2,2)=(vbles(1,2,3,2)-2d0*vbles(1,2,2,2)+vbles(1,2,1,2))/dx**2
   phi(3,3)=(vbles(1,2,2,3)-2d0*vbles(1,2,2,2)+vbles(1,2,2,1))/dx**2

   phi(1,2) = 0.25*(vbles(1,3,3,2)+vbles(1,1,1,2)-vbles(1,3,1,2)-vbles(1,1,3,2))/dx**2
   phi(1,3) = 0.25*(vbles(1,3,2,3)+vbles(1,1,2,1)-vbles(1,3,2,1)-vbles(1,1,2,3))/dx**2
   phi(2,3) = 0.25*(vbles(1,2,3,3)+vbles(1,2,1,1)-vbles(1,2,3,1)-vbles(1,2,1,3))/dx**2
   phi(2,1) = phi(1,2)
   phi(3,1) = phi(1,3)
   phi(3,2) = phi(2,3)
   x = sqrt(phi_(1)**2 +phi_(2)**2 + phi_(3)**2 + 1d-9)
   Lap = 0d0
   do i=1,3
      xx(i) = phi_(i)/x
      Lap = Lap + phi(i,i)
   enddo
   AA = 1.-3.*epsilon + 4.*epsilon*(xx(1)**4 + xx(2)**4 + xx(3)**4)
   do i = 1,3
      A_(i) = 16.*epsilon*xx(i)**3
      A(i,i)= 48.*epsilon*xx(i)**2
   enddo
   LapS = 0d0
   H = 0d0
   do i=1,3
      LapS = LapS + 2.*(A_(i)**2+AA*A(i,i))-4.*xx(i)*AA*A_(i)
      do j=1,3
         LapS = LapS - 2*xx(i)*xx(j)*A_(i)*A_(j)
         H = H + phi(i,j)*xx(i)*xx(j)
      enddo
   enddo
   H = Lap - H
   
   GradientEnergyPCB = - AA**2*Lap - 0.5 * LapS*H



      
end function GradientEnergyPCB

double precision function GradientEnergy(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3,3)
integer, intent (in):: lp
!double precision ::phix,phiy,X,Y,u,v,A
!double precision :: G11,G22,G12,M11,M12,g,GG
double precision :: LapPhi
double precision :: A,B,phix,phiy,phiz,phi11,phi22,phi33,phi12,phi13,phi23,G11,G22,G33,G12,G13,G23,B11,B22,B33,B12,B13,B23,B1,B2,B3,d
!call Calc_GradPhi(LapPhi,vbles,dx)
!!$
!!$phix = 0.5*(vbles(1,3,2) - vbles(1,1,2))/dx 
!!$phiy = 0.5*(vbles(1,2,3) - vbles(1,2,1))/dx 
!!$u = phix*phix
!!$v = phiy*phiy
!!$g  = u + v
!!$if(g.gt.1d-20)then
!!$
!!$   X = u/g
!!$   Y= v/g
!!$   A = A_0*(1d0+epsilon_tilde*(X**2+Y**2))
!!$
!!$   GG = 16*A_0*epsilon_tilde*X*Y*(A_0*epsilon_tilde*(X-Y)**2+A)
!!$   G11 = Y*(GG +4*epsilon_tilde*A*A_0*(X-Y))+ A**2
!!$   G22 = X*(GG +4*epsilon_tilde*A*A_0*(Y-X))+ A**2
!!$   G12 = -phix*phiy*GG/g
!!$
!!$
!!$
!!$   M11 = 0.5*(vbles(1,3,2)+vbles(1,1,2)-vbles(1,2,3)-vbles(1,2,1))/dx**2  !(phi_xx - phi_yy)/2
!!$   M12 = 0.25*(vbles(1,3,3)+vbles(1,1,1)-vbles(1,1,3)-vbles(1,3,1))/dx**2 !phi_xy
!!$
!!$   GradientEnergy = 0.5*LapPhi(1)*(G11+G22)+M11*(G11-G22)+2*M12*G12
!!$
!!$
!!$else
!   GradientEnergy = (A_0*(1d0+epsilon_tilde))**2*LapPhi

   d=1d-2/g_lambda

   phix=0.5*(vbles(1,3,2,2)-vbles(1,1,2,2))/dx
   phiy=0.5*(vbles(1,2,3,2)-vbles(1,2,1,2))/dx
   phiz=0.5*(vbles(1,2,2,3)-vbles(1,2,2,1))/dx
   phi11=(vbles(1,3,2,2)-2d0*vbles(1,2,2,2)+vbles(1,1,2,2))/dx**2
   phi22=(vbles(1,2,3,2)-2d0*vbles(1,2,2,2)+vbles(1,2,1,2))/dx**2
   phi33=(vbles(1,2,2,3)-2d0*vbles(1,2,2,2)+vbles(1,2,2,1))/dx**2

   phi12 = 0.25*(vbles(1,3,3,2)+vbles(1,1,1,2)-vbles(1,3,1,2)-vbles(1,1,3,2))/dx**2
   phi13 = 0.25*(vbles(1,3,2,3)+vbles(1,1,2,1)-vbles(1,3,2,1)-vbles(1,1,2,3))/dx**2
   phi23 = 0.25*(vbles(1,2,3,3)+vbles(1,2,1,1)-vbles(1,2,3,1)-vbles(1,2,1,3))/dx**2
   
   A = 0.5*(d**2+phix**2+phiy**2+phiz**2)
   B = 4*(phix**4+phiy**4+phiz**4)/A**2
   !!==============
   B1 = 16*phix**3/A**2-8*(phix**4+phiy**4+phiz**4)/A**3*phix
   B2 = 16*phiy**3/A**2-8*(phix**4+phiy**4+phiz**4)/A**3*phiy
   B3 = 16*phiz**3/A**2-8*(phix**4+phiy**4+phiz**4)/A**3*phiz
   B11 = 48*phix**2/A**2-64*phix**4/A**3+24*(phix**4+phiy**4+phiz**4)/A**4*phix**2-8*(phix**4+phiy**4+phiz**4)/A**3
   B22 = 48*phiy**2/A**2-64*phiy**4/A**3+24*(phix**4+phiy**4+phiz**4)/A**4*phiy**2-8*(phix**4+phiy**4+phiz**4)/A**3
   B33 = 48*phiz**2/A**2-64*phiz**4/A**3+24*(phix**4+phiy**4+phiz**4)/A**4*phiz**2-8*(phix**4+phiy**4+phiz**4)/A**3
   B12 = -32*phix**3/A**3*phiy-32*phiy**3/A**3*phix+24*(phix**4+phiy**4+phiz**4)/A**4*phix*phiy
   B13 = -32*phix**3/A**3*phiz-32*phiz**3/A**3*phix+24*(phix**4+phiy**4+phiz**4)/A**4*phix*phiz
   B23 = -32*phiy**3/A**3*phiz-32*phiz**3/A**3*phiy+24*(phix**4+phiy**4+phiz**4)/A**4*phiy*phiz
   !!==============
   G11 = (epsilon*B-3*epsilon+1)**2+4*phix*(epsilon*B-3*epsilon+1)*epsilon*B1+2*A*epsilon**2*B1**2+2*A*(epsilon*B-3*epsilon+1)*epsilon*B11
   G22 = (epsilon*B-3*epsilon+1)**2+4*phiy*(epsilon*B-3*epsilon+1)*epsilon*B2+2*A*epsilon**2*B2**2+2*A*(epsilon*B-3*epsilon+1)*epsilon*B22
   G33 = (epsilon*B-3*epsilon+1)**2+4*phiz*(epsilon*B-3*epsilon+1)*epsilon*B3+2*A*epsilon**2*B3**2+2*A*(epsilon*B-3*epsilon+1)*epsilon*B33
   G12 = 2*phix*(epsilon*B-3*epsilon+1)*epsilon*B2+2*phiy*(epsilon*B-3*epsilon+1)*epsilon*B1+2*A*epsilon**2*B2*B1+2*A*(epsilon*B-3*epsilon+1)*epsilon*B12
   G13 = 2*phix*(epsilon*B-3*epsilon+1)*epsilon*B3+2*phiz*(epsilon*B-3*epsilon+1)*epsilon*B1+2*A*epsilon**2*B3*B1+2*A*(epsilon*B-3*epsilon+1)*epsilon*B13
   G23 = 2*phiy*(epsilon*B-3*epsilon+1)*epsilon*B3+2*phiz*(epsilon*B-3*epsilon+1)*epsilon*B2+2*A*epsilon**2*B3*B2+2*A*(epsilon*B-3*epsilon+1)*epsilon*B23
   

GradientEnergy = phi11*G11+phi22*G22+phi33*G33+2d0*(phi12*G12+phi13*G13+phi23*G23)

GradientEnergy = - GradientEnergy

!
end function GradientEnergy








!T_K converts non-dimensional parameter T to dimensioned T_K(T)
double precision function T_K(T)
use solution_parameters
implicit none
double precision, intent (in)::T


 
 T_K = g_T0*(1.+T)

 end function T_K
 double precision function TemperatureRHS(vbles,dx,c,T,phi,c_dot,phi_dot)
use solution_parameters
implicit none
double precision, intent (in):: dx,c,T,phi,c_dot,phi_dot,vbles(N_unk+1,3,3,3)

double precision, dimension(3,3)::stencil =reshape((/0.,1.,0.,1.,-4.,1.,0.,1.,0./),(/3,3/))
double precision, dimension(3,3,3)::stencil3d =reshape((/0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 1., -6., 1., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0./),(/3,3,3/))
double precision, dimension(3)::stencil_x =(/-.5,.0,.5/)
double precision:: FreeEnergy,GradientEnergy,potential,x,a,xa,dFdc(3,3,3),d_H(3,3,3) !x coordinate
double precision:: Cp= 28.5!29.9758394     !J/mol/K
double precision :: D_heat_c 
double precision :: T_,c_,phi_,MolPerVol,X_energy,Y_energy,div,Grad_DW
integer :: ii,jj,kk
if(.not.thermal)then
   write(*,*)"error in TemperatureRHS: not thermal"
   stop
endif



  MolPerVol = g_MolPerVol(1)*(1.-c)+g_MolPerVol(2)*c  !mol/m^3



! to help with floating point arithmetic
  X_energy=g_R*g_T0  !J/mol



  TemperatureRHS=0.

  TemperatureRHS = FreeEnergy(c,T,phi,c_dot,phi_dot,N_phases+2) !J K/mol
  Cp = g_C  !J/K/mol

  Y_energy=Cp*g_T0  !J/mol  
  TemperatureRHS = TemperatureRHS/Y_energy  !no dimension

  
  do ii=1,3
     do jj=1,3
        do kk=1,3
        T_=vbles(N_phases+2,ii,jj,kk)*stencil3d(ii,jj,kk)
        TemperatureRHS = TemperatureRHS + g_D_tilde*T_/(dx*dx)
     enddo
  enddo
  enddo
! div(D_heat_c grad(dFdc))
  if(heat_c_term)then
     do ii=1,3
        do jj=1,3
           do kk=1,3
              T_=vbles(2,ii,jj,kk)
              c_=vbles(3,ii,jj,kk)
              phi_=vbles(1,ii,jj,kk)
              dFdc(ii,jj,kk)=(FreeEnergy(c_,T_,Phi_,c_dot,phi_dot,3)/X_energy)**2
              D_heat_c =  g_D(1)*phi+(1d0-phi)*g_D(2)
              d_H(ii,jj,kk) = 0.5*c_*(1.-c_)*(D_heat_c/g_Dch)*(X_energy/Y_energy)   
           enddo
        enddo
     enddo
     div=0d0
!z=0
     div = div  + (dFdc(3,2,2)-dFdc(2,2,2))*(d_H(3,2,2)+d_H(2,2,2)) 
     div = div  + (dFdc(1,2,2)-dFdc(2,2,2))*(d_H(1,2,2)+d_H(2,2,2)) 
     div = div  + (dFdc(2,3,2)-dFdc(2,2,2))*(d_H(2,3,2)+d_H(2,2,2)) 
     div = div  + (dFdc(2,1,2)-dFdc(2,2,2))*(d_H(2,1,2)+d_H(2,2,2))
!y=0
     div = div  + (dFdc(2,2,3)-dFdc(2,2,2))*(d_H(2,2,3)+d_H(2,2,2)) 
     div = div  + (dFdc(2,2,1)-dFdc(2,2,2))*(d_H(2,2,1)+d_H(2,2,2)) 
     div = div  + (dFdc(3,2,2)-dFdc(2,2,2))*(d_H(3,2,2)+d_H(2,2,2)) 
     div = div  + (dFdc(1,2,2)-dFdc(2,2,2))*(d_H(1,2,2)+d_H(2,2,2))
!x=0
     div = div  + (dFdc(2,3,2)-dFdc(2,2,2))*(d_H(2,3,2)+d_H(2,2,2)) 
     div = div  + (dFdc(2,1,2)-dFdc(2,2,2))*(d_H(2,1,2)+d_H(2,2,2)) 
     div = div  + (dFdc(2,2,3)-dFdc(2,2,2))*(d_H(2,2,3)+d_H(2,2,2)) 
     div = div  + (dFdc(2,2,1)-dFdc(2,2,2))*(d_H(2,2,1)+d_H(2,2,2))

     
     div =0.5*div/(dx*dx) ! 0.5 comes from the average of d_H. 3D is identical to 2D

     heat_size(3)=max(heat_size(3),abs(div*Y_energy)) !how big is this term
     Q_heat(3)=div*Y_energy
     TemperatureRHS = TemperatureRHS + div
  endif
! grad double well term in phase equation
  if(grad_DW_term)then
     Grad_DW = (g_lambda*GradientEnergy(vbles,1,dx)+potential(Phi,c,1)/g_lambda)*g_W(1)/MolPerVol*phi_dot !J/mol
     heat_size(4)=max(heat_size(4),abs(Grad_DW))
     Q_heat(4)=-Grad_DW
    TemperatureRHS = TemperatureRHS - Grad_DW/(Cp*g_T0)   !no dimension  
  endif


end function TemperatureRHS

double precision function SoluteRHS(vbles,dx,c_c)
! this uses volume element method. This guarantees conservation of solute. The main source of solute is the variation of free energy with phase: d (d F/dc) = ...+ (d^2 F/dc d phi) d phi + ..., which is revealed by a plot of c_dot.
use solution_parameters
implicit none
double precision, intent (in):: dx,c_c,vbles(N_unk+1,3,3,3)
integer ii,jj,ip,ih,jv,kd
double precision c,T,Phi,FreeEnergy,D
double precision :: potential,D_c,f_c,MolPerVol
double precision :: c_dot=0d0,phi_dot=0d0,beta,phi_

if(AlNi.eq.8)approx=.true.

SoluteRHS=0.

!z=0
!x-coordinate transform values
!! Right
 ih=3
 jv=2
 kd=2
 include 'c_stencil.f90'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Left
 ih=1
 jv=2
 kd=2
 include 'c_stencil.f90'


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!up
 ih=2
 jv=3
 kd=2
 include 'c_stencil.f90'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !down
 ih=2
 jv=1
 kd=2
 include 'c_stencil.f90'

!x=0
!! Right
 ih=2!3
 jv=3!2
 kd=2!2
 include 'c_stencil.f90'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Left
 ih=2!1
 jv=1!2
 kd=2!2
 include 'c_stencil.f90'


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!up
 ih=2!2
 jv=2!3
 kd=3!2
 include 'c_stencil.f90'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !down
 ih=2!2
 jv=2!1
 kd=1!2
 include 'c_stencil.f90'
!y=0
!! Right
 ih=2!2!3
 jv=2!3!2
 kd=3!2!2
 include 'c_stencil.f90'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Left
 ih=2!2!1
 jv=2!1!2
 kd=1!2!2
 include 'c_stencil.f90'


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!up
 ih=3!2!2
 jv=2!2!3
 kd=2!3!2
 include 'c_stencil.f90'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !down
 ih=1!2!2
 jv=2!2!1
 kd=2!1!2
 include 'c_stencil.f90'


 g_c_force = SoluteRHS/(g_Dch*dx*dx*g_R*g_T0)
 SoluteRHS = g_c_force

return

end function SoluteRHS


double precision function potential(phi,c,lp)
use solution_parameters
implicit none

double precision, intent (in)::phi,c
integer, intent(in):: lp
integer :: ip,jp

potential=0d0
if(lp.eq.1)then
   potential = 2d0*phi*(1d0-phi)*(1d0-2d0*phi)
elseif(lp.ne.1)then
   potential =  0d0
else
   stop
endif

return


end function potential




double precision function PhaseRHS(vbles,dx,lp,Phi,c,T)
use solution_parameters
implicit none
double precision, intent(in)::dx,c,T,Phi,vbles(N_unk+1,3,3,3)
integer, intent(in):: lp
double precision :: LapPhi
double precision GradientEnergy,potential,FreeEnergy, GradientEnergyPCB,GradientEnergyPCB2,GradientEnergyPCB3, GradientEnergyHexPrism, GradientEnergyHexPrismHess
double precision ::M_tilde 
double precision :: c_dot=0,phi_dot=0
double precision :: MolPerVol,beta,phix,phiy,phiz

if(AlNi.eq.8)then
   if(g_crystal.eq.1)then
      approx = .true.
      M_tilde = 2d1!debug
      phaseRHS=- M_tilde*(GradientEnergyPCB2(vbles,lp,dx)+potential(Phi,c,lp)/g_lambda**2&
           +FreeEnergy(c,T,Phi,0,0,1)/(g_L(1)*g_lambda**2))
   elseif(g_crystal.eq.2)then
      approx = .true.
      M_tilde = 1d0!debug                                                                                                                                  
      phaseRHS=- M_tilde*(GradientEnergyPCB3(vbles,lp,dx)+potential(Phi,c,lp)/g_lambda**2&
           +FreeEnergy(c,T,Phi,0,0,1)/(g_L(1)*g_lambda**2))
   elseif(g_crystal.eq.3)then !hexagonal prism
      approx = .true.
      M_tilde = 1d0!debug
      if(g_kinetic)then
         phix = 0.5*(vbles(1,3,2,2) - vbles(1,1,2,2))/dx 
         phiy = 0.5*(vbles(1,2,3,2) - vbles(1,2,1,2))/dx 
         phiz = 0.5*(vbles(1,2,2,3) - vbles(1,2,2,1))/dx 
         ! This is for needle growth. A_factor can be any positive number
         M_tilde = M_tilde / (1+a_factor*epsilon)*(1 + a_factor*epsilon*abs(phix)/sqrt(phix**2+phiy**2+phiz**2+0.0001))
      endif
      phaseRHS=- M_tilde*( GradientEnergyHexPrism(vbles,lp,dx)+potential(Phi,c,lp)/g_lambda**2&
           +g_driveforce*FreeEnergy(c,T,Phi,0,0,1)/(g_L(1)*g_lambda**2))
   endif
else

   MolPerVol=(1.-c)*g_MolPerVol(1)+c*g_MolPerVol(2)
   
   !assemble and non-dimensionlise phaseRHS 
   
   M_tilde = (1.-c)+c*g_M(2)/g_M(1)
   phaseRHS=- M_tilde*(GradientEnergy(vbles,lp,dx)+potential(Phi,c,lp)/g_lambda**2&
        +FreeEnergy(c,T,Phi,c_dot,phi_dot,1)*MolPerVol/(g_W(1)*g_lambda**2))
endif


end function PhaseRHS






!


subroutine get_RHS_Fi_2d(i,j,k,lb,dx,Fi)
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use physicaldata
  use tree
  implicit none

integer, intent(in)::i,j,k,lb
double precision, intent(in)::dx
double precision, intent(out):: Fi(N_unk)
double precision :: c,T,phi,c_dot,phi_dot
double precision phi_rhs,rfactor,rf4,rf5,rf6,S
double precision TemperatureRHS,SoluteRHS,T_k,phi_
double precision vbles(N_unk+1,3,3,3),unk_(M_unk),MolPerVol
logical Not_number
integer ip,n,ii,jj,kk

  Fi=0.

  do ii=1,3
     do jj=1,3
        do kk=1,3
           
           do ip=1,N_unk
              vbles(ip,ii,jj,kk)=unk(1+(ip-1)*nunkvbles,i+ii-2,j+jj-2,k+kk-2,lb)
           enddo
        enddo
     enddo
  enddo

  
  ip=N_unk+1
  do ii=1,3
     do jj=1,3
        do kk=1,3
           if(mode_1.eq.1) then
              vbles(ip,ii,jj,kk)=(unk(1,i+ii-2,j+jj-2,k+kk-2,lb)-&
                   unk(2,i+ii-2,j+jj-2,k+kk-2,lb))/dt
           else
              rfactor = dt/dtold
              rf4 = rfactor*rfactor/(rfactor+1.0)
              rf5 = (rfactor+1.0)
              rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
              vbles(ip,ii,jj,kk)=(rf4*unk(5,i+ii-2,j+jj-2,k+kk-2,lb)-&
                   rf5*unk(4,i+ii-2,j+jj-2,k+kk-2,lb)+&
                   rf6*unk(1,i+ii-2,j+jj-2,k+kk-2,lb))/dt
           endif
           if(AlNi.ne.8)then
              c=vbles(N_unk,ii,jj,kk)
              vbles(ip,ii,jj,kk)=-vbles(ip,ii,jj,kk)*g_W(1)/((1.-c)+c*g_M(2)/g_M(1))
              MolPerVol=g_MolPerVol(1)*(1.-c)+g_MolPerVol(2)*c
              vbles(ip,ii,jj,kk)=vbles(ip,ii,jj,kk)/MolPerVol*g_lambda**2
           endif
        enddo
     enddo
  enddo
!(1.-c)+c*g_M(2)/g_M(1)
!  if(neigh(1,2,lb).eq.bc_cond_sym.and. i.eq.9)then
!     ii = 3
!     do jj=1,3
!        vbles(2,ii,jj)=delta
!     enddo
!  endif



  do ii=1,M_unk
     unk_(ii) = unk(ii,i,j,k,lb)
  enddo


  call get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_)

!
  if(thermal)then
     unk(1+N_unk*6,i,j,k,lb)=Q_heat(1)
     unk(2+N_unk*6,i,j,k,lb)=Q_heat(2)
     unk(3+N_unk*6,i,j,k,lb)=Q_heat(3)
     unk(4+N_unk*6,i,j,k,lb)=Q_heat(4)
  endif

 
end subroutine get_RHS_Fi_2d

subroutine get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_)

  use solution_parameters

  implicit none

!
integer, intent(in)::mode_1,nunkvbles,lnblocks
double precision, intent(in)::dx,dt,dtold,unk_(M_unk)
double precision, intent(out):: Fi(N_unk)
double precision, intent(in):: vbles(N_unk+1,3,3,3)
double precision :: c,T,phi,c_dot,phi_dot
double precision rfactor,rf4,rf5,rf6,S,GradientEnergy,potential,vbless(N_unk+1,3,3,3)
double precision TemperatureRHS,SoluteRHS,PhaseRHS,T_k,phi_,FreeEnergy,MolPerVol
logical Not_number
integer ip,n,i,j,k

vbless=0d0


Fi=0d0

  if (mode_1.eq.1) then
     phi_dot = (vbles(1,2,2,2)-unk_(2))/dt
     n=nunkvbles
     c_dot = (vbles(2,2,2,2)-unk_(2+n))/dt
  elseif(mode_1.eq.2)then
     rfactor = dt/dtold

     rf4 = rfactor*rfactor/(rfactor+1.0)
     rf5 = (rfactor+1.0)
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
     
     phi_dot = (rf4*unk_(5)-rf5*unk_(4)+rf6*vbles(1,2,2,2))/dt

     n=nunkvbles
     c_dot = (rf4*unk_(5+n)-rf5*unk_(4+n)+rf6*vbles(2,2,2,2))/dt
  else
     write(*,*)"mode_1 not = 1 or 2"
     stop
  endif


  phi=vbles(1,2,2,2)
  if(thermal)  T = vbles(3,2,2,2)
  c = vbles(2,2,2,2)

  



  Fi(1) = PhaseRHS(vbles,dx,1,Phi,c,T)

  if(thermal)then
     Fi(3) = TemperatureRHS(vbles,dx,c,T,phi,c_dot,phi_dot)
  endif

  Fi(2) = SoluteRHS(vbles,dx,c)


  do ip=1,N_unk
     n=(ip-1)*nunkvbles
     if(mode_1.le.1)then
        Fi(ip)=vbles(ip,2,2,2)-dt*Fi(ip)-unk_(n+2)
     else if(mode_1.eq.2)then
        Fi(ip)=vbles(ip,2,2,2)-dt*Fi(ip)/rf6-unk_(n+2)
     end if
  enddo

!  Fi(4)=vbles(4,2,2)-Fi(4)



end subroutine get_RHS_Fi

subroutine Calc_GradPhi(LapPhi,vbles,dx)
use solution_parameters
implicit none
double precision, intent(in)::dx,vbles(N_unk+1,3,3,3)
double precision, intent(out)::LapPhi
double precision, dimension(3,3,3)::stencil3d =reshape((/0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 1., -6., 1., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0./),(/3,3,3/))
integer ip,ii,jj,kk



 LapPhi=0.
 do ii=1,3
    do jj=1,3
       do kk=1,3
          
          LapPhi=LapPhi+stencil3d(ii,jj,kk)*vbles(1,ii,jj,kk)/(dx*dx)
          
       end do
    end do
 end do


 
  end subroutine Calc_GradPhi


double precision function FreeEnergy(c,T,phi,c_dot,phi_dot,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi,c_dot,phi_dot
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(2),hb(2),hr(2)
double precision :: Tk,Cp,R
double precision :: FE
double precision Gibbs_FE_liq,Gibbs_FE_sol,Gibbs_FE_l,Gibbs_FE_e,T_K,g_phi,EntropyMixing,AlFeFreeEnergy !functions
integer :: i,j


if(AlNi.eq.1)then
   if(lp.eq.0)then !free energy proper J/mol
      
      FreeEnergy =              Gibbs_FE_liq(c,T,0,0)*g_phi(phi,0)
      FreeEnergy = FreeEnergy + Gibbs_FE_sol(c,T,0,0)*g_phi(1d0-phi,0)
      
      return
   elseif(lp.eq.1)then !phase
      
      FreeEnergy =              Gibbs_FE_liq(c,T,0,0)*g_phi(phi,1)
      FreeEnergy = FreeEnergy - Gibbs_FE_sol(c,T,0,0)*g_phi(phi,1)
      
      return
   elseif(lp.eq.2)then !Solute
      
      FreeEnergy = FreeEnergy + Gibbs_FE_liq(c,T,1,0)*g_phi(phi,0)
      FreeEnergy = FreeEnergy + Gibbs_FE_sol(c,T,1,0)*(1.-g_phi(phi,0))
      
      return
   elseif(lp.eq.3)then !Temperature Keep in J/mol
      Tk=T_k(T)
      
      !-(1-T*d/dT)dF/dphi
      if(heat_phi_term)then
         !  T (d/dT)dF/dphi
         FreeEnergy =  Tk*(Gibbs_FE_liq(c,T,0,1) - Gibbs_FE_sol(c,T,0,1))
         ! -Df/Dphi
         FreeEnergy = FreeEnergy - Gibbs_FE_liq(c,T,0,0) + Gibbs_FE_sol(c,T,0,0)
         
         FreeEnergy = FreeEnergy*phi_dot*g_phi(phi,1)
         
         Q_heat(1) = FreeEnergy
      endif !end phi_dot_term
      !   heat capacity = -T (d/dT) dF/dT
      g_C =  -Tk*(g_phi(phi,0)*Gibbs_FE_liq(c,T,0,2)+g_phi(1d0-phi,0)*Gibbs_FE_sol(c,T,0,2))
      return
   else
      write(*,*)"bad current_var"
      stop
   endif
elseif(AlNi.eq.2)then
   if(lp.eq.0)then !free energy proper J/mol
      
      FreeEnergy =              Gibbs_FE_l(c,T,0,0)*g_phi(phi,0)
      FreeEnergy = FreeEnergy + Gibbs_FE_e(c,T,0,0)*g_phi(1d0-phi,0)
      
      return
   elseif(lp.eq.1)then !phase
      
      FreeEnergy =              Gibbs_FE_l(c,T,0,0)*g_phi(phi,1)
      FreeEnergy = FreeEnergy - Gibbs_FE_e(c,T,0,0)*g_phi(phi,1)
      
      return
   elseif(lp.eq.2)then !Solute
      
      FreeEnergy = FreeEnergy + Gibbs_FE_l(c,T,1,0)*g_phi(phi,0)
      FreeEnergy = FreeEnergy + Gibbs_FE_e(c,T,1,0)*(1.-g_phi(phi,0))
      
      return
   elseif(lp.eq.3)then !Temperature Keep in J/mol
      Tk=T_k(T)
      
      !-(1-T*d/dT)dF/dphi
      if(heat_phi_term)then
         !  T (d/dT)dF/dphi
         FreeEnergy =  Tk*(Gibbs_FE_l(c,T,0,1) - Gibbs_FE_e(c,T,0,1))
         ! -Df/Dphi
         FreeEnergy = FreeEnergy - Gibbs_FE_l(c,T,0,0) + Gibbs_FE_e(c,T,0,0)
         
         FreeEnergy = FreeEnergy*phi_dot*g_phi(phi,1)
         
         Q_heat(1) = FreeEnergy
      endif !end phi_dot_term
      !   heat capacity = -T (d/dT) dF/dT
      g_C =  -Tk*(g_phi(phi,0)*Gibbs_FE_l(c,T,0,2)+g_phi(1d0-phi,0)*Gibbs_FE_e(c,T,0,2))
      return
   else
      write(*,*)"bad current_var"
      stop
   endif
elseif(AlNi.eq.8)then

   if(approx)then
      if(lp.eq.0)then
         FreeEnergy =   (0.5*g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))**2+g_AlFe_c(2,2)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,3) )*g_phi(phi,0)&
              +  (0.5*g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))**2+g_AlFe_c(1,2)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,3) )*g_phi(1-phi,0)
      elseif(lp.eq.1)then
         FreeEnergy =   (0.5*g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))**2+g_AlFe_c(2,2)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,3) )*g_phi(phi,1)&
              -  (0.5*g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))**2+g_AlFe_c(1,2)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,3) )*g_phi(phi,1)
      elseif(lp.eq.2)then
         FreeEnergy =   (g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,2) )*g_phi(phi,0)&
              +  (g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,2) )*g_phi(1-phi,0)
      elseif(lp.eq.-2)then !f_cc
         FreeEnergy = (1 - g_phi(phi,0))*g_AlFe_c(1,4) + g_phi(phi,0)*g_AlFe_c(2,4)

      elseif(lp.eq.-1)then !f_c_phi
         FreeEnergy = (g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,2) )*g_phi(phi,1)&
              -  (g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,2) )*g_phi(phi,1)

      else
         write(*,*)"lp = ",lp
         stop
      endif

   else

      FreeEnergy =    AlFeFreeEnergy(c,T,phi,lp)
   endif
   return

else
   write(*,*)"models: 1,2 "
   stop
endif
end function FreeEnergy


double precision function AlFeFreeEnergy(c,T,phi,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi(N_phases)
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(2),hb(2),hr(2)
double precision :: Tk,Cp,R,S
double precision :: FE,dh=1d-6
double precision Gibbs_FE_l_AlFe,Gibbs_FE_e_AlFe,Gibbs_FE_e2_AlFe,g_phi,dAlFeFreeEnergy
integer :: i,j
logical :: g_Tflag = .false.


AlFeFreeEnergy=0d0


if(lp.eq.0)then !free energy proper J/mol
!!$   AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,T,0,0)*g_phi(phi(1),0)
!!$   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e_AlFe(c,T,0,0)*g_phi(phi(2),0)
!!$   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e2_AlFe(c,T,0,0)*g_phi(phi(3),0)

   AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,T,0,0)*g_phi(phi(1),0)
   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e2_AlFe(c,T,0,0)*(1-g_phi(phi(1),0))
   return
elseif(lp.le.1)then
   if(g_Tflag)then
      write(*,*)g_Tflag
      stop
   else
      AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,T,0,0)*g_phi(phi(1),1)
      AlFeFreeEnergy = AlFeFreeEnergy -   Gibbs_FE_e2_AlFe(c,T,0,0)*(g_phi(phi(1),1))
      return
   endif
elseif(lp.eq.2)then !Solute


   AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,T,1,0)*g_phi(phi(1),0)
   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e2_AlFe(c,T,1,0)*(1-g_phi(phi(1),0))

elseif(lp.eq.3)then !Temperature Keep in J/mol
   AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,T,0,1)*g_phi(phi(1),0)
   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e_AlFe(c,T,0,1)*g_phi(phi(2),0)
   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e2_AlFe(c,T,0,1)*g_phi(phi(3),0)
   
   S = g_phi(phi(1),0)+g_phi(phi(2),0)+g_phi(phi(3),0)
   
   AlFeFreeEnergy = AlFeFreeEnergy/S 
elseif(lp.eq.-1)then !for anti-trapping
write(*,*)"not coded"
stop
elseif(lp.eq.-2)then !for anti-trapping f_cc
write(*,*)"not coded"
stop
else
   write(*,*)"bad current_var"
   stop
endif

end function

double precision function Gibbs_FE_l_AlFe(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: EntropyMixing,T_k,RK(5)
  double precision :: RKl,ent,gl,glA,glB,Tk,Gibbs_T,Rla
  double precision :: alfehlA(8),alfehlB(8)

! AlFe below
Gibbs_FE_l_AlFe=0

Tk=T_k(T)
if (Tk .gt. 1811.0)then
  alfehlA=AlFehlA3 
  alfehlB=AlFehlB2 
elseif (Tk .gt. 933.47)then
  alfehlA=AlFehlA3 
  alfehlB=AlFehlB1 
elseif (Tk .gt. 700.0)then
  alfehlA=AlFehlA2 
  alfehlB=AlFehlB1 
else
  alfehlA=AlFehlA1 
  alfehlB=AlFehlB1 
endif


if(d_c.eq.0.and.d_T.eq.0)then
   glA=gibbs_T(alfehlA,Tk,0)
   glB=gibbs_T(alfehlB,Tk,0)
   call MixingVector(RK,c,0)
   RKl=       (AlFehrkl0(1)+AlFehrkl0(2)*Tk)*RK(1)
   RKl=RKl + (AlFehrkl1(1)+AlFehrkl1(2)*Tk)*RK(2)
   RKl=RKL + (AlFehrkl2(1)+AlFehrkl2(2)*Tk)*RK(3)

   ent = g_R*Tk*EntropyMixing(c,0)
   gl = glA*(1d0-c)+glB*c+ent+RKl
   Gibbs_FE_l_AlFe=gl
   return
elseif(d_c.eq.0.and.d_T.eq.1)then
   glA=gibbs_T(alfehlA,Tk,1)
   glB=gibbs_T(alfehlB,Tk,1)
   call MixingVector(RK,c,0)
   RKl=      AlFehrkl0(2)*RK(1)
   RKl=RKl + AlFehrkl1(2)*RK(2)
   RKl=RKL + AlFehrkl2(2)*RK(3)
   ent = g_R*EntropyMixing(c,0)
   gl = glA*(1d0-c)+glB*c+ent+RKl
   Gibbs_FE_l_AlFe=gl
   return
elseif(d_c.eq.0.and.d_T.eq.2)then
elseif(d_c.eq.1.and.d_T.eq.0)then
   glA=gibbs_T(alfehlA,Tk,0)
   glB=gibbs_T(alfehlB,Tk,0)
   call MixingVector(RK,c,1)
   RKl=       (AlFehrkl0(1)+AlFehrkl0(2)*Tk)*RK(1)
   RKl=RKl + (AlFehrkl1(1)+AlFehrkl1(2)*Tk)*RK(2)
   RKl=RKL + (AlFehrkl2(1)+AlFehrkl2(2)*Tk)*RK(3)

   ent = g_R*Tk*EntropyMixing(c,1)
   gl = -glA + glB+ent+RKl
   Gibbs_FE_l_AlFe=gl

elseif(d_c.eq.1.and.d_T.eq.1)then
elseif(d_c.eq.1.and.d_T.eq.2)then
elseif(d_c.eq.2.and.d_T.eq.0)then

   call MixingVector(RK,c,2)
   RKl=       (AlFehrkl0(1)+AlFehrkl0(2)*Tk)*RK(1)
   RKl=RKl + (AlFehrkl1(1)+AlFehrkl1(2)*Tk)*RK(2)
   RKl=RKL + (AlFehrkl2(1)+AlFehrkl2(2)*Tk)*RK(3)

   ent = g_R*Tk*EntropyMixing(c,2)
   gl = ent+RKl
   Gibbs_FE_l_AlFe=gl
else
   write(*,*)"unexpected d_c,d_T combination",d_c,d_T
   stop
endif

end function Gibbs_FE_l_AlFe

double precision function AlFey1(c,d)
integer, intent(in)::d
double precision, intent(in) :: c
if(c.gt.0.27246.or.c.lt.0.2350)then
   write(*,*)"c out of bound in AlFey1()"
   stop
endif
if(d.eq.0)then
      AlFey1=(0.2350/c - 0.8625)/0.1375 
      return
elseif(d.eq.1)then
      AlFey1=-(0.2350/c**2)/0.1375 
      return
else
   write(*,*)"d =0 or 1"
endif



end function AlFey1


double precision function Gibbs_FE_e_AlFe(c,T,d_c,d_T)
  use solution_parameters
  use time_dep_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: EntropyMixing,T_k,RK(5)
  double precision :: RKa,ent,ga,gaA,gaB,Tk,Gibbs_T
  double precision :: haA(8),haB(8)




! AlFe below
Tk=T_k(T)
if (Tk .gt. 1811.0)then
  haA=AlFehaA3 
  haB=AlFehaB2 
elseif (Tk .gt. 933.47)then
  haA=AlFehaA3 
  haB=AlFehaB1 
elseif (Tk .gt. 700.0)then
  haA=AlFehaA2 
  haB=AlFehaB1 
else
  haA=AlFehaA1 
  haB=AlFehaB1 
endif







if(d_c.eq.0.and.d_T.eq.0)then
   gaA=gibbs_T(haA,Tk,0)
   gaB=gibbs_T(haB,Tk,0)
   call MixingVector(RK,c,0)
   RKa=      (AlFehrka0(1)+AlFehrka0(2)*Tk)*RK(1)
   RKa=RKa + (AlFehrka1(1)+AlFehrka1(2)*Tk)*RK(2)
   ent = g_R*Tk*EntropyMixing(c,0)
   ga = gaA*(1d0-c)+gaB*c+ent+RKa
   Gibbs_FE_e_AlFe=ga
elseif(d_c.eq.0.and.d_T.eq.1)then
   gaA=gibbs_T(haA,Tk,1)
   gaB=gibbs_T(haB,Tk,1)
   call MixingVector(RK,c,0)
   RKa=      AlFehrka0(2)*RK(1)
   RKa=RKa + AlFehrka1(2)*RK(2)
   ent = g_R*EntropyMixing(c,0)
   ga = gaA*(1d0-c)+gaB*c+ent+RKa
   Gibbs_FE_e_AlFe=ga
elseif(d_c.eq.0.and.d_T.eq.2)then
elseif(d_c.eq.1.and.d_T.eq.0)then
   gaA=gibbs_T(haA,Tk,0)
   gaB=gibbs_T(haB,Tk,0)
   call MixingVector(RK,c,1)
   RKa=      (AlFehrka0(1)+AlFehrka0(2)*Tk)*RK(1)
   RKa=RKa + (AlFehrka1(1)+AlFehrka1(2)*Tk)*RK(2)
   ent = g_R*Tk*EntropyMixing(c,1)
   ga = -gaA + gaB+ent+RKa
   Gibbs_FE_e_AlFe=ga
elseif(d_c.eq.1.and.d_T.eq.1)then
elseif(d_c.eq.1.and.d_T.eq.2)then
elseif(d_c.eq.2.and.d_T.eq.0)then
else
   write(*,*)"unexpected d_c,d_T combination",d_c,d_T
   stop
endif

end function Gibbs_FE_e_AlFe

!!!!!!!!!!!!!!!!!!!!!!!!!1
!Al (of the AlFe)
double precision function Gibbs_FE_e2_AlFe(c,T,d_c,d_T)
  use solution_parameters
  use time_dep_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: EntropyMixing,T_k,RK(5)
  double precision :: Tk,Gibbs_T,entSLb,dentSLb,AlFey1,y1,dy1
  double precision :: alfeh1(8),alfeh2(8),sLb1,dsLb1,sLb2,dsLb2,gb1,gb,dgb1,dgb



! AlFe below
Gibbs_FE_e2_AlFe=0



Tk=T_k(T)
if (Tk .gt. 1811.0)then
  alfeh1=AlFehslb1_4 
  alfeh2=AlFehslb2_4 
elseif (Tk .gt. 933.47)then
  alfeh1=AlFehslb1_3 
  alfeh2=AlFehslb2_3 
elseif (Tk .gt. 700.0)then
  alfeh1=AlFehslb1_2 
  alfeh2=AlFehslb2_2 
else
  alfeh1=AlFehslb1_1 
  alfeh2=AlFehslb2_1 
endif

if(d_c.eq.0.and.d_T.eq.0)then
   sLb1=gibbs_T(alfeh1,Tk,0)
   sLb2=gibbs_T(alfeh2,Tk,0)
   y1 = AlFey1(max(min(0.27246,c),0.2350),0)
   entSLb=0.1375*g_R*Tk*EntropyMixing(y1,0)
   gb1=SLb1*y1 + SLb2*(1d0-y1) 
   gb=(gb1 + entSLb)/(0.8625 + 0.1375*y1)
   Gibbs_FE_e2_AlFe=gb
   return


elseif(d_c.eq.0.and.d_T.eq.1)then
   sLb1=gibbs_T(alfeh1,Tk,1)
   sLb2=gibbs_T(alfeh2,Tk,1)
   y1 = AlFey1(max(min(0.27246,c),0.2350),0)
   entSLb=0.1375*g_R*EntropyMixing(y1,0)
   gb1=SLb1*y1 + SLb2*(1d0-y1) 
   gb=(gb1 + entSLb)/(0.8625 + 0.1375*y1)
   Gibbs_FE_e2_AlFe=gb
   return
elseif(d_c.eq.0.and.d_T.eq.2)then
elseif(d_c.eq.1.and.d_T.eq.0)then
   sLb1=gibbs_T(alfeh1,Tk,0)
   sLb2=gibbs_T(alfeh2,Tk,0)
!   dsLb1=0d0
!   dsLb2=0d0
   y1 = AlFey1(max(min(0.27246,c),0.2350),0)
   dy1 = AlFey1(max(min(0.27246,c),0.2350),1)
   entSLb = 0.1375*g_R*Tk*EntropyMixing(y1,0)
   dentSLb = 0.1375*g_R*Tk*EntropyMixing(y1,1)
   gb1 =SLb1*y1 + SLb2*(1d0-y1) 
   dgb1 = (SLb1 - SLb2)
   dgb = -0.1375*(gb1 + entSLb)/(0.8625 + 0.1375*y1)**2 + (dgb1 + dentSLb)/(0.8625 + 0.1375*y1)
   Gibbs_FE_e2_AlFe = dgb*dy1
   return
elseif(d_c.eq.1.and.d_T.eq.1)then
elseif(d_c.eq.1.and.d_T.eq.2)then
elseif(d_c.eq.2.and.d_T.eq.0)then
else
   write(*,*)"unexpected d_c,d_T combination",d_c,d_T
   stop
endif

end function Gibbs_FE_e2_AlFe


double precision function EntropyMixing(c,d)
double precision, intent (in)::c
integer, intent (in)::d

if(c.eq.1d0.or.c.eq.0d0)then
   write(*,*)"c=",c, "in EntropyMixing"
   stop
endif
 
 if(d.eq.1)then
  EntropyMixing = log(c/(1.-c))
 else
  EntropyMixing = c*log(c)+(1.-c)*log(1.-c)
 endif
end function EntropyMixing


subroutine GibbsVector(V,T,d)
double precision, intent (out)::V(8)
double precision, intent (in)::T
integer, intent (in)::d
if(d.eq.0)then
 V(1)=1.0
 V(2)=T 
 V(3)=T*log(T)
 V(4)=T**2
 V(5)=T**3
 V(6)=1./T
 V(7)=T**7
 V(8)=1.0/T**9
else if(d.eq.1)then
 V(1)=0.
 V(2)=1.
 V(3)=1.+log(T)
 V(4)=2.*T
 V(5)=3.*T**2
 V(6)=-1./T**2
 V(7)=7.*T**6
 V(8)=-9.0/T**10
elseif(d.eq.2)then
 V(1)=0.
 V(2)=0.
 V(3)=1./T
 V(4)=2.
 V(5)=6.*T
 V(6)=2./T**3
 V(7)=42.*T**5
 V(8)=90./T**11
else
 write(*,*)"bad diff order in GibbsVector"
stop
endif
end subroutine GibbsVector
subroutine MixingVector(RK,c,d)
double precision, intent (out)::RK(5)
double precision, intent (in)::c
integer, intent (in)::d
integer i
if(d.eq.0)then
 RK(1)=c*(1.-c)
 RK(2)=RK(1)*(1.-2.*c)
 RK(3)=RK(2)*(1.-2.*c)
 RK(4)=RK(3)*(1.-2.*c)
 RK(5)=RK(4)*(1.-2.*c)
elseif(d.eq.1)then
 RK(1)=1.-2.*c
 RK(2)=1.+(-6.+6.*c)*c
 RK(3)=1.+(-10.+(24.-16.*c)*c)*c
 RK(4)=1.+(-14.+(54.+(-80.+40.*c)*c)*c)*c
 RK(5)=1.+(-18.+(96.+(-224.+(240.-96.*c)*c)*c)*c)*c
elseif(d.eq.2)then
 RK(1)=-2.
 RK(2)=12.*c-6. 
 RK(3)=-10.+(48.-48.*c)*c
 RK(4)=-14.+(108.+(-240.+160.*c)*c)*c
 RK(5)=-18.+(192.+(-672.+(960.-480.*c)*c)*c)*c
elseif(d.eq.3)then
 RK(1)=0.
 RK(2)=12.
 RK(3)=-96.*c+48.
 RK(4)=108.+(-480.+480.*c)*c
 RK(5)=192.+(-1344.+(2880.-1920.*c)*c)*c
else
 write(*,*)"bad diff order in MixingVector"
 stop
endif
end subroutine MixingVector


double precision function g_phi(phi,d)
  double precision, intent (in):: phi
  integer, intent (in):: d
  if(d.eq.1)then
   g_phi = 6.*phi*(1.-phi)
!    g_phi = 30.*phi**4-60.*phi**3+30.*phi**2
  elseif(d.eq.2)then
   g_phi = 6.-12.*phi
!    g_phi = 120.*phi**3-180.*phi**2+60.*phi
  else !assume 0
   g_phi = 3.*phi**2-2.*phi**3
!    g_phi = 6.*phi**5-15.*phi**4+10.*phi**3
  endif
end function g_phi


logical function Not_number(x)
double precision, intent(in):: x
if(x.eq.x)then
Not_number=.false.
else
Not_number=.true.
endif
end function Not_number

double precision function f_c(x,df)
double precision, intent (in)::x
integer, intent(in)::df

f_c=0.
if(df.eq.0)then
if(x*(1-x).gt.0.)f_c=x*(1.-x)
elseif(df.eq.1)then
 if(x*(1-x).gt.0.)f_c=1.-2.*x
endif
end function f_c
 
subroutine get_y1_y2(y,d,y1_y2)
integer, intent(in)::d
double precision, intent(in) :: y(2)
double precision, intent(out) :: y1_y2(6)



y1_y2=0d0
if(d.eq.0)then
   y1_y2(1)=y(1)*y(2)
   y1_y2(2)=(1d0-y(1))*y(2)
   y1_y2(3)=(1d0-y(2))*y(1)
   y1_y2(4)=(1d0-y(1))*(1d0-y(2))
   y1_y2(5)=y(1)*(1d0-y(1))
   y1_y2(6)=y(2)*(1d0-y(2))
elseif(d.eq.1)then
   y1_y2(1)=y(2)
   y1_y2(2)=-y(2)
   y1_y2(3)=(1d0-y(2))
   y1_y2(4)=-(1d0-y(2))
   y1_y2(5)=1d0-2d0*y(1)
   y1_y2(6)=0d0
elseif(d.eq.2)then
   y1_y2(1)=y(1)
   y1_y2(2)=(1d0-y(1))
   y1_y2(3)=-y(1)
   y1_y2(4)=-(1d0-y(1))
   y1_y2(5)=0d0
   y1_y2(6)=1d0-2d0*y(2)
else
   write(*,*)"d=0 to 2"
   stop
endif
end subroutine get_y1_y2




double precision function y1(c,d)
integer, intent(in)::d
double precision, intent(in) :: c
if(c.gt.1d0-1d-9.or.c.lt.1d-9)then
   write(*,*)"c out of bound in y1()"
   stop
endif
if(d.eq.0)then
   if(c.ge.0.4d0)then
      y1=1d0-1d-5
      return
   else
      y1=2.5*c
      return
   endif
elseif(d.eq.1)then
   if(c.ge.0.4d0)then
      y1=1d-5
      return
   else
      y1=2.5
      return
!   else
!      write(*,*)"d y1/dc not determined at c=",c
!      stop
   endif
else
   write(*,*)"d =0 or 1"
endif
end function y1

double precision function y2(c,d)
integer, intent(in)::d
double precision, intent(in) :: c
if(c.gt.1d0-1d-9.or.c.lt.1d-9)then
   write(*,*)"c out of bound in y2()"
   stop
endif

if(d.eq.0)then
   if(c.ge.0.5d0)then
      y2=1d-5
      return
   elseif(c.gt.0.4d0)then
      y2=1d0-(2d0-5d0*c)/(c-1d0)
      return
   else
      y2=1d0-1d-5
      return
   endif
elseif(d.eq.1)then
!   if(c.eq.0.5.or.c.eq.0.4)then
!       write(*,*)"d y2/dc not determined at c =",c
!       stop
   if(c.gt.0.5d0)then
      y2=1d-5
      return
   elseif(c.gt.0.4d0)then
      y2=-3d0/(1d0-c)**2
      return
   else
      y2=1d-5
      return
   endif
else
   write(*,*)"d =0 or 1"
endif
end function y2





double precision function y3(c,d)
integer, intent(in)::d
double precision, intent(in) :: c
if(c.gt.1d0-1d-9.or.c.lt.1d-9)then
   write(*,*)"c out of bound in y3()"
   stop
endif
if(d.eq.0)then
   if(c.gt.0.5d0)then
      y3=2d0*(1d0-c)
      return
   else
      y3=1d0-1d-5
      return
   endif
elseif(d.eq.1)then
   if(c.gt.0.5d0)then
      y3=-2d0
      return
   else
      y3=0d0
      return
!   else
!      write(*,*)"d y3/dc not determined at c=",c
!      stop
   endif
else
   write(*,*)"d =0 or 1"
endif
end function y3


double precision function y4(c,d)
integer, intent(in)::d
double precision, intent(in) :: c
if(c.gt.1d0-1d-9.or.c.lt.1d-9)then
   write(*,*)"c out of bound in y4()",c
   stop
endif
if(d.eq.0)then
   if(c.ge.0.5d0)then
      y4=1d0-1d-5
      return
   else
      y4=c/(1d0-c)
      return
   endif
elseif(d.eq.1)then
   if(c.ge.0.5d0)then
      y4=0d0
      return
   else
      y4=1d0/(1d0-c)**2
      return
!   else
!      write(*,*)"d y4/dc not determined at c=",c
!      stop
   endif
else
   write(*,*)"d =0 or 1"
endif
end function y4



double precision function gibbs_T(h,T,d)
double precision, intent(in) :: T,h(8)
integer, intent(in) :: d
double precision V(8)
integer i
gibbs_T=0d0
call GibbsVector(V,T,d)
do i= 1,8
   gibbs_T=gibbs_T+h(i)*V(i)
enddo
end function gibbs_T


double precision function Gibbs_FE_l(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5)
  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gl_,Tk,Gibbs_T
  double precision :: hlA_(8),hlB_(8)
  Tk=T_k(T)
  if (Tk > 1728.0) then
     hlA_=hlA3 
     hlB_=hlB2 
  elseif (Tk > 933.47)then
     hlA_=hlA3 
     hlB_=hlB1 
  elseif (Tk > 700.0)then
     hlA_=hlA2 
     hlB_=hlB1 
  else
     hlA_=hlA1 
     hlB_=hlB1 
  endif
  if(d_c.eq.0.and.d_T.eq.0)then
     Tk=T_k(T)
     glA=gibbs_T(hlA_,Tk,0)
     glB=gibbs_T(hlB_,Tk,0)
     call MixingVector(RK,c,0)
     RKl_=         (hrkl0(1)+hrkl0(2)*Tk)*RK(1)
     RKl_=RKl_ + (hrkl1(1)+hrkl1(2)*Tk)*RK(2)
     RKl_=RKl_ + (hrkl2(1)+hrkl2(2)*Tk)*RK(3)
     RKl_=RKl_ + (hrkl3(1)+hrkl3(2)*Tk)*RK(4)
     RKl_=RKl_ + (hrkl4(1)+hrkl4(2)*Tk)*RK(5)
     ent_=g_R*Tk*EntropyMixing(c,0)
     gl_=glA*(1d0-c)+glB*c+ent_+RKl_
     Gibbs_FE_l=gl_
  elseif(d_c.eq.0.and.d_T.eq.1)then
     Tk=T_k(T)
     glA=gibbs_T(hlA_,Tk,1)
     glB=gibbs_T(hlB_,Tk,1)
     call MixingVector(RK,c,0)
     RKl_=         hrkl0(2)*RK(1)
     RKl_=RKl_ +   hrkl1(2)*RK(2)
     RKl_=RKl_ +   hrkl2(2)*RK(3)
     RKl_=RKl_ +   hrkl3(2)*RK(4)
     RKl_=RKl_ +   hrkl4(2)*RK(5)
     ent_=g_R*EntropyMixing(c,0)
     gl_=glA*(1d0-c)+glB*c+ent_+RKl_
     Gibbs_FE_l=gl_
  elseif(d_c.eq.0.and.d_T.eq.2)then
     Tk=T_k(T)
     glA=gibbs_T(hlA_,Tk,2)
     glB=gibbs_T(hlB_,Tk,2)
     call MixingVector(RK,c,0)
     gl_=glA*(1d0-c)+glB*c
     Gibbs_FE_l=gl_
  elseif(d_c.eq.1.and.d_T.eq.0)then
     Tk=T_k(T)
     glA=gibbs_T(hlA_,Tk,0)
     glB=gibbs_T(hlB_,Tk,0)
     call MixingVector(RK,c,1)
     RKl_=         (hrkl0(1)+hrkl0(2)*Tk)*RK(1)
     RKl_=RKl_ + (hrkl1(1)+hrkl1(2)*Tk)*RK(2)
     RKl_=RKl_ + (hrkl2(1)+hrkl2(2)*Tk)*RK(3)
     RKl_=RKl_ + (hrkl3(1)+hrkl3(2)*Tk)*RK(4)
     RKl_=RKl_ + (hrkl4(1)+hrkl4(2)*Tk)*RK(5)
     ent_=g_R*Tk*EntropyMixing(c,1)
     gl_=-glA+glB+ent_+RKl_
     Gibbs_FE_l=gl_
  elseif(d_c.eq.1.and.d_T.eq.1)then
     Tk=T_k(T)
     glA=gibbs_T(hlA_,Tk,1)
     glB=gibbs_T(hlB_,Tk,1)
     call MixingVector(RK,c,0)
     RKl_=         hrkl0(2)*RK(1)
     RKl_=RKl_ +   hrkl1(2)*RK(2)
     RKl_=RKl_ +   hrkl2(2)*RK(3)
     RKl_=RKl_ +   hrkl3(2)*RK(4)
     RKl_=RKl_ +   hrkl4(2)*RK(5)
     ent_=g_R*EntropyMixing(c,1)
     gl_=-glA+glB+ent_+RKl_
     Gibbs_FE_l=gl_
  elseif(d_c.eq.1.and.d_T.eq.2)then
     Tk=T_k(T)
     glA=gibbs_T(hlA_,Tk,2)
     glB=gibbs_T(hlB_,Tk,2)
     call MixingVector(RK,c,0)
     gl_=-glA+glB
     Gibbs_FE_l=gl_


  else
     write(*,*)"unexpected d_c,d_T combination",d_c,d_T
     stop
  endif

  end function Gibbs_FE_l


double precision function Gibbs_FE_e(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in)::d_c,d_T
  double precision :: SLe(6),yy(2),yy_d(2),y1_y2(6),y1_y2_d1(6),y1_y2_d2(6)
  double precision :: y3,y4,EntropyMixing,T_k
  double precision ::entSL_,entSL_d,ge1_,ge1_d1,ge1_d2,ge2_,ge2_d1,ge2_d2,ge_,ge_d,Tk,Gibbs_T
  double precision :: hlA_(8),hlB_(8),haA_(8),haB_(8),hslb1_(8),hslb2_(8),hslb3_(8),hslb4_(8),hc_(8),hsle1_(8),hsle2_(8),hsle3_(8),hsle4_(8)
  
 Tk=T_k(T)
  if (Tk > 1728.0) then
  hlA_=hlA3 
  hlB_=hlB2 
  haA_=haA3 
  haB_=haB2 
  hslb1_=hslb1_4 
  hslb2_=hslb2_3 
  hslb3_=hslb3_4 
  hslb4_=hslb4_4 
  hc_=hc4 
  hsle1_=hsle1_4 
  hsle2_=hsle2_3 
  hsle3_=hsle3_2 
  hsle4_=hsle4_2 
elseif (Tk > 933.47)then
  hlA_=hlA3 
  hlB_=hlB1 
  haA_=haA3 
  haB_=haB1 
  hslb1_=hslb1_3 
  hslb2_=hslb2_3 
  hslb3_=hslb3_3 
  hslb4_=hslb4_3 
  hc_=hc3 
  hsle1_=hsle1_3 
  hsle2_=hsle2_3 
  hsle3_=hsle3_1 
  hsle4_=hsle4_1 
elseif (Tk > 700.0)then
  hlA_=hlA2 
  hlB_=hlB1 
  haA_=haA2 
  haB_=haB1 
  hslb1_=hslb1_2 
  hslb2_=hslb2_2 
  hslb3_=hslb3_2 
  hslb4_=hslb4_2 
  hc_=hc2 
  hsle1_=hsle1_2 
  hsle2_=hsle2_2 
  hsle3_=hsle3_1 
  hsle4_=hsle4_1 
else
  hlA_=hlA1 
  hlB_=hlB1 
  haA_=haA1 
  haB_=haB1 
  hslb1_=hslb1_1 
  hslb2_=hslb2_1 
  hslb3_=hslb3_1 
  hslb4_=hslb4_1 
  hc_=hc1 
  hsle1_=hsle1_1 
  hsle2_=hsle2_1 
  hsle3_=hsle3_1 
  hsle4_=hsle4_1 
endif

  if(d_c.eq.0.and.d_T.eq.0)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(hsle1_,Tk,0)
     SLe(2)=gibbs_T(hsle2_,Tk,0)
     SLe(3)=gibbs_T(hsle3_,Tk,0)
     SLe(4)=gibbs_T(hsle4_,Tk,0)
     SLe(5)=hsle5(1)+hsle5(2)*Tk
     SLe(6)=hsle6(1)+hsle6(2)*Tk

     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     entSL_=g_R*Tk*(EntropyMixing(yy(1),0)+EntropyMixing(yy(2),0))
     call get_y1_y2(yy,0,y1_y2)
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge2_=SLe(5)*y1_y2(5)+SLe(6)*y1_y2(6)
     ge_=(ge1_+ge2_+entSL_)/(1d0+yy(1))
     Gibbs_FE_e=ge_
  elseif(d_c.eq.1.and.d_T.eq.0)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(hsle1_,Tk,0)
     SLe(2)=gibbs_T(hsle2_,Tk,0)
     SLe(3)=gibbs_T(hsle3_,Tk,0)
     SLe(4)=gibbs_T(hsle4_,Tk,0)
     SLe(5)=hsle5(1)+hsle5(2)*Tk
     SLe(6)=hsle6(1)+hsle6(2)*Tk

     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     yy_d(1)=y4(c,1)
     yy_d(2)=y3(c,1)
     entSL_=g_R*Tk*(EntropyMixing(yy(1),0)+EntropyMixing(yy(2),0))
     entSL_d=g_R*Tk*(EntropyMixing(yy(1),1)*yy_d(1)+EntropyMixing(yy(2),1)*yy_d(2))


     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1)
     call get_y1_y2(yy,2,y1_y2_d2)
     
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge2_=SLe(5)*y1_y2(5)+SLe(6)*y1_y2(6)
     ge1_d1=(SLe(1)*y1_y2_d1(1)+SLe(2)*y1_y2_d1(2)+SLe(3)*y1_y2_d1(3)+SLe(4)*y1_y2_d1(4))*yy_d(1)
     ge1_d2=(SLe(1)*y1_y2_d2(1)+SLe(2)*y1_y2_d2(2)+SLe(3)*y1_y2_d2(3)+SLe(4)*y1_y2_d2(4))*yy_d(2)
     ge2_d1=(SLe(5)*y1_y2_d1(5)+SLe(6)*y1_y2_d1(6))*yy_d(1)
     ge2_d2=(SLe(5)*y1_y2_d2(5)+SLe(6)*y1_y2_d2(6))*yy_d(2)
     ge_=(ge1_+ge2_+entSL_)/(1d0+yy(1))
     ge_d=-ge_/(1d0+yy(1))*yy_d(1)+(ge1_d1+ge1_d2+ge2_d1+ge2_d2+entSL_d)/(1d0+yy(1))
     Gibbs_FE_e=ge_d
  elseif(d_c.eq.1.and.d_T.eq.1)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(hsle1_,Tk,1)
     SLe(2)=gibbs_T(hsle2_,Tk,1)
     SLe(3)=gibbs_T(hsle3_,Tk,1)
     SLe(4)=gibbs_T(hsle4_,Tk,1)
     SLe(5)=hsle5(2)
     SLe(6)=hsle6(2)

     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     yy_d(1)=y4(c,1)
     yy_d(2)=y3(c,1)
     entSL_=g_R*(EntropyMixing(yy(1),0)+EntropyMixing(yy(2),0))
     entSL_d=g_R*(EntropyMixing(yy(1),1)*yy_d(1)+EntropyMixing(yy(2),1)*yy_d(2))


     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1)
     call get_y1_y2(yy,2,y1_y2_d2)
     
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge2_=SLe(5)*y1_y2(5)+SLe(6)*y1_y2(6)
     ge1_d1=(SLe(1)*y1_y2_d1(1)+SLe(2)*y1_y2_d1(2)+SLe(3)*y1_y2_d1(3)+SLe(4)*y1_y2_d1(4))*yy_d(1)
     ge1_d2=(SLe(1)*y1_y2_d2(1)+SLe(2)*y1_y2_d2(2)+SLe(3)*y1_y2_d2(3)+SLe(4)*y1_y2_d2(4))*yy_d(2)
     ge2_d1=(SLe(5)*y1_y2_d1(5)+SLe(6)*y1_y2_d1(6))*yy_d(1)
     ge2_d2=(SLe(5)*y1_y2_d2(5)+SLe(6)*y1_y2_d2(6))*yy_d(2)
     ge_=(ge1_+ge2_+entSL_)/(1d0+yy(1))
     ge_d=-ge_/(1d0+yy(1))*yy_d(1)+(ge1_d1+ge1_d2+ge2_d1+ge2_d2+entSL_d)/(1d0+yy(1))
     Gibbs_FE_e=ge_d
  elseif(d_c.eq.1.and.d_T.eq.2)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(hsle1_,Tk,2)
     SLe(2)=gibbs_T(hsle2_,Tk,2)
     SLe(3)=gibbs_T(hsle3_,Tk,2)
     SLe(4)=gibbs_T(hsle4_,Tk,2)


     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     yy_d(1)=y4(c,1)
     yy_d(2)=y3(c,1)



     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1)
     call get_y1_y2(yy,2,y1_y2_d2)
     
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge1_d1=(SLe(1)*y1_y2_d1(1)+SLe(2)*y1_y2_d1(2)+SLe(3)*y1_y2_d1(3)+SLe(4)*y1_y2_d1(4))*yy_d(1)
     ge1_d2=(SLe(1)*y1_y2_d2(1)+SLe(2)*y1_y2_d2(2)+SLe(3)*y1_y2_d2(3)+SLe(4)*y1_y2_d2(4))*yy_d(2)
     ge_=ge1_/(1d0+yy(1))
     ge_d=-ge_/(1d0+yy(1))*yy_d(1)+(ge1_d1+ge1_d2)/(1d0+yy(1))
     Gibbs_FE_e=ge_d
  elseif(d_c.eq.0.and.d_T.eq.1)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(hsle1_,Tk,1)
     SLe(2)=gibbs_T(hsle2_,Tk,1)
     SLe(3)=gibbs_T(hsle3_,Tk,1)
     SLe(4)=gibbs_T(hsle4_,Tk,1)
     SLe(5)=hsle5(2)
     SLe(6)=hsle6(2)

     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     entSL_=g_R*(EntropyMixing(yy(1),0)+EntropyMixing(yy(2),0))
     call get_y1_y2(yy,0,y1_y2)
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge2_=SLe(5)*y1_y2(5)+SLe(6)*y1_y2(6)
     ge_=(ge1_+ge2_+entSL_)/(1d0+yy(1))
     Gibbs_FE_e=ge_
  elseif(d_c.eq.0.and.d_T.eq.2)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(hsle1_,Tk,2)
     SLe(2)=gibbs_T(hsle2_,Tk,2)
     SLe(3)=gibbs_T(hsle3_,Tk,2)
     SLe(4)=gibbs_T(hsle4_,Tk,2)
     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     call get_y1_y2(yy,0,y1_y2)
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge_=ge1_/(1d0+yy(1))
     Gibbs_FE_e=ge_
  endif
end function Gibbs_FE_e

double precision function Gibbs_FE_liq(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5)
  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gl_,Tk,Gibbs_T

  
  if(d_c.eq.0.and.d_T.eq.0)then
     Tk=T_k(T)
     glA=gibbs_T(hlA,Tk,0)
     glB=gibbs_T(hlB,Tk,0)
     call MixingVector(RK,c,0)
     RKl_=         (hrkl0(1)+hrkl0(2)*Tk)*RK(1)
     RKl_=RKl_ + (hrkl1(1)+hrkl1(2)*Tk)*RK(2)
     RKl_=RKl_ + (hrkl2(1)+hrkl2(2)*Tk)*RK(3)
     RKl_=RKl_ + (hrkl3(1)+hrkl3(2)*Tk)*RK(4)
     RKl_=RKl_ + (hrkl4(1)+hrkl4(2)*Tk)*RK(5)
     ent_=g_R*Tk*EntropyMixing(c,0)
     gl_=glA*(1d0-c)+glB*c+ent_+RKl_
     Gibbs_FE_liq=gl_
  elseif(d_c.eq.0.and.d_T.eq.1)then
     Tk=T_k(T)
     glA=gibbs_T(hlA,Tk,1)
     glB=gibbs_T(hlB,Tk,1)
     call MixingVector(RK,c,0)
     RKl_=         hrkl0(2)*RK(1)
     RKl_=RKl_ +   hrkl1(2)*RK(2)
     RKl_=RKl_ +   hrkl2(2)*RK(3)
     RKl_=RKl_ +   hrkl3(2)*RK(4)
     RKl_=RKl_ +   hrkl4(2)*RK(5)
     ent_=g_R*EntropyMixing(c,0)
     gl_=glA*(1d0-c)+glB*c+ent_+RKl_
     Gibbs_FE_liq=gl_
  elseif(d_c.eq.0.and.d_T.eq.2)then
     Tk=T_k(T)
     glA=gibbs_T(hlA,Tk,2)
     glB=gibbs_T(hlB,Tk,2)
     call MixingVector(RK,c,0)
     gl_=glA*(1d0-c)+glB*c
     Gibbs_FE_liq=gl_
  elseif(d_c.eq.1.and.d_T.eq.0)then
     Tk=T_k(T)
     glA=gibbs_T(hlA,Tk,0)
     glB=gibbs_T(hlB,Tk,0)
     call MixingVector(RK,c,1)
     RKl_=         (hrkl0(1)+hrkl0(2)*Tk)*RK(1)
     RKl_=RKl_ + (hrkl1(1)+hrkl1(2)*Tk)*RK(2)
     RKl_=RKl_ + (hrkl2(1)+hrkl2(2)*Tk)*RK(3)
     RKl_=RKl_ + (hrkl3(1)+hrkl3(2)*Tk)*RK(4)
     RKl_=RKl_ + (hrkl4(1)+hrkl4(2)*Tk)*RK(5)
     ent_=g_R*Tk*EntropyMixing(c,1)
     gl_=-glA+glB+ent_+RKl_
     Gibbs_FE_liq=gl_
  elseif(d_c.eq.1.and.d_T.eq.1)then
     Tk=T_k(T)
     glA=gibbs_T(hlA,Tk,1)
     glB=gibbs_T(hlB,Tk,1)
     call MixingVector(RK,c,0)
     RKl_=         hrkl0(2)*RK(1)
     RKl_=RKl_ +   hrkl1(2)*RK(2)
     RKl_=RKl_ +   hrkl2(2)*RK(3)
     RKl_=RKl_ +   hrkl3(2)*RK(4)
     RKl_=RKl_ +   hrkl4(2)*RK(5)
     ent_=g_R*EntropyMixing(c,1)
     gl_=-glA+glB+ent_+RKl_
     Gibbs_FE_liq=gl_
  elseif(d_c.eq.1.and.d_T.eq.2)then
     Tk=T_k(T)
     glA=gibbs_T(hlA,Tk,2)
     glB=gibbs_T(hlB,Tk,2)
     call MixingVector(RK,c,0)
     gl_=-glA+glB
     Gibbs_FE_liq=gl_
  else
     write(*,*)"unexpected d_c,d_T combination",d_c,d_T
     stop
  endif

  end function Gibbs_FE_liq




double precision function Gibbs_FE_sol(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in)::d_c,d_T
  double precision :: SLb(6),yy(2),yy_d(2),y1_y2(6),y1_y2_d1(6),y1_y2_d2(6)
  double precision :: y1,y2,EntropyMixing,T_k
  double precision ::entSL_,entSL_d,gb1_,gb1_d1,gb1_d2,gb2_,gb2_d1,gb2_d2,gb_,gb_d,Tk,Gibbs_T
  
  if(d_c.eq.0.and.d_T.eq.0)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(hslb1,Tk,0)
     SLb(2)=gibbs_T(hslb2,Tk,0)
     SLb(3)=gibbs_T(hslb3,Tk,0)
     SLb(4)=gibbs_T(hslb4,Tk,0)
     SLb(5)=hslb5(1)+hslb5(2)*Tk
     SLb(6)=hslb6(1)+hslb6(2)*Tk
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     entSL_=g_R*Tk*(EntropyMixing(yy(1),0)/3.+EntropyMixing(yy(2),0)/6.)
     call get_y1_y2(yy,0,y1_y2)
     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)
     gb2_=SLb(5)*y1_y2(5)+SLb(6)*y1_y2(6)
     gb_=(gb1_+gb2_+entSL_)/(1-yy(2)/6d0)
     Gibbs_FE_sol=gb_
  elseif(d_c.eq.0.and.d_T.eq.1)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(hslb1,Tk,1)
     SLb(2)=gibbs_T(hslb2,Tk,1)
     SLb(3)=gibbs_T(hslb3,Tk,1)
     SLb(4)=gibbs_T(hslb4,Tk,1)
     SLb(5)=hslb5(2)
     SLb(6)=hslb6(2)
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     entSL_=g_R*(EntropyMixing(yy(1),0)/3.+EntropyMixing(yy(2),0)/6.)
     call get_y1_y2(yy,0,y1_y2)
     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)
     gb2_=SLb(5)*y1_y2(5)+SLb(6)*y1_y2(6)
     gb_=(gb1_+gb2_+entSL_)/(1-yy(2)/6d0)
     Gibbs_FE_sol=gb_
  elseif(d_c.eq.0.and.d_T.eq.2)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(hslb1,Tk,2)
     SLb(2)=gibbs_T(hslb2,Tk,2)
     SLb(3)=gibbs_T(hslb3,Tk,2)
     SLb(4)=gibbs_T(hslb4,Tk,2)
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     call get_y1_y2(yy,0,y1_y2)
     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)
     gb_=gb1_/(1-yy(2)/6d0)
     Gibbs_FE_sol=gb_
  elseif(d_c.eq.1.and.d_T.eq.0)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(hslb1,Tk,0)
     SLb(2)=gibbs_T(hslb2,Tk,0)
     SLb(3)=gibbs_T(hslb3,Tk,0)
     SLb(4)=gibbs_T(hslb4,Tk,0)
     SLb(5)=hslb5(1)+hslb5(2)*Tk
     SLb(6)=hslb6(1)+hslb6(2)*Tk
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     yy_d(1)=y1(c,1)
     yy_d(2)=y2(c,1)

     entSL_=g_R*Tk*(EntropyMixing(yy(1),0)/3.+EntropyMixing(yy(2),0)/6.)
     entSL_d=g_R*Tk*(yy_d(1)*EntropyMixing(yy(1),1)/3.+yy_d(2)*EntropyMixing(yy(2),1)/6.)


     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1) !d(y1-y2)/d(y1)
     call get_y1_y2(yy,2,y1_y2_d2) !d(y1-y2)/d(y2)

     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)

     gb1_d1=SLb(1)*y1_y2_d1(1)+SLb(2)*y1_y2_d1(2)+SLb(3)*y1_y2_d1(3)+SLb(4)*y1_y2_d1(4)
     gb1_d1=gb1_d1*yy_d(1) !d(y1)/d(c)*d(gb1)/d(y1)

     gb1_d2=SLb(1)*y1_y2_d2(1)+SLb(2)*y1_y2_d2(2)+SLb(3)*y1_y2_d2(3)+SLb(4)*y1_y2_d2(4)
     gb1_d2=gb1_d2*yy_d(2)!d(y2)/d(c)*d(gb1)/d(y2)


     gb2_=SLb(5)*y1_y2(5)+SLb(6)*y1_y2(6)

     gb2_d1=SLb(5)*y1_y2_d1(5)+SLb(6)*y1_y2_d1(6)
     gb2_d1=gb2_d1*yy_d(1) !d(y1)/d(c)*d(gb2)/d(y1) 

     gb2_d2=SLb(5)*y1_y2_d2(5)+SLb(6)*y1_y2_d2(6)
     gb2_d2=gb2_d2*yy_d(2)!d(y2)/d(c)*d(gb2)/d(y2) 


     gb_=(gb1_+gb2_+entSL_)/(1-yy(2)/6d0)
     gb_d=(gb1_d1+gb1_d2+gb2_d1+gb2_d2+entSL_d)/(1-yy(2)/6d0)-(gb_)/(1-yy(2)/6d0)*(-yy_d(2)/6d0)
     Gibbs_FE_sol=gb_d
  elseif(d_c.eq.1.and.d_T.eq.1)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(hslb1,Tk,1)
     SLb(2)=gibbs_T(hslb2,Tk,1)
     SLb(3)=gibbs_T(hslb3,Tk,1)
     SLb(4)=gibbs_T(hslb4,Tk,1)
     SLb(5)=hslb5(2)
     SLb(6)=hslb6(2)
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     yy_d(1)=y1(c,1)
     yy_d(2)=y2(c,1)

     entSL_=g_R*(EntropyMixing(yy(1),0)/3.+EntropyMixing(yy(2),0)/6.)
     entSL_d=g_R*(yy_d(1)*EntropyMixing(yy(1),1)/3.+yy_d(2)*EntropyMixing(yy(2),1)/6.)


     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1) !d(y1-y2)/d(y1)
     call get_y1_y2(yy,2,y1_y2_d2) !d(y1-y2)/d(y2)

     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)

     gb1_d1=SLb(1)*y1_y2_d1(1)+SLb(2)*y1_y2_d1(2)+SLb(3)*y1_y2_d1(3)+SLb(4)*y1_y2_d1(4)
     gb1_d1=gb1_d1*yy_d(1) !d(y1)/d(c)*d(gb1)/d(y1)

     gb1_d2=SLb(1)*y1_y2_d2(1)+SLb(2)*y1_y2_d2(2)+SLb(3)*y1_y2_d2(3)+SLb(4)*y1_y2_d2(4)
     gb1_d2=gb1_d2*yy_d(2)!d(y2)/d(c)*d(gb1)/d(y2)


     gb2_=SLb(5)*y1_y2(5)+SLb(6)*y1_y2(6)

     gb2_d1=SLb(5)*y1_y2_d1(5)+SLb(6)*y1_y2_d1(6)
     gb2_d1=gb2_d1*yy_d(1) !d(y1)/d(c)*d(gb2)/d(y1) 

     gb2_d2=SLb(5)*y1_y2_d2(5)+SLb(6)*y1_y2_d2(6)
     gb2_d2=gb2_d2*yy_d(2)!d(y2)/d(c)*d(gb2)/d(y2) 


     gb_=(gb1_+gb2_+entSL_)/(1-yy(2)/6d0)
     gb_d=(gb1_d1+gb1_d2+gb2_d1+gb2_d2+entSL_d)/(1-yy(2)/6d0)-(gb_)/(1-yy(2)/6d0)*(-yy_d(2)/6d0)
     Gibbs_FE_sol=gb_d
  elseif(d_c.eq.1.and.d_T.eq.2)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(hslb1,Tk,2)
     SLb(2)=gibbs_T(hslb2,Tk,2)
     SLb(3)=gibbs_T(hslb3,Tk,2)
     SLb(4)=gibbs_T(hslb4,Tk,2)
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     yy_d(1)=y1(c,1)
     yy_d(2)=y2(c,1)

     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1) !d(y1-y2)/d(y1)
     call get_y1_y2(yy,2,y1_y2_d2) !d(y1-y2)/d(y2)

     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)

     gb1_d1=SLb(1)*y1_y2_d1(1)+SLb(2)*y1_y2_d1(2)+SLb(3)*y1_y2_d1(3)+SLb(4)*y1_y2_d1(4)
     gb1_d1=gb1_d1*yy_d(1) !d(y1)/d(c)*d(gb1)/d(y1)

     gb1_d2=SLb(1)*y1_y2_d2(1)+SLb(2)*y1_y2_d2(2)+SLb(3)*y1_y2_d2(3)+SLb(4)*y1_y2_d2(4)
     gb1_d2=gb1_d2*yy_d(2)!d(y2)/d(c)*d(gb1)/d(y2)

     gb_=gb1_/(1-yy(2)/6d0)
     gb_d=gb1_d1/(1-yy(2)/6d0)-(gb_)/(1-yy(2)/6d0)*(-yy_d(2)/6d0)
     Gibbs_FE_sol=gb_d

  else
     write(*,*)"unexpected d_c d_T combination",d_c,d_T
  endif
  
end function Gibbs_FE_sol
