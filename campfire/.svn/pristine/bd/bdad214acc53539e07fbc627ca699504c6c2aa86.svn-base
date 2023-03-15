!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  Utilities for phase field stencils in 2-d and 3-d
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
maximum = amaximum(p,n)
return
maximum = p(1)
j=1
do i = 2, n
   if(p(i).gt.p(j))then
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
double precision :: s,m,maximum
amaximum = 0d0
s=0.5
do i = 1, n
   amaximum = amaximum + ((((s+(1.0-s)*p(i))**8)**8)**8)
enddo
amaximum=exp(log(amaximum)/512)

amaximum=(amaximum-s)/(1.0-s)



end function amaximum

double precision function SurfaceEnergy(p,x,n)
use solution_parameters
implicit none
double precision, intent (in):: p(2,100),x(2)
integer, intent (in) :: n
integer i,j
double precision :: amaximum,maximum,q(100),qmax,xx,xhat(2)

!current extimate of max |grad(phi)| = 0.3
!xx = sqrt(x(1)**2 + x(2)**2 + 1d-16)

!xhat = x/xx

q=0d0
do j=1,n
   do i=1,2
      q(j) = q(j)  + p(i,j)*x(i)
   enddo
enddo

qmax=(epsilon+maximum(q,n))/(1d0+epsilon)

SurfaceEnergy =0.5*qmax**2
end function SurfaceEnergy

double precision function GradientEnergyHexPrismVol3(vbles,lp,dx)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision :: phi_(2),phi(2,2),A_(4),A(2,2),AA(2,2)
double precision :: yy(2),y,xx(2),x,d=4d-2,h
double precision :: p(2,6),pi=3.141592654,q(6),dh(2),pi6=.5235987758
integer :: i,j,np=6,k,l,iSurfaceEnergy
double precision :: SurfaceEnergy,AnisEnergy,UL(2),UR(2),DL(2),DR(2),ML(2),MR(2),MU(2),MD(2)
double precision, dimension(4,4)::stencil =reshape((/1.,-8.,8.,-1.,-8.,64.,-64.,8.,8.,-64.,64.,-8.,-1.,8.,-8.,1./),(/4,4/))
double precision, dimension(2,2)::stencil2 =reshape((/1.,-1.,-1.,1./),(/2,2/))
double precision :: cc(2)=(/-1,1/)


p= g_phex

GradientEnergyHexPrismVol3=0

   xx(1) = 0.5*(vbles(1,3,2)-vbles(1,1,2))/dx
   xx(2) = 0.5*(vbles(1,2,3)-vbles(1,2,1))/dx
   x = sqrt(xx(1)**2 + xx(2)**2)

   H = 0
   if(x.gt.1d-9)then
      xx=xx/x
      H = AnisEnergy(p,xx,6)**2*(vbles(1,1,2)+vbles(1,2,1)+vbles(1,3,2)+vbles(1,2,3)-4*vbles(1,2,2))/dx**2
   endif



   GradientEnergyHexPrismVol3 = - H





end function GradientEnergyHexPrismVol3
double precision function GradientEnergyHexPrismVol2(vbles,lp,dx)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision :: phi_(2),phi(2,2),A_(4),A(2,2),AA(2,2)
double precision :: yy(2),y,xx(2),zz(2),ww(2),x,d=4d-2,h
double precision :: p(2,6),pi=3.141592654,q(6),dh(2),pi6=.5235987758
integer :: i,j,np=6,k,l,iSurfaceEnergy
double precision :: SurfaceEnergy,AnisEnergy,UL(2),UR(2),DL(2),DR(2),ML(2),MR(2),MU(2),MD(2)
double precision, dimension(4,4)::stencil =reshape((/1.,-8.,8.,-1.,-8.,64.,-64.,8.,8.,-64.,64.,-8.,-1.,8.,-8.,1./),(/4,4/))
double precision, dimension(2,2)::stencil2 =reshape((/1.,-1.,-1.,1./),(/2,2/))
double precision :: cc(2)=(/-1,1/)
!!$do i=1,6
!!$p(1,i)=cos(pi*(i-1)/3d0)
!!$p(2,i)=sin(pi*(i-1)/3d0)
!!$enddo
!!$p(1,3)=cos(pi*0.8)
!!$p(2,3)=sin(0.8*pi)

p= g_phex

GradientEnergyHexPrismVol2=0


!   d=1d-5

!AnisEnergy(p,x,n)

   A_ = 0d0

   xx(1) = 0.5*(vbles(1,3,2)-vbles(1,1,2))/dx
   xx(2) = 0.5*(vbles(1,2,3)-vbles(1,2,1))/dx
   x = sqrt(xx(1)**2 + xx(2)**2)


   if(x.gt.1d-8)then
   d = 1d-3


   
!UR
   UR(1) = 0.5*(vbles(1,3,3)-vbles(1,2,3)+vbles(1,3,2)-vbles(1,2,2))/dx
   UR(2) = 0.5*(vbles(1,3,3)-vbles(1,3,2)+vbles(1,2,3)-vbles(1,2,2))/dx
   DR(1) = 0.5*(vbles(1,3,2)-vbles(1,2,2)+vbles(1,3,1)-vbles(1,2,1))/dx
   DR(2) = 0.5*(vbles(1,3,2)-vbles(1,3,1)+vbles(1,2,2)-vbles(1,2,1))/dx
   UL(1) = 0.5*(vbles(1,2,3)-vbles(1,1,3)+vbles(1,2,2)-vbles(1,1,2))/dx
   UL(2) = 0.5*(vbles(1,2,3)-vbles(1,2,2)+vbles(1,1,3)-vbles(1,1,2))/dx
   DL(1) = 0.5*(vbles(1,2,2)-vbles(1,1,2)+vbles(1,2,1)-vbles(1,1,1))/dx
   DL(2) = 0.5*(vbles(1,2,2)-vbles(1,2,1)+vbles(1,1,2)-vbles(1,1,1))/dx

   MR(1) = (vbles(1,3,2)-vbles(1,2,2))/dx
   MR(2) =  0.25*(vbles(1,3,3)-vbles(1,3,1)+vbles(1,2,3)-vbles(1,2,1))/dx
   ML(1) = (vbles(1,2,2)-vbles(1,1,2))/dx
   ML(2) =  0.25*(vbles(1,2,3)-vbles(1,2,1)+vbles(1,1,3)-vbles(1,1,1))/dx
   MU(2) = (vbles(1,2,3)-vbles(1,2,2))/dx
   MU(1) =  0.25*(vbles(1,3,3)-vbles(1,1,3)+vbles(1,3,2)-vbles(1,1,2))/dx
   MD(2) = (vbles(1,2,2)-vbles(1,2,1))/dx
   MD(1) =  0.25*(vbles(1,3,2)-vbles(1,1,2)+vbles(1,3,1)-vbles(1,1,1))/dx

   yy(1) = UR(1) + d
   yy(2) = UR(2)
   xx(1) = UR(1) - d
   xx(2) = UR(2)
   zz(1) = UR(1) +2*d
   zz(2) = UR(2)
   ww(1) = UR(1) -2*d
   ww(2) = UR(2)
   A_(1) = A_(1) + AnisEnergy(p,UR,np)*0.5/6.*(AnisEnergy(p,ww,np)-8* AnisEnergy(p,xx,np)+8* AnisEnergy(p,yy,np) -AnisEnergy(p,zz,np))/d/12d0 
!   A_(1) = A_(1) + AnisEnergy(p,UR,np)*0.5/6.*(AnisEnergy(p,yy,np) - AnisEnergy(p,xx,np))/d 


   yy(1) = MR(1) + d
   yy(2) = MR(2)
   xx(1) = MR(1) - d
   xx(2) = MR(2)
   zz(1) = UR(1) +2*d
   zz(2) = UR(2)
   ww(1) = UR(1) -2*d
   ww(2) = UR(2)

   A_(1) = A_(1) + AnisEnergy(p,MR,np)*1./3.*(AnisEnergy(p,ww,np)-8* AnisEnergy(p,xx,np)+8* AnisEnergy(p,yy,np) -AnisEnergy(p,zz,np))/d/12d0 
!LR
   yy(1) = DR(1) + d
   yy(2) = DR(2) 
   xx(1) = DR(1) - d
   xx(2) = DR(2)
   zz(1) = UR(1) +2*d
   zz(2) = UR(2)
   ww(1) = UR(1) -2*d
   ww(2) = UR(2)
 
   A_(1) = A_(1) + AnisEnergy(p,DR,np)*0.5/6.*(AnisEnergy(p,ww,np)-8* AnisEnergy(p,xx,np)+8* AnisEnergy(p,yy,np) -AnisEnergy(p,zz,np))/d/12d0 

!Aright done
!UL

   yy(1) = UL(1) + d
   yy(2) = UL(2) 
   xx(1) = UL(1) - d
   xx(2) = UL(2)
   zz(1) = UR(1) +2*d
   zz(2) = UR(2)
   ww(1) = UR(1) -2*d
   ww(2) = UR(2) 
   A_(2) = A_(2) + AnisEnergy(p,UL,np)*0.5/6.*(AnisEnergy(p,ww,np)-8* AnisEnergy(p,xx,np)+8* AnisEnergy(p,yy,np) -AnisEnergy(p,zz,np))/d/12d0 


   yy(1) = ML(1) + d
   yy(2) = ML(2) 
   xx(1) = ML(1) - d
   xx(2) = ML(2) 
   zz(1) = UR(1) +2*d
   zz(2) = UR(2)
   ww(1) = UR(1) -2*d
   ww(2) = UR(2) 
   A_(2) = A_(2) + AnisEnergy(p,ML,np)*1./3.*(AnisEnergy(p,ww,np)-8* AnisEnergy(p,xx,np)+8* AnisEnergy(p,yy,np) -AnisEnergy(p,zz,np))/d/12d0 
!DL

   yy(1) = DL(1) + d
   yy(2) = DL(2) 
   xx(1) = DL(1) - d
   xx(2) = DL(2)
   zz(1) = UR(1) +2*d
   zz(2) = UR(2)
   ww(1) = UR(1) -2*d
   ww(2) = UR(2)
 
   A_(2) = A_(2) + AnisEnergy(p,DL,np)*0.5/6.*(AnisEnergy(p,ww,np)-8* AnisEnergy(p,xx,np)+8* AnisEnergy(p,yy,np) -AnisEnergy(p,zz,np))/d/12d0 
!Aleft done
!UR
   yy(1) = UR(1) 
   yy(2) = UR(2) + d
   xx(1) = UR(1) 
   xx(2) = UR(2) - d
   zz(1) = UR(1) 
   zz(2) = UR(2) + 2*d
   ww(1) = UR(1) 
   ww(2) = UR(2) - 2*d

   A_(3) = A_(3) + AnisEnergy(p,UR,np)*0.5/6.*(AnisEnergy(p,ww,np)-8* AnisEnergy(p,xx,np)+8* AnisEnergy(p,yy,np) -AnisEnergy(p,zz,np))/d/12d0 


   yy(1) = MU(1)
   yy(2) = MU(2) + d
   xx(1) = MU(1)
   xx(2) = MU(2) - d
   zz(1) = UR(1) 
   zz(2) = UR(2) + 2*d
   ww(1) = UR(1) 
   ww(2) = UR(2) - 2*d

   A_(3) = A_(3) + AnisEnergy(p,MU,np)*1./3.*(AnisEnergy(p,ww,np)-8* AnisEnergy(p,xx,np)+8* AnisEnergy(p,yy,np) -AnisEnergy(p,zz,np))/d/12d0 
!UL
   yy(1) = UL(1)
   yy(2) = UL(2) + d 
   xx(1) = UL(1)
   xx(2) = UL(2) - d 
   zz(1) = UR(1) 
   zz(2) = UR(2) + 2*d
   ww(1) = UR(1) 
   ww(2) = UR(2) - 2*d

   A_(3) = A_(3) + AnisEnergy(p,UL,np)*0.5/6.*(AnisEnergy(p,ww,np)-8* AnisEnergy(p,xx,np)+8* AnisEnergy(p,yy,np) -AnisEnergy(p,zz,np))/d/12d0 

!Aup done
   yy(1) = DR(1) 
   yy(2) = DR(2) + d
   xx(1) = DR(1) 
   xx(2) = DR(2) - d
   zz(1) = UR(1) 
   zz(2) = UR(2) + 2*d
   ww(1) = UR(1) 
   ww(2) = UR(2) - 2*d
   A_(4) = A_(4) + AnisEnergy(p,DR,np)*0.5/6.*(AnisEnergy(p,ww,np)-8* AnisEnergy(p,xx,np)+8* AnisEnergy(p,yy,np) -AnisEnergy(p,zz,np))/d/12d0 


   yy(1) = MD(1)
   yy(2) = MD(2) + d
   xx(1) = MD(1)
   xx(2) = MD(2) - d
   zz(1) = UR(1) 
   zz(2) = UR(2) + 2*d
   ww(1) = UR(1) 
   ww(2) = UR(2) - 2*d
   A_(4) = A_(4) + AnisEnergy(p,MD,np)*1./3.*(AnisEnergy(p,ww,np)-8* AnisEnergy(p,xx,np)+8* AnisEnergy(p,yy,np) -AnisEnergy(p,zz,np))/d/12d0 


   yy(1) = DL(1)
   yy(2) = DL(2) + d 
   xx(1) = DL(1)
   xx(2) = DL(2) - d 
   zz(1) = UR(1) 
   zz(2) = UR(2) + 2*d
   ww(1) = UR(1) 
   ww(2) = UR(2) - 2*d

   A_(4) = A_(4) + AnisEnergy(p,DL,np)*0.5/6.*(AnisEnergy(p,ww,np)-8* AnisEnergy(p,xx,np)+8* AnisEnergy(p,yy,np) -AnisEnergy(p,zz,np))/d/12d0 
!Adown done




   H = (A_(1) - A_(2) +A_(3) - A_(4))/dx


endif

   GradientEnergyHexPrismVol2 = - H





end function GradientEnergyHexPrismVol2




double precision function AnisEnergy(p,x,n)
use solution_parameters
implicit none
double precision, intent (in):: p(2,100),x(2)
integer, intent (in) :: n
integer i,j
double precision :: maximum,q(100),qmax



q=0d0
do j=1,n
   do i=1,2
      q(j) = q(j)+  p(i,j)*x(i)
   enddo
enddo

qmax=(sqrt(x(1)**2+x(2)**2)+epsilon*maximum(q,n))/(1d0+epsilon)

AnisEnergy =qmax
end function AnisEnergy



!!!!!!!!

double precision function GradientEnergyHexPrism(vbles,lp,dx)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision :: phi_(2),phi(2,2),A_(2),A(2,2),AA(2,2)
double precision :: yy(2),y,xx(2),x,d=4d-2,h
double precision :: p(2,6),pi=3.141592654,q(6),dh(2),pi6=.5235987758
integer :: i,j,np=6,k,l,iSurfaceEnergy
double precision :: SurfaceEnergy
double precision, dimension(4,4)::stencil =reshape((/1.,-8.,8.,-1.,-8.,64.,-64.,8.,8.,-64.,64.,-8.,-1.,8.,-8.,1./),(/4,4/))
double precision, dimension(2,2)::stencil2 =reshape((/1.,-1.,-1.,1./),(/2,2/))
double precision :: cc(2)=(/-1,1/)
do i=1,6
p(1,i)=cos(pi*(i-1)/3d0)
p(2,i)=sin(pi*(i-1)/3d0)
enddo
p(1,3)=cos(0.8*pi)
p(2,3)=sin(0.8*pi)


GradientEnergyHexPrism=0
   phi_(1)=0.5*(vbles(1,3,2)-vbles(1,1,2))/dx
   phi_(2)=0.5*(vbles(1,2,3)-vbles(1,2,1))/dx
   phi(1,1)=(vbles(1,3,2)-2d0*vbles(1,2,2)+vbles(1,1,2))/dx**2
   phi(2,2)=(vbles(1,2,3)-2d0*vbles(1,2,2)+vbles(1,2,1))/dx**2


   phi(1,2) = 0.25*(vbles(1,3,3)+vbles(1,1,1)-vbles(1,3,1)-vbles(1,1,3))/dx**2
   phi(2,1) = phi(1,2)

   do i=1,2
      xx(i) = phi_(i)
   enddo


   x = sqrt(xx(1)**2+xx(2)**2+1d-9)

!   xx=xx/x


   d=1d-1*x



   A=0d0
   do i=1,2
!!$      dh=0d0
!!$      dh(i)=d
!!$      yy=xx-2*dh
!!$      A(i,i) =         -SurfaceEnergy(p,yy,np)
!!$      yy=xx-dh
!!$      A(i,i) = A(i,i) + 16*SurfaceEnergy(p,yy,np)
!!$      A(i,i) = A(i,i) - 30*SurfaceEnergy(p,xx,np)
!!$      yy=xx+dh
!!$      A(i,i) = A(i,i) + 16*SurfaceEnergy(p,yy,np)
!!$      yy=xx+2*dh
!!$      A(i,i) = A(i,i) - SurfaceEnergy(p,yy,np)

      dh=0d0
      dh(i)=d
      yy=xx-dh
      A(i,i) = SurfaceEnergy(p,yy,np)
      A(i,i) = A(i,i) -2d0*SurfaceEnergy(p,xx,np)
      yy=xx+dh
      A(i,i) = A(i,i) + SurfaceEnergy(p,yy,np)
      A(i,i) = A(i,i)/d**2
   enddo

!!$
!!$   i=1
!!$   j=2
!!$   dh=0d0
!!$   dh(i)=d
!!$   dh(j)=d
!!$   yy=xx + dh
!!$   A(i,j) = SurfaceEnergy(p,yy,np)
!!$   
!!$   yy=xx - dh
!!$   A(i,j) = A(i,j) + SurfaceEnergy(p,yy,np)
!!$   
!!$   dh=0d0
!!$   dh(i)=d
!!$   dh(j)=-d
!!$   yy=xx+dh
!!$   A(i,j) = A(i,j) - SurfaceEnergy(p,yy,np)
!!$   
!!$   yy=xx-dh
!!$   A(i,j) = A(i,j) - SurfaceEnergy(p,yy,np)
!!$   
!!$   A(i,j)=A(i,j)*0.25/d**2
!!$   A(j,i)=A(i,j)
   
!!!!!!!!!!!!!
   A(1,2)=0d0
   do i=1,2
      do j=1,2
         yy(1)=xx(1)+cc(i)*d
         yy(2)=xx(2)+cc(j)*d
         A(1,2) = A(1,2) + stencil2(i,j)*SurfaceEnergy(p,yy,np)
      enddo
   enddo
   A(1,2)=A(1,2)*0.25d0/d**2
   A(2,1)=A(1,2)

   

   H = 0d0
   do i=1,2
      do j=1,2
         H = H + phi(i,j)*A(i,j) 
      enddo
   enddo


   GradientEnergyHexPrism = - H





end function GradientEnergyHexPrism
!!!

!!!!!!!!

double precision function GradientEnergyHexPrismVol(vbles,lp,dx)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision :: phi_(2),phi(2,2),A_(4),A(2,2),AA(2,2)
double precision :: yy(2),y,xx(2),x,d=4d-2,h
double precision :: p(2,6),pi=3.141592654,q(6),dh(2),pi6=.5235987758
integer :: i,j,np=6,k,l,iSurfaceEnergy
double precision :: SurfaceEnergy,AnisEnergy
double precision, dimension(4,4)::stencil =reshape((/1.,-8.,8.,-1.,-8.,64.,-64.,8.,8.,-64.,64.,-8.,-1.,8.,-8.,1./),(/4,4/))
double precision, dimension(2,2)::stencil2 =reshape((/1.,-1.,-1.,1./),(/2,2/))
double precision :: cc(2)=(/-1,1/)
!!$do i=1,6
!!$p(1,i)=cos(pi*(i-1)/3d0)
!!$p(2,i)=sin(pi*(i-1)/3d0)
!!$enddo
!!$p(1,3)=cos(pi*0.8)
!!$p(2,3)=sin(0.8*pi)

p= g_phex

GradientEnergyHexPrismVol=0


!   d=1d-5

!AnisEnergy(p,x,n)

   A_ = 0d0

   d = 0.003


   xx(1) = (vbles(1,3,2)-vbles(1,2,2))/dx
   xx(2) =  0.25*(vbles(1,3,3)-vbles(1,3,1)+vbles(1,2,3)-vbles(1,2,1))/dx
   
!   x = sqrt(xx(1)**2+xx(2)**2+1d-16)


   yy(1) = xx(1) + d
   yy(2) = xx(2) 
   A_(1) = A_(1) + (SurfaceEnergy(p,yy,np) - SurfaceEnergy(p,xx,np))/d 



   xx(1) = (vbles(1,2,2)-vbles(1,1,2))/dx
   xx(2) =  0.25*(vbles(1,1,3)-vbles(1,1,1)+vbles(1,2,3)-vbles(1,2,1))/dx
   yy(1) = xx(1) + d
   yy(2) = xx(2) 
   A_(2) = A_(2) + (SurfaceEnergy(p,yy,np) - SurfaceEnergy(p,xx,np))/d 



   xx(2) = (vbles(1,2,3)-vbles(1,2,2))/dx
   xx(1) =  0.25*(vbles(1,3,3)-vbles(1,1,3)+vbles(1,3,2)-vbles(1,1,2))/dx
   yy(1) = xx(1) 
   yy(2) = xx(2) + d 
   A_(3) = A_(3) + (SurfaceEnergy(p,yy,np) - SurfaceEnergy(p,xx,np))/d 




   xx(2) = (vbles(1,2,2)-vbles(1,2,1))/dx
   xx(1) =  0.25*(vbles(1,3,1)-vbles(1,1,1)+vbles(1,3,2)-vbles(1,1,2))/dx
   yy(1) = xx(1) 
   yy(2) = xx(2) + d
   A_(4) = A_(4) + (SurfaceEnergy(p,yy,np) - SurfaceEnergy(p,xx,np))/d 



   H = (A_(1) - A_(2) +A_(3) - A_(4))/dx




   GradientEnergyHexPrismVol = - H





end function GradientEnergyHexPrismVol
!!!




!Facet anisotropy
double precision function GradientEnergy_circ(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision ::phix,phiy
double precision :: M11,M12,M22,g,q,qq,kappa,dirac,ddx
double precision :: g0,h0,dh,gx,hx,gxx,hxx,gy,hy,gyy,hyy,hxy,phixx,phiyy,phixy,f,fy,fx,fxx,fxy,fyy,f0,p,px,pxx
double precision :: X1,X2,theta,pi,A,Ax,Ay,A2,W
integer :: i,imodel=1





phixx = (vbles(1,3,2) - 2d0*vbles(1,2,2) + vbles(1,1,2))/dx**2
phiyy = (vbles(1,2,3) - 2d0*vbles(1,2,2) + vbles(1,2,1))/dx**2







if(g_kinetic.and. (.not. g_needle))then

   GradientEnergy_circ = -phixx - phiyy
   return
else

   
   phix = 0.5*(vbles(1,3,2) - vbles(1,1,2))/dx 
   phiy = 0.5*(vbles(1,2,3) - vbles(1,2,1))/dx 
   phixy = 0.25*(vbles(1,3,3)+vbles(1,1,1)-vbles(1,3,1)-vbles(1,1,3))/dx**2
   
   
   g0  = sqrt(phix*phix +phiy*phiy)
   if(g0.gt.1d-8)then
      
      
      
      
      X1 = max(phix/g0,1d-9)
      X2 = max(phiy/g0,1d-9)
      
      
      
      
      theta = atan2(X2,X1)
      pi = 4.*atan2(1.,1.)
      
      
      !      A = (1d0 + epsilon*abs(4*X1**3-3*X1))**2
!!$      if(abs(theta).le.pi/6.)then
!!$         A = 1./(1.+epsilon)+2.*epsilon*X1**2./(1.+epsilon)
!!$         A2 = 4.*epsilon*X2**2./(1.+epsilon)-4.*epsilon*X1**2./(1.+epsilon)
!!$      else
!!$         A = (2.*epsilon*sqrt(3.)*X1*X2-2.*epsilon*X1**2+3*epsilon+2.)/(2.+2.*epsilon)
!!$         A2 = -4.*epsilon*(sqrt(3.)*X1*X2-X1**2+0.5)/(1+epsilon)
!!$      endif
      if(imodel.eq. 3)then
!         A = 1+2*sin(theta)*cos(theta)
!         A2 = -8*sin(theta)*cos(theta)

         A = sin(2*theta) + 1
         A2 = -4*sin(2*theta)

         kappa = (phixx*X2**2 - 2.* phixy*X1*X2 + phiyy*X1**2)
         h0 = 0.01
         dirac = 0d0

         if(abs(theta).lt.h0)then
            dh = theta
            dirac = 15./16.*(1-dh/h0)**2*(1+dh/h0)**2/h0*4.
         elseif(abs(theta-pi/2.).lt.h0)then
            dh = theta-pi/2.
            dirac = 15./16.*(1-dh/h0)**2*(1+dh/h0)**2/h0*4.
         endif

         A2 = A2 + dirac
 
      elseif(imodel.eq.1)then
!hexagon
         if(abs(theta).le.pi/6.)then
!            A = epsilon + 2*X1**2
!            A2 =-4*X1**2+4*X2**2

            A = 2*cos(theta)**2
            A2 = 4-8*cos(theta)**2

         else
            A = 2*cos(theta-pi/3.)**2
            A2 = 4-8*cos(theta-pi/3.)**2

!            A = epsilon + 1.5-X1**2+sqrt(3.)*X1*X2
!            A2 = 2*X1**2-2*X2**2-4*sqrt(3.)*X1*X2
         endif
         
         
         kappa = (phixx*X2**2 - 2.* phixy*X1*X2 + phiyy*X1**2)
         h0 = 0.02
         dirac = 0d0
         
         if(abs(theta-0.5*pi).lt.h0)then
            dh = theta-0.5*pi
            dirac = 15./16.*(1-dh/h0)**2*(1+dh/h0)**2/h0*2.*sqrt(3.)
         elseif(abs(theta-pi/6.).lt.h0)then
            dh = theta-pi/6.
            dirac = 15./16.*(1-dh/h0)**2*(1+dh/h0)**2/h0*2.*sqrt(3.)   
         endif
   
         A2 = A2 + dirac
      elseif(imodel.eq.2)then
         epsilon=0.125
         A = (1d0 + epsilon*abs(4*X1**3-3*X1))**2
         kappa = (phixx*X2**2 - 2.* phixy*X1*X2 + phiyy*X1**2)
         h0 = 0.02
         dirac = 0d0
         
         if(abs(theta-0.5*pi).lt.h0)then
            dh = theta-0.5*pi
            dirac = 15./16.*(1-dh/h0)**2*(1+dh/h0)**2/h0*12d0*epsilon
           !   A = A + dh**6/(32*h0**5)-5*dh**4/(32*h0**3)+15*dh**2/(32*h0) -11*h0**2/32
         elseif(abs(theta-pi/6.).lt.h0)then
            dh = theta-pi/6.
            dirac = 15./16.*(1-dh/h0)**2*(1+dh/h0)**2/h0*12d0*epsilon
         endif
         
         A2= 18*epsilon*(epsilon*((1d0 - 18*X1**2 + 48*X1**4 - 32*X1**6)) - abs(4*X1**3-3*X1)) + dirac
         A2 = A2 + dirac
      else
         write(*,*)"what anisotropy model"
         stop
      endif
      GradientEnergy_circ = (phixx+phiyy)*A + 0.5*A2*kappa 
      
   else
      GradientEnergy_circ = phixx + phiyy
   endif
   
endif






GradientEnergy_circ = - GradientEnergy_circ




end function GradientEnergy_circ

!Facet anisotropy
double precision function GradientEnergy(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision ::phix,phiy
double precision :: M11,M12,M22,g,q,qq
double precision :: g0,h0,gx,hx,gxx,hxx,gy,hy,gyy,hyy,hxy,phixx,phiyy,phixy,f,fy,fx,fxx,fxy,fyy,f0,p,px,pxx
double precision :: v,vx,vxx,w,wx,wy,wxx,wyy,wxy





phixx = (vbles(1,3,2) - 2d0*vbles(1,2,2) + vbles(1,1,2))/dx**2
phiyy = (vbles(1,2,3) - 2d0*vbles(1,2,2) + vbles(1,2,1))/dx**2








if(g_kinetic.and. (.not. g_needle))then

   GradientEnergy = -phixx - phiyy
   return
else
   q=0.5/g_lambda!g_q0
   
   phix = 0.5*(vbles(1,3,2) - vbles(1,1,2))/dx 
   phiy = 0.5*(vbles(1,2,3) - vbles(1,2,1))/dx 
   phixy = 0.25*(vbles(1,3,3)+vbles(1,1,1)-vbles(1,3,1)-vbles(1,1,3))/dx**2
   
   g0  = phix*phix +phiy*phiy +q**2
   gx = 2d0*phix
   gy = 2d0*phiy
   gxx=2d0
   gyy=2d0




   f   = phix**3-3*phix*phiy**2
   fx  = 3*phix**2 - 3*phiy**2
   fxx = 6*phix
   fy  = -6*phix*phiy
   fyy = -6*phix
   fxy = -6*phiy






   h0 = f**2 + q**6
   hx = 2*f*fx
   hxx = 2*fx**2+2*f*fxx
   hy = 2*f*fy
   hyy = 2*fy**2+2*f*fyy
   hxy = 2*fy*fx+2*f*fxy

   
   m11 = -(-2*h0*g0**2*hxx+g0**2*hx**2+2*h0*g0*hx*gx+2*h0**2*g0*gxx-3*h0**2*gx**2)/h0/g0**3/sqrt(h0/g0)*epsilon/3  
   m11 = m11 + gxx/2
   m22 = -(-2*h0*g0**2*hyy+g0**2*hy**2+2*h0*g0*gy*hy+2*h0**2*g0*gyy-3*h0**2*gy**2)/h0/g0**3/sqrt(h0/g0)*epsilon/3
   m22 = m22 + gyy/2
   m12 = -(-2*h0*g0**2*hxy+g0**2*hx*hy+h0*g0*gx*hy+h0*g0*gy*hx-3*h0**2*gx*gy)/h0/g0**3/sqrt(h0/g0)*epsilon/3 




   GradientEnergy = phixx*m11+2d0*phixy*m12+phiyy*m22
endif






GradientEnergy = - GradientEnergy




end function GradientEnergy


double precision function GradientEnergyx(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision ::phix,phiy,X,Y,u,v,A
double precision :: G11,G22,G12,M11,M12,g,GG
double precision :: LapPhi(1)




call Calc_GradPhi(LapPhi(1),vbles,dx)

phix = 0.5*(vbles(1,3,2) - vbles(1,1,2))/dx 
phiy = 0.5*(vbles(1,2,3) - vbles(1,2,1))/dx 
u = phix*phix
v = phiy*phiy
g  = u + v
if(g.gt.1d-20)then

   X = u/g
   Y= v/g
   A = A_0*(1d0+epsilon_tilde*(X**2+Y**2))

   GG = 16*A_0*epsilon_tilde*X*Y*(A_0*epsilon_tilde*(X-Y)**2+A)
   G11 = Y*(GG +4*epsilon_tilde*A*A_0*(X-Y))+ A**2
   G22 = X*(GG +4*epsilon_tilde*A*A_0*(Y-X))+ A**2
   G12 = -phix*phiy*GG/g



   M11 = 0.5*(vbles(1,3,2)+vbles(1,1,2)-vbles(1,2,3)-vbles(1,2,1))/dx**2  !(phi_xx - phi_yy)/2
   M12 = 0.25*(vbles(1,3,3)+vbles(1,1,1)-vbles(1,1,3)-vbles(1,3,1))/dx**2 !phi_xy

   GradientEnergyx = 0.5*LapPhi(1)*(G11+G22)+M11*(G11-G22)+2*M12*G12


else
   GradientEnergyx = (A_0*(1d0+epsilon_tilde))**2*LapPhi(1)
endif

GradientEnergyx = - GradientEnergyx

!
end function GradientEnergyx

double precision function GradientEnergy_2theta(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision ::phix,phiy,X,Y,u,v,A,Ax,Ay,Axy,Axx,Ayy,phixy,phixx,phiyy,fxx,fxy,fyy
double precision :: G11,G22,G12,M11,M12,g,GG,q=0.1
double precision :: LapPhi(1), s1,s2,phi,potential,Omega,c




!call Calc_GradPhi(LapPhi(1),vbles,dx)
phi = vbles(1,2,2)
c = vbles(3,2,2)
phix  = 0.5*(vbles(1,3,2) - vbles(1,1,2))/dx 
phix=phix*g_r_factor
phiy  = 0.5*(vbles(1,2,3) - vbles(1,2,1))/dx 
phixy = 0.25*(vbles(1,3,3)+vbles(1,1,1)-vbles(1,1,3)-vbles(1,3,1))/dx**2 
phixy=phixy*g_r_factor
phixx = (vbles(1,1,2)-2*vbles(1,2,2)+vbles(1,3,2))/dx**2
phixx=phixx*g_r_factor**2
phiyy = (vbles(1,2,1)-2*vbles(1,2,2)+vbles(1,2,3))/dx**2

LapPhi(1)=phixx+phiyy
u = phix*phix
v = phiy*phiy
g  = u + v
if(g.gt.1d-20)then


!!$   A   = 1d0+epsilon*(2*u/g-1)
!!$   Ax  = epsilon*4*phix*phiy**2/g**2
!!$   Ay  = -epsilon*4*phix**2*phiy/g**2
!!$   Axy = epsilon*8*phix*phiy*(phix**2-phiy**2)/g**3
!!$   Axx = -epsilon*4*phiy**2*(3*phix**2-phiy**2)/g**3
!!$   Ayy = -epsilon*4*phix**2*(phix**2-3*phiy**2)/g**3

A  =  1+epsilon*sqrt(((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**3)
Ax = epsilon/sqrt(((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**3)*(2*(phix**3-3*phix*phiy**2)*(3*phix**2-3*phiy**2)/(phix**2+phiy**2+q**2)**3-6*((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**4*phix)/2
Ay = epsilon/sqrt(((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**3)*(-12*(phix**3-3*phix*phiy**2)*phix*phiy/(phix**2+phiy**2+q**2)**3-6*((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**4*phiy)/2
Axx = -epsilon/sqrt(((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**3)**3*(2*(phix**3-3*phix*phiy**2)*(3*phix**2-3*phiy**2)/(phix**2+phiy**2+q**2)**3-6*((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**4*phix)**2/4+epsilon/sqrt(((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**3)*(2*(3*phix**2-3*phiy**2)**2/(phix**2+phiy**2+q**2)**3+12*(phix**3-3*phix*phiy**2)*phix/(phix**2+phiy**2+q**2)**3-24*(phix**3-3*phix*phiy**2)*(3*phix**2-3*phiy**2)/(phix**2+phiy**2+q**2)**4*phix+48*((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**5*phix**2-6*((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**4)/2
s1 = -epsilon/sqrt(((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**3)**3*(2*(phix**3-3*phix*phiy**2)*(3*phix**2-3*phiy**2)/(phix**2+phiy**2+q**2)**3-6*((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**4*phix)*(-12*(phix**3-3*phix*phiy**2)*phix*phiy/(phix**2+phiy**2+q**2)**3-6*((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**4*phiy)/4
s2 = epsilon/sqrt(((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**3)*(-12*phix*phiy*(3*phix**2-3*phiy**2)/(phix**2+phiy**2+q**2)**3-12*(phix**3-3*phix*phiy**2)*phiy/(phix**2+phiy**2+q**2)**3-12*(phix**3-3*phix*phiy**2)*(3*phix**2-3*phiy**2)/(phix**2+phiy**2+q**2)**4*phiy+72*(phix**3-3*phix*phiy**2)*phix**2*phiy/(phix**2+phiy**2+q**2)**4+48*((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**5*phix*phiy)/2
Axy  = s1+s2
Ayy =  -epsilon/sqrt(((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**3)**3*(-12*(phix**3-3*phix*phiy**2)*phix*phiy/(phix**2+phiy**2+q**2)**3-6*((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**4*phiy)**2/4+epsilon/sqrt(((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**3)*(72*phix**2*phiy**2/(phix**2+phiy**2+q**2)**3-12*(phix**3-3*phix*phiy**2)*phix/(phix**2+phiy**2+q**2)**3+144*(phix**3-3*phix*phiy**2)*phix*phiy**2/(phix**2+phiy**2+q**2)**4+48*((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**5*phiy**2-6*((phix**3-3*phix*phiy**2)**2+q**6)/(phix**2+phiy**2+q**2)**4)/2


if(g_equal_delta)then
   write(*,*)"g_equal_delta"
   stop
   Omega  = (phi**2*(1-phi)**2)/g_lambda**2
   fxx = (Ax*Ax+A*Axx)*(g+2*Omega)+2*A*Ax*phix+2*A*Ax*phix ! + A**2
   fyy = (Ay*Ay+A*Ayy)*(g+2*Omega)+2*A*Ay*phiy+2*A*Ay*phiy
   fxy = (Ax*Ay+A*Axy)*(g+2*Omega)+2*A*Ax*phiy+2*A*Ay*phix ! + A**2
else
   fxx = (Ax*Ax+A*Axx)*g+2*A*Ax*phix+2*A*Ax*phix ! + A**2
   fyy = (Ay*Ay+A*Ayy)*g+2*A*Ay*phiy+2*A*Ay*phiy
   fxy = (Ax*Ay+A*Axy)*g+2*A*Ax*phiy+2*A*Ay*phix ! + A**2
endif

   GradientEnergy_2theta = phixx*fxx+2*phixy*fxy +phiyy*fyy + A**2*LapPhi(1) - A**2*phi*(1-phi)*(1-2*phi)/g_lambda**2

else
   GradientEnergy_2theta = LapPhi(1)
endif

GradientEnergy_2theta = - GradientEnergy_2theta

!
end function GradientEnergy_2theta


double precision function GradientEnergy_ellipse(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision ::phix,phiy,X,Y,u,v,A,Ax,Ay,Axy,Axx,Ayy,phixy,phixx,phiyy,fxx,fxy,fyy
double precision :: G11,G22,G12,M11,M12,g,GG
double precision :: LapPhi(1)




call Calc_GradPhi(LapPhi(1),vbles,dx)

phix  = 0.5*(vbles(1,3,2) - vbles(1,1,2))/dx 
phiy  = 0.5*(vbles(1,2,3) - vbles(1,2,1))/dx 
phixy = 0.25*(vbles(1,3,3)+vbles(1,1,1)-vbles(1,1,3)-vbles(1,3,1))/dx**2 
phixx = (vbles(1,1,2)-2*vbles(1,2,2)+vbles(1,3,2))/dx**2
phiyy = (vbles(1,2,1)-2*vbles(1,2,2)+vbles(1,2,3))/dx**2

   GradientEnergy_ellipse = -phixx/(1d0-epsilon)**2 - phiyy  




!
end function GradientEnergy_ellipse







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
double precision, intent (in):: dx,c,T,phi(N_phases),c_dot,phi_dot(N_phases),vbles(N_unk+1,3,3)

double precision, dimension(3,3)::stencil =reshape((/0.,1.,0.,1.,-4.,1.,0.,1.,0./),(/3,3/))
double precision, dimension(3)::stencil_x =(/-.5,.0,.5/)
double precision:: FreeEnergy,GradientEnergyx,potential,x,a,xa,dFdc(3,3),d_H(3,3) !x coordinate
double precision:: Cp= 28.5!29.9758394     !J/mol/K
double precision :: D_heat_c 
double precision :: T_,c_,phi_(N_phases),MolPerVol,X_energy,Y_energy,div,Grad_DW
integer :: ii,jj,kk


if(g_pure)then
  MolPerVol = g_MolPerVol(1)
else
  MolPerVol = g_MolPerVol(1)*(1.-c)+g_MolPerVol(2)*c  !mol/m^3
endif


! to help with floating point arithmetic
  X_energy=g_R*g_T0  !J/mol



  TemperatureRHS=0.

!  TemperatureRHS = FreeEnergy(c,T,phi,c_dot,phi_dot,N_phases+1,.false.) !J K/mol
if(g_pure)then
  TemperatureRHS = FreeEnergy(g_c0,T,phi,0d0,phi_dot,N_phases+1,g_approx) !J K/mol  Using quad approximation
else
  TemperatureRHS = FreeEnergy(c,T,phi,c_dot,phi_dot,N_phases+1,g_approx) !J K/mol  Using quad approximation
endif



  Cp = g_C  !J/K/mol

  Y_energy=Cp*g_T0  !J/mol  
  TemperatureRHS = TemperatureRHS/Y_energy  !no dimension

  
  do ii=1,3
     do jj=1,3
        T_=vbles(N_phases+1,ii,jj)*stencil(ii,jj)
        TemperatureRHS = TemperatureRHS + g_D_tilde*T_/(dx*dx)
     enddo
  enddo
  
! div(D_heat_c grad(dFdc))
if(.not.g_pure)then
  if(heat_c_term)then
     do ii=1,3
        do jj=1,3
           T_=vbles(2,ii,jj)
           c_=vbles(3,ii,jj)
           phi_(1)=vbles(1,ii,jj)
           dFdc(ii,jj)=(FreeEnergy(c_,T_,Phi_,c_dot,phi_dot,3,g_approx)/X_energy)**2

           !using quad approximation and tables
!           dFdc(ii,jj)=(FreeEnergy(c_,T_,Phi_,c_dot,phi_dot,3,.false.)/X_energy)**2
           D_heat_c =  g_D(1)*phi(1)+(1d0-phi(1))*g_D(2)
           d_H(ii,jj) = 0.5*c_*(1.-c_)*(D_heat_c/g_Dch)*(X_energy/Y_energy)   
        enddo
     enddo
     div=0d0
     div = div  + (dFdc(3,2)-dFdc(2,2))*(d_H(3,2)+d_H(2,2)) 
     div = div  + (dFdc(1,2)-dFdc(2,2))*(d_H(1,2)+d_H(2,2)) 
     div = div  + (dFdc(2,3)-dFdc(2,2))*(d_H(2,3)+d_H(2,2)) 
     div = div  + (dFdc(2,1)-dFdc(2,2))*(d_H(2,1)+d_H(2,2)) 
     
     div =0.5*div/(dx*dx)

     heat_size(3)=max(heat_size(3),abs(div*Y_energy)) !how big is this term
     Q_heat(3)=div*Y_energy
     TemperatureRHS = TemperatureRHS + div
  endif
endif
! grad double well term in phase equation

  if(grad_DW_term)then
     Grad_DW = (g_lambda*GradientEnergyx(vbles,1,dx)+potential(Phi,c,1)/g_lambda)*g_W(1)/MolPerVol*phi_dot(1) !J/mol
!     heat_size(4)=max(heat_size(4),abs(Grad_DW))
!     Q_heat(4)=-Grad_DW
    TemperatureRHS = TemperatureRHS - Grad_DW/(Cp*g_T0)   !no dimension  
  endif


end function TemperatureRHS

double precision function SoluteRHS(vbles,dx,c_c,a_trap)
! this uses volume element method. This guarantees conservation of solute. The main source of solute is the variation of free energy with phase: d (d F/dc) = ...+ (d^2 F/dc d phi) d phi + ..., which is revealed by a plot of c_dot.
use solution_parameters
implicit none
double precision, intent (in):: dx,c_c,vbles(N_unk+1,3,3),a_trap
integer ii,jj,ip,ih,jv
double precision c,T,Phi(N_phases),FreeEnergy,D
double precision, dimension(3,3)::stencil11=reshape((/0.,6.,0.,6.,-24.,6.,0.,6.,0./),(/3,3/))
double precision, dimension(3)::stencil1=(/-.5,0.,.5/)
double precision :: potential,D_c,f_c,MolPerVol
double precision :: c_dot=0d0,phi_dot(1)=0d0,beta,phi_,PhiSoluteRHS=0d0
double precision :: phidot, abs_grad_phi,u,FE

!!!


if(g_pure)then
write(*,*)"in SoluteRHS"
stop
endif

SoluteRHS=0.


!!
!x-coordinate transform values
!! Right
 ih=3
 jv=2
 include 'c_stencil.f90'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Left
 ih=1
 jv=2
 include 'c_stencil.f90'


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!up
 ih=2
 jv=3
 include 'c_stencil.f90'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !down
 ih=2
 jv=1
 include 'c_stencil.f90'
!!$
!!$ PhiSoluteRHS=0d0
!!$ if(cross_term)then
!!$    !!
!!$    !x-coordinate transform values
!!$    !! Right
!!$    ih=3
!!$    jv=2
!!$    
!!$    include 'phi_stencil.f90'
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    !!Left
!!$    ih=1
!!$    jv=2
!!$    
!!$    include 'phi_stencil.f90'
!!$    
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    !up
!!$    ih=2
!!$    jv=3
!!$    
!!$    include 'phi_stencil.f90'
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    !down
!!$    ih=2
!!$    jv=1
!!$    
!!$    include 'phi_stencil.f90'
!!$ endif
!!$ 
!!$ g_anti_trap = PhiSoluteRHS/(g_Dch*dx*dx*g_R*g_T0)!
 g_c_force = SoluteRHS/(g_Dch*dx*dx*g_R*g_T0)
! SoluteRHS = (SoluteRHS+PhiSoluteRHS)/(g_Dch*dx*dx*g_R*g_T0) 
 SoluteRHS = SoluteRHS/(g_Dch*dx*dx*g_R*g_T0) 

return

end function SoluteRHS
!!$
!!$double precision function PhiSoluteRHS(vbles,dx,c_c)
!!$! this uses volume element method. This guarantees conservation of solute. The main source of solute is the variation of free energy with phase: d (d F/dc) = ...+ (d^2 F/dc d phi) d phi + ..., which is revealed by a plot of c_dot.
!!$use solution_parameters
!!$implicit none
!!$double precision, intent (in):: dx,c_c,vbles(N_unk+1,3,3)
!!$integer ii,jj,ip,ih,jv
!!$double precision c,T,Phi(N_phases),FreeEnergy,D
!!$double precision, dimension(3,3)::stencil11=reshape((/0.,6.,0.,6.,-24.,6.,0.,6.,0./),(/3,3/))
!!$double precision, dimension(3)::stencil1=(/-.5,0.,.5/)
!!$double precision :: potential,D_c,f_c,SoluteRHS
!!$double precision :: c_dot=0d0,phi_dot(1)=0d0, a_trap
!!$
!!$
!!$
!!$SoluteRHS=0.
!!$
!!$
!!$!!
!!$!x-coordinate transform values
!!$!! Right
!!$ ih=3
!!$ jv=2
!!$ include 'c_stencil.f90'
!!$
!!$ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!Left
!!$ ih=1
!!$ jv=2
!!$ include 'c_stencil.f90'
!!$
!!$
!!$ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!up
!!$ ih=2
!!$ jv=3
!!$ include 'c_stencil.f90'
!!$
!!$ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$ !down
!!$ ih=2
!!$ jv=1
!!$ include 'c_stencil.f90'
!!$
!!$ PhiSoluteRHS = SoluteRHS/(g_Dch*dx*dx)
!!$
!!$return
!!$
!!$end function PhiSoluteRHS
double precision function dw(x,d)
double precision, intent(in)::x
integer, intent(in):: d
if(d.eq.0)then
dw = x**2*(1d0 - x)**2
elseif(d.eq.1)then
dw = 2d0*x*(1d0-x)*(1d0-2d0*x)
endif
end function dw

double precision function potential(phi,c,lp)
use solution_parameters
implicit none

double precision, intent (in)::phi(N_phases),c
integer, intent(in):: lp


potential=0d0
if(lp.eq.1)then
   potential =  2d0*phi(1)*(1d0-phi(1))*(1d0-2d0*phi(1))
elseif(lp.eq.3)then
   potential = 0d0
else
   stop
endif

return


end function potential

double precision function potentialx(phi,c,lp)
use solution_parameters
implicit none

double precision, intent (in)::phi(N_phases),c
integer, intent(in):: lp
integer :: ip,jp
double precision :: u,Omega,Omega1,n ,dw,h_phi,u_phi,g,dg,p,dp

potentialx=0d0
if(g_equal_delta.and. lp.eq.1)return
if(g_skew_potential)then
   if(phi(1).lt.1d-6)return
   u=0.5*(log(2d0-phi(1))-log(phi(1)))
   if(lp.eq.1)then
      Omega1 = u*phi(1)*(phi(1)-2d0)*(2d0*u*(phi(1)-1d0)+1d0)
      if(g_pure)then
         potentialx = ((1.-c)+c*g_W(2)/g_W(1))* Omega1
      else
         potentialx = Omega1
      endif
   elseif(lp.eq.3)then
      Omega = 0.5*(phi(1)*(2d0-phi(1))*u)**2
      potentialx = (g_W(2)/g_W(1)-1d0)* Omega
   else
      stop
   endif
elseif(g_power_skew)then
   n = g_power
   if(lp.eq.1)then
      Omega1 = ((n*phi(1)+2d0*phi(1)-2d0)*(n+2d0)**2*phi(1)*(-(-1+phi(1))*(n+2d0)/n)**n/(-1d0+phi(1))/64d0)
      if(g_pure)then
         potentialx = Omega1
      else
         potentialx = ((1.-c)+c*g_W(2)/g_W(1))*Omega1
      endif
   elseif(lp.eq.3)then
      Omega = ((2d0+n)*phi(1)/8d0)**2*((1d0-phi(1))*(1+2d0/n))**n 
      potentialx = (g_W(2)/g_W(1)-1d0)* Omega
   else
      stop
   endif
elseif(g_s_skew)then
   if(g_pure)then
      if(lp.eq.0)then
         !AMM
         if(g_AMM.eq.1)then
            potentialx = 0.5D2*phi(1)**2*(-1.D0+phi(1))**2/(0.2035214409D0+0.2041189996D1*phi(1))**2         
            
            potentialx = potentialx/0.4673293612476353D2
         elseif(g_AMM.eq.2)then
            potentialx = phi(1)**2*(-1.D0+phi(1))**2
         else
            write(*,*)"g_AKK 1 or 2 only"
         endif
      elseif(lp.eq.1)then
         !AMM
         if(g_AMM.eq.1)then
            potentialx = 15625.D0*phi(1)*(-1.D0+phi(1))*(-0.2035214409D2+0.4070428819D2*phi(1)+0.2041189996D3*phi(1)**2)/(0.5088036022D1+0.510297499D2*phi(1))**3
            potentialx = potentialx/0.4673293612476353D2
         elseif(g_AMM.eq.2)then
            potentialx = 0.9346587224952705D2*phi(1)-0.2803976167485812D3*phi(1)**2+0.1869317444990541D3*phi(1)**3
            potentialx = potentialx/0.4673293612476353D2
         else
            write(*,*)"g_AKK 1 or 2 only"
         endif
      endif
   else !solute with a skew
      if(lp.eq.0)then
         !AMM
         if(g_AMM.eq.1)then
            potentialx = 0.5D2*phi(1)**2*(-1.D0+phi(1))**2/(0.2035214409D0+0.2041189996D1*phi(1))**2         
            
            potentialx = potentialx/0.4673293612476353D2
         elseif(g_AMM.eq.2)then
            potentialx = phi(1)**2*(-1.D0+phi(1))**2
         else
            write(*,*)"g_AMM 1 or 2 only"
         endif
         
         potentialx = potentialx * ((1.-c)+c*g_W(2)/g_W(1))
      elseif(lp.eq.1)then
         !AMM
         if(g_AMM.eq.1)then
            potentialx = 15625.D0*phi(1)*(-1.D0+phi(1))*(-0.2035214409D2+0.4070428819D2*phi(1)+0.2041189996D3*phi(1)**2)/(0.5088036022D1+0.510297499D2*phi(1))**3
            potentialx = potentialx/0.4673293612476353D2
         elseif(g_AMM.eq.2)then
            potentialx = 0.9346587224952705D2*phi(1)-0.2803976167485812D3*phi(1)**2+0.1869317444990541D3*phi(1)**3
            potentialx = potentialx/0.4673293612476353D2
         else
            write(*,*)"g_AMM 1 or 2 only"
         endif
         potentialx = potentialx * ((1.-c)+c*g_W(2)/g_W(1))
      elseif(lp.eq.3)then
         !AMM
         if(g_AMM.eq.1)then
            potentialx = 0.5D2*phi(1)**2*(-1.D0+phi(1))**2/(0.2035214409D0+0.2041189996D1*phi(1))**2         
            
            potentialx = potentialx/0.4673293612476353D2
         elseif(g_AMM.eq.2)then
            potentialx = phi(1)**2*(-1.D0+phi(1))**2
         else
            write(*,*)"g_AMM 1 or 2 only"
         endif
         potentialx = potentialx * (g_W(2)/g_W(1)-1d0)
      endif
   endif
else
   if(lp.eq.1)then
      if(g_pure)then
         potentialx =  2d0*phi(1)*(1d0-phi(1))*(1d0-2d0*phi(1))
      else
         potentialx =  ((1.-c)+c*g_W(2)/g_W(1))* 2d0*phi(1)*(1d0-phi(1))*(1d0-2d0*phi(1))
      endif
   elseif(lp.eq.3)then
      potentialx = (g_W(2)/g_W(1)-1d0)* phi(1)**2*(1d0-phi(1))**2
   else
      stop
   endif
endif
return


end function potentialx



double precision function PhaseRHS(vbles,dx,lp,Phi,c,T)
use solution_parameters
implicit none
double precision, intent(in)::dx,c,T,Phi(N_phases),vbles(N_unk+1,3,3)
integer, intent(in):: lp
double precision :: LapPhi(N_phases)
double precision GradientEnergy,GradientEnergyx,GradientEnergy_2theta,gradientEnergy_ellips,potential,FreeEnergy,PhiSoluteRHS,GradientEnergyHexPrism,SurfaceEnergy
double precision GradientEnergy_circ,GradientEnergyHexPrismVol,GradientEnergyHexPrismVol2,GradientEnergyHexPrismVol3
double precision ::M_tilde 
double precision :: c_dot=0,phi_dot=0
double precision :: MolPerVol,beta
double precision :: phix,phiy, xx(2)
double precision :: u,v,g,A,X,Y,a_factor=1d1
logical :: new_formulation=.true.
integer :: n=6


if(g_pure)then
   MolPerVol=g_MolPerVol(1)

!!$   phix = 0.5*(vbles(1,3,2) - vbles(1,1,2))/dx 
!!$   phiy = 0.5*(vbles(1,2,3) - vbles(1,2,1))/dx 
!!$   M_tilde = 1d0+epsilon* sqrt((phix**3-3*phix*phiy**2)**2/(phix**2+phiy**2+0.0001)**3)
!!$   M_tilde = M_tilde/(1d0+epsilon)
   M_tilde = 1d0
elseif(AlNi.eq.8)then
   M_tilde = 1d0
else
   MolPerVol=(1.-c)*g_MolPerVol(1)+c*g_MolPerVol(2)
!   M_tilde = (1.-c)+c*g_M(2)/g_M(1)
   M_tilde = g_mu(1)*((1.-c)+c*g_M(2)/g_M(1))
endif

 !assemble and non-dimensionlise phaseRHS 


if(g_kinetic)then
   if(hex_facet)then
      if(g_needle)then
         phix = 0.5*(vbles(1,3,2) - vbles(1,1,2))/dx 
         phiy = 0.5*(vbles(1,2,3) - vbles(1,2,1))/dx 
! This is for needle growth. A_factor can be any positive number
         M_tilde = M_tilde / (1+a_factor*epsilon)*(1 + a_factor*epsilon*abs(phix)/sqrt(phix**2+phiy**2+0.0001))

         


      else
         phix = 0.5*(vbles(1,3,2) - vbles(1,1,2))/dx
         phiy = 0.5*(vbles(1,2,3) - vbles(1,2,1))/dx
         u = phix*phix
         v = phiy*phiy
         g  = u + v
         if(g.gt.1d-20)then

            X = u/g
            Y = v/g
            A = (A_0*(1d0-epsilon_tilde/15*(X**2+Y**2)))**2
         else 
            A=1d0
         endif
         M_tilde=M_tilde*A

      endif
   else
      phix = 0.5*(vbles(1,3,2) - vbles(1,1,2))/dx
      phiy = 0.5*(vbles(1,2,3) - vbles(1,2,1))/dx
      u = phix*phix
      v = phiy*phiy
      g  = u + v
      if(g.gt.1d-20)then
         
         X = u/g
         Y = v/g
         A = (A_0*(1d0+epsilon_tilde*(X**2+Y**2)))**2
      else
         A=1d0
      endif
      M_tilde=M_tilde*A
   endif
endif

if(hex_facet)then
   if(AlNi.eq.8)then
!SurfaceEnergy(p,x,n)

!      phaseRHS=- M_tilde*(GradientEnergy_circ(vbles,lp,dx)+potential(Phi,c,lp)/g_lambda**2&
!           +FreeEnergy(c,T,Phi,0,0,1,.true.)/(g_L(1)*g_lambda**2))
!!$      phix = 0.5*(vbles(1,3,2) - vbles(1,1,2))/dx
!!$      phiy = 0.5*(vbles(1,2,3) - vbles(1,2,1))/dx
!!$      g = sqrt(phix**2 + phiy**2)
!!$      if(g.gt.1d-9)then
!!$         xx(1) = phix/g
!!$         xx(2) = phiy/g
!!$         u =   SurfaceEnergy(g_phex,xx,n)
!!$         M_tilde= M_tilde/u
!!$      endif
      phaseRHS=- M_tilde*(GradientEnergyHexPrismVol2(vbles,lp,dx)+potential(Phi,c,lp)/g_lambda**2&
           +FreeEnergy(c,T,Phi,0,0,1,.true.)/(g_L(1)*g_lambda**2))


   elseif(new_formulation)then
      phaseRHS=- M_tilde*(GradientEnergy_circ(vbles,lp,dx)+potential(Phi,c,lp)/g_lambda**2&
           +FreeEnergy(c,T,Phi,c_dot,phi_dot,1,g_approx)*MolPerVol/(g_W(1)*g_lambda**2))
   else
      phaseRHS=- M_tilde*(GradientEnergy(vbles,lp,dx)+potential(Phi,c,lp)/g_lambda**2&
           +FreeEnergy(c,T,Phi,c_dot,phi_dot,1,g_approx)*MolPerVol/(g_W(1)*g_lambda**2))
   endif
elseif(g_ellipse)then
   phaseRHS=- M_tilde*(GradientEnergy_2theta(vbles,lp,dx)+potential(Phi,c,lp)/g_lambda**2&
        +FreeEnergy(c,T,Phi,c_dot,phi_dot,1,g_approx)*MolPerVol/(g_W(1)*g_lambda**2))
else !standard: 1+epsilon*cos(4 theta)
   phaseRHS=- M_tilde*(GradientEnergyx(vbles,lp,dx)+potential(Phi,c,lp)/g_lambda**2&
        +FreeEnergy(c,T,Phi,c_dot,phi_dot,1,g_approx)*MolPerVol/(g_W(1)*g_lambda**2))
endif
! if(cross_term)then
!    beta = g_beta/g_lambda !
!    phaseRHS = phaseRHS + g_sym* beta*PhiSoluteRHS(vbles,dx,c)/(g_R*g_T0) 
! endif
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
double precision :: c,T,phi(N_phases),c_dot,phi_dot(N_phases)
double precision phi_rhs(N_phases),rfactor,rf4,rf5,rf6,S
double precision TemperatureRHS,SoluteRHS,T_k,phi_
double precision vbles(N_unk+1,3,3),unk_(M_unk),MolPerVol,vbles10(N_unk+1,10,10)
logical Not_number
integer ip,n,ii,jj,kk


  do ii=1,10
     do jj=1,10
        do ip=1,N_unk
           vbles10(ip,ii,jj)=unk(1+(ip-1)*nunkvbles,ii,jj,1,lb)
        enddo
     enddo
  enddo


  Fi=0.



  do ii=1,3
     do jj=1,3
        do ip=1,N_unk
           vbles(ip,ii,jj)=unk(1+(ip-1)*nunkvbles,i+ii-2,j+jj-2,1,lb)
        enddo
     enddo
  enddo



  ip=N_unk+1
  do ii=1,3
     do jj=1,3
        if(mode_1.eq.1) then
           vbles(ip,ii,jj)=(unk(1,i+ii-2,j+jj-2,1,lb)-&
                unk(2,i+ii-2,j+jj-2,1,lb))/dt
        else
           rfactor = dt/dtold
           rf4 = rfactor*rfactor/(rfactor+1.0)
           rf5 = (rfactor+1.0)
           rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
           vbles(ip,ii,jj)=(rf4*unk(5,i+ii-2,j+jj-2,1,lb)-&
                rf5*unk(4,i+ii-2,j+jj-2,1,lb)+&
                rf6*unk(1,i+ii-2,j+jj-2,1,lb))/dt
        endif

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


  call get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_,vbles10,i,j)

!  unk(2+3*nunkvbles,i,j,1,lb)=g_anti_trap
!  unk(3+3*nunkvbles,i,j,1,lb)=g_c_force
!
if(.not.g_pure)then
  unk(1+N_unk*6,i,j,1,lb)=Q_heat(1)
  unk(2+N_unk*6,i,j,1,lb)=Q_heat(2)
  unk(3+N_unk*6,i,j,1,lb)=Q_heat(3)
  unk(4+N_unk*6,i,j,1,lb)=Q_heat(4)
endif

 
end subroutine get_RHS_Fi_2d

subroutine get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_,vbles10,ix,iy)

  use solution_parameters

  implicit none

!
integer, intent(in)::mode_1,nunkvbles,lnblocks
double precision, intent(in)::dx,dt,dtold,unk_(M_unk),ix,iy
double precision, intent(out):: Fi(N_unk)
double precision, intent(in):: vbles(N_unk+1,3,3),vbles10(N_unk+1,10,10)
double precision :: c,T,phi(N_phases),c_dot,phi_dot(N_phases)
double precision rfactor,rf4,rf5,rf6,S,GradientEnergy,potential
double precision TemperatureRHS,SoluteRHS,PhaseRHS,T_k,phi_,FreeEnergy,MolPerVol,a_trap,grad_phi,grad_c,sum
logical Not_number
integer ip,n,i,j,k


Fi=0d0


  if (mode_1.eq.1) then
!     phi_dot(1) = (vbles(1,2,2)-unk_(2))/dt
     phi_dot(1)=vbles(N_unk+1,2,2)
     if(.not.g_pure)then
        n=2*nunkvbles
        c_dot = (vbles(3,2,2)-unk_(2+n))/dt
     else 
        c_dot=0d0
     endif
  elseif(mode_1.eq.2)then
     rfactor = dt/dtold

     rf4 = rfactor*rfactor/(rfactor+1.0)
     rf5 = (rfactor+1.0)
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
     
!!     phi_dot(1) = (rf4*unk_(5)-rf5*unk_(4)+rf6*vbles(1,2,2))/dt
     phi_dot(1) = vbles(N_unk+1,2,2)

     if(.not.g_pure)then
        n=2*nunkvbles
        c_dot = (rf4*unk_(5+n)-rf5*unk_(4+n)+rf6*vbles(3,2,2))/dt
     else
        c_dot=0d0
     endif
  else
     write(*,*)"mode_1 not = 1 or 2"
     stop
  endif


  phi(1)=vbles(1,2,2)
  T = vbles(2,2,2)!
  if(.not.g_pure)then
     c = vbles(3,2,2)
  else
     c=0d0
     c_dot=0d0
  endif
  


!

  Fi(1) = PhaseRHS(vbles,dx,1,Phi,c,T)

  if(thermal)then
     Fi(2) = TemperatureRHS(vbles,dx,c,T,phi,c_dot,phi_dot)
   else
     Fi(2)=0d0
  endif

  a_trap=0d0
  if(g_alpha.ne.0d0)then
     a_trap =g_alpha*sqrt(8d0)*g_lambda*g_ratio*(g_Dch/g_D(1))



!!$   equals g_alpha * d * (g_Dch/g_D(1))
  endif
!!$     if(abs(vbles(1,2,2)-0.5).lt.0.1)then
!!$        write(*,*)dx/(0.5*sqrt((vbles(1,3,2)-vbles(1,1,2))**2+(vbles(1,2,3)-vbles(1,2,1))**2))
!!$        write(*,*)sqrt(8d0)*g_lambda*g_ratio,dx
!!$     endif



!
  if(.not.g_pure)  Fi(3) = SoluteRHS(vbles,dx,c,a_trap)

  if(grad_DW_term.and. thermal) then
     write(*,*)"grad_dw term needs looking at"
     stop
     MolPerVol=(1.-c)*g_MolPerVol(1)+c*g_MolPerVol(2)
     Fi(4) = (g_lambda**2*GradientEnergy(vbles,1,dx)+potential(phi,c,1))/MolPerVol*g_W(1) +FreeEnergy(c,T,Phi,c_dot,phi_dot,1,g_approx)
  else
     if(g_pure)then
        Fi(3)=0d0
     else
        Fi(4)=0d0
     endif
  endif

  do ip=1,N_unk
     n=(ip-1)*nunkvbles
     if(mode_1.le.1)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)-unk_(n+2)
     else if(mode_1.eq.2)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)/rf6-unk_(n+2)
     end if
  enddo

!  Fi(4)=vbles(4,2,2)-Fi(4)
!


end subroutine get_RHS_Fi

subroutine Calc_GradPhi(LapPhi,vbles,dx)
use solution_parameters
implicit none
double precision, intent(in)::dx,vbles(N_unk+1,3,3)
double precision, intent(out)::LapPhi(N_phases)
double precision, dimension(3,3)::stencil=reshape((/1.,4.,1.,4.,-20.,4.,1.,4.,1./),(/3,3/))
double precision, dimension(3,3):: stencil5=reshape((/0.,1.,0.,1.,-4.,1.,0.,1.,0./),(/3,3/))
integer ip,ii,jj




 LapPhi=0.
 do ip=1,N_phases
    do ii=1,3
       do jj=1,3
          LapPhi(ip)=LapPhi(ip)+stencil(ii,jj)*vbles(ip,ii,jj)/(6.*dx*dx)
       end do
    end do
 end do
 return
 
end subroutine Calc_GradPhi


double precision function FreeEnergy(c,T,phi,c_dot,phi_dot,lp,approx)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi(N_phases),c_dot,phi_dot(N_phases)
integer, intent(in)::lp
logical, intent(in)::approx
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(2),hb(2),hr(2)
double precision :: Tk,Cp,R,dc,df,dfdc
double precision :: FE,dh=1d-6
double precision Gibbs_FE_liq,Gibbs_FE_sol,Gibbs_FE_l,Gibbs_FE_e,T_K,g_phi,EntropyMixing,Gibbs_SiGeliq,Gibbs_SiGesol !functions
double precision :: Gibbs_FE_l_AlSi,Gibbs_FE_e_AlSi,cc,Gibbs_Fe_l_PbSn,Gibbs_Fe_e_PbSn,Gibbs_Fe_l_NiCu,Gibbs_Fe_e_NiCu
double precision :: Gibbs_Fe_La,Gibbs_Fe_a,AlFeFreeEnergy
integer :: i,j

if(g_pure.and. lp.gt. 2)then
   write(*,*)"lp> 2:",lp
   stop
endif


if(AlNi.eq.1)then
   if(lp.eq.0)then !free energy proper J/mol
      
      FreeEnergy =              Gibbs_FE_liq(c,T,0,0)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_FE_sol(c,T,0,0)*g_phi(1d0-phi(1),0)
      
      return
   elseif(lp.eq.1)then !phase
      
      FreeEnergy =              Gibbs_FE_liq(c,T,0,0)*g_phi(phi(1),1)
      FreeEnergy = FreeEnergy - Gibbs_FE_sol(c,T,0,0)*g_phi(phi(1),1)
      
      return
   elseif(lp.eq.3)then !Solute
      
      FreeEnergy = FreeEnergy + Gibbs_FE_liq(c,T,1,0)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_FE_sol(c,T,1,0)*(1.-g_phi(phi(1),0))
      
      return
   elseif(lp.eq.2)then !Temperature Keep in J/mol
      Tk=T_k(T)
      
      !-(1-T*d/dT)dF/dphi
!      if(heat_phi_term)then
         !  T (d/dT)dF/dphi
         FreeEnergy =  Tk*(Gibbs_FE_liq(c,T,0,1) - Gibbs_FE_sol(c,T,0,1))
         ! -Df/Dphi
         FreeEnergy = FreeEnergy - Gibbs_FE_liq(c,T,0,0) + Gibbs_FE_sol(c,T,0,0)
         g_latent_heat = FreeEnergy         
         FreeEnergy = FreeEnergy*phi_dot(1)*g_phi(phi(1),1)
         
         Q_heat(1) = FreeEnergy
!      endif !end phi_dot_term
      if(heat_c_term)then 
         !   -(1- T(d/dT))dF/dc 
         FE=(Tk*Gibbs_FE_liq(c,T,1,1)*g_phi(phi(1),0)- Gibbs_FE_liq(c,T,1,0))*g_phi(phi(1),0)
         FE=FE+(Tk*Gibbs_FE_sol(c,T,1,1)- Gibbs_FE_sol(c,T,1,0))*(1.-g_phi(phi(1),0))
         FE=FE*c_dot
         FreeEnergy = FreeEnergy + FE
         Q_heat(2)= FE
         !      FreeEnergy = FreeEnergy + Tk*Gibbs_FE_liq(c,T,1,1)*g_phi(phi(1),0)*c_dot
         !      FreeEnergy = FreeEnergy + Tk*Gibbs_FE_sol(c,T,1,1)*(1.-g_phi(phi(1),0))*c_dot
         !      FreeEnergy = FreeEnergy - Gibbs_FE_liq(c,T,1,0)*g_phi(phi(1),0)*c_dot
         !      FreeEnergy = FreeEnergy - Gibbs_FE_sol(c,T,1,0)*(1.-g_phi(phi(1),0))*c_dot
      endif !c_dot_term
      !   heat capacity = -T (d/dT) dF/dT
      g_C =  -Tk*(g_phi(phi(1),0)*Gibbs_FE_liq(c,T,0,2)+g_phi(1d0-phi(1),0)*Gibbs_FE_sol(c,T,0,2))
      return
   elseif(lp.eq.-1)then !for anti-trapping
      
      FreeEnergy = FreeEnergy + Gibbs_FE_l(c,T,1,0)*g_phi(phi(1),1)
      FreeEnergy = FreeEnergy - Gibbs_FE_e(c,T,1,0)*g_phi(phi(1),1)
      
      return
   elseif(lp.eq.-2)then !for anti-trapping f_cc
      
      !resorting to approximation for gibbs_fe_e as it looks too complicated
      FreeEnergy =              0.5*(Gibbs_FE_liq(c+dh,T,1,0)-Gibbs_FE_liq(c-dh,T,1,0))/dh*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + 0.5*(Gibbs_FE_sol(c+dh,T,1,0)-Gibbs_FE_sol(c-dh,T,1,0))/dh*(1d0-g_phi(phi(1),0))
      
      return
   else
      write(*,*)"bad current_var"
      stop
   endif
elseif(AlNi.eq.2)then

!!$   write(*,*)c,T,phi(1)
!!$   write(*,*)lp
!!$   write(*,*)Gibbs_FE_l(c,T,0,0),Gibbs_FE_e(c,T,0,0)
!!$   write(*,*)Gibbs_FE_l(c,0d0,0,0),Gibbs_FE_e(c,0d0,0,0)
!!$   stop

   if(lp.eq.0)then !free energy proper J/mol
      
      FreeEnergy =              Gibbs_FE_l(c,T,0,0)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_FE_e(c,T,0,0)*g_phi(1d0-phi(1),0)
      
      return
   elseif(lp.eq.1)then !phase

      
      FreeEnergy =              Gibbs_FE_l(c,T,0,0)*g_phi(phi(1),1)
      FreeEnergy = FreeEnergy - Gibbs_FE_e(c,T,0,0)*g_phi(phi(1),1)
      
      return
   elseif(lp.eq.3)then !Solute
      
      FreeEnergy = FreeEnergy + Gibbs_FE_l(c,T,1,0)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_FE_e(c,T,1,0)*(1.-g_phi(phi(1),0))
      
      return
   elseif(lp.eq.-1)then !for anti-trapping f_c_phi
      FreeEnergy =              (Gibbs_FE_l(c,T,1,0) - Gibbs_FE_e(c,T,1,0))*g_phi(phi(1),1)
      return
   elseif(lp.eq.-2)then !for anti-trapping f_cc

      !resorting to approximation for gibbs_fe_e as it looks too complicated
      FreeEnergy =              0.5*(Gibbs_FE_l(c+dh,T,1,0)-Gibbs_FE_l(c-dh,T,1,0))/dh*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + 0.5*(Gibbs_FE_e(c+dh,T,1,0)-Gibbs_FE_e(c-dh,T,1,0))/dh*(1d0-g_phi(phi(1),0))

      return
   elseif(lp.eq.2)then !Temperature Keep in J/mol
      Tk=T_k(T)
      
      !-(1-T*d/dT)dF/dphi
!      if(heat_phi_term)then
         !  T (d/dT)dF/dphi
         FreeEnergy =  Tk*(Gibbs_FE_l(c,T,0,1) - Gibbs_FE_e(c,T,0,1))
         ! -Df/Dphi
         FreeEnergy = FreeEnergy - Gibbs_FE_l(c,T,0,0) + Gibbs_FE_e(c,T,0,0)
         g_latent_heat = FreeEnergy         
         FreeEnergy = FreeEnergy*phi_dot(1)*g_phi(phi(1),1)
         
         Q_heat(1) = FreeEnergy
!      endif !end phi_dot_term
!      if(heat_c_term)then 
         !   -(1- T(d/dT))dF/dc 
         FE=(Tk*Gibbs_FE_l(c,T,1,1)*g_phi(phi(1),0)- Gibbs_FE_l(c,T,1,0))*g_phi(phi(1),0)
         FE=FE+(Tk*Gibbs_FE_e(c,T,1,1)- Gibbs_FE_e(c,T,1,0))*(1.-g_phi(phi(1),0))
         FE=FE*c_dot
         FreeEnergy = FreeEnergy + FE
         Q_heat(2)= FE
         !      FreeEnergy = FreeEnergy + Tk*Gibbs_FE_l(c,T,1,1)*g_phi(phi(1),0)*c_dot
         !      FreeEnergy = FreeEnergy + Tk*Gibbs_FE_e(c,T,1,1)*(1.-g_phi(phi(1),0))*c_dot
         !      FreeEnergy = FreeEnergy - Gibbs_FE_l(c,T,1,0)*g_phi(phi(1),0)*c_dot
         !      FreeEnergy = FreeEnergy - Gibbs_FE_e(c,T,1,0)*(1.-g_phi(phi(1),0))*c_dot
!      endif !c_dot_term
      !   heat capacity = -T (d/dT) dF/dT
      g_C =  -Tk*(g_phi(phi(1),0)*Gibbs_FE_l(c,T,0,2)+g_phi(1d0-phi(1),0)*Gibbs_FE_e(c,T,0,2))
      return
   else
      write(*,*)"bad current_var"
      stop
   endif

elseif(AlNi.eq.3)then
   if(lp.eq.0)then !free energy proper J/mol
      
      FreeEnergy =              Gibbs_SiGeliq(c,T,0,0)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_SiGesol(c,T,0,0)*g_phi(1d0-phi(1),0)
      
      return
   elseif(lp.eq.1)then !phase
      
      FreeEnergy =              Gibbs_SiGeliq(c,T,0,0)*g_phi(phi(1),1)
      FreeEnergy = FreeEnergy - Gibbs_SiGesol(c,T,0,0)*g_phi(phi(1),1)
      
      return
   elseif(lp.eq.3)then !Solute
      
      FreeEnergy = FreeEnergy + Gibbs_SiGeliq(c,T,1,0)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_SiGesol(c,T,1,0)*(1.-g_phi(phi(1),0))
      
      return
   elseif(lp.eq.-1)then !for anti-trapping
      
      FreeEnergy = FreeEnergy + Gibbs_SiGeliq(c,T,1,0)*g_phi(phi(1),1)
      FreeEnergy = FreeEnergy + Gibbs_SiGesol(c,T,1,0)*(1.-g_phi(phi(1),1))
      
      return
   elseif(lp.eq.2)then !Temperature Keep in J/mol
      Tk=T_k(T)
      
      !-(1-T*d/dT)dF/dphi
 !     if(heat_phi_term)then
         !  T (d/dT)dF/dphi
         FreeEnergy =  Tk*(Gibbs_SiGeliq(c,T,0,1) - Gibbs_SiGesol(c,T,0,1))
         ! -Df/Dphi
         FreeEnergy = FreeEnergy - Gibbs_SiGeliq(c,T,0,0) + Gibbs_SiGesol(c,T,0,0)
         g_latent_heat = FreeEnergy         
         FreeEnergy = FreeEnergy*phi_dot(1)*g_phi(phi(1),1)
         
         Q_heat(1) = FreeEnergy
 !     endif !end phi_dot_term
 !     if(heat_c_term)then 
         !   -(1- T(d/dT))dF/dc 
         FE=(Tk*Gibbs_SiGeliq(c,T,1,1)*g_phi(phi(1),0)- Gibbs_SiGeliq(c,T,1,0))*g_phi(phi(1),0)
         FE=FE+(Tk*Gibbs_SiGesol(c,T,1,1)- Gibbs_SiGesol(c,T,1,0))*(1.-g_phi(phi(1),0))
         FE=FE*c_dot
         FreeEnergy = FreeEnergy + FE
         Q_heat(2)= FE
         !      FreeEnergy = FreeEnergy + Tk*Gibbs_SiGeliq(c,T,1,1)*g_phi(phi(1),0)*c_dot
         !      FreeEnergy = FreeEnergy + Tk*Gibbs_SiGesol(c,T,1,1)*(1.-g_phi(phi(1),0))*c_dot
         !      FreeEnergy = FreeEnergy - Gibbs_SiGeliq(c,T,1,0)*g_phi(phi(1),0)*c_dot
         !      FreeEnergy = FreeEnergy - Gibbs_SiGesol(c,T,1,0)*(1.-g_phi(phi(1),0))*c_dot
 !     endif !c_dot_term
      !   heat capacity = -T (d/dT) dF/dT
      g_C =  -Tk*(g_phi(phi(1),0)*Gibbs_SiGeliq(c,T,0,2)+g_phi(1d0-phi(1),0)*Gibbs_SiGesol(c,T,0,2))
      return
   else
      write(*,*)"bad current_var"
      stop
   endif
elseif(AlNi.eq.4)then !AlSi
   if(lp.eq.0)then !free energy proper J/mol  
      FreeEnergy =              Gibbs_Fe_l_AlSi(c,T,0,0,approx)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_FE_e_AlSi(c,T,0,0,approx)*g_phi(1d0-phi(1),0)
      
      return
   elseif(lp.eq.1)then !phase
      
      FreeEnergy =              Gibbs_Fe_l_AlSi(c,T,0,0,approx)*g_phi(phi(1),1)
      FreeEnergy = FreeEnergy - Gibbs_FE_e_AlSi(c,T,0,0,approx)*g_phi(phi(1),1)
      
      return
   elseif(lp.eq.3)then !Solute
      
      FreeEnergy =  Gibbs_Fe_l_AlSi(c,T,1,0,approx)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_FE_e_AlSi(c,T,1,0,approx)*(1.-g_phi(phi(1),0))
      
      return
   elseif(lp.eq.2)then !Temperature Keep in J/mol

      Tk=T_k(T)
      
      !-(1-T*d/dT)dF/dphi
      !     if(heat_phi_term)then
      !  T (d/dT)dF/dphi
      
      FreeEnergy =  Tk*(Gibbs_Fe_l_AlSi(c,T,0,1,approx) - Gibbs_FE_e_AlSi(c,T,0,1,approx))

         ! -Df/Dphi
         FreeEnergy = FreeEnergy - Gibbs_Fe_l_AlSi(c,T,0,0,approx) + Gibbs_FE_e_AlSi(c,T,0,0,approx)
         g_latent_heat = FreeEnergy         
         FreeEnergy = FreeEnergy*phi_dot(1)*g_phi(phi(1),1)
         
         Q_heat(1) = FreeEnergy
 !     endif !end phi_dot_term
      if(heat_c_term)then 
         !   -(1- T(d/dT))dF/dc 
         FE=(Tk*Gibbs_Fe_l_AlSi(c,T,1,1,approx)*g_phi(phi(1),0)- Gibbs_Fe_l_AlSi(c,T,1,0,approx))*g_phi(phi(1),0)
         FE=FE+(Tk*Gibbs_FE_e_AlSi(c,T,1,1,approx)- Gibbs_FE_e_AlSi(c,T,1,0,approx))*(1.-g_phi(phi(1),0))
         FE=FE*c_dot
         FreeEnergy = FreeEnergy + FE
         Q_heat(2)= FE

      endif !c_dot_term
      !   heat capacity = -T (d/dT) dF/dT
      g_C =  -Tk*(g_phi(phi(1),0)*Gibbs_Fe_l_AlSi(c,T,0,2,approx)+g_phi(1d0-phi(1),0)*Gibbs_FE_e_AlSi(c,T,0,2,approx))

      return
 !     else
 !        write(*,*)approx
 !        stop
 !     endif
   else
      write(*,*)"bad current_var"
      stop
   endif
elseif(AlNi.eq.5)then !PbSn
   if(lp.eq.0)then !free energy proper J/mol  
      FreeEnergy =              Gibbs_Fe_l_PbSn(c,T,0,0,approx)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_FE_e_PbSn(c,T,0,0,approx)*g_phi(1d0-phi(1),0)

      return
   elseif(lp.eq.1)then !phase
      if(g_plapp)then
         dc = g_c_max-g_c_min
         df = g_f_c_max - g_f_c_min
         dfdc =       Gibbs_Fe_l_PbSn(c+(1d0-g_phi(phi(1),0))*dc,T,1,0,approx)*g_phi(phi(1),0)
         dfdc = dfdc +Gibbs_Fe_e_PbSn(c - g_phi(phi(1),0)*dc,T,1,0,approx)*(1.-g_phi(phi(1),0))


!!$         if(abs(phi(1)-0.5).lt.0.1)then
!!$            FreeEnergy =              Gibbs_Fe_l_PbSn(c,T,0,0,approx)
!!$            FreeEnergy = FreeEnergy - Gibbs_FE_e_PbSn(c,T,0,0,approx)
!!$            FreeEnergy = FreeEnergy * g_phi(phi(1),1)
!!$            write(*,*)"1",FreeEnergy,phi(1),c
!!$            FreeEnergy =               Gibbs_Fe_l_PbSn(c+(1d0-g_phi(phi(1),0))*dc,T,0,0,approx)
!!$            FreeEnergy = FreeEnergy -  Gibbs_Fe_e_PbSn(c - g_phi(phi(1),0)*dc,T,0,0,approx)  - dc * dfdc
!!$            FreeEnergy = FreeEnergy * g_phi(phi(1),1)
!!$            write(*,*)"2",FreeEnergy,phi(1),c
!!$
!!$         else
            FreeEnergy =               Gibbs_Fe_l_PbSn(c+(1d0-g_phi(phi(1),0))*dc,T,0,0,approx)
            FreeEnergy = FreeEnergy -  Gibbs_Fe_e_PbSn(c - g_phi(phi(1),0)*dc,T,0,0,approx)  -dc * dfdc
            FreeEnergy = FreeEnergy * g_phi(phi(1),1)
!         endif


      else
         FreeEnergy =              Gibbs_Fe_l_PbSn(c,T,0,0,approx)
         FreeEnergy = FreeEnergy - Gibbs_FE_e_PbSn(c,T,0,0,approx)
         FreeEnergy = FreeEnergy * g_phi(phi(1),1)
      endif
      return
   elseif(lp.eq.3)then !Solute
      if(g_plapp)then      
         dc = g_c_max-g_c_min
         df = g_f_c_max - g_f_c_min


         if(.false.)then
            FreeEnergy =              Gibbs_Fe_l_PbSn(c,T,1,0,approx)*g_phi(phi(1),0)
            FreeEnergy = FreeEnergy + Gibbs_FE_e_PbSn(c,T,1,0,approx)*(1.-g_phi(phi(1),0))
            write(*,*)"1",FreeEnergy,phi(1),c
            
            
            FreeEnergy =              Gibbs_Fe_l_PbSn(c+(1d0-g_phi(phi(1),0))*dc,T,1,0,approx)*g_phi(phi(1),0)
            FreeEnergy = FreeEnergy + Gibbs_Fe_e_PbSn(c - g_phi(phi(1),0)*dc,T,1,0,approx)*(1.-g_phi(phi(1),0))
            write(*,'(7F14.5)')FreeEnergy,phi(1),c+(1d0-g_phi(phi(1),0))*dc,c - g_phi(phi(1),0)*dc,dc,g_c_min,g_c_max

         else
            FreeEnergy =              Gibbs_Fe_l_PbSn(c+(1d0-g_phi(phi(1),0))*dc,T,1,0,approx)*g_phi(phi(1),0)
            FreeEnergy = FreeEnergy + Gibbs_Fe_e_PbSn(c - g_phi(phi(1),0)*dc,T,1,0,approx)*(1.-g_phi(phi(1),0))
         endif
         

      else
         FreeEnergy =              Gibbs_Fe_l_PbSn(c,T,1,0,approx)*g_phi(phi(1),0)
         FreeEnergy = FreeEnergy + Gibbs_FE_e_PbSn(c,T,1,0,approx)*(1.-g_phi(phi(1),0))
      endif
      return
   elseif(lp.eq.2)then !Temperature Keep in J/mol

      Tk=T_k(T)
      
      !-(1-T*d/dT)dF/dphi
      !     if(heat_phi_term)then
      !  T (d/dT)dF/dphi
      
      FreeEnergy =  Tk*(Gibbs_Fe_l_PbSn(c,T,0,1,approx) - Gibbs_FE_e_PbSn(c,T,0,1,approx))

         ! -Df/Dphi
         FreeEnergy = FreeEnergy - Gibbs_Fe_l_PbSn(c,T,0,0,approx) + Gibbs_FE_e_PbSn(c,T,0,0,approx)
         g_latent_heat = FreeEnergy         
         FreeEnergy = FreeEnergy*phi_dot(1)*g_phi(phi(1),1)
         
         Q_heat(1) = FreeEnergy
 !     endif !end phi_dot_term
      if(heat_c_term)then 
         !   -(1- T(d/dT))dF/dc 
         FE=(Tk*Gibbs_Fe_l_PbSn(c,T,1,1,approx)*g_phi(phi(1),0)- Gibbs_Fe_l_PbSn(c,T,1,0,approx))*g_phi(phi(1),0)
         FE=FE+(Tk*Gibbs_FE_e_PbSn(c,T,1,1,approx)- Gibbs_FE_e_PbSn(c,T,1,0,approx))*(1.-g_phi(phi(1),0))
         FE=FE*c_dot
         FreeEnergy = FreeEnergy + FE
         Q_heat(2)= FE

      endif !c_dot_term
      !   heat capacity = -T (d/dT) dF/dT
      g_C =  -Tk*(g_phi(phi(1),0)*Gibbs_Fe_l_PbSn(c,T,0,2,approx)+g_phi(1d0-phi(1),0)*Gibbs_FE_e_PbSn(c,T,0,2,approx))

      return
 !     else
 !        write(*,*)approx
 !        stop
 !     endif
   elseif(lp.eq.-1)then !for anti-trapping f_c_phi
      FreeEnergy =              (Gibbs_FE_l_PbSn(c,T,1,0,.false.) - Gibbs_FE_e_PbSn(c,T,1,0,.false.))*g_phi(phi(1),1)
      return
   elseif(lp.eq.-2)then !for anti-trapping f_cc

      !resorting to approximation for gibbs_fe_e as it looks too complicated
      FreeEnergy =              0.5*(Gibbs_FE_l_PbSn(c+dh,T,1,0,.false.)-Gibbs_FE_l_PbSn(c-dh,T,1,0,.false.))/dh*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + 0.5*(Gibbs_FE_e_PbSn(c+dh,T,1,0,.false.)-Gibbs_FE_e_PbSn(c-dh,T,1,0,.false.))/dh*(1d0-g_phi(phi(1),0))

      return
   else
      write(*,*)"bad current_var",lp
      stop
   endif
elseif(AlNi.eq.6)then !NiCu                                                                                                                                                                                                                  
   if(lp.eq.0)then !free energy proper J/mol                                                                                                                                                                                                 
      FreeEnergy =              Gibbs_Fe_l_NiCu(c,T,0,0,approx)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_FE_e_NiCu(c,T,0,0,approx)*g_phi(1d0-phi(1),0)

      FreeEnergy = FreeEnergy*2d0/3d0

      return
   elseif(lp.eq.1)then !phase                                                                                                                                                                                                                




      FreeEnergy =              Gibbs_Fe_l_NiCu(c,T,0,0,approx)*g_phi(phi(1),1)
      FreeEnergy = FreeEnergy - Gibbs_FE_e_NiCu(c,T,0,0,approx)*g_phi(phi(1),1)
      FreeEnergy = FreeEnergy*2d0/3d0
      return
   elseif(lp.eq.3)then !Solute                                                                                                                                                                                                               

      FreeEnergy = FreeEnergy + Gibbs_Fe_l_NiCu(c,T,1,0,approx)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_FE_e_NiCu(c,T,1,0,approx)*(1.-g_phi(phi(1),0))
      FreeEnergy = FreeEnergy*2d0/3d0
      return
   elseif(lp.eq.-1)then !for anti-trapping f_c_phi
      FreeEnergy =              (Gibbs_FE_l_NiCu(c,T,1,0,.false.) - Gibbs_FE_e_NiCu(c,T,1,0,.false.))*g_phi(phi(1),1)
      FreeEnergy = FreeEnergy*2d0/3d0
      return
   elseif(lp.eq.-2)then !for anti-trapping f_cc

      !resorting to approximation for gibbs_fe_e as it looks too complicated
      FreeEnergy =              0.5*(Gibbs_FE_l_NiCu(c+dh,T,1,0,.false.)-Gibbs_FE_l_NiCu(c-dh,T,1,0,.false.))/dh*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + 0.5*(Gibbs_FE_e_NiCu(c+dh,T,1,0,.false.)-Gibbs_FE_e_NiCu(c-dh,T,1,0,.false.))/dh*(1d0-g_phi(phi(1),0))
      FreeEnergy = FreeEnergy*2d0/3d0

      return
   endif
elseif(AlNi.eq.7)then !pure Aluminium
   if(lp.eq.0)then !free energy proper J/mol
      FreeEnergy =              Gibbs_Fe_La(g_c0,T,0,0,approx)*g_phi(phi(1),0)
      FreeEnergy = FreeEnergy + Gibbs_FE_a(g_c0,T,0,0,approx)*g_phi(1d0-phi(1),0)
      return
   elseif(lp.eq.1)then
      FreeEnergy =              Gibbs_Fe_La(g_c0,T,0,0,approx)*g_phi(phi(1),1)
      FreeEnergy = FreeEnergy - Gibbs_FE_a(g_c0,T,0,0,approx)*g_phi(phi(1),1)
   elseif(lp.eq.2)then !Temperature Keep in J/mol

      Tk=T_k(T)
      
      !-(1-T*d/dT)dF/dphi
      !     if(heat_phi_term)then
      !  T (d/dT)dF/dphi
      
      FreeEnergy =  Tk*(Gibbs_Fe_la(g_c0,T,0,1,approx) - Gibbs_FE_a(g_c0,T,0,1,approx))


      ! -Df/Dphi
      FreeEnergy = FreeEnergy - Gibbs_Fe_la(g_c0,T,0,0,approx) + Gibbs_FE_a(g_c0,T,0,0,approx)
      g_latent_heat = FreeEnergy         
      FreeEnergy = FreeEnergy*phi_dot(1)*g_phi(phi(1),1)
      
      Q_heat(1) = FreeEnergy
      !     endif !end phi_dot_term
      !   heat capacity = -T (d/dT) dF/dT
!      g_C =  -Tk*(g_phi(phi(1),0)*Gibbs_Fe_la(g_c0,T,0,2,approx)+g_phi(1d0-phi(1),0)*Gibbs_FE_a(g_c0,T,0,2,approx))
      g_C=28.0

      return
   endif
elseif(AlNi.eq.8)then

   if(approx)then
      if(lp.eq.0)then
         FreeEnergy =   (0.5*g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))**2+g_AlFe_c(2,2)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,3) )*g_phi(phi(1),0)&
              +  (0.5*g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))**2+g_AlFe_c(1,2)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,3) )*g_phi(1-phi(1),0)
      elseif(lp.eq.1)then
         FreeEnergy =   (0.5*g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))**2+g_AlFe_c(2,2)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,3) )*g_phi(phi(1),1)&
              -  (0.5*g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))**2+g_AlFe_c(1,2)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,3) )*g_phi(phi(1),1)
      elseif(lp.eq.3)then
         FreeEnergy =   (g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,2) )*g_phi(phi(1),0)&
              +  (g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,2) )*g_phi(1-phi(1),0)
      elseif(lp.eq.-2)then !f_cc
         FreeEnergy = (1 - g_phi(phi(1),0))*g_AlFe_c(1,4) + g_phi(phi(1),0)*g_AlFe_c(2,4)

      elseif(lp.eq.-1)then !f_c_phi
         FreeEnergy = (g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,2) )*g_phi(phi(1),1)&
              -  (g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,2) )*g_phi(phi(1),1)

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
elseif(lp.eq.3)then !Solute


   AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,T,1,0)*g_phi(phi(1),0)
   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e2_AlFe(c,T,1,0)*(1-g_phi(phi(1),0))

elseif(lp.eq.2)then !Temperature Keep in J/mol
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

!end AlFe


!!!!!!!!!!!!!!!
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




double precision function EntropyMixing(c,d)
double precision, intent (in)::c
integer, intent (in)::d

if(c.eq.1d0.or.c.eq.0d0)then
   write(*,*)"c=",c, "in EntropyMixing"
   stop
endif
 
if(d.eq.1)then
   EntropyMixing = log(c/(1.-c))
elseif(d.eq.2)then
   EntropyMixing = 1d0/c + 1d0/(1d0-c)
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

double precision function h_phi(phi,d)
  double precision, intent (in):: phi
  integer, intent (in):: d
  if(d.eq.0)then
     h_phi = phi-phi**2
  elseif(d.eq.1)then
     h_phi =1d0 - 2d0*phi
  endif
end function h_phi

double precision function u_phi(phi,gamma)
double precision, intent (in):: phi,gamma
double precision h_phi
u_phi = phi + gamma * h_phi(phi,0)**2
end function u_phi

double precision function g_phi(phi,d)
use solution_parameters
  double precision, intent (in):: phi
  integer, intent (in):: d
  double precision :: n,h_phi,u_phi,u,h
  g_phi=0d0
  if(phi.lt.1d-16)return
  if(g_skew_potential)then
     if(d.eq.1)then
        g_phi =  ((-1/(2-phi)/2-1/phi/2)*(phi**2-phi**3/3)+(log(2-phi)/2-log(phi)/2)*(-phi**2+2*phi)-1.D0/3.D0+phi/3+2.D0/3.D0/(2-phi))&
             /(2.D0/3.D0*log(2.D0)-1.D0/6.D0)
     elseif(d.eq.2)then
        g_phi = ((-1/((2-phi)**2)/2+1/phi**2/2)*(phi**2-phi**3/3)+2*(-1/(2-phi)/2-1/phi/2)*(-phi**2+2*phi)&
             +(log(2-phi)/2-log(phi)/2)*(-2*phi+2)+1.D0/3.D0+2.D0/3.D0/(2-phi)**2)/(2.D0/3.D0*log(2.D0)-1.D0/6.D0)
     else !assume 0
        g_phi = ((log(2-phi)/2-log(phi)/2)*(phi**2-phi**3/3)-phi/3&
             +phi**2/6-2.D0/3.D0*log(2-phi)+2.D0/3.D0*log(2.D0))/(2.D0/3.D0*log(2.D0)-1.D0/6.D0)
     endif
  elseif(g_power_skew)then
     n = g_power 
     if(d.eq.1)then
        g_phi = 0.25*phi*(n+4)*(n+2)/(1-phi)*(1-phi)**(0.5*n+1) 
     elseif(d.eq.2)then
        g_phi = -0.125*(n*phi+2*phi-2)*(n+4)*(n+2)/(1-phi)**2*(1-phi)**(0.5*n+1)
     else !assume 0
        g_phi = 1-0.5*((2+n)*phi+2)*(1-phi)**(0.5*n+1)
     endif
!
  elseif(g_s_skew)then
     if(d.eq.0)then

!AMM
        if(g_AMM.eq.1)then
           g_phi = -0.1485934166D1*phi**2+0.3268185149D1*phi-0.3258617531D0*log(2544018011.D0+0.2551487495D11*phi)+0.7057191429D1
        elseif(g_AMM.eq.2)then
           g_phi = -0.2D1*phi**3+0.3D1*phi**2
        else
           write(*,*)"g_AMM is 1 or 2"
        endif


     elseif(d.eq.1)then 
!AMM

        if(g_AMM.eq.1)then
           g_phi = -0.1516536976D3*phi*(-1.D0+phi)/(0.5088036022D1+0.510297499D2*phi)
        elseif(g_AMM.eq.2)then
           g_phi = -0.6D1*phi*(-1.D0+phi)
        else
           write(*,*)"g_AMM is 1 or 2"
        endif



     endif

  else
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
!      write*,*)"d y1/dc not determined at c=",c
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
   write(*,*)"c out of bound in y4()"
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
  elseif(d_c.eq.2.and.d_T.eq.0)then
     Tk=T_k(T)
     glA=gibbs_T(hlA_,Tk,0)
     glB=gibbs_T(hlB_,Tk,0)
     call MixingVector(RK,c,2)
     RKl_=         (hrkl0(1)+hrkl0(2)*Tk)*RK(1)
     RKl_=RKl_ + (hrkl1(1)+hrkl1(2)*Tk)*RK(2)
     RKl_=RKl_ + (hrkl2(1)+hrkl2(2)*Tk)*RK(3)
     RKl_=RKl_ + (hrkl3(1)+hrkl3(2)*Tk)*RK(4)
     RKl_=RKl_ + (hrkl4(1)+hrkl4(2)*Tk)*RK(5)
     ent_=g_R*Tk*EntropyMixing(c,2)
     gl_=-glA+glB+ent_+RKl_
     Gibbs_FE_l=gl_
  else
     write(*,*)"unexpected d_c,d_T combination",d_c,d_T
     stop
  endif

  end function Gibbs_FE_l


double precision function Gibbs_FE_La(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: T_k
  double precision :: Tk,Gibbs_T
  double precision :: hlA_(8)
  Tk=T_k(T)
  if (Tk > 1728.0) then
     hlA_=hlA3 
  elseif (Tk > 933.47)then
     hlA_=hlA3 
  elseif (Tk > 700.0)then
     hlA_=hlA2 
  else
     hlA_=hlA1 
  endif

     Gibbs_FE_La=gibbs_T(hlA_,Tk,d_T)

!!$  if(d_T.eq.0)then
!!$     Gibbs_FE_La=gibbs_T(hlA_,Tk,0)
!!$  elseif(d_T.eq.1)then
!!$     Gibbs_FE_La=gibbs_T(hlA_,Tk,1)
!!$  elseif(d_T.eq.2)then
!!$     Gibbs_FE_La=gibbs_T(hlA_,Tk,2)
!!$  else
!!$     write(*,*)"unexpected d_c,d_T combination",d_c,d_T
!!$     stop
!!$  endif

  end function Gibbs_FE_La


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



double precision function Gibbs_FE_a(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in)::d_c,d_T
  double precision :: T_k
  double precision ::Tk,Gibbs_T
  double precision :: haA_(8)

 Tk=T_k(T)
  if (Tk > 1728.0) then
  haA_=haA3 
elseif (Tk > 933.47)then
  haA_=haA3 
elseif (Tk > 700.0)then
  haA_=haA2 
else
  haA_=haA1 
endif

     Tk=T_k(T)
     Gibbs_FE_a=gibbs_T(haA_,Tk,d_T)
!!$  if(d_T.eq.0)then
!!$     Tk=T_k(T)
!!$     Gibbs_FE_a=gibbs_T(haA_,Tk,0)
!!$
!!$  elseif(d_T.eq.1)then
!!$     Tk=T_k(T)
!!$     Gibbs_FE_a=gibbs_T(haA_,Tk,1)
!!$
!!$  elseif(d_T.eq.2)then
!!$     Tk=T_k(T)
!!$     Gibbs_FE_a=gibbs_T(haA_,Tk,2)
!!$  endif
end function Gibbs_FE_a
!begin AlSi

!y:=(x,y1,y2,dy1)->(1-x)*((1-x)*y1+x*(dy1/2+y1))+x*(x*y2+(1-x)*(dy1/2+y1));
!Y:=(c,x1,x2,y1,y2,dy1)->y((c-x1)/(x2-x1),y1,y2,dy1*(x2-x1));

double precision function Bezier_1(x,y1,y2,dy1,d)
double precision, intent(in):: x,y1,y2,dy1
integer, intent(in) ::d
if(d.eq.0)then
   Bezier_1 = (1d0-x)*((1-x)*y1+x*(dy1/2+y1))+x*(x*y2+(1-x)*(dy1/2+y1))
else
   Bezier_1 = (1d0-x)*dy1+x*(2d0*(y2-y1)-dy1)
endif
end function Bezier_1

double Precision function Bezier_2(c,x1,x2,y1,y2,dy1,d)
double precision, intent(in)::c,x1,x2,y1,y2,dy1
integer, intent(in) ::d
double precision Bezier_1
if(d.eq.0)then
   Bezier_2 = Bezier_1((c-x1)/(x2-x1),y1,y2,dy1*(x2-x1),0)
else
   Bezier_2 = Bezier_1((c-x1)/(x2-x1),y1,y2,dy1*(x2-x1),1)/(x2-x1)
endif
end function Bezier_2
double precision function Quad(c,c0,f,fc,fcc,d)
double precision, intent(in)::c,c0,f,fc,fcc
integer, intent(in)::d
if(d.eq.0)then
   Quad = 0.5*fcc*(c-c0)**2 + fc*(c-c0) + f
else
   Quad = fcc*(c-c0) + fc
endif
end function Quad

double precision function Gibbs_FE_l_AlSi(c,T,d_c,d_T,approx)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  logical, intent(in):: approx
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5), gs1A,gs1B, gs2A,gs2B
  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gl_,Tk,Gibbs_T
  double precision :: hlA_(8),hlB_(8),hs1A_(8),hs1B_(8),hs2A_(8),hs2B_(8)
  double precision :: Bezier_2,Quad
  integer :: T_low,T_high
  double precision :: T_weight


! AlSi below
if(approx)then
   write(*,*)"approx routine not finishes yet:Gibbs_FE_l_AlSi"
   stop
   Tk=T_k(T)   
   if(d_c.eq.0.and.d_T.eq.0)then
      Gibbs_FE_l_AlSi =  Bezier_2(c,g_c_min,g_c_mid,g_FE(1,1)+T*g_dfdT0(1,1)+0.5*T**2*g_d2fdT0(1,1),g_FE(1,2)+T*g_dfdT0(1,2)+0.5*T**2*g_d2fdT0(1,2),g_f_c,0)
   elseif(d_c.eq.1.and.d_T.eq.0)then
      Gibbs_FE_l_AlSi =  Bezier_2(c,g_c_min,g_c_mid,g_FE(1,1)+T*g_dfdT0(1,1)+0.5*T**2*g_d2fdT0(1,1),g_FE(1,2)+T*g_dfdT0(1,2)+0.5*T**2*g_d2fdT0(1,2),g_f_c,1)
   elseif(d_c.eq.0.and.d_T.eq.1)then
      Gibbs_FE_l_AlSi =   (c-g_c_min)/(g_c_mid-g_c_min)*g_dfdT0(1,1) + (g_c_mid - c)/(g_c_mid-g_c_min)* g_dfdT0(1,2)
   elseif(d_c.eq.1.and.d_T.eq.1)then
      Gibbs_FE_l_AlSi =   (c-g_c_min)/(g_c_mid-g_c_min)*g_dfdcT0(1,1) + (g_c_mid - c)/(g_c_mid-g_c_min)* g_dfdcT0(1,2)
   elseif(d_c.eq.0.and.d_T.eq.2)then
      Gibbs_FE_l_AlSi =      (c-g_c_min)/(g_c_mid-g_c_min)*g_d2fdT0(1,1) + (g_c_mid - c)/(g_c_mid-g_c_min)* g_d2fdT0(1,2)
   else
      write(*,*)"in g_approx for T"
      stop
      Gibbs_FE_l_AlSi =  0d0
   endif
else


Tk=T_k(T)
if (Tk .gt. 1687.0) then
  hlA_=AlSihlA3 
  hlB_=AlSihlB2 
elseif(Tk .gt. 933.47) then
  hlA_=AlSihlA3 
  hlB_=AlSihlB1 
elseif (Tk .gt. 700.0) then
  hlA_=AlSihlA2 
  hlB_=AlSihlB1 
else
  hlA_=AlSihlA1 
  hlB_=AlSihlB1 
endif

if(d_c.eq.0.and.d_T.eq.0)then
   glA=gibbs_T(hlA_,Tk,0)
   glB=gibbs_T(hlB_,Tk,0)
   call MixingVector(RK,c,0)
   RKl_=         (AlSihrkl0(1)+AlSihrkl0(2)*Tk)*RK(1)
   RKl_=RKl_ + (AlSihrkl1(1)+AlSihrkl1(2)*Tk)*RK(2)
   RKl_=RKl_ + (AlSihrkl2(1)+AlSihrkl2(2)*Tk)*RK(3)
   
   ent_=g_R*Tk*EntropyMixing(c,0)
   gl_=glA*(1d0-c)+glB*c+ent_+RKl_
   Gibbs_FE_l_AlSi=gl_
elseif(d_c.eq.0.and.d_T.eq.1)then
   glA=gibbs_T(hlA_,Tk,1)
   glB=gibbs_T(hlB_,Tk,1)
   call MixingVector(RK,c,0)
   RKl_=         AlSihrkl0(2)*RK(1)
   RKl_=RKl_ + AlSihrkl1(2)*RK(2)
   RKl_=RKl_ + AlSihrkl2(2)*RK(3)
   
   ent_=g_R*EntropyMixing(c,0)
   gl_=glA*(1d0-c)+glB*c+ent_+RKl_
   Gibbs_FE_l_AlSi=gl_
elseif(d_c.eq.0.and.d_T.eq.2)then
   glA=gibbs_T(hlA_,Tk,2)
   glB=gibbs_T(hlB_,Tk,2)
   gl_=glA*(1d0-c)+glB*c
   Gibbs_FE_l_AlSi=gl_
elseif(d_c.eq.1.and.d_T.eq.0)then
   glA=gibbs_T(hlA_,Tk,0)
   glB=gibbs_T(hlB_,Tk,0)
   call MixingVector(RK,c,1)
   RKl_=         (AlSihrkl0(1)+AlSihrkl0(2)*Tk)*RK(1)
   RKl_=RKl_ + (AlSihrkl1(1)+AlSihrkl1(2)*Tk)*RK(2)
   RKl_=RKl_ + (AlSihrkl2(1)+AlSihrkl2(2)*Tk)*RK(3)
   
   ent_=g_R*Tk*EntropyMixing(c,1)
   gl_=-glA+glB+ent_+RKl_
   Gibbs_FE_l_AlSi=gl_
elseif(d_c.eq.1.and.d_T.eq.1)then
   glA=gibbs_T(hlA_,Tk,1)
   glB=gibbs_T(hlB_,Tk,1)
   call MixingVector(RK,c,1)
   RKl_=       AlSihrkl0(2)*RK(1)
   RKl_=RKl_ + AlSihrkl1(2)*RK(2)
   RKl_=RKl_ + AlSihrkl2(2)*RK(3)
   
   ent_=g_R*EntropyMixing(c,1)
   gl_=-glA+glB+ent_+RKl_
   Gibbs_FE_l_AlSi=gl_
elseif(d_c.eq.1.and.d_T.eq.2)then
   glA=gibbs_T(hlA_,Tk,2)
   glB=gibbs_T(hlB_,Tk,2)
   gl_=-glA+glB
   Gibbs_FE_l_AlSi=gl_
elseif(d_c.eq.2.and.d_T.eq.0)then
   call MixingVector(RK,c,2)
   RKl_=         (AlSihrkl0(1)+AlSihrkl0(2)*Tk)*RK(1)
   RKl_=RKl_ + (AlSihrkl1(1)+AlSihrkl1(2)*Tk)*RK(2)
   RKl_=RKl_ + (AlSihrkl2(1)+AlSihrkl2(2)*Tk)*RK(3)
   
   ent_=g_R*Tk*EntropyMixing(c,2)
   Gibbs_FE_l_AlSi=ent_+RKl_
else
   write(*,*)"unexpected d_c,d_T combination",d_c,d_T
   stop
endif
endif
end function Gibbs_FE_l_AlSi


double precision function Gibbs_FE_e_AlSi(c,T,d_c,d_T,approx)
  use solution_parameters
  use time_dep_parameters
  implicit none
  double precision, intent(in) :: c,T
  logical, intent(in)::approx
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5), gs1A,gs1B, gs2A,gs2B

  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gs2_,Tk,Gibbs_T
  double precision :: hlA_(8),hlB_(8),hs1A_(8),hs1B_(8),hs2A_(8),hs2B_(8)
  double precision :: t_factor
  double precision :: Bezier_2,Quad
  integer :: T_low,T_high
  double precision :: T_weight


! AlSi below
if(approx)then


   Tk=T_k(T)   
   if(d_c.eq.0.and.d_T.eq.0)then
      Gibbs_FE_e_AlSi =  Bezier_2(c,g_c_max,g_c_mid,g_FE(2,1)+T*g_dfdT0(2,1)+0.5*T**2*g_d2fdT0(2,1),g_FE(2,2)+T*g_dfdT0(2,2)+0.5*T**2*g_d2fdT0(2,2),g_f_c,0)
   elseif(d_c.eq.1.and.d_T.eq.0)then
      Gibbs_FE_e_AlSi =  Bezier_2(c,g_c_max,g_c_mid,g_FE(2,1)+T*g_dfdT0(2,1)+0.5*T**2*g_d2fdT0(2,1),g_FE(2,2)+T*g_dfdT0(2,2)+0.5*T**2*g_d2fdT0(2,2),g_f_c,1)
   elseif(d_c.eq.0.and.d_T.eq.1)then
      Gibbs_FE_e_AlSi =   (c-g_c_max)/(g_c_mid-g_c_max)*g_dfdT0(2,1) + (g_c_mid - c)/(g_c_mid-g_c_max)* g_dfdT0(2,2)
   elseif(d_c.eq.1.and.d_T.eq.1)then
      Gibbs_FE_e_AlSi =   (c-g_c_max)/(g_c_mid-g_c_max)*g_dfdcT0(2,1) + (g_c_mid - c)/(g_c_mid-g_c_max)* g_dfdcT0(2,2)
   elseif(d_c.eq.0.and.d_T.eq.2)then
      Gibbs_FE_e_AlSi =   (c-g_c_max)/(g_c_mid-g_c_max)*g_d2fdT0(2,1) + (g_c_mid - c)/(g_c_mid-g_c_max)* g_d2fdT0(2,2)
else
   write(*,*)"in Gibbs_FE_e_AlSi",d_c,d_T
   stop
   Gibbs_FE_e_AlSi =  0d0
endif
else



Tk=T_k(T)
if (Tk .gt. 1687.0) then
  hs2A_=AlSihs2A3 
  hs2B_=AlSihs2B2 
elseif(Tk .gt. 933.47) then
  hs2A_=AlSihs2A3 
  hs2B_=AlSihs2B1 
elseif (Tk .gt. 700.0) then
  hs2A_=AlSihs2A2 
  hs2B_=AlSihs2B1 
else
  hs2A_=AlSihs2A1 
  hs2B_=AlSihs2B1 
endif

if(d_c.eq.0.and.d_T.eq.0)then
   gs2A=gibbs_T(hs2A_,Tk,0)
   gs2B=gibbs_T(hs2B_,Tk,0)
   call MixingVector(RK,c,0)
   RKl_=       (AlSihrks20(1)+AlSihrks20(2)*Tk)*RK(1)
   RKl_=RKl_ + (AlSihrks21(1)+AlSihrks21(2)*Tk)*RK(2)
   RKl_=RKl_ + (AlSihrks22(1)+AlSihrks22(2)*Tk)*RK(3)
   ent_=g_R*Tk*EntropyMixing(c,0)
   gs2_=gs2A*(1d0-c)+gs2B*c+ent_+RKl_
   Gibbs_FE_e_AlSi=gs2_
elseif(d_c.eq.0.and.d_T.eq.1)then
   gs2A=gibbs_T(hs2A_,Tk,1)
   gs2B=gibbs_T(hs2B_,Tk,1)
   call MixingVector(RK,c,0)
   RKl_=       AlSihrks20(2)*RK(1)
   RKl_=RKl_ + AlSihrks21(2)*RK(2)
   RKl_=RKl_ + AlSihrks22(2)*RK(3)
   ent_=g_R*EntropyMixing(c,0)
   gs2_=gs2A*(1d0-c)+gs2B*c+ent_+RKl_
   Gibbs_FE_e_AlSi=gs2_
elseif(d_c.eq.0.and.d_T.eq.2)then
   gs2A=gibbs_T(hs2A_,Tk,2)
   gs2B=gibbs_T(hs2B_,Tk,2)
   gs2_=gs2A*(1d0-c)+gs2B*c
   Gibbs_FE_e_AlSi=gs2_
elseif(d_c.eq.1.and.d_T.eq.0)then
   gs2A=gibbs_T(hs2A_,Tk,0)
   gs2B=gibbs_T(hs2B_,Tk,0)
   call MixingVector(RK,c,1)
   RKl_=       (AlSihrks20(1)+AlSihrks20(2)*Tk)*RK(1)
   RKl_=RKl_ + (AlSihrks21(1)+AlSihrks21(2)*Tk)*RK(2)
   RKl_=RKl_ + (AlSihrks22(1)+AlSihrks22(2)*Tk)*RK(3)
   ent_=g_R*Tk*EntropyMixing(c,1)
   gs2_=-gs2A+gs2B+ent_+RKl_
   Gibbs_FE_e_AlSi=gs2_
elseif(d_c.eq.1.and.d_T.eq.1)then
   gs2A=gibbs_T(hs2A_,Tk,1)
   gs2B=gibbs_T(hs2B_,Tk,1)
   call MixingVector(RK,c,1)
   RKl_=       AlSihrks20(2)*RK(1)
   RKl_=RKl_ + AlSihrks21(2)*RK(2)
   RKl_=RKl_ + AlSihrks22(2)*RK(3)
   ent_=g_R*EntropyMixing(c,1)
   gs2_=-gs2A+gs2B+ent_+RKl_
   Gibbs_FE_e_AlSi=gs2_
elseif(d_c.eq.1.and.d_T.eq.2)then
   gs2A=gibbs_T(hs2A_,Tk,2)
   gs2B=gibbs_T(hs2B_,Tk,2)
   gs2_=-gs2A+gs2B
   Gibbs_FE_e_AlSi=gs2_
elseif(d_c.eq.2.and.d_T.eq.0)then
   call MixingVector(RK,c,2)
   RKl_=       (AlSihrks20(1)+AlSihrks20(2)*Tk)*RK(1)
   RKl_=RKl_ + (AlSihrks21(1)+AlSihrks21(2)*Tk)*RK(2)
   RKl_=RKl_ + (AlSihrks22(1)+AlSihrks22(2)*Tk)*RK(3)
   ent_=g_R*Tk*EntropyMixing(c,2)
   Gibbs_FE_e_AlSi=ent_+RKl_
else
   write(*,*)"unexpected d_c,d_T combination",d_c,d_T
   stop
endif
endif
end function Gibbs_FE_e_AlSi

!end AlSi


double precision function Gibbs_FE_l_NiCu(c,T,d_c,d_T,approx)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  logical, intent(in):: approx
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5), gs1A,gs1B, gs2A,gs2B
  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gl_,Tk,Gibbs_T

Tk=T_k(T)
if(d_c.eq.0.and.d_T.eq.0)then

   Gibbs_FE_l_NiCu = 0d0 

   Gibbs_FE_l_NiCu = Gibbs_FE_l_NiCu + g_R*Tk*EntropyMixing(c,0)
elseif(d_c.eq.1.and.d_T.eq.0)then
   Gibbs_FE_l_NiCu = Gibbs_FE_l_NiCu + g_R*Tk*EntropyMixing(c,1)
else


endif
end function Gibbs_FE_l_NiCu



double precision function Gibbs_FE_e_NiCu(c,T,d_c,d_T,approx)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  logical, intent(in):: approx
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5), gs1A,gs1B, gs2A,gs2B
  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gl_,Tk,Gibbs_T

Tk=T_k(T)
if(d_c.eq.0.and.d_T.eq.0)then
   Gibbs_FE_e_NiCu = g_W(1)/ g_MolPerVol(1)*(1-c)*1.6*(Tk-g_T(1))/g_T(1)+g_W(2)/ g_MolPerVol(2)*c*2.0*(Tk-g_T(2))/g_T(2)
   Gibbs_FE_e_NiCu = Gibbs_FE_e_NiCu*2d0*g_lambda + g_R*Tk*EntropyMixing(c,0)
elseif(d_c.eq.1.and.d_T.eq.0)then
   Gibbs_FE_e_NiCu = g_W(1)/ g_MolPerVol(1)*(-1)*1.6*(Tk-g_T(1))/g_T(1)+g_W(2)/ g_MolPerVol(2)*1*2.0*(Tk-g_T(2))/g_T(2)
   Gibbs_FE_e_NiCu = Gibbs_FE_e_NiCu*2d0*g_lambda + g_R*Tk*EntropyMixing(c,1)
else
   write(*,*)"d_c,d_T=",d_c,d_T
   stop
endif
end function Gibbs_FE_e_NiCu


double precision function Gibbs_FE_l_PbSn(c,T,d_c,d_T,approx)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  logical, intent(in):: approx
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5), gs1A,gs1B, gs2A,gs2B
  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gl_,Tk,Gibbs_T


  if(g_quad_extension.and.c.ge.g_c_max)then
     if(d_c.eq.0)then
        Gibbs_FE_l_PbSn = g_f_c_max + g_f1_c_max*(c-g_c_max)+0.5*g_f2_c_max*(c-g_c_max)**2
     elseif(d_c.eq.1)then
        Gibbs_FE_l_PbSn = g_f1_c_max+g_f2_c_max*(c-g_c_max)
     elseif(d_c.eq.2)then
        Gibbs_FE_l_PbSn = g_f2_c_max
     endif
     return
  endif



! PbSn
if(approx)then
   write(*,*)"approx routine not finished yet:Gibbs_FE_l_PbSn"
   stop
endif


Tk=T_k(T)

!!!!!!!!!!
!!$  call GibbsVector(V,Tk,0)
!!$  call MixingVector(RK,c,0)
!!$   FreeEnergy=0.
!!$   ha=0.
!!$   hb=0.
!!$   hr=0.
!!$   do i=1,n_phases
!!$      do j = 1,8
!!$         ha(i) = ha(i) + hsa(i,j)*V(j)
!!$         hb(i) = hb(i) + hsb(i,j)*V(j)
!!$      enddo
!!$   enddo
!!$   do i=1,n_phases
!!$      do j = 1,4
!!$         hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
!!$      enddo
!!$   enddo
!!$   S=0d0
!!$   do i=1,n_phases
!!$      S=S+g_phi(phi(i),0)
!!$   enddo
!!$   do i=1,n_phases
!!$      FreeEnergy = FreeEnergy + (c*hb(i)+(1.-c)*ha(i)+hr(i))*g_phi(phi(i),0)/S
!!$   enddo
!!$   FreeEnergy = FreeEnergy + R*Tk*EntropyMixing(c,0)
!!$   return
!!!!!!!




if(d_c.eq.0.and.d_T.eq.0)then
   glA=gibbs_T(PbSnhsa1,Tk,0)
   glB=gibbs_T(PbSnhsb1,Tk,0)
   call MixingVector(RK,c,0)
   RKl_=         (PbSnhsrk1(1)+PbSnhsrk1(2)*Tk)*RK(1)
   RKl_=  RKl_ + (PbSnhsrk1(3)+PbSnhsrk1(4)*Tk)*RK(2)
   ent_=g_R*Tk*EntropyMixing(c,0)
   gl_=glA*(1d0-c)+glB*c+ent_+RKl_
   Gibbs_FE_l_PbSn=gl_
elseif(d_c.eq.0.and.d_T.eq.1)then
   glA=gibbs_T(PbSnhsa1,Tk,1)
   glB=gibbs_T(PbSnhsb1,Tk,1)
   call MixingVector(RK,c,0)
   RKl_=         PbSnhsrk1(2)*RK(1)
   RKl_=RKl_ + PbSnhsrk1(2)*RK(2)
   ent_=g_R*EntropyMixing(c,0)
   gl_=glA*(1d0-c)+glB*c+ent_+RKl_
   Gibbs_FE_l_PbSn=gl_
elseif(d_c.eq.0.and.d_T.eq.2)then
   glA=gibbs_T(PbSnhsa1,Tk,2)
   glB=gibbs_T(PbSnhsb1,Tk,2)
   gl_=glA*(1d0-c)+glB*c
   Gibbs_FE_l_PbSn=gl_
elseif(d_c.eq.1.and.d_T.eq.0)then
   glA=gibbs_T(PbSnhsa1,Tk,0)
   glB=gibbs_T(PbSnhsb1,Tk,0)
   call MixingVector(RK,c,1)
   RKl_=         (PbSnhsrk1(1)+PbSnhsrk1(2)*Tk)*RK(1)
   RKl_=  RKl_ + (PbSnhsrk1(3)+PbSnhsrk1(4)*Tk)*RK(2)
   ent_=g_R*Tk*EntropyMixing(c,1)
   gl_=-glA+glB+ent_+RKl_
   Gibbs_FE_l_PbSn=gl_
elseif(d_c.eq.1.and.d_T.eq.1)then
   glA=gibbs_T(PbSnhsa1,Tk,1)
   glB=gibbs_T(PbSnhsb1,Tk,1)
   call MixingVector(RK,c,1)
   RKl_=       PbSnhsrk1(2)*RK(1)
   RKl_=RKl_ + PbSnhsrk1(2)*RK(2)
   ent_=g_R*EntropyMixing(c,1)
   gl_=-glA+glB+ent_+RKl_
   Gibbs_FE_l_PbSn=gl_
elseif(d_c.eq.1.and.d_T.eq.2)then
   glA=gibbs_T(PbSnhsa1,Tk,2)
   glB=gibbs_T(PbSnhsb1,Tk,2)
   gl_=-glA+glB
   Gibbs_FE_l_PbSn=gl_
elseif(d_c.eq.2.and.d_T.eq.0)then
   call MixingVector(RK,c,2)
   RKl_=         (PbSnhsrk1(1)+PbSnhsrk1(2)*Tk)*RK(1)
   RKl_=  RKl_ + (PbSnhsrk1(3)+PbSnhsrk1(4)*Tk)*RK(2)
   ent_=g_R*Tk*EntropyMixing(c,2)
   Gibbs_FE_l_PbSn=ent_+RKl_
else
   write(*,*)"unexpected d_c,d_T combination",d_c,d_T
   stop
endif

end function Gibbs_FE_l_PbSn


double precision function Gibbs_FE_e_PbSn(c,T,d_c,d_T,approx)
  use solution_parameters
  use time_dep_parameters
  implicit none
  double precision, intent(in) :: c,T
  logical, intent(in)::approx
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5), gs1A,gs1B, gs2A,gs2B

  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gs2_,Tk,Gibbs_T



  if(g_quad_extension.and.c.le.g_c_min)then
     if(d_c.eq.0)then
        Gibbs_FE_e_PbSn = g_f_c_min + g_f1_c_min*(c-g_c_min)+0.5*g_f2_c_min*(c-g_c_min)**2
     elseif(d_c.eq.1)then
        Gibbs_FE_e_PbSn = g_f1_c_min+g_f2_c_min*(c-g_c_min)
     elseif(d_c.eq.2)then
        Gibbs_FE_e_PbSn = g_f2_c_min
     endif
     return
  endif
  
! AlSi below
if(approx)then
   write(*,*)"in approx, which is not coded for PbSn"
   stop

endif




Tk=T_k(T)


if(d_c.eq.0.and.d_T.eq.0)then
   gs2A=gibbs_T(PbSnhsa2,Tk,0)
   gs2B=gibbs_T(PbSnhsb2,Tk,0)
   call MixingVector(RK,c,0)
   RKl_= (PbSnhsrk2(1)+PbSnhsrk2(2)*Tk)*RK(1)

   ent_=g_R*Tk*EntropyMixing(c,0)
   gs2_=gs2A*(1d0-c)+gs2B*c+ent_+RKl_
   Gibbs_FE_e_PbSn=gs2_

elseif(d_c.eq.0.and.d_T.eq.1)then
   gs2A=gibbs_T(PbSnhsa2,Tk,1)
   gs2B=gibbs_T(PbSnhsb2,Tk,1)
   call MixingVector(RK,c,0)
   RKl_=       PbSnhsrk2(2)*RK(1)
   ent_=g_R*EntropyMixing(c,0)
   gs2_=gs2A*(1d0-c)+gs2B*c+ent_+RKl_
   Gibbs_FE_e_PbSn=gs2_
elseif(d_c.eq.0.and.d_T.eq.2)then
   gs2A=gibbs_T(PbSnhsa2,Tk,2)
   gs2B=gibbs_T(PbSnhsb2,Tk,2)
   gs2_=gs2A*(1d0-c)+gs2B*c
   Gibbs_FE_e_PbSn=gs2_
elseif(d_c.eq.1.and.d_T.eq.0)then
   gs2A=gibbs_T(PbSnhsa2,Tk,0)
   gs2B=gibbs_T(PbSnhsb2,Tk,0)
   call MixingVector(RK,c,1)
   RKl_= (PbSnhsrk2(1)+PbSnhsrk2(2)*Tk)*RK(1)
   ent_=g_R*Tk*EntropyMixing(c,1)
   gs2_=-gs2A+gs2B+ent_+RKl_
   Gibbs_FE_e_PbSn=gs2_
elseif(d_c.eq.1.and.d_T.eq.1)then
   gs2A=gibbs_T(PbSnhsa2,Tk,1)
   gs2B=gibbs_T(PbSnhsb2,Tk,1)
   call MixingVector(RK,c,1)
   RKl_= PbSnhsrk2(2)*RK(1)
   ent_=g_R*EntropyMixing(c,1)
   gs2_=-gs2A+gs2B+ent_+RKl_
   Gibbs_FE_e_PbSn=gs2_
elseif(d_c.eq.1.and.d_T.eq.2)then
   gs2A=gibbs_T(PbSnhsa2,Tk,2)
   gs2B=gibbs_T(PbSnhsb2,Tk,2)
   gs2_=-gs2A+gs2B
   Gibbs_FE_e_PbSn=gs2_
elseif(d_c.eq.2.and.d_T.eq.0)then
   call MixingVector(RK,c,2)
   RKl_= (PbSnhsrk2(1)+PbSnhsrk2(2)*Tk)*RK(1)
   ent_=g_R*Tk*EntropyMixing(c,2)
   Gibbs_FE_e_PbSn=ent_+RKl_
else
   write(*,*)"unexpected d_c,d_T combination",d_c,d_T
   stop
endif

end function Gibbs_FE_e_PbSn




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

double precision function Gibbs_SiGeliq(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: glSi,glGe,RK(5)
  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gl_,Tk,Gibbs_T
  double precision :: hlSi(8),hlGe(8)


  Tk=T_k(T)
  if (Tk > 1687.0)then
     hlSi = hlSi2 ;
     hlGe = hlGe3 ;
  elseif (Tk > 1211.4)then
     hlSi = hlSi1 ;
     hlGe = hlGe3 ;
  elseif (Tk > 900.0)then
     hlSi = hlSi1 ;
     hlGe = hlGe2 ;
  else
     hlSi = hlSi1 ;
     hlGe = hlGe1 ;
  endif
  
  
  if(d_c.eq.0.and.d_T.eq.0)then
     glSi=gibbs_T(hlSi,Tk,0)
     glGe=gibbs_T(hlGe,Tk,0)
     call MixingVector(RK,c,0)
     RKl_=         (hrklSiGe0(1)+hrklSiGe0(2)*Tk)*RK(1)
     RKl_=RKl_ + (hrklSiGe1(1)+hrklSiGe1(2)*Tk)*RK(2)
     RKl_=RKl_ + (hrklSiGe2(1)+hrklSiGe2(2)*Tk)*RK(3)
     ent_=g_R*Tk*EntropyMixing(c,0)
     gl_=glSi*(1d0-c)+glGe*c+ent_+RKl_
     Gibbs_SiGeliq=gl_
  elseif(d_c.eq.0.and.d_T.eq.1)then
     Tk=T_k(T)
     glSi=gibbs_T(hlSi,Tk,1)
     glGe=gibbs_T(hlGe,Tk,1)
     call MixingVector(RK,c,0)
     RKl_=         hrklSiGe0(2)*RK(1)
     RKl_=RKl_ +   hrklSiGe1(2)*RK(2)
     RKl_=RKl_ +   hrklSiGe2(2)*RK(3)
     ent_=g_R*EntropyMixing(c,0)
     gl_=glSi*(1d0-c)+glGe*c+ent_+RKl_
     Gibbs_SiGeliq=gl_
  elseif(d_c.eq.0.and.d_T.eq.2)then
     Tk=T_k(T)
     glSi=gibbs_T(hlSi,Tk,2)
     glGe=gibbs_T(hlGe,Tk,2)
     call MixingVector(RK,c,0)
     gl_=glSi*(1d0-c)+glGe*c
     Gibbs_SiGeliq=gl_
  elseif(d_c.eq.1.and.d_T.eq.0)then
     Tk=T_k(T)
     glSi=gibbs_T(hlSi,Tk,0)
     glGe=gibbs_T(hlGe,Tk,0)
     call MixingVector(RK,c,1)
     RKl_=         (hrklSiGe0(1)+hrklSiGe0(2)*Tk)*RK(1)
     RKl_=RKl_ + (hrklSiGe1(1)+hrklSiGe1(2)*Tk)*RK(2)
     RKl_=RKl_ + (hrklSiGe2(1)+hrklSiGe2(2)*Tk)*RK(3)
     ent_=g_R*Tk*EntropyMixing(c,1)
     gl_=-glSi+glGe+ent_+RKl_
     Gibbs_SiGeliq=gl_
  elseif(d_c.eq.1.and.d_T.eq.1)then
     Tk=T_k(T)
     glSi=gibbs_T(hlSi,Tk,1)
     glGe=gibbs_T(hlGe,Tk,1)
     call MixingVector(RK,c,0)
     RKl_=         hrklSiGe0(2)*RK(1)
     RKl_=RKl_ +   hrklSiGe1(2)*RK(2)
     RKl_=RKl_ +   hrklSiGe2(2)*RK(3)
     ent_=g_R*EntropyMixing(c,1)
     gl_=-glSi+glGe+ent_+RKl_
     Gibbs_SiGeliq=gl_
  elseif(d_c.eq.1.and.d_T.eq.2)then
     Tk=T_k(T)
     glSi=gibbs_T(hlSi,Tk,2)
     glGe=gibbs_T(hlGe,Tk,2)
     call MixingVector(RK,c,0)
     gl_=-glSi+glGe
     Gibbs_SiGeliq=gl_
  else
     write(*,*)"unexpected d_c,d_T combination",d_c,d_T
     stop
  endif

  end function Gibbs_SiGeliq

double precision function Gibbs_SiGesol(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: gsSi,gsGe,RK(5)
  double precision :: EntropyMixing,T_k
  double precision :: RKs_,ent_,gs_,Tk,Gibbs_T
  double precision :: hsSi(8),hsGe(8)



  Tk=T_k(T)
  if (Tk > 1687.0)then
     hsSi = hsSi2 ;
     hsGe = hsGe3 ;
  elseif (Tk > 1211.4)then
     hsSi = hsSi1 ;
     hsGe = hsGe3 ;
  elseif (Tk > 900.0)then
     hsSi = hsSi1 ;
     hsGe = hsGe2 ;
  else
     hsSi = hsSi1 ;
     hsGe = hsGe1 ;
  endif

  
  if(d_c.eq.0.and.d_T.eq.0)then
     gsSi=gibbs_T(hsSi,Tk,0)
     gsGe=gibbs_T(hsGe,Tk,0)
     call MixingVector(RK,c,0)
     RKs_=         (hrksSiGe0(1)+hrksSiGe0(2)*Tk)*RK(1)
     RKs_=RKs_ + (hrksSiGe1(1)+hrksSiGe1(2)*Tk)*RK(2)
     RKs_=RKs_ + (hrksSiGe2(1)+hrksSiGe2(2)*Tk)*RK(3)
     ent_=g_R*Tk*EntropyMixing(c,0)
     gs_=gsSi*(1d0-c)+gsGe*c+ent_+RKs_
     Gibbs_SiGesol=gs_
  elseif(d_c.eq.0.and.d_T.eq.1)then
     Tk=T_k(T)
     gsSi=gibbs_T(hsSi,Tk,1)
     gsGe=gibbs_T(hsGe,Tk,1)
     call MixingVector(RK,c,0)
     RKs_=         hrksSiGe0(2)*RK(1)
     RKs_=RKs_ +   hrksSiGe1(2)*RK(2)
     RKs_=RKs_ +   hrksSiGe2(2)*RK(3)
     ent_=g_R*EntropyMixing(c,0)
     gs_=gsSi*(1d0-c)+gsGe*c+ent_+RKs_
     Gibbs_SiGesol=gs_
  elseif(d_c.eq.0.and.d_T.eq.2)then
     Tk=T_k(T)
     gsSi=gibbs_T(hsSi,Tk,2)
     gsGe=gibbs_T(hsGe,Tk,2)
     call MixingVector(RK,c,0)
     gs_=gsSi*(1d0-c)+gsGe*c
     Gibbs_SiGesol=gs_
  elseif(d_c.eq.1.and.d_T.eq.0)then
     Tk=T_k(T)
     gsSi=gibbs_T(hsSi,Tk,0)
     gsGe=gibbs_T(hsGe,Tk,0)
     call MixingVector(RK,c,1)
     RKs_=         (hrksSiGe0(1)+hrksSiGe0(2)*Tk)*RK(1)
     RKs_=RKs_ + (hrksSiGe1(1)+hrksSiGe1(2)*Tk)*RK(2)
     RKs_=RKs_ + (hrksSiGe2(1)+hrksSiGe2(2)*Tk)*RK(3)
     ent_=g_R*Tk*EntropyMixing(c,1)
     gs_=-gsSi+gsGe+ent_+RKs_
     Gibbs_SiGesol=gs_
  elseif(d_c.eq.1.and.d_T.eq.1)then
     Tk=T_k(T)
     gsSi=gibbs_T(hsSi,Tk,1)
     gsGe=gibbs_T(hsGe,Tk,1)
     call MixingVector(RK,c,0)
     RKs_=         hrksSiGe0(2)*RK(1)
     RKs_=RKs_ +   hrksSiGe1(2)*RK(2)
     RKs_=RKs_ +   hrksSiGe2(2)*RK(3)
     ent_=g_R*EntropyMixing(c,1)
     gs_=-gsSi+gsGe+ent_+RKs_
     Gibbs_SiGesol=gs_
  elseif(d_c.eq.1.and.d_T.eq.2)then
     Tk=T_k(T)
     gsSi=gibbs_T(hsSi,Tk,2)
     gsGe=gibbs_T(hsGe,Tk,2)
     call MixingVector(RK,c,0)
     gs_=-gsSi+gsGe
     Gibbs_SiGesol=gs_
  else
     write(*,*)"unexpected d_c,d_T combination",d_c,d_T
     stop
  endif

  end function Gibbs_SiGesol




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
