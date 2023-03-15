!!!!!!!freee!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!! functions for the model
!!!!  

!!!!!!!!

!!$
!!$module ExampleFuncs
!!$
!!$contains
!!$
!!$double precision function f1(c,phi)
!!$implicit none
!!$double precision, intent(in):: c,phi
!!$!FreeEnergy(c,T,phi,lp,approx)
!!$double precision :: FreeEnergy
!!$
!!$f1 = FreeEnergy(x,0d0,phi,2,.false.)!/ abs(FreeEnergy(0.25,0d0,0d0,0,.false.) )
!!$
!!$end function f1
!!$
!!$double precision function f2(c,phi)
!!$implicit none
!!$double precision, intent(in):: c,phi
!!$double precision :: FreeEnergy
!!$
!!$f2 = FreeEnergy(c,0d0,phi,2,.false.)!/ abs(FreeEnergy(0.25,0d0,0d0,0,.false.)  )
!!$
!!$
!!$end function f2
!!$
!!$double precision function f3(x)
!!$implicit none
!!$double precision, intent(in):: x
!!$!FreeEnergy(c,T,phi,lp,approx)                                                                                                                                                 
!!$double precision :: FreeEnergy
!!$
!!$f3 = FreeEnergy(x,0d0,0d0,0,.false.)!/ abs(FreeEnergy(0.25,0d0,0d0,0,.false.))
!!$
!!$end function f3
!!$double precision function f4(x)
!!$implicit none
!!$double precision, intent(in):: x
!!$!FreeEnergy(c,T,phi,lp,approx)                                                                                                                                                 
!!$double precision :: FreeEnergy
!!$
!!$f4 = FreeEnergy(x,0d0,1d0,0,.false.)!/ abs(FreeEnergy(0.25,0d0,0d0,0,.false.) )
!!$
!!$end function f4
!!$
!!$
!!$end module ExampleFuncs 


double precision function Cphi_f(Cphi,x,n,d)
implicit none
integer, intent(in)::n,d
double precision, intent(in)::Cphi(-2:100002,6),x
integer :: i
double precision :: w

i=int(x*n)
w=i+1-x*n

Cphi_f = w*Cphi(i,d)+(1.-w)*Cphi(i+1,d)

end function Cphi_f

double precision function Y_func(phi,cs,cL,Ans)
implicit none
double precision, intent(in):: phi,cs,cl,Ans(3)
double precision g_phi,a4
a4 = cL-cS-ans(1)-ans(2)-ans(3)
Y_func = cS + Ans(1)*g_phi(phi,0)+Ans(2)*g_phi(phi,0)**2+Ans(3)*g_phi(phi,0)**3+a4*g_phi(phi,0)**4
end function Y_func

double precision function  RootFinderBiSection(iflag,phi,df)
!use ExampleFuncs 
implicit none
logical, intent (in):: iflag
double precision, intent(in)::df,phi
double precision :: a1,a2,d,a3,u3
double precision :: FreeEnergy
integer :: kk,k

a1=0.
a2=1.
d=1e-11
kk=200
do k =1, kk 
   if(abs(a1-a2) > d)then
      a3=(a1+a2)*0.5
      u3 = FreeEnergy(a3,0d0,phi,2,iflag)-df
      if(u3.ge. 0.)then
         a2=a3
      else
         a1=a3
      endif
   else
      exit
   endif
enddo
RootFinderBiSection = a3
end function RootFinderBiSection

double precision function  RootFinderBiSection2(iflag)
  use solution_parameters
  implicit none
  logical, intent (in):: iflag
  double precision :: a1,a2,d,a3,u3,cL,cS
  double precision :: FreeEnergy,RootFinderBiSection
  integer :: kk,k
  
  a1=0.
  a2=10000.
  d=1e-11
  kk=200
  do k =1, kk 
     if(abs(a1-a2) > d)then
        a3=(a1+a2)*0.5
        g_T0 = a3
        cS=RootFinderBiSection(.false.,0d0,0.)
        cL=RootFinderBiSection(.false.,1d0,0.)
        u3 = (FreeEnergy(cS,0.,0d0,0,.false.)-FreeEnergy(cL,0.,1d0,0,.false.))/(cS-cL)
!        write(*,*)cL,cS,u3,g_t0
        if(u3.le. 0.)then
           a2=a3
        else
           a1=a3
        endif
     else
        exit
     endif
  enddo
  RootFinderBiSection2 = a3
end function RootFinderBiSection2





double precision function GradientEnergyHex(dx,vbles10,ix,iy)
  use paramesh_dimensions 
  use solution_parameters
  use multigrid_parameters  ! for min_dx, total_vars
  implicit none
  double precision, intent (in)::dx,vbles10(N_unk,10,10)
  integer, intent (in):: ix,iy
  double precision :: phi_(2)=0.,phi(2,2)=0.,A_(2),A(2,2),AA(2,2),vbles(N_unk,3,3)
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
!!$p(1,3)=cos(0.8*pi)
!!$p(2,3)=sin(0.8*pi)
  
  
  GradientEnergyHex=0
if(g_linear)then
   do i=1,3
      do j=1,3
         vbles(1,i,j)=vbles10(1,i+ix-2,j+iy-2)
      enddo
   enddo

   phi_(1)=0.5*(vbles(1,3,2)-vbles(1,1,2))/dx
   phi_(2)=0.5*(vbles(1,2,3)-vbles(1,2,1))/dx
   phi(1,1)=(vbles(1,3,2)-2d0*vbles(1,2,2)+vbles(1,1,2))/dx**2
   phi(2,2)=(vbles(1,2,3)-2d0*vbles(1,2,2)+vbles(1,2,1))/dx**2
   phi(1,2) = 0.25*(vbles(1,3,3)+vbles(1,1,1)-vbles(1,3,1)-vbles(1,1,3))/dx**2
   phi(2,1) = phi(1,2)

else  
!!$   phi_=0d0
!!$   phi=0d0
!!$  do i=1,g_order
!!$     phi_(1)  = phi_(1) + vbles10(1,i,iy)*g_diff1(i,ix)
!!$     phi_(2)  = phi_(2) + vbles10(1,ix,i)*g_diff1(i,iy)
!!$     phi(1,1) = phi(1,1) + vbles10(1,i,iy)*g_diff2(i,ix)
!!$     phi(2,2) = phi(2,2) + vbles10(1,ix,i)*g_diff2(i,iy)
!!$  enddo
!!$  phi_(1)  = phi_(1)*(1.+g_y10(g_order-1))/(nxb*dx)
!!$  phi_(2)  = phi_(2)*(1.+g_y10(g_order-1))/(nxb*dx)
!!$  phi(1,1) = phi(1,1)*((1.+g_y10(g_order-1))/(nxb*dx))**2
!!$  phi(2,2) = phi(2,2)*((1.+g_y10(g_order-1))/(nxb*dx))**2
!!$  do i=1,g_order
!!$     do j=1,g_order
!!$        phi(1,2) = phi(1,2) + vbles10(1,i,j)*g_diff1(i,ix)*g_diff1(j,iy)
!!$     enddo
!!$  enddo
!!$  phi(1,2) = phi(1,2)*((1.+g_y10(g_order-1))/(nxb*dx))**2
!!$  phi(2,1)=phi(1,2)
endif  
  

  
  do i=1,2
     xx(i) = phi_(i)
  enddo
  
  
!  x = sqrt(xx(1)**2+xx(2)**2+1d-9)
  
  !   xx=xx/x
  
  
!  d=1d-1*x
  d = 1d-3
  
  
  
  A=0d0
  do i=1,2
     
     dh=0d0
     dh(i)=d
     yy=xx-dh
     A(i,i) = SurfaceEnergy(p,yy,np)
     A(i,i) = A(i,i) -2d0*SurfaceEnergy(p,xx,np)
     yy=xx+dh
     A(i,i) = A(i,i) + SurfaceEnergy(p,yy,np)
     A(i,i) = A(i,i)/d**2
  enddo
  
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
  
  
  GradientEnergyHex = - H
  
  
  
  
  
end function GradientEnergyHex
!!!

double precision function GradientEnergySquare(dx,vbles10,ix,iy)
  use time_dep_parameters
  use paramesh_dimensions 
  use solution_parameters
  use multigrid_parameters  ! for min_dx, total_vars
  use physicaldata
  implicit none
  double precision, intent (in)::dx,vbles10(N_unk,10,10)
  integer, intent (in):: ix,iy
  double precision :: phi_(2)=0.,phi(2,2)=0.,A_(2),A(2,2),AA(2,2),vbles(N_unk,3,3)
  double precision :: yy(2),y,xx(2),x,d=4d-2,h
  double precision :: p(2,6),pi=3.141592654,q(6),dh(2),pi6=.5235987758
  integer :: i,j,np=4,k,l,iSurfaceEnergy
  double precision :: SurfaceEnergy
  double precision, dimension(4,4)::stencil =reshape((/1.,-8.,8.,-1.,-8.,64.,-64.,8.,8.,-64.,64.,-8.,-1.,8.,-8.,1./),(/4,4/))
  double precision, dimension(2,2)::stencil2 =reshape((/1.,-1.,-1.,1./),(/2,2/))
  double precision :: cc(2)=(/-1,1/)
  

  do i=1,np
     p(1,i)=cos(pi*(i-1)/2d0+0.25*pi)
     p(2,i)=sin(pi*(i-1)/2d0+0.25*pi)
  enddo
!!$p(1,3)=cos(0.8*pi)
!!$p(2,3)=sin(0.8*pi)
  
  
  GradientEnergySquare=0
if(g_linear)then
!!$   do i=1,3
!!$      do j=1,3
!!$         vbles(1,i,j)=vbles10(1,i+ix-2,j+iy-2)
!!$      enddo
!!$   enddo

   do i=1,3
      do j=1,3
         vbles(1,i,j)=g_gamma*unk(4,i+ix-2,j+iy-2,1,g_block) + (1-g_gamma)*unk(1,i+ix-2,j+iy-2,1,g_block)
      enddo
   enddo
   
   phi_(1)=0.5*(vbles(1,3,2)-vbles(1,1,2))/dx
   phi_(2)=0.5*(vbles(1,2,3)-vbles(1,2,1))/dx
   phi(1,1)=(vbles(1,3,2)-2d0*vbles(1,2,2)+vbles(1,1,2))/dx**2
   phi(2,2)=(vbles(1,2,3)-2d0*vbles(1,2,2)+vbles(1,2,1))/dx**2
   phi(1,2) = 0.25*(vbles(1,3,3)+vbles(1,1,1)-vbles(1,3,1)-vbles(1,1,3))/dx**2
   phi(2,1) = phi(1,2)

!   phi = phi*(1d0 +4d-2*sin(time*1d6*ix*iy))
!write(*,*)sin(time*1d6),"time random"
else  
endif  
  

  
  do i=1,2
     xx(i) = phi_(i)
  enddo
  
  
!  x = sqrt(xx(1)**2+xx(2)**2+1d-9)
  
  !   xx=xx/x
  
  
!  d=1d-1*x
  d = 1d-3
  
  
  
  A=0d0
  do i=1,2
     
     dh=0d0
     dh(i)=d
     yy=xx-dh
     A(i,i) = SurfaceEnergy(p,yy,np)
     A(i,i) = A(i,i) -2d0*SurfaceEnergy(p,xx,np)
     yy=xx+dh
     A(i,i) = A(i,i) + SurfaceEnergy(p,yy,np)
     A(i,i) = A(i,i)/d**2
  enddo
  
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
  
  
  if(g_axis_symmetric)then
     dh=0d0
     dh(2)=d
     yy = xx+dh
     A(2,2) =  SurfaceEnergy(p,yy,np)
     yy = xx-dh
     A(2,2) = A(2,2) - SurfaceEnergy(p,yy,np)
     A(2,2) = 0.5*A(2,2)/d
     H = H + A(2,2)/g_xy_position(2)
  endif


  GradientEnergySquare = - H
  
  
  
  
  
end function GradientEnergySquare
!!!


double precision function GradientEnergyHexPrismVol2(vbles10,lp,dx,ix,iy)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
implicit none
double precision, intent (in)::dx,vbles10(N_unk,10,10)
integer, intent (in):: lp,ix,iy
double precision :: phi_(2),phi(2,2),A_(4),A(2,2),AA(2,2),vbles(N_unk,3,3)
double precision :: yy(2),y,xx(2),zz(2),ww(2),x,d=4d-2,h
double precision :: p(2,6),pi=3.141592654,q(6),dh(2),pi6=.5235987758
integer :: i,j,np=6,k,l,iSurfaceEnergy
double precision :: SurfaceEnergy,AnisEnergy,UL(2),UR(2),DL(2),DR(2),ML(2),MR(2),MU(2),MD(2)
double precision, dimension(4,4)::stencil =reshape((/1.,-8.,8.,-1.,-8.,64.,-64.,8.,8.,-64.,64.,-8.,-1.,8.,-8.,1./),(/4,4/))
double precision, dimension(2,2)::stencil2 =reshape((/1.,-1.,-1.,1./),(/2,2/))
double precision :: cc(2)=(/-1,1/)
do i=1,6
   p(1,i)=cos(pi*(i-1)/3d0)
   p(2,i)=sin(pi*(i-1)/3d0)
enddo
!!$p(1,3)=cos(pi*0.8)
!!$p(2,3)=sin(0.8*pi)

if(g_linear)then
   do i=1,3
      do j=1,3
         vbles(1,i,j)=vbles10(1,i+ix-2,j+iy-2)
      enddo
   enddo
else
endif


!p= g_phex

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

!qmax=(sqrt(x(1)**2+x(2)**2)+epsilon*maximum(q,n))/(1d0+epsilon)
qmax = maximum(q,n)
AnisEnergy =qmax
end function AnisEnergy




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
!   amaximum = amaximum + ((((s+(1.0-s)*p(i))**8)**8))
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
xx = sqrt(x(1)**2 + x(2)**2)

!xhat = x/xx

q=0d0
do j=1,n
   do i=1,2
      q(j) = q(j)  + p(i,j)*x(i)
   enddo
enddo

qmax=epsilon*xx+maximum(q,n)

!SurfaceEnergy =0.5*qmax**2 + 0.5*epsilon*xx
SurfaceEnergy =0.5*qmax**2
end function SurfaceEnergy


double precision function GradientEnergyIsotropic(dx,vbles10,ix,iy)
  use solution_parameters
  use paramesh_dimensions 
  implicit none
  double precision, intent (in)::dx,vbles10(N_unk,10,10)
  integer, intent (in):: ix,iy
  integer :: i,j
  double precision :: vbles(N_unk,3,3),dphi0=0d0,dphi1=0d0,dphi2=0d0



!if(ix.le.2.or.ix.ge.9.or.iy.le.2.or.iy.ge.9 )then
if(g_linear)then
   do i=1,3
      do j=1,3
         vbles(1,i,j)=vbles10(1,i+ix-2,j+iy-2)
      enddo
   enddo
   dphi2 = vbles(1,1,2)+vbles(1,3,2)+vbles(1,2,1)+vbles(1,2,3)-4.*vbles(1,2,2)
   GradientEnergyIsotropic = -dphi2/dx**2   
else
!!$   GradientEnergyIsotropic = 0d0
!!$   do i=1,g_order
!!$      do j=1,g_order
!!$         GradientEnergyIsotropic = GradientEnergyIsotropic + vbles10(1,i,j)*(g_diff2(i,ix)*g_diff0(j,iy)+g_diff0(i,ix)*g_diff2(j,iy))
!!$      enddo
!!$   enddo
!!$   do i=1,g_order 
!!$      GradientEnergyIsotropic = GradientEnergyIsotropic + vbles10(1,i,iy)*g_diff2(i,ix) + vbles10(1,ix,i)*g_diff2(i,iy)
!!$   enddo
!!$
!!$
!!$   GradientEnergyIsotropic = -GradientEnergyIsotropic*((1.+g_y10(g_order-1))/(nxb*dx))**2
endif

!!$
!!$   do i=1,3
!!$      do j=1,3
!!$         vbles(1,i,j)=vbles10(1,i+ix-2,j+iy-2)
!!$      enddo
!!$   enddo
!!$if(abs(vbles(1,2,2)-0.5).lt.0.1)then
!!$   dphi2 = vbles(1,1,2)+vbles(1,3,2)+vbles(1,2,1)+vbles(1,2,3)-4.*vbles(1,2,2)
!!$   write(*,*)-dphi2/dx**2,GradientEnergyIsotropic
!!$endif


!!$if(abs(vbles10(1,ix,iy)-0.5).lt.0.1.and.ix.eq.5)then
!!$   dphi0 =0d0
!!$   dphi1 =0d0
!!$   dphi2 =0d0
!!$   do i=1,10
!!$      do j=1,10
!!$         dphi0=dphi0+vbles10(1,i,j)*g_diff0(i,ix)*g_diff0(j,iy)
!!$         dphi1=dphi1+vbles10(1,i,j)*g_diff1(i,ix)*g_diff0(j,iy)
!!$         dphi2=dphi2+vbles10(1,i,j)*g_diff2(i,ix)*g_diff0(j,iy)
!!$      enddo
!!$   enddo
!!$   write(*,*)dphi0,vbles(1,2,2),dx
!!$   write(*,*)dphi1/(4*dx),0.5*(vbles(1,3,2)-vbles(1,1,2))/dx,dx
!!$   write(*,*)dphi2/(4*dx)**2,(vbles(1,3,2)-2.*vbles(1,2,2)+vbles(1,1,2))/dx**2,dx
!!$   write(*,*)
!!$endif


!GradientEnergyIsotropic = (vbles(1,1,2)+vbles(1,3,2)+vbles(1,2,1)+vbles(1,2,3)-4*vbles(1,2,2))/dx**2
!GradientEnergyIsotropic = -GradientEnergyIsotropic/(dx)**2



!!$if(abs(dphi2-GradientEnergyIsotropic).gt.0.001.and.abs(vbles(1,2,2)-0.5).lt.0.1.and.ix.ge.5.and.ix.le.6.and. iy.ge.5.and.iy.le.6)then
!!$   write(*,*)dphi2,GradientEnergyIsotropic,ix,iy
!!$   write(*,*)dx
!!$   do i=1,10
!!$      write(*,'("[",9(F12.5,","),F12.5,"]")')vbles10(1,i,:)
!!$   enddo
!!$   write(*,*)
!!$   do i=1,3
!!$      write(*,*)vbles(1,i,:)
!!$   enddo
!!$   stop
!!$endif






end function  GradientEnergyIsotropic



!T_K converts non-dimensional parameter T to dimensioned T_K(T)
double precision function T_K(T)
use solution_parameters
implicit none
double precision, intent (in)::T


 
 T_K = g_T0*(1.+T)

 end function T_K


double precision function SoluteRHS(dx,c_c,a_trap,vbles10,ix,iy)
use solution_parameters
use physicaldata
implicit none
double precision, intent (in):: dx,c_c,a_trap,vbles10(N_unk,10,10)
integer,intent (in) :: ix,iy
integer ii,jj,ip,ih,jv,i,j
double precision c,T,Phi,FreeEnergy,D
double precision :: potential,D_c,f_c,MolPerVol,radial_factor
double precision :: c_dot=0d0,phi_dot=0d0,beta,phi_,PhiSoluteRHS=0d0
double precision :: phidot, abs_grad_phi,u,FE,vbles(N_unk,3,3),Y_phi,gh_phi
double precision :: AlFey1
!

do i=1,3
   do j=1,3
      do ip=1,N_unk
         if(ip.eq.1)then
            vbles(ip,i,j) = (g_gamma*unk(4+(ip-1)*6,i+ix-2,j+iy-2,1,g_block) +  (1-g_gamma)*vbles10(ip,i+ix-2,j+iy-2))
         else
            vbles(ip,i,j) = (g_gamma*unk(4+(ip-1)*6,i+ix-2,j+iy-2,1,g_block) +  (1-g_gamma)*vbles10(ip,i+ix-2,j+iy-2))
         endif
      enddo
   enddo
enddo

!!!
if(New_potential)then
do i=1,3
   do j=1,3
!         vbles(1,i,j) = Y_phi(gh_phi(vbles10(1,i+ix-2,j+iy-2),0),0)
!         vbles(1,i,j) = Y_phi(gh_phi(unk(4,i+ix-2,j+iy-2,1,g_block),0),0)
      phi = g_gamma*unk(4,i+ix-2,j+iy-2,1,g_block) +  (1-g_gamma)*vbles10(1,i+ix-2,j+iy-2)
         vbles(1,i,j) = Y_phi(gh_phi(phi,0),0)

   enddo
enddo
endif




SoluteRHS=0.


!!
!x-coordinate transform values
!! Right
 ih=3
 jv=2
 include 'c_stencil.f90'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
! g_c_force = SoluteRHS/(g_Dch*dx*dx*g_R*g_T0)
 SoluteRHS = SoluteRHS/(g_Dch*dx*dx*g_R*g_T0) 

return

end function SoluteRHS

double precision function SoluteRHSHP(dx,c_c,a_trap,vbles10,ix,iy)
use solution_parameters
  use paramesh_dimensions
implicit none
double precision, intent (in):: dx,c_c,a_trap,vbles10(N_unk,10,10)
integer,intent (in) :: ix,iy
integer ii,jj,ip,ih,jv,i,j
double precision c,T,Phi,FreeEnergy,D
double precision :: potential,D_c,f_c,MolPerVol
double precision :: c_dot=0d0,phi_dot=0d0,beta,phi_,PhiSoluteRHS=0d0
double precision :: phidot, abs_grad_phi,u,FE,vbles(N_unk,3,3)
double precision :: x1,x2,y1,y2
!!$do i=1,3
!!$   do j=1,3
!!$      vbles(:,i,j) = vbles10(:,i+ix-2,j+iy-2)
!!$   enddo
!!$enddo

!!$   GradientEnergyIsotropic = 0d0
!!$   do i=1,10
!!$      do j=1,10
!!$         GradientEnergyIsotropic = GradientEnergyIsotropic + vbles10(1,i,j)*(g_diff2(i,ix)*g_diff0(j,iy)+g_diff0(i,ix)*g_diff2(j,iy))
!!$      enddo
!!$   enddo



SoluteRHSHP=0.

D_c = 0d0
do i=1,g_order
   do j=1,g_order
      SoluteRHSHP = SoluteRHSHP + FreeEnergy(vbles10(2,i,j),0d0,vbles10(1,i,j),2,g_approx)*(g_diff2(i,ix)*g_diff0(j,iy)+g_diff0(i,ix)*g_diff2(j,iy))
      D_c = D_c + (vbles10(1,i,j)*g_D(1)+(1.-vbles10(1,i,j))*g_D(2))*f_c(vbles10(2,i,j),0)*g_diff0(i,ix)*g_diff0(j,iy)
   enddo
enddo
!D_c = (vbles10(1,ix,iy)*g_D(1)+(1.-vbles10(1,ix,iy))*g_D(2))*f_c(vbles10(2,ix,iy),0)
SoluteRHSHP = SoluteRHSHP * D_c

x1=0d0
x2=0d0
y1=0d0
y2=0d0
!!$do i=1,g_order
!!$   do j=1,g_order
!!$      x1 = x1 + (vbles10(1,i,j)*g_D(1)+(1.-vbles10(1,i,j))*g_D(2))*f_c(vbles10(2,i,j),0)*g_diff1(i,ix)*g_diff0(j,iy)
!!$      x2 = x2 + FreeEnergy(vbles10(2,i,j),0d0,vbles10(1,i,j),2,g_approx)                *g_diff1(i,ix)*g_diff0(j,iy)
!!$      y1 = y1 + (vbles10(1,i,j)*g_D(1)+(1.-vbles10(1,i,j))*g_D(2))*f_c(vbles10(2,i,j),0)*g_diff0(i,ix)*g_diff1(j,iy)
!!$      y2 = y2 + FreeEnergy(vbles10(2,i,j),0d0,vbles10(1,i,j),2,g_approx)                *g_diff0(i,ix)*g_diff1(j,iy)
!!$enddo
!!$enddo


do i=1,g_order
      x1 = x1 + (vbles10(1,i,iy)*g_D(1)+(1.-vbles10(1,i,iy))*g_D(2))*f_c(vbles10(2,i,iy),0)*g_diff1(i,ix)
      x2 = x2 + FreeEnergy(vbles10(2,i,iy),0d0,vbles10(1,i,iy),2,g_approx)                *g_diff1(i,ix)
enddo
do j=1,g_order
      y1 = y1 + (vbles10(1,ix,j)*g_D(1)+(1.-vbles10(1,ix,j))*g_D(2))*f_c(vbles10(2,ix,j),0)*g_diff1(j,iy)
      y2 = y2 + FreeEnergy(vbles10(2,ix,j),0d0,vbles10(1,ix,j),2,g_approx)                *g_diff1(j,iy)
enddo

SoluteRHSHP = SoluteRHSHP + x1*x2 + y1*y2

!SoluteRHSHP = SoluteRHSHP/(g_Dch*0.25*nxb*nxb*dx*dx*g_y10(g_order-1)*g_y10(g_order-1)*g_R*g_T0) 
SoluteRHSHP = SoluteRHSHP/(g_Dch*g_R*g_T0*(nxb*dx/(1.+g_y10(g_order-1)))**2)

return

end function SoluteRHSHP

double precision function Y_phi(phi,d)
use solution_parameters
implicit none

double precision, intent (in)::phi
integer, intent (in)::d
double precision :: cL,cS,Z_c,g_phi 

cS = g_AlFe_c(1,1)
cL = g_AlFe_c(2,1)
if(d.eq.0)then
Y_phi = Z_c(cS + phi *(cL - cS),0)
else
Y_phi = Z_c(cS + phi *(cL - cS),1)*(cL - cS)
endif

end function Y_phi



double precision function Z_c(c,d)
use solution_parameters
implicit none

double precision, intent (in)::c
integer, intent(in) :: d
double precision :: FreeEnergy,fL,fS,Y,dfL,dFS 

if(d.eq.0)then
   fL= FreeEnergy(c,0d0,1.,2,.true.)
   fS= FreeEnergy(c,0d0,0.,2,.true.)
   
   Z_c = (g_df- fS)/(fL - fS) 
else !1st derivative
   fL= FreeEnergy(c,0d0,1.,2,.true.)
   fS= FreeEnergy(c,0d0,0.,2,.true.)
   dfL = 0.5*(FreeEnergy(c+1d-8,0d0,1.,2,.true.) - FreeEnergy(c-1d-8,0d0,1.,2,.true.))/1d-8
   dfS = 0.5*(FreeEnergy(c+1d-8,0d0,0.,2,.true.) - FreeEnergy(c-1d-8,0d0,0.,2,.true.))/1d-8



!!$   dfL = FreeEnergy(c,0d0,1.,3,.true.)
!!$   dfS = FreeEnergy(c,0d0,0.,3,.true.)

   Z_c = (-dfS*(fL - fS) - (g_df- fS)*(dfL - dfS))/(fL - fS)**2
end if
end function Z_c

double precision function piecewiseQ(x,xx)
implicit none
double precision, intent (in):: x,xx

if(x.lt.xx)then
   piecewiseQ = -0.1875000000*x**2/xx**2 + 0.3750000000*x/xx
elseif(x.lt.0.5)then
   piecewiseQ = -0.7500000000*x**2/(4.*xx**2 - 4.*X + 1) + 1.500000000*xx*x/(4.*xx**2 - 4.*X + 1) - 0.1875000000*(4.*xx - 1.)/(4.*xx**2 - 4.*X + 1)
endif
end function piecewiseQ

double precision function potential(phi,lp)
use solution_parameters
implicit none

double precision, intent (in)::phi
integer, intent(in)::lp
double precision :: Cphi_f,FreeEnergy,c,dc,gh_phi,cS,cL,Y_phi,phi_,dphi_,height=0.15,left=.1,piecewiseQ
!FreeEnergy(c,T,phi,lp,approx)
! Cphi_f(Cphi,x,n,d)



if(lp.eq.1)then
   if(g_obstacle.le.0d0)then
      potential =  2d0*phi*(1d0-phi)*(1d0-2d0*phi)
!!$   elseif(g_obstacle.lt.0.25)then
!!$      if(phi.lt.g_obstacle)then
!!$         potential = -(3.*phi**2)/(4.*g_obstacle**2) + (3.*phi)/(4.*g_obstacle)
!!$      elseif(phi.lt.1.0-g_obstacle)then
!!$         potential = 0d0
!!$      else
!!$         potential = 3.*(-1. + phi)*(-1. + phi + g_obstacle)/(4.*g_obstacle**2)
!!$      endif
!!$   elseif(g_obstacle.lt.1.0)then
!!$      if(phi.le.0.5)then
!!$         potential = piecewiseQ(phi,left)
!!$      else
!!$         potential = -piecewiseQ(1d0-phi,left)
!!$      endif
   elseif(g_obstacle.lt.1.0)then
      if(phi.lt.left)then
         potential = phi/left*g_obstacle
      elseif(phi.lt.1.0-left)then
         potential = (0.5-phi)/(0.5-left)*g_obstacle
      else
         potential = (phi-1.)/left*g_obstacle
      endif
   else
      write(*,*)"potential: g_obstacle out of bounds",g_obstacle
      stop
   endif

   if(New_potential)then
      phi_ = gh_phi(phi,0)
      dphi_ = gh_phi(phi,1)
 !     c = Cphi_f(g_Cphi,phi,g_Ncphi,2)
!      dc = Cphi_f(g_Cphi,phi,g_Ncphi,3)
      cS = g_AlFe_c(1,1)
      cL = g_AlFe_c(2,1)
      !   potential = potential +(- FreeEnergy(c,0d0,phi,1,.true.)-FreeEnergy(c,0d0,phi,2,.true.)*dc + g_df*(g_AlFe_c(2,1)-g_AlFe_c(1,1))*gh_phi(phi,1))/g_L(1)
      !   potential = potential +(- FreeEnergy(c,0d0,phi,1,.true.) + g_df*(gh_phi(phi,1) - dc))/g_L(1)
!      potential = potential +(- FreeEnergy(c,0d0,phi,1,.true.) -g_df*dc  + g_df*(g_AlFe_c(2,1)-g_AlFe_c(1,1))*gh_phi(phi,1))/g_L(1)
      potential = potential - FreeEnergy(cS+phi_*(cL-cS),0d0,Y_phi(phi_,0),1,.true.)*dphi_*Y_phi(phi_,1)/g_L(1)
   end if
else
   potential = phi**2*(1.-phi)**2
   if(New_potential)then
  !    c = Cphi_f(g_Cphi,phi,g_Ncphi,2)
  !    dc = Cphi_f(g_Cphi,phi,g_Ncphi,3)
      cS = g_AlFe_c(1,1)
      cL = g_AlFe_c(2,1)

      !   potential = potential +(- FreeEnergy(c,0d0,phi,1,.true.)-FreeEnergy(c,0d0,phi,2,.true.)*dc + g_df*(g_AlFe_c(2,1)-g_AlFe_c(1,1))*gh_phi(phi,1))/g_L(1)
      !   potential = potential +(- FreeEnergy(c,0d0,phi,1,.true.) + g_df*(gh_phi(phi,1) - dc))/g_L(1)
!      potential = potential +(- FreeEnergy(c,0d0,phi,0,.true.) + g_AlFe_c(1,3) + g_df*(g_AlFe_c(2,1)-g_AlFe_c(1,1))*gh_phi(phi,0))/g_L(1)
!      potential =potential +(- FreeEnergy(cS+phi*(cL-cS),0d0,Y_phi(phi,0),0,.true.) +( g_AlFe_c(1,3) + phi*(g_AlFe_c(2,3)-g_AlFe_c(1,3))))/g_L(1)



!!!!!!!!
      phi_ = gh_phi(phi,0)

      if(g_quad)then
         potential =potential +(- FreeEnergy(cS+phi_*(cL-cS),0d0,Y_phi(phi_,0),0,.true.) +( g_AlFe_c(1,3) + phi_*(g_AlFe_c(2,3)-g_AlFe_c(1,3))))/g_L(1)
      else
         potential = (- FreeEnergy(cS+phi_*(cL-cS),0d0,Y_phi(phi_,0),0,.true.) +( g_AlFe_c(1,4) + phi_*(g_AlFe_c(2,4)-g_AlFe_c(1,4))))/g_L(1)
      endif

   end if


endif

end function potential


double precision function PhaseRHS(dx,vbles10,ix,iy)
use solution_parameters
use physicaldata
implicit none
double precision, intent(in)::dx,vbles10(N_unk,10,10)
integer, intent(in):: ix,iy
double precision :: vbles(N_unk,3,3)
double precision potential,FreeEnergy,AnisEnergy
double precision :: vertices(2,4),xx(2)
double precision GradientEnergyIsotropic,GradientEnergyHexPrismVol2,GradientEnergyHex,GradientEnergySquare
double precision ::M_tildec,phi,c,T , Y_phi,gh_phi,phi_,dphi_,Sum,phix,phiy,grad_phi,Mobility=0.1
integer :: i,j,ip

Sum=0d0
do i=1,3
   do j=1,3
      do ip=1,N_unk
         vbles(ip,i,j)=vbles10(ip,i+ix-2,j+iy-2)
      enddo
      Sum = Sum + vbles(1,i,j)
   enddo
enddo

phi = g_gamma*unk(4,ix,iy,1,g_block) + (1-g_gamma)*unk(1,ix,iy,1,g_block)!vbles(1,2,2)
c   = g_gamma*unk(4+6,ix,iy,1,g_block)+ (1-g_gamma)*unk(7,ix,iy,1,g_block)!vbles(2,2,2)
T=0d0
if(New_potential)then
   phi_  = gh_phi(phi,0)
   dphi_ = gh_phi(phi,1)

   if(g_Isotropic)then
      phaseRHS=- g_mu_kinetic*(GradientEnergyIsotropic(dx,vbles10,ix,iy)+potential(Phi,1)/g_lambda**2&
           +FreeEnergy(c,T,Y_phi(Phi_,0),1,.true.)*Y_phi(phi_,1)*dphi_/(g_L(1)*g_lambda**2))
   else
      if(g_mu_kinetic.gt.1d0)then
!!$         vertices(1,1)=1d0
!!$         vertices(2,1)=0d0
!!$         vertices(1,2)=-1d0
!!$         vertices(2,2)=0d0

         phix=0.5*(vbles(1,3,2)-vbles(1,1,2))/dx
         phiy=0.5*(vbles(1,2,3)-vbles(1,2,1))/dx
         grad_phi = sqrt(phix**2 + phiy**2+1d-9)

         xx(1) = phix/grad_phi
         xx(2) = phiy/grad_phi
         mobility = abs(xx(1))  + abs(xx(2))/g_mu_kinetic

      endif
!!$      phaseRHS = -g_mu_kinetic*(GradientEnergyHex(dx,vbles10,ix,iy)+potential(Phi,1)/g_lambda**2&
!!$           +FreeEnergy(c,T,Y_phi(Phi_,0),1,.true.)*Y_phi(phi_,1)*dphi_/(g_L(1)*g_lambda**2))
      phaseRHS = -mobility*(GradientEnergySquare(dx,vbles10,ix,iy)+potential(phi,1)/g_lambda**2&
           +FreeEnergy(c,T,Y_phi(Phi_,0),1,.true.)*Y_phi(phi_,1)*dphi_/(g_L(1)*g_lambda**2))

   end if
else
   if(g_Isotropic)then
      phaseRHS=- g_mu_kinetic*(GradientEnergyIsotropic(dx,vbles10,ix,iy)+potential(Phi,1)/g_lambda**2&
           +FreeEnergy(c,T,Phi,1,.true.)/(g_L(1)*g_lambda**2))
   else
      phaseRHS = -g_mu_kinetic*(GradientEnergyHex(dx,vbles10,ix,iy)+potential(Phi,1)/g_lambda**2&
           +FreeEnergy(c,T,Phi,1,.true.)/(g_L(1)*g_lambda**2))
   end if
end if

!phaseRHS=- ( GradientEnergyHexPrismVol2(vbles10,1,dx,ix,iy)  + potential(Phi)/g_lambda**2&
!           +FreeEnergy(c,T,Phi,1,.true.)/(g_L(1)*g_lambda**2))





end function PhaseRHS

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
  double precision rfactor,rf4,rf5,rf6,S
  double precision unk_(M_unk),MolPerVol,vbles10(N_unk,10,10)
  integer ip,n,ii,jj,kk





  do ii=1,g_order
     do jj=1,g_order
        do ip=1,N_unk
           vbles10(ip,ii,jj)=unk(1+(ip-1)*nunkvbles,ii,jj,1,lb)
        enddo
     enddo
  enddo
  


  call get_RHS_Fi(dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,vbles10,i,j,lb)

end subroutine get_RHS_Fi_2d

subroutine get_position(lb,dx,i,j)
use tree
use solution_parameters
implicit none
double precision, intent(in)::dx
integer, intent(in)::lb,i,j
g_xy_position(1) =bnd_box(1,1,lb)+(i-1.5)*dx
g_xy_position(2) =bnd_box(1,2,lb)+(j-1.5)*dx
end subroutine get_position


subroutine get_RHS_Fi(dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,vbles10,ix,iy,lb)
  use solution_parameters
  use physicaldata
  implicit none
  integer, intent(in)::mode_1,nunkvbles,lnblocks,ix,iy,lb
  double precision, intent(in)::dx,dt,dtold
  double precision, intent(inout):: Fi(N_unk)
  double precision, intent(in):: vbles10(N_unk,10,10)
  double precision :: c,T,phi,c_dot,phi_dot,vbles(N_unk,3,3),weight=0.9
  double precision rfactor,rf4,rf5,rf6,S,GradientEnergy,potential
  double precision SoluteRHS,SoluteRHSHP,PhaseRHS,T_k,phi_,FreeEnergy,MolPerVol,a_trap,grad_phi,grad_c,sum
  logical Not_number
  integer ip,n,i,j,k


g_block = lb
  
!put position in g_xy_position
call get_position(lb,dx,ix,iy)
  do i=1,3
     do j=1,3
        do k=1,n_unk
           vbles(k,i,j)=vbles10(k,i+ix-2,j+iy-2)
        enddo
     enddo
  enddo




  if (mode_1.eq.1) then
     g_phi_dot = (unk(1,ix,iy,1,lb)-unk(2,ix,iy,1,lb))/dt
  elseif(mode_1.eq.2)then
     rfactor = dt/dtold
     rf4 = rfactor*rfactor/(rfactor+1.0)
     rf5 = (rfactor+1.0)
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
     g_phi_dot = (rf4*unk(5,ix,iy,1,lb)-rf5*unk(4,ix,iy,1,lb)+rf6*unk(1,ix,iy,1,lb))/dt
  endif




  Fi=0d0

  
  if(mode_1.eq.2)then
     rfactor = dt/dtold
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
    endif
  
  
  
  phi = g_gamma*unk(4,ix,iy,1,lb) + (1-g_gamma)* unk(1,ix,iy,1,lb)!vbles(1,2,2)
  c   = g_gamma*unk(10,ix,iy,1,lb) + (1-g_gamma)* unk(7,ix,iy,1,lb)!vbles(2,2,2)

  
  Fi(1) = PhaseRHS(dx,vbles10,ix,iy)


  if(g_quad.and.(.not. new_potential).and.g_vel.gt.0)then
     a_trap=0d0
     if(g_alpha.ne.0d0)then
        a_trap =g_alpha*sqrt(8d0)*g_lambda*g_ratio*(g_Dch/g_D(1))
     endif
  endif

  if(g_linear)then
     Fi(2) = SoluteRHS(dx,c,a_trap,vbles10,ix,iy)
  else
     Fi(2) = SoluteRHSHP(dx,c,a_trap,vbles10,ix,iy)
  endif




  do ip=1,N_unk
     n=1+(ip-1)*nunkvbles
     if(mode_1.le.1)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)!-unk_(n+1)
     else if(mode_1.eq.2)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)/rf6!-unk_(n+1)
     end if
  enddo




end subroutine get_RHS_Fi



double precision function FreeEnergy(c,T,phi,lp,approx)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi
integer, intent(in)::lp
logical, intent(in)::approx
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(2),hb(2),hr(2)
double precision :: Tk,Cp,R,dc,df,dfdc
double precision :: FE,dh=1d-6
double precision Gibbs_FE_liq,Gibbs_FE_sol,Gibbs_FE_l,Gibbs_FE_e,T_K,g_phi,EntropyMixing,Gibbs_SiGeliq,Gibbs_SiGesol !functions
double precision :: Gibbs_FE_l_AlSi,Gibbs_FE_e_AlSi,cc,Gibbs_Fe_l_PbSn,Gibbs_Fe_e_PbSn,Gibbs_Fe_l_NiCu,Gibbs_Fe_e_NiCu
double precision :: Gibbs_Fe_La,Gibbs_Fe_a,AlFeFreeEnergy,Cphi_f,U_CurveS,U_CurveL
integer :: i,j


!(g_AlFe_LUT,xx0,g_Ncphi,2)
if(approx)then
   if(g_LUT)then
      if(lp.eq.0)then
         FreeEnergy = Cphi_f(g_AlFe_LUT,c,g_Ncphi,4)*g_phi(phi,0) +  Cphi_f(g_AlFe_LUT,c,g_Ncphi,1)*(1d0 -g_phi(phi,0))
      elseif(lp.eq.1)then
         FreeEnergy = Cphi_f(g_AlFe_LUT,c,g_Ncphi,4)*g_phi(phi,1) -  Cphi_f(g_AlFe_LUT,c,g_Ncphi,1)*g_phi(phi,1)
      elseif(lp.eq.2)then
         FreeEnergy =  Cphi_f(g_AlFe_LUT,c,g_Ncphi,5)*g_phi(phi,0) +  Cphi_f(g_AlFe_LUT,c,g_Ncphi,2)*(1d0 -g_phi(phi,0))
      elseif(lp.eq.3)then
         FreeEnergy =  Cphi_f(g_AlFe_LUT,c,g_Ncphi,6)*g_phi(phi,0) +  Cphi_f(g_AlFe_LUT,c,g_Ncphi,3)*(1d0 -g_phi(phi,0))
      else
         write(*,*)"FreeEnergy"
         stop
      endif
   else
      if(g_Quad)then
         if(lp.eq.0)then
            FreeEnergy =   (0.5*g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))**2+g_AlFe_c(2,2)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,3) )*g_phi(phi,0)&
                 +  (0.5*g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))**2+g_AlFe_c(1,2)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,3) )*g_phi(1-phi,0)
            
         elseif(lp.eq.1)then
            FreeEnergy =   (0.5*g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))**2+g_AlFe_c(2,2)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,3) )*g_phi(phi,1)&
                 -  (0.5*g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))**2+g_AlFe_c(1,2)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,3) )*g_phi(phi,1)
         elseif(lp.eq.2)then
            FreeEnergy =   (g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,2) )*g_phi(phi,0)&
                 +  (g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,2) )*g_phi(1-phi,0)
         elseif(lp.eq.3)then
            FreeEnergy =   g_AlFe_c(2,4)*g_phi(phi,0) + g_AlFe_c(1,4)*g_phi(1-phi,0)
         elseif(lp.eq.-2)then !f_cc
            FreeEnergy = (1 - g_phi(phi,0))*g_AlFe_c(1,4) + g_phi(phi,0)*g_AlFe_c(2,4)
            
         elseif(lp.eq.-1)then !f_c_phi
            FreeEnergy = (g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,2) )*g_phi(phi,1)&
                 -  (g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,2) )*g_phi(phi,1)
            
         else
            write(*,*)"lp = ",lp
            stop
         endif
      else !quartic
         if(lp.eq.0)then
            FreeEnergy = U_CurveS(c,0)*(1-g_phi(phi,0))+U_CurveL(c,0)*g_phi(phi,0)
         elseif(lp.eq.1)then
            FreeEnergy = - U_CurveS(c,0)*g_phi(phi,1)+U_CurveL(c,0)*g_phi(phi,1)
         elseif(lp.eq.2)then
            FreeEnergy = U_CurveS(c,1)*(1-g_phi(phi,0))+U_CurveL(c,1)*g_phi(phi,0)
         elseif(lp.eq.3)then
            FreeEnergy = U_CurveS(c,2)*(1-g_phi(phi,0))+U_CurveL(c,2)*g_phi(phi,0)
         else
            write(*,*)"FreeEnergy in mp_stencil_terms.f90, and lp =",lp
            stop
         endif
      endif !end Quadratic or quartic
   endif ! end LUT of polynomials
else
!thermodynamic data base   
   FreeEnergy =    AlFeFreeEnergy(c,T,phi,lp)

endif !end approx or thermo base
return

end function FreeEnergy

double precision function AlSiFreeEnergy(c,T,phi,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(2),hb(2),hr(2)
double precision :: Tk,Cp,R,S
double precision :: FE,dh=1d-6
double precision Gibbs_FE_l_AlSi,Gibbs_FE_e_AlSi,g_phi,dAlSiFreeEnergy
integer :: i,j
logical :: g_Tflag = .false.


AlSiFreeEnergy=0d0


if(lp.eq.0)then !free energy proper J/mol
!Gibbs_FE_e_AlSi
   AlSiFreeEnergy =                    Gibbs_FE_l_AlSi(c,T,0,0)*g_phi(phi,0)
   AlSiFreeEnergy = AlSiFreeEnergy +   Gibbs_FE_e_AlSi(c,T,0,0)*(1-g_phi(phi,0))
   return
elseif(lp.le.1)then
   if(g_Tflag)then
      write(*,*)g_Tflag
      stop
   else
      AlSiFreeEnergy =                    Gibbs_FE_l_AlSi(c,T,0,0)*g_phi(phi,1)
      AlSiFreeEnergy = AlSiFreeEnergy -   Gibbs_FE_e_AlSi(c,T,0,0)*(g_phi(phi,1))
      return
   endif
elseif(lp.eq.2)then !Solute


   AlSiFreeEnergy =                    Gibbs_FE_l_AlSi(c,T,1,0)*g_phi(phi,0)
   AlSiFreeEnergy = AlSiFreeEnergy +   Gibbs_FE_e_AlSi(c,T,1,0)*(1-g_phi(phi,0))

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

end function AlSiFreeEnergy

double precision function AlFeFreeEnergy(c,T,phi,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi
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

   AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,0.,0,0)*g_phi(phi,0)
   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e2_AlFe(c,0.,0,0)*(1-g_phi(phi,0))
   return
elseif(lp.le.1)then
   if(g_Tflag)then
      write(*,*)g_Tflag
      stop
   else
      AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,0.,0,0)*g_phi(phi,1)
      AlFeFreeEnergy = AlFeFreeEnergy -   Gibbs_FE_e2_AlFe(c,0.,0,0)*(g_phi(phi,1))
      return
   endif
elseif(lp.eq.2)then !Solute


   AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,0.,1,0)*g_phi(phi,0)
   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e2_AlFe(c,0.,1,0)*(1-g_phi(phi,0))

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
  use solution_parameters
  integer, intent(in)::d
  double precision, intent(in) :: c
  !double precision :: cL=0.235,cR=0.2724637681 
  double precision :: cL,cR 
  
  cL = g_left_limit
  cR = g_right_limit
  
  if(d.eq.0)then
     !      AlFey1=(0.2350/c - 0.8625)/0.1375 
     
     AlFey1 = (c-cR)/(cL-cR)
     
     AlFey1 = min(0.9999,max(0.0001,Alfey1))
     return
  elseif(d.eq.1)then
!!$   if(c.gt.cR.or.c.lt.cL)then
!!$      write(*,*)"c out of bound in AlFey1()"
!!$      stop
!!$   endif
     
     !      AlFey1=-(0.2350/c**2)/0.1375 
     AlFey1 = 1/(cL-cR)
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
!c=g_right_limit37681 - 0.0374637681*cc


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
   
   if(c.lt.g_left_limit)then
      y1 = AlFey1(g_left_limit,0)
      entSLb=0.1375*g_R*Tk*EntropyMixing(y1,0)
      gb1=SLb1*y1 + SLb2*(1d0-y1) 
      gb=(gb1 + entSLb)/(0.8625 + 0.1375*y1)
      Gibbs_FE_e2_AlFe=gb + (c-g_left_limit)*g_left_tan + 0.5*(c-g_left_limit)**2*g_left_acc 
   elseif(c.gt.g_right_limit)then
      y1 = AlFey1(g_right_limit,0)
      entSLb=0.1375*g_R*Tk*EntropyMixing(y1,0)
      gb1=SLb1*y1 + SLb2*(1d0-y1) 
      gb=(gb1 + entSLb)/(0.8625 + 0.1375*y1)
      Gibbs_FE_e2_AlFe=gb + (c-g_right_limit)*g_right_tan +0.5*(c-g_right_limit)**2*g_right_acc+1./6.*(c-g_right_limit)**3*g_right_3
   else
      y1 = AlFey1(c,0)
      entSLb=0.1375*g_R*Tk*EntropyMixing(y1,0)
      gb1=SLb1*y1 + SLb2*(1d0-y1) 
      gb=(gb1 + entSLb)/(0.8625 + 0.1375*y1)
      Gibbs_FE_e2_AlFe=gb
   endif
   return


elseif(d_c.eq.0.and.d_T.eq.1)then
   sLb1=gibbs_T(alfeh1,Tk,1)
   sLb2=gibbs_T(alfeh2,Tk,1)
   y1 = AlFey1(max(min(g_right_limit,c),g_left_limit),0)
   entSLb=0.1375*g_R*EntropyMixing(y1,0)
   gb1=SLb1*y1 + SLb2*(1d0-y1) 
   gb=(gb1 + entSLb)/(0.8625 + 0.1375*y1)
   return
elseif(d_c.eq.0.and.d_T.eq.2)then
elseif(d_c.eq.1.and.d_T.eq.0)then
   sLb1=gibbs_T(alfeh1,Tk,0)
   sLb2=gibbs_T(alfeh2,Tk,0)
!   dsLb1=0d0
!   dsLb2=0d0
   y1 = AlFey1(max(min(g_right_limit,c),g_left_limit),0)
   dy1 = AlFey1(max(min(g_right_limit,c),g_left_limit),1)
   entSLb = 0.1375*g_R*Tk*EntropyMixing(y1,0)
   dentSLb = 0.1375*g_R*Tk*EntropyMixing(y1,1)
   gb1 =SLb1*y1 + SLb2*(1d0-y1) 
   dgb1 = (SLb1 - SLb2)
   dgb = -0.1375*(gb1 + entSLb)/(0.8625 + 0.1375*y1)**2 + (dgb1 + dentSLb)/(0.8625 + 0.1375*y1)
   Gibbs_FE_e2_AlFe = dgb*dy1
   if(c.lt.g_left_limit)then
      Gibbs_FE_e2_AlFe =Gibbs_FE_e2_AlFe  + (c-g_left_limit)*g_left_acc 
   elseif(c.gt.g_right_limit)then
      Gibbs_FE_e2_AlFe = Gibbs_FE_e2_AlFe + (c-g_right_limit)*g_right_acc + 0.5*(c-g_right_limit)**2*g_right_3
   endif



   


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
   EntropyMixing = log(c)-log(1.-c)
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

double precision function gg_phi(phi,d)
use solution_parameters
  double precision, intent (in):: phi
  integer, intent (in):: d
  double precision :: n,h_phi,u_phi,u,h,g_phi
  gg_phi = g_phi(phi,d)
  return
  gg_phi=0d0
  if(phi.lt.1d-16)return
  if(g_skew_potential)then
     if(d.eq.1)then
        gg_phi =  ((-1/(2-phi)/2-1/phi/2)*(phi**2-phi**3/3)+(log(2-phi)/2-log(phi)/2)*(-phi**2+2*phi)-1.D0/3.D0+phi/3+2.D0/3.D0/(2-phi))&
             /(2.D0/3.D0*log(2.D0)-1.D0/6.D0)
     elseif(d.eq.2)then
        gg_phi = ((-1/((2-phi)**2)/2+1/phi**2/2)*(phi**2-phi**3/3)+2*(-1/(2-phi)/2-1/phi/2)*(-phi**2+2*phi)&
             +(log(2-phi)/2-log(phi)/2)*(-2*phi+2)+1.D0/3.D0+2.D0/3.D0/(2-phi)**2)/(2.D0/3.D0*log(2.D0)-1.D0/6.D0)
     else !assume 0
        gg_phi = ((log(2-phi)/2-log(phi)/2)*(phi**2-phi**3/3)-phi/3&
             +phi**2/6-2.D0/3.D0*log(2-phi)+2.D0/3.D0*log(2.D0))/(2.D0/3.D0*log(2.D0)-1.D0/6.D0)
     endif
  elseif(g_power_skew)then
     n = g_power 
     if(d.eq.1)then
        gg_phi = 0.25*phi*(n+4)*(n+2)/(1-phi)*(1-phi)**(0.5*n+1) 
     elseif(d.eq.2)then
        gg_phi = -0.125*(n*phi+2*phi-2)*(n+4)*(n+2)/(1-phi)**2*(1-phi)**(0.5*n+1)
     else !assume 0
        gg_phi = 1-0.5*((2+n)*phi+2)*(1-phi)**(0.5*n+1)
     endif
!
  elseif(g_s_skew)then
     if(d.eq.0)then

!AMM
        if(g_AMM.eq.1)then
           gg_phi = -0.1485934166D1*phi**2+0.3268185149D1*phi-0.3258617531D0*log(2544018011.D0+0.2551487495D11*phi)+0.7057191429D1
        elseif(g_AMM.eq.2)then
           gg_phi = -0.2D1*phi**3+0.3D1*phi**2
        else
           write(*,*)"g_AMM is 1 or 2"
        endif


     elseif(d.eq.1)then 
!AMM

        if(g_AMM.eq.1)then
           gg_phi = -0.1516536976D3*phi*(-1.D0+phi)/(0.5088036022D1+0.510297499D2*phi)
        elseif(g_AMM.eq.2)then
           gg_phi = -0.6D1*phi*(-1.D0+phi)
        else
           write(*,*)"g_AMM is 1 or 2"
        endif



     endif

  else
     if(d.eq.1)then
        gg_phi = 6.*phi*(1.-phi)
!        gg_phi = 30.*phi**4-60.*phi**3+30.*phi**2
     elseif(d.eq.2)then
                gg_phi = 6.-12.*phi
!    gg_phi = 120.*phi**3-180.*phi**2+60.*phi
 else !assume 0
            gg_phi = 3.*phi**2-2.*phi**3
!    gg_phi = 6.*phi**5-15.*phi**4+10.*phi**3
     endif
  endif
end function gg_phi


double precision function gh_phi(phi,d)
use solution_parameters
  double precision, intent (in):: phi
  integer, intent (in):: d
  gh_phi=0d0
  if(d.eq.1)then
             gh_phi = 6.*phi*(1.-phi)
!     gh_phi = 30.*phi**4-60.*phi**3+30.*phi**2
  elseif(d.eq.2)then
                    gh_phi = 6.-12.*phi
 !    gh_phi = 120.*phi**3-180.*phi**2+60.*phi
  else !assume 0
               gh_phi = 3.*phi**2-2.*phi**3
!     gh_phi = 6.*phi**5-15.*phi**4+10.*phi**3
  endif

end function gh_phi


double precision function g_phi(phi,d)
use solution_parameters
  double precision, intent (in):: phi
  integer, intent (in):: d
 double precision :: gh_phi
!!$ if(d.eq.1)then
!!$    if(phi.lt.1d-7.or.phi.gt.1.-1d-7)then
!!$       g_phi=0d0
!!$       return
!!$    endif
!!$ endif
!!$ g_phi = gh_phi(phi,d)
!!$ return
  if(New_potential)then
     
     g_phi=0d0
     if(d.eq.1)then
        g_phi = 1d0
     !   !             g_phi = 6.*phi*(1.-phi)
        !g_phi = 30.*phi**4-60.*phi**3+30.*phi**2
     elseif(d.eq.2)then
        g_phi = 0d0
        !     g_phi = 6.-12.*phi
        !g_phi = 120.*phi**3-180.*phi**2+60.*phi
     else !assume 0
        g_phi = phi
        ! g_phi = 3.*phi**2-2.*phi**3
        !g_phi = 6.*phi**5-15.*phi**4+10.*phi**3
     endif
     
  else
     g_phi=gh_phi(phi,d)
  end if
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

double precision function Gibbs_FE_l_AlSi(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5), gs1A,gs1B, gs2A,gs2B
  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gl_,Tk,Gibbs_T
  double precision :: hlA_(8),hlB_(8),hs1A_(8),hs1B_(8),hs2A_(8),hs2B_(8)
  double precision :: Bezier_2,Quad
  integer :: T_low,T_high
  double precision :: T_weight
  logical :: approx=.false.

! AlSi below



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

end function Gibbs_FE_l_AlSi


double precision function Gibbs_FE_e_AlSi(c,T,d_c,d_T)
  use solution_parameters
  use time_dep_parameters
  implicit none
  double precision, intent(in) :: c,T
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
