!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!! functions for the model
!!!!  

!!!!!!!!
subroutine set_mollifier
  use solution_parameters
implicit none

double precision :: area,phi,d,phi1,phi2,q=20.
integer :: i
!
d = g_mol_d
g_mol_epsilon = 2.*d/g_mol_n
do i = 1, g_mol_n-1
   g_mollifier(i,1) = -d+g_mol_epsilon*i
enddo
!!$
!!$area = 0.
!!$do i = 1, g_mol_n-1
!!$   phi1=exp(-q*(g_mollifier(i,1)/d+0.5))
!!$   phi2=exp(q*(g_mollifier(i,1)/d-0.5))
!!$   phi = phi1 + phi2
!!$   area = area + 1/(1+phi)
!!$enddo
!!$do i = 1, g_mol_n-1
!!$   phi1=exp(-q*(g_mollifier(i,1)/d+0.5))
!!$   phi2=exp(q*(g_mollifier(i,1)/d-0.5))
!!$   phi = phi1 + phi2
!!$
!!$   g_mollifier(i,2) = 1./(1.+phi)/area
!!$enddo
!!$area=0
!!$do i = 1, g_mol_n-1
!!$   phi1=exp(-q*(g_mollifier(i,1)/d+0.5))
!!$   phi2=exp(q*(g_mollifier(i,1)/d-0.5))
!!$   phi = phi1 + phi2
!!$   g_mollifier(i,3) = -1./(1.+phi)**2*(-q*phi1+q*phi2)/d
!!$   area = area - g_mollifier(i,3)*g_mollifier(i,1)
!!$enddo
!!$do i = 1, g_mol_n-1
!!$   g_mollifier(i,3) = g_mollifier(i,3)/area
!!$enddo
area = 0.
do i = 1, g_mol_n-1 
   phi = (g_mollifier(i,1)/d-1.)**2* (g_mollifier(i,1)/d+1.)**2
   area = area + phi
enddo
do i = 1,  g_mol_n-1
   g_mollifier(i,2) = (g_mollifier(i,1)/d-1.)**2* (g_mollifier(i,1)/d+1.)**2/area 
enddo
area = 0.
do i = 1,  g_mol_n-1
   g_mollifier(i,3) =-4./(d**4)*(g_mollifier(i,1)/d-1.)* (g_mollifier(i,1)/d+1.)*g_mollifier(i,1)
   area = area - g_mollifier(i,3)*g_mollifier(i,1)
enddo

do i = 1,  g_mol_n-1
   g_mollifier(i,3) =g_mollifier(i,3)/area 
enddo



end subroutine set_mollifier

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
!!!!!!!

double precision function potential3D(phi,lp,n)
use solution_parameters
use multigrid_parameters
use physicaldata
implicit none

double precision, intent (in)::phi,n(3)
integer, intent(in)::lp
double precision :: Cphi_f,FreeEnergy,c,dc,gh_phi,cS,cL,Y_phi,phi_,dphi_,p=0.01,left=0.1,h
integer :: i
!FreeEnergy(c,T,phi,lp,approx)
! Cphi_f(Cphi,x,n,d)

if(lp.eq.1)then
   if(g_mesh_dependence)then
      
      h = g_lambda/min_dx
      potential3D = -6.*phi
      do i=1,3
         potential3D = potential3D + phi/(phi+(1.-phi)*exp(n(i)/h))+ phi/(phi+(1.-phi)*exp(-n(i)/h))
      enddo
!!$      potential3D = (-1 + phi)*(2*phi - 1)*phi + 1/12*(-1 + phi)*(24*phi**3 - 36*phi**2 + 14*phi - 1)*phi*(n(1)**4 + n(2)**4 + n(3)**4)/h**2
      potential3D = potential3D*h*h
      return
   endif

   if(g_obstacle.le.0d0)then
      potential3D =  phi*(1d0-phi)*(1d0-2d0*phi)
      if(g_pot_h)then
         h = g_lambda/min_dx
         potential3D = potential3D*(2*exp(h) - exp(2*h) - 1)/(((-1 + phi)*exp(h) - phi)*h**2*(phi*exp(h) - phi + 1))
      endif

      return
   elseif(g_obstacle.lt.1.0)then
      if(phi.lt.left)then
         potential3D = phi/left*g_obstacle
      elseif(phi.lt.1.0-left)then
         potential3D = (0.5-phi)/(0.5-left)*g_obstacle
      else
         potential3D = (phi-1.)/left*g_obstacle
      endif
      
   else
      write(*,*)"potential: g_obstacle out of bounds",g_obstacle
      stop
   endif
   return
   potential3D = potential3D*g_barrier


!!$
!!$   if(phi.lt.p)then
!!$      potential3D = 0.2*phi/p
!!$   elseif(phi.lt.1.0-p)then
!!$      potential3D = 0.2*(0.5-phi)/(0.5-p)
!!$   else
!!$      potential3D = (phi-1.)*0.2/p
!!$   endif
   if(New_potential)then
      phi_ = gh_phi(phi,0)
      dphi_ = gh_phi(phi,1)
      cS = g_AlFe_c(1,1)
      cL = g_AlFe_c(2,1)
      potential3D = potential3D - FreeEnergy(cS+phi_*(cL-cS),0d0,Y_phi(phi_,0),1,.true.)*dphi_*Y_phi(phi_,1)/g_L(1)
   end if
else
   potential3D = phi**2*(1.-phi)**2
   if(New_potential)then
      cS = g_AlFe_c(1,1)
      cL = g_AlFe_c(2,1)
      phi_ = gh_phi(phi,0)
      potential3D =potential3D +(- FreeEnergy(cS+phi_*(cL-cS),0d0,Y_phi(phi_,0),0,.true.) +( g_AlFe_c(1,3) + phi_*(g_AlFe_c(2,3)-g_AlFe_c(1,3))))/g_L(1)
   end if
endif

end function potential3D


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

  g_phi=0d0
  if(d.eq.1)then
             g_phi = 6.*phi*(1.-phi)
  elseif(d.eq.2)then
                    g_phi = 6.-12.*phi
  else !assume 0
               g_phi = 3.*phi**2-2.*phi**3
  endif

!!$  if(New_potential)then
!!$     
!!$     g_phi=0d0
!!$     if(d.eq.1)then
!!$        g_phi = 1d0
!!$     elseif(d.eq.2)then
!!$        g_phi = 0d0
!!$     else !assume 0
!!$        g_phi = phi
!!$     endif
!!$  else
!!$     g_phi=gh_phi(phi,d)
!!$  end if
end function g_phi


double precision function SurfaceEnergym(p,x,n,der)
use solution_parameters
implicit none
double precision, intent (in):: p(3,100),x(3)
integer, intent (in) :: n,der
integer :: i,j,k,m=128
double precision :: amaximum,q(100),qmax,xx,tmp,S0,p3(100,3,3),radius,p1(100),p2(100,3)
double precision :: ID(3,3)=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/)),d



d = g_d_facets





do i = 1,n
   do j=1,3
      do k=1,3
         p3(i,j,k) = p(j,i)*p(k,i) + id(j,k)*d!*exp(-g_time/g_time_end)
      enddo
   enddo
enddo


p1 = 0.
p2 = 0.
do i = 1,n
   do j=1,3
      do k=1,3
         p1(i) = p1(i) + 0.5*x(j)*p3(i,j,k)*x(k)
         p2(i,k) = p2(i,k) + x(j)*p3(i,j,k)
      enddo
   enddo
enddo
S0 = 0.

do i = 1,n
   if(p1(i).gt.0.)then
!      S0 = S0 + p1(i)**m

      S0 = S0 + exp(m*log(p1(i)))
   endif
enddo

if(der.eq.0)then
!   SurfaceEnergym = S0**(1./m)
   SurfaceEnergym = exp(1./m*log(S0))

else
   SurfaceEnergym = 0.
   do i = 1,n
      if(p1(i).gt.0.)then
!         SurfaceEnergym = SurfaceEnergym  +S0**(1./m-1.)*p1(i)**(m-1)*p2(i,der)     
        SurfaceEnergym = SurfaceEnergym  + exp((1./m-1.)*log(S0) + (m-1.)*log(p1(i)))*p2(i,der)
      endif
   enddo
endif


return

end function SurfaceEnergym

double precision function SurfaceEnergy0(p,x,n)
use solution_parameters
implicit none
double precision, intent (in):: p(3,100),x(3)
integer, intent (in) :: n
integer :: i,j,k

double precision :: ID(3,3)=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/)),p3(100,3,3),d,p1(100),S0

d = g_d_facets
do i = 1,n
   do j=1,3
      do k=1,3
         p3(i,j,k) = (p(j,i)*p(k,i) + id(j,k)*d)/(1+d)
      enddo
   enddo
enddo


p1 = 0.

do i = 1,n
   do j=1,3
      do k=1,3
         p1(i) = p1(i) + 0.5*x(j)*p3(i,j,k)*x(k)
      enddo
   enddo
enddo

S0 = 0.
do i = 1,n
   if(p1(i).gt.S0)S0 = p1(i)
enddo
SurfaceEnergy0 = S0

end function SurfaceEnergy0


double precision function SurfaceEnergy(p,x,n,der)
use solution_parameters
implicit none
double precision, intent (in):: p(3,100),x(3)
integer, intent (in) :: n,der
integer :: i,j,k,n_mol
double precision :: y(3),t(3),u(3),v(3),SurfaceEnergy0

SurfaceEnergy = 0.
k=der
y = x
do n_mol = 1,g_mol_n - 1
   y(k) = x(k) - g_mollifier(n_mol,1)
   SurfaceEnergy = SurfaceEnergy  + g_mollifier(n_mol,3)*SurfaceEnergy0(g_HexCol,y,12)
enddo



if(.false.)then
if(der.eq.1)then
   do i=1,g_mol_n -1
      do j=1,g_mol_n-1
         do k = 1,g_mol_n - 1
            y(1) = x(1) - g_mollifier(i,1)
            y(2) = x(2) - g_mollifier(j,1)
            y(3) = x(3) - g_mollifier(k,1)
            SurfaceEnergy = SurfaceEnergy  + g_mollifier(i,3)*g_mollifier(j,2)*g_mollifier(k,2)*SurfaceEnergy0(g_HexCol,y,12)
         enddo
      enddo
   enddo
elseif(der.eq.2)then
   do i=1,g_mol_n -1
      do j=1,g_mol_n-1
         do k = 1,g_mol_n - 1
            y(1) = x(1) - g_mollifier(i,1)
            y(2) = x(2) - g_mollifier(j,1)
            y(3) = x(3) - g_mollifier(k,1)
            SurfaceEnergy = SurfaceEnergy  + g_mollifier(i,2)*g_mollifier(j,3)*g_mollifier(k,2)*SurfaceEnergy0(g_HexCol,y,12)
         enddo
      enddo
   enddo
else
   do i=1,g_mol_n -1
      do j=1,g_mol_n-1
         do k = 1,g_mol_n - 1
            y(1) = x(1) - g_mollifier(i,1)
            y(2) = x(2) - g_mollifier(j,1)
            y(3) = x(3) - g_mollifier(k,1)
            SurfaceEnergy = SurfaceEnergy  + g_mollifier(i,2)*g_mollifier(j,2)*g_mollifier(k,3)*SurfaceEnergy0(g_HexCol,y,12)
         enddo
      enddo
   enddo
endif
endif
end function SurfaceEnergy



double precision function AnisPowEnergy(p,x,n,der)
use solution_parameters
implicit none
double precision, intent (in):: p(3,100),x(3)
integer, intent (in) :: n,der
integer :: i,j,k,m=32
double precision :: amaximum,q(100),qmax,xx,tmp,S0,p3(100,3,3),radius,p1(100),p2(100,3),a0
double precision :: ID(3,3)=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/)),d



d = g_d_facets

p1=0.
do i = 1,n
   p1(i)=p1(i)+p(j,i)*x(j)
enddo
S0 = 0.
do i = 1,n
   if(p1(i).gt.0.)then
      S0 = S0 + exp(m*log(p1(i)))
   endif
enddo

if(der.eq.0)then

   AnisPowEnergy = exp(1./m*log(S0))

else

   AnisPowEnergy = 0.

   do i = 1,n
      if(p1(i).gt.0.)then
!!$        AnisPowEnergy = AnisPowEnergy  + exp((1./m-1.)*log(S0) + (m-1.)*log(p1(i)))*p(der,i)
!        AnisPowEnergy = AnisPowEnergy  + S0**(1./m-1.)*p1(i)**(m-1.)*p(der,i)
!!$
        AnisPowEnergy = AnisPowEnergy  + exp((1./m-1.)*log(S0)+(m-1.)*log(p1(i)))*p(der,i)

      endif
   enddo
endif








return

end function AnisPowEnergy


double precision function GradientEnergyHexColumnx(dx,vbles10,ix,iy,iz)!(vbles,lp,dx)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
use physicaldata
implicit none
double precision, intent (in)::dx,vbles10(N_unk,10,10,10)
integer :: ix,iy,iz
double precision :: phi_(3),phi(3,3),phi_tmp(3,3),A_(3),A(3,3),AA(3,3),fgf
double precision :: yy(3),y,xx(3),x,px,LapS,H,Lap,d=4d-2,c(6)=(/1,2,3,1,2,3/),rr(3)=(/1,1,2/)
double precision :: p(3,20),pi=3.141592654,q(12),dh(3),pi6=.5235987758,t1(3),t2(3),t3(3)
integer :: i,j,np=12,k,l,iSurfaceEnergy
double precision :: SurfaceEnergy,maxA,psi=1.618033988749895,vbles(N_unk,3,3,3)
double precision :: ID(3,3)=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))
!!$
!!$do i=1,3
!!$   do j=1,3
!!$      do k=1,3
!!$         vbles(1,i,j,k)=vbles10(1,i+ix-2,j+iy-2,k+iz-2)
!!$      enddo
!!$   enddo
!!$enddo


   do i=1,3
      do j=1,3
         do k=1,3
            vbles(1,i,j,k)=g_gamma*unk(4,i+ix-2,j+iy-2,k+iz-2,g_block) + (1-g_gamma)*unk(1,i+ix-2,j+iy-2,k+iz-2,g_block)
         enddo
      enddo
   enddo
   


GradientEnergyHexColumnx=0
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


   x = sqrt(xx(1)**2+xx(2)**2+xx(3)**2+1d-6)

   xx=xx/x


  d=g_d_facets
  np=12

   do i=1,3
      dh=0d0
      dh(i)=d
      yy=xx-dh
      A(i,i) = SurfaceEnergy(g_Hexcol,yy,np)
      A(i,i) = A(i,i) -2d0*SurfaceEnergy(g_Hexcol,xx,np)
      yy=xx+dh
      A(i,i) = A(i,i) + SurfaceEnergy(g_Hexcol,yy,np)
      A(i,i) = A(i,i)/d**2
   enddo

   do i=1,2
      do j=i+1,3

         dh=0d0
         dh(i)=d
         dh(j)=d
         yy=xx + dh
         A(i,j) = SurfaceEnergy(g_Hexcol,yy,np)

         yy=xx - dh
         A(i,j) = A(i,j) + SurfaceEnergy(g_Hexcol,yy,np)

         dh=0d0
         dh(i)=d
         dh(j)=-d
         yy=xx+dh
         A(i,j) = A(i,j) - SurfaceEnergy(g_Hexcol,yy,np)

         yy=xx-dh
         A(i,j) = A(i,j) - SurfaceEnergy(g_Hexcol,yy,np)

         A(i,j)=A(i,j)*0.25/d**2
         A(j,i)=A(i,j)
      enddo
   enddo



   


   do i=1,3
      do j=1,3
         H = H + phi(i,j)*A(i,j) 
      enddo
   enddo


   GradientEnergyHexColumnx = - H





end function GradientEnergyHexColumnx

!!!!!!!!!!!!



double precision function GradientEnergyCube(dx,vbles10,ix,iy,iz)!(vbles,lp,dx)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
use physicaldata
implicit none
double precision, intent (in)::dx,vbles10(N_unk,10,10,10)
integer :: ix,iy,iz
double precision :: phi_p(3),phi_m(3),phi_(3),phi(3,3),phi_tmp(3,3),A0,A1(3),A2(3,3),AA(3,3),fgf,tmp
double precision :: yy(3),y,xx(3),x,px,LapS,H,Lap,d=4d-2,c(6)=(/1,2,3,1,2,3/),rr(3)=(/1,1,2/)
double precision :: p(3,20),pi=3.141592654,q(12),dh(3),pi6=.5235987758,t1(3),t2(3),t3(3)
integer :: i,j,np=12,k,l,iSurfaceEnergy,SS
double precision :: SurfaceEnergy,maxA,psi=1.618033988749895,vbles(N_unk,3,3,3),B1(0:3),B2(0:3,1:3),B3(0:3,1:3,1:3)
double precision :: ID(3,3)=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/)),abs_


GradientEnergyCube=0
   do i=1,3
      do j=1,3
         do k=1,3
            tmp = g_gamma*unk(4,i+ix-2,j+iy-2,k+iz-2,g_block) + (1-g_gamma)*unk(1,i+ix-2,j+iy-2,k+iz-2,g_block)
            vbles(1,i,j,k) = tmp 
         enddo
      enddo
   enddo
   




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


   x = sqrt(xx(1)**2+xx(2)**2+xx(3)**2+1d-12)
   xx=xx/x !debug 7th Oct
  d = g_d_facets
  np=8
  A2 = 0d0
  do i=1,3
     A2(i,i)=1d0
  enddo
  if(.true. .and. x.gt.1d-6)then
        do i=1,3
!           B1(i)=sqrt(xx(i)**2+d**2*x**2)
           B1(i)=sqrt(xx(i)**2+d**2)
        enddo
        do j=1,3
           do i =1,3
              B2(j,i)=(id(j,i)+d**2)*xx(i)/B1(j)
           enddo
        enddo
        do j=1,3     
           do i=1,3
              do k=1,3
                 B3(j,i,k)=1d0/B1(j)*((id(j,i)+d**2)*id(j,k)-B2(j,i)*B2(j,k))
              enddo
           enddo
        enddo
        A0 = 0d0
        A1 = 0d0
        A2 = 0d0
        do j=1,3
           A0 = A0 + B1(j)/(1+d)
           do i = 1,3
              A1(i)=A1(i)+B2(j,i)
              do k = 1,3
                 A2(i,k)=A2(i,k)+B3(j,i,k)
              enddo
           enddo
        enddo
        do i = 1,3
           do j=1,3
              AA(i,j) = (A0*A2(i,j)+A1(i)*A1(j))/(1+d)**2 
           enddo
        enddo
     
     H = 0d0
     do i=1,3
        do j=1,3
           H = H + phi(i,j)*AA(i,j)
        enddo
     enddo
       
     GradientEnergyCube = - H/3.
  else
   
     H = 0d0
     do i=1,3
        do j=1,3
           H = H + phi(i,j)*A2(i,j) 
        enddo
     enddo
     GradientEnergyCube = - H
  endif
  return


  if(.false.)then
     A0 = abs_(xx(1),d,0)+abs_(xx(2),d,0)+abs_(xx(3),d,0)  + epsilon
     A2 = 0d0
     do i=1,3  
        A1(i)   = abs_(xx(i),d,1)  + epsilon * xx(i)
        A2(i,i) = abs_(xx(i),d,2)  + epsilon  
     enddo
     
     do i = 1,3
        do j=1,3
           A2(i,j)=A2(i,j) - epsilon * xx(i)*xx(j)
           AA(i,j) = A0*A2(i,j)+A1(i)*A1(j) 
        enddo
     enddo

     H = 0d0
     do i=1,3
        do j=1,3
           H = H + phi(i,j)*AA(i,j)
        enddo
     enddo

     GradientEnergyCube = - H
     return
  endif

  if(.false.)then
  call AnisotropyPower(g_cube,phi_,np,A2,epsilon)
   H = 0d0
   do i=1,3
      do j=1,3
         H = H + phi(i,j)*A2(i,j) 
      enddo
   enddo
   if(H.gt.0d0.or.H.le.0d0)then
      GradientEnergyCube = - H
   else
      GradientEnergyCube = 0d0
   endif
   return
endif
if(.false.)then
   do i=1,3
      dh=0d0
      dh(i)=d
      yy=xx-dh

      A2(i,i) = SurfaceEnergy(g_cube,yy,np)
      A2(i,i) = A2(i,i) -2d0*SurfaceEnergy(g_cube,xx,np)
      yy=xx+dh
      A2(i,i) = A2(i,i) + SurfaceEnergy(g_cube,yy,np)
      A2(i,i) = A2(i,i)/d**2
   enddo

   do i=1,2
      do j=i+1,3

         dh=0d0
         dh(i)=d
         dh(j)=d
         yy=xx + dh
         A2(i,j) = SurfaceEnergy(g_cube,yy,np)

         yy=xx - dh
         A2(i,j) = A2(i,j) + SurfaceEnergy(g_cube,yy,np)

         dh=0d0
         dh(i)=d
         dh(j)=-d
         yy=xx+dh
         A2(i,j) = A2(i,j) - SurfaceEnergy(g_cube,yy,np)

         yy=xx-dh
         A2(i,j) = A2(i,j) - SurfaceEnergy(g_cube,yy,np)

         A2(i,j)=A2(i,j)*0.25/d**2
         A2(j,i)=A2(i,j)
      enddo
   enddo
endif
   
   H = 0d0
   do i=1,3
      do j=1,3
         H = H + phi(i,j)*A2(i,j) 
      enddo
   enddo
   GradientEnergyCube = - H

   return
! using div d/dgrad(phi) 1/2 A^2: does not work well

  if(.false.)then
     A0 = abs_(phi_(1),d,0)+abs_(phi_(2),d,0)+abs_(phi_(3),d,0)
     A2 = 0d0
     do i=1,3  
        A1(i)   = abs_(phi_(i),d,1)
        A2(i,i) = abs_(phi_(i),d,2)
     enddo
     do i = 1,3
        do j=1,3
           AA(i,j) = A0*A2(i,j)+A1(i)*A1(j)
           if(i.eq.j)AA(i,i)=AA(i,i)+epsilon**2
        enddo
     enddo
     return
  elseif(.true.)then
     A0 = abs_(phi_(1),d,0)+abs_(phi_(2),d,0)+abs_(phi_(3),d,0)  + epsilon*x
     A2 = 0d0
     do i=1,3  
        A1(i)   = abs_(phi_(i),d,1)  + epsilon * phi_(i)/x
        A2(i,i) = abs_(phi_(i),d,2)  + epsilon * 1/x 
     enddo
     
     do i = 1,3
        do j=1,3
           A2(i,j)=A2(i,j) - epsilon * phi_(i)*phi_(j)/x**3
           AA(i,j) = A0*A2(i,j)+A1(i)*A1(j) 
        enddo
     enddo

     H = 0d0
     do i=1,3
        do j=1,3
           H = H + phi(i,j)*AA(i,j)
        enddo
     enddo

     GradientEnergyCube = - H
     return
  else
     if(x.gt.1d-12)then
        B1(0)=d*x
        do i=1,3
           B1(i)=sqrt(phi_(i)**2+d**2*x**2)
        enddo
        do j=0,3
           do i =1,3
              B2(j,i)=(id(j,i)+d**2)*phi_(i)/B1(j)
           enddo
        enddo
        do j=0,3     
           do i=1,3
              do k=1,3
                 B3(j,i,k)=1d0/B1(j)*((id(j,i)+d**2)*id(j,k)-B2(j,i)*B2(j,k))
              enddo
           enddo
        enddo
        A0 = 0d0
        A1 = 0d0
        A2 = 0d0
        do j=1,3
           A0 = A0 + B1(j)
           do i = 1,3
              A1(i)=A1(i)+B2(j,i)
              do k = 1,3
                 A2(i,k)=A2(i,k)+B3(j,i,k)
              enddo
           enddo
        enddo
        do i = 1,3
           do j=1,3
              AA(i,j) = A0*A2(i,j)+A1(i)*A1(j) 
           enddo
        enddo
     else
        AA=0d0
     endif
  endif


   H = 0d0
   do i=1,3
      do j=1,3
         H = H + phi(i,j)*AA(i,j)
      enddo
   enddo







   GradientEnergyCube = - H



   




end function GradientEnergyCube


double precision function GradientEnergyOct(dx,vbles10,ix,iy,iz)!(vbles,lp,dx)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
use physicaldata
implicit none
double precision, intent (in)::dx,vbles10(N_unk,10,10,10)
integer :: ix,iy,iz,m=16
double precision :: phi_p(3),phi_m(3),phi_(3),phi(3,3),phi_tmp(3,3),A0,A1(3),A2(3,3),AA(3,3),fgf,tmp
double precision :: yy(3),y,xx(3),x,px,LapS,H,Lap,d=4d-2,c(6)=(/1,2,3,1,2,3/),rr(3)=(/1,1,2/)
double precision :: p(3,20),pi=3.141592654,q(12),dh(3),pi6=.5235987758,t1(3),t2(3),t3(3),S0,S1(3),S2(3,3)
integer :: i,j,np=12,k,l,iSurfaceEnergy,SS
double precision :: SurfaceEnergy,maxA,psi=1.618033988749895,vbles(N_unk,3,3,3),B1(1:4),B2(1:4,1:3),B3(1:4,1:3,1:3)
double precision :: ID(3,3)=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/)),abs_
double precision :: qq(4,3)=reshape((/1.,1.,1.,1.,1.,-1.,1.,-1.,1.,1.,-1.,-1./),(/4,3/))

!!$
!!$do i=1,3
!!$   do j=1,3
!!$      do k=1,3
!!$         vbles(1,i,j,k)=vbles10(1,i+ix-2,j+iy-2,k+iz-2)
!!$      enddo
!!$   enddo
!!$enddo
GradientEnergyOct=0
   do i=1,3
      do j=1,3
         do k=1,3
            tmp = g_gamma*unk(4,i+ix-2,j+iy-2,k+iz-2,g_block) + (1-g_gamma)*unk(1,i+ix-2,j+iy-2,k+iz-2,g_block)
            vbles(1,i,j,k) = tmp 
         enddo
      enddo
   enddo
   




   phi_(1)=0.5*(vbles(1,3,2,2)-vbles(1,1,2,2))/dx
   phi_(2)=0.5*(vbles(1,2,3,2)-vbles(1,2,1,2))/dx
   phi_(3)=0.5*(vbles(1,2,2,3)-vbles(1,2,2,1))/dx

!!$   phi_p(1)=(vbles(1,3,2,2)-vbles(1,2,2,2))/dx
!!$   phi_p(2)=(vbles(1,2,3,2)-vbles(1,2,2,2))/dx
!!$   phi_p(3)=(vbles(1,2,2,3)-vbles(1,2,2,2))/dx
!!$
!!$   phi_m(1)=-(vbles(1,1,2,2)-vbles(1,2,2,2))/dx
!!$   phi_m(2)=-(vbles(1,2,1,2)-vbles(1,2,2,2))/dx
!!$   phi_m(3)=-(vbles(1,2,2,1)-vbles(1,2,2,2))/dx



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
         xx(i) = abs(phi_(i))
   enddo


   x = sqrt(xx(1)**2+xx(2)**2+xx(3)**2+1d-12)

   xx = xx/x

  d = g_d_facets

     S0 = xx(1)**m+xx(2)**m+xx(3)**m +1d-12
     do i = 1,3
        S1(i) = m*xx(i)**(m-1)
     enddo
     S2=0d0
     do i=1,3
        S2(i,i)=m*(m-1)*xx(i)**(m-2)
     enddo
     
     do j=1,3
        do k=1,3
           AA(j,k)=(1./m)*(2./m-1.)*S0**(2./m-2.)*S1(j)*S1(k)+(1./m)*S0**(2./m-1.)*S2(j,k)
        enddo
     enddo

  do i=1,3
     AA(i,i) = AA(i,i) +d
  enddo


!  tmp = (xx(1)*xx(1))**m+(xx(2)*xx(2))**m+(xx(3)*xx(3))**m +(d*d)**m

!!$  do i=1,3
!!$     do j=1,3
!!$!        AA(i,j) = -4.*xx(i)**(2*m-1)*xx(j)**(2*m-1)*exp((1./m-2.)*log(tmp))*(m-1)
!!$        AA(i,j) = tmp**(1/m-2)*xx(i)**(2*m-1)*xx(j)**(2*m-1)*1./m*(1./m-1.)*4.*m**2

!        AA(i,j)=4.*tmp**(d-2.)*(xx(i)*xx(j))**(2/d-1.)*(1.-1./d)
!!$     enddo
!!$  enddo
!!$  do i=1,3
!!$     AA(i,i)=AA(i,i) + tmp**(1/m-1)*xx(i)**(2*m-2)*(1./m)*2*m*(2*m-1.)
!     AA(i,i)=AA(i,i)+2*xx(i)**(2*m-2)*exp((1./m-1)*log(tmp))*(2*m-1) 
!     AA(i,i)=AA(i,i) +(4./d-2.)*(xx(i)*xx(i))**(1./d-2.)*tmp**(d-1.)
!!$  enddo 
!!$  A0 = log(exp((xx(1)/d)**2+(xx(2)/d)**2+(xx(3)/d)**2))*d**2
!  tmp = exp((xx(1)/d)**2+(xx(2)/d)**2+(xx(3)/d)**2)
!!$  do i=1,3
!!$     A1(i)=2*xx(i)*exp((xx(i)/d)**2)/tmp
!!$  enddo

!!$
!!$  do i=1,2
!!$     do j=i+1,3
!!$        if(i.ne.j)then
!!$           AA(i,j)=- 4.*xx(i)*xx(j)/d**2*exp((xx(i)/d)**2+(xx(j)/d)**2)/tmp**2
!!$           AA(j,i)=AA(i,j)
!!$        endif
!!$     enddo
!!$  enddo
!!$  do i=1,3
!!$     AA(i,i)=(2.+4.*(xx(i)/d)**2)*exp((xx(i)/d)**2)/tmp - 4.*(xx(i)/d)**2*exp(2*(xx(i)/d)**2)/tmp**2
!!$  enddo
!!$
!!$  do i=1,3
!!$     do j=1,3
!!$           AA(i,j)=- 4.*xx(i)*xx(j)/d**2*exp((xx(i)/d)**2+(xx(j)/d)**2)/tmp**2
!!$     enddo
!!$  enddo
!!$  do i=1,3
!!$     AA(i,i) = AA(i,i) +(2.+4.*(xx(i)/d)**2)*exp((xx(i)/d)**2)/tmp
!!$  enddo


  H = 0d0

  do i=1,3
     do j=1,3
        H = H + phi(i,j)*AA(i,j)
     enddo
  enddo


  GradientEnergyOct = - H*0.5
  return
  
  




end function GradientEnergyOct




double precision function GradientEnergyOct2(dx,vbles10,ix,iy,iz)!(vbles,lp,dx)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
use physicaldata
implicit none
double precision, intent (in)::dx,vbles10(N_unk,10,10,10)
integer :: ix,iy,iz,m=16
double precision :: phi_p(3),phi_m(3),phi_(3),phi(3,3),phi_tmp(3,3),A0,A1(3),A2(3,3),AA(3,3),fgf,tmp,minn
double precision :: yy(3),y,xx(3),x,px,LapS,H,Lap,d=4d-2,c(6)=(/1,2,3,1,2,3/),rr(3)=(/1,1,2/)
double precision :: p(3,20),pi=3.141592654,q(12),dh(3),pi6=.5235987758,t1(3),t2(3),t3(3),S0,S1(3),S2(3,3)
integer :: i,j,ii,np=12,k,l,iSurfaceEnergy,SS
double precision :: SurfaceEnergy,maxA,psi=1.618033988749895,vbles(N_unk,3,3,3),B1(1:4),B2(1:4,1:3),B3(1:4,1:3,1:3)
double precision :: ID(3,3)=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/)),abs_
double precision :: qq(4,3)=reshape((/1.,1.,1.,1.,1.,-1.,1.,-1.,1.,1.,-1.,-1./),(/4,3/))


GradientEnergyOct2=0
   do i=1,3
      do j=1,3
         do k=1,3
            tmp = g_gamma*unk(4,i+ix-2,j+iy-2,k+iz-2,g_block) + (1-g_gamma)*unk(1,i+ix-2,j+iy-2,k+iz-2,g_block)
            vbles(1,i,j,k) = tmp 
         enddo
      enddo
   enddo
   




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
         xx(i) = abs(phi_(i))
   enddo


   x = sqrt(xx(1)**2+xx(2)**2+xx(3)**2+1d-12)

   xx = abs(xx/x)

   ii=1

   minn = 1d9
   do i=1,3
      if(xx(i).lt.minn)then
         ii=i
         minn=xx(i)
      endif
   enddo


   j=1
   do i=1,3
      if(i.ne.ii)then
         yy(j)=xx(i)
         j=j+1
      endif
   enddo

   xx = yy

   yy(1)=0.5*(xx(1)+xx(2))
   yy(2)=0.5*abs(xx(1)-xx(2))


  d = g_d_facets

  np=8
  A2 = 0d0
  do i=1,3
     A2(i,i)=1d0
  enddo

  do i=1,3
     B1(i)=sqrt(xx(i)**2+d**2*x**2)
  enddo
  do j=1,3
     do i =1,3
        B2(j,i)=(id(j,i)+d**2)*xx(i)/B1(j)
     enddo
  enddo
  do j=1,3     
     do i=1,3
        do k=1,3
           B3(j,i,k)=1d0/B1(j)*((id(j,i)+d**2)*id(j,k)-B2(j,i)*B2(j,k))
        enddo
     enddo
  enddo
  A0 = 0d0
  A1 = 0d0
  A2 = 0d0
  do j=1,3
     A0 = A0 + B1(j)/(1+d)
     do i = 1,3
        A1(i)=A1(i)+B2(j,i)
        do k = 1,3
           A2(i,k)=A2(i,k)+B3(j,i,k)
        enddo
     enddo
  enddo
  do i = 1,3
     do j=1,3
        AA(i,j) = (A0*A2(i,j)+A1(i)*A1(j))/(1+d)**2 
     enddo
  enddo
  
  
  
  H = 0d0
  
  do i=1,3
     do j=1,3
        H = H + phi(i,j)*AA(i,j)
     enddo
  enddo
  

  GradientEnergyOct2 = - H*0.5

  
  




end function GradientEnergyOct2




double precision function GradientEnergyDodec(dx,vbles10,ix,iy,iz)!(vbles,lp,dx)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
implicit none
double precision, intent (in)::dx,vbles10(N_unk,10,10,10)
integer :: ix,iy,iz
double precision :: phi_(3),phi(3,3),phi_tmp(3,3),A_(3),A(3,3),AA(3,3),fgf
double precision :: yy(3),y,xx(3),x,px,LapS,H,Lap,d=4d-2,c(6)=(/1,2,3,1,2,3/),rr(3)=(/1,1,2/)
double precision :: p(3,20),pi=3.141592654,q(12),dh(3),pi6=.5235987758,t1(3),t2(3),t3(3)
integer :: i,j,np=12,k,l,iSurfaceEnergy
double precision :: SurfaceEnergy,maxA,psi=1.618033988749895,vbles(N_unk,3,3,3)
double precision :: ID(3,3)=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))

do i=1,3
   do j=1,3
      do k=1,3
         vbles(1,i,j,k)=vbles10(1,i+ix-2,j+iy-2,k+iz-2)
      enddo
   enddo
enddo


GradientEnergyDodec=0
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


!   xx = sqrt(xx(1)**2+xx(2)**2+xx(3)**2)

!   xx=xx/x


  d=1d-6
  np=20

   do i=1,3
      dh=0d0
      dh(i)=d
      yy=xx-dh
      A(i,i) = SurfaceEnergy(g_Dodec,yy,np)
      A(i,i) = A(i,i) -2d0*SurfaceEnergy(g_Dodec,xx,np)
      yy=xx+dh
      A(i,i) = A(i,i) + SurfaceEnergy(g_Dodec,yy,np)
      A(i,i) = A(i,i)/d**2
   enddo

   do i=1,2
      do j=i+1,3

         dh=0d0
         dh(i)=d
         dh(j)=d
         yy=xx + dh
         A(i,j) = SurfaceEnergy(g_Dodec,yy,np)

         yy=xx - dh
         A(i,j) = A(i,j) + SurfaceEnergy(g_Dodec,yy,np)

         dh=0d0
         dh(i)=d
         dh(j)=-d
         yy=xx+dh
         A(i,j) = A(i,j) - SurfaceEnergy(g_Dodec,yy,np)

         yy=xx-dh
         A(i,j) = A(i,j) - SurfaceEnergy(g_Dodec,yy,np)

         A(i,j)=A(i,j)*0.25/d**2
         A(j,i)=A(i,j)
      enddo
   enddo



   

   H = 0d0
   do i=1,3
      do j=1,3
         H = H + phi(i,j)*A(i,j) 
      enddo
   enddo


   GradientEnergyDodec = - H





end function GradientEnergyDodec






double precision function GradientEnergyHexColy(dx,vbles10,ix,iy,iz)!(vbles,lp,dx)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
implicit none
double precision, intent (in)::dx,vbles10(N_unk,10,10,10)
integer :: ix,iy,iz
double precision :: phi_(3),phi(3,3),phi_tmp(3,3),A_(3),A(3,3),AA(3,3),fgf,S0,S1(3),S2(3,3),tmp
double precision :: yy(3),y,xx(3),x,px,LapS,H,Lap,d=4d-2,c(6)=(/1,2,3,1,2,3/),rr(3)=(/1,1,2/)
double precision :: p(3,20),pi=3.141592654,q(12),dh(3),pi6=.5235987758,t1(3),t2(3),t3(3),radius
integer :: i,j,np=12,k,l,iSurfaceEnergy,m=16
double precision :: SurfaceEnergy,maxA,psi=1.618033988749895,vbles(N_unk,3,3,3)
double precision :: ID(3,3)=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))

do i=1,3
   do j=1,3
      do k=1,3
         vbles(1,i,j,k)=vbles10(1,i+ix-2,j+iy-2,k+iz-2)
      enddo
   enddo
enddo


GradientEnergyHexColy=0
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


   x = sqrt(xx(1)**2+xx(2)**2+xx(3)**2)
   AA = 0d0
   if(x.gt.1d-6)then

      xx = xx/x
      
      d = g_d_facets
      
      radius = sqrt(g_HexCol(1,1)**2+g_HexCol(2,1)**2+g_HexCol(3,1)**2)
      g_HexCol=g_HexCol/radius
      
      S0 = 0.
      S1 = 0.
      S2 = 0.
      do i=1,12
         tmp=0.
         do l=1,3
            tmp = tmp + xx(l)*g_HexCol(l,i) 
         enddo
         if(tmp.gt.0.)then
            S0 = S0 + tmp**m
            do j=1,3
               S1(j)=S1(j)+m*g_HexCol(j,i)*tmp**(m-1) 
            enddo
            do j=1,3
               do k=1,3
                  S2(j,k)=S2(j,k) + m*(m-1)*g_HexCol(j,i)*g_HexCol(k,i)*tmp**(m-2)
               enddo
            enddo
         endif
      enddo
      S0 = S0 + 1d-12
      
      
      
      do j=1,3
         do k=1,3
            AA(j,k)=(1./m)*(2./m-1.)*S0**(2./m-2.)*S1(j)*S1(k)+(1./m)*S0**(2./m-1.)*S2(j,k)
         enddo
      enddo
   endif
   do i=1,3
      AA(i,i) = AA(i,i) +d
   enddo
   
   

   H = 0d0
   do i=1,3
      do j=1,3
         H = H + phi(i,j)*AA(i,j) 
      enddo
   enddo


   GradientEnergyHexColy = - H
   return
end function GradientEnergyHexColy


double precision function GradientEnergyHexCol(dx,vbles10,ix,iy,iz)!(vbles,lp,dx)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
implicit none
double precision, intent (in)::dx,vbles10(N_unk,10,10,10)
integer :: ix,iy,iz
double precision :: phi_(3),phi(3,3),phi_tmp(3,3),A_(3),A(3,3),AA(3,3),fgf,S0,S1(3),S2(3,3),tmp,f0,f1(3),f2(3,3)
double precision :: y,xx(3),yy(3),x,px,LapS,H,Lap,d=4d-2,c(6)=(/1,2,3,1,2,3/),rr(3)=(/1,1,2/)
double precision :: p(3,20),pi=3.141592654,q(12),dh(3),pi6=.5235987758,t1(3),t2(3),t3(3),radius,p3(12,3,3),p1(12),p2(12,3)
integer :: i,j,np=12,k,l,ii,iSurfaceEnergy,m=8,n_mol,m_mol
double precision :: SurfaceEnergy,maxA,psi=1.618033988749895,vbles(N_unk,3,3,3),AnisPowEnergy,AnisEnergy,SurfaceEnergy0
double precision :: ID(3,3)=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))


H=0.
do i=1,3
   do j=1,3
      do k=1,3
         vbles(1,i,j,k)=vbles10(1,i+ix-2,j+iy-2,k+iz-2)
      enddo
   enddo
enddo


GradientEnergyHexCol=0
   phi_(1)=   0.5*(vbles(1,3,2,2)-vbles(1,1,2,2))/dx

   phi_(2)=   0.5*(vbles(1,2,3,2)-vbles(1,2,1,2))/dx

   phi_(3)=   0.5*(vbles(1,2,2,3)-vbles(1,2,2,1))/dx

   phi(1,1)=         g_hex*(vbles(1,3,2,2)-2d0*vbles(1,2,2,2)+vbles(1,1,2,2))/dx**2
   phi(1,1)=phi(1,1)+0.5*(1-g_hex)*(vbles(1,3,3,2)-2d0*vbles(1,2,3,2)+vbles(1,1,3,2))/dx**2
   phi(1,1)=phi(1,1)+0.5*(1-g_hex)*(vbles(1,3,1,2)-2d0*vbles(1,2,1,2)+vbles(1,1,1,2))/dx**2

   phi(2,2)=         g_hex*(vbles(1,2,3,2)-2d0*vbles(1,2,2,2)+vbles(1,2,1,2))/dx**2
   phi(2,2)=phi(2,2)+0.5*(1-g_hex)*(vbles(1,3,3,2)-2d0*vbles(1,3,2,2)+vbles(1,3,1,2))/dx**2
   phi(2,2)=phi(2,2)+0.5*(1-g_hex)*(vbles(1,1,3,2)-2d0*vbles(1,1,2,2)+vbles(1,1,1,2))/dx**2

   phi(3,3)=(vbles(1,2,2,3)-2d0*vbles(1,2,2,2)+vbles(1,2,2,1))/dx**2

   phi(1,2) = 0.25*(vbles(1,3,3,2)+vbles(1,1,1,2)-vbles(1,3,1,2)-vbles(1,1,3,2))/dx**2
   phi(1,3) = 0.25*(vbles(1,3,2,3)+vbles(1,1,2,1)-vbles(1,3,2,1)-vbles(1,1,2,3))/dx**2
   phi(2,3) = 0.25*(vbles(1,2,3,3)+vbles(1,2,1,1)-vbles(1,2,3,1)-vbles(1,2,1,3))/dx**2
   phi(2,1) = phi(1,2)
   phi(3,1) = phi(1,3)
   phi(3,2) = phi(2,3)

   if(.false.)then
      ii=0
      do k=1,3
         do j=1,3
            do i=1,3 !i needs to be faster for the 1D stencil
               ii=ii+1
               phi(1,1) = phi(1,1)+vbles(1,i,j,k)*g_Sx(ii)/dx**2
               phi(2,2) = phi(2,2)+vbles(1,i,j,k)*g_Sy(ii)/dx**2
               phi(3,3) = phi(3,3)+vbles(1,i,j,k)*g_Sz(ii)/dx**2
            enddo
         enddo
      enddo
   endif
   do i=1,3
         xx(i) = phi_(i)
   enddo


   x = sqrt(xx(1)**2+xx(2)**2+xx(3)**2)
   S2 = 0.
   S1 = 0.
   if(x.gt.1d-6)then
      
      !xx = xx/x
      xx = xx*g_lambda
      if(.false.)then ! uses approx of abs and two maximums. Does not work
         do j=1,3
            do k=1,3
               S2(j,k)=AnisEnergy(g_HexCol,xx,12,j,0)*AnisEnergy(g_HexCol,xx,12,k,0) + AnisEnergy(g_HexCol,xx,12,0,0)*AnisEnergy(g_HexCol,xx,12,j,k)
            enddo
         enddo
      endif
      

      if(.true.)then
         yy = xx
         do k=1,3
            do j=1,3
               do n_mol = 1,g_mol_n - 1
                  do m_mol = 1,g_mol_n - 1
                     yy(k) = xx(k) - g_mollifier(n_mol,1)
                     if(k.eq.j)then
                        yy(k) = yy(k) - g_mollifier(m_mol,1)
                     else
                        yy(j) = xx(j) - g_mollifier(m_mol,1)
                     endif
                     S2(j,k) = S2(j,k) + g_mollifier(n_mol,3)*g_mollifier(m_mol,3)*SurfaceEnergy0(g_HexCol,yy,12)
                  enddo
               enddo
            enddo
         enddo
      endif

      if(.false.)then
         yy = xx
         do k=1,3
            do n_mol = 1,g_mol_n - 1
               yy(k) = xx(k) - g_mollifier(n_mol,1)
               do j=1,3
                  S2(j,k) = S2(j,k) + g_mollifier(n_mol,3)*SurfaceEnergy(g_HexCol,yy,12,j)
               enddo
            enddo
            yy(k)=xx(k)
         enddo
      endif
      
      if(.false.)then
         do j=1,3
            do k=1,3
               S2(j,k)=AnisPowEnergy(g_HexCol,xx,12,j)*AnisPowEnergy(g_HexCol,xx,12,k)
            enddo
         enddo
      endif
      if(.false.)then
         do i = 1,3
            do j =1,3
               phi_tmp(i,j)=S2(i,j)
               do k = 1,3
                  phi_tmp(i,j) = phi_tmp(i,j) - S2(i,k)*xx(k)*xx(j)-S2(k,j)*xx(k)*xx(i)
                  do l = 1,3
                     phi_tmp(i,j) = phi_tmp(i,j) + 2.*S2(k,l)*xx(k)*xx(l)*xx(i)*xx(j)
                  enddo
               enddo
            enddo
         enddo
         
         S2 = phi_tmp
      endif

      
      H = 0d0
      do i=1,3
         do j=1,3
            H = H + phi(i,j)*S2(i,j) 
         enddo
      enddo
   endif

   GradientEnergyHexCol = - H

   return
end function GradientEnergyHexCol



double precision function GradientEnergyHexColAug(dx,vbles10,ix,iy,iz)!(vbles,lp,dx)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
implicit none
double precision, intent (in)::dx,vbles10(N_unk,10,10,10)
integer :: ix,iy,iz
double precision :: phi_(3),phi(3,3),phi_tmp(3,3),A_(3),A(3,3),AA(3,3),fgf,S0,S1(3),S2(3,3),tmp,f0,f1(3),f2(3,3)
double precision :: yy(3),y,xx(3),x,px,LapS,H,Lap,d=4d-2,c(6)=(/1,2,3,1,2,3/),rr(3)=(/1,1,2/)
double precision :: p(3,20),pi=3.141592654,q(12),dh(3),pi6=.5235987758,t1(3),t2(3),t3(3),radius,p3(12,3,3),p1(12),p2(12,3)
integer :: i,j,np=12,k,l,iSurfaceEnergy,m=8
double precision :: SurfaceEnergy,maxA,psi=1.618033988749895,vbles(N_unk,3,3,3)
double precision :: ID(3,3)=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))

do i=1,3
   do j=1,3
      do k=1,3
         vbles(1,i,j,k)=vbles10(1,i+ix-2,j+iy-2,k+iz-2)
      enddo
   enddo
enddo


GradientEnergyHexColAug=0
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


   x = sqrt(xx(1)**2+xx(2)**2+xx(3)**2)
   f2 = 0d0
   if(x.gt.1d-6)then
      
      xx = xx/x
      
      
      
      
      d = g_d_facets
      
      radius = sqrt(g_HexCol(1,1)**2+g_HexCol(2,1)**2+g_HexCol(3,1)**2)
      g_HexCol=g_HexCol/radius
      
      do i = 1,12
         do j=1,3
            do k=1,3
               p3(i,j,k) = g_HexCol(j,i)*g_HexCol(k,i) + id(j,k)*d
            enddo
         enddo
      enddo

      
      
      if(.true.)then
      p1 = 0.
      p2 = 0.
      do i = 1,12
         do j=1,3
            do k=1,3
               p1(i) = p1(i) + 0.5*xx(j)*p3(i,j,k)*xx(k)
               p2(i,k) = p2(i,k) + xx(j)*p3(i,j,k)
            enddo
         enddo
      enddo
         S0 = 0.
         S1 = 0.
         S2 = 0.
         do i = 1,12
            if(p1(i).gt.0.)then
               S0 = S0 + p1(i)**m
               do j=1,3
                  S1(j) = S1(j) + m*p1(i)**(m-1)*p2(i,j)
                  do k=1,3
                     S2(j,k) = S2(j,k) + m*(m-1)*p1(i)**(m-2)*p2(i,j)*p2(i,k) + m*p1(i)**(m-1)*p3(i,j,k)
                  enddo
               enddo
            endif
         enddo
         !      f0 = S0**(1./m) !not needed apert from calculating f1 and f2
         do j=1,3
            !         f1(j) = (1./m)*S0**(1./m -1.)*S1(j) ! not needed apart fro finding f2
            do k=1,3
               f2(j,k)=(1./m)*(1./m -1.)*S0**(1./m -2.)*S1(j)*S1(k) + (1./m)*S0**(1./m -1.)*S2(j,k)
            enddo
         enddo
         

      else
         A = 0.
         np=12
         
         do i=1,3
            dh=0d0
            dh(i)=d
            yy=xx-dh
            A(i,i) = SurfaceEnergy(g_Hexcol,yy,np)
            A(i,i) = A(i,i) -2d0*SurfaceEnergy(g_Hexcol,xx,np)
            yy=xx+dh
            A(i,i) = A(i,i) + SurfaceEnergy(g_Hexcol,yy,np)
            A(i,i) = A(i,i)/d**2
         enddo
         
         do i=1,2
            do j=i+1,3
               
               dh=0d0
               dh(i)=d
               dh(j)=d
               yy=xx + dh
               A(i,j) = SurfaceEnergy(g_Hexcol,yy,np)
               
               yy=xx - dh
               A(i,j) = A(i,j) + SurfaceEnergy(g_Hexcol,yy,np)
               
               dh=0d0
               dh(i)=d
               dh(j)=-d
               yy=xx+dh
               A(i,j) = A(i,j) - SurfaceEnergy(g_Hexcol,yy,np)
               
               yy=xx-dh
               A(i,j) = A(i,j) - SurfaceEnergy(g_Hexcol,yy,np)
               
               A(i,j)=A(i,j)*0.25/d**2
               A(j,i)=A(i,j)
            enddo
         enddo
         f2 = A    
         
      endif
   endif



   
   

   H = 0d0
   do i=1,3
      do j=1,3
         H = H + phi(i,j)*f2(i,j) 
      enddo
   enddo


   GradientEnergyHexColAug = - H
   return
end function GradientEnergyHexColAug



double precision function GradientEnergyTrunOct(dx,vbles10,ix,iy,iz)!(vbles,lp,dx)
use solution_parameters
use multigrid_parameters  ! for min_dx, total_vars
implicit none
double precision, intent (in)::dx,vbles10(N_unk,10,10,10)
integer :: ix,iy,iz
double precision :: phi_(3),phi(3,3),phi_tmp(3,3),A_(3),A(3,3),AA(3,3),fgf,S0,S1(3),S2(3,3),tmp
double precision :: yy(3),y,xx(3),x,px,LapS,H,Lap,d=4d-2,c(6)=(/1,2,3,1,2,3/),rr(3)=(/1,1,2/)
double precision :: p(3,20),pi=3.141592654,q(12),dh(3),pi6=.5235987758,t1(3),t2(3),t3(3)
integer :: i,j,np=12,k,l,iSurfaceEnergy,m=32
double precision :: SurfaceEnergy,maxA,psi=1.618033988749895,vbles(N_unk,3,3,3)
double precision :: ID(3,3)=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))

do i=1,3
   do j=1,3
      do k=1,3
         vbles(1,i,j,k)=vbles10(1,i+ix-2,j+iy-2,k+iz-2)
      enddo
   enddo
enddo


GradientEnergyTrunOct=0
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


   x = sqrt(xx(1)**2+xx(2)**2+xx(3)**2+1d-12)
   
   xx = xx/x
   
   d = g_d_facets
   
   
   S0 = 0.
   S1 = 0.
   S2 = 0.
   do i=1,24
      tmp=0.
      do l=1,3
         tmp = tmp + xx(l)*g_TrunOct(l,i) 
      enddo
      if(tmp.gt.0.)then
         S0 = S0 + tmp**m
         do j=1,3
            S1(j)=S1(j)+m*g_TrunOct(j,i)*tmp**(m-1) 
         enddo
         do j=1,3
            do k=1,3
               S2(j,k)=S2(j,k) + m*(m-1)*g_TrunOct(j,i)*g_TrunOct(k,i)*tmp**(m-2)
            enddo
         enddo
      endif
   enddo
   S0 = S0 + 1d-12


!!$   S1=0.
!!$   do j=1,3
!!$      do i=1,24
!!$         tmp=0.
!!$         do k=1,3
!!$            tmp=tmp + xx(k)*g_TrunOct(k,i)
!!$         enddo
!!$         if(tmp.gt.0.)then
!!$            S1(j)=S1(j)+m*g_TrunOct(j,i)*tmp**(m-1) 
!!$         endif
!!$      enddo
!!$   enddo
!!$   S2 =0.
!!$   do j=1,3
!!$      do k=1,3
!!$         do i=1,24
!!$            tmp = 0.
!!$            do l = 1,3
!!$               tmp = tmp + xx(l)*g_TrunOct(l,i)
!!$            enddo
!!$            if(tmp.gt.0.)then
!!$               S2(j,k)=S2(j,k) + m*(m-1)*g_TrunOct(j,i)*g_TrunOct(k,i)*tmp**(m-2)
!!$            endif
!!$         enddo
!!$      enddo
!!$   enddo
!!$
!!$     S0 = xx(1)**m+xx(2)**m+xx(3)**m +1d-12
!!$     do i = 1,3
!!$        S1(i) = m*xx(i)**(m-1)
!!$     enddo
!!$     S2=0d0
!!$     do i=1,3
!!$        S2(i,i)=m*(m-1)*xx(i)**(m-2)
!!$     enddo
     
   do j=1,3
      do k=1,3
         AA(j,k)=(1./m)*(2./m-1.)*S0**(2./m-2.)*S1(j)*S1(k)+(1./m)*S0**(2./m-1.)*S2(j,k)
      enddo
   enddo
   
   do i=1,3
      AA(i,i) = AA(i,i) +d
   enddo
   H = 0d0
   do i=1,3
      do j=1,3
         H = H + phi(i,j)*AA(i,j) 
      enddo
   enddo


   GradientEnergyTrunOct = - H
   return

!   xx = sqrt(xx(1)**2+xx(2)**2+xx(3)**2)

!   xx=xx/x


  d=1d-4
  np=24

   do i=1,3
      dh=0d0
      dh(i)=d
      yy=xx-dh
      A(i,i) = SurfaceEnergy(g_TrunOct,yy,np)
      A(i,i) = A(i,i) -2d0*SurfaceEnergy(g_TrunOct,xx,np)
      yy=xx+dh
      A(i,i) = A(i,i) + SurfaceEnergy(g_TrunOct,yy,np)
      A(i,i) = A(i,i)/d**2
   enddo

   do i=1,2
      do j=i+1,3

         dh=0d0
         dh(i)=d
         dh(j)=d
         yy=xx + dh
         A(i,j) = SurfaceEnergy(g_TrunOct,yy,np)

         yy=xx - dh
         A(i,j) = A(i,j) + SurfaceEnergy(g_TrunOct,yy,np)

         dh=0d0
         dh(i)=d
         dh(j)=-d
         yy=xx+dh
         A(i,j) = A(i,j) - SurfaceEnergy(g_TrunOct,yy,np)

         yy=xx-dh
         A(i,j) = A(i,j) - SurfaceEnergy(g_TrunOct,yy,np)

         A(i,j)=A(i,j)*0.25/d**2
         A(j,i)=A(i,j)
      enddo
   enddo



   

   H = 0d0
   do i=1,3
      do j=1,3
         H = H + phi(i,j)*A(i,j) 
      enddo
   enddo


   GradientEnergyTrunOct = - H





end function GradientEnergyTrunOct


double precision function AnisEnergy(p,x,n,d1,d2)
use solution_parameters
implicit none
double precision, intent (in):: p(3,100),x(3)
integer, intent (in) :: n,d1,d2
integer i,j, imaximum
double precision :: maximum,q(100),qi,dt(2),t,dtt



q=0d0
do j=1,n
   do i=1,3
      q(j) = q(j)+  p(i,j)*x(i)
   enddo
enddo

!qmax=(epsilon*sqrt(x(1)**2+x(2)**2+x(3)**2)+maximum(q,n))/(1d0+epsilon)

AnisEnergy = maximum(q,n)
return
i=imaximum(q,n)
qi=q(i)
q(i)=0d0
j=imaximum(q,n)
q(i)=qi


if(d1.eq.0)then
   AnisEnergy = 0.5*(q(i)+q(j))+0.5*sqrt((q(i)-q(j))**2+g_d_facets**2*(x(1)**2+x(2)**2+x(3)**2))
elseif(d2.eq.0)then
   t = (q(i)-q(j))**2+g_d_facets**2*(x(1)**2+x(2)**2+x(3)**2)
   dt(1)= 2.*(q(i)-q(j))*(p(d1,i)-p(d1,j))+2.*g_d_facets**2*x(d1)
   AnisEnergy = 0.5*(p(d1,i)+p(d1,j))+0.5*0.5/sqrt(t)*dt(1)
else
   t =(q(i)-q(j))**2+g_d_facets**2*(x(1)**2+x(2)**2+x(3)**2)
   dt(1)=2.*(q(i)-q(j))*(p(d1,i)-p(d1,j))+2.*g_d_facets**2*x(d1)
   dt(2)=2.*(q(i)-q(j))*(p(d2,i)-p(d2,j))+2.*g_d_facets**2*x(d2)
   dtt=2.*(p(d2,i)-p(d2,j))*(p(d1,i)-p(d1,j))
   if(d1.eq.d2)then
      dtt = dtt + 2.*g_d_facets**2
   endif
   AnisEnergy = -0.5*0.5*0.5/sqrt(t)**3*dt(1)*dt(2)+0.5*0.5/sqrt(t)*dtt
endif


end function AnisEnergy




integer function imaximum(p,n)
implicit none
double precision, intent (in):: p(100)
integer, intent (in) :: n
integer i,j

j=1
do i = 2, n
if(p(i).gt.p(j))then
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
double precision:: q(100),amaximum,tmp




maximum = p(1)
j=1
do i = 2, n
   if(p(i).gt.p(j))then
      j=i
   endif
enddo
maximum = p(j)




end function maximum

double precision function amaximum(p,n)
use solution_parameters
implicit none
double precision, intent (in):: p(100)
integer, intent (in) :: n
integer i,j
double precision :: s,m,maximum
amaximum = 0d0

!s=g_amax


do i = 1, n
!   amaximum = amaximum + ((((s+(1.0-s)*p(i))**8)**8)**8)
   amaximum = amaximum + ((((((p(i)*p(i))**2)**2)**2)**2)**2)
enddo

amaximum=sqrt(sqrt(sqrt(sqrt(sqrt(sqrt(amaximum)))))) 



end function amaximum


double precision function Abs_(x,h,d)
double precision, intent (in):: x,h
integer, intent (in):: d
if(d.eq.0)then
   Abs_ = sqrt(x**2+h**2) - h
elseif(d.eq.1)then
   Abs_ = x/sqrt(x**2+h**2)
else
   Abs_ = -x**2/sqrt(x**2+h**2)**3 + 1d0/sqrt(x**2+h**2)
endif

end function Abs_



          
subroutine AnisotropyPower(p,x,n,A_,epsilon)
implicit none
double precision, intent (in):: p(3,8),x(3),epsilon
integer, intent (in) :: n
double precision, intent (out)::A_(3,3)
integer :: i,j,k,nn,m=512
double precision :: B(100),xx,A=0d0, A0,A1(3),A2(3,3),tmp,amaximum,phi(3,3),d=1d-4

!current extimate of max |grad(phi)| = 0.3
xx = sqrt(x(1)**2 + x(2)**2 + x(3)**2) + 1d-12
A0=0d0
A1=0d0
A2=0d0
A = 0d0
B=1d-12
A_=0d0
do i=1,3
   A_(i,i)=1d0
enddo
if(xx.gt.1d-12)then
   nn = 0
   do i = 1,n
      tmp = 0d0
      do j=1,3
         tmp = tmp + p(j,i)*x(j)/xx
      enddo
      B(i) = tmp
   enddo

   do i=1,n
!      A = A + exp(log(B(i))*m)
!       A = A + ((((((B(i)**2)**2)**2)**2)**2)**2)
       A = A + (((((((((B(i)*B(i))**2)**2)**2)**2)**2)**2)**2)**2)
   enddo
!   A0 = exp(log(A)/m)
!   A0 = sqrt(sqrt(sqrt(sqrt(sqrt(sqrt(A))))))
   A0 =  sqrt(sqrt(sqrt(sqrt(sqrt(sqrt(sqrt(sqrt(sqrt(A)))))))))
   do j=1,n
      B(j) = B(j)/A0
   enddo
   do i = 1,3
      do j=1,n
!         A1(i) = A1(i) + exp(log(B(j))*(m-1d0))*p(i,j)
!         if(abs(B(j)).gt.0d0) A1(i) = A1(i) + (((B(j)**3)**3)**7)*p(i,j)
         if(abs(B(j)).gt.0d0) A1(i) = A1(i) + (B(j)*B(j))**255*B(j)*p(i,j)
      enddo
      
  enddo
   do i = 1,3
      do k = 1,3
         do j=1,n
!            A2(i,k) = A2(i,k) +  (m-1d0)/A0*exp(log(B(j))*(m-2d0))*p(i,j)*(p(k,j)-B(j)*A1(k))
            A2(i,k) = A2(i,k) +  (m-1d0)/A0*(B(j)*B(j))**255*p(i,j)*(p(k,j)-B(j)*A1(k))
         enddo
      enddo
   enddo
   
   A0 =  A0 + epsilon
   A1 =  A1 + epsilon * x/xx
   do i=1,3
      A2(i,i) = A2(i,i) + epsilon
   enddo
   do i =1,3
      do j=1,3
         A2(i,j) = A2(i,j)-epsilon*x(i)*x(j)/xx**2
      enddo
   enddo
   
   do i=1,3
      do j=1,3
         A_(i,j) = A0*A2(i,j) + A1(i)*A1(j) 

      enddo
   enddo
endif


end subroutine AnisotropyPower


double precision function GradientEnergyIsotropic(dx,vbles10,ix,iy,iz)
  use solution_parameters
  use paramesh_dimensions 
  implicit none
  double precision, intent (in)::dx,vbles10(N_unk,10,10,10)
  integer, intent (in):: ix,iy,iz
  integer :: i,j,k
  double precision :: vbles(N_unk,3,3,3),dphi0=0d0,dphi1=0d0,dphi2=0d0




   do i=1,3
      do j=1,3
         do k=1,3
            vbles(1,i,j,k)=vbles10(1,i+ix-2,j+iy-2,k+iz-2)
         enddo
      enddo
   enddo
   dphi2 = vbles(1,1,2,2)+vbles(1,3,2,2)+vbles(1,2,1,2)+vbles(1,2,3,2)+vbles(1,2,2,3)+vbles(1,2,2,1)-6.*vbles(1,2,2,2)
   GradientEnergyIsotropic = -dphi2/dx**2   

end function  GradientEnergyIsotropic



!T_K converts non-dimensional parameter T to dimensioned T_K(T)
double precision function T_K(T)
use solution_parameters
implicit none
double precision, intent (in)::T


 
 T_K = g_T0*(1.+T)

 end function T_K

!! sharpensphi for use in c eqn
double precision function sharpen(phi,n)
double precision, intent (in):: phi
integer, intent(in) :: n
double precision :: U

     if(phi.gt.0.5+0.5d0/n)then
        sharpen = 1d0
     elseif(phi.lt.0.5-0.5/n)then
        sharpen=0d0
     else
        sharpen = (phi-0.5)*n+0.5
     endif

!!$
!!$
!!$if(phi.lt.0.5)then
!!$U = phi/(1-phi)**n
!!$sharpen = U/(1d0+U)
!!$else
!!$U = (1-phi)/(phi)**n
!!$sharpen = 1/(1d0+U)
!!$endif
!
end function sharpen

double precision function SoluteRHS(dx,c_c,a_trap,vbles10,ix,iy,iz,mode_1,dt,rf4,rf5,rf6,lb)
use solution_parameters
use physicaldata
implicit none
double precision, intent (in):: dx,c_c,a_trap,vbles10(N_unk,10,10,10),dt,rf4,rf5,rf6
integer,intent (in) :: ix,iy,iz,mode_1,lb
integer ii,jj,ip,ih,jv,kd,i,j,k
double precision c,T,Phi,FreeEnergy,D,norm,X_phi !
double precision :: potential,D_c,f_c,MolPerVol,sharpen
double precision :: c_dot=0d0,phi_dot=0d0,beta,phi_,PhiSoluteRHS=0d0,phaseRHS,phaseRHSbulk
double precision :: phidot, abs_grad_phi,u,FE,vbles(N_unk,3,3,3)
double precision :: Y_phi,gh_phi,AlFey1
double precision :: c_neigh, c_centr, y_neigh, y_centr, phi_centr,phi_neigh!
double precision :: phi_t,X_t,vbles10fine(N_unk,10,10,10)




do i=1,3
   do j=1,3
      do k=1,3
         do ip=1,N_unk
!            vbles(ip,i,j,k) = vbles10(ip,i+ix-2,j+iy-2,k+iz-2)
            vbles(ip,i,j,k) = g_gamma*unk(4+(ip-1)*6,i+ix-2,j+iy-2,k+iz-2,g_block) +  (1-g_gamma)*unk(1+(ip-1)*6,i+ix-2,j+iy-2,k+iz-2,g_block)
         enddo
      enddo
   enddo
enddo


!!$if(New_potential)then!!
!!$do i=1,3
!!$   do j=1,3
!!$      do k=1,3
!!$
!!$         phi =  g_gamma*unk(4,i+ix-2,j+iy-2,k+iz-2,g_block) +  (1-g_gamma)*unk(1,i+ix-2,j+iy-2,k+iz-2,g_block)
!!$
!!$         vbles(1,i,j,k) = Y_phi(gh_phi(phi,0),0)
!!$      enddo
!!$   enddo
!!$enddo
!!$endif




SoluteRHS=0.


!!
!x-coordinate transform values
!! Right
 ih=3
 jv=2
 kd=2
 include 'c_stencil.f90'
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Left!!
 ih=1
 jv=2
 kd=2
 include 'c_stencil.f90'
 
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

 !!!!!!!!!!!!!!!!!!
 !front
 ih=2
 jv=2
 kd=1
 include 'c_stencil.f90'
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !back
 ih=2
 jv=2
 kd=3
 include 'c_stencil.f90'




! SoluteRHS = SoluteRHS/(g_Dch*dx*dx*g_R*g_T0) 
 SoluteRHS = SoluteRHS/(g_Dch*dx*dx)


 
 if(g_plapp)then !!

    SoluteRHS = g_accel_factor*SoluteRHS - g_accel_factor*(g_xx1-g_xx0)*gh_phi(vbles(1,2,2,2),1)*PhaseRHS(dx,vbles10,ix,iy,iz)
!!$    SoluteRHS = SoluteRHS - g_accel_factor*(g_xx1-g_xx0)*gh_phi(vbles(1,2,2,2),1)*PhaseRHS(dx,vbles10,ix,iy,iz)

!!$ ! Remember to compensate for g_accel when choosing g_Dch. 
 endif

return

end function SoluteRHS


double precision function potential(phi)
use solution_parameters

implicit none

double precision, intent (in)::phi
double precision :: c=0.01
potential =  g_barrier*2d0*phi*(1d0-phi)*(1d0-2d0*phi)

!!$if(phi.lt.c)then
!!$   potential = (1.-2.*c)/c*phi
!!$elseif(phi.lt.1.0-c)then
!!$   potential = 1d0 - 2.*phi
!!$else
!!$   potential = (1.-2.*c)/c*(1.-phi)
!!$endif
end function potential


double precision function PhaseRHS(dx,vbles10,ix,iy,iz)
use solution_parameters
use physicaldata
implicit none
double precision, intent(in)::dx,vbles10(N_unk,10,10,10)
integer, intent(in):: ix,iy,iz
double precision :: vbles(N_unk,3,3,3)
double precision potential,FreeEnergy,potential3D
double precision GradientEnergyIsotropic,GradientEnergyHexPrismVol2,GradientEnergyDodec,GradientEnergyTrunOct,GradientEnergyHexColumn,GradientEnergyCube,GradientEnergyOct,GradientEnergyHexCol
double precision ::M_tildec,phi,c,T ,phi_,dphi_,gh_phi,Y_phi,phiz,phix,phiy,gradphi, M_,xx(3),anisenergy,hexcol(3,12)
integer :: i,j,k,ip

do i=1,3
   do j=1,3
      do k=1,3
         do ip=1,N_unk
!            vbles(ip,i,j,k)=vbles10(ip,i+ix-2,j+iy-2,k+iz-2)
            vbles(ip,i,j,k)=g_gamma*unk(4+(ip-1)*6,i+ix-2,j+iy-2,k+iz-2,g_block) + (1-g_gamma)*unk(1+(ip-1)*6,i+ix-2,j+iy-2,k+iz-2,g_block)
         enddo
      enddo
   enddo
enddo


!!$xx(1) = phix
!!$xx(2) = phiy
!!$xx(3) = phiz

!hexcol = g_hexcol
!!$do i=1,12
!!$hexcol(1,i) = g_hexcol(1,i)*(1d0-kinetic_force)
!!$hexcol(2,i) = g_hexcol(2,i)*(1d0-kinetic_force)
!!$hexcol(3,i) = g_hexcol(3,i)*kinetic_force
!!$enddo


!M_ = (1.0-kinetic_force) + kinetic_force*abs(xx(3))

!phiz = 0.1+abs(0.5*(vbles(1,2,2,3)-vbles(1,2,2,1)))/sqrt(phix**2+phiy**2+phiz**2+1d-6)


phi = vbles(1,2,2,2)
c=vbles(2,2,2,2)
T=0d0

phi = g_gamma*unk(4,ix,iy,iz,g_block) + (1-g_gamma)*unk(1,ix,iy,iz,g_block)!vbles(1,2,2)
c   = g_gamma*unk(4+6,ix,iy,iz,g_block)+ (1-g_gamma)*unk(7,ix,iy,iz,g_block)!vbles(2,2,2)
T=0d0


!!$phaseRHS = -(GradientEnergyDodec(dx,vbles10,ix,iy,iz)+potential(Phi)/g_lambda**2&
!!$           +FreeEnergy(c,T,Phi,1,.true.)/(g_L(1)*g_lambda**2))

!!$if(New_potential)then
!!$
!!$
!!$   phi_  = gh_phi(phi,0)
!!$   dphi_ = gh_phi(phi,1)
!!$   phaseRHS = -(GradientEnergyTrunOct(dx,vbles10,ix,iy,iz)+potential3D(Phi,1)/g_lambda**2&
!!$        +FreeEnergy(c,T,Y_phi(Phi_,0),1,.true.)*Y_phi(phi_,1)*dphi_/(g_L(1)*g_lambda**2))

!!$   phaseRHS = -M_*(GradientEnergyHexColumn(dx,vbles10,ix,iy,iz)+potential3D(Phi,1)/g_lambda**2&
!!$        +FreeEnergy(c,T,Y_phi(Phi_,0),1,.true.)*Y_phi(phi_,1)*dphi_/(g_L(1)*g_lambda**2))
!!$else

!!$   phaseRHS = -M_*(GradientEnergyHexColumn(dx,vbles10,ix,iy,iz)+potential3D(Phi,1)/g_lambda**2&
!!$        +FreeEnergy(c,T,Phi,1,.true.)/(g_L(1)*g_lambda**2))



!!$
!M_ = 1d0

if(g_shape.eq.1)then
   if(aspect_ratio.ne.-1.)then
      phix = 0.5*(vbles(1,3,2,2)-vbles(1,1,2,2))
      phiy = 0.5*(vbles(1,2,3,2)-vbles(1,2,1,2))
      phiz = 0.5*(vbles(1,2,2,3)-vbles(1,2,2,1))
      gradphi = sqrt(phix**2+phiy**2+phiz**2+1d-6)
      
      xx(1) = phix/gradphi
      xx(2) = phiy/gradphi
      xx(3) = phiz/gradphi
!      M_ = anisenergy(g_cube,xx,8)/1.280789275
      M_ = abs(xx(1))
      
   else
      M_ = 1
   endif

   phaseRHS = -M_*(GradientEnergyCube(dx,vbles10,ix,iy,iz)+potential3D(Phi,1)/g_lambda**2&
        +FreeEnergy(c,T,Phi,1,.true.)/(g_L(1)*g_lambda**2))
else
   phix = 0.5*(vbles(1,3,2,2)-vbles(1,1,2,2))
   phiy = 0.5*(vbles(1,2,3,2)-vbles(1,2,1,2))
   phiz = 0.5*(vbles(1,2,2,3)-vbles(1,2,2,1))
   gradphi = sqrt(phix**2+phiy**2+phiz**2+1d-6)
   
   xx(1) = phix/gradphi
   xx(2) = phiy/gradphi
   xx(3) = phiz/gradphi
   M_ = anisenergy(g_hexCol,xx,12)
   phaseRHS = -M_*(GradientEnergyHexCol(dx,vbles10,ix,iy,iz)+potential3D(Phi,1,xx)/g_lambda**2&
        +FreeEnergy(c,T,Phi,1,.true.)/(g_L(1)*g_lambda**2))
endif

!!$   phaseRHS = -M_*(GradientEnergyOct(dx,vbles10,ix,iy,iz)+potential3D(Phi,1)/g_lambda**2&
!!$        +FreeEnergy(c,T,Phi,1,.true.)/(g_L(1)*g_lambda**2))
!!$M_ = 1d0
!!$phaseRHS = -M_*(GradientEnergyIsotropic(dx,vbles10,ix,iy,iz)+potential3D(Phi,1)/g_lambda**2&
!!$           +FreeEnergy(c,T,Phi,1,.true.)/(g_L(1)*g_lambda**2))
!

!!$   phaseRHS = -M_*(GradientEnergyTrunOct(dx,vbles10,ix,iy,iz)+potential3D(Phi,1)/g_lambda**2&
!!$        +FreeEnergy(c,T,Phi,1,.true.)/(g_L(1)*g_lambda**2))


!!!$endif




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
  double precision unk_(M_unk),MolPerVol,vbles10(N_unk,10,10,10)
  integer ip,n,ii,jj,kk





  do ii=1,g_order
     do jj=1,g_order
        do kk=1,g_order
           do ip=1,N_unk
              vbles10(ip,ii,jj,kk)=unk(1+(ip-1)*nunkvbles,ii,jj,kk,lb)
           enddo
        enddo
     enddo
  enddo


  call get_RHS_Fi(dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,vbles10,i,j,k,lb)

end subroutine get_RHS_Fi_2d

subroutine get_RHS_Fi(dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,vbles10,ix,iy,iz,lb)
  use solution_parameters
  use physicaldata
  implicit none
  integer, intent(in)::mode_1,nunkvbles,lnblocks,ix,iy,iz,lb
  double precision, intent(in)::dx,dt,dtold
  double precision, intent(inout):: Fi(N_unk)
  double precision, intent(in):: vbles10(N_unk,10,10,10)
  double precision :: c,T,phi,c_dot,phi_dot,vbles(N_unk,3,3,3),weight=0.9
  double precision rfactor,rf4,rf5,rf6,S,GradientEnergy,potential
  double precision SoluteRHS,SoluteRHSHP,PhaseRHS,T_k,phi_,FreeEnergy,MolPerVol,a_trap,grad_phi,grad_c,sum
  logical Not_number
  integer ip,n,i,j,k
  double precision :: phi_t,phix,phiy,phiz,norm,norm0
  
  g_dt = dt
  g_dtold = dtold !for velocity calculation in mp_output
  g_block = lb


  do i=1,3
     do j=1,3
        do k=1,3
           do ip=1,n_unk
              vbles(ip,i,j,k)=vbles10(ip,i+ix-2,j+iy-2,k+iz-2)
           enddo
        enddo
     enddo
  enddo


  Fi=0d0
!!$  
  
  if(mode_1.eq.2)then
     rfactor = dt/dtold
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
    endif
  
  
!!$  
!!$  phi = vbles(1,2,2,2)
!!$  c   = vbles(2,2,2,2)

  phi = g_gamma*unk(4,ix,iy,iz,lb) + (1-g_gamma)* unk(1,ix,iy,iz,lb)
  c   = g_gamma*unk(10,ix,iy,iz,lb) + (1-g_gamma)* unk(7,ix,iy,iz,lb)



  Fi(1) = PhaseRHS(dx,vbles10,ix,iy,iz)


!!$i=2
!!$j=2
!!$k=2
!!$
!!$ phix = 0.5*(vbles(1,i+1,j,k)-vbles(1,i-1,j,k))/dx
!!$ phiy = 0.5*(vbles(1,i,j+1,k)-vbles(1,i,j-1,k))/dx
!!$ phiz = 0.5*(vbles(1,i,j,k+1)-vbles(1,i,j,k-1))/dx
!!$ 
!!$ norm  = 1/g_width!max(1d-6,vbles(1,2,2,2)*(1d0-vbles(1,2,2,2))/g_lambda)
!!$ norm0 = max(1d-6,sqrt(phix**2 + phiy**2+phiz**2))
 
 if (mode_1.eq.1) then
!    phi_t = (unk(1,i,j,k,lb)-unk(2,i,j,k,lb))/dt
 elseif(mode_1.eq.2)then
    rfactor = dt/dtold
    rf4 = rfactor*rfactor/(rfactor+1.0)
    rf5 = (rfactor+1.0)
    rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
!    phi_t = (rf4*unk(5,i,j,k,lb)-rf5*unk(4,i,j,k,lb)+rf6*unk(1,i,j,k,lb))/dt
 endif


if(g_antitrapping)a_trap = (g_AlFe_c(2,1)-g_AlFe_c(1,1)) !!!! = (Delta c)(phi_dot)/grad(phi)!


  Fi(2) = SoluteRHS(dx,c,a_trap,vbles10,ix,iy,iz,mode_1,dt,rf4,rf5,rf6,lb)
 
!!$  a_trap=0d0
!!$  if(g_alpha.ne.0d0)then
!!$     a_trap =g_alpha*sqrt(8d0)*g_lambda*g_ratio*(g_Dch/g_D(1))
!!$  endif







  do ip=1,N_unk
     n=1+(ip-1)*nunkvbles
     if(mode_1.le.1)then
        Fi(ip)=vbles(ip,2,2,2)-dt*Fi(ip) 
     elseif(mode_1.eq.2)then
        Fi(ip)=vbles(ip,2,2,2)-dt*Fi(ip)/rf6
     endif
  enddo




end subroutine get_RHS_Fi


double precision function X_phi(phi)
implicit none
double precision, intent (in):: phi
double precision :: L=1d-9,R=1-1d-9 
X_phi = 0
if(phi.lt.L)then
X_phi = log(L)-log(1-L)
elseif(phi.gt.R)then
X_phi = log(R)-log(1-R)
else
X_phi = log(phi)-log(1-phi)
endif
end function X_phi


double precision function FreeEnergy(c,T,phi,lp,approx)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi
integer, intent(in)::lp
logical, intent(in)::approx
double precision :: a0,a1,c1,c0,x1,x0,dx1dc,dx0dc,dx1dphi,dx0dphi,h,f_sol_0,f_liq_0,f1,f0,gphi,dgphi,df1dx1,df0dx0,gh_phi
FreeEnergy =0d0
      c0 = g_AlFe_c(1,1) !solid common tang
      c1 = g_AlFe_c(2,1) !liquid  common tang
      h =  g_AlFe_c(1,2) ! slope (both curves)
      f_sol_0 = g_AlFe_c(1,3)
      f_liq_0 = g_AlFe_c(2,3)
      a0 = g_AlFe_c(1,4) !solid acceleration
      a1 = g_AlFe_c(2,4) !liquid acceleration
      gphi = gh_phi(phi,0)
!      dgphi = gh_phi(phi,1)
      if(g_Wheeler)then
         x1 = c
         x0 = c
      else
!!$         x1 = ((a0 * c0 - a1 * c1) * gphi + (c - c0) * a0 + a1 * c1) / ((a0 - a1) * gphi + a1)
!!$         x0 = ((a0 * c0 - a1 * c1) * gphi + a1 * c) / ((a0 - a1) * gphi + a1)

         x1 = ((a0 * c0 - a1 * c1) * phi + (c - c0) * a0 + a1 * c1) / ((a0 - a1) * phi + a1)
         x0 = ((a0 * c0 - a1 * c1) * phi + a1 * c) / ((a0 - a1) * phi + a1)
      endif
      f1 = 0.5*a1*(x1-c1)**2+h*(x1-c1)+f_liq_0
      f0 = 0.5*a0*(x0-c0)**2+h*(x0-c0)+f_sol_0
if(lp.eq.0)then
      FreeEnergy = gphi*f1 + (1-gphi)*f0
elseif(lp.eq.1)then
   dgphi = gh_phi(phi,1)
   if(g_Plapp)then !
      FreeEnergy = dgphi*(h-c)*(c1-c0)
      return
   endif
   
   if(g_Wheeler)then
      FreeEnergy = dgphi*f1 -dgphi*f0
   else
!!$         dx1dphi = (-a0 * ((c - c0) * a0 - a1 * (c - c1)) / (gphi * a0 - (-1 + gphi) * a1) ** 2)
!!$         dx0dphi = (-a1 * ((-c + c1) * a1 + (c - c0) * a0) / ((1 - gphi) * a1 + gphi * a0) ** 2)
!!$         df1dx1 = a1*(x1-c1)+h
!!$         df0dx0 = a0*(x0-c0)+h
!!$      
!!$         FreeEnergy = dgphi*f1 -dgphi*f0 + gphi*df1dx1*dx1dphi*dgphi + (1-gphi)*df0dx0*dx0dphi*dgphi





!!$         dx1dphi = (-a0 * ((c - c0) * a0 - a1 * (c - c1)) / (phi * a0 - (-1 + phi) * a1) ** 2)
!!$         dx0dphi = (-a1 * ((-c + c1) * a1 + (c - c0) * a0) / ((1 - phi) * a1 + phi * a0) ** 2)

      dx1dphi = c0-c1
      dx0dphi = c0-c1
      df1dx1 = a1*(x1-c1)+h
      df0dx0 = a0*(x0-c0)+h
      
      FreeEnergy = dgphi*f1 -dgphi*f0 + gphi*df1dx1*dx1dphi + (1-gphi)*df0dx0*dx0dphi
      
      
      !         FreeEnergy = FreeEnergy + 2.6*phi*(-1 + phi)*(-1 + 2*phi)
   endif
elseif(lp.eq.2)then
   if(g_plapp)then
      write(*,*)"no c term in plapp"
      stop
   endif
   df1dx1 = a1*(x1-c1)+h 
   df0dx0 = a0*(x0-c0)+h 
   if(g_Wheeler)then
      FreeEnergy = gphi*df1dx1 + (1-gphi)*df0dx0
   else
!!$            dx1dc = a0 / ((a0 - a1) * gphi + a1)
!!$            dx0dc = a1 / ((a0 - a1) * gphi + a1)
!!$
!!$
!!$            FreeEnergy = gphi*df1dx1*dx1dc + (1-gphi)*df0dx0*dx0dc 
      
      
      dx1dc = a0 / ((a0 - a1) * phi + a1)
      dx0dc = a1 / ((a0 - a1) * phi + a1)
      
      
      FreeEnergy = gphi*df1dx1*dx1dc + (1-gphi)*df0dx0*dx0dc
   endif
else
   write(*,*)"lp>2 not coded"
   stop
endif


end function FreeEnergy

double precision function FreeEnergy_9Feb(c,T,phi,lp,approx)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi
integer, intent(in)::lp
logical, intent(in)::approx
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(2),hb(2),hr(2)
double precision :: Tk,Cp,R,dc,df,dfdc
double precision :: FE,dh=1d-6,phi_,c_
double precision Gibbs_FE_liq,Gibbs_FE_sol,Gibbs_FE_l,Gibbs_FE_e,T_K,g_phi,EntropyMixing,Gibbs_SiGeliq,Gibbs_SiGesol !functions
double precision :: Gibbs_FE_l_AlSi,Gibbs_FE_e_AlSi,cc,Gibbs_Fe_l_PbSn,Gibbs_Fe_e_PbSn,Gibbs_Fe_l_NiCu,Gibbs_Fe_e_NiCu,gh_phi, t_param
double precision :: Gibbs_Fe_La,Gibbs_Fe_a,AlFeFreeEnergy,Cphi_f,U_CurveS,U_CurveL,c_liq,c_sol,s0,s1,s2,ds1,dc_liq,dc_sol,d_c_c_liq,d_c_c_sol
integer :: i,j


!(g_AlFe_LUT,xx0,g_Ncphi,2)
if(approx)then
   
   !t_param = 1d0 - exp(-g_time/g_tinf)
   
!   t_param = min(g_time/g_tinf,1d0)
   
!   g_AlFe_c(1,3) = t_param*(g_f0(1)) + (1d0-t_param)*g_f0(2)
   
!!$if(t_param.gt.0.99)then
!!$write(*,*)g_time,g_f0
!!$write(*,*)g_AlFe_c(:,3)
!!$stop
!!$endif

   
   
   if(lp.eq.0)then
      
      
      
!!$      
!!$      s0 = g_AlFe_c(1, 4)*g_AlFe_c(1, 1) - g_AlFe_c(2, 4)*g_AlFe_c(2, 1) - g_AlFe_c(1, 2) + g_AlFe_c(2, 2)
!!$      s1 = (g_AlFe_c(1, 4) - g_AlFe_c(2, 4))*gh_phi(phi,0) + g_AlFe_c(2, 4)
!!$
!!$      c_sol = (S0*gh_phi(phi,0)+g_AlFe_c(2,4)*c)/s1
!!$      c_liq = (S0*(gh_phi(phi,0)-1)+g_AlFe_c(1,4)*c)/s1

!modification 30 sep lp=0

      s1 = (g_AlFe_c(1, 4)*g_AlFe_c(1, 1) - g_AlFe_c(2, 4)*g_AlFe_c(2, 1))*gh_phi(phi,0) + g_AlFe_c(2, 4)*c
      s0 = (g_AlFe_c(1, 4) - g_AlFe_c(2, 4))*gh_phi(phi,0) + g_AlFe_c(2, 4)
      c_sol = s1/s0
      s2 = (g_AlFe_c(1, 4)*g_AlFe_c(1, 1) - g_AlFe_c(2, 4)*g_AlFe_c(2, 1))*(gh_phi(phi,0)- 1d0) + g_AlFe_c(1, 4)*c
      c_liq = s2/s0

!!!! end of mod 30 sep
      
if(g_Wheeler)then      
            c_liq = c
            c_sol = c
endif
      FreeEnergy_9Feb =   (0.5*g_AlFe_c(2,4)*(c_liq-g_AlFe_c(2,1))**2+g_AlFe_c(2,2)*(c_liq-g_AlFe_c(2,1))+g_AlFe_c(2,3) )*gh_phi(phi,0)&
           +  (0.5*g_AlFe_c(1,4)*(c_sol-g_AlFe_c(1,1))**2+g_AlFe_c(1,2)*(c_sol-g_AlFe_c(1,1))+g_AlFe_c(1,3) )*gh_phi(1-phi,0)
      
   elseif(lp.eq.1)then
!!$      s0 = g_AlFe_c(1, 4)*g_AlFe_c(1, 1) - g_AlFe_c(2, 4)*g_AlFe_c(2, 1) - g_AlFe_c(1, 2) + g_AlFe_c(2, 2)
!!$      s1 = (g_AlFe_c(1, 4) - g_AlFe_c(2, 4))*gh_phi(phi,0) + g_AlFe_c(2, 4)
!!$      c_sol = (S0*gh_phi(phi,0)+g_AlFe_c(2,4)*c)/s1
!!$      c_liq = (S0*(gh_phi(phi,0)-1)+g_AlFe_c(1,4)*c)/s1
     !modification 30 sep lp = 1

      s1 = (g_AlFe_c(1, 4)*g_AlFe_c(1, 1) - g_AlFe_c(2, 4)*g_AlFe_c(2, 1))*gh_phi(phi,0) + g_AlFe_c(2, 4)*c
      s0 = (g_AlFe_c(1, 4) - g_AlFe_c(2, 4))*gh_phi(phi,0) + g_AlFe_c(2, 4)
      c_sol = s1/s0
      s2 = (g_AlFe_c(1, 4)*g_AlFe_c(1, 1) - g_AlFe_c(2, 4)*g_AlFe_c(2, 1))*(gh_phi(phi,0)- 1d0) + g_AlFe_c(1, 4)*c
      c_liq = s2/s0

!!!! end of mod 30 sep

if(g_Wheeler)then      
            c_liq = c
            c_sol = c
endif 
!!$
!!$            c_liq = c
!!$            c_sol = c
      FreeEnergy_9Feb =   (0.5*g_AlFe_c(2,4)*(c_liq-g_AlFe_c(2,1))**2+g_AlFe_c(2,2)*(c_liq-g_AlFe_c(2,1))+g_AlFe_c(2,3) )*gh_phi(phi,0)&
           +  (0.5*g_AlFe_c(1,4)*(c_sol-g_AlFe_c(1,1))**2+g_AlFe_c(1,2)*(c_sol-g_AlFe_c(1,1))+g_AlFe_c(1,3) )*gh_phi(1-phi,0)
      if(phi.gt.0.5)then
         dh=-1d-8
      else
         dh=1d-8
      endif
      phi_ = phi + dh
      
!!$      s0 = g_AlFe_c(1, 4)*g_AlFe_c(1, 1) - g_AlFe_c(2, 4)*g_AlFe_c(2, 1) - g_AlFe_c(1, 2) + g_AlFe_c(2, 2)
!!$      s1 = (g_AlFe_c(1, 4) - g_AlFe_c(2, 4))*gh_phi(phi_,0) + g_AlFe_c(2, 4)
!!$      c_sol = (S0*gh_phi(phi_,0)+g_AlFe_c(2,4)*c)/s1
!!$      c_liq = (S0*(gh_phi(phi_,0)-1)+g_AlFe_c(1,4)*c)/s1
!modification 30 sep lp=1

      s1 = (g_AlFe_c(1, 4)*g_AlFe_c(1, 1) - g_AlFe_c(2, 4)*g_AlFe_c(2, 1))*gh_phi(phi_,0) + g_AlFe_c(2, 4)*c
      s0 = (g_AlFe_c(1, 4) - g_AlFe_c(2, 4))*gh_phi(phi_,0) + g_AlFe_c(2, 4)
      c_sol = s1/s0
      s2 = (g_AlFe_c(1, 4)*g_AlFe_c(1, 1) - g_AlFe_c(2, 4)*g_AlFe_c(2, 1))*(gh_phi(phi_,0)- 1d0) + g_AlFe_c(1, 4)*c
      c_liq = s2/s0

!!!! end of mod 30 sep      


      if(g_Wheeler)then      
         c_liq = c
         c_sol = c
      endif
      FreeEnergy_9Feb = ((0.5*g_AlFe_c(2,4)*(c_liq-g_AlFe_c(2,1))**2+g_AlFe_c(2,2)*(c_liq-g_AlFe_c(2,1))+g_AlFe_c(2,3) )*gh_phi(phi_,0)&
           +  (0.5*g_AlFe_c(1,4)*(c_sol-g_AlFe_c(1,1))**2+g_AlFe_c(1,2)*(c_sol-g_AlFe_c(1,1))+g_AlFe_c(1,3) )*gh_phi(1-phi_,0) - FreeEnergy_9Feb)/dh
      
      
   elseif(lp.eq.2)then
      
!!$      s0 = g_AlFe_c(1, 4)*g_AlFe_c(1, 1) - g_AlFe_c(2, 4)*g_AlFe_c(2, 1) - g_AlFe_c(1, 2) + g_AlFe_c(2, 2)
!!$      s1 = (g_AlFe_c(1, 4) - g_AlFe_c(2, 4))*gh_phi(phi,0) + g_AlFe_c(2, 4)
!!$      c_sol = (S0*gh_phi(phi,0)+g_AlFe_c(2,4)*c)/s1
!!$      c_liq = (S0*(gh_phi(phi,0)-1)+g_AlFe_c(1,4)*c)/s1
      
!modification 30 sep lp =2

      s1 = (g_AlFe_c(1, 4)*g_AlFe_c(1, 1) - g_AlFe_c(2, 4)*g_AlFe_c(2, 1))*gh_phi(phi,0) + g_AlFe_c(2, 4)*c
      s0 = (g_AlFe_c(1, 4) - g_AlFe_c(2, 4))*gh_phi(phi,0) + g_AlFe_c(2, 4)
      c_sol = s1/s0
      s2 = (g_AlFe_c(1, 4)*g_AlFe_c(1, 1) - g_AlFe_c(2, 4)*g_AlFe_c(2, 1))*(gh_phi(phi,0)- 1d0) + g_AlFe_c(1, 4)*c
      c_liq = s2/s0

!!!! end of mod 30 sep      

      if(g_Wheeler)then      
         c_liq = c
         c_sol = c
      endif
!!$
!!$            c_liq = c
!!$            c_sol = c
      
      FreeEnergy_9Feb =   (0.5*g_AlFe_c(2,4)*(c_liq-g_AlFe_c(2,1))**2+g_AlFe_c(2,2)*(c_liq-g_AlFe_c(2,1))+g_AlFe_c(2,3) )*gh_phi(phi,0)&
           +  (0.5*g_AlFe_c(1,4)*(c_sol-g_AlFe_c(1,1))**2+g_AlFe_c(1,2)*(c_sol-g_AlFe_c(1,1))+g_AlFe_c(1,3) )*gh_phi(1-phi,0)
      if(c.gt.0.5)then
         dh=-1d-8
      else
         dh=1d-8
      endif
      c_ = c + dh
      
!!$      s0 = g_AlFe_c(1, 4)*g_AlFe_c(1, 1) - g_AlFe_c(2, 4)*g_AlFe_c(2, 1) - g_AlFe_c(1, 2) + g_AlFe_c(2, 2)
!!$      s1 = (g_AlFe_c(1, 4) - g_AlFe_c(2, 4))*gh_phi(phi,0) + g_AlFe_c(2, 4)
!!$      c_sol = (S0*gh_phi(phi,0)+g_AlFe_c(2,4)*c_)/s1
!!$      c_liq = (S0*(gh_phi(phi,0)-1)+g_AlFe_c(1,4)*c_)/s1
      
!modification 30 sep lp =2

      s1 = (g_AlFe_c(1, 4)*g_AlFe_c(1, 1) - g_AlFe_c(2, 4)*g_AlFe_c(2, 1))*gh_phi(phi,0) + g_AlFe_c(2, 4)*c_
      s0 = (g_AlFe_c(1, 4) - g_AlFe_c(2, 4))*gh_phi(phi,0) + g_AlFe_c(2, 4)
      c_sol = s1/s0
      s2 = (g_AlFe_c(1, 4)*g_AlFe_c(1, 1) - g_AlFe_c(2, 4)*g_AlFe_c(2, 1))*(gh_phi(phi,0)- 1d0) + g_AlFe_c(1, 4)*c_
      c_liq = s2/s0

!!!! end of mod 30 sep      

      if(g_Wheeler)then      
         c_liq = c_
         c_sol = c_
      endif

!!$            c_liq = c_
!!$            c_sol = c_
      FreeEnergy_9Feb = ((0.5*g_AlFe_c(2,4)*(c_liq-g_AlFe_c(2,1))**2+g_AlFe_c(2,2)*(c_liq-g_AlFe_c(2,1))+g_AlFe_c(2,3) )*gh_phi(phi,0)&
           +  (0.5*g_AlFe_c(1,4)*(c_sol-g_AlFe_c(1,1))**2+g_AlFe_c(1,2)*(c_sol-g_AlFe_c(1,1))+g_AlFe_c(1,3) )*gh_phi(1-phi,0) - FreeEnergy_9Feb)/dh
      
      
   else
      write(*,*)"lp = ",lp
      stop
   endif
else
   !thermodynamic data base   
   FreeEnergy_9Feb =    AlFeFreeEnergy(c,T,phi,lp)
   write(*,*)".true in free energy"
   stop
endif !end approx or thermo base
return

end function FreeEnergy_9Feb

double precision function FreeEnergyWBM(c,T,phi,lp,approx)
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
double precision :: Gibbs_Fe_La,Gibbs_Fe_a,AlFeFreeEnergy
integer :: i,j
if(approx)then
   if(lp.eq.0)then
      FreeEnergyWBM =   (0.5*g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))**2+g_AlFe_c(2,2)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,3) )*g_phi(phi,0)&
           +  (0.5*g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))**2+g_AlFe_c(1,2)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,3) )*g_phi(1-phi,0)
   elseif(lp.eq.1)then
      FreeEnergyWBM =   (0.5*g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))**2+g_AlFe_c(2,2)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,3) )*g_phi(phi,1)&
           -  (0.5*g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))**2+g_AlFe_c(1,2)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,3) )*g_phi(phi,1)
   elseif(lp.eq.2)then
      FreeEnergyWBM =   (g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,2) )*g_phi(phi,0)&
           +  (g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,2) )*g_phi(1-phi,0)
   elseif(lp.eq.-2)then !f_cc
      FreeEnergyWBM = (1 - g_phi(phi,0))*g_AlFe_c(1,4) + g_phi(phi,0)*g_AlFe_c(2,4)
      
   elseif(lp.eq.-1)then !f_c_phi
      FreeEnergyWBM = (g_AlFe_c(2,4)*(c-g_AlFe_c(2,1))+g_AlFe_c(2,2) )*g_phi(phi,1)&
           -  (g_AlFe_c(1,4)*(c-g_AlFe_c(1,1))+g_AlFe_c(1,2) )*g_phi(phi,1)
      
   else
      write(*,*)"lp = ",lp
      stop
   endif
   
else
   
   FreeEnergyWBM =    AlFeFreeEnergy(c,T,phi,lp)
endif
return

end function FreeEnergyWBM

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
   AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,T,0,0)*g_phi(phi,0)
   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e2_AlFe(c,T,0,0)*(1-g_phi(phi,0))
!   return
elseif(lp.le.1)then
   if(g_Tflag)then
      write(*,*)g_Tflag
      stop
   else
      AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,T,0,0)*g_phi(phi,1)
      AlFeFreeEnergy = AlFeFreeEnergy -   Gibbs_FE_e2_AlFe(c,T,0,0)*(g_phi(phi,1))
 !     return
   endif
elseif(lp.eq.2)then !Solute


   AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,T,1,0)*g_phi(phi,0)
   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e2_AlFe(c,T,1,0)*(1-g_phi(phi,0))

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

AlFeFreeEnergy = AlFeFreeEnergy/(g_R*g_T0)

end function

!
double precision function AlFey1(c,d)
integer, intent(in)::d
double precision, intent(in) :: c



if(d.eq.0)then
      AlFey1=(0.2350/c - 0.8625)/0.1375 
      AlFey1 = min(0.9999,max(0.0001,AlFey1))
      return
elseif(d.eq.1)then
if(c.gt.0.27246.or.c.lt.0.2350)then
   write(*,*)"c out of bound in AlFey1()"
   stop
endif
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

double precision function g_phi_8Jul2020(phi,d)
use solution_parameters
  double precision, intent (in):: phi
  integer, intent (in):: d
  double precision :: n,h_phi,u_phi,u,h
  g_phi_8Jul2020=0d0
  if(phi.lt.1d-16)return
  if(g_skew_potential)then
     if(d.eq.1)then
        g_phi_8Jul2020 =  ((-1/(2-phi)/2-1/phi/2)*(phi**2-phi**3/3)+(log(2-phi)/2-log(phi)/2)*(-phi**2+2*phi)-1.D0/3.D0+phi/3+2.D0/3.D0/(2-phi))&
             /(2.D0/3.D0*log(2.D0)-1.D0/6.D0)
     elseif(d.eq.2)then
        g_phi_8Jul2020 = ((-1/((2-phi)**2)/2+1/phi**2/2)*(phi**2-phi**3/3)+2*(-1/(2-phi)/2-1/phi/2)*(-phi**2+2*phi)&
             +(log(2-phi)/2-log(phi)/2)*(-2*phi+2)+1.D0/3.D0+2.D0/3.D0/(2-phi)**2)/(2.D0/3.D0*log(2.D0)-1.D0/6.D0)
     else !assume 0
        g_phi_8Jul2020 = ((log(2-phi)/2-log(phi)/2)*(phi**2-phi**3/3)-phi/3&
             +phi**2/6-2.D0/3.D0*log(2-phi)+2.D0/3.D0*log(2.D0))/(2.D0/3.D0*log(2.D0)-1.D0/6.D0)
     endif
  elseif(g_power_skew)then
     n = g_power 
     if(d.eq.1)then
        g_phi_8Jul2020 = 0.25*phi*(n+4)*(n+2)/(1-phi)*(1-phi)**(0.5*n+1) 
     elseif(d.eq.2)then
        g_phi_8Jul2020 = -0.125*(n*phi+2*phi-2)*(n+4)*(n+2)/(1-phi)**2*(1-phi)**(0.5*n+1)
     else !assume 0
        g_phi_8Jul2020 = 1-0.5*((2+n)*phi+2)*(1-phi)**(0.5*n+1)
     endif
!
  elseif(g_s_skew)then
     if(d.eq.0)then

!AMM
        if(g_AMM.eq.1)then
           g_phi_8Jul2020 = -0.1485934166D1*phi**2+0.3268185149D1*phi-0.3258617531D0*log(2544018011.D0+0.2551487495D11*phi)+0.7057191429D1
        elseif(g_AMM.eq.2)then
           g_phi_8Jul2020 = -0.2D1*phi**3+0.3D1*phi**2
        else
           write(*,*)"g_AMM is 1 or 2"
        endif


     elseif(d.eq.1)then 
!AMM

        if(g_AMM.eq.1)then
           g_phi_8Jul2020 = -0.1516536976D3*phi*(-1.D0+phi)/(0.5088036022D1+0.510297499D2*phi)
        elseif(g_AMM.eq.2)then
           g_phi_8Jul2020 = -0.6D1*phi*(-1.D0+phi)
        else
           write(*,*)"g_AMM is 1 or 2"
        endif



     endif

  else
     if(d.eq.1)then
        g_phi_8Jul2020 = 6.*phi*(1.-phi)
        !    g_phi_8Jul2020 = 30.*phi**4-60.*phi**3+30.*phi**2
     elseif(d.eq.2)then
        g_phi_8Jul2020 = 6.-12.*phi
        !    g_phi_8Jul2020 = 120.*phi**3-180.*phi**2+60.*phi
     else !assume 0
        g_phi_8Jul2020 = 3.*phi**2-2.*phi**3
        !    g_phi_8Jul2020 = 6.*phi**5-15.*phi**4+10.*phi**3
     endif
  endif
end function g_phi_8Jul2020


logical function Not_number(x)
double precision, intent(in):: x
if(x.eq.x)then
Not_number=.false.
else
Not_number=.true.
endif
end function Not_number
!
double precision function f_c(x,df)
double precision, intent (in)::x
integer, intent(in)::df
f_c=1.
if(df.eq.1)f_c=0.
return
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
