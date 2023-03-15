!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Warren Boetinger Loginova model    !!
!!Coded by Peter Bollada 17 Sep 2015 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!T_K converts non-dimensional parameter T to dimensioned T_K(T)
double precision function T_K(T)
implicit none
double precision, intent (in)::T
double precision :: T_eut= 456.14
  
 T_eut = 456.14
 
 T_K = T_eut*(1.+T)

 end function T_K
 double precision function TemperatureRHS(vbles,dx)
implicit none
double precision, intent (in):: dx,vbles(3,3,3)

  TemperatureRHS=0.
  return
 

 end function TemperatureRHS

!

double precision function SoluteRHS(vbles,dx)
use solution_parameters
implicit none
double precision, intent (in):: vbles(3,3,3),dx


SoluteRHS=0.
return


end function SoluteRHS

double precision function potential(phi,c)
use solution_parameters
implicit none

double precision, intent (in)::phi,c

potential=0.
return

end function potential


double precision function PhaseRHS(vbles,dx)
use solution_parameters
implicit none

double precision, intent(in)::dx,vbles(3,3,3)

 PhaseRHS=0.
return

end function PhaseRHS





subroutine get_RHS_Fi_2d(i,j,k,lb,dx,Fi)
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none

integer, intent(in)::i,j,k,lb
double precision, intent(in)::dx
double precision, intent(out):: Fi(3)
double precision :: c,T,phi,c_dot,phi_dot
double precision phi_rhs,rfactor,rf4,rf5,rf6,S
double precision TemperatureRHS,SoluteRHS,T_k,phi_
double precision vbles(3,3,3),unk_(18)
logical Not_number
integer ip,n,ii,jj,kk




  Fi=0.

  do ii=1,3
     do jj=1,3
        do ip=1,3
           vbles(ip,ii,jj)=unk(1+(ip-1)*nunkvbles,i+ii-2,j+jj-2,k,lb)
        enddo
     enddo
  enddo




  if(neigh(1,2,lb).eq.bc_cond_sym.and. i.eq.9)then
     ii = 3
     do jj=1,3
        ip=2
        vbles(ip,ii,jj)=delta
     enddo
  endif


  do ii=1,18
     unk_(ii) = unk(ii,i,j,k,lb)
  enddo

  
  call get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_)
  do ip=1,3
     unk(1+(ip-1)*nunkvbles,i,j,k,lb)=vbles(ip,2,2)
  enddo



 
end subroutine get_RHS_Fi_2d

subroutine get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_)
use solution_parameters
  implicit none


integer, intent(in)::mode_1,nunkvbles,lnblocks
double precision, intent(in)::dx,dt,dtold,unk_(18)
double precision, intent(out):: Fi(3)
double precision, intent(in):: vbles(3,3,3)
double precision :: c,T,phi,c_dot,phi_dot
double precision phi_rhs,dphi_rhs,rfactor,rf4,rf5,rf6,S
double precision TemperatureRHS,SoluteRHS,T_k,phi_
logical Not_number
integer ip,n,ii,jj,kk

Fi=0d0
   
  if (mode_1.eq.1) then
     


     phi_dot = (vbles(1,2,2)-unk_(2))/dt

     n=(2)*nunkvbles
     c_dot = (vbles(3,2,2)-unk_(2+n))/dt
  elseif(mode_1.eq.2)then
     rfactor = dt/dtold


     rf4 = rfactor*rfactor/(rfactor+1.0)
     rf5 = (rfactor+1.0)
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
     
     phi_dot = (rf4*unk_(5)-rf5*unk_(4)+rf6*vbles(1,2,2))/dt

     n=2*nunkvbles
     c_dot = (rf4*unk_(5+n)-rf5*unk_(4+n)+rf6*vbles(3,2,2))/dt
  else
     write(*,*)"mode_1 not = 1 or 2"
     stop
  endif

  phi = vbles(1,2,2)
  T   = vbles(2,2,2)
  c   = vbles(3,2,2)

  Fi(1)=0d0
  Fi(2) = TemperatureRHS(vbles,dx)
  Fi(3) = SoluteRHS(vbles,dx)

  do ip=1,3
     n=(ip-1)*nunkvbles
     if(mode_1.le.1)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)-unk_(n+2)
     else if(mode_1.eq.2)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)/rf6-unk_(n+2)
     end if
  enddo


end subroutine get_RHS_Fi




double precision function FreeEnergy(c,T,phi)
use solution_parameters
implicit none


double precision, intent (in):: T,c,phi
FreeEnergy=0d0
return
end function FreeEnergy

double precision function EntropyMixing(c,d)
double precision, intent (in)::c
integer, intent (in)::d
 
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
if(x.gt.1e20)then
Not_number=.true.
elseif(x.lt.-1e20)then
Not_number=.true.
elseif(x.le.1e20.or.x.ge.-1e20)then
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
double precision function f_cc(x,df)
double precision, intent (in)::x
integer, intent(in)::df
if(x*(1.-x).le.0.)then
f_cc=0.
else
 if(df.eq.0)then
  f_cc=x*(1.-x)*(log(x)-log(1.-x))
 elseif(df.eq.1)then
  f_cc=log(x)-log(1.-x)-2.*x*log(x)+2.*x*log(1.-x)+1.
 endif
endif
end function f_cc


