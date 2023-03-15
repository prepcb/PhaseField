!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  Utilities for phase field stencils in 2-d and 3-d
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




double precision function GradientEnergy(LapPhi,vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::LapPhi(3),dx,vbles(5,3,3)
integer, intent (in):: lp
double precision ::phi_x(3),phi_y(3),phi(3),div_q
integer :: ii

!isotropic                                                                                                

  phi(1) = vbles(1,2,2)
  phi(2) = vbles(2,2,2)
  phi(3) = vbles(3,2,2)

  phi_x(1) =  (vbles(1,3,2)-vbles(1,1,2))*0.5/dx
  phi_x(2) =  (vbles(2,3,2)-vbles(2,1,2))*0.5/dx
  phi_x(3) =  (vbles(3,3,2)-vbles(3,1,2))*0.5/dx

  phi_y(1) =  (vbles(1,2,3)-vbles(1,2,1))*0.5/dx
  phi_y(2) =  (vbles(2,2,3)-vbles(2,2,1))*0.5/dx
  phi_y(3) =  (vbles(3,2,3)-vbles(3,2,1))*0.5/dx
  
  GradientEnergy = 0d0
  do ii=1,N_phases
     if(ii.ne.lp)then
        GradientEnergy = GradientEnergy - eps22(lp,ii)*(2.*((phi(lp)*phi_x(ii)-phi(ii)*phi_x(lp))*phi_x(ii)&
             + (phi(lp)*phi_y(ii)-phi(ii)*phi_y(lp))*phi_y(ii)) -&
             (phi(ii)*LapPhi(lp)-phi(lp)*LapPhi(ii))*phi(ii) )
     endif
  enddo
  
  
end function GradientEnergy










!T_K converts non-dimensional parameter T to dimensioned T_K(T)
double precision function T_K(T)
implicit none
double precision, intent (in)::T
double precision :: T_eut= 456.14
!Alternate non-dimensionalisation
! double precision T_b,C_b,L_b, 
!  T_b = 505.07
!  C_b = 210.
!  L_b = 7029.
!  T_K = L_b/C_b*T+T_b
  
 T_eut = 456.14
 
 T_K = T_eut*(1.+T)

 end function T_K
 double precision function TemperatureRHS(vbles,dx,c,T,phi,c_dot,phi_dot)
!This is not the most efficient routine. It takes Div (D*c*(1-c)Grad(h)) by
!examining its neighboring values. Quicker is to evaluate all points in
!one sweep then use them for this routine. But that approach would have to be !done outside this routine.
implicit none
double precision, intent (in):: dx,c,T,phi(3),c_dot,phi_dot(3),vbles(5,3,3)
!integer, intent (in)::df
double precision, dimension(3,3)::stencil =reshape((/0.,1.,0.,1.,-4.,1.,0.,1.,0./),(/3,3/))
double precision, dimension(3)::stencil_x =(/-.5,.0,.5/)
double precision:: FreeEnergy,x,a,xa,dFdc(3,3) !x coordinate
double precision:: Cp= 28.5!29.9758394     !J/mol/K
! double precision:: D=2.110795e-5      !m^2/s
double precision:: T_eut=456.14  !K
! double precision ::D_c = 2.976894361E-9 !m^2/s characteristic diffusivity
! double precision ::D_ch = 1.126392460619016E-011 !hunt Jackson
double precision ::D_ch = 6.7d-10,R = 8.31 !J/(K mol) !hunt Jackson
double precision, dimension(3):: D1=(/6.7e-10,6.7e-12,6.7e-12/)!Jackson Hunt
! double precision :: D_tilde = 7090.594237 !dimensionless diffusivity
double precision :: D_tilde = 2630.6602369, D_heat_c !dimensionless debug Le=10.
double precision :: T_,c_,phi_(3)
logical :: heat_c_term=.false.
integer :: ii,jj,kk, N_phases=3
double precision FF(3)

!  rhoa(1)=10660. !liquid Lead kg/m^3
!  rhoa(2)=11340. !Solid Lead
!  rhoa(3)=11340. !The same (for want of a better option)
!  rhob(1)=6999.  !liquid Tin
!  rhob(2)=7365.  !Solid Tin (White)
!  rhob(3)=7365.  !The same

!  rho=0.25*(11340.+ 7365.+6999.+10660.)!kg/m^3 keep it simple
!  kappa =35.3 !J/(s m)
  TemperatureRHS=0.


  D_tilde = 100.
  !    return !delete for isothermal
  
  TemperatureRHS = FreeEnergy(c,T,phi,c_dot,phi_dot,N_phases+1) !J K/mol
  !  TemperatureRHS = TemperatureRHS/Cp !means LHS is T_dot only
  ! heat capacity is in FreeEnergy function now Cp = -T(d/dT)dF/dT
  TemperatureRHS = TemperatureRHS/T_eut !non dimensional T
  
  do ii=1,3
     do jj=1,3
        !  T_=unk(1+N_phases*nunkvbles,i-2+ii,j-2+jj,k,lb)*stencil(ii,jj)
        T_=vbles(4,ii,jj)*stencil(ii,jj)
        TemperatureRHS = TemperatureRHS + D_tilde*T_/(dx*dx)
     enddo
  enddo
  
  ! replaces c_dot_term an heat_flux_term which does not appear to operate well
  if(heat_c_term)then
     do ii=1,3
        do jj=1,3
           !   T_=unk(1+N_phases*nunkvbles,i-2+ii,j-2+jj,k,lb)
           T_=vbles(4,ii,jj)
           !   c_=unk(1+(N_phases+1)*nunkvbles,i-2+ii,j-2+jj,k,lb)
           c_=vbles(5,ii,jj)
           do kk=1,N_phases
              !    phi_(kk)=unk(1+(kk-1)*nunkvbles,i-2+ii,j-2+jj,k,lb)
              phi_(kk)=vbles(kk,ii,jj)
           enddo
           dFdc(ii,jj)=(1.+T_)*(log(c_)-log(1.-c_))+FreeEnergy(c_,T_,Phi_,c_dot,phi_dot,N_phases+2)
        enddo
     enddo
     D_heat_c = 0.
     do ii=1,N_phases
        D_heat_c = D_heat_c + D1(ii)*phi(ii)
     enddo
     D_heat_c = (c*(1.-c)*D_heat_c/D_ch)/(Cp*T_eut)
     
     !  D_heat_c* grad(dfdc) dot grad(dFdc)
     TemperatureRHS = TemperatureRHS + D_heat_c*((dFdc(3,2)-dFdc(1,2))**2 + (dFdc(2,3)-dFdc(2,1))**2)/(4.*dx*dx)
     
  endif
  !  TemperatureRHS = TemperatureRHS + v_mesh*(unk(1+3*nunkvbles,i+1,j,k,lb)-unk(1+3*nunkvbles,i,j,k,lb))/(dx)
  
 

 end function TemperatureRHS

!

double precision function SoluteRHS(vbles,dx,c_dot,phi_dot,c_c)
! this uses volume element method. This guarantees conservation of solute. The main source of solute is the variation of free energy with phase: d (d F/dc) = ...+ (d^2 F/dc d phi) d phi + ..., which is revealed by a plot of c_dot.
use solution_parameters
implicit none
double precision, intent (in):: dx,c_dot,c_c,vbles(5,3,3)
integer ii,jj,ip,ih,jv
double precision c,T,Phi(3),FreeEnergy,D,phi_dot(3)
double precision, dimension(3,3)::stencil11=reshape((/0.,6.,0.,6.,-24.,6.,0.,6.,0./),(/3,3/))
double precision, dimension(3)::stencil1=(/-.5,0.,.5/)
double precision, dimension(3):: D1=(/6.7e-10,6.7e-12,6.7e-12/)!Jackson Hunt
double precision :: potential,D_, D_c,D_cc,f_c,f_cc
double precision ::D_ch =  6.7d-10! (m^2/s)    characteristic diffusivity delta2/tau=M*eps2(1)


SoluteRHS=0.
!
!if(df.eq.0)then
!x-coordinate transform values
!! Right
 ih=3
 jv=2
 include 'c_stencil.h'
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Left
 ih=1
 jv=2
 include 'c_stencil.h'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!up
 ih=2
 jv=3
 include 'c_stencil.h'
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !down
 ih=2
 jv=1
 include 'c_stencil.h'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SoluteRHS = SoluteRHS/D_ch/(dx*dx) 
! SoluteRHS = SoluteRHS+v_mesh*(unk(1+4*nunkvbles,i+1,j,k,lb)-unk(1+4*nunkvbles,i,j,k,lb))/(dx)
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

return

end function SoluteRHS

double precision function potential(phi,c,lp)
use solution_parameters
implicit none

double precision, intent (in)::phi(3),c
integer, intent(in):: lp
integer :: ip,jp



potential=0.

if(lp.le.N_phases)then 

   do ip=1,N_phases
      potential = potential + eps22(lp,ip)*2d0*phi(lp)*phi(ip)**2 
   enddo



endif
return


end function potential


double precision function PhaseRHS(vbles,dx,lp,Phi,c,T,c_dot,phi_dot)
use solution_parameters
implicit none
!integer, intent(in)::i,j,k,lb
double precision, intent(in)::dx,c,T,c_dot,phi_dot(3),Phi(3),vbles(5,3,3)
integer, intent(in):: lp
double precision :: LapPhi(3)
double precision GradientEnergy,potential,FreeEnergy
double precision rf4,rf5,rf6,rfactor
double precision ::M_a= 0.948190 ,M_b= 0.5680384!m^3/(J s) mobility
double precision ::t_a = 600.61,t_b = 505.07
integer ii,jj
double precision ::W_a(3,3)
double precision ::M_tilde !s characteristic time
double precision :: delta2 = 4E-18!decrease to sharpen. This i.c. must increase by similar factor and 
double precision Cross_term
double precision :: damping=1d0
integer ip,n

!!!!! Calculate characteristic time and space (squared)

! tau_a=1./(M_a*W_a)
 PhaseRHS=0.
! if(lp.eq.3)return
! w_a(1,1)=7e7 ! see notes = 0.07 / delta = 0.07/1e-9 
 w_a(1,1)=7.06d7*0.5094748161 ! = (3/sqrt(2))(sigma_Pb)/delta = (3/sqrt(2))(33.3 ergs/cm^2)/(2e-9m)
 
 call Calc_GradPhi(LapPhi,vbles,dx)
 !assemble and non-dimensionlise phaseRHS 
 
 M_tilde = 1.0-0.159*c ! see notes Anisotropic2d3d.tex
 
 phaseRHS=M_tilde*(GradientEnergy(LapPhi,vbles,lp,dx)-potential(Phi,c,lp)*scale_factor**2/lambda**2&
      -scale_factor**2/lambda*FreeEnergy(c,T,Phi,c_dot,phi_dot,lp)/W_a(1,1))


end function PhaseRHS



subroutine Project_Phase(vbles,dx,phi_rhs,Phi,c,T,c_dot,phi_dot)
use solution_parameters
implicit none
!integer, intent(in)::df
double precision, intent(in)::dx,Phi(3),c,T,c_dot,phi_dot(3),vbles(5,3,3)
double precision, intent(out) :: phi_rhs(3)
double precision PhaseRHS,Wk(3)
integer ip,jp,kp,lp,izero
double precision, dimension(3,3)::I2 = reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))
double precision, dimension(3,3)::Projector = 0.,M=reshape((/1.,1.,1.,1.,1.,1.,1.,1.,1./),(/3,3/))



  phi_rhs=0.

  !assemble Projector
  izero=0
  do ip=1,N_phases
     if(1d0-Phi(ip).lt.1d-16)izero=ip
  enddo
  
  if(izero.eq.0)then
     
     ! maple generated A
     if(.false.)then
   Projector(1,1) = MOBILITY(1,2) * phi(1) * phi(2) / (1 - phi(1)) / (1 - phi(2))&
      + MOBILITY(1,3) * phi(1) * phi(3) / (1 - phi(1)) / (1 - phi(3))
      Projector(2,1) = -MOBILITY(1,2) * phi(1) * phi(2) / (1 - phi(1)) / (1 - phi(2)&
     )
      Projector(3,1) = -MOBILITY(1,3) * phi(1) * phi(3) / (1 - phi(1)) / (1 - phi(3)&
     )
      Projector(1,2) = -MOBILITY(1,2) * phi(1) * phi(2) / (1 - phi(1)) / (1 - phi(2)&
           )
      Projector(2,2) = MOBILITY(1,2) * phi(1) * phi(2) / (1 - phi(1)) / (1 - phi(2))&
      + MOBILITY(2,3) * phi(2) * phi(3) / (1 - phi(2)) / (1 - phi(3))
      Projector(3,2) = -MOBILITY(2,3) * phi(2) * phi(3) / (1 - phi(2)) / (1 - phi(3)&
     )
      Projector(1,3) = -MOBILITY(1,3) * phi(1) * phi(3) / (1 - phi(1)) / (1 - phi(3)&
     )
      Projector(2,3) = -MOBILITY(2,3) * phi(2) * phi(3) / (1 - phi(2)) / (1 - phi(3)&
    )
      Projector(3,3) = MOBILITY(1,3) * phi(1) * phi(3) / (1 - phi(1)) / (1 - phi(3))&
     + MOBILITY(2,3) * phi(2) * phi(3) / (1 - phi(2)) / (1 - phi(3))
   end if

   Projector(1,1) = 1 + phi(1) ** 2 / (phi(1) + phi(3)) ** 2 + phi(1) ** 2 / (phi(1) + phi(2)) ** 2
   Projector(1,2) = -phi(2) / (phi(2) + phi(3)) - phi(1) / (phi(1) + phi(3)) + phi(1) / (phi(1) + phi(2)) ** 2 * phi(2)
   Projector(1,3) = -phi(3) / (phi(2) + phi(3)) + phi(1) / (phi(1) + phi(3)) ** 2 * phi(3) - phi(1) / (phi(1) + phi(2))
   Projector(2,1) = -phi(2) / (phi(2) + phi(3)) - phi(1) / (phi(1) + phi(3)) + phi(1) / (phi(1) + phi(2)) ** 2 * phi(2)
   Projector(2,2) = phi(2) ** 2 / (phi(2) + phi(3)) ** 2 + 1 + phi(2) ** 2 / (phi(1) + phi(2)) ** 2
   Projector(2,3) = phi(2) / (phi(2) + phi(3)) ** 2 * phi(3) - phi(3) / (phi(1) + phi(3)) - phi(2) / (phi(1) + phi(2))
   Projector(3,1) = -phi(3) / (phi(2) + phi(3)) + phi(1) / (phi(1) + phi(3))** 2 * phi(3) - phi(1) / (phi(1) + phi(2))
   Projector(3,2) = phi(2) / (phi(2) + phi(3)) ** 2 * phi(3) - phi(3) / (phi(1) + phi(3)) - phi(2) / (phi(1) + phi(2))
   Projector(3,3) = phi(3) ** 2 / (phi(2) + phi(3)) ** 2 + phi(3) ** 2 / (phi(1) + phi(3)) ** 2 + 1




     
  else
!     if(.false.)then
        Projector(1,1) = MOBILITY(1,2)+MOBILITY(1,3)
        Projector(2,2) = MOBILITY(2,1)+MOBILITY(2,3)
        Projector(3,3) = MOBILITY(3,2)+MOBILITY(3,1)
        Projector(1,2) = -MOBILITY(1,2)
        Projector(1,3) = -MOBILITY(1,3)
        Projector(2,3) = -MOBILITY(2,3)
        Projector(2,1) = -MOBILITY(1,2)
        Projector(3,1) = -MOBILITY(1,3)
        Projector(3,2) = -MOBILITY(2,3)
        
        Projector = Projector/2.0
!     end if
!     Projector=0.

  endif
  
  
  
  do ip=1,N_phases
     Wk(ip) =  PhaseRHS(vbles,dx,ip,Phi,c,T,c_dot,phi_dot) 
  enddo
  
  
  !Project phi_vector
  do ip=1,N_phases
     phi_rhs(ip)=0.!q_factor(ip)*tip*dt
     do jp=1,N_phases
        phi_rhs(ip)=phi_rhs(ip)+Projector(ip,jp)*Wk(jp)
     enddo
  enddo

  
end subroutine Project_Phase





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
double precision, intent(out):: Fi(5)
double precision :: c,T,phi(3),c_dot,phi_dot(3)
double precision phi_rhs(3),rfactor,rf4,rf5,rf6,S
double precision TemperatureRHS,SoluteRHS,T_k,phi_
double precision vbles(5,3,3),unk_(30)
logical Not_number
integer ip,n,ii,jj,kk
  Fi=0.

  do ii=1,3
     do jj=1,3
        do ip=1,5
           vbles(ip,ii,jj)=unk(1+(ip-1)*nunkvbles,i+ii-2,j+jj-2,k,lb)
        enddo
     enddo
  enddo




  if(neigh(1,2,lb).eq.bc_cond_sym.and. i.eq.9)then
     ii = 3
     do jj=1,3
        ip=4
!        vbles(ip,ii,jj)=delta  !derichlet
     enddo
  endif

!  if(neigh(1,1,lb).eq.bc_cond_sym.and. i.eq.2)then
!     ii = 1
!     do jj=1,3
!        ip=4
!        vbles(ip,ii,jj)=delta
!     enddo
!  endif

  do ii=1,30
     unk_(ii) = unk(ii,i,j,k,lb)
  enddo

  
  call get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_)
  do ip=1,5
     unk(1+(ip-1)*nunkvbles,i,j,k,lb)=vbles(ip,2,2)
  enddo



 
end subroutine get_RHS_Fi_2d

subroutine get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_)
use solution_parameters
  implicit none


integer, intent(in)::mode_1,nunkvbles,lnblocks
double precision, intent(in)::dx,dt,dtold,unk_(30)
double precision, intent(out):: Fi(5)
double precision, intent(in):: vbles(5,3,3)
double precision :: c,T,phi(3),c_dot,phi_dot(3)
double precision phi_rhs(3),dphi_rhs(3),rfactor,rf4,rf5,rf6,S
double precision TemperatureRHS,SoluteRHS,T_k,phi_
logical Not_number
integer ip,n,ii,jj,kk

Fi=0d0
   
  if (mode_1.eq.1) then
     
     do ip=1,N_phases
        n=(ip-1)*nunkvbles
        phi_dot(ip) = (vbles(ip,2,2)-unk_(n+2))/dt
     enddo
     n=(N_phases+1)*nunkvbles
     c_dot = (vbles(5,2,2)-unk_(2+n))/dt
  elseif(mode_1.eq.2)then
     rfactor = dt/dtold

!     rf1 = (rfactor+1.0)*dt/(2.0*rfactor+1.0)
!     rf2 = (rfactor+1.0)*(rfactor+1.0)/(2.0*rfactor+1.0)
!     rf3 = rfactor*rfactor/(2.0*rfactor+1.0)

     rf4 = rfactor*rfactor/(rfactor+1.0)
     rf5 = (rfactor+1.0)
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
     
     do ip=1,N_phases
        n=(ip-1)*nunkvbles
        phi_dot(ip) = (rf4*unk_(5+n)-rf5*unk_(4+n)+rf6*vbles(ip,2,2))/dt
     enddo
     n=(N_phases+1)*nunkvbles
     c_dot = (rf4*unk_(5+n)-rf5*unk_(4+n)+rf6*vbles(5,2,2))/dt
  else
     write(*,*)"mode_1 not = 1 or 2"
     stop
  endif

  do ip=1,N_phases
     phi(ip)=vbles(ip,2,2)
  enddo
  
  
  T = vbles(4,2,2)!debug pcb
  c = vbles(5,2,2)

  
  call Project_Phase( vbles,dx,phi_rhs,phi,c,T,c_dot,phi_dot)
  do ip=1,N_phases
     Fi(ip)=phi_rhs(ip)
  enddo
  
  Fi(N_phases+1) = TemperatureRHS(vbles,dx,c,T,phi,c_dot,phi_dot)
  Fi(N_phases+2) = SoluteRHS(vbles,dx,c_dot,phi_dot,c)

  do ip=1,n_phases+2
     n=(ip-1)*nunkvbles
     if(mode_1.le.1)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)-unk_(n+2)
     else if(mode_1.eq.2)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)/rf6-unk_(n+2)
     end if
  enddo


end subroutine get_RHS_Fi


! for now only calculates LapPhi
subroutine Calc_GradPhi(LapPhi,vbles,dx)

implicit none
double precision, intent(in)::dx,vbles(5,3,3)
double precision, intent(out)::LapPhi(3)
double precision, dimension(3,3)::stencil=reshape((/1.,4.,1.,4.,-20.,4.,1.,4.,1./),(/3,3/))
double precision, dimension(3,3):: stencil5=reshape((/0.,1.,0.,1.,-4.,1.,0.,1.,0./),(/3,3/))
integer ip,ii,jj
integer :: N_phases=3
 
 LapPhi=0.
 do ip=1,N_phases
    do ii=1,3
       do jj=1,3
          LapPhi(ip)=LapPhi(ip)+stencil5(ii,jj)*vbles(ip,ii,jj)/(dx*dx)
       end do
    end do
 end do
 return
 
end subroutine Calc_GradPhi



double precision function FreeEnergy(c,T,phi,c_dot,phi_dot,lp)
use solution_parameters
implicit none

double precision ta,tb,L_a,L_b,R,mola,molb,rhoa(3),rhob(3)
!double precision,dimension (3,8)::hsa,hsb,hsrk
double precision, intent (in):: T,c,phi(3),c_dot,phi_dot(3)
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha,hb,hr,h(3),f_phi,EntropyMixing,t_a,E_const
double precision :: Tk,T_K,T_eut=456.14,Cp,RE
double precision :: MolPerVol=54185.!=(rhoa(1)+rhob(1))/(mola+molb)
double precision :: nu_Pb=54730.,nu_Sn=61505.

! double precision :: MolPerVol=1.347708895e5 !debug
integer :: i,j,TT
logical :: c_dot_term=.false. ! for extra term for PRL paper
!non dimension T to dimensioned Tk for gibbs expansion etc.

MolPerVol= nu_Pb*(1.-c)+nu_Sn*c
Tk=T_K(T)
R = 8.31 !J/(K mol)
! NU = 7.42e-6 !m^3 !taken from phase_8.08/init_eut_Pb_Sn_8.h
! Latent heat
L_a = 4773.925095 !J/mol 
L_b = 7029.016991
! Setting Lead diffusivity =1 then implies characteristic length
!delta_a= 1.489764340*10^(-7)*m
! and characteristic time
!tau_a = 9.809908287 10^(-10) s,

! not used: reference only
!!!!!!!!!!!!!!!!!!!!!!!
! rhoa = 10660. !kg/m^3 liquid Lead at melting point
! rhob = 6999.  !              Tin
! ta = 600.61 !Lead
! tb = 505.07 !Tin
!!!!!!!!!!!!!!!!!!!!!!!!
!  mola=.2072     !kg/mol
!  molb= 0.1187  
! rhoa(1)=10660. !liquid Lead kg/m^3
! rhob(1)=6999.  !liquid Tin
! rhoa(2)=11340. !Solid Lead
! rhoa(3)=11340. !The same (for want of a better option)
! rhob(2)=7365.  !Solid Tin (White)
! rhob(3)=7365.  !The same
! t_a=600.61 
! E_const=(R*rhoa(1)*t_a)/mola !~ 2.5 E8 J/m^3
! write(*,*)E_const
! stop
!note: but perhaps Delta T, undercooling is better than t_a
! E_const= 2.56779906e8 !J/m^3 for non-dimensionlising D F/D c

!convert Latent heat yo J/m^3 used in temperature equation.
!This approximation will be superceeded
!Not needed as can work with J/mol in temperature equation
! do i=1,3
! L_a(i)=L_a(i)*rhoa(i)/mola !J/m^3
! L_b(i)=L_b(i)*rhob(i)/molb
! enddo




if(lp.eq.0)then !free energy proper
 call GibbsVector(V,Tk,0)
 call MixingVector(RK,c,0)
 FreeEnergy=0.
 do i=1,3
 ha=0.
 hb=0.
 hr=0.
 do j = 1,8
  ha = ha + hsa(i,j)*V(j)
  hb = hb + hsb(i,j)*V(j)
 enddo
 do j = 1,4
   hr = hr + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
 enddo
!Convert to J/m^3 
!  ha=ha*MolPerVol
!  hb=hb*MolPerVol 
!  hr=hr*MolPerVol
 
 FreeEnergy = FreeEnergy + (c*hb+(1.-c)*ha+hr)*f_phi(phi,i,0,0)
 enddo
  FreeEnergy = FreeEnergy + R*Tk*EntropyMixing(c,0)

elseif(lp.le.N_phases)then !phase = current_var
 FreeEnergy = 0.
 do i=1,N_phases
 call GibbsVector(V,Tk,0)
 call MixingVector(RK,c,0)

 ha=0.
 hb=0.
 hr=0.
 do j = 1,8
  ha = ha + hsa(i,j)*V(j)
  hb = hb + hsb(i,j)*V(j)
 enddo
 do j = 1,4
   hr = hr + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
 enddo
!Convert to J/m^3 
 ha=ha*MolPerVol
 hb=hb*MolPerVol 
 hr=hr*MolPerVol
 FreeEnergy = FreeEnergy +(c*hb+(1.-c)*ha+hr)*f_phi(phi,i,lp,0)
 enddo
 return
elseif(lp.eq.N_phases+2)then !Solute
! the below term does not discretise well. So this term is combined
!with the diffusivity (dependent on c and phi_i and T)
!  FreeEnergy = R*Tk*EntropyMixing(c,1) !dimension (J/mol)
 FreeEnergy = 0.!debug
 call GibbsVector(V,Tk,0)
 call MixingVector(RK,c,1)

 do i= 1,N_phases
  ha=0.
  hb=0.
  hr=0.
  do j = 1,8
   ha = ha + hsa(i,j)*V(j)
   hb = hb + hsb(i,j)*V(j)
  enddo
  do j = 1,4
   hr = hr + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
  enddo
! add in other factor and non dimensionalise
! Folch Plapp
  FreeEnergy = FreeEnergy+(hb-ha+hr)*f_phi(phi,i,0,0)
 enddo
! add in other factor and non dimensionalise

 
 
 ! non-dimensionalise Free energy
  FreeEnergy = FreeEnergy/(R*T_eut)
elseif(lp.eq.N_phases+1)then !Temperature Keep in J/mol
!  T (d/dT)dF/dphi
  call GibbsVector(V,Tk,1) !first derivative of free energy
  call MixingVector(RK,c,0)
  FreeEnergy = 0.
  do i= 1,N_phases
   ha=0.
   hb=0.
   hr=0.
   do j = 1,8
    ha = ha + hsa(i,j)*V(j)
    hb = hb + hsb(i,j)*V(j)
   enddo
   do j = 1,4
    hr = hr + hsrk(i, 2*j)*RK(j)
   enddo
   do j=1,N_phases
    FreeEnergy = FreeEnergy + Tk*(c*hb+(1.-c)*ha+hr)*f_phi(phi,i,j,0)*phi_dot(j)
   enddo
  enddo
!    pure -dF/dphi term  
  call GibbsVector(V,Tk,0)
  call MixingVector(RK,c,0)
  do i= 1,N_phases
   ha=0.
   hb=0.
   hr=0.
   do j = 1,8
    ha = ha + hsa(i,j)*V(j)
    hb = hb + hsb(i,j)*V(j)
   enddo
   do j = 1,4
    hr = hr + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
   enddo
   do j=1,N_phases
    FreeEnergy = FreeEnergy - (c*hb+(1.-c)*ha+hr)*f_phi(phi,i,j,0)*phi_dot(j)
   enddo
  enddo
 if(c_dot_term)then !there is a c_dot term see PRL small for large Le
!    T(d/dT)dF/dc =-T (ds/dc) ! see AT14.pdf for details
! setting c_dot_term = .false. results in faster growth  
  call GibbsVector(V,Tk,1)
  call MixingVector(RK,c,1)
  RE = R*EntropyMixing(c,1)
  do i= 1,N_phases
   ha=0.
   hb=0.
   hr=0.
   do j = 1,8
    ha = ha + hsa(i,j)*V(j)
    hb = hb + hsb(i,j)*V(j)
   enddo
   do j = 1,4
    hr = hr +hsrk(i, 2*j)*RK(j)
   enddo
   FreeEnergy = FreeEnergy + Tk*(RE+(hb-ha+hr)*f_phi(phi,i,0,0))*c_dot
  enddo
 endif !c_dot_term
!   heat capacity = -T (d/dT) dF/dT
  call GibbsVector(V,Tk,2)
  Cp=0.
  do i= 1,N_phases
   ha=0.
   hb=0.
   do j = 1,8
    ha = ha + hsa(i,j)*V(j)
    hb = hb + hsb(i,j)*V(j)
   enddo
   Cp = Cp - (c*hb+(1.-c)*ha)*f_phi(phi,i,0,0)
  enddo
   Cp = Cp*Tk
  FreeEnergy = FreeEnergy/Cp
else
 write(*,*)"bad current_var"
 stop
endif
end function FreeEnergy

double precision function EntropyMixing(c,d)
double precision, intent (in)::c
integer, intent (in)::d

! EntropyMixing =0.
! if((1.-c).le.1e-1)then
! if(d.eq.1) EntropyMixing = 2.197224577
! return
! endif
 
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

! double precision function g_phi(phi,d)
!   double precision, intent (in):: phi
!   integer, intent (in):: d
!   if(d.eq.1)then
!    g_phi = 30.*phi**2*(1.0-phi)**2
!   elseif(d.eq.2)then
!    g_phi = 60.0*phi*(2.0*phi-1.0)*(phi-1.0)
!   else !assume 0
!    g_phi = 6.0*phi**5 - 15.0*phi**4 + 10.0*phi**3
!   endif
! end function g_phi
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
double precision function gg_phi(phi,i,j,k)
  double precision, intent (in):: phi(3)
  integer, intent (in):: i,j,k
  double precision g_phi 
   
   if(k.eq.0)then
    if(j.eq.0)then
     gg_phi =  g_phi(phi(i),0)
    elseif(j.eq.i)then
     gg_phi =  g_phi(phi(i),1)
    else
     gg_phi = 0.
    endif
   elseif(k.eq.i.and.i.eq.j)then
     gg_phi =  g_phi(phi(i),2)
   else
    gg_phi = 0.
   endif
      
end function gg_phi
double precision function f_phi(phi,i,j,k)
  double precision, intent (in):: phi(3)
  double precision g_phi,s,gg_phi
  integer, intent (in):: i,j,k
  s=g_phi(phi(1),0)+g_phi(phi(2),0)+g_phi(phi(3),0)
  if(j.eq.0)then
   f_phi = g_phi(phi(i),0)/s
  elseif(k.eq.0)then
   f_phi = gg_phi(phi,i,j,0)/s-g_phi(phi(i),0)*gg_phi(phi,j,j,0)/s**2
  else
   f_phi = gg_phi(phi,i,j,k)/s - gg_phi(phi,i,j,0)*gg_phi(phi,k,k,0)/s**2&
        - gg_phi(phi,i,k,0)*gg_phi(phi,j,j,0)/s**2 - g_phi(phi(i),0)*gg_phi(phi,j,j,k)/s**2 &
        + 2.*gg_phi(phi,j,j,0)*g_phi(phi(i),0)*gg_phi(phi,k,k,0)/s**3
  endif
end function f_phi

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



double precision function flatten(x)
double precision, intent (in)::x
double precision:: h=0.01 !decrease to flatten more
! flatten=h*tanh(x/h)
flatten=max(-h,min(h,x))
end function flatten
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


double precision function phi_(phi,ip)
double precision, intent (in)::phi(3)
integer,intent (in):: ip
! phi_=(3.*phi(ip)-phi(1)-phi(2)-phi(3)+1.)/3.
phi_=phi(ip)/(phi(1)+phi(2)+phi(3))
end function phi_
