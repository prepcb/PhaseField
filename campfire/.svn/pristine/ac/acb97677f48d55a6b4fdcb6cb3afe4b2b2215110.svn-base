!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  solution_parameters
!!!!  * Module containing application specific parameters
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Latest version: Chris Goodyer, May 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module solution_parameters
  ! Module for storing all the parameters for the problem being solved.
  ! End user should keep track of which are and aren't used
  ! The only critical ones are 
  ! t0: start time. For most cases this is 0 but if you need to use a non-zero 
  ! start time then look no further. This is here in case start time depends on
  ! other solution parameters.
  ! nucleate_x/y/z are the co-ordinates of the nucleation point
  ! nuc_radius is the nucleation radius
  ! anti_trapping_mod is a multiplier for the antitrapping terms, normally 1.0
  ! solute sets whether code is phase/thermal (.false.) or phase/solute (.true.) IFF only 2 real vars
  ! For thermal/phase le still NOT used to set D_therm
  double precision, parameter :: a_2=0.626663849
  double precision, parameter :: a_1=5.0*sqrt(2.0)/8.0

  integer :: eighthdmn

  double precision :: D_therm, A_0, lambda, delta, scale_factor,epsilon, epsilon_tilde
  double precision :: D_solute, ke, le, mcinf,T_scale,c_scale,v_mesh

  double precision :: nucleate_x, nucleate_y, nucleate_z, nuc_radius
  double precision :: anti_trapping_mod=1.0,maxsol=0.0
  double precision :: eps2(3),eps22(3,3),mobility(3,3)=1d0
  logical :: init_c=.true.
  logical :: solute
  integer N_phases
  double precision :: max_x,c_min,c_max
  double precision :: hsa(3,8),hsb(3,8),hsrk(3,8)
  double precision :: F_time=0d0,dF_time=0d0,inv_time=0d0
  double precision :: width = 10d0
  integer :: rrblock(12,5),rrarray(4,2),rrrblock(64,5),rrrarray(16,2)
  double precision :: lamda(3,3)=0d0,c12(3),c23(3),c31(3),cc12(3),cc23(3),cc31(3)
  
contains
!  subroutine get_h_psi(psi,h)
!    implicit none
!    double precision, intent(in) :: psi
!    double precision, intent(out) :: h
!    h=psi
!  end subroutine get_h_psi
  subroutine get_dh_dpsi(psi,dh)
    implicit none
    double precision, intent(in) :: psi
    double precision, intent(out) :: dh
    dh=1.0
  end subroutine get_dh_dpsi
!  subroutine get_d2h_dpsi2(psi,d2h)
!    implicit none
!    double precision, intent(in) :: psi
!    double precision, intent(out) :: d2h
!    d2h=0.0
!  end subroutine get_d2h_dpsi2
end module solution_parameters

