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

  double precision :: D_therm, A_0,  delta, scale_factor,epsilon, epsilon_tilde
  double precision :: D_solute, ke, le, mcinf,T_scale,c_scale,v_mesh

  double precision :: nuc_radius
  double precision :: anti_trapping_mod=1.0,maxsol=0.0
  logical :: init_c=.true.
  logical :: solute
  integer N_phases
  double precision :: max_x,c_min,c_max
  double precision :: hsa(3,8),hsb(3,8),hsrk(3,8)
  double precision :: F_time=0d0,dF_time=0d0,inv_time=0d0
  double precision :: width = 10d0
  integer :: rrblock(12,5),rrarray(4,2),rrrblock(64,5),rrrarray(16,2)
  integer :: N_unk,M_unk,NN_unk
! parameters for Ni-Cu
  double precision :: g_M(2),g_sigma(2)=(/0.37,0.29/),g_T(2)=(/1728.,1358./),g_W(2)
  double precision :: g_L(2)=(/2.35e9,1.728e9/),g_delta(2)=(/2.2e-9 ,2.2e-9 /),g_mu(2)=(/0.0033,0.0039/),g_epsilon(2)
  double precision :: g_C(2)=(/5.42e6,3.96e6/),g_MolPerVol(2)=(/1.517d5,1.406d5/),g_D(2)=(/1e-9,1e-13/) !liquid, solid
  double precision :: g_tau,g_lengthSQ,g_Dch, g_R = 8.31,g_T0=1549.15,g_d0=5d-10,g_lambda=1.0 !J/(K mol)
  logical :: heat_c_term=.true.,heat_phi_term=.true.,grad_DW_term=.true.,cross_term=.true.,thermal=.true.,g_1D=.false.
  integer :: g_sym= -1,g_max_lvl=1,g_output_rate=100,g_nbx=1,g_nby=1,g_nbz=1
  double precision :: heat_size(4)=0d0,g_c0=0.5,g_D_tilde=1d2,Q_heat(4),sum_heat(4)=0d0,g_beta=0.0,g_Le=100.0,g_max_time=1d4
  double precision :: g_grad(3,3),g_max_dt=1.0,g_anti_trap,g_c_force, g_grid
  double precision :: d_10_10(10,10)
! end parameters for Ni-Cu

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

