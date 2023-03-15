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
  integer :: rrblock(12,5),rrarray(4,2),rrrblock(64,5),rrrarray(16,2),rrblock3d(32,6),rrarray3d(8,3)
  integer :: N_unk,M_unk,NN_unk
! parameters for Al-Ni (modified in mp_problem.f90)
  double precision :: g_M(2),g_sigma(2)=(/0.37,0.29/),g_T(2)=(/1728.,1358./),g_W(2)
  double precision :: g_L(2)=(/2.35e9,1.728e9/),g_delta(2)=(/2.2e-9 ,2.2e-9 /),g_mu(2)=(/0.0033,0.0039/),g_epsilon(2)
  double precision :: g_C,g_MolPerVol(2)=(/1.517d5,1.406d5/),g_D(2)=(/1e-9,1e-13/) !liquid, solid
  double precision :: g_tau,g_lengthSQ,g_Dch, g_R = 8.31,g_T0=1549.15,g_d0=5d-10,g_lambda=1.0 !J/(K mol)
  logical :: heat_c_term=.true.,heat_phi_term=.true.,grad_DW_term=.true.,cross_term=.true.,thermal=.true.,g_1D=.false.
  integer :: g_sym= -1,g_max_lvl=1,g_output_rate=100,g_nbx=1,g_nby=1,g_nbz=1,g_grid,AlNi=2
  double precision :: heat_size(4)=0d0,g_c0=0.5,g_D_tilde=1d2,Q_heat(4),sum_heat(4)=0d0,g_beta=0.0,g_Le=100.0,g_max_time=1d4
  double precision :: g_grad(3,3),g_max_dt=1.0,g_anti_trap,g_c_force
  double precision ::g_s1,g_s2,g_s3

  double precision :: hlA(8)=(/-795.996, 177.430178, -31.748192, 0, 0, 0, 0, 0/)
  double precision :: hlB(8)=(/11235.527, 108.457, -22.096, -0.0048047, 0, 0, -3.82318e-21, 0/)
  double precision :: hslb1(8)=(/-354240.204,  848.531459, -139.436576, -0.0096814, 0, 0, 0, -3.691572e+28/)
  double precision :: hslb2(8)=(/24023.11,916.355765,-158.74096,0,0,0,0,-6.15262e+28/)
  double precision :: hslb3(8)=(/-420170.259,973.724709,-161.532576,-0.0145221,0,0,0,-3.691572e+28/)
  double precision :: hslb4(8)=(/-41906.943,1041.549015,-180.83696,-0.0048407,0,0,0,-6.15262e+28/)
  double precision :: hrkl0(2)=(/-207109.28,41.31501/),hrkl1(2)=(/-10185.79,5.8714/),hrkl2(2)=(/81204.81,-31.95713/)
  double precision :: hrkl3(2)=(/4365.35,-2.51632/),hrkl4(2)=(/-22101.64,13.16341/)



  double precision ::hrka0(2)=(/-162407.75, 16.212965/)
  double precision ::hrka1(2)=(/73417.798, -34.914168/)
  double precision ::hrka2(2)=(/33471.014, -9.8373558/)
  double precision ::hrka3(2)=(/-30758.01, 10.25267/)

  double precision :: hslb5(2)=(/-193484.18,131.79/),hslb6(2)=(/-22001.7,7.0332/) 
  double precision ::hlA1(8)=(/3028.879,  125.251171, -24.3671976,-0.001884662,-8.77664e-7,74092,7.9337e-20,0/) 
  double precision ::hlA2(8)=(/-271.21, 211.206579,-38.5844296,0.018531982,-5.764227e-06,74092,7.9337e-20,0/) 
  double precision ::hlA3(8)=(/-795.996,177.430178,-31.748192,0,0,0,0,0/) 
  double precision ::hlB1(8)=(/11235.527,108.457,-22.096,-0.0048047,0,0,-3.82318e-21,0/) 
  double precision ::hlB2(8)=(/-9549.817,268.597977,-43.1,0,0,0,0,0/) 
  double precision ::haA1(8)=(/-7976.15,137.093038,-24.3671976,-0.001884662,-8.77664e-7,74092.0,0,0/) 
  double precision ::haA2(8)=(/-11276.24,223.048446,-38.5844296,0.018531982,-5.764227e-6,74092.0,0,0/)
  double precision ::haA3(8)=(/-11278.378,188.684153,-31.748192,0,0,0,0,-1.230524e+28/) 
  double precision ::haB1(8)=(/-5179.159,117.854,-22.096,-0.0048407,0,0,0,0/) 
  double precision ::haB2(8)=(/-27840.62,279.134977,-43.1,0,0,0,0,1.12754e+31/) 
  double precision ::hc1(8)=(/-55760.63225, 144.5824085,-23.7993982,-0.0026236715,-6.58248e-7,55569,0,0/) 
  double precision ::hc2(8)=(/-58235.69975, 209.0489645,-34.4623222, 0.0126888115,-4.32317025e-6, 55569,0,0/) 
  double precision ::hc3(8)=(/-58237.30325, 183.27574475,-29.335144,-0.001210175,0,0,0,-9.22893e+27/) 
  double precision ::hc4(8)=(/-63902.6685, 223.595989,-34.586144,0,0,0,0,2.80962107e+30/) 
  
  double precision ::hslb1_1(8)=(/-344333.52, 693.758114,-117.2935928,-0.015335386,-2.632992e-6, 222276,0, 0/) 
  double precision ::hslb1_2(8)=(/-354233.79,951.624338,-159.9452888,0.045914546,-1.7292681e-5,222276,0,0/) 
  double precision ::hslb1_3(8)=(/-354240.204,848.531459,-139.436576,-0.0096814,0,0,0,-3.691572e+28/) 
  double precision ::hslb1_4(8)=(/-399563.126,1171.093413,-181.444576,0,0,0,0,2.251388428e+31/) 
  double precision ::hslb2_1(8)=(/40534.25,658.40019,-121.835988,-0.00942331,-4.38832e-6,370460,0,0/) 
  double precision ::hslb2_2(8)=(/24033.8,1088.17723,-192.922148,0.09265991,-2.8821135e-5,370460,0,0/) 
  double precision ::hslb2_3(8)=(/24023.11,916.355765,-158.74096,0,0,0,0,-6.15262e+28/) 
  double precision ::hslb3_1(8)=(/-410263.575,818.951364,-139.3895928,-0.020176086,-2.632992e-6,222276,0,0/) 
  double precision ::hslb3_2(8)=(/-420163.845,1076.817588,-182.0412888,0.041073846,-1.7292681e-05,222276,0,0/) 
  double precision ::hslb3_3(8)=(/-420170.259,973.724709,-161.532576,-0.0145221,0,0,0,-3.691572e+28/) 
  double precision ::hslb3_4(8)=(/-488154.642,1457.56764,-224.544576,0,0,0,0,3.378928428e+31/) 
  double precision ::hslb4_1(8)=(/-25395.803,783.59344,-143.931988,-0.01426401,-4.38832e-6,370460,0,0/) 
  double precision ::hslb4_2(8)=(/-41896.253,1213.37048,-215.018148,0.08781921,-2.8821135e-05,370460,0,0/) 
  double precision ::hslb4_3(8)=(/-41906.943,1041.549015,-180.83696,-0.0048407,0,0,0,-6.15262e+28/) 
  double precision ::hslb4_4(8)=(/-64568.404,1202.829992,-201.84096,0,0,0,0,1.121387380e+31/) 
  
  
  
  double precision ::hsle1_1(8)=(/-146754.525,272.983788,-46.4631976,-0.006725362,-8.77664e-7,74092,0,0/) 
  double precision ::hsle1_2(8)=(/-150054.615,358.939196,-60.6804296,0.013691282,-5.764227e-6,74092,0,0/) 
  double precision ::hsle1_3(8)=(/-150056.753,324.574903,-53.8441920,-0.004840700,0,0,0,-1.230524e+28/) 
  double precision ::hsle1_4(8)=(/-172718.214,485.855880,-74.8481920,0,0,0,0,1.126309476e+31/) 
  double precision ::hsle2_1(8)=(/12106.85,131.280038,-24.3671976,-0.001884662,-8.776640e-7,74092,0,0/) 
  double precision ::hsle2_2(8)=(/8806.760,217.235446,-38.5844296,0.018531982,-5.764227e-6,74092,0,0/) 
  double precision ::hsle2_3(8)=(/8804.622,182.871153,-31.7481920,0,0,0,0,-1.230524e+28/) 
  double precision ::hsle3_1(8)=(/7071.850,228.596000,-44.1920000,-0.0096814,0,0,0,0/) 
  double precision ::hsle3_2(8)=(/-38251.072,551.157954,-86.2,0,0,0,0,2.25508e+31/) 
  double precision ::hsle4_1(8)=(/165933.225,86.89225,-22.096,-0.0048407,0,0,0,0/) 
  double precision ::hsle4_2(8)=(/143271.764,248.173227,-43.10,0,0,0,0,1.12754e+31/) 
  
  double precision ::hsle5(2)=(/-64024.38,26.49419/) 
  double precision ::hsle6(2)=(/-52440.88,11.30117/) 
  



! end parameters for Al-Ni

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

