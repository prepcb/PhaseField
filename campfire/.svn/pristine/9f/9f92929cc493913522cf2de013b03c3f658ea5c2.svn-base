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

  integer :: eighthdmn, g_block
  double precision :: g_phase,g_gamma=0.8,g_xy_position(2), g_slider=0d0
  double precision :: g_Diff0(10,10),g_Diff1(10,10),g_Diff2(10,10)
  double precision :: D_therm, A_0,  delta, scale_factor,epsilon, epsilon_tilde,g_x0=0,g_y0=0,g_obstacle=0d0
  double precision :: D_solute, ke, le, mcinf,T_scale,c_scale,v_mesh,g_phex(2,6),g_phi_dot

  double precision :: nuc_radius,g_drive
  double precision :: anti_trapping_mod=1.0,maxsol=0.0
  logical :: init_c=.true.,g_linear=.true.,g_Isotropic=.false.,g_axis_symmetric=.false.,g_Halley=.false.
  logical :: solute,g_diagonal=.false.,g_triangular=.false.
  logical :: g_curve = .false.
  integer :: N_phases,g_AMM=2,g_order=10
!  double precision :: g_left_limit=0.235d0,g_right_limit=0.272463768d0
  double precision :: g_left_limit=0.2350d0,g_right_limit=0.261463,g_left_tan=0d0,g_left_acc=0d0,g_right_tan=0d0,g_right_acc=0d0,g_right_3=0d0
  double precision :: hsa(3,8),hsb(3,8),hsrk(3,8),max_x
  double precision :: F_time=0d0,dF_time=0d0,inv_time=0d0,g_time=0d0
  double precision :: width = 10d0, g_cinfty=0.5
  integer :: N_unk,M_unk,NN_unk
! parameters for Al-Ni (modified in mp_problem.f90)
  double precision :: g_M(2),g_sigma(2)=(/0.37,0.29/),g_T(2)=(/1728.,1358./),g_W(2)
  double precision :: g_L(2)=(/2.35e9,1.728e9/),g_delta(2)=(/2.2e-9 ,2.2e-9 /),g_mu(2)=(/0.0033,0.0039/),g_epsilon(2),g_weight,g_mu_kinetic=1d0
  double precision :: g_C,g_MolPerVol(2)=(/1.517d5,1.406d5/),g_D(2)=(/1e-9,1e-13/) !liquid, solid
  double precision :: g_tau,g_lengthSQ,g_Dch, g_R = 8.31,g_T0=0d0,g_T_i=0d0,g_d0=5d-10,g_lambda=1.0 !J/(K mol)
  logical :: heat_c_term=.true.,heat_phi_term=.true.,grad_DW_term=.true.,cross_term=.true.,thermal=.true.,g_1D=.false.,hex_facet=.false.
  integer :: g_sym= -1,g_max_lvl=1,g_output_rate=100,g_nbx=1,g_nby=1,g_nbz=1,g_grid,AlNi=2,g_npoints=0
  double precision :: heat_size(4)=0d0,g_c0=0.5,g_D_tilde=1d2,Q_heat(4),sum_heat(4)=0d0,g_beta=0.0,g_Le=100.0,g_max_time=1d6
  double precision :: g_grad(3,3),g_max_dt,g_anti_trap,g_c_force,g_alpha,g_at_factor=1d0,g_vel=0.02,g_vel_max,g_ratio=1d0,g_ratio0=1d0,g_alpha0,g_kurv=0d0
  double precision ::g_s1,g_s2,g_s3, g_latent_heat, g_interface
  double precision :: g_c_min,g_c_max,g_f_c_min,g_f_c_max,g_f1_c_min,g_f2_c_min,g_f1_c_max,g_f2_c_max,g_c_mid,g_FE(2,2),g_f_c,g_f1,g_f1c,g_f1cc,g_f2,g_f2c,g_f2cc,g_c_par(4)
  double precision :: g_c_min_V(0:80),g_c_max_V(0:80),g_c_mid_V(0:80),g_FE_V(2,2,0:80),g_f_c_V(0:80),g_dfdcT0(2,2),g_dfdT0(2,2),g_d2fdT0(2,2),g_q0=0.2,g_power=2d0,g_s=0d0,g_r_factor=1d0
  logical :: g_approx = .false.,g_kinetic=.false.,g_k_mult=.true.,g_skew_potential=.false., g_power_skew=.false.,g_s_skew=.false.,g_dw_skew,g_plapp=.false.
  logical :: g_pure = .false.,g_quad_extension=.false.,g_ellipse=.false.,g_equal_delta=.false.,g_needle=.false.,g_AlFe_flag=.true.,new_potential=.false.,g_LUT=.false.,g_Quad=.true.
  doubleprecision :: g_AlFe_c(2,6)=0d0,g_Cphi(-2:100002,6),g_df,g_AlFe_LUT(-2:100002,6)
  integer :: g_Ncphi = 100000
!  double precision :: g_y10(10)=(/-1., -0.7777777778, -0.5555555556, -0.3333333333, -0.1111111111, 0.1111111111, 0.3333333333, 0.5555555556, 0.7777777778, 1./)
!  double precision :: g_y10(10)=(/-1., -0.875, -0.625, -0.375, -0.125, 0.125, 0.375, 0.625, 0.875, 1./)

  doubleprecision :: g_y10(10)=(/-1.,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1./)
!
!AlFe

double precision ::AlFeR=8.31 

double precision ::AlFehlA1(8)=(/3028.879, 125.251171, -24.3671976, -0.001884662, -8.77664e-7, 74092d0, 7.9337e-20, 0d0/) 
double precision ::AlFehlA2(8)=(/-271.21, 211.206579, -38.5844296, 0.018531982, -5.764227e-6, 74092d0, 7.9337e-20, 0d0/) 
double precision ::AlFehlA3(8)=(/-795.996, 177.430178, -31.748192, 0d0, 0d0, 0d0, 0d0, 0d0/) 
double precision ::AlFehlB1(8)=(/13265.87, 117.57557, -23.5143, -0.00439752, -5.8927e-8, 77359d0, -3.67516e-21, 0d0/) 
double precision ::AlFehlB2(8)=(/-10838.83, 291.302, -46.0, 0d0, 0d0, 0d0, 0d0, 0d0/) 
double precision ::AlFehaA1(8)=(/-7976.15, 137.093038, -24.3671976, -0.001884662, -8.77664e-7, 74092.0, 0d0, 0d0/) 
double precision ::AlFehaA2(8)=(/-11276.24, 223.048446, -38.5844296, 0.018531982, -5.764227e-6, 74092.0, 0d0, 0d0/)
double precision ::AlFehaA3(8)=(/-11278.378, 188.684153, -31.748192, 0d0, 0d0, 0d0, 0d0, -1.230524e+28/) 
double precision ::AlFehaB1(8)=(/-236.7, 132.416, -24.6643, -0.00375752, -5.8927e-8, 77359d0, 0d0, 0d0/) 
double precision ::AlFehaB2(8)=(/-27097.3963, 300.252559, -46.0, 0d0, 0d0, 0d0, 0d0, 2.78854e+31/) 

double precision ::AlFehslb1_1(8)=(/-36528.11525, 141.48766407, -24.166766664, -0.00247518363, -6.85260805e-7, 74859.745, 0d0, 0d0/) 
double precision ::AlFehslb1_2(8)=(/-39052.6841, 207.24355119, -35.042949144, 0.01314354903, -4.4234815e-6, 74859.745, 0d0, 0d0/) 
double precision ::AlFehslb1_3(8)=(/-39054.31967, 180.954867045, -29.81322738, -0.0010334172, -1.3847845e-8, 18179.365, 0d0, -9.4135086e+27/) 
double precision ::AlFehslb1_4(8)=(/-45307.500705, 222.121826295, -35.09736688, 0d0, 0d0, 0d0, 0d0, 5.386256991e+30/) 
double precision ::AlFehslb2_1(8)=(/-32498.294625, 122.453971345, -20.816276994, -0.002216042605, -5.64582005e-07, 64672.095, 0d0, 0d0/) 
double precision ::AlFehslb2_2(8)=(/-34569.1011, 176.390989865, -29.737590074, 0.010595401505, -3.6309002875e-6, 64672.095, 0d0, 0d0/) 
double precision ::AlFehslb2_3(8)=(/-34570.442695, 154.827396007, -25.44785098, -0.0010334172, -1.3847845e-8, 18179.365, 0d0, -7.7215381e+27/) 
double precision ::AlFehslb2_4(8)=(/-40823.62373, 195.994355257, -30.73199048, 0d0, 0d0, 0d0, 0d0, 5.387948962e+30/)


double precision ::AlFehrkl0(2)=(/-91976.5, 22.1314/)
double precision ::AlFehrkl1(2)=(/-5672.58, 4.8728/)
double precision ::AlFehrkl2(2)=(/121.9, 0d0/)
double precision ::AlFehrka0(2)=(/-76066.1, 18.6758/)
double precision ::AlFehrka1(2)=(/21167.4, 1.3398/)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision :: hlA(8)=(/-795.996, 177.430178, -31.748192, 0d0, 0d0, 0d0, 0d0, 0d0/)
  double precision :: hlB(8)=(/11235.527, 108.457, -22.096, -0.0048047, 0d0, 0d0, -3.82318e-21, 0d0/)
  double precision :: hslb1(8)=(/-354240.204,  848.531459, -139.436576, -0.0096814, 0d0, 0d0, 0d0, -3.691572e+28/)
  double precision :: hslb2(8)=(/24023.11,916.355765,-158.74096,0d0,0d0,0d0,0d0,-6.15262e+28/)
  double precision :: hslb3(8)=(/-420170.259,973.724709,-161.532576,-0.0145221,0d0,0d0,0d0,-3.691572e+28/)
  double precision :: hslb4(8)=(/-41906.943,1041.549015,-180.83696,-0.0048407,0d0,0d0,0d0,-6.15262e+28/)
  double precision :: hrkl0(2)=(/-207109.28,41.31501/),hrkl1(2)=(/-10185.79,5.8714/),hrkl2(2)=(/81204.81,-31.95713/)
  double precision :: hrkl3(2)=(/4365.35,-2.51632/),hrkl4(2)=(/-22101.64,13.16341/)



  double precision ::hrka0(2)=(/-162407.75, 16.212965/)
  double precision ::hrka1(2)=(/73417.798, -34.914168/)
  double precision ::hrka2(2)=(/33471.014, -9.8373558/)
  double precision ::hrka3(2)=(/-30758.01, 10.25267/)

  double precision :: hslb5(2)=(/-193484.18,131.79/),hslb6(2)=(/-22001.7,7.0332/) 
  double precision ::hlA1(8)=(/3028.879,  125.251171, -24.3671976,-0.001884662,-8.77664e-7,74092.,7.9337d-20,0d0/) 
  double precision ::hlA2(8)=(/-271.21, 211.206579,-38.5844296,0.018531982,-5.764227e-06,74092.,7.9337e-20,0d0/) 
  double precision ::hlA3(8)=(/-795.996,177.430178,-31.748192,0d0,0d0,0d0,0d0,0d0/) 
  double precision ::hlB1(8)=(/11235.527,108.457,-22.096,-0.0048047,0d0,0d0,-3.82318e-21,0d0/) 
  double precision ::hlB2(8)=(/-9549.817,268.597977,-43.1,0d0,0d0,0d0,0d0,0d0/) 
  double precision ::haA1(8)=(/-7976.15,137.093038,-24.3671976,-0.001884662,-8.77664e-7,74092.0d0,0d0,0d0/) 
  double precision ::haA2(8)=(/-11276.24,223.048446,-38.5844296,0.018531982,-5.764227e-6,74092.0d0,0d0,0d0/)
  double precision ::haA3(8)=(/-11278.378,188.684153,-31.748192,0d0,0d0,0d0,0d0,-1.230524e+28/) 
  double precision ::haB1(8)=(/-5179.159,117.854,-22.096,-0.0048407,0d0,0d0,0d0,0d0/) 
  double precision ::haB2(8)=(/-27840.62,279.134977,-43.1,0d0,0d0,0d0,0d0,1.12754e+31/) 
  double precision ::hc1(8)=(/-55760.63225, 144.5824085,-23.7993982,-0.0026236715,-6.58248e-7,55569.,0d0,0d0/) 
  double precision ::hc2(8)=(/-58235.69975, 209.0489645,-34.4623222, 0.0126888115,-4.32317025e-6, 55569.,0d0,0d0/) 
  double precision ::hc3(8)=(/-58237.30325, 183.27574475,-29.335144,-0.001210175,0d0,0d0,0d0,-9.22893e+27/) 
  double precision ::hc4(8)=(/-63902.6685, 223.595989,-34.586144,0d0,0d0,0d0,0d0,2.80962107d+30/) 
  
  double precision ::hslb1_1(8)=(/-344333.52, 693.758114,-117.2935928,-0.015335386,-2.632992e-6, 222276.,0d0, 0d0/) 
  double precision ::hslb1_2(8)=(/-354233.79,951.624338,-159.9452888,0.045914546,-1.7292681e-5,222276.,0d0,0d0/) 
  double precision ::hslb1_3(8)=(/-354240.204,848.531459,-139.436576,-0.0096814,0d0,0d0,0d0,-3.691572e+28/) 
  double precision ::hslb1_4(8)=(/-399563.126,1171.093413,-181.444576,0d0,0d0,0d0,0d0,2.251388428e+31/) 
  double precision ::hslb2_1(8)=(/40534.25,658.40019,-121.835988,-0.00942331,-4.38832e-6,370460d0,0d0,0d0/) 
  double precision ::hslb2_2(8)=(/24033.8,1088.17723,-192.922148,0.09265991,-2.8821135e-5,370460d0,0d0,0d0/) 
  double precision ::hslb2_3(8)=(/24023.11,916.355765,-158.74096,0d0,0d0,0d0,0d0,-6.15262e+28/) 
  double precision ::hslb3_1(8)=(/-410263.575,818.951364,-139.3895928,-0.020176086,-2.632992e-6,222276.,0d0,0d0/) 
  double precision ::hslb3_2(8)=(/-420163.845,1076.817588,-182.0412888,0.041073846,-1.7292681e-05,222276.,0d0,0d0/) 
  double precision ::hslb3_3(8)=(/-420170.259,973.724709,-161.532576,-0.0145221,0d0,0d0,0d0,-3.691572e+28/) 
  double precision ::hslb3_4(8)=(/-488154.642,1457.56764,-224.544576,0d0,0d0,0d0,0d0,3.378928428e+31/) 
  double precision ::hslb4_1(8)=(/-25395.803,783.59344,-143.931988,-0.01426401,-4.38832e-6,370460d0,0d0,0d0/) 
  double precision ::hslb4_2(8)=(/-41896.253,1213.37048,-215.018148,0.08781921,-2.8821135e-05,370460d0,0d0,0d0/) 
  double precision ::hslb4_3(8)=(/-41906.943,1041.549015,-180.83696,-0.0048407,0d0,0d0,0d0,-6.15262e+28/) 
  double precision ::hslb4_4(8)=(/-64568.404,1202.829992,-201.84096,0d0,0d0,0d0,0d0,1.121387380e+31/) 
  
  
  
  double precision ::hsle1_1(8)=(/-146754.525,272.983788,-46.4631976,-0.006725362,-8.77664e-7,74092.,0d0,0d0/) 
  double precision ::hsle1_2(8)=(/-150054.615,358.939196,-60.6804296,0.013691282,-5.764227e-6,74092.,0d0,0d0/) 
  double precision ::hsle1_3(8)=(/-150056.753,324.574903,-53.8441920d0,-0.004840700d0,0d0,0d0,0d0,-1.230524e+28/) 
  double precision ::hsle1_4(8)=(/-172718.214,485.855880d0,-74.8481920d0,0d0,0d0,0d0,0d0,1.126309476e+31/) 
  double precision ::hsle2_1(8)=(/12106.85,131.280038,-24.3671976,-0.001884662,-8.776640e-7,74092.,0d0,0d0/) 
  double precision ::hsle2_2(8)=(/8806.760d0,217.235446,-38.5844296,0.018531982,-5.764227e-6,74092.,0d0,0d0/) 
  double precision ::hsle2_3(8)=(/8804.622,182.871153,-31.7481920d0,0d0,0d0,0d0,0d0,-1.230524e+28/) 
  double precision ::hsle3_1(8)=(/7071.850d0,228.596000d0,-44.1920000d0,-0.0096814,0d0,0d0,0d0,0d0/) 
  double precision ::hsle3_2(8)=(/-38251.072,551.157954,-86.2,0d0,0d0,0d0,0d0,2.25508e+31/) 
  double precision ::hsle4_1(8)=(/165933.225,86.89225,-22.096,-0.0048407,0d0,0d0,0d0,0d0/) 
  double precision ::hsle4_2(8)=(/143271.764,248.173227,-43.10d0,0d0,0d0,0d0,0d0,1.12754e+31/) 
  
  double precision ::hsle5(2)=(/-64024.38,26.49419/) 
  double precision ::hsle6(2)=(/-52440.88,11.30117/) 
  
! end parameters for Al-Ni

!begin parameters for Si-Ge ;


double precision ::hrklSiGe0(2)=(/6610.0d0,-0.354/)
double precision ::hrklSiGe1(2)=(/0d0,0d0/)
double precision ::hrklSiGe2(2)=(/0d0,0d0/)

double precision ::hrksSiGe0(2)=(/3500.0d0,0d0/)
double precision ::hrksSiGe1(2)=(/0d0,0d0/)
double precision ::hrksSiGe2(2)=(/0d0,0d0/)


double precision ::hlSi1(8)=(/42533.751,107.13742,-22.8317533,-0.001912904,-3.552e-9,176667.,2.09307e-21,0d0/)
double precision ::hlSi2(8)=(/40370.523,137.722298,-27.196,0d0,0d0,0d0,0d0,0d0/)

double precision ::hlGe1(8)=(/27655.337,134.94853,-29.5337682,0.005568297,-1.51369400000e-6,163298.,8.56632e-21,0d0/)
double precision ::hlGe2(8)=(/31452.25,72.173826,-19.8536239,-0.003672527,0d0,0d0,8.56632e-21,0d0/)
double precision ::hlGe3(8)=(/27243.473,126.324186,-27.6144,0d0,0d0,0d0,0d0,0d0/)

double precision ::hsSi1(8)=(/-8162.609,137.2368590d0,-22.8317533,-0.001912904,-3.552e-9,176667.,0d0,0d0/)
double precision ::hsSi2(8)=(/-9457.642,167.281367,-27.196,0d0,0d0,0d0,0d0,-4.20369d30/)

double precision ::hsGe1(8)=(/-9486.153,165.635573,-29.5337682,0.005568297,-1.513694e-6,163298.,0d0,0d0/)
double precision ::hsGe2(8)=(/-5689.239,102.86087,-19.8536239,-0.003672527,0d0,0d0,0d0,0d0/)
double precision ::hsGe3(8)=(/-9548.204,156.708024,-27.6144,0d0,0d0,0d0,0d0,-8.59809e+28/)



!end parameters for Si-Gi

!begin parameters for AlSi
!AlSi

double precision ::AlSihlA1(8)=(/3028.879,125.251171,-24.3671976,-0.001884662,-8.77664e-7,74092.,7.9337e-20,0d0/)
double precision ::AlSihlA2(8)=(/-271.21,211.206579,-38.5844296,0.018531982,-5.764227e-6,74092.,7.9337e-20,0d0/)
double precision ::AlSihlA3(8)=(/-795.996,177.430178,-31.748192,0d0,0d0,0d0,0d0,0d0/)

double precision ::AlSihlB1(8)=(/42533.751,107.13742,-22.8317533,-0.001912904,-3.552e-9,176667.,2.09307e-21,0d0/)
double precision ::AlSihlB2(8)=(/40370.523,137.722298,-27.196,0d0,0d0,0d0,0d0,0d0/)

double precision ::AlSihs1A1(8)=(/-7976.15,137.093038,-24.3671976,-0.001884662,-8.77664e-7,74092.,0d0,0d0/)
double precision ::AlSihs1A2(8)=(/-11276.24,223.048446,-38.5844296,0.018531982,-5.764227e-6,74092.,0d0,0d0/)
double precision ::AlSihs1A3(8)=(/-11278.378,188.684153,-31.748192,0d0,0d0,0d0,0d0,-1.230524e28/)

double precision ::AlSihs1B1(8)=(/42837.391,115.436859,-22.8317533,-0.001912904,-3.552e-9,176667.,0d0,0d0/)
double precision ::AlSihs1B2(8)=(/41542.358,145.481367,-27.196,0d0,0d0,0d0,0d0,-4.20369e30/)

double precision ::AlSihs2A1(8)=(/-7976.15,167.093038,-24.3671976,-0.001884662,-8.77664e-7,74092.,0d0,0d0/)
double precision ::AlSihs2A2(8)=(/-11276.24,253.048446,-38.58442960d0,0.018531982,-5.764227e-6,74092.,0d0,0d0/)
double precision ::AlSihs2A3(8)=(/-11278.378,218.684153,-31.748192,0d0,0d0,0d0,0d0,-1.230524e28/)

double precision ::AlSihs2B1(8)=(/-8162.609,137.2368590d0,-22.8317533,-0.001912904,-3.552e-9,176667.,0d0,0d0/)
double precision ::AlSihs2B2(8)=(/-9457.642,167.281367,-27.196,0d0,0d0,0d0,0d0,-4.20369e30/)

double precision ::AlSihrkl0(2)=(/-11340.1,-1.23394/)
double precision ::AlSihrkl1(2)=(/-3530.93,1.35993/)
double precision ::AlSihrkl2(2)=(/2265.39,0d0/)

double precision ::AlSihrks10(2)=(/-3143.78, 0.39297/) !hrks10=[-3143.78  0.39297] ;
double precision ::AlSihrks11(2)=(/0d0, 0d0/)
double precision ::AlSihrks12(2)=(/0d0, 0d0/)

double precision ::AlSihrks20(2)=(/113246.16, -47.55509/)!,h_factor=1d0  hrks20=[113246.16 -47.55509] 
double precision ::AlSihrks21(2)=(/0d0,0d0/)
double precision ::AlSihrks22(2)=(/0d0,0d0/)
!for PbSn

double precision ::PbSnhsa1(8)=(/-2977.961,93.949561,-24.5242231,-0.365895e-2,-0.24395d-6,0d0,-0.6019d-18,0d0/)
double precision ::PbSnhsa2(8)=(/-7650.085,101.700244,-24.5242231,-0.365895e-2, -0.24395e-6,0d0,0d0,0d0/)
double precision ::PbSnhsa3(8)=(/-7161.085,105.220244,-24.5242231,-0.365895e-2, -0.24395e-6,0d0,0d0,0d0/)
double precision ::PbSnhsb1(8)=(/1247.957,51.355548,-15.961,-0.188702e-1,0.3121167e-5,-61960.0,0.147031e-17,0d0/)
double precision ::PbSnhsb2(8)=(/-345.135,56.983315,-15.961,-0.188702e-1,0.3121167e-5,-61960.0,0d0,0d0/)
double precision ::PbSnhsb3(8)=(/-5855.135,65.443315,-15.961,-0.188702e-1, 0.3121167e-5,-61960.0,0d0,0d0/)
double precision ::PbSnhsrk1(8)=(/6204.5,-.67981,791.7,-1.5219,0d0,0d0,0d0,0d0/)
double precision ::PbSnhsrk2(8)=(/7145.3,-2.30237,0d0,0d0,0d0,0d0,0d0,0d0/)
double precision ::PbSnhsrk3(8)=(/19700.0,-15.89,0d0,0d0,0d0,0d0,0d0,0d0/)

!============
!!$  PbSnhsa(1, 1) = -2977.961 ! J/mol 
!!$  PbSnhsa(1, 2) = 93.949561 
!!$  PbSnhsa(1, 3) = -24.5242231 
!!$  PbSnhsa(1, 4) = -0.365895e-2 
!!$  PbSnhsa(1, 5) = -0.24395e-6 
!!$  PbSnhsa(1, 6) = 0. 
!!$  PbSnhsa(1, 7) = -0.6019e-18 
!!$  PbSnhsa(1, 8) = 0. 

!!$  PbSnhsa(2, 1) = -7650.085 
!!$  PbSnhsa(2, 2) = 101.700244 
!!$  PbSnhsa(2, 3) = -24.5242231 
!!$  PbSnhsa(2, 4) = -0.365895e-2 
!!$  PbSnhsa(2, 5) = -0.24395e-6 
!!$  PbSnhsa(2, 6) = 0. 
!!$  PbSnhsa(2, 7) = 0. 
!!$  PbSnhsa(2, 8) = 0. 
!!$
!!$  PbSnhsa(3, 1) = -7161.085 
!!$  PbSnhsa(3, 2) = 105.220244 
!!$  PbSnhsa(3, 3) = -24.5242231 
!!$  PbSnhsa(3, 4) = -0.365895e-2 
!!$  PbSnhsa(3, 5) = -0.24395e-6 
!!$  PbSnhsa(3, 6) = 0. 
!!$  PbSnhsa(3, 7) = 0. 
!!$  PbSnhsa(3, 8) = 0. 

!!$  PbSnhsb(1, 1) = 1247.957
!!$  PbSnhsb(1, 2) = 51.355548 
!!$  PbSnhsb(1, 3) = -15.961 
!!$  PbSnhsb(1, 4) = -0.188702e-1 
!!$  PbSnhsb(1, 5) = 0.3121167e-5 
!!$  PbSnhsb(1, 6) = -61960.0 
!!$  PbSnhsb(1, 7) = 0.147031e-17 
!!$  PbSnhsb(1, 8) = 0. 
!!$
!!$  PbSnhsb(2, 1) = -345.135 
!!$  PbSnhsb(2, 2) = 56.983315 
!!$  PbSnhsb(2, 3) = -15.961 
!!$  PbSnhsb(2, 4) = -0.188702e-1 
!!$  PbSnhsb(2, 5) = 0.3121167e-5 
!!$  PbSnhsb(2, 6) = -61960.0 
!!$  PbSnhsb(2, 7) = 0. 
!!$  PbSnhsb(2, 8) = 0. 

!!$  PbSnhsb(3, 1) = -5855.135 
!!$  PbSnhsb(3, 2) = 65.443315 
!!$  PbSnhsb(3, 3) = -15.961 
!!$  PbSnhsb(3, 4) = -0.188702e-1 
!!$  PbSnhsb(3, 5) = 0.3121167e-5 
!!$  PbSnhsb(3, 6) = -61960.0 
!!$  PbSnhsb(3, 7) = 0. 
!!$  PbSnhsb(3, 8) = 0. 

!!$  PbSnhsrk(1, 1) = 6204.5 
!!$  PbSnhsrk(1, 2) = -.67981 
!!$  PbSnhsrk(1, 3) = 791.7 
!!$  PbSnhsrk(1, 4) = -1.5219 
!!$  PbSnhsrk(1, 5) = 0. 
!!$  PbSnhsrk(1, 6) = 0. 
!!$  PbSnhsrk(1, 7) = 0. 
!!$  PbSnhsrk(1, 8) = 0. 

!!$
!!$  PbSnhsrk(2, 1) = 7145.3 
!!$  PbSnhsrk(2, 2) = -2.30237 
!!$  PbSnhsrk(2, 3) = 0. 
!!$  PbSnhsrk(2, 4) = 0. 
!!$  PbSnhsrk(2, 5) = 0. 
!!$  PbSnhsrk(2, 6) = 0. 
!!$  PbSnhsrk(2, 7) = 0. 
!!$  PbSnhsrk(2, 8) = 0. 
!!$
!!$  PbSnhsrk(3, 1) =6e3!19700.0 !6e3 
!!$  PbSnhsrk(3, 2) =0.! -15.89 
!!$  PbSnhsrk(3, 3) = 0. 
!!$  PbSnhsrk(3, 4) = 0. 
!!$  PbSnhsrk(3, 5) = 0. 
!!$  PbSnhsrk(3, 6) = 0. 
!!$  PbSnhsrk(3, 7) = 6e3 
!!$  PbSnhsrk(3, 8) = 0. 


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

