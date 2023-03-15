!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_read_checkpoint                                          REQUIRED
!!!!  * Calls amr_checkpoint_re with application specific extra parameters
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_pretimeloop_output                                       REQUIRED
!!!!  * Any processing or output needed before first timestep
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_close_files                                              REQUIRED
!!!!  * Closes any files opened in application section
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Latest version: Chris Goodyer, May 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_read_checkpoint(pid, noprocs)
  use paramesh_interfaces
  use time_dep_parameters
  use checkpoint_parameters
  implicit none
  integer,intent(in) :: pid, noprocs

  if(pid.eq.0) write(6,*) 'Checkpoint on : ', chk_chkpt
  call amr_checkpoint_re(chk_chkpt, user_attr_1=chk_dt, user_attr_2=chk_t, user_attr_3=chk_dtold)
  if(pid.eq.0) write(6,*) 'Checkpoint read with ', chk_t, chk_dt, chk_dtold

end subroutine app_read_checkpoint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_pretimeloop_output(pid, noprocs)
  use multigrid_parameters  ! for min_dx, total_vars
  use solution_parameters   ! for most of the parameters
  use tree                  ! for grid_xmax
  use paramesh_dimensions   ! for nxb
  use checkpoint_parameters ! for chk_chkpt
  implicit none
  integer,intent(in) :: pid, noprocs

  !CEG Add README output
  if(pid.eq.0.and.chk_chkpt.eq.1) then
     open (unit=61,file="README")
     write(61,*) 'Parameter list: Pb Sn'
     write(61,'(A13,F12.5)') 'anisotropy  ',  epsilon
     write(61,'(A13,F12.5,A5)') 'T           ',  g_T0+delta*g_T0,"K    "
     write(61,'(A13,F12.5,A5)') 'T_M A       ',  g_T(1),"K    "
     write(61,'(A13,F12.5,A5)') 'T_M B       ',  g_T(2),"K    "
     write(61,'(A13,F12.5)')    'c ave       ',  g_c0
     write(61,'(A13,F12.5,A5)') 'char t      ',  g_tau*1d9,"ns   "
     write(61,'(A13,F12.5,A5)') 'char L      ',  sqrt(g_lengthSQ)*1e9,"nm   "
     write(61,'(A13,F12.5,A5)') 'char V      ',  sqrt(g_lengthSQ)/g_tau,"m/s   "
     write(61,'(A13,E12.5,A5)')'char diffuse',  g_Dch,"m^2/s"
     write(61,'(A13,E12.5,A5)') 'D_c liq     ',  g_D(1),"m^2/s"
     write(61,'(A13,E12.5,A5)') 'D_c sol     ',  g_D(2),"m^2/s"
     write(61,'(A13,F12.5)')    'D_c/D_phi   ', g_D(1)/g_Dch
     write(61,'(A13,F12.5)')    'Le= D_T/D_s ', g_D_tilde/g_D(1)*g_Dch
     write(61,'(A13,F12.5,A5)') 'Delta A     ',  g_delta(1)*1d9,"nm   "
     write(61,'(A13,F12.5,A5)') 'Delta B     ',  g_delta(2)*1d9,"nm   "
     write(61,'(A13,F12.5)')    'delta/d0    ', g_delta(1)/5d-10
     write(61,'(A13,F12.5,A5)') 'mu A        ',  g_mu(1),"m/K/s"
     write(61,'(A13,F12.5,A5)') 'mu B        ',  g_mu(2),"m/K/s"
     write(61,'(A13,F12.5,A5)') 'sigma A     ',  g_sigma(1),"J/m^2"
     write(61,'(A13,F12.5,A5)') 'sigma B     ',  g_sigma(2),"J/m^2"
     write(61,'(A13,F12.5)')    'nucleus rad' , nuc_radius
     write(61,'(A13,F12.5)')    'Domain size ', grid_xmax
     write(61,'(A13,F12.5)')    'Eutect width', width
     write(61,'(A13,F12.5,A5)') 'Finest dx   ', min_dx
     write(61,'(A13,F12.5,A5)') 'alpha - psi ', g_alpha
     write(61,'(A13,F12.5,A5)') 'Defect tol  ', defect_tol_min*1d7,"E-7"
     write(61,'(A13,I12,A5)')   'Block size  ', nxb
     write(61,'(A13,I12,A5)')   'total_vars  ', total_vars
     write(61,'(A13,L)')        'phi heating ',  heat_phi_term
     write(61,'(A13,L)')        'sol  heating',  heat_c_term
     write(61,'(A13,L)')        'grad heating',  Grad_DW_term
     write(61,'(A13,L)')        'cross term  ',  cross_term
     write(61,'(A13,L)')        'thermal run ',  thermal
     write(61,'(A13,L)')        'Folch-Plapp ',  folch
     write(61,'(A13,L)')        '1D          ',  g_1D
     write(61,'(A13,I2)')        'sym:1,-1,0  ',  g_sym
     write(61,'(A13,I2)')        'max level   ',  g_max_lvl
     write(61,'(A13,F12.5)')    'AntiTrapping',  g_beta
     write(61,'(A13,F12.5)')    'max time    ',  g_max_time
     write(61,'(A13,I2)')    'smooths     ',  smooth_count_max
     write(61,'(A13,I2)')    'output rate ',  g_output_rate
     write(61,'(A13,F12.5)')    'max dt    ',  g_max_dt
     write(61,*) ' '
     flush(61)
  endif

!   if(pid.eq.0) then
!      print *,"Lambda: ",lambda,"Dtherm",D_therm,"Dsolute",D_solute
!      open (unit=125,file="tip_rad.txt")
!      open (unit=126,file="tip_loc.txt")
!      open (unit=133,file="real_time.txt")
!      !open (unit=301,file="Vcycnt.txt")
!   endif
  if(pid.eq.0) then
     if (chk_chkpt.eq.1) then
!        open (unit=125,file="tip_rad.txt")
        open (unit=126,file="tip_loc.txt")
!        open (unit=133,file="real_time.txt")
        !open (unit=301,file="Vcycnt.txt")
     else
!        open (unit=125,file="tip_rad.txt",position='append')
        open (unit=126,file="tip_loc.txt",position='append')
!        open (unit=133,file="real_time.txt",position='append')
     endif
  endif

end subroutine app_pretimeloop_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_close_files(pid, noprocs)
  implicit none
  integer,intent(in) :: pid, noprocs

  if(pid.eq.0) then
     close(61)
     close(125)
     close(126)
     close(133)
  endif

end subroutine app_close_files  

