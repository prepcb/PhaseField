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
     write(61,*) 'Parameter list:'
     write(61,'(A13,F12.5)') 'anisotropy  ',  epsilon
     write(61,'(A13,F12.5,A5)') 'Delta T     ',  delta*g_T0,"K    "
     write(61,'(A13,F12.5,A5)') 'Ref Temp    ',  g_T0,"K    "
     write(61,'(A13,F12.5,A5)') 'T_M Ni      ',  g_T(1),"K    "
     write(61,'(A13,F12.5,A5)') 'T_M Cu      ',  g_T(2),"K    "
     write(61,'(A13,F12.5)')    'c ave       ',  g_c0
     write(61,'(A13,F12.5,A5)') 'char t      ',  g_tau*1d9,"ns   "
     write(61,'(A13,F12.5,A5)') 'char L      ',  sqrt(g_lengthSQ)*1e9,"nm   "
     write(61,'(A13,E12.5,A5)')'char diffuse',  g_Dch,"m^2/s"
     write(61,'(A13,E12.5,A5)') 'D_c liq     ',  g_D(1),"m^2/s"
     write(61,'(A13,E12.5,A5)') 'D_c sol     ',  g_D(2),"m^2/s"
     write(61,'(A13,F12.5)')    'D_c/d_phi   ', g_D(1)/g_Dch
     write(61,'(A13,F12.5)')    'Le= D_T/D_s ', g_D_tilde/g_D(1)*g_Dch
     write(61,'(A13,F12.5,A5)') 'Delta Ni    ',  g_delta(1)*1d9,"nm   "
     write(61,'(A13,F12.5,A5)') 'Delta Cu    ',  g_delta(2)*1d9,"nm   "
     write(61,'(A13,F12.5)')    'delta/d0    ', g_delta(1)/5d-10
     write(61,'(A13,F12.5,A5)') 'beta Ni     ',  g_beta(1),"m/K/s"
     write(61,'(A13,F12.5,A5)') 'beta Cu     ',  g_beta(2),"m/K/s"
     write(61,'(A13,F12.5)')    'nucleus rad' , nuc_radius
     write(61,'(A13,F12.5)')    'Domain size ', grid_xmax
     write(61,'(A13,F12.5,A5)') 'Finest dx   ', min_dx
     write(61,'(A13,I12,A5)')   'Block size  ', nxb
     write(61,'(A13,I12,A5)')   'total_vars  ', total_vars
     write(61,'(A13,L)')        'phi heating ',  heat_phi_term
     write(61,'(A13,L)')        'sol  heating',  heat_c_term
     write(61,'(A13,L)')        'grad heating',  Grad_DW_term
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

