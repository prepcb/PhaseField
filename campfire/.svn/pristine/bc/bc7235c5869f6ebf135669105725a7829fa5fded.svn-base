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
     write(61,*) 'Mc-Inf      ',  mcinf
     write(61,*) 'kappa_E     ',  ke
     write(61,*) 'epsilon     ',  epsilon
     write(61,*) 'lambda      ',  lambda
     write(61,*) 'D_solute    ', D_solute
     write(61,*) 'D_therm     ', D_therm
     write(61,*) 'Lewis no.   ',  le
     write(61,*) 'Delta       ',  delta
     write(61,*) 'nucleus rad ', nuc_radius
!     write(61,*) 'End time    ',  end_time
     write(61,*) 'Domain size ', grid_xmax
     write(61,*) 'Finest dx   ', min_dx
     write(61,*) 'Block size  ', nxb
     write(61,*) 'total_vars  ', total_vars
     write(61,*) 'Solute on?  ', solute
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
     print *,"Lambda: ",lambda,"Dtherm",D_therm,"Dsolute",D_solute
     if (chk_chkpt.eq.1) then
        open (unit=125,file="tip_rad.txt")
        open (unit=126,file="tip_loc.txt")
        open (unit=133,file="real_time.txt")
        !open (unit=301,file="Vcycnt.txt")
     else
        open (unit=125,file="tip_rad.txt",position='append')
        open (unit=126,file="tip_loc.txt",position='append')
        open (unit=133,file="real_time.txt",position='append')
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

