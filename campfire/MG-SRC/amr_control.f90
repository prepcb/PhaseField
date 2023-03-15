!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!!  amr_mg_control
!!!!  * Shell for calling main time-stepping code
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!!  amr_mg_time_dep_control
!!!!  * Main time-stepping loop
!!!!  * Calls various application specific functions:
!!!!       -> app_pretimeloop_output
!!!!       -> app_output_perstep
!!!!       -> app_output_occasional
!!!!       -> app_close_files
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!!  get_input_parameters
!!!!  * Processes command line inputs to read checkpoint data
!!!!  * Calls application specific function
!!!!          app_read_checkpoint
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!!  timestep_adjust
!!!!  * Changes timestep size
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!!  nodebound
!!!!  * [unused] constraint on very small values in unk() arrays
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Chris Goodyer, May 2012                                              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine amr_mg_control(pid,noprocs)
! Wrapper for controlling whole solver process
  integer, intent(in) :: pid, noprocs
  call amr_mg_time_dep_control(pid,noprocs)

end subroutine amr_mg_control

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_mg_time_dep_control(pid,noprocs)
  use workspace
  use tree
  use physicaldata
  use paramesh_dimensions
  use time_dep_parameters
  use multigrid_parameters
  use paramesh_interfaces
  use checkpoint_parameters
  use generic_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, noprocs
  integer :: back_step_ct,time_step,i_step,i_step_total,ierr,level,iopt, lb, npts
  double precision defect,end_time,old_defect_1, rmsres(total_vars)
  integer failer, finegridblocks(2)
  double precision :: timeprev, timenew
  double precision local_r
  double precision :: wtime, wtimeprev,time_total
  integer :: i,j,k, g, G_lnblocks, G_count_nodetype, count_nodetype

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! C2F feature, Tim Yang, Jan 2014
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: C2F_turning_point(C2F_desired_level-C2F_start_level+1)
  integer :: C2F_count, C2F_notice
  C2F_turning_point(:) = -1
  C2F_count = 1
  C2F_notice = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  failer=0
  rmsres(:) = 0.0
  npts = 0

  end_time=time+simulation_time
  incset=0
  if(pid.eq.0.and.verbose.ge.2) write(6,*) 'end time is ',end_time

! Pre timestep output
  
  call app_pretimeloop_output(pid, noprocs)

  if (verbose.ge.2.and.chk_chkpt.gt.1) then
     if (pid.eq.0) write(6,*) '----------Loaded solution--------------'
     call app_output_perstep(pid, noprocs, chk_chkpt-1)
     if (pid.eq.0) write(6,*) '---------------------------------------'
  endif

  if(pid.eq.0.and.verbose.ge.2)print *,"Entered multigrid control"
  iopt=1
  back_step_ct=0
  allocate(has_children(1:maxblocks),stat=time_step)!Using time_step as dud var
  allocate(nodetype_copy(1:maxblocks),stat=time_step)!Using time_step as dud var
!  call amr_guardcell(pid,iopt,nguard)
  !has_children(:) = .false.
  
!  call pf_guardcell(pid,iopt,mg_max_level)
  call amr_mg_init()
  !!--tim modified, THIS IS WRONG when parallel, checkpoint file and adaptive grid are used.
  if(local_adaptation.eq.0)then
    call amr_multigrid_block_types(pid,noprocs,finegridblocks)
  endif

  !!--tim added for more information
  call MPI_Reduce(lnblocks, G_lnblocks, 2, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
  if(pid.eq.0) write(*,*) "total blocks from checkpoint file:", G_lnblocks
  count_nodetype = 0
  G_count_nodetype = 0
  do lb=1,lnblocks
    if(nodetype(lb).eq.1)then
      count_nodetype = count_nodetype + 1
    endif
  enddo
  call MPI_Reduce(count_nodetype, G_count_nodetype, 2, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
  if(pid.eq.0) write(*,*) "the number of nodetype 1 blocks are:", G_count_nodetype
  
!  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!  call amr_guardcell(pid,iopt,nguard)
!  call pf_guardcell(pid,iopt,mg_max_level)

  if(pid.eq.0.and.verbose.ge.1)print *,"Number of 'real' variables:",total_vars
  
  i_step_total=0
  if(pid.eq.0.and.verbose.ge.2)print *,"Starting timestepping"
!  call ceg_block_balancing(pid, noprocs)
!  call amr_multigrid_block_types(pid,noprocs,finegridblocks)
!  stop
  wtimeprev = MPI_Wtime()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! C2F feature, Tim Yang, Jan 2014
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(C2F_switch.eq.1)lrefine_max = C2F_start_level
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(output_rate.ge.100.and. multiphase .and. chk_chkpt.gt.1)chk_chkpt=chk_chkpt*100-99

  do time_step=chk_chkpt,max_time_iter+back_step_ct






     if(pid.eq.0.and.verbose.eq.1)write(*,*)"secs passed.....", MPI_Wtime()-start_secs




   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! C2F feature, Tim Yang, Jan 2014
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(C2F_switch.eq.1)then
      !!--tim
      !!-- user edit
      C2F_turning_point(1) = 5 !! after how many time_steps from start should we do C2F?
      if(C2F_desired_level-C2F_start_level.eq.2)C2F_turning_point(2) = 10

      !!--tim
      !!-- no need to change here
      if(pid.eq.0.and.verbose.gt.1)write(*,*)"C2F is enabled"
      if(verbose.gt.1)then
        do i=1,C2F_desired_level-C2F_start_level
          if(pid.eq.0)write(*,*) "C2F turning point",i, "is at time-step",chk_chkpt+C2F_turning_point(i)
        enddo
      endif
      if(time_step-back_step_ct.eq.chk_chkpt+C2F_turning_point(C2F_count))then
        lrefine_max = lrefine_max + 1
        if(pid.eq.0.and.verbose.gt.1)write(*,*)"C2F current refinement level is", lrefine_max
        if(lrefine_max.gt.C2F_desired_level)then
          if(pid.eq.0)write(*,*) "FATAL ERROR, current refinement is greater than desired"
          stop
        endif
        C2F_count = C2F_count + 1
        C2F_notice = 1
      endif
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! Little bit of code to ensure the last time step ends on simulation_time
     !     if((time+dt).gt.end_time)then
     !        dt=end_time-time
     !     end if
    dt= min(dt,max_dt)

    time=time+dt
     if(time.ge.end_time)exit

     if(pid.eq.0.and.verbose.ge.2)print *,"max_time_iter, back_step_ct:: ", max_time_iter,back_step_ct
     nodetype_copy(1:lnblocks)=nodetype(1:lnblocks)
     call amr_multigrid_child_set()
     if(pid.eq.0.and.verbose.ge.1)write(*,*)!Blank line between timesteps
     if(pid.eq.0.and.verbose.ge.0)write(*,*)"Starting time step no.",time_step-back_step_ct,".  Time:",time-dt," to ",time
     if(pid.eq.0.and.verbose.ge.1)write(*,*)"Current dt:",dt

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!! C2F feature, Tim Yang, Jan 2014
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     if(C2F_switch.eq.1.and.C2F_notice.eq.1)then
       if(pid.eq.0)write(*,*)"C2F notice: ------------------------------------------------------"
       if(pid.eq.0)write(*,*)"            refinement level has been increased to:", lrefine_max
       if(pid.eq.0)write(*,*)"            remaining refinement level to increase:", C2F_desired_level-lrefine_max
       if(pid.eq.0)write(*,*)"------------------------------------------------------------------"
       C2F_notice = 0
     endif
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

     call amr_mg_get_max_refinement(pid,noprocs)
     call reset_nodetype(noprocs,pid,mg_max_level)
    
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Adaptivity tests
     chk_restart = 0
     if(local_adaptation.ge.1.and.failer.eq.0.and.chk_restart.ne.1)then

         ! DO NOT DO THIS CALL IF grow=1
        !call amr_guardcell(pid,iopt,nguard)

        ! Just do main variable guardcells
        if(explicit)then
           gcell_on_cc(:)=.false.
           do current_var=1,total_vars
              current_unk=1+(current_var-1)*nunkvbles
              gcell_on_cc(current_unk)=.true.
              current_unk=2+(current_var-1)*nunkvbles
              gcell_on_cc(current_unk)=.true.
              current_unk=3+(current_var-1)*nunkvbles
              gcell_on_cc(current_unk)=.true.
           end do
        else
           gcell_on_cc(:)=.false.
           do current_var=1,total_vars
              current_unk=1+(current_var-1)*nunkvbles
              gcell_on_cc(current_unk)=.true.
              current_unk=4+(current_var-1)*nunkvbles
              gcell_on_cc(current_unk)=.true.
              current_unk=5+(current_var-1)*nunkvbles
              gcell_on_cc(current_unk)=.true.
           end do
        end if
        ! Loop over grids to match guardcells
        do g = lrefine_max, 2, -1
           call reset_nodetype(noprocs,pid,mg_max_level)
!           call reset_nodetype(noprocs,pid,g)
!           call pf_guardcell(pid,1,g)

           call pf_restrict(pid, 1, 0, .true., g)

           call reset_nodetype(noprocs, pid, g-1)
           call pf_guardcell(pid,iopt,g-1)
        end do
        call reset_nodetype(noprocs, pid, mg_max_level)

        if(pid.eq.0.and.verbose.ge.1)print *,pid,"Testing refinement"
        if(local_adaptation.eq.1) then
           ! Default adaptation
           call amr_test_refinement(pid, lrefine_min, lrefine_max)
        else
           ! User supplied adaptation
           call app_test_refinement(pid, lrefine_min, lrefine_max)
        endif
        if(pid.eq.0.and.verbose.ge.2) call system('date +"%k:%M:%S:%N"')
        if(pid.eq.0) wtimeprev = MPI_Wtime()
        if(pid.eq.0) then
           call cpu_time(timenew)
           timeprev=timenew
        endif
        
        call amr_refine_derefine
                        
        if(pid.eq.0) then
           wtime = MPI_Wtime()
           call cpu_time(timenew)
           if (verbose.ge.2) print*,'Sys/Real adapt time was ', timenew-timeprev, wtime-wtimeprev
           timeprev=timenew
           wtimeprev=wtime
        endif
        if(pid.eq.0.and.verbose.ge.2) call system('date +"%k:%M:%S:%N"')

!        call amr_multigrid_block_types(pid,noprocs,finegridblocks)
        call amr_prolong(pid,iopt,nguard)
     
        ! CEG try domain growing
        if (grow.eq.1) call amr_grow_domain(pid, noprocs, time_step-back_step_ct)


        !call amr_multigrid_block_types(pid,noprocs)
        if(pid.eq.0.and.verbose.ge.2)print *,"Reinitialising multigrid"
        call amr_mg_init() !NOTE FOR CHRIS, THIS IS NEW THING WHICH STOPPED SEG FAULT
        nodetype_copy(1:lnblocks)=nodetype(1:lnblocks)
        call amr_mg_get_max_refinement(pid,noprocs)
        call amr_multigrid_child_set()
        call reset_nodetype(noprocs,pid,mg_max_level)
        !call amr_mg_init() !NOTE FOR CHRIS, THIS IS OLD THING WHICH NEEDED TO BE MOVED UP
        call amr_multigrid_block_types(pid,noprocs,finegridblocks)
     else if (verbose.ge.3) then
        call amr_multigrid_block_types(pid,noprocs,finegridblocks)
     end if  ! end if test ref

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Copy existing soln(s) to rhs
     !!--tim added for more information
     call MPI_Reduce(lnblocks, G_lnblocks, 2, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
     if(pid.eq.0) write(*,*) "the number of total blocks are:", G_lnblocks
!!--tim     count_nodetype = 0
!!--tim     G_count_nodetype = 0
!!--tim     do lb=1,lnblocks
!!--tim       if(nodetype(lb).eq.1)then
!!--tim         count_nodetype = count_nodetype + 1
!!--tim       endif
!!--tim     enddo
!!--tim     call MPI_Reduce(count_nodetype, G_count_nodetype, 2, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
!!--tim     if(pid.eq.0) write(*,*) "the number of nodetype 1 blocks are:", G_count_nodetype

!     call amr_guardcell(pid,iopt,nguard)
     call pf_guardcell(pid,iopt,mg_max_level)
     if(explicit)then
        call amr_mg_runtime_rhs_explicit()
     else
        call amr_mg_runtime_rhs()
     end if
     if(pid.eq.0.and.verbose.ge.2)print *,"Refinement/guardcell update/RHS setting done, start v-cycles"

     defect=0.0
     old_defect_1=0.0



     do i_step=1,max_v_cycle_count
        i_step_total=i_step_total+1
        g_n_vcycle = g_n_vcycle + 1
        if(pid.eq.0)write(*,*)"g v cycle",g_n_vcycle 


        ! Back up copies of old defects
        old_defect_1=defect
        defect=0.0

        if (multigrid_on.ge.1) then
           call amr_multigrid_v_cycle(pid,noprocs,mg_max_level,defect, rmsres, npts)
        else
           call amr_non_linear_solve(pid,noprocs,mg_max_level,defect, rmsres, npts)
        endif

        call amr_vcycle_summary(pid, noprocs, i_step, defect, old_defect_1, npts, rmsres)

        if (shampine_test.eq.1) then 
           call shampine_convergence(i_step, failer, pid)
           if (defect.lt.defect_tol_min) exit
           if (failer.eq.1) exit           
        else
           ! Final check, if converged, stop cycles
           if(defect.lt.defect_tol_min)exit
        end if
  
        ! Start testing defect to see whether it's growing/too big
        if(defect.gt.defect_tol_max)then
           if(pid.eq.0.and.verbose.ge.1)print *,"Defect too large, aborting"
           exit
        end if
         if(i_step.ne.1.and.old_defect_1/defect.lt.1.01)then
         if(pid.eq.0.and.verbose.ge.1)print *,"Convergence rate",old_defect_1/defect," too small, aborting"
         exit
         endif
        !CEG changed 2 to 20
!        if(i_step.gt.20)then
           ! If defects grow, abort
!           if(defect.gt.old_defect_1)then
!              if(pid.eq.0.and.verbose.ge.1)print *,"Defects growing, aborting v-cycles"
!              exit
!           end if
!        end if
        if (i_step.eq.max_v_cycle_count)then
           ! Don't waste time on cycles going nowhere
           if(pid.eq.0.and.verbose.ge.1) print *,"Maximum V-cycles reached" 
        endif
     end do

!     if (time_step-back_step_ct.gt.2) then
!        call ceg_local_timestep(pid, noprocs, local_r)
!     else
        local_r = -1.0
!     endif

     ! Timer end code
     if(pid.eq.0) then
        wtime = MPI_Wtime()
        call cpu_time(timenew)
        if (verbose.ge.2) write(6,*)'Sys/Real solve time was ', timenew-timeprev, wtime-wtimeprev
!        write (133,*)time, wtime-wtimeprev, finegridblocks(1),finegridblocks(2), i_step
        timeprev=timenew
        wtimeprev=wtime
!        call flush(133)
        if (verbose.ge.4) write(*,*)"Vcycle: ",i_step,"Defect=",defect
        if (verbose.ge.1) then
           write(*,*)"Total V-cycles used:",i_step_total,"Avg V-cycles per timestep: ",(i_step_total*1.0)/(1.0*(time_step-back_step_ct))
        endif
        if (verbose.ge.1) write(*,*)''
!        write (unit=301,fmt=*)time,i_step,(i_step_total*1.0)/(1.0*(time_step-back_step_ct))
     endif

     failer=0


     call timestep_adjust(defect, local_r, failer, i_step, pid,old_defect_1/defect)
     if (failer.eq.1) then
        i_step_total=i_step_total-i_step
        back_step_ct=back_step_ct+1
        if (dt.lt.min_dt) then
           if(pid.eq.0) print*,'Tiny dt - aborting'
           call MPI_Abort(MPI_COMM_WORLD, failer, ierr)
        endif
     else ! i.e. if (failer.eq.0) then
     !!--tim output for successful time-step
!!--tim     call tim_visual_paraview(time_step)
            
        ! Output on every step
        time_total = MPI_WTIME()-start_secs
!        if(mod(time_step-back_step_ct,10).eq.0)then
           call app_output_perstep(pid, noprocs, time_step-back_step_ct)
!        endif
        if(mod(time_step-back_step_ct,output_rate).eq.0.or. time_total.gt.max_time) then
           !Output every output_rate steps
           if (time_total.gt.max_time)max_time=max_time+144d2
           call app_output_occasional(pid, noprocs, time_step-back_step_ct)
           if(time_total.gt.max_time)max_time=max_time + 144d2
           call cpu_time(timenew)
           if (pid.eq.0) then
              wtime = MPI_Wtime()
              if (verbose.ge.2) print*,'Output took ', timenew-timeprev, wtime-wtimeprev
              timeprev=timenew
           endif
           if(time_total.gt.max_time)exit
        endif






        ! Set to BDF2 now we have passed the first timestep
        if(explicit)then
           mode_1 = 1
        else
           mode_1 = 2
        endif
!        call nodebound
     endif
     chk_restart = 0

     if (time_step-back_step_ct.eq.multigrid_on.and.multigrid_on.gt.1) multigrid_on = 0



! #WARNING DEBUGGING get outs on
! DEBUGGING get outs
!     if (mode_1.eq.1.and.dt.gt.0.02) dt = 0.02
!     if (dt.gt.0.1) dt = 0.1

!     if (mode_1.eq.2) return
!     if (mode_1.eq.2) exit
!     if (time_step-back_step_ct.eq.90) return!585
!     if (time_step-back_step_ct.eq.100) exit!585
!     if (time_step-back_step_ct.eq.10000) return!585


!max_dt adjust
!!$  if(mod(time_step-back_step_ct,10).eq.0)then
!!$     if(pid.eq.0)write(*,*)"vcycle ratio =", g_n_vcycle
!!$!     if(dt.eq.max_dt.and.1d2/real(g_n_vcycle).lt.0.7)then
!!$     if(1d2/real(g_n_vcycle).lt.0.5)then
!!$        max_dt=max_dt*0.9
!!$     elseif(dt.eq.max_dt.and.1d2/real(g_n_vcycle).gt.0.70)then
!!$        max_dt=max_dt*1.2
!!$     endif
!!$     g_n_vcycle=0
!!$  endif
    
  end do !end main time loop

! Close all output text files

  call app_close_files(pid, noprocs)

!  if(pid.eq.0)print *,"Mode: ",mode_1
!  if(pid.eq.0)print *,"Min level: ",mg_min_lvl
!  call amr_multigrid_block_types(pid,noprocs,finegridblocks)

  deallocate(has_children)
  deallocate(nodetype_copy)

! print out the final grid state
!  call amr_plotfile_chombo(time_step-back_step_ct)
!  call amr_checkpoint_wr(time_step-back_step_ct, user_attr_1=dt, user_attr_2=time, user_attr_3=dtold)
!  write(61,*) 'Step, time', time_step-back_step_ct, time

!101 format(e19.13)

end subroutine amr_mg_time_dep_control

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_input_parameters(pid, noprocs)
  use multigrid_parameters
  use time_dep_parameters  
  use generic_parameters
  use tree                 
  use paramesh_dimensions  
  use paramesh_interfaces
  use checkpoint_parameters

  implicit none
  include 'mpif.h'
  integer,intent(in) :: pid, noprocs
  integer :: ierr
  !BUFFER holds the command line arguments
  CHARACTER *100 BUFFER  
  integer, external     :: iargc
  logical :: there
  integer lb
  if (IArgC().gt.0) then            ! Do we have command line arguments?
     if (IArgC().eq.2) then         ! Are we expecting to specify which checkpoint to load
        call getarg(1,BUFFER)
        read(BUFFER,*) chk_restart
        if (chk_restart.eq.1) then  ! Have we done it right?
           call getarg(2,BUFFER)
           read(BUFFER,*) chk_chkpt
        else
           if (pid.eq.0.and.verbose.ge.0) then
              write(6,*) 'ERROR: Command line usage is:'
              write(6,*) './executable          for new run'
              write(6,*) './executable 1 <n>    to load checkpoint <n>'
              write(6,*) './executable 2        to load latest checkpoint, saved in CHK.out'
           endif
           stop
        endif
     else if (IArgC().eq.1) then
        if(pid.eq.0) then
           inquire(file="CHK.out",exist=there)
           if (there) then
              open (unit=99,file="CHK.out")
              read(99,*) chk_chkpt
              close(99)
           else
              if (verbose.ge.1) write(6,*) 'File CHK.out not found.  Starting new simulation.'
              chk_chkpt=0
           endif
        endif
        call MPI_Bcast(chk_chkpt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     else
        if (pid.eq.0.and.verbose.ge.0) then
           write(6,*) 'ERROR: Command line usage is:'
           write(6,*) './executable          for new run'
           write(6,*) './executable 1 <n>    to load checkpoint <n>'
           write(6,*) './executable 2        to load latest checkpoint, saved in CHK.out'
        endif
        stop
     endif

     if (chk_chkpt.gt.0) then
        chk_checkf = "hdf5"

        call app_read_checkpoint(pid, noprocs)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! C2F feature, Tim Yang, Jan 2014
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(C2F_switch.eq.1)then
          if(pid.eq.0)write(*,*) "C2F notice: ------------------------------------------------------"
          if(pid.eq.0)write(*,*) "current refinement level from input file: ", lrefine(lnblocks)
          if(pid.eq.0)write(*,*) "desired refinement level to reach       : ", C2F_desired_level
          if(pid.eq.0)write(*,*) "level difference is : ", C2F_desired_level-C2F_start_level
          if(pid.eq.0)write(*,*) "------------------------------------------------------------------"
          if(pid.eq.0.and.C2F_desired_level-C2F_start_level.gt.1)then 
            if(pid.eq.0)write(*,*) "Error, so far the code only does C2F upto 1 level."
            stop
          endif
          call amr_relocate_meshing(pid, noprocs)
!!--tim          call amr_mg_init()
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        call amr_guardcell(pid,1,1)

!        call amr_refine_derefine
!        if(pid.eq.0) write(6,*) 'Checkpoint load reordered again'
        time = chk_t
        dt = chk_dt
        dtold = chk_dtold
        if (verbose.ge.3) call graphviz_output(pid, noprocs, chk_chkpt)

        if (pid.eq.0.and.verbose.ge.2 ) then  
           write(6,'(A)') 'Global domain limits'
           write(6,'(F9.2,X,F9.2,X,F9.2,X,F9.2,X,F9.2,X,F9.2,X)') grid_xmin,grid_xmax,grid_ymin,grid_ymax,grid_zmin,grid_zmax
        endif

!        call amr_multigrid_block_types(pid, noprocs, finegridblocks)
!        call amr_multigrid_block_timer(pid,noprocs)
        chk_chkpt=chk_chkpt+1

        mode_1 = 2
     endif
  endif

  call mpi_amr_global_domain_limits() ! This routine sets grid_xmax, grid_xmin, ...
  call set_domain_limits()            ! This is a local routine to turn output from mpi_amr_global_domain_limits into boundary_box array

  ! Initialise everything
  refine(1:lnblocks) = .false.
  call amr_refine_derefine


  if (pid.eq.0.and.verbose.ge.2) then  
     write(6,'(A)') 'Global domain limits'
     write(6,'(F9.2,X,F9.2,X,F9.2,X,F9.2,X,F9.2,X,F9.2,X)') grid_xmin,grid_xmax,grid_ymin,grid_ymax,grid_zmin,grid_zmax
  endif

  ! For new jobs do initial adaptation
  if (chk_chkpt.eq.0) then
     if(pid.eq.0.and.verbose.ge.1) write(6,*) ' starting new job'
     chk_chkpt=1
     call amr_initial_meshing(pid, noprocs)
  endif

end subroutine get_input_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine timestep_adjust(defect, local_r, failer, i_step, pid,conv_rate)
  use generic_parameters
  use workspace
  use tree
  use physicaldata
  use paramesh_dimensions
  use time_dep_parameters
  use multigrid_parameters
  use paramesh_interfaces!, only :amr_guardcell, & 
       !&                         amr_prolong, & 
       !&                         amr_restrict
  implicit none
  integer, intent(in) :: i_step, pid
  integer, intent(inout) :: failer
  double precision, intent(inout) :: defect, local_r,conv_rate

! Not done Local error test
  if (local_r.lt.0.0) then
     ! James's method:
     !   In all failure cases dt->.75*dt and timestep is restarted (failure is lack of convergence or NANs)
     !   If converged in small (<4) cycles dt->10/9*dt then next cycle
     !   If converged in large (>8) cycles dt->9/10*dt then next cycle (note, not reset since it converged)
     if(defect.eq.0.0)then
        ! If it's = 0, then in all likelyhood either everything's melted, or
        ! the grid's full of NANs
        if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[37;1;48m"
        if(pid.eq.0.and.verbose.ge.1)print *,"Zero defect - let's ignore for now - treat as tiny"
        if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[0m"
        if (MOD(incset,stepsize_increase_frequency).eq.0) then
           if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[38;1;47m"
           if(pid.eq.0.and.verbose.ge.1)print *,"Low v-cycle count, increasing dt"
           if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[0m"
           dtold=dt
           dt=tstep_increase_factor*dt
           if(dt.gt.max_dt)then
              dt=max_dt
           end if
!           incset=1
           incset=1
        else
           incset = incset+1
        endif

!        if(pid.eq.0)print *,"Zero defect, asumming failure"
!        call amr_mg_time_step_reset()
!        back_step_ct=back_step_ct+1
!        time=time-dt
!        dt=dt*.75
!        incset = 0
!        failer=1
     else if(defect.gt.defect_too_big)then
        ! Pretty obvious really
        if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[39;1;46m"
        if(pid.eq.0.and.verbose.ge.1)print *,"Large defect, asumming failure"
        if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[0m"
        call amr_mg_time_step_reset()
        time=time-dt
        dt=dt*tstep_decrease_factor
        incset = 0
        failer=1
     else if(i_step.le.low_vcycle_count.and.(defect.lt.defect_tol_min))then

!!!???
        if (MOD(incset,stepsize_increase_frequency).eq.0) then
           if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[40;1;45m"
           if(pid.eq.0.and.verbose.ge.1)print *,"Low v-cycle count, increasing dt"
           if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[0m"
           dtold=dt
           if(conv_rate.gt.3d1.or.conv_rate.eq.0d0)dt=tstep_increase_factor*dt
           if(pid.eq.0.and.verbose.ge.1)print *,"dt=",dt
           if(dt.gt.max_dt)then
              dt=max_dt
           end if
           incset = 1
        else
           incset = incset+1
        endif
     else if((i_step.ge.max_v_cycle_count).or.(defect.gt.defect_tol_min.and.shampine_test.eq.0))then
!     else if((i_step.ge.10).or.(defect.gt.defect_tol))then
!         if(pid.eq.0)print *,"Timestep decreased removed"
        if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[41;1;44m"
        if(pid.eq.0.and.verbose.ge.1)print *,"Excessive v-cycle count/lack of convergence, decreasing dt and restarting time step"
        if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[0m"
        call amr_mg_time_step_reset()
        time=time-dt
        dt=dt*tstep_decrease_factor
        !max_dt=dt
        incset = 0
        failer=1
     else if(i_step.gt.high_vcycle_count)then
        if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[42;1;43m"
        if(pid.eq.0.and.verbose.ge.1)print *,"High v-cycle count, decreasing dt"   
        if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[0m"
        dtold=dt
        dt=tstep_decrease_factor*dt
        incset = 0
     end if
  else     ! Local error estimation techniques
     if(defect.gt.defect_too_big.or.(i_step.ge.max_v_cycle_count).or.(defect.gt.defect_too_big.and.shampine_test.eq.0))then
        ! Pretty obvious really
        if(pid.eq.0.and.verbose.ge.1)print *,"LEE: Large defect, asumming failure", defect
        call amr_mg_time_step_reset()
        time=time-dt
        dt=dt*0.5
        incset = 0
        failer=1
     else if (local_r.lt.lee_tol_failing) then
        if(pid.eq.0.and.verbose.ge.1)print *,"LEE: Tiny local_r, asumming failure", local_r
        call amr_mg_time_step_reset()
        time=time-dt
        dt=dt*0.5
        incset = 0
        failer=1
     else if (local_r.lt.lee_tol_halving) then
        if(pid.eq.0.and.verbose.ge.1)print *,"LEE: Small local_r, halving stepsize", local_r
        dtold=dt
        dt=dt*0.5
        incset = 0
     else if (local_r.lt.lee_tol_decrease) then
        if(pid.eq.0.and.verbose.ge.1)print *,"LEE: local_r requests decreased timestep", local_r
        dtold=dt
        dt=dt*local_r
        incset = 0
     else if (local_r.gt.lee_tol_increase.and.dt.lt.max_dt) then
        if (MOD(incset,stepsize_increase_frequency).eq.0) then
           if(pid.eq.0.and.verbose.ge.1)print *,"LEE: local_r requests increased timestep", local_r
           dtold=dt
           dt=dt*2.0
           if (dt.gt.max_dt) dt = max_dt
           incset = 1
        else
           incset = incset + 1
        endif
     endif
  endif

end subroutine timestep_adjust


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nodebound
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none

  integer i, j, k, lb, v


!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb,v,k,j,i)
!!!$OMP DO SCHEDULE(DYNAMIC,8) 
  do lb=1,lnblocks
     if(nodetype(lb).eq.1) then
        do k=kl_bnd+1,ku_bnd-1
           do j=jl_bnd+1,ju_bnd-1
              do i=il_bnd+1,iu_bnd-1
                 do v = 1, nvar
                    if (ABS(unk(v,i,j,k,lb)).lt.1.0e-15) then
!                       print *,v,i,j,k,lb,unk(v,i,j,k,lb), ABS(unk(v,i,j,k,lb)
                       unk(v,i,j,k,lb)=0.0
                    endif
                 end do
              end do
           end do
        end do
     end if
   end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 

end subroutine nodebound

