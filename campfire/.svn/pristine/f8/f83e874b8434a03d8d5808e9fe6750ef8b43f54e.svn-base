!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_mg_runtime_rhs
!!!!  * Sets up right hand side variable per unknown using previous timesteps
!!!!  * Also predicts initial guess at solution for new timestep
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_mg_time_step_reset
!!!!  * Rolls back to start of previous timestep if last attempt failed
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_multigrid_v_cycle
!!!!  * Recursive 2-grid V-cycle call
!!!!  * Calls user supplied function
!!!!             -> app_mg_smooth
!!!! 
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_mg_get_defect
!!!!  * Calculates defect contribution and the maximum defect and RMS residual
!!!!  * Calls user supplied function
!!!!             -> app_mg_get_defect_var
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_multigrid_child_set
!!!!  * Finds out which blocks have children and sets has_children() array
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  reset_nodetype
!!!!  * This resets nodetypes so that MLAT can work alongside guardcell call
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_mg_get_max_refinement
!!!!  * Finds maximum refinement level on a processor
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_mg_get_defect
!!!!  * This is the alternative call to amr_multigrid_v_cycle if multigrid
!!!!    is not being used
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Chris Goodyer, September 2012                                        !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine mp_convection_MG(i,j,k,lb,phi_est,c_est)
  use time_dep_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree

  implicit none

integer, intent(in) :: i,j,k, lb
double precision, intent(out):: phi_est, c_est
double precision :: dx,xi,xj,xk
double precision :: phix,phiy,phiz,norm,norm0,phi_t 
double precision ::rfactor,rf4,rf5,rf6
 dx = bsize(1,lb)/real(nxb)
 phix = 0.5*(unk(1+0*nunkvbles,i+1,j,k,lb) - unk(1+0*nunkvbles,i-1,j,k,lb))/dx
 phiy = 0.5*(unk(1+0*nunkvbles,i,j+1,k,lb) - unk(1+0*nunkvbles,i,j-1,k,lb))/dx
 phiz = 0.5*(unk(1+0*nunkvbles,i,j,k+1,lb) - unk(1+0*nunkvbles,i,j,k-1,lb))/dx
 norm = phix**2 + phiy**2+phiz**2+1d-3
 norm0 = sqrt(phix**2 + phiy**2+phiz**2)

 if (mode_1.eq.1) then
    phi_t = (unk(1,i,j,1,lb)-unk(2,i,j,1,lb))/dt
 elseif(mode_1.eq.2)then
    rfactor = dt/dtold
    rf4 = rfactor*rfactor/(rfactor+1.0)
    rf5 = (rfactor+1.0)
    rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
    phi_t = (rf4*unk(5,i,j,1,lb)-rf5*unk(4,i,j,1,lb)+rf6*unk(1,i,j,1,lb))/dt
 endif

 xi = phi_t*dt/dx*phix/norm
 xj = phi_t*dt/dx*phiy/norm
 xk = phi_t*dt/dx*phiz/norm

 phi_est = (1d0 - abs(xi))*(1d0 - abs(xj))*(1d0 - abs(xk))*unk(1+0*nunkvbles,i,j,k,lb)
 phi_est = phi_est + (abs(xi))*(1d0 - abs(xj))*(1d0 - abs(xk))*unk(1+0*nunkvbles,i+int(sign(1d0,xi)),j,k,lb)
 phi_est = phi_est + (abs(xj))*(1d0 - abs(xk))*(1d0 - abs(xi))*unk(1+0*nunkvbles,i,j+int(sign(1d0,xj)),k,lb)
 phi_est = phi_est + (abs(xk))*(1d0 - abs(xi))*(1d0 - abs(xj))*unk(1+0*nunkvbles,i,j,k+int(sign(1d0,xk)),lb)
 phi_est = phi_est + abs(xi)*abs(xj)*(1-abs(xk))*unk(1+0*nunkvbles,i+int(sign(1d0,xi)),j+int(sign(1d0,xj)),k,lb)
 phi_est = phi_est + abs(xj)*abs(xk)*(1-abs(xi))*unk(1+0*nunkvbles,i,j+int(sign(1d0,xj)),k+int(sign(1d0,xk)),lb)
 phi_est = phi_est + abs(xk)*abs(xi)*(1-abs(xj))*unk(1+0*nunkvbles,i+int(sign(1d0,xi)),j,k+int(sign(1d0,xk)),lb)
 phi_est = phi_est + abs(xk)*abs(xi)*abs(xj)*unk(1+0*nunkvbles,i+int(sign(1d0,xi)),j+int(sign(1d0,xj)),k+int(sign(1d0,xk)),lb)
 
 
 c_est = (1d0 - abs(xi))*(1d0 - abs(xj))*(1d0 - abs(xk))*unk(1+1*nunkvbles,i,j,k,lb)
 c_est = c_est + (abs(xi))*(1d0 - abs(xj))*(1d0 - abs(xk))*unk(1+1*nunkvbles,i+int(sign(1d0,xi)),j,k,lb)
 c_est = c_est + (abs(xj))*(1d0 - abs(xk))*(1d0 - abs(xi))*unk(1+1*nunkvbles,i,j+int(sign(1d0,xj)),k,lb)
 c_est = c_est + (abs(xk))*(1d0 - abs(xi))*(1d0 - abs(xj))*unk(1+1*nunkvbles,i,j,k+int(sign(1d0,xk)),lb)
 c_est = c_est + abs(xi)*abs(xj)*(1-abs(xk))*unk(1+1*nunkvbles,i+int(sign(1d0,xi)),j+int(sign(1d0,xj)),k,lb)
 c_est = c_est + abs(xj)*abs(xk)*(1-abs(xi))*unk(1+1*nunkvbles,i,j+int(sign(1d0,xj)),k+int(sign(1d0,xk)),lb)
 c_est = c_est + abs(xk)*abs(xi)*(1-abs(xj))*unk(1+1*nunkvbles,i+int(sign(1d0,xi)),j,k+int(sign(1d0,xk)),lb)
 c_est = c_est + abs(xk)*abs(xi)*abs(xj)*unk(1+1*nunkvbles,i+int(sign(1d0,xi)),j+int(sign(1d0,xj)),k+int(sign(1d0,xk)),lb)


end subroutine mp_convection_MG


subroutine amr_mg_runtime_rhs()
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  use multigrid_parameters 
  implicit none
  integer :: lb,i,j,k,v,u
  double precision :: rfactor, rf1, rf2, rf3
  double precision :: phi_est,c_est
 
  rfactor = dt/dtold
  rf1 = (rfactor+1.0)*dt/(2.0*rfactor+1.0)
  rf2 = (rfactor+1.0)*(rfactor+1.0)/(2.0*rfactor+1.0)
  rf3 = rfactor*rfactor/(2.0*rfactor+1.0)
!  rf1=dt
!  rf2=1.0
!  rf3=0.0

  ! Note this could be done more simply, however this is left as is for additional functionality
  if(lnblocks.gt.0) then
!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb,i,j,k,u,v)
!!!$OMP DO SCHEDULE(DYNAMIC,1)
     do lb=1,lnblocks
        !if(nodetype(lb).eq.1) then
!           if(nvar.gt.0) then
              do k=kl_bnd,ku_bnd
                 do j=jl_bnd,ju_bnd
                    do i=il_bnd,iu_bnd
                       
                       if(convection_predictor)call mp_convection_MG(i,j,k,lb,phi_est,c_est)
                       do v=1,total_vars
                     
                          ! set the start point for the real variable and its workspace data
                          u=1+(v-1)*nunkvbles
                         
                          if(mode_1.eq.1)then
                                     
                             unk(u+1,i,j,k,lb)=unk(u,i,j,k,lb)
! back up previous timesteps
                             unk(u+4,i,j,k,lb)=unk(u+3,i,j,k,lb)
                             unk(u+3,i,j,k,lb)=unk(u,i,j,k,lb)
! set initial guess as first order prediction
                             unk(u,i,j,k,lb)=unk(u+3,i,j,k,lb) + &
                                rfactor*(unk(u+3,i,j,k,lb)-unk(u+4,i,j,k,lb))
                          else if(mode_1.eq.2) then
! back up previous timesteps
                             unk(u+4,i,j,k,lb)=unk(u+3,i,j,k,lb)
                             unk(u+3,i,j,k,lb)=unk(u,i,j,k,lb)
! set up RHS based on previous solutions
!                             unk(u+1,i,j,k,lb)=unk(u,i,j,k,lb)
                             unk(u+1,i,j,k,lb)= rf2*unk(u+3,i,j,k,lb) - rf3*unk(u+4,i,j,k,lb)
                             if(convection_predictor) then
                                if(v.eq.1)then
                                   unk(u+1,i,j,k,lb) = phi_est
                                elseif(v.eq.2)then
                                   unk(u+1,i,j,k,lb) = c_est
                                endif
                             else
! set initial guess as first order prediction
                             unk(u,i,j,k,lb)=unk(u+3,i,j,k,lb) + &
                                rfactor*(unk(u+3,i,j,k,lb)-unk(u+4,i,j,k,lb))
                             endif
!                             unk(u,i,j,k,lb)=unk(u+3,i,j,k,lb) 
!                             unk(u,i,j,k,lb)=0.8*unk(u+3,i,j,k,lb)+0.2*unk(u+4,i,j,k,lb)

!                          else
!                             unk(u+1,i,j,k,lb)=unk(u,i,j,k,lb)/dt
                          
                          end if

                       end do
                    end do
                 end do
              end do
           !end if
        !end if
     end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 
  end if

  if (shampine_test.eq.1) call shampine_store()


end subroutine amr_mg_runtime_rhs

subroutine amr_mg_runtime_rhs_explicit()
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  implicit none
  integer :: lb,i,j,k,v,u



  if(lnblocks.gt.0) then
     do lb=1,lnblocks
              do k=kl_bnd,ku_bnd
                 do j=jl_bnd,ju_bnd
                    do i=il_bnd,iu_bnd
                       do v=1,total_vars
                          u=1+(v-1)*nunkvbles
                          unk(u+2,i,j,k,lb)=unk(u+1,i,j,k,lb)
                          unk(u+1,i,j,k,lb)=unk(u,i,j,k,lb)
                       end do
                    end do
                 end do
              end do
     end do
  end if




end subroutine amr_mg_runtime_rhs_explicit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_mg_time_step_reset()
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use tree
!  use solution_parameters
  implicit none
  integer :: lb,i,j,k,u,v
  ! Note this could be done simpler, however this is left as is for additional functionality
  if(lnblocks.gt.0) then
!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb,i,j,k,u,v)
!!!$OMP DO SCHEDULE(DYNAMIC,1)
     do lb=1,lnblocks
        !if(nvar.gt.0) then
           do k=kl_bnd,ku_bnd
              do j=jl_bnd,ju_bnd
                 do i=il_bnd,iu_bnd
                    do v=1,total_vars
                       ! set the start point for the real variable
                       u=1+(v-1)*nunkvbles
                       if(mode_1.eq.1)then
                          unk(u,i,j,k,lb)=unk(u+1,i,j,k,lb)
                       else if(mode_1.eq.2)then
                          unk(u,i,j,k,lb)=unk(u+3,i,j,k,lb)
                          unk(u+3,i,j,k,lb)=unk(u+4,i,j,k,lb)
                       end if
                    end do
                 end do
              end do
           end do
        !end if
     end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 
  end if

end subroutine amr_mg_time_step_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine amr_multigrid_v_cycle(pid,noprocs,level,global_max_defect, rmsres, npts)
!
  use paramesh_interfaces
  use multigrid_parameters
  use generic_parameters
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use tree
  use workspace
  use paramesh_comm_data ! Needed for definition of amr_mpi_real
  implicit none  
  include 'mpif.h'
  integer, intent(in) :: pid,noprocs, level
  integer, intent(inout) :: npts
  double precision, intent(inout) :: global_max_defect, rmsres(total_vars)
  integer :: smooth_count,lb,nlayers,iopt,ierr
  double precision :: local_max_defect, global_total, smooth_count_max_
  integer :: i,j,k,new_level
  integer :: lbs
  double precision :: timing_beg, timing_end, timing_beg2
  double precision :: timing_presmooth,timing_presmoothloop, timing_restrict,timing_loopingcode,timing_loopingcode2
  double precision :: timing_postsmooth,timing_postsmoothloop





  ! Plan:
  ! Need to cycle over all the leaf (node) blocks smoothing
  ! Next job is to find the defect
  ! Now do the FAS bit...
  ! Restrict if leaf nodes are not the coarsest level
  ! Possibly have to swap lots of data about here, 
  !    probably put this in it's own sub for my sanity
  ! Call solver again (recursive)
  ! More swapping? If so own sub for sanity
  ! Refine back
  ! Reverse FAS bit
  ! Smooth again
  ! Find defect/residual
  ! Setup a couple of parameters
  iopt = 1

  nlayers = nguard

  !!! ------------- added by MHC ------------------
!!$  if (level.eq.lrefine_max) then
!!$     smooth_count_max = smooth_count_max_finest
!!$  else
!!$     smooth_count_max = smooth_count_max_coarse
!!$  endif

  !!! ---------------------------------------------
  ! Calculate initial residual if wanted for output purposes only


  if (verbose.ge.1.and.level.eq.lrefine_max) then


    do current_var=1,total_vars
        local_max_defect=0.0
        call amr_mg_get_defect(pid, level, local_max_defect, rmsres, npts)
    end do
  endif



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Multigrid bit
  ! Note to self, this isn't as simple as suggested on the website
  ! this restriction (and prolongation) only operates on 'work'
  if(level.gt.mg_min_lvl)then

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Pre-smooth
     do smooth_count=1, smooth_count_max

        call app_mg_smooth(pid, level)



        call MPI_Allreduce(mg_update_total, global_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                           MPI_COMM_WORLD, ierr)
        if (global_total.le.0.0) then
           if (pid.eq.0.and.verbose.ge.4) write(*,*) '(A,I,A)', 'Pre Update get-out after ',smooth_count, ' iterations'
           exit
        endif
     end do


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! CEG - i.e. coarsening

     ! CEG In order to fix the MG coarsening and defect calculations we need to coarsen 
     !     the solutions first of ALL variables, then calculate the residuals on both 
     !     fine and coarse grids without the additional RHS contribution, before adding 
     !     the defect into unk(2) etc
     
     ! * Coarsen all solutions
     gcell_on_cc(:)=.false.
     do current_var=1,total_vars
        current_unk=1+(current_var-1)*nunkvbles
        gcell_on_cc(current_unk)=.true.
        gcell_on_cc(current_unk+1)=.true.
        gcell_on_cc(current_unk+3)=.true.
        gcell_on_cc(current_unk+4)=.true.
     end do

     call pf_restrict(pid, 1, 0, .false., level)
     call reset_nodetype(noprocs, pid, level-1)
     call pf_guardcell(pid,iopt,level-1)





!if(pid==0) write(*,*) 'level', level, 'after restrict','in',pid

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Start looping code block 1
     ! Note to self, within the "looping code" only make guardcell calls on work array
     ! put the full guardcell call at the end


     do current_var=1,total_vars
     
        ! set the start point for the real variable and it's workspace data
        current_unk=1+(current_var-1)*nunkvbles
        current_work=1
        ! Reset to current ref level, this only really needs to be done for current var>1
        ! Defect calculation
        if (pid.eq.0.and.verbose.ge.4.and.current_var.eq.1)write(*,*) 'Residual to coarsening'
        call zero_level_work_array(level)
        local_max_defect=0.0

        call amr_mg_get_defect(pid, level, local_max_defect, rmsres, npts)


!            call pf_restrict(pid, 1, 0, .false., level)
!!$            call reset_nodetype(noprocs, pid, level-1)
!!$            call pf_guardcell(pid,iopt,level-1)



       do lbs=block_starts(level), block_starts(level+1)-1
           if(morton_flag) then
              lb = block_list(lbs)
           else
              lb = lbs
           endif
           !work2= dh
           do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
                 do i=il_bnd+nguard,iu_bnd-nguard
                    work(i,j,k,lb,2)=unk(current_unk+5,i,j,k,lb)
                 enddo
              enddo
           enddo
        enddo


        ! Need guardcells of defect setting up, defect is in work(2) so call guardcell with iopt=3
        call reset_nodetype(noprocs,pid,level)
        call pf_guardcell(pid,current_work+2,level)

        ! CEG Restrict work(2)
        call pf_mg_restrict(noprocs,pid,level)
        call reset_nodetype(noprocs,pid,level-1)


!!$        do lbs=block_starts(level), block_starts(level+1)-1
!!$           if(morton_flag) then
!!$              lb = block_list(lbs)
!!$           else
!!$              lb = lbs
!!$           endif
!!$           work(:,:,:,lb,1)=work(:,:,:,lb,2)
!!$        enddo


        ! Need guardcells of defect setting up, defect is in work(2) so call guardcell with iopt=3
        ! presumable iopt=1 is unk, iopt=2 is work(1) -pcb

 

        ! CEG Restrict work(2)
!        call pf_mg_restrict(noprocs,pid,level)
!!$            gcell_on_cc(:)=.false.
!!$            gcell_on_cc(current_unk+5)=.true.
!!$        call pf_restrict(pid, 1, 0, .false., level)
!!$        call pf_guardcell(pid,1,level)
!!$        call pf_restrict(pid, 3, 0, .false., level)
!!$        call reset_nodetype(noprocs,pid,level)
!!$        call pf_guardcell(pid,3,level)
!!$        call reset_nodetype(noprocs,pid,level-1)

 

       ! CEG Now calculate RHS defect
        if (pid.eq.0.and.verbose.ge.4.and.current_var.eq.1)write(*,*) 'Residual of coarsened solution'
        local_max_defect=0.0
        call amr_mg_get_defect(pid, level-1,local_max_defect, rmsres, npts)

       ! CEG  Finally add FAS RHS to solution

!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb,i,j,k)
!!!$OMP DO SCHEDULE(DYNAMIC,1)
        do lbs=block_starts(level-1), block_starts(level)-1
           if(morton_flag) then
              lb = block_list(lbs)
           else
              lb = lbs
           endif
           do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
                 do i=il_bnd+nguard,iu_bnd-nguard
                    if(has_children(lb)) then
                       ! RHS <-- previous timesteps - defect of coarsened solution
!                       if(work(i,j,k,lb,1).gt.1d-10)  write(*,*)unk(current_unk+1,i,j,k,lb),work(i,j,k,lb,2), work(i,j,k,lb,1)

                       ! new f2h = restricted(dh)+A(restricted(u))
                       !         = restricted(dh)+(f2h -d2h)
                       !         = work(1)  + unk(curr+1) - work(2)
                       ! res(fh)==f2h
                       !f2h = res(fh)- d2h + res(dh)
                       !    = unk(curr+1) - work2 + work1
                       unk(current_unk+1,i,j,k,lb)=unk(current_unk+1,i,j,k,lb) - unk(current_var*nunkvbles,i,j,k,lb) + work(i,j,k,lb,2)
!                       unk(current_unk+1,i,j,k,lb)=unk(current_unk+1,i,j,k,lb) - work(i,j,k,lb,2) + work(i,j,k,lb,1)
                       unk(current_unk+2,i,j,k,lb)=unk(current_unk,i,j,k,lb)   ! Store coarsened value as no longer done via work
                    else
                       unk(current_unk+2,i,j,k,lb)=unk(current_unk,i,j,k,lb)
                    end if
                 end do
              end do
           end do
        end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 
        gcell_on_cc(current_unk:current_unk+4)=.false.  ! Preparatory for setting at end of loop for each variable
        gcell_on_cc(current_unk)=.true.
        gcell_on_cc(current_unk+1)=.true.
!        gcell_on_cc(current_unk+2)=.true.
        gcell_on_cc(current_unk+3)=.true.
        gcell_on_cc(current_unk+4)=.true.

     end do !! enddo of current_var=1,total_vars
     call pf_guardcell(pid,iopt,level-1)
     gcell_on_cc(:)=.true.

     ! End looping code block 1
!if(pid==0) write(*,*) 'level', level, 'after end of code block1','in',pid
             
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Ok Think everything is now assembled let's try that recursive call
     new_level=level-1
     ! CEG - i.e. coarsened
     call reset_nodetype(noprocs, pid, level-1)
        
     
     call amr_multigrid_v_cycle(pid, noprocs, new_level, global_max_defect, rmsres, npts)

     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     ! Start looping code block 2
     ! Store coarse grid correction work(...1)
     ! Now we need to subtract the restricted version of v from the current (unk(3))

     do current_var=1,total_vars,1
        call reset_nodetype(noprocs,pid,level-1)
        ! set the start point for the real variable and its workspace data
        current_unk=1+(current_var-1)*nunkvbles
        current_work=1
!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb,i,j,k)
!!!$OMP DO SCHEDULE(DYNAMIC,1)
!if(pid==0) write(*,*) 'level', level, 'before work= subtr','in',pid
        do lbs=block_starts(level-1), block_starts(level)-1
           if(morton_flag) then
              lb = block_list(lbs)
           else
              lb = lbs
           endif
           ! Put result straight into work ready for prolongation
           do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
                 do i=il_bnd+nguard,iu_bnd-nguard
                    if(has_children(lb))then
                       ! work(...,1) is change in solution over the solve
                       work(i,j,k,lb,current_work) = unk(current_unk,i,j,k,lb) - unk(current_unk+2,i,j,k,lb)
                    else
                       work(i,j,k,lb,current_work) = 0.0
                    endif
                 end do
              end do
           end do
        end do

!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 
!if(pid==0) write(*,*) 'level', level, 'before guardcell 1',' in',pid

        call pf_guardcell(pid,1+current_work,level)
        ! Reset to current ref level, this only really needs to be done for current var>1
!if(pid==0) write(*,*) 'level', level, 'after guardcell 1',' in',pid

        call amr_mg_prolong(noprocs,pid,level)  ! Prolongs work array
!if(pid==0) write(*,*) 'level', level, 'after prolong','in',pid
        call reset_nodetype(noprocs, pid, level)
!if(pid==0) write(*,*) 'level', level, 'reset nodetype1','in',pid
        call pf_guardcell(pid,1+current_work,level+1)
!if(pid==0) write(*,*) 'level', level, 'guardcell 2','in',pid
        call reset_nodetype(noprocs,pid,level)
!if(pid==0) write(*,*) 'level', level, 'after prolong reset nodetype','in',pid

!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb,i,j,k)
!!!$OMP DO SCHEDULE(DYNAMIC,1)
        do lbs=block_starts(level), block_starts(level+1)-1
           if(morton_flag) then
              lb = block_list(lbs)
           else
              lb = lbs
           endif
           do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
                 do i=il_bnd+nguard,iu_bnd-nguard
                    ! Fine grid solution <-- pre-coarsened solution + coarse grid correction
                    unk(current_unk,i,j,k,lb) = unk(current_unk,i,j,k,lb) + work(i,j,k,lb,current_work)
                 end do
              end do
           end do
        end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 
        gcell_on_cc(:)=.false.
        gcell_on_cc(current_unk)=.true.
        gcell_on_cc(current_unk+3)=.true.
        gcell_on_cc(current_unk+4)=.true.
        call pf_guardcell(pid,iopt,level)
        gcell_on_cc(:)=.true.
     end do
     ! End looping code block 2
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!if(pid==0) write(*,*) 'level', level, 'after end of code block2','in',pid

     ! Post Smooth
!!$     if(level.ge.lrefine_max-1.and.ndim.eq.3)then
!!$        smooth_count_max_=4*smooth_count_max
!!$     else
!!$        smooth_count_max_=smooth_count_max
!!$     endif
     do smooth_count=1,smooth_count_max
        call app_mg_smooth(pid,level)
        call MPI_Allreduce(mg_update_total, global_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                           MPI_COMM_WORLD, ierr)
        if (global_total.le.0.0) then
           if (pid.eq.0.and.verbose.ge.4) write(*, *) 'Post Update get out after ',smooth_count, ' iterations'
           exit
        endif
     end do
!if(pid==0) write(*,*) 'level', level, 'after end of postsmoothing','in',pid

  else
     ! CEG Coarse grid solve
     global_max_defect = 1.0
     do smooth_count=1,solve_count_max
        if(global_max_defect.lt.defect_tol_min*defect_tol_min)then
           exit
        end if
           
        call app_mg_smooth(pid,level)
        call MPI_Allreduce(mg_update_total, global_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                           MPI_COMM_WORLD, ierr)
        if (global_total.le.0.0) then
           if (pid.eq.0.and.verbose.ge.4) write(*, *) 'CG Update get out after ',smooth_count, ' iterations'
           exit
        endif
     end do
!if(pid==0) write(*,*) 'level', level, 'coarsest solve: after solve ','in',pid

     local_max_defect=0.0
     do current_var=1,total_vars
        ! set the start point for the real variable and its workspace data
        current_unk=1+(current_var-1)*nunkvbles
        current_work=1
        if (pid.eq.0.and.verbose.ge.4.and.current_var.eq.1)write(*,*) 'Coarse grid solve defect'
        call amr_mg_get_defect(pid, level,local_max_defect, rmsres, npts)
!if(pid==0) write(*,*) 'level', level, 'coarsest solve: in loop getdefect ', current_var

     end do
     call MPI_AllReduce(local_max_defect,global_max_defect,&
                        1,amr_mpi_real,MPI_MAX,MPI_COMM_WORLD,ierr)
!if(pid==0) write(*,*) 'level', level, 'coarsest solve: after getdefect ','in',pid

  end if

  ! Defect calculation
  local_max_defect=0.0
  do current_var=1,total_vars
     current_unk=1+(current_var-1)*nunkvbles
     current_work=1
     if (pid.eq.0.and.verbose.ge.4.and.current_var.eq.1)write(*,*) 'Post smooths defect'
     call amr_mg_get_defect(pid, level,local_max_defect, rmsres, npts)
  end do


  ! NOTE: max of defects only needed at the end of v-cycle
  call MPI_AllReduce(local_max_defect,global_max_defect,&
       1,amr_mpi_real,MPI_MAX,MPI_COMM_WORLD,ierr)

end subroutine amr_multigrid_v_cycle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_mg_get_defect(pid, level, local_max_defect, rms, npts)
  use multigrid_parameters
  use generic_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, level
  double precision, intent(inout) :: local_max_defect
  integer :: lb, npts, Gnpts, ierr, i
  double precision, intent(inout) ::  rms(total_vars)
  double precision :: global_max_defect, outres, Grms(total_vars), locrms
  integer :: lbs
  rms(current_var) = 0.0
  npts = 0
  if(current_var.gt.1)npts=1
  locrms = 0.0
!$OMP PARALLEL SHARED(level, local_max_defect, npts, locrms) PRIVATE(lb)
! REDUCTION(MAX:locrms)
!  if(total_vars.le.3.or.current_var.eq.1)then
  if (multigrid_on.ge.1) then


!$OMP DO SCHEDULE(DYNAMIC,1) REDUCTION(+:locrms) REDUCTION(+:npts) REDUCTION(MAX:local_max_defect)
     do lbs = block_starts(level), block_starts(level+1)-1
           if(morton_flag) then
             lb = block_list(lbs)
           else
              lb = lbs
           endif

        call app_mg_get_defect_var(lb,local_max_defect, locrms, npts)


     end do
!$OMP END DO NOWAIT
  else 
!$OMP DO SCHEDULE(DYNAMIC,1) REDUCTION(+:locrms) REDUCTION(+:npts) REDUCTION(MAX:local_max_defect)
     do lb = 1, lnblocks
        if (.not.has_children(lb)) call app_mg_get_defect_var(lb,local_max_defect, locrms, npts)
     end do
!$OMP END DO NOWAIT
  endif
!  endif
  rms(current_var) = locrms
!  if(npts.eq.0)then
!     write(*,*)npts,current_var,locrms
!     stop
!  endif
!$OMP END PARALLEL 
!  print *, locrms, rms(current_var), npts
  ! Output if verbose >= 3
  if ((verbose.eq.2.and.current_var.eq.1).or.verbose.ge.4) then
     call MPI_Reduce(local_max_defect,global_max_defect,              &
                        1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
     call MPI_Reduce(rms, Grms, total_vars, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     call MPI_Reduce(npts, Gnpts, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     outres = 0.0
     if (pid.eq.0) then
        do i = 1, total_vars
           Grms(i) = Grms(i)/real(Gnpts)
           outres = outres + Grms(i)
           Grms(i) = sqrt(Grms(i))
        end do
     endif

     if (pid.eq.0) write(*, '(A,I2,X,A,X,I2,X,A,E10.3,4X,A,E10.3)') 'Grid ', &
                level, 'Vble', current_var, 'Max defect: ', global_max_defect, 'RMS:', Grms(current_var)
  end if

  return

end subroutine amr_mg_get_defect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine zero_level_work_array(level)
  use multigrid_parameters
  use generic_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  include 'mpif.h'
  integer, intent(in) :: level
  integer :: lb,lbs

!!!$OMP PARALLEL SHARED(level) PRIVATE(lb)
!!!$OMP DO SCHEDULE(DYNAMIC,1)
  do lbs=block_starts(level), block_starts(level+1)-1
           if(morton_flag) then
              lb = block_list(lbs)
           else
              lb = lbs
           endif
     work(:,:,:,lb,:) = 0.0
  end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 

  return

end subroutine zero_level_work_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_multigrid_child_set()
  ! Sweeps through all (local) blocks and checks if they're parents or not
  ! This needs to be called once at start and once after all adaptivity events
  use paramesh_dimensions
  use physicaldata
  use tree
  use multigrid_parameters
  implicit none
  include 'mpif.h'
  integer lb
  ! Only operating on blocks not on finest possible as these should never have chldren set true
!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb)
!!!$OMP DO SCHEDULE(DYNAMIC,1)
  do lb=1,lnblocks
!  do lb=1,block_starts(lrefine_max)-1
     if(child(1,1,lb).eq.-1)then
        has_children(lb)=.false.
     else if(child(1,1,lb).gt.0)then
        has_children(lb)=.true.
     else
        write(*,*)"Odd value in child array for block ", lb, "Value: ", child(1,lb,1)
     end if
  end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 
end subroutine amr_multigrid_child_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine reset_nodetype(noprocs,pid,level)
  use paramesh_dimensions
  use physicaldata
  use tree
  use multigrid_parameters
  implicit none
  integer, intent(in) :: noprocs,pid,level
  integer :: lb
  ! This resets nodetypes so that mlat can work alongside guardcell call
  ! Basically it's a fudge and only really sorts which blocks are nodetype 1
  ! It does have one obvious bonus, that is nodetype becomes redundant outside of GC
  ! calls and this can reset it if needed
  ! Conditions for nodetypes
  ! if a block is more refined than the current "level" do nothing
  ! if a block is same as current level set to 1
  ! if a block is less refined than current level AND has no children, set to 1
  ! if a block is less refined than current level AND has children, leave alone
  ! <snip>
  ! James: Guess what, I have a NEW way of doing this now
!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb)
!!!$OMP DO SCHEDULE(DYNAMIC,1)
  do lb=1,lnblocks
     nodetype(lb)=nodetype_copy(lb)
     if(lrefine(lb).gt.level)then
        nodetype(lb) = -1
     else if(lrefine(lb).eq.level)then
        nodetype(lb) = 1
     else if((lrefine(lb).eq.level-1).and.(nodetype(lb).ne.1))then
        nodetype(lb) = 2
     end if
  end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 

  call amr_get_new_nodetypes (noprocs, pid, level)

end subroutine reset_nodetype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_mg_get_max_refinement(pid,noprocs)
  ! This subroutine sweeps through all blocks and finds the max ref level on the current processor
  ! Then it passes the max to all procs
  ! Why is this needed? Well multigrid needs to know what level to start at
  ! Should be called after all refinement events
  use tree
  use multigrid_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid,noprocs
  integer :: lb,local_max
  local_max=0
!!!$OMP PARALLEL SHARED(lnblocks,local_max) PRIVATE(lb)
!!!$OMP DO SCHEDULE(DYNAMIC,1) REDUCTION (MAX:local_max)
  do lb=1,lnblocks
     if(lrefine(lb).gt.local_max)then
        local_max=lrefine(lb)
     end if
  end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 

  call MPI_AllReduce(local_max,mg_max_level,&
       1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,lb)

end subroutine amr_mg_get_max_refinement

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_non_linear_solve(pid, noprocs, level, global_max_defect, rmsres, npts)
!
  use paramesh_interfaces
  use multigrid_parameters
  use generic_parameters
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use tree
  use workspace
  use paramesh_comm_data ! Needed for definition of amr_mpi_real
  implicit none  
  include 'mpif.h'
  integer, intent(in) :: pid,noprocs, level
  integer, intent(inout) :: npts
  double precision, intent(inout) :: global_max_defect, rmsres(total_vars)
  integer :: smooth_count,lb,nlayers,iopt,ierr
  double precision :: local_max_defect, global_total
  integer :: i,j,k,new_level,lvl

  ! Plan:
  ! Need to cycle over all the leaf (node) blocks smoothing
  ! Next job is to find the defect
  ! Now do the FAS bit...
  ! Restrict if leaf nodes are not the coarsest level
  ! Possibly have to swap lots of data about here, 
  !    probably put this in it's own sub for my sanity
  ! Call solver again (recursive)
  ! More swapping? If so own sub for sanity
  ! Refine back
  ! Reverse FAS bit
  ! Smooth again
  ! Find defect/residual
  ! Setup a couple of parameters

  iopt = 1

  nlayers = nguard
  if(explicit)then
     gcell_on_cc(:) = .false.


     do current_var = 1,total_vars
        current_unk = 1+(current_var-1)*nunkvbles
        gcell_on_cc(current_unk) = .true.
     end do
     
     call amr_guardcell(pid,iopt,nguard)
     call app_mg_smooth(pid, -1)

     return
  endif
    
  ! Defect calculation
  local_max_defect = 0.0
  do lvl=1,level
  do current_var = 1, total_vars
     current_unk = 1 + (current_var-1)*nunkvbles
     current_work = 1
     if (pid.eq.0.and.verbose.ge.2.and.current_var.eq.1) write(*,*) 'Post smooths defect'
     if(.not.has_children(lb))call amr_mg_get_defect(pid, lvl, local_max_defect, rmsres, npts)
  end do
  enddo

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! smooth
  
  
  ! * Do guardcells on all solutions
  
  gcell_on_cc(:) = .true.
  do current_var = 1,total_vars
     current_unk = 1+(current_var-1)*nunkvbles
     gcell_on_cc(current_unk) = .true.
     gcell_on_cc(current_unk+3) = .true.
     gcell_on_cc(current_unk+4) = .true.
  end do



  
  do smooth_count = 1, smooth_count_max

!     call amr_guardcell(pid,iopt,nguard)
     call app_mg_smooth(pid, -1)
     call amr_guardcell(pid,iopt,nguard)

  end do
  
  ! Defect calculation
  local_max_defect = 0.0
  do lvl=1,level
  do current_var = 1, total_vars
     current_unk = 1 + (current_var-1)*nunkvbles
     current_work = 1
     if (pid.eq.0.and.verbose.ge.2.and.current_var.eq.1) write(*,*) 'Post smooths defect'
     if(.not.has_children(lb))call amr_mg_get_defect(pid, lvl, local_max_defect, rmsres, npts)
  end do
  enddo
  ! NOTE: max of defects only needed at the end of v-cycle
  call MPI_AllReduce(local_max_defect, global_max_defect, 1, amr_mpi_real, MPI_MAX, MPI_COMM_WORLD, ierr)

       
end subroutine amr_non_linear_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

