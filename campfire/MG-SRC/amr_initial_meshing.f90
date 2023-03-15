!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!! amr_initial_meshing
!!!! * Performs global and local refinement on initial solution
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Chris Goodyer, May 2012                                              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_initial_meshing(pid, noprocs)
  use paramesh_interfaces
  use paramesh_dimensions
  use tree
  use time_dep_parameters
  use multigrid_parameters
  use generic_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, noprocs
  integer ::  numlocalref, tmpI(2), loop_count, iopt, nlayers, ierr

  if (local_adaptation.eq.0) numglobalref = lrefine_max - 1
!     numglobalref = lrefine_max - 1 !debug for eutectic
  if (numglobalref.gt.lrefine_max-1) numglobalref = lrefine_max-1

  if (lrefine_max-1.lt.numglobalref) then
     numglobalref = lrefine_max - 1
     numlocalref = 1
  else
     numlocalref = lrefine_max - numglobalref - 1
     if (numlocalref.eq.0) numlocalref = 1
  endif

  iopt = 1
  nlayers = nguard
  
  ! Initialise everything
  refine(1:lnblocks) = .false.
  call amr_refine_derefine

  ! Global refinement
!  print *,'do refinement'
  if (numglobalref.gt.0) then  
     do loop_count=1, numglobalref
        if(pid.eq.0.and.verbose.ge.2) print *,'Pre refinement', loop_count, lnblocks
        call amr_initial_soln
        refine(1:lnblocks) = .true.
     
        ! refine grid
        call amr_refine_derefine
        call amr_prolong(pid,iopt,nlayers)
        call pf_guardcell(pid,iopt,numglobalref+loop_count+1)
     enddo
  end if

  ! Local refinement
  if (numlocalref.gt.0) then
     do loop_count=1,numlocalref
        if(pid.eq.0.and.verbose.ge.2) print *,'Initial refinement', lrefine_min, lrefine_max, loop_count
        call amr_initial_soln
        call amr_test_refinement(pid, lrefine_min, lrefine_max, 1)
        call amr_refine_derefine
        call amr_prolong(pid,iopt,nlayers)
        call pf_guardcell(pid,iopt,numglobalref+loop_count+1)
     enddo
  else
     call amr_initial_soln
     call amr_prolong(pid,iopt,nlayers)
     call pf_guardcell(pid,iopt,numglobalref+loop_count+1)
  endif
    
end subroutine amr_initial_meshing


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!! amr_relocate_meshing
!!!! * essentially required for C2F feature
!!!! * Performs global and local refinement on a check point input
!!!! * this subroutine only refines blocks on desired finer refinement level,
!!!!   without changing values of current data.
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Tim Yang, Jan 2014. Modified from Chris Goodyer's amr_initial_meshing!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_relocate_meshing(pid, noprocs)
  use paramesh_interfaces
  use paramesh_dimensions
  use tree
  use time_dep_parameters
  use multigrid_parameters
  use generic_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, noprocs
  integer ::  numlocalref, tmpI(2), loop_count, iopt, nlayers, ierr, lb

write(*,*)"in Tims"
stop


  if (local_adaptation.eq.0) numglobalref = lrefine_max - 1
!     numglobalref = lrefine_max - 1 !debug for eutectic
  if (numglobalref.gt.lrefine_max-1) numglobalref = lrefine_max-1

  if (lrefine_max-1.lt.numglobalref) then
     numglobalref = lrefine_max - 1
     numlocalref = 1
  else
     numlocalref = lrefine_max - numglobalref - 1
     if (numlocalref.eq.0) numlocalref = 1
  endif

  iopt = 1
  nlayers = nguard
  
  ! Initialise everything
  refine(1:lnblocks) = .false.
  call amr_refine_derefine

  ! Global refinement
!  print *,'do refinement'
!!  if (numglobalref.gt.0) then  
!!     do loop_count=1, numglobalref
!!        if(pid.eq.0.and.verbose.ge.2) print *,'Pre refinement', loop_count, lnblocks
!!--tim        call amr_initial_soln
!!        refine(1:lnblocks) = .true.

!1        if(pid.eq.0)then
          do lb=1,lnblocks
            if(lrefine(lb).eq.C2F_start_level)then
              refine(lb)=.true.
 !!             write(*,*) lnblocks, parent(1,lb), lrefine(lnblocks), "-----------------"
              exit
            endif
          enddo
 !1       else
!1          refine(1:lnblocks)=.false.
!1        endif
!1        write(*,*) "----", pid, maxblocks_tr

        
     
        ! refine grid
        call amr_refine_derefine
        call amr_prolong(pid,iopt,nlayers)
        call pf_guardcell(pid,iopt,C2F_start_level+1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!! two jumps
        if(C2F_desired_level-C2F_start_level.gt.1)then
          do lb=1,lnblocks
            if(lrefine(lb).eq.C2F_start_level+1)then
              refine(lb)=.true.
 !!             write(*,*) lnblocks, parent(1,lb), lrefine(lnblocks), "-----------------"
              exit
            endif
          enddo
        call amr_refine_derefine
        call amr_prolong(pid,iopt,nlayers)
        call pf_guardcell(pid,iopt,C2F_start_level+2)         

        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!     enddo
!!  end if

!        call pf_mg_init()
        call amr_mg_init()
!! write(*,*) "should not get called"
!!        do lb=1,lnblocks
!!          if(lrefine(lb).eq.C2F_desired_level)then
!!            refine(lb)=.false.
!!            derefine(lb)=.true.
!!          endif
!!        enddo
!!        call amr_refine_derefine
!!        call amr_prolong(pid,iopt,nlayers)
        call pf_guardcell(pid,iopt,C2F_desired_level)

  ! Local refinement
  if (numlocalref.gt.0) then
!1     do loop_count=1,numlocalref
        if(pid.eq.0.and.verbose.ge.2) print *,'Initial refinement', lrefine_min, lrefine_max, loop_count
!!--tim        call amr_initial_soln
!!--        call amr_test_refinement(pid, lrefine_min, lrefine_max, 1)
!!--        call amr_refine_derefine
!!--        call amr_prolong(pid,iopt,nlayers)
!!--        call pf_guardcell(pid,iopt,numglobalref+loop_count+1)
!!     enddo
  else
!!--tim     call amr_initial_soln
!!--     call amr_prolong(pid,iopt,nlayers)
!!--     call pf_guardcell(pid,iopt,numglobalref+loop_count+1)
  endif
    
end subroutine amr_relocate_meshing
