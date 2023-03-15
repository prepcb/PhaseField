!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!  
!!!!  amr_test_refinement
!!!!  * Default mesh refinement testing routine
!!!!  * Looks at gradients in cell to decide if needs to be refined
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Chris Goodyer, May 2012                                              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

! PCB 9/3/16
! Modification at end using mg_min_lvl allows control over the coarsest uniform grid. It seems that
! paramesh does not allow level 1 to be the coarsest uniform grid and ofetn not 2. This means that
!  mg_min_lvl=3 and mg_min_lvl=1 have the same refinement but differ in the multigrid level. But
!  mg_min_lvl=4 has a different grid (more elements)

subroutine amr_test_refinement(pid, llrefine_min, llrefine_max, initialadapt)

  use paramesh_dimensions
  use physicaldata
  use tree
  use multigrid_parameters
  use refinement_parameters
  implicit none
  integer, intent(in) :: pid, llrefine_min, llrefine_max
  integer, optional :: initialadapt
  double precision :: dx,error,e1,e2,e3,error_max(total_vars),error_mod(total_vars), dist(8), dist2(4), errnum, errdenom
  integer :: ndel,nlayers,i,j,k,lb,iopt, inadapter, u, v
  logical :: allsolid,unkgt=.false.
!  uniform_mesh=.false.

  ndel = (nguard_work - nguard)*npgs
  ! Re-initialize the refinement and derefinement flag arrays
  if(uniform_mesh)then
     refine(:)   = .true.
     derefine(:) = .false.
     return
  endif

 




  refine(:)   = .false.
  derefine(:) = .false.

  error_max(:) = 0.0
  error_mod(:)=1d0

!  error_mod(1)=1d0
  iopt=1
  nlayers=1
  inadapter=0
  if(present(initialadapt)) inadapter = initialadapt
  
  
!!$  if (allsolid_test.eq.1) then
!!$     allsolid = .true.
!!$  else
!!$     allsolid = .false.
!!$  endif
  

  if(lnblocks.gt.0) then
     do lb=1,lnblocks
        unkgt=.false.
        if (inadapter.eq.1) then
           if (ndim.eq.2) then
              if(total_vars.le.3)then
                 if(uniformstart.or.sqrt(bnd_box(1,1,lb)**2+bnd_box(1,2,lb)**2).lt.2d1)refine(lb) = .true.
              else !multiphase eutectic
!                 if(abs(bnd_box(1,1,lb))+abs(bnd_box(1,2,lb)).lt.1d1)refine(lb) = .true.
                 if(abs(bnd_box(1,1,lb)).lt.1d1)refine(lb) = .true.
!                 write(*,*)lb,bnd_box(1,1,lb),bnd_box(1,2,lb)
!                 refine(lb) = .true. 
              endif
           else           
              dist(1) = abs(bnd_box(1,1,lb))+abs(bnd_box(1,2,lb))+abs(bnd_box(1,3,lb))
              dist(2) = abs(bnd_box(1,1,lb))+abs(bnd_box(1,2,lb))+abs(bnd_box(2,3,lb))
              dist(3) = abs(bnd_box(1,1,lb))+abs(bnd_box(2,2,lb))+abs(bnd_box(1,3,lb))
              dist(4) = abs(bnd_box(1,1,lb))+abs(bnd_box(2,2,lb))+abs(bnd_box(2,3,lb))
              dist(5) = abs(bnd_box(2,1,lb))+abs(bnd_box(1,2,lb))+abs(bnd_box(1,3,lb))
              dist(6) = abs(bnd_box(2,1,lb))+abs(bnd_box(1,2,lb))+abs(bnd_box(2,3,lb))
              dist(7) = abs(bnd_box(2,1,lb))+abs(bnd_box(2,2,lb))+abs(bnd_box(1,3,lb))
              dist(8) = abs(bnd_box(2,1,lb))+abs(bnd_box(2,2,lb))+abs(bnd_box(2,3,lb))
              
              if (MINVAL(dist).lt.1.0) refine(lb) = .true.
              
           endif
           
           
        else if(nodetype(lb).eq.1) then
           if(g_gpure.and. .false.)then
              write(*,*)"in test refinement"
              stop
              do k=kl_bnd+nguard*k3d, ku_bnd-nguard*k3d 
                 do j=jl_bnd+nguard*k2d, ju_bnd-nguard*k2d 
                    do i=il_bnd+nguard, iu_bnd-nguard
                       if(1d0-unk(1,i,j,k,lb).gt.1d-4 .and. unk(1,i,j,k,lb).gt.1d-4)then
                          refine(lb) = .true.
                          derefine(lb) = .false.
                       else
                          refine(lb) = .false.
                          derefine(lb) = .true.
                          
                       endif
                    enddo
                 enddo
              enddo
           else
              
              
              
              if(total_vars.le.5)then
                 do v = 2,2!!!total_vars
                    u = 1 + (v-1)*nunkvbles
                    error_max(v) = 0.0
                    do k=kl_bnd+nguard*k3d, ku_bnd-nguard*k3d 
                       do j=jl_bnd+nguard*k2d, ju_bnd-nguard*k2d 
                          do i=il_bnd+nguard, iu_bnd-nguard
                             
                             
                             
                             
!!$                             e1 = (unk(u,i+1,j,k,lb)-unk(u,i,j,k,lb))**2
!!$                             e1 = e1 +(unk(u,i-1,j,k,lb)-unk(u,i,j,k,lb))**2
!!$                             e1 = e1 +(unk(u,i,j+1,k,lb)-unk(u,i,j,k,lb))**2
!!$                             e1 = e1 +(unk(u,i,j-1,k,lb)-unk(u,i,j,k,lb))**2

                             e1 = unk(u,i+1,j,k,lb)-4*unk(u,i,j,k,lb)+unk(u,i-1,j,k,lb)
                             e1 = e1 +unk(u,i,j-1,k,lb) +unk(u,i,j+1,k,lb)
                             if(ndim.eq.2)e1 = abs(e1)
                             


!                             if(g_gpure) e1=e1*(1d0-unk(1,i,j,k,lb))**4
                             if (ndim.eq.3) then
!                                e1 = e1 +(unk(u,i,j,k+1,lb)-unk(u,i,j,k,lb))**2
!                                e1 = e1 +(unk(u,i,j,k-1,lb)-unk(u,i,j,k,lb))**2
                            
                                e1 = e1 + unk(u,i,j,k+1,lb)+unk(u,i,j,k-1,lb) -2d0* unk(u,i,j,k,lb)
!!$                                e1 = 0d0
!!$                                if(0.5-abs(unk(1,i,j,k,lb)-0.5).gt.0.1)e1 =1d0 

                                error = 1d3*abs(e1)
                                

                             else
                                error = 1d2*e1*error_mod(v)
!!$                                if(1d0-unk(1,i,j,k,lb).gt.2d-1 .and. unk(1,i,j,k,lb).gt.1d-4.and. g_gpure)then
!!$                                   error = error*1d2
!!$                                endif
                             endif
                             
                             if (error.gt.error_max(v)) error_max(v)=error*error_mod(v)
                             if(total_vars.eq.3)then
                                if (unk(1,i,j,k,lb).lt.0.999) allsolid = .false.
                             else
                                if (unk(1,i,j,k,lb).gt.0.001) allsolid = .false.
                             endif
                          enddo
                       enddo
                    enddo
                    
                    
                 end do
                 
                 
                 if (allsolid) then
!                    print *, 'Derefining block ', lb
                    refine(lb) = .false.
                    derefine(lb) = .true.
!!$              elseif(abs(1d0-unk(1,i,j,k,lb)).gt.1d-12)then
!!$                 refine(lb)=.true.
!!$                 derefine(lb) = .false.
!!$              elseif(unkgt)then
!!$                 refine(lb) = .true.
!!$                 derefine(lb) = .false.
                    
                 else
                    errdenom = total_vars
                    errnum = 0.0
                    
                    do v = 1, total_vars
                       errnum = errnum + error_max(v)
                       
                    enddo
                    
                    error=errnum/errdenom
                    if(nodetype(lb).eq.1)then
                       if (( error .ge. ctore.and.lrefine(lb).lt.llrefine_max).or.lrefine(lb).lt.mg_min_lvl)then
                          refine(lb) = .true.
                       else
                          if ( error .lt. 0.5*ctode.and. lrefine(lb).gt.mg_min_lvl) then
                             derefine(lb) = .true.
                          endif
                       endif
                    endif
                    

                 endif
              endif
           endif
        endif
     end do
  endif
end subroutine amr_test_refinement

