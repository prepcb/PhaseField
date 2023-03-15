!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_vcycle_summary
!!!!  * Prints info on convergence of individual V-cycles
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_multigrid_block_types
!!!!  * Establishes how many blocks on each level and prints out
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Chris Goodyer, September 2012                                        !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_vcycle_summary(pid, noprocs, i_step, defect, old_defect, npts, rms)
  use multigrid_parameters
  use paramesh_dimensions
  use generic_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, noprocs, i_step, npts
  double precision, intent(inout) :: defect, old_defect, rms(total_vars)
  double precision :: convrate, Grms(total_vars), outres
  integer :: Gnpts, ierr, i

  call MPI_Reduce(rms, Grms, total_vars, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(npts, Gnpts, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  outres = 0.0
  if (pid.eq.0) then
     do i = 1, total_vars
        Grms(i) = Grms(i)/real(Gnpts)
        outres = outres + Grms(i)
        rms(i) = sqrt(Grms(i))
     end do
!     print *, 'RMS = ', Grms(1:total_vars),rms(1:total_vars)
  endif
  outres = sqrt(outres)

  if(pid.eq.0.and.verbose.ne.0)then
     if(i_step.eq.1) then
        if(pid.eq.0)write(*,2001) i_step, defect, outres
     else
        convrate = old_defect/defect
        if(pid.eq.0)write(*,2000) i_step, defect, outres, convrate
     end if
   end if

 return

 2001   format ("Vcycle: ",I2," Max Defect=",E14.7," RMS res", E14.7) 
 2000   format ("Vcycle: ",I2," Max Defect=",E14.7," RMS res", E14.7," rate=",F10.3) 
end subroutine amr_vcycle_summary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_multigrid_block_types(pid,noprocs,finegridblocks)
!
!    Checks all blocks and sees what ref levels/node types they all are  
!
  use paramesh_dimensions
  use physicaldata
  use tree
  use generic_parameters
  use multigrid_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid,noprocs
  integer, intent(inout) :: finegridblocks(2)
  integer :: blocks_type(-1:lrefine_max),blocks_ref(-1:lrefine_max)
  integer :: temp(-1:lrefine_max),recv_buff(1:(lrefine_max+2)*(noprocs))
  integer :: lb,ierr,proc, i
  do lb=-1,lrefine_max
     temp(lb)=0
     blocks_ref(lb)=0
  end do
  do lb=1,lnblocks
     ! Put all data into temp on each proc, then gather it on proc 1 into the final data structs
     if(.not.has_children(lb))temp(lrefine(lb))=temp(lrefine(lb))+1
  end do
  call MPI_Gather(temp(-1),lrefine_max+2,MPI_INTEGER, &
       recv_buff(1),lrefine_max+2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if(pid.eq.0)then
     if (verbose.ge.2) call coutput(lrefine_max, noprocs, recv_buff)
!        print*,'F90 version', noprocs
!        do i=1, lrefine_max+2
!           print *, i, recv_buff(i)
!        end do

     do lb=-1,lrefine_max
        do proc=0,noprocs-1
         blocks_ref(lb)=blocks_ref(lb)+recv_buff(lb+2+proc*(lrefine_max+2))
        end do
     end do    
  end if
  ! Done ref level, now do nodetype
  ! Assuming that nodetype has the same range as lrefine
  do lb=-1,lrefine_max
     temp(lb)=0
     blocks_type(lb)=0
  end do
  do lb=1,lnblocks
     temp(nodetype(lb))=temp(nodetype(lb))+1
  end do
  finegridblocks(2) = 0
  call MPI_Gather(temp(-1),lrefine_max+2,MPI_INTEGER, &
       recv_buff(1),lrefine_max+2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if(pid.eq.0)then
     do lb=-1,lrefine_max
        do proc=0,noprocs-1
           blocks_type(lb)=blocks_type(lb)+recv_buff(lb+2+proc*(lrefine_max+2))
        end do
     end do
     if (verbose.ge.1) print *,"Ref/Type No. | Blocks by type | Blocks by ref"
     do lb=-1,lrefine_max
        finegridblocks(2) = finegridblocks(2) + blocks_ref(lb)
        if (verbose.ge.1) print *,lb,blocks_type(lb),blocks_ref(lb)
     end do
     finegridblocks(1) = blocks_ref(lrefine_max)
  end if

end subroutine amr_multigrid_block_types
