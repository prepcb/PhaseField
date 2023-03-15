!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!! Routines to balance the computational load per multigrid grid.
!!!! Replaces Morton ordering within Paramesh.
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!! ceg_block_balancing_link
!!!! * Stub calling ceg_block_balancing
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!! ceg_block_balancing
!!!! * Establishes how many blocks are needed on each processor 
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!! ceg_setup_newloc
!!!! * Fills the new_loc array which tells Paramesh which blocks to move where
!!!! * Also sorts out the ordering of blocks per processor so levels are
!!!!   contiguous in memory per processor
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Chris Goodyer, May 2012                                              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This file contains routines for balancing the distribution of block for MLAT versions of PARAMESH

!subroutine ceg_block_balancing(pid, noprocs)
subroutine ceg_block_balancing_link(pid, noprocs)

  use paramesh_dimensions
  use physicaldata
  use tree
  use paramesh_interfaces
  integer, intent(in) :: pid, noprocs
  integer, dimension(2, maxblocks_tr) :: new_loc

interface
  subroutine ceg_block_balancing(pid, noprocs, new_loc)
  use paramesh_dimensions
  use physicaldata
  use tree
  use paramesh_interfaces
    integer, intent(in) :: pid, noprocs
!    integer, intent(inout) :: new_loc(:,:)
    integer, dimension(2, maxblocks_tr), intent(inout) :: new_loc
    end subroutine ceg_block_balancing
end interface

  call ceg_block_balancing(pid, noprocs, new_loc)

!  print*,'out of ceg_block_balancing_link'

end subroutine ceg_block_balancing_link

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ceg_block_balancing(pid, noprocs, new_loc)

  use paramesh_dimensions
  use physicaldata
  use tree
  use generic_parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Use paramesh_interfaces, only : amr_compute_morton,              & 
                                  amr_sort_morton             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, noprocs
!  logical, intent(inout) :: mor_flag  !!! added by MCH
  integer :: g, b, bb,bbp,bbl,p, pp, proc, procp, T, ierr, mvd, wanted, total
  integer, dimension(2, maxblocks_tr), intent(inout) :: new_loc
  integer :: i, j, k, bl,blk,i_tot,j_tot,k_tot,k1,k2
!!!!!!!!!!!!! for morton ordering !!!!!!!!!!!!!!!
  Integer,allocatable :: mort_no(:,:)
  Integer,save :: num_called1
!!!!!!!!!!!!! for levelmorton ordering !!!!!!!!!!
  Integer :: displs(noprocs+1), lnblocks_each_proc(noprocs)
  Integer :: tgt(noprocs)
  Integer,allocatable :: mort_no_row6(:), mort_no_row6_all(:)
  Integer,allocatable :: lrefine_all(:), pid_from(:),nodetype_all(:)
  Integer,allocatable :: pid_from_level(:), pid_to_level(:)
  Integer,allocatable :: GridGlobalID(:)
  Integer,allocatable :: temp_list_blocks(:)
  double precision,allocatable ::  bsize_all(:)
  Integer :: blocks_proc_level(lrefine_max*noprocs)
  Integer :: proc_divider_level(noprocs+1)
  Integer :: counter_proc(noprocs),one(10000),two(10000),both(10000),NT_two(10000),NT_both(10000),pid_two(10000),pid_one(10000),pid_both(10000),two_tot(10000),both_tot(10000)
 
  Integer :: total_no_blocks, start, ending
  Integer :: level_total, totalnum_level
  integer :: aaatmp
  integer,allocatable :: aaapid(:),aaabl(:),aaaNT(:),aaaGID(:)
  integer :: list_nodetype(10000),parent_i,leaf_i,one_i(10000),two_i(10000),leaf_i_tot,both_i
!!!!!!!!!!!!! end of setup enviroment !!!!!!!!!!
  num_called1 = num_called1 + 1

  !!!!!!!!!!!
  ! Stage 1 : Discover how many blocks each processor has on each level
  !!!!!!!!!!!
!write(*,*) 'alternative morton',num_called1
  new_loc(:,:)=-1 
!  ltemp = .False.
!  If (present(reorder_grid)) ltemp = reorder_grid

  allocate(mort_no(6,2*maxblocks_tr))

  mort_no = 0

  Call amr_compute_morton (mort_no)






  allocate(mort_no_row6(lnblocks))
  mort_no_row6 = mort_no(6,1:lnblocks)

  call MPI_AllGather(lnblocks, 1, MPI_INTEGER, &
       lnblocks_each_proc, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
  total_no_blocks = sum(lnblocks_each_proc)


  displs(1) = 0
  do i = 2,noprocs+1
     displs(i) = displs(i-1)+ lnblocks_each_proc(i-1)
  enddo 


  allocate(mort_no_row6_all(total_no_blocks), &
           lrefine_all(total_no_blocks), &
           nodetype_all(total_no_blocks), &
           GridGlobalID(total_no_blocks), &
!           bsize_all(total_no_blocks), &
           temp_list_blocks(total_no_blocks))
  allocate(pid_from(total_no_blocks), pid_from_level(total_no_blocks), &
       pid_to_level(total_no_blocks))


  call MPI_AllGatherv(mort_no_row6, lnblocks, MPI_INTEGER, &
       mort_no_row6_all ,lnblocks_each_proc, & 
       displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)

  call MPI_AllGatherv(lrefine, lnblocks, MPI_INTEGER, &
       lrefine_all ,lnblocks_each_proc , &
       displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)

  do b = 1, total_no_blocks
     GridGlobalID(b) = b
  enddo

  pid_from = 0 
  do proc = 0, noprocs-1
     start = displs(proc+1)+1
     ending = displs(proc+2)
     do b = start, ending
        pid_from(b)=proc
     enddo
  enddo
  call MPI_AllGatherv(nodetype, lnblocks, MPI_INTEGER, &
       nodetype_all ,lnblocks_each_proc , &
       displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)


  call sort(mort_no_row6_all, lrefine_all, GridglobalID, nodetype_all,pid_from, &
       total_no_blocks, num_called1, pid)



  if(PCB_order)then
     do g = lrefine_max-1,1,-1
        i=0
        j=0
        do bl = 1, total_no_blocks
           if(lrefine_all(bl).eq.g)then
              if(nodetype_all(bl).ne.1)then
                 i = i + 1
                 two(i) = gridglobalID(bl)
                 pid_two(i) = pid_from(bl)
              elseif(nodetype_all(bl).eq.1)then
                 j = j + 1
                 one(j) = gridglobalID(bl)
                 pid_one(j) = pid_from(bl)
              endif
           endif
        enddo
        i_tot=i
        j_tot=j
        k = i + j
        k_tot=i + j
        
        
        
        
        do bl = noprocs,1,-1
           two_tot(bl)=i/bl
           i=i-two_tot(bl)
        enddo
        do bl=noprocs,1,-1
           both_tot(bl)=k/bl
           k=k-both_tot(bl)
        enddo
        
        
        
        k=0
        k2=0
        k1=0
        
        
        
        do i=1,noprocs
           do j=1,two_tot(i)
              k=k+1
              k2=k2+1
              both(k)=two(k2)
              pid_both(k)=pid_two(k2)
              NT_both(k)=nodetype_all(both(k))
           enddo
           do j=1,both_tot(i)-two_tot(i)
              k=k+1
              k1=k1+1
              both(k)=one(k1)
              pid_both(k)=pid_one(k1)
              NT_both(k)=nodetype_all(both(k))
           enddo
        enddo
        
        bb = 0
        do bl = 1, total_no_blocks
           if(lrefine_all(bl).eq.g)then
              bb=bb+1  
              gridglobalID(bl) = both(bb)
              pid_from(bl) = pid_both(bb)
              nodetype_all(bl)=NT_both(bb)
           endif
        enddo
     enddo
  endif






  counter_proc = 0
  blocks_proc_level = 0
  do g = 1,lrefine_max
     level_total = 0
     temp_list_blocks = 0
     pid_from_level = 0
     pid_to_level = 0
     
     do bl = 1, total_no_blocks
        if(lrefine_all(bl)==g) then
           level_total = level_total + 1
           temp_list_blocks(level_total) = GridglobalID(bl)
           pid_from_level(level_total) = pid_from(bl)
        endif
     enddo
     totalnum_level = level_total

     do proc = noprocs,1,-1
        tgt(proc) = level_total /proc
        level_total = level_total - tgt(proc)
     end do



     proc_divider_level(1)=0
     do i = 1,noprocs
        proc_divider_level(i+1) = proc_divider_level(i)+ tgt(i)
     enddo

     do proc = 1, noprocs
        start = proc_divider_level(proc)+1
        ending = proc_divider_level(proc+1)
        do b = start, ending
              pid_to_level(b)= proc -1
        enddo
     enddo


     
 
     do bl = 1, totalnum_level
        counter_proc(pid_to_level(bl)+1)=counter_proc(pid_to_level(bl)+1)+1

        if(pid_from_level(bl)==pid) then
           new_loc(1, temp_list_blocks(bl)-displs(pid_from_level(bl)+1)) = counter_proc(pid_to_level(bl)+1)
           new_loc(2, temp_list_blocks(bl)-displs(pid_from_level(bl)+1)) = pid_to_level(bl)
        endif        
     enddo

     
  enddo



!if(num_called1==3) stop
  return
end subroutine ceg_block_balancing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ceg_set_blockstarts(pid, noprocs)

  use paramesh_dimensions
  use physicaldata
  use tree
  use paramesh_interfaces
  integer, intent(in) :: pid, noprocs
  integer :: locblocks(lrefine_max)
  integer :: g, b, count


  locblocks(:) = 0
  do g = 1, lrefine_max
     do b = 1, lnblocks
        if(lrefine(b)==g) then
           locblocks(g)=locblocks(g)+1
        endif
     enddo
  enddo

  block_starts(1) = 1
  do g = 1, lrefine_max
     block_starts(g+1) =  block_starts(g)+locblocks(g)
  end do

  do g = 1, lrefine_max
     count = 0
     do b = 1, lnblocks
        if(lrefine(b)==g) then
           block_list(count+block_starts(g)) = b
           count = count+1           
        endif
     enddo
  enddo

!if(pid==0) write(*,*) 'block_list',block_list,'in',pid
!if(pid==0) write(*,*) 'lrefine',lrefine(1:lnblocks),'in',pid


end subroutine ceg_set_blockstarts
! --------------------------------------------------------------------

! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

   SUBROUTINE  FindMinimum(x, Start, End, Location, num_called)
      IMPLICIT  NONE
      INTEGER, INTENT(IN)                :: Start, End, num_called
!      DOUBLE PRECISION, DIMENSION(1:End), INTENT(IN) :: x
      INTEGER, DIMENSION(1:End), INTENT(IN) :: x
      INTEGER, INTENT(INOUT)             :: Location
      REAL                               :: Minimum
      INTEGER                            :: i
      Minimum  = x(Start)          ! assume the first is the min
      Location = Start             ! record its position
      DO i = Start+1, End          ! start with next elements
         IF (x(i) < Minimum) THEN  !   if x(i) less than the min?
            Minimum  = x(i)        !      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
!      FindMinimum = Location            ! return the position
   END SUBROUTINE  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------
   SUBROUTINE  Swap(a, b, ai1, bi1, ai2, bi2, ai3, bi3,ai4,bi4)
      IMPLICIT  NONE
      INTEGER, INTENT(INOUT) :: a, b
      INTEGER, INTENT(INOUT) :: ai1, bi1, ai2, bi2, ai3, bi3, ai4, bi4
      INTEGER             :: Tempi


      Tempi = a
      a    = b
      b    = Tempi

      Tempi = ai1
      ai1    = bi1
      bi1    = Tempi

      Tempi = ai2
      ai2    = bi2
      bi2    = Tempi

      Tempi = ai3
      ai3    = bi3
      bi3    = Tempi

      Tempi = ai4
      ai4    = bi4
      bi4    = Tempi


   END SUBROUTINE  Swap
! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(x,y,z,w,v,Size, num_called, pid)
      IMPLICIT  NONE
      INTEGER, INTENT(IN)                   :: Size, num_called, pid
!      DOUBLE PRECISION, DIMENSION(1:SIZE), INTENT(INOUT) :: x
      INTEGER, DIMENSION(1:SIZE), INTENT(INOUT) :: x
      INTEGER, DIMENSION(1:SIZE), INTENT(INOUT) :: y,z,w,v
!      double precision,  DIMENSION(1:SIZE), INTENT(INOUT) :: u
      INTEGER                               :: i
      INTEGER                               :: Location
      DO i = 1, Size-1             ! except for the last
         CALL FindMinimum(x, i, Size, location,num_called)  ! find min from this to last
         CALL Swap(x(i),x(Location), y(i),y(Location), z(i),z(Location),w(i),w(Location),v(i),v(Location),w(Location))  ! swap this and the minimum
      END DO
   END SUBROUTINE  Sort
! --------------------------------------------------------------------
