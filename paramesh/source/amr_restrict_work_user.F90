!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_restrict_work_genorder
!! NAME
!!
!!   amr_restrict_work_genorder
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_work_genorder(datainw,dataoutw,order)
!!   Call amr_restrict_work_genorder(real array, real array, integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(in)    :: datain(:,:,:,:)  data to restrict                         
!!   Real,    Intent(inout) :: dataout(:,:,:,:)  data which is restricted and returned   
!!   Integer, Intent(in)    :: order order of interpolating polynomial
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   workspace
!!
!! CALLS
!!
!! RETURNS
!!
!!   Restricted data returned in array 'dataoutw'.
!!
!! DESCRIPTION
!!
!!   This routine performs interpolation for the restriction operation on
!!   cell centered data stored in 'work'.  It uses a lagrange polynomial
!!   interpolation scheme.  Also, the interpolation stencil is automatically
!!   shifted to avoid using data in guardcells. CAUTION: you must realize that
!!   use of this routine with 'order' set to values higher than 1 MAY
!!   cause asymmetric interpolation and your results may lose symmetry.
!!
!!   Data is passed in in the array 'datainw' and returned in the array
!!   'dataoutw'.  The order of the interpolating polynomial is also passed
!!   in the variable 'order' and can take on value ranging from 1 to 5.
!!   The last argument 'ivar' specifies which variable in 'work' to apply
!!   the interpolation to.
!!
!! AUTHORS
!!
!!   Written :     Kevin Olson          March 2004
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_restrict_work_user(datainw,dataoutw,iopt,order,i_ind,j_ind)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use workspace

      Implicit None

!-----Input/Output arguments.
      Real,    Intent(in)    :: datainw(:,:,:)
      Real,    Intent(inout) :: dataoutw(:,:,:)
      Integer, Intent(in)    :: iopt, order,i_ind,j_ind




!-----Local arrays and variables.
      Integer :: i,j,k
      Real :: parent(10,10),childLD(10,10),childLU(10,10),childRD(10,10),childRU(10,10)
      Real :: s


      if(i_ind.ne.1.and.i_ind.ne.2)then
         write(*,*)"amr_restrict_work_user",i_ind,j_ind
         stop
      endif

!-----Being executable code.
      if(i_ind.eq.1.and.j_ind.eq.1)childLD(:,:)=datainw(:,:,1)
      if(i_ind.eq.1.and.j_ind.eq.2)childLU(:,:)=datainw(:,:,1)
      if(i_ind.eq.2.and.j_ind.eq.1)childRD(:,:)=datainw(:,:,1)
      if(i_ind.eq.2.and.j_ind.eq.2)childRU(:,:)=datainw(:,:,1)


      parent=0d0
      call  pcb_restrict(parent,childLD,childLU,childRD,childRU,10,i_ind,j_ind)

      dataoutw(:,:,:)=0d0
      do i=1,4
         do j=1,4
            dataoutw(2*i,2*j,1)=parent(i+1+(i_ind-1)*4,j+1+(j_ind-1)*4)
         enddo
      enddo

      Return
      End Subroutine amr_restrict_work_user
