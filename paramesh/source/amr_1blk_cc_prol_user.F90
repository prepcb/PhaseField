!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_cc_prol_user
!! NAME
!!  
!!   amr_1blk_cc_prol_user
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_cc_prol_user(recv,ia,ib,ja,jb,ka,kb,
!!                                idest,ioff,joff,koff,mype,ivar)
!!   Call amr_1blk_cc_prol_user(real,
!!                                integer, integer, integer, integer,
!!                                integer, integer, integer, integer,
!!                                integer, integer, integer, integer)
!!
!!
!! ARGUMENTS
!!
!!  Real,    intent(inout) :: recv(:,:,:,:)
!!    Data array holding the data extracted from unk which will be prolonged
!!    and placed into the unk1 array.
!!
!!  Integer, Intent(in) :: ia,ib,ja,jb,ka,kb
!!    Integers which control the limits into unk1 where the prolonged data
!!    will be placed.
!!
!!  Integer, Intent(in) :: idest,ioff,joff,koff,mype
!!    Idest controls which 'layer' into which the prolonged data will be
!!    placed in unk1.  ioff, joff and koff are offsets.  mype is is the
!!    local processor id.
!!
!!  Integer, intent(in) :: ivar
!!    ivar is the variable number in unk which is prolonged.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!! 
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   prolong_arrays
!!   
!! CALLS
!! 
!!   amr_abort
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit the data in recv is interpolated
!!   and placed into the unk1 array.
!!
!! DESCRIPTION
!!
!!   This routine takes data from the array recv, originally extracted 
!!   from the solution array unk, and performs a prolongation operation 
!!   on it, between the bounds ranges ia to ib, ja to jb, and ka to kb. 
!!   The data in recv is from a parent block and the
!!   result of the prolongation operation is written directly into one
!!   layer of the working block array unk1(...,idest).
!!   The position of the child within the parent block is specified by 
!!   the ioff, joff and koff arguments.
!!
!!   This particular prolongation is simple user interpolation. It can
!!   be used for blocks with an even or odd number of grid cells.
!!
!!   It is applied to all UNK variables whose corresponding element
!!   of interp_mask is set to 1.
!!
!!   Conservative prolongation. Special treatment for the  cells immediately
!!   adjacent to a boundary 
!!   (ie i=nguard,nguard+1,iu_bnd1-nguard,iu_bnd1-nguard+1
!!   and likewise for j and k indeces) 
!!   if using an even number of grid cells per block along that axis. 
!!   No special treatment is required when the number of cells is odd.
!!
!!   Note: before using this routine in your program, make sure that the
!!   routine prolong_unk_fun_init has been called.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          January 1997
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
    
#include "paramesh_preprocessor.fh"

      Subroutine amr_1blk_cc_prol_user               & 
        (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,  & 
         mype,ivar,order)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use prolong_arrays

      Implicit None

!-----Input/Output Variables
      Real,    Intent(inout) :: recv(:,:,:,:)
      Integer, Intent(in)    :: ia,ib,ja,jb,ka,kb
      Integer, Intent(in)    :: idest,ioff,joff,koff,mype
      Integer, Intent(in)    :: ivar,order

!-----Local variables
      Real    :: dx,dy,dz,cx,cy,cz
      Integer :: icl,icu,jcl,jcu,kcl,kcu,i_ind,j_ind,k_ind
      Integer :: i,j,k,i1,j1,k1,i1p,j1p,k1p

      Real :: childLD(10,10),childLU(10,10),childRD(10,10),childRU(10,10)
      Real :: pcb_parent(10,10)
      Logical :: iboundary=.false.




childLD=0.
childLU=0.
childRD=0.
childRU=0.

!-----Begin exectable code

!      write(*,*)"in amr_1blk_cc_prol_user "
!      stop

      If (prol_init.ne.100) Then
       Write(*,*) 'PARAMESH ERROR !'
       Write(*,*) 'Error : prolong_gen_unk_fun. ',  & 
             'You must call amr_initialize ',       & 
             'before you can use this routine!'
       Call amr_abort
      End If  ! End If (prol_init.ne.100)


! Set the bounds on the loop controlling the interpolation.



      icl=ia
      icu=ib
      jcl=ja
      jcu=jb
      kcl=ka
      kcu=kb

      i_ind = 1
      j_ind = 1
      k_ind = 1
      If (ioff.gt.0) i_ind = 2
      If (joff.gt.0) j_ind = 2
      If (koff.gt.0) k_ind = 2

      do i=1,order
         do j=1,order
            pcb_parent(i,j) = recv(ivar,i,j,1)
         enddo
      enddo
      

      call pcb_prolong(pcb_parent,childLD,childLU,childRD,childRU,order,i_ind,j_ind,ia,ib,ja,jb)


!-----Interpolation loop.


      if(i_ind.eq.1.and.j_ind.eq.1)then
         Do j=ja,jb
            Do i=ia,ib
               unk1(ivar,i,j,1,idest) =childLD(i,j)
            enddo
         enddo
      elseif(i_ind.eq.1.and.j_ind.eq.2)then
         Do j=ja,jb
            Do i=ia,ib
               unk1(ivar,i,j,1,idest) =childLU(i,j)
            enddo
         enddo
      elseif(i_ind.eq.2.and.j_ind.eq.1)then
         Do j=ja,jb
            Do i=ia,ib
               unk1(ivar,i,j,1,idest) =childRD(i,j)
            enddo
         enddo
      else
         Do j=ja,jb
            Do i=ia,ib
               unk1(ivar,i,j,1,idest) =childRU(i,j)
            enddo
         enddo
      End if




      Return
      End Subroutine amr_1blk_cc_prol_user

subroutine pcb_prolong(parent,childLD,childLU,childRD,childRU,order,i_ind,j_ind,ia,ib,ja,jb)

implicit none
Real, intent (out) :: childLD(10,10),childLU(10,10),childRD(10,10),childRU(10,10)
Real, intent (in):: parent(10,10)
Real :: f(10,10)
integer, intent (in) :: order,i_ind,j_ind,ia,ib,ja,jb
integer :: i,j,ii,jj
Real :: y(10),xx,yy
Real :: y4(10)=(/-1.,-0.447213595499957939282,0.447213595499957939282,1.,0.,0.,0.,0.,0.,0./)
Real :: y6(10)=(/-1.,-0.765055323929464692851,-0.2852315164806450963142,0.2852315164806450963142,0.765055323929464692851,1.,0.,0.,0.,0./)
Real :: y8(10)=(/-1.,-0.8717401485096066153375,-0.5917001814331423021445,-0.2092992179024788687687,0.2092992179024788687687,0.5917001814331423021445,0.8717401485096066153375  ,1.,0.,0./)
Real :: y10(10)=(/-1.,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1./)
Real :: yy10(10)=(/-1.,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1./)
Real :: pcb_interpolate,pcb_LagrangePoly
Real :: Child(10,10),Child1(10,10)

if(order.eq.4)then
   y=y4
elseif(order.eq.6)then
   y=y6
elseif(order.eq.8)then
   y=y8
elseif(order.eq.10)then
   y=yy10
else
   write(*,*)"order out of range in pcb_LagrangePoly: 4,6,8,10 only"
   stop
endif

f = parent

   if(i_ind.eq.1.and.j_ind.eq.1)then
      do i=ia,ib
         xx = 0.5*y(i) + 0.25*(y(1)+y(2))
         do j=ja,jb
            yy = 0.5*y(j) + 0.25*(y(1)+y(2))
            childLD(i,j) =  pcb_interpolate(order,xx,yy,f)
         enddo
      enddo


      return
   elseif(i_ind.eq.1.and.j_ind.eq.2)then
      do i=ia,ib
         xx = 0.5*y(i) + 0.25*(y(1)+y(2))
         do j=ja,jb
           yy = 0.5*y(j) + 0.25*(y(order-1)+y(order))
            childLU(i,j) =  pcb_interpolate(order,xx,yy,f)
         enddo
      enddo
      return
   elseif(i_ind.eq.2.and.j_ind.eq.1)then
      do i=ia,ib
         xx = 0.5*y(i) + 0.25*(y(order-1)+y(order))
         do j=ja,jb
            yy = 0.5*y(j) + 0.25*(y(1)+y(2))
            childRD(i,j) =  pcb_interpolate(order,xx,yy,f)
         enddo
      enddo
      return
   else
      do i=ia,ib
         xx = 0.5*y(i) + 0.25*(y(order-1)+y(order))
         do j=ja,jb
            yy =   0.5*y(j) + 0.25*(y(order-1)+y(order))
            childRU(i,j) =  pcb_interpolate(order,xx,yy,f)
         enddo
      enddo
   endif





end subroutine pcb_prolong


subroutine pcb_prolong_X(parent,childLD,childLU,childRD,childRU,order,i_ind,j_ind,ia,ib,ja,jb)

implicit none
Real, intent (out) :: childLD(10,10),childLU(10,10),childRD(10,10),childRU(10,10)
Real, intent (in):: parent(10,10)
Real :: f(10,10)
integer, intent (in) :: order,i_ind,j_ind,ia,ib,ja,jb
integer :: i,j,ii,jj
Real :: y(10),xx,yy
Real :: y4(10)=(/-1.,-0.447213595499957939282,0.447213595499957939282,1.,0.,0.,0.,0.,0.,0./)
Real :: y6(10)=(/-1.,-0.765055323929464692851,-0.2852315164806450963142,0.2852315164806450963142,0.765055323929464692851,1.,0.,0.,0.,0./)
Real :: y8(10)=(/-1.,-0.8717401485096066153375,-0.5917001814331423021445,-0.2092992179024788687687,0.2092992179024788687687,0.5917001814331423021445,0.8717401485096066153375  ,1.,0.,0./)
Real :: y10(10)=(/-1.,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1./)
Real :: yy10(10)=(/-1.,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1./)
Real :: pcb_interpolate,pcb_LagrangePoly
Real :: Child(10,10),Child1(10,10)
!Real :: MatrxInt(8,10)=reshape((/ 0.0120621,     0.0120621,     0.0120621,     0.0120621,     0.0120621,     0.0120621,     0.0120621,     0.0120621,    -0.0352092,    -0.0352092,    -0.0352092,    -0.0352092,    -0.0352092,    -0.0352092,    -0.0352092,    -0.0352092,     0.0803518,     0.0803518,     0.0803518,     0.0803518,     0.0803518,     0.0803518,     0.0803518,     0.0803518,     0.9904206,     0.9904206,     0.9904206,     0.9904206,     0.9904206,     0.9904206,     0.9904206,     0.9904206,    -0.0691747,    -0.0691747,    -0.0691747,    -0.0691747,    -0.0691747,    -0.0691747,    -0.0691747,    -0.0691747,     0.0348038,     0.0348038,     0.0348038,     0.0348038,     0.0348038,     0.0348038,     0.0348038,     0.0348038,    -0.0223571,    -0.0223571,    -0.0223571,    -0.0223571,    -0.0223571,    -0.0223571,    -0.0223571,    -0.0223571,     0.0154878,     0.0154878,     0.0154878,     0.0154878,     0.0154878,     0.0154878,     0.0154878,     0.0154878,    -0.0104058,    -0.0104058,    -0.0104058,    -0.0104058,    -0.0104058,    -0.0104058,    -0.0104058,    -0.0104058,     0.0040207,     0.0040207,     0.0040207,     0.0040207,     0.0040207,     0.0040207,     0.0040207,     0.0040207/),(/8,10/))


if(order.eq.10)then
   yy10(order)=2.-y10(order)
   yy10(1)=-yy10(order)
elseif(order.eq.6)then
   yy10(order)=2.-y6(order)
   yy10(1)=-yy10(order)
elseif(order.eq.4)then
   yy10(order)=2.-y4(order)
   yy10(1)=-yy10(order)
else
   write(*,*)" in pcb_prolong. order=",order 
   stop
endif
      
      



if(order.eq.4)then
   y=y4
elseif(order.eq.6)then
   y=y6
elseif(order.eq.8)then
   y=y8
elseif(order.eq.10)then
   y=yy10
else
   write(*,*)"order out of range in pcb_LagrangePoly: 4,6,8,10 only"
   stop
endif

f = parent
i=1
do j=2,order-1
   f(i,j) = 0.5*(parent(i,j)+parent(i+1,j))
   f(j,i) = 0.5*(parent(j,i)+parent(j,i+1))
enddo
i=order
do j=2,order-1
   f(i,j) = 0.5*(parent(i,j)+parent(i-1,j))
   f(j,i) = 0.5*(parent(j,i)+parent(j,i-1))
enddo
f(1,1) =         0.25*(parent(1,1)+        parent(1,2)+          parent(2,1)+           parent(2,2))
f(1,order)     = 0.25*(parent(1,order)+    parent(1,order-1)+    parent(2,order)+       parent(2,order-1))
f(order,1)     = 0.25*(parent(order,1)+    parent(order-1,1)+    parent(order,2)+       parent(order-1,2))
f(order,order) = 0.25*(parent(order,order)+parent(order,order-1)+parent(order-1,order)+ parent(order-1,order-1))



   if(i_ind.eq.1.and.j_ind.eq.1)then
      do i=ia,ib
         xx = (y(i)-1.)*0.5 
         do j=ja,jb
            yy =  (y(j)-1.)*0.5 
            childLD(i,j) =  pcb_interpolate(order,xx,yy,f)
         enddo
      enddo


      return
   elseif(i_ind.eq.1.and.j_ind.eq.2)then
      do i=ia,ib
         xx = (y(i)-1)*0.5 
         do j=ja,jb
            yy =  (y(j)+1.)*0.5 
            childLU(i,j) =  pcb_interpolate(order,xx,yy,f)
         enddo
      enddo
      return
   elseif(i_ind.eq.2.and.j_ind.eq.1)then
      do i=ia,ib
         xx = (y(i)+1)*0.5 
         do j=ja,jb
            yy =  (y(j)-1.)*0.5
               childRD(i,j) =  pcb_interpolate(order,xx,yy,f)
         enddo
      enddo
      return
   else
      do i=ia,ib
         xx = (y(i)+1)*0.5 
         do j=ja,jb
            yy =  (y(j)+1.)*0.5 
            childRU(i,j) =  pcb_interpolate(order,xx,yy,f)
         enddo
      enddo
   endif





end subroutine pcb_prolong_X


