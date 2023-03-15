!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_restrict_unk_genorder
!! NAME
!!
!!   amr_restrict_unk_genorder
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_unk_genorder (datain, dataout, order, ivar)
!!   Call amr_restrict_unk_genorder (real array, real array, integer, integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(in)    :: datain(:,:,:,:)  data to restrict
!!   Real,    Intent(inout) :: dataout(:,:,:,:)  data which is restricted and returned
!!   Integer, Intent(in)    :: order order of interpolating polynomial to use
!!   Integer, Intent(in)    :: ivar  variable number in unk to restrict
!!
!! INCLUDES
!! 
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!
!! CALLS
!! 
!! RETURNS
!!
!!   Restricted data returned in array 'dataout'.  
!!
!! DESCRIPTION
!!
!!   This routine performs interpolation for the restriction operation on
!!   cell centered data stored in 'unk'.  It uses a lagrange polynomial
!!   interpolation scheme.  Also, the interpolation stencil is automatically
!!   shifted to avoid using data in guardcells. CAUTION: you must realize that
!!   use of this routine with 'order' set to values higher than 1 MAY 
!!   cause asymmetric interpolation and your results may lose symmetry.
!!
!!   Data is passed in via the array 'datain' and returned in the array 
!!   'dataout'.  The order of the interpolating polynomial is also passed 
!!   in the variable 'order' and can take on value ranging from 1 to 5.  
!!   The last argument 'ivar' specifies which variable in 'unk' to apply 
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

      Subroutine amr_restrict_unk_user(datain,dataout,order,ivar,i_ind,j_ind) 

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata

      Implicit None

!-----Input/Output arguments.
      Real,    Intent(in)    :: datain(:,:,:,:)
      Real,    Intent(inout) :: dataout(:,:,:,:)
      Integer, Intent(in)    :: order, ivar,i_ind,j_ind

!-----Local arrays and variables.
      Integer :: i,j,k,ihalf
      Real :: parent(10,10),childLD(10,10),childLU(10,10),childRD(10,10),childRU(10,10)




!-----Being executable code.

      if(i_ind.eq.1.and.j_ind.eq.1)then
         do i=1,order
            do j=1,order
               childLD(i,j)=datain(ivar,i,j,1)
            enddo
         enddo
      elseif(i_ind.eq.1.and.j_ind.eq.2)then
         do i=1,order
            do j=1,order
         childLU(i,j)=datain(ivar,i,j,1)
            enddo
         enddo
      elseif(i_ind.eq.2.and.j_ind.eq.1)then
         do i=1,order
            do j=1,order
               childRD(i,j)=datain(ivar,i,j,1)
            enddo
         enddo
      else
         do i=1,order
            do j=1,order
               childRU(i,j)=datain(ivar,i,j,1)
            enddo
         enddo
      endif
      call  pcb_restrict(parent,childLD,childLU,childRD,childRU,order,i_ind,j_ind)
      dataout(ivar,:,:,:)=0d0
      ihalf=order/2-1
      do i=1,ihalf
         do j=1,ihalf
            dataout(ivar,2*i,2*j,1)=parent(i+1+(i_ind-1)*4,j+1+(j_ind-1)*4)
         enddo
      enddo


      return

      End Subroutine amr_restrict_unk_user
subroutine pcb_restrict(parent,childLD,childLU,childRD,childRU,order,i_ind,j_ind)

implicit none
Real, intent (in) :: childLD(10,10),childLU(10,10),childRD(10,10),childRU(10,10)
Real, intent (out):: parent(10,10)
integer, intent (in):: i_ind,j_ind
Real :: f(10,10)
integer, intent (in) :: order
integer :: i,j,ii,jj
Real :: y(10),xx,yy
Real :: y4(10)=(/-1.,-0.447213595499957939282,0.447213595499957939282,1.,0.,0.,0.,0.,0.,0./)
Real :: y6(10)=(/-1.,-0.765055323929464692851,-0.2852315164806450963142,0.2852315164806450963142,0.765055323929464692851,1.,0.,0.,0.,0./)
Real :: y8(10)=(/-1.,-0.8717401485096066153375,-0.5917001814331423021445,-0.2092992179024788687687,0.2092992179024788687687,0.5917001814331423021445,0.8717401485096066153375  ,1.,0.,0./)
Real :: y10(10)=(/-1.,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1./)
Real :: pcb_interpolate,pcb_LagrangePoly




if (i_ind.gt.2)then
   write(*,*)"first amr_restrict_unk_user", i_ind,j_ind
   stop
endif

if(order.eq.4)then
   y=y4
elseif(order.eq.6)then
   y=y6
elseif(order.eq.8)then
   y=y8
elseif(order.eq.10)then
   y=y10
else
   write(*,*)"order out of range in pcb_LagrangePoly: 4,6,8,10 only"
   stop
endif


parent=0.
!LD
if(i_ind.eq.1.and.j_ind.eq.1)then
do i =1, order
   do j=1,order
      f(i,j)=childLD(i,j)
   enddo
enddo
do i = 2,order/2
   xx = y(i)*2. -0.5*(y(1)+y(2))
   do j=2,order/2
      yy = y(j)*2. -0.5*(y(1)+y(2))
      parent(i,j)=pcb_interpolate(order,xx,yy,f)
   enddo
enddo

!LU
elseif(i_ind.eq.1.and.j_ind.eq.2)then

do i =1, order
   do j=1,order
      f(i,j)=childLU(i,j)
   enddo
enddo

do i = 2,order/2
   xx = y(i)*2. -0.5*(y(1)+y(2))
   do j=order/2+1,order-1
      yy = y(j)*2. -0.5*(y(order-1)+y(order))
      parent(i,j)=pcb_interpolate(order,xx,yy,f)
   enddo
enddo
!RD
elseif(i_ind.eq.2.and.j_ind.eq.1)then

do i =1, order
   do j=1,order
      f(i,j)=childRD(i,j)
   enddo
enddo
do i = order/2+1,order-1
   xx = y(i)*2. -0.5*(y(order-1)+y(order))
   do j=  2,order/2
      yy = y(j)*2. -0.5*(y(1)+y(2))
      parent(i,j)=pcb_interpolate(order,xx,yy,f)
   enddo
enddo
!RU
elseif(i_ind.eq.2.and.j_ind.eq.2)then

do i =1, order
   do j=1,order
      f(i,j)=childRU(i,j)
   enddo
enddo
do i = order/2+1,order-1
   xx = y(i)*2. -0.5*(y(order-1)+y(order))
   do j=  order/2+1,order-1
      yy = y(j)*2. -0.5*(y(order-1)+y(order))
      parent(i,j)=pcb_interpolate(order,xx,yy,f)
   enddo
enddo
else
write(*,*)"error in restrict",i_ind,j_ind
stop
endif


end subroutine pcb_restrict


subroutine pcb_restrict_X(parent,childLD,childLU,childRD,childRU,order,i_ind,j_ind)

implicit none
Real, intent (in) :: childLD(10,10),childLU(10,10),childRD(10,10),childRU(10,10)
Real, intent (out):: parent(10,10)
integer, intent (in):: i_ind,j_ind
Real :: f(10,10)
integer, intent (in) :: order
integer :: i,j,ii,jj
Real :: y(10),xx,yy
Real :: y4(10)=(/-1.,-0.447213595499957939282,0.447213595499957939282,1.,0.,0.,0.,0.,0.,0./)
Real :: y6(10)=(/-1.,-0.765055323929464692851,-0.2852315164806450963142,0.2852315164806450963142,0.765055323929464692851,1.,0.,0.,0.,0./)
Real :: y8(10)=(/-1.,-0.8717401485096066153375,-0.5917001814331423021445,-0.2092992179024788687687,0.2092992179024788687687,0.5917001814331423021445,0.8717401485096066153375  ,1.,0.,0./)
Real :: y10(10)=(/-1.,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1./)
Real :: pcb_interpolate,pcb_LagrangePoly
Real :: MatrxInt(4,10)=reshape((/ -0.1258726,    -0.0001992,     0.0239232,     0.0112625,     0.6165841,     0.0005766,    -0.0634861,    -0.0289813,     0.6425251,    -0.0012672,     0.1014969,     0.0424743,    -0.2033408,     0.9999972,    -0.1734516,    -0.0594062,     0.1154220,     0.0012793,     0.4579133,     0.0864772,    -0.0774335,    -0.0006215,     0.7917217,    -0.1431766,     0.0557597,     0.0003948,    -0.2087598,     0.3558988,    -0.0408415,    -0.0002721,     0.1143992,     0.8627066,     0.0282123,     0.0001824,    -0.0698899,    -0.1841350,    -0.0110148,    -0.0000704,     0.0261332,     0.0568797/),(/4,10/))



if (i_ind.gt.2)then
   write(*,*)"first amr_restrict_unk_user", i_ind,j_ind
   stop
endif

if(order.eq.4)then
   y=y4
elseif(order.eq.6)then
   y=y6
elseif(order.eq.8)then
   y=y8
elseif(order.eq.10)then
   y=y10
else
   write(*,*)"order out of range in pcb_LagrangePoly: 4,6,8,10 only"
   stop
endif


parent=0.
!LD
if(i_ind.eq.1.and.j_ind.eq.1)then
do i =1, order
   do j=1,order
      f(i,j)=childLD(i,j)
   enddo
enddo
do i = 2,order/2
   xx = 1.0 + y(i)*2.
   do j=2,order/2
      yy = 1.0 + y(j)*2.
      parent(i,j)=pcb_interpolate(order,xx,yy,f)
   enddo
enddo

!LU
elseif(i_ind.eq.1.and.j_ind.eq.2)then

do i =1, order
   do j=1,order
      f(i,j)=childLU(i,j)
   enddo
enddo

do i = 2,order/2
   xx = 1 + y(i)*2.
   do j=order/2+1,order-1
      yy = -1.0 + y(j)*2.
      parent(i,j)=pcb_interpolate(order,xx,yy,f)
   enddo
enddo
!RD
elseif(i_ind.eq.2.and.j_ind.eq.1)then

do i =1, order
   do j=1,order
      f(i,j)=childRD(i,j)
   enddo
enddo
do i = order/2+1,order-1
   xx = -1.0 + y(i)*2.
   do j=  2,order/2
      yy = 1.0 + y(j)*2.
      parent(i,j)=pcb_interpolate(order,xx,yy,f)
   enddo
enddo
!RU
elseif(i_ind.eq.2.and.j_ind.eq.2)then

do i =1, order
   do j=1,order
      f(i,j)=childRU(i,j)
   enddo
enddo
do i = order/2+1,order-1
   xx = -1.0 + y(i)*2.
   do j=  order/2+1,order-1
      yy = -1.0 + y(j)*2.
      parent(i,j)=pcb_interpolate(order,xx,yy,f)
   enddo
enddo
else
write(*,*)"error in restrict",i_ind,j_ind
stop
endif


end subroutine pcb_restrict_X



Real function pcb_interpolate(order,x,y,f)
implicit none
integer, intent (in) :: order
Real, intent (in) :: x,y,f(10,10)
integer :: i,j
Real :: pcb_LagrangePoly

if(order.eq.4)then
elseif(order.eq.6)then
elseif(order.eq.8)then
elseif(order.eq.10)then
else
   write(*,*)"order out of range in pcb_LagrangePoly: 4,6,8,10 only"
   stop   
endif
pcb_interpolate=0d0
do i=1,order
   do j=1,order
      pcb_interpolate = pcb_interpolate + f(i,j)*pcb_LagrangePoly(order,x,i,0)*pcb_LagrangePoly(order,y,j,0) 
   enddo
enddo



end function pcb_interpolate



Real function pcb_LagrangePoly0(order,x,j,k,l)

implicit none
integer, intent (in) :: order,j,k,l
Real, intent (in) :: x
integer :: i
Real :: y(10)
Real :: y4(10)=(/-1.,-0.447213595499957939282,0.447213595499957939282,1.,0.,0.,0.,0.,0.,0./)
Real :: y6(10)=(/-1.,-0.765055323929464692851,-0.2852315164806450963142,0.2852315164806450963142,0.765055323929464692851,1.,0.,0.,0.,0./)
Real :: y8(10)=(/-1.,-0.8717401485096066153375,-0.5917001814331423021445,-0.2092992179024788687687,0.2092992179024788687687,0.5917001814331423021445,0.8717401485096066153375  ,1.,0.,0./)
Real :: y10(10)=(/-1.,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1./)

if(order.eq.4)then
   y=y4
elseif(order.eq.6)then
   y=y6
elseif(order.eq.8)then
   y=y8
elseif(order.eq.10)then
   y=y10
else
   write(*,*)"order out of range in pcb_LagrangePoly: 4,6,8,10 only"
   stop
endif


pcb_LagrangePoly0=1d0
do i = 1,order
   if(i.ne.j.and.i.ne.k .and. i.ne.l)pcb_LagrangePoly0 = pcb_LagrangePoly0*(x-y(i))
enddo

end function pcb_LagrangePoly0

Real function pcb_LagrangePoly(order,x,j)

implicit none
integer, intent (in) :: order,j
Real, intent (in) :: x
integer :: i
Real :: y(10),pcb_LagrangePoly0
Real :: y4(10)=(/-1.,-0.447213595499957939282,0.447213595499957939282,1.,0.,0.,0.,0.,0.,0./)
Real :: y6(10)=(/-1.,-0.765055323929464692851,-0.2852315164806450963142,0.2852315164806450963142,0.765055323929464692851,1.,0.,0.,0.,0./)
Real :: y8(10)=(/-1.,-0.8717401485096066153375,-0.5917001814331423021445,-0.2092992179024788687687,0.2092992179024788687687,0.5917001814331423021445,0.8717401485096066153375  ,1.,0.,0./)
Real :: y10(10)=(/-1.,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1./)

if(order.eq.4)then
   y=y4
elseif(order.eq.6)then
   y=y6
elseif(order.eq.8)then
   y=y8
elseif(order.eq.10)then
   y=y10
else
   write(*,*)"order out of range in pcb_LagrangePoly: 4,6,8,10 only"
   stop
endif
pcb_LagrangePoly = pcb_LagrangePoly0(order,x,j,0,0)/pcb_LagrangePoly0(order,y(j),j,0,0)
end function pcb_LagrangePoly
