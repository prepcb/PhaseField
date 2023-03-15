!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_mg_smooth                                                 REQUIRED
!!!!  * Main smoothing call
!!!!  * Loops over real variables calling pf_mg_smooth_var
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  pf_mg_smooth_var
!!!!  * Smoothing of each variable
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_mg_get_defect_var                                         REQUIRED
!!!!  * Calculates FAS MG defect, max defect and RMS residual
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Latest version: Chris Goodyer, May 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_mg_smooth(pid,level)
  use paramesh_dimensions
  use tree
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use paramesh_interfaces
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid,level
  integer lb, iopt
  iopt=1
  ! Cycle over all blocks
  do current_var=1,total_vars
     ! set the start point for the real variable and it's workspace data
     current_unk=1+(current_var-1)*nunkvbles
     current_work=1
!     do lb=1,lnblocks
      do lb=block_starts(level), block_starts(level+1)-1
         call pf_mg_smooth_var(lb)
     end do
     ! Update gcells for "current_unk" only
     gcell_on_cc(:)=.false.
     gcell_on_cc(current_unk)=.true.
!     call amr_guardcell(pid,iopt,nguard)
     call pf_guardcell(pid,iopt,level)
     gcell_on_cc(:)=.true.
  end do
end subroutine app_mg_smooth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pf_mg_smooth_var(lb)
  ! Smoother for node centered variables
  ! lb is the block number
  ! n is the variable number being smoothed
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  integer,intent(in) :: lb
  integer :: i,j,k
  double precision :: new_soln(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
  double precision :: dx,Fi,dFi_dvi
  double precision :: LapMult,tau_val
  double precision :: rfactor, rf1, rf2, rf3, dx2inv

  rfactor = dt/dtold
  rf1 = (rfactor+1.0)*dt/(2.0*rfactor+1.0)
  rf2 = (rfactor+1.0)*(rfactor+1.0)/(2.0*rfactor+1.0)
  rf3 = rfactor*rfactor/(2.0*rfactor+1.0)

  dx = bsize(1,lb)/real(nxb)
  dx2inv = 1.0/(dx*dx)

  do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
     do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
        do i=il_bnd+nguard,iu_bnd-nguard
           if (ndim.eq.3) then
              call get_stencil(i,j,k,dx,lb,1,Fi,LapMult)
              call get_stencil_diff(i,j,k,dx,lb,1,dFi_dvi,LapMult)
           else
              call get_stencil_2d(i,j,dx,lb,1,Fi,dx2inv,LapMult,tau_val)
              call get_stencil_diff_2d(i,j,dx,lb,1,dFi_dvi,LapMult,tau_val)
           endif

           if(mode_1.eq.1)then
              Fi = unk(current_unk,i,j,k,lb)-dt*Fi-unk(current_unk+1,i,j,k,lb)
              dFi_dvi = 1.0-dt*dFi_dvi
           else if(mode_1.eq.2)then
              Fi = unk(current_unk,i,j,k,lb)-rf1*Fi-unk(current_unk+1,i,j,k,lb)
! is actually
!              Fi=unk(current_unk,i,j,k,lb)-rf1*Fi-(rf2*unk(current_unk+3,i,j,k,lb) - rf3*unk(current_unk+4,i,j,k,lb)) 
              dFi_dvi = 1.0-rf1*dFi_dvi
           end if
           new_soln(i,j,k) = unk(current_unk,i,j,k,lb)-weight*Fi/dFi_dvi
        end do
     end do
  end do
  unk(current_unk,:,:,:,lb) = new_soln(:,:,:)
end subroutine pf_mg_smooth_var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_mg_get_defect_var(lb, max_defect, rms, npts)
  ! Defect for node centered variables
  ! lb is the block number
  ! n is the variable number whose defect is being calculated
  ! defect will be put into work(i,j,k,lb,n+1)
  use multigrid_parameters
  use paramesh_dimensions
  use time_dep_parameters
  use solution_parameters
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: lb
  integer, intent(inout) :: npts
  double precision, intent(inout) :: max_defect, rms
  integer :: i,j,k
  double precision :: dx,Fi, dx2inv
  double precision :: rfactor, rf1, rf2, rf3
  double precision :: LapMult,tau_val

  rfactor = dt/dtold
  rf1 = (rfactor+1.0)*dt/(2.0*rfactor+1.0)
  rf2 = (rfactor+1.0)*(rfactor+1.0)/(2.0*rfactor+1.0)
  rf3 = rfactor*rfactor/(2.0*rfactor+1.0)
!  rf1=dt
!  rf2=1.0
!  rf3=0.0

  dx = bsize(1,lb)/real(nxb)
  dx2inv = 1.0/(dx*dx)

  do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
     do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
        do i=il_bnd+nguard,iu_bnd-nguard
           if (ndim.eq.3) then
              call get_stencil(i,j,k,dx,lb,1,Fi,LapMult)
           else
              call get_stencil_2d(i,j,dx,lb,1,Fi,dx2inv,LapMult,tau_val)
           endif
           if(mode_1.eq.1)then
              Fi=unk(current_unk,i,j,k,lb)-dt*Fi-unk(current_unk+1,i,j,k,lb)
           else if(mode_1.eq.2)then
              Fi=unk(current_unk,i,j,k,lb)-rf1*Fi-unk(current_unk+1,i,j,k,lb)
           end if

           ! Add new residual to work(2).  
           !   Will be zero at start of coarsening.  
           !   This then sets the value that gets coarsened and then updates with defect
           work(i,j,k,lb,current_work+1) = -work(i,j,k,lb,current_work+1)-Fi


           if(Fi.lt.0.0)then
              Fi=-Fi
           end if
           if(Fi.gt.max_defect)then
              max_defect=Fi
           end if
           if (.not.has_children(lb)) then
              rms = rms + Fi*Fi
              npts = npts + 1
           endif
        end do
     end do
  end do

  return
end subroutine app_mg_get_defect_var

