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
!!!!  Latest version: Chris Goodyer, September 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_mg_smooth(pid,level)
  use paramesh_dimensions
  use tree
  use multigrid_parameters
  use generic_parameters
  use paramesh_dimensions
  use physicaldata
  use paramesh_interfaces
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, level
  integer lb, iopt
  double precision :: AA(3,3),BB(3)
  iopt=1

  mg_update_total = 0.0
  AA=0.
  BB=0.
  ! Cycle over all blocks
  do current_var=1,total_vars
     ! set the start point for the real variable and it's workspace data
     current_unk=1+(current_var-1)*nunkvbles
     current_work=1
!$OMP PARALLEL SHARED(level) PRIVATE(lb)
     if (multigrid_on.ge.1) then
!$OMP DO SCHEDULE(DYNAMIC,16)
        do lb=block_starts(level), block_starts(level+1)-1
            call pf_mg_smooth_var(lb,AA,BB)
        end do
!$OMP END DO  NOWAIT
     else
!$OMP DO SCHEDULE(DYNAMIC,16)
        do lb=1, lnblocks
            if (.not.has_children(lb)) call pf_mg_smooth_var(lb)
        end do
!$OMP END DO  NOWAIT
     endif
!$OMP END PARALLEL 

     ! Update gcells for "current_unk" only
     gcell_on_cc(:)=.false.
     gcell_on_cc(current_unk)=.true.
     if (multigrid_on.ge.1) then
        call pf_guardcell(pid,iopt,level)
     else
        call amr_guardcell(pid,iopt,nguard)
     endif
     gcell_on_cc(:)=.true.
  end do
end subroutine app_mg_smooth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pf_mg_smooth_var(lb,AA,BB)
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
  double precision, intent(inout):: AA(3,3),BB(3)
  integer,intent(in) :: lb
  integer :: i,j,k,ii
  double precision :: new_soln(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,3)
  double precision :: dx,Fi,dFi_dvi
  double precision :: LapMult,tau_val
  double precision :: rfactor, rf1, rf2, rf3, dx2inv
  double precision :: CC(3)

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
              dFi_dvi = exp(-dt*dFi_dvi)
           else if(mode_1.eq.2)then
              Fi = unk(current_unk,i,j,k,lb)-rf1*Fi-unk(current_unk+1,i,j,k,lb)
              dFi_dvi = exp(-rf1*dFi_dvi)
           end if
!            AA(current_var,current_var)=dFi_dvi
!            BB(current_var)=Fi
           
!            new_soln(i,j,k) = unk(current_unk,i,j,k,lb)-weight*Fi/dFi_dvi
!            if(current_var.eq.3)then
!            
!             call TLinear(AA,BB,CC)
!             do ii=1,3
!             new_soln(i,j,k,ii) = unk(1+(ii-1)*nunkvbles,i,j,k,lb) - weight*CC(ii)
!             enddo
!            endif
           new_soln(i,j,k,current_var) = unk(1+(current_var-1)*nunkvbles,i,j,k,lb) - weight*Fi/dFi_dvi
           mg_update_total = mg_update_total + ABS(Fi/dFi_dvi)
        end do
     end do
  end do
!   if(current_var.eq.3)then
  do ii=1,3
    unk(1+(ii-1)*nunkvbles,:,:,:,lb) = new_soln(:,:,:,ii)
  enddo
!   endif
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

           rms = rms + Fi*Fi
           npts = npts + 1
        end do
     end do
  end do

  return
end subroutine app_mg_get_defect_var
!*********************************************************
!*   SOLVE A LINEAR SYSTEM BY TRIANGULARIZATION METHOD   *
!* ----------------------------------------------------- *
!* SAMPLE RUN:                                           *
!*                                                       *
!* Linear system AX = B                                  *
!*                                                       *
!* Matrix A:                                             *
!* 2.00    1.00   -2.00                                  *
!* 3.00   -1.00    1.00                                  *
!* 7.00    5.00   -3.00                                  *
!*                                                       *
!* Right side:                                           *
!* B(1) =    1.00                                        *
!* B(2) =    0.00                                        *
!* B(3) =    0.00                                        *
!*                                                       *
!* Solution of AX=B:                                     *
!* X(1) =    0.0625                                      *
!* X(2) =   -0.0500                                      *
!* X(3) =   -0.6875                                      *
!*                                                       *
!* ----------------------------------------------------- *
!* Reference: "Méthodes de calcul numérique - Tome 1 By  *
!*             Claude Nowakowski, PS1 1981" [BIBLI 07].  *
!*                                                       *
!*                    F90 Release By J-P Moreau, Paris.  *
!*********************************************************
subroutine Tlinear(AA,BB,X)
double precision,intent(inout) :: AA(3,3),BB(3)
double precision,intent(out) :: X(3)
double precision ::  A(3,3),B(3)
INTEGER, PARAMETER :: N = 3
do i=1,N
B(i)=BB(i)
do j=1,N
A(i,j)=AA(i,j)
enddo
enddo


! Transform A into triangular matrix
DO K = 1, N-1
  DO I = K + 1, N
    B(I) = B(I) - A(I, K) / A(K, K) * B(K)
    DO J = K + 1, N
      A(I, J) = A(I, J) - A(I, K) / A(K, K) * A(K, J)
    END DO
  END DO
END DO

! Solve triangular system
X(N) = B(N) / A(N, N)
DO I = N - 1, 1, -1
  S = 0.d0
  DO K = I + 1, N
    S = S + A(I, K) * X(K)
  END DO
  X(I) = (B(I) - S) / A(I, I)
END DO
X(4)=BB(4)/AA(4,4)
X(5)=BB(5)/AA(5,5)

END
