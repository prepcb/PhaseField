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
  use time_dep_parameters
  use solution_parameters
  use  generic_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid,level
  integer lb,lbs, iopt,i,j,k
  integer :: ii,jj,kk
  logical Not_number
  double precision:: new_soln(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,N_unk),dx
  double precision:: block_defect





  if(level.gt.-1)then
     iopt=1
     
     
     
     
     g_time=time
     
     
     ! Cycle over all blocks
     do lbs=block_starts(level), block_starts(level+1)-1
        if(morton_flag) then
           lb = block_list(lbs)
        else
           lb = lbs
        endif
           dx = bsize(1,lb)/real(nxb)
           
           
           
           call pf_mg_smooth_var(lb,level,new_soln,pid)

           do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d            
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
                 do i=il_bnd+nguard,iu_bnd-nguard
                    if(new_soln(i,j,k,1).le.0d0.or.new_soln(i,j,k,1).ge.0d0)then
                    else
                       write(*,*)"in smooth",i,j,lb,dx
                       stop
                    endif
                    if(abs(new_soln(i,j,k,1)).lt.1d-10)new_soln(i,j,k,1)=0d0
                    if(abs(1d0-new_soln(i,j,k,1)).lt.1d-10)new_soln(i,j,k,1)=1d0
                    if(.not.g_plapp)new_soln(i,j,k,2) = max(0d0,min(1d0,new_soln(i,j,k,2)))
                    do current_var=1,N_unk
                       unk(1+(current_var-1)*nunkvbles,i,j,k,lb)=new_soln(i,j,k,current_var)
                    enddo
                 enddo
              enddo
           enddo
     end do
  else
!!$     write(*,*)"mp_smooth.f90"
!!$     stop
     do lbs=1,lnblocks
        if(morton_flag) then
           lb = block_list(lbs)
        else
           lb = lbs
        endif


        if (.not.has_children(lb))then

!           call pf_mg_smooth_var(lb,level,new_soln,pid)
!           write(*,*)"g_max_lvl",g_max_lvl
           call pf_mg_smooth_var(lb,g_max_lvl,new_soln,pid)
!           write(*,*)"after smooth",lrefine(lb),lb

           do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d            
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
                 do i=il_bnd+nguard,iu_bnd-nguard
                    if(new_soln(i,j,k,1).le.0d0.or.new_soln(i,j,k,1).ge.0d0)then
                    else
                       write(*,*)"in smooth",i,j,lb,dx
                       stop
                    endif
                    
                    if(abs(new_soln(i,j,k,1)).lt.1d-10)new_soln(i,j,k,1)=0d0
                    if(abs(1d0-new_soln(i,j,k,1)).lt.1d-10)new_soln(i,j,k,1)=1d0
                    if(.not.g_plapp)new_soln(i,j,k,2) = max(0d0,min(1d0,new_soln(i,j,k,2)))
                    do current_var=1,N_unk
                       unk(1+(current_var-1)*nunkvbles,i,j,k,lb)=new_soln(i,j,k,current_var)
                    enddo
                 enddo
              enddo
           enddo


!!$           do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d            
!!$              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
!!$                 do i=il_bnd+nguard,iu_bnd-nguard
!!$                    if(new_soln(i,j,k,1).le.0d0.or.new_soln(i,j,k,1).ge.0d0)then
!!$                    else
!!$                       write(*,*)"in smooth",i,j,k,lb,dx
!!$                       stop
!!$                    endif
!!$                    
!!$                    if(abs(new_soln(i,j,k,1)).lt.1d-10)new_soln(i,j,k,1)=0d0
!!$                    if(abs(1d0-new_soln(i,j,k,1)).lt.1d-10)new_soln(i,j,k,1)=1d0
!!$                    new_soln(i,j,k,2) = max(0d0,min(1d0,new_soln(i,j,k,2)))
!!$                    do current_var=1,N_unk
!!$                       !                    if(.not.Not_number(new_soln(i,j,1,current_var)))unk(1+(current_var-1)*nunkvbles,i,j,1,lb)=new_soln(i,j,1,current_var)
!!$                       if(.not.Not_number(new_soln(i,j,k,current_var)))unk(1+(current_var-1)*nunkvbles,i,j,k,lb)=new_soln(i,j,k,current_var)
!!$                    enddo
!!$                 enddo
!!$              enddo
!!$           enddo

        endif
     enddo
  endif


     if (multigrid_on.ge.1) then
        call pf_guardcell(pid,iopt,level)
     else

!!$        do i = 1,g_max_lvl
!!$           call pf_guardcell(pid,iopt,i)
!!$        enddo
 !       write(*,*)"before guardcell",nguard
!        call amr_guardcell(pid,iopt,nguard)
 !       write(*,*)"after gaurdcell"
     endif
 
 
!  current_work=1


end subroutine app_mg_smooth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pf_mg_smooth_var(lb,level,new_soln,pid)
  ! Smoother for node centered  variables
  ! lb is the block number
  ! n is the variable number being smoothed
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  integer,intent(in) :: lb,level,pid
  integer :: i,j,k,aa,bb,ii,jj,kk
  double precision, intent(out) :: new_soln(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,N_unk)
  double precision :: dx,Fi(N_unk,nxb,nxb,nxb),D_Fi(N_unk,nxb,nxb,nxb),Wk(N_unk,nxb,nxb,nxb)
  integer:: ifid
  logical :: mode_1_log
  logical :: Not_number,smooth_flag=.true.



  if (mode_1 == 0) then
     mode_1_log = .false.
  else
     mode_1_log = .true.
  end if
  dx = bsize(1,lb)/real(nxb)
  k=1
  d_Fi=0d0
  !Form Jacobi matrix
  


  call Get_Jacobi(lb, dx,Fi,d_Fi,smooth_flag)
  if(smooth_flag)then

     do k=1,nxb
        do j=1,nxb
           do i=1,nxb
           do aa=1,N_unk
              Fi(aa,i,j,k)=Fi(aa,i,j,k)-unk(2+(aa-1)*nunkvbles,i+1,j+1,k+1,lb) 
           enddo
        enddo
     enddo
  enddo
  if(g_diagonal)then
     call TLinearDiag(d_Fi,Fi,Wk)
  else
     write(*,*)"pf_mg_smooth_var. Not set up for block solve"
     stop
     call TLinear(d_Fi,Fi,Wk)
  endif
  do ii=1,N_unk
     do i=1,nxb
        do j=1,nxb
           do k=1,nxb
              new_soln(i+1,j+1,k+1,ii)=unk(1+(ii-1)*nunkvbles,i+1,j+1,k+1,lb)-Wk(ii,i,j,k)
           enddo
        enddo
     enddo
  enddo
endif
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
  integer :: i,j,k,ii,n
  double precision :: dx,Fi(N_unk)
  double precision :: rfactor, rf1, rf2, rf3
  logical :: mode_1_log,Not_number,iflag=.false.







  dx = bsize(1,lb)/real(nxb)


  if(current_var.eq.1)then
     do k=2,nxb+1
        do j=2,nxb+1
           do i=2,nxb+1
              
              call get_RHS_Fi_2d(i,j,k,lb,dx,Fi)
              
              do ii=1,N_unk
                 n=1+(ii-1)*nunkvbles
                 if(Fi(ii).le.0d0.or.Fi(ii).ge.0d0)then
                 else
                    write(*,*)i,j,ii,"in app_mg_get_defect_var",Fi
                    stop
                 endif
                 !d = f - A(u)
                 unk(n+5,i,j,k,lb) = unk(n+1,i,j,k,lb)  -Fi(ii)
                 Fi(ii)=unk(n+5,i,j,k,lb) 
                 Fi(ii)=abs(Fi(ii)) !apply weighting
                 rms=rms+Fi(ii)*Fi(ii)
                 npts = npts + 1
                 if(Fi(ii).gt.max_defect)then
                    max_defect=Fi(ii) !debug
                 end if
                 
              enddo
           end do
        enddo
     end do
  endif

end subroutine app_mg_get_defect_var


subroutine Get_Jacobi(lb, dx,Fi,d_Fi,smooth_flag)
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
 
integer, intent(in)::lb
double precision, intent(in)::dx
double precision, intent(out) :: Fi(N_unk,nxb,nxb,nxb),d_Fi(N_unk,nxb,nxb,nxb)
logical, intent(inout) :: smooth_flag
double precision :: hh=1d-10
integer ip,ifid,n,ii,jj,aa,kk,ll
double precision :: Fid(N_unk,nxb,nxb,nxb), unkM(M_unk,10,10,10)
double precision :: c,rfactor,rf4,rf5,rf6,vbles10(N_unk,10,10,10),vbles10d(N_unk,10,10,10),v_sum=0d0

  v_sum = 0d0
  do ii=1,g_order
     do jj=1,g_order
        do kk=1,g_order
           do ip=1,N_unk
              vbles10(ip,ii,jj,kk)=unk(1+(ip-1)*nunkvbles,ii,jj,kk,lb)
              v_sum = v_sum + abs(vbles10(ip,ii,jj,kk))
           enddo
        enddo
     enddo
  enddo

  if(v_sum.lt.1d-16)then
     smooth_flag = .false.
     return
  endif




! Call RHS, which includes time derivative
  Fi = 0d0
  do ii=1,nxb
     do jj=1,nxb
        do kk=1,nxb
           call get_RHS_Fi(dx,Fi(:,ii,jj,kk),mode_1,nunkvbles,dt,dtold,lnblocks,vbles10,ii+1,jj+1,kk+1,lb)
        enddo
     enddo
  enddo
  

  vbles10d = vbles10
!  if(g_diagonal)then
  D_fi=0d0
  do kk=1,nxb
     do jj=1,nxb
        do ii=1,nxb
           do ifid = 1,N_unk
              hh=abs(hh)
              if(vbles10(ifid,nxb/2,nxb/2,nxb/2).gt.0.5)then
                 hh=-abs(hh)
              endif
              vbles10d(ifid,ii+1,jj+1,kk+1)=vbles10(ifid,ii+1,jj+1,kk+1)+hh
              call get_RHS_Fi(dx,Fid(:,ii,jj,kk),mode_1,nunkvbles,dt,dtold,lnblocks,vbles10d,ii+1,jj+1,kk+1,lb)
              D_Fi(ifid,ii,jj,kk)=(Fid(ifid,ii,jj,kk)-Fi(ifid,ii,jj,kk))/hh
              vbles10d(ifid,ii+1,jj+1,kk+1)=vbles10(ifid,ii+1,jj+1,kk+1)
           enddo
        enddo
     enddo
  enddo


end subroutine Get_Jacobi




!*********************************************************
!*   SOLVE A LINEAR SYSTEM BY TRIANGULARIZATION METHOD   *
!*********************************************************
subroutine TlinearDiag(AA,BB,X)
use solution_parameters
double precision,intent(inout) :: AA(NN_unk),BB(NN_unk)
double precision,intent(out) :: X(NN_unk)
double precision ::  A(NN_unk,NN_unk),B(NN_unk)
INTEGER :: N,i 
N=NN_unk
do i = 1,N
   X(i) = BB(i)/AA(i)
enddo

end subroutine TlinearDiag

subroutine Tlinear(AA,BB,X)
use solution_parameters
double precision,intent(inout) :: AA(NN_unk,NN_unk),BB(NN_unk)
double precision,intent(out) :: X(NN_unk)
double precision ::  A(NN_unk,NN_unk),B(NN_unk)
INTEGER :: N 
N=NN_unk 
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

END subroutine Tlinear

subroutine Triangular(AA,BB,X)
use solution_parameters
double precision,intent(inout) :: AA(NN_unk,NN_unk),BB(NN_unk)
double precision,intent(out) :: X(NN_unk)
double precision ::  A(NN_unk,NN_unk),B(NN_unk)
INTEGER :: N 
N=NN_unk 
do i=1,N
B(i)=BB(i)
do j=1,N
A(i,j)=AA(i,j)
enddo
enddo

! Solve triangular system
X(N) = B(N) / A(N, N)
DO I = N - 1, 1, -1
  S = 0.d0
  DO K = I + 1, N
    S = S + A(I, K) * X(K)
  END DO
  X(I) = (B(I) - S) / A(I, I)
END DO

END subroutine Triangular
