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
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid,level
  integer lb, iopt,i,j,k
  integer :: ii,jj
  logical Not_number
  double precision:: new_soln(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,N_unk),phi(N_phases),dx
  iopt=1

  ! Cycle over all blocks

  do lb=block_starts(level), block_starts(level+1)-1
     dx = bsize(1,lb)/real(nxb)
     call pf_mg_smooth_var(lb,level,new_soln,pid)
     do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
        do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
           do i=il_bnd+nguard,iu_bnd-nguard
              new_soln(i,j,k,1) = max(0d0,min(1d0,new_soln(i,j,k,1)))
              new_soln(i,j,k,3) = max(0d0,min(1d0,new_soln(i,j,k,3)))
              do current_var=1,N_unk
                 if(.not.Not_number(new_soln(i,j,k,current_var)))unk(1+(current_var-1)*nunkvbles,i,j,k,lb)=new_soln(i,j,k,current_var)
              enddo
           enddo
        enddo
     enddo
  end do
  

  


  
  current_work=1
  
  gcell_on_cc(:)=.true.
  
  call pf_guardcell(pid,iopt,level)

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
  integer :: i,j,k,ii,aa
  double precision, intent(out) :: new_soln(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,N_unk)
  double precision :: dx,Fi(NN_unk),D_Fi(NN_unk,NN_unk),Wk(NN_unk)
  integer:: ifid
  logical :: mode_1_log
  logical Not_number

  if (mode_1 == 0) then
  mode_1_log = .false.
  else
  mode_1_log = .true.
  end if





  dx = bsize(1,lb)/real(nxb)

  d_Fi=0d0
  do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d-1,2 
     do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d-1,2 
        do i=il_bnd+nguard,iu_bnd-nguard-1,2
           
           
           
           !Form Jacobi matrix        
           do ifid=1,N_unk
              call GET_RHS_FI_2D_D(i, j, k, lb, dx,ifid,Fi,d_Fi)
           enddo
           
           
           
           call Tlinear(d_Fi,Fi,Wk)
           
           
           do ii=1,N_unk
              do aa=1,4
                 new_soln(i+rrarray3d(aa,1),j+rrarray3d(aa,2),k+rrarray3d(aa,3),ii)=unk(1+(ii-1)*nunkvbles,i+rrarray3d(aa,1),j+rrarray3d(aa,2),k+rrarray3d(aa,3),lb)-weight*Wk(ii+(aa-1)*N_unk)
              enddo
           enddo
        end do
     end do
  enddo
  
  

  
     
  
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
  double precision :: wt(N_unk)
  logical mode_1_log,Not_number




  wt =1d0


  rfactor = dt/dtold
  rf1 = (rfactor+1.0)*dt/(2.0*rfactor+1.0)
  rf2 = (rfactor+1.0)*(rfactor+1.0)/(2.0*rfactor+1.0)
  rf3 = rfactor*rfactor/(2.0*rfactor+1.0)
  dx = bsize(1,lb)/real(nxb)


  do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
     do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
        do i=il_bnd+nguard,iu_bnd-nguard
         if(current_var.eq.1)then
           call get_RHS_Fi_2d(i,j,k,lb,dx,Fi)
           do ii=1,N_unk
            n=(ii-1)*nunkvbles
            unk(n+nunkvbles,i,j,k,lb)=-Fi(ii)
            Fi(ii)=abs(Fi(ii)*wt(ii)) !apply weighting
            npts = npts + 1
            if(Fi(ii).gt.max_defect)then
              max_defect=Fi(ii) !debug
            end if
           enddo
         endif !current_var=1
         work(i,j,k,lb,2) = -work(i,j,k,lb,2)+ unk(current_var*nunkvbles,i,j,k,lb)
        end do
     end do
  end do
  return
end subroutine app_mg_get_defect_var


subroutine GET_RHS_FI_2D_D(i, j, k, lb, dx,ifid,Fi,d_Fi)
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
 
integer, intent(in)::i,j,k,lb,ifid
double precision, intent(in)::dx
double precision, intent(inout) :: Fi(NN_unk),d_Fi(NN_unk,NN_unk)
double precision :: unk_(M_unk),vblesd(N_unk+1,3,3,3),hh=1d-10
integer ip,n,ii,jj,kk,aa
double precision :: FiM(N_unk,4),FidM(N_unk,4), unkM(M_unk,4)
double precision :: vblesM(N_unk+1,3,3,3,4),c,rfactor,rf4,rf5,rf6,MolPerVol



do ii=1,3
   do jj=1,3
      do kk=1,3
         do ip=1,N_unk
            do aa = 1,4
               vblesM(ip,ii,jj,kk,aa)=unk(1+(ip-1)*nunkvbles,i+ii-2+rrarray3d(aa,1),j+jj-2+rrarray3d(aa,2),k+kk-2+rrarray3d(aa,3),lb)
            enddo
         enddo
      enddo
   enddo
enddo

  ip=N_unk+1
  do ii=1,3
     do jj=1,3
        do kk=1,3
           do aa = 1,4
              if(mode_1.eq.1) then
                 vblesM(ip,ii,jj,kk,aa)=(unk(1,i+ii-2+rrarray3d(aa,1),j+jj-2+rrarray3d(aa,2),k+kk-2+rrarray3d(aa,3),lb)-&
                      unk(2,i+ii-2+rrarray3d(aa,1),j+jj-2+rrarray3d(aa,2),k+kk-2+rrarray3d(aa,3),lb))/dt
              else
                 rfactor = dt/dtold
                 rf4 = rfactor*rfactor/(rfactor+1.0)
                 rf5 = (rfactor+1.0)
                 rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
                 vblesM(ip,ii,jj,kk,aa)=(rf4*unk(5,i+ii-2+rrarray3d(aa,1),j+jj-2+rrarray3d(aa,2),k+kk-2+rrarray3d(aa,3),lb)-&
                      rf5*unk(4,i+ii-2+rrarray3d(aa,1),j+jj-2+rrarray3d(aa,2),k+kk-2+rrarray3d(aa,3),lb)+&
                      rf6*unk(1,i+ii-2+rrarray3d(aa,1),j+jj-2+rrarray3d(aa,2),k+kk-2+rrarray3d(aa,3),lb))/dt
              endif
              c=vblesM(N_unk,ii,jj,kk,aa)
              vblesM(ip,ii,jj,kk,aa)=-vblesM(ip,ii,jj,kk,aa)*g_W(1)/((1.-c)+c*g_M(2)/g_M(1))
              MolPerVol=g_MolPerVol(1)*(1.-c)+g_MolPerVol(2)*c
              vblesM(ip,ii,jj,kk,aa)=vblesM(ip,ii,jj,kk,aa)/MolPerVol*g_lambda**2
           enddo
        enddo
     enddo
  enddo
  ! boundary condition. If Jacobi geaometry changes these lines must change too  
!  if(neigh(1,2,lb).eq.bc_cond_sym.and. i.eq.9)then
!     ii = 3
!     do jj=1,3
!        ip=N_phases+1
!        vblesM(ip,ii,jj,4)=delta  ! 
!        vblesM(ip,ii,jj,2)=delta  !
!     enddo
!  endif



  do aa=1,4
     do ii=1,M_unk
        unkM(ii,aa)= unk(ii,i+rrarray3d(aa,1),j+rrarray3d(aa,2),k+kk-2+rrarray3d(aa,3),lb)
     enddo
enddo

! Call RHS, which includes time derivative
  do aa=1,4
     call get_RHS_Fi(vblesM(:,:,:,:,aa),dx,FiM(:,aa),mode_1,nunkvbles,dt,dtold,lnblocks,unkM(:,aa))
  enddo

!Numerical Jacobi (2 x 2 x 5)^2
!set up correct sign of "h" for numerical differentiation



  hh=abs(hh)
  if(ifid.eq.N_phases+1)then
     continue
  elseif(vblesM(ifid,2,2,2,1).gt.0.5)then
     hh=-abs(hh)
  endif

!all this is 2D of course, so see the 2D code...
! Compute Jacobi, 12 out of 16 possible blocks corresponding to a 2 x 2 area within each block (of normally 8 x 8)
! the 1,1 and 2,2 and 3,3 and 4,4 are the usual diagonal terms for the pointwise Newton.
! Matrix rrblock(12,6) contains the relations between different elements in the 2 x 2 
!--------------------------------------------
!          |          |          |          |
!    1     |    5     |     7    |          |
!          |          |          |          |
!--------------------------------------------
!          |          |          |          |
!    6     |    2     |          |     9    |
!          |          |          |          |
!--------------------------------------------
!          |          |          |          |
!    8     |          |    3     |     11   |
!          |          |          |          |
!--------------------------------------------
!          |          |          |          |
!          |    10    |    12    |     4    |
!          |          |          |          |
!-------------------------------------------

! Each number in the above contain a 5 x 5 sub block for each of the 5 variables
! blocks 1 to 4 are the usual pointwise elements at location: (i,j); (i+1,j);(i,j+1) and (i+1,j+1)
! block 5 and 6 is the coupling between (i,j) and (i+1,j)
! block 7 and 8 is the coupling between (i,j) and (i,j+1)
! block 9 and 10 is the coupling between (i+1,j) and (i+1,j+1)
! block 11 and 12 is the coupling between (i,j+1) and (i+1,j+1)
! The missing sub blocks do not contribute (much) being diagonal elements in the mesh
!
! All this can be changed by changing rrblock, in mp_problem.f90  to allow, for example,  line Jacobi
! of 4 elements in a line (i,j);(i+1,j);(i+2,j);(i+3,j). The difficulty is getting the elements of rrblock() correct!
!
! the first 2 lines of rrblock: rrblock(:,1) and rrblock(:,2) contains position in 20 x 20 Jacobi matrix
! the next 2 lines rrblock(j,3)  refers to position in 2 x 2 elements, e.g. rrblock(1,3)=3 is (i,j+1)
! The final 2 lines rrblock(j,4) and rrblock(j,5) also refer to position (i,j) where (rrblock(7,4)=2,rrblock(7,5)=3) is (i,j+1)
! thus addressing the correct element in the 3 x 3 stencil, (2,2) being the central node.




  do jj=1,32
     vblesd=vblesM(:,:,:,:,rrblock3d(jj,3))
     vblesd(ifid,rrblock3d(jj,4),rrblock3d(jj,5),rrblock3d(jj,6))=vblesM(ifid,rrblock3d(jj,4),rrblock3d(jj,5),rrblock3d(jj,6),rrblock3d(jj,3))+hh
     
     
     call get_RHS_Fi(vblesd,dx,FidM(:,rrblock3d(jj,3)),mode_1,nunkvbles,dt,dtold,lnblocks,unkM(:,rrblock3d(jj,3)))
     
     do ip=1,total_vars
        D_Fi(ip+rrblock3d(jj,1),ifid+rrblock3d(jj,2))=(fidM(ip,rrblock3d(jj,3))-FiM(ip,rrblock3d(jj,3)))/hh
        Fi(ip+rrblock3d(jj,1)) = FiM(ip,rrblock3d(jj,3))
     enddo
  enddo




 






end subroutine get_RHS_Fi_2d_D


!*********************************************************
!*   SOLVE A LINEAR SYSTEM BY TRIANGULARIZATION METHOD   *
!*********************************************************
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
