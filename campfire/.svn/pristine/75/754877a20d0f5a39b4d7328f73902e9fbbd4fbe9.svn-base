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
  integer :: ii,jj,ip
  logical Not_number
  double precision:: new_soln(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,3)
  double precision:: vbles(4,3,3),phaseRHS_c,c,phi,dx


  iopt=1


  


  ! Cycle over all blocks

     if(level.gt.-1)then
        do lb=block_starts(level), block_starts(level+1)-1
           call pf_mg_smooth_var(lb,level,new_soln)
           do current_var=1,total_vars
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
                 do i=il_bnd+nguard,iu_bnd-nguard
                    if(.not.Not_number(new_soln(i,j,1,current_var))) unk(1+(current_var-1)*nunkvbles,i,j,1,lb)=new_soln(i,j,1,current_var)
                 enddo
              enddo
           enddo
        end do
     else

        do lb=1, lnblocks
           current_unk=1+(current_var-1)*nunkvbles
           if (.not.has_children(lb))then
              call pf_mg_smooth_var(lb,level,new_soln)
              do current_var=1,total_vars
                 do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
                    do i=il_bnd+nguard,iu_bnd-nguard
                       if(.not.Not_number(new_soln(i,j,1,current_var)))unk(1+(current_var-1)*nunkvbles,i,j,1,lb)=new_soln(i,j,1,current_var)
                    enddo
                 enddo
              enddo
           endif
        enddo
     endif


  current_work=1
  do current_var=1,total_vars
     ! set the start point for the real variable and it's workspace data
     current_unk=1+(current_var-1)*nunkvbles
     ! Update gcells for "current_unk" only
     gcell_on_cc(:)=.false.
     gcell_on_cc(current_unk)=.true.
     call pf_guardcell(pid,iopt,level)
     gcell_on_cc(:)=.true.
  end do






end subroutine app_mg_smooth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pf_mg_smooth_var(lb,level,new_soln)
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
  integer,intent(in) :: lb,level
  integer :: i,j,k,ii,jj,kk
  double precision, intent(out) :: new_soln(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,3)
  real :: dx,Fi(nvar),D_Fi(3,3),sum_new_soln,FFi(3),Wk(3),Wj(3),W_max,diag(3)
  double precision :: rfactor, rf1, rf2, rf3
  double precision :: identity(3,3)=0d0,fid(nvar)
  integer:: iflg=0,n,i2,j2,k2,ifid
  logical :: mode_1_log
  logical Not_number
  
  if (mode_1 == 0) then
  mode_1_log = .false.
  else
  mode_1_log = .true.
  end if



  do ii=1,3
   identity(ii,ii)=1d0
  enddo  
  d_Fi=0.
  rfactor = dt/dtold
  rf1 = (rfactor+1.0)*dt/(2.0*rfactor+1.0)
 
  rf2 = (rfactor+1.0)*(rfactor+1.0)/(2.0*rfactor+1.0)
  rf3 = rfactor*rfactor/(2.0*rfactor+1.0)
  dx = bsize(1,lb)/real(nxb)

  
  k=1



  do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
     do i=il_bnd+nguard,iu_bnd-nguard

        call get_RHS_Fi_2d(i, j, k, lb, dx, Fi)


        if(time.gt.1d0)then
        do ifid=1,total_vars

           call GET_RHS_FI_2D_D(i, j, k, lb, dx,  fi, fid, ifid)

           do ii=1,total_vars
              D_Fi(ii,ifid)=fid(ii)
           enddo
        enddo
        endif


        do ii=1,total_vars
           n=(ii-1)*nunkvbles
           if(mode_1.le.1)then
              FFi(ii)=unk(1+n,i,j,k,lb)-dt*Fi(ii)-unk(2+n,i,j,k,lb)
              do jj=1,total_vars
                 d_Fi(ii,jj)=identity(ii,jj)-dt*D_Fi(ii,jj)
              enddo
           else if(mode_1.eq.2)then
              FFi(ii)=unk(1+n,i,j,k,lb)-rf1*Fi(ii)-unk(2+n,i,j,k,lb)
              do jj=1,total_vars
                 d_Fi(ii,jj)=identity(ii,jj)-rf1*D_Fi(ii,jj)
              enddo
           end if
        enddo

        call Tlinear(d_Fi,FFi,Wk)
        do ii=1,total_vars
           new_soln(i,j,k,ii)=unk(1+(ii-1)*nunkvbles,i,j,k,lb)-weight*Wk(ii)
!           new_soln(i,j,k,ii)=unk(1+(ii-1)*nunkvbles,i,j,k,lb)-FFi(ii)/d_Fi(ii,ii)
        enddo
     end do
  end do


end subroutine pf_mg_smooth_var


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine app_mg_get_defect_var(lb, max_defect, rms, npts)
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
  integer :: i,j,k,ii,jj,n,kk,ip
  double precision :: dx,Fi(3), dx2inv,x0,y0,fexp
  double precision :: rfactor, rf1, rf2, rf3
  double precision :: LapMult,tau_val
  double precision :: c,T,phi,c_dot,phi_dot,G4,S_f,S_p,S_defect
  double precision :: wt(2)=(/1e0,1e0/)
  double precision :: vbles(3,3,3),phaseRHS_c
  logical mode_1_log,Not_number



  rfactor = dt/dtold
  rf1 = (rfactor+1.0)*dt/(2.0*rfactor+1.0)
  rf2 = (rfactor+1.0)*(rfactor+1.0)/(2.0*rfactor+1.0)
  rf3 = rfactor*rfactor/(2.0*rfactor+1.0)
  
  dx = bsize(1,lb)/real(nxb)
  dx2inv = 1.0/(dx*dx)





  do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
     do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
        do i=il_bnd+nguard,iu_bnd-nguard
         if(current_var.eq.1)then
           call get_RHS_Fi_2d(i,j,k,lb,dx,Fi)

           do ii=1,total_vars
            n=(ii-1)*nunkvbles
            if(mode_1.le.1)then
              Fi(ii)=unk(n+1,i,j,k,lb)-dt*Fi(ii)-unk(n+2,i,j,k,lb)
            else if(mode_1.eq.2)then
              Fi(ii)=unk(n+1,i,j,k,lb)-rf1*Fi(ii)-unk(n+2,i,j,k,lb)
            end if
            unk(n+nunkvbles,i,j,k,lb)=-Fi(ii)
            Fi(ii)=abs(Fi(ii)*wt(ii)) !apply weighting
            npts = npts + 1
            if(Fi(ii).gt.max_defect)then
              max_defect=Fi(ii) 
            end if
           enddo
         endif !current_var=1
         work(i,j,k,lb,2) = -work(i,j,k,lb,2)+ unk(current_var*nunkvbles,i,j,k,lb)
        end do
     end do
  end do

  return
end subroutine app_mg_get_defect_var


subroutine GET_RHS_FI_2D_D(i, j, k, lb, dx, fi,  fid, ifid)
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
double precision, intent(out):: Fid(3)
double precision, intent(out):: Fi(3)
double precision vbles(3,3,3),unk_(18),phid,cd,psid,c,phi,psi,PhaseRHS_c
double precision :: vblesd(3,3,3)=0.,rf4,rf5,rf6,rfactor,x
integer ip,n,ii,jj




  Fid=0d0
  phid=0d0
  cd=0d0
  psid=0d0
  if(ifid.eq.1)phid=1d0
  if(ifid.eq.2)cd=1d0
  if(ifid.eq.3)psid=1d0
  
  do ii=1,3
     do jj=1,3
        do ip=1,total_vars
           vbles(ip,ii,jj)=unk(1+(ip-1)*nunkvbles,i+ii-2,j+jj-2,k,lb)
        enddo
     enddo
  enddo




  phi=vbles(1,2,2)
  c=vbles(2,2,2)
  psi=vbles(3,2,2)
  
  do ii=1,nvar
     unk_(ii) = unk(ii,i,j,k,lb)
  enddo
  x = bnd_box(1,1,lb)+(i-1.5)*dx

  call GET_RHS_FI_D(vbles, dx, fi, fid, mode_1, nunkvbles, dt, dtold&
       & , lnblocks, unk_, phi, phid, c, cd, psi, psid, time, x)


  
  do ip=1,total_vars
     unk(1+(ip-1)*nunkvbles,i,j,k,lb)=vbles(ip,2,2)
  enddo
end subroutine get_RHS_Fi_2d_D
