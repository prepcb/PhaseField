!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_output_perstep                          REQUIRED
!!!!  * Post V-cycles anything to be done every step
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_output_occasional                       REQUIRED
!!!!  * Post V-cycles anything to be done every 'output_frequency' steps
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  phase_field_check
!!!!  * Main routine for obtaining tip location and radius
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  phase_loc_rad
!!!!  * Computes tip location and radius using quadratic interpolation
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  get_quadratic_interp
!!!!  * Does quadratic interpolation of supplied data
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  pf_output_2d
!!!!  * Prints points showing 2-d phase field boundary
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Latest version: Chris Goodyer, May 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Generic routines for output

subroutine app_output_perstep(pid, noprocs, inumber)
  implicit none
  integer, intent(in) :: inumber,pid, noprocs

  call phase_field_check(pid, inumber)

  return
end subroutine app_output_perstep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_output_occasional(pid, noprocs, inumber)
  use multigrid_parameters
  use solution_parameters
  use paramesh_dimensions
  use time_dep_parameters
  use paramesh_interfaces
  use checkpoint_parameters
  implicit none
  integer, intent(in) :: inumber, pid, noprocs
  logical :: chombo = .false.
  logical chk_gd
  integer inumber_100
  inumber_100=inumber
  if(g_output_rate.ge.100)inumber_100=inumber/100

  ! 2-d points around front
  if (ndim.eq.2) call pf_output_2d(pid, inumber_100)

  ! Chombo output
  if (chombo)then
     if(pid.eq.0)then
        print *,"Output Chombo file"
     end if
     ! Write file to disk
     call amr_plotfile_chombo(inumber_100)
     write(61,*) 'Step, time', inumber_100, time
     flush(61) 
  end if

  ! Checkpointing output
  if(pid.eq.0) print *,"Checkpoint!"
  chk_gd = .FALSE.
  chk_checkf = 'hdf5'
  call amr_checkpoint_wr(inumber_100, user_attr_1=dt, user_attr_2=time, user_attr_3=dtold)
  if (pid.eq.0) then
     open (unit=99,file="CHK.out")
     write(99,*) inumber_100
     close(99)
  endif



  return
end subroutine app_output_occasional

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine phase_field_check(pid, inumber_100)
  ! Not really a user edit but since it's specific to the problem
  ! being simulated it's sat here.
  use paramesh_dimensions
  use physicaldata
  use tree
  use paramesh_comm_data ! Needed for definition of amr_mpi_real
  use paramesh_interfaces
  use solution_parameters
  use time_dep_parameters
  use multigrid_parameters
  use refinement_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: inumber_100,pid
  integer :: ierr,lb,i,j,k,mid_x,mid_y,mid_z, found
  integer :: ii,jj,n_leaf
  double precision :: phase,umin,umax,phimin,phimax,u,phi(n_phases),phi_liq,c,T,temp,gradphi,gradphi2,theta1,theta2
  double precision :: dx,dist,phi_val_1,phi_val_2,phi_val_3,phi_val_4,solute_max,solute_min,s1,s2,s3,s4,s5,s6,s7,s8
  double precision :: new_loc, new_rad,new_DW, tempr,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,sol_max,sol_min,heat(4)
  double precision :: xl, xm, xr, phil, phim, phir, thetal, thetam, thetar, sigmastartemp,heats(3)
  logical :: on_axisx=.false.,on_axisy=.false.
  double precision :: FreeEnergy,TotalEnergy0,TotalEnergy,vbles(N_unk+1,3,3),GradientEnergy,MolPerVol, potential
  double precision :: sum_c,sum_c0,sum_liq,sum_liq0
!  character (len=19) :: filename
  !print *,"Guardcell axis node count:",axis_node_count
  !axis_node_count=0
  
  vbles=0d0


  phase=0.0
  umin=1e30
  umax=-1e30
  phimin=1e30
  phimax=-1e30
  mid_x=il_bnd+nxb/2+nguard
  mid_y=jl_bnd+nyb/2+nguard
  if (ndim.eq.3) mid_z=kl_bnd+nzb/2+nguard
  found=0
  if(pid.eq.0)print *,"Phase field summary for time step: ",inumber_100
! CEG removed profiles
!  if(pid.gt.0)print *,"Warning, solute profiler not yet coded for parallelism"
!  if(mod(inumber_100,100).eq.0)then
!     write (filename,1001) "solute_profile_",inumber_100
!     if(pid.eq.0)open (unit=136,file=filename)
!     write (filename,1001) "phase_profile_",inumber_100
!     if(pid.eq.0)open (unit=137,file=filename)
!     write (filename,1001) "thermal_profile_",inumber_100
!     if(pid.eq.0)open (unit=138,file=filename)
!  end if
  

  dist = 0.0
  sigmastartemp = 0.0
  sol_max=-1d9
  sol_min= 1d9
  s1=0d0
  s2=0d0
  s3=0d0
  s4=0d0
  s5=0d0
  s6=0d0
  s7=0d0
  s8=0d0
  TotalEnergy=0d0
!   if (nvar.eq.2.and.solute) found=1
  n_leaf=0
  do lb=1,lnblocks
!!$     if(.not.has_children(lb)) then
!!$        do i=il_bnd+1,iu_bnd-1
!!$           do j=jl_bnd+1,ju_bnd-1
!!$              phi_liq = unk(1,i,j,1,lb)
!!$              c   = unk(1+(N_unk-1)*nunkvbles,i,j,1,lb)
!!$              s5 = s5 + c*phi_liq**2*(1d0 - phi_liq)**2
!!$              s6 = s6 + phi_liq**2*(1d0 - phi_liq)**2
!!$           enddo
!!$        enddo
!!$     endif
!     if(nodetype(lb).eq.1) then
        dx = bsize(1,lb)/real(nxb)


        ! First find out if we're "on axis"
        ! how do we know what the "axis" is?
        ! Simple the nucleate_(x,y,z) parameter holds the co-ords
        ! To be on an axis 2 out of 3 of the block faces must match these co-ords
        ! Only testing lower co-ords so for corner nucleation the nucleus MUST BE at min co-ords for each dimension
        on_axisx=.false.
        on_axisy=.false.
        
        if(bnd_box(1,2,lb).lt.1d-15)then
           ! lower y bnd matches
           if(ndim.eq.2.or.bnd_box(1,3,lb).lt.1d-15)then
              ! lower z bnd matches
              ! y and z co-ords match so we're on the x axis
              on_axisx=.true.
              !axis_node_count=axis_node_count+1
           end if
        end if

        if(bnd_box(1,1,lb).lt.1d-15)then
           ! lower y bnd matches
           if(ndim.eq.2.or.bnd_box(1,3,lb).lt.1d-15)then
              ! lower z bnd matches
              ! y and z co-ords match so we're on the x axis
              on_axisy=.true.
              !axis_node_count=axis_node_count+1
           end if
        end if
!        j=jl_bnd+nguard
!        k=kl_bnd+nguard*k3d
        if(on_axisy)then
           if(dx.lt.32.0)call phase_loc_rad_y(lb, s5,s6,s7,s8)
        endif
        if(on_axisx)then
           if(dx.lt.32.0)call phase_loc_rad(lb, s1,s2,s3,s4)
        end if

!     end if
!     if (found.eq.2) exit
     
     if(.not.has_children(lb).and.bnd_box(1,2,lb).lt.1d-15) then
        do i=il_bnd+1,iu_bnd-1
           if(unk(1,i,2,1,lb).gt.1d-3.and. unk(1,i,2,1,lb).lt.1d0-1d-3)then
              sol_max=max(unk(1+(N_unk-1)*nunkvbles,i,2,1,lb),sol_max)
              sol_min=min(unk(1+(N_unk-1)*nunkvbles,i,2,1,lb),sol_min)
           endif
        enddo
     endif
     if(thermal)then
        if(.not.has_children(lb)) then
           do j=jl_bnd+1,ju_bnd-1
              do i=il_bnd+1,iu_bnd-1
                 sum_heat(1)=sum_heat(1) + unk(1+N_unk,i,j,1,lb)*dx*dx
                 sum_heat(2)=sum_heat(2) + unk(2+N_unk,i,j,1,lb)*dx*dx+ unk(3+18,i,j,1,lb)*dx*dx
                 sum_heat(3)=sum_heat(4) + unk(4+N_unk,i,j,1,lb)*dx*dx
              end do
           end do
        endif
     else
        sum_heat=0d0
     endif




     if(.not.has_children(lb)) then
     n_leaf=n_leaf+1
        do j=jl_bnd+1,ju_bnd-1
           do i=il_bnd+1,iu_bnd-1
              phi(1) = unk(1,i,j,1,lb)
              phi(2) = unk(1+nunkvbles,i,j,1,lb)
              phi(3) = unk(1+2*nunkvbles,i,j,1,lb)
              do k=1,3
                 vbles(k,2,2)=unk(1+(k-1)*nunkvbles,i,j,1,lb)
                 vbles(k,1,2)=unk(1+(k-1)*nunkvbles,i-1,j,1,lb)
                 vbles(k,2,1)=unk(1+(k-1)*nunkvbles,i,j-1,1,lb)
                 vbles(k,3,2)=unk(1+(k-1)*nunkvbles,i+1,j,1,lb)
                 vbles(k,2,3)=unk(1+(k-1)*nunkvbles,i,j+1,1,lb)
                 vbles(k,1,1)=unk(1+(k-1)*nunkvbles,i-1,j-1,1,lb)
                 vbles(k,1,3)=unk(1+(k-1)*nunkvbles,i-1,j+1,1,lb)
                 vbles(k,3,1)=unk(1+(k-1)*nunkvbles,i+1,j-1,1,lb)
                 vbles(k,3,3)=unk(1+(k-1)*nunkvbles,i+1,j+1,1,lb)
              enddo
              T      = unk(1+3*nunkvbles,i,j,1,lb)
              c      = unk(1+4*nunkvbles,i,j,1,lb)
              MolPerVol=(1.-c)*g_MolPerVol(1)+c*g_MolPerVol(2)
              TotalEnergy = TotalEnergy + (g_lambda**2*GradientEnergy(vbles,0,dx)&
                   +FreeEnergy(c,T,phi,0d0,phi,0)*MolPerVol/g_W(1)&
                   +potential(Phi,c,0))*dx*dx


           enddo
        enddo
     endif


     if(.not.has_children(lb)) then

        sum_c = 0d0
        sum_liq = 0d0
        do j=jl_bnd+1,ju_bnd-1
           do i=il_bnd+1,iu_bnd-1
              if(dx.eq.g_min_dx)then
                 c      = unk(1+(N_unk-1)*nunkvbles,i,j,1,lb)
                 phi(1) = unk(1,i,j,1,lb)
                 sum_c = sum_c + c*phi(1)*dx*dx
                 sum_liq = sum_liq + phi(1)*dx*dx
              endif
           enddo
        enddo
     endif
  end do

  call MPI_REDUCE(sum_c,sum_c0,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(sum_liq,sum_liq0,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  call MPI_REDUCE(TotalEnergy,TotalEnergy0,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  !s1=s1+dw_phi            !double well
  !s2=s2+dw_phi*phi_x      !phi'(x) *dw
  !s3=s3+dw_phi*x0         !x *dw
  !s4=s4+dw_phi*phi_yy     !phi''(y)*dw
  !s5=s5 + c*phi_liq**2*(1d0 - phi_liq)**2 !get concentration at tip
  !s6=s6 + phi_liq**2*(1d0 - phi_liq)**2   ! global DW different to s1




  
  call MPI_REDUCE(sol_max,solute_max,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(sol_min,solute_min,1,amr_mpi_real,MPI_MIN,0,MPI_COMM_WORLD,ierr)

  call MPI_REDUCE(sum_heat(1),heats(1),1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(sum_heat(2),heats(2),1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(sum_heat(3),heats(3),1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!  print *,"Phase_field_check axis node count:",axis_node_count

  call MPI_REDUCE(s1,temp1,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(s2,temp2,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(s3,temp3,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(s4,temp4,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(s6,temp6,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(s8,temp8,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  g_nuc_radius = temp6/temp8

  call MPI_Bcast(g_nuc_radius, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  if(g_1D)temp3=0d0
  sum_heat=0d0
!  if(temp3/temp1.lt.nuc_radius*0.5)then
!     write(*,*)"melting: tip< nuc_radius/2"
!  endif


!  call MPI_REDUCE(temp_min,temp4,1,amr_mpi_real,MPI_MIN,0,MPI_COMM_WORLD,ierr)
  if(pid.eq.0)write (6,'(A,A)') 27, "[36;1;47m"

  !s1=s1+dw_phi            !double well
  !s2=s2+dw_phi*phi_x      !phi'(x) *dw
  !s3=s3+dw_phi*x0         !x *dw
  !s4=s4+dw_phi*phi_yy     !phi''(y)*dw
  !s5=


  if(pid.eq.0)write(*,'(6e12.5)')sum_c0/sum_liq0,solute_max,solute_min,solute_min/solute_max  
  if(pid.eq.0)print *,"Tip position and radius and tip width and facet curv"
  if(pid.eq.0)print *,temp3/temp1,temp2/temp4,0.8*temp1/temp2,temp6/temp8
  if(pid.eq.0)print *,"dt",dt
!  if(pid.eq.0)write(126,'(8e15.7)') time, temp3/temp1,temp2/temp4,solute_min,solute_max,solute_min/solute_max,temp1/temp2*0.8,temp6/temp8
  if(pid.eq.0)write(126,'(8e15.7)') time, temp3/temp1,temp2/temp4,solute_min,solute_max,solute_min/solute_max,temp1/temp2*0.8,sum_c0/sum_liq0

!  call MPI_REDUCE(phase,temp,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!  if(pid.eq.0)print *,"Amount of solid",temp

  call MPI_REDUCE(thetam,temp,1,amr_mpi_real,MPI_MIN,0,MPI_COMM_WORLD,ierr)
!   if ((nvar.eq.3.or.(.not.solute)).and.pid.eq.0)print *,"Sigma star temperature", temp

  if(pid.eq.0)write (6,'(A,A)') 27, "[0m"
  if(pid.eq.0)call flush(125)
  if(pid.eq.0)call flush(126)
! CEG removed profiles and p2, p3, p4, out1, out2, out3
1001 format(A15,I4) 
  if(time.gt.g_max_time)then
     write(*,*)"time greater than",g_max_time
     stop
  end if


end subroutine phase_field_check
!!!!!!!!!!
! PCB's routine for radius etc
subroutine  phase_loc_rad(lb, s1,s2,s3,s4)
  use paramesh_dimensions
  use physicaldata
  use tree

  implicit none
  double precision, intent(inout):: s1,s2,s3,s4
  integer i, j, k, lb
  double precision cur_loc, cur_rad, x1a, x1b, x2a, x2b, a5, b1, b3, b5, dx
  double precision x0, x1, x2,x3, xm1,phi,phi_x,phi_yy,phi_xx,radius,y1,xloc,L1,L2,L3
  double precision dw_phi
 integer isng
  dx = bsize(1,lb)/real(nxb)


!

  j = jl_bnd+nguard
!  y1 = bnd_box(1,2,lb)+(j-1.5)*dx
  k = kl_bnd+nguard*k3d
  cur_loc=1.0

  do i=il_bnd+nguard,iu_bnd-nguard
!     if (unk(1,i,j,k,lb).gt.1d-5.and.unk(1,i,j,k,lb).lt.1d0-1d-5) then
         x0 = bnd_box(1,1,lb)+(i-1.5)*dx
         phi_x  = 0.5*(unk(1,i+1,j,k,lb)-unk(1,i-1,j,k,lb))/dx
!         if(ndim.eq.3)then
!         phi_yy = 0.5*(-unk(1,i,j,k,lb)+unk(1,i,j+1,k+1,lb))/(dx*dx)
!         else
         phi_yy = (-unk(1,i,j,k,lb)+unk(1,i,j+1,k,lb))/(dx*dx)
!         phi_yy = (-15d0*unk(1,i,j,k,lb)+16d0*unk(1,i,j+1,k,lb)-unk(1,i,j+2,k,lb))/(6d0*dx*dx)
!         endif
!         b5 =  dx*atanh(unk(1,i,j,k,lb))/(atanh(unk(1,i,j,k,lb))-atanh(unk(1,i+1,j,k,lb)))
!         b5 =  dx*(unk(1,i,j,k,lb)-0.5)/((unk(1,i,j,k,lb))-(unk(1,i+1,j,k,lb)))
         phi=unk(1,i,j,k,lb)
         dw_phi=phi**2*(1d0-phi)**2
         s1=s1+dw_phi
         s2=s2+dw_phi*phi_x
         s3=s3+dw_phi*x0
         s4=s4+dw_phi*phi_yy
!     endif
  enddo
  
!  cur_loc=x0+b5

!  cur_rad = abs(phi_x/phi_yy)+b5  
!  cur_rad=abs(s2/s1)
!  cur_loc=s3/s1

end subroutine phase_loc_rad
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine  phase_loc_rad_y(lb, s1,s2,s3,s4)
  use paramesh_dimensions
  use physicaldata
  use tree

  implicit none
  double precision, intent(inout):: s1,s2,s3,s4
  integer i, j, k, lb
  double precision cur_loc, cur_rad, x1a, x1b, x2a, x2b, a5, b1, b3, b5, dx
  double precision x0, x1, x2,x3, xm1,phi,phi_x,phi_yy,phi_xx,radius,y1,xloc,L1,L2,L3
  double precision dw_phi
 integer isng
  dx = bsize(1,lb)/real(nxb)


!

  j = il_bnd+nguard
!  y1 = bnd_box(1,2,lb)+(j-1.5)*dx
  k = kl_bnd+nguard*k3d
  cur_loc=1.0

  do i=jl_bnd+nguard,ju_bnd-nguard
!     if (unk(1,i,j,k,lb).gt.1d-5.and.unk(1,i,j,k,lb).lt.1d0-1d-5) then
         x0 = bnd_box(1,2,lb)+(i-1.5)*dx
!         y0 =  bnd_box(1,2,lb)+(j-1.5)*dx
         phi_x  = 0.5*(unk(1,j,i+1,k,lb)-unk(1,j,i-1,k,lb))/dx
!         if(ndim.eq.3)then
!         phi_yy = 0.5*(-unk(1,i,j,k,lb)+unk(1,i,j+1,k+1,lb))/(dx*dx)
!         else
         phi_yy = (-unk(1,j,i,k,lb)+unk(1,j+1,i,k,lb))/(dx*dx)
!         phi_yy = (-15d0*unk(1,i,j,k,lb)+16d0*unk(1,i,j+1,k,lb)-unk(1,i,j+2,k,lb))/(6d0*dx*dx)
!         endif
!         b5 =  dx*atanh(unk(1,i,j,k,lb))/(atanh(unk(1,i,j,k,lb))-atanh(unk(1,i+1,j,k,lb)))
!         b5 =  dx*(unk(1,i,j,k,lb)-0.5)/((unk(1,i,j,k,lb))-(unk(1,i+1,j,k,lb)))
         phi=unk(1,j,i,k,lb)
         dw_phi=phi**2*(1d0-phi)**2
         s1=s1+dw_phi
         s2=s2+dw_phi*phi_x
         s3=s3+dw_phi*x0
         s4=s4+dw_phi*phi_yy
!     endif
  enddo
  
!  cur_loc=x0+b5

!  cur_rad = abs(phi_x/phi_yy)+b5  
!  cur_rad=abs(s2/s1)
!  cur_loc=s3/s1

end subroutine phase_loc_rad_y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! CEG routine for calculating tip location and radius using a series of quadratics

subroutine  phase_loc_rad_x(lb, cur_loc, cur_rad)
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  integer i, j, k, lb
  double precision cur_loc, cur_rad, x1a, x1b, x2a, x2b, a5, b1, b3, b5, dx
  double precision x0, x1, x2, xm1

 

  dx = bsize(1,lb)/real(nxb)

  j = jl_bnd+nguard
  k = kl_bnd+nguard*k3d
  do i=il_bnd+nguard,iu_bnd-nguard
     if (unk(1,i,j,k,lb).lt.0.5 .and. unk(1,i+1,j,k,lb).ge.0.5) then

        x1 = bnd_box(1,1,lb)+(i-1.5)*dx
        x0 = x1-dx
        xm1 = x0-dx
        x2 = x1+dx

        call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j,k,lb), unk(1,i,j,k,lb), unk(1,i+1,j,k,lb), x1a)
        if (unk(1,i,j+1,k,lb).lt.0.5.or.i.le.2) then
           call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j+1,k,lb), unk(1,i,j+1,k,lb), unk(1,i+1,j+1,k,lb), x1b)
        else
           call get_quadratic_interp(xm1, x0, x1, unk(1,i-2,j+1,k,lb), unk(1,i-1,j+1,k,lb), unk(1,i,j+1,k,lb), x1b)
        endif
        if (ndim.eq.3) then
           if (unk(1,i,j,k+1,lb).gt.0.0.or.i.le.2) then
              call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j,k+1,lb), unk(1,i,j,k+1,lb), unk(1,i+1,j,k+1,lb), x2a)
           else
              call get_quadratic_interp(xm1, x0, x1, unk(1,i-2,j,k+1,lb), unk(1,i-1,j,k+1,lb), unk(1,i,j,k+1,lb), x2a)
           endif
           if (unk(1,i,j+1,k+1,lb).gt.0.0.or.i.le.2) then
              call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j+1,k+1,lb), unk(1,i,j+1,k+1,lb), unk(1,i+1,j+1,k+1,lb), x2b)
           else
              call get_quadratic_interp(xm1, x0, x1, unk(1,i-2,j+1,k+1,lb), unk(1,i-1,j+1,k+1,lb), unk(1,i,j+1,k+1,lb), x2b)
           endif
        else
           x2a = x1a
           x2b = x1b
        endif
        exit
     endif
  enddo

  b1 = 0.125*(9.0*x1a-x2a)
  b3 = 0.125*(9.0*x1b-x2b)
  b5 = 0.125*(9.0*b1-b3)

  a5 = 0.5*(b3-b1)/(dx*dx)

  cur_loc=b5
  
! Wolfram-Alpha gives radius of curvature of y=f(x) as R=( (1+(y')^2)^1.5 )/|y''|
  cur_rad = 1.0/abs(2.0*a5)
end subroutine phase_loc_rad_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_quadratic_interp(x0, x1, x2, y0, y1, y2, xstar)
  double precision x0, x1, x2, y0, y1, y2, xstar, a, b, c
  double precision yy0,yy1,yy2
! next 3 lines PCB for phi=0 to 1
  yy0= 2.*(1.-y0)-1.
  yy1= 2.*(1.-y1)-1.
  yy2= 2.*(1.-y2)-1.
  b = ((yy0-yy2)*(x1*x1-x2*x2)-(yy1-yy2)*(x0*x0-x2*x2))/((x0-x2)*(x1*x1-x2*x2)-(x1-x2)*(x0*x0-x2*x2))
  a = (yy0-b*(x0-x2)-yy2) / (x0*x0-x2*x2)
  c = yy2-a*x2*x2-b*x2

! Choose -b MINUS (...)/2a as want left sided root
  xstar = (-b-sqrt(b*b-4.0*a*c))/(2.0*a)

end subroutine get_quadratic_interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function Gradoutput(vbles,lp,dx,gi,gj,lb)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp,gi,gj,lb
double precision :: phix(3),phiy(3),phixx(3),phiyy(3),phixy(3),LapPhi(5),eps,d_eps,gradientEnergy0,phi(3),s,phi_(3),U(3),dU(3)
integer :: i,j,k,ctr(6)=(/1,2,3,1,2,3/)
double precision :: d(3)=(/1,1,1/) 
double precision :: g0,X1,X2,theta,nu,A,A2,pi,kappa,kappa0,dh,dirac,epsilon0=1d0,x0,y0,xcoord,ycoord

  d(1) = 0.5*(g_eps22(1,2)+g_eps22(1,3)-g_eps22(2,3))
  d(2) = 0.5*(g_eps22(2,3)+g_eps22(2,1)-g_eps22(3,1))
  d(3) = 0.5*(g_eps22(3,1)+g_eps22(3,2)-g_eps22(1,2))
  
  d(2) = 1.0
  d(1) = g_d1
  d(3) = 1.0
  
  
  pi = 4.*atan2(1.,1.)

  phi(1)=vbles(1,2,2)
  phi(2)=vbles(2,2,2)
  phi(3)=vbles(3,2,2)
  if(lp.le.n_phases)then
     do i=1,n_phases
        call Calc_GradPhi(LapPhi(i),vbles(i,:,:),dx)
     enddo
     do i=1,n_phases
        phix(i) = 0.5*(vbles(i,3,2)-vbles(i,1,2))/dx
        phiy(i) = 0.5*(vbles(i,2,3)-vbles(i,2,1))/dx
        phixx(i)= (vbles(i,3,2)-2d0*vbles(i,2,2)+vbles(i,1,2))/dx**2
        phiyy(i)= (vbles(i,2,3)-2d0*vbles(i,2,2)+vbles(i,2,1))/dx**2
        phixy(i)= 0.25*(vbles(i,3,3)+vbles(i,1,1)-vbles(i,3,1)-vbles(i,1,3))/dx**2
     enddo
     U(1)=phi(1)+3*phi(2)*phi(3)
     U(2)=phi(2)+3*phi(3)*phi(1)
     U(3)=phi(3)+3*phi(1)*phi(2)
     Gradoutput = 0d0
     if(g_model.eq.4)then
       if(lp.eq.0)then
           do i=1,n_phases
              Gradoutput=Gradoutput + 0.5*d(i)*U(i)*(phix(i)**2+phiy(i)**2) 
           enddo
           return
        endif 


        g0  = sqrt(phix(lp)**2 +phiy(lp)**2)
        if(g0.gt.1d-9)then
           X1=max((phix(lp))/g0,1d-9)
           X2=max((phiy(lp))/g0,0d0)
           theta=atan2(abs(X2),abs(X1))
           if(abs(theta).le.pi/6.)then
              A =2*X1**2
              A2 =(-4*X1**2+4*X2**2)
           else
              A = (1.5-X1**2+sqrt(3.)*X1*X2)
              A2 = (2*X1**2-2*X2**2-4*sqrt(3.)*X1*X2)
           endif
           
           kappa =phixx(lp)*X2**2 - 2.* phixy(lp)*X1*X2 + phiyy(lp)*X1**2
!           kappa0 = min(0d0,kappa)
           
!           h0= 0.02

!           h0=min(0.1,abs(4./g_nuc_radius))
           dirac = 0d0
           if(abs(theta-0.5*pi).lt.h0)then

              dh = theta-0.5*pi
!              A =  ((11*h0**6-15*h0**4*dh**2+5*h0**2*dh**4-dh**6)*(-2*cos(Pi/6.-h0)*sin(Pi/6.-h0))+16*cos(pi/6.-h0)**2*h0**5)/(16*h0**5)
!              dirac = 15./16.*(1-dh/h0)**2*(1+dh/h0)**2/h0*2.*sqrt(3.)
!              dirac = 2.*sqrt(3.)*g0/dx
             dirac = 15./16.*(1-dh/h0)**2*(1+dh/h0)**2/h0*sqrt(3.)

!              dirac =(-30*(dh)**4+60*(dh)**2*h0**2-30*h0**4)*(-2*cos(Pi/6.-h0)*sin(Pi/6.-h0))/(8*h0**5)
!              A2 = A2 + dirac
           elseif(abs(theta-pi/6.).lt.h0)then
              dh = theta-pi/6.
!              A =  ((11*h0**6-15*h0**4*dh**2+5*h0**2*dh**4-dh**6)*(-2*cos(Pi/6.-h0)*sin(Pi/6.-h0))+16*cos(pi/6.-h0)**2*h0**5)/(16*h0**5)
!              dirac = 15./16.*(1-dh/h0)**2*(1+dh/h0)**2/h0*
             dirac = 15./16.*(1-dh/h0)**2*(1+dh/h0)**2/h0*sqrt(3.)

!              dirac = 2.*sqrt(3.)*g0/dx
!              dirac =(-30*(dh)**4+60*(dh)**2*h0**2-30*h0**4)*(-2*cos(Pi/6.-h0)*sin(Pi/6.-h0))/(8*h0**5)
!              A2 = A2 +dirac
           endif


        else
           A = 1d0
           A2 = 0
        endif



        dU(lp)=1d0
        dU(ctr(lp+1))=3d0*phi(ctr(lp+2))
        dU(ctr(lp+2))=3d0*phi(ctr(lp+1))

           Gradoutput = d(lp)*U(lp)*(LapPhi(lp)*A+0.5*A2*kappa+0.5*dirac*kappa)&
                - d(ctr(lp+0))*dU(ctr(lp+0))*0.5*(phix(ctr(lp+0))**2+phiy(ctr(lp+0))**2)*A&
                - d(ctr(lp+1))*dU(ctr(lp+1))*0.5*(phix(ctr(lp+1))**2+phiy(ctr(lp+1))**2)*A&
                - d(ctr(lp+2))*dU(ctr(lp+2))*0.5*(phix(ctr(lp+2))**2+phiy(ctr(lp+2))**2)*A





     else
        if(lp.eq.0)then
           do i=1,n_phases
              Gradoutput=Gradoutput + 0.5*d(i)*U(i)*(phix(i)**2+phiy(i)**2) 
           enddo
           return
        endif
        Gradoutput = d(ctr(lp))*U(ctr(lp))*LapPhi(ctr(lp))&
             + 3*d(ctr(lp))*phi(ctr(lp+2))*(phix(ctr(lp))*phix(ctr(lp+1)) + phiy(ctr(lp))*phiy(ctr(lp+1))) &
             + 3*d(ctr(lp))*phi(ctr(lp+1))*(phix(ctr(lp))*phix(ctr(lp+2)) + phiy(ctr(lp))*phiy(ctr(lp+2))) &
             +0.5*d(ctr(lp))*(phix(ctr(lp))*phix(ctr(lp)) + phiy(ctr(lp))*phiy(ctr(lp))) &
             -1.5*d(ctr(lp+1))*phi(ctr(lp+2))*(phix(ctr(lp+1))*phix(ctr(lp+1)) + phiy(ctr(lp+1))*phiy(ctr(lp+1))) &
             -1.5*d(ctr(lp+2))*phi(ctr(lp+1))*(phix(ctr(lp+2))*phix(ctr(lp+2)) + phiy(ctr(lp+2))*phiy(ctr(lp+2))) 
     endif


     Gradoutput = -Gradoutput

  else
     write(*,*)"in Gradoutput: lp > n_phases"
  endif

end function Gradoutput




! Rewritten by Peter Bollada

double precision function xcoord(i,lb)



  use paramesh_interfaces
  use refinement_parameters



  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use paramesh_comm_data ! Needed for definition of amr_mpi_real
  use tree
  use io
  use workspace
  implicit none
  include 'mpif.h'
integer, intent (in)::i,lb
double precision :: dx
  
  dx = bsize(1,lb)/real(nxb)
  xcoord =  bnd_box(1,1,lb)+(i-1.5)*dx
end function xcoord

double precision function ycoord(i,lb)
  use paramesh_interfaces
  use refinement_parameters
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use paramesh_comm_data ! Needed for definition of amr_mpi_real
  use tree
  use io
  use workspace
  implicit none
  include 'mpif.h'
integer, intent (in)::i,lb
double precision :: dx
  
  dx = bsize(1,lb)/real(nxb)
  ycoord =  bnd_box(1,2,lb)+(i-1.5)*dx
end function ycoord


subroutine pf_output_2d(pid, fnum)

  use paramesh_interfaces
  use refinement_parameters
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use paramesh_comm_data ! Needed for definition of amr_mpi_real
  use tree
  use io
  use workspace
  implicit none
  include 'mpif.h'
  integer pid,ierr
  integer i, j, k, v, lb, fnum,ip,ii,jj
  character*74 fname,gname
  character (len=2)  :: proc_string
  character (len=2)  :: grid_string
  character (len=6)  :: fnum_string
  double precision dx
  double precision x0,y0,p0,n0,xposition, FreeEnergy,potential,Fi(N_unk),s12(2),s23(2),s31(2),a12(2),aa12(2),a23(2),aa23(2),a31(2),aa31(2),s,hump,pi,theta1,theta2,theta3,theta4,t1,t2,t3,tt1,tt2,tt3
  double precision :: b12(2),b23(2),b31(2),bb12(2),bb23(2),bb31(2),dp(3,2),p(3),c12(2),c23(2),c31(2),cc12(2),cc23(2),cc31(2),d12(2),d23(2),d31(2)
  double precision:: t12,t23,t31,tt12,tt23,tt31,dp1(2),ddp1(2),dp2(2),ddp2(2),dp3(2),ddp3(2),dd12(2),dd23(2),dd31(2)
  double precision :: e12,e23,e31,ee12,ee23,ee31,f12,f23,f31,ff12,ff23,ff31,ee,ff,q12,q23,q31
  double precision :: a1(2),a2(2),a3(2),aa1(2),aa2(2),aa3(2),b,d,a
  double precision ::c,T,phi(n_phases),c_dot,phi_dot(3),rfactor,rf4,rf5,rf6,K_heat(4)
  double precision :: vbles(N_unk+1,3,3),G1,Gradoutput,xcoord,ycoord
  integer, dimension(16) :: Hdum= (/ 1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1,  0, -1 , 1, -1 /) 
  double precision, dimension(2) :: dH
  integer, dimension(2,8) :: H

  integer ::summ=0
 
  H = RESHAPE(Hdum,shape=(/2,8/))
     
   Write (proc_string, '(i2.2)') pid
   Write (grid_string, '(i2.2)') lrefine_max
   Write (fnum_string, '(i6.6)') fnum
    rfactor = dt/dtold
   rf4 = rfactor*rfactor/(rfactor+1.0)
   rf5 = (rfactor+1.0)
   rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
  do lb=1,lnblocks
     if(nodetype(lb).eq.1) then
        dx = bsize(1,lb)/real(nxb)

        do j=jl_bnd+1,ju_bnd-1
           do i=il_bnd+1,iu_bnd-1
              do v = 1, 3*nunkvbles
                 if (ABS(unk(v,i,j,1,lb)).lt.1.0e-20) unk(v,i,j,1,lb)=0.0
              end do
           end do
        end do
     end if
  end do
  !  close(140)
  
  fname = trim(output_dir)// '2dcontour2_' // grid_string // "_" // proc_string // "_" // fnum_string
!   gname = trim(output_dir)// 'bulk2_' // grid_string // "_" // proc_string // "_"
   gname = trim(output_dir)// 'bulk_' // grid_string // "_" // proc_string // "_"  // fnum_string
  write(*,*)fname,gname

  phi_dot=0.
  c_dot=0.
   open (unit=150,file=gname)
   b12=0
   b23=0
   b31=0
   d12=0
   d23=0
   d31=0
   a12=0
   a23=0
   a31=0
   c12=0
   c23=0
   c31=0

  do lb=1,lnblocks
    if(.not.has_children(lb)) then
      dx = bsize(1,lb)/real(nxb)
        do j=jl_bnd+1,ju_bnd-1
           do i=il_bnd+1,iu_bnd-1
!!$
!!$
!!$              do ip = 1,3
!!$                 p(ip)=unk(1+(ip-1)*nunkvbles,i,j,1,lb)
!!$                 dp(ip,1)=0.5*(unk(1+(ip-1)*nunkvbles,i+1,j,1,lb)-unk(1+(ip-1)*nunkvbles,i-1,j,1,lb))/dx
!!$                 dp(ip,2)=0.5*(unk(1+(ip-1)*nunkvbles,i,j+1,1,lb)-unk(1+(ip-1)*nunkvbles,i,j-1,1,lb))/dx
!!$              enddo
!!$
!!$
!!$              s12 = p(1)*p(2)*dp(3,:)
!!$              s23 = p(2)*p(3)*dp(1,:)
!!$              s31 = p(3)*p(1)*dp(2,:)
!!$
!!$
!!$
!!$
!!$              hump = dx*dx
!!$              
!!$              c12 = c12 + s12*hump
!!$              c23 = c23 + s23*hump
!!$              c31 = c31 + s31*hump
!!$              
!!$
!!$
!!$
!!$              hump = (p(1)*p(2)*p(3))*dx*dx
!!$
!!$              s12 = p(1)*dp(2,:)-p(2)*dp(1,:)
!!$              s23 = p(2)*dp(3,:)-p(3)*dp(2,:)
!!$              s31 = p(3)*dp(1,:)-p(1)*dp(3,:)
!!$
!!$              a12 = a12 + s12*hump
!!$              a23 = a23 + s23*hump
!!$              a31 = a31 + s31*hump
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$              hump =        dp(1,1)*dp(2,2)-dp(2,1)*dp(1,2)
!!$              hump = hump + dp(2,1)*dp(3,2)-dp(3,1)*dp(2,2)
!!$              hump = hump + dp(3,1)*dp(1,2)-dp(1,1)*dp(3,2)
!!$              hump = abs(hump)*dx*dx
!!$
!!$              d12 = d12 + s12*hump
!!$              d23 = d23 + s23*hump
!!$              d31 = d31 + s31*hump
!!$
!!$              hump=0
!!$              if(p(1)*p(2)*p(3).gt.0.1/27.)hump=1d0
!!$
!!$
!!$
!!$
!!$              b12 = b12 + s12*hump*dx*dx
!!$              b23 = b23 + s23*hump*dx*dx
!!$              b31 = b31 + s31*hump*dx*dx
!!$
!!$
!!$





             x0 =  bnd_box(1,1,lb)+(i-1.5)*dx
             y0 =  bnd_box(1,2,lb)+(j-1.5)*dx
             if(n_phases.eq.2)then
                phi(1)=unk(1+0*nunkvbles,i,j,1,lb)
                T=unk(1+(N_unk-2)*nunkvbles,i,j,1,lb)
                c=unk(1+(N_unk-1)*nunkvbles,i,j,1,lb)
                write(150,'(5(E16.5,X))')x0,y0,phi(1),T,c
             elseif(n_phases.eq.3)then

                phi(1)=unk(1+0*nunkvbles,i,j,1,lb)
                phi(2)=unk(1+1*nunkvbles,i,j,1,lb)
                phi(3)=unk(1+2*nunkvbles,i,j,1,lb)
                T=unk(1+(N_unk-2)*nunkvbles,i,j,1,lb)
                c=unk(1+(N_unk-1)*nunkvbles,i,j,1,lb)

                do ii=1,3
                   do jj=1,3
                      do ip=1,N_unk
                         vbles(ip,ii,jj)=unk(1+(ip-1)*nunkvbles,i+ii-2,j+jj-2,1,lb)
                      enddo
                   enddo
                enddo

                g1 = Gradoutput(vbles,1,dx,i,j,lb)


                write(150,'(7(E16.7,X))')x0,y0,phi(1),phi(2),phi(3),c,c


             endif
           end do
        end do
    endif
  end do

  close(150)
!!$  pi=180/3.14159
!!$
!!$
!!$
!!$  call MPI_REDUCE(b12,bb12,2,amr_mpi_real,MPI_sum,0,MPI_COMM_WORLD,ierr)
!!$  call MPI_REDUCE(b23,bb23,2,amr_mpi_real,MPI_sum,0,MPI_COMM_WORLD,ierr)
!!$  call MPI_REDUCE(b31,bb31,2,amr_mpi_real,MPI_sum,0,MPI_COMM_WORLD,ierr)
!!$
!!$
!!$
!!$  call MPI_REDUCE(d12,dd12,2,amr_mpi_real,MPI_sum,0,MPI_COMM_WORLD,ierr)
!!$  call MPI_REDUCE(d23,dd23,2,amr_mpi_real,MPI_sum,0,MPI_COMM_WORLD,ierr)
!!$  call MPI_REDUCE(d31,dd31,2,amr_mpi_real,MPI_sum,0,MPI_COMM_WORLD,ierr)
!!$
!!$  call MPI_REDUCE(a12,aa12,2,amr_mpi_real,MPI_sum,0,MPI_COMM_WORLD,ierr)
!!$  call MPI_REDUCE(a23,aa23,2,amr_mpi_real,MPI_sum,0,MPI_COMM_WORLD,ierr)
!!$  call MPI_REDUCE(a31,aa31,2,amr_mpi_real,MPI_sum,0,MPI_COMM_WORLD,ierr)
!!$
!!$  call MPI_REDUCE(c12,cc12,2,amr_mpi_real,MPI_sum,0,MPI_COMM_WORLD,ierr)
!!$  call MPI_REDUCE(c23,cc23,2,amr_mpi_real,MPI_sum,0,MPI_COMM_WORLD,ierr)
!!$  call MPI_REDUCE(c31,cc31,2,amr_mpi_real,MPI_sum,0,MPI_COMM_WORLD,ierr)
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$ if(pid.eq.0)then
!!$
!!$      bb12 = bb12/ sqrt((bb12(1))**2+(bb12(2))**2)
!!$      bb23 = bb23/ sqrt((bb23(1))**2+(bb23(2))**2)
!!$      bb31 = bb31/ sqrt((bb31(1))**2+(bb31(2))**2)
!!$
!!$      dd12 = dd12/ sqrt((dd12(1))**2+(dd12(2))**2)
!!$      dd23 = dd23/ sqrt((dd23(1))**2+(dd23(2))**2)
!!$      dd31 = dd31/ sqrt((dd31(1))**2+(dd31(2))**2)
!!$
!!$      aa12 = aa12/ sqrt((aa12(1))**2+(aa12(2))**2)
!!$      aa23 = aa23/ sqrt((aa23(1))**2+(aa23(2))**2)
!!$      aa31 = aa31/ sqrt((aa31(1))**2+(aa31(2))**2)
!!$
!!$      cc12 = cc12/ sqrt((cc12(1))**2+(cc12(2))**2)
!!$      cc23 = cc23/ sqrt((cc23(1))**2+(cc23(2))**2)
!!$      cc31 = cc31/ sqrt((cc31(1))**2+(cc31(2))**2)
!!$
!!$theta1 = acos(bb12(1)*bb31(1)+bb12(2)*bb31(2))*180/3.14159
!!$theta2 = acos(dd12(1)*dd31(1)+dd12(2)*dd31(2))*180/3.14159
!!$theta3 = acos(aa12(1)*aa31(1)+aa12(2)*aa31(2))*180/3.14159
!!$theta4 = acos(cc12(1)*cc31(1)+cc12(2)*cc31(2))*180/3.14159
!!$
!!$
!!$
!!$
!!$    open (unit=155,file="angle.txt",position = 'append')
!!$    write(155,'(5(F6.1,X))') theta1,theta2,theta3,theta4, 180-180./3.14159*2*asin(1/sqrt(2+2*g_d1**2))
!!$    close(155)
!!$ endif

end subroutine pf_output_2d

double precision function R2Z(x)
double precision, intent (in)::x
R2Z=x
if(abs(x).lt.1e-90)R2Z=0E0
end  function R2Z

