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
  use paramesh_dimensions
  use time_dep_parameters
  use paramesh_interfaces
  use checkpoint_parameters
  use solution_parameters
  implicit none
  integer, intent(in) :: inumber, pid, noprocs
  logical :: chombo = .false.
  logical chk_gd

  ! 2-d points around front
  if (ndim.eq.2) call pf_output_2d(pid, inumber)
  if (ndim.eq.3.and.paraview_output) then
     call pf_output_3d(pid, inumber)!PCB output
!     call pf_output_3d_block(pid, inumber)!PCB output
     return
  end if
  ! Chombo output
  if (chombo)then
     if(pid.eq.0)then
        print *,"Output Chombo file"
     end if
     ! Write file to disk
     call amr_plotfile_chombo(inumber)
     write(61,*) 'Step, time', inumber, time
     flush(61) 
  end if

  ! Checkpointing output
  if(pid.eq.0) print *,"Checkpoint!"
  chk_gd = .FALSE.
  chk_checkf = "hdf5"
!  call amr_checkpoint_wr(inumber, user_attr_1=dt, user_attr_2=time, user_attr_3=dtold, check_format=chk_checkf)
  call amr_checkpoint_wr(inumber, user_attr_1=dt, user_attr_2=time, user_attr_3=dtold)
  if (pid.eq.0) then
     open (unit=99,file="CHK.out")
     write(99,*) inumber
     close(99)
  endif

  return
end subroutine app_output_occasional

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine phase_field_check(pid, inumber)
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
  implicit none
  include 'mpif.h'
  integer, intent(in) :: inumber,pid
  integer :: ierr,lb,i,j,k,mid_x,mid_y,mid_z, found
  double precision :: phase,umin,umax,phimin,phimax,u,phi,temp,gradphi,gradphi2,theta1,theta2
  double precision :: dx,dist,phi_val_1,phi_val_2,phi_val_3,phi_val_4
  double precision :: new_loc, new_rad, temp2
  double precision :: xl, xm, xr, phil, phim, phir, thetal, thetam, thetar, sigmastartemp
  logical :: on_axisx=.false.
  double precision :: G_energy,energy
!  character (len=19) :: filename
  !print *,"Guardcell axis node count:",axis_node_count
  !axis_node_count=0
  phase=0.0
  umin=1e30
  umax=-1e30
  phimin=1e30
  phimax=-1e30
  mid_x=il_bnd+nxb/2+nguard
  mid_y=jl_bnd+nyb/2+nguard
  if (ndim.eq.3) mid_z=kl_bnd+nzb/2+nguard
  found=0
  if(pid.eq.0)print *,"Phase field summary for time step: ",inumber
! CEG removed profiles
!  if(pid.gt.0)print *,"Warning, solute profiler not yet coded for parallelism"
!  if(mod(inumber,100).eq.0)then
!     write (filename,1001) "solute_profile_",inumber
!     if(pid.eq.0)open (unit=136,file=filename)
!     write (filename,1001) "phase_profile_",inumber
!     if(pid.eq.0)open (unit=137,file=filename)
!     write (filename,1001) "thermal_profile_",inumber
!     if(pid.eq.0)open (unit=138,file=filename)
!  end if
  

  dist = 0.0
  sigmastartemp = 0.0
  if (nvar.eq.2.and.solute) found=1
  do lb=1,lnblocks
     if(nodetype(lb).eq.1) then
        dx = bsize(1,lb)/real(nxb)
        ! First find out if we're "on axis"
        ! how do we know what the "axis" is?
        ! Simple the nucleate_(x,y,z) parameter holds the co-ords
        ! To be on an axis 2 out of 3 of the block faces must match these co-ords
        ! Only testing lower co-ords so for corner nucleation the nucleus MUST BE at min co-ords for each dimension
        on_axisx=.false.
        
        if(bnd_box(1,2,lb).eq.nucleate_y)then
           ! lower y bnd matches
           if(ndim.eq.2.or.bnd_box(1,3,lb).eq.nucleate_z)then
              ! lower z bnd matches
              ! y and z co-ords match so we're on the x axis
              on_axisx=.true.
              !axis_node_count=axis_node_count+1
           end if
        end if
        j=jl_bnd+nguard
        k=kl_bnd+nguard*k3d

        if(on_axisx)then
           ! On x axis now test whether the interface is within this block
           ! Note testing from il_bnd to iu_bnd
           ! This ensure that if interface is between 2 blocks we still find it
           ! Basically we're checking guardcells as well as computational domain
           if(unk(1,il_bnd+nguard,j,k,lb).gt.0.0.and.&
              unk(1,iu_bnd-nguard+1,j,k,lb).le.0.0)then
!              print *, lb, unk(1,il_bnd+nguard,j,k,lb),unk(1,iu_bnd-nguard+1,j,k,lb)
              call phase_loc_rad(lb, new_loc, new_rad)
              found=found+1
              

!              print *,'fnd', unk(1,iu_bnd,j,k,lb)

           end if

           if (unk(1,il_bnd+nguard,j,k,lb).gt.-0.9.and.&
               unk(1,iu_bnd-nguard+1,j,k,lb).le.-0.9)then
              if (nvar.eq.3.or.(.not.solute)) then
                 do i = il_bnd+nguard, iu_bnd-nguard
                    if (unk(1,i+1,j,k,lb).le.-0.9) exit
                 end do
                 ! Set up local variables to do this double variable interpolation, first in X for xm then in theta for thetam
                 xl = bnd_box(1,1,lb)+(i-1.5)*dx
                 xr = bnd_box(1,1,lb)+(i+1-1.5)*dx
                 phil = unk(1,i,j,k,lb)
                 phim = -0.9
                 phir = unk(1,i+1,j,k,lb)
                 thetal = unk(1+nunkvbles,i,j,k,lb)
                 thetar = unk(1+nunkvbles,i+1,j,k,lb)

                 xm = xl + (phim-phil)*(xr-xl)/(phir-phil)
                 thetam = thetal + (xm-xl)*(thetar-thetal)/(xr-xl)
                 found=found+1
              endif
           end if
        end if
     end if
!     if (found.eq.2) exit
  end do

! Energy compuation - should decrease
  call amr_mg_get_Energy(pid,Energy)
  call MPI_Reduce(Energy, G_Energy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
  if (pid.eq.0)write(*,*) 'Post smooths Energy Calculation',G_energy
  
!  print *,"Phase_field_check axis node count:",axis_node_count
  call MPI_REDUCE(new_rad,temp2,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  new_rad=temp2
  call MPI_REDUCE(new_loc,temp2,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  new_loc=temp2
  if(pid.eq.0)write (6,'(A,A)') 27, "[36;1;47m"
  if(pid.eq.0)print *,"Tip Radius",new_rad
  if(pid.eq.0.and.new_loc.gt.2e0)write(unit=125,fmt=*) time, new_rad, G_energy
  call MPI_REDUCE(dist,temp,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  
  if(pid.eq.0)print *,"Tip position",new_loc
  if(pid.eq.0.and.new_loc.gt.2e0)write(unit=126,fmt=*) time, new_loc

  call MPI_REDUCE(thetam,temp,1,amr_mpi_real,MPI_MIN,0,MPI_COMM_WORLD,ierr)
  if ((nvar.eq.3.or.(.not.solute)).and.pid.eq.0)print *,"Sigma star temperature", temp

  if(pid.eq.0)write (6,'(A,A)') 27, "[0m"
  if(pid.eq.0)call flush(125)
  if(pid.eq.0)call flush(126)

  ! Get out for failing code: only if tip_location is tiny we must have gone wrong
  if (pid.eq.0.and.abs(temp2).lt.-1.0e6) call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
! CEG removed profiles and p2, p3, p4, out1, out2, out3
1001 format(A15,I4) 
end subroutine phase_field_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! CEG routine for calculating tip location and radius using a series of quadratics

subroutine  phase_loc_rad(lb, cur_loc, cur_rad)
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  integer i, j, k, lb
  double precision cur_loc, cur_rad, x1a, x1b, x2a, x2b, a5, b1, b3, b5, dx
  double precision x0, x1, x2,x3, xm1,phi,phi_x,phi_yy,phi_xx,radius,y1,xloc,L1,L2,L3
  integer isng
  dx = bsize(1,lb)/real(nxb)

  j = jl_bnd+nguard
  y1 = bnd_box(1,2,lb)+(j-1.5)*dx
  k = kl_bnd+nguard*k3d
  cur_loc=1.0
  do i=il_bnd+nguard,iu_bnd-nguard
     if (unk(1,i+1,j,k,lb).lt.0.0.and.unk(1,i,j,k,lb).gt.0.0) then
         x0 = bnd_box(1,1,lb)+(i-1.5)*dx
         phi_x  = 0.5*(unk(1,i+1,j,k,lb)-unk(1,i-1,j,k,lb))/dx
         if(ndim.eq.3)then
         phi_yy = 0.5*(-(unk(1,i,j,k,lb))+(unk(1,i,j+1,k+1,lb)))/(dx*dx)
         else
         phi_yy = (-(unk(1,i,j,k,lb))+(unk(1,i,j+1,k,lb)))/(dx*dx)
         endif
!         b5 =  dx*atanh(unk(1,i,j,k,lb))/(atanh(unk(1,i,j,k,lb))-atanh(unk(1,i+1,j,k,lb)))
         b5 =  dx*(unk(1,i,j,k,lb))/((unk(1,i,j,k,lb))-(unk(1,i+1,j,k,lb)))
     endif
  enddo
  
  cur_loc=x0+b5
  cur_rad = abs(phi_x/phi_yy)+b5  
end subroutine phase_loc_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_quadratic_interp(x0, x1, x2, y0, y1, y2, xstar)
  double precision x0, x1, x2, y0, y1, y2, xstar, a, b, c

  b = ((y0-y2)*(x1*x1-x2*x2)-(y1-y2)*(x0*x0-x2*x2))/((x0-x2)*(x1*x1-x2*x2)-(x1-x2)*(x0*x0-x2*x2))
  a = (y0-b*(x0-x2)-y2) / (x0*x0-x2*x2)
  c = y2-a*x2*x2-b*x2

! Choose -b MINUS (...)/2a as want left sided root
  xstar = (-b-sqrt(b*b-4.0*a*c))/(2.0*a)

end subroutine get_quadratic_interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Rewritten by Peter Bollada

subroutine pf_output_2d(pid, fnum)
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use io

  implicit none

  integer pid
  integer i, j, k, v, lb, fnum
  character*80 fname
  character (len=2)  :: proc_string
  character (len=2)  :: grid_string
  character (len=5)  :: fnum_string
  double precision dx
  double precision x0,y0,p0,n0
  integer, dimension(16) :: Hdum= (/ 1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1,  0, -1 , 1, -1 /) 
  double precision, dimension(2) :: dH
  integer, dimension(2,8) :: H

  H = RESHAPE(Hdum,shape=(/2,8/))
     
!  write(fname,'(A,I0.2)'),"/tmp/2dout_",lrefine_max,"_",fnum 
   Write (proc_string, '(i2.2)') pid
   Write (grid_string, '(i2.2)') lrefine_max
   Write (fnum_string, '(i5.5)') fnum
!  fname = trim(output_dir) // '2dout_' // grid_string // "_" // fnum_string
!  open (unit=140,file=fname)
 
  do lb=1,lnblocks
     if(nodetype(lb).eq.1) then
        dx = bsize(1,lb)/real(nxb)

        do j=jl_bnd+1,ju_bnd-1
           do i=il_bnd+1,iu_bnd-1
              do v = 1, 15
                 if (ABS(unk(v,i,j,1,lb)).lt.1.0e-15) unk(v,i,j,1,lb)=0.0
              end do
!              write(140,'(E12.5,X,E12.5,X,E12.5,X,E12.5,X,E12.5,X)') &
!                 bnd_box(1,1,lb)+(i-1.5)*dx, bnd_box(1,2,lb)+(j-1.5)*dx, unk(1,i,j,1,lb), unk(6,i,j,1,lb), unk(11,i,j,1,lb)
           end do
        end do
     end if
  end do
  !  close(140)
  
  fname = trim(output_dir)// '/2dcontour_' // grid_string // "_" // proc_string // "_" // fnum_string
!  fname = 'output/'// '2dcontour_' // grid_string // "_" // !proc_string // "_" // fnum_string
  
  !  write(fname,'(A,I0.2)'),"/tmp/2dcontour_",lrefine_max,"_",fnum
  open (unit=140,file=fname)
  
  
  do lb=1,lnblocks
     if(nodetype(lb).eq.1) then
        
        dx = bsize(1,lb)/real(nxb)
        
        do j=jl_bnd+1,ju_bnd-1
           do i=il_bnd+1,iu_bnd-1
                 if (unk(1,i,j,1,lb).ge.0.0.and.&
             (   unk(1,i+1,j  ,1,lb).lt.0.0.or.&
                 unk(1,i+1,j+1,1,lb).lt.0.0.or.&
                 unk(1,i  ,j+1,1,lb).lt.0.0.or.&
                 unk(1,i-1,j+1,1,lb).lt.0.0.or.&
                 unk(1,i-1,j  ,1,lb).lt.0.0.or.&
                 unk(1,i-1,j-1,1,lb).lt.0.0.or.&
                 unk(1,i  ,j-1,1,lb).lt.0.0.or.&
                 unk(1,i+1,j-1,1,lb).lt.0.0)&
                    )then 
                x0 =  bnd_box(1,1,lb)+(i-1.5)*dx
                y0 =  bnd_box(1,2,lb)+(j-1.5)*dx
                n0=0.0
                p0= unk(1,i,j,1,lb)
                do k=1,8
                   if( unk(1,i+H(1,k),j+H(2,k),1,lb)<n0)then
                       n0=unk(1,i+H(1,k),j+H(2,k)  ,1,lb)
                       dH=H(:,k)
                   end if
                end do
                  
                write(140,'(5(E12.5,X))') &
                (-n0*x0+p0*(x0+dH(1)*dx))/(p0-n0),&
                (-n0*y0+p0*(y0+dH(2)*dx))/(p0-n0),&
                unk(1,i  ,j  ,1,lb),&    
                unk(6,i,j,1,lb), unk(11,i,j,1,lb)
                 
              endif
           end do
        end do
     end if
  end do
  
end subroutine pf_output_2d

! Rewritten by Peter Bollada

subroutine pf_output_3d_old(pid, fnum)
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use io

  implicit none

  integer pid
  integer i, j, k, v, lb, fnum,n_count
  character*80 fname
  character (len=3)  :: proc_string
  character (len=2)  :: grid_string
  character (len=5)  :: fnum_string
  double precision dx,dy,dz,x_reduce
  double precision x0,y0,p0,n0,z0
  integer, dimension(16) :: Hdum= (/ 1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1,  0, -1 , 1, -1 /) 
  double precision, dimension(2) :: dH
  integer, dimension(2,8) :: H
  double precision :: u(3)
  
  H = RESHAPE(Hdum,shape=(/2,8/))
  
  !  write(fname,'(A,I0.2)'),"/tmp/2dout_",lrefine_max,"_",fnum 
  Write (proc_string, '(i3.3)') pid
  Write (grid_string, '(i2.2)') lrefine_max
  Write (fnum_string, '(i5.5)') fnum
!!!!!!!!!!!!! 1D output x-axis output !!!!!!!!!!!!!!!!   
  fname = trim(output_dir) // '2dout_' // grid_string // "_" // proc_string // "_" // fnum_string
  open (unit=140,file=fname)
  do lb=1,lnblocks
     if (.not.has_children(lb).and.bnd_box(1,2,lb).eq.0.and.bnd_box(1,3,lb).eq.0)then
        dx = bsize(1,lb)/real(nxb)
        do i=il_bnd+1,iu_bnd-1
           u(1)=0.25*(unk(1,i,1,1,lb)+unk(1,i,1,2,lb)+unk(1,i,2,1,lb)+unk(1,i,2,2,lb))
           u(2)=0.25*(unk(6,i,1,1,lb)+unk(6,i,1,2,lb)+unk(6,i,2,1,lb)+unk(6,i,2,2,lb))
           u(3)=0.25*(unk(11,i,1,1,lb)+unk(11,i,1,2,lb)+unk(11,i,2,1,lb)+unk(11,i,2,2,lb))
           write(140,'(E12.5,X,E12.5,X,E12.5,X,E12.5,X)') &
                bnd_box(1,1,lb)+(i-1.5)*dx,  u
        end do
     end if
  end do
  close(140)

!!!!!!!!!!!!!! 3D interface output !!!!!!!!!!!!!!!  
  fname = trim(output_dir)// '/3dcontour_' // grid_string // "_" // proc_string // "_" // fnum_string
  open (unit=140,file=fname)
  do lb=1,lnblocks
     if(.not.has_children(lb)) then
        dx = bsize(1,lb)/real(nxb)
        do k=kl_bnd, ku_bnd-nguard*k3d 
           do j=jl_bnd+1,ju_bnd-nguard*k2d
              do i=il_bnd+1,iu_bnd-nguard
                 if (unk(1,i,j,k,lb).ge.-0.2.and.unk(1,i,j,k,lb).le.0.2)then
                    x0 =  bnd_box(1,1,lb)+(i-1.5)*dx
                    y0 =  bnd_box(1,2,lb)+(j-1.5)*dx
                    z0 =  bnd_box(1,3,lb)+(k-1.5)*dx
                    write(140,'(3(F12.5,X))') x0,y0,z0
                 endif
              end do
           end do
        enddo
     end if
  end do
  close(140)
end subroutine pf_output_3d_old


subroutine pf_output_3d(pid, fnum)
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use io

  implicit none

  integer pid
  integer i, j, k, v, lb, fnum,n_count,ii
  character*80 fname
  character (len=3)  :: proc_string
  character (len=2)  :: grid_string
  character (len=5)  :: fnum_string
  double precision dx,dy,dz,x_reduce
  double precision x0,y0,p0,n0,z0,phix,phiy,phiz,gradphisq
  integer, dimension(16) :: Hdum= (/ 1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1,  0, -1 , 1, -1 /) 
  double precision, dimension(2) :: dH
  integer, dimension(2,8) :: H
  double precision :: u(3), cell(8)
  double precision, dimension(3,8) :: xx
  
  H = RESHAPE(Hdum,shape=(/2,8/))
  write(*,*)"writing..................................."
  !  write(fname,'(A,I0.2)'),"/tmp/2dout_",lrefine_max,"_",fnum 
  Write (proc_string, '(i3.3)') pid
  Write (grid_string, '(i2.2)') lrefine_max
  Write (fnum_string, '(i5.5)') fnum


!!!!!!!!!!!!!! 3D interface output !!!!!!!!!!!!!!!  
 fname=trim(output_dir)// '/iso_' // grid_string // "_" // fnum_string // "_" //proc_string//".vtk"
 write(*,*)fname
  open (unit=140,file=fname)
  
  write(140,'(A)')"# vtk DataFile Version 2.0"
  write(140,'(A)')"Unstructured Grid Example"
  write(140,'(A)')"ASCII"
  write(140,*)
  write(140,'(A)')"DATASET UNSTRUCTURED_GRID"
  n_count=0 
  do lb=1,lnblocks
     if(.not.has_children(lb)) then
        dx = bsize(1,lb)/real(nxb)
        do k=1,9 
           do j=1,9
              do i=1,9
                 cell(1)=unk(1,i,j,k,lb)
                 cell(2)=unk(1,i+1,j,k,lb)
                 cell(3)=unk(1,i,j+1,k,lb)
                 cell(4)=unk(1,i+1,j+1,k,lb)
                 cell(5)=unk(1,i,j,k+1,lb)
                 cell(6)=unk(1,i+1,j,k+1,lb)
                 cell(7)=unk(1,i,j+1,k+1,lb)
                 cell(8)=unk(1,i+1,j+1,k+1,lb)

                 if (minval(cell).le.0d0.and. maxval(cell).gt.0d0)then ! found cell
                    n_count=n_count+1
                 endif
              end do
           end do
        enddo
     end if
  end do
  write(140,'(A,I10,A)')"POINTS", n_count*8, "  FLOAT"
  do lb=1,lnblocks
     if(.not.has_children(lb)) then
        dx = bsize(1,lb)/real(nxb)
        do k=1,9 
           do j=1,9
              do i=1,9
                 cell(1)=unk(1,i,j,k,lb)
                 cell(2)=unk(1,i+1,j,k,lb)
                 cell(3)=unk(1,i,j+1,k,lb)
                 cell(4)=unk(1,i+1,j+1,k,lb)
                 cell(5)=unk(1,i,j,k+1,lb)
                 cell(6)=unk(1,i+1,j,k+1,lb)
                 cell(7)=unk(1,i,j+1,k+1,lb)
                 cell(8)=unk(1,i+1,j+1,k+1,lb)

                 if (minval(cell).le.0d0.and. maxval(cell).gt.0d0)then ! found cell

                    xx(1,1) =  bnd_box(1,1,lb)+(i-1.5)*dx
                    xx(2,1) =  bnd_box(1,2,lb)+(j-1.5)*dx
                    xx(3,1) =  bnd_box(1,3,lb)+(k-1.5)*dx
                    xx(1,2) =  bnd_box(1,1,lb)+(i+1-1.5)*dx
                    xx(2,2) =  bnd_box(1,2,lb)+(j-1.5)*dx
                    xx(3,2) =  bnd_box(1,3,lb)+(k-1.5)*dx
                    xx(1,3) =  bnd_box(1,1,lb)+(i-1.5)*dx
                    xx(2,3) =  bnd_box(1,2,lb)+(j+1-1.5)*dx
                    xx(3,3) =  bnd_box(1,3,lb)+(k-1.5)*dx
                    xx(1,4) =  bnd_box(1,1,lb)+(i+1-1.5)*dx
                    xx(2,4) =  bnd_box(1,2,lb)+(j+1-1.5)*dx
                    xx(3,4) =  bnd_box(1,3,lb)+(k-1.5)*dx
                    xx(1,5) =  bnd_box(1,1,lb)+(i-1.5)*dx
                    xx(2,5) =  bnd_box(1,2,lb)+(j-1.5)*dx
                    xx(3,5) =  bnd_box(1,3,lb)+(k+1-1.5)*dx
                    xx(1,6) =  bnd_box(1,1,lb)+(i+1-1.5)*dx
                    xx(2,6) =  bnd_box(1,2,lb)+(j-1.5)*dx
                    xx(3,6) =  bnd_box(1,3,lb)+(k+1-1.5)*dx
                    xx(1,7) =  bnd_box(1,1,lb)+(i-1.5)*dx
                    xx(2,7) =  bnd_box(1,2,lb)+(j+1-1.5)*dx
                    xx(3,7) =  bnd_box(1,3,lb)+(k+1-1.5)*dx
                    xx(1,8) =  bnd_box(1,1,lb)+(i+1-1.5)*dx
                    xx(2,8) =  bnd_box(1,2,lb)+(j+1-1.5)*dx
                    xx(3,8) =  bnd_box(1,3,lb)+(k+1-1.5)*dx
                    write(140,'(3(F12.5,X))') xx

                 endif
              end do
           end do
        enddo
     end if
  end do
  write(140,'(A,2I10)')"CELLS", n_count,9*n_count
  do i =1,n_count
     write(140,'(9(I10,X))') 8,(ii+(i-1)*8,ii=0,7)
  end do
  write(140,*)
    write(140,'(A,I10)')"CELL_TYPES", n_count
  do i =1, n_count
     write(140,'(I10)')11
  enddo
  write(140,'(A,I10)')"POINT_DATA",n_count*8
  write(140,'(A)')"SCALARS phi float 1"
  write(140,'(A)')"LOOKUP_TABLE default"
 do lb=1,lnblocks
     if(.not.has_children(lb)) then
        dx = bsize(1,lb)/real(nxb)
        do k=1,9
           do j=1,9
              do i=1,9
                 cell(1)=unk(1,i,j,k,lb)
                 cell(2)=unk(1,i+1,j,k,lb)
                 cell(3)=unk(1,i,j+1,k,lb)
                 cell(4)=unk(1,i+1,j+1,k,lb)
                 cell(5)=unk(1,i,j,k+1,lb)
                 cell(6)=unk(1,i+1,j,k+1,lb)
                 cell(7)=unk(1,i,j+1,k+1,lb)
                 cell(8)=unk(1,i+1,j+1,k+1,lb)

                 if (minval(cell).le.0d0.and. maxval(cell).gt.0d0)then
                    write(140,'(8(F12.5,X))') cell
                 endif
              end do
           end do
        enddo
     end if
  end do


  close(140)
end subroutine pf_output_3d


!!!!! -----------------------------------------------------------------
subroutine pf_output_3d_block(pid, fnum)
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use io

  implicit none

  integer pid
  integer i, j, k, v, lb, fnum,n_count, ii, jj, kk
  character*80 fname
  character (len=3)  :: proc_string
  character (len=2)  :: grid_string
  character (len=5)  :: fnum_string
  character (len=6)  :: blevel_string
  double precision dx,dy,dz,x_reduce
  double precision x0,y0,p0,n0,z0,phix,phiy,phiz,gradphisq
  integer, dimension(16) :: Hdum= (/ 1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1,  0, -1 , 1, -1 /) 
  double precision, dimension(2) :: dH
  integer, dimension(2,8) :: H
  double precision :: u(3), cell(8)
  double precision, dimension(3,8) :: xx
  H = RESHAPE(Hdum,shape=(/2,8/))
  write(*,*)"writing..................................."
  !  write(fname,'(A,I0.2)'),"/tmp/2dout_",lrefine_max,"_",fnum 
  Write (proc_string, '(i3.3)') pid
  Write (grid_string, '(i2.2)') lrefine_max
  Write (fnum_string, '(i5.5)') fnum

!!!!!!!!!!!!!! 3D interface output !!!!!!!!!!!!!!!  
  do ii = 7, 7
     Write (blevel_string, '(i2.2)') ii
     fname=trim(output_dir)// '/isoceg_' // grid_string // "_" // fnum_string // "_" //proc_string//"_"//blevel_string
     open (unit=ii+139,file=fname)
  enddo
  n_count=0 
  do lb=1,lnblocks
     if(.not.has_children(lb)) then
        dx = bsize(1,lb)/real(nxb)
        do k=1,9 
           do j=1,9
              do i=1,9
                 cell(1)=unk(1,i,j,k,lb)
                 cell(2)=unk(1,i+1,j,k,lb)
                 cell(3)=unk(1,i,j+1,k,lb)
                 cell(4)=unk(1,i+1,j+1,k,lb)
                 cell(5)=unk(1,i,j,k+1,lb)
                 cell(6)=unk(1,i+1,j,k+1,lb)
                 cell(7)=unk(1,i,j+1,k+1,lb)
                 cell(8)=unk(1,i+1,j+1,k+1,lb)
                 if (minval(cell).le.0d0.and. maxval(cell).gt.0d0)then ! found cell
                    n_count=n_count+1
                 endif
              end do
           end do
        enddo
     end if
  end do
  do lb=1,lnblocks
     if(.not.has_children(lb)) then
        dx = bsize(1,lb)/real(nxb)
        do k=5,5 
           do j=5,5
              do i=5,5
!                 cell(1)=unk(1,i,j,k,lb)
!                 cell(2)=unk(1,i+1,j,k,lb)
!                 cell(3)=unk(1,i,j+1,k,lb)
!                 cell(4)=unk(1,i+1,j+1,k,lb)
!                 cell(5)=unk(1,i,j,k+1,lb)
!                 cell(6)=unk(1,i+1,j,k+1,lb)
!                 cell(7)=unk(1,i,j+1,k+1,lb)
!                 cell(8)=unk(1,i+1,j+1,k+1,lb)

!                 if (minval(cell).le.0d0.and. maxval(cell).gt.0d0)then ! found cell

                    xx(1,1) =  bnd_box(1,1,lb)
                    xx(2,1) =  bnd_box(1,2,lb)
                    xx(3,1) =  bnd_box(1,3,lb)
!                    xx(1,1) =  bnd_box(1,1,lb)+(i-1.5)*dx
!                    xx(2,1) =  bnd_box(1,2,lb)+(j-1.5)*dx
!                    xx(3,1) =  bnd_box(1,3,lb)+(k-1.5)*dx
!                    xx(1,2) =  bnd_box(1,1,lb)+(i+1-1.5)*dx
!                    xx(2,2) =  bnd_box(1,2,lb)+(j-1.5)*dx
!                    xx(3,2) =  bnd_box(1,3,lb)+(k-1.5)*dx
!                    xx(1,3) =  bnd_box(1,1,lb)+(i-1.5)*dx
!                    xx(2,3) =  bnd_box(1,2,lb)+(j+1-1.5)*dx
!                    xx(3,3) =  bnd_box(1,3,lb)+(k-1.5)*dx
!                    xx(1,4) =  bnd_box(1,1,lb)+(i+1-1.5)*dx
!                    xx(2,4) =  bnd_box(1,2,lb)+(j+1-1.5)*dx
!                    xx(3,4) =  bnd_box(1,3,lb)+(k-1.5)*dx
!                    xx(1,5) =  bnd_box(1,1,lb)+(i-1.5)*dx
!                    xx(2,5) =  bnd_box(1,2,lb)+(j-1.5)*dx
!                    xx(3,5) =  bnd_box(1,3,lb)+(k+1-1.5)*dx
!                    xx(1,6) =  bnd_box(1,1,lb)+(i+1-1.5)*dx
!                    xx(2,6) =  bnd_box(1,2,lb)+(j-1.5)*dx
!                    xx(3,6) =  bnd_box(1,3,lb)+(k+1-1.5)*dx
!                    xx(1,7) =  bnd_box(1,1,lb)+(i-1.5)*dx
!                    xx(2,7) =  bnd_box(1,2,lb)+(j+1-1.5)*dx
!                    xx(3,7) =  bnd_box(1,3,lb)+(k+1-1.5)*dx
!                    xx(1,8) =  bnd_box(1,1,lb)+(i+1-1.5)*dx
!                    xx(2,8) =  bnd_box(1,2,lb)+(j+1-1.5)*dx
!                    xx(3,8) =  bnd_box(1,3,lb)+(k+1-1.5)*dx

                    if(lrefine(lb)==1) then
!                       write(140,'(4(E12.5,X))') xx(1,1),xx(2,1),xx(3,1),real(pid*1.0d0)
                    elseif(lrefine(lb)==2) then
!                       write(141,'(4(E12.5,X))') xx(1,1),xx(2,1),xx(3,1),real(pid*1.0d0)
                    elseif(lrefine(lb)==3) then
!                       write(142,'(4(E12.5,X))') xx(1,1),xx(2,1),xx(3,1),real(pid*1.0d0)
                    elseif(lrefine(lb)==4) then 
!                       write(143,'(4(E12.5,X))') xx(1,1),xx(2,1),xx(3,1),real(pid*1.0d0)
                    elseif(lrefine(lb)==5) then
!                       write(144,'(4(E12.5,X))') xx(1,1),xx(2,1),xx(3,1),real(pid*1.0d0)
                    elseif(lrefine(lb)==6) then
!                       write(145,'(4(E12.5,X))') xx(1,1),xx(2,1),xx(3,1),real(pid*1.0d0)
                    elseif(lrefine(lb)==7) then
                       write(146,'(4(E12.5,X))') xx(1,1),xx(2,1),xx(3,1),real(pid*1.0d0)
                    endif
               !  endif
              end do
           end do
        enddo
     end if
  end do


  close(140)
end subroutine pf_output_3d_block
