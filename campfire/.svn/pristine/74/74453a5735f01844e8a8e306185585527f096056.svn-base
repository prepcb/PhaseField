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
  implicit none
  integer, intent(in) :: inumber, pid, noprocs
  logical :: chombo = .false.
  logical chk_gd

  ! 2-d points around front
  if (ndim.eq.2) call pf_output_2d(pid, inumber)

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
  chk_checkf = 'hdf5'
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
  use refinement_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: inumber,pid
  integer :: ierr,lb,i,j,k,mid_x,mid_y,mid_z, found
  integer :: ii,jj
  double precision :: phase,umin,umax,phimin,phimax,u,phi,temp,gradphi,gradphi2,theta1,theta2
  double precision :: dx,dist,phi_val_1,phi_val_2,phi_val_3,phi_val_4
  double precision :: new_loc, new_rad, tempr,temp2,temp3,temp4,temp_max,temp_min
  double precision :: xl, xm, xr, phil, phim, phir, thetal, thetam, thetar, sigmastartemp
  logical :: on_axisx=.false.
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
  temp_max=-1d1
  temp_min= 1d1

!   if (nvar.eq.2.and.solute) found=1
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
                   
           if(unk(1,il_bnd+nguard,j,k,lb).lt.0.5.and.&
              unk(1,iu_bnd-nguard+1,j,k,lb).ge.0.5)then
!              print *, lb, unk(1,il_bnd+nguard,j,k,lb),unk(1,iu_bnd-nguard+1,j,k,lb)
              call phase_loc_rad(lb, new_loc, new_rad)
              found=found+1

!              print *,'fnd', unk(1,iu_bnd,j,k,lb)

           end if
 
           if (unk(1,il_bnd+nguard,j,k,lb).gt.-0.9.and.&
               unk(1,iu_bnd-nguard+1,j,k,lb).le.-0.9)then
              if (nvar.eq.3.or.(.not.solute)) then
                 do i = il_bnd+nguard, iu_bnd-nguard
                    if (unk(1,i+1,j,k,lb).ge.0.05) exit
                 end do
                 ! Set up local variables to do this double variable interpolation, first in X for xm then in theta for thetam
                 xl = bnd_box(1,1,lb)+(i-1.5)*dx
                 xr = bnd_box(1,1,lb)+(i+1-1.5)*dx
                 phil = unk(1,i,j,k,lb)
                 phim = -0.9
                 phir = 1.-2.*unk(1,i+1,j,k,lb)
                 thetal = 1.-2.*unk(1+nunkvbles,i,j,k,lb)
                 thetar = 1.-2.*unk(1+nunkvbles,i+1,j,k,lb)

                 xm = xl + (phim-phil)*(xr-xl)/(phir-phil)
                 thetam = thetal + (xm-xl)*(thetar-thetal)/(xr-xl)
                 found=found+1
              endif
                 
           end if
 
        end if
    
     end if
!     if (found.eq.2) exit
     
!     if(.not.has_children(lb).and.on_axisx) then
        do j=jl_bnd+1,ju_bnd-1
           do i=il_bnd+1,iu_bnd-1
              temp_max=max(unk(1+3*nunkvbles,i,j,1,lb),temp_max)
!              temp_min=min(unk(1+3*nunkvbles,i,j,1,lb),temp_min)
           enddo
        enddo
!     endif
  end do
!  print *,"Phase_field_check axis node count:",axis_node_count
!  call MPI_REDUCE(new_rad,tempr,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(temp_max,temp3,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
!  call MPI_REDUCE(temp_min,temp4,1,amr_mpi_real,MPI_MIN,0,MPI_COMM_WORLD,ierr)
  if(pid.eq.0)write (6,'(A,A)') 27, "[36;1;47m"
!  if(pid.eq.0)print *,"Tip Radius",temp3
!  if(pid.eq.0)write(*,*)temp3,temp4,temp3/temp4
!  if(pid.eq.0)write(unit=125,fmt=*) time, temp2, temp3

!  call MPI_REDUCE(dist,temp,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(new_loc,temp2,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  
  if(pid.eq.0)print *,"Tip position",temp2
  if(pid.eq.0)print *,"dt",dt
  if(pid.eq.0)write(unit=126,fmt=*) time, temp2,temp3
!  call MPI_REDUCE(phase,temp,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!  if(pid.eq.0)print *,"Amount of solid",temp

  call MPI_REDUCE(thetam,temp,1,amr_mpi_real,MPI_MIN,0,MPI_COMM_WORLD,ierr)
!   if ((nvar.eq.3.or.(.not.solute)).and.pid.eq.0)print *,"Sigma star temperature", temp

  if(pid.eq.0)write (6,'(A,A)') 27, "[0m"
  if(pid.eq.0)call flush(125)
  if(pid.eq.0)call flush(126)
! CEG removed profiles and p2, p3, p4, out1, out2, out3
1001 format(A15,I4) 



end subroutine phase_field_check
!!!!!!!!!!
! PCB's routine for radius etc
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


!

  j = jl_bnd+nguard
  y1 = bnd_box(1,2,lb)+(j-1.5)*dx
  k = kl_bnd+nguard*k3d
  cur_loc=1.0
  do i=il_bnd+nguard,iu_bnd-nguard
     if (unk(1,i+1,j,k,lb).ge.0.5.and.unk(1,i,j,k,lb).lt.0.5) then
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
! Linearly interpolated points
        ! Get x1a - x-position of interface on line y=h/2,z=h/2 using x-y plane
!        x1a=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j,k,lb)/(unk(1,i+1,j,k,lb)-unk(1,i,j,k,lb)))
        ! Get x1b - x-position of interface on line y=3h/2,z=h/2 using x-y plane
!        x1b=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j+1,k,lb)/(unk(1,i+1,j+1,k,lb)-unk(1,i,j+1,k,lb)))
        ! Get x2a - x-position of interface on line y=h/2,z=3h/2 using x-y plane
!        x2a=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j,k+1,lb)/(unk(1,i+1,j,k+1,lb)-unk(1,i,j,k+1,lb)))
        ! Get x2b - x-position of interface on line y=3h/2,z=3h/2 using x-y plane
!        x2b=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j+1,k+1,lb)/(unk(1,i+1,j+1,k+1,lb)-unk(1,i,j+1,k+1,lb)))

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

! Rewritten by Peter Bollada

subroutine pf_output_2d(pid, fnum)
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use io
  use workspace
  implicit none

  integer pid
  integer i, j, k, v, lb, fnum
  character*74 fname,gname
  character (len=2)  :: proc_string
  character (len=2)  :: grid_string
  character (len=6)  :: fnum_string
  double precision dx
  double precision x0,y0,p0,n0,xposition, FreeEnergy,potential
  double precision ::c,T,phi(3),c_dot,phi_dot(3),rfactor,rf4,rf5,rf6
  integer, dimension(16) :: Hdum= (/ 1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1,  0, -1 , 1, -1 /) 
  double precision, dimension(2) :: dH
  integer, dimension(2,8) :: H
  double precision :: cc,TT,pp(3),R2Z
  double precision :: Phi11(3,2),Phi12(3,2,2),LapPhi(3),GradientEnergy,psi,lappsi
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
!   gname = trim(output_dir)// 'bulk2_' // grid_string // "_" // proc_string // "_" // fnum_string
   gname = trim(output_dir)// 'bulk2_' // grid_string // "_" // proc_string // "_" 
  write(*,*)fname,gname

  phi_dot=0.
  c_dot=0.
   open (unit=150,file=gname)
  do lb=1,lnblocks
    if(.not.has_children(lb)) then
      dx = bsize(1,lb)/real(nxb)
        do j=jl_bnd+1,ju_bnd-1
           do i=il_bnd+1,iu_bnd-1
             x0 =  bnd_box(1,1,lb)+(i-1.5)*dx
             y0 =  bnd_box(1,2,lb)+(j-1.5)*dx
             phi(1)=unk(1+0*nunkvbles,i,j,1,lb)
             phi(2)=unk(1+1*nunkvbles,i,j,1,lb)
             phi(3)=unk(1+2*nunkvbles,i,j,1,lb)
             c =unk(1+4*nunkvbles,i,j,1,lb)
!             psi = unk(1+5*nunkvbles,i,j,1,lb)
!             lappsi = (unk(1+5*nunkvbles,i-1,j,1,lb)+unk(1+5*nunkvbles,i,j-1,1,lb)+unk(1+5*nunkvbles,i+1,j,1,lb)&
!                  +unk(1+5*nunkvbles,i,j+1,1,lb)-4d0*psi)/(dx*dx)
!             if(unk(1+3*nunkvbles,i,j,1,lb).gt.delta)then
                write(150,'(6(E16.5,X))')x0,y0,phi(2),phi(3),unk(1+3*nunkvbles,i,j,1,lb)/T_scale,c
!             endif
           end do
        end do
    endif
  end do
  close(150)

end subroutine pf_output_2d

double precision function R2Z(x)
double precision, intent (in)::x
R2Z=x
if(abs(x).lt.1e-90)R2Z=0E0
end  function R2Z

