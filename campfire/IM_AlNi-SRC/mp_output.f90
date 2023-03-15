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
  double precision :: T_

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



  if(g_quad.or. (.not. new_potential))then
     if(g_alpha.ne.0d0)then
        
        g_ratio = g_ratio0
        if(g_vel.gt.0d0.or. .true.)then
           if(AlNi.eq.8)call find_antitrapping_AlFe
           if(pid.eq.0)then 
              write(*,*)"antitrapping =", g_alpha,g_alpha*g_ratio
              open (unit=270,file="Antitrapping.txt")
              write(270,'(3F14.7)')g_alpha,g_c_min,g_c_max
              close(270)
           endif
        else
           g_alpha=1d-6
        endif
     endif
     
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
  double precision :: phase,umin,umax,phimin,phimax,u,phi,temp,gradphi,gradphi2,theta1,theta2,phi_x,phi_y
  double precision :: dx,dist,phi_val_1,phi_val_2,phi_val_3,phi_val_4,solute_max,solute_min,s1,s2,s3,s4,s5,s6,s7,sy1,sy2,sy3,sy4
  double precision :: new_loc, new_rad,new_DW, tempr,temp1,temp2,temp3,temp4,temp5,sol_max,sol_min,heat(4),tempy1,tempy2,tempy3,tempy4,tempy5,temp6,temp7

  double precision :: xl, xm, xr, phil, phim, phir, thetal, thetam, thetar, sigmastartemp,heats(3),x0,y0,kappa,area,kappa_g,area_g
  logical :: on_axisx=.false.,on_axisy=.false.
  double precision :: Tk
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

  

  dist = 0.0
  sigmastartemp = 0.0
  sol_max=-1d9
  sol_min= 1d9
  s1=0d0
  s2=0d0
  s3=0d0
  s4=0d0
  s5=0d0
  s6 = 0d0
  s7 = 0d0
  sy1=0d0
  sy2=0d0
  sy3=0d0
  sy4=0d0
  kappa=0d0
  area = 0d0


  do lb=1,lnblocks

  
  
     if(nodetype(lb).eq.1) then
        dx = bsize(1,lb)/real(nxb)
        on_axisx=.false.

        
        if(bnd_box(1,2,lb).lt.1d-15)then
           if(ndim.eq.2.or.bnd_box(1,3,lb).lt.1d-15)then
              on_axisx=.true.
           end if
        end if




        if(on_axisx)then
           call phase_loc_rad(lb, s1,s2,s3,s4,s5,s6,s7)
        end if
        

! this code is to calculate position on the y axis (we do not need radius)
! for when we have hexagonal facets. Then the ratio x0/y0 is constant for a pure hexagon
! but varies in dendritic growth
        if(hex_facet.or.g_ellipse)then
           on_axisy=.false.
           if(bnd_box(1,1,lb).lt.1d-15)then
              if(ndim.eq.2.or.bnd_box(1,3,lb).lt.1d-15)then
                 on_axisy=.true.
              end if
           end if
           if(on_axisy)then
              call phase_loc_rady(lb, sy1,sy2,sy3,sy4)
           end if
        endif
        




        if(.not.has_children(lb).and.bnd_box(1,1,lb).lt.1d-15) then
           do i=il_bnd+1,iu_bnd-1
              if(unk(1,i,2,1,lb).gt.1d-3.and. unk(1,i,2,1,lb).lt.1d0-1d-3)then
                 sol_max=max(unk(1+nunkvbles,i,2,1,lb),sol_max)
                 sol_min=min(unk(1+nunkvbles,i,2,1,lb),sol_min)
              endif
           enddo
        endif

!        this section computes the normalised length of the dendrite surface. 2 * int(int phi^2(1-phi)^2       
!                                                                               ----------------------------        
!                                                                               Pi *    int x * phi^2(1-phi)^2
        if(.not.has_children(lb)) then
           do i=il_bnd+1,iu_bnd-1
              do j=jl_bnd+1,ju_bnd-1
                 phi = unk(1,i,j,1,lb)
                 phi_x = 0.5*(unk(1,i+1,j,1,lb) - unk(1,i-1,j,1,lb))
                 phi_y = 0.5*(unk(1,i,j+1,1,lb) - unk(1,i,j-1,1,lb))
                 kappa = kappa + sqrt(phi_x**2+phi_y**2)*dx
                 area = area + (1d0-phi)*dx*dx
!!$                 sum_heat(1)=sum_heat(1) + unk(1+18,i,j,1,lb)*dx*dx
!!$                 sum_heat(2)=sum_heat(2) + unk(2+18,i,j,1,lb)*dx*dx+ unk(3+18,i,j,1,lb)*dx*dx
!!$                 sum_heat(3)=sum_heat(4) + unk(4+18,i,j,1,lb)*dx*dx
              end do
           end do
        endif
        
     endif
  end do

  call MPI_REDUCE(kappa,kappa_g,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(area,area_g,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)


  !  s4 phi''(y)
  !  s3 x  
  !  s2 phi'(x)
  !  s1 normaliser phi^2*(1-phi)^2
if(.not.g_pure)then
  call MPI_REDUCE(sol_max,solute_max,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(sol_min,solute_min,1,amr_mpi_real,MPI_MIN,0,MPI_COMM_WORLD,ierr)
endif
!!$  call MPI_REDUCE(sum_heat(1),heats(1),1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!$  call MPI_REDUCE(sum_heat(2),heats(2),1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!$  call MPI_REDUCE(sum_heat(3),heats(3),1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)




  call MPI_ALLREDUCE(s1,temp1,1,amr_mpi_real,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(s2,temp2,1,amr_mpi_real,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(s3,temp3,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(s4,temp4,1,amr_mpi_real,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(s5,temp5,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(s6,temp6,1,amr_mpi_real,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(s7,temp7,1,amr_mpi_real,MPI_SUM,MPI_COMM_WORLD,ierr)


!  call MPI_AllReduce(temp6, g_vel, 1, amr_mpi_real, MPI_MAX, MPI_COMM_WORLD, ierr)


  g_kurv = temp4/temp2
  g_vel = temp6
  g_T_i = temp7

  g_ratio0 = 1d0/(sqrt(8d0)*g_lambda*temp1)


  if(hex_facet.or. g_ellipse)then
     call MPI_REDUCE(sy4,tempy4,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call MPI_REDUCE(sy2,tempy2,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call MPI_REDUCE(sy1,tempy1,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call MPI_REDUCE(sy3,tempy3,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  endif

  if(g_1D)temp3=0d0
  sum_heat=0d0
  heats=0d0


!  call MPI_REDUCE(temp_min,temp4,1,amr_mpi_real,MPI_MIN,0,MPI_COMM_WORLD,ierr)
  if(pid.eq.0)write (6,'(A,A)') 27, "[36;1;47m"
!  if(pid.eq.0)print *,"Tip Radius",temp3
!  if(pid.eq.0)write(*,*)temp3,temp4,temp3/temp4
!  if(pid.eq.0)write(unit=125,fmt=*) time, temp2, temp3

!  call MPI_REDUCE(dist,temp,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
!  call MPI_REDUCE(new_loc,temp2,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)



  !  s1 normaliser phi^2*(1-phi)^2
  !  s2 phi'(x)
  !  s3 x 
  !  s4 phi''(y)

 
!  if((solute_max- solute_min)/(g_c_max-g_c_min).lt.1.0)then
!     if(time.lt.20)g_alpha = 2d0*g_alpha0
!  else
!     g_alpha = g_alpha0
!  endif



  if(pid.eq.0.and..not.g_pure)write(*,'(6e12.5)')solute_max,solute_min,solute_min/solute_max  
  if(hex_facet)then
!     g_interface=temp1/temp2*0.8
     
     if(pid.eq.0)then
        print *,"Tip position, radius,  delta, curvature on facet, surface length"
!        write(*,'(5F14.5)') temp3/temp1,temp2/temp4,temp3/temp1*tempy1/tempy3,tempy4/tempy2,kappa_g/sqrt(area_g)*0.5641895835
        write(*,'(5F14.5)') temp3,temp2/temp4,1d0/temp1,tempy4/tempy2,kappa_g/sqrt(area_g)*0.5641895835

     endif
  else
     if(pid.eq.0)write(*,*)"      x             rho           delta         delta/delta0"
     if(pid.eq.0) write(*,'(4F14.6)')temp3,temp2/temp4,  1/temp1, g_ratio0
     if(.not.g_pure)then
        if(pid.eq.0)write(*,*)"      Delta c       k_E           V             Tip T"
        if(pid.eq.0) write(*,'(4F14.6)')(solute_max- solute_min)/(g_c_max-g_c_min),solute_min/solute_max*g_c_max/g_c_min,temp6,g_T0*(1d0+temp7)
     else
        if(pid.eq.0)write(*,*)"     Tip T"
        if(pid.eq.0) write(*,'(1F14.6)') g_T0*(1d0+temp7)
     endif
  endif
  if(pid.eq.0)print *,"dt",dt, g_alpha
  if(hex_facet)then
! time, position, radius, min solute, max solute , solute ratio, length/sqrt(area)* normalising factor 
!     
     if(pid.eq.0.and.g_needle)then
        write(126,'(7e14.5)') time, temp3,1d0/temp4,tempy4/tempy2,solute_min,solute_max,temp3/tempy3
     else
        if(pid.eq.0)write(126,'(7e14.5)') time, temp3,temp2/temp4,tempy4/tempy2,solute_min,solute_max,kappa_g/sqrt(area_g)*0.5641895835
     endif
  elseif(g_pure)then
     if(pid.eq.0)write(126,'(7e14.5)') time, temp3,temp2/temp4,g_T0*(1d0+temp7),g_T0*(1d0+temp7),g_T0*(1d0+temp7),g_T0*(1d0+temp7)
  elseif(g_needle)then
     if(pid.eq.0)write(126,'(7e14.5)') time, temp3,temp2/temp4,solute_min,solute_max,(solute_max- solute_min)/(g_c_max-g_c_min),temp3/tempy3
  else
     if(pid.eq.0)write(126,'(7e14.5)') time, temp3,temp2/temp4,solute_min,solute_max,(solute_max- solute_min)/(g_c_max-g_c_min),1d0/temp7
  endif

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

subroutine  phase_loc_rad(lb, s1,s2,s3,s4,s5,s6,s7)
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use time_dep_parameters

  implicit none
  double precision, intent(inout):: s1,s2,s3,s4,s5,s6,s7
  integer i, j, k, lb
  double precision cur_loc, cur_rad, x1a, x1b, x2a, x2b, a5, b1, b3, b5, dx
  double precision x0,y0, x1, x2,x3, xm1,phi,phi_x,phi_yy,phi_xx,radius,y1,xloc,L1,L2,L3
  double precision dw_phi,potential
  double precision rf4,rf5,rf6,rfactor
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
!         phi_yy = -(4d0*unk(1,i,j,k,lb)+unk(1,i,j+1,k,lb)-5d0*unk(1,i,j+2,k,lb))/(dx*dx)/14d0

!         phi_yy = -1d0/42d0*(5d0*unk(1,i,j,k,lb)+3d0*unk(1,i,j+1,k,lb)-unk(1,i,j+2,k,lb)-7d0*unk(1,i,j+3,k,lb))/(dx*dx)

         phi_xx = (unk(1,i-1,j,k,lb)-2d0*unk(1,i,j,k,lb)+unk(1,i+1,j,k,lb))/(dx*dx)

!         phi_yy = (-15d0*unk(1,i,j,k,lb)+16d0*unk(1,i,j+1,k,lb)-unk(1,i,j+2,k,lb))/(6d0*dx*dx)
!         endif
!         b5 =  dx*atanh(unk(1,i,j,k,lb))/(atanh(unk(1,i,j,k,lb))-atanh(unk(1,i+1,j,k,lb)))
!         b5 =  dx*(unk(1,i,j,k,lb)-0.5)/((unk(1,i,j,k,lb))-(unk(1,i+1,j,k,lb)))
         phi=unk(1,i,j,k,lb)
         dw_phi=2d0*phi*(phi-1d0)*(2d0*phi-1d0)
!         dw_phi=abs(phi_x)


         s2 = s2+phi_x * dx
!         s2 = s2+phi_x * dx
!         s2 = s2 + potential(phi,0d0,0)*phi_x*dx
!         s2 = s2+phi_x * dx * exp(-phi*4)
         s3 = s3+(1d0-phi)*dx
!         s1 = s1+dw_phi*dx
! here using K *d^2 =  int (x-R)^2 phi_x dx = int (x^2 - 2 R x + R^2) phi_ x dx
         if(phi_xx.gt.0d0)then
            s1 = s1 + phi_xx * dx
         else
            s5 = s5 - phi_xx * dx
         endif
         
!         s4 = s4+potential(phi,0d0,0)*phi_yy*dx
         s4 = s4+phi_yy*dx



         rfactor = dt/dtold
         rf4 = rfactor*rfactor/(rfactor+1.0)
         rf5 = (rfactor+1.0)
         rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
         s6 = s6 - dx*((rf4*unk(5,i,j,1,lb)-&
              rf5*unk(4,i,j,1,lb)+&
              rf6*unk(1,i,j,1,lb))/dt)
         
         s7 = s7 + phi_x*unk(1+nunkvbles,i,j,1,lb)
!     endif

  enddo
  
!  cur_loc=x0+b5

!  cur_rad = abs(phi_x/phi_yy)+b5  
!  cur_rad=abs(s2/s1)
!  cur_loc=s3/s1

end subroutine phase_loc_rad
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  phase_loc_rady(lb, s1,s2,s3,s4)
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

  i = il_bnd+nguard

  k = kl_bnd+nguard*k3d
  cur_loc=1.0

  do j=jl_bnd+nguard,ju_bnd-nguard
         x0 = bnd_box(1,2,lb)+(j-1.5)*dx !=y0
         phi_x  = 0.5*(unk(1,i,j+1,k,lb)-unk(1,i,j-1,k,lb))/dx !phi'(y)
         phi_yy = (-unk(1,i,j,k,lb)+unk(1,i+1,j,k,lb))/(dx*dx) !phi''(x)
         phi=unk(1,i,j,k,lb)
!         dw_phi=phi**2*(1d0-phi)**2
         dw_phi= abs(phi_x)
         s1=s1+dw_phi*dx
         s3=s3+(1d0-phi)*dx
         s2=s2+dw_phi*phi_x
         s4=s4+dw_phi*phi_yy
  enddo


end subroutine phase_loc_rady
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
  integer i, j, k, v, lb, fnum,ii,jj
  character*75 fname,gname
  character (len=2)  :: proc_string
  character (len=2)  :: grid_string
  character (len=7)  :: fnum_string
  double precision dx, gphi(2), AnisEnergy,vbles(N_unk,3,3)
  double precision x0,y0,p0,n0,xposition, FreeEnergy,potential,Gibbs_FE_e,Gibbs_FE_l
  double precision ::c,T,phi,rfactor,rf4,rf5,rf6,K_heat(4)
  integer, dimension(16) :: Hdum= (/ 1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1,  0, -1 , 1, -1 /) 
  double precision, dimension(2) :: dH
  integer, dimension(2,8) :: H
  double precision :: pcb_LagrangePoly,xx,yy,pcb_interpolate,unkM(2,10,10),unk2(2,10,10)
  integer ::summ=0,npoints=4
 
  H = RESHAPE(Hdum,shape=(/2,8/))
     
   Write (proc_string, '(i2.2)') pid
   Write (grid_string, '(i2.2)') lrefine_max
   Write (fnum_string, '(i7.7)') fnum
    rfactor = dt/dtold
   rf4 = rfactor*rfactor/(rfactor+1.0)
   rf5 = (rfactor+1.0)
   rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)


   do lb=1,lnblocks
      if(nodetype(lb).eq.1) then
         dx = bsize(1,lb)/real(nxb)
         
         do j=jl_bnd+1,ju_bnd-1
            do i=il_bnd+1,iu_bnd-1
               do v = 1, n_unk*nunkvbles
                  if (ABS(unk(v,i,j,1,lb)).lt.1.0e-20) unk(v,i,j,1,lb)=0.0
               end do
            end do
         end do
      end if
   end do
   

  
  fname = trim(output_dir)// '2dcontour2_' // grid_string // "_" // proc_string // "_" // fnum_string
!   gname = trim(output_dir)// 'bulk2_' // grid_string // "_" // proc_string // "_"
   gname = trim(output_dir)// 'bulk_' // grid_string // "_" // proc_string // "_"  // fnum_string




   open (unit=150,file=gname)

!   if(g_linear)then
      do lb=1,lnblocks
         if(.not.has_children(lb)) then
            dx = bsize(1,lb)/real(nxb)
            do j=jl_bnd+1,ju_bnd-1
               do i=il_bnd+1,iu_bnd-1
                  x0 =  bnd_box(1,1,lb)+(i-1.5)*dx
                  y0 =  bnd_box(1,2,lb)+(j-1.5)*dx
                  phi = Unk(1,i,j,1,lb)
                  c   = Unk(2+nunkvbles,i,j,1,lb)
                  write(150,'(3(E16.5,","),E16.5)')x0/4d1,y0/4d1,1-phi,c*10
               end do
            end do
         endif
      end do
      close(150)
      open (unit=150, file="One_dimensional.txt"//proc_string)
      do lb=1,lnblocks
         !         if(.not.has_children(lb)) then
         
         if(nodetype(lb).eq.1) then
            if(bnd_box(1,2,lb).lt.1d-15)then 
               dx = bsize(1,lb)/real(nxb)
               j = jl_bnd+nguard
               do i=il_bnd+1,iu_bnd-1
                  x0 =  bnd_box(1,1,lb)+(i-1.5)*dx
                  phi = Unk(1,i,j,1,lb)
                  c   = Unk(2+nunkvbles,i,j,1,lb)
                  write(150,'(3E16.5)')x0,1-phi,c
                  
               end do
            endif
         endif
      end do
      close(150)

!   else
      
!!$      do lb=1,lnblocks
!!$         if(.not.has_children(lb)) then
!!$
!!$            dx = bsize(1,lb)/real(nxb)
!!$            do i=1,g_order
!!$               do j=1,g_order
!!$                  unkM(1,i,j)=unk(1,i,j,1,lb)
!!$                  unkM(2,i,j)=unk(1+nunkvbles,i,j,1,lb)
!!$               enddo
!!$            enddo
!!$            
!!$            if(g_npoints.le.4)then
!!$               do j=2,nxb+1
!!$                  y0 = bnd_box(1,2,lb)+(0.5+g_y10(j)/(1.+g_y10(g_order-1)))*nxb*dx
!!$                  do i=2,nxb+1
!!$                     x0 = bnd_box(1,1,lb)+(0.5+g_y10(i)/(1.+g_y10(g_order-1)))*nxb*dx
!!$                     phi = Unk(1,i,j,1,lb)
!!$                     c   = Unk(1+nunkvbles,i,j,1,lb)
!!$                     write(150,'(3(E16.5,","),E16.5)')x0/4d1,y0/4d1,1-phi,c*10
!!$                  end do
!!$               end do
!!$ 
!!$            else
!!$               do j=0,g_npoints
!!$                  yy=(1.+g_y10(g_order-1))/g_npoints*j-1. +0.5*(1.-g_y10(g_order-1))
!!$                  y0 =  bnd_box(1,2,lb)+(0.5+yy/(1.+g_y10(g_order-1)))*nxb*dx
!!$                  do i=1,g_npoints
!!$                     xx=(1.+g_y10(g_order-1))/g_npoints*i-1. +0.5*(1.-g_y10(g_order-1))
!!$                     x0 =  bnd_box(1,1,lb)+(0.5+xx/(1.+g_y10(g_order-1)))*nxb*dx
!!$                     phi = pcb_interpolate(nxb+2,xx,yy,UnkM(1,:,:))
!!$                     c   = pcb_interpolate(nxb+2,xx,yy,UnkM(2,:,:))
!!$                     write(150,'(3(E16.5,","),E16.5)')x0/4d1,y0/4d1,1-phi,c*10
!!$                  end do
!!$               end do
!!$
!!$            endif
!!$         endif
!!$      end do
 !  endif



end subroutine pf_output_2d

double precision function R2Z(x)
double precision, intent (in)::x
R2Z=x
if(abs(x).lt.1e-90)R2Z=0E0
end  function R2Z

