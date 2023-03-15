!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_initial_soln_blk                                         REQUIRED
!!!!  * Sets the initial solution on a block
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_initial_soln_end                                         REQUIRED
!!!!  * Called at the end of the initial value setting
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_test_refinement                                          REQUIRED
!!!!  * Application specific refinement tests - unused in phase field code
!!!!  * Only called if local_adaptation=2
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Latest version: Chris Goodyer, May 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_initial_soln_blk(lb)
  
!--------------------------------------------------------------
! include files 
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use solution_parameters
  implicit none
  include 'mpif.h'
  
  integer, intent(in) :: lb
  integer :: nguard0, i, j, k
  double precision :: radius, pi, dx, dy, dz, xi, yi, zi
!--------------------------------------------------------------

  pi=4.0*datan(1.0)
  nguard0 = nguard*npgs
! set values for unk
  dx = bsize(1,lb)/real(nxb)
  do k=kl_bnd+nguard0*k3d,ku_bnd-nguard0*k3d
     zi = bnd_box(1,3,lb) + dx*(real(k-nguard0)-.5) - nucleate_z
     if (ndim.eq.2) zi = 0.0
     do j=jl_bnd+nguard0*k2d,ju_bnd-nguard0*k2d
        yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5) - nucleate_y
        do i=il_bnd+nguard0,iu_bnd-nguard0
           unk(1,i,j,k,lb) = 0.0
           xi = bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5) - nucleate_x
           
           ! Phase
           unk(1,i,j,k,lb) = -tanh(0.6*((xi*xi+yi*yi+zi*zi)**0.5-nuc_radius))
           unk(2,i,j,k,lb) = 0.0
           unk(3,i,j,k,lb) = 0.0
           unk(4,i,j,k,lb) = -tanh(0.6*((xi*xi+yi*yi+zi*zi)**0.5-nuc_radius))
           unk(5,i,j,k,lb) = -tanh(0.6*((xi*xi+yi*yi+zi*zi)**0.5-nuc_radius))

           ! Temperature
           unk(1+nunkvbles,i,j,k,lb) = delta+abs(delta)*0.5*(unk(1,i,j,k,lb)+1.0)
           !unk(1+nunkvbles,i,j,k,lb) = delta+abs(delta)*0.5*(-tanh(0.05*((xi*xi+yi*yi+zi*zi)-nuc_radius))+1.0)

           unk(2+nunkvbles,i,j,k,lb) = 0.0
           unk(3+nunkvbles,i,j,k,lb) = 0.0
           unk(4+nunkvbles,i,j,k,lb) = unk(1+nunkvbles,i,j,k,lb)
           unk(5+nunkvbles,i,j,k,lb) = unk(1+nunkvbles,i,j,k,lb)

           ! Solute concentration
           if(nvar.eq.3*nunkvbles) then
              ! Solute and thermal
              unk(1+2*nunkvbles,i,j,k,lb) = 0.0
              unk(2+2*nunkvbles,i,j,k,lb) = 0.0
              unk(3+2*nunkvbles,i,j,k,lb) = 0.0
              unk(4+2*nunkvbles,i,j,k,lb) = 0.0
              unk(5+2*nunkvbles,i,j,k,lb) = 0.0
           else if (solute) then
              ! Solute only - no thermal
              unk(1+nunkvbles,i,j,k,lb) = 0.0
              unk(2+nunkvbles,i,j,k,lb) = 0.0
              unk(3+nunkvbles,i,j,k,lb) = 0.0
              unk(4+nunkvbles,i,j,k,lb) = 0.0
              unk(5+nunkvbles,i,j,k,lb) = 0.0
           endif

        enddo
     enddo
  enddo

  return

end subroutine app_initial_soln_blk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_initial_soln_end                                   
  return
end subroutine app_initial_soln_end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_test_refinement
  return
end subroutine app_test_refinement
