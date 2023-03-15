!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_initial_soln_blk                                         REQUIRED
!!!!  * Sets the initial solution on a block
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
  integer :: nguard0, i, j, k,ip,ii
  double precision :: radius, pi, dx, dy, dz, xi, yi, zi,u1
  double precision c_ave,ci,Gibbs_FE_e,d0




  if (ndim.eq.2) then
     write(*,*)"only 3D. Set ndim =3"
     stop
  endif


  c_ave=g_c0


  nguard0 = nguard*npgs
  unk(:,:,:,:,lb) =0.

! set values for unk
  dx = bsize(1,lb)/real(nxb)
  do k=kl_bnd+nguard0*k3d,ku_bnd-nguard0*k3d
     zi = bnd_box(1,3,lb) + dx*(real(k-nguard0)-.5) 
     do j=jl_bnd+nguard0*k2d,ju_bnd-nguard0*k2d
        yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5)
        do i=il_bnd+nguard0,iu_bnd-nguard0
           xi = bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5)
           do ii=1,nunkvbles-1
              !phase (liquid)
!              unk(ii,i,j,k,lb)=0.5 +0.5*tanh((sqrt(xi**2+yi**2+zi**2)-nuc_radius)/d0)
              unk(ii,i,j,k,lb)= 1d0/(1.+exp(-(sqrt(xi**2+yi**2+zi**2)-nuc_radius)/(g_lambda*0.5)))
                ! Temperature            
              unk(ii+2*nunkvbles,i,j,k,lb) =delta
              ! solute
              unk(ii+nunkvbles,i,j,k,lb) =c_ave
           enddo
        enddo
     enddo
  enddo


  
end subroutine app_initial_soln_blk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_initial_soln_end                                   
  return
end subroutine app_initial_soln_end


subroutine app_test_refinement
  return
end subroutine app_test_refinement


