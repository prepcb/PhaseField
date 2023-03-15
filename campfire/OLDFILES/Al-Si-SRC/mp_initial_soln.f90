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
  double precision c_ave,ci,Gibbs_FE_e,aa,bb,cc
  double precision :: FreeEnergy,c_dot(1)=0d0,phi_dot(1)=0d0,phi(1),phiX

!FreeEnergy(c,T,phi,c_dot,phi_dot,lp)

  c_ave=g_c0



  if(.false.)then

     open (unit=99,file="FE.txt")
     j=10000
     do i=1,j-1
        u1=(1.*i)/j
        phi(1)=0d0
        aa = FreeEnergy(u1,0d0,phi,c_dot,phi_dot,0)
        phi(1)=1d0
        bb = FreeEnergy(u1,0d0,phi,c_dot,phi_dot,0)
        write(99,*)aa,bb
     enddo
     close(99)
     
     write(*,*)"done in FE.txt"
     
     stop
  endif
  


  
  nguard0 = nguard*npgs
  unk(:,:,:,:,lb) =0.
! set values for unk
  dx = bsize(1,lb)/real(nxb)
  do k=kl_bnd+nguard0*k3d,ku_bnd-nguard0*k3d
     zi = bnd_box(1,3,lb) + dx*(real(k-nguard0)-.5) 
     if (ndim.eq.2) zi = 0.0
     do j=jl_bnd+nguard0*k2d,ju_bnd-nguard0*k2d
        yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5)
        do i=il_bnd+nguard0,iu_bnd-nguard0
           xi = bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5)
           do ii=1,nunkvbles-1
!phase (liquid)
              if(g_1D)then
                 unk(ii,i,j,k,lb)=0.5 +0.5*tanh(xi-nuc_radius)
              else
                 unk(ii,i,j,k,lb)=0.5 +0.5*tanh((xi**2+yi**2)**0.5-nuc_radius)
              endif
! Temperature            
              unk(ii+nunkvbles,i,j,k,lb) =delta
! solute
              unk(ii+2*nunkvbles,i,j,k,lb) =(1d0-unk(ii,i,j,k,lb))*0.99+unk(ii,i,j,k,lb)*c_ave


           enddo
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


subroutine app_test_refinement
  return
end subroutine app_test_refinement

