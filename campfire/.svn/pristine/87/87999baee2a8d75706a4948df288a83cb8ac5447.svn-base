double precision function A_grad(phi,grad_phi)
use solution_parameters
implicit none
double precision, intent (in)::phi(3),grad_phi(3,2)
double precision :: r3(3,3,2),A3(3,3,2),r2(3,3)
integer :: i,j,a
do i=1,3
do j=1,3
do a=1,2
   r3(i,j,a) = phi(i)*grad_phi(j,a)-phi(j)*grad_phi(i,a)
enddo
enddo
enddo

do i=1,3
do j=1,3
   r2=r3(i,j,1)**2+r3(i,j,2)**2
enddo
enddo

do i=1,3
do j=1,3
do a=1,2
   A3(i,j,a) = (1d0+epsilon_tilde*(r3(i,j,1)**4/r2(i,j)**2 + r3(i,j,2)**4/r2(i,j)**2))*r3(i,j,a) 
enddo
enddo
enddo

A_grad = 0d0

do i=2,3
do j=1,i-1
do a=1,2
   A_grad = A_grad + 0.5*A3(i,j,a)**2 
enddo
enddo
enddo
end function A_grad
