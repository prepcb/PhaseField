double precision function A_grad(grad_phi)
use solution_parameters
implicit none
double precision, intent (in)::grad_phi(2)
double precision :: A1(2),r,XX,YY
integer :: a





  XX=grad_phi(1)*grad_phi(1)
  YY=grad_phi(2)*grad_phi(2)
  
  r=XX+YY
  
  
  
  do a=1,2
     if(r.gt.1d-20)then
        A1(a) = (1d0+epsilon_tilde*((XX/r)**2 + (YY/r)**2))*grad_phi(a) 
     else
        A1(a) = (1d0+epsilon_tilde)*grad_phi(a)
     end if
  enddo
  
  A_grad = eps22(1,2)*0.5*(A1(2)**2 +A1(2)**2 )
  
  
end function A_grad
DOUBLE PRECISION FUNCTION A_GRAD_D(grad_phi, grad_phid)
  USE SOLUTION_PARAMETERS
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: grad_phi(2)
  DOUBLE PRECISION, INTENT(IN) :: grad_phid(2)
  DOUBLE PRECISION :: a1(2), r, xx, yy
  DOUBLE PRECISION :: a1d(2), rd, xxd, yyd
  INTEGER :: a


  xxd = grad_phid(1)*grad_phi(1) + grad_phi(1)*grad_phid(1)
  xx = grad_phi(1)*grad_phi(1)
  yyd = grad_phid(2)*grad_phi(2) + grad_phi(2)*grad_phid(2)
  yy = grad_phi(2)*grad_phi(2)
  rd = xxd + yyd
  r = xx + yy
  a1d = 0.D0
  DO a=1,2
    IF (r .GT. 1d-20) THEN
      a1d(a) = epsilon_tilde*(2*xx*(xxd*r-xx*rd)/r**3+2*yy*(yyd*r-yy*rd)&
&       /r**3)*grad_phi(a) + (1d0+epsilon_tilde*((xx/r)**2+(yy/r)**2))*&
&       grad_phid(a)
      a1(a) = (1d0+epsilon_tilde*((xx/r)**2+(yy/r)**2))*grad_phi(a)
    ELSE
      a1d(a) = (1d0+epsilon_tilde)*grad_phid(a)
      a1(a) = (1d0+epsilon_tilde)*grad_phi(a)
    END IF
  END DO

  a_grad_d = eps22(1,2)*2*a1(2)*a1d(2)

END FUNCTION A_GRAD_D
