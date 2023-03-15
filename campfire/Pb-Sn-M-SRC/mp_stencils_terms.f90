!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  Utilities for phase field stencils in 2-d and 3-d
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



double precision function EntropyMixing(c,d)
double precision, intent (in)::c
integer, intent (in)::d
double precision:: cc,a=1d-6
b = 1.-a

!!
!!$if(c.gt.1d0-1d-16.or.c.lt.1d-16)then
!!$   write(*,*)"c=",c, "in EntropyMixing"
!!$   stop
!!$
!!$endif
!quadratic extension
if(c.lt.a)then
   if(d.eq.0) EntropyMixing = 0.5*(1d0/a + 1d0/(1d0-a))*(c-a)**2 + log(a/(1d0-a))*(c-a)+a*log(a)+(1.-a)*log(1.-a)
   if(d.eq.1) EntropyMixing =    (1d0/a + 1d0/(1d0-a))*(c-a) + log(a/(1d0-a))
   if(d.eq.2) EntropyMixing =     (1d0/a + 1d0/(1d0-a))

else if(c.gt.b)then
   if(d.eq.0) EntropyMixing = 0.5*(1d0/b + 1d0/(1d0-b))*(c-b)**2 + log(b/(1d0-b))*(c-b)+b*log(b)+(1.-b)*log(1.-b)
   if(d.eq.1) EntropyMixing =     (1d0/b + 1d0/(1d0-b))*(c-b) + log(b/(1d0-b))
   if(d.eq.2) EntropyMixing =     (1d0/b + 1d0/(1d0-b))
else
   if(d.eq.1)then
      EntropyMixing = log(c/(1d0-c))
   elseif(d.eq.2)then
      EntropyMixing = 1d0/c + 1d0/(1d0-c)
   else
      EntropyMixing = c*log(c)+(1.-c)*log(1.-c)
   endif
endif
end function EntropyMixing

 
subroutine get_y1_y2(y,d,y1_y2)
integer, intent(in)::d
double precision, intent(in) :: y(2)
double precision, intent(out) :: y1_y2(6)



y1_y2=0d0
if(d.eq.0)then
   y1_y2(1)=y(1)*y(2)
   y1_y2(2)=(1d0-y(1))*y(2)
   y1_y2(3)=(1d0-y(2))*y(1)
   y1_y2(4)=(1d0-y(1))*(1d0-y(2))
   y1_y2(5)=y(1)*(1d0-y(1))
   y1_y2(6)=y(2)*(1d0-y(2))
elseif(d.eq.1)then
   y1_y2(1)=y(2)
   y1_y2(2)=-y(2)
   y1_y2(3)=(1d0-y(2))
   y1_y2(4)=-(1d0-y(2))
   y1_y2(5)=1d0-2d0*y(1)
   y1_y2(6)=0d0
elseif(d.eq.2)then
   y1_y2(1)=y(1)
   y1_y2(2)=(1d0-y(1))
   y1_y2(3)=-y(1)
   y1_y2(4)=-(1d0-y(1))
   y1_y2(5)=0d0
   y1_y2(6)=1d0-2d0*y(2)
else
   write(*,*)"d=0 to 2"
   stop
endif
end subroutine get_y1_y2


!!$or i = 1:N, 
!!$  if (c(i) > 0.27246)
!!$    y1(i) = NaN ;
!!$  elseif (c(i) > 0.2350)
!!$    y1(i) = (0.2350/c(i) - 0.8625)/0.1375 ;
!!$  else
!!$    y1(i) = NaN ;
!!$  end
!!$end

double precision function AlFey1(c,d)
integer, intent(in)::d
double precision, intent(in) :: c
if(c.gt.0.27246.or.c.lt.0.2350)then
   write(*,*)"c out of bound in AlFey1()"
   stop
endif
if(d.eq.0)then
      AlFey1=(0.2350/c - 0.8625)/0.1375 
      return
elseif(d.eq.1)then
      AlFey1=-(0.2350/c**2)/0.1375 
      return
else
   write(*,*)"d =0 or 1"
endif



end function AlFey1


double precision function y1(c,d)
integer, intent(in)::d
double precision, intent(in) :: c
if(c.gt.1d0-1d-9.or.c.lt.1d-9)then
   write(*,*)"c out of bound in y1()"
   stop
endif
if(d.eq.0)then
   if(c.ge.0.4d0)then
      y1=1d0-1d-5
      return
   else
      y1=2.5*c
      return
   endif
elseif(d.eq.1)then
   if(c.ge.0.4d0)then
      y1=1d-5
      return
   else
      y1=2.5
      return
!   else
!      write(*,*)"d y1/dc not determined at c=",c
!      stop
   endif
else
   write(*,*)"d =0 or 1"
endif
end function y1

double precision function y2(c,d)
integer, intent(in)::d
double precision, intent(in) :: c
if(c.gt.1d0-1d-9.or.c.lt.1d-9)then
   write(*,*)"c out of bound in y2()"
   stop
endif

if(d.eq.0)then
   if(c.ge.0.5d0)then
      y2=1d-5
      return
   elseif(c.gt.0.4d0)then
      y2=1d0-(2d0-5d0*c)/(c-1d0)
      return
   else
      y2=1d0-1d-5
      return
   endif
elseif(d.eq.1)then
!   if(c.eq.0.5.or.c.eq.0.4)then
!       write(*,*)"d y2/dc not determined at c =",c
!       stop
   if(c.gt.0.5d0)then
      y2=1d-5
      return
   elseif(c.gt.0.4d0)then
      y2=-3d0/(1d0-c)**2
      return
   else
      y2=1d-5
      return
   endif
else
   write(*,*)"d =0 or 1"
endif
end function y2





double precision function y3(c,d)
integer, intent(in)::d
double precision, intent(in) :: c
if(c.gt.1d0-1d-9.or.c.lt.1d-9)then
   write(*,*)"c out of bound in y3()"
   stop
endif
if(d.eq.0)then
   if(c.gt.0.5d0)then
      y3=2d0*(1d0-c)
      return
   else
      y3=1d0-1d-5
      return
   endif
elseif(d.eq.1)then
   if(c.gt.0.5d0)then
      y3=-2d0
      return
   else
      y3=0d0
      return
!   else
!      write(*,*)"d y3/dc not determined at c=",c
!      stop
   endif
else
   write(*,*)"d =0 or 1"
endif
end function y3


double precision function y4(c,d)
integer, intent(in)::d
double precision, intent(in) :: c
if(c.gt.1d0-1d-9.or.c.lt.1d-9)then
   write(*,*)"c out of bound in y4()"
   stop
endif
if(d.eq.0)then
   if(c.ge.0.5d0)then
      y4=1d0-1d-5
      return
   else
      y4=c/(1d0-c)
      return
   endif
elseif(d.eq.1)then
   if(c.ge.0.5d0)then
      y4=0d0
      return
   else
      y4=1d0/(1d0-c)**2
      return
!   else
!      write(*,*)"d y4/dc not determined at c=",c
!      stop
   endif
else
   write(*,*)"d =0 or 1"
endif
end function y4



double precision function gibbs_T(h,T,d)
double precision, intent(in) :: T,h(8)
integer, intent(in) :: d
double precision V(8)
integer i
gibbs_T=0d0
call GibbsVector(V,T,d)
do i= 1,8
   gibbs_T=gibbs_T+h(i)*V(i)
enddo
end function gibbs_T


double precision function Gibbs_FE_l(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5)
  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gl_,Tk,Gibbs_T
  double precision :: hlA_(8),hlB_(8)


  Tk=T_k(T)
  if (Tk > 1728.0) then
     hlA_=AlNihlA3 
     hlB_=AlNihlB2 
  elseif (Tk > 933.47)then
     hlA_=AlNihlA3 
     hlB_=AlNihlB1 
  elseif (Tk > 700.0)then
     hlA_=AlNihlA2 
     hlB_=AlNihlB1 
  else
     hlA_=AlNihlA1 
     hlB_=AlNihlB1 
  endif
  if(d_c.eq.0.and.d_T.eq.0)then
     Tk=T_k(T)
     glA=gibbs_T(hlA_,Tk,0)
     glB=gibbs_T(hlB_,Tk,0)
     call MixingVector(RK,c,0)
     RKl_=         (AlNihrkl0(1)+AlNihrkl0(2)*Tk)*RK(1)
     RKl_=RKl_ + (AlNihrkl1(1)+AlNihrkl1(2)*Tk)*RK(2)
     RKl_=RKl_ + (AlNihrkl2(1)+AlNihrkl2(2)*Tk)*RK(3)
     RKl_=RKl_ + (AlNihrkl3(1)+AlNihrkl3(2)*Tk)*RK(4)
     RKl_=RKl_ + (AlNihrkl4(1)+AlNihrkl4(2)*Tk)*RK(5)
     ent_=g_R*Tk*EntropyMixing(c,0)
     gl_=glA*(1d0-c)+glB*c+ent_+RKl_
     Gibbs_FE_l=gl_
  elseif(d_c.eq.0.and.d_T.eq.1)then
     Tk=T_k(T)
     glA=gibbs_T(hlA_,Tk,1)
     glB=gibbs_T(hlB_,Tk,1)
     call MixingVector(RK,c,0)
     RKl_=         AlNihrkl0(2)*RK(1)
     RKl_=RKl_ +   AlNihrkl1(2)*RK(2)
     RKl_=RKl_ +   AlNihrkl2(2)*RK(3)
     RKl_=RKl_ +   AlNihrkl3(2)*RK(4)
     RKl_=RKl_ +   AlNihrkl4(2)*RK(5)
     ent_=g_R*EntropyMixing(c,0)
     gl_=glA*(1d0-c)+glB*c+ent_+RKl_
     Gibbs_FE_l=gl_
  elseif(d_c.eq.0.and.d_T.eq.2)then
     Tk=T_k(T)
     glA=gibbs_T(hlA_,Tk,2)
     glB=gibbs_T(hlB_,Tk,2)
     call MixingVector(RK,c,0)
     gl_=glA*(1d0-c)+glB*c
     Gibbs_FE_l=gl_
  elseif(d_c.eq.1.and.d_T.eq.0)then
     Tk=T_k(T)
     glA=gibbs_T(hlA_,Tk,0)
     glB=gibbs_T(hlB_,Tk,0)
     call MixingVector(RK,c,1)
     RKl_=         (AlNihrkl0(1)+AlNihrkl0(2)*Tk)*RK(1)
     RKl_=RKl_ + (AlNihrkl1(1)+AlNihrkl1(2)*Tk)*RK(2)
     RKl_=RKl_ + (AlNihrkl2(1)+AlNihrkl2(2)*Tk)*RK(3)
     RKl_=RKl_ + (AlNihrkl3(1)+AlNihrkl3(2)*Tk)*RK(4)
     RKl_=RKl_ + (AlNihrkl4(1)+AlNihrkl4(2)*Tk)*RK(5)
     ent_=g_R*Tk*EntropyMixing(c,1)
     gl_=-glA+glB+ent_+RKl_
     Gibbs_FE_l=gl_
  elseif(d_c.eq.1.and.d_T.eq.1)then
     Tk=T_k(T)
     glA=gibbs_T(hlA_,Tk,1)
     glB=gibbs_T(hlB_,Tk,1)
     call MixingVector(RK,c,0)
     RKl_=         AlNihrkl0(2)*RK(1)
     RKl_=RKl_ +   AlNihrkl1(2)*RK(2)
     RKl_=RKl_ +   AlNihrkl2(2)*RK(3)
     RKl_=RKl_ +   AlNihrkl3(2)*RK(4)
     RKl_=RKl_ +   AlNihrkl4(2)*RK(5)
     ent_=g_R*EntropyMixing(c,1)
     gl_=-glA+glB+ent_+RKl_
     Gibbs_FE_l=gl_
  elseif(d_c.eq.1.and.d_T.eq.2)then
     Tk=T_k(T)
     glA=gibbs_T(hlA_,Tk,2)
     glB=gibbs_T(hlB_,Tk,2)
     call MixingVector(RK,c,0)
     gl_=-glA+glB
     Gibbs_FE_l=gl_
  elseif(d_c.eq.2.and.d_T.eq.0)then
     Tk=T_k(T)
     glA=gibbs_T(hlA_,Tk,0)
     glB=gibbs_T(hlB_,Tk,0)
     call MixingVector(RK,c,2)
     RKl_=         (AlNihrkl0(1)+AlNihrkl0(2)*Tk)*RK(1)
     RKl_=RKl_ + (AlNihrkl1(1)+AlNihrkl1(2)*Tk)*RK(2)
     RKl_=RKl_ + (AlNihrkl2(1)+AlNihrkl2(2)*Tk)*RK(3)
     RKl_=RKl_ + (AlNihrkl3(1)+AlNihrkl3(2)*Tk)*RK(4)
     RKl_=RKl_ + (AlNihrkl4(1)+AlNihrkl4(2)*Tk)*RK(5)
     ent_=g_R*Tk*EntropyMixing(c,2)
     gl_=-glA+glB+ent_+RKl_
     Gibbs_FE_l=gl_
  else
     write(*,*)"unexpected d_c,d_T combination",d_c,d_T
     stop
  endif

  end function Gibbs_FE_l


double precision function Gibbs_FE_e(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in)::d_c,d_T
  double precision :: SLe(6),yy(2),yy_d(2),y1_y2(6),y1_y2_d1(6),y1_y2_d2(6)
  double precision :: y3,y4,EntropyMixing,T_k
  double precision ::entSL_,entSL_d,ge1_,ge1_d1,ge1_d2,ge2_,ge2_d1,ge2_d2,ge_,ge_d,Tk,Gibbs_T
  double precision :: AlNihlA_(8),AlNihlB_(8),AlNihaA_(8),AlNihaB_(8),AlNihslb1_(8),AlNihslb2_(8),AlNihslb3_(8),AlNihslb4_(8),AlNihc_(8),AlNihsle1_(8),AlNihsle2_(8),AlNihsle3_(8),AlNihsle4_(8)
  
 Tk=T_k(T)
  if (Tk > 1728.0) then
  AlNihlA_=AlNihlA3 
  AlNihlB_=AlNihlB2 
  AlNihaA_=AlNihaA3 
  AlNihaB_=AlNihaB2 
  AlNihslb1_=AlNihslb1_4 
  AlNihslb2_=AlNihslb2_3 
  AlNihslb3_=AlNihslb3_4 
  AlNihslb4_=AlNihslb4_4 
  AlNihc_=AlNihc4 
  AlNihsle1_=AlNihsle1_4 
  AlNihsle2_=AlNihsle2_3 
  AlNihsle3_=AlNihsle3_2 
  AlNihsle4_=AlNihsle4_2 
elseif (Tk > 933.47)then
  AlNihlA_=AlNihlA3 
  AlNihlB_=AlNihlB1 
  AlNihaA_=AlNihaA3 
  AlNihaB_=AlNihaB1 
  AlNihslb1_=AlNihslb1_3 
  AlNihslb2_=AlNihslb2_3 
  AlNihslb3_=AlNihslb3_3 
  AlNihslb4_=AlNihslb4_3 
  AlNihc_=AlNihc3 
  AlNihsle1_=AlNihsle1_3 
  AlNihsle2_=AlNihsle2_3 
  AlNihsle3_=AlNihsle3_1 
  AlNihsle4_=AlNihsle4_1 
elseif (Tk > 700.0)then
  AlNihlA_=AlNihlA2 
  AlNihlB_=AlNihlB1 
  AlNihaA_=AlNihaA2 
  AlNihaB_=AlNihaB1 
  AlNihslb1_=AlNihslb1_2 
  AlNihslb2_=AlNihslb2_2 
  AlNihslb3_=AlNihslb3_2 
  AlNihslb4_=AlNihslb4_2 
  AlNihc_=AlNihc2 
  AlNihsle1_=AlNihsle1_2 
  AlNihsle2_=AlNihsle2_2 
  AlNihsle3_=AlNihsle3_1 
  AlNihsle4_=AlNihsle4_1 
else
  AlNihlA_=AlNihlA1 
  AlNihlB_=AlNihlB1 
  AlNihaA_=AlNihaA1 
  AlNihaB_=AlNihaB1 
  AlNihslb1_=AlNihslb1_1 
  AlNihslb2_=AlNihslb2_1 
  AlNihslb3_=AlNihslb3_1 
  AlNihslb4_=AlNihslb4_1 
  AlNihc_=AlNihc1 
  AlNihsle1_=AlNihsle1_1 
  AlNihsle2_=AlNihsle2_1 
  AlNihsle3_=AlNihsle3_1 
  AlNihsle4_=AlNihsle4_1 
endif

  if(d_c.eq.0.and.d_T.eq.0)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(AlNihsle1_,Tk,0)
     SLe(2)=gibbs_T(AlNihsle2_,Tk,0)
     SLe(3)=gibbs_T(AlNihsle3_,Tk,0)
     SLe(4)=gibbs_T(AlNihsle4_,Tk,0)
     SLe(5)=AlNihsle5(1)+AlNihsle5(2)*Tk
     SLe(6)=AlNihsle6(1)+AlNihsle6(2)*Tk

     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     entSL_=g_R*Tk*(EntropyMixing(yy(1),0)+EntropyMixing(yy(2),0))
     call get_y1_y2(yy,0,y1_y2)
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge2_=SLe(5)*y1_y2(5)+SLe(6)*y1_y2(6)
     ge_=(ge1_+ge2_+entSL_)/(1d0+yy(1))
     Gibbs_FE_e=ge_
  elseif(d_c.eq.1.and.d_T.eq.0)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(AlNihsle1_,Tk,0)
     SLe(2)=gibbs_T(AlNihsle2_,Tk,0)
     SLe(3)=gibbs_T(AlNihsle3_,Tk,0)
     SLe(4)=gibbs_T(AlNihsle4_,Tk,0)
     SLe(5)=AlNihsle5(1)+AlNihsle5(2)*Tk
     SLe(6)=AlNihsle6(1)+AlNihsle6(2)*Tk

     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     yy_d(1)=y4(c,1)
     yy_d(2)=y3(c,1)
     entSL_=g_R*Tk*(EntropyMixing(yy(1),0)+EntropyMixing(yy(2),0))
     entSL_d=g_R*Tk*(EntropyMixing(yy(1),1)*yy_d(1)+EntropyMixing(yy(2),1)*yy_d(2))


     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1)
     call get_y1_y2(yy,2,y1_y2_d2)
     
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge2_=SLe(5)*y1_y2(5)+SLe(6)*y1_y2(6)
     ge1_d1=(SLe(1)*y1_y2_d1(1)+SLe(2)*y1_y2_d1(2)+SLe(3)*y1_y2_d1(3)+SLe(4)*y1_y2_d1(4))*yy_d(1)
     ge1_d2=(SLe(1)*y1_y2_d2(1)+SLe(2)*y1_y2_d2(2)+SLe(3)*y1_y2_d2(3)+SLe(4)*y1_y2_d2(4))*yy_d(2)
     ge2_d1=(SLe(5)*y1_y2_d1(5)+SLe(6)*y1_y2_d1(6))*yy_d(1)
     ge2_d2=(SLe(5)*y1_y2_d2(5)+SLe(6)*y1_y2_d2(6))*yy_d(2)
     ge_=(ge1_+ge2_+entSL_)/(1d0+yy(1))
     ge_d=-ge_/(1d0+yy(1))*yy_d(1)+(ge1_d1+ge1_d2+ge2_d1+ge2_d2+entSL_d)/(1d0+yy(1))
     Gibbs_FE_e=ge_d

  elseif(d_c.eq.1.and.d_T.eq.1)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(AlNihsle1_,Tk,1)
     SLe(2)=gibbs_T(AlNihsle2_,Tk,1)
     SLe(3)=gibbs_T(AlNihsle3_,Tk,1)
     SLe(4)=gibbs_T(AlNihsle4_,Tk,1)
     SLe(5)=AlNihsle5(2)
     SLe(6)=AlNihsle6(2)

     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     yy_d(1)=y4(c,1)
     yy_d(2)=y3(c,1)
     entSL_=g_R*(EntropyMixing(yy(1),0)+EntropyMixing(yy(2),0))
     entSL_d=g_R*(EntropyMixing(yy(1),1)*yy_d(1)+EntropyMixing(yy(2),1)*yy_d(2))


     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1)
     call get_y1_y2(yy,2,y1_y2_d2)
     
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge2_=SLe(5)*y1_y2(5)+SLe(6)*y1_y2(6)
     ge1_d1=(SLe(1)*y1_y2_d1(1)+SLe(2)*y1_y2_d1(2)+SLe(3)*y1_y2_d1(3)+SLe(4)*y1_y2_d1(4))*yy_d(1)
     ge1_d2=(SLe(1)*y1_y2_d2(1)+SLe(2)*y1_y2_d2(2)+SLe(3)*y1_y2_d2(3)+SLe(4)*y1_y2_d2(4))*yy_d(2)
     ge2_d1=(SLe(5)*y1_y2_d1(5)+SLe(6)*y1_y2_d1(6))*yy_d(1)
     ge2_d2=(SLe(5)*y1_y2_d2(5)+SLe(6)*y1_y2_d2(6))*yy_d(2)
     ge_=(ge1_+ge2_+entSL_)/(1d0+yy(1))
     ge_d=-ge_/(1d0+yy(1))*yy_d(1)+(ge1_d1+ge1_d2+ge2_d1+ge2_d2+entSL_d)/(1d0+yy(1))
     Gibbs_FE_e=ge_d
  elseif(d_c.eq.1.and.d_T.eq.2)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(AlNihsle1_,Tk,2)
     SLe(2)=gibbs_T(AlNihsle2_,Tk,2)
     SLe(3)=gibbs_T(AlNihsle3_,Tk,2)
     SLe(4)=gibbs_T(AlNihsle4_,Tk,2)


     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     yy_d(1)=y4(c,1)
     yy_d(2)=y3(c,1)



     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1)
     call get_y1_y2(yy,2,y1_y2_d2)
     
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge1_d1=(SLe(1)*y1_y2_d1(1)+SLe(2)*y1_y2_d1(2)+SLe(3)*y1_y2_d1(3)+SLe(4)*y1_y2_d1(4))*yy_d(1)
     ge1_d2=(SLe(1)*y1_y2_d2(1)+SLe(2)*y1_y2_d2(2)+SLe(3)*y1_y2_d2(3)+SLe(4)*y1_y2_d2(4))*yy_d(2)
     ge_=ge1_/(1d0+yy(1))
     ge_d=-ge_/(1d0+yy(1))*yy_d(1)+(ge1_d1+ge1_d2)/(1d0+yy(1))
     Gibbs_FE_e=ge_d
  elseif(d_c.eq.0.and.d_T.eq.1)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(AlNihsle1_,Tk,1)
     SLe(2)=gibbs_T(AlNihsle2_,Tk,1)
     SLe(3)=gibbs_T(AlNihsle3_,Tk,1)
     SLe(4)=gibbs_T(AlNihsle4_,Tk,1)
     SLe(5)=AlNihsle5(2)
     SLe(6)=AlNihsle6(2)

     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     entSL_=g_R*(EntropyMixing(yy(1),0)+EntropyMixing(yy(2),0))
     call get_y1_y2(yy,0,y1_y2)
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge2_=SLe(5)*y1_y2(5)+SLe(6)*y1_y2(6)
     ge_=(ge1_+ge2_+entSL_)/(1d0+yy(1))
     Gibbs_FE_e=ge_
  elseif(d_c.eq.0.and.d_T.eq.2)then
     Tk=T_k(T)
     SLe(1)=gibbs_T(AlNihsle1_,Tk,2)
     SLe(2)=gibbs_T(AlNihsle2_,Tk,2)
     SLe(3)=gibbs_T(AlNihsle3_,Tk,2)
     SLe(4)=gibbs_T(AlNihsle4_,Tk,2)
     yy(1)=y4(c,0)
     yy(2)=y3(c,0)
     call get_y1_y2(yy,0,y1_y2)
     ge1_=SLe(1)*y1_y2(1)+SLe(2)*y1_y2(2)+SLe(3)*y1_y2(3)+SLe(4)*y1_y2(4)
     ge_=ge1_/(1d0+yy(1))
     Gibbs_FE_e=ge_
  endif
end function Gibbs_FE_e
!begin AlSi

!y:=(x,y1,y2,dy1)->(1-x)*((1-x)*y1+x*(dy1/2+y1))+x*(x*y2+(1-x)*(dy1/2+y1));
!Y:=(c,x1,x2,y1,y2,dy1)->y((c-x1)/(x2-x1),y1,y2,dy1*(x2-x1));

double precision function Bezier_1(x,y1,y2,dy1,d)
double precision, intent(in):: x,y1,y2,dy1
integer, intent(in) ::d
if(d.eq.0)then
   Bezier_1 = (1d0-x)*((1-x)*y1+x*(dy1/2+y1))+x*(x*y2+(1-x)*(dy1/2+y1))
else
   Bezier_1 = (1d0-x)*dy1+x*(2d0*(y2-y1)-dy1)
endif
end function Bezier_1

double Precision function Bezier_2(c,x1,x2,y1,y2,dy1,d)
double precision, intent(in)::c,x1,x2,y1,y2,dy1
integer, intent(in) ::d
double precision Bezier_1
if(d.eq.0)then
   Bezier_2 = Bezier_1((c-x1)/(x2-x1),y1,y2,dy1*(x2-x1),0)
else
   Bezier_2 = Bezier_1((c-x1)/(x2-x1),y1,y2,dy1*(x2-x1),1)/(x2-x1)
endif
end function Bezier_2
double precision function Quad(c,c0,f,fc,fcc,d)
double precision, intent(in)::c,c0,f,fc,fcc
integer, intent(in)::d
if(d.eq.0)then
   Quad = 0.5*fcc*(c-c0)**2 + fc*(c-c0) + f
else
   Quad = fcc*(c-c0) + fc
endif
end function Quad





double precision function Gibbs_FE_sol(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in)::d_c,d_T
  double precision :: SLb(6),yy(2),yy_d(2),y1_y2(6),y1_y2_d1(6),y1_y2_d2(6)
  double precision :: y1,y2,EntropyMixing,T_k
  double precision ::entSL_,entSL_d,gb1_,gb1_d1,gb1_d2,gb2_,gb2_d1,gb2_d2,gb_,gb_d,Tk,Gibbs_T
  
  if(d_c.eq.0.and.d_T.eq.0)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(AlNihslb1,Tk,0)
     SLb(2)=gibbs_T(AlNihslb2,Tk,0)
     SLb(3)=gibbs_T(AlNihslb3,Tk,0)
     SLb(4)=gibbs_T(AlNihslb4,Tk,0)
     SLb(5)=AlNihslb5(1)+AlNihslb5(2)*Tk
     SLb(6)=AlNihslb6(1)+AlNihslb6(2)*Tk
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     entSL_=g_R*Tk*(EntropyMixing(yy(1),0)/3.+EntropyMixing(yy(2),0)/6.)
     call get_y1_y2(yy,0,y1_y2)
     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)
     gb2_=SLb(5)*y1_y2(5)+SLb(6)*y1_y2(6)
     gb_=(gb1_+gb2_+entSL_)/(1-yy(2)/6d0)
     Gibbs_FE_sol=gb_
  elseif(d_c.eq.0.and.d_T.eq.1)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(AlNihslb1,Tk,1)
     SLb(2)=gibbs_T(AlNihslb2,Tk,1)
     SLb(3)=gibbs_T(AlNihslb3,Tk,1)
     SLb(4)=gibbs_T(AlNihslb4,Tk,1)
     SLb(5)=AlNihslb5(2)
     SLb(6)=AlNihslb6(2)
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     entSL_=g_R*(EntropyMixing(yy(1),0)/3.+EntropyMixing(yy(2),0)/6.)
     call get_y1_y2(yy,0,y1_y2)
     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)
     gb2_=SLb(5)*y1_y2(5)+SLb(6)*y1_y2(6)
     gb_=(gb1_+gb2_+entSL_)/(1-yy(2)/6d0)
     Gibbs_FE_sol=gb_
  elseif(d_c.eq.0.and.d_T.eq.2)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(AlNihslb1,Tk,2)
     SLb(2)=gibbs_T(AlNihslb2,Tk,2)
     SLb(3)=gibbs_T(AlNihslb3,Tk,2)
     SLb(4)=gibbs_T(AlNihslb4,Tk,2)
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     call get_y1_y2(yy,0,y1_y2)
     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)
     gb_=gb1_/(1-yy(2)/6d0)
     Gibbs_FE_sol=gb_
  elseif(d_c.eq.1.and.d_T.eq.0)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(AlNihslb1,Tk,0)
     SLb(2)=gibbs_T(AlNihslb2,Tk,0)
     SLb(3)=gibbs_T(AlNihslb3,Tk,0)
     SLb(4)=gibbs_T(AlNihslb4,Tk,0)
     SLb(5)=AlNihslb5(1)+AlNihslb5(2)*Tk
     SLb(6)=AlNihslb6(1)+AlNihslb6(2)*Tk
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     yy_d(1)=y1(c,1)
     yy_d(2)=y2(c,1)

     entSL_=g_R*Tk*(EntropyMixing(yy(1),0)/3.+EntropyMixing(yy(2),0)/6.)
     entSL_d=g_R*Tk*(yy_d(1)*EntropyMixing(yy(1),1)/3.+yy_d(2)*EntropyMixing(yy(2),1)/6.)


     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1) !d(y1-y2)/d(y1)
     call get_y1_y2(yy,2,y1_y2_d2) !d(y1-y2)/d(y2)

     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)

     gb1_d1=SLb(1)*y1_y2_d1(1)+SLb(2)*y1_y2_d1(2)+SLb(3)*y1_y2_d1(3)+SLb(4)*y1_y2_d1(4)
     gb1_d1=gb1_d1*yy_d(1) !d(y1)/d(c)*d(gb1)/d(y1)

     gb1_d2=SLb(1)*y1_y2_d2(1)+SLb(2)*y1_y2_d2(2)+SLb(3)*y1_y2_d2(3)+SLb(4)*y1_y2_d2(4)
     gb1_d2=gb1_d2*yy_d(2)!d(y2)/d(c)*d(gb1)/d(y2)


     gb2_=SLb(5)*y1_y2(5)+SLb(6)*y1_y2(6)

     gb2_d1=SLb(5)*y1_y2_d1(5)+SLb(6)*y1_y2_d1(6)
     gb2_d1=gb2_d1*yy_d(1) !d(y1)/d(c)*d(gb2)/d(y1) 

     gb2_d2=SLb(5)*y1_y2_d2(5)+SLb(6)*y1_y2_d2(6)
     gb2_d2=gb2_d2*yy_d(2)!d(y2)/d(c)*d(gb2)/d(y2) 


     gb_=(gb1_+gb2_+entSL_)/(1-yy(2)/6d0)
     gb_d=(gb1_d1+gb1_d2+gb2_d1+gb2_d2+entSL_d)/(1-yy(2)/6d0)-(gb_)/(1-yy(2)/6d0)*(-yy_d(2)/6d0)
     Gibbs_FE_sol=gb_d
  elseif(d_c.eq.1.and.d_T.eq.1)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(AlNihslb1,Tk,1)
     SLb(2)=gibbs_T(AlNihslb2,Tk,1)
     SLb(3)=gibbs_T(AlNihslb3,Tk,1)
     SLb(4)=gibbs_T(AlNihslb4,Tk,1)
     SLb(5)=AlNihslb5(2)
     SLb(6)=AlNihslb6(2)
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     yy_d(1)=y1(c,1)
     yy_d(2)=y2(c,1)

     entSL_=g_R*(EntropyMixing(yy(1),0)/3.+EntropyMixing(yy(2),0)/6.)
     entSL_d=g_R*(yy_d(1)*EntropyMixing(yy(1),1)/3.+yy_d(2)*EntropyMixing(yy(2),1)/6.)


     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1) !d(y1-y2)/d(y1)
     call get_y1_y2(yy,2,y1_y2_d2) !d(y1-y2)/d(y2)

     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)

     gb1_d1=SLb(1)*y1_y2_d1(1)+SLb(2)*y1_y2_d1(2)+SLb(3)*y1_y2_d1(3)+SLb(4)*y1_y2_d1(4)
     gb1_d1=gb1_d1*yy_d(1) !d(y1)/d(c)*d(gb1)/d(y1)

     gb1_d2=SLb(1)*y1_y2_d2(1)+SLb(2)*y1_y2_d2(2)+SLb(3)*y1_y2_d2(3)+SLb(4)*y1_y2_d2(4)
     gb1_d2=gb1_d2*yy_d(2)!d(y2)/d(c)*d(gb1)/d(y2)


     gb2_=SLb(5)*y1_y2(5)+SLb(6)*y1_y2(6)

     gb2_d1=SLb(5)*y1_y2_d1(5)+SLb(6)*y1_y2_d1(6)
     gb2_d1=gb2_d1*yy_d(1) !d(y1)/d(c)*d(gb2)/d(y1) 

     gb2_d2=SLb(5)*y1_y2_d2(5)+SLb(6)*y1_y2_d2(6)
     gb2_d2=gb2_d2*yy_d(2)!d(y2)/d(c)*d(gb2)/d(y2) 


     gb_=(gb1_+gb2_+entSL_)/(1-yy(2)/6d0)
     gb_d=(gb1_d1+gb1_d2+gb2_d1+gb2_d2+entSL_d)/(1-yy(2)/6d0)-(gb_)/(1-yy(2)/6d0)*(-yy_d(2)/6d0)
     Gibbs_FE_sol=gb_d
  elseif(d_c.eq.1.and.d_T.eq.2)then
     Tk=T_k(T)
     SLb(1)=gibbs_T(AlNihslb1,Tk,2)
     SLb(2)=gibbs_T(AlNihslb2,Tk,2)
     SLb(3)=gibbs_T(AlNihslb3,Tk,2)
     SLb(4)=gibbs_T(AlNihslb4,Tk,2)
     SLb=SLb/6d0
     yy(1)=y1(c,0)
     yy(2)=y2(c,0)
     yy_d(1)=y1(c,1)
     yy_d(2)=y2(c,1)

     call get_y1_y2(yy,0,y1_y2)
     call get_y1_y2(yy,1,y1_y2_d1) !d(y1-y2)/d(y1)
     call get_y1_y2(yy,2,y1_y2_d2) !d(y1-y2)/d(y2)

     gb1_=SLb(1)*y1_y2(1)+SLb(2)*y1_y2(2)+SLb(3)*y1_y2(3)+SLb(4)*y1_y2(4)

     gb1_d1=SLb(1)*y1_y2_d1(1)+SLb(2)*y1_y2_d1(2)+SLb(3)*y1_y2_d1(3)+SLb(4)*y1_y2_d1(4)
     gb1_d1=gb1_d1*yy_d(1) !d(y1)/d(c)*d(gb1)/d(y1)

     gb1_d2=SLb(1)*y1_y2_d2(1)+SLb(2)*y1_y2_d2(2)+SLb(3)*y1_y2_d2(3)+SLb(4)*y1_y2_d2(4)
     gb1_d2=gb1_d2*yy_d(2)!d(y2)/d(c)*d(gb1)/d(y2)

     gb_=gb1_/(1-yy(2)/6d0)
     gb_d=gb1_d1/(1-yy(2)/6d0)-(gb_)/(1-yy(2)/6d0)*(-yy_d(2)/6d0)
     Gibbs_FE_sol=gb_d

  else
     write(*,*)"unexpected d_c d_T combination",d_c,d_T
  endif
  
end function Gibbs_FE_sol




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function AlNiFreeEnergy(c,T,phi,c_dot,phi_dot,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi(N_phases),c_dot,phi_dot(N_phases)
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(2),hb(2),hr(2)
double precision :: Tk,Cp,R
double precision :: FE,dh=1d-6
double precision Gibbs_FE_sol,Gibbs_FE_l,Gibbs_FE_e,T_K,g_phi,EntropyMixing,Gibbs_SiGeliq,Gibbs_SiGesol !functions
double precision :: Gibbs_FE_l_AlSi,Gibbs_FE_e_AlSi,cc,Gibbs_Fe_l_PbSn,Gibbs_Fe_e_PbSn,S,dS
integer :: i,j



AlNiFreeEnergy=0d0


if(lp.eq.0)then !free energy proper J/mol
   
   !      AlNiFreeEnergy =              Gibbs_FE_liq(c,T,0,0)*g_phi(phi(1),0)
   !      AlNiFreeEnergy = AlNiFreeEnergy + Gibbs_FE_sol(c,T,0,0)*g_phi(1d0-phi(1),0)
   
   AlNiFreeEnergy =                    Gibbs_FE_l(c,T,0,0)*g_phi(phi(1),0)
   AlNiFreeEnergy = AlNiFreeEnergy + Gibbs_FE_sol(c,T,0,0)*g_phi(phi(2),0)
   AlNiFreeEnergy = AlNiFreeEnergy +   Gibbs_FE_e(c,T,0,0)*g_phi(phi(3),0)
   
   S = 1d0/(g_phi(phi(1),0)+g_phi(phi(2),0)+g_phi(phi(3),0))
   
   AlNiFreeEnergy = AlNiFreeEnergy*S 
   return
elseif(lp.le.n_phases)then

   S = 1d0/(g_phi(phi(1),0)+g_phi(phi(2),0)+g_phi(phi(3),0))
   
   dS = -1d0/S**2*g_phi(phi(lp),1)
   AlNiFreeEnergy =                    Gibbs_FE_l(c,T,0,0)*g_phi(phi(1),0)
   AlNiFreeEnergy = AlNiFreeEnergy + Gibbs_FE_sol(c,T,0,0)*g_phi(phi(2),0)
   AlNiFreeEnergy = AlNiFreeEnergy +   Gibbs_FE_e(c,T,0,0)*g_phi(phi(3),0)
   
   AlNiFreeEnergy = AlNiFreeEnergy*dS
   if(lp.eq.1)then
      AlNiFreeEnergy = AlNiFreeEnergy +   Gibbs_FE_l(c,T,0,0)*g_phi(phi(1),1)*S
   elseif(lp.eq.2)then
      AlNiFreeEnergy = AlNiFreeEnergy + Gibbs_FE_sol(c,T,0,0)*g_phi(phi(2),1)*S
   else
      AlNiFreeEnergy = AlNiFreeEnergy +   Gibbs_FE_e(c,T,0,0)*g_phi(phi(3),1)*S
   endif
   return
elseif(lp.eq.n_phases+2)then !Solute
   
   
   AlNiFreeEnergy =                    Gibbs_FE_l(c,T,1,0)*g_phi(phi(1),0)
   AlNiFreeEnergy = AlNiFreeEnergy + Gibbs_FE_sol(c,T,1,0)*g_phi(phi(2),0)
   AlNiFreeEnergy = AlNiFreeEnergy +   Gibbs_FE_e(c,T,1,0)*g_phi(phi(3),0)
   
   S = 1d0/(g_phi(phi(1),0)+g_phi(phi(2),0)+g_phi(phi(3),0))

   AlNiFreeEnergy = AlNiFreeEnergy*S
   return
elseif(lp.eq.n_phases+1)then !Temperature Keep in J/mol
   Tk=T_k(T)

      !-(1-T*d/dT)dF/dphi  
   g_latent_heat(1) =  Tk*  Gibbs_FE_l(c,T,0,1)-  Gibbs_FE_l(c,T,0,0)
   g_latent_heat(2) =  Tk*Gibbs_FE_sol(c,T,0,1)-Gibbs_FE_sol(c,T,0,0)
   g_latent_heat(3) =  Tk*  Gibbs_FE_e(c,T,0,1)-  Gibbs_FE_e(c,T,0,0)
   


   !-(1-T*d/dT)dF/dphi
   if(heat_phi_term)then
      

      
      AlNiFreeEnergy = g_latent_heat(1)*phi_dot(1)*g_phi(phi(1),1) +g_latent_heat(2)*phi_dot(2)*g_phi(phi(2),1) +g_latent_heat(3)*phi_dot(3)*g_phi(phi(3),1)
      
      
      
      Q_heat(1) = AlNiFreeEnergy
      !      endif !end phi_dot_term
!!$      if(heat_c_term)then 
!!$         write(*,*)"heat_c_term not coded yet"
!!$         stop
!!$         !   -(1- T(d/dT))dF/dc 
!!$         FE=(Tk*Gibbs_FE_liq(c,T,1,1)*g_phi(phi(1),0)- Gibbs_FE_liq(c,T,1,0))*g_phi(phi(1),0)
!!$         FE=FE+(Tk*Gibbs_FE_sol(c,T,1,1)- Gibbs_FE_sol(c,T,1,0))*(1.-g_phi(phi(1),0))
!!$         FE=FE*c_dot
!!$         AlNiFreeEnergy = AlNiFreeEnergy + FE
!!$         Q_heat(2)= FE
!!$         !      AlNiFreeEnergy = AlNiFreeEnergy + Tk*Gibbs_FE_liq(c,T,1,1)*g_phi(phi(1),0)*c_dot
!!$         !      AlNiFreeEnergy = AlNiFreeEnergy + Tk*Gibbs_FE_sol(c,T,1,1)*(1.-g_phi(phi(1),0))*c_dot
!!$         !      AlNiFreeEnergy = AlNiFreeEnergy - Gibbs_FE_liq(c,T,1,0)*g_phi(phi(1),0)*c_dot
!!$         !      AlNiFreeEnergy = AlNiFreeEnergy - Gibbs_FE_sol(c,T,1,0)*(1.-g_phi(phi(1),0))*c_dot
!!$      endif !c_dot_term
      
      
      !   heat capacity = -T (d/dT) dF/dT
      !      g_C =  -Tk*(g_phi(phi(1),0)*Gibbs_FE_liq(c,T,0,2)+g_phi(1d0-phi(1),0)*Gibbs_FE_sol(c,T,0,2))
      g_C = -Tk*(Gibbs_FE_l(c,T,0,2)*g_phi(phi(1),0)+Gibbs_FE_sol(c,T,0,2)*g_phi(phi(2),0)+Gibbs_FE_e(c,T,0,2)*g_phi(phi(3),0))
      
      return
   endif
elseif(lp.eq.-1)then !for anti-trapping
   write(*,*)"lp=-1 in AlNiFreeEnergy"
   stop
   AlNiFreeEnergy = AlNiFreeEnergy + Gibbs_FE_l(c,T,1,0)*g_phi(phi(1),1)
   AlNiFreeEnergy = AlNiFreeEnergy - Gibbs_FE_e(c,T,1,0)*g_phi(phi(1),1)
   
   return
elseif(lp.eq.-2)then !for anti-trapping f_cc
   write(*,*)"lp=-2 in AlNiFreeEnergy"
   stop      
   !resorting to approximation for gibbs_fe_e as it looks too complicated
   AlNiFreeEnergy =              0.5*(Gibbs_FE_l(c+dh,T,1,0)-Gibbs_FE_l(c-dh,T,1,0))/dh*g_phi(phi(1),0)
   AlNiFreeEnergy = AlNiFreeEnergy + 0.5*(Gibbs_FE_sol(c+dh,T,1,0)-Gibbs_FE_sol(c-dh,T,1,0))/dh*(1d0-g_phi(phi(1),0))
   
   return
else
   write(*,*)"bad current_var"
   stop
endif

end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function AlFeFreeEnergy(c,T,phi,c_dot,phi_dot,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi(N_phases),c_dot,phi_dot(N_phases)
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(2),hb(2),hr(2)
double precision :: Tk,Cp,R,S
double precision :: FE,dh=1d-6
double precision Gibbs_FE_l_AlFe,Gibbs_FE_e_AlFe,Gibbs_FE_e2_AlFe,g_phi,dAlFeFreeEnergy
integer :: i,j



AlFeFreeEnergy=0d0


if(lp.eq.0)then !free energy proper J/mol
   AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,T,0,0)*g_phi(phi(1),0)
   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e_AlFe(c,T,0,0)*g_phi(phi(2),0)
   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e2_AlFe(c,T,0,0)*g_phi(phi(3),0)
   
   S = g_phi(phi(1),0)+g_phi(phi(2),0)+g_phi(phi(3),0)
   
   AlFeFreeEnergy = AlFeFreeEnergy/S 
   return
elseif(lp.le.n_phases)then
   if(g_Tflag)then
      write(*,*)g_Tflag
      stop
      AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,T,0,1)*g_phi(phi(1),0)
      AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e_AlFe(c,T,0,1)*g_phi(phi(2),0)
      AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e2_AlFe(c,T,0,1)*g_phi(phi(3),0)
      
      S = g_phi(phi(1),0)+g_phi(phi(2),0)+g_phi(phi(3),0)
      dAlFeFreeEnergy =   Gibbs_FE_l_AlFe(c,T,0,0)*g_phi(phi(lp),1)
      AlFeFreeEnergy = dAlFeFreeEnergy/S - AlFeFreeEnergy/S**2
      return
   else
      AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,T,0,0)*g_phi(phi(1),0)
      AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e_AlFe(c,T,0,0)*g_phi(phi(2),0)
      AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e2_AlFe(c,T,0,0)*g_phi(phi(3),0)
      
      S = g_phi(phi(1),0)+g_phi(phi(2),0)+g_phi(phi(3),0)
      if(lp.eq.1)then
         dAlFeFreeEnergy =   Gibbs_FE_l_AlFe(c,T,0,0)
      elseif(lp.eq.2)then
         dAlFeFreeEnergy =   Gibbs_FE_e_AlFe(c,T,0,0)
      else
         dAlFeFreeEnergy =   Gibbs_FE_e2_AlFe(c,T,0,0)
      endif
      AlFeFreeEnergy = (dAlFeFreeEnergy/S - AlFeFreeEnergy/S**2)*g_phi(phi(lp),1)
      return
   endif
elseif(lp.eq.n_phases+2)then !Solute
   AlFeFreeEnergy =                     Gibbs_FE_l_AlFe(c,T,1,0)*g_phi(phi(1),0)
   AlFeFreeEnergy = AlFeFreeEnergy +    Gibbs_FE_e_AlFe(c,T,1,0)*g_phi(phi(2),0)
   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e2_AlFe(c,T,1,0)*g_phi(phi(3),0)
   
   S = g_phi(phi(1),0)+g_phi(phi(2),0)+g_phi(phi(3),0)
   
   AlFeFreeEnergy = AlFeFreeEnergy/S 
elseif(lp.eq.n_phases+1)then !Temperature Keep in J/mol
   AlFeFreeEnergy =                    Gibbs_FE_l_AlFe(c,T,0,1)*g_phi(phi(1),0)
   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e_AlFe(c,T,0,1)*g_phi(phi(2),0)
   AlFeFreeEnergy = AlFeFreeEnergy +   Gibbs_FE_e2_AlFe(c,T,0,1)*g_phi(phi(3),0)
   
   S = g_phi(phi(1),0)+g_phi(phi(2),0)+g_phi(phi(3),0)
   
   AlFeFreeEnergy = AlFeFreeEnergy/S 
elseif(lp.eq.-1)then !for anti-trapping
write(*,*)"not coded"
stop
elseif(lp.eq.-2)then !for anti-trapping f_cc
write(*,*)"not coded"
stop
else
   write(*,*)"bad current_var"
   stop
endif

end function



subroutine Project_Phase(vbles,dx,phi_rhs,Phi_,c,T)
use solution_parameters
implicit none
!integer, intent(in)::df
double precision, intent(in)::dx,c,T,phi_(n_phases),vbles(N_unk+1,3,3)
double precision, intent(out) :: phi_rhs(n_phases)
double precision PhaseRHS,Wk(5)
integer ip,jp,izero
double precision, dimension(5,5)::Projector = 0.
double precision :: q=1d-9,phi(5),S
double precision :: Anisotropy(3,3),normal(3,3,2),g(3,3),eps = 0.1

!!$do ip= 1,n_phases
!!$   do jp = 1 ,n_phases
!!$      normal(ip,jp,1) = vbles(jp,2,2)*0.5*(vbles(ip,3,2) - vbles(ip,1,2))/dx - vbles(ip,2,2)*0.5*(vbles(jp,3,2) - vbles(jp,1,2))/dx
!!$      normal(ip,jp,2) = vbles(jp,2,2)*0.5*(vbles(ip,2,3) - vbles(ip,2,1))/dx - vbles(ip,2,2)*0.5*(vbles(jp,2,3) - vbles(jp,2,1))/dx
!!$      g(ip,jp) = sqrt(normal(ip,jp,1)**2 + normal(ip,jp,2)**2 + 1d-16)
!!$   enddo
!!$enddo
!!$do ip=1,n_phases
!!$   do jp=1,n_phases
!!$      normal(ip,jp,1) = normal(ip,jp,1)/g(ip,jp)
!!$      normal(ip,jp,2) = normal(ip,jp,2)/g(ip,jp)
!!$   enddo
!!$enddo

!!$
!!$Anisotropy = 0d0
!!$Anisotropy(2,3) = (1d0 + epsilon*abs(normal(2,3,2)**3-3*normal(2,3,2)*normal(2,3,1)**2))/(1d0+epsilon)
!!$Anisotropy(3,2) = Anisotropy(2,3)
!!$Anisotropy(1,2) = (1d0 + epsilon*abs(normal(1,2,2)**3-3*normal(1,2,2)*normal(1,2,1)**2))/(1d0+epsilon)
!!$Anisotropy(2,1) = Anisotropy(1,2)
!!$Anisotropy(1,3) = 1d0
!!$Anisotropy(3,1) = Anisotropy(1,3)
!!$

!cos^6 + 6 cos^5sin + 15 cos^4sin^2 + 20 cos^3sin^3 +15 cos^2sin^4 + 6 cos sin^5 + sin^6 
!cos^6 +            - 15 cos^4sin^2 +               +15 cos^2sin^4 +             - sin^6 
!!$
!!$Anisotropy = 0d0
!!$Anisotropy(2,3) = (1d0 + 0.5*epsilon*(1+normal(2,3,2)**6 - 15*normal(2,3,2)**4*normal(2,3,1)**2 + 15*normal(2,3,2)**2*normal(2,3,1)**4 - normal(2,3,1)**6))/(1d0+epsilon)
!!$Anisotropy(3,2) = Anisotropy(2,3)
!!$Anisotropy(1,2) = (1d0 + 0.5*epsilon*(1+normal(1,2,2)**6 - 15*normal(1,2,2)**4*normal(1,2,1)**2 + 15*normal(1,2,2)**2*normal(1,2,1)**4 - normal(1,2,1)**6))/(1d0+epsilon)
!!$Anisotropy(2,1) = Anisotropy(1,2)
!!$Anisotropy(1,3) = 1d0
!!$Anisotropy(3,1) = Anisotropy(1,3)
!!$
!!$if(n_phase.gt.3)then
!!$   write(*,*)"in projectphase, N_phases > 3"
!!$endif

  S=0d0
  do ip = 1,n_phases
!     phi(ip) = max(1d-2,min(1d0-1d-2,phi_(ip)))
     phi(ip) = max(0d0,min(1d0,phi_(ip)))
     S=S+phi(ip)
  enddo
  phi=phi/S

  phi_rhs=0.

  !assemble Projector
        
  Projector=0d0
  do ip=2,n_phases
     do jp=1,ip-1
        Projector(ip,jp)=-(phi(ip)+q)/(1d0-phi(ip)+2*q)*(phi(jp)+q)/(1d0-phi(jp)+2*q)
!                      if(g_model.eq.1)Projector(ip,jp)=Projector(ip,jp)*g_Lij(ip,jp)
!        Projector(ip,jp) = Projector(ip,jp)*Anisotropy(ip,jp)
     enddo
  enddo
  do ip=2,n_phases
     do jp=1,ip-1
        Projector(jp,ip)=Projector(ip,jp)
     enddo
  enddo
  do ip=1,N_phases
     do jp=1,N_phases
        if(jp.ne.ip)Projector(ip,ip)=Projector(ip,ip)-Projector(ip,jp)
     enddo
  enddo
!debug
!!$Projector =0d0
!!$Projector(1,1)=1
!!$Projector(1,3)=-1
!!$Projector(3,1)=-1
!!$Projector(3,3)=1
  



  
  do ip=1,N_phases
     Wk(ip) =  PhaseRHS(vbles,dx,ip,Phi,c,T)
  enddo        
  
  
  !Project phi_vector
  do ip=1,N_phases
     phi_rhs(ip)=0d0
     do jp=1,N_phases
        phi_rhs(ip)=phi_rhs(ip)+Projector(ip,jp)*Wk(jp)
     enddo
  enddo

  
end subroutine Project_Phase

double precision function GradientEnergy(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision :: phix(3),phiy(3),LapPhi(5),eps,d_eps,gradientEnergy0,phi(3),s,phi_(3)
integer :: i,j,k


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
     enddo
     GradientEnergy = 0d0
     if(lp.eq.0)then
        do i=1,n_phases
           do j=1,n_phases
              GradientEnergy=GradientEnergy + g_eps22(i,j)*phi(j)*(phix(i)**2+phiy(i)**2) 
           enddo
        enddo
        return
     endif
     GradientEnergy = 0d0
     do i=1,n_phases
        
        GradientEnergy = GradientEnergy - g_eps22(lp,i)*(2d0*(phix(lp)*phix(i)+phiy(lp)*phiy(i))+2d0*phi(i)*LapPhi(lp)-(phix(i)**2+phiy(i)**2))

     enddo

  else
     write(*,*)"in GradientEnergy: lp > n_phases"
  endif

end function GradientEnergy

double precision function GradientEnergyDelta(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision :: phix(3),phiy(3),phixx(3),phiyy(3),phixy(3),LapPhi(5),eps,d_eps,gradientEnergy0,phi(3),s,phi_(3),U(3),dU(3)
integer :: i,j,k,ctr(6)=(/1,2,3,1,2,3/)
double precision :: d(3)=(/1,1,1/) 
double precision :: g0,X1,X2,theta,A,A2,pi,kappa,kappa0,dh,dirac,epsilon0=1d0

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
     GradientEnergyDelta = 0d0

!!!!!!model 4 !!!!!!!!1


     if(g_model.eq.4)then
        if(lp.eq.0)then
           do i=1,n_phases
              GradientEnergyDelta=GradientEnergyDelta + 0.5*d(i)*U(i)*(phix(i)**2+phiy(i)**2) 
           enddo
           return
        endif
        
        
        g0  = sqrt(phix(lp)**2 +phiy(lp)**2)
        if(g0.gt.1d-9)then
           X1=max((phix(lp))/g0,1d-9)
           X2=max((phiy(lp))/g0,0d0)
           theta=atan2(abs(X2),abs(X1))
           if(abs(theta).le.pi/6.)then
              !              A =2*X1**2
              !              A2 =(-4*X1**2+4*X2**2)
              A = cos(theta)**2
              A2=2.*cos(theta)*sin(theta)
           else
              !              A = (1.5-X1**2+sqrt(3.)*X1*X2)
              !              A2 = (2*X1**2-2*X2**2-4*sqrt(3.)*X1*X2)
              A = cos(theta-pi/3.)**2
              A2 = 2.*cos(theta-pi/3.)*sin(theta-pi/3.)
           endif
           
           kappa = phixx(lp)*X2**2 - 2.* phixy(lp)*X1*X2 + phiyy(lp)*X1**2
           !           kappa0=min(g0*1d-1,kappa)
           
           !           h0=min(0.1,abs(4./g_nuc_radius))
           
           
           
           dirac = 0d0
           if(abs(theta-0.5*pi).lt.h0)then
              
              dh = theta-0.5*pi
              !              A =  ((11*h0**6-15*h0**4*dh**2+5*h0**2*dh**4-dh**6)*(-2*cos(Pi/6.-h0)*sin(Pi/6.-h0))+16*cos(pi/6.-h0)**2*h0**5)/(16*h0**5)
              dirac = 15./16.*(1-dh/h0)**2*(1+dh/h0)**2/h0*sqrt(3.)
              !              dirac = sqrt(3.)*g0/dx
              !              dirac =(-30*(dh)**4+60*(dh)**2*h0**2-30*h0**4)*(-2*cos(Pi/6.-h0)*sin(Pi/6.-h0))/(8*h0**5)
              !              A2 = A2 +dirac
           elseif(abs(theta-pi/6.).lt.h0)then
              dh = theta-pi/6.
              !              A =  ((11*h0**6-15*h0**4*dh**2+5*h0**2*dh**4-dh**6)*(-2*cos(Pi/6.-h0)*sin(Pi/6.-h0))+16*cos(pi/6.-h0)**2*h0**5)/(16*h0**5)
              dirac = 15./16.*(1-dh/h0)**2*(1+dh/h0)**2/h0*sqrt(3.)
              !              dirac = sqrt(3.)*g0/dx
              !              dirac =(-30*(dh)**4+60*(dh)**2*h0**2-30*h0**4)*(-2*cos(Pi/6.-h0)*sin(Pi/6.-h0))/(8*h0**5)
              !              A2 = A2 + dirac
           endif
           
           
        else
           A = 1d0
           A2 = 0
        endif
        
        
        
        dU(lp)=1d0
        dU(ctr(lp+1))=3d0*phi(ctr(lp+2))
        dU(ctr(lp+2))=3d0*phi(ctr(lp+1))
        
        GradientEnergyDelta = d(lp)*U(lp)*(LapPhi(lp)*A+0.5*A2*kappa+0.5*dirac*kappa)&
             - d(ctr(lp+0))*dU(ctr(lp+0))*0.5*(phix(ctr(lp+0))**2+phiy(ctr(lp+0))**2)*A&
             - d(ctr(lp+1))*dU(ctr(lp+1))*0.5*(phix(ctr(lp+1))**2+phiy(ctr(lp+1))**2)*A&
             - d(ctr(lp+2))*dU(ctr(lp+2))*0.5*(phix(ctr(lp+2))**2+phiy(ctr(lp+2))**2)*A
        
        
!!$
!!$        if(lp.eq.3)then
!!$           GradientEnergyDelta = d(lp)*U(lp)*(LapPhi(lp)*(A-0.5)+0.5*A2*kappa)&
!!$                - d(ctr(lp+0))*dU(ctr(lp+0))*0.5*(phix(ctr(lp+0))**2+phiy(ctr(lp+0))**2)*(A-0.5)&
!!$                - d(ctr(lp+1))*dU(ctr(lp+1))*0.5*(phix(ctr(lp+1))**2+phiy(ctr(lp+1))**2)*0.5&
!!$                - d(ctr(lp+2))*dU(ctr(lp+2))*0.5*(phix(ctr(lp+2))**2+phiy(ctr(lp+2))**2)*0.5
!!$        elseif(lp.eq.1)then
!!$           GradientEnergyDelta = d(lp)*U(lp)*(LapPhi(lp)*0.5)&
!!$                - d(ctr(lp+0))*dU(ctr(lp+0))*0.5*(phix(ctr(lp+0))**2+phiy(ctr(lp+0))**2)*0.5&
!!$                - d(ctr(lp+1))*dU(ctr(lp+1))*0.5*(phix(ctr(lp+1))**2+phiy(ctr(lp+1))**2)*0.5&
!!$                - d(ctr(lp+2))*dU(ctr(lp+2))*0.5*(phix(ctr(lp+2))**2+phiy(ctr(lp+2))**2)*(A-0.5)
!!$        else
!!$           GradientEnergyDelta = d(lp)*U(lp)*(LapPhi(lp)*0.5)&
!!$                - d(ctr(lp+0))*dU(ctr(lp+0))*0.5*(phix(ctr(lp+0))**2+phiy(ctr(lp+0))**2)*0.5&
!!$                - d(ctr(lp+1))*dU(ctr(lp+1))*0.5*(phix(ctr(lp+1))**2+phiy(ctr(lp+1))**2)*(A-0.5)&
!!$                - d(ctr(lp+2))*dU(ctr(lp+2))*0.5*(phix(ctr(lp+2))**2+phiy(ctr(lp+2))**2)*0.5
!!$        endif
!!!!!!! end of model 4 !!!!!!!!!!!!!
     else
        if(lp.eq.0)then
           do i=1,n_phases
              GradientEnergyDelta=GradientEnergyDelta + 0.5*d(i)*U(i)*(phix(i)**2+phiy(i)**2) 
           enddo
           return
        endif
        GradientEnergyDelta = d(ctr(lp))*U(ctr(lp))*LapPhi(ctr(lp))&
             + 3*d(ctr(lp))*phi(ctr(lp+2))*(phix(ctr(lp))*phix(ctr(lp+1)) + phiy(ctr(lp))*phiy(ctr(lp+1))) &
             + 3*d(ctr(lp))*phi(ctr(lp+1))*(phix(ctr(lp))*phix(ctr(lp+2)) + phiy(ctr(lp))*phiy(ctr(lp+2))) &
             +0.5*d(ctr(lp))*(phix(ctr(lp))*phix(ctr(lp)) + phiy(ctr(lp))*phiy(ctr(lp))) &
             -1.5*d(ctr(lp+1))*phi(ctr(lp+2))*(phix(ctr(lp+1))*phix(ctr(lp+1)) + phiy(ctr(lp+1))*phiy(ctr(lp+1))) &
             -1.5*d(ctr(lp+2))*phi(ctr(lp+1))*(phix(ctr(lp+2))*phix(ctr(lp+2)) + phiy(ctr(lp+2))*phiy(ctr(lp+2))) 
     endif
!#3*delta[1]^2*dp[1]*dp[2]*p[3]+3*delta[1]^2*dp[1]*dp[3]*p[2]+(1/2)*delta[1]^2*dp[1]^2+3*ddp[1]*delta[1]^2*p[2]*p[3]+ddp[1]*delta[1]^2*p[1]-(3/2*(delta[2]^2))*p[3]*dp[2#]^2-(3/2*(delta[3]^2))*p[2]*dp[3]^2
        


     GradientEnergyDelta = -GradientEnergyDelta

  else
     write(*,*)"in GradientEnergyDelta: lp > n_phases"
  endif

end function GradientEnergyDelta




!!!!!!!!!!!!
double precision function GradientEnergyPCB(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision :: phix(3),phiy(3),LapPhi(5),eps,d_eps,gradientEnergy0,phi(3),s,phi_(3)
integer :: i,j,k


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
     enddo
     GradientEnergyPCB = 0d0
     if(lp.eq.0)then
        do i=1,n_phases
           do j=1,n_phases
              GradientEnergyPCB=GradientEnergyPCB - 0.5*g_eps22(i,j)*(phi(i)+phi(j))*(4d0-3d0*(phi(i)+phi(j))*(phix(i)*phix(j) + phiy(i)*phiy(j)))
           enddo
        enddo
        return
     endif
     GradientEnergyPCB = 0d0
     do i=1,n_phases
        
        GradientEnergyPCB = GradientEnergyPCB 

     enddo

  else
     write(*,*)"in GradientEnergy: lp > n_phases"
  endif

end function GradientEnergyPCB









!!!!!!!!!!!!






double precision function GradientEnergy_toth(vbles,lp,dx)
use solution_parameters
implicit none
double precision, intent (in)::dx,vbles(N_unk+1,3,3)
integer, intent (in):: lp
double precision :: phix,phiy,LapPhi(5),eps,d_eps,GradientEnergy_toth0,phi(3),s,phi_(3)
integer :: i


    if(lp.le.n_phases)then
       do i=1,n_phases
          call Calc_GradPhi(LapPhi(i),vbles(i,:,:),dx)
       enddo
       GradientEnergy_toth0 = 0d0
       do i=1,n_phases
          phix = 0.5*(vbles(i,3,2)-vbles(i,1,2))/dx
          phiy = 0.5*(vbles(i,2,3)-vbles(i,2,1))/dx
          GradientEnergy_toth0 = GradientEnergy_toth0 + 0.5*(phix**2 +phiy**2)
          phi(i) = vbles(i,2,2) 
       enddo

      


      phi(1)=vbles(1,2,2)
      phi(2)=vbles(2,2,2)
      phi(3)=vbles(3,2,2)
      
      
      
      phi_=phi+1d-3
      eps=0d0
      s=(phi_(1)*phi_(2))**2+(phi_(2)*phi_(3))**2+(phi_(3)*phi_(1))**2
      eps=eps+g_eps22(1,2)*(phi_(1)*phi_(2))**2
      eps=eps+g_eps22(2,3)*(phi_(2)*phi_(3))**2
      eps=eps+g_eps22(3,1)*(phi_(3)*phi_(1))**2
      eps=eps/s

      if(lp.eq.0)then
         GradientEnergy_toth=eps*GradientEnergy_toth0
         return
      endif

      if(lp.eq.1)d_eps=(2*phi_(1)*phi_(2)**2*phi_(3)**2*(g_eps22(1, 2)*phi_(2)**2-g_eps22(2, 3)*phi_(2)**2-g_eps22(2, 3)*phi_(3)**2+g_eps22(3, 1)*phi_(3)**2))/s**2
      if(lp.eq.2)d_eps=(2*phi_(2)*phi_(3)**2*phi_(1)**2*(g_eps22(2, 3)*phi_(3)**2-g_eps22(3, 1)*phi_(3)**2-g_eps22(3, 1)*phi_(1)**2+g_eps22(1, 2)*phi_(1)**2))/s**2
      if(lp.eq.3)d_eps=(2*phi_(3)*phi_(1)**2*phi_(2)**2*(g_eps22(3, 1)*phi_(1)**2-g_eps22(1, 2)*phi_(1)**2-g_eps22(1, 2)*phi_(2)**2+g_eps22(2, 3)*phi_(2)**2))/s**2
      
     GradientEnergy_toth = -eps*LapPhi(lp) + d_eps*GradientEnergy_toth0


!!$      GradientEnergy = 0d0
!!$      do i=1,n_phases
!!$         GradientEnergy = GradientEnergy - g_eps22(i,lp)*(LapPhi(lp)-LapPhi(i))
!!$      enddo

   else
      write(*,*)"in GradientEnergy: lp > n_phases"
   endif







end function GradientEnergy_toth








!T_K converts non-dimensional parameter T to dimensioned T_K(T)
double precision function T_K(T)
use solution_parameters
implicit none
double precision, intent (in)::T


 
 T_K = g_T0*(1.+T)

 end function T_K


double precision function SoluteRHS(vbles,dx,c_c)!
! this uses volume element method. This guarantees conservation of solute. The main source of solute is the variation of free energy with phase: d (d F/dc) = ...+ (d^2 F/dc d phi) d phi + ..., which is revealed by a plot of c_dot.
use solution_parameters
implicit none
double precision, intent (in):: dx,c_c,vbles(N_unk+1,3,3)
integer ii,jj,ip,ih,jv
double precision c,T,Phi(n_phases),FreeEnergy,D
double precision, dimension(3,3)::stencil11=reshape((/0.,6.,0.,6.,-24.,6.,0.,6.,0./),(/3,3/))
double precision, dimension(3)::stencil1=(/-.5,0.,.5/)
double precision :: potential,D_c,f_c,MolPerVol
double precision :: c_dot=0d0,phi_dot=0d0,beta,phi_,s1,s2!



SoluteRHS=0.


!x-coordinate transform values
!! Right
 ih=3
 jv=2
 include 'c_stencil.f90'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Left
 ih=1
 jv=2
 include 'c_stencil.f90'


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!up
 ih=2
 jv=3
 include 'c_stencil.f90'

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !down
 ih=2
 jv=1
 include 'c_stencil.f90'


if(g_model.eq.6.or.g_model.eq.0.or.g_model.eq.5)then
 g_c_force = SoluteRHS/(g_Dch*dx*dx)
 SoluteRHS = SoluteRHS/(g_Dch*dx*dx) 
else
 g_c_force = SoluteRHS/(g_Dch*dx*dx*g_R*g_T0)
 SoluteRHS = SoluteRHS/(g_Dch*dx*dx*g_R*g_T0) 
endif

return

end function SoluteRHS

double precision function potential(phi,c,lp)
use solution_parameters
implicit none

double precision, intent (in)::phi(n_phases),c
integer, intent(in):: lp
integer :: i,j
double precision :: S,potential0,eps,d_eps,phi_(3), q=0d0
potential = 0d0
if(folch)then
   if(lp.le.n_phases)then
      potential0 = 0d0
      do i=1,n_phases
         potential0 = potential0 + phi(i)**2*(1d0-phi(i))**2*0.5
      enddo
      
      
      
      if(lp.eq.0)then
         potential = potential0
         return
      endif

      potential = phi(lp)*(1-phi(lp))*(1-2*phi(lp))
      
   else
      write(*,*)"in potential: lp > n_phases?"
      stop
   endif
elseif(PCB_potential)then

   if(lp.eq.0)then
         do i=2,n_phases
            do j=1,i-1
               potential = potential   -(31/2)*phi(i)**4*phi(j)-(35/2)*phi(i)**3*phi(j)**2&
                    -(35/2)*phi(i)**2*phi(j)**3-(31/2)*phi(i)*phi(j)**4+(31/2)*phi(i)**3*phi(j)&
                    +18*phi(i)**2*phi(j)**2+(31/2)*phi(i)*phi(j)**3
            enddo
         enddo
         potential = potential/16
      elseif(lp.le.n_phases)then
         do j=1,n_phases
            if(j.ne.lp)potential = potential -62*phi(lp)**3*phi(j)-(105/2)*phi(lp)**2*phi(j)**2-35*phi(lp)*phi(j)**3&
                 -(31/2)*phi(j)**4+(93/2)*phi(lp)**2*phi(j)+36*phi(lp)*phi(j)**2+(31/2)*phi(j)**3
         enddo
         potential = potential/16         
      else
         write(*,*)"lp>N"
         stop
      endif

elseif(BJM)then
   if(c_delta)then

      if(lp.eq.0)then
         do i=2,n_phases
            do j=1,i-1
               potential = potential + phi(j)**2*phi(i)**2*(2d0-phi(i)-phi(j))**4
            enddo
         enddo
      elseif(lp.le.n_phases)then
         do j=1,n_phases
            if(j.ne.lp)potential = potential + phi(lp)*phi(j)**2*(2-phi(lp)-phi(j))**3*(8-12*phi(lp)-4*phi(j))
         enddo
      else
         write(*,*)"lp>N"
         stop
      endif
   else
      if(lp.eq.0)then
         do i=2,n_phases
            do j=1,i-1
               potential = potential + g_eps22(i,j)*phi(j)**2*phi(i)**2*(2d0-phi(i)-phi(j))**4
            enddo
         enddo
      elseif(lp.le.n_phases)then
         do j=1,n_phases
            if(j.ne.lp)potential = potential + g_eps22(j,lp)*phi(lp)*phi(j)**2*(2-phi(lp)-phi(j))**3*(8-12*phi(lp)-4*phi(j))
         enddo
      else
         write(*,*)"lp>N"
         stop
      endif
   endif
endif
!!$
!!$
!!$!4th order
!!$elseif(.true.)then
!!$   if(lp.eq.0)then
!!$      do i=1,n_phases
!!$         do j=1,n_phases
!!$            potential = potential + g_eps22(i,j)*phi(j)*phi(i)**2*(1d0-phi(i))
!!$         enddo
!!$      enddo
!!$   else
!!$      
!!$      do i=1,n_phases
!!$         potential = potential + g_eps22(i,lp)*(phi(i)**2*(1-phi(i)) + phi(i)*(2*phi(lp)-3*phi(lp)**2))
!!$         !            if(phi(lp).gt.0)potential = potential - phi(1)*phi(2)*phi(3)/phi(lp)
!!$      enddo
!!$   endif
!!$   ! 5th order
!!$else
!!$   if(lp.eq.0)then
!!$      do i=1,n_phases
!!$         do j=1,n_phases
!!$            potential = potential + g_eps22(i,j)*phi(j)*phi(i)**2*(1d0-phi(i))**2
!!$         enddo
!!$      enddo
!!$   else
!!$      
!!$      do i=1,n_phases
!!$         potential = potential + g_eps22(i,lp)*(phi(i)**2*(1-phi(i))**2 + 2*phi(i)*phi(lp)*(1-phi(lp))*(1-2*phi(lp)))
!!$         !            if(phi(lp).gt.0)potential = potential - phi(1)*phi(2)*phi(3)/phi(lp)
!!$      enddo
!!$   endif
!!$else
!!$   if(lp.le.n_phases)then
!!$      potential = 0d0
!!$      if(lp.eq.0)then
!!$         do i=1,n_phases
!!$            do j=1,n_phases
!!$               potential = potential + g_eps22(i,j)*phi(i)*0.5*phi(i)**2*(1-phi(i)) 
!!$               potential = potential + phi(1)*phi(2)*phi(3)
!!$            enddo
!!$         enddo
!!$      else
!!$         do i = 1, n_phases
!!$            potential = potential + g_eps22(i,lp)*(0.5*phi(i)**2*(1-phi(i))+phi(i)*(phi(lp)-1.5*phi(lp)**2)) 
!!$            potential = potential + phi(1)*phi(2)*phi(3)/abs(phi(lp)+1d-10)
!!$         enddo
!!$      endif
!!$   endif
!!$endif
!!$else !toth et al potential
!!$   if(lp.le.n_phases)then
!!$      potential0 = 1d0/12d0
!!$      do i=1,n_phases
!!$         potential0 = potential0 + phi(i)**4/4d0-phi(i)**3/3d0
!!$      enddo
!!$      do i=2,n_phases
!!$         do j=1,i-1
!!$            potential0 = potential0 + 0.5*phi(i)**2*phi(j)**2
!!$         enddo
!!$      enddo
!!$      phi_=phi+1d-3
!!$      eps=0d0
!!$      s=(phi_(1)*phi_(2))**2+(phi_(2)*phi_(3))**2+(phi_(3)*phi_(1))**2
!!$      eps=eps+g_eps22(1,2)*(phi_(1)*phi_(2))**2
!!$      eps=eps+g_eps22(2,3)*(phi_(2)*phi_(3))**2
!!$      eps=eps+g_eps22(3,1)*(phi_(3)*phi_(1))**2
!!$      eps=eps/s
!!$
!!$      if(lp.eq.0)then
!!$         potential = eps*potential0
!!$         return
!!$      endif
!!$      
!!$      potential = phi(lp)**3-phi(lp)**2
!!$      S=0d0
!!$      do i=1,n_phases
!!$         S = S + phi(i)**2
!!$      enddo
!!$      potential =  potential + phi(lp)*(S-phi(lp))
!!$      
!!$      
!!$   
!!$      if(lp.eq.1)d_eps=(2*phi_(1)*phi_(2)**2*phi_(3)**2*(g_eps22(1, 2)*phi_(2)**2-g_eps22(2, 3)*phi_(2)**2-g_eps22(2, 3)*phi_(3)**2+g_eps22(3, 1)*phi_(3)**2))/s**2
!!$      if(lp.eq.2)d_eps=(2*phi_(2)*phi_(3)**2*phi_(1)**2*(g_eps22(2, 3)*phi_(3)**2-g_eps22(3, 1)*phi_(3)**2-g_eps22(3, 1)*phi_(1)**2+g_eps22(1, 2)*phi_(1)**2))/s**2
!!$      if(lp.eq.3)d_eps=(2*phi_(3)*phi_(1)**2*phi_(2)**2*(g_eps22(3, 1)*phi_(1)**2-g_eps22(1, 2)*phi_(1)**2-g_eps22(1, 2)*phi_(2)**2+g_eps22(2, 3)*phi_(2)**2))/s**2
!!$      
!!$      potential=eps*potential + d_eps*potential0
!!$      
!!$      
!!$   else
!!$      write(*,*)"in potential: lp > n_phases?"
!!$      stop
!!$   endif
!!$   endif
end function potential

!test


double precision function PhaseRHS(vbles,dx,lp,Phi,c,T)
use solution_parameters
implicit none
double precision, intent(in)::dx,c,T,Phi(n_phases),vbles(N_unk+1,3,3)
integer, intent(in):: lp
double precision GradientEnergy,potential,FreeEnergy,PhiSoluteRHS,GradientEnergy_toth,GradientEnergyDelta
double precision ::M_tilde 
double precision :: c_dot=0,phi_dot(n_phases)
double precision :: MolPerVol,beta

phi_dot=0d0
!MolPerVol=(1.-c)*g_MolPerVol(1)+c*g_MolPerVol(2)

 !assemble and non-dimensionlise phaseRHS 

!M_tilde = (1.-c)+c*g_M(2)/g_M(1)
M_tilde =  1d0!is inverse propto Length scale
if(.true.)then
   phaseRHS=- M_tilde*(GradientEnergyDelta(vbles,lp,dx)+potential(Phi,c,lp)&
      +FreeEnergy(c,T,Phi,c_dot,phi_dot,lp)/delta_f)
elseif(folch)then
   phaseRHS=- M_tilde*(GradientEnergy_toth(vbles,lp,dx)+potential(Phi,c,lp)/g_lambda**2&
      +FreeEnergy(c,T,Phi,c_dot,phi_dot,lp)*MolPerVol/(g_W(1)*g_lambda**2))
else
   phaseRHS=- M_tilde*(GradientEnergy(vbles,lp,dx)+potential(Phi,c,lp)/g_lambda**2&
      +FreeEnergy(c,T,Phi,c_dot,phi_dot,lp)*MolPerVol/(g_W(1)*g_lambda**2))
endif

end function PhaseRHS






!


subroutine get_RHS_Fi_2d(i,j,k,lb,dx,Fi)
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use physicaldata
  use tree
  implicit none

integer, intent(in)::i,j,k,lb
double precision, intent(in)::dx
double precision, intent(out):: Fi(N_unk)
double precision :: c,T,phi(n_phases),c_dot,phi_dot(n_phases)
double precision phi_rhs(N_phases),rfactor,rf4,rf5,rf6,S
double precision TemperatureRHS,SoluteRHS,T_k,phi_
double precision vbles(N_unk+1,3,3),unk_(M_unk),MolPerVol,temp_vbles(N_unk+1,3,3)
logical Not_number
integer ip,n,ii,jj,kk

  Fi=0.

 do ii=1,3
     do jj=1,3
        do ip=1,N_unk
           vbles(ip,ii,jj)=unk(1+(ip-1)*nunkvbles,i+ii-2,j+jj-2,1,lb)
        enddo
     enddo
  enddo


  do ii=1,M_unk
     unk_(ii) = unk(ii,i,j,k,lb)
  enddo


  call get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_)
 
end subroutine get_RHS_Fi_2d

subroutine get_RHS_Fi(vbles,dx,Fi,mode_1,nunkvbles,dt,dtold,lnblocks,unk_)

  use solution_parameters

  implicit none

!
integer, intent(in)::mode_1,nunkvbles,lnblocks
double precision, intent(in)::dx,dt,dtold,unk_(M_unk)
double precision, intent(out):: Fi(N_unk)
double precision, intent(inout):: vbles(N_unk+1,3,3)
double precision :: c,T,phi(n_phases),c_dot,phi_dot(n_phases),phi_rhs(n_phases)
double precision rfactor,rf4,rf5,rf6,S,GradientEnergy,potential,Projector(n_phases,n_phases)
double precision TemperatureRHS,SoluteRHS,PhaseRHS,T_k,phi_,FreeEnergy,MolPerVol,Sum
logical Not_number
integer ip,n,i,j


Fi=0d0

  if (mode_1.eq.1) then
     do ip=1,n_phases
        n=(ip-1)*nunkvbles
        phi_dot(ip) = (vbles(ip,2,2)-unk_(n+2))/dt
     enddo

     n=(N_unk-1)*nunkvbles !puts pointer at beginning of solute part
     c_dot = (vbles(N_unk,2,2)-unk_(2+n))/dt
  elseif(mode_1.eq.2)then
     rfactor = dt/dtold

     rf4 = rfactor*rfactor/(rfactor+1.0)
     rf5 = (rfactor+1.0)
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
     do ip=1,n_phases
        n=(ip-1)*nunkvbles
        phi_dot(ip) = (rf4*unk_(n+5)-rf5*unk_(n+4)+rf6*vbles(ip,2,2))/dt
     enddo

     n=(N_unk-1)*nunkvbles
     c_dot = (rf4*unk_(5+n)-rf5*unk_(4+n)+rf6*vbles(N_unk,2,2))/dt
  else
     write(*,*)"mode_1 not = 1 or 2"
     stop
  endif

  sum=0d0
  do ip=1,n_phases
     phi(ip)=max(0d0,min(1d0,vbles(ip,2,2)))
     sum=sum+phi(ip)
  enddo
  do ip=1,n_phases
     phi(ip)=phi(ip)/sum
  enddo


  T = 0d0!vbles(N_unk-1,2,2)

  c = vbles(N_unk,2,2)

  


  call Project_Phase(vbles,dx,phi_rhs,Phi,c,T)
  sum=0d0
  do ip=1,n_phases
     Fi(ip) = phi_rhs(ip)
!     sum =sum + Fi(ip)
  enddo



  Fi(N_unk) = SoluteRHS(vbles,dx,c)

  do ip=1,N_unk
     n=(ip-1)*nunkvbles
     if(mode_1.le.1)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)!-unk_(n+2)
     else if(mode_1.eq.2)then
        Fi(ip)=vbles(ip,2,2)-dt*Fi(ip)/rf6!-unk_(n+2)
     end if
  enddo





end subroutine get_RHS_Fi

subroutine Calc_GradPhi(LapPhi,dphi,dx)
use solution_parameters
implicit none
double precision, intent(in)::dx,dphi(3,3)
double precision, intent(out)::LapPhi
double precision, dimension(3,3)::stencil=reshape((/1.,4.,1.,4.,-20.,4.,1.,4.,1./),(/3,3/))
double precision, dimension(3,3):: stencil5=reshape((/0.,1.,0.,1.,-4.,1.,0.,1.,0./),(/3,3/))
integer ip,ii,jj

 LapPhi=0.
 do ii=1,3
    do jj=1,3
       LapPhi=LapPhi+stencil(ii,jj)*dphi(ii,jj)/(6.*dx*dx)
    end do
 end do
  
end subroutine Calc_GradPhi!



!!$double precision function EntropyMixing(c,d)
!!$double precision, intent (in)::c
!!$integer, intent (in)::d
!!$
!!$if(c.ge.1d0.or.c.le.0d0)then
!!$   write(*,*)"c=",c, "in EntropyMixing"
!!$   stop
!!$endif
!!$ 
!!$ if(d.eq.1)then
!!$  EntropyMixing = log(c/(1.-c))
!!$ else
!!$  EntropyMixing = c*log(c)+(1.-c)*log(1.-c)
!!$ endif
!!$end function EntropyMixing


subroutine GibbsVector(V,T,d)
double precision, intent (out)::V(8)
double precision, intent (in)::T
integer, intent (in)::d
if(d.eq.0)then
 V(1)=1.0
 V(2)=T
 V(3)=T*log(T)
 V(4)=T**2
 V(5)=T**3
 V(6)=1./T
 V(7)=T**7
 V(8)=1.0/T**9
else if(d.eq.1)then
 V(1)=0.
 V(2)=1.
 V(3)=1.+log(T)
 V(4)=2.*T
 V(5)=3.*T**2
 V(6)=-1./T**2
 V(7)=7.*T**6
 V(8)=-9.0/T**10
elseif(d.eq.2)then
 V(1)=0.
 V(2)=0.
 V(3)=1./T
 V(4)=2.
 V(5)=6.*T
 V(6)=2./T**3
 V(7)=42.*T**5
 V(8)=90./T**11
else
 write(*,*)"bad diff order in GibbsVector"
stop
endif
end subroutine GibbsVector
subroutine MixingVector(RK,c,d)
double precision, intent (out)::RK(5)
double precision, intent (in)::c
integer, intent (in)::d
integer i
if(d.eq.0)then
 RK(1)=c*(1.-c)
 RK(2)=RK(1)*(1.-2.*c)
 RK(3)=RK(2)*(1.-2.*c)
 RK(4)=RK(3)*(1.-2.*c)
 RK(5)=RK(4)*(1.-2.*c)
elseif(d.eq.1)then
 RK(1)=1.-2.*c
 RK(2)=1.+(-6.+6.*c)*c
 RK(3)=1.+(-10.+(24.-16.*c)*c)*c
 RK(4)=1.+(-14.+(54.+(-80.+40.*c)*c)*c)*c
 RK(5)=1.+(-18.+(96.+(-224.+(240.-96.*c)*c)*c)*c)*c
elseif(d.eq.2)then
 RK(1)=-2.
 RK(2)=12.*c-6. 
 RK(3)=-10.+(48.-48.*c)*c
 RK(4)=-14.+(108.+(-240.+160.*c)*c)*c
 RK(5)=-18.+(192.+(-672.+(960.-480.*c)*c)*c)*c
elseif(d.eq.3)then
 RK(1)=0.
 RK(2)=12.
 RK(3)=-96.*c+48.
 RK(4)=108.+(-480.+480.*c)*c
 RK(5)=192.+(-1344.+(2880.-1920.*c)*c)*c
else
 write(*,*)"bad diff order in MixingVector"
 stop
endif
end subroutine MixingVector


double precision function g_phi(phi,d)
  double precision, intent (in):: phi
  integer, intent (in):: d
  if(d.eq.1)then
   g_phi = 6.*phi*(1.-phi)
!    g_phi = 30.*phi**4-60.*phi**3+30.*phi**2
  elseif(d.eq.2)then
   g_phi = 6.-12.*phi
!    g_phi = 120.*phi**3-180.*phi**2+60.*phi
  else !assume 0
   g_phi = 3.*phi**2-2.*phi**3
!    g_phi = 6.*phi**5-15.*phi**4+10.*phi**3
  endif
end function g_phi


logical function Not_number(x)
double precision, intent(in):: x
if(x.eq.x)then
Not_number=.false.
else
Not_number=.true.
endif
end function Not_number

double precision function f_c(x,df)
double precision, intent (in)::x
integer, intent(in)::df

f_c=0.
if(df.eq.0)then
if(x*(1-x).gt.0.)f_c=x*(1.-x)
elseif(df.eq.1)then
   write(*,*)"in f_c with df=1"
   stop
! if(x*(1-x).gt.0.)
f_c=1.-2.*x
endif
end function f_c

double precision function FreeEnergy(c,T,phi,c_dot,phi_dot,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi(n_phases),c_dot,phi_dot(n_phases)
integer, intent(in)::lp
double precision :: PbSnFreeEnergy,AlNiFreeEnergy,AlSiFreeEnergy,ToyFreeEnergy,AlFeFreeEnergy
double precision :: PbSnFreeEnergyKKS,PbSnFreeEnergyWBM ,PbSnFreeEnergyBJM
if(g_model.eq.0)then
   FreeEnergy = PbSnFreeEnergyWBM(c,T,phi,c_dot,phi_dot,lp)
elseif(g_model.eq.1)then
   FreeEnergy = AlNiFreeEnergy(c,T,phi,c_dot,phi_dot,lp)
elseif(g_model.eq.2)then
   FreeEnergy = AlSiFreeEnergy(c,T,phi,c_dot,phi_dot,lp)
elseif(g_model.eq.3)then
   FreeEnergy = ToyFreeEnergy(c,T,phi,c_dot,phi_dot,lp)
elseif(g_model.eq.4)then
   FreeEnergy = AlFeFreeEnergy(c,T,phi,c_dot,phi_dot,lp)
elseif(g_model.eq.5)then
   FreeEnergy = PbSnFreeEnergyKKS(c,T,phi,c_dot,phi_dot,lp)
elseif(g_model.eq.6)then
   FreeEnergy = PbSnFreeEnergyBJM(c,T,phi,c_dot,phi_dot,lp)
   return
else
   write(*,*)"g_model is currently 0 to 5"
   stop
endif

end function FreeEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!1

double precision function Gibbs_FE_e_AlSi(c,T,d_c,d_T)
  use solution_parameters
  use time_dep_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5), gs1A,gs1B, gs2A,gs2B

  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gs2_,Tk,Gibbs_T
  double precision :: hlA_(8),hlB_(8),hs1A_(8),hs1B_(8),hs2A_(8),hs2B_(8)
  double precision :: t_factor
  double precision :: Bezier_2,Quad
  integer :: T_low,T_high
  double precision :: T_weight


! AlSi below
!Si of the AlSi



Tk=T_k(T)
if (Tk .gt. 1687.0) then
  hs2A_=AlSihs2A3 
  hs2B_=AlSihs2B2 
elseif(Tk .gt. 933.47) then
  hs2A_=AlSihs2A3 
  hs2B_=AlSihs2B1 
elseif (Tk .gt. 700.0) then
  hs2A_=AlSihs2A2 
  hs2B_=AlSihs2B1 
else
  hs2A_=AlSihs2A1 
  hs2B_=AlSihs2B1 
endif

if(d_c.eq.0.and.d_T.eq.0)then
   gs2A=gibbs_T(hs2A_,Tk,0)
   gs2B=gibbs_T(hs2B_,Tk,0)
   call MixingVector(RK,c,0)
   RKl_=       (AlSihrks20(1)+AlSihrks20(2)*Tk)*RK(1)
   RKl_=RKl_ + (AlSihrks21(1)+AlSihrks21(2)*Tk)*RK(2)
   RKl_=RKl_ + (AlSihrks22(1)+AlSihrks22(2)*Tk)*RK(3)
   ent_=g_R*Tk*EntropyMixing(c,0)
   gs2_=gs2A*(1d0-c)+gs2B*c+ent_+RKl_
   Gibbs_FE_e_AlSi=gs2_
elseif(d_c.eq.0.and.d_T.eq.1)then
   gs2A=gibbs_T(hs2A_,Tk,1)
   gs2B=gibbs_T(hs2B_,Tk,1)
   call MixingVector(RK,c,0)
   RKl_=       AlSihrks20(2)*RK(1)
   RKl_=RKl_ + AlSihrks21(2)*RK(2)
   RKl_=RKl_ + AlSihrks22(2)*RK(3)
   ent_=g_R*EntropyMixing(c,0)
   gs2_=gs2A*(1d0-c)+gs2B*c+ent_+RKl_
   Gibbs_FE_e_AlSi=gs2_
elseif(d_c.eq.0.and.d_T.eq.2)then
   gs2A=gibbs_T(hs2A_,Tk,2)
   gs2B=gibbs_T(hs2B_,Tk,2)
   gs2_=gs2A*(1d0-c)+gs2B*c
   Gibbs_FE_e_AlSi=gs2_
elseif(d_c.eq.1.and.d_T.eq.0)then
   gs2A=gibbs_T(hs2A_,Tk,0)
   gs2B=gibbs_T(hs2B_,Tk,0)
   call MixingVector(RK,c,1)
   RKl_=       (AlSihrks20(1)+AlSihrks20(2)*Tk)*RK(1)
   RKl_=RKl_ + (AlSihrks21(1)+AlSihrks21(2)*Tk)*RK(2)
   RKl_=RKl_ + (AlSihrks22(1)+AlSihrks22(2)*Tk)*RK(3)
   ent_=g_R*Tk*EntropyMixing(c,1)
   gs2_=-gs2A+gs2B+ent_+RKl_
   Gibbs_FE_e_AlSi=gs2_
elseif(d_c.eq.1.and.d_T.eq.1)then
   gs2A=gibbs_T(hs2A_,Tk,1)
   gs2B=gibbs_T(hs2B_,Tk,1)
   call MixingVector(RK,c,1)
   RKl_=       AlSihrks20(2)*RK(1)
   RKl_=RKl_ + AlSihrks21(2)*RK(2)
   RKl_=RKl_ + AlSihrks22(2)*RK(3)
   ent_=g_R*EntropyMixing(c,1)
   gs2_=-gs2A+gs2B+ent_+RKl_
   Gibbs_FE_e_AlSi=gs2_

elseif(d_c.eq.1.and.d_T.eq.2)then
   gs2A=gibbs_T(hs2A_,Tk,2)
   gs2B=gibbs_T(hs2B_,Tk,2)
   gs2_=-gs2A+gs2B
   Gibbs_FE_e_AlSi=gs2_
elseif(d_c.eq.2.and.d_T.eq.0)then
   call MixingVector(RK,c,2)
   RKl_=       (AlSihrks20(1)+AlSihrks20(2)*Tk)*RK(1)
   RKl_=RKl_ + (AlSihrks21(1)+AlSihrks21(2)*Tk)*RK(2)
   RKl_=RKl_ + (AlSihrks22(1)+AlSihrks22(2)*Tk)*RK(3)
   ent_=g_R*Tk*EntropyMixing(c,2)
   Gibbs_FE_e_AlSi=ent_+RKl_
else
   write(*,*)"unexpected d_c,d_T combination",d_c,d_T
   stop
endif

end function Gibbs_FE_e_AlSi

!!!!!!!!!!!!!!!!!!!!!!!!!1
!Al (of the alsi)
double precision function Gibbs_FE_e2_AlSi(c,T,d_c,d_T)
  use solution_parameters
  use time_dep_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5), gs1A,gs1B, gs2A,gs2B

  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gs2_,Tk,Gibbs_T
  double precision :: hlA_(8),hlB_(8),hs1A_(8),hs1B_(8),hs2A_(8),hs2B_(8)
  double precision :: t_factor
  double precision :: Bezier_2,Quad
  integer :: T_low,T_high
  double precision :: T_weight


! AlSi below




Tk=T_k(T)
if (Tk .gt. 1687.0) then
  hs2A_=AlSihs1A3 
  hs2B_=AlSihs1B2 
elseif(Tk .gt. 933.47) then
  hs2A_=AlSihs1A3 
  hs2B_=AlSihs1B1 
elseif (Tk .gt. 700.0) then
  hs2A_=AlSihs1A2 
  hs2B_=AlSihs1B1 
else
  hs2A_=AlSihs1A1 
  hs2B_=AlSihs1B1 
endif

if(d_c.eq.0.and.d_T.eq.0)then
   gs2A=gibbs_T(hs2A_,Tk,0)
   gs2B=gibbs_T(hs2B_,Tk,0)
   call MixingVector(RK,c,0)
   RKl_=       (AlSihrks10(1)+AlSihrks10(2)*Tk)*RK(1)
   RKl_=RKl_ + (AlSihrks11(1)+AlSihrks11(2)*Tk)*RK(2)
   RKl_=RKl_ + (AlSihrks12(1)+AlSihrks12(2)*Tk)*RK(3)
   ent_=g_R*Tk*EntropyMixing(c,0)
   gs2_=gs2A*(1d0-c)+gs2B*c+ent_+RKl_
   Gibbs_FE_e2_AlSi=gs2_
elseif(d_c.eq.0.and.d_T.eq.1)then
   gs2A=gibbs_T(hs2A_,Tk,1)
   gs2B=gibbs_T(hs2B_,Tk,1)
   call MixingVector(RK,c,0)
   RKl_=       AlSihrks10(2)*RK(1)
   RKl_=RKl_ + AlSihrks11(2)*RK(2)
   RKl_=RKl_ + AlSihrks12(2)*RK(3)
   ent_=g_R*EntropyMixing(c,0)
   gs2_=gs2A*(1d0-c)+gs2B*c+ent_+RKl_
   Gibbs_FE_e2_AlSi=gs2_
elseif(d_c.eq.0.and.d_T.eq.2)then
   gs2A=gibbs_T(hs2A_,Tk,2)
   gs2B=gibbs_T(hs2B_,Tk,2)
   gs2_=gs2A*(1d0-c)+gs2B*c
   Gibbs_FE_e2_AlSi=gs2_
elseif(d_c.eq.1.and.d_T.eq.0)then
   gs2A=gibbs_T(hs2A_,Tk,0)
   gs2B=gibbs_T(hs2B_,Tk,0)
   call MixingVector(RK,c,1)
   RKl_=       (AlSihrks10(1)+AlSihrks10(2)*Tk)*RK(1)
   RKl_=RKl_ + (AlSihrks11(1)+AlSihrks11(2)*Tk)*RK(2)
   RKl_=RKl_ + (AlSihrks12(1)+AlSihrks12(2)*Tk)*RK(3)
   ent_=g_R*Tk*EntropyMixing(c,1)
   gs2_=-gs2A+gs2B+ent_+RKl_
   Gibbs_FE_e2_AlSi=gs2_
elseif(d_c.eq.1.and.d_T.eq.1)then
   gs2A=gibbs_T(hs2A_,Tk,1)
   gs2B=gibbs_T(hs2B_,Tk,1)
   call MixingVector(RK,c,1)
   RKl_=       AlSihrks10(2)*RK(1)
   RKl_=RKl_ + AlSihrks11(2)*RK(2)
   RKl_=RKl_ + AlSihrks12(2)*RK(3)
   ent_=g_R*EntropyMixing(c,1)
   gs2_=-gs2A+gs2B+ent_+RKl_
   Gibbs_FE_e2_AlSi=gs2_
elseif(d_c.eq.1.and.d_T.eq.2)then
   gs2A=gibbs_T(hs2A_,Tk,2)
   gs2B=gibbs_T(hs2B_,Tk,2)
   gs2_=-gs2A+gs2B
   Gibbs_FE_e2_AlSi=gs2_
elseif(d_c.eq.2.and.d_T.eq.0)then
   call MixingVector(RK,c,2)
   RKl_=       (AlSihrks10(1)+AlSihrks10(2)*Tk)*RK(1)
   RKl_=RKl_ + (AlSihrks11(1)+AlSihrks11(2)*Tk)*RK(2)
   RKl_=RKl_ + (AlSihrks12(1)+AlSihrks12(2)*Tk)*RK(3)
   ent_=g_R*Tk*EntropyMixing(c,2)
   Gibbs_FE_e2_AlSi=ent_+RKl_
else
   write(*,*)"unexpected d_c,d_T combination",d_c,d_T
   stop
endif

end function Gibbs_FE_e2_AlSi

!end AlSi


!!!!!!!!!!!!!!!
double precision function Gibbs_FE_l_AlSi(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: glA,glB,RK(5), gs1A,gs1B, gs2A,gs2B
  double precision :: EntropyMixing,T_k
  double precision :: RKl_,ent_,gl_,Tk,Gibbs_T
  double precision :: hlA_(8),hlB_(8),hs1A_(8),hs1B_(8),hs2A_(8),hs2B_(8)
  double precision :: Bezier_2,Quad
  integer :: T_low,T_high
  double precision :: T_weight


! AlSi below


Tk=T_k(T)


if (Tk .gt. 1687.0) then
  hlA_=AlSihlA3 
  hlB_=AlSihlB2 
elseif(Tk .gt. 933.47) then
  hlA_=AlSihlA3 
  hlB_=AlSihlB1 
elseif (Tk .gt. 700.0) then
  hlA_=AlSihlA2 
  hlB_=AlSihlB1 
else
  hlA_=AlSihlA1 
  hlB_=AlSihlB1 
endif

if(d_c.eq.0.and.d_T.eq.0)then
   glA=gibbs_T(hlA_,Tk,0)
   glB=gibbs_T(hlB_,Tk,0)
   call MixingVector(RK,c,0)
   RKl_=         (AlSihrkl0(1)+AlSihrkl0(2)*Tk)*RK(1)
   RKl_=RKl_ + (AlSihrkl1(1)+AlSihrkl1(2)*Tk)*RK(2)
   RKl_=RKl_ + (AlSihrkl2(1)+AlSihrkl2(2)*Tk)*RK(3)
   
   ent_=g_R*Tk*EntropyMixing(c,0)
   gl_=glA*(1d0-c)+glB*c+ent_+RKl_
   Gibbs_FE_l_AlSi=gl_
elseif(d_c.eq.0.and.d_T.eq.1)then
   glA=gibbs_T(hlA_,Tk,1)
   glB=gibbs_T(hlB_,Tk,1)
   call MixingVector(RK,c,0)
   RKl_=         AlSihrkl0(2)*RK(1)
   RKl_=RKl_ + AlSihrkl1(2)*RK(2)
   RKl_=RKl_ + AlSihrkl2(2)*RK(3)
   
   ent_=g_R*EntropyMixing(c,0)
   gl_=glA*(1d0-c)+glB*c+ent_+RKl_

   Gibbs_FE_l_AlSi=gl_
elseif(d_c.eq.0.and.d_T.eq.2)then
   glA=gibbs_T(hlA_,Tk,2)
   glB=gibbs_T(hlB_,Tk,2)
   gl_=glA*(1d0-c)+glB*c
   Gibbs_FE_l_AlSi=gl_
elseif(d_c.eq.1.and.d_T.eq.0)then
   glA=gibbs_T(hlA_,Tk,0)
   glB=gibbs_T(hlB_,Tk,0)
   call MixingVector(RK,c,1)
   RKl_=         (AlSihrkl0(1)+AlSihrkl0(2)*Tk)*RK(1)
   RKl_=RKl_ + (AlSihrkl1(1)+AlSihrkl1(2)*Tk)*RK(2)
   RKl_=RKl_ + (AlSihrkl2(1)+AlSihrkl2(2)*Tk)*RK(3)
   
   ent_=g_R*Tk*EntropyMixing(c,1)
   gl_=-glA+glB+ent_+RKl_
   Gibbs_FE_l_AlSi=gl_
elseif(d_c.eq.1.and.d_T.eq.1)then
   glA=gibbs_T(hlA_,Tk,1)
   glB=gibbs_T(hlB_,Tk,1)
   call MixingVector(RK,c,1)
   RKl_=       AlSihrkl0(2)*RK(1)
   RKl_=RKl_ + AlSihrkl1(2)*RK(2)
   RKl_=RKl_ + AlSihrkl2(2)*RK(3)
   
   ent_=g_R*EntropyMixing(c,1)
   gl_=-glA+glB+ent_+RKl_
   Gibbs_FE_l_AlSi=gl_
elseif(d_c.eq.1.and.d_T.eq.2)then
   glA=gibbs_T(hlA_,Tk,2)
   glB=gibbs_T(hlB_,Tk,2)
   gl_=-glA+glB
   Gibbs_FE_l_AlSi=gl_
elseif(d_c.eq.2.and.d_T.eq.0)then
   call MixingVector(RK,c,2)
   RKl_=         (AlSihrkl0(1)+AlSihrkl0(2)*Tk)*RK(1)
   RKl_=RKl_ + (AlSihrkl1(1)+AlSihrkl1(2)*Tk)*RK(2)
   RKl_=RKl_ + (AlSihrkl2(1)+AlSihrkl2(2)*Tk)*RK(3)
   
   ent_=g_R*Tk*EntropyMixing(c,2)
   Gibbs_FE_l_AlSi=ent_+RKl_
else
   write(*,*)"unexpected d_c,d_T combination",d_c,d_T
   stop
endif

end function Gibbs_FE_l_AlSi

!!!!!!!!!!!!AlFe
!!!!!!!!!!!!!!!!!!!!!!!!!1

double precision function Gibbs_FE_e_AlFe(c,T,d_c,d_T)
  use solution_parameters
  use time_dep_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: EntropyMixing,T_k,RK(5)
  double precision :: RKa,ent,ga,gaA,gaB,Tk,Gibbs_T
  double precision :: haA(8),haB(8)




! AlFe below
Tk=T_k(T)
if (Tk .gt. 1811.0)then
  haA=AlFehaA3 
  haB=AlFehaB2 
elseif (Tk .gt. 933.47)then
  haA=AlFehaA3 
  haB=AlFehaB1 
elseif (Tk .gt. 700.0)then
  haA=AlFehaA2 
  haB=AlFehaB1 
else
  haA=AlFehaA1 
  haB=AlFehaB1 
endif







if(d_c.eq.0.and.d_T.eq.0)then
   gaA=gibbs_T(haA,Tk,0)
   gaB=gibbs_T(haB,Tk,0)
   call MixingVector(RK,c,0)
   RKa=      (AlFehrka0(1)+AlFehrka0(2)*Tk)*RK(1)
   RKa=RKa + (AlFehrka1(1)+AlFehrka1(2)*Tk)*RK(2)
   ent = g_R*Tk*EntropyMixing(c,0)
   ga = gaA*(1d0-c)+gaB*c+ent+RKa
   Gibbs_FE_e_AlFe=ga
elseif(d_c.eq.0.and.d_T.eq.1)then
   gaA=gibbs_T(haA,Tk,1)
   gaB=gibbs_T(haB,Tk,1)
   call MixingVector(RK,c,0)
   RKa=      AlFehrka0(2)*RK(1)
   RKa=RKa + AlFehrka1(2)*RK(2)
   ent = g_R*EntropyMixing(c,0)
   ga = gaA*(1d0-c)+gaB*c+ent+RKa
   Gibbs_FE_e_AlFe=ga
elseif(d_c.eq.0.and.d_T.eq.2)then
elseif(d_c.eq.1.and.d_T.eq.0)then
   gaA=gibbs_T(haA,Tk,0)
   gaB=gibbs_T(haB,Tk,0)
   call MixingVector(RK,c,1)
   RKa=      (AlFehrka0(1)+AlFehrka0(2)*Tk)*RK(1)
   RKa=RKa + (AlFehrka1(1)+AlFehrka1(2)*Tk)*RK(2)
   ent = g_R*Tk*EntropyMixing(c,1)
   ga = -gaA + gaB+ent+RKa
   Gibbs_FE_e_AlFe=ga
elseif(d_c.eq.1.and.d_T.eq.1)then
elseif(d_c.eq.1.and.d_T.eq.2)then
elseif(d_c.eq.2.and.d_T.eq.0)then
else
   write(*,*)"unexpected d_c,d_T combination",d_c,d_T
   stop
endif

end function Gibbs_FE_e_AlFe

!!!!!!!!!!!!!!!!!!!!!!!!!1
!Al (of the AlFe)
double precision function Gibbs_FE_e2_AlFe(c,T,d_c,d_T)
  use solution_parameters
  use time_dep_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: EntropyMixing,T_k,RK(5)
  double precision :: Tk,Gibbs_T,entSLb,dentSLb,AlFey1,y1,dy1
  double precision :: hslb1(8),hslb2(8),sLb1,dsLb1,sLb2,dsLb2,gb1,gb,dgb1,dgb


! AlFe below
Gibbs_FE_e2_AlFe=0



Tk=T_k(T)
if (Tk .gt. 1811.0)then
  hslb1=AlFehslb1_4 
  hslb2=AlFehslb2_4 
elseif (Tk .gt. 933.47)then
  hslb1=AlFehslb1_3 
  hslb2=AlFehslb2_3 
elseif (Tk .gt. 700.0)then
  hslb1=AlFehslb1_2 
  hslb2=AlFehslb2_2 
else
  hslb1=AlFehslb1_1 
  hslb2=AlFehslb2_1 
endif

if(d_c.eq.0.and.d_T.eq.0)then
   sLb1=gibbs_T(hslb1,Tk,0)
   sLb2=gibbs_T(hslb2,Tk,0)
   y1 = AlFey1(max(min(0.27246,c),0.2350),0)
   entSLb=0.1375*g_R*Tk*EntropyMixing(y1,0)
   gb1=SLb1*y1 + SLb2*(1d0-y1) 
   gb=(gb1 + entSLb)/(0.8625 + 0.1375*y1)
   Gibbs_FE_e2_AlFe=gb
   return


elseif(d_c.eq.0.and.d_T.eq.1)then
   sLb1=gibbs_T(hslb1,Tk,1)
   sLb2=gibbs_T(hslb2,Tk,1)
   y1 = AlFey1(max(min(0.27246,c),0.2350),0)
   entSLb=0.1375*g_R*EntropyMixing(y1,0)
   gb1=SLb1*y1 + SLb2*(1d0-y1) 
   gb=(gb1 + entSLb)/(0.8625 + 0.1375*y1)
   Gibbs_FE_e2_AlFe=gb
   return
elseif(d_c.eq.0.and.d_T.eq.2)then
elseif(d_c.eq.1.and.d_T.eq.0)then
   sLb1=gibbs_T(hslb1,Tk,0)
   sLb2=gibbs_T(hslb2,Tk,0)
!   dsLb1=0d0
!   dsLb2=0d0
   y1 = AlFey1(max(min(0.27246,c),0.2350),0)
   dy1 = AlFey1(max(min(0.27246,c),0.2350),1)
   entSLb = 0.1375*g_R*Tk*EntropyMixing(y1,0)
   dentSLb = 0.1375*g_R*Tk*EntropyMixing(y1,1)
   gb1 =SLb1*y1 + SLb2*(1d0-y1) 
   dgb1 = (SLb1 - SLb2)
   dgb = -0.1375*(gb1 + entSLb)/(0.8625 + 0.1375*y1)**2 + (dgb1 + dentSLb)/(0.8625 + 0.1375*y1)
   Gibbs_FE_e2_AlFe = dgb*dy1
   return
elseif(d_c.eq.1.and.d_T.eq.1)then
elseif(d_c.eq.1.and.d_T.eq.2)then
elseif(d_c.eq.2.and.d_T.eq.0)then
else
   write(*,*)"unexpected d_c,d_T combination",d_c,d_T
   stop
endif

end function Gibbs_FE_e2_AlFe

!end AlFe


!!!!!!!!!!!!!!!
double precision function Gibbs_FE_l_AlFe(c,T,d_c,d_T)
  use solution_parameters
  implicit none
  double precision, intent(in) :: c,T
  integer,intent(in):: d_c,d_T
  double precision :: EntropyMixing,T_k,RK(5)
  double precision :: RKl,ent,gl,glA,glB,Tk,Gibbs_T,Rla
  double precision :: hlA(8),hlB(8)

! AlFe below
Gibbs_FE_l_AlFe=0

Tk=T_k(T)
if (Tk .gt. 1811.0)then
  hlA=AlFehlA3 
  hlB=AlFehlB2 
elseif (Tk .gt. 933.47)then
  hlA=AlFehlA3 
  hlB=AlFehlB1 
elseif (Tk .gt. 700.0)then
  hlA=AlFehlA2 
  hlB=AlFehlB1 
else
  hlA=AlFehlA1 
  hlB=AlFehlB1 
endif


if(d_c.eq.0.and.d_T.eq.0)then
   glA=gibbs_T(hlA,Tk,0)
   glB=gibbs_T(hlB,Tk,0)
   call MixingVector(RK,c,0)
   RKl=       (AlFehrkl0(1)+AlFehrkl0(2)*Tk)*RK(1)
   RKl=RKl + (AlFehrkl1(1)+AlFehrkl1(2)*Tk)*RK(2)
   RKl=RKL + (AlFehrkl2(1)+AlFehrkl2(2)*Tk)*RK(3)

   ent = g_R*Tk*EntropyMixing(c,0)
   gl = glA*(1d0-c)+glB*c+ent+RKl
   Gibbs_FE_l_AlFe=gl
   return
elseif(d_c.eq.0.and.d_T.eq.1)then
   glA=gibbs_T(hlA,Tk,1)
   glB=gibbs_T(hlB,Tk,1)
   call MixingVector(RK,c,0)
   RKl=      AlFehrkl0(2)*RK(1)
   RKl=RKl + AlFehrkl1(2)*RK(2)
   RKl=RKL + AlFehrkl2(2)*RK(3)
   ent = g_R*EntropyMixing(c,0)
   gl = glA*(1d0-c)+glB*c+ent+RKl
   Gibbs_FE_l_AlFe=gl
   return
elseif(d_c.eq.0.and.d_T.eq.2)then
elseif(d_c.eq.1.and.d_T.eq.0)then
   glA=gibbs_T(hlA,Tk,0)
   glB=gibbs_T(hlB,Tk,0)
   call MixingVector(RK,c,1)
   RKl=       (AlFehrkl0(1)+AlFehrkl0(2)*Tk)*RK(1)
   RKl=RKl + (AlFehrkl1(1)+AlFehrkl1(2)*Tk)*RK(2)
   RKl=RKL + (AlFehrkl2(1)+AlFehrkl2(2)*Tk)*RK(3)

   ent = g_R*Tk*EntropyMixing(c,1)
   gl = -glA + glB+ent+RKl
   Gibbs_FE_l_AlFe=gl

elseif(d_c.eq.1.and.d_T.eq.1)then
elseif(d_c.eq.1.and.d_T.eq.2)then
elseif(d_c.eq.2.and.d_T.eq.0)then
else
   write(*,*)"unexpected d_c,d_T combination",d_c,d_T
   stop
endif

end function Gibbs_FE_l_AlFe


!!!!!!!!!!!!!!!
double precision function ToyFreeEnergyPhase(c,phi,lp)
use solution_parameters
implicit none
double precision, intent (in):: c,phi(n_phases)
integer, intent(in)::lp
double precision :: S,ff1,ff2,ff3,g_phi,c1,c2,c3
integer :: i

c1=0.5
c2=0.25
c3=0.75

   S=0d0
   do i=1,n_phases
      S=S+g_phi(phi(i),0)
   enddo
   
   ff1 = (c-c1 +(c1-c2)*g_phi(phi(2),0)/S+(c1-c3)*g_phi(phi(3),0)/S)**2+g_phi(phi(1),0)*g_T0
   ff2 = (c-c2 +(c2-c3)*g_phi(phi(3),0)/S+(c2-c1)*g_phi(phi(1),0)/S)**2+g_phi(phi(1),0)*g_T0
   ff3 = (c-c3 +(c3-c1)*g_phi(phi(1),0)/S+(c3-c2)*g_phi(phi(2),0)/S)**2+g_phi(phi(1),0)*g_T0
!!$
!!$      ff1 = (c + 0.25 *g_phi(phi(2),0)/S&
!!$           - 0.25*g_phi(phi(3),0)/S  - 0.5)**2 + g_T0&
!!$           - g_T0*g_phi(phi(2),0)/S - g_T0*g_phi(phi(3),0)/S
!!$      
!!$      ff1 = (
!!$
!!$      ff2 = (c - 0.50*g_phi(phi(3),0)/S&
!!$           - 0.25*g_phi(phi(1),0)/S - 0.25)**2&
!!$           + g_T0 *g_phi(phi(1),0)/S
!!$      
!!$      ff3 = (c + 0.25*g_phi(phi(1),0)/S&
!!$           + 0.50*g_phi(phi(2),0)/S - 0.75)**2&
!!$           + g_T0*g_phi(phi(1),0)/S
      
      ToyFreeEnergyPhase =                  ff1*g_phi(phi(1),0)/S
      ToyFreeEnergyPhase = ToyFreeEnergyPhase +  ff2*g_phi(phi(2),0)/S
      ToyFreeEnergyPhase = ToyFreeEnergyPhase +  ff3*g_phi(phi(3),0)/S

   


end function ToyFreeEnergyPhase

double precision function ToyFreeEnergy     (c,T,phi,c_dot,phi_dot,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi(n_phases),c_dot,phi_dot(n_phases)
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(n_phases),hb(n_phases),hr(n_phases),f_phi,g_phi,EntropyMixing
double precision :: Tk,T_K,Cp,R
double precision :: FE,S,dS,Gibbs_Fe_l_AlSi,Gibbs_FE_e_AlSi,Gibbs_FE_e2_AlSi,ff1,ff2,ff3,dff1,dff2,dff3
double precision :: dff11,dff12,dff13,dff21,dff22,dff23,dff31,dff32,dff33
double precision :: ToyFreeEnergyPhase,h=1d-6,phi_(3),c_
integer :: i,j
logical :: plapp_PCB=.true.
ToyFreeEnergy=0d0




if(lp.eq.0)then !free energy proper J/mol 
if(plapp_PCB)then

   ToyFreeEnergy = ToyFreeEnergyPhase(c,phi,lp)
else
   S=0d0
   do i=1,n_phases
      S=S+g_phi(phi(i),0)
   enddo
   ToyFreeEnergy =                  ((c-0.5)**2+g_T0)*g_phi(phi(1),0)/S
   ToyFreeEnergy = ToyFreeEnergy + (c-0.25)**2*g_phi(phi(2),0)/S
   ToyFreeEnergy = ToyFreeEnergy + (c-0.75)**2*g_phi(phi(3),0)/S 

endif   

elseif(lp.le.n_phases)then !phase = current_var

if(plapp_PCB)then
   if(phi(lp).gt.0.5)then
      h=-abs(h)
   else
      h=abs(h)
   endif
   phi_ = phi
   phi_(lp)=phi(lp)+h
   ToyFreeEnergy = (ToyFreeEnergyPhase(c,phi_,lp) - ToyFreeEnergyPhase(c,phi,lp))/h

else
  
  S=0d0
   do i=1,n_phases
      S=S+g_phi(phi(i),0)
   enddo
   if(lp.eq.1)then
      ToyFreeEnergy =  ((c-0.5)**2+g_T0)*g_phi(phi(lp),1)/S
   elseif(lp.eq.2)then
      ToyFreeEnergy =  (c-0.25)**2*g_phi(phi(lp),1)/S
   elseif(lp.eq.3)then
      ToyFreeEnergy =  (c-0.75)**2*g_phi(phi(lp),1)/S
   endif
   
   
   dS = -1d0/S**2*g_phi(phi(lp),1)
   ToyFreeEnergy = ToyFreeEnergy+ ((c-0.5)**2+g_T0)*g_phi(phi(1),0)*dS
   ToyFreeEnergy = ToyFreeEnergy+       (c-0.25)**2*g_phi(phi(2),0)*dS
   ToyFreeEnergy = ToyFreeEnergy+       (c-0.75)**2*g_phi(phi(3),0)*dS
   
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




   return



elseif(lp.eq.N_unk)then !Solute
if(plapp_PCB)then
   if(c.gt.0.5)then
      h=-abs(h)
   else
      h=abs(h)
   endif
   c_ = c + h
   ToyFreeEnergy = (ToyFreeEnergyPhase(c_,phi,lp) - ToyFreeEnergyPhase(c,phi,lp))/h

else
   S=0d0
   do i=1,n_phases
      S=S+g_phi(phi(i),0)
   enddo
   ToyFreeEnergy =                 2*(c-0.5) *g_phi(phi(1),0)/S
   ToyFreeEnergy = ToyFreeEnergy + 2*(c-0.25)*g_phi(phi(2),0)/S
   ToyFreeEnergy = ToyFreeEnergy + 2*(c-0.75)*g_phi(phi(3),0)/S

endif

return
elseif(lp.eq.N_unk-1)then !Temperature Keep in J/mol

   ToyFreeEnergy =  0
   return
else
   write(*,*)"bad current_var"
   stop
endif
end function ToyFreeEnergy

!!!!!!!!!!!!!!!

double precision function AlSiFreeEnergy     (c,T,phi,c_dot,phi_dot,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi(n_phases),c_dot,phi_dot(n_phases)
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(n_phases),hb(n_phases),hr(n_phases),f_phi,g_phi,EntropyMixing
double precision :: Tk,T_K,Cp,R
double precision :: FE,S,dS,Gibbs_Fe_l_AlSi,Gibbs_FE_e_AlSi,Gibbs_FE_e2_AlSi
integer :: i,j

AlSiFreeEnergy=0d0




if(lp.eq.0)then !free energy proper J/mol 

   S=0d0
   do i=1,n_phases
      S=S+g_phi(phi(i),0)
   enddo
   AlSiFreeEnergy =                  Gibbs_Fe_l_AlSi(c,T,0,0)*g_phi(phi(1),0)/S
   AlSiFreeEnergy = AlSiFreeEnergy + Gibbs_FE_e_AlSi(c,T,0,0)*g_phi(phi(2),0)/S
   AlSiFreeEnergy = AlSiFreeEnergy + Gibbs_FE_e2_AlSi(c,T,0,0)*g_phi(phi(3),0)/S !3rd phase debug
   
   return

elseif(lp.le.n_phases)then !phase = current_var

  
  S=0d0
   do i=1,n_phases
      S=S+g_phi(phi(i),0)
   enddo
   if(lp.eq.1)then
      AlSiFreeEnergy =  Gibbs_Fe_l_AlSi(c,T,0,0)*g_phi(phi(lp),1)/S
   elseif(lp.eq.2)then
      AlSiFreeEnergy =  Gibbs_FE_e_AlSi(c,T,0,0)*g_phi(phi(lp),1)/S
   elseif(lp.eq.3)then
      AlSiFreeEnergy =  Gibbs_FE_e2_AlSi(c,T,0,0)*g_phi(phi(lp),1)/S !3rd phase to implement debug
   endif
   
   
   dS = -1d0/S**2*g_phi(phi(lp),1)
   AlSiFreeEnergy = AlSiFreeEnergy+ Gibbs_Fe_l_AlSi(c,T,0,0)*g_phi(phi(1),0)*dS
   AlSiFreeEnergy = AlSiFreeEnergy+ Gibbs_Fe_e_AlSi(c,T,0,0)*g_phi(phi(2),0)*dS
   AlSiFreeEnergy = AlSiFreeEnergy+ Gibbs_Fe_e2_AlSi(c,T,0,0)*g_phi(phi(3),0)*dS


   return



elseif(lp.eq.N_unk)then !Solute
   S=0d0
   do i=1,n_phases
      S=S+g_phi(phi(i),0)
   enddo
   AlSiFreeEnergy =                  Gibbs_Fe_l_AlSi(c,T,1,0)*g_phi(phi(1),0)/S
   AlSiFreeEnergy = AlSiFreeEnergy + Gibbs_FE_e_AlSi(c,T,1,0)*g_phi(phi(2),0)/S
   AlSiFreeEnergy = AlSiFreeEnergy + Gibbs_FE_e2_AlSi(c,T,1,0)*g_phi(phi(3),0)/S !should be the 3rd phase debug
   return

elseif(lp.eq.N_unk-1)then !Temperature Keep in J/mol

   AlSiFreeEnergy =  0
   return
      Tk=T_k(T)
      
      !-(1-T*d/dT)dF/dphi
      !     if(heat_phi_term)then
      !  T (d/dT)dF/dphi
      
      AlSiFreeEnergy =  Tk*(Gibbs_Fe_l_AlSi(c,T,0,1) - Gibbs_FE_e_AlSi(c,T,0,1))

         ! -Df/Dphi
         AlSiFreeEnergy = AlSiFreeEnergy - Gibbs_Fe_l_AlSi(c,T,0,0) + Gibbs_FE_e_AlSi(c,T,0,0)
         g_latent_heat = AlSiFreeEnergy         
         AlSiFreeEnergy = AlSiFreeEnergy*phi_dot(1)*g_phi(phi(1),1)
         
         Q_heat(1) = AlSiFreeEnergy
 !     endif !end phi_dot_term
      if(heat_c_term)then 
         !   -(1- T(d/dT))dF/dc 
         FE=(Tk*Gibbs_Fe_l_AlSi(c,T,1,1)*g_phi(phi(1),0)- Gibbs_Fe_l_AlSi(c,T,1,0))*g_phi(phi(1),0)
         FE=FE+(Tk*Gibbs_FE_e_AlSi(c,T,1,1)- Gibbs_FE_e_AlSi(c,T,1,0))*(g_phi(phi(2),0))
         FE=FE*c_dot
         AlSiFreeEnergy = AlSiFreeEnergy + FE
         Q_heat(2)= FE

      endif !c_dot_term
      !   heat capacity = -T (d/dT) dF/dT
      g_C =  -Tk*(g_phi(phi(1),0)*Gibbs_Fe_l_AlSi(c,T,0,2)+g_phi(phi(2),0)*Gibbs_FE_e_AlSi(c,T,0,2))

      return

else
   write(*,*)"bad current_var"
   stop
endif
end function AlSiFreeEnergy


!!!!!!!!!!!!!!!!!!!

double precision function PbSnFreeEnergy(c,T,phi,c_dot,phi_dot,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi(n_phases),c_dot,phi_dot(n_phases)
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(n_phases),hb(n_phases),hr(n_phases),f_phi,g_phi,EntropyMixing
double precision :: Tk,T_K,Cp,R
double precision :: FE,S,dS,f(3)
integer :: i,j
double precision:: PbSn_i,a0,a1,a2,h=0.0001,c0=0.01,c1=.99
PbSnFreeEnergy=0d0



R=g_R
Tk=T_K(T)

if(lp.eq.0)then !free energy proper J/mol





   PbSnFreeEnergy=0.

   S=0d0
   do i=1,n_phases
      S=S+g_phi(phi(i),0)
   enddo
   do i=1,n_phases
      f(i) = PbSn_i(c,i,0)
      PbSnFreeEnergy = PbSnFreeEnergy + f(i)*g_phi(phi(i),0)/S
   enddo




   return
elseif(lp.le.n_phases)then !phase = current_var

   if(zero_start)return
   call GibbsVector(V,Tk,0)
   call MixingVector(RK,c,0)
   PbSnFreeEnergy=0.
   ha=0.
   hb=0.
   hr=0.
   do i=1,n_phases
      do j = 1,8
         ha(i) = ha(i) + hsa(i,j)*V(j)
         hb(i) = hb(i) + hsb(i,j)*V(j)
      enddo
   enddo
   do i=1,n_phases
      do j = 1,4
         hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
      enddo
   enddo

   S=0d0
   do i=1,n_phases
      S=S+g_phi(phi(i),0)
   enddo
   PbSnFreeEnergy = (c*hb(lp)+(1.-c)*ha(lp)+hr(lp))*g_phi(phi(lp),1)/S
   do i = 1, n_phases
      dS = -1d0/S**2*g_phi(phi(lp),1)
      PbSnFreeEnergy = PbSnFreeEnergy+ (c*hb(i)+(1.-c)*ha(i)+hr(i))*g_phi(phi(i),0)*dS
   enddo
  
 return
elseif(lp.eq.N_unk)then !Solute

   PbSnFreeEnergy=R*Tk*EntropyMixing(c,1)
   if(zero_start)return
   call GibbsVector(V,Tk,0)
   call MixingVector(RK,c,1)

   ha=0.
   hb=0.
   hr=0.
   do i=1,n_phases
      do j = 1,8
         ha(i) = ha(i) + hsa(i,j)*V(j)
         hb(i) = hb(i) + hsb(i,j)*V(j)
      enddo
   enddo
   do i=1,n_phases
      do j = 1,4
         hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
      enddo
   enddo
   S=0d0
   do i=1,n_phases
      S=S+g_phi(phi(i),0)
   enddo
   do i=1,n_phases
      PbSnFreeEnergy = PbSnFreeEnergy + (hb(i) - ha(i)+hr(i))*g_phi(phi(i),0)/S
   enddo

   return

elseif(lp.eq.N_unk-1)then !Temperature Keep in J/mol


   PbSnFreeEnergy = 0d0
   FE = 0d0
   if(heat_phi_term)then

!  T (d/dT)dF/dphi
      call GibbsVector(V,Tk,1) !first derivative of free energy
      call MixingVector(RK,c,0)      
      ha=0.
      hb=0.
      hr=0.
      do i=1,n_phases
         do j = 1,8
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      do i=1,n_phases
         do j = 1,4
            hr(i) = hr(i) + hsrk(i, 2*j)*RK(j)
         enddo!
      enddo

      do i=1,n_phases

         FE = FE + (c*hb(i)+(1.-c)*ha(i)+hr(i))*g_phi(phi(i),1)*phi_dot(i)*Tk

      enddo

      
      
      !    pure -dF/dphi term  
      call GibbsVector(V,Tk,0)
      
      
      ha=0.
      hb=0.
      hr=0.
      do i=1,n_phases
         do j = 1,8
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      do i=1,n_phases
         do j = 1,4
            hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
         enddo
      enddo
      do i=1,n_phases
         FE = FE - (c*hb(i)+(1.-c)*ha(i)+hr(i))*g_phi(phi(i),1)*phi_dot(i)

      enddo
      Q_heat(1)=FE
      heat_size(1)=max(heat_size(1),abs(FE))
      PbSnFreeEnergy=FE
   endif !end phi_dot_term
   if(heat_c_term)then 
      !    T(d/dT)dF/dc =-T (ds/dc) ! see AT14.pdf for details

      call GibbsVector(V,Tk,1)
      call MixingVector(RK,c,1)      
      FE=0d0
      ha=0.
      hb=0.
      hr=0.
      do i=1,n_phases
         do j = 1,8
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      do i=1,n_phases
         do j = 1,4
            hr(i) = hr(i) + hsrk(i, 2*j)*RK(j)
         enddo
      enddo
      do i=1,n_phases
         FE = FE + (hb(i)-ha(i)+hr(i))*g_phi(phi(i),0)*Tk*c_dot 
      enddo


      call GibbsVector(V,Tk,0)
      ha=0.
      hb=0.
      hr=0.
      do i=1,n_phases
         do j = 1,8
            ha(i) = ha(i) + hsa(i,j)*V(j)
            hb(i) = hb(i) + hsb(i,j)*V(j)
         enddo
      enddo
      do i=1,n_phases
         do j = 1,4
            hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
         enddo
      enddo
      do i=1,n_phases
         FE = FE - (hb(i)-ha(i)+hr(i))*g_phi(phi(i),0)*c_dot
      enddo
      
      heat_size(2)=max(heat_size(2),abs(FE))
      Q_heat(2)=FE
      PbSnFreeEnergy=PbSnFreeEnergy + FE
   endif !c_dot_term
!   heat capacity = -T (d/dT) dF/dT
   call GibbsVector(V,Tk,2)

   g_C=0d0
   ha=0.
   hb=0.
   do j = 1,8
      do i=1,n_phases
         ha(i) = ha(i) + hsa(i,j)*V(j)
         hb(i) = hb(i) + hsb(i,j)*V(j)
      enddo
   enddo
!   g_C(1) = g_C(1) - ha(i)*g_phi(phi,0) - ha(j)*(1.-g_phi(phi,0))
!   g_C(2) = g_C(2) - hb(i)*g_phi(phi,0) - hb(j)*(1.-g_phi(phi,0))
   do i=1,n_phases
      g_C(1) = g_C(1) - ha(i)*g_phi(phi(i),0)
      g_C(2) = g_C(2) - hb(i)*g_phi(phi(i),0)
   enddo
   g_C=g_C*Tk

   return
else
   write(*,*)"bad current_var"
   stop
endif
end function PbSnFreeEnergy

!!!!!!!!!!!!!!!!!!!

double precision function PbSnFreeEnergyKKS(c,T,phi,c_dot,phi_dot,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi(n_phases),c_dot,phi_dot(n_phases)
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(n_phases),hb(n_phases),hr(n_phases),f_phi,g_phi,EntropyMixing,f(3),df(3)
double precision :: Tk,T_K,Cp,R,tt
double precision :: FE,S,dS,cc(3),PbSn_i,a0,a1,a2,h=0.001,S1,S2,S3,aleft,accel=0.0
integer :: i,j



S=0d0
do i=1,n_phases
   S=S+g_phi(phi(i),0)
enddo

do i=1,n_phases
   cc(i) = c + g_aa(i)
   do j=1,n_phases
      cc(i)=cc(i)-g_phi(phi(j),0)*g_aa(j)/S
!      cc(i)=cc(i)-phi(j)*g_aa(j)
   enddo
enddo
R=g_R
Tk=T_K(T)

if(lp.eq.0)then !free energy proper J/mol

   PbSnFreeEnergyKKS=0.


   do i=1,n_phases
      if(cc(i).lt.g_bb(i))then
         f(i) = PbSn_i(g_bb(i),i,0)+PbSn_i(g_bb(i),i,1)*(cc(i)-g_bb(i)) + 0.5*PbSn_i(g_bb(i),i,2)*(cc(i)-g_bb(i))**2*accel
      elseif(cc(i).gt.g_aa(3))then
         f(i) = PbSn_i(g_aa(3),i,0)+PbSn_i(g_aa(3),i,1)*(cc(i)-g_aa(3))+0.5*PbSn_i(g_aa(3),i,2)*(cc(i)-g_aa(3))**2*accel
      else
         f(i)=PbSn_i(cc(i),i,0)
      endif
   enddo

   do i=1,n_phases
      PbSnFreeEnergyKKS = PbSnFreeEnergyKKS + f(i)*g_phi(phi(i),0)/S
   enddo

   return
elseif(lp.le.n_phases)then !phase = current_var
!
! d/d phi_k (g(phi_i)/S * f_i(c_i)) = g'(phi_k)/S^2(f_k - sum_i g(phi_i)f_i/S) - g_aa(k) sum_ig(phi_i) f_i'(c_i)/S
!                                   = S2         (f(k) - S1         ) - g_aa(k)* S3/S

! n.b. d/d c(c_i) = 1
!      d/d phi_k((c_i)= - g_aa(k)
   

   do i=1,n_phases
      if(cc(i).lt.g_bb(i))then
         f(i) = PbSn_i(g_bb(i),i,0)+PbSn_i(g_bb(i),i,1)*(cc(i)-g_bb(i)) + 0.5*PbSn_i(g_bb(i),i,2)*(cc(i)-g_bb(i))**2*accel
      elseif(cc(i).gt.g_aa(3))then
         f(i) = PbSn_i(g_aa(3),i,0)+PbSn_i(g_aa(3),i,1)*(cc(i)-g_aa(3)) + 0.5*PbSn_i(g_aa(3),i,2)*(cc(i)-g_aa(3))**2*accel
      else
         f(i)=PbSn_i(cc(i),i,0)
      endif
   enddo
   S1 =0.
   do i=1,n_phases
      S1 = S1 + f(i)*g_phi(phi(i),0)
   enddo
   
   S2 = g_phi(phi(lp),1)/S
   
   do i=1,n_phases
      if(cc(i).lt.g_bb(i))then
         df(i) = PbSn_i(g_bb(i),i,1) + PbSn_i(g_bb(i),i,2)*(cc(i)-g_bb(i))*accel
      elseif(cc(i).gt.g_aa(3))then
         df(i) = PbSn_i(g_aa(3),i,1) + PbSn_i(g_aa(3),i,2)*(cc(i)-g_aa(3))*accel
      else
         df(i)=PbSn_i(cc(i),i,1)
      endif
   enddo

   S3 = 0.
   do i=1,n_phases
      S3 = S3 +df(i)*g_phi(phi(i),0)
   enddo
   

   PbSnFreeEnergyKKS = S2*(f(lp) - S1/S + S3/S*(c-cc(lp)))

!!$   PbSnFreeEnergyKKS = S2*(f(lp) - S1/S) - g_aa(lp)*S3/S
!!$   do i=1,n_phases
!!$      S3 = S3 + df(i)*g_phi(phi(i),0)
!!$   enddo

 return
elseif(lp.eq.N_unk)then !Solute




   do i=1,n_phases
      if(cc(i).lt.g_bb(i))then
         df(i) = PbSn_i(g_bb(i),i,1)+PbSn_i(g_bb(i),i,2)*(cc(i)-g_bb(i))*accel
      elseif(cc(i).gt.g_aa(3))then
         df(i) = PbSn_i(g_aa(3),i,1) + PbSn_i(g_aa(3),i,2)*(cc(i)-g_aa(3))*accel
      else
         df(i)=PbSn_i(cc(i),i,1)
      endif
   enddo

   PbSnFreeEnergyKKS=0.
   do i=1,n_phases
      PbSnFreeEnergyKKS = PbSnFreeEnergyKKS + df(i)*g_phi(phi(i),0)/S
   enddo

   return
else
   write(*,*)"bad current_var"
   stop
endif
end function PbSnFreeEnergyKKS


!!!!!!!!!!!!!!!!!!!

double precision function PbSnFreeEnergyWBM(c,T,phi,c_dot,phi_dot,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi(n_phases),c_dot,phi_dot(n_phases)
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(n_phases),hb(n_phases),hr(n_phases),f_phi,g_phi,EntropyMixing,f(3),df(3)
double precision :: Tk,T_K,Cp,R,tt
double precision :: FE,S,dS,cc(3),PbSn_i,a0,a1,a2,h=0.001,S1,S2,S3,aleft,accel=0.0
integer :: i,j



S=0d0
do i=1,n_phases
   S=S+g_phi(phi(i),0)
enddo

do i=1,n_phases
   cc(i) = c ! temporary fix. We really need to change "cc(i)" to "c" 
enddo


R=g_R
Tk=T_K(T)

if(lp.eq.0)then !free energy proper J/mol

   PbSnFreeEnergyWBM=0.


   do i=1,n_phases
      if(cc(i).lt.g_bb(i))then
         f(i) = PbSn_i(g_bb(i),i,0)+PbSn_i(g_bb(i),i,1)*(cc(i)-g_bb(i)) + 0.5*PbSn_i(g_bb(i),i,2)*(cc(i)-g_bb(i))**2*accel
      elseif(cc(i).gt.g_aa(3))then
         f(i) = PbSn_i(g_aa(3),i,0)+PbSn_i(g_aa(3),i,1)*(cc(i)-g_aa(3))+0.5*PbSn_i(g_aa(3),i,2)*(cc(i)-g_aa(3))**2*accel
      else
         f(i)=PbSn_i(cc(i),i,0)
      endif
   enddo

   do i=1,n_phases
      PbSnFreeEnergyWBM = PbSnFreeEnergyWBM + f(i)*g_phi(phi(i),0)/S
   enddo

   return
elseif(lp.le.n_phases)then !phase = current_var
!
! d/d phi_k (g(phi_i)/S * f_i(c_i)) = g'(phi_k)/S^2(f_k - sum_i g(phi_i)f_i/S) - g_aa(k) sum_ig(phi_i) f_i'(c_i)/S
!                                   = S2         (f(k) - S1         ) - g_aa(k)* S3/S

! n.b. d/d c(c_i) = 1
!      d/d phi_k((c_i)= - g_aa(k)
   

   do i=1,n_phases
      if(cc(i).lt.g_bb(i))then
         f(i) = PbSn_i(g_bb(i),i,0)+PbSn_i(g_bb(i),i,1)*(cc(i)-g_bb(i)) + 0.5*PbSn_i(g_bb(i),i,2)*(cc(i)-g_bb(i))**2*accel
      elseif(cc(i).gt.g_aa(3))then
         f(i) = PbSn_i(g_aa(3),i,0)+PbSn_i(g_aa(3),i,1)*(cc(i)-g_aa(3)) + 0.5*PbSn_i(g_aa(3),i,2)*(cc(i)-g_aa(3))**2*accel
      else
         f(i)=PbSn_i(cc(i),i,0)
      endif
   enddo
   S1 =0.
   do i=1,n_phases
      S1 = S1 + f(i)*g_phi(phi(i),0)
   enddo
   
   S2 = g_phi(phi(lp),1)/S
   
   do i=1,n_phases
      if(cc(i).lt.g_bb(i))then
         df(i) = PbSn_i(g_bb(i),i,1) + PbSn_i(g_bb(i),i,2)*(cc(i)-g_bb(i))*accel
      elseif(cc(i).gt.g_aa(3))then
         df(i) = PbSn_i(g_aa(3),i,1) + PbSn_i(g_aa(3),i,2)*(cc(i)-g_aa(3))*accel
      else
         df(i)=PbSn_i(cc(i),i,1)
      endif
   enddo

   S3 = 0.
   do i=1,n_phases
      S3 = S3 +df(i)*g_phi(phi(i),0)
   enddo
   

   PbSnFreeEnergyWBM = S2*(f(lp) - S1/S + S3/S*(c-cc(lp)))

!!$   PbSnFreeEnergyWBM = S2*(f(lp) - S1/S) - g_aa(lp)*S3/S
!!$   do i=1,n_phases
!!$      S3 = S3 + df(i)*g_phi(phi(i),0)
!!$   enddo

 return
elseif(lp.eq.N_unk)then !Solute




   do i=1,n_phases
      if(cc(i).lt.g_bb(i))then
         df(i) = PbSn_i(g_bb(i),i,1)+PbSn_i(g_bb(i),i,2)*(cc(i)-g_bb(i))*accel
      elseif(cc(i).gt.g_aa(3))then
         df(i) = PbSn_i(g_aa(3),i,1) + PbSn_i(g_aa(3),i,2)*(cc(i)-g_aa(3))*accel
      else
         df(i)=PbSn_i(cc(i),i,1)
      endif
   enddo

   PbSnFreeEnergyWBM=0.
   do i=1,n_phases
      PbSnFreeEnergyWBM = PbSnFreeEnergyWBM + df(i)*g_phi(phi(i),0)/S
   enddo

   return
else
   write(*,*)"bad current_var"
   stop
endif
end function PbSnFreeEnergyWBM



!!!!!!!!!!!!!!!!!!!


! function change from c to c1,c2, c3.
! d=0 the basic function
! d ne. 0 . j=0 is c derivative and j-1,2,3 is derivatives wrt phi(3)
double precision function c_interpolant(c,phi,i,j,d)
double precision, intent (in)::c,phi(3)
integer, intent (in)::d,i,j
double precision ::  delta(3) = (/0d0,2.462193664,-2.462193664/),dd,dd_(3),g_phi,g(3),s,h=1d-6,phi_(3),ddh



s = g_phi(phi(1),0)+g_phi(phi(2),0)+g_phi(phi(3),0)

g(1)=g_phi(phi(1),0)/s
g(2)=g_phi(phi(2),0)/s
g(3)=g_phi(phi(3),0)/s


dd = g(1)*delta(1)+g(2)*delta(2)+g(3)*delta(3)


dd_(1) =exp(delta(1)-dd)
dd_(2) =exp(delta(2)-dd)
dd_(3) =exp(delta(3)-dd)



if(d.eq.0)then
   c_interpolant = c/((1-c)*dd_(i)+c)
elseif(j.eq.0)then !c derivative
   c_interpolant = dd_(i)/((1-c)*dd_(i)+c)**2
else !phi derivative is j index
   c_interpolant = (1-c)*c*delta(j)*dd_(i)/((1-c)*dd_(i)+c)**2
endif





end function c_interpolant


! function change from c to c1,c2, c3.
! d=0 the basic function
! d ne. 0 . j=0 is c derivative and j-1,2,3 is derivatives wrt phi(3)
double precision function c_interpolantg(c,phi,i,j,d)
double precision, intent (in)::c,phi(3)
integer, intent (in)::d,i,j
double precision ::  delta(3) = (/0d0,2.462193664,-2.462193664/),dd,dd_(3)
dd = phi(1)*delta(1)+phi(2)*delta(2)+phi(3)*delta(3)
dd_(1) =exp(delta(1)-dd)
dd_(2) =exp(delta(2)-dd)
dd_(3) =exp(delta(3)-dd)


if(d.eq.0)then
   c_interpolantg = c/((1-c)*dd_(i)+c)
elseif(j.eq.0)then !c derivative
   c_interpolantg = dd_(i)/((1-c)*dd_(i)+c)**2
else !phi derivative is j index
   c_interpolantg = (1-c)*c*delta(j)*dd_(i)/((1-c)*dd_(i)+c)**2
endif




end function c_interpolantg

double precision function hc_slope(c,d,i)
double precision, intent (in)::c
integer, intent (in)::d,i
![1.692915096*c + 0.2703149202, 1.692915096*c, 1.692915096*c]
if(d.eq.0)then
   hc_slope= 1.692915096*c
   if(i.eq.1) hc_slope = hc_slope + 0.2703149202
elseif(d.eq.1)then
   hc_slope= 1.692915096
else
   write(*,*)"hc_slope. d>1"
   stop
endif
end function hc_slope


double precision function PbSnFreeEnergyBJM(c,T,phi,c_dot,phi_dot,lp)
use solution_parameters
implicit none
double precision, intent (in):: T,c,phi(n_phases),c_dot,phi_dot(n_phases)
integer, intent(in)::lp
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(n_phases),hb(n_phases),hr(n_phases),f_phi,g_phi,EntropyMixing,f(3),df(3)
double precision :: Tk,T_K,Cp,R,tt
double precision :: FE,S,dS,cc(3),dcc(3),cc_phi(3),PbSn_i,a0,a1,a2,h=0.001,S1,S2,S3,S4,aleft,accel=0.0,phi_(3),hh=1d-6
integer :: i,j
double precision :: c_interpolant, hc_slope,hc(3),dhc(3)


!work
!c_interpolant(c,phi,i,j,d),hc_slope (c,d,i)

S=0d0
do i=1,n_phases
   S=S+g_phi(phi(i),0)
enddo
!work
phi_ = phi
phi_(lp)=phi(lp)+h
do i=1,n_phases
   cc(i)  = c_interpolant(c,phi,i,0,0)
   dcc(i) = c_interpolant(c,phi,i,0,1)
!   cc_phi(i) = c_interpolant(c,phi,i,lp,1)
   cc_phi(i) = (c_interpolant(c,phi_,i,lp,0) - c_interpolant(c,phi,i,0,0))/h
enddo


R=g_R
Tk=T_K(T)

if(lp.eq.0)then !free energy proper J/mol

   PbSnFreeEnergyBJM=0.


   do i=1,n_phases
!!$      if(cc(i).lt.g_bb(i))then
!!$         f(i) = PbSn_i(g_bb(i),i,0)+PbSn_i(g_bb(i),i,1)*(cc(i)-g_bb(i)) + 0.5*PbSn_i(g_bb(i),i,2)*(cc(i)-g_bb(i))**2*accel
!!$      elseif(cc(i).gt.g_aa(3))then
!!$         f(i) = PbSn_i(g_aa(3),i,0)+PbSn_i(g_aa(3),i,1)*(cc(i)-g_aa(3))+0.5*PbSn_i(g_aa(3),i,2)*(cc(i)-g_aa(3))**2*accel
!!$      else
         f(i)=PbSn_i(cc(i),i,0)
!      endif
      f(i) = f(i)  - hc_slope (cc(i),0,i)  + hc_slope(c,0,i)
   enddo

   do i=1,n_phases
      PbSnFreeEnergyBJM = PbSnFreeEnergyBJM + f(i)*g_phi(phi(i),0)/S
   enddo



   return
elseif(lp.le.n_phases)then !phase = current_var


   do i=1,n_phases
!!$      if(cc(i).lt.g_bb(i))then
!!$         f(i) = PbSn_i(g_bb(i),i,0)+PbSn_i(g_bb(i),i,1)*(cc(i)-g_bb(i)) + 0.5*PbSn_i(g_bb(i),i,2)*(cc(i)-g_bb(i))**2*accel
!!$      elseif(cc(i).gt.g_aa(3))then
!!$         f(i) = PbSn_i(g_aa(3),i,0)+PbSn_i(g_aa(3),i,1)*(cc(i)-g_aa(3)) + 0.5*PbSn_i(g_aa(3),i,2)*(cc(i)-g_aa(3))**2*accel
!!$      else
         f(i)=PbSn_i(cc(i),i,0)
!      endif
      f(i) = f(i) - hc_slope(cc(i),0,i)  + hc_slope(c,0,i)
   enddo

   S1 =0.
   do i=1,n_phases
      S1 = S1 + f(i)*g_phi(phi(i),0)
   enddo
   
   S2 = g_phi(phi(lp),1)/S
   
   do i=1,n_phases
      if(cc(i).lt.g_bb(i))then
         df(i) = PbSn_i(g_bb(i),i,1) + PbSn_i(g_bb(i),i,2)*(cc(i)-g_bb(i))*accel
      elseif(cc(i).gt.g_aa(3))then
         df(i) = PbSn_i(g_aa(3),i,1) + PbSn_i(g_aa(3),i,2)*(cc(i)-g_aa(3))*accel
      else
         df(i)=PbSn_i(cc(i),i,1)
      endif
      df(i) = df(i)*cc_phi(i)  - hc_slope(cc(i),1,i)*cc_phi(i)
   enddo

   S3 = 0.
   do i=1,n_phases
      S3 = S3 +df(i)*g_phi(phi(i),0)
   enddo


   

   PbSnFreeEnergyBJM = S2*(f(lp) - S1/S + S3/S)  




 return
elseif(lp.eq.N_unk)then !Solute


   do i=1,n_phases
!!$      if(cc(i).lt.g_bb(i))then
!!$         df(i) = PbSn_i(g_bb(i),i,1)+PbSn_i(g_bb(i),i,2)*(cc(i)-g_bb(i))*accel
!!$      elseif(cc(i).gt.g_aa(3))then
!!$         df(i) = PbSn_i(g_aa(3),i,1) + PbSn_i(g_aa(3),i,2)*(cc(i)-g_aa(3))*accel
!!$      else
         df(i)=PbSn_i(cc(i),i,1)
!      endif
      df(i) = df(i)*dcc(i) - hc_slope(cc(i),1,i)*dcc(i) + hc_slope(c,1,i)
   enddo

   PbSnFreeEnergyBJM=0.
   do i=1,n_phases
      PbSnFreeEnergyBJM = PbSnFreeEnergyBJM + df(i)*g_phi(phi(i),0)/S
   enddo

   return
else
   write(*,*)"bad current_var"
   stop
endif
end function PbSnFreeEnergyBJM

double precision function PbSn_i(c,i,d)
use solution_parameters
implicit none
integer, intent(in):: i,d
double precision, intent(in)::c
integer:: j
double precision:: T_K,R,Tk
double precision, dimension (8):: V
double precision, dimension (5):: RK
double precision ha(n_phases),hb(n_phases),hr(n_phases),EntropyMixing
R=g_R
Tk=T_K(0.)

if(d.eq.0)then

   
   call GibbsVector(V,Tk,0)
   call MixingVector(RK,c,0)
   
   ha=0.
   hb=0.
   hr=0.
   
   do j = 1,8
      ha(i) = ha(i) + hsa(i,j)*V(j)
      hb(i) = hb(i) + hsb(i,j)*V(j)
   enddo
   do j = 1,4
      hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
   enddo
   PbSn_i = (c*hb(i)+(1.-c)*ha(i)+hr(i))+ R*Tk*EntropyMixing(c,0)
elseif(d.eq.1)then

   
   call GibbsVector(V,Tk,0)
   call MixingVector(RK,c,1)
   
   ha=0.
   hb=0.
   hr=0.
   
   do j = 1,8
      ha(i) = ha(i) + hsa(i,j)*V(j)
      hb(i) = hb(i) + hsb(i,j)*V(j)
   enddo
   do j = 1,4
      hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
   enddo
   PbSn_i = hb(i)-ha(i)+hr(i)+ R*Tk*EntropyMixing(c,1)
elseif(d.eq.2)then



   call MixingVector(RK,c,2)


   hr=0.

   do j = 1,4
      hr(i) = hr(i) + (hsrk(i,2*j-1)+hsrk(i, 2*j)*Tk)*RK(j)
   enddo
   PbSn_i = hr(i)+ R*Tk*EntropyMixing(c,2)
else
   write(*,*)"PbSn_i: derivative > 1 not coded"
   stop
endif

PbSn_i = PbSn_i/(g_R*g_T0)

end function PbSn_i
