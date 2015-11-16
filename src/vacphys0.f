!##############################################################################
! module vacphys0 - mhd

!=============================================================================
subroutine conserve(ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w)

! Transform primitive variables into conservative ones

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)
!-----------------------------------------------------------------------------


! Calculate total energy from pressure, kinetic and magnetic energy

w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,e_)=w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,p_)/(eqpar(gamma_)-1)+half*((w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_))&
   *(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,v1_)**2+w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,v2_)**2+w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,v3_)**2)+((w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b1_))&
   **2+(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b2_))**2&
   +(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_))**2))&
   +((w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b1_)*w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,bg1_))+(w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,b2_)*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg2_))&
   +(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_)*w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,bg3_)))


! Convert velocity to momentum
w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m1_)=(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,rhob_))*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,v1_)
w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m2_)=(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,rhob_))*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,v2_)
w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m3_)=(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,rhob_))*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,v3_);


return
end

!=============================================================================
subroutine primitive(ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w)

! Transform conservative variables into primitive ones

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)
!-----------------------------------------------------------------------------


! Calculate pressure

call getpthermal(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,tmp)

w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,p_)=tmp(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3)

! Convert momentum to velocity
w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,v1_)=w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,m1_)/(w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_))
w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,v2_)=w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,m2_)/(w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_))
w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,v3_)=w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,m3_)/(w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_));

return
end

!=============================================================================
subroutine getv(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,v)

! Calculate v_idim=m_idim/rho within ix

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw),&
   v(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
!-----------------------------------------------------------------------------

oktest=index(teststr,'getv')>=1
if(oktest)write(*,*)'GetV w:',w(ixtest1,ixtest2,ixtest3,iwtest)

v(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,m0_+idim)/(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)&
   +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_))

if(oktest)write(*,*)'GetV v:',v(ixtest1,ixtest2,ixtest3)

return 
end


!=============================================================================
subroutine getcmax(new_cmax,w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,&
   cmax)

! Calculate cmax_idim=cfast_i+abs(v_idim) within ix^L
! where cfast_i=sqrt(0.5*(cf**2+sqrt(cf**4-4*cs**2*b_i**2/rho)))
! and cf**2=b**2/rho+cs**2/rho is the square of the speed of the fast wave 
! perpendicular to the magnetic field, and cs is the sound speed.

include 'vacdef.f'

logical:: new_cmax
integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw),&
   cmax(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
double precision:: csound2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   cfast2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
save csound2,cfast2
!-----------------------------------------------------------------------------
oktest=index(teststr,'getcmax')>=1

!Direction independent part of getcmax:
if(new_cmax)then
   new_cmax=.false.
   call getcsound2(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,csound2)
   if(oktest)write(*,*)'csound2:',csound2(ixtest1,ixtest2,ixtest3)
   cfast2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=((w(ixmin1:ixmax1,&
      ixmin2:ixmax2,ixmin3:ixmax3,b1_)+w(ixmin1:ixmax1,ixmin2:ixmax2,&
      ixmin3:ixmax3,bg1_))**2+(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
      b2_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg2_))**2&
      +(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_)+w(ixmin1:ixmax1,&
      ixmin2:ixmax2,ixmin3:ixmax3,bg3_))**2 )/(w(ixmin1:ixmax1,ixmin2:ixmax2,&
      ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_))&
      +csound2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)
end if
if(oktest)write(*,*)'cfast2:',cfast2(ixtest1,ixtest2,ixtest3)

cmax(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=sqrt(half*(cfast2&
   (ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)+ sqrt(cfast2(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3)**2-4*csound2(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3)* ((w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b0_+idim)&
   +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg0_+idim))**2)&
   /(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,rhob_))))) +abs(w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,m0_+idim)/(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)&
   +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_)))

if(oktest) write(*,*)'cmax:',cmax(ixtest1,ixtest2,ixtest3)


return 
end

!=============================================================================
subroutine getcsound2prim(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,csound2)

! Calculate the square of the thermal sound speed csound2 within ix^L
! from the primitive variables in w.
! csound2=gamma*p/rho

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw),&
   csound2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)call die&
   ('Correct Getcsound2prim for NONIDEAL gas in vacphys.t.mhd')

csound2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=eqpar(gamma_)&
   *(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,p_)+(eqpar(gamma_)&
   -one)*(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,eb_)-half&
   *( (w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg1_))**2&
   +(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg2_))**2+(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,bg3_))**2 )))/(w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_))

return 
end

!=============================================================================
subroutine getcsound2(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,csound2)

! Calculate the square of the thermal sound speed csound2 within ix^L.
! csound2=gamma*p/rho

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw),&
   csound2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)call die&
   ('Correct Getcsound2 for NONIDEAL gas in vacphys.t.mhd')

oktest=index(teststr,'getcsound2')>=1
if(oktest) write(*,*)'Getcsound2'

call getpthermal(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,csound2)
if(oktest) write(*,*)'p(ixtest)=',csound2(ixtest1,ixtest2,ixtest3)
csound2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=eqpar(gamma_)&
   *(csound2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)+(eqpar(gamma_)&
   -one)*(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,eb_)-half&
   *( (w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg1_))**2&
   +(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg2_))**2+(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,bg3_))**2 )))/(w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_))

return 
end

!=============================================================================
subroutine getpthermal(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,p)

!!! This subroutine should not use tmp,tmp2


include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw),&
   p(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
!-----------------------------------------------------------------------------


p(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=half*( w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,m1_)**2+w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,m2_)**2+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m3_)&
   **2 )/(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)&
   +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_))

p(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=p(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3)+ half*( (w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b1_)&
   **2)+(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b2_)**2)&
   +(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_)**2) )&
   +( (w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b1_)*w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,bg1_))+(w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,b2_)*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg2_))&
   +(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_)*w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,bg3_)) )

p(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(eqpar(gamma_)&
   -one)*(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,e_)-p(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3))


return 
end

!=============================================================================
subroutine getptotal(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,p)

include 'vacdef.f'

double precision::  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw),&
   p(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),gamma
integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)call die&
   ('Correct GetPtotal for NONIDEAL gas in vacphys.t.mhd')

gamma=eqpar(gamma_)

p(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(gamma-two)*(( (w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,b1_)*w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,bg1_))+(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b2_)&
   *w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg2_))+(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,b3_)*w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,bg3_)) )+ half*( (w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,b1_))**2.d0+(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
   b2_))**2.d0+(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_))**2.d0 ))

p(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(gamma-one)*(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,e_)-half*( w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3,m1_)**2.d0+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m2_)&
   **2.d0+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m3_)**2.d0 )&
   /(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3,rhob_)))-p(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixmin3:ixmax3)


return 
end

!=============================================================================
subroutine getptotal_bg(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,p)

include 'vacdef.f'

double precision::  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw),&
   p(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),gamma
integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)call die&
   ('Correct GetPtotal for NONIDEAL gas in vacphys.t.mhd')

gamma=eqpar(gamma_)

p(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(eqpar(gamma_)&
   -one)*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,eb_)-half*(eqpar(gamma_)&
   -two)*( (w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg1_)**2.d0)&
   +(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg2_)**2.d0)&
   +(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg3_)**2.d0) )    

return 
end

