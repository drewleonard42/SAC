  !##############################################################################
  ! module vacusr - uniform spiral ! setvac -d=22 -g=204,204 -p=hdadiab -u=sim1
  
  
!==============================================================================
!
!    THE FOLLOWING SUBROUTINES ADD GRAVITATIONAL SOURCE TERMS, SET GRAVITY
!
!------------------------------------------------------------------------------
!    See vacusr.t.gravity and vacusrpar.t.gravity for an example of usage
!
!    Gravitational force is added to the momentum equation:
!
!    d m_i/dt += rho*eqpar(grav0_+i)
!
!    Gravitational work is added to the energy equation (if present):
!
!    de/dt += Sum_i m_i*eqpar(grav0_+i)
!
!    The eqpar(grav1_),eqpar(grav2_),... coefficients are the components of 
!    the gravitational acceleration in each dimension. Set them to 0 for no
!    gravity in that direction. 
!    The !!! comments show how a grav array could be used for a spatially
!    (and maybe temporally) varying gravitational field.
!    The setgrav subroutine has to be completed then.
!
!============================================================================
subroutine addsource_grav(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iws,qtC,w,qt,wnew)

! Add gravity source calculated from w to wnew within ixO for all variables 
! in iws. w is at time qtC, wnew is advanced from qt to qt+qdt.

include 'vacdef.f'

integer::          ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iws(niw_)
double precision:: qdt,qtC,qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw),&
   wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)
integer:: iiw,iw,idim
!!! ! For a spatially varying gravity define the common grav array
!!! double precision:: grav(ixG^T,ndim)
!!! common /gravity/ grav

!-----------------------------------------------------------------------------

!!! ! If grav needs to be calculated only once do it for the whole grid
!!! if(it==itmin)call setgrav(w,ixG^L,ixG^L,grav)
!!! ! Otherwise call setgrav in every time step
!!! call setgrav(w,ixI^L,ixO^L,grav)

! add sources from gravity
do iiw=1,iws(niw_); iw=iws(iiw)
   select case(iw)
   case(m1_,m2_,m3_)
      ! dm_i/dt= +rho*g_i
      idim=iw-m0_
      if(abs(eqpar(grav0_+idim))>smalldouble) wnew(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,m0_+idim)=wnew(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,m0_+idim)+ qdt*eqpar(grav0_&
         +idim)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_))

!          wnew(ixO^S,m0_+idim)=wnew(ixO^S,m0_+idim)+ &
!              qdt*eqpar(grav0_+idim)*(w(ixO^S,rho_)+w(ixO^S,rhob_))

      !!! ! For a spatially varying gravity use instead of the above lines
      !!! wnew(ixO^S,m0_+idim)=wnew(ixO^S,m0_+idim)+ &
      !!!    qdt*grav(ixO^S,idim)*(w(ixO^S,rho_)+w(ixO^S,rhob_))

   case(e_)
      ! de/dt= +g_i*m_i
      do idim=1,ndim
         if(abs(eqpar(grav0_+idim))>smalldouble) wnew(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,ee_)=wnew(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,ee_)+ qdt*eqpar(grav0_&
            +idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)&
            *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,m0_&
            +idim)/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)&
            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rhob_))

!            wnew(ixO^S,ee_)=wnew(ixO^S,ee_)+ &
!               qdt*eqpar(grav0_+idim)*w(ixO^S,m0_+idim)

         !!! ! For a spatially varying gravity use instead of the above lines
         !!! wnew(ixO^S,ee_)=wnew(ixO^S,ee_)+ &
         !!!    qdt*grav(ixO^S,idim)*w(ixO^S,m0_+idim)

      end do
   end select ! iw
end do        ! iiw

return
end
!=============================================================================
!!! subroutine setgrav(w,ixI^L,ixO^L,grav)

! Set the gravitational acceleration within ixO based on x(ixI,ndim) 
! and/or w(ixI,nw)

!!! include 'vacdef.f'

!!! double precision:: w(ixG^T,nw),grav(ixG^T,ndim)
!!! integer:: ixI^L,ixO^L
!----------------------------------------------------------------------------
!!! return
!!! end
!=============================================================================

subroutine getdt_grav(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3)

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)
integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim
double precision:: dtgrav
save dtgrav

!!! ! For spatially varying gravity you need a common grav array
!!! double precision:: grav(ixG^T,ndim)
!!! common/gravity/grav

!----------------------------------------------------------------------------

oktest=index(teststr,'getdt')>=1

if(it==itmin)then
   ! If gravity is descibed by the equation parameters, use this:
   dtgrav=bigdouble
   do idim=1,ndim
      if(abs(eqpar(grav0_+idim))>zero)dtgrav=min(dtgrav,one&
         /sqrt(maxval(abs(eqpar(grav0_+idim))/dx(ixMmin1:ixMmax1,&
         ixMmin2:ixMmax2,ixMmin3:ixMmax3,1:ndim))))
   enddo
   !!! ! For spatially varying gravity use this instead of the lines above:
   !!! call setgrav(w,ixG^L,ixM^L,grav)
   !!! ! If gravity does not change with time, calculate dtgrav here:
   !!! dtgrav=one/sqrt(maxval(abs(grav(ixM^S,1:ndim))/dx(ixM^S,1:ndim)))
endif

!!! ! If gravity changes with time, calculate dtgrav here:
!!! dtgrav=one/sqrt(maxval(abs(grav(ixM^S,1:ndim))/dx(ixM^S,1:ndim)))

      call mpiallreduce(dtgrav,MPI_MIN)

! limit the time step
dt=min(dt,dtgrav)
if(oktest)write(*,*)'Gravity limit for dt:',dtgrav

return
end

!=============================================================================

!==============================================================================
subroutine addsource_visc(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iws,qtC,w,qt,wnew)

! Add viscosity source to wnew within ixO 

include 'vacdef.f'

integer::          ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iws(niw_)
double precision:: qdt,qtC,qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw),&
   wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)

integer:: ix,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,idir,jdir,iiw,iw
double precision:: tmp2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   nushk(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,ndim)


double precision:: tmprhoL(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
    tmprhoR(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3), tmprhoC(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,ixGlo3:ixGhi3)
double precision:: tmpVL(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
    tmpVR(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3), tmpVC(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,ixGlo3:ixGhi3)
double precision:: tmpBL(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
    tmpBR(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3), tmpBC(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,ixGlo3:ixGhi3)

double precision:: tmpL(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   tmpR(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3), tmpC(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,ixGlo3:ixGhi3)

double precision:: nuL(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   nuR(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)

integer:: jxmin1,jxmin2,jxmin3,jxmax1,jxmax2,jxmax3,hxmin1,hxmin2,hxmin3,&
   hxmax1,hxmax2,hxmax3, hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3

double precision:: c_ene,c_shk

integer:: i,j,k,l,m,ii0,ii1,t00

double precision:: sB

!-----------------------------------------------------------------------------

! Calculating viscosity sources 
! involves second derivatives, two extra layers
call ensurebound(2,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qtC,w)
ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmin3=ixOmin3-1;ixmax1=ixOmax1+1
ixmax2=ixOmax2+1;ixmax3=ixOmax3+1;

!sehr wichtig
call setnushk(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,nushk)

do idim=1,ndim
      tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,rho_)
      call setnu(w,rho_,idim,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         nuR,nuL)      
      call gradient1L(tmp,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,tmp2)
      tmpL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
         =(nuL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
         +nushk(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim))&
         *tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)         
      call gradient1R(tmp,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,tmp2)
      tmpR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
         =(nuR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
         +nushk(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim))&
         *tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
      wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,rho_)&
         =wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,rho_)&
         +(tmpR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
         -tmpL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))&
         /dx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim)*qdt
enddo
   

do idim=1,ndim
      tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,e_)-half*((w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,b1_)**2+w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,b2_)**2+w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,b3_)**2)+(w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,m1_)**2+w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,m2_)**2+w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,m3_)**2)/(w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,rho_)+w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,rhob_)))
      call setnu(w,173,idim,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         nuR,nuL)      
      call gradient1L(tmp,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,tmp2)
      tmpL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
         =(nuL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
         +nushk(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim))&
         *tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)      
      call gradient1R(tmp,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,tmp2)
      tmpR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
         =(nuR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
         +nushk(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim))&
         *tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
      wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)&
         =wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)&
         +(tmpR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
         -tmpL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))&
         /dx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim)*qdt
enddo




      tmprhoC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
         =w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,rho_)&
         +w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,rhob_)



do k=1,ndim
        jxmin1=ixmin1+kr(k,1);jxmin2=ixmin2+kr(k,2);jxmin3=ixmin3+kr(k,3)
        jxmax1=ixmax1+kr(k,1);jxmax2=ixmax2+kr(k,2);jxmax3=ixmax3+kr(k,3); 
        hxmin1=ixmin1-kr(k,1);hxmin2=ixmin2-kr(k,2);hxmin3=ixmin3-kr(k,3)
        hxmax1=ixmax1-kr(k,1);hxmax2=ixmax2-kr(k,2);hxmax3=ixmax3-kr(k,3);
        tmprhoL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=((w(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3,rhob_))+(w(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,&
           rho_)+w(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,rhob_)))/two
        tmprhoR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=((w(jxmin1:jxmax1,&
           jxmin2:jxmax2,jxmin3:jxmax3,rho_)+w(jxmin1:jxmax1,jxmin2:jxmax2,&
           jxmin3:jxmax3,rhob_))+(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_)))/two

   do l=1,ndim
        call setnu(w,l+m0_,k,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
           nuR,nuL)      
        tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
           =w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,m0_&
           +l)/(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,rho_)&
           +w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,rhob_))


      do ii1=0,1
        if (ii1 .eq. 0) then
                           i=k
                           ii0=l
                        else
                           i=l
                           ii0=k
        endif



        if (i .eq. k) then 
           tmpVL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)&
              =(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m0_+ii0)&
              +w(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,m0_+ii0))/two
           tmpVR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)&
              =(w(jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3,m0_+ii0)&
              +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m0_+ii0))/two

           call gradient1L(tmp,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,k,&
              tmp2)
           tmpL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              =(nuL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              +nushk(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,k))&
              *tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
           call gradient1R(tmp,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,k,&
              tmp2)
           tmpR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              =(nuR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              +nushk(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,k))&
              *tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) 

           tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              =(tmprhoR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              *tmpR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              -tmprhoL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              *tmpL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))&
              /dx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,k)/two
 
           wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,m0_&
              +ii0)=wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,m0_&
              +ii0)+tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)*qdt

           tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              =(tmpVR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              *tmpR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              -tmpVL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              *tmpL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))&
              /dx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,k)/two

           wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)&
              =wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)&
              +tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)*qdt
        endif




        if (i .ne. k) then
           call gradient1(tmp,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,k,&
              tmp2)
           tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              =tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              *(nuL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              +nuR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              +two*nushk(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,k))&
              /two/two

           tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              =tmprhoC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              *tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
           call gradient1(tmp,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,i,&
              tmpC)

           wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,m0_&
              +ii0)=wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,m0_&
              +ii0)+tmpC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)*qdt

           tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              =w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,m0_&
              +ii0)*tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
           call gradient1(tmp,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,i,&
              tmpC)

           wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)&
              =wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)&
              +tmpC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)*qdt
        endif

     enddo
    enddo
   enddo





do k=1,ndim
 do l=1,ndim

   if (k .ne. l) then

    call setnu(w,b0_+l,k,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,nuR,&
       nuL)

    do ii1=0,1

      if (ii1 .eq. 0) then
              ii0=k
              m=l
              sB=-1.d0
              j=k
      endif

      if (ii1 .eq. 1) then 
              ii0=l    !ii0 is index B
              m=k      !first derivative
              sB=1.d0  !sign B
              j=l      !first B in energy
      endif



!print*,'k,l,m,j,ii0,ii1=',k,l,m,j,ii0,ii1



      if (m .eq. k) then

           jxmin1=ixmin1+kr(m,1);jxmin2=ixmin2+kr(m,2);jxmin3=ixmin3+kr(m,3)
           jxmax1=ixmax1+kr(m,1);jxmax2=ixmax2+kr(m,2);jxmax3=ixmax3+kr(m,3); 
           hxmin1=ixmin1-kr(m,1);hxmin2=ixmin2-kr(m,2);hxmin3=ixmin3-kr(m,3)
           hxmax1=ixmax1-kr(m,1);hxmax2=ixmax2-kr(m,2);hxmax3=ixmax3-kr(m,3);
           tmpBL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)&
              =(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b0_+j)&
              +w(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,b0_+j))/two
           tmpBR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)&
              =(w(jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3,b0_+j)&
              +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b0_+j))/two

           tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              =w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,b0_+l)

           call gradient1L(tmp,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,k,&
              tmp2)
           tmpL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              =(nuL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))&
              *tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
           call gradient1R(tmp,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,k,&
              tmp2)
           tmpR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              =(nuR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))&
              *tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) 

           wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,b0_&
              +ii0)=wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,b0_&
              +ii0)+sB*(tmpR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              -tmpL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))&
              /dx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,k)*qdt

           wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)&
              =wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)&
              +sB*(tmpR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              *tmpBR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              -tmpL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              *tmpBL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))&
              /dx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,k)*qdt


      endif



      if (m .ne. k) then

           tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              =w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,b0_+l)

           call gradient1(tmp,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,k,&
              tmp2)

           tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              =tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              *(nuL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              +nuR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))/two

           call gradient1(tmp2,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,m,&
              tmpC)

           wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,b0_&
              +ii0)=wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,b0_&
              +ii0)+sB*tmpC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              *qdt

           tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              =tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
              *w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,b0_+j)

           call gradient1(tmp2,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,m,&
              tmpC)

           wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)&
              =wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)&
              +sB*tmpC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)*qdt

      endif


      enddo
   endif
 enddo
enddo




return
end

!=============================================================================
subroutine setnu(w,iw,idim,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,nuR,nuL)

! Set the viscosity coefficient nu within ixO based on w(ixI). 

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)
double precision:: d1R(ixGlo1:ixGhi1+1,ixGlo2:ixGhi2+1,ixGlo3:ixGhi3&
   +1),d1L(ixGlo1:ixGhi1+1,ixGlo2:ixGhi2+1,ixGlo3:ixGhi3+1)
double precision:: d3R(ixGlo1:ixGhi1+1,ixGlo2:ixGhi2+1,ixGlo3:ixGhi3&
   +1),d3L(ixGlo1:ixGhi1+1,ixGlo2:ixGhi2+1,ixGlo3:ixGhi3+1)
double precision:: md3R(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   md3L(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
double precision:: md1R(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   md1L(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
double precision:: nuR(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   nuL(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)

double precision:: c_tot, c_hyp,cmax(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3), tmp_nu(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim, iw
integer:: kxmin1,kxmin2,kxmin3,kxmax1,kxmax2,kxmax3,jxmin1,jxmin2,jxmin3,&
   jxmax1,jxmax2,jxmax3,hxmin1,hxmin2,hxmin3,hxmax1,hxmax2,hxmax3,gxmin1,&
   gxmin2,gxmin3,gxmax1,gxmax2,gxmax3,ixFFmin1,ixFFmin2,ixFFmin3,ixFFmax1,&
   ixFFmax2,ixFFmax3,jxFFmin1,jxFFmin2,jxFFmin3,jxFFmax1,jxFFmax2,jxFFmax3,&
   hxFFmin1,hxFFmin2,hxFFmin3,hxFFmax1,hxFFmax2,hxFFmax3
integer:: ix_1,ix_2,ix_3

integer:: ixFlo1,ixFlo2,ixFlo3,ixFhi1,ixFhi2,ixFhi3,ixFmin1,ixFmin2,ixFmin3,&
   ixFmax1,ixFmax2,ixFmax3,ixYlo1,ixYlo2,ixYlo3,ixYhi1,ixYhi2,ixYhi3

logical:: new_cmax

double precision:: tmp_nuI(ixGlo1:ixGhi1+2,ixGlo2:ixGhi2+2,ixGlo3:ixGhi3+2)

integer:: k,iwc

integer:: ix,ixe

      

integer :: nmpirequest, mpirequests(2)
integer :: mpistatus(MPI_STATUS_SIZE,2)
common /mpirecv/ nmpirequest,mpirequests,mpistatus


integer:: hpe,jpe

double precision:: tgtbufferR1(1,ixGlo2:ixGhi2+2,ixGlo3:ixGhi3&
   +2),tgtbufferR2(ixGlo1:ixGhi1+2,1,ixGlo3:ixGhi3+2),tgtbufferR3&
   (ixGlo1:ixGhi1+2,ixGlo2:ixGhi2+2,1)
double precision:: tgtbufferL1(1,ixGlo2:ixGhi2+2,ixGlo3:ixGhi3&
   +2),tgtbufferL2(ixGlo1:ixGhi1+2,1,ixGlo3:ixGhi3+2),tgtbufferL3&
   (ixGlo1:ixGhi1+2,ixGlo2:ixGhi2+2,1)
double precision:: srcbufferR1(1,ixGlo2:ixGhi2+2,ixGlo3:ixGhi3&
   +2),srcbufferR2(ixGlo1:ixGhi1+2,1,ixGlo3:ixGhi3+2),srcbufferR3&
   (ixGlo1:ixGhi1+2,ixGlo2:ixGhi2+2,1)
double precision:: srcbufferL1(1,ixGlo2:ixGhi2+2,ixGlo3:ixGhi3&
   +2),srcbufferL2(ixGlo1:ixGhi1+2,1,ixGlo3:ixGhi3+2),srcbufferL3&
   (ixGlo1:ixGhi1+2,ixGlo2:ixGhi2+2,1)

integer:: n



!----------------------------------------------------------------------------

new_cmax=.true.

call getcmax(new_cmax,w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,cmax)
c_tot=maxval(cmax(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))

      call mpiallreduce(c_tot,MPI_MAX)



c_hyp=0.4d0 ! 1.4d0 ! 0.6

if (iw.eq.b1_.or.iw.eq.b2_.or.iw.eq.b3_) c_hyp=0.04d0 ! 2d0

if (iw .eq. rho_) c_hyp=0.04d0 !5d0

if (iw .eq. 173) c_hyp=0.04d0 !2d0



!if (iw.eq.b^D_|.or.) c_hyp=0.02d0

!if (iw .eq. rho_) c_hyp=0.02d0

!if (iw .eq. 173) c_hyp=0.02d0

        
if (iw .ne. 173) then     
        tmp_nu(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)=w(ixGlo1:ixGhi1,&
           ixGlo2:ixGhi2,ixGlo3:ixGhi3,iw)
        if (iw.eq.m1_.or.iw.eq.m2_.or.iw.eq.m3_) tmp_nu(ixGlo1:ixGhi1,&
           ixGlo2:ixGhi2,ixGlo3:ixGhi3)=w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
           ixGlo3:ixGhi3,iw)/(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
           rho_)+w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,rhob_))
endif

if (iw .eq. 173) tmp_nu(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)&
   =w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,e_)-half*((w(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,ixGlo3:ixGhi3,b1_)**2+w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3,b2_)**2+w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,b3_)**2)&
   +(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,m1_)**2+w(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,ixGlo3:ixGhi3,m2_)**2+w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3,m3_)**2)/(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,rho_)&
   +w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,rhob_)))


ixYlo1=ixmin1-2;ixYlo2=ixmin2-2;ixYlo3=ixmin3-2;ixYhi1=ixmax1+2
ixYhi2=ixmax2+2;ixYhi3=ixmax3+2;

ixFlo1=ixYlo1+1;ixFlo2=ixYlo2+1;ixFlo3=ixYlo3+1;ixFhi1=ixYhi1+1
ixFhi2=ixYhi2+1;ixFhi3=ixYhi3+1;

tmp_nuI(ixFlo1:ixFhi1,ixFlo2:ixFhi2,ixFlo3:ixFhi3)=tmp_nu(ixYlo1:ixYhi1,&
   ixYlo2:ixYhi2,ixYlo3:ixYhi3)


     

call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)

n = (ixFhi1-ixFlo1+1)*(ixFhi2-ixFlo2+1)*(ixFhi3-ixFlo3+1)   

select case(idim)
  case(1)

 n=n/(ixFhi1-ixFlo1+1)

  case(2)

 n=n/(ixFhi2-ixFlo2+1)

  case(3)

 n=n/(ixFhi3-ixFlo3+1)


end select



select case(idim)
   case(1)

 if(npe1>1)then

   nmpirequest =0
   mpirequests(1:2) = MPI_REQUEST_NULL


!source
   srcbufferL1(1,ixFlo2:ixFhi2,ixFlo3:ixFhi3)=tmp_nuI(ixFlo1&
      +4,ixFlo2:ixFhi2,ixFlo3:ixFhi3) !left, lower

   srcbufferR1(1,ixFlo2:ixFhi2,ixFlo3:ixFhi3)=tmp_nuI(ixFhi1&
      -4,ixFlo2:ixFhi2,ixFlo3:ixFhi3) !right, upper

   call mpineighbors(1,hpe,jpe)


   if (mpiupperB(1)) nmpirequest=nmpirequest+1
   if (mpiupperB(1)) call MPI_IRECV(tgtbufferR1(1),n,MPI_DOUBLE_PRECISION,&
       jpe,10*jpe+0,MPI_COMM_WORLD, mpirequests(nmpirequest),ierrmpi)

   if (mpilowerB(1)) nmpirequest=nmpirequest+1
   if (mpilowerB(1)) call MPI_IRECV(tgtbufferL1(1),n,MPI_DOUBLE_PRECISION,&
       hpe,10*hpe+1,MPI_COMM_WORLD, mpirequests(nmpirequest),ierrmpi)

   call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)

   if (mpiupperB(1)) call MPI_RSEND(srcbufferR1(1),n,MPI_DOUBLE_PRECISION,&
       jpe,10*ipe+1,MPI_COMM_WORLD,ierrmpi)

   if (mpilowerB(1)) call MPI_RSEND(srcbufferL1(1),n,MPI_DOUBLE_PRECISION,&
       hpe,10*ipe+0,MPI_COMM_WORLD,ierrmpi)

   call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)

!target
   tmp_nuI(ixFhi1+1,ixFlo2:ixFhi2,ixFlo3:ixFhi3)=tgtbufferR1(1,ixFlo2:ixFhi2,&
      ixFlo3:ixFhi3) !right, upper R

   tmp_nuI(ixFlo1-1,ixFlo2:ixFhi2,ixFlo3:ixFhi3)=tgtbufferL1(1,ixFlo2:ixFhi2,&
      ixFlo3:ixFhi3) !left, lower  L


  endif
    case(2)

 if(npe2>1)then

   nmpirequest =0
   mpirequests(1:2) = MPI_REQUEST_NULL


!source
   srcbufferL2(ixFlo1:ixFhi1,1,ixFlo3:ixFhi3)=tmp_nuI(ixFlo1:ixFhi1,ixFlo2&
      +4,ixFlo3:ixFhi3) !left, lower

   srcbufferR2(ixFlo1:ixFhi1,1,ixFlo3:ixFhi3)=tmp_nuI(ixFlo1:ixFhi1,ixFhi2&
      -4,ixFlo3:ixFhi3) !right, upper

   call mpineighbors(2,hpe,jpe)


   if (mpiupperB(2)) nmpirequest=nmpirequest+1
   if (mpiupperB(2)) call MPI_IRECV(tgtbufferR2(1),n,MPI_DOUBLE_PRECISION,&
       jpe,10*jpe+0,MPI_COMM_WORLD, mpirequests(nmpirequest),ierrmpi)

   if (mpilowerB(2)) nmpirequest=nmpirequest+1
   if (mpilowerB(2)) call MPI_IRECV(tgtbufferL2(1),n,MPI_DOUBLE_PRECISION,&
       hpe,10*hpe+1,MPI_COMM_WORLD, mpirequests(nmpirequest),ierrmpi)

   call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)

   if (mpiupperB(2)) call MPI_RSEND(srcbufferR2(1),n,MPI_DOUBLE_PRECISION,&
       jpe,10*ipe+1,MPI_COMM_WORLD,ierrmpi)

   if (mpilowerB(2)) call MPI_RSEND(srcbufferL2(1),n,MPI_DOUBLE_PRECISION,&
       hpe,10*ipe+0,MPI_COMM_WORLD,ierrmpi)

   call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)

!target
   tmp_nuI(ixFlo1:ixFhi1,ixFhi2+1,ixFlo3:ixFhi3)=tgtbufferR2(ixFlo1:ixFhi1,1,&
      ixFlo3:ixFhi3) !right, upper R

   tmp_nuI(ixFlo1:ixFhi1,ixFlo2-1,ixFlo3:ixFhi3)=tgtbufferL2(ixFlo1:ixFhi1,1,&
      ixFlo3:ixFhi3) !left, lower  L


  endif
    case(3)

 if(npe3>1)then

   nmpirequest =0
   mpirequests(1:2) = MPI_REQUEST_NULL


!source
   srcbufferL3(ixFlo1:ixFhi1,ixFlo2:ixFhi2,1)=tmp_nuI(ixFlo1:ixFhi1,&
      ixFlo2:ixFhi2,ixFlo3+4) !left, lower

   srcbufferR3(ixFlo1:ixFhi1,ixFlo2:ixFhi2,1)=tmp_nuI(ixFlo1:ixFhi1,&
      ixFlo2:ixFhi2,ixFhi3-4) !right, upper

   call mpineighbors(3,hpe,jpe)


   if (mpiupperB(3)) nmpirequest=nmpirequest+1
   if (mpiupperB(3)) call MPI_IRECV(tgtbufferR3(1),n,MPI_DOUBLE_PRECISION,&
       jpe,10*jpe+0,MPI_COMM_WORLD, mpirequests(nmpirequest),ierrmpi)

   if (mpilowerB(3)) nmpirequest=nmpirequest+1
   if (mpilowerB(3)) call MPI_IRECV(tgtbufferL3(1),n,MPI_DOUBLE_PRECISION,&
       hpe,10*hpe+1,MPI_COMM_WORLD, mpirequests(nmpirequest),ierrmpi)

   call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)

   if (mpiupperB(3)) call MPI_RSEND(srcbufferR3(1),n,MPI_DOUBLE_PRECISION,&
       jpe,10*ipe+1,MPI_COMM_WORLD,ierrmpi)

   if (mpilowerB(3)) call MPI_RSEND(srcbufferL3(1),n,MPI_DOUBLE_PRECISION,&
       hpe,10*ipe+0,MPI_COMM_WORLD,ierrmpi)

   call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)

!target
   tmp_nuI(ixFlo1:ixFhi1,ixFlo2:ixFhi2,ixFhi3+1)=tgtbufferR3(ixFlo1:ixFhi1,&
      ixFlo2:ixFhi2,1) !right, upper R

   tmp_nuI(ixFlo1:ixFhi1,ixFlo2:ixFhi2,ixFlo3-1)=tgtbufferL3(ixFlo1:ixFhi1,&
      ixFlo2:ixFhi2,1) !left, lower  L


  endif

end select

call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)



if (iw .eq. 173) then 
   iwc=e_ 
else 
   iwc=iw
endif

do k=0,1  !left-right bc

if (typeB(iwc,2*idim-1+k) .ne. 'mpi') then
      if (upperB(2*idim-1+k)) then

          select case(idim)
               case(1)
                     tmp_nuI(ixFhi1+1,ixFlo2:ixFhi2,ixFlo3:ixFhi3)&
                        =tmp_nuI(ixFhi1-5,ixFlo2:ixFhi2,ixFlo3:ixFhi3)
                case(2)
                     tmp_nuI(ixFlo1:ixFhi1,ixFhi2+1,ixFlo3:ixFhi3)&
                        =tmp_nuI(ixFlo1:ixFhi1,ixFhi2-5,ixFlo3:ixFhi3)
                case(3)
                     tmp_nuI(ixFlo1:ixFhi1,ixFlo2:ixFhi2,ixFhi3&
                        +1)=tmp_nuI(ixFlo1:ixFhi1,ixFlo2:ixFhi2,ixFhi3-5)
            
          end select

      else

          select case(idim)
               case(1)
                     tmp_nuI(ixFlo1-1,ixFlo2:ixFhi2,ixFlo3:ixFhi3)&
                        =tmp_nuI(ixFlo1+5,ixFlo2:ixFhi2,ixFlo3:ixFhi3)
                case(2)
                     tmp_nuI(ixFlo1:ixFhi1,ixFlo2-1,ixFlo3:ixFhi3)&
                        =tmp_nuI(ixFlo1:ixFhi1,ixFlo2+5,ixFlo3:ixFhi3)
                case(3)
                     tmp_nuI(ixFlo1:ixFhi1,ixFlo2:ixFhi2,ixFlo3&
                        -1)=tmp_nuI(ixFlo1:ixFhi1,ixFlo2:ixFhi2,ixFlo3+5)
            
          end select

      endif
endif

enddo 

        ixFmin1=ixFlo1+1;ixFmin2=ixFlo2+1;ixFmin3=ixFlo3+1;ixFmax1=ixFhi1-1
        ixFmax2=ixFhi2-1;ixFmax3=ixFhi3-1; 

        kxmin1=ixFmin1+2*kr(idim,1);kxmin2=ixFmin2+2*kr(idim,2)
        kxmin3=ixFmin3+2*kr(idim,3);kxmax1=ixFmax1+2*kr(idim,1)
        kxmax2=ixFmax2+2*kr(idim,2);kxmax3=ixFmax3+2*kr(idim,3); !5:66
        jxmin1=ixFmin1+kr(idim,1);jxmin2=ixFmin2+kr(idim,2)
        jxmin3=ixFmin3+kr(idim,3);jxmax1=ixFmax1+kr(idim,1)
        jxmax2=ixFmax2+kr(idim,2);jxmax3=ixFmax3+kr(idim,3); !4:65
        hxmin1=ixFmin1-kr(idim,1);hxmin2=ixFmin2-kr(idim,2)
        hxmin3=ixFmin3-kr(idim,3);hxmax1=ixFmax1-kr(idim,1)
        hxmax2=ixFmax2-kr(idim,2);hxmax3=ixFmax3-kr(idim,3); !2:63
        gxmin1=ixFmin1-2*kr(idim,1);gxmin2=ixFmin2-2*kr(idim,2)
        gxmin3=ixFmin3-2*kr(idim,3);gxmax1=ixFmax1-2*kr(idim,1)
        gxmax2=ixFmax2-2*kr(idim,2);gxmax3=ixFmax3-2*kr(idim,3); !1:62

        ixFFmin1=ixFlo1;ixFFmin2=ixFlo2;ixFFmin3=ixFlo3;ixFFmax1=ixFhi1
        ixFFmax2=ixFhi2;ixFFmax3=ixFhi3; !2:65
        jxFFmin1=ixFlo1+kr(idim,1);jxFFmin2=ixFlo2+kr(idim,2)
        jxFFmin3=ixFlo3+kr(idim,3);jxFFmax1=ixFhi1+kr(idim,1)
        jxFFmax2=ixFhi2+kr(idim,2);jxFFmax3=ixFhi3+kr(idim,3); !3:66
        hxFFmin1=ixFlo1-kr(idim,1);hxFFmin2=ixFlo2-kr(idim,2)
        hxFFmin3=ixFlo3-kr(idim,3);hxFFmax1=ixFhi1-kr(idim,1)
        hxFFmax2=ixFhi2-kr(idim,2);hxFFmax3=ixFhi3-kr(idim,3); !1:64
        
        d3R(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ixFmin3:ixFmax3)=abs(3.d0&
           *(tmp_nuI(jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3)&
           -tmp_nuI(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ixFmin3:ixFmax3))&
           -(tmp_nuI(kxmin1:kxmax1,kxmin2:kxmax2,kxmin3:kxmax3)&
           -tmp_nuI(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3))) !3:64
        d1R(ixFFmin1:ixFFmax1,ixFFmin2:ixFFmax2,ixFFmin3:ixFFmax3)&
           =abs(tmp_nuI(jxFFmin1:jxFFmax1,jxFFmin2:jxFFmax2,&
           jxFFmin3:jxFFmax3)-tmp_nuI(ixFFmin1:ixFFmax1,ixFFmin2:ixFFmax2,&
           ixFFmin3:ixFFmax3)) !2:65

        do ix_1=ixmin1,ixmax1
        do ix_2=ixmin2,ixmax2
        do ix_3=ixmin3,ixmax3    !3:62  +1=4:63

          md3R(ix_1,ix_2,ix_3)=maxval(d3R(ix_1+1-kr(idim,1):ix_1+1&
             +kr(idim,1),ix_2+1-kr(idim,2):ix_2+1+kr(idim,2),ix_3&
             +1-kr(idim,3):ix_3+1+kr(idim,3)))
          md1R(ix_1,ix_2,ix_3)=maxval(d1R(ix_1+1-2*kr(idim,1):ix_1+1&
             +2*kr(idim,1),ix_2+1-2*kr(idim,2):ix_2+1+2*kr(idim,2),ix_3&
             +1-2*kr(idim,3):ix_3+1+2*kr(idim,3)))

        enddo
        enddo
        enddo

        WHERE (md1R(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3).gt.0.d0)
          nuR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=c_tot*c_hyp&
             *md3R(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)&
             /md1R(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)&
             *dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idim)
        ELSEWHERE 
          nuR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
        END WHERE
        
        maxviscoef=max(maxval(nuR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)),&
            maxviscoef)


!************

        d3L(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ixFmin3:ixFmax3)=abs(3.d0&
           *(tmp_nuI(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ixFmin3:ixFmax3)&
           -tmp_nuI(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3))&
           -(tmp_nuI(jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3)&
           -tmp_nuI(gxmin1:gxmax1,gxmin2:gxmax2,gxmin3:gxmax3)))
        d1L(ixFFmin1:ixFFmax1,ixFFmin2:ixFFmax2,ixFFmin3:ixFFmax3)&
           =abs(tmp_nuI(ixFFmin1:ixFFmax1,ixFFmin2:ixFFmax2,&
           ixFFmin3:ixFFmax3)-tmp_nuI(hxFFmin1:hxFFmax1,hxFFmin2:hxFFmax2,&
           hxFFmin3:hxFFmax3))    

        do ix_1=ixmin1,ixmax1
        do ix_2=ixmin2,ixmax2
        do ix_3=ixmin3,ixmax3

          md3L(ix_1,ix_2,ix_3)=maxval(d3L(ix_1+1-kr(idim,1):ix_1+1&
             +kr(idim,1),ix_2+1-kr(idim,2):ix_2+1+kr(idim,2),ix_3&
             +1-kr(idim,3):ix_3+1+kr(idim,3)))
          md1L(ix_1,ix_2,ix_3)=maxval(d1L(ix_1+1-2*kr(idim,1):ix_1+1&
             +2*kr(idim,1),ix_2+1-2*kr(idim,2):ix_2+1+2*kr(idim,2),ix_3&
             +1-2*kr(idim,3):ix_3+1+2*kr(idim,3)))

        enddo
        enddo
        enddo

        WHERE (md1L(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3).gt.0.d0)
          nuL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=c_tot*c_hyp&
             *md3L(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)&
             /md1L(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)&
             *dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idim)
        ELSEWHERE 
          nuL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0  
        END WHERE

        maxviscoef=max(maxval(nuL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)),&
            maxviscoef)

              call mpiallreduce(maxviscoef,MPI_MAX)

return
end


!=============================================================================
!=============================================================================
SUBROUTINE setnushk(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,nushk)

  INCLUDE 'vacdef.f'

  DOUBLE PRECISION:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw),&
     tmp2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),nushk(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2,ixGlo3:ixGhi3,ndim)
  DOUBLE PRECISION:: c_shk
  
  DOUBLE PRECISION:: tmp3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
  
  INTEGER:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim, iw,i
  
  INTEGER:: ix_1,ix_2
  
  DO idim=1,ndim
     nushk(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idim)=0.d0
  END DO
  
  
!!$go to 100
!!$c_shk=0.5d0
!!$
!!$tmp3(ix^S)=0.d0
!!$
!!$!**************************BEGIN shock viscosity*******************************
!!$      do idim=1,ndim
!!$         tmp(ix^S)=w(ix^S,m0_+idim)/(w(ix^S,rho_)+w(ix^S,rhob_))
!!$         call gradient1(tmp,ix^L,idim,tmp2)
!!$         tmp3(ix^S)=tmp3(ix^S)+tmp2(ix^S)
!!$       enddo
!!$      do idim=1,ndim
!!$        nushk(ix^S,idim)=tmp3(ix^S)*(dx(ix^S,idim)**2.d0)*c_shk
!!$	WHERE (tmp3(ix^S) .ge. 0.d0)
!!$!	  nushk(ix^S,idim)=0.d0
!!$	END WHERE
!!$	nushk(ix^S,idim)=abs(nushk(ix^S,idim))
!!$      enddo
!!$!****************************END shock viscosity*******************************

!!$100 continue


  RETURN
END SUBROUTINE setnushk



!=============================================================================
subroutine getdt_visc(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3)

! Check diffusion time limit for dt < dtdiffpar * dx**2 / (nu/rho)

! Based on Hirsch volume 2, p.631, eq.23.2.17

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw),dtdiff_visc
integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim, ix_1,ix_2

integer:: aa

! For spatially varying nu you need a common nu array
 double precision::tmpdt(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
     nuL(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),nuR(ixGlo1:ixGhi1,&
    ixGlo2:ixGhi2,ixGlo3:ixGhi3), nushk(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
    ixGlo3:ixGhi3,ndim)
 common/visc/nuL
 common/visc/nuR
!-----------------------------------------------------------------------------

 call setnushk(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,nushk)

 dtdiffpar=0.25d0

 do idim=1,ndim
   tmpdt(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(maxviscoef&
      +nushk(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idim)) !/(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_))   ! 1/dt
   dtdiff_visc=dtdiffpar/maxval(tmpdt(ixmin1:ixmax1,ixmin2:ixmax2,&
      ixmin3:ixmax3)/(dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idim)**2))
         call mpiallreduce(dtdiff_visc,MPI_MIN)
   dt=min(dt,dtdiff_visc)
 end do
 
 maxviscoef=0.d0

return
end


!***** 2-point central finite difference gradient******

subroutine gradient1(q,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,gradq)
include 'vacdef.f'
integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim
double precision:: q(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   gradq(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
integer:: hxmin1,hxmin2,hxmin3,hxmax1,hxmax2,hxmax3,kxmin1,kxmin2,kxmin3,&
   kxmax1,kxmax2,kxmax3
integer:: minx11,minx12,minx13,maxx11,maxx12,maxx13,k
!-----------------------------------------------------------------------------

hxmin1=ixmin1-kr(idim,1);hxmin2=ixmin2-kr(idim,2);hxmin3=ixmin3-kr(idim,3)
hxmax1=ixmax1-kr(idim,1);hxmax2=ixmax2-kr(idim,2);hxmax3=ixmax3-kr(idim,3);
kxmin1=ixmin1+kr(idim,1);kxmin2=ixmin2+kr(idim,2);kxmin3=ixmin3+kr(idim,3)
kxmax1=ixmax1+kr(idim,1);kxmax2=ixmax2+kr(idim,2);kxmax3=ixmax3+kr(idim,3);
gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(q(kxmin1:kxmax1,&
   kxmin2:kxmax2,kxmin3:kxmax3)-q(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3))&
   /dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idim)/two

 minx11=ixmin1+kr(idim,1);minx12=ixmin2+kr(idim,2);minx13=ixmin3+kr(idim,3);
 maxx11=ixmax1-kr(idim,1);maxx12=ixmax2-kr(idim,2);maxx13=ixmax3-kr(idim,3);
 
 do k=0,1  !left-right bc
 
 if (typeB(1,2*idim-1+k) .ne. 'mpi') then
 if (upperB(2*idim-1+k)) then
 select case(idim)
    case(1)
 gradq(ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
 gradq(maxx11,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
    case(2)
 gradq(ixmin1:ixmax1,ixmax2,ixmin3:ixmax3)=0.d0
 gradq(ixmin1:ixmax1,maxx12,ixmin3:ixmax3)=0.d0
    case(3)
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmax3)=0.d0
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,maxx13)=0.d0

 end select
 else
 select case(idim)
    case(1)
 gradq(ixmin1,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
 gradq(minx11,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
    case(2)
 gradq(ixmin1:ixmax1,ixmin2,ixmin3:ixmax3)=0.d0
 gradq(ixmin1:ixmax1,minx12,ixmin3:ixmax3)=0.d0
    case(3)
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3)=0.d0
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,minx13)=0.d0

 end select
 endif
 endif
 enddo


return
end

!=============================================================================


!*****left upwind forward 2-point non-central finite difference gradient******

subroutine gradient1L(q,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,gradq)
include 'vacdef.f'
integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim
double precision:: q(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   gradq(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
integer:: hxmin1,hxmin2,hxmin3,hxmax1,hxmax2,hxmax3
integer:: minx11,minx12,minx13,maxx11,maxx12,maxx13,k
!-----------------------------------------------------------------------------

hxmin1=ixmin1-kr(idim,1);hxmin2=ixmin2-kr(idim,2);hxmin3=ixmin3-kr(idim,3)
hxmax1=ixmax1-kr(idim,1);hxmax2=ixmax2-kr(idim,2);hxmax3=ixmax3-kr(idim,3);
gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(q(ixmin1:ixmax1,&
   ixmin2:ixmax2,ixmin3:ixmax3)-q(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3))&
   /dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idim)

 minx11=ixmin1+kr(idim,1);minx12=ixmin2+kr(idim,2);minx13=ixmin3+kr(idim,3);
 maxx11=ixmax1-kr(idim,1);maxx12=ixmax2-kr(idim,2);maxx13=ixmax3-kr(idim,3);
 
 do k=0,1  !left-right bc
 
 if (typeB(1,2*idim-1+k) .ne. 'mpi') then
 if (upperB(2*idim-1+k)) then
 select case(idim)
    case(1)
 gradq(ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
 gradq(maxx11,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
    case(2)
 gradq(ixmin1:ixmax1,ixmax2,ixmin3:ixmax3)=0.d0
 gradq(ixmin1:ixmax1,maxx12,ixmin3:ixmax3)=0.d0
    case(3)
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmax3)=0.d0
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,maxx13)=0.d0

 end select
 else
 select case(idim)
    case(1)
 gradq(ixmin1,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
 gradq(minx11,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
    case(2)
 gradq(ixmin1:ixmax1,ixmin2,ixmin3:ixmax3)=0.d0
 gradq(ixmin1:ixmax1,minx12,ixmin3:ixmax3)=0.d0
    case(3)
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3)=0.d0
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,minx13)=0.d0

 end select
 endif
 endif
 enddo


return
end

!=============================================================================

!*****right upwind forward 2-point non-central finite difference gradient*****

subroutine gradient1R(q,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,gradq)
include 'vacdef.f'
integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim
double precision:: q(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   gradq(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
integer:: hxmin1,hxmin2,hxmin3,hxmax1,hxmax2,hxmax3
integer:: minx11,minx12,minx13,maxx11,maxx12,maxx13,k
!-----------------------------------------------------------------------------

hxmin1=ixmin1+kr(idim,1);hxmin2=ixmin2+kr(idim,2);hxmin3=ixmin3+kr(idim,3)
hxmax1=ixmax1+kr(idim,1);hxmax2=ixmax2+kr(idim,2);hxmax3=ixmax3+kr(idim,3);
gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(q(hxmin1:hxmax1,&
   hxmin2:hxmax2,hxmin3:hxmax3)-q(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))&
   /dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idim)

 minx11=ixmin1+kr(idim,1);minx12=ixmin2+kr(idim,2);minx13=ixmin3+kr(idim,3);
 maxx11=ixmax1-kr(idim,1);maxx12=ixmax2-kr(idim,2);maxx13=ixmax3-kr(idim,3);
 
 do k=0,1  !left-right bc
 
 if (typeB(1,2*idim-1+k) .ne. 'mpi') then
 if (upperB(2*idim-1+k)) then
 select case(idim)
    case(1)
 gradq(ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
 gradq(maxx11,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
    case(2)
 gradq(ixmin1:ixmax1,ixmax2,ixmin3:ixmax3)=0.d0
 gradq(ixmin1:ixmax1,maxx12,ixmin3:ixmax3)=0.d0
    case(3)
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmax3)=0.d0
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,maxx13)=0.d0

 end select
 else
 select case(idim)
    case(1)
 gradq(ixmin1,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
 gradq(minx11,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
    case(2)
 gradq(ixmin1:ixmax1,ixmin2,ixmin3:ixmax3)=0.d0
 gradq(ixmin1:ixmax1,minx12,ixmin3:ixmax3)=0.d0
    case(3)
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3)=0.d0
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,minx13)=0.d0

 end select
 endif
 endif
 enddo


return
end

  !=============================================================================
  SUBROUTINE specialini(ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w)

    INCLUDE 'vacdef.f'

    INTEGER:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    DOUBLE PRECISION:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)

    RETURN
  END SUBROUTINE specialini


  !=============================================================================
  SUBROUTINE specialsource(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iws,qtC,wCT,qt,w)

    INCLUDE 'vacdef.f'

    INTEGER:: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
       ixOmin3,ixOmax1,ixOmax2,ixOmax3,iws(niw_)
    DOUBLE PRECISION:: qdt,qtC,qt,wCT(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,nw),w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw) !Intent(IN)

    INTEGER:: iw
    INTEGER:: ix_1,ix_2,ix_3

    DOUBLE PRECISION:: s_period, xc1, xc2, xc3, xxmax, yymax, zzmax
    DOUBLE PRECISION:: xc1Mm, xc2Mm, xc3Mm
    DOUBLE PRECISION:: xx, yy, zz
    DOUBLE PRECISION:: vvx(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
        vvy(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3), vvz(ixGlo1:ixGhi1,&
       ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    DOUBLE PRECISION:: AA, B, ux, uy, phi, theta
    DOUBLE PRECISION:: delta_z, delta_x, delta_y, exp_x, exp_y, exp_z,&
        exp_xyz, tdep
    !-----------------------------------------------------------------------------

    eqpar(eta_)=0.d0
    eqpar(nu_)=1.0d0
    
    CALL addsource_grav(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iws,qtC,wCT,qt,w)
    
    IF(ABS(eqpar(nu_))>smalldouble)CALL addsource_visc(qdt,ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,iws,qtC,wCT,qt,w)

    vvx(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3) = 0.d0
    vvy(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3) = 0.d0
    vvz(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3) = 0.d0


    xc1Mm=0.1   !Mm        z axis
    xc2Mm=1.0 !0.99d0  !Mm        x axis
    xc3Mm=1.0 !0.99d0  !Mm        y axis


    xc1=xc1Mm*1.0d6  !m        z axis
    xc2=xc2Mm*1.0d6  !m        x axis
    xc3=xc3Mm*1.0d6  !m        y axis

    xxmax=2.d6
    yymax=2.d6
    zzmax=1.6d6
    
    delta_z=0.05d6
    delta_x=0.1d6
    delta_y=0.1d6

 !   AA=80000.d0 ! ~2km/s hor/vert

    !AA=0.005 !~2.4k
    AA=10 !0.3k/s tor
    B = 0.05

    s_period=120.d0

    tdep=SIN(qt*2.d0*pi/s_period)!*EXP(-(qt/(4.d0 * s_period))**2.d0) 
          IF (ipe.EQ.0) PRINT*, "tdep==",tdep
    !tdep=sin(qt*2.d0*pi/s_period) 

    !tdep=1.d0 ! continious rotation

    DO ix_1=ixImin1,ixImax1
       DO ix_2=ixImin2,ixImax2
          DO ix_3=ixImin3,ixImax3


             xx=x(ix_1,ix_2,ix_3,2)-xc2
             yy=x(ix_1,ix_2,ix_3,3)-xc3
             zz=x(ix_1,ix_2,ix_3,1)-xc1  
             
             exp_z = EXP(-zz**2.d0/(delta_z**2.d0))
             exp_x = EXP(-xx**2.d0/(delta_x**2.d0))
             exp_y = EXP(-yy**2.d0/(delta_y**2.d0))

             exp_xyz=exp_x*exp_y*exp_z
             
             theta = ATAN2(yy,xx)
             phi = ATAN2(1.d0,B)
             

             ux = COS(theta + phi)
             uy =  SIN(theta + phi)
             
             vvx(ix_1,ix_2,ix_3) = (ux / SQRT(ux**2 + uy**2)) * exp_xyz &
                * tdep * AA
             vvy(ix_1,ix_2,ix_3) = (uy / SQRT(ux**2 + uy**2)) * exp_xyz &
                * tdep * AA

          ENDDO
       ENDDO
    ENDDO


    DO ix_1=ixImin1,ixImax1
       DO ix_2=ixImin2,ixImax2
          DO ix_3=ixImin3,ixImax3
             
            !w(ix_1,ix_2,ix_3,m1_)=w(ix_1,ix_2,ix_3,m1_)+(w(ix_1,ix_2,ix_3,rho_)+w(ix_1,ix_2,ix_3,rhob_))*vvz(ix_1,ix_2,ix_3)*qdt
             
             w(ix_1,ix_2,ix_3,m2_)=w(ix_1,ix_2,ix_3,m2_)+(w(ix_1,ix_2,ix_3,&
                rho_)+w(ix_1,ix_2,ix_3,rhob_))*vvx(ix_1,ix_2,ix_3)*qdt
             
             w(ix_1,ix_2,ix_3,m3_)= w(ix_1,ix_2,ix_3,m3_)+(w(ix_1,ix_2,ix_3,&
                rho_)+w(ix_1,ix_2,ix_3,rhob_))*vvy(ix_1,ix_2,ix_3)*qdt

             
             w(ix_1,ix_2,ix_3,e_)=w(ix_1,ix_2,ix_3,e_)+(w(ix_1,ix_2,ix_3,&
                rho_)+w(ix_1,ix_2,ix_3,rhob_))*(vvx(ix_1,ix_2,ix_3)**2.d0 &
                + vvy(ix_1,ix_2,ix_3)**2.d0 + vvz(ix_1,ix_2,ix_3)**2.d0)*qdt&
                /2.d0


          ENDDO
       ENDDO
    ENDDO


          IF (ipe.EQ.0) WRITE(*,*) '***time=',qt

  END SUBROUTINE specialsource



!=============================================================================
SUBROUTINE specialbound(qt,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,iw,iB,w)
  INCLUDE 'vacdef.f'

  INTEGER:: ix_1,ix_2

  INTEGER:: iwmin,iwmax,idimmin,idimmax
  DOUBLE PRECISION:: qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)
  INTEGER:: ix,ix1,ix2,ix3,ixe,ixf,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,&
     ixpairmin1,ixpairmin2,ixpairmin3,ixpairmax1,ixpairmax2,ixpairmax3,idim,&
     iw,iB

  CALL die('not defined')
  RETURN
END SUBROUTINE specialbound

!=============================================================================
SUBROUTINE getdt_special(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3)

  ! If the Coriolis force is made very strong it may require time step limiting,
  ! but this is not implemented here.

  INCLUDE 'vacdef.f'
  DOUBLE PRECISION:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)
  INTEGER:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
  !----------------------------------------------------------------------------

  !call getdt_diff(w,ix^L)

  IF(ABS(eqpar(nu_))>smalldouble)CALL getdt_visc(w,ixmin1,ixmin2,ixmin3,&
     ixmax1,ixmax2,ixmax3)

  CALL getdt_grav(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3)

  RETURN
END SUBROUTINE getdt_special

!=============================================================================
SUBROUTINE specialeta(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idirmin)

  INCLUDE 'vacdef.f'

  DOUBLE PRECISION:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)
  INTEGER:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idirmin
!------------------------------------------------------------------------------

  STOP 'specialeta is not defined'
END SUBROUTINE specialeta

!------------------------------------------------------------------------------
SUBROUTINE savefilelog_special(qunit,w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
   ixmax3)
  
  ! This is a save log file routine to calculate and save out Vpar Vperp and Vaz
  ! It mimics savefileout_bin to mantain compatibility with usual readin routines
  
  INCLUDE 'vacdef.f'

  INTEGER:: qunit,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
  DOUBLE PRECISION :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)

  CALL die("log save is not defined")

END SUBROUTINE savefilelog_special

