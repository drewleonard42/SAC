!##############################################################################
! module vacphys.mhd0 - common subroutines for mhd and mhdiso

INCLUDE:vacphys.mhdroe.t^IFTVD
INCLUDE:vacproc.projectb.t^NOONED^IFPOISSON
INCLUDE:vacproc.constrainb.t^NOONED^IFCT
INCLUDE:vacphys.mhdres.t^IFRES
!=============================================================================

SUBROUTINE physini

  ! Tell VAC which variables are vectors, set default entropy coefficients
  USE constants
  USE common_variables

  INTEGER:: il
  !-----------------------------------------------------------------------------

  iw_vector(1)=m0_; iw_vector(2)=b0_

  ! The values of the constants are taken from Ryu & Jones ApJ 442, 228
  DO il=1,nw
     SELECT CASE(il)
     CASE(fastRW_,fastLW_,slowRW_,slowLW_)
        entropycoef(il)=0.2
     CASE(alfvRW_,alfvLW_)
        entropycoef(il)=0.4
     CASE default
        entropycoef(il)= -one
     END SELECT
  END DO

  RETURN
END SUBROUTINE physini

!=============================================================================
SUBROUTINE process(count,idim^LIM,w)

  ! Process w before it is advected in directions idim^LIM, or before save
  ! count=1 and 2 for first and second (half step) processing during advection
  ! count=ifile+2 for saving results into the file indexed by ifile

  USE constants
  USE common_variables

  INTEGER:: count,idim^LIM
  DOUBLE PRECISION:: w(ixG^T,nw)

  LOGICAL:: oktime
  DOUBLE PRECISION:: cputime,time1,timeproc
  DATA timeproc /0.D0/

  ! The processing should eliminate divergence of B.
  !-----------------------------------------------------------------------------

  oktest=INDEX(teststr,'process')>=1
  oktime=INDEX(teststr,'timeproc')>=1

  IF(oktest)WRITE(*,*)'Process it,idim^LIM,count',it,idim^LIM,count

  IF(oktime)time1=cputime()

  {^NOONED
  IF(count==0)THEN
     IF(divbconstrain)THEN
        {^IFCT CALL constrainb(w)
        IF(.FALSE.)}CALL die('CT module is OFF: setvac -on=ct; make vac')
     ENDIF
  ELSE
     ! Use the projection scheme 
     {^IFPOISSON CALL projectb(w)
     IF(.FALSE.)}CALL die('Poisson module is OFF: setvac -on=poisson;make vac')
  ENDIF
  }

  IF(oktime)THEN
     time1=cputime()-time1
     timeproc=timeproc+time1
     WRITE(*,*)'Time.Proc:',time1,timeproc
  ENDIF

  RETURN
END SUBROUTINE process

!=============================================================================
SUBROUTINE getdt(w,ix^L)

  ! If resistivity is  not zero, check diffusion time limit for dt

  USE constants
  USE common_variables

  DOUBLE PRECISION:: w(ixG^T,nw)
  INTEGER:: ix^L
  !-----------------------------------------------------------------------------

  oktest=INDEX(teststr,'getdt')>=1
  IF(oktest)WRITE(*,*)'GetDt'

  IF(eqpar(eta_)==zero)RETURN

  {^IFRES CALL getdt_res(w,ix^L)}
  {^NORES WRITE(*,*)'Error: Resistive MHD module is OFF'
  CALL die('Recompile with setvac -on=resist or set eqpar(eta_)=0')}

  RETURN
END SUBROUTINE getdt

!=============================================================================
SUBROUTINE getdivb(w,ixO^L,divb)

  ! Calculate div B within ixO

  USE constants
  USE common_variables

  INTEGER::          ixO^L,ix^L,idim
  DOUBLE PRECISION:: w(ixG^T,nw),divb(ixG^T)
  !-----------------------------------------------------------------------------

  oktest=INDEX(teststr,'getdivb')>=1
  IF(oktest)WRITE(*,*)'getdivb ixO=',ixO^L

  IF(fourthorder)THEN
     ix^L=ixO^L^LADD2;
  ELSE
     ix^L=ixO^L^LADD1;
  ENDIF
  divb(ixO^S)=zero
  DO idim=1,ndim
     tmp(ix^S)=w(ix^S,b0_+idim)+w(ix^S,bg0_+idim)

     CALL gradient4(.FALSE.,tmp,ixO^L,idim,tmp2)

     divb(ixO^S)=divb(ixO^S)+tmp2(ixO^S)
  ENDDO

  IF(oktest)THEN
     WRITE(*,*)'divb:',divb(ixtest^D)
     !   write(*,*)'bx=',w(ixtest1-1:ixtest1+1,ixtest2,b1_)
     !   write(*,*)'by=',w(ixtest1,ixtest2-1:ixtest2+1,b2_)
     !   write(*,*)'x =',x(ixtest1-1:ixtest1+1,ixtest2,1)
     !   write(*,*)'y =',x(ixtest1,ixtest2-1:ixtest2+1,2)
     !   write(*,*)'dx=',dx(ixtest1,ixtest2,1)
     !   write(*,*)'dy=',dx(ixtest1,ixtest2,2)
  ENDIF

  RETURN
END SUBROUTINE getdivb

!=============================================================================
SUBROUTINE getflux(w,ix^L,iw,idim,f,transport)

  ! Calculate non-transport flux f_idim[iw] within ix^L.
  ! Set transport=.true. if a transport flux should be added

  USE constants
  USE common_variables

  INTEGER::          ix^L,iw,idim
  DOUBLE PRECISION:: w(ixG^T,nw),f(ixG^T), fb(ixG^T)
  LOGICAL::          transport
  !-----------------------------------------------------------------------------

  oktest= INDEX(teststr,'getflux')>=1
  IF(oktest.AND.iw==iwtest)WRITE(*,*)'Getflux idim,w:',&
       idim,w(ixtest^D,iwtest)

  transport=.TRUE.


  SELECT CASE(iw)
  CASE(rho_)
     f(ix^S)=w(ix^S,rhob_)*w(ix^S,m0_+idim)/(w(ix^S,rho_)+w(ix^S,rhob_))

     {CASE(m^C_)
     IF(idim==^C)THEN

        CALL getptotal(w,ix^L,f)
        CALL getptotal_bg(w,ix^L,fb)
        fb(ix^S)=0.d0

        f(ix^S)=f(ix^S)+fb(ix^S)-(w(ix^S,b^C_)*w(ix^S,bg0_+idim)+w(ix^S,b0_+idim)*w(ix^S,bg^C_))-&
             w(ix^S,b^C_)*w(ix^S,b0_+idim) !-&
        !  w(ix^S,bg0_+idim)*w(ix^S,bg^C_)    !remove for perturbed
     ELSE
        f(ix^S)=-(w(ix^S,b^C_)*w(ix^S,bg0_+idim)+w(ix^S,b0_+idim)*w(ix^S,bg^C_))-&
             w(ix^S,b^C_)*w(ix^S,b0_+idim) !-&
        !   -w(ix^S,bg0_+idim)*w(ix^S,bg^C_)  !remove for perturbed

     ENDIF\}

  CASE(e_)

     CALL getptotal(w,ix^L,f)
     CALL getptotal_bg(w,ix^L,fb)      
     fb(ix^S)=0.d0      

     f(ix^S)=(w(ix^S,m0_+idim)*(f(ix^S)+fb(ix^S))-&
          w(ix^S,b0_+idim)*( ^C&(w(ix^S,bg^C_))*w(ix^S,m^C_)+ )-&
          w(ix^S,bg0_+idim)*( ^C&(w(ix^S,b^C_))*w(ix^S,m^C_)+ ))/(w(ix^S,rho_)+w(ix^S,rhob_))+&	      
          w(ix^S,eb_)*w(ix^S,m0_+idim)/(w(ix^S,rho_)+w(ix^S,rhob_))-&
          w(ix^S,b0_+idim)*( ^C&(w(ix^S,b^C_))*w(ix^S,m^C_)+ )/(w(ix^S,rho_)+w(ix^S,rhob_)) !&
     !  -w(ix^S,bg0_+idim)*( ^C&(w(ix^S,bg^C_))*w(ix^S,m^C_)+ )/(w(ix^S,rho_)+w(ix^S,rhob_))  ! remove for perturbed




     {
  CASE(b^C_)
     IF(idim==^C) THEN
        f(ix^S)= zero
        transport=.FALSE.
     ELSE

        f(ix^S)= -w(ix^S,m^C_)/(w(ix^S,rho_)+w(ix^S,rhob_))*(w(ix^S,b0_+idim)+w(ix^S,bg0_+idim))+ &
             w(ix^S,m0_+idim)/(w(ix^S,rho_)+w(ix^S,rhob_))*w(ix^S,bg^C_)

     ENDIF
\}

  CASE default
     print*, iw
     CALL die('Error in getflux: unknown flow variable')
  END SELECT

  IF(oktest.AND.iw==iwtest)WRITE(*,*)'transport,f:',&
       transport,f(ixtest^D)

  RETURN
END SUBROUTINE getflux

!=============================================================================
SUBROUTINE addsource(qdt,ixI^L,ixO^L,iws,qtC,w,qt,wnew)

  ! Add sources from resistivity and Powell solver

  USE constants
  USE common_variables

  INTEGER::          ixI^L,ixO^L,iws(niw_)
  DOUBLE PRECISION:: qdt,qtC,qt,w(ixG^T,nw),wnew(ixG^T,nw)
  !-----------------------------------------------------------------------------

  oktest=INDEX(teststr,'addsource')>=1
  IF(oktest)WRITE(*,*)'Addsource, compactres,divbfix:',compactres,divbfix
  IF(oktest)WRITE(*,*)'Before adding source:',wnew(ixtest^D,iwtest)

  ! Sources for resistivity in eqs. for e, B1, B2 and B3
  IF(ABS(eqpar(eta_))>smalldouble)THEN
     {^IFRES
     IF(compactres)THEN
        CALL addsource_res1(qdt,ixI^L,ixO^L,iws,qtC,w,qt,wnew)
     ELSE
        CALL addsource_res2(qdt,ixI^L,ixO^L,iws,qtC,w,qt,wnew)
     ENDIF
     IF(oktest)WRITE(*,*)'With resistive source:',wnew(ixtest^D,iwtest)
     }
     {^NORES WRITE(*,*)'Error: Resistive MHD module is OFF'
     CALL die('Recompile with setvac -on=resist or set eqpar(eta_)=0')}
  ENDIF


  ! Sources related to div B in the Powell solver
  {^NOONED IF(divbfix) CALL addsource_divb(qdt,ixI^L,ixO^L,iws,qtC,w,qt,wnew)}

  IF(oktest)WRITE(*,*)'After adding source:',wnew(ixtest^D,iwtest)

  RETURN
END SUBROUTINE addsource

!=============================================================================
SUBROUTINE addsource_divb(qdt,ixI^L,ixO^L,iws,qtC,w,qt,wnew)

  ! Add Powell's divB related sources to wnew within ixO if possible, 
  ! otherwise shrink ixO

  USE constants
  USE common_variables

  INTEGER::          ixI^L,ixO^L,iws(niw_),iiw,iw
  DOUBLE PRECISION:: qdt,qtC,qt,w(ixG^T,nw),wnew(ixG^T,nw)
  DOUBLE PRECISION:: divb(ixG^T)
  !-----------------------------------------------------------------------------

  ! Calculating div B involves first derivatives
  CALL ensurebound(1,ixI^L,ixO^L,qtC,w)

  ! We calculate now div B
  CALL getdivb(w,ixO^L,divb)
  divb(ixO^S)=qdt*divb(ixO^S)

  DO iiw=1,iws(niw_); iw=iws(iiw)
     SELECT CASE(iw)
        {CASE(m^C_)
        wnew(ixO^S,iw)=wnew(ixO^S,iw)-(w(ixO^S,b^C_)+w(ixO^S,bg^C_))*divb(ixO^S)
        }
        {CASE(b^C_)
        wnew(ixO^S,iw)=wnew(ixO^S,iw)-w(ixO^S,m^C_)/(w(ixO^S,rho_)+w(ixO^S,rhob_))*divb(ixO^S)
        }
     CASE(e_)
        wnew(ixO^S,iw)=wnew(ixO^S,iw)-&
             (^C&w(ixO^S,m^C_)*(w(ixO^S,b^C_)+w(ixO^S,bg^C_))+ )/(w(ixO^S,rho_)+w(ixO^S,rhob_))*divb(ixO^S)
     END SELECT
  END DO

  RETURN
END SUBROUTINE addsource_divb


!=============================================================================
! end module vacphys.mhd0
!##############################################################################
