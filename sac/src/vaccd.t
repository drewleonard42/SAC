!#############################################################################
! module vaccd
! Centered difference scheme
!=============================================================================

!=============================================================================
SUBROUTINE centdiff4(qdt,ixI^L,ixO^L,iws,idim^LIM,qtC,wCT,qt,w)

  ! Advance the iws flow variables from t to t+qdt within ixO^L by 
  ! fourth order centered  differencing in space the dw/dt+dF_i(w)/dx_i=S 
  ! type equation. 
  ! wCT contains the time centered variables at time qtC for flux and source.
  ! w is the old value at qt on input and the new value at qt+qdt on output.

  USE constants
  USE common_variables

  DOUBLE PRECISION:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)
  INTEGER:: ixI^L,ixO^L,iws(niw_),idim^LIM
  LOGICAL :: transport

  DOUBLE PRECISION:: v(ixG^T),f(ixG^T), fb(ixG^T)
  INTEGER:: iiw,iw,ix^L,idim,idir
  !-----------------------------------------------------------------------------


  ! Two extra layers are needed in each direction for which fluxes are added.
  ix^L=ixO^L;
  DO idim= idim^LIM
     ix^L=ix^L^LADD2*kr(idim,^D);
  ENDDO
  IF(ixI^L^LTix^L|.OR.|.OR.) CALL die( &
       'Error in CentDiff4: Non-conforming input limits')

  ! Add fluxes to w
  DO idim= idim^LIM
     ix^L=ixO^L^LADD2*kr(idim,^D);

     CALL getv(wCT,ix^L,idim,v)

     DO iiw=1,iws(niw_); iw=iws(iiw)
        !   print*,'iiw', iiw,idim,idir
        ! Get non-transported flux
        CALL getflux(wCT,ix^L,iw,idim,f,transport)

        ! Add transport flux
        IF(transport)f(ix^S)=f(ix^S)+v(ix^S)*wCT(ix^S,iw)

        ! Add divergence of flux
        CALL gradient4(.FALSE.,f,ixO^L,idim,tmp)
        w(ix^S,iw)=w(ix^S,iw)-qdt*tmp(ix^S)

        SELECT CASE(iw)

        CASE(e_)

           CALL gradient4(.FALSE.,v,ixO^L,idim,tmp)   
           CALL getptotal_bg(w,ix^L,fb)

           w(ix^S,iw)=w(ix^S,iw)-qdt*tmp(ix^S)*fb(ix^S)

           DO idir= idim^LIM 
              CALL gradient4(.FALSE.,v,ixO^L,idir,tmp)   
              w(ix^S,iw)=w(ix^S,iw)+qdt*w(ix^S,bg0_+idir)*w(ix^S,bg0_+idim)*tmp(ix^S)
           ENDDO

        END SELECT

     END DO    !next iw
  END DO       !next idim


  IF(sourceunsplit) &
       CALL addsource2(qdt*(idimmax-idimmin+one)/ndim, &
       ixI^L,ixO^L,iws,qtC,wCT,qt,w)


END SUBROUTINE centdiff4

!=============================================================================
! end module vaccd
!#############################################################################
