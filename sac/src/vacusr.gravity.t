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
SUBROUTINE addsource_grav(qdt,ixI^L,ixO^L,iws,qtC,w,qt,wnew)

  ! Add gravity source calculated from w to wnew within ixO for all variables 
  ! in iws. w is at time qtC, wnew is advanced from qt to qt+qdt.

  USE constants
  USE common_variables

  INTEGER::          ixI^L,ixO^L,iws(niw_)
  DOUBLE PRECISION:: qdt,qtC,qt,w(ixG^T,nw),wnew(ixG^T,nw)
  INTEGER:: iiw,iw,idim
!!! ! For a spatially varying gravity define the common grav array
!!! double precision:: grav(ixG^T,ndim)
!!! common /gravity/ grav

  !-----------------------------------------------------------------------------

!!! ! If grav needs to be calculated only once do it for the whole grid
!!! if(it==itmin)call setgrav(w,ixG^L,ixG^L,grav)
!!! ! Otherwise call setgrav in every time step
!!! call setgrav(w,ixI^L,ixO^L,grav)

  ! add sources from gravity
  DO iiw=1,iws(niw_); iw=iws(iiw)
     SELECT CASE(iw)
     CASE(m^D_)
        ! dm_i/dt= +rho*g_i
        idim=iw-m0_
        IF(ABS(eqpar(grav0_+idim))>smalldouble) &
             wnew(ixO^S,m0_+idim)=wnew(ixO^S,m0_+idim)+ &
             qdt*eqpar(grav0_+idim)*(w(ixO^S,rho_))

        !          wnew(ixO^S,m0_+idim)=wnew(ixO^S,m0_+idim)+ &
        !              qdt*eqpar(grav0_+idim)*(w(ixO^S,rho_)+w(ixO^S,rhob_))

!!! ! For a spatially varying gravity use instead of the above lines
!!! wnew(ixO^S,m0_+idim)=wnew(ixO^S,m0_+idim)+ &
!!!    qdt*grav(ixO^S,idim)*(w(ixO^S,rho_)+w(ixO^S,rhob_))

     CASE(e_)
        ! de/dt= +g_i*m_i
        DO idim=1,ndim
           IF(ABS(eqpar(grav0_+idim))>smalldouble) &
                wnew(ixO^S,ee_)=wnew(ixO^S,ee_)+ &
                qdt*eqpar(grav0_+idim)*w(ixO^S,rho_)*w(ixO^S,m0_+idim)/(w(ixO^S,rho_)+w(ixO^S,rhob_))

           !            wnew(ixO^S,ee_)=wnew(ixO^S,ee_)+ &
           !               qdt*eqpar(grav0_+idim)*w(ixO^S,m0_+idim)

!!! ! For a spatially varying gravity use instead of the above lines
!!! wnew(ixO^S,ee_)=wnew(ixO^S,ee_)+ &
!!!    qdt*grav(ixO^S,idim)*w(ixO^S,m0_+idim)

        END DO
     END SELECT ! iw
  END DO        ! iiw

  RETURN
END SUBROUTINE addsource_grav
!=============================================================================
!!! subroutine setgrav(w,ixI^L,ixO^L,grav)

! Set the gravitational acceleration within ixO based on x(ixI,ndim) 
! and/or w(ixI,nw)

!!! include 'vacdef.f90'

!!! double precision:: w(ixG^T,nw),grav(ixG^T,ndim)
!!! integer:: ixI^L,ixO^L
!----------------------------------------------------------------------------
!!! return
!!! end
!=============================================================================

SUBROUTINE getdt_grav(w,ix^L)

  USE constants
  USE common_variables

  DOUBLE PRECISION :: w(ixG^T,nw)
  INTEGER :: ix^L,idim
  DOUBLE PRECISION, SAVE :: dtgrav
 
!!! ! For spatially varying gravity you need a common grav array
!!! double precision:: grav(ixG^T,ndim)
!!! common/gravity/grav

  !----------------------------------------------------------------------------

  oktest=INDEX(teststr,'getdt')>=1

  IF(it==itmin)THEN
     ! If gravity is descibed by the equation parameters, use this:
     dtgrav=bigdouble
     DO idim=1,ndim
        IF(ABS(eqpar(grav0_+idim))>zero)&
             dtgrav=MIN(dtgrav,&
             one/SQRT(MAXVAL(ABS(eqpar(grav0_+idim))/dx(ixM^S,1:ndim))))
     ENDDO
!!! ! For spatially varying gravity use this instead of the lines above:
!!! call setgrav(w,ixG^L,ixM^L,grav)
!!! ! If gravity does not change with time, calculate dtgrav here:
!!! dtgrav=one/sqrt(maxval(abs(grav(ixM^S,1:ndim))/dx(ixM^S,1:ndim)))
  ENDIF

!!! ! If gravity changes with time, calculate dtgrav here:
!!! dtgrav=one/sqrt(maxval(abs(grav(ixM^S,1:ndim))/dx(ixM^S,1:ndim)))

  {^IFMPI CALL mpiallreduce(dtgrav,MPI_MIN)}

  ! limit the time step
  dt=MIN(dt,dtgrav)
  IF(oktest)WRITE(*,*)'Gravity limit for dt:',dtgrav

  RETURN
END SUBROUTINE getdt_grav

!=============================================================================
