!##############################################################################
! module vacusr - nul
!##############################################################################

!=============================================================================
SUBROUTINE specialini(ix^L,w)

! Initialize w for VACINI, user-defined

INCLUDE 'vacdef.f'

INTEGER:: ix^L
DOUBLE PRECISION:: w(ixG^T,nw)
!-----------------------------------------------------------------------------

CALL die('Special initial condition is not defined')
END SUBROUTINE specialini
!=============================================================================

!=============================================================================
SUBROUTINE specialbound(qt,ix^L,iw,iB,w)

! Calculates the boundary values in the iB-th boundary segment, user-defined

INCLUDE 'vacdef.f'

INTEGER:: ix^L,iw,iB
DOUBLE PRECISION:: qt,w(ixG^T,nw)
!-----------------------------------------------------------------------------

CALL die('Special boundary is not defined')
END SUBROUTINE specialbound
!=============================================================================

!=============================================================================
! specialsource -- for sources other than resistivity
! getdt_special -- for time step conditions other than CFL or resistivity
! specialeta    -- for non-constant resistivity with eqpar(eta_)<zero
!=============================================================================
SUBROUTINE specialsource(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

! Calculate w(iw)=w(iw)+qdt*SOURCE(wCT,iw) within ixO for all indices
! iw=iws(iiw) where the indirect index iiw=1..iws(niw_).
! wCT is at time qCT, while w is at time qt on input and qt+qdt on output.
!
! You may want to use second order accurate time integration if 
! "sourcesplit=T" and/or the "typefull='tvd'" method is used in the par-file.
!
! If the source needs wCT values outside ixO, ensurebound should be used:
!
! call ensurebound(dix,ixI^L,ixO^L,qtC,wCT)
!
! where "dix" is the number of extra layers needed, typically 1 or 2.

INCLUDE 'vacdef.f'

INTEGER:: ixI^L,ixO^L,iws(niw_)
DOUBLE PRECISION:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)

RETURN
END SUBROUTINE specialsource

!=============================================================================
SUBROUTINE getdt_special(w,ix^L)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the VACPHYS 
! module have already been called. 

INCLUDE 'vacdef.f'

DOUBLE PRECISION:: w(ixG^T,nw)
INTEGER:: ix^L
!-----------------------------------------------------------------------------

! Something like this
!
! dt_mine=...
! dt=min(dt,dt_mine)

RETURN
END SUBROUTINE getdt_special

!=============================================================================
SUBROUTINE specialeta(w,ix^L,idirmin)

! Set the common "eta" array for resistive MHD based on w and the common
! "current" variable which has components between idirmin and 3.
! If resistivity should be treated implicitly, the "gradeta" array 
! has to be set. It may use the closest neighbours of "w" (jx,ix,hx), 
! and/or local values (ix) of "current" only, i.e. it has to be compact.

INCLUDE 'vacdef.f'

DOUBLE PRECISION:: w(ixG^T,nw)
INTEGER:: ix^L,idirmin

CALL die('specialeta is not defined')
END SUBROUTINE specialeta
!=============================================================================

!=============================================================================
SUBROUTINE readfileini_special(w)

! Reads from unitini,filenameini in user-defined format.
! Check readfileini_asc and readfileini_bin in vacio.t on what should be done.

INCLUDE 'vacdef.f'

DOUBLE PRECISION:: w(ixG^T,nw)
!-----------------------------------------------------------------------------

CALL die('Special readfileini is not defined')
END SUBROUTINE readfileini_special
!=============================================================================
SUBROUTINE savefileout_special(qunit,w,ix^L)

! Save current results into filenameout in user-defined format.
! Check savefileout_asc and savefileout_bin in vacio.t on what should be done.

INCLUDE 'vacdef.f'

INTEGER:: qunit,ix^L
DOUBLE PRECISION:: w(ixG^T,nw)
!-----------------------------------------------------------------------------

CALL die('Special savefileout is not defined')
END SUBROUTINE savefileout_special
!=============================================================================
SUBROUTINE savefilelog_special(qunit,w,ix^L)

! Save user-defined log data into filename(filelog_) in user-defined format.
! Check savefilelog_default on opening the file etc.

INCLUDE 'vacdef.f'

INTEGER:: qunit,ix^L
DOUBLE PRECISION:: w(ixG^T,nw)
!-----------------------------------------------------------------------------

CALL die('Special savefilelog is not defined')
END SUBROUTINE savefilelog_special
!=============================================================================

!##############################################################################
! end module vacusr - nul
!##############################################################################
