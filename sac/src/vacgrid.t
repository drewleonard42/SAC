!#############################################################################
! module vacgrid
! Subroutines for boundaries, grid, divergence of flux, gradients
! Also the limiter functions for TVD, TVDLF, TVDMU schemes

{INCLUDE:vacgrid.gencoord.t ^IFGEN}
!INCLUDE:vacgrid.setnozzle.t
!=============================================================================
SUBROUTINE boundsetup

  ! Variables describing boundaries
  !
  !    iB=1..nB                                    - boundary ID
  !    ixBmin(idim,iB):ixBmax(idim,iD)             - location of boundary
  !    idimB(iB)                                   - direction orthogonal to the 
  !                                                  boundary layer
  !    upperB(iB)                                  - .true. if boundary is at a
  !                                                  high index of the phys. grid
  !    typeB(iw,iB)                                - boundary type string

  USE constants
  USE common_variables

  INTEGER:: ix^L,iB,jB,iw,idim,idm,ixG^LIM(ndim),ixM^LIM(ndim)
  !-----------------------------------------------------------------------------

  oktest=INDEX(teststr,'boundsetup')>=1
  IF(oktest)WRITE(*,*)'BoundSetup'

  ! If ixBmax(1,1)==0, the user did not specify boundary regions. Setup default.
  ! A schematic picture for ndim=2, where the numbers are iB for each region, and
  ! *-s are the extra layers for qx. 
  !
  !                     ************
  !                     *1444444424*    ixGmax2
  !                     *4144444442*
  !                     *11      22*    ixMmax2
  !                     *11      22*
  !                     *11      22*    ixMmin2
  !                     *1333333323*
  !                     *3133333332*    ixGmin2
  !                     ************
  !                      i i    i i
  !                      x x    x x
  !                      G M    M G
  !                      m m    m m
  !                      i i    a a
  !                      n n    x x
  !                      1 1    1 1

  ^D&ixG^LIM(^D)=ixG^DL;
  ^D&ixM^LIM(^D)=ixM^DL;

  IF(ixBmax(1,1)==0)THEN
     nB=2*ndim
     DO iB=1,nB
        idim=(iB+1)/2
        idimB(iB) = idim
        upperB(iB)= 2*idim==iB
        DO idm=1,ndim
           ixB^LIM(idm,iB)=ixG^LIM(idm);
        END DO
        IF(upperB(iB))THEN
           ixBmin(idim,iB)=ixMmax(idim)+1
           ixBmax(idim,iB)=ixGmax(idim)
        ELSE
           ixBmin(idim,iB)=ixGmin(idim)
           ixBmax(idim,iB)=ixMmin(idim)-1
        END IF
     END DO
  ELSE
     ! Check if boundary regions are inside grid and in between mesh and grid
     DO iB=1,nB
        ix^L=ixB^LIM(^D,iB);
        {^IFMPI
        ! Convert global limits into local grid limits 
        ix^L=ix^L-ixpemin^D+1;
        ! Cut off extra cells at MPI boundaries
        ^D&IF(ipe^D>0)ixmin^D=MAX(ixmin^D,ixGmin^D)\
        ^D&IF(ipe^D<npe^D-1)ixmax^D=MIN(ixmax^D,ixGmax^D)\
        ! Change the global ixB^LIM(^D,iB) to the local index limits
        ^D&ixBmin(^D,iB)=ixmin^D;
        ^D&ixBmax(^D,iB)=ixmax^D;
        ! Check if local boundary region is empty
        IF(ixmin^D>ixmax^D|.OR.)THEN
           ! The limits do not contain cells, change the boundary type to 'mpi'
           typeB(1:nw,iB)='mpi'
           CYCLE
        ENDIF
        |\}
        IF(ixmin^D<ixGmin^D.OR.ixmax^D>ixGmax^D|.OR.) THEN
           WRITE(*,*)'Error for boundary region iB,ixL=',iB,ix^L
           CALL die('Error in BoundSetup: Boundary region is outside grid')
        ENDIF
        SELECT CASE(idimB(iB))
           {CASE(^D)
           IF(upperB(iB))THEN
              IF(ixmin^D-1/=ixMmax^D.OR.ixmax^D/=ixGmax^D)WRITE(*,*)&
                   'Warning in BoundSetup: Boundary does not fit, iB:',iB
           ELSE
              IF(ixmax^D+1/=ixMmin^D.OR.ixmin^D/=ixGmin^D)WRITE(*,*)&
                   'Warning in BoundSetup: Boundary does not fit, iB:',iB
           ENDIF\}
        END SELECT
     END DO
  END IF

  ! Identify the periodic pairs if they are not defined in boundlist
  ! Check type, direction, orientation, and size before matching pairs
  DO iB=1,nB
     IF(typeB(1,iB)=='periodic'.AND.ipairB(iB)==0)THEN
        DO jB=iB+1,nB
           IF(typeB(1,jB)=='periodic'.AND.ipairB(jB)==0.AND.&
                idimB(iB)==idimB(jB).AND.(upperB(iB).NEQV.upperB(jB)).AND.&
                {ixBmax(^D,iB)-ixBmin(^D,iB)==ixBmax(^D,jB)-ixBmin(^D,jB)|.AND.})THEN
              ipairB(iB)=jB; ipairB(jB)=iB
           ENDIF
        END DO
        IF(ipairB(iB)==0)CALL die('Error in BoundSetup: No periodic pair')
     END IF
  END DO

  {^IFMPI
  DO iB=1,nB
     ! Change boundary type if processor is not at the edge or if this is a
     ! periodic boundary and there are more than 1 processors in this direction
     SELECT CASE(idimB(iB))
        {CASE(^D)
        IF(typeB(1,iB)=='periodic' .AND. npe^D>1) THEN
           typeB(1:nw,iB)='mpiperiod'
        ELSE IF(upperB(iB) .AND. ipe^D<npe^D-1 .OR. &
             .NOT.upperB(iB) .AND. ipe^D>0) THEN
           typeB(1:nw,iB)='mpi'
        ENDIF\}
     END SELECT
  END DO
  }


  ! symm0 means zero orthogonal flux via the boundary 
  DO iw=1,nw
     DO iB=1,nB
        IF(typeB(iw,iB)=='symm0') nofluxB(iw,idimB(iB))=.TRUE.
     ENDDO
  ENDDO

  IF(oktest)THEN
     DO iB=1,nB
        {WRITE(unitterm,*)'ipe=',ipe ^IFMPI}
        WRITE(unitterm,*)'iB,idimB,upperB,typeB:',iB,idimB(iB),upperB(iB),&
             ' ',typeB(1,iB)
        WRITE(unitterm,*)'ixBmin:',(ixBmin(idm,iB),idm=1,ndim)
        WRITE(unitterm,*)'ixBmax:',(ixBmax(idm,iB),idm=1,ndim)
     END DO
     DO iB=1,nB
        IF(typeB(1,iB)=='periodic')WRITE(*,*)'Pairs:',iB,ipairB(iB)
     END DO
  END IF

  RETURN
END SUBROUTINE boundsetup

!=============================================================================
SUBROUTINE ensurebound(dix,ixI^L,ixO^L,qt,w)

  ! Check if there is enough information for calculating derivatives.
  ! The requirement is that ixI is wider than ixO by dix. 
  ! Adjust ixI and ixO. Call getboundary if needed.

  USE constants
  USE common_variables

  INTEGER:: dix,ixI^L,ixO^L
  DOUBLE PRECISION:: qt,w(ixG^T,nw)
  !-----------------------------------------------------------------------------

  oktest=INDEX(teststr,'ensurebound')>0
  IF(oktest)WRITE(*,*)'EnsureBound dix,ixI,ixO:',dix,',',ixI^L,',',ixO^L

  ! Check wether ixO+dix is within the grid
  IF(ixG^L^LTixO^L^LADDdix|.OR.|.OR.)THEN
     ixO^L=ixM^L;
  ENDIF
  ! Check whether ixI is wider than ixO by at least dix otherwise getboundary
  IF(ixI^L^LTixO^L^LADDdix|.OR.|.OR.)THEN
     ixI^L=ixG^L;
     CALL getboundary(qt,1,nw,1,ndim,w)
     IF(oktest)WRITE(*,*)'calling getboundary'
  END IF

  IF(oktest)WRITE(*,*)'Final       dix,ixI,ixO:',dix,',',ixI^L,',',ixO^L

  RETURN
END SUBROUTINE ensurebound

!=============================================================================
SUBROUTINE getboundary(qt,iw^LIM,idim^LIM,w)

  USE constants
  USE common_variables

  INTEGER:: iw^LIM,idim^LIM
  DOUBLE PRECISION:: qt,w(ixG^T,1:nw)
  INTEGER:: ix,ix^D,ixe,ixf,ix^L,ixpair^L,idim,iw,iB
  INTEGER:: iwv,jdim
  DOUBLE PRECISION:: coeffnormal,coefftransv
  LOGICAL:: initialized


  INTEGER:: ixee
  INTEGER:: ixb^L

  initialized = .FALSE.
  !-----------------------------------------------------------------------------


  ixb^L=ixG^LL;

  oktest=INDEX(teststr,'getboundary')>=1
  IF(oktest)WRITE(*,*)'GetBoundary it,step:',it,step
  IF(oktest)WRITE(*,*)'GetBoundary wold:',w(ixtest^D,iwtest)

  IF(extraB)CALL specialbound(qt,ixM^L,0,0,w)
  IF(smallfix)CALL keeppositive(ixM^L,w)

  {^IFMPI
  ! Get boundaries from other PE-s
  CALL mpibound(nw,w)
  }

  iB=0
  DO
     iB=iB+1
     IF(oktest)WRITE(*,*)'iB  :',iB
     IF(iB>nB) EXIT
     idim=idimB(iB)
     ! Only boundary segments in the required direction(s) are filled in
     IF(idimmin>idim.OR.idimmax<idim) CYCLE

     ix^L=ixB^LIM(^D,iB);

     ! The possibly shifted coordinates parallel to the boundary layer 
     ! are defined by the PAIR of the periodic/overlapping boundary.
     ! Put the location of the source of information into ixpairL.
     IF(ipairB(iB)>0)THEN
        ixpair^L=ixB^LIM(^D,ipairB(iB));
        SELECT CASE(idim)
           {CASE(^D)
           IF(upperB(iB))THEN
              ixpair^LIM^D=ixpair^LIM^D+dixB^LIM^D;
           ELSE
              ixpair^LIM^D=ixpair^LIM^D-dixB^LIM^D;
           ENDIF
           \}
        END SELECT
     ENDIF

     DO iw= iw^LIM
        IF(oktest)WRITE(*,*)'  iw:',iw
        IF(oktest)WRITE(*,*)'typeB(iw,iB):',typeB(iw,iB)
        SELECT CASE (typeB(iw,iB))



        CASE('contCD4')
           SELECT CASE(idim)
              {   CASE(^D)

              IF(upperB(iB))THEN


                 ixe=ixmin^D-2
                 ixee=ixmin^D-3

!HPF$ INDEPENDENT
                 ix= ixmin^D
                 w(ix^D%ix^S,iw)=w(ixe^D%ix^S,iw)

                 ix= ixmax^D
                 w(ix^D%ix^S,iw)=w(ixee^D%ix^S,iw)

              ELSE
                 ixe=ixmax^D+3
                 ixee=ixmax^D+2

!HPF$ INDEPENDENT
                 ix= ixmin^D
                 w(ix^D%ix^S,iw)=w(ixe^D%ix^S,iw)

                 ix= ixmax^D
                 w(ix^D%ix^S,iw)=w(ixee^D%ix^S,iw)

              ENDIF


              }
           END SELECT

        CASE('zero')
           SELECT CASE(idim)
              {   CASE(^D)

              IF(upperB(iB))THEN

                 CALL primitive(ixG^LL,w)

                 ixe=ixmin^D-2
                 ixee=ixmin^D-3

                 ixbmax^D=ixmax^D
                 ixbmin^D=ixee


!HPF$ INDEPENDENT
                 ix= ixmin^D

                 w(ix^D%ix^S,iw)=w(ixe^D%ix^S,iw)

                 ix= ixmax^D
                 w(ix^D%ix^S,iw)=w(ixee^D%ix^S,iw)

                 CALL conserve(ixG^LL,w)

              ELSE

                 CALL primitive(ixG^LL,w)

                 ixe=ixmax^D+3
                 ixee=ixmax^D+2


!HPF$ INDEPENDENT
                 ix= ixmin^D
                 w(ix^D%ix^S,iw)=w(ixe^D%ix^S,iw)

                 ix= ixmax^D
                 w(ix^D%ix^S,iw)=w(ixee^D%ix^S,iw)


                 CALL conserve(ixG^LL,w)

              ENDIF



              }



           END SELECT



        CASE('cont','fixed')
           ! For 'cont' copy w at the edge into the whole boundary region.
           ! Fot 'fixed' copy w first, later use values stored in fixB.
           ! For fullgridini=T store the w values read from the file in fixB.
           SELECT CASE(idim)
              {CASE(^D)
              IF(upperB(iB))THEN
                 ixe=ixmin^D-1
              ELSE
                 ixe=ixmax^D+1
              ENDIF
              IF(fixedB(iw,iB))THEN
!HPF$ INDEPENDENT
                 {DO ix^DD=ixmin^DD,ixmax^DD\}
                 w(ix^DD,iw)=fixB^D(ix^D-ixe^D%ix^DD,iw)
                 {ENDDO^DD&\} 
              ELSE IF(typeB(iw,iB)=='cont' .OR. .NOT.fullgridini) THEN
!HPF$ INDEPENDENT
                 DO ix= ix^DL
                    w(ix^D%ix^S,iw)=w(ixe^D%ix^S,iw)
                 END DO
              END IF\}
           END SELECT
        CASE ('cont1','fixed1','grad1')
           ! First order extrapolation from edge and inner edge to the boundary
           ! 'cont1'  extrapolates in every time step, can be numericly unstable.
           ! 'fixed1' extrapolates first, stores VALUES into fixB, then restores.
           ! 'grad1'  extrapolates first, stores DIFFERENCES into fixB, then 
           !          adds the stored differences to the current edge value.
           SELECT CASE(idim)
              {CASE(^D)
              IF(upperB(iB))THEN
                 ixe=ixmin^D-1; ixf=ixe-1
              ELSE
                 ixe=ixmax^D+1; ixf=ixe+1
              ENDIF
              IF(fixedB(iw,iB))THEN
                 IF(typeB(iw,iB)=='grad1')THEN
!HPF$ INDEPENDENT
                    {DO ix^DD=ixmin^DD,ixmax^DD\}
                    w(ix^DD,iw)=fixB^D(ix^D-ixe^D%ix^DD,iw)+w(ixe^D%ix^DD,iw)
                    {ENDDO^DD&\}
                 ELSE
!HPF$ INDEPENDENT
                    {DO ix^DD=ixmin^DD,ixmax^DD\}
                    w(ix^DD,iw)=fixB^D(ix^D-ixe^D%ix^DD,iw)
                    {ENDDO^DD&\}
                 ENDIF
              ELSE IF(typeB(iw,iB)=='cont1'.OR. .NOT.fullgridini)THEN !HPF_
                 !HPF_ endif
                 !HPF_ if(.not.fixedB(iw,iB).and.&
                 !HPF_    (typeB(iw,iB)=='cont1'.or. .not.fullgridini))then
!HPF$ INDEPENDENT
                 DO ix= ix^DL
                    w(ix^D%ix^S,iw)=&
                         (ABS(ix-ixe)+1)*w(ixe^D%ix^S,iw)- &
                         ABS(ix-ixe)   *w(ixf^D%ix^S,iw)  
                 END DO
              END IF\}
           END SELECT
        CASE('periodic')
           ! Update boundary by translation of w by width of mesh (and shift)
           w(ix^S,iw)=w(ixpair^S,iw)
        CASE('symm','symm0','asymm')
           ! Reflect w into the boundary region, multiply by -1 for "asymm"
           ! In generalized coordinates take into account the other vector
           ! components for vector variables. The symmetry of the transverse 
           ! component is based on typeB for the jdim=idim+1 -th component.
           ! ixe is used for the reflection, normal vectors are taken at ixf+1/2
           IF(gencoord.AND.vectoriw(iw)>=0)THEN
              ! Determine direction of vector component, symmetry coefficients
              ! for normal and transverse vector components
              iwv=vectoriw(iw); jdim=idim+1-(idim/ndim)*ndim 
              coeffnormal=1; IF(typeB(iwv+idim,iB)=='asymm') coeffnormal=-1
              coefftransv=1; IF(typeB(iwv+jdim,iB)=='asymm') coefftransv=-1
           ENDIF
           SELECT CASE(idim)
              {CASE(^D)
              IF(upperB(iB))THEN
                 ixe=2*ixmin^D-1; ixf=ixmin^D-1
              ELSE
                 ixe=2*ixmax^D+1; ixf=ixmax^D
              ENDIF
              IF(gencoord.AND.vectoriw(iw)>=0)THEN
!HPF$ INDEPENDENT
                 DO ix= ix^DL
                    w(ix^D%ix^S,iw)=zero
                    DO jdim=1,ndim
                       w(ix^D%ix^S,iw)=w(ix^D%ix^S,iw)+&
                            normalC(ixf^D%ix^S,idim,jdim)*w(ixe-ix^D%ix^S,iwv+jdim)
                    END DO
                    w(ix^D%ix^S,iw)=w(ix^D%ix^S,iw)*&
                         normalC(ixf^D%ix^S,idim,iw-iwv)*(coeffnormal-coefftransv)&
                         +w(ixe-ix^D%ix^S,iw)*coefftransv
                 END DO
              ELSE
!HPF$ INDEPENDENT
                 DO ix= ix^DL
                    w(ix^D%ix^S,iw)=w(ixe-ix^D%ix^S,iw)
                 END DO
                 IF(typeB(iw,iB)=='asymm') w(ix^S,iw)=-w(ix^S,iw)
              ENDIF
              \}
           END SELECT
        CASE('special')
           ! Skip special now, we do it after normal boundary type variables
           !HPF_ if(.false.)write(*,*)'Avoiding xlhpf90 compiler bug'
           {^IFMPI
        CASE('mpi','mpiperiod')
           ! This boundary is handled by MPI \}
        CASE default
           WRITE(uniterr,*)'Error in GetBoundary: boundary type(', &
                iw,iB,')=',typeB(iw,iB),' is not implemented!'
           CALL die(' ')
        END SELECT ! typeB(iw,iB)
     END DO ! next iw
     ! Do special boundaries
     DO iw= iw^LIM
        IF(oktest)WRITE(*,*)'special, iw:',iw
        IF(typeB(iw,iB).EQ.'special')CALL specialbound(qt,ix^L,iw,iB,w)
     END DO ! next iw
  END DO ! next iB

  ! Fixed boundaries (fixed,fixed1) or fixed gradients (grad1) are stored into
  ! fixB after all boundaries have been updated.
  ! This needs to be done in the very first time step only.
  IF(.NOT.initialized)THEN
     initialized=.TRUE.
     DO iB= 1,nB
        ix^L=ixB^LIM(^D,iB);
        DO iw= iw^LIM
           IF( (typeB(iw,iB)=='fixed'.OR.typeB(iw,iB)=='fixed1'.OR.&
                typeB(iw,iB)=='grad1') .AND. .NOT.fixedB(iw,iB))THEN
              fixedB(iw,iB)=.TRUE.
              SELECT CASE(idimB(iB))
                 {CASE(^D)
                 IF(upperB(iB))THEN
                    ixe=ixmin^D-1
                 ELSE
                    ixe=ixmax^D+1
                 ENDIF
                 IF(typeB(iw,iB)=='grad1')THEN
!HPF$ INDEPENDENT
                    {DO ix^DD= ixmin^DD,ixmax^DD\}
                    fixB^D(ix^D-ixe^D%ix^DD,iw)=w(ix^DD,iw)-w(ixe^D%ix^DD,iw)
                    {ENDDO^DD&\}
                 ELSE
!HPF$ INDEPENDENT
                    {DO ix^DD= ixmin^DD,ixmax^DD\}
                    fixB^D(ix^D-ixe^D%ix^DD,iw)=w(ix^DD,iw)
                    {ENDDO^DD&\}
                 ENDIF
                 \}
              END SELECT
           END IF
        END DO ! iw
     END DO    ! iB
  END IF

  IF(oktest)WRITE(*,*)'GetBoundary wnew:',w(ixtest^D,iwtest)

  RETURN 
END SUBROUTINE getboundary

!=============================================================================
SUBROUTINE setnoflux(iw,idim,ix^L,fRC,ixR^L,fLC,ixL^L)

  ! Set flux in direction idim to zero for variable iw if typeB is 'symm0'
  ! in a boundary region

  USE constants
  USE common_variables

  INTEGER:: iw,idim,ix^L,ixL^L,ixR^L
  DOUBLE PRECISION:: fRC(ixG^T), fLC(ixG^T)
  INTEGER:: iB,ixe,ixB^L

  !-----------------------------------------------------------------------------

  oktest=INDEX(teststr,'setnoflux')>=1

  DO iB=1,nB
     IF(typeB(iw,iB)=='symm0'.AND.idimB(iB)==idim)THEN
        ixB^L=ixB^LIM(^D,iB);
        ! Calculate edge index and set the appropriate flux to zero
        SELECT CASE(idim)
           {CASE(^D)
           IF(upperB(iB))THEN
              ixe=ixBmin^D-1+ixRmin^D-ixmin^D
              IF(oktest)WRITE(*,*)'Setnoflux it,idim,iw,ixe:', &
                   it,idim,iw,ixe,fRC(ixe^D%ixtest^DD)
              fRC(ixe^D%ixB^S)=zero
           ELSE
              ixe=ixBmax^D+1+ixLmin^D-ixmin^D
              IF(oktest)WRITE(*,*)'Setnoflux it,idim,iw,ixe:',&
                   it,idim,iw,ixe,fLC(ixe^D%ixtest^DD)
              fLC(ixe^D%ixB^S)=zero
           ENDIF
           \}
        END SELECT
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE setnoflux

!=============================================================================
SUBROUTINE gridsetup1

  ! Cartesian or polar grid. Determine x at the boundaries.
  ! Determine often needed combinations of x, such as dx or dvolume.
  ! Determine variables for axial symmetry
  !
  ! ixe          - edge coordinate of the grid touching the boundary region
  ! ixf          - coordinate inside of ixe
  ! qx           - x with an extended index range for calculation of dx

  USE constants
  USE common_variables

  INTEGER:: ix^L,hx^L,jx^L
  INTEGER:: ix,ixe,ixf,idim,jdim
  DOUBLE PRECISION:: qx(IXG^LL^LADD1:,ndim)
  DOUBLE PRECISION:: r(IXGLO1-1:IXGHI1+1),rC(IXGLO1-1:IXGHI1+1)

  !-----------------------------------------------------------------------------

  oktest=INDEX(teststr,'gridsetup')>=1
  IF(oktest)WRITE(*,*)'GridSetup1'

  ! Calculate qx in the boundary layers by linearly extrapolating x from
  ! the touching edge of the grid (ixe) and the inner edge (ixf).

  {^IFMPI
  ! Fill in ghost cells for x from MPI neighbors
  CALL mpibound(ndim,x)
  }
  qx(ixG^LL^LADD1:,1:ndim)=zero
  qx(ixG^S,1:ndim) = x(ixG^S,1:ndim)
  DO idim=1,ndim
     ix^L=ixG^L^LADD1;
     SELECT CASE(idim)
        {CASE(^D)
        ! First do the upper layers
        ixmax^D=ixGmax^D+1; ixmin^D=ixMmax^D+1
        IF(fullgridini .OR.mpiupperB(^D)^IFMPI) ixmin^D=ixGmax^D+1
        ixe=ixmin^D-1; ixf=ixe-1
!!!forall replaced by do for sake of sequentializing (Adaptor)
        DO jdim=1,ndim
           DO ix= ix^DL
              qx(ix^D%ix^S,jdim)=(1+ABS(ixe-ix))*qx(ixe^D%ix^S,jdim)- &
                   ABS(ixe-ix) *qx(ixf^D%ix^S,jdim)
           END DO
        END DO
        ! Next do the lower layers
        ixmin^D=ixGmin^D-1; ixmax^D=ixMmin^D-1
        IF(fullgridini .OR.mpilowerB(^D)^IFMPI) ixmax^D=ixGmin^D-1
        ixe=ixmax^D+1; ixf=ixe+1
        DO jdim=1,ndim
           DO ix= ix^DL
              qx(ix^D%ix^S,jdim)=(1+ABS(ixe-ix))*qx(ixe^D%ix^S,jdim)- &
                   ABS(ixe-ix) *qx(ixf^D%ix^S,jdim)
           END DO
        END DO\}
     END SELECT
  ENDDO

  x(ixG^S,1:ndim)=qx(ixG^S,1:ndim)

  ! Loop with ^D instead of idim to avoid an SP2 xlphf90 compiler error
  ! if qx is distributed. But it should not be
  !{
  !^D&jx^L=ixG^L+kr(^D,^DD); hx^L=ixG^L-kr(^D,^DD);
  !   dx(ixG^S,^D)=half*(qx(jx^S,^D)-qx(hx^S,^D))
  !   if(oktest)write(*,*)'dx,qxj,qxh:',dx(ixtest^DD,^D),&
  !      qx(ixtest^DD+kr(^D,^DD),^D),qx(ixtest^DD-kr(^D,^DD),^D)
  !\}

  DO idim=1,ndim
     jx^L=ixG^L+kr(idim,^D); hx^L=ixG^L-kr(idim,^D);
     dx(ixG^S,idim)=half*(qx(jx^S,idim)-qx(hx^S,idim))
  END DO

  ! Calculate geometrical factors for axial symmetry based on Boris FCT.
  ! Fluxes are multiplied by areaC. The cell volume is proportional to areadx.
  ! Gradient type source terms are multiplied by areaside = darea/areadx.

  IF(oktest)WRITE(*,*)'Start calculating geometrical factors'



  IF(oktest)WRITE(*,*)'Start calculating cell volumes'

  ! Calculate volume of cells and total volume of mesh
  IF(typeaxial=='slab')THEN
     dvolume(ixG^S)= ^D&dx(ixG^S,^D)*
  ELSE
     FORALL(ix= ixG^LIM1:) dvolume(ix,ixG^SE)= ^D&areadx(ix)^%1dx(ix,ixG^SE,^D)*
  ENDIF
  volume=SUM(dvolume(ixM^S))
  {^IFMPI
  ! Add up volumes
  CALL mpiallreduce(volume,MPI_SUM)
  ! Correct volumes in 2nd ghost cells from neighboring processors
  CALL mpibound(1,dvolume)
  }

  ! For polar grid dx_phi=r*dphi. 
  IF(polargrid)dx(ixG^S,pphi_)=x(ixG^S,r_)*dx(ixG^S,pphi_)

  IF(oktest)WRITE(*,*)'Finish GridSetup1'

  RETURN 
END SUBROUTINE gridsetup1

!=============================================================================

SUBROUTINE gradient4(realgrad,q,ix^L,idim,gradq)

  USE constants
  USE common_variables

  LOGICAL:: realgrad
  INTEGER:: ix^L,idim
  DOUBLE PRECISION:: q(ixG^T),gradq(ixG^T)
  INTEGER:: kx^L,jx^L,hx^L,gx^L
  INTEGER:: minx1^D,maxx1^D,k
  !-----------------------------------------------------------------------------

  !SHIFT
  kx^L=ix^L+2*kr(idim,^D);
  !SHIFT MORE
  jx^L=ix^L+kr(idim,^D);
  !SHIFT MORE
  hx^L=ix^L-kr(idim,^D);
  !SHIFT MORE
  gx^L=ix^L-2*kr(idim,^D);

  !SHIFT BEGIN
  gradq(ix^S)=-(q(kx^S)-8.D0*(q(jx^S)-q(hx^S))-q(gx^S))/dx(ix^S,idim)/12.D0
  !SHIFT END

  minx1^D=ixmin^D+kr(idim,^D);
  maxx1^D=ixmax^D-kr(idim,^D);

  DO k=0,1  !left-right bc

     IF (typeB(1,2*idim-1+k) .NE. 'mpi') THEN
        IF (upperB(2*idim-1+k)) THEN
           SELECT CASE(idim)
              {   CASE(^D)
              gradq(ixmax^D^D%ix^S)=0.d0
              gradq(maxx1^D^D%ix^S)=0.d0
              }
           END SELECT
        ELSE
           SELECT CASE(idim)
              {   CASE(^D)
              gradq(ixmin^D^D%ix^S)=0.d0
              gradq(minx1^D^D%ix^S)=0.d0
              }
           END SELECT
        ENDIF
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE gradient4


!=============================================================================
SUBROUTINE laplace4(q,ix^L,laplaceq)

!!! This subroutine should not use tmp or tmp2

  ! Calculate 4th order laplace of q within ixL in Cartesian direction idim

!!! We assume uniform Cartesian grid in slab symmetry for now

  USE constants
  USE common_variables

  INTEGER:: ix^L
  DOUBLE PRECISION:: q(ixG^T),laplaceq(ixG^T)

  INTEGER:: idim,kx^L,jx^L,hx^L,gx^L
  !-----------------------------------------------------------------------------

  IF(gencoord)CALL die('Error: laplace4 does not work for gen.coords!')
  IF(typeaxial/='slab')&
       CALL die('Error: laplace4 does not work for axial symmetry yet!')

  oktest= INDEX(teststr,'lpalace')>=1

  laplaceq(ix^S)=zero

  DO idim=1,ndim
     !SHIFT
     kx^L=ix^L+2*kr(idim,^D); 
     !SHIFT MORE
     jx^L=ix^L+kr(idim,^D); 
     !SHIFT MORE
     hx^L=ix^L-kr(idim,^D);
     !SHIFT MORE
     gx^L=ix^L-2*kr(idim,^D);

     !SHIFT BEGIN
     laplaceq(ix^S)=laplaceq(ix^S)+&
          (q(kx^S)+q(gx^S)+30*q(ix^S)-16*(q(jx^S)+q(hx^S)))/dx(ix^S,idim)**2/12
     !SHIFT END

     IF(oktest)WRITE(*,*)'idim,q(kx,jx,ix,hx,gx):',idim,&
          q(ixtest^D+2*kr(idim,^D)),q(ixtest^D+kr(idim,^D)),q(ixtest^D),&
          q(ixtest^D-kr(idim,^D)),q(ixtest^D-2*kr(idim,^D))
  ENDDO

  IF(oktest)WRITE(*,*)'laplaceq:',laplaceq(ixtest^D)

  RETURN
END SUBROUTINE laplace4

!=============================================================================
! end module vacgrid
!##############################################################################



