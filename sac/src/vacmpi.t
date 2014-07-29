!=============================================================================
SUBROUTINE mpiinit

  ! Initialize MPI variables

  USE constants
  USE common_variables
  !----------------------------------------------------------------------------
  CALL MPI_INIT(ierrmpi)
  CALL MPI_COMM_RANK (MPI_COMM_WORLD, ipe, ierrmpi)
  CALL MPI_COMM_SIZE (MPI_COMM_WORLD, npe, ierrmpi)

  ! unset values for directional processor numbers
  npe^D=-1;
  ! default value for test processor
  ipetest=0

  RETURN
END SUBROUTINE mpiinit

!==============================================================================
SUBROUTINE mpifinalize

  USE constants
  USE common_variables

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
  CALL MPI_FINALIZE(ierrmpi)

  RETURN
END SUBROUTINE mpifinalize

!==============================================================================
SUBROUTINE ipe2ipeD(qipe,qipe^D)

  ! Convert serial processor index to directional processor indexes

  USE constants
  USE common_variables

  INTEGER:: qipe^D, qipe
  !-----------------------------------------------------------------------------
  qipe1 = qipe - npe1*(qipe/npe1)
  {qipe2 = qipe/npe1 - npe2*(qipe/(npe1*npe2)) ^NOONED}
  {qipe3 = qipe/(npe1*npe2)                    ^IFTHREED}

  RETURN
END SUBROUTINE ipe2ipeD

!==============================================================================
SUBROUTINE ipeD2ipe(qipe^D,qipe)

  ! Convert directional processor indexes to serial processor index

  USE constants
  USE common_variables

  INTEGER:: qipe^D, qipe
  !-----------------------------------------------------------------------------
  qipe = qipe1 {^NOONED + npe1*qipe2} {^IFTHREED + npe1*npe2*qipe3}

  RETURN
END SUBROUTINE ipeD2ipe

!==============================================================================
SUBROUTINE mpisetnpeDipeD(name)

  ! Set directional processor numbers and indexes based on a filename.
  ! The filename contains _np followed by np^D written with 2 digit integers.
  ! For example _np0203 means np1=2, np2=3 for 2D.

  USE constants
  USE common_variables
  CHARACTER(^LENNAME) :: name, nametail
  INTEGER:: i,qnpe^D
  LOGICAL:: npeDknown,npeDinname
  !-----------------------------------------------------------------------------

  oktest = INDEX(teststr,'mpisetnpeDipeD')>0
  IF(oktest)WRITE(*,*)'mpisetnpeDipeD ipe,name=',ipe,name

  ! Check if npe^D is already known
  npeDknown  = npe1>0
  IF(npedknown .AND. npe^D* /= npe)THEN
     WRITE(*,*)'npe=',npe,' /= product of npe^D=',npe^D
     CALL mpistop('ERROR in setnpeDipeD')
  ENDIF

  ! Check if npe^D is given in the name
  i=INDEX(name,'_np')+3
  npeDinname = i>3

  IF(.NOT.(npeDknown.OR.npeDinname))CALL mpistop( &
       'ERROR in setnpeDipeD: npeD is neither known nor given in name='//name)

  IF(npeDinname)THEN
     ! read npe^D from name
     READ(name(i:i+5),'(3i2)') qnpe^D
     i=i+2*^ND
     nametail=name(i:^LENNAME)
  ENDIF

  IF( npeDknown .AND. npeDinname )THEN
     ! Check agreement
     IF( qnpe^D/=npe^D|.OR. )THEN
        WRITE(*,*)'npe^D=',npe^D,' /= qnpe^D=',qnpe^D,' read from filename=',name
        CALL mpistop('ERROR in mpisetnpeDipeD')
     ENDIF
  ENDIF

  IF(npeDinname .AND. .NOT.npeDknown)THEN
     ! set npe^D based on name
     npe^D=qnpe^D;
     IF( npe^D* /= npe)THEN
        WRITE(*,*)'npe=',npe,' /= product of npe^D=',npe^D,&
             ' read from filename=',name
        CALL mpistop('ERROR in setnpeDipeD')
     ENDIF
  ENDIF

  ! Get directional processor indexes
  CALL ipe2ipeD(ipe,ipe^D)

  IF(npeDknown .AND. .NOT.npeDinname)THEN
     ! insert npe^D into name
     i=INDEX(name,'.')
     nametail=name(i:^LENNAME)
     WRITE(name(i:^LENNAME),"('_np',3i2.2)") npe^D
     i = i+3+2*^ND
  ENDIF

  ! insert ipe number into the filename
  WRITE(name(i:^LENNAME),"('_',i3.3,a)") ipe,nametail(1:^LENNAME-i-4)

  ! Set logicals about MPI boundaries for this processor
  {^DLOOP
  mpiupperB(^D)=ipe^D<npe^D-1
  mpilowerB(^D)=ipe^D>0 \}

  IF(oktest)WRITE(*,*)'mpisetnpeDipeD: ipe,npeD,ipeD,name=',ipe,npe^D,ipe^D,name

  RETURN
END SUBROUTINE mpisetnpeDipeD

!==============================================================================
SUBROUTINE mpineighbors(idir,hpe,jpe)

  ! Find the hpe and jpe processors on the left and right side of this processor 
  ! in direction idir. The processor cube is taken to be periodic in every
  ! direction.

  USE constants
  USE common_variables

  INTEGER :: idir,hpe,jpe,hpe^D,jpe^D
  !-----------------------------------------------------------------------------
  hpe^D=ipe^D-kr(^D,idir);
  jpe^D=ipe^D+kr(^D,idir);
  {^DLOOP
  IF(hpe^D<0)hpe^D=npe^D-1
  IF(jpe^D>=npe^D)jpe^D=0\}

  CALL ipeD2ipe(hpe^D,hpe)
  CALL ipeD2ipe(jpe^D,jpe)

  RETURN
END SUBROUTINE mpineighbors
!==============================================================================
SUBROUTINE mpigridsetup

  ! Distribute a grid of size nxall^D onto PE-s arranged in a cube of size npe^D

  USE constants
  USE common_variables
  !-----------------------------------------------------------------------------
!!!write(*,*)'nxall,npe=',nxall^D,npe^D

  ! Grid sizes on the processors (except for the last ones in some direction)
  ! This formula optimizes the load balance when nx^D is not a multiple of npe^D
  nxpe^D=(nxall^D-1)/npe^D+1; 

  ! Global grid indexes of the first grid point stored on this PE
  ixPEmin^D=ipe^D*nxpe^D+1;

  ! The last processors in a direction may have smaller grid sizes than nxpe
  {^DLOOP 
  IF(ipe^D < npe^D-1)THEN
     nx^D = nxpe^D
  ELSE
     nx^D = nxall^D - ixpemin^D + 1
  ENDIF
  \}

  ! Global grid index of the last grid point stored on this PE
  ixPEmax^D=ixPEmin^D+nx^D-1;

  RETURN
END SUBROUTINE mpigridsetup

!=============================================================================
SUBROUTINE mpireduce(a,mpifunc)

  ! reduce input for one PE 0 using mpifunc

  USE constants

  DOUBLE PRECISION :: a, alocal
  INTEGER          :: mpifunc, ierrmpi
  !----------------------------------------------------------------------------
  alocal = a
  CALL MPI_REDUCE(alocal,a,1,MPI_DOUBLE_PRECISION,mpifunc,&
       0,MPI_COMM_WORLD,ierrmpi)

  RETURN
END SUBROUTINE mpireduce

!==============================================================================
SUBROUTINE mpiallreduce(a,mpifunc)

  ! reduce input onto all PE-s using mpifunc

  USE constants

  DOUBLE PRECISION :: a, alocal
  INTEGER          :: mpifunc, ierrmpi
  !-----------------------------------------------------------------------------
  alocal = a
  CALL MPI_ALLREDUCE(alocal,a,1,MPI_DOUBLE_PRECISION,mpifunc,&
       MPI_COMM_WORLD,ierrmpi)

  RETURN
END SUBROUTINE mpiallreduce

!==============================================================================
SUBROUTINE mpiix(ix^D,jpe)

  ! Convert ix^D physical cell index on the global grid to local indexes 
  ! and set the processor number jpe to the processor that contains the cell

  USE constants
  USE common_variables
  INTEGER :: ix^D, jpe, jpe^D
  !-----------------------------------------------------------------------------

  ! Directional processor indexes
  jpe^D=(ix^D-ixMmin^D)/nxpe^D;

  ! Conversion to local index
  ix^D=ix^D-jpe^D*nxpe^D;

  ! Get MPI processor index
  CALL ipeD2ipe(jpe^D,jpe)

  RETURN
END SUBROUTINE mpiix

!==============================================================================
SUBROUTINE mpiixlimits(ix^L)

  ! Convert global index limits to local index limits for this PE

  USE constants
  USE common_variables
  INTEGER :: ix^L
  !-----------------------------------------------------------------------------
  {^DLOOP
  IF(ixmin^D > ixPEmax^D)THEN
     ixmin^D = nx^D
     ixmax^D = nx^D-1
  ELSEIF(ixmax^D < ixPEmin^D)THEN
     ixmax^D = 0
     ixmin^D = 1
  ELSE
     ixmin^D = MAX(ixmin^D,ixPEmin^D) - ixPEmin^D + 1
     ixmax^D = MIN(ixmax^D,ixPEmax^D) - ixPEmin^D + 1
  ENDIF
  \}

  RETURN
END SUBROUTINE mpiixlimits
!==============================================================================

SUBROUTINE mpistop(message)

  ! Stop MPI run in an orderly fashion

  USE constants
  USE common_variables

  CHARACTER(*) :: message
  INTEGER :: nerrmpi

  !------------------------------------------------------------------------------
  WRITE(*,*)'ERROR for processor',ipe,':'
  WRITE(*,*)message
  CALL MPI_abort(MPI_COMM_WORLD, nerrmpi, ierrmpi)

  STOP
END SUBROUTINE mpistop

!==============================================================================
SUBROUTINE mpibound(nvar,var)

  ! Fill in ghost cells of var(ixG,nvar) from other processors

  USE constants
  USE common_variables

  INTEGER :: nvar
  DOUBLE PRECISION :: var(ixG^T,nvar)

  ! processor indexes for left and right neighbors
  INTEGER :: hpe,jpe
  ! index limits for the left and right side mesh and ghost cells 
  INTEGER :: ixLM^L, ixRM^L, ixLG^L, ixRG^L
  LOGICAL :: periodic

  ! There can be at most 2 receives in any direction for each PE
  INTEGER :: nmpirequest, mpirequests(2)
  INTEGER :: mpistatus(MPI_STATUS_SIZE,2)
  COMMON /mpirecv/ nmpirequest,mpirequests,mpistatus
  !-----------------------------------------------------------------------------
  oktest=INDEX(teststr,'mpibound')>0
  IF(oktest)WRITE(*,*)'mpibound ipe,nvar,varold=',&
       ipe,nvar,var(ixtest^D,MIN(nvar,iwtest))

  {^DLOOP
  IF(npe^D>1)THEN
     nmpirequest =0
     mpirequests(1:2) = MPI_REQUEST_NULL

     periodic=typeB(1,2*^D)=='mpiperiod'

     ! Left and right side ghost cell regions (target)
     ixLG^L=ixG^L; ixLGmax^D=ixMmin^D-1;
     ixRG^L=ixG^L; ixRGmin^D=ixMmax^D+1;

     ! Left and right side mesh cell regions (source)
     ixLM^L=ixG^L; ixLMmin^D=ixMmin^D; ixLMmax^D=ixMmin^D+dixBmin^D-1;
     ixRM^L=ixG^L; ixRMmax^D=ixMmax^D; ixRMmin^D=ixMmax^D-dixBmax^D+1;

     ! Obtain left and right neighbor processors for this direction
     CALL mpineighbors(^D,hpe,jpe)

     IF(oktest)THEN
        WRITE(*,*)'mpibound ipe,idir=',ipe,^D
        WRITE(*,*)'mpibound ipe,ixLG=',ipe,ixLG^L
        WRITE(*,*)'mpibound ipe,ixRG=',ipe,ixRG^L
        WRITE(*,*)'mpibound ipe,ixLM=',ipe,ixLM^L
        WRITE(*,*)'mpibound ipe,ixRM=',ipe,ixRM^L
        WRITE(*,*)'mpibound ipe,hpe,jpe=',ipe,hpe,jpe
     ENDIF

     ! receive right (2) boundary from left neighbor hpe
     IF(mpilowerB(^D) .OR. periodic)CALL mpirecvbuffer(nvar,ixRM^L,hpe,2)
     ! receive left (1) boundary from right neighbor jpe
     IF(mpiupperB(^D) .OR. periodic)CALL mpirecvbuffer(nvar,ixLM^L,jpe,1)
     ! Wait for all receives to be posted
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
     ! Ready send left (1) boundary to left neighbor hpe
     IF(mpilowerB(^D) .OR. periodic)CALL mpisend(nvar,var,ixLM^L,hpe,1)
     ! Ready send right (2) boundary to right neighbor
     IF(mpiupperB(^D) .OR. periodic)CALL mpisend(nvar,var,ixRM^L,jpe,2)
     ! Wait for messages to arrive
     CALL MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)
     ! Copy buffer received from right (2) physical cells into left ghost cells
     IF(mpilowerB(^D) .OR. periodic)CALL mpibuffer2var(2,nvar,var,ixLG^L)
     ! Copy buffer received from left (1) physical cells into right ghost cells
     IF(mpiupperB(^D) .OR. periodic)CALL mpibuffer2var(1,nvar,var,ixRG^L)
  ENDIF
  \}

  IF(oktest)WRITE(*,*)'mpibound ipe,varnew=',ipe,var(ixtest^D,MIN(nvar,iwtest))

  RETURN
END SUBROUTINE mpibound

!==============================================================================
SUBROUTINE mpisend(nvar,var,ix^L,qipe,iside)

  ! Send var(ix^L,1:nvar) to processor qipe.
  ! jside is 0 for min and 1 for max side of the grid for the sending PE

  USE constants
  USE common_variables

  INTEGER :: nvar
  DOUBLE PRECISION :: var(ixG^T,nvar)
  INTEGER :: ix^L, qipe, iside, n, ix^D, ivar
  !----------------------------------------------------------------------------
  oktest = INDEX(teststr,'mpisend')>0

  n=0
  DO ivar=1,nvar
     {DO ix^DB=ixmin^DB,ixmax^DB;}
     n=n+1
     sendbuffer(n)=var(ix^D,ivar)
     {ENDDO^DLOOP\}
  END DO

  IF(oktest)THEN
     WRITE(*,*)'mpisend ipe-->qipe,iside,itag',ipe,qipe,iside,10*ipe+iside
     WRITE(*,*)'mpisend ipe,ix^L,var=',ipe,ix^L,var(ixtest^D,MIN(iwtest,nvar))
  ENDIF

  CALL MPI_RSEND(sendbuffer(1),n,MPI_DOUBLE_PRECISION,qipe,10*ipe+iside,&
       MPI_COMM_WORLD,ierrmpi)

  RETURN
END SUBROUTINE mpisend

!==============================================================================
SUBROUTINE mpirecvbuffer(nvar,ix^L,qipe,iside)

  ! receive buffer for a ghost cell region of size ix^L sent from processor qipe
  ! and sent from side iside of the grid

  USE constants
  USE common_variables

  INTEGER:: nvar, ix^L, qipe, iside, n

  INTEGER :: nmpirequest, mpirequests(2)
  INTEGER :: mpistatus(MPI_STATUS_SIZE,2)
  COMMON /mpirecv/ nmpirequest,mpirequests,mpistatus
  !----------------------------------------------------------------------------

  oktest = INDEX(teststr,'mpirecv')>0

  n = nvar* ^D&(ixmax^D-ixmin^D+1)*

  IF(oktest)WRITE(*,*)'mpirecv ipe<--qipe,iside,itag,n',&
       ipe,qipe,iside,10*qipe+iside,n

  nmpirequest = nmpirequest + 1
  CALL MPI_IRECV(recvbuffer(1,iside),n,MPI_DOUBLE_PRECISION,qipe,10*qipe+iside,&
       MPI_COMM_WORLD,mpirequests(nmpirequest),ierrmpi)

  RETURN
END SUBROUTINE mpirecvbuffer

!==============================================================================
SUBROUTINE mpibuffer2var(iside,nvar,var,ix^L)

  ! Copy mpibuffer(:,iside) into var(ix^L,1:nvar)
  USE constants
  USE common_variables

  INTEGER :: nvar
  DOUBLE PRECISION:: var(ixG^T,nvar)
  INTEGER:: ix^L,iside,n,ix^D,ivar
  !-----------------------------------------------------------------------------
  oktest = INDEX(teststr,'buffer2var')>0

  n=0
  DO ivar=1,nvar
     {DO ix^DB=ixmin^DB,ixmax^DB;}
     n=n+1
     var(ix^D,ivar)=recvbuffer(n,iside)
     {ENDDO^DLOOP\}
  END DO

  IF(oktest)WRITE(*,*)'buffer2var: ipe,iside,ix^L,var',&
       ipe,iside,ix^L,var(ixtest^D,MIN(iwtest,nvar))

  RETURN
END SUBROUTINE mpibuffer2var
