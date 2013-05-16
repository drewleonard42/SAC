!=============================================================================
subroutine mpiinit

! Initialize MPI variables

include 'vacdef.f'
!----------------------------------------------------------------------------
call MPI_INIT(ierrmpi)
call MPI_COMM_RANK (MPI_COMM_WORLD, ipe, ierrmpi)
call MPI_COMM_SIZE (MPI_COMM_WORLD, npe, ierrmpi)

! unset values for directional processor numbers
npe1=-1;npe2=-1;npe3=-1;
! default value for test processor
ipetest=0

return
end

!==============================================================================
subroutine mpifinalize

include 'vacdef.f'

call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
call MPI_FINALIZE(ierrmpi)

return
end

!==============================================================================
subroutine ipe2ipeD(qipe,qipe1,qipe2,qipe3)

! Convert serial processor index to directional processor indexes

include 'vacdef.f'

integer:: qipe1,qipe2,qipe3, qipe
!-----------------------------------------------------------------------------
qipe1 = qipe - npe1*(qipe/npe1)
qipe2 = qipe/npe1 - npe2*(qipe/(npe1*npe2)) 
qipe3 = qipe/(npe1*npe2)                    

return
end

!==============================================================================
subroutine ipeD2ipe(qipe1,qipe2,qipe3,qipe)

! Convert directional processor indexes to serial processor index

include 'vacdef.f'

integer:: qipe1,qipe2,qipe3, qipe
!-----------------------------------------------------------------------------
qipe = qipe1  + npe1*qipe2  + npe1*npe2*qipe3

return
end

!==============================================================================
subroutine mpisetnpeDipeD(name)

! Set directional processor numbers and indexes based on a filename.
! The filename contains _np followed by np^D written with 2 digit integers.
! For example _np0203 means np1=2, np2=3 for 2D.

include 'vacdef.f'
character*79 :: name, nametail
integer:: i,qnpe1,qnpe2,qnpe3
logical:: npeDknown,npeDinname
!-----------------------------------------------------------------------------

oktest = index(teststr,'mpisetnpeDipeD')>0
if(oktest)write(*,*)'mpisetnpeDipeD ipe,name=',ipe,name

! Check if npe^D is already known
npeDknown  = npe1>0
if(npedknown .and. npe1*npe2*npe3 /= npe)then
   write(*,*)'npe=',npe,' /= product of npe^D=',npe1,npe2,npe3
   call mpistop('ERROR in setnpeDipeD')
endif

! Check if npe^D is given in the name
i=index(name,'_np')+3
npeDinname = i>3

if(.not.(npeDknown.or.npeDinname))call mpistop( &
   'ERROR in setnpeDipeD: npeD is neither known nor given in name='//name)

if(npeDinname)then
   ! read npe^D from name
   read(name(i:i+5),'(3i2)') qnpe1,qnpe2,qnpe3
   i=i+2*3
   nametail=name(i:79)
endif

if( npeDknown .and. npeDinname )then
   ! Check agreement
   if( qnpe1/=npe1.or.qnpe2/=npe2.or.qnpe3/=npe3 )then
      write(*,*)'npe^D=',npe1,npe2,npe3,' /= qnpe^D=',qnpe1,qnpe2,qnpe3,&
         ' read from filename=',name
      call mpistop('ERROR in mpisetnpeDipeD')
   endif
endif

if(npeDinname .and. .not.npeDknown)then
   ! set npe^D based on name
   npe1=qnpe1;npe2=qnpe2;npe3=qnpe3;
   if( npe1*npe2*npe3 /= npe)then
      write(*,*)'npe=',npe,' /= product of npe^D=',npe1,npe2,npe3,&
         ' read from filename=',name
      call mpistop('ERROR in setnpeDipeD')
   endif
endif

! Get directional processor indexes
call ipe2ipeD(ipe,ipe1,ipe2,ipe3)

if(npeDknown .and. .not.npeDinname)then
   ! insert npe^D into name
   i=index(name,'.')
   nametail=name(i:79)
   write(name(i:79),"('_np',3i2.2)") npe1,npe2,npe3
   i = i+3+2*3
endif

! insert ipe number into the filename
write(name(i:79),"('_',i3.3,a)") ipe,nametail(1:79-i-4)

! Set logicals about MPI boundaries for this processor

mpiupperB(1)=ipe1<npe1-1
mpilowerB(1)=ipe1>0 

mpiupperB(2)=ipe2<npe2-1
mpilowerB(2)=ipe2>0 

mpiupperB(3)=ipe3<npe3-1
mpilowerB(3)=ipe3>0 

if(oktest)write(*,*)'mpisetnpeDipeD: ipe,npeD,ipeD,name=',ipe,npe1,npe2,npe3,&
   ipe1,ipe2,ipe3,name

return
end

!==============================================================================
subroutine mpineighbors(idir,hpe,jpe)

! Find the hpe and jpe processors on the left and right side of this processor 
! in direction idir. The processor cube is taken to be periodic in every
! direction.

include 'vacdef.f'

integer :: idir,hpe,jpe,hpe1,hpe2,hpe3,jpe1,jpe2,jpe3
!-----------------------------------------------------------------------------
hpe1=ipe1-kr(1,idir);hpe2=ipe2-kr(2,idir);hpe3=ipe3-kr(3,idir);
jpe1=ipe1+kr(1,idir);jpe2=ipe2+kr(2,idir);jpe3=ipe3+kr(3,idir);

if(hpe1<0)hpe1=npe1-1
if(jpe1>=npe1)jpe1=0

if(hpe2<0)hpe2=npe2-1
if(jpe2>=npe2)jpe2=0

if(hpe3<0)hpe3=npe3-1
if(jpe3>=npe3)jpe3=0

call ipeD2ipe(hpe1,hpe2,hpe3,hpe)
call ipeD2ipe(jpe1,jpe2,jpe3,jpe)

return
end
!==============================================================================
subroutine mpigridsetup

! Distribute a grid of size nxall^D onto PE-s arranged in a cube of size npe^D

include 'vacdef.f'
!-----------------------------------------------------------------------------
!!!write(*,*)'nxall,npe=',nxall^D,npe^D

! Grid sizes on the processors (except for the last ones in some direction)
! This formula optimizes the load balance when nx^D is not a multiple of npe^D
nxpe1=(nxall1-1)/npe1+1;nxpe2=(nxall2-1)/npe2+1;nxpe3=(nxall3-1)/npe3+1; 

! Global grid indexes of the first grid point stored on this PE
ixPEmin1=ipe1*nxpe1+1;ixPEmin2=ipe2*nxpe2+1;ixPEmin3=ipe3*nxpe3+1;

! The last processors in a direction may have smaller grid sizes than nxpe
 
if(ipe1 < npe1-1)then
   nx1 = nxpe1
else
   nx1 = nxall1 - ixpemin1 + 1
endif

  
if(ipe2 < npe2-1)then
   nx2 = nxpe2
else
   nx2 = nxall2 - ixpemin2 + 1
endif

  
if(ipe3 < npe3-1)then
   nx3 = nxpe3
else
   nx3 = nxall3 - ixpemin3 + 1
endif


! Global grid index of the last grid point stored on this PE
ixPEmax1=ixPEmin1+nx1-1;ixPEmax2=ixPEmin2+nx2-1;ixPEmax3=ixPEmin3+nx3-1;

return
end

!=============================================================================
subroutine mpireduce(a,mpifunc)

! reduce input for one PE 0 using mpifunc

include 'mpif.h'

double precision :: a, alocal
integer          :: mpifunc, ierrmpi
!----------------------------------------------------------------------------
alocal = a
call MPI_REDUCE(alocal,a,1,MPI_DOUBLE_PRECISION,mpifunc,0,MPI_COMM_WORLD,&
   ierrmpi)

return
end

!==============================================================================
subroutine mpiallreduce(a,mpifunc)

! reduce input onto all PE-s using mpifunc

include 'mpif.h'

double precision :: a, alocal
integer          :: mpifunc, ierrmpi
!-----------------------------------------------------------------------------
alocal = a
call MPI_ALLREDUCE(alocal,a,1,MPI_DOUBLE_PRECISION,mpifunc,MPI_COMM_WORLD,&
   ierrmpi)

return
end

!==============================================================================
subroutine mpiix(ix1,ix2,ix3,jpe)

! Convert ix^D physical cell index on the global grid to local indexes 
! and set the processor number jpe to the processor that contains the cell

include 'vacdef.f'
integer :: ix1,ix2,ix3, jpe, jpe1,jpe2,jpe3
!-----------------------------------------------------------------------------

! Directional processor indexes
jpe1=(ix1-ixMmin1)/nxpe1;jpe2=(ix2-ixMmin2)/nxpe2;jpe3=(ix3-ixMmin3)/nxpe3;

! Conversion to local index
ix1=ix1-jpe1*nxpe1;ix2=ix2-jpe2*nxpe2;ix3=ix3-jpe3*nxpe3;

! Get MPI processor index
call ipeD2ipe(jpe1,jpe2,jpe3,jpe)

return
end

!==============================================================================
subroutine mpiixlimits(ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3)

! Convert global index limits to local index limits for this PE

include 'vacdef.f'
integer :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
!-----------------------------------------------------------------------------

if(ixmin1 > ixPEmax1)then
   ixmin1 = nx1
   ixmax1 = nx1-1
elseif(ixmax1 < ixPEmin1)then
   ixmax1 = 0
   ixmin1 = 1
else
   ixmin1 = max(ixmin1,ixPEmin1) - ixPEmin1 + 1
   ixmax1 = min(ixmax1,ixPEmax1) - ixPEmin1 + 1
endif


if(ixmin2 > ixPEmax2)then
   ixmin2 = nx2
   ixmax2 = nx2-1
elseif(ixmax2 < ixPEmin2)then
   ixmax2 = 0
   ixmin2 = 1
else
   ixmin2 = max(ixmin2,ixPEmin2) - ixPEmin2 + 1
   ixmax2 = min(ixmax2,ixPEmax2) - ixPEmin2 + 1
endif


if(ixmin3 > ixPEmax3)then
   ixmin3 = nx3
   ixmax3 = nx3-1
elseif(ixmax3 < ixPEmin3)then
   ixmax3 = 0
   ixmin3 = 1
else
   ixmin3 = max(ixmin3,ixPEmin3) - ixPEmin3 + 1
   ixmax3 = min(ixmax3,ixPEmax3) - ixPEmin3 + 1
endif


return
end
!==============================================================================

subroutine mpistop(message)

! Stop MPI run in an orderly fashion

include 'vacdef.f'

character(*) :: message
integer :: nerrmpi

!------------------------------------------------------------------------------
write(*,*)'ERROR for processor',ipe,':'
write(*,*)message
call MPI_abort(MPI_COMM_WORLD, nerrmpi, ierrmpi)

stop
end

!==============================================================================
subroutine mpibound(nvar,var)

! Fill in ghost cells of var(ixG,nvar) from other processors

include 'vacdef.f'

integer :: nvar
double precision :: var(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nvar)

! processor indexes for left and right neighbors
integer :: hpe,jpe
! index limits for the left and right side mesh and ghost cells 
integer :: ixLMmin1,ixLMmin2,ixLMmin3,ixLMmax1,ixLMmax2,ixLMmax3, ixRMmin1,&
   ixRMmin2,ixRMmin3,ixRMmax1,ixRMmax2,ixRMmax3, ixLGmin1,ixLGmin2,ixLGmin3,&
   ixLGmax1,ixLGmax2,ixLGmax3, ixRGmin1,ixRGmin2,ixRGmin3,ixRGmax1,ixRGmax2,&
   ixRGmax3
logical :: periodic

! There can be at most 2 receives in any direction for each PE
integer :: nmpirequest, mpirequests(2)
integer :: mpistatus(MPI_STATUS_SIZE,2)
common /mpirecv/ nmpirequest,mpirequests,mpistatus
!-----------------------------------------------------------------------------
oktest=index(teststr,'mpibound')>0
if(oktest)write(*,*)'mpibound ipe,nvar,varold=',ipe,nvar,var(ixtest1,ixtest2,&
   ixtest3,min(nvar,iwtest))


if(npe1>1)then
   nmpirequest =0
   mpirequests(1:2) = MPI_REQUEST_NULL

   periodic=typeB(1,2*1)=='mpiperiod'

   ! Left and right side ghost cell regions (target)
   ixLGmin1=ixGmin1;ixLGmin2=ixGmin2;ixLGmin3=ixGmin3;ixLGmax1=ixGmax1
   ixLGmax2=ixGmax2;ixLGmax3=ixGmax3; ixLGmax1=ixMmin1-1;
   ixRGmin1=ixGmin1;ixRGmin2=ixGmin2;ixRGmin3=ixGmin3;ixRGmax1=ixGmax1
   ixRGmax2=ixGmax2;ixRGmax3=ixGmax3; ixRGmin1=ixMmax1+1;

   ! Left and right side mesh cell regions (source)
   ixLMmin1=ixGmin1;ixLMmin2=ixGmin2;ixLMmin3=ixGmin3;ixLMmax1=ixGmax1
   ixLMmax2=ixGmax2;ixLMmax3=ixGmax3; ixLMmin1=ixMmin1
   ixLMmax1=ixMmin1+dixBmin1-1;
   ixRMmin1=ixGmin1;ixRMmin2=ixGmin2;ixRMmin3=ixGmin3;ixRMmax1=ixGmax1
   ixRMmax2=ixGmax2;ixRMmax3=ixGmax3; ixRMmax1=ixMmax1
   ixRMmin1=ixMmax1-dixBmax1+1;

   ! Obtain left and right neighbor processors for this direction
   call mpineighbors(1,hpe,jpe)

   if(oktest)then
      write(*,*)'mpibound ipe,idir=',ipe,1
      write(*,*)'mpibound ipe,ixLG=',ipe,ixLGmin1,ixLGmin2,ixLGmin3,ixLGmax1,&
         ixLGmax2,ixLGmax3
      write(*,*)'mpibound ipe,ixRG=',ipe,ixRGmin1,ixRGmin2,ixRGmin3,ixRGmax1,&
         ixRGmax2,ixRGmax3
      write(*,*)'mpibound ipe,ixLM=',ipe,ixLMmin1,ixLMmin2,ixLMmin3,ixLMmax1,&
         ixLMmax2,ixLMmax3
      write(*,*)'mpibound ipe,ixRM=',ipe,ixRMmin1,ixRMmin2,ixRMmin3,ixRMmax1,&
         ixRMmax2,ixRMmax3
      write(*,*)'mpibound ipe,hpe,jpe=',ipe,hpe,jpe
   endif

   ! receive right (2) boundary from left neighbor hpe
   if(mpilowerB(1) .or. periodic)call mpirecvbuffer(nvar,ixRMmin1,ixRMmin2,&
      ixRMmin3,ixRMmax1,ixRMmax2,ixRMmax3,hpe,2)
   ! receive left (1) boundary from right neighbor jpe
   if(mpiupperB(1) .or. periodic)call mpirecvbuffer(nvar,ixLMmin1,ixLMmin2,&
      ixLMmin3,ixLMmax1,ixLMmax2,ixLMmax3,jpe,1)
   ! Wait for all receives to be posted
   call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
   ! Ready send left (1) boundary to left neighbor hpe
   if(mpilowerB(1) .or. periodic)call mpisend(nvar,var,ixLMmin1,ixLMmin2,&
      ixLMmin3,ixLMmax1,ixLMmax2,ixLMmax3,hpe,1)
   ! Ready send right (2) boundary to right neighbor
   if(mpiupperB(1) .or. periodic)call mpisend(nvar,var,ixRMmin1,ixRMmin2,&
      ixRMmin3,ixRMmax1,ixRMmax2,ixRMmax3,jpe,2)
   ! Wait for messages to arrive
   call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)
   ! Copy buffer received from right (2) physical cells into left ghost cells
   if(mpilowerB(1) .or. periodic)call mpibuffer2var(2,nvar,var,ixLGmin1,&
      ixLGmin2,ixLGmin3,ixLGmax1,ixLGmax2,ixLGmax3)
   ! Copy buffer received from left (1) physical cells into right ghost cells
   if(mpiupperB(1) .or. periodic)call mpibuffer2var(1,nvar,var,ixRGmin1,&
      ixRGmin2,ixRGmin3,ixRGmax1,ixRGmax2,ixRGmax3)
endif


if(npe2>1)then
   nmpirequest =0
   mpirequests(1:2) = MPI_REQUEST_NULL

   periodic=typeB(1,2*2)=='mpiperiod'

   ! Left and right side ghost cell regions (target)
   ixLGmin1=ixGmin1;ixLGmin2=ixGmin2;ixLGmin3=ixGmin3;ixLGmax1=ixGmax1
   ixLGmax2=ixGmax2;ixLGmax3=ixGmax3; ixLGmax2=ixMmin2-1;
   ixRGmin1=ixGmin1;ixRGmin2=ixGmin2;ixRGmin3=ixGmin3;ixRGmax1=ixGmax1
   ixRGmax2=ixGmax2;ixRGmax3=ixGmax3; ixRGmin2=ixMmax2+1;

   ! Left and right side mesh cell regions (source)
   ixLMmin1=ixGmin1;ixLMmin2=ixGmin2;ixLMmin3=ixGmin3;ixLMmax1=ixGmax1
   ixLMmax2=ixGmax2;ixLMmax3=ixGmax3; ixLMmin2=ixMmin2
   ixLMmax2=ixMmin2+dixBmin2-1;
   ixRMmin1=ixGmin1;ixRMmin2=ixGmin2;ixRMmin3=ixGmin3;ixRMmax1=ixGmax1
   ixRMmax2=ixGmax2;ixRMmax3=ixGmax3; ixRMmax2=ixMmax2
   ixRMmin2=ixMmax2-dixBmax2+1;

   ! Obtain left and right neighbor processors for this direction
   call mpineighbors(2,hpe,jpe)

   if(oktest)then
      write(*,*)'mpibound ipe,idir=',ipe,2
      write(*,*)'mpibound ipe,ixLG=',ipe,ixLGmin1,ixLGmin2,ixLGmin3,ixLGmax1,&
         ixLGmax2,ixLGmax3
      write(*,*)'mpibound ipe,ixRG=',ipe,ixRGmin1,ixRGmin2,ixRGmin3,ixRGmax1,&
         ixRGmax2,ixRGmax3
      write(*,*)'mpibound ipe,ixLM=',ipe,ixLMmin1,ixLMmin2,ixLMmin3,ixLMmax1,&
         ixLMmax2,ixLMmax3
      write(*,*)'mpibound ipe,ixRM=',ipe,ixRMmin1,ixRMmin2,ixRMmin3,ixRMmax1,&
         ixRMmax2,ixRMmax3
      write(*,*)'mpibound ipe,hpe,jpe=',ipe,hpe,jpe
   endif

   ! receive right (2) boundary from left neighbor hpe
   if(mpilowerB(2) .or. periodic)call mpirecvbuffer(nvar,ixRMmin1,ixRMmin2,&
      ixRMmin3,ixRMmax1,ixRMmax2,ixRMmax3,hpe,2)
   ! receive left (1) boundary from right neighbor jpe
   if(mpiupperB(2) .or. periodic)call mpirecvbuffer(nvar,ixLMmin1,ixLMmin2,&
      ixLMmin3,ixLMmax1,ixLMmax2,ixLMmax3,jpe,1)
   ! Wait for all receives to be posted
   call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
   ! Ready send left (1) boundary to left neighbor hpe
   if(mpilowerB(2) .or. periodic)call mpisend(nvar,var,ixLMmin1,ixLMmin2,&
      ixLMmin3,ixLMmax1,ixLMmax2,ixLMmax3,hpe,1)
   ! Ready send right (2) boundary to right neighbor
   if(mpiupperB(2) .or. periodic)call mpisend(nvar,var,ixRMmin1,ixRMmin2,&
      ixRMmin3,ixRMmax1,ixRMmax2,ixRMmax3,jpe,2)
   ! Wait for messages to arrive
   call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)
   ! Copy buffer received from right (2) physical cells into left ghost cells
   if(mpilowerB(2) .or. periodic)call mpibuffer2var(2,nvar,var,ixLGmin1,&
      ixLGmin2,ixLGmin3,ixLGmax1,ixLGmax2,ixLGmax3)
   ! Copy buffer received from left (1) physical cells into right ghost cells
   if(mpiupperB(2) .or. periodic)call mpibuffer2var(1,nvar,var,ixRGmin1,&
      ixRGmin2,ixRGmin3,ixRGmax1,ixRGmax2,ixRGmax3)
endif


if(npe3>1)then
   nmpirequest =0
   mpirequests(1:2) = MPI_REQUEST_NULL

   periodic=typeB(1,2*3)=='mpiperiod'

   ! Left and right side ghost cell regions (target)
   ixLGmin1=ixGmin1;ixLGmin2=ixGmin2;ixLGmin3=ixGmin3;ixLGmax1=ixGmax1
   ixLGmax2=ixGmax2;ixLGmax3=ixGmax3; ixLGmax3=ixMmin3-1;
   ixRGmin1=ixGmin1;ixRGmin2=ixGmin2;ixRGmin3=ixGmin3;ixRGmax1=ixGmax1
   ixRGmax2=ixGmax2;ixRGmax3=ixGmax3; ixRGmin3=ixMmax3+1;

   ! Left and right side mesh cell regions (source)
   ixLMmin1=ixGmin1;ixLMmin2=ixGmin2;ixLMmin3=ixGmin3;ixLMmax1=ixGmax1
   ixLMmax2=ixGmax2;ixLMmax3=ixGmax3; ixLMmin3=ixMmin3
   ixLMmax3=ixMmin3+dixBmin3-1;
   ixRMmin1=ixGmin1;ixRMmin2=ixGmin2;ixRMmin3=ixGmin3;ixRMmax1=ixGmax1
   ixRMmax2=ixGmax2;ixRMmax3=ixGmax3; ixRMmax3=ixMmax3
   ixRMmin3=ixMmax3-dixBmax3+1;

   ! Obtain left and right neighbor processors for this direction
   call mpineighbors(3,hpe,jpe)

   if(oktest)then
      write(*,*)'mpibound ipe,idir=',ipe,3
      write(*,*)'mpibound ipe,ixLG=',ipe,ixLGmin1,ixLGmin2,ixLGmin3,ixLGmax1,&
         ixLGmax2,ixLGmax3
      write(*,*)'mpibound ipe,ixRG=',ipe,ixRGmin1,ixRGmin2,ixRGmin3,ixRGmax1,&
         ixRGmax2,ixRGmax3
      write(*,*)'mpibound ipe,ixLM=',ipe,ixLMmin1,ixLMmin2,ixLMmin3,ixLMmax1,&
         ixLMmax2,ixLMmax3
      write(*,*)'mpibound ipe,ixRM=',ipe,ixRMmin1,ixRMmin2,ixRMmin3,ixRMmax1,&
         ixRMmax2,ixRMmax3
      write(*,*)'mpibound ipe,hpe,jpe=',ipe,hpe,jpe
   endif

   ! receive right (2) boundary from left neighbor hpe
   if(mpilowerB(3) .or. periodic)call mpirecvbuffer(nvar,ixRMmin1,ixRMmin2,&
      ixRMmin3,ixRMmax1,ixRMmax2,ixRMmax3,hpe,2)
   ! receive left (1) boundary from right neighbor jpe
   if(mpiupperB(3) .or. periodic)call mpirecvbuffer(nvar,ixLMmin1,ixLMmin2,&
      ixLMmin3,ixLMmax1,ixLMmax2,ixLMmax3,jpe,1)
   ! Wait for all receives to be posted
   call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
   ! Ready send left (1) boundary to left neighbor hpe
   if(mpilowerB(3) .or. periodic)call mpisend(nvar,var,ixLMmin1,ixLMmin2,&
      ixLMmin3,ixLMmax1,ixLMmax2,ixLMmax3,hpe,1)
   ! Ready send right (2) boundary to right neighbor
   if(mpiupperB(3) .or. periodic)call mpisend(nvar,var,ixRMmin1,ixRMmin2,&
      ixRMmin3,ixRMmax1,ixRMmax2,ixRMmax3,jpe,2)
   ! Wait for messages to arrive
   call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)
   ! Copy buffer received from right (2) physical cells into left ghost cells
   if(mpilowerB(3) .or. periodic)call mpibuffer2var(2,nvar,var,ixLGmin1,&
      ixLGmin2,ixLGmin3,ixLGmax1,ixLGmax2,ixLGmax3)
   ! Copy buffer received from left (1) physical cells into right ghost cells
   if(mpiupperB(3) .or. periodic)call mpibuffer2var(1,nvar,var,ixRGmin1,&
      ixRGmin2,ixRGmin3,ixRGmax1,ixRGmax2,ixRGmax3)
endif


if(oktest)write(*,*)'mpibound ipe,varnew=',ipe,var(ixtest1,ixtest2,ixtest3,&
   min(nvar,iwtest))

return
end

!==============================================================================
subroutine mpisend(nvar,var,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,qipe,&
   iside)

! Send var(ix^L,1:nvar) to processor qipe.
! jside is 0 for min and 1 for max side of the grid for the sending PE

include 'vacdef.f'

integer :: nvar
double precision :: var(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nvar)
integer :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, qipe, iside, n, ix1,ix2,&
   ix3, ivar
!----------------------------------------------------------------------------
oktest = index(teststr,'mpisend')>0

n=0
do ivar=1,nvar
   do ix3=ixmin3,ixmax3;do ix2=ixmin2,ixmax2;do ix1=ixmin1,ixmax1;
      n=n+1
      sendbuffer(n)=var(ix1,ix2,ix3,ivar)
   enddo
   enddo
   enddo
end do

if(oktest)then
   write(*,*)'mpisend ipe-->qipe,iside,itag',ipe,qipe,iside,10*ipe+iside
   write(*,*)'mpisend ipe,ix^L,var=',ipe,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
      ixmax3,var(ixtest1,ixtest2,ixtest3,min(iwtest,nvar))
endif

call MPI_RSEND(sendbuffer(1),n,MPI_DOUBLE_PRECISION,qipe,10*ipe&
   +iside,MPI_COMM_WORLD,ierrmpi)

return
end

!==============================================================================
subroutine mpirecvbuffer(nvar,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,qipe,&
   iside)

! receive buffer for a ghost cell region of size ix^L sent from processor qipe
! and sent from side iside of the grid

include 'vacdef.f'

integer:: nvar, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, qipe, iside, n

integer :: nmpirequest, mpirequests(2)
integer :: mpistatus(MPI_STATUS_SIZE,2)
common /mpirecv/ nmpirequest,mpirequests,mpistatus
!----------------------------------------------------------------------------

oktest = index(teststr,'mpirecv')>0

n = nvar* (ixmax1-ixmin1+1)*(ixmax2-ixmin2+1)*(ixmax3-ixmin3+1)

if(oktest)write(*,*)'mpirecv ipe<--qipe,iside,itag,n',ipe,qipe,iside,10*qipe&
   +iside,n

nmpirequest = nmpirequest + 1
call MPI_IRECV(recvbuffer(1,iside),n,MPI_DOUBLE_PRECISION,qipe,10*qipe&
   +iside,MPI_COMM_WORLD,mpirequests(nmpirequest),ierrmpi)

return
end

!==============================================================================
subroutine mpibuffer2var(iside,nvar,var,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
   ixmax3)

! Copy mpibuffer(:,iside) into var(ix^L,1:nvar)
include 'vacdef.f'

integer :: nvar
double precision:: var(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nvar)
integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,iside,n,ix1,ix2,ix3,ivar
!-----------------------------------------------------------------------------
oktest = index(teststr,'buffer2var')>0

n=0
do ivar=1,nvar
   do ix3=ixmin3,ixmax3;do ix2=ixmin2,ixmax2;do ix1=ixmin1,ixmax1;
      n=n+1
      var(ix1,ix2,ix3,ivar)=recvbuffer(n,iside)
   enddo
   enddo
   enddo
end do

if(oktest)write(*,*)'buffer2var: ipe,iside,ix^L,var',ipe,iside,ixmin1,ixmin2,&
   ixmin3,ixmax1,ixmax2,ixmax3,var(ixtest1,ixtest2,ixtest3,min(iwtest,nvar))

return
end
