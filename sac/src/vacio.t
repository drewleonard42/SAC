!##############################################################################
! module vacio

!=============================================================================
SUBROUTINE readparameters(w)

  ! This subroutine sets or reads all default parameters from par/DEFAULT,
  ! then it reads the par/PROBLEM parameter file through standard input,
  ! and the initial data from data/PROBLEM.ini as soon as the filename is read 
  ! from the parameter file.

  USE constants
  USE common_varibles

  DOUBLE PRECISION:: w(ixG^T,nw)

  CHARACTER(^LENTYPE):: typepred(nw),typefull(nw),typeimpl(nw),typefilter(nw)
  DOUBLE PRECISION:: muscleta
  INTEGER:: i,j,k,iw,idim,iB,ifile,isave
  LOGICAL:: implmrpc,globalixtest^IFMPI

  ! The use of NAMELIST is not recommended by Fortran 90. It could be replaced
  ! by some string manipulations, but that is difficult in Fortran 77/Adaptor.
  ! I use NAMELIST, since it is simple and convenient.

  ! SM: NAMELIST is a defined part of the F90+ specification.

  NAMELIST /testlist/   teststr,ixtest1,ixtest2,ixtest3,&
       iwtest,idimtest,ipetest^IFMPI
  NAMELIST /filelist/   filenameini,filename,varnames{,^IFMPI npe^D}, &
       typefileini,typefileout,typefilelog,&
       snapshotini,snapshotout,fullgridini,fullgridout,dixB^L
  NAMELIST /savelist/   tsave,itsave,dtsave,ditsave
  NAMELIST /stoplist/   itmax,tmax,tmaxexact,dtmin,residmin,residmax,t,it,&
       cputimemax
  NAMELIST /methodlist/ wnames,fileheadout,eqpar,&
       typeadvance,typefull,typepred,typeimpl,typefilter,&
       typelimiter,typeentropy,entropycoef,artcomp,&
       typelimited,typefct,typetvd,typeaxial,&
       typepoisson,typeconstrain,constraincoef,&
       useprimitive,muscleta,musclomega,&
       acmwidth,acmnolim,acmcoef,acmexpo,fourthorder,&
       implmrpc,&
       dimsplit,typedimsplit,sourcesplit,typesourcesplit,&
       sourceunsplit,&
       divbfix,divbwave,divbconstrain,angmomfix,compactres,&
       smallfix,smallp,smallpcoeff,smallrho,smallrhocoeff,&
       vacuumrho,nproc,procpar
  NAMELIST /boundlist/  nB,typeB,ixB^LIM,idimB,upperB,extraB,typeBscalar,ipairB
  NAMELIST /paramlist/  courantpar,dtpar,dtdiffpar,dtcantgrow,slowsteps,&
       implmrpcpar,&
       implpar,impldiffpar,implerror,implrelax,impldwlimit,&
       implrestart,implrestart2,impliter,impliternr,&
       typeimplinit,typeimpliter,typeimplmat,&
       implnewton,implconserv,implcentered,&
       implnewmat,implpred,impl3level,impljacfast,implsource
  !-----------------------------------------------------------------------------

!!! Set new scalars for sake of unaltered par/DEFAULT
  constraincoef=one
  cputimemax=bigdouble
!!!

  ! Set default values for arrays (except the ones read from the inifile)

  DO ifile=1,nfile
     DO isave=1,nsavehi
        tsave(isave,ifile)=bigdouble   ! t  of saves into the output files 
        itsave(isave,ifile)=biginteger ! it of saves into the output files 
     END DO
     dtsave(ifile)=bigdouble           ! time between saves
     ditsave(ifile)=biginteger         ! timesteps between saves
     isavet(ifile)=1                   ! index for saves by t
     isaveit(ifile)=1                  ! index for saves by it
  END DO

  DO iw=1,nw
     typepred(iw)='default'     ! Predictor scheme (will be adjusted later)
     typefull(iw)='tvdlf'       ! Full step scheme
     typeimpl(iw)='nul'         ! Implicit step scheme
     typefilter(iw)='nul'       ! Filter scheme
     typelimiter(iw)='minmod'   ! Limiter type for flow variables/characteristics
     typeentropy(iw)='nul'      ! Entropy fix type
     artcomp(iw)=.FALSE.        ! No artificial compression for Harten type TVD
  END DO

  DO iw=1,nw
     acmcoef(iw)=-one       ! Coefficients (0,1) for the dissipative fluxes
  ENDDO                     ! negative value means no coefficient is used

  DO i=1,nfile+2            ! Elements define processing for the fullstep,
     nproc(i)=0             ! halfstep and the nfile output files. If the value
  END DO                    ! is 0, no processing. For nproc(1) and nproc(2)
  ! the value defines the proc. frequency. Negative
  ! value results in a call at every sweep. Positive
  ! value N results in a call at every N-th step before
  ! the first sweep. For nproc(ifile+2) the nonzero 
  ! values cause processing for that file.
  DO i=1,nprocpar
     procpar(i)=-one        ! Parameters for processing
  END DO
  smallp=-one               ! Default small pressure, redefined in getpthermal
  smallrho=-one             ! Default small density, redefined in keeppositive
  vacuumrho=-one            ! Density representing vacuum

  nB=2*ndim                 ! If nB is not specified by the user, gridsetup 
  DO iB=1,nhiB              ! will create 2*ndim boundary regions by default.
     DO iw=1,nw             
        typeB(iw,iB) ='cont'  ! Default boundary type
        fixedB(iw,iB)=.FALSE. ! Fixed boundaries are not extrapolated into yet
     END DO
     ipairB(iB)=0           ! periodic pair is unknown, but can be set or guessed
  END DO
  DO iw=1,nw
     DO idim=1,ndim
        nofluxB(iw,idim)=.FALSE.  ! No zero flux condition for variables
     ENDDO
  END DO
  dixB^L=2;                 ! Default width of boundary regions
  ixBmax(1,1)=0             ! An impossible value if user specifies boundaries.

  ! Read scalar parameters from the par/DEFAULT file

  unitpar=unitini-1
  OPEN(unitpar,file='par/DEFAULT',status='old')

  READ(unitpar,testlist)
  READ(unitpar,filelist)
  READ(unitpar,savelist)
  READ(unitpar,stoplist)
  READ(unitpar,methodlist)
  READ(unitpar,boundlist)
  READ(unitpar,paramlist)

  CLOSE(unitpar)

  ! end defaults

  ! Initialize Kronecker delta, and Levi-Civita tensor
  DO i=1,3
     DO j=1,3
        IF(i==j)THEN
           kr(i,j)=1
        ELSE
           kr(i,j)=0
        ENDIF
        DO k=1,3
           IF(i==j.OR.j==k.OR.k==i)THEN
              lvc(i,j,k)=0
           ELSE IF(i+1==j.OR.i-2==j)THEN
              lvc(i,j,k)=1
           ELSE
              lvc(i,j,k)=-1
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  ! Initialize error conunters and equation parameters
  DO i=1,nerrcode
     nerror(i)=0
  END DO
  DO i=1,neqpar+nspecialpar
     eqpar(i)=zero
  END DO

  ! read from STDIN
  unitpar=unitstdin
  {^IFMPI
  ! MPI reads from a file
  unitpar=unitini-1
  OPEN(unitpar,file='vac.par',status='old')
  }

  ! Start reading parameters from standard input, i.e. "< par/PROBLEM"
  {ipetest = -1 ^IFMPI}
  READ(unitpar,testlist)
  {^IFMPI
  ! ixtest^D is given for the full grid if ipetest was not set explicitly
  globalixtest = ipetest < 0
  ipetest      = MAX(ipetest,0)
  ! Erase test string for other processors unless ipetest >= npe (test all PEs)
  IF(ipetest<npe.AND.ipetest/=ipe)teststr=''
  }
  oktest=INDEX(teststr,'readparameters')>=1
  IF(oktest) WRITE(unitterm,*)'ReadParameters'
  IF(oktest) WRITE(unitterm,testlist)

  varnames='default'
  READ(unitpar,filelist)
  {^IFMPI
  ! Extract and check the directional processor numbers and indexes
  ! and concat the PE number to the input and output filenames
  CALL mpisetnpeDipeD(filenameini)
  CALL mpisetnpeDipeD(filename(fileout_))
  }

  IF(oktest) THEN 
     {^IFMPI WRITE(unitterm,*)'npe^D=',npe^D}
     WRITE(unitterm,*)filenameini
     DO ifile=1,nfile
        WRITE(unitterm,*)filename(ifile)
     ENDDO
     WRITE(unitterm,*)'Type of ini/out and log files:',&
          typefileini,typefileout,typefilelog
     IF(varnames/='default')WRITE(unitterm,*)'Varnames:',varnames
     IF(snapshotini>0)WRITE(unitterm,*)'Snapshotini:',snapshotini
     IF(snapshotout>0)WRITE(unitterm,*)'Snapshotout:',snapshotout
     WRITE(unitterm,*)'Fullgridini,out:',fullgridini,fullgridout
  ENDIF

  CALL readfileini(w)

  ! Default for output header line
  fileheadout=fileheadini

  {^IFMPI
  ! Reset global test cell indexes to local ones
  IF(globalixtest)CALL mpiix(ixtest^D,ipetest)
  }

  READ(unitpar,savelist)
  DO ifile=1,nfile
     IF(dtsave(ifile)<bigdouble/2.AND.oktest) &
          WRITE(unitterm,'(" DTSAVE  for file",i2," =",g10.5)') &
          ifile,dtsave(ifile)
     IF(ditsave(ifile)<biginteger.AND.oktest) &
          WRITE(unitterm,'(" DITSAVE for file",i2," =",i10)') &
          ifile,ditsave(ifile)
     IF(tsave(1,ifile)==bigdouble.AND.itsave(1,ifile)==biginteger.AND. &
          dtsave(ifile)==bigdouble.AND.ditsave(ifile)==biginteger.AND.verbose) &
          WRITE(uniterr,*)'Warning in ReadParameters: ', &
          'No save condition for file ',ifile
  ENDDO

  READ(unitpar,stoplist)
  IF(oktest)THEN
     IF(tmax<bigdouble)         WRITE(unitterm,*) 'TMAX= ',tmax
     IF(tmaxexact.AND.oktest)   WRITE(unitterm,*) 'TMAXEXACT=',tmaxexact
     IF(itmax<biginteger)       WRITE(unitterm,*) 'ITMAX=',itmax
     IF(dtmin>smalldouble)      WRITE(unitterm,*) 'DTMIN=',dtmin
     IF(residmin>smalldouble)   WRITE(unitterm,*) 'RESIDMIN=',residmin
     IF(residmax<bigdouble)     WRITE(unitterm,*) 'RESIDMAX=',residmax
     IF(cputimemax<bigdouble) WRITE(unitterm,*) 'CPUTIMEMAX=',cputimemax
  ENDIF

  IF(tmax==bigdouble.AND.itmax==biginteger) WRITE(uniterr,*) &
       'Warning in ReadParameters: Neither tmax nor itmax are given!'

  READ(unitpar,methodlist)


  typefull1='nul'
  typepred1='nul'
  typeimpl1='nul'
  typefilter1='nul'
  iw_full(niw_)=0
  iw_impl(niw_)=0
  iw_semi(niw_)=0
  iw_filter(niw_)=0
  DO iw=1,nw
     IF(typefull(iw)/='nul')THEN
        typefull1=typefull(iw)
        iw_full(niw_)=iw_full(niw_)+1
        iw_full(iw_full(niw_))=iw
     END IF
     IF(typepred(iw)/='nul')THEN
        typepred1=typepred(iw)
     ENDIF
     IF(typeimpl(iw)/='nul')THEN
        typeimpl1=typeimpl(iw)
        iw_impl(niw_)=iw_impl(niw_)+1
        iw_impl(iw_impl(niw_))=iw
     END IF
     IF(typefull(iw)/='nul'.AND.typeimpl(iw)=='nul')THEN
        iw_semi(niw_)=iw_semi(niw_)+1
        iw_semi(iw_semi(niw_))=iw
     END IF
     IF(typefilter(iw)/='nul')THEN
        typefilter1=typefilter(iw)
        iw_filter(niw_)=iw_filter(niw_)+1
        iw_filter(iw_filter(niw_))=iw
     END IF
  END DO

  IF(typefull1=='source'.AND.dimsplit.AND.sourceunsplit) &
       WRITE(*,*)'Warning: dimensional splitting for source terms!'

  IF(typepred1=='default')THEN
     SELECT CASE(typefull1)
     CASE('fct')
        typepred1='fct'
     CASE('tvdlf','tvdmu')
        typepred1='hancock'
     CASE('source')
        typepred1='source'
     CASE('tvdlf1','tvdmu1','tvd','tvd1','tvdmc','cd','cd4','mc','nul')
        typepred1='nul'
     CASE default
        CALL die('No default predictor for full step='//typefull1)
     END SELECT
  ENDIF

  muscleta1=(one-muscleta)/4
  muscleta2=(one+muscleta)/4

  polargrid=.FALSE.
  IF(typeaxial=='cylinder')THEN
     IF(phi_>ndir)CALL die(&
          'Error: typeaxial=cylinder with phi>ndir is impossible!')
     IF(phi_<=ndim.AND.phi_>=2)THEN
        polargrid=.TRUE.
        WRITE(*,*) 'Using polar coordinates...'
        IF(ndir==3.AND.(z_<2.OR.z_>ndir.OR.z_==phi_))CALL die( &
             'z direction is not set correctly! Use setup -z=? and make vac')
        IF(gencoord)WRITE(*,*)'generalized coordinates in the r,z plane...'
     ENDIF
  ENDIF

  IF(angmomfix)THEN
     WRITE(*,*)'Angular momentum conservation is on (angmomfix=T) for iw=',mphi_
     IF(typeaxial/='cylinder') CALL die(&
          'angmomfix works in cylindrical symmetry only!'//&
          ' Set typeaxial in par file')
     IF(phi_<2.OR.phi_>ndir) CALL die(&
          'phi direction does not satisfy 1<phi<=ndir! Use setvac -phi!')
     IF(gencoord)THEN
        WRITE(*,*)'angmomfix=T in generalized coordinates...'
        IF(polargrid.AND.phi_/=ndir) CALL die(&
             'phi=ndir is required for angmomfix in gen. coordinates!')
     ENDIF
  ENDIF

  ! Artificial compression requires Harten type TVD
  DO iw=1,nw
     IF(artcomp(iw))typetvd='harten'
  ENDDO

  ! Harmonize the parameters for dimensional splitting and source splitting
  IF(typedimsplit   =='default'.AND.     dimsplit)   typedimsplit='xyyx'
  IF(typedimsplit   =='default'.AND..NOT.dimsplit)   typedimsplit='unsplit'
  IF(typesourcesplit=='default'.AND.     sourcesplit)typesourcesplit='sfs'
  IF(typesourcesplit=='default'.AND..NOT.sourcesplit)typesourcesplit='unsplit'
  dimsplit   = typedimsplit   /='unsplit'
  sourcesplit= typesourcesplit/='unsplit'
  IF(sourcesplit)sourceunsplit=.FALSE.

  {^NOONED
  IF(.NOT.divbconstrain.OR.b0_==0)typeconstrain='nul'
  IF(typeconstrain=='nul')divbconstrain=.FALSE.

  IF(typeconstrain/='nul')THEN
     divbfix=.FALSE.
     !???divbwave=.false.???
     nproc(1:nfile+2)=0
     IF(typeaxial/='slab')WRITE(*,*)&
          'constrainb=T keeps div B = 0 ',&
          'only approximately in cylindrical symm.'
  ENDIF

  IF(divbfix.AND.(typephys=='mhd'.OR.typephys=='mhdiso').AND.&
       (.NOT.sourcesplit).AND.(.NOT.sourceunsplit))&
       CALL die('divbfix=T requires unsplitsource=T or splitsource=T !')

  IF(typefilter1=='tvdmu'.OR.typefilter1=='tvdlf')typelimited='predictor'
  IF((.NOT.dimsplit.OR.typeimpl1=='tvdlf1'.OR.typeimpl1=='tvdmu1')&
       .AND.typelimited=='original')typelimited='previous'
  }

  READ(unitpar,boundlist)
  IF(nB>nhiB) CALL die('Error in ReadParameters: Too many boundary regions')
  IF(ixBmax(1,1)==0.AND.oktest)THEN
     DO iB=1,nB
        WRITE(unitterm,'(" TYPEB(",i2,")       = ",100a10)') &
             iB,(typeB(iw,iB),iw=1,nw)
     END DO
     WRITE(unitterm,*)' EXTRAB = ',extraB
  END IF
  IF(ixBmax(1,1)/=0.AND.oktest) WRITE(unitterm,boundlist)
  IF(dixB^L>dixB^LLIM|.OR.|.OR.) CALL die( &
       'Error in ReadParameters: adjust dixBlo,dixBhi')

  ! Initialize implicit parameters as needed
  implsource=sourceunsplit
  IF(typeimpl1=='nul')THEN
     ! Explicit time integration
     implpar=-one
     IF(implmrpc)WRITE(*,*)'Warning: implmrpc=T but typeimpl is not set!'
     typeimplinit = 'unused'
     typeimplmat  = 'unused'
  ELSE IF(implmrpc)THEN
     ! MR-PC time integration
     implpar=one
     typeimpliter='vac_mrpc'
     typeimplmat ='free'
     impliter=       1
     implrestart=    5
     implerror=      1.D-3
     IF(residmin>0)THEN
        ! MR-PC for steady state
        impl3level=     .FALSE.
        typeimplinit=   'nul'
     ELSE
        ! MR-PC for time accurate
        impl3level=     .TRUE.
        typeimplinit=   'explicit2'
     ENDIF
  ELSE IF(typeimpl1=='source')THEN
     ! Semi-implicit time integration for sources

     IF(.NOT.sourceunsplit)CALL die('Implicit sources must be unsplit!')

     IF(residmin>zero)THEN
        WRITE(*,*)'Warning: Semi-implicit sources for steady state!'
        implerror=residmin
     ELSE
        implerror=1.D-5
     ENDIF

     ! Explicit fluxes, implicit sources --> Trapezoidal scheme
     impl3level=.FALSE.
     implpar=half
     typeimplinit='explicit'

     ! Sources only --> Matrix free approach
     typeimplmat='free'

  ELSE
     ! Implicit time integration
     IF(residmin>zero)THEN
        !Steady state calculation --> Backward Euler scheme
        impl3level=.FALSE.
        implpar=one
        typeimplinit='nul'
        implerror=residmin
        IF(iw_semi(niw_)>0)WRITE(*,*)&
             'Warning: some variables are explicit for steady state!'
     ELSE
        ! Implicit fluxes (and sources) --> BDF2 scheme
        impl3level=.TRUE.
        ! First step is Backward Euler, but it could be trapezoidal !!!
        implpar=one
        typeimplinit='nul'
        implerror=1.D-3
     ENDIF
     {^IFONED 
     typeimpliter='tridiag'
     typeimplmat='with'
     }
     {^NOONED
     typeimpliter='vac_bicg'
     typeimplmat='prec'
     }
     IF((typeimpl1=='tvdlf1'.OR.typeimpl1=='cd').AND..NOT.sourceunsplit)THEN
        impljacfast=.TRUE.
     ENDIF
  ENDIF

  READ(unitpar,paramlist)
  {^IFMPI
  CLOSE(unitpar)
  }

  IF(typeimpl1/='nul')THEN
     ! Check and/or set parameters for implicit calculations

     IF(implpar<=zero)CALL die('Error: implpar<=0 for implicit calculation!')
     IF(implpar>one)CALL die('Error: implpar>1 for implicit calculation!')

     IF(implmrpc.AND.typeimpliter/='vac_mrpc')&
          CALL die('Error: implmrpc=T requires typeimpliter=vac_mrpc')

     IF(implpred)THEN
        implconserv=.FALSE.
        typeadvance='onestep'
        implpar=one
     ENDIF

     IF(impl3level)THEN
        implpred=.FALSE.
        implconserv=.FALSE.
     ENDIF

     IF(typeimpliter=='tridiag'.AND.ndim==1)typeimplmat='with'
     IF(typeimplmat=='prec'.AND.ndim==1)typeimpliter='tridiag'

     ! steady state warnings
     IF(residmin>zero)THEN
        IF(implpar/=one)WRITE(*,*)'Warning: implpar<1 for steady state!'
        IF(impl3level)WRITE(*,*)'Warning: 3 level for steady state!'
     ENDIF

     IF(impljacfast.AND.typeimpl1/='tvdlf1'.AND.typeimpl1/='cd')&
          CALL die('Error: impljacfast=T works with typeimpl=tvdlf1 or cd!')

     IF(typeimpliter=='vac_mrpc')THEN
        IF(implnewton)CALL die('Error: MRPC with Newton-Raphson iteration!')
        IF(residmin>zero.AND.typeimplinit/='nul')&
             WRITE(*,*)'Warning: MRPC for steady state',&
             ' should not use typeimplinit=',typeimplinit
     ENDIF

  ENDIF

  ! If all predictor steps are 'nul' then 'twostep' method reduces to 'onestep'
  IF(typeadvance=='twostep'.AND.typepred1=='nul')typeadvance='onestep'

  IF(typefull1=='tvd'.AND.typeadvance/='onestep')&
       CALL die('tvd method should only be used as a onestep method')

  IF(typefull1=='tvd' .AND. .NOT.dimsplit)&
       WRITE(*,*)'Warning: One step TVD without dimensional splitting !!!'

  SELECT CASE(typeadvance)
  CASE('onestep','adams2')
     nstep=1
  CASE('twostep')
     nstep=2
  CASE('threestep')
     nstep=3
  CASE('fourstep','sterck','jameson')
     nstep=4
  CASE default
     CALL die('Unknown typeadvance='//typeadvance)
  END SELECT

  IF(oktest) WRITE(unitterm,methodlist)

  IF(oktest) WRITE(unitterm,paramlist)

  CALL setheaderstrings

  RETURN 
END SUBROUTINE readparameters

!=============================================================================
SUBROUTINE readfileini(w)

  ! Reads from unitini named filenameini in typefilini format.
  !
  ! The file may contain more than one snapshots, in which case the last set is 
  ! read. The compatibility of initial data with internal parameters is checked.
  
  USE constants
  USE common_varibles

  DOUBLE PRECISION:: w(ixG^T,nw)

  LOGICAL:: fileexist
  CHARACTER(91):: fhead
  !-----------------------------------------------------------------------------

  IF(typefileini=='auto')THEN
     INQUIRE(FILE=filenameini,EXIST=fileexist)
     IF(.NOT.fileexist) CALL die('Stop: file does not exist, filenameini='//&
          filenameini)
     OPEN(unitini,FILE=filenameini,STATUS='old')
     READ(unitini,'(a91)')fhead
     CLOSE(unitini)

     IF(ICHAR(fhead(1:1))/=0.AND.ICHAR(fhead(2:2))/=0.AND. &
          ICHAR(fhead(3:3))/=0.AND.ICHAR(fhead(4:4))/=0)THEN
        typefileini='ascii'
     ELSE IF(ICHAR(fhead(89:89))==0.AND.ICHAR(fhead(90:90))==0.AND. &
          (ICHAR(fhead(88:88))==24.OR.ICHAR(fhead(91:91))==24))THEN
        typefileini='binary'
     ELSE
        typefileini='special'
     ENDIF
     IF(verbose)THEN
        WRITE(*,*)'Auto typefileini=',typefileini
        IF(typefileini=='special'.AND. &
             ICHAR(fhead(89:89))==0.AND.ICHAR(fhead(90:90))==0.AND. &
             (ICHAR(fhead(88:88))==20.OR.ICHAR(fhead(91:91))==20)) &
             WRITE(*,*)'Looks like a real*4 file'
     ENDIF
  ENDIF
  IF(typefileout=='auto')THEN
     typefileout=typefileini
     IF(verbose)WRITE(*,*)'Auto typefileout=',typefileout
  ENDIF

  SELECT CASE(typefileini)
  CASE('ascii')
     CALL readfileini_asc(w)
  CASE('binary')
     CALL readfileini_bin(w)
  CASE default
     CALL die('Error in VAC: Unknown typefileini='//typefileini)
  END SELECT

  RETURN
END SUBROUTINE readfileini

!=============================================================================
SUBROUTINE readfileini_asc(w)

  ! Reads from unitini, filenameini in ASCII format.
  !
  ! The file may contain more than one snapshots, in which case the last set is 
  ! read. The compatibility of the initial data with internal parameters checked.
  !
  ! Variables in the order they are read from the file:
  !
  !   fileheadini - a header identifying the input file
  !   it,t        - the initial timestep and time
  !   ndimini     - dimensionality of grid,   Test: ==ndim
  !   neqparini   - number of eq. parameters, Test: <=neqpar+nspecialpar
  !   nwini       - number of flow variables, Test: ==nw
  !   nx          - the grid dimensions,      Test: <=ixGhi-dixBmax-ixMmin+1
  !   eqpar       - equation parameters from filenameini
  !   varnamesini - names of the coordinates, variables, equation parameters
  !                 eg. 'x y rho mx my e bx by  gamma eta'
  !   x           - the (initial) coordinates
  !   w           - the initial flow variables

  USE constants
  USE common_varibles

  DOUBLE PRECISION:: w(ixG^T,nw)

  LOGICAL:: fileexist
  INTEGER:: ios                           ! 0 if not EOF, -1 if EOF, >0 if error
  INTEGER:: ndimini,neqparini,neqparin,nwini,nwin ! values describing input data
  INTEGER:: ix^L,ix^D,idim,iw,ieqpar,snapshot
  DOUBLE PRECISION:: eqparextra,wextra
  CHARACTER(^LENNAME) :: varnamesini
  !-----------------------------------------------------------------------------

  oktest=INDEX(teststr,'readfileini')>=1

  IF(oktest) WRITE(unitterm,*)'ReadFileIni'

  INQUIRE(file=filenameini,exist=fileexist)
  IF(.NOT.fileexist) CALL die('Stop: file does not exist, filenameini='//&
       filenameini)
  OPEN(unitini,file=filenameini,status='old')

  snapshot=0
  DO
     READ(unitini,'(a)',iostat=ios,END=100)fileheadini
     IF(ios<0)EXIT                ! Cycle until the last recorded state
     IF(oktest) WRITE(unitterm,*)'fileheadini=',fileheadini(1:30)
     READ(unitini,*,iostat=ios)it,t,ndimini,neqparini,nwini
     IF(oktest) WRITE(unitterm, &
          "('it=',i7,' t=',g10.3,' ndim=',i3,' neqpar=',i3,' nw=',i3)")&
          it,t,ndimini,neqparini,nwini
     gencoord= ndimini<0
     CALL checkNdimNeqparNw(ndimini,neqparini,nwini,neqparin,nwin)
     READ(unitini,*,iostat=ios)nx
     IF(oktest) WRITE(unitterm,"('nx =',3i4)")nx
     CALL setixGixMix(ix^L)
     READ(unitini,*,iostat=ios)(eqpar(ieqpar),ieqpar=1,neqparin),&
          (eqparextra,ieqpar=neqparin+1,neqparini)
     IF(oktest) WRITE(unitterm,*)eqpar
     READ(unitini,'(a)',iostat=ios)varnamesini
     IF(varnames=='default')varnames=varnamesini
     IF(oktest) WRITE(unitterm,*)varnames

     {DO ix^DB=ixmin^DB,ixmax^DB;}
     READ(unitini,*,iostat=ios)(x(ix^D,idim),idim=1,ndim),&
          (w(ix^D,iw),iw=1,nwin),(wextra,iw=nwin+1,nwini)
     {END DO^D&\}
     IF(ios/=0)THEN
        WRITE(uniterr,*)'Stop: iostat=',ios
        CALL die('Error in reading file')
     END IF
     snapshot=snapshot+1
     IF(snapshot==snapshotini)EXIT
  END DO

100 CONTINUE

  CLOSE(unitini)

  IF(oktest) WRITE(*,*)'x,w:',&
       x(ixtest^D,idimtest),w(ixtest^D,iwtest)
  IF(oktest) WRITE(*,*)'x,w:',&
       x(ixtest^D,idimtest),w(ixtest^D,1:nw)

  RETURN
END SUBROUTINE readfileini_asc

!=============================================================================
SUBROUTINE readfileini_bin(w)

  ! Reads from unitini,filenameini in binary format.
  !
  ! The file may contain more than one snapshots, in which case the last set is 
  ! read unless snapshotini is set. 
  ! The compatibility of initial data with internal parameters is checked.

  USE constants
  USE common_varibles

  DOUBLE PRECISION:: w(ixG^T,nw)

  LOGICAL:: fileexist
  INTEGER:: ios                           ! 0 if not EOF, -1 if EOF, >0 if error
  INTEGER:: ndimini,neqparini,neqparin,nwini,nwin ! values describing input data
  INTEGER:: ix^L,idim,iw,ieqpar,snapshot
  DOUBLE PRECISION:: eqparextra
  CHARACTER(^LENNAME) :: varnamesini
  !-----------------------------------------------------------------------------

  oktest=INDEX(teststr,'readfileini')>=1

  IF(oktest) WRITE(unitterm,*)'ReadFileIni'

  INQUIRE(file=filenameini,exist=fileexist)
  IF(.NOT.fileexist) CALL die('Stop: file does not exist, filenameini='//&
       filenameini)
  OPEN(unitini,file=filenameini,status='old',form='unformatted')

  snapshot=0
  DO
     ! Read filehead
     READ(unitini,iostat=ios) fileheadini !END=100

     IF(ios<0)EXIT                ! Cycle until the last recorded state
     IF(oktest) WRITE(unitterm,*)'fileheadini=',fileheadini(1:30)

     ! Read params
     READ(unitini,iostat=ios)it,t,ndimini,neqparini,nwini
     IF(oktest) WRITE(unitterm, &
          "('it=',i7,' t=',g10.3,' ndim=',i3,' neqpar=',i3,' nw=',i3)")&
          it,t,ndimini,neqparini,nwini
     gencoord= ndimini<0
     ! Validate parameters?
     CALL checkNdimNeqparNw(ndimini,neqparini,nwini,neqparin,nwin)

     ! Read nx
     READ(unitini,iostat=ios)nx
     IF(oktest) WRITE(unitterm,"('nx =',3i4)")nx
     ! This set's up the global indicies based on nx and also
     ! deals with the MPI indicies etc.
     CALL setixGixMix(ix^L) 

     ! Read eqpar
     READ(unitini,iostat=ios)(eqpar(ieqpar),ieqpar=1,neqparin),&
          (eqparextra,ieqpar=neqparin+1,neqparini)
     IF(oktest) WRITE(unitterm,*)eqpar

     ! Read varnamesini
     READ(unitini,iostat=ios)varnamesini
     IF(varnames=='default')varnames=varnamesini
     IF(oktest) WRITE(unitterm,*)varnames

     ! Read x array
     READ(unitini,iostat=ios)(x(ix^S,idim),idim=1,ndim)

     ! Read w array
     ! To conform savefileout_bin we use loop for iw
     DO iw=1,nwin
        READ(unitini,iostat=ios)w(ix^S,iw)
     END DO
     IF(ios/=0)THEN
        WRITE(uniterr,*)'Error in ReadFileIni: iostat=',ios
        CALL die('Error in reading file')
     END IF
     snapshot=snapshot+1
     IF(snapshot==snapshotini)EXIT
  END DO

!100 CONTINUE

  CLOSE(unitini)

  IF(oktest) WRITE(*,*)'x,w:',&
       x(ixtest^D,idimtest),w(ixtest^D,iwtest)
  IF(oktest) WRITE(*,*)'x,w:',&
       x(ixtest^D,idimtest),w(ixtest^D,1:nw)

  RETURN
END SUBROUTINE readfileini_bin

!=============================================================================
SUBROUTINE checkNdimNeqparNw(ndimini,neqparini,nwini,neqparin,nwin)

  USE constants
  USE common_varibles

  INTEGER:: ndimini,neqparini,nwini,neqparin,nwin
  !-----------------------------------------------------------------------------

  IF(ndim/=ABS(ndimini))THEN
     WRITE(*,*)'Error in ReadFileini: ndimini=',ndimini
     CALL die('Incompatible dimensionalities')
  ENDIF

  IF(neqpar+nspecialpar/=neqparini)WRITE(*,"(a,i3,a,i3)")&
       'Warning in ReadFileini: number of eq.params=',neqpar,&
       ' /= neqparini=',neqparini

  IF(nw/=nwini)WRITE(*,"(a,i3,a,i3)")&
       'Warning in ReadFileini: number of variables nw=',nw,&
       ' /= nwini=',nwini

  IF((neqpar+nspecialpar/=neqparini.OR.nw/=nwini).AND.varnames=='default')&
       CALL die('Define varnames (in &filelist for VAC, in 3rd line for VACINI)!')

  ! The number of equation parameters and variables to be read
  neqparin=MIN(neqparini,neqpar+nspecialpar)
  nwin=MIN(nwini,nw)

  RETURN
END SUBROUTINE checkNdimNeqparNw


!=============================================================================
SUBROUTINE setixGixMix(ix^L)

  USE constants
  USE common_varibles

  INTEGER:: ix^L,qnx^IFMPI
  !-----------------------------------------------------------------------------

  ixGmin^D=ixGlo^D;
  ixMmin^D=ixGmin^D+dixBmin^D;

  ! Shave off ghost cells from nx
  IF(fullgridini)THEN
     {^DLOOP
     IF(ipe^D==0)^IFMPI       nx(^D)=nx(^D)-dixBmin^D
     IF(ipe^D==npe^D-1)^IFMPI nx(^D)=nx(^D)-dixBmax^D
     \}
  ENDIF

  ! Calculate mesh and grid sizes
  ixMmax^D=ixMmin^D+nx(^D)-1;
  ixGmax^D=ixMmax^D+dixBmax^D;

  ! Set the index range for this grid
  IF(fullgridini)THEN
     ix^L=ixG^L;
     {^IFMPI
     ! Set index range to mesh value if the boundary is not an outer bounary
     ^D&IF(ipe^D>0)ixmin^D=ixMmin^D\
     ^D&IF(ipe^D<npe^D-1)ixmax^D=ixMmax^D\
     }
  ELSE
     ix^L=ixM^L;
  ENDIF

  IF(ixGmax^D>ixGhi^D|.OR.)THEN
     WRITE(uniterr,*)'Stop: nxhi=',ixGhi^D-dixBmax^D-ixMmin^D+1
     CALL die('Error in SetixGixMix')
  END IF

  nx^D=nx(^D);

  {^IFMPI
  ! set global grid size by adding up nx and 
  ! dividing by the number of processors in the orthogonal plane
  {^DLOOP
  CALL MPI_allreduce(nx^D,nxall^D,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
  nxall^D=nxall^D*npe^D/npe;
  \}

  ! Distribute global grid onto processor cube
  CALL mpigridsetup
  }

  RETURN
END SUBROUTINE setixGixMix

!=============================================================================
SUBROUTINE setheaderstrings

  ! Check and/or put physics and equation parameter names into file header

  USE constants
  USE common_varibles

  INTEGER:: i
  CHARACTER(^LENTYPE) :: physics
  !-----------------------------------------------------------------------------

  ! Check physics or add _typephysNDIMNDIR

  WRITE(physics,'(a,i1,i1)')typephys,^ND,^NC

  i=INDEX(fileheadout,'_')
  IF(i>=1)THEN
     IF(physics/=fileheadout(i+1:i+^LENTYPE))THEN
        WRITE(*,*)'This code is configured to ',physics
        CALL die('Error:  physics in file is '//fileheadout(i+1:i+^LENTYPE))
     ENDIF
  ELSE
     i=^LENNAME-LEN(typephys)-3
     DO
        IF(fileheadout(i:i)/=' ' .OR. i==1)EXIT
        i=i-1
     ENDDO
     fileheadout=fileheadout(1:i)//'_'//physics
     WRITE(*,*)'Warning: incomplete input headline.',&
          ' Added to output headline _',physics
  ENDIF

  ! Check for equation parameter names in varnames, add them if missing

  IF(varnames/='default' .AND. INDEX(varnames,eqparname)<=0)THEN
     i=^LENNAME-LEN(eqparname)-3
     DO
        IF(varnames(i:i)/=' ' .OR. i==1)EXIT
        i=i-1
     ENDDO
     varnames=varnames(1:i)//'   '//eqparname
  ENDIF

  ! Check for special equation parameter names in varnames, add them if missing

  IF(varnames/='default' .AND. INDEX(varnames,specialparname)<=0)THEN
     i=^LENNAME-LEN(specialparname)-3
     DO
        IF(varnames(i:i)/=' ' .OR. i==1)EXIT
        i=i-1
     ENDDO
     varnames=varnames(1:i)//'   '//specialparname
  ENDIF

  RETURN
END SUBROUTINE setheaderstrings

!=============================================================================
SUBROUTINE savefile(ifile,w)

  USE constants
  USE common_varibles

  INTEGER:: ifile,ix^L
  DOUBLE PRECISION:: w(ixG^T,nw)
  CHARACTER(10):: itstring
  !-----------------------------------------------------------------------------

  !if(nproc(ifile+2)>0) call process(ifile+2,1,ndim,w)

  ! In most cases the mesh should be saved
  ix^L=ixM^L;

  SELECT CASE(ifile)
  CASE(fileout_)
     ! Produce the output file name
     filenameout=filename(fileout_)
     IF(snapshotout>0.AND.isaveout>0)THEN
        IF(isaveout==snapshotout*(isaveout/snapshotout))THEN
           CLOSE(unitini+ifile)
           WRITE(itstring,'(i10)')isaveout+1
           filenameout=filenameout(1:INDEX(filenameout,' ')-1)//'_'// &
                itstring(10-INT(alog10(isaveout+1.5)):10) 
        ENDIF
     ENDIF
     isaveout=isaveout+1
     IF(fullgridout)THEN
        ix^L=ixG^L;
        {^IFMPI
        ! Set index range to mesh value if not an outer boundary
        ^D&IF(ipe^D>0)ixmin^D=ixMmin^D\
        ^D&IF(ipe^D<npe^D-1)ixmax^D=ixMmax^D\
        }
     END IF
     SELECT CASE(typefileout)
     CASE('ascii')
        CALL savefileout_asc(unitini+ifile,w,ix^L)
     CASE('binary')
        CALL savefileout_bin(unitini+ifile,w,ix^L)
     CASE default
        CALL die('Error in SaveFile: Unknown typefileout:'//typefileout)
     END SELECT
  CASE(filelog_)
     SELECT CASE(typefilelog)
     CASE('default')
        CALL savefilelog_default(unitini+ifile,w,ix^L)
     CASE default
        CALL die('Error in SaveFile: Unknown typefilelog:'//typefilelog)
     END SELECT
  CASE default
     WRITE(*,*) 'No save method is defined for ifile=',ifile
     CALL die(' ')
  END SELECT

  RETURN 
END SUBROUTINE savefile

!=============================================================================
SUBROUTINE savefileout_asc(qunit,w,ix^L)

  ! This version saves into filename(fileout_) ASCII data at every save time
  ! in full accordance with the ReadFileini subroutine, except that the first
  ! line is fileheadout and not fileheadini.

  USE constants
  USE common_varibles

  INTEGER:: qunit,ix^L,ix^D,iw,idim,ndimout
  DOUBLE PRECISION:: w(ixG^T,nw),qw(nw)
  LOGICAL:: fileopen
  !-----------------------------------------------------------------------------

  INQUIRE(qunit,opened=fileopen)

  IF(.NOT.fileopen)OPEN(qunit,file=filenameout,status='unknown')

  IF(gencoord)THEN
     ndimout= -ndim
  ELSE
     ndimout= ndim
  ENDIF

  WRITE(qunit,"(a)")fileheadout
  WRITE(qunit,"(i7,1pe13.5,3i3)")it,t,ndimout,neqpar+nspecialpar,nw
  WRITE(qunit,"(3i4)") ixmax^D-ixmin^D+1
  WRITE(qunit,"(100(1pe13.5))")eqpar
  WRITE(qunit,"(a)")varnames
  {DO ix^DB= ixmin^DB,ixmax^DB \}
  ! Values with magnitude less than smalldouble are written as 0d0
  WHERE(ABS(w(ix^D,1:nw))>5.0d-16)
     qw(1:nw)=w(ix^D,1:nw)
  ELSEWHERE
     qw(1:nw)=0d0
  END WHERE
  WRITE(qunit,"(100(1pe18.10))")x(ix^D,1:ndim),qw(1:nw)
ENDDO^D&;

CALL flushunit(qunit)

RETURN 
END SUBROUTINE savefileout_asc

!=============================================================================
SUBROUTINE savefileout_bin(qunit,w,ix^L)

  ! This version saves into filename(fileout_) binary data at every save time
  ! in full accordance with the ReadFileini subroutine, except that the first
  ! line is fileheadout and not fileheadini.

  USE constants
  USE common_varibles

  INTEGER:: qunit,ix^L,idim,iw,ndimout
  DOUBLE PRECISION:: w(ixG^T,nw)
  LOGICAL:: fileopen

  !**************** slice
  INTEGER:: s_ixmax^D, prom_ixmax^D, s_ixmin^D, prom_ixmin^D  
  !**************** endslice


  !-----------------------------------------------------------------------------

  INQUIRE(qunit,opened=fileopen)
  IF(.NOT.fileopen)&
       OPEN(qunit,file=filenameout,status='unknown',form='unformatted')

  IF(gencoord)THEN
     ndimout= -ndim
  ELSE
     ndimout= ndim
  ENDIF

  WRITE(qunit)fileheadout
  WRITE(qunit)it,t,ndimout,neqpar+nspecialpar,nw
  WRITE(qunit) ixmax^D-ixmin^D+1
  WRITE(qunit)eqpar
  WRITE(qunit)varnames
  WRITE(qunit)(x(ix^S,idim),idim=1,ndim)

  ! write(qunit)w(ix^S,1:nw) produces segmentation fault on Alpha, thus loop

  DO iw=1,nw
     WRITE(qunit)w(ix^S,iw)
  END DO

  CALL flushunit(qunit)


  !**************** slice *********************************

  !inquire(qunit,opened=fileopen)
  !if(.not.fileopen)&
  !   open(qunit,file=filenameout,status='unknown',form='unformatted')

  !if(gencoord)then
  !   ndimout= -ndim
  !else
  !   ndimout= ndim
  !endif

  !prom_ixmax^D=ixmax^D
  !prom_ixmin^D=ixmin^D
  !ixmax1=2
  !ixmin1=0

  !write(qunit)fileheadout
  !write(qunit)it,t,ndimout,neqpar+nspecialpar,nw
  !write(qunit) ixmax^D-ixmin^D+1
  !write(qunit)eqpar
  !write(qunit)varnames
  !write(qunit)(x(ix^S,idim),idim=1,ndim)

  ! write(qunit)w(ix^S,1:nw) produces segmentation fault on Alpha, thus loop
  !
  !do iw=1,nw
  !   write(qunit)w(ix^S,iw)
  !end do

  !call flushunit(qunit)


  !************* end slice ********************************** 

  RETURN 
END SUBROUTINE savefileout_bin

!=============================================================================
SUBROUTINE savefilelog_default(qunit,w,ix^L)

  ! This version saves into filename(filelog_) the following formatted data:
  !
  !   fileheadout
  !   STRING_DESCRIBING_COLUMNS
  !   it t dt wmean(1) wmean(2) ... wmean(nw) residual
  !   it t dt wmean(1) wmean(2) ... wmean(nw) residual
  !   it t dt wmean(1) wmean(2) ... wmean(nw) residual
  !   etc.
  !
  ! at every save time. wmean is the volume averaged w, residual is saved 
  ! if residmin>0 is set in the parfile.

  USE constants
  USE common_varibles

  INTEGER:: qunit,ix^L
  DOUBLE PRECISION:: w(ixG^T,nw)
  INTEGER:: iw
  LOGICAL:: fileopen
  DOUBLE PRECISION:: wmean(nw)
  !-----------------------------------------------------------------------------

  IF(ipe==0)THEN^IFMPI

     INQUIRE(qunit,opened=fileopen)
     IF(.NOT.fileopen)THEN
        OPEN(qunit,file=filename(filelog_),status='unknown')
        WRITE(qunit,'(a)')fileheadout
        IF(residmin>zero.OR.residmax<bigdouble)THEN
           WRITE(qunit,'(a15,a55,a9)')'it   t   dt    ',wnames,' residual'
        ELSE
           WRITE(qunit,'(a15,a64)')   'it   t   dt    ',wnames
        ENDIF
     ENDIF

  ENDIF^IFMPI

  DO iw=1,nw 
     wmean(iw)=SUM(dvolume(ix^S)*w(ix^S,iw))/volume
     {^IFMPI CALL mpireduce(wmean(iw),MPI_SUM)}
  END DO

  IF(ipe==0)THEN^IFMPI

     IF(residmin>zero.OR.residmax<bigdouble)THEN
        WRITE(qunit,'(i7,100(1pe13.5))')it,t,dt,wmean,residual
     ELSE
        WRITE(qunit,'(i7,100(1pe13.5))')it,t,dt,wmean
     ENDIF
     CALL flushunit(qunit)

  ENDIF^IFMPI

  RETURN 
END SUBROUTINE savefilelog_default

!=============================================================================
SUBROUTINE die(message)

  CHARACTER(*) :: message
  !-----------------------------------------------------------------------------
  {^IFMPI CALL mpistop(message)}
  WRITE(*,*)message

  STOP
END SUBROUTINE die
!=============================================================================
SUBROUTINE flushunit(qunit)

  !USE F90_UNIX_IO,ONLY : flush ! F90=f95 (NAG)
  IMPLICIT NONE

  INTEGER::qunit, ierror

  !CALL FLUSH(qunit) ! F90=f95 (NAG)
  !call flush(qunit)   ! OS=Linux, SunOS, UNICOS, T3E, Fujitsu
  !call flush_(qunit)  ! OS=AIX, F90=xlf

  ! no flush on Linux IA64 with Intel compiler

  RETURN
END SUBROUTINE flushunit
!=============================================================================
! end module vacio
!##############################################################################



