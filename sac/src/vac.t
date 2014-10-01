!##############################################################################
! module vac

!=============================================================================
PROGRAM vac

  ! SAC - Sheffield Advanced Code, 2010 S. Shelag, V Fedun and R. Erdelyi
  ! Based on the
  ! Versatile Advection Code, (c) Gabor Toth. Started on Nov 8, 1994.

  ! Pulled upto FORTRAN 2008 by Stuart Mumford 2013

  USE constants
  USE common_variables

  INTEGER:: ifile,ierrcode,iw
  REAL(kind=8):: w(ixG^T,nw),wnrm2,dtold,time0,time1

  ! functions
  LOGICAL:: timetofinish,timetosave
  REAL(kind=8):: cputime
  !-----------------------------------------------------------------------------
  {CALL mpiinit ^IFMPI}

  verbose=.TRUE. .AND.ipe==0^IFMPI
  IF(verbose)THEN
     WRITE(*,'(a)')'VAC 4.52 configured to'
     WRITE(*,'(a)')'  -d=22 -phi=0 -z=0 -g=100,100 -p=mhd -u=default'
     WRITE(*,'(a)')'  -on=cd,rk'
     WRITE(*,'(a)')'  -off=mc,fct,tvdlf,tvd,impl,poisson,ct,gencoord,resist,mpi'
     {^IFMPI WRITE(*,'(a,i3,a)')'Running on ',npe,' processors'}
  ENDIF

  IF(ipe==0)^IFMPI  time0=cputime()

  CALL physini                  ! Initialize physics dependent variables
  CALL readparameters(w)        ! Read filenames and parameters for advection
  ! Read initial data, set ixM,ixG,gencoord

  {^NOGEN 
  IF(gencoord)THEN
     WRITE(*,*) 'Error: input file contains general grid'
     WRITE(*,*) 'Recompile vac after setvac -on=gencoord is set.'
  ENDIF
  }
  CALL boundsetup               ! Initialize boundary data
  ! Initialize grid geometry
  IF(gencoord)THEN
     {^IFGEN CALL gridsetup2}
     {^NOGEN CALL die('Error: gencoord module is off')}
  ELSE
     CALL gridsetup1
  ENDIF

  CALL startup                  ! Initialize it, t, headers etc.

  IF(verbose)WRITE(*,'(a,f10.3,a)')'Start Advance  ',cputime()-time0,' sec'

  CALL getboundary(t,1,nw,1,ndim,w)
  DO
     DO ifile=1,nfile
        IF(timetosave(ifile)) CALL savefile(ifile,w)
     END DO

     ! Determine time step
     IF(dtpar>zero)THEN
        dt=dtpar
     ELSE
        IF(courantpar>zero)CALL getdt_courant(w,ixM^L)
        CALL getdt(w,ixM^L)
        CALL getdt_special(w,ixM^L)

        IF(dtcantgrow.AND.it>itmin)dt=MIN(dt,dtold)
        dtold=dt
     ENDIF
     IF(dtmrpc>zero)dt=MIN(dt,dtmrpc)

     IF (timetofinish(time0)) EXIT

     ! For slowsteps == 1, use dtpar in the first time step ONLY
     IF(slowsteps==1.AND.it==itmin)dtpar=-one

     ! For slowsteps > 1, reduce dt for the first few steps 
     IF(slowsteps>it-itmin+1) dt=dt*(one-(one-(it-itmin+one)/slowsteps)**2)

     IF(tmaxexact)dt=MIN(dt,tmax-t+smalldouble)

     ! Store w into wold for residual calculations and 
     ! for TVD limiting based on the previous time step.
     wold(ixG^S,1:nw)=w(ixG^S,1:nw)

     ! Advance w (except variables with typefull='nul') by dt in the full grid
     CALL advance(iw_full,w)

     IF(residmin>zero .OR. residmax<bigdouble)THEN
        ! calculate true residual ||w_n+1-w_n|| for steady-state calculations
        residual=zero
        DO iw=1,nw
           wnrm2=SUM(w(ixM^S,iw)**2)
           {^IFMPI CALL mpiallreduce(wnrm2,MPI_SUM)}
           IF(wnrm2<smalldouble)wnrm2=one
           residual = residual + SUM((w(ixM^S,iw)-wold(ixM^S,iw))**2)/wnrm2
        ENDDO
        {^IFMPI CALL mpiallreduce(residual,MPI_SUM)}
        residual=SQRT(residual/nw)
     ENDIF

     it=it+1
     t=t+dt
     PRINT*, "****t=", t
  END DO

  time1=cputime()-time0


  DO ifile=1,nfile
     IF(itsavelast(ifile)<it)CALL savefile(ifile,w)
     CLOSE(unitini+ifile)
  ENDDO

  IF(verbose)WRITE(*,'(a,f10.3,a)')'Finish Advance ',time1,' sec'

  IF(ipe==0)THEN^IFMPI

     IF(dt<dtmin)WRITE(unitterm,*)'Warning: dt<dtmin !!!'
     IF(time1>cputimemax)WRITE(unitterm,*)'Warning: cputimemax exceeded !!!'
     DO ierrcode=1,nerrcode
        IF(nerror(ierrcode)>0)THEN
           WRITE(*,"(a,i2,a,i5,a)")'Error (code=',&
                ierrcode,') occurred ',nerror(ierrcode),' times !!!'
           SELECT CASE(ierrcode)
           CASE(toosmallp_)
              WRITE(*,"(a)")'Error description: Pressure below psmall'
           CASE(toosmallr_)
              WRITE(*,"(a)")'Error description: Density below rhosmall'
           CASE(couranterr_)
              WRITE(*,"(a)")'Error description: Courant number above 1'
           CASE(poissonerr_)
              WRITE(*,"(a)")'Error description: Poisson solver failed'
           END SELECT
        ENDIF
     END DO
  ENDIF^IFMPI

  IF(verbose)THEN
     IF(implpar>zero)THEN
        WRITE(*,*)'Number of explicit evaluations:',nexpl
        WRITE(*,*)'Number of Newton iterations   :',nnewton
        WRITE(*,*)'Number of linear iterations   :',niter
        WRITE(*,*)'Number of MatVecs             :',nmatvec
     ENDIF

     IF(residmin>zero .OR. residmax<bigdouble)THEN
        WRITE(*,*)'Number of time steps          :',it-itmin
        WRITE(*,*)'Residual and residmin         :',residual,residmin
     ENDIF

     WRITE(*,'(a,f10.3,a)')'Finished VAC   ',cputime()-time0,' sec'
  ENDIF

  {CALL mpifinalize ^IFMPI}

END PROGRAM vac

!=============================================================================
SUBROUTINE startup

  USE constants
  USE common_variables

  INTEGER:: ifile,iw,ivector,idim,qnvector
  !-----------------------------------------------------------------------------

  ! Initialize dtmrpc which will be calculated by MRPC
  dtmrpc=-one

  ! Initialize dtcourant, which will be calculated by TVD, TVD-MUSCL or TVDLF
  DO idim=1,ndim
     dtcourant(idim)=bigdouble
  ENDDO

  ! If dtpar is set, and not only for the first time step (slowsteps/=1)
  ! then set courantpar<0, so that dtcourant is not calculated at all
  IF(dtpar>zero.AND.slowsteps/=1)courantpar= -one

  itmin=it
  DO ifile=1,nfile
     tsavelast(ifile)=t
     itsavelast(ifile)=it
  END DO
  isaveout=0

  ! Initialize vectoriw based on iw_vector. vectoriw=-1 for scalar variables,
  DO iw=1,nw
     vectoriw(iw)=-1
  END DO
  ! It points to the 0-th component (iw_vector=m0_,b0_,...) for vector variables.
  ! Only the first ndim components of the vector variables are rotated in
  ! generalized coordinates. 
  ! qnvector is only used to avoid compiler warning when nvector=0

  qnvector=nvector
  DO ivector=1,qnvector
     DO idim=1,ndim
        vectoriw(iw_vector(ivector)+idim)=iw_vector(ivector)
     END DO
  END DO

  {^IFPHI
  !For 3D polar (r,z,phi) grid m_phi behaves as a scalar
  !and we take advantage of that to conserve angular momentum
  IF(angmomfix.AND.polargrid.AND.gencoord)vectoriw(mphi_)=-1
  }

  ! Initial value for residual and counters
  residual=bigdouble
  nexpl=0
  nnewton=0
  nmatvec=0
  niter=0

  {^IFCT
  ! fstore should be zero for fluxCD, fluxCT, ryuCT schemes
  fstore(ixG^S,1:ndim)=zero
  }

  RETURN
END SUBROUTINE startup

!=============================================================================
SUBROUTINE advance(iws,w)

  ! w(iws,t) -> w(iws,t+qdt) based on typedimsplit and typesourcesplit
  !
  ! Add split sources and fluxes with unsplit sources

  USE constants
  USE common_variables

  INTEGER:: iws(niw_)
  REAL(kind=8):: w(ixG^T,nw), w1(ixG^T,nw)
  !-----------------------------------------------------------------------------

  ! Add split sources berforehand if this is required
  IF(sourcesplit)THEN
     w1(ixG^S,1:nw)=w(ixG^S,1:nw)
     SELECT CASE(typesourcesplit)
     CASE('sf')
        CALL addsource2(dt,ixG^L,ixM^L,iws,t,w1,t,w)
     CASE('sfs')
        CALL addsource2(dt/2,ixG^L,ixM^L,iws,t,w1,t,w)
     CASE('ssf')
        CALL addsource2(dt/2,ixG^L,ixG^L,iws,t,w,t,w1)
        CALL addsource2(dt,ixG^L,ixM^L,iws,t,w1,t,w)
     CASE('ssfss')
        CALL addsource2(dt/4,ixG^L,ixG^L,iws,t,w,t,w1)
        CALL addsource2(dt/2,ixG^L,ixM^L,iws,t,w1,t,w)
     CASE default
        CALL die('Error: Unknown typesourcesplit='//typesourcesplit)
     END SELECT
     CALL getboundary(t,1,nw,1,ndim,w)
  ENDIF

  ! Add fluxes and unsplit sources explicitly or implicitly
  IF(typeimpl1=='nul')THEN
     CALL advance_expl(typefull1,ixG^L,iws,w1,w)
  ELSE
     {^IFIMPL
     CALL advance_impl(ixG^L,w1,w)
     IF(.FALSE.)} CALL die('IMPL module is switched off')
  ENDIF

  ! Add split sources afterwards if this is required
  IF(sourcesplit)THEN
     SELECT CASE(typesourcesplit)
     CASE('sfs')
        w1(ixG^S,1:nw)=w(ixG^S,1:nw)
        CALL addsource2(dt/2,ixG^L,ixM^L,iws,t+dt,w1,t+dt,w)
        CALL getboundary(t+dt,1,nw,1,ndim,w)
     CASE('ssfss')
        w1(ixG^S,1:nw)=w(ixG^S,1:nw)
        CALL addsource2(dt/4,ixG^L,ixG^L,iws,t+dt,w ,t+dt,w1)
        CALL addsource2(dt/2,ixG^L,ixM^L,iws,t+dt,w1,t+dt, w)
        CALL getboundary(t+dt,1,nw,1,ndim,w)
     END SELECT
  ENDIF

  RETURN
END SUBROUTINE advance

!=============================================================================
SUBROUTINE advance_expl(method,ix^L,iws,w1,w)

  ! w(t) -> w(t+qdt) within ix^L based on typedimsplit, typesourcesplit, nproc
  !
  ! Add fluxes and unsplit sources, possibly with dimensional splitting
  ! Boundaries should be kept updated by addsource2 and advect
  !
  ! w1 can be ised freely.

  USE constants
  USE common_variables

  CHARACTER(^LENTYPE) :: method
  INTEGER :: ix^L,iws(niw_)
  REAL(kind=8) :: w(ixG^T,nw),w1(ixG^T,nw)

  LOGICAL :: firstsweep,lastsweep
  !-----------------------------------------------------------------------------

  oktest=INDEX(teststr,'advance')>=1
  IF(oktest)WRITE(*,*)'Advance method,it,w:',method,' ',it,w(ixtest^D,iwtest)

  IF(ix^L/=ixG^L|.OR.|.OR.)&
       CALL die('Error in Advance: No subgrids implemented yet...')

  nexpl=nexpl+1
  firstsweep=.TRUE.
  IF(dimsplit)THEN
     IF((it/2)*2.EQ.it .OR. typedimsplit=='xy')THEN
        !If typedimsplit='xy', always do the sweeps in order of increasing idim,
        !otherwise for even parity of "it" only, and reverse order for odd. 
        DO idimsplit=1,ndim
           lastsweep= idimsplit==ndim
           CALL advect(method,ix^L,iws,idimsplit,idimsplit,w1,w,firstsweep,lastsweep)
        ENDDO
     ELSE
        ! If the parity of "it" is odd and typedimsplit=xyyx, do sweeps backwards
        DO idimsplit=ndim,1,-1
           lastsweep= idimsplit==1
           CALL advect(method,ix^L,iws,idimsplit,idimsplit,w1,w,firstsweep,lastsweep)
        ENDDO
     ENDIF
  ELSE
     ! Add fluxes from all directions at once
     lastsweep= .TRUE.
     CALL advect(method,ix^L,iws,1,ndim,w1,w,firstsweep,lastsweep)
  ENDIF

  IF(typefilter1/='nul')THEN
     ! We use updated w for the filter fluxes
     w1(ix^S,1:nw)=w(ix^S,1:nw)

     ! Filter according to typefilter1
     SELECT CASE(typefilter1)
        {^IFTVD
     CASE('tvd1')
        ! Call tvdlimit with tvd1 (2nd order Lax-Wendroff terms off)
        CALL tvdlimit(typefilter1,dt,ix^L,ix^L^LSUB2,iw_filter,1,ndim,w1,t+dt,w)
        }
        {^IFTVDLF
     CASE('tvdlf','tvdmu','tvdlf1','tvdmu1')
        ! Call tvdmusclf with filter method and physical fluxes off
        CALL tvdmusclf(.FALSE.,typefilter1,dt,ix^L,ix^L^LSUB2,iw_filter,1,ndim,&
             t+dt,w1,t+dt,w)
        }
     CASE default
        CALL die('Error in Advance: typefilter='&
             //typefilter1//' is unknown or module is switched off!')
     endselect
     CALL getboundary(t+dt,1,nw,1,ndim,w)
  ENDIF

  CALL process(0,1,ndim,w)

  IF(oktest)WRITE(*,*)'Advance new w:',w(ixtest^D,iwtest)

  RETURN
END SUBROUTINE advance_expl

!=============================================================================
SUBROUTINE advect(method,ix^L,iws,idim^LIM,w1,w,firstsweep,lastsweep) 
  ! Process w if nproc/=0:   		call process
  ! Add fluxes and unsplit sources in 
  ! directions idim=idimmin..idimmax:	call advect1
  !
  ! Depending on typeadvance and implpar call advect1 several times

  USE constants
  USE common_variables

  CHARACTER(^LENTYPE):: method
  INTEGER:: ix^L,iws(niw_),idim^LIM
  REAL(kind=8):: w1(ixG^T,nw),w(ixG^T,nw)

  ! For most Runge-Kutta type schemes one more full array is needed
  ! For classical RK4 another array is needed
  {^ANDIFRK 
  REAL(kind=8):: w2(ixG^T,nw),w3(ixG^T,nw)
  }

  !!!MEMORY Needed for typeadvance='adams2' only
  {^IFIMPL{^ANDIFRK SAVE w2}}

  LOGICAL, INTENT(INOUT) :: firstsweep, lastsweep
  !-----------------------------------------------------------------------------

  oktest=INDEX(teststr,'advect')>=1
  IF(oktest)WRITE(*,*)'Advect method w:',method,' ',w(ixtest^D,iwtest)

  ! For negative "nproc(1)" call process, if positive check whether this is the
  ! first sweep and if "it-itmin" is an integer multiple of "nproc(1)" 
  ! (the frequency of processing before the whole timestep)
  ! Processing is done in advance_impl for implicit methods
  IF(nproc(1)/=0.AND.implpar<=zero)THEN
     IF(nproc(1)<0.OR.(firstsweep.AND.it-itmin==((it-itmin)/nproc(1))*nproc(1)))&
          CALL process(1,idim^LIM,w)
  END IF

  ! Typically use "method" and at least one extra variable w1
  w1(ix^S,1:nw)=w(ix^S,1:nw)

  istep=0
  SELECT CASE(typeadvance)
  CASE('onestep')
     CALL advect1(method,dt,ix^L,iws,idim^LIM,t,w1,t,w,firstsweep,lastsweep)
     {^IFIMPL{^ANDIFRK
  CASE('adams2')
     ! w=w+dt*R+dt/2*[R-(dt/dtold)*R_n-1]
     ! Use w1=w+dt*R and w2=R_n-1/dtold
     IF(it==itmin)THEN
        CALL advect1(method,dt,ix^L,iws,idim^LIM,t,w1,t,w,firstsweep,lastsweep)
        w2(ix^S,1:nw)=(w(ix^S,1:nw)-w1(ix^S,1:nw))/dt**2
     ELSE
        CALL advect1(method,dt,ix^L,iws,idim^LIM,t,w,t,w1,firstsweep,lastsweep)
        w1(ix^S,1:nw)=w1(ix^S,1:nw)-w(ix^S,1:nw)
        w(ix^S,1:nw)=w(ix^S,1:nw)+w1(ix^S,1:nw)&
             +half*(w1(ix^S,1:nw)-dt**2*w2(ix^S,1:nw))
        w2(ix^S,1:nw)=w1(ix^S,1:nw)/dt**2
     ENDIF
     }}
  CASE('twostep')
     ! do predictor step with typepred method to calculate w1 from w, then
     ! full step with method. Fluxes and unsplit sources are taken at w1.
     CALL advect1(typepred1,dt/2,ix^L,iws,idim^LIM,t     ,w,t,w1,firstsweep,lastsweep)
     CALL advect1(method   ,dt  ,ix^L,iws,idim^LIM,t+dt/2,w1,t,w,firstsweep,lastsweep)
     {^ANDIFRK
  CASE('threestep')
     ! Shu-s third order method based on eq 2.15b of Yee II with signs corrected 
     CALL advect1(method,dt      ,ix^L,iws,idim^LIM,t     ,w ,t     ,w1,firstsweep,lastsweep)
     w2(ix^S,1:nw)=3*quarter*w(ix^S,1:nw)+quarter*w1(ix^S,1:nw)
     CALL advect1(method,dt/4    ,ix^L,iws,idim^LIM,t+dt  ,w1,t+dt/4,w2,firstsweep,lastsweep)
     w(ix^S,1:nw)=w(ix^S,1:nw)/3+w2(ix^S,1:nw)*(two/3)
     CALL advect1(method,dt*two/3,ix^L,iws,idim^LIM,t+dt/2,w2,t+dt/3,w ,firstsweep,lastsweep)
     }
     {^ANDIFRK
  CASE('fourstep')
     ! Classical four step Runge-Kutta
     ! w1=w+Dt/2*k1
     CALL advect1(method,dt/2    ,ix^L,iws,idim^LIM,t     ,w ,t,w1,firstsweep,lastsweep)
     ! w2=w+Dt/2*k2
     w2(ix^S,1:nw)=w(ix^S,1:nw)
     CALL advect1(method,dt/2    ,ix^L,iws,idim^LIM,t+dt/2,w1,t,w2,firstsweep,lastsweep)
     ! w3=w+dt*k3
     w3(ix^S,1:nw)=w(ix^S,1:nw)   
     CALL advect1(method,dt      ,ix^L,iws,idim^LIM,t+dt/2,w2,t,w3,firstsweep,lastsweep)
     ! w1=(w1+2*w2+w3)/3=Dt*(k1+2*k2+2*k3)/6
     w1(ix^S,1:nw)=(w1(ix^S,1:nw)+2*w2(ix^S,1:nw)+w3(ix^S,1:nw)-4*w(ix^S,1:nw))/3
     ! w=w+Dt*k4/6
     CALL advect1(method,dt/6    ,ix^L,iws,idim^LIM,t+dt  ,w3,t,w,firstsweep,lastsweep)
     ! w=w+Dt*(k1+2*k2+2*k3+k4)/6
     w(ix^S,1:nw)=w(ix^S,1:nw)+w1(ix^S,1:nw)
     }
     {^ANDIFRK
  CASE('sterck')
     ! H. Sterck has this fourstep time integration, w2 is needed
     CALL advect1(method,dt*0.12,ix^L,iws,idim^LIM,t        ,w ,t,w1,firstsweep,lastsweep)
     w2(ix^S,1:nw)=w(ix^S,1:nw)
     CALL advect1(method,dt/4   ,ix^L,iws,idim^LIM,t+dt*0.12,w1,t,w2,firstsweep,lastsweep)
     w1(ix^S,1:nw)=w(ix^S,1:nw)
     CALL advect1(method,dt/2   ,ix^L,iws,idim^LIM,t+dt/4   ,w2,t,w1,firstsweep,lastsweep)
     CALL advect1(method,dt     ,ix^L,iws,idim^LIM,t+dt/2   ,w1,t,w,firstsweep,lastsweep)
     }
     {^ANDIFRK
  CASE('jameson')
     ! Based on eq.2.15 of Yee II
     CALL advect1(method,dt/4,ix^L,iws,idim^LIM,t     ,w ,t,w1,firstsweep,lastsweep)
     w2(ix^S,1:nw)=w(ix^S,1:nw)
     CALL advect1(method,dt/3,ix^L,iws,idim^LIM,t+dt/4,w1,t,w2,firstsweep,lastsweep)
     w1(ix^S,1:nw)=w(ix^S,1:nw)
     CALL advect1(method,dt/2,ix^L,iws,idim^LIM,t+dt/3,w2,t,w1,firstsweep,lastsweep)
     CALL advect1(method,dt  ,ix^L,iws,idim^LIM,t+dt/2,w1,t,w,firstsweep,lastsweep)
     }
  CASE default
     WRITE(*,*)'typeadvance=',typeadvance
     WRITE(*,*)'Error in Advect: Unknown time integration method or RK is off'
     CALL die('Correct typeadvance or: cd src; setvac -on=rk; make vac')
  END SELECT

  IF(oktest)WRITE(*,*)'Advect final w:',w(ixtest^D,iwtest)

  firstsweep=.FALSE.

  RETURN
END SUBROUTINE advect

!=============================================================================
SUBROUTINE advect1(method,qdt,ixI^L,iws,idim^LIM,qtC,wCT,qt,w,firstsweep,lastsweep)

  ! Process if not first advection and nproc<0 is set
  ! Advect w to w+qdt*dF(wCT)_idim/dx_idim+qdt*((idimmax-idimmin+1)/ndim)*S(wCT)
  ! getboundaries

  USE constants
  USE common_variables

  CHARACTER(^LENTYPE) :: method
  INTEGER:: ixI^L,ixO^L,iws(niw_),idim^LIM,idim
  REAL(kind=8):: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)

  LOGICAL, INTENT(INOUT) :: firstsweep,lastsweep
  !-----------------------------------------------------------------------------

  istep=istep+1

  IF(INDEX(teststr,'saveadvect1')>=1) CALL savefile(fileout_,wCT)

  oktest=INDEX(teststr,'advect1')>=1
  IF(oktest)WRITE(*,*)'Advect1 istep,wCT,w:',&
       istep,wCT(ixtest^D,iwtest),w(ixtest^D,iwtest)

  ! In the first step wCT=w thus wCT is already processed if there is processing.
  ! Otherwise for negative "nproc(2)" call process, if positive check whether 
  ! this is the first sweep and if "it-itmin" is an integer multiple of 
  ! "nproc(2)" (the frequency of processing before intermediate steps)
  ! No processing here for implicit methods
  IF(istep>1.AND.nproc(2)/=0.AND.implpar<=zero)THEN
     IF(nproc(2)<0.OR.(firstsweep.AND.it-itmin==((it-itmin)/nproc(2))*nproc(2)))&
          CALL process(2,idim^LIM,w)
  END IF

  ! Shrink ixO^L in all directions by 2
  ixO^L=ixI^L^LSUB2;

  SELECT CASE(method)
     {^ANDIFCD

  CASE('cd4')
     CALL centdiff4(qdt,ixI^L,ixO^L,iws,idim^LIM,qtC,wCT,qt,w)
     }

  CASE('source')
     IF(sourceunsplit)CALL addsource2(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)
  CASE('nul')
     ! There is nothing to do
     !HPF_ if(.false.)write(*,*)'This avoids an xlhpf compiler bug'
  CASE default
     WRITE(*,*)'Error in Advect1:',method,' is unknown or switched off!'
     CALL die('Error in Advect1:'//method//' is unknown or switched off!')
  END SELECT

  CALL getboundary(qt+qdt,1,nw,1,ndim,w)

  IF(oktest)WRITE(*,*)'Advect1 final w:',w(ixtest^D,iwtest)

  RETURN
END SUBROUTINE advect1

!=============================================================================
SUBROUTINE addsource2(qdt,ixII^L,ixOO^L,iws,qtC,wCT,qt,w)

  ! Add source within ixOO for iws: w=w+qdt*S[wCT]

  USE constants
  USE common_variables

  INTEGER:: ixI^L,ixO^L,ixII^L,ixOO^L,iws(niw_)
  REAL(kind=8):: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)
  !-----------------------------------------------------------------------------

  oktest=INDEX(teststr,'addsource')>=1
  IF(oktest)WRITE(*,*)'Add Source qdt,wCT,w:',&
       qdt,wCT(ixtest^D,iwtest),w(ixtest^D,iwtest)

  ! AddSource and SpecialSource may shrink ixO or expand ixI for derivatives 
  ixI^L=ixII^L; ixO^L=ixOO^L;

  CALL specialsource(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

  CALL     addsource(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

  IF(oktest)WRITE(*,*)'wnew:',w(ixtest^D,iwtest)

  ! If AddSource/SpecialSource shrunk ixO, getboundary is needed.
  IF(ixO^L^LTixOO^L|.OR.|.OR.)THEN
     CALL getboundary(qt+qdt,1,nw,1,ndim,w)
     IF(oktest)WRITE(*,*)'wnew after getboundary:',w(ixtest^D,iwtest)
  END IF

  RETURN
END SUBROUTINE addsource2

!=============================================================================
LOGICAL FUNCTION timetofinish(time0)

  ! Finish when it or t reached its maximum expected value, or dt is too small,
  ! or residual is small enough. Other conditions may be included.

  USE constants
  USE common_variables

  REAL(kind=8):: time0, cputime
  LOGICAL:: okfinish
  !-----------------------------------------------------------------------------

  okfinish = it>=itmax .OR. t>=tmax .OR. dt<dtmin .OR. &
       (it>itmin.AND.(residual<residmin .OR. residual>residmax))

  IF(cputimemax < bigdouble .AND. .NOT.okfinish) &
       okfinish= cputimemax <= cputime()-time0

  timetofinish=okfinish

  RETURN
END FUNCTION timetofinish

!=============================================================================
LOGICAL FUNCTION timetosave(ifile)

  ! Save times are defined by either tsave(isavet(ifile),ifile) or 
  ! itsave(isaveit(ifile),ifile) or dtsave(ifile) or ditsave(ifile)
  ! Other conditions may be included.

  USE constants
  USE common_variables

  INTEGER:: ifile
  LOGICAL:: oksave
  !-----------------------------------------------------------------------------

  oksave=.FALSE.
  IF(t>=tsave(isavet(ifile),ifile))THEN
     oksave=.TRUE.
     isavet(ifile)=isavet(ifile)+1
  END IF
  IF(it==itsave(isaveit(ifile),ifile))THEN
     oksave=.TRUE.
     isaveit(ifile)=isaveit(ifile)+1
  END IF
  IF(it==itsavelast(ifile)+ditsave(ifile)) oksave=.TRUE.
  IF(t >=tsavelast(ifile) +dtsave(ifile) ) oksave=.TRUE.
  IF(oksave)THEN
     tsavelast(ifile) =t
     itsavelast(ifile)=it
  END IF
  timetosave=oksave

  RETURN
END FUNCTION timetosave

!=============================================================================
SUBROUTINE getdt_courant(w,ix^L)

  ! Ensure that the courant conditions is met
  ! Calculate the time for the  maximum propagation speed cmax_i to cross dx_i
  ! in each i directions then take minimum for all grid points in the mesh and 
  ! for all i directions, finally multiply by courantpar.
  !
  ! If TVD or TDLF provides dtcourant(idim) we take the minimum of those.
  ! In case of generalized coordinates dtcourant(idim) is correct due to the
  ! rotations while the value calculated here does not use a rotation.

  USE constants
  USE common_variables

  REAL(kind=8):: w(ixG^T,nw),cmax(ixG^T),courantmax,dtold
  INTEGER:: ix^L,idim
  LOGICAL:: new_cmax
  !-----------------------------------------------------------------------------

  oktest=INDEX(teststr,'getdt')>=1

  IF(oktest) WRITE(*,*)'getdt_courant'

  dtold=dt
  dt=bigdouble
  courantmax=zero
  new_cmax=.TRUE.
  DO idim=1,ndim
     IF(dtcourant(idim)<bigdouble)THEN
        ! If dtcourant(idim) is calculated, use it
!!!      if(it==itmin+1)write(*,*)'second order correction in dt_courant!!!'
!!!      dt=min(dt,dtcourant(idim),dtcourant(idim)**2/dtold,1.1*dtold)
        dt=MIN(dt,dtcourant(idim))
        IF(oktest) WRITE(*,*)'idim,dtcourant(idim)',idim,dtcourant(idim)
     ELSE
        ! dx>0, but cmax>=0 may actually be 0, thus we calculate 
        ! max(cmax/dx) rather than min(dx/cmax).

        CALL getcmax(new_cmax,w,ix^L,idim,cmax)
        courantmax=MAX(courantmax,MAXVAL(cmax(ix^S)/dx(ix^S,idim)))

        IF(gencoord.AND.it==itmin+1.AND.verbose)WRITE(*,*)&
             'Warning in GetDtCourant: for gencoord approx. only',&
             ', better use TVD-type methods'
        IF(oktest) WRITE(*,*)'idim,cmax:',idim,cmax(ixtest^D)
        IF(oktest) WRITE(*,*)'max(c/dx)',MAXVAL(cmax(ix^S)/dx(ix^S,idim))
     ENDIF
  END DO
  {^IFMPI CALL mpiallreduce(courantmax,MPI_MAX)}
  IF(INDEX(teststr,'dtdecline')<1)THEN
     DO idim=1,ndim
        dtcourant(idim)=bigdouble
     ENDDO
  ENDIF
  IF(courantmax>smalldouble) dt=MIN(dt,courantpar/courantmax)

  IF(oktest) WRITE(*,*)'GetDtCourant dt=',dt

  RETURN 
END SUBROUTINE getdt_courant

!=============================================================================
REAL(kind=8) FUNCTION cputime()

  ! Return cputime in seconds as a double precision number.
  ! For g77 compiler replace F77_ with F77_ everywhere in this function
  ! so that f90tof77 does not touch the system_clock function.

  INTEGER:: clock,clockrate,count_max !HPF_ !F77_
  !F77_ real:: etime,total,tarray(2)
  !F77_ external etime
  !HPF_ real:: timef
  !-----------------------------------------------------------------------------

  cputime=-1.D0                             ! No timing

  CALL SYSTEM_CLOCK(clock,clockrate,count_max) !HPF_ !F77_
  cputime=clock*(1.D0/clockrate)               !HPF_ !F77_
  !F77_ total = etime(tarray)
  !F77_ cputime=tarray(1)
  !HPF_ cputime=timef()/1.0D3

  !cputime=second()   ! Cray CF77 or F90 (total CPU time for more CPU-s)
  !cputime=secondr()  ! Cray CF77 or F90 (elapsed time for more CPU-s)

  RETURN
END FUNCTION cputime
!=============================================================================
! end module vac
!##############################################################################
