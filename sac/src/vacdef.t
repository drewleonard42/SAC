MODULE phys_constants
  IMPLICIT NONE
  SAVE
  
  INTEGER, PARAMETER :: biginteger=10000000
  REAL(kind=8), PARAMETER :: pi= 3.1415926535897932384626433832795
  real(kind=8), parameter :: mu0=4.d-7*pi
  REAL(kind=8), PARAMETER :: smalldouble=1.d-99, bigdouble=1.d+99
  REAL(kind=8), PARAMETER :: zero=0d0, one=1d0, two=2d0, half=0.5d0, quarter=0.25d0
  
END MODULE phys_constants

MODULE code_constants
  IMPLICIT NONE
  SAVE

  ! DEFINITIONS OF GLOBAL PARAMETERS
  ! Indices for cylindrical coordinates FOR TESTS, negative value when not used:
  INTEGER, PARAMETER :: r_=1, phi_=^PHI, z_=^Z

  ! Indices for cylindrical coordinates FOR INDEXING, always positive
  INTEGER, PARAMETER :: pphi_=^PPHI, zz_=^ZZ

  INTEGER, PARAMETER :: ixGlo^D=1

  ! The next line is edited by SETVAC
  INTEGER, PARAMETER :: ixGhi1=100,ixGhi2=100,ixGhimin=100,ixGhimax=100
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER, PARAMETER :: ndim=^ND, ndir=^NC

  INTEGER, PARAMETER :: dixBlo=2,dixBhi=2

  INTEGER, PARAMETER :: nhiB=10           ! maximum No. boundary sections

  INTEGER, PARAMETER :: nsavehi=100       ! maximum No. saves into outputfiles
  ! defined by arrays of tsave or itsave


  INTEGER, PARAMETER :: filelog_=1,fileout_=2,nfile=2 ! outputfiles

  INTEGER, PARAMETER :: unitstdin=5,unitterm=6,uniterr=6,unitini=10 ! Unit names. 
  ! Outputfiles use unitini+1..initini+nfile
  ! Default parfiles uses unitini-1

  INTEGER, PARAMETER :: toosmallp_=1, toosmallr_=2, couranterr_=3, poissonerr_=4
  INTEGER, PARAMETER :: nerrcode=4

  {^IFPOISSON
  INTEGER, PARAMETER :: nwrk=7 !Size of work array for VACPOISSON
  } 

END MODULE code_constants

MODULE constants
  USE phys_constants
  USE code_constants
  {^IFMPI
  ! MPI header file
  USE mpi !include 'mpif.h'
  }
  IMPLICIT NONE
  SAVE

  !Old imports really could do with making into modules
  include 'vacpar.f90'
  include 'vacusrpar.f90'

  INTEGER, PARAMETER:: nhi=nw*{ixGhi^D*} ! Maximum number of unknowns for VACIMPL

  INTEGER, PARAMETER :: niw_=nw+1         !Indexname for size of iw index array

  !##############################################################################
  ! Define MPI only things 
  !##############################################################################
  {^IFMPI
  ! Size of work array for VACIMPL
  INTEGER, PARAMETER :: nwork=(7+nw*(2*ndim+1))*{ixGhi^D*}*nw 
  ! Max. number of implicit variables
  INTEGER, PARAMETER :: nwimplhi=nw 

  ! Buffer size for sending and receiving ghost cells.
  ! Both w and x are sent with MPI.

  ! The next line is edited by SETVAC
  INTEGER, PARAMETER :: maxndimnw = nw
  INTEGER, PARAMETER :: nmpibuffer = dixBhi*maxndimnw*{ixGhi^D*} / ixGhimin
  }

END MODULE constants


module common_variables
  USE constants
  SAVE

  {^IFMPI 
  INTEGER:: ipe, ipe^D, npe, npe^D, nxall^D, nxpe^D, ierrmpi
  INTEGER:: ixPEmin^D, ixPEmax^D
  LOGICAL:: mpiupperB(ndim),mpilowerB(ndim)
  REAL(kind=8):: sendbuffer(nmpibuffer)
  REAL(kind=8):: recvbuffer(nmpibuffer,2)
  }

  ! Unit for reading input parameters.
  INTEGER :: unitpar

  ! Logical to set verbosity. For MPI parallel run only PE 0 is verbose
  LOGICAL :: verbose

  ! General temporary arrays, any subroutine call may change them 
  ! except for subroutines which say the opposite in their header
  REAL(kind=8):: tmp(ixG^T),tmp2(ixG^T)

  ! Number of errors during calculation
  INTEGER:: nerror(nerrcode)

  !Kronecker delta and Levi-Civita tensors
  INTEGER:: kr(3,3),lvc(3,3,3)

  !Grid parameters
  INTEGER:: ixM^L,ixG^L,nx^D,nx(ndim)
  INTEGER:: dixB^L

  ! Global x array edges
  real(kind=8), dimension(3) :: x_left_edge = (/ 0, 0, 0 /), x_right_edge = (/ 1, 1, 1 /)
  ! x and dx are local for HPF
  REAL(kind=8):: x(IXG^T,ndim),dx(IXG^T,ndim)
  REAL(kind=8):: volume,dvolume(IXG^T)
  REAL(kind=8):: area(IXGLO1:IXGHI1),areaC(IXGLO1:IXGHI1)
  REAL(kind=8):: areadx(IXGLO1:IXGHI1),areaside(IXGLO1:IXGHI1)

  ! Variables for generalized coordinates and polargrid
  LOGICAL::          gencoord, polargrid
  {^IFGEN REAL(kind=8):: surfaceC(IXG^T,ndim),normalC(IXG^T,ndim,ndim)}
  {^NOGEN REAL(kind=8):: surfaceC(2^D&,ndim), normalC(2^D&,ndim,ndim)}

  !Boundary region parameters
  REAL(kind=8):: fixB^D(-dixBlo:dixBhi^D%ixGLO^DD:ixGHI^DD,nw)
  INTEGER:: nB,ixB^LIM(ndim,nhiB),idimB(nhiB),ipairB(nhiB)
  LOGICAL:: upperB(nhiB),fixedB(nw,nhiB),nofluxB(nw,ndim),extraB
  CHARACTER(^LENTYPE) :: typeB(nw,nhiB),typeBscalar(nhiB)

  !Equation and method parameters
  real(kind=8), target :: eqpar(neqpar+nspecialpar),procpar(nprocpar)

  ! Time step control parameters
  REAL(kind=8):: courantpar,dtpar,dtdiffpar,dtcourant(ndim),dtmrpc
  LOGICAL:: dtcantgrow
  INTEGER:: slowsteps

  ! Parameters for the implicit techniques
  {^IFPOISSON REAL(kind=8):: wrk(ixG^T,nwrk) } 
  {^IFIMPL REAL(kind=8):: work(nwork) }
  INTEGER:: nwimpl,nimpl
  REAL(kind=8):: implpar,impldiffpar,implerror,implrelax,impldwlimit
  INTEGER:: implrestart,implrestart2,impliter,impliternr,implmrpcpar
  CHARACTER(^LENTYPE) :: typeimplinit,typeimpliter,typeimplmat
  LOGICAL:: implconserv,implnewton,implcentered,implnewmat
  LOGICAL:: implpred,impl3level,impljacfast,implsource

  !Method switches
  INTEGER:: iw_full(niw_),iw_semi(niw_),iw_impl(niw_),iw_filter(niw_)
  INTEGER:: iw_vector(nvector+1),vectoriw(nw)
  ! The upper bound+1 in iw_vector avoids F77 warnings when nvector=0
  CHARACTER(^LENTYPE) :: typefull1,typepred1,typeimpl1,typefilter1
  CHARACTER(^LENTYPE) :: typelimited,typefct,typetvd,typeaxial
  CHARACTER(^LENTYPE) :: typepoisson, typeconstrain
  CHARACTER(^LENTYPE) :: typelimiter(nw),typeentropy(nw)
  CHARACTER(^LENTYPE) :: typeadvance, typedimsplit, typesourcesplit
  LOGICAL:: dimsplit,sourcesplit,sourceunsplit,artcomp(nw),useprimitive
  LOGICAL:: divbfix,divbwave,divbconstrain,angmomfix,compactres,smallfix
  INTEGER:: idimsplit
  INTEGER:: nproc(nfile+2)
  REAL(kind=8):: entropycoef(nw),constraincoef
  REAL(kind=8):: smallp,smallpcoeff,smallrho,smallrhocoeff,vacuumrho
  REAL(kind=8):: muscleta1,muscleta2,musclomega,acmcoef(nw),acmexpo
  LOGICAL:: acmnolim, fourthorder
  INTEGER:: acmwidth

  !Previous time step and residuals
  REAL(kind=8):: wold(ixG^T,nw),residual,residmin,residmax

  ! Flux storage for flux-CT and flux-CD methods !!! for MHD only !!! 
  {^IFCT REAL(kind=8):: fstore(ixG^T,ndim) }

  !Time parameters
  INTEGER:: step,istep,nstep,it,itmin,itmax,nexpl,nnewton,niter,nmatvec
  REAL(kind=8):: t,tmax,dt,dtmin,cputimemax
  LOGICAL:: tmaxexact
  REAL(kind=8):: tsave(nsavehi,nfile),tsavelast(nfile),dtsave(nfile)
  INTEGER:: itsave(nsavehi,nfile),itsavelast(nfile),ditsave(nfile)
  INTEGER:: isavet(nfile),isaveit(nfile)

  !File parameters
  CHARACTER(^LENNAME) :: filenameini,filenameout,filename(nfile)
  CHARACTER(^LENNAME) :: fileheadini,fileheadout,varnames,wnames
  CHARACTER(^LENTYPE) :: typefileini,typefileout,typefilelog
  LOGICAL::             fullgridini,fullgridout
  INTEGER::             snapshotini,snapshotout,isaveout

  !Test parameters
  CHARACTER(^LENNAME) :: teststr
  INTEGER:: ixtest1,ixtest2,ixtest3,iwtest,idimtest,ipetest^IFMPI
  LOGICAL:: oktest    !This is a local variable for all subroutines and functions

  REAL(kind=8):: maxviscoef

end module common_variables
