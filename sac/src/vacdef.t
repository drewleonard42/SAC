module phys_constants
  implicit none
  save
  
  integer, parameter :: biginteger=10000000
  double precision, parameter :: pi= 3.1415926535897932384626433832795
  double precision, parameter :: smalldouble=1.d-99, bigdouble=1.d+99
  double precision, parameter :: zero=0d0, one=1d0, two=2d0, half=0.5d0, quarter=0.25d0
  
end module phys_constants

module code_constants
  implicit none
  save

  ! DEFINITIONS OF GLOBAL PARAMETERS
  ! Indices for cylindrical coordinates FOR TESTS, negative value when not used:
  integer, parameter :: r_=1, phi_=^PHI, z_=^Z

  ! Indices for cylindrical coordinates FOR INDEXING, always positive
  integer, parameter :: pphi_=^PPHI, zz_=^ZZ

  integer, parameter :: ixGlo^D=1

  ! The next line is edited by SETVAC
  integer, parameter :: ixGhi1=104,ixGhi2=104,ixGhimin=104,ixGhimax=104
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, parameter :: ndim=^ND, ndir=^NC

  integer, parameter :: dixBlo=2,dixBhi=2

  integer, parameter :: nhiB=10           ! maximum No. boundary sections

  integer, parameter :: nsavehi=100       ! maximum No. saves into outputfiles
  ! defined by arrays of tsave or itsave


  integer, parameter :: filelog_=1,fileout_=2,nfile=2 ! outputfiles

  integer, parameter :: unitstdin=5,unitterm=6,uniterr=6,unitini=10 ! Unit names. 
  ! Outputfiles use unitini+1..initini+nfile
  ! Default parfiles uses unitini-1

  integer, parameter :: toosmallp_=1, toosmallr_=2, couranterr_=3, poissonerr_=4
  integer, parameter :: nerrcode=4

  {^IFPOISSON
  integer, parameter :: nwrk=7 !Size of work array for VACPOISSON
  } 

end module code_constants

module constants
  use phys_constants
  use code_constants
  {^IFMPI
  ! MPI header file
  use mpi !include 'mpif.h'
  }
  implicit none
  save

  !Old imports really could do with making into modules
  include 'vacpar.f90'
  include 'vacusrpar.f90'

  integer, parameter:: nhi=nw*{ixGhi^D*} ! Maximum number of unknowns for VACIMPL

  integer, parameter :: niw_=nw+1         !Indexname for size of iw index array

  !##############################################################################
  ! Define MPI only things 
  !##############################################################################
  {^IFMPI
  ! Size of work array for VACIMPL
  integer, parameter :: nwork=(7+nw*(2*ndim+1))*{ixGhi^D*}*nw 
  ! Max. number of implicit variables
  integer, parameter :: nwimplhi=nw 

  ! Buffer size for sending and receiving ghost cells.
  ! Both w and x are sent with MPI.

  ! The next line is edited by SETVAC
  integer, parameter :: maxndimnw = nw
  integer, parameter :: nmpibuffer = dixBhi*maxndimnw*{ixGhi^D*} / ixGhimin
  }

end module constants


module common_varibles
  use constants
  save

  {^IFMPI 
  integer:: ipe, ipe^D, npe, npe^D, nxall^D, nxpe^D, ierrmpi
  integer:: ixPEmin^D, ixPEmax^D
  logical:: mpiupperB(ndim),mpilowerB(ndim)
  double precision:: sendbuffer(nmpibuffer)
  double precision:: recvbuffer(nmpibuffer,2)
  }

  ! Unit for reading input parameters.
  integer :: unitpar

  ! Logical to set verbosity. For MPI parallel run only PE 0 is verbose
  logical :: verbose

  ! General temporary arrays, any subroutine call may change them 
  ! except for subroutines which say the opposite in their header
  double precision:: tmp(ixG^T),tmp2(ixG^T)

  ! Number of errors during calculation
  integer:: nerror(nerrcode)

  !Kronecker delta and Levi-Civita tensors
  integer:: kr(3,3),lvc(3,3,3)

  !Grid parameters
  integer:: ixM^L,ixG^L,nx^D,nx(ndim)
  integer:: dixB^L
  ! x and dx are local for HPF
  double precision:: x(IXG^T,ndim),dx(IXG^T,ndim)
  double precision:: volume,dvolume(IXG^T)
  double precision:: area(IXGLO1:IXGHI1),areaC(IXGLO1:IXGHI1)
  double precision:: areadx(IXGLO1:IXGHI1),areaside(IXGLO1:IXGHI1)

  ! Variables for generalized coordinates and polargrid
  logical::          gencoord, polargrid
  {^IFGEN double precision:: surfaceC(IXG^T,ndim),normalC(IXG^T,ndim,ndim)}
  {^NOGEN double precision:: surfaceC(2^D&,ndim), normalC(2^D&,ndim,ndim)}

  !Boundary region parameters
  double precision:: fixB^D(-dixBlo:dixBhi^D%ixGLO^DD:ixGHI^DD,nw)
  integer:: nB,ixB^LIM(ndim,nhiB),idimB(nhiB),ipairB(nhiB)
  logical:: upperB(nhiB),fixedB(nw,nhiB),nofluxB(nw,ndim),extraB
  character(^LENTYPE) :: typeB(nw,nhiB),typeBscalar(nhiB)

  !Equation and method parameters
  double precision:: eqpar(neqpar+nspecialpar),procpar(nprocpar)

  ! Time step control parameters
  double precision:: courantpar,dtpar,dtdiffpar,dtcourant(ndim),dtmrpc
  logical:: dtcantgrow
  integer:: slowsteps

  ! Parameters for the implicit techniques
  {^IFPOISSON double precision:: wrk(ixG^T,nwrk) } 
  {^IFIMPL double precision:: work(nwork) }
  integer:: nwimpl,nimpl
  double precision:: implpar,impldiffpar,implerror,implrelax,impldwlimit
  integer:: implrestart,implrestart2,impliter,impliternr,implmrpcpar
  character(^LENTYPE) :: typeimplinit,typeimpliter,typeimplmat
  logical:: implconserv,implnewton,implcentered,implnewmat
  logical:: implpred,impl3level,impljacfast,implsource

  !Method switches
  integer:: iw_full(niw_),iw_semi(niw_),iw_impl(niw_),iw_filter(niw_)
  integer:: iw_vector(nvector+1),vectoriw(nw)
  ! The upper bound+1 in iw_vector avoids F77 warnings when nvector=0
  character(^LENTYPE) :: typefull1,typepred1,typeimpl1,typefilter1
  character(^LENTYPE) :: typelimited,typefct,typetvd,typeaxial
  character(^LENTYPE) :: typepoisson, typeconstrain
  character(^LENTYPE) :: typelimiter(nw),typeentropy(nw)
  character(^LENTYPE) :: typeadvance, typedimsplit, typesourcesplit
  logical:: dimsplit,sourcesplit,sourceunsplit,artcomp(nw),useprimitive
  logical:: divbfix,divbwave,divbconstrain,angmomfix,compactres,smallfix
  integer:: idimsplit
  integer:: nproc(nfile+2)
  double precision:: entropycoef(nw),constraincoef
  double precision:: smallp,smallpcoeff,smallrho,smallrhocoeff,vacuumrho
  double precision:: muscleta1,muscleta2,musclomega,acmcoef(nw),acmexpo
  logical:: acmnolim, fourthorder
  integer:: acmwidth

  !Previous time step and residuals
  double precision:: wold(ixG^T,nw),residual,residmin,residmax

  ! Flux storage for flux-CT and flux-CD methods !!! for MHD only !!! 
  {^IFCT double precision:: fstore(ixG^T,ndim) }

  !Time parameters
  integer:: step,istep,nstep,it,itmin,itmax,nexpl,nnewton,niter,nmatvec
  double precision:: t,tmax,dt,dtmin,cputimemax
  logical:: tmaxexact
  double precision:: tsave(nsavehi,nfile),tsavelast(nfile),dtsave(nfile)
  integer:: itsave(nsavehi,nfile),itsavelast(nfile),ditsave(nfile)
  integer:: isavet(nfile),isaveit(nfile)

  !File parameters
  character(^LENNAME) :: filenameini,filenameout,filename(nfile)
  character(^LENNAME) :: fileheadini,fileheadout,varnames,wnames
  character(^LENTYPE) :: typefileini,typefileout,typefilelog
  logical::             fullgridini,fullgridout
  integer::             snapshotini,snapshotout,isaveout

  !Test parameters
  character(^LENNAME) :: teststr
  integer:: ixtest1,ixtest2,ixtest3,iwtest,idimtest,ipetest^IFMPI
  logical:: oktest    !This is a local variable for all subroutines and functions

  double precision:: maxviscoef

end module common_varibles
