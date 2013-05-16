!#############################################################################
! module vacgrid
! Subroutines for boundaries, grid, divergence of flux, gradients
! Also the limiter functions for TVD, TVDLF, TVDMU schemes


!INCLUDE:vacgrid.setnozzle.t
!=============================================================================
subroutine boundsetup

! Variables describing boundaries
!
!    iB=1..nB                                    - boundary ID
!    ixBmin(idim,iB):ixBmax(idim,iD)             - location of boundary
!    idimB(iB)                                   - direction orthogonal to the 
!                                                  boundary layer
!    upperB(iB)                                  - .true. if boundary is at a
!                                                  high index of the phys. grid
!    typeB(iw,iB)                                - boundary type string

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,iB,jB,iw,idim,idm,&
   ixGmin(ndim),ixGmax(ndim),ixMmin(ndim),ixMmax(ndim)
!-----------------------------------------------------------------------------

oktest=index(teststr,'boundsetup')>=1
if(oktest)write(*,*)'BoundSetup'

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

ixGmin(1)=ixGmin1;ixGmax(1)=ixGmax1;ixGmin(2)=ixGmin2;ixGmax(2)=ixGmax2
ixGmin(3)=ixGmin3;ixGmax(3)=ixGmax3;
ixMmin(1)=ixMmin1;ixMmax(1)=ixMmax1;ixMmin(2)=ixMmin2;ixMmax(2)=ixMmax2
ixMmin(3)=ixMmin3;ixMmax(3)=ixMmax3;

if(ixBmax(1,1)==0)then
   nB=2*ndim
   do iB=1,nB
      idim=(iB+1)/2
      idimB(iB) = idim
      upperB(iB)= 2*idim==iB
      do idm=1,ndim
         ixBmin(idm,iB)=ixGmin(idm);ixBmax(idm,iB)=ixGmax(idm);
      end do
      if(upperB(iB))then
         ixBmin(idim,iB)=ixMmax(idim)+1
         ixBmax(idim,iB)=ixGmax(idim)
      else
         ixBmin(idim,iB)=ixGmin(idim)
         ixBmax(idim,iB)=ixMmin(idim)-1
      end if
   end do
else
   ! Check if boundary regions are inside grid and in between mesh and grid
   do iB=1,nB
      ixmin1=ixBmin(1,iB);ixmin2=ixBmin(2,iB);ixmin3=ixBmin(3,iB)
      ixmax1=ixBmax(1,iB);ixmax2=ixBmax(2,iB);ixmax3=ixBmax(3,iB);
           
      ! Convert global limits into local grid limits 
      ixmin1=ixmin1-ixpemin1+1;ixmin2=ixmin2-ixpemin2+1
      ixmin3=ixmin3-ixpemin3+1;ixmax1=ixmax1-ixpemin1+1
      ixmax2=ixmax2-ixpemin2+1;ixmax3=ixmax3-ixpemin3+1;
      ! Cut off extra cells at MPI boundaries
      if(ipe1>0)ixmin1=max(ixmin1,ixGmin1)
           if(ipe2>0)ixmin2=max(ixmin2,ixGmin2)
           if(ipe3>0)ixmin3=max(ixmin3,ixGmin3)
      if(ipe1<npe1-1)ixmax1=min(ixmax1,ixGmax1)
           if(ipe2<npe2-1)ixmax2=min(ixmax2,ixGmax2)
           if(ipe3<npe3-1)ixmax3=min(ixmax3,ixGmax3)
 !Change the global ixBmin(1,2,3,iB),ixBmax(1,2,3,iB) to the local index limits
      ixBmin(1,iB)=ixmin1;ixBmin(2,iB)=ixmin2;ixBmin(3,iB)=ixmin3;
      ixBmax(1,iB)=ixmax1;ixBmax(2,iB)=ixmax2;ixBmax(3,iB)=ixmax3;
      ! Check if local boundary region is empty
      if(ixmin1>ixmax1.or.ixmin2>ixmax2.or.ixmin3>ixmax3)then
          ! The limits do not contain cells, change the boundary type to 'mpi'
          typeB(1:nw,iB)='mpi'
          CYCLE
      endif
      
      if(ixmin1<ixGmin1.or.ixmax1>ixGmax1.or.ixmin2<ixGmin2&
         .or.ixmax2>ixGmax2.or.ixmin3<ixGmin3.or.ixmax3>ixGmax3) then
         write(*,*)'Error for boundary region iB,ixL=',iB,ixmin1,ixmin2,&
            ixmin3,ixmax1,ixmax2,ixmax3
         call die('Error in BoundSetup: Boundary region is outside grid')
      endif
      select case(idimB(iB))
      case(1)
         if(upperB(iB))then
            if(ixmin1-1/=ixMmax1.or.ixmax1/=ixGmax1)write(*,*)&
               'Warning in BoundSetup: Boundary does not fit, iB:',iB
         else
            if(ixmax1+1/=ixMmin1.or.ixmin1/=ixGmin1)write(*,*)&
               'Warning in BoundSetup: Boundary does not fit, iB:',iB
         endif 
      case(2)
         if(upperB(iB))then
            if(ixmin2-1/=ixMmax2.or.ixmax2/=ixGmax2)write(*,*)&
               'Warning in BoundSetup: Boundary does not fit, iB:',iB
         else
            if(ixmax2+1/=ixMmin2.or.ixmin2/=ixGmin2)write(*,*)&
               'Warning in BoundSetup: Boundary does not fit, iB:',iB
         endif 
      case(3)
         if(upperB(iB))then
            if(ixmin3-1/=ixMmax3.or.ixmax3/=ixGmax3)write(*,*)&
               'Warning in BoundSetup: Boundary does not fit, iB:',iB
         else
            if(ixmax3+1/=ixMmin3.or.ixmin3/=ixGmin3)write(*,*)&
               'Warning in BoundSetup: Boundary does not fit, iB:',iB
         endif 
      end select
   end do
end if

! Identify the periodic pairs if they are not defined in boundlist
! Check type, direction, orientation, and size before matching pairs
do iB=1,nB
   if(typeB(1,iB)=='periodic'.and.ipairB(iB)==0)then
      do jB=iB+1,nB
         if(typeB(1,jB)=='periodic'.and.ipairB(jB)==0.and.idimB(iB)&
            ==idimB(jB).and.(upperB(iB).neqv.upperB(jB)).and.ixBmax(1,iB)&
            -ixBmin(1,iB)==ixBmax(1,jB)-ixBmin(1,jB).and.ixBmax(2,iB)&
            -ixBmin(2,iB)==ixBmax(2,jB)-ixBmin(2,jB).and.ixBmax(3,iB)&
            -ixBmin(3,iB)==ixBmax(3,jB)-ixBmin(3,jB))then
            ipairB(iB)=jB; ipairB(jB)=iB
         endif
      end do
      if(ipairB(iB)==0)call die('Error in BoundSetup: No periodic pair')
   end if
end do

     
do iB=1,nB
   ! Change boundary type if processor is not at the edge or if this is a
   ! periodic boundary and there are more than 1 processors in this direction
   select case(idimB(iB))
   case(1)
      if(typeB(1,iB)=='periodic' .and. npe1>1) then
          typeB(1:nw,iB)='mpiperiod'
      else if(upperB(iB) .and. ipe1<npe1-1 .or. &
              .not.upperB(iB) .and. ipe1>0) then
          typeB(1:nw,iB)='mpi'
      endif 
     case(2)
      if(typeB(1,iB)=='periodic' .and. npe2>1) then
          typeB(1:nw,iB)='mpiperiod'
      else if(upperB(iB) .and. ipe2<npe2-1 .or. &
              .not.upperB(iB) .and. ipe2>0) then
          typeB(1:nw,iB)='mpi'
      endif 
     case(3)
      if(typeB(1,iB)=='periodic' .and. npe3>1) then
          typeB(1:nw,iB)='mpiperiod'
      else if(upperB(iB) .and. ipe3<npe3-1 .or. &
              .not.upperB(iB) .and. ipe3>0) then
          typeB(1:nw,iB)='mpi'
      endif 
   end select
end do



! symm0 means zero orthogonal flux via the boundary 
do iw=1,nw
   do iB=1,nB
     if(typeB(iw,iB)=='symm0') nofluxB(iw,idimB(iB))=.true.
   enddo
enddo

if(oktest)then
   do iB=1,nB
      write(unitterm,*)'ipe=',ipe      
      write(unitterm,*)'iB,idimB,upperB,typeB:',iB,idimB(iB),upperB(iB),' ',&
         typeB(1,iB)
      write(unitterm,*)'ixBmin:',(ixBmin(idm,iB),idm=1,ndim)
      write(unitterm,*)'ixBmax:',(ixBmax(idm,iB),idm=1,ndim)
   end do
   do iB=1,nB
      if(typeB(1,iB)=='periodic')write(*,*)'Pairs:',iB,ipairB(iB)
   end do
end if

return
end

!=============================================================================
subroutine ensurebound(dix,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w)

! Check if there is enough information for calculating derivatives.
! The requirement is that ixI is wider than ixO by dix. 
! Adjust ixI and ixO. Call getboundary if needed.

include 'vacdef.f'

integer:: dix,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
   ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision:: qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)
!-----------------------------------------------------------------------------

oktest=index(teststr,'ensurebound')>0
if(oktest)write(*,*)'EnsureBound dix,ixI,ixO:',dix,',',ixImin1,ixImin2,&
   ixImin3,ixImax1,ixImax2,ixImax3,',',ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
   ixOmax2,ixOmax3

! Check wether ixO+dix is within the grid
if(ixGmin1>ixOmin1-dix.or.ixGmin2>ixOmin2-dix.or.ixGmin3>ixOmin3&
   -dix.or.ixGmax1<ixOmax1+dix.or.ixGmax2<ixOmax2+dix.or.ixGmax3<ixOmax3&
   +dix)then
   ixOmin1=ixMmin1;ixOmin2=ixMmin2;ixOmin3=ixMmin3;ixOmax1=ixMmax1
   ixOmax2=ixMmax2;ixOmax3=ixMmax3;
endif
! Check whether ixI is wider than ixO by at least dix otherwise getboundary
if(ixImin1>ixOmin1-dix.or.ixImin2>ixOmin2-dix.or.ixImin3>ixOmin3&
   -dix.or.ixImax1<ixOmax1+dix.or.ixImax2<ixOmax2+dix.or.ixImax3<ixOmax3&
   +dix)then
   ixImin1=ixGmin1;ixImin2=ixGmin2;ixImin3=ixGmin3;ixImax1=ixGmax1
   ixImax2=ixGmax2;ixImax3=ixGmax3;
   call getboundary(qt,1,nw,1,ndim,w)
   if(oktest)write(*,*)'calling getboundary'
end if

if(oktest)write(*,*)'Final       dix,ixI,ixO:',dix,',',ixImin1,ixImin2,&
   ixImin3,ixImax1,ixImax2,ixImax3,',',ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
   ixOmax2,ixOmax3

return
end

!=============================================================================
subroutine getboundary(qt,iwmin,iwmax,idimmin,idimmax,w)

include 'vacdef.f'

integer:: iwmin,iwmax,idimmin,idimmax
double precision:: qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)
integer:: ix,ix1,ix2,ix3,ixe,ixf,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,&
   ixpairmin1,ixpairmin2,ixpairmin3,ixpairmax1,ixpairmax2,ixpairmax3,idim,iw,&
   iB
integer:: iwv,jdim
double precision:: coeffnormal,coefftransv
logical:: initialized

integer:: ixee
integer:: ixbmin1,ixbmin2,ixbmin3,ixbmax1,ixbmax2,ixbmax3

data initialized/.false./
!-----------------------------------------------------------------------------


ixbmin1=ixGlo1;ixbmin2=ixGlo2;ixbmin3=ixGlo3;ixbmax1=ixGhi1;ixbmax2=ixGhi2
ixbmax3=ixGhi3;

oktest=index(teststr,'getboundary')>=1
if(oktest)write(*,*)'GetBoundary it,step:',it,step
if(oktest)write(*,*)'GetBoundary wold:',w(ixtest1,ixtest2,ixtest3,iwtest)

if(extraB)call specialbound(qt,ixMmin1,ixMmin2,ixMmin3,ixMmax1,ixMmax2,&
   ixMmax3,0,0,w)
if(smallfix)call keeppositive(ixMmin1,ixMmin2,ixMmin3,ixMmax1,ixMmax2,ixMmax3,&
   w)

     
! Get boundaries from other PE-s
call mpibound(nw,w)


iB=0
do
   iB=iB+1
   if(oktest)write(*,*)'iB  :',iB
   if(iB>nB) exit
   idim=idimB(iB)
   ! Only boundary segments in the required direction(s) are filled in
   if(idimmin>idim.or.idimmax<idim) cycle

   ixmin1=ixBmin(1,iB);ixmin2=ixBmin(2,iB);ixmin3=ixBmin(3,iB)
   ixmax1=ixBmax(1,iB);ixmax2=ixBmax(2,iB);ixmax3=ixBmax(3,iB);

   ! The possibly shifted coordinates parallel to the boundary layer 
   ! are defined by the PAIR of the periodic/overlapping boundary.
   ! Put the location of the source of information into ixpairL.
   if(ipairB(iB)>0)then
      ixpairmin1=ixBmin(1,ipairB(iB));ixpairmin2=ixBmin(2,ipairB(iB))
      ixpairmin3=ixBmin(3,ipairB(iB));ixpairmax1=ixBmax(1,ipairB(iB))
      ixpairmax2=ixBmax(2,ipairB(iB));ixpairmax3=ixBmax(3,ipairB(iB));
      select case(idim)
      case(1)
         if(upperB(iB))then
            ixpairmin1=ixpairmin1+dixBmin1;ixpairmax1=ixpairmax1+dixBmax1;
         else
            ixpairmin1=ixpairmin1-dixBmin1;ixpairmax1=ixpairmax1-dixBmax1;
         endif
      
      case(2)
         if(upperB(iB))then
            ixpairmin2=ixpairmin2+dixBmin2;ixpairmax2=ixpairmax2+dixBmax2;
         else
            ixpairmin2=ixpairmin2-dixBmin2;ixpairmax2=ixpairmax2-dixBmax2;
         endif
      
      case(3)
         if(upperB(iB))then
            ixpairmin3=ixpairmin3+dixBmin3;ixpairmax3=ixpairmax3+dixBmax3;
         else
            ixpairmin3=ixpairmin3-dixBmin3;ixpairmax3=ixpairmax3-dixBmax3;
         endif
      
      end select
   endif

   do iw= iwmin,iwmax
      if(oktest)write(*,*)'  iw:',iw
      if(oktest)write(*,*)'typeB(iw,iB):',typeB(iw,iB)
      select case (typeB(iw,iB))



      case('contCD4')
         select case(idim)
         case(1)

            if(upperB(iB))then

      	    
               ixe=ixmin1-2
	       ixee=ixmin1-3

               !HPF$ INDEPENDENT
               ix= ixmin1
	          w(ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)=w(ixe,ixmin2:ixmax2,&
   ixmin3:ixmax3,iw)

               ix= ixmax1
		  w(ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)=w(ixee,ixmin2:ixmax2,ixmin3:ixmax3,&
   iw)
	       
            else
               ixe=ixmax1+3
	       ixee=ixmax1+2

               !HPF$ INDEPENDENT
               ix= ixmin1
                  w(ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)=w(ixe,ixmin2:ixmax2,&
                     ixmin3:ixmax3,iw)

               ix= ixmax1
		  w(ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)=w(ixee,ixmin2:ixmax2,ixmin3:ixmax3,&
   iw)
	       
            endif


             case(2)

            if(upperB(iB))then

      	    
               ixe=ixmin2-2
	       ixee=ixmin2-3

               !HPF$ INDEPENDENT
               ix= ixmin2
	          w(ixmin1:ixmax1,ix,ixmin3:ixmax3,iw)=w(ixmin1:ixmax1,ixe,&
   ixmin3:ixmax3,iw)

               ix= ixmax2
		  w(ixmin1:ixmax1,ix,ixmin3:ixmax3,iw)=w(ixmin1:ixmax1,ixee,ixmin3:ixmax3,&
   iw)
	       
            else
               ixe=ixmax2+3
	       ixee=ixmax2+2

               !HPF$ INDEPENDENT
               ix= ixmin2
                  w(ixmin1:ixmax1,ix,ixmin3:ixmax3,iw)=w(ixmin1:ixmax1,ixe,&
                     ixmin3:ixmax3,iw)

               ix= ixmax2
		  w(ixmin1:ixmax1,ix,ixmin3:ixmax3,iw)=w(ixmin1:ixmax1,ixee,ixmin3:ixmax3,&
   iw)
	       
            endif


             case(3)

            if(upperB(iB))then

      	    
               ixe=ixmin3-2
	       ixee=ixmin3-3

               !HPF$ INDEPENDENT
               ix= ixmin3
	          w(ixmin1:ixmax1,ixmin2:ixmax2,ix,iw)=w(ixmin1:ixmax1,ixmin2:ixmax2,&
   ixe,iw)

               ix= ixmax3
		  w(ixmin1:ixmax1,ixmin2:ixmax2,ix,iw)=w(ixmin1:ixmax1,ixmin2:ixmax2,ixee,&
   iw)
	       
            else
               ixe=ixmax3+3
	       ixee=ixmax3+2

               !HPF$ INDEPENDENT
               ix= ixmin3
                  w(ixmin1:ixmax1,ixmin2:ixmax2,ix,iw)=w(ixmin1:ixmax1,&
                     ixmin2:ixmax2,ixe,iw)

               ix= ixmax3
		  w(ixmin1:ixmax1,ixmin2:ixmax2,ix,iw)=w(ixmin1:ixmax1,ixmin2:ixmax2,ixee,&
   iw)
	       
            endif


         
         end select

      case('zero')
         select case(idim)
         case(1)

            if(upperB(iB))then
              
	      call primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,w)
	      
               ixe=ixmin1-2
	       ixee=ixmin1-3

	       ixbmax1=ixmax1
	       ixbmin1=ixee
	       
	                     
               !HPF$ INDEPENDENT
               ix= ixmin1
	       
                  w(ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)=w(ixe,ixmin2:ixmax2,&
                     ixmin3:ixmax3,iw)

               ix= ixmax1
		  w(ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)=w(ixee,ixmin2:ixmax2,ixmin3:ixmax3,&
   iw)
	       
	      call conserve(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,w)	         
	       
           else
	       
	      call primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,w)
               
               ixe=ixmax1+3
	       ixee=ixmax1+2
	       
       
               !HPF$ INDEPENDENT
               ix= ixmin1
                  w(ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)=w(ixe,ixmin2:ixmax2,&
                     ixmin3:ixmax3,iw)

               ix= ixmax1
		  w(ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)=w(ixee,ixmin2:ixmax2,ixmin3:ixmax3,&
   iw)
	

	      call conserve(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
   w)	         	       
	       
            endif
	
   	

             case(2)

            if(upperB(iB))then
              
	      call primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,w)
	      
               ixe=ixmin2-2
	       ixee=ixmin2-3

	       ixbmax2=ixmax2
	       ixbmin2=ixee
	       
	                     
               !HPF$ INDEPENDENT
               ix= ixmin2
	       
                  w(ixmin1:ixmax1,ix,ixmin3:ixmax3,iw)=w(ixmin1:ixmax1,ixe,&
                     ixmin3:ixmax3,iw)

               ix= ixmax2
		  w(ixmin1:ixmax1,ix,ixmin3:ixmax3,iw)=w(ixmin1:ixmax1,ixee,ixmin3:ixmax3,&
   iw)
	       
	      call conserve(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,w)	         
	       
           else
	       
	      call primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,w)
               
               ixe=ixmax2+3
	       ixee=ixmax2+2
	       
       
               !HPF$ INDEPENDENT
               ix= ixmin2
                  w(ixmin1:ixmax1,ix,ixmin3:ixmax3,iw)=w(ixmin1:ixmax1,ixe,&
                     ixmin3:ixmax3,iw)

               ix= ixmax2
		  w(ixmin1:ixmax1,ix,ixmin3:ixmax3,iw)=w(ixmin1:ixmax1,ixee,ixmin3:ixmax3,&
   iw)
	

	      call conserve(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
   w)	         	       
	       
            endif
	
   	

             case(3)

            if(upperB(iB))then
              
	      call primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,w)
	      
               ixe=ixmin3-2
	       ixee=ixmin3-3

	       ixbmax3=ixmax3
	       ixbmin3=ixee
	       
	                     
               !HPF$ INDEPENDENT
               ix= ixmin3
	       
                  w(ixmin1:ixmax1,ixmin2:ixmax2,ix,iw)=w(ixmin1:ixmax1,&
                     ixmin2:ixmax2,ixe,iw)

               ix= ixmax3
		  w(ixmin1:ixmax1,ixmin2:ixmax2,ix,iw)=w(ixmin1:ixmax1,ixmin2:ixmax2,ixee,&
   iw)
	       
	      call conserve(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,w)	         
	       
           else
	       
	      call primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,w)
               
               ixe=ixmax3+3
	       ixee=ixmax3+2
	       
       
               !HPF$ INDEPENDENT
               ix= ixmin3
                  w(ixmin1:ixmax1,ixmin2:ixmax2,ix,iw)=w(ixmin1:ixmax1,&
                     ixmin2:ixmax2,ixe,iw)

               ix= ixmax3
		  w(ixmin1:ixmax1,ixmin2:ixmax2,ix,iw)=w(ixmin1:ixmax1,ixmin2:ixmax2,ixee,&
   iw)
	

	      call conserve(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
   w)	         	       
	       
            endif
	
   	

         
          
          
          
         end select
     


      case('cont','fixed')
         ! For 'cont' copy w at the edge into the whole boundary region.
         ! Fot 'fixed' copy w first, later use values stored in fixB.
         ! For fullgridini=T store the w values read from the file in fixB.
         select case(idim)
         case(1)
            if(upperB(iB))then
               ixe=ixmin1-1
            else
               ixe=ixmax1+1
            endif
            if(fixedB(iw,iB))then
               !HPF$ INDEPENDENT
               do ix1=ixmin1,ixmax1
         do ix2=ixmin2,ixmax2
         do ix3=ixmin3,ixmax3
                  w(ix1,ix2,ix3,iw)=fixB1(ix1-ixe,ix2,ix3,iw)
               enddo
         enddo
         enddo 
            else if(typeB(iw,iB)=='cont' .or. .not.fullgridini) then
               !HPF$ INDEPENDENT
                do ix= ixmin1,ixmax1
                  w(ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)=w(ixe,ixmin2:ixmax2,&
                     ixmin3:ixmax3,iw)
               end do
            end if 
         case(2)
            if(upperB(iB))then
               ixe=ixmin2-1
            else
               ixe=ixmax2+1
            endif
            if(fixedB(iw,iB))then
               !HPF$ INDEPENDENT
               do ix1=ixmin1,ixmax1
         do ix2=ixmin2,ixmax2
         do ix3=ixmin3,ixmax3
                  w(ix1,ix2,ix3,iw)=fixB2(ix1,ix2-ixe,ix3,iw)
               enddo
         enddo
         enddo 
            else if(typeB(iw,iB)=='cont' .or. .not.fullgridini) then
               !HPF$ INDEPENDENT
                do ix= ixmin2,ixmax2
                  w(ixmin1:ixmax1,ix,ixmin3:ixmax3,iw)=w(ixmin1:ixmax1,ixe,&
                     ixmin3:ixmax3,iw)
               end do
            end if 
         case(3)
            if(upperB(iB))then
               ixe=ixmin3-1
            else
               ixe=ixmax3+1
            endif
            if(fixedB(iw,iB))then
               !HPF$ INDEPENDENT
               do ix1=ixmin1,ixmax1
         do ix2=ixmin2,ixmax2
         do ix3=ixmin3,ixmax3
                  w(ix1,ix2,ix3,iw)=fixB3(ix1,ix2,ix3-ixe,iw)
               enddo
         enddo
         enddo 
            else if(typeB(iw,iB)=='cont' .or. .not.fullgridini) then
               !HPF$ INDEPENDENT
                do ix= ixmin3,ixmax3
                  w(ixmin1:ixmax1,ixmin2:ixmax2,ix,iw)=w(ixmin1:ixmax1,&
                     ixmin2:ixmax2,ixe,iw)
               end do
            end if 
         end select
      case ('cont1','fixed1','grad1')
         ! First order extrapolation from edge and inner edge to the boundary
         ! 'cont1'  extrapolates in every time step, can be numericly unstable.
         ! 'fixed1' extrapolates first, stores VALUES into fixB, then restores.
         ! 'grad1'  extrapolates first, stores DIFFERENCES into fixB, then 
         !          adds the stored differences to the current edge value.
         select case(idim)
         case(1)
            if(upperB(iB))then
               ixe=ixmin1-1; ixf=ixe-1
            else
               ixe=ixmax1+1; ixf=ixe+1
            endif
            if(fixedB(iw,iB))then
               if(typeB(iw,iB)=='grad1')then
                  !HPF$ INDEPENDENT
                  do ix1=ixmin1,ixmax1
         do ix2=ixmin2,ixmax2
         do ix3=ixmin3,ixmax3
                     w(ix1,ix2,ix3,iw)=fixB1(ix1-ixe,ix2,ix3,iw)&
                        +w(ixe,ix2,ix3,iw)
                  enddo
         enddo
         enddo
               else
                  !HPF$ INDEPENDENT
                  do ix1=ixmin1,ixmax1
         do ix2=ixmin2,ixmax2
         do ix3=ixmin3,ixmax3
                     w(ix1,ix2,ix3,iw)=fixB1(ix1-ixe,ix2,ix3,iw)
                  enddo
         enddo
         enddo
               endif
            else if(typeB(iw,iB)=='cont1'.or. .not.fullgridini)then !HPF_
            !HPF_ endif
            !HPF_ if(.not.fixedB(iw,iB).and.&
            !HPF_    (typeB(iw,iB)=='cont1'.or. .not.fullgridini))then
               !HPF$ INDEPENDENT
               do ix= ixmin1,ixmax1
                  w(ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)=&
                     (abs(ix-ixe)+1)*w(ixe,ixmin2:ixmax2,ixmin3:ixmax3,iw)- &
                      abs(ix-ixe)   *w(ixf,ixmin2:ixmax2,ixmin3:ixmax3,iw)  
               end do
            end if 
         case(2)
            if(upperB(iB))then
               ixe=ixmin2-1; ixf=ixe-1
            else
               ixe=ixmax2+1; ixf=ixe+1
            endif
            if(fixedB(iw,iB))then
               if(typeB(iw,iB)=='grad1')then
                  !HPF$ INDEPENDENT
                  do ix1=ixmin1,ixmax1
         do ix2=ixmin2,ixmax2
         do ix3=ixmin3,ixmax3
                     w(ix1,ix2,ix3,iw)=fixB2(ix1,ix2-ixe,ix3,iw)&
                        +w(ix1,ixe,ix3,iw)
                  enddo
         enddo
         enddo
               else
                  !HPF$ INDEPENDENT
                  do ix1=ixmin1,ixmax1
         do ix2=ixmin2,ixmax2
         do ix3=ixmin3,ixmax3
                     w(ix1,ix2,ix3,iw)=fixB2(ix1,ix2-ixe,ix3,iw)
                  enddo
         enddo
         enddo
               endif
            else if(typeB(iw,iB)=='cont1'.or. .not.fullgridini)then !HPF_
            !HPF_ endif
            !HPF_ if(.not.fixedB(iw,iB).and.&
            !HPF_    (typeB(iw,iB)=='cont1'.or. .not.fullgridini))then
               !HPF$ INDEPENDENT
               do ix= ixmin2,ixmax2
                  w(ixmin1:ixmax1,ix,ixmin3:ixmax3,iw)=&
                     (abs(ix-ixe)+1)*w(ixmin1:ixmax1,ixe,ixmin3:ixmax3,iw)- &
                      abs(ix-ixe)   *w(ixmin1:ixmax1,ixf,ixmin3:ixmax3,iw)  
               end do
            end if 
         case(3)
            if(upperB(iB))then
               ixe=ixmin3-1; ixf=ixe-1
            else
               ixe=ixmax3+1; ixf=ixe+1
            endif
            if(fixedB(iw,iB))then
               if(typeB(iw,iB)=='grad1')then
                  !HPF$ INDEPENDENT
                  do ix1=ixmin1,ixmax1
         do ix2=ixmin2,ixmax2
         do ix3=ixmin3,ixmax3
                     w(ix1,ix2,ix3,iw)=fixB3(ix1,ix2,ix3-ixe,iw)&
                        +w(ix1,ix2,ixe,iw)
                  enddo
         enddo
         enddo
               else
                  !HPF$ INDEPENDENT
                  do ix1=ixmin1,ixmax1
         do ix2=ixmin2,ixmax2
         do ix3=ixmin3,ixmax3
                     w(ix1,ix2,ix3,iw)=fixB3(ix1,ix2,ix3-ixe,iw)
                  enddo
         enddo
         enddo
               endif
            else if(typeB(iw,iB)=='cont1'.or. .not.fullgridini)then !HPF_
            !HPF_ endif
            !HPF_ if(.not.fixedB(iw,iB).and.&
            !HPF_    (typeB(iw,iB)=='cont1'.or. .not.fullgridini))then
               !HPF$ INDEPENDENT
               do ix= ixmin3,ixmax3
                  w(ixmin1:ixmax1,ixmin2:ixmax2,ix,iw)=&
                     (abs(ix-ixe)+1)*w(ixmin1:ixmax1,ixmin2:ixmax2,ixe,iw)- &
                      abs(ix-ixe)   *w(ixmin1:ixmax1,ixmin2:ixmax2,ixf,iw)  
               end do
            end if 
         end select
      case('periodic')
         ! Update boundary by translation of w by width of mesh (and shift)
         w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,iw)=w&
            (ixpairmin1:ixpairmax1,ixpairmin2:ixpairmax2,&
            ixpairmin3:ixpairmax3,iw)
      case('symm','symm0','asymm')
         ! Reflect w into the boundary region, multiply by -1 for "asymm"
         ! In generalized coordinates take into account the other vector
         ! components for vector variables. The symmetry of the transverse 
         ! component is based on typeB for the jdim=idim+1 -th component.
         ! ixe is used for the reflection, normal vectors are taken at ixf+1/2
         if(gencoord.and.vectoriw(iw)>=0)then
            ! Determine direction of vector component, symmetry coefficients
            ! for normal and transverse vector components
            iwv=vectoriw(iw); jdim=idim+1-(idim/ndim)*ndim 
            coeffnormal=1; if(typeB(iwv+idim,iB)=='asymm') coeffnormal=-1
            coefftransv=1; if(typeB(iwv+jdim,iB)=='asymm') coefftransv=-1
         endif
         select case(idim)
         case(1)
            if(upperB(iB))then
               ixe=2*ixmin1-1; ixf=ixmin1-1
            else
               ixe=2*ixmax1+1; ixf=ixmax1
            endif
            if(gencoord.and.vectoriw(iw)>=0)then
               !HPF$ INDEPENDENT
               do ix= ixmin1,ixmax1
                  w(ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)=zero
                  do jdim=1,ndim
                     w(ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)=w(ix,ixmin2:ixmax2,&
                        ixmin3:ixmax3,iw)+&
                        normalC(ixf,ixmin2:ixmax2,ixmin3:ixmax3,idim,jdim)&
                           *w(ixe-ix,ixmin2:ixmax2,ixmin3:ixmax3,iwv+jdim)
                  end do
                  w(ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)=w(ix,ixmin2:ixmax2,&
                     ixmin3:ixmax3,iw)*&
                     normalC(ixf,ixmin2:ixmax2,ixmin3:ixmax3,idim,iw&
                        -iwv)*(coeffnormal-coefftransv)&
                     +w(ixe-ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)*coefftransv
               end do
            else
               !HPF$ INDEPENDENT
               do ix= ixmin1,ixmax1
                  w(ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)=w(ixe&
                     -ix,ixmin2:ixmax2,ixmin3:ixmax3,iw)
               end do
               if(typeB(iw,iB)=='asymm') w(ixmin1:ixmax1,ixmin2:ixmax2,&
                  ixmin3:ixmax3,iw)=-w(ixmin1:ixmax1,ixmin2:ixmax2,&
                  ixmin3:ixmax3,iw)
            endif  
         
         case(2)
            if(upperB(iB))then
               ixe=2*ixmin2-1; ixf=ixmin2-1
            else
               ixe=2*ixmax2+1; ixf=ixmax2
            endif
            if(gencoord.and.vectoriw(iw)>=0)then
               !HPF$ INDEPENDENT
               do ix= ixmin2,ixmax2
                  w(ixmin1:ixmax1,ix,ixmin3:ixmax3,iw)=zero
                  do jdim=1,ndim
                     w(ixmin1:ixmax1,ix,ixmin3:ixmax3,iw)=w(ixmin1:ixmax1,ix,&
                        ixmin3:ixmax3,iw)+&
                        normalC(ixmin1:ixmax1,ixf,ixmin3:ixmax3,idim,jdim)&
                           *w(ixmin1:ixmax1,ixe-ix,ixmin3:ixmax3,iwv+jdim)
                  end do
                  w(ixmin1:ixmax1,ix,ixmin3:ixmax3,iw)=w(ixmin1:ixmax1,ix,&
                     ixmin3:ixmax3,iw)*&
                     normalC(ixmin1:ixmax1,ixf,ixmin3:ixmax3,idim,iw&
                        -iwv)*(coeffnormal-coefftransv)&
                     +w(ixmin1:ixmax1,ixe-ix,ixmin3:ixmax3,iw)*coefftransv
               end do
            else
               !HPF$ INDEPENDENT
               do ix= ixmin2,ixmax2
                  w(ixmin1:ixmax1,ix,ixmin3:ixmax3,iw)=w(ixmin1:ixmax1,ixe&
                     -ix,ixmin3:ixmax3,iw)
               end do
               if(typeB(iw,iB)=='asymm') w(ixmin1:ixmax1,ixmin2:ixmax2,&
                  ixmin3:ixmax3,iw)=-w(ixmin1:ixmax1,ixmin2:ixmax2,&
                  ixmin3:ixmax3,iw)
            endif  
         
         case(3)
            if(upperB(iB))then
               ixe=2*ixmin3-1; ixf=ixmin3-1
            else
               ixe=2*ixmax3+1; ixf=ixmax3
            endif
            if(gencoord.and.vectoriw(iw)>=0)then
               !HPF$ INDEPENDENT
               do ix= ixmin3,ixmax3
                  w(ixmin1:ixmax1,ixmin2:ixmax2,ix,iw)=zero
                  do jdim=1,ndim
                     w(ixmin1:ixmax1,ixmin2:ixmax2,ix,iw)=w(ixmin1:ixmax1,&
                        ixmin2:ixmax2,ix,iw)+&
                        normalC(ixmin1:ixmax1,ixmin2:ixmax2,ixf,idim,jdim)&
                           *w(ixmin1:ixmax1,ixmin2:ixmax2,ixe-ix,iwv+jdim)
                  end do
                  w(ixmin1:ixmax1,ixmin2:ixmax2,ix,iw)=w(ixmin1:ixmax1,&
                     ixmin2:ixmax2,ix,iw)*&
                     normalC(ixmin1:ixmax1,ixmin2:ixmax2,ixf,idim,iw&
                        -iwv)*(coeffnormal-coefftransv)&
                     +w(ixmin1:ixmax1,ixmin2:ixmax2,ixe-ix,iw)*coefftransv
               end do
            else
               !HPF$ INDEPENDENT
               do ix= ixmin3,ixmax3
                  w(ixmin1:ixmax1,ixmin2:ixmax2,ix,iw)=w(ixmin1:ixmax1,&
                     ixmin2:ixmax2,ixe-ix,iw)
               end do
               if(typeB(iw,iB)=='asymm') w(ixmin1:ixmax1,ixmin2:ixmax2,&
                  ixmin3:ixmax3,iw)=-w(ixmin1:ixmax1,ixmin2:ixmax2,&
                  ixmin3:ixmax3,iw)
            endif  
         
         end select
      case('special')
            ! Skip special now, we do it after normal boundary type variables
            !HPF_ if(.false.)write(*,*)'Avoiding xlhpf90 compiler bug'
           
      case('mpi','mpiperiod')
            ! This boundary is handled by MPI 
      case default
         write(uniterr,*)'Error in GetBoundary: boundary type(', iw,iB,')=',&
            typeB(iw,iB),' is not implemented!'
         call die(' ')
      end select ! typeB(iw,iB)
   end do ! next iw
   ! Do special boundaries
   do iw= iwmin,iwmax
      if(oktest)write(*,*)'special, iw:',iw
      if(typeB(iw,iB).eq.'special')call specialbound(qt,ixmin1,ixmin2,ixmin3,&
         ixmax1,ixmax2,ixmax3,iw,iB,w)
   end do ! next iw
end do ! next iB

! Fixed boundaries (fixed,fixed1) or fixed gradients (grad1) are stored into
! fixB after all boundaries have been updated.
! This needs to be done in the very first time step only.
if(.not.initialized)then
   initialized=.true.
   do iB= 1,nB
      ixmin1=ixBmin(1,iB);ixmin2=ixBmin(2,iB);ixmin3=ixBmin(3,iB)
      ixmax1=ixBmax(1,iB);ixmax2=ixBmax(2,iB);ixmax3=ixBmax(3,iB);
      do iw= iwmin,iwmax
         if( (typeB(iw,iB)=='fixed'.or.typeB(iw,iB)=='fixed1'&
            .or.typeB(iw,iB)=='grad1') .and. .not.fixedB(iw,iB))then
            fixedB(iw,iB)=.true.
            select case(idimB(iB))
               case(1)
                  if(upperB(iB))then
                     ixe=ixmin1-1
                  else
                     ixe=ixmax1+1
                  endif
                  if(typeB(iw,iB)=='grad1')then
                     !HPF$ INDEPENDENT
                     do ix1= ixmin1,ixmax1
               do ix2= ixmin2,ixmax2
               do ix3= ixmin3,ixmax3
                      fixB1(ix1-ixe,ix2,ix3,iw)=w(ix1,ix2,ix3,iw)&
                         -w(ixe,ix2,ix3,iw)
                     enddo
               enddo
               enddo
                  else
                     !HPF$ INDEPENDENT
                     do ix1= ixmin1,ixmax1
               do ix2= ixmin2,ixmax2
               do ix3= ixmin3,ixmax3
                      fixB1(ix1-ixe,ix2,ix3,iw)=w(ix1,ix2,ix3,iw)
                     enddo
               enddo
               enddo
                  endif
               
               case(2)
                  if(upperB(iB))then
                     ixe=ixmin2-1
                  else
                     ixe=ixmax2+1
                  endif
                  if(typeB(iw,iB)=='grad1')then
                     !HPF$ INDEPENDENT
                     do ix1= ixmin1,ixmax1
               do ix2= ixmin2,ixmax2
               do ix3= ixmin3,ixmax3
                      fixB2(ix1,ix2-ixe,ix3,iw)=w(ix1,ix2,ix3,iw)&
                         -w(ix1,ixe,ix3,iw)
                     enddo
               enddo
               enddo
                  else
                     !HPF$ INDEPENDENT
                     do ix1= ixmin1,ixmax1
               do ix2= ixmin2,ixmax2
               do ix3= ixmin3,ixmax3
                      fixB2(ix1,ix2-ixe,ix3,iw)=w(ix1,ix2,ix3,iw)
                     enddo
               enddo
               enddo
                  endif
               
               case(3)
                  if(upperB(iB))then
                     ixe=ixmin3-1
                  else
                     ixe=ixmax3+1
                  endif
                  if(typeB(iw,iB)=='grad1')then
                     !HPF$ INDEPENDENT
                     do ix1= ixmin1,ixmax1
               do ix2= ixmin2,ixmax2
               do ix3= ixmin3,ixmax3
                      fixB3(ix1,ix2,ix3-ixe,iw)=w(ix1,ix2,ix3,iw)&
                         -w(ix1,ix2,ixe,iw)
                     enddo
               enddo
               enddo
                  else
                     !HPF$ INDEPENDENT
                     do ix1= ixmin1,ixmax1
               do ix2= ixmin2,ixmax2
               do ix3= ixmin3,ixmax3
                      fixB3(ix1,ix2,ix3-ixe,iw)=w(ix1,ix2,ix3,iw)
                     enddo
               enddo
               enddo
                  endif
               
            end select
         end if 
      end do ! iw
   end do    ! iB
end if

if(oktest)write(*,*)'GetBoundary wnew:',w(ixtest1,ixtest2,ixtest3,iwtest)

return 
end

!=============================================================================
subroutine setnoflux(iw,idim,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,fRC,&
   ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,fLC,ixLmin1,ixLmin2,&
   ixLmin3,ixLmax1,ixLmax2,ixLmax3)

! Set flux in direction idim to zero for variable iw if typeB is 'symm0'
! in a boundary region

include 'vacdef.f'

integer:: iw,idim,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,ixLmin1,ixLmin2,&
   ixLmin3,ixLmax1,ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
   ixRmax3
double precision:: fRC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
    fLC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
integer:: iB,ixe,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3

!-----------------------------------------------------------------------------

oktest=index(teststr,'setnoflux')>=1

do iB=1,nB
   if(typeB(iw,iB)=='symm0'.and.idimB(iB)==idim)then
      ixBmin1=ixBmin(1,iB);ixBmin2=ixBmin(2,iB);ixBmin3=ixBmin(3,iB)
      ixBmax1=ixBmax(1,iB);ixBmax2=ixBmax(2,iB);ixBmax3=ixBmax(3,iB);
      ! Calculate edge index and set the appropriate flux to zero
      select case(idim)
      case(1)
      if(upperB(iB))then
         ixe=ixBmin1-1+ixRmin1-ixmin1
         if(oktest)write(*,*)'Setnoflux it,idim,iw,ixe:', &
             it,idim,iw,ixe,fRC(ixe,ixtest2,ixtest3)
         fRC(ixe,ixBmin2:ixBmax2,ixBmin3:ixBmax3)=zero
      else
         ixe=ixBmax1+1+ixLmin1-ixmin1
         if(oktest)write(*,*)'Setnoflux it,idim,iw,ixe:',&
             it,idim,iw,ixe,fLC(ixe,ixtest2,ixtest3)
         fLC(ixe,ixBmin2:ixBmax2,ixBmin3:ixBmax3)=zero
      endif
      
      case(2)
      if(upperB(iB))then
         ixe=ixBmin2-1+ixRmin2-ixmin2
         if(oktest)write(*,*)'Setnoflux it,idim,iw,ixe:', &
             it,idim,iw,ixe,fRC(ixtest1,ixe,ixtest3)
         fRC(ixBmin1:ixBmax1,ixe,ixBmin3:ixBmax3)=zero
      else
         ixe=ixBmax2+1+ixLmin2-ixmin2
         if(oktest)write(*,*)'Setnoflux it,idim,iw,ixe:',&
             it,idim,iw,ixe,fLC(ixtest1,ixe,ixtest3)
         fLC(ixBmin1:ixBmax1,ixe,ixBmin3:ixBmax3)=zero
      endif
      
      case(3)
      if(upperB(iB))then
         ixe=ixBmin3-1+ixRmin3-ixmin3
         if(oktest)write(*,*)'Setnoflux it,idim,iw,ixe:', &
             it,idim,iw,ixe,fRC(ixtest1,ixtest2,ixe)
         fRC(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixe)=zero
      else
         ixe=ixBmax3+1+ixLmin3-ixmin3
         if(oktest)write(*,*)'Setnoflux it,idim,iw,ixe:',&
             it,idim,iw,ixe,fLC(ixtest1,ixtest2,ixe)
         fLC(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixe)=zero
      endif
      
      end select
   endif
enddo

return
end

!=============================================================================
subroutine gridsetup1

! Cartesian or polar grid. Determine x at the boundaries.
! Determine often needed combinations of x, such as dx or dvolume.
! Determine variables for axial symmetry
!
! ixe          - edge coordinate of the grid touching the boundary region
! ixf          - coordinate inside of ixe
! qx           - x with an extended index range for calculation of dx

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,hxmin1,hxmin2,hxmin3,&
   hxmax1,hxmax2,hxmax3,jxmin1,jxmin2,jxmin3,jxmax1,jxmax2,jxmax3
integer:: ix,ixe,ixf,idim,jdim
double precision:: qx(IXGlo1-1:IXGhi1+1,IXGlo2-1:IXGhi2+1,IXGlo3-1:IXGhi3&
   +1,ndim)
double precision:: r(IXGLO1-1:IXGHI1+1),rC(IXGLO1-1:IXGHI1+1)

!-----------------------------------------------------------------------------

oktest=index(teststr,'gridsetup')>=1
if(oktest)write(*,*)'GridSetup1'

! Calculate qx in the boundary layers by linearly extrapolating x from
! the touching edge of the grid (ixe) and the inner edge (ixf).

     
! Fill in ghost cells for x from MPI neighbors
call mpibound(ndim,x)

qx(ixGlo1-1:ixGhi1+1,ixGlo2-1:ixGhi2+1,ixGlo3-1:ixGhi3+1,1:ndim)=zero
qx(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim) &
   = x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim)
do idim=1,ndim
   ixmin1=ixGmin1-1;ixmin2=ixGmin2-1;ixmin3=ixGmin3-1;ixmax1=ixGmax1+1
   ixmax2=ixGmax2+1;ixmax3=ixGmax3+1;
   select case(idim)
      case(1)
         ! First do the upper layers
         ixmax1=ixGmax1+1; ixmin1=ixMmax1+1
         if(fullgridini .or.mpiupperB(1)     ) ixmin1=ixGmax1+1
         ixe=ixmin1-1; ixf=ixe-1
         !!!forall replaced by do for sake of sequentializing (Adaptor)
	 do jdim=1,ndim
            do ix= ixmin1,ixmax1
               qx(ix,ixmin2:ixmax2,ixmin3:ixmax3,jdim)=(1+abs(ixe&
                  -ix))*qx(ixe,ixmin2:ixmax2,ixmin3:ixmax3,jdim)- &
                                     abs(ixe-ix) *qx(ixf,ixmin2:ixmax2,&
                                        ixmin3:ixmax3,jdim)
            end do
         end do
         ! Next do the lower layers
         ixmin1=ixGmin1-1; ixmax1=ixMmin1-1
         if(fullgridini .or.mpilowerB(1)     ) ixmax1=ixGmin1-1
         ixe=ixmax1+1; ixf=ixe+1
	 do jdim=1,ndim
            do ix= ixmin1,ixmax1
                qx(ix,ixmin2:ixmax2,ixmin3:ixmax3,jdim)=(1+abs(ixe&
                   -ix))*qx(ixe,ixmin2:ixmax2,ixmin3:ixmax3,jdim)- &
                                      abs(ixe-ix) *qx(ixf,ixmin2:ixmax2,&
                                         ixmin3:ixmax3,jdim)
            end do
         end do 
      case(2)
         ! First do the upper layers
         ixmax2=ixGmax2+1; ixmin2=ixMmax2+1
         if(fullgridini .or.mpiupperB(2)     ) ixmin2=ixGmax2+1
         ixe=ixmin2-1; ixf=ixe-1
         !!!forall replaced by do for sake of sequentializing (Adaptor)
	 do jdim=1,ndim
            do ix= ixmin2,ixmax2
               qx(ixmin1:ixmax1,ix,ixmin3:ixmax3,jdim)=(1+abs(ixe&
                  -ix))*qx(ixmin1:ixmax1,ixe,ixmin3:ixmax3,jdim)- &
                                     abs(ixe-ix) *qx(ixmin1:ixmax1,ixf,&
                                        ixmin3:ixmax3,jdim)
            end do
         end do
         ! Next do the lower layers
         ixmin2=ixGmin2-1; ixmax2=ixMmin2-1
         if(fullgridini .or.mpilowerB(2)     ) ixmax2=ixGmin2-1
         ixe=ixmax2+1; ixf=ixe+1
	 do jdim=1,ndim
            do ix= ixmin2,ixmax2
                qx(ixmin1:ixmax1,ix,ixmin3:ixmax3,jdim)=(1+abs(ixe&
                   -ix))*qx(ixmin1:ixmax1,ixe,ixmin3:ixmax3,jdim)- &
                                      abs(ixe-ix) *qx(ixmin1:ixmax1,ixf,&
                                         ixmin3:ixmax3,jdim)
            end do
         end do 
      case(3)
         ! First do the upper layers
         ixmax3=ixGmax3+1; ixmin3=ixMmax3+1
         if(fullgridini .or.mpiupperB(3)     ) ixmin3=ixGmax3+1
         ixe=ixmin3-1; ixf=ixe-1
         !!!forall replaced by do for sake of sequentializing (Adaptor)
	 do jdim=1,ndim
            do ix= ixmin3,ixmax3
               qx(ixmin1:ixmax1,ixmin2:ixmax2,ix,jdim)=(1+abs(ixe&
                  -ix))*qx(ixmin1:ixmax1,ixmin2:ixmax2,ixe,jdim)- &
                                     abs(ixe-ix) *qx(ixmin1:ixmax1,&
                                        ixmin2:ixmax2,ixf,jdim)
            end do
         end do
         ! Next do the lower layers
         ixmin3=ixGmin3-1; ixmax3=ixMmin3-1
         if(fullgridini .or.mpilowerB(3)     ) ixmax3=ixGmin3-1
         ixe=ixmax3+1; ixf=ixe+1
	 do jdim=1,ndim
            do ix= ixmin3,ixmax3
                qx(ixmin1:ixmax1,ixmin2:ixmax2,ix,jdim)=(1+abs(ixe&
                   -ix))*qx(ixmin1:ixmax1,ixmin2:ixmax2,ixe,jdim)- &
                                      abs(ixe-ix) *qx(ixmin1:ixmax1,&
                                         ixmin2:ixmax2,ixf,jdim)
            end do
         end do 
   end select
enddo

x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim)=qx(ixGmin1:ixGmax1,&
   ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim)

! Loop with ^D instead of idim to avoid an SP2 xlphf90 compiler error
! if qx is distributed. But it should not be
!{
!^D&jx^L=ixG^L+kr(^D,^DD); hx^L=ixG^L-kr(^D,^DD);
!   dx(ixG^S,^D)=half*(qx(jx^S,^D)-qx(hx^S,^D))
!   if(oktest)write(*,*)'dx,qxj,qxh:',dx(ixtest^DD,^D),&
!      qx(ixtest^DD+kr(^D,^DD),^D),qx(ixtest^DD-kr(^D,^DD),^D)
!\}

do idim=1,ndim
   jxmin1=ixGmin1+kr(idim,1);jxmin2=ixGmin2+kr(idim,2)
   jxmin3=ixGmin3+kr(idim,3);jxmax1=ixGmax1+kr(idim,1)
   jxmax2=ixGmax2+kr(idim,2);jxmax3=ixGmax3+kr(idim,3)
   hxmin1=ixGmin1-kr(idim,1);hxmin2=ixGmin2-kr(idim,2)
   hxmin3=ixGmin3-kr(idim,3);hxmax1=ixGmax1-kr(idim,1)
   hxmax2=ixGmax2-kr(idim,2);hxmax3=ixGmax3-kr(idim,3);
   dx(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,idim)=half&
      *(qx(jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3,idim)-qx(hxmin1:hxmax1,&
      hxmin2:hxmax2,hxmin3:hxmax3,idim))
end do

! Calculate geometrical factors for axial symmetry based on Boris FCT.
! Fluxes are multiplied by areaC. The cell volume is proportional to areadx.
! Gradient type source terms are multiplied by areaside = darea/areadx.

if(oktest)write(*,*)'Start calculating geometrical factors'



if(oktest)write(*,*)'Start calculating cell volumes'

! Calculate volume of cells and total volume of mesh
if(typeaxial=='slab')then
   dvolume(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3)&
      = dx(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1)&
      *dx(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,2)&
      *dx(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,3)
else
   forall(ix= ixGmin1:ixGmax1) dvolume(ix,ixGmin2:ixGmax2,ixGmin3:ixGmax3)&
      = areadx(ix)*dx(ix,ixGmin2:ixGmax2,ixGmin3:ixGmax3,2)&
      *dx(ix,ixGmin2:ixGmax2,ixGmin3:ixGmax3,3)
endif
volume=sum(dvolume(ixMmin1:ixMmax1,ixMmin2:ixMmax2,ixMmin3:ixMmax3))
     
! Add up volumes
call mpiallreduce(volume,MPI_SUM)
! Correct volumes in 2nd ghost cells from neighboring processors
call mpibound(1,dvolume)


! For polar grid dx_phi=r*dphi. 
if(polargrid)dx(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,pphi_)&
   =x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,r_)*dx(ixGmin1:ixGmax1,&
   ixGmin2:ixGmax2,ixGmin3:ixGmax3,pphi_)

if(oktest)write(*,*)'Finish GridSetup1'

return 
end

!=============================================================================
 
 subroutine gradient4(realgrad,q,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,&
    idim,gradq)
 include 'vacdef.f'
 logical:: realgrad
 integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim
 double precision:: q(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
    gradq(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
 integer:: kxmin1,kxmin2,kxmin3,kxmax1,kxmax2,kxmax3,jxmin1,jxmin2,jxmin3,&
    jxmax1,jxmax2,jxmax3,hxmin1,hxmin2,hxmin3,hxmax1,hxmax2,hxmax3,gxmin1,&
    gxmin2,gxmin3,gxmax1,gxmax2,gxmax3
 integer:: minx11,minx12,minx13,maxx11,maxx12,maxx13,k
 !-----------------------------------------------------------------------------
 
 !SHIFT
 kxmin1=ixmin1+2*kr(idim,1);kxmin2=ixmin2+2*kr(idim,2)
 kxmin3=ixmin3+2*kr(idim,3);kxmax1=ixmax1+2*kr(idim,1)
 kxmax2=ixmax2+2*kr(idim,2);kxmax3=ixmax3+2*kr(idim,3);
 !SHIFT MORE
 jxmin1=ixmin1+kr(idim,1);jxmin2=ixmin2+kr(idim,2);jxmin3=ixmin3+kr(idim,3)
 jxmax1=ixmax1+kr(idim,1);jxmax2=ixmax2+kr(idim,2);jxmax3=ixmax3+kr(idim,3);
 !SHIFT MORE
 hxmin1=ixmin1-kr(idim,1);hxmin2=ixmin2-kr(idim,2);hxmin3=ixmin3-kr(idim,3)
 hxmax1=ixmax1-kr(idim,1);hxmax2=ixmax2-kr(idim,2);hxmax3=ixmax3-kr(idim,3);
 !SHIFT MORE
 gxmin1=ixmin1-2*kr(idim,1);gxmin2=ixmin2-2*kr(idim,2)
 gxmin3=ixmin3-2*kr(idim,3);gxmax1=ixmax1-2*kr(idim,1)
 gxmax2=ixmax2-2*kr(idim,2);gxmax3=ixmax3-2*kr(idim,3);
 
 !SHIFT BEGIN
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=-(q(kxmin1:kxmax1,&
    kxmin2:kxmax2,kxmin3:kxmax3)-8.D0*(q(jxmin1:jxmax1,jxmin2:jxmax2,&
    jxmin3:jxmax3)-q(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3))&
    -q(gxmin1:gxmax1,gxmin2:gxmax2,gxmin3:gxmax3))/dx(ixmin1:ixmax1,&
    ixmin2:ixmax2,ixmin3:ixmax3,idim)/12.D0
 !SHIFT END
 
 minx11=ixmin1+kr(idim,1);minx12=ixmin2+kr(idim,2);minx13=ixmin3+kr(idim,3);
 maxx11=ixmax1-kr(idim,1);maxx12=ixmax2-kr(idim,2);maxx13=ixmax3-kr(idim,3);
 
 do k=0,1  !left-right bc
 
 if (typeB(1,2*idim-1+k) .ne. 'mpi') then
 if (upperB(2*idim-1+k)) then
 select case(idim)
    case(1)
 gradq(ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
 gradq(maxx11,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
    case(2)
 gradq(ixmin1:ixmax1,ixmax2,ixmin3:ixmax3)=0.d0
 gradq(ixmin1:ixmax1,maxx12,ixmin3:ixmax3)=0.d0
    case(3)
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmax3)=0.d0
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,maxx13)=0.d0

 end select
 else
 select case(idim)
    case(1)
 gradq(ixmin1,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
 gradq(minx11,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
    case(2)
 gradq(ixmin1:ixmax1,ixmin2,ixmin3:ixmax3)=0.d0
 gradq(ixmin1:ixmax1,minx12,ixmin3:ixmax3)=0.d0
    case(3)
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3)=0.d0
 gradq(ixmin1:ixmax1,ixmin2:ixmax2,minx13)=0.d0

 end select
 endif
 endif
 enddo
 
 return
 end


!=============================================================================
subroutine laplace4(q,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,laplaceq)

!!! This subroutine should not use tmp or tmp2

! Calculate 4th order laplace of q within ixL in Cartesian direction idim

!!! We assume uniform Cartesian grid in slab symmetry for now

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
double precision:: q(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
   laplaceq(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)

integer:: idim,kxmin1,kxmin2,kxmin3,kxmax1,kxmax2,kxmax3,jxmin1,jxmin2,jxmin3,&
   jxmax1,jxmax2,jxmax3,hxmin1,hxmin2,hxmin3,hxmax1,hxmax2,hxmax3,gxmin1,&
   gxmin2,gxmin3,gxmax1,gxmax2,gxmax3
!-----------------------------------------------------------------------------

if(gencoord)call die('Error: laplace4 does not work for gen.coords!')
if(typeaxial/='slab')call die&
   ('Error: laplace4 does not work for axial symmetry yet!')

oktest= index(teststr,'lpalace')>=1

laplaceq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=zero

do idim=1,ndim
   !SHIFT
   kxmin1=ixmin1+2*kr(idim,1);kxmin2=ixmin2+2*kr(idim,2)
   kxmin3=ixmin3+2*kr(idim,3);kxmax1=ixmax1+2*kr(idim,1)
   kxmax2=ixmax2+2*kr(idim,2);kxmax3=ixmax3+2*kr(idim,3); 
   !SHIFT MORE
   jxmin1=ixmin1+kr(idim,1);jxmin2=ixmin2+kr(idim,2);jxmin3=ixmin3+kr(idim,3)
   jxmax1=ixmax1+kr(idim,1);jxmax2=ixmax2+kr(idim,2);jxmax3=ixmax3+kr(idim,3)
   
   !SHIFT MORE
   hxmin1=ixmin1-kr(idim,1);hxmin2=ixmin2-kr(idim,2);hxmin3=ixmin3-kr(idim,3)
   hxmax1=ixmax1-kr(idim,1);hxmax2=ixmax2-kr(idim,2);hxmax3=ixmax3-kr(idim,3);
   !SHIFT MORE
   gxmin1=ixmin1-2*kr(idim,1);gxmin2=ixmin2-2*kr(idim,2)
   gxmin3=ixmin3-2*kr(idim,3);gxmax1=ixmax1-2*kr(idim,1)
   gxmax2=ixmax2-2*kr(idim,2);gxmax3=ixmax3-2*kr(idim,3);

   !SHIFT BEGIN
   laplaceq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=laplaceq(ixmin1:ixmax1,&
      ixmin2:ixmax2,ixmin3:ixmax3)+(q(kxmin1:kxmax1,kxmin2:kxmax2,&
      kxmin3:kxmax3)+q(gxmin1:gxmax1,gxmin2:gxmax2,gxmin3:gxmax3)&
      +30*q(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)-16*(q(jxmin1:jxmax1,&
      jxmin2:jxmax2,jxmin3:jxmax3)+q(hxmin1:hxmax1,hxmin2:hxmax2,&
      hxmin3:hxmax3)))/dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idim)**2&
      /12
   !SHIFT END

   if(oktest)write(*,*)'idim,q(kx,jx,ix,hx,gx):',idim,q(ixtest1&
      +2*kr(idim,1),ixtest2+2*kr(idim,2),ixtest3+2*kr(idim,3)),q(ixtest1&
      +kr(idim,1),ixtest2+kr(idim,2),ixtest3+kr(idim,3)),q(ixtest1,ixtest2,&
      ixtest3),q(ixtest1-kr(idim,1),ixtest2-kr(idim,2),ixtest3&
      -kr(idim,3)),q(ixtest1-2*kr(idim,1),ixtest2-2*kr(idim,2),ixtest3&
      -2*kr(idim,3))
enddo

if(oktest)write(*,*)'laplaceq:',laplaceq(ixtest1,ixtest2,ixtest3)

return
end

!=============================================================================
! end module vacgrid
!##############################################################################



