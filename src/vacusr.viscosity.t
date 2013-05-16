
!==============================================================================
subroutine addsource_visc(qdt,ixI^L,ixO^L,iws,qtC,w,qt,wnew)

! Add viscosity source to wnew within ixO 

include 'vacdef.f'

integer::          ixI^L,ixO^L,iws(niw_)
double precision:: qdt,qtC,qt,w(ixG^T,nw),wnew(ixG^T,nw)

integer:: ix,ix^L,idim,idir,jdir,iiw,iw
double precision:: tmp2(ixG^T),nushk(ixG^T,ndim)


double precision:: tmprhoL(ixG^T), tmprhoR(ixG^T), tmprhoC(ixG^T)
double precision:: tmpVL(ixG^T), tmpVR(ixG^T), tmpVC(ixG^T)
double precision:: tmpBL(ixG^T), tmpBR(ixG^T), tmpBC(ixG^T)

double precision:: tmpL(ixG^T),tmpR(ixG^T), tmpC(ixG^T)

double precision:: nuL(ixG^T),nuR(ixG^T)

integer:: jx^L,hx^L, hxO^L

double precision:: c_ene,c_shk

integer:: i,j,k,l,m,ii0,ii1,t00

double precision:: sB

!-----------------------------------------------------------------------------

! Calculating viscosity sources 
! involves second derivatives, two extra layers
call ensurebound(2,ixI^L,ixO^L,qtC,w)
ix^L=ixO^L^LADD1;

!sehr wichtig
call setnushk(w,ix^L,nushk)

do idim=1,ndim
      tmp(ixI^S)=w(ixI^S,rho_)
      call setnu(w,rho_,idim,ixO^L,nuR,nuL)      
      call gradient1L(tmp,ix^L,idim,tmp2)
      tmpL(ixI^S)=(nuL(ixI^S)+nushk(ixI^S,idim))*tmp2(ixI^S)	     
      call gradient1R(tmp,ix^L,idim,tmp2)
      tmpR(ixI^S)=(nuR(ixI^S)+nushk(ixI^S,idim))*tmp2(ixI^S)
      wnew(ixI^S,rho_)=wnew(ixI^S,rho_)+(tmpR(ixI^S)-tmpL(ixI^S))/dx(ixI^S,idim)*qdt
enddo
   

do idim=1,ndim
      tmp(ixI^S)=w(ixI^S,e_)-half*((^C&w(ixI^S,b^C_)**2+)+(^C&w(ixI^S,m^C_)**2+)/(w(ixI^S,rho_)+w(ixI^S,rhob_)))
      call setnu(w,173,idim,ixO^L,nuR,nuL)      
      call gradient1L(tmp,ix^L,idim,tmp2)
      tmpL(ixI^S)=(nuL(ixI^S)+nushk(ixI^S,idim))*tmp2(ixI^S)      
      call gradient1R(tmp,ix^L,idim,tmp2)
      tmpR(ixI^S)=(nuR(ixI^S)+nushk(ixI^S,idim))*tmp2(ixI^S)
      wnew(ixI^S,e_)=wnew(ixI^S,e_)+(tmpR(ixI^S)-tmpL(ixI^S))/dx(ixI^S,idim)*qdt
enddo




      tmprhoC(ixI^S)=w(ixI^S,rho_)+w(ixI^S,rhob_)



do k=1,ndim
        jx^L=ix^L+kr(k,^D); 
        hx^L=ix^L-kr(k,^D);
	tmprhoL(ix^S)=((w(ix^S,rho_)+w(ix^S,rhob_))+(w(hx^S,rho_)+w(hx^S,rhob_)))/two
	tmprhoR(ix^S)=((w(jx^S,rho_)+w(jx^S,rhob_))+(w(ix^S,rho_)+w(ix^S,rhob_)))/two

   do l=1,ndim
	call setnu(w,l+m0_,k,ixO^L,nuR,nuL)      
	tmp(ixI^S)=w(ixI^S,m0_+l)/(w(ixI^S,rho_)+w(ixI^S,rhob_))


      do ii1=0,1
        if (ii1 .eq. 0) then
                           i=k
                           ii0=l
                        else
                           i=l
                           ii0=k
        endif



        if (i .eq. k) then 
           tmpVL(ix^S)=(w(ix^S,m0_+ii0)+w(hx^S,m0_+ii0))/two
	   tmpVR(ix^S)=(w(jx^S,m0_+ii0)+w(ix^S,m0_+ii0))/two

           call gradient1L(tmp,ix^L,k,tmp2)
           tmpL(ixI^S)=(nuL(ixI^S)+nushk(ixI^S,k))*tmp2(ixI^S)
           call gradient1R(tmp,ix^L,k,tmp2)
           tmpR(ixI^S)=(nuR(ixI^S)+nushk(ixI^S,k))*tmp2(ixI^S) 

           tmp2(ixI^S)=(tmprhoR(ixI^S)*tmpR(ixI^S)-tmprhoL(ixI^S)*tmpL(ixI^S))/dx(ixI^S,k)/two
 
	   wnew(ixI^S,m0_+ii0)=wnew(ixI^S,m0_+ii0)+tmp2(ixI^S)*qdt

           tmp2(ixI^S)=(tmpVR(ixI^S)*tmpR(ixI^S)-tmpVL(ixI^S)*tmpL(ixI^S))/dx(ixI^S,k)/two

	   wnew(ixI^S,e_)=wnew(ixI^S,e_)+tmp2(ixI^S)*qdt
	endif




	if (i .ne. k) then
           call gradient1(tmp,ix^L,k,tmp2)
           tmp2(ixI^S)=tmp2(ixI^S)*(nuL(ixI^S)+nuR(ixI^S)+two*nushk(ixI^S,k))/two/two

           tmp(ixI^S)=tmprhoC(ixI^S)*tmp2(ixI^S)
           call gradient1(tmp,ix^L,i,tmpC)

	   wnew(ixI^S,m0_+ii0)=wnew(ixI^S,m0_+ii0)+tmpC(ixI^S)*qdt

           tmp(ixI^S)=w(ixI^S,m0_+ii0)*tmp2(ixI^S)
           call gradient1(tmp,ix^L,i,tmpC)

	   wnew(ixI^S,e_)=wnew(ixI^S,e_)+tmpC(ixI^S)*qdt
	endif

     enddo
    enddo
   enddo





do k=1,ndim
 do l=1,ndim

   if (k .ne. l) then

    call setnu(w,b0_+l,k,ixO^L,nuR,nuL)

    do ii1=0,1

      if (ii1 .eq. 0) then
              ii0=k
              m=l
              sB=-1.d0
              j=k
      endif

      if (ii1 .eq. 1) then 
              ii0=l    !ii0 is index B
              m=k      !first derivative
              sB=1.d0  !sign B
              j=l      !first B in energy
      endif



!print*,'k,l,m,j,ii0,ii1=',k,l,m,j,ii0,ii1



      if (m .eq. k) then

           jx^L=ix^L+kr(m,^D); 
           hx^L=ix^L-kr(m,^D);
           tmpBL(ix^S)=(w(ix^S,b0_+j)+w(hx^S,b0_+j))/two
	   tmpBR(ix^S)=(w(jx^S,b0_+j)+w(ix^S,b0_+j))/two

           tmp(ixI^S)=w(ixI^S,b0_+l)

           call gradient1L(tmp,ix^L,k,tmp2)
           tmpL(ixI^S)=(nuL(ixI^S))*tmp2(ixI^S)
           call gradient1R(tmp,ix^L,k,tmp2)
           tmpR(ixI^S)=(nuR(ixI^S))*tmp2(ixI^S) 

           wnew(ixI^S,b0_+ii0)=wnew(ixI^S,b0_+ii0)+sB*(tmpR(ixI^S)-tmpL(ixI^S))/dx(ixI^S,k)*qdt

           wnew(ixI^S,e_)=wnew(ixI^S,e_)+sB*(tmpR(ixI^S)*tmpBR(ixI^S)-tmpL(ixI^S)*tmpBL(ixI^S))/dx(ixI^S,k)*qdt


      endif



      if (m .ne. k) then

           tmp(ixI^S)=w(ixI^S,b0_+l)

           call gradient1(tmp,ix^L,k,tmp2)

           tmp2(ixI^S)=tmp2(ixI^S)*(nuL(ixI^S)+nuR(ixI^S))/two

           call gradient1(tmp2,ix^L,m,tmpC)

           wnew(ixI^S,b0_+ii0)=wnew(ixI^S,b0_+ii0)+sB*tmpC(ixI^S)*qdt

           tmp2(ixI^S)=tmp2(ixI^S)*w(ixI^S,b0_+j)

           call gradient1(tmp2,ix^L,m,tmpC)

           wnew(ixI^S,e_)=wnew(ixI^S,e_)+sB*tmpC(ixI^S)*qdt

      endif


      enddo
   endif
 enddo
enddo




return
end

!=============================================================================
subroutine setnu(w,iw,idim,ix^L,nuR,nuL)

! Set the viscosity coefficient nu within ixO based on w(ixI). 

include 'vacdef.f'

double precision:: w(ixG^T,nw)
double precision:: d1R(^SIDEADO),d1L(^SIDEADO)
double precision:: d3R(^SIDEADO),d3L(^SIDEADO)
double precision:: md3R(ixG^T),md3L(ixG^T)
double precision:: md1R(ixG^T),md1L(ixG^T)
double precision:: nuR(ixG^T),nuL(ixG^T)

double precision:: c_tot, c_hyp,cmax(ixG^T), tmp_nu(ixG^T)
integer:: ix^L,idim, iw
integer:: kx^L,jx^L,hx^L,gx^L,ixFF^L,jxFF^L,hxFF^L
integer:: ix_1,ix_2,ix_3

integer:: ixF^LL,ixF^L,ixY^LL

logical:: new_cmax

double precision:: tmp_nuI(^SIDEADD)

integer:: k,iwc

integer:: ix,ixe

{^IFMPI 

integer :: nmpirequest, mpirequests(2)
integer :: mpistatus(MPI_STATUS_SIZE,2)
common /mpirecv/ nmpirequest,mpirequests,mpistatus


integer:: hpe,jpe

double precision:: tgtbufferR^D(1^D%^LM)
double precision:: tgtbufferL^D(1^D%^LM)
double precision:: srcbufferR^D(1^D%^LM)
double precision:: srcbufferL^D(1^D%^LM)

integer:: n

}

!----------------------------------------------------------------------------

new_cmax=.true.

call getcmax(new_cmax,w,ix^L,idim,cmax)
c_tot=maxval(cmax(ix^S))

{^IFMPI call mpiallreduce(c_tot,MPI_MAX)}



c_hyp=0.4d0 ! 1.4d0 ! 0.6

if (iw.eq.b^D_|.or.) c_hyp=0.2d0 ! 2d0

if (iw .eq. rho_) c_hyp=0.45d0 !5d0

if (iw .eq. 173) c_hyp=0.2d0 !2d0



!if (iw.eq.b^D_|.or.) c_hyp=0.02d0

!if (iw .eq. rho_) c_hyp=0.02d0

!if (iw .eq. 173) c_hyp=0.02d0

        
if (iw .ne. 173) then     
	tmp_nu(ixG^T)=w(ixG^T,iw)
        if (iw.eq.m^D_|.or.) tmp_nu(ixG^T)=w(ixG^T,iw)/(w(ixG^T,rho_)+w(ixG^T,rhob_))
endif

if (iw .eq. 173) tmp_nu(ixG^T)=w(ixG^T,e_)-half*((^C&w(ixG^T,b^C_)**2+)+(^C&w(ixG^T,m^C_)**2+)/(w(ixG^T,rho_)+w(ixG^T,rhob_)))


ixY^LL=ix^L^LADD2;

ixF^LL=ixY^LL+1;

tmp_nuI(ixF^T)=tmp_nu(ixY^T)


{^IFMPI

call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)

n = ^D&(ixFhi^D-ixFlo^D+1)*   

select case(idim)
 { case(^D)

 n=n/(ixFhi^D-ixFlo^D+1)

 }
end select



select case(idim)
{   case(^D)

 if(npe^D>1)then

   nmpirequest =0
   mpirequests(1:2) = MPI_REQUEST_NULL


!source
   srcbufferL^D(1^D%ixF^T)=tmp_nuI(ixFlo^D+4^D%ixF^T) !left, lower

   srcbufferR^D(1^D%ixF^T)=tmp_nuI(ixFhi^D-4^D%ixF^T) !right, upper

   call mpineighbors(^D,hpe,jpe)


   if (mpiupperB(^D)) nmpirequest=nmpirequest+1
   if (mpiupperB(^D)) call MPI_IRECV(tgtbufferR^D(1),n,MPI_DOUBLE_PRECISION, jpe,10*jpe+0,MPI_COMM_WORLD, mpirequests(nmpirequest),ierrmpi)

   if (mpilowerB(^D)) nmpirequest=nmpirequest+1
   if (mpilowerB(^D)) call MPI_IRECV(tgtbufferL^D(1),n,MPI_DOUBLE_PRECISION, hpe,10*hpe+1,MPI_COMM_WORLD, mpirequests(nmpirequest),ierrmpi)

   call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)

   if (mpiupperB(^D)) call MPI_RSEND(srcbufferR^D(1),n,MPI_DOUBLE_PRECISION, jpe,10*ipe+1,MPI_COMM_WORLD,ierrmpi)

   if (mpilowerB(^D)) call MPI_RSEND(srcbufferL^D(1),n,MPI_DOUBLE_PRECISION, hpe,10*ipe+0,MPI_COMM_WORLD,ierrmpi)

   call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)

!target
   tmp_nuI(ixFhi^D+1^D%ixF^T)=tgtbufferR^D(1^D%ixF^T) !right, upper R

   tmp_nuI(ixFlo^D-1^D%ixF^T)=tgtbufferL^D(1^D%ixF^T) !left, lower  L


  endif 
}
end select

call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
}


if (iw .eq. 173) then 
   iwc=e_ 
else 
   iwc=iw
endif

do k=0,1  !left-right bc

if (typeB(iwc,2*idim-1+k) .ne. 'mpi') then
      if (upperB(2*idim-1+k)) then

          select case(idim)
            {   case(^D)
                     tmp_nuI(ixFhi^D+1^D%ixF^T)=tmp_nuI(ixFhi^D-5^D%ixF^T)
             }
          end select

      else

          select case(idim)
            {   case(^D)
                     tmp_nuI(ixFlo^D-1^D%ixF^T)=tmp_nuI(ixFlo^D+5^D%ixF^T)
             }
          end select

      endif
endif

enddo 

        ixF^L=ixF^LL^LSUB1; 

        kx^L=ixF^L+2*kr(idim,^D);  !5:66
        jx^L=ixF^L+kr(idim,^D);  !4:65
        hx^L=ixF^L-kr(idim,^D);  !2:63
        gx^L=ixF^L-2*kr(idim,^D);  !1:62

	ixFF^L=ixF^LL;   !2:65
        jxFF^L=ixF^LL+kr(idim,^D);  !3:66
	hxFF^L=ixF^LL-kr(idim,^D);  !1:64
	
        d3R(ixF^S)=abs(3.d0*(tmp_nuI(jx^S)-tmp_nuI(ixF^S))-(tmp_nuI(kx^S)-tmp_nuI(hx^S))) !3:64
        d1R(ixFF^S)=abs(tmp_nuI(jxFF^S)-tmp_nuI(ixFF^S)) !2:65

	{do ix_^D=ixmin^D,ixmax^D\}    !3:62  +1=4:63

	  md3R(ix_^D)=maxval(d3R(ix_^D+1-kr(idim,^D):ix_^D+1+kr(idim,^D)))
	  md1R(ix_^D)=maxval(d1R(ix_^D+1-2*kr(idim,^D):ix_^D+1+2*kr(idim,^D)))

	{enddo\}

	WHERE (md1R(ix^S).gt.0.d0)
	  nuR(ix^S)=c_tot*c_hyp*md3R(ix^S)/md1R(ix^S)*dx(ix^S,idim)
	ELSEWHERE 
	  nuR(ix^S)=0.d0
	END WHERE
        
        maxviscoef=max(maxval(nuR(ix^S)), maxviscoef)


!************

        d3L(ixF^S)=abs(3.d0*(tmp_nuI(ixF^S)-tmp_nuI(hx^S))-(tmp_nuI(jx^S)-tmp_nuI(gx^S)))
        d1L(ixFF^S)=abs(tmp_nuI(ixFF^S)-tmp_nuI(hxFF^S))    

	{do ix_^D=ixmin^D,ixmax^D\}

	  md3L(ix_^D)=maxval(d3L(ix_^D+1-kr(idim,^D):ix_^D+1+kr(idim,^D)))
	  md1L(ix_^D)=maxval(d1L(ix_^D+1-2*kr(idim,^D):ix_^D+1+2*kr(idim,^D)))

	{enddo\}

	WHERE (md1L(ix^S).gt.0.d0)
	  nuL(ix^S)=c_tot*c_hyp*md3L(ix^S)/md1L(ix^S)*dx(ix^S,idim)
	ELSEWHERE 
	  nuL(ix^S)=0.d0  
	END WHERE

        maxviscoef=max(maxval(nuL(ix^S)), maxviscoef)

        {^IFMPI call mpiallreduce(maxviscoef,MPI_MAX)}

return
end


!=============================================================================
!=============================================================================
subroutine setnushk(w,ix^L,nushk)

include 'vacdef.f'

double precision:: w(ixG^T,nw),tmp2(ixG^T),nushk(ixG^T,ndim)
double precision:: c_shk

double precision:: tmp3(ixG^T)

integer:: ix^L,idim, iw,i

integer:: ix_1,ix_2

do idim=1,ndim
nushk(ix^S,idim)=0.d0
enddo


go to 100
c_shk=0.5d0

tmp3(ix^S)=0.d0

!**************************BEGIN shock viscosity*******************************
      do idim=1,ndim
         tmp(ix^S)=w(ix^S,m0_+idim)/(w(ix^S,rho_)+w(ix^S,rhob_))
         call gradient1(tmp,ix^L,idim,tmp2)
         tmp3(ix^S)=tmp3(ix^S)+tmp2(ix^S)
       enddo
      do idim=1,ndim
        nushk(ix^S,idim)=tmp3(ix^S)*(dx(ix^S,idim)**2.d0)*c_shk
	WHERE (tmp3(ix^S) .ge. 0.d0)
!	  nushk(ix^S,idim)=0.d0
	END WHERE
	nushk(ix^S,idim)=abs(nushk(ix^S,idim))
      enddo
!****************************END shock viscosity*******************************

100 continue


return
end



!=============================================================================
subroutine getdt_visc(w,ix^L)

! Check diffusion time limit for dt < dtdiffpar * dx**2 / (nu/rho)

! Based on Hirsch volume 2, p.631, eq.23.2.17

include 'vacdef.f'

double precision:: w(ixG^T,nw),dtdiff_visc
integer:: ix^L,idim, ix_1,ix_2

integer:: aa

! For spatially varying nu you need a common nu array
 double precision::tmpdt(ixG^T), nuL(ixG^T),nuR(ixG^T), nushk(ixG^T,ndim)
 common/visc/nuL
 common/visc/nuR
!-----------------------------------------------------------------------------

 call setnushk(w,ix^L,nushk)

 dtdiffpar=0.25d0

 do idim=1,ndim
   tmpdt(ix^S)=(maxviscoef+nushk(ix^S,idim))       !/(w(ix^S,rho_)+w(ix^S,rhob_))   ! ~1/dt
   dtdiff_visc=dtdiffpar/maxval(tmpdt(ix^S)/(dx(ix^S,idim)**2))
   {^IFMPI call mpiallreduce(dtdiff_visc,MPI_MIN)}
   dt=min(dt,dtdiff_visc)
 end do
 
 maxviscoef=0.d0

return
end


!***** 2-point central finite difference gradient******

subroutine gradient1(q,ix^L,idim,gradq)
include 'vacdef.f'
integer:: ix^L,idim
double precision:: q(ixG^T),gradq(ixG^T)
integer:: hx^L,kx^L
integer:: minx1^D,maxx1^D,k
!-----------------------------------------------------------------------------

hx^L=ix^L-kr(idim,^D);
kx^L=ix^L+kr(idim,^D);
gradq(ix^S)=(q(kx^S)-q(hx^S))/dx(ix^S,idim)/two

 minx1^D=ixmin^D+kr(idim,^D);
 maxx1^D=ixmax^D-kr(idim,^D);
 
 do k=0,1  !left-right bc
 
 if (typeB(1,2*idim-1+k) .ne. 'mpi') then
 if (upperB(2*idim-1+k)) then
 select case(idim)
 {   case(^D)
 gradq(ixmax^D^D%ix^S)=0.d0
 gradq(maxx1^D^D%ix^S)=0.d0
 }
 end select
 else
 select case(idim)
 {   case(^D)
 gradq(ixmin^D^D%ix^S)=0.d0
 gradq(minx1^D^D%ix^S)=0.d0
 }
 end select
 endif
 endif
 enddo


return
end

!=============================================================================


!*****left upwind forward 2-point non-central finite difference gradient******

subroutine gradient1L(q,ix^L,idim,gradq)
include 'vacdef.f'
integer:: ix^L,idim
double precision:: q(ixG^T),gradq(ixG^T)
integer:: hx^L
integer:: minx1^D,maxx1^D,k
!-----------------------------------------------------------------------------

hx^L=ix^L-kr(idim,^D);
gradq(ix^S)=(q(ix^S)-q(hx^S))/dx(ix^S,idim)

 minx1^D=ixmin^D+kr(idim,^D);
 maxx1^D=ixmax^D-kr(idim,^D);
 
 do k=0,1  !left-right bc
 
 if (typeB(1,2*idim-1+k) .ne. 'mpi') then
 if (upperB(2*idim-1+k)) then
 select case(idim)
 {   case(^D)
 gradq(ixmax^D^D%ix^S)=0.d0
 gradq(maxx1^D^D%ix^S)=0.d0
 }
 end select
 else
 select case(idim)
 {   case(^D)
 gradq(ixmin^D^D%ix^S)=0.d0
 gradq(minx1^D^D%ix^S)=0.d0
 }
 end select
 endif
 endif
 enddo


return
end

!=============================================================================

!*****right upwind forward 2-point non-central finite difference gradient*****

subroutine gradient1R(q,ix^L,idim,gradq)
include 'vacdef.f'
integer:: ix^L,idim
double precision:: q(ixG^T),gradq(ixG^T)
integer:: hx^L
integer:: minx1^D,maxx1^D,k
!-----------------------------------------------------------------------------

hx^L=ix^L+kr(idim,^D);
gradq(ix^S)=(q(hx^S)-q(ix^S))/dx(ix^S,idim)

 minx1^D=ixmin^D+kr(idim,^D);
 maxx1^D=ixmax^D-kr(idim,^D);
 
 do k=0,1  !left-right bc
 
 if (typeB(1,2*idim-1+k) .ne. 'mpi') then
 if (upperB(2*idim-1+k)) then
 select case(idim)
 {   case(^D)
 gradq(ixmax^D^D%ix^S)=0.d0
 gradq(maxx1^D^D%ix^S)=0.d0
 }
 end select
 else
 select case(idim)
 {   case(^D)
 gradq(ixmin^D^D%ix^S)=0.d0
 gradq(minx1^D^D%ix^S)=0.d0
 }
 end select
 endif
 endif
 enddo


return
end
