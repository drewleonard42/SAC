!##############################################################################
! module vacphys - mhd

!##############################################################################
! module vacphys.mhd0 - common subroutines for mhd and mhdiso





!=============================================================================
subroutine physini

! Tell VAC which variables are vectors, set default entropy coefficients

include 'vacdef.f'
integer:: il
!-----------------------------------------------------------------------------

iw_vector(1)=m0_; iw_vector(2)=b0_

! The values of the constants are taken from Ryu & Jones ApJ 442, 228
do il=1,nw
   select case(il)
   case(fastRW_,fastLW_,slowRW_,slowLW_)
      entropycoef(il)=0.2
   case(alfvRW_,alfvLW_)
      entropycoef(il)=0.4
   case default
      entropycoef(il)= -one
   end select
end do

return
end

!=============================================================================
subroutine process(count,idimmin,idimmax,w)

! Process w before it is advected in directions idim^LIM, or before save
! count=1 and 2 for first and second (half step) processing during advection
! count=ifile+2 for saving results into the file indexed by ifile

include 'vacdef.f'

integer:: count,idimmin,idimmax
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)

logical:: oktime
double precision:: cputime,time1,timeproc
data timeproc /0.D0/

! The processing should eliminate divergence of B.
!-----------------------------------------------------------------------------

oktest=index(teststr,'process')>=1
oktime=index(teststr,'timeproc')>=1

if(oktest)write(*,*)'Process it,idim^LIM,count',it,idimmin,idimmax,count

if(oktime)time1=cputime()


if(count==0)then
   if(divbconstrain)then
      call die('CT module is OFF: setvac -on=ct; make vac')
   endif
else
   ! Use the projection scheme 
   call die('Poisson module is OFF: setvac -on=poisson;make vac')
endif


if(oktime)then
   time1=cputime()-time1
   timeproc=timeproc+time1
   write(*,*)'Time.Proc:',time1,timeproc
endif

return
end

!=============================================================================
subroutine getdt(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3)

! If resistivity is  not zero, check diffusion time limit for dt

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)
integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
!-----------------------------------------------------------------------------

oktest=index(teststr,'getdt')>=1
if(oktest)write(*,*)'GetDt'

if(eqpar(eta_)==zero)return


       write(*,*)'Error: Resistive MHD module is OFF'
call die('Recompile with setvac -on=resist or set eqpar(eta_)=0')

return
end

!=============================================================================
subroutine getdivb(w,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)

! Calculate div B within ixO

include 'vacdef.f'

integer::          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ixmin1,&
   ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw),&
   divb(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
!-----------------------------------------------------------------------------

oktest=index(teststr,'getdivb')>=1
if(oktest)write(*,*)'getdivb ixO=',ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
   ixOmax3

if(fourthorder)then
   ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmin3=ixOmin3-2;ixmax1=ixOmax1+2
   ixmax2=ixOmax2+2;ixmax3=ixOmax3+2;
else
   ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmin3=ixOmin3-1;ixmax1=ixOmax1+1
   ixmax2=ixOmax2+1;ixmax3=ixOmax3+1;
endif
divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
do idim=1,ndim
   tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=w(ixmin1:ixmax1,&
      ixmin2:ixmax2,ixmin3:ixmax3,b0_+idim)+w(ixmin1:ixmax1,ixmin2:ixmax2,&
      ixmin3:ixmax3,bg0_+idim)

      call gradient4(.false.,tmp,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3,idim,tmp2)

   divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=divb(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      ixOmin3:ixOmax3)
enddo

if(oktest)then
   write(*,*)'divb:',divb(ixtest1,ixtest2,ixtest3)
!   write(*,*)'bx=',w(ixtest1-1:ixtest1+1,ixtest2,b1_)
!   write(*,*)'by=',w(ixtest1,ixtest2-1:ixtest2+1,b2_)
!   write(*,*)'x =',x(ixtest1-1:ixtest1+1,ixtest2,1)
!   write(*,*)'y =',x(ixtest1,ixtest2-1:ixtest2+1,2)
!   write(*,*)'dx=',dx(ixtest1,ixtest2,1)
!   write(*,*)'dy=',dx(ixtest1,ixtest2,2)
endif

return
end

!=============================================================================
subroutine getflux(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,iw,idim,f,&
   transport)

! Calculate non-transport flux f_idim[iw] within ix^L.
! Set transport=.true. if a transport flux should be added

include 'vacdef.f'

integer::          ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,iw,idim
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw),&
   f(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3), fb(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,ixGlo3:ixGhi3)
logical::          transport
!-----------------------------------------------------------------------------

oktest= index(teststr,'getflux')>=1
if(oktest.and.iw==iwtest)write(*,*)'Getflux idim,w:',idim,w(ixtest1,ixtest2,&
   ixtest3,iwtest)

transport=.true.

     
select case(iw)
   case(rho_)
      f(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=w(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3,rhob_)*w(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3,m0_+idim)/(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
         rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_))
      
   case(m1_)
      if(idim==1)then
            
          call getptotal(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,f)
	  call getptotal_bg(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,fb)
          fb(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
	  
          f(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=f(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)+fb(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)-(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b1_)&
             *w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg0_+idim)&
             +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b0_&
             +idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg1_))-&
                           w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b1_)&
                              *w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                              b0_+idim) !-&
			 !w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg1_)    !remove for perturbed
      else
          f(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=-(w(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,b1_)*w(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,bg0_+idim)+w(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,b0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,bg1_))-&
                    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b1_)&
                       *w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b0_+idim) !-&
		 !-w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg1_)  !remove for perturbed
 
      endif 
   case(m2_)
      if(idim==2)then
            
          call getptotal(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,f)
	  call getptotal_bg(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,fb)
          fb(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
	  
          f(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=f(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)+fb(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)-(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b2_)&
             *w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg0_+idim)&
             +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b0_&
             +idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg2_))-&
                           w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b2_)&
                              *w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                              b0_+idim) !-&
			 !w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg2_)    !remove for perturbed
      else
          f(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=-(w(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,b2_)*w(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,bg0_+idim)+w(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,b0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,bg2_))-&
                    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b2_)&
                       *w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b0_+idim) !-&
		 !-w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg2_)  !remove for perturbed
 
      endif 
   case(m3_)
      if(idim==3)then
            
          call getptotal(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,f)
	  call getptotal_bg(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,fb)
          fb(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0
	  
          f(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=f(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)+fb(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)-(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_)&
             *w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg0_+idim)&
             +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b0_&
             +idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg3_))-&
                           w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_)&
                              *w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                              b0_+idim) !-&
			 !w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg3_)    !remove for perturbed
      else
          f(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=-(w(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,b3_)*w(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,bg0_+idim)+w(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,b0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,bg3_))-&
                    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_)&
                       *w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b0_+idim) !-&
		 !-w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg3_)  !remove for perturbed
 
      endif 

   case(e_)
  
      call getptotal(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,f)
      call getptotal_bg(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,fb)      
      fb(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=0.d0      

      f(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(w(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3,m0_+idim)*(f(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3)+fb(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))&
         -w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b0_+idim)&
         *( (w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg1_))&
         *w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m1_)+(w(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3,bg2_))*w(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3,m2_)+(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
         bg3_))*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m3_) )&
         -w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg0_+idim)&
         *( (w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b1_))&
         *w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m1_)+(w(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3,b2_))*w(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3,m2_)+(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
         b3_))*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m3_) ))&
         /(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)&
         +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_))+&        
              !  -w(ix^S,bg0_+idim)*( ^C&(w(ix^S,bg^C_))*w(ix^S,m^C_)+ )/(w(ix^S,rho_)+w(ix^S,rhob_))  ! remove for perturbed
               w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,eb_)&
                  *w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m0_&
                  +idim)/(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)&
                  +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_))&
                  -              w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                  b0_+idim)*( (w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                  b1_))*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m1_)&
                  +(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b2_))&
                  *w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m2_)&
                  +(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_))&
                  *w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m3_) )&
                  /(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)&
                  +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_)) !
              


   case(b1_)
      if(idim==1) then
         f(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)= zero
         transport=.false.
      else
      
         f(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)= -w(ixmin1:ixmax1,&
            ixmin2:ixmax2,ixmin3:ixmax3,m1_)/(w(ixmin1:ixmax1,ixmin2:ixmax2,&
            ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
            rhob_))*(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b0_+idim)&
            +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg0_+idim))+ &
                  w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m0_&
                     +idim)/(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                     rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                     rhob_))*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg1_)
		  
      endif  
   case(b2_)
      if(idim==2) then
         f(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)= zero
         transport=.false.
      else
      
         f(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)= -w(ixmin1:ixmax1,&
            ixmin2:ixmax2,ixmin3:ixmax3,m2_)/(w(ixmin1:ixmax1,ixmin2:ixmax2,&
            ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
            rhob_))*(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b0_+idim)&
            +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg0_+idim))+ &
                  w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m0_&
                     +idim)/(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                     rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                     rhob_))*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg2_)
		  
      endif  
   case(b3_)
      if(idim==3) then
         f(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)= zero
         transport=.false.
      else
      
         f(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)= -w(ixmin1:ixmax1,&
            ixmin2:ixmax2,ixmin3:ixmax3,m3_)/(w(ixmin1:ixmax1,ixmin2:ixmax2,&
            ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
            rhob_))*(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b0_+idim)&
            +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg0_+idim))+ &
                  w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m0_&
                     +idim)/(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                     rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                     rhob_))*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg3_)
		  
      endif  

   case default
      call die('Error in getflux: unknown flow variable')
end select

if(oktest.and.iw==iwtest)write(*,*)'transport,f:',transport,f(ixtest1,ixtest2,&
   ixtest3)

return
end

!=============================================================================
subroutine addsource(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iws,qtC,w,qt,wnew)

! Add sources from resistivity and Powell solver

include 'vacdef.f'

integer::          ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iws(niw_)
double precision:: qdt,qtC,qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw),&
   wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)
!-----------------------------------------------------------------------------

oktest=index(teststr,'addsource')>=1
if(oktest)write(*,*)'Addsource, compactres,divbfix:',compactres,divbfix
if(oktest)write(*,*)'Before adding source:',wnew(ixtest1,ixtest2,ixtest3,&
   iwtest)

! Sources for resistivity in eqs. for e, B1, B2 and B3
if(abs(eqpar(eta_))>smalldouble)then
   
          write(*,*)'Error: Resistive MHD module is OFF'
   call die('Recompile with setvac -on=resist or set eqpar(eta_)=0')
endif


! Sources related to div B in the Powell solver
 if(divbfix) call addsource_divb(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
    ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iws,qtC,w,qt,wnew)

if(oktest)write(*,*)'After adding source:',wnew(ixtest1,ixtest2,ixtest3,&
   iwtest)

return
end

!=============================================================================
subroutine addsource_divb(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iws,qtC,w,qt,wnew)

! Add Powell's divB related sources to wnew within ixO if possible, 
! otherwise shrink ixO

include 'vacdef.f'

integer::          ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iws(niw_),iiw,iw
double precision:: qdt,qtC,qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw),&
   wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)
double precision:: divb(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
!-----------------------------------------------------------------------------

! Calculating div B involves first derivatives
call ensurebound(1,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qtC,w)

! We calculate now div B
call getdivb(w,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)
divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=qdt*divb&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

do iiw=1,iws(niw_); iw=iws(iiw)
   select case(iw)
      case(m1_)
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            -(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_)&
            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,bg1_))&
            *divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      case(m2_)
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            -(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b2_)&
            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,bg2_))&
            *divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      case(m3_)
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            -(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b3_)&
            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,bg3_))&
            *divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
     
      case(b1_)
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,m1_)&
            /(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)&
            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rhob_))&
            *divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      case(b2_)
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,m2_)&
            /(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)&
            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rhob_))&
            *divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      case(b3_)
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,m3_)&
            /(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)&
            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rhob_))&
            *divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
     
      case(e_)
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            =wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)&
            -(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,m1_)&
            *(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_)&
            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,bg1_))&
            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,m2_)&
            *(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b2_)&
            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,bg2_))&
            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,m3_)&
            *(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b3_)&
            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,bg3_)) )&
            /(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)&
            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rhob_))&
            *divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
   end select
end do

return
end


!=============================================================================
! end module vacphys.mhd0
!##############################################################################

!=============================================================================

subroutine keeppositive(ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w)

! Keep pressure and density positive (following D.Ryu)

include 'vacdef.f'

integer::          ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)
logical:: toosmallp
!-----------------------------------------------------------------------------


   ! Where rho is small use vacuum state: rho=vacuumrho, v=0, p=smallp, same B
   where((w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)&
      +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_))<smallrho)
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m1_)=zero
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m2_)=zero
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m3_)=zero;
!!!      ^C&w(ix^S,m^C_)=w(ix^S,m^C_)/w(ix^S,rho_)*vacuumrho;
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)=vacuumrho&
         -w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_)
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,e_)=smallp/(eqpar(gamma_)&
         -one)+half*(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b1_)**2&
         +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b2_)**2&
         +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_)**2)&
         -w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,eb_)
   endwhere


! Calculate pressure without clipping toosmall values (.false.)
call getpthermal(w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,tmp)

toosmallp=any(tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)<max(zero,smallp))

if(toosmallp)then
   nerror(toosmallp_)=nerror(toosmallp_)+1
   if(nerror(toosmallp_)==1)then
      write(*,'(a,i2,a,i7)')'Too small pressure (code=',toosmallp_,') at it=',&
         it
      write(*,*)'Value < smallp: ',minval(tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3)),smallp
!     write(*,*)'Location: ',minloc(tmp(ix^S)) !F77_  
      write(*,*)'Processor:',ipe      
   endif
   if(smallp>zero)w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,e_)&
      =max(tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3),smallp)&
      /(eqpar(gamma_)-1)+half*((w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
      m1_)**2+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m2_)**2&
      +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m3_)**2)&
      /w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)+(w(ixmin1:ixmax1,&
      ixmin2:ixmax2,ixmin3:ixmax3,b1_)**2+w(ixmin1:ixmax1,ixmin2:ixmax2,&
      ixmin3:ixmax3,b2_)**2+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_)&
      **2))-w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,eb_)
endif

return
end

!=============================================================================
! end module vacphys - mhd
!##############################################################################
