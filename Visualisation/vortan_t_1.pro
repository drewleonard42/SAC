pro curl,vx,vy,vz,wx,wy,wz,xax,yax

nptx=n_elements(xax)
npty=n_elements(yax)

dvzdy=dblarr(nptx,npty,nptx)
dvydz=dblarr(nptx,npty,nptx)
dvxdz=dblarr(nptx,npty,nptx)
dvzdx=dblarr(nptx,npty,nptx)
dvydx=dblarr(nptx,npty,nptx)
dvxdy=dblarr(nptx,npty,nptx)


for i=0,nptx-1 do begin
 for j=0,nptx-1 do begin
  dvzdy(i,*,j)=deriv(yax,vz(i,*,j))
  dvxdy(i,*,j)=deriv(yax,vx(i,*,j))
 endfor
endfor

for i=0,npty-1 do begin
 for j=0,nptx-1 do begin
  dvydz(j,i,*)=deriv(xax,vy(j,i,*))
  dvxdz(j,i,*)=deriv(xax,vx(j,i,*))
  dvzdx(*,i,j)=deriv(xax,vz(*,i,j))
  dvydx(*,i,j)=deriv(xax,vy(*,i,j))
 endfor
endfor

wx=dvzdy-dvydz
wy=dvxdz-dvzdx
wz=dvydx-dvxdy

end




pro cross,ax,ay,az,bx,by,bz,cx,cy,cz
cx=ay*bz-az*by
cy=az*bx-ax*bz
cz=ax*by-ay*bx
end





pro div,ax,ay,az,res,xax,yax

nptx=n_elements(xax)
npty=n_elements(yax)
dvxdx=dblarr(nptx,npty,nptx)
dvydy=dblarr(nptx,npty,nptx)
dvzdz=dblarr(nptx,npty,nptx)

for i=0,nptx-1 do begin
 for j=0,nptx-1 do begin
  dvydy(i,*,j)=deriv(yax,ay(i,*,j))
 endfor
endfor

for i=0,npty-1 do begin
 for j=0,nptx-1 do begin
  dvzdz(j,i,*)=deriv(xax,az(j,i,*))
  dvxdx(*,i,j)=deriv(xax,ax(*,i,j))
 endfor
endfor

res=dvxdx+dvydy+dvzdz
end


pro grad,field,ax,ay,az,xax,yax
nptx=n_elements(xax)
npty=n_elements(yax)

for i=0,nptx-1 do begin
 for j=0,nptx-1 do begin
  ay(i,*,j)=deriv(yax,field(i,*,j))
 endfor
endfor

for i=0,npty-1 do begin
 for j=0,nptx-1 do begin
  az(j,i,*)=deriv(xax,field(j,i,*))
  ax(*,i,j)=deriv(xax,field(*,i,j))
 endfor
endfor

end




pro t1,wx,wy,wz,vx,vy,vz,t1x,t1y,t1z,xax,yax
nptx=n_elements(xax)
npty=n_elements(yax)
dvxdx=dblarr(nptx,npty,nptx)
dvydx=dblarr(nptx,npty,nptx)
dvzdx=dblarr(nptx,npty,nptx)
dvxdy=dblarr(nptx,npty,nptx)
dvydy=dblarr(nptx,npty,nptx)
dvzdy=dblarr(nptx,npty,nptx)
dvxdz=dblarr(nptx,npty,nptx)
dvydz=dblarr(nptx,npty,nptx)
dvzdz=dblarr(nptx,npty,nptx)

for i=0,nptx-1 do begin
 for j=0,nptx-1 do begin
  dvzdy(i,*,j)=deriv(yax,vz(i,*,j))
  dvydy(i,*,j)=deriv(yax,vy(i,*,j))
  dvxdy(i,*,j)=deriv(yax,vx(i,*,j))  
 endfor
endfor

for i=0,npty-1 do begin
 for j=0,nptx-1 do begin
  dvzdz(j,i,*)=deriv(xax,vz(j,i,*))
  dvxdx(*,i,j)=deriv(xax,vx(*,i,j))
  dvxdz(j,i,*)=deriv(xax,vx(j,i,*))
  dvzdx(*,i,j)=deriv(xax,vz(*,i,j))
  dvydz(j,i,*)=deriv(xax,vy(j,i,*))
  dvydx(*,i,j)=deriv(xax,vy(*,i,j))
 endfor
endfor

t1x=wx*dvxdx+wy*dvxdy+wz*dvxdz
t1y=wx*dvydx+wy*dvydy+wz*dvydz
t1z=wx*dvzdx+wy*dvzdy+wz*dvzdz

end


pro t2,pre,rho,t2x,t2y,t2z,xax,yax

nptx=n_elements(xax)
npty=n_elements(yax)

gxpre=dblarr(nptx,npty,nptx)
gypre=dblarr(nptx,npty,nptx)
gzpre=dblarr(nptx,npty,nptx)
gxrho=dblarr(nptx,npty,nptx)
gyrho=dblarr(nptx,npty,nptx)
gzrho=dblarr(nptx,npty,nptx)

grad,pre,gxpre,gypre,gzpre,xax,yax


grad,1.d0/rho,gxrho,gyrho,gzrho,xax,yax

cross,gxrho,gyrho,gzrho,gxpre,gypre,gzpre,t2x,t2y,t2z


end




pro t3,bx,by,bz,rho,t3x,t3y,t3z,xax,yax

nptx=n_elements(xax)
npty=n_elements(yax)

bnbx=dblarr(nptx,npty,nptx)
bnby=dblarr(nptx,npty,nptx)
bnbz=dblarr(nptx,npty,nptx)

gbpx=dblarr(nptx,npty,nptx)
gbpy=dblarr(nptx,npty,nptx)
gbpz=dblarr(nptx,npty,nptx)

ax=dblarr(nptx,npty,nptx)
ay=dblarr(nptx,npty,nptx)
az=dblarr(nptx,npty,nptx)

grad,1.d0/rho,ax,ay,az,xax,yax


t1,bx,by,bz,bx,by,bz,bnbx,bnby,bnbz,xax,yax

bpr=(bx^2.d0+by^2.d0+bz^2.d0)/2.d0

grad,bpr,gbpx,gbpy,gbpz,xax,yax

bnbx=gbpx-bnbx
bnby=gbpy-bnby
bnbz=gbpz-bnbz

cross,ax,ay,az,bnbx,bnby,bnbz,t3x,t3y,t3z

end




pro t4,bx,by,bz,rho,t4x,t4y,t4z,xax,yax

nptx=n_elements(xax)
npty=n_elements(yax)

bnbx=dblarr(nptx,npty,nptx)
bnby=dblarr(nptx,npty,nptx)
bnbz=dblarr(nptx,npty,nptx)

t1,bx,by,bz,bx,by,bz,bnbx,bnby,bnbz,xax,yax

curl,bnbx,bnby,bnbz,t4x,t4y,t4z,xax,yax

t4x=t4x/rho
t4y=t4y/rho
t4z=t4z/rho


end











;; put your numbers
;; npt=480 horizontal directions
;; nptz=100 vertical direction


;; axes
;; xax=dindgen(npt)/npt*12.0d8 horizontal axis
;; zax=dindgen(nptz)/nptz*1.4d8 vertical axis

;; coef=sqrt(4*3.14159265)



;; for i=0,nmod-1 do begin

;; here you get the model snapshot

;; vxx=vxx/rho
;; vyy=vyy/rho
;; vzz=vzz/rho

;; bxx=bxx*coef
;; byy=byy*coef
;; bzz=bzz*coef

;; wx=dblarr(npt,nptz,npt)
;; wy=dblarr(npt,nptz,npt)
;; wz=dblarr(npt,nptz,npt)

;; t1x=dblarr(npt,nptz,npt)
;; t1y=dblarr(npt,nptz,npt)
;; t1z=dblarr(npt,nptz,npt)

;; t2x=dblarr(npt,nptz,npt)
;; t2y=dblarr(npt,nptz,npt)
;; t2z=dblarr(npt,nptz,npt)

;; t3x=dblarr(npt,nptz,npt)
;; t3y=dblarr(npt,nptz,npt)
;; t3z=dblarr(npt,nptz,npt)

;; t4x=dblarr(npt,nptz,npt)
;; t4y=dblarr(npt,nptz,npt)
;; t4z=dblarr(npt,nptz,npt)



;; curl,vxx,vyy,vzz,wx,wy,wz,xax,zax

;; t1,wx,wy,wz,vxx,vyy,vzz,t1x,t1y,t1z,xax,zax

;; t2,pre,rho,t2x,t2y,t2z,xax,zax

;; t3,bxx,byy,bzz,rho,t3x,t3y,t3z,xax,zax

;; t4,bxx,byy,bzz,rho,t4x,t4y,t4z,xax,zax




;; endfor





;; end


