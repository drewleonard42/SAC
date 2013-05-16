function XYZ_coord, xyzmt,nnp,npp,i
  xyz=dblarr(1)
  for j=0, nnp-1 do begin       ; kolichestvo peremescheniy
     if j eq 0 then begin
        xyz=xyzmt(i)
     endif else begin
        xyz=[xyz,xyzmt(j*npp+i)]
     endelse  
  endfor

  RETURN,xyz
end  

function B_x,x_i,y_i,z_i
  RETURN,y_i
end  
function B_y,x_i,y_i,z_i
  RETURN,x_i
end  

function B_z,x_i,y_i,z_i
  RETURN,z_i
end  

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW

   ; Create some windows.
   Window, Title='IDL - Data Window', xsize=600,ysize=600, $
        /FREE, XPOS = 800, YPOS = 400
   dataWin = !D.Window


   Window, 1, Title='IDL - Trajectory', xsize=900,ysize=300, $
           XPOS = 800, YPOS = 400

   Window, 2, Title='IDL - Vz', xsize=350,ysize=330, $
           XPOS = 800, YPOS = 750


print, 'data or function ? (data "1", function "2")'
;read,aa

aa=1

;**********  begin data from file

ni=0
nj=0
mu=4.0*!PI/1.0e7
ijj=0.d0
ijjp=0.d0
ij=0


i_files=1
n_files=1



xyz_color_point=dblarr(1)

x0n=dblarr(1)
y0n=dblarr(1)
z0n=dblarr(1)

xmt=dblarr(1)
ymt=dblarr(1)
zmt=dblarr(1)


cxmi=0.7d0
cymi=0.7d0
czmi=0.7d0

cxma=1.3d0
cyma=1.3d0
czma=1.3d0

while i_files le n_files do begin

close,1
;openr,1,'/fastdata/smq11sjm/3D_data/3D_tube128_vert_p30_C_A90000_400.out',/f77_unf
;openr,1,'/fastdata/smq11sjm/3D_data/3D_tube128_vert_p30_C_A90000_400.out',/f77_unf
;openr,1,'/fastdata/smq11sjm/3D_data/3D_tube128_vert_p30_A90000_100S.out',/f77_unf
;openr,1,'/fastdata/smq11sjm/3D_data/3D_tube128_p30_A90000_horiz.out',/f77_unf


;***********************  files path *************************

dr='/home/smq11sjm/sac_3D/Visualisation/pngt'



;; traj=dr+'/ends_moves_three_full/traj_three/'
;; traj_proj=dr+'/ends_moves_three_full/traj_proj_three/'
;; vz_slice=dr+'/ends_moves_three_full/vz_slice/'

traj=dr+'/traj_three/'
traj_proj=dr+'/traj_proj_three/'
vz_slice=dr+'/vz_slice/'

tarr=dblarr(1)
headline='                                                                               '
it=long(1)
ndim=long(1)
neqpar=long(1)
nw=long(1)
varname='                                                                               '
time=double(1)
dum=long(1)
dumd=long(1)
nn=0

while not(eof(1)) do begin
readu,1,headline
print,'head = ', headline
readu,1,it,time,ndim,neqpar,nw
print,'it, time', it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
tarr=[tarr,time]
ndim=abs(ndim)
print, '*',ndim 
nx=lonarr(ndim)
readu,1,nx
print, 'nx = ', nx
eqpar=dblarr(neqpar)
readu,1,eqpar
print,eqpar
readu,1,varname
print,varname

n1=nx(0) 
n2=nx(1)
n3=nx(2)

xx=dblarr(n1,n2,n3,ndim)

print, '****************'

readu,1,xx

wi=dblarr(n1,n2,n3)
w=dblarr(n1,n2,n3,nw)

for iw=0,nw-1 do begin
 readu,1,wi
  w(*,*,*,iw)=wi
  print, iw
endfor

if nj eq 0 then dt=time
nj=1

;close,1

;***********************************************
nzs=0
;nz=190

nxs=0
nxe=127

nys=0
nye=127

;nzs=0
nz=120

;nxs=0
;nxe=99

;nys=0
;nye=99
;stop
;***********************************************

bx=dblarr(nz-nzs+1,nxe-nxs+1,nye-nys+1)
by=dblarr(nz-nzs+1,nxe-nxs+1,nye-nys+1)
bz=dblarr(nz-nzs+1,nxe-nxs+1,nye-nys+1)

vx=dblarr(nz-nzs+1,nxe-nxs+1,nye-nys+1)
vy=dblarr(nz-nzs+1,nxe-nxs+1,nye-nys+1)
vz=dblarr(nz-nzs+1,nxe-nxs+1,nye-nys+1)

rho=dblarr(nz-nzs+1,nxe-nxs+1,nye-nys+1)

if aa eq 1 then begin

;************************** B E G I N   F I L E *******************************

bx(0:nz-nzs,0:nxe-nxs,0:nye-nys)=w(nzs:nz,nxs:nxe,nys:nye,6)+w(nzs:nz,nxs:nxe,nys:nye,11)
by(0:nz-nzs,0:nxe-nxs,0:nye-nys)=w(nzs:nz,nxs:nxe,nys:nye,7)+w(nzs:nz,nxs:nxe,nys:nye,12)
bz(0:nz-nzs,0:nxe-nxs,0:nye-nys)=w(nzs:nz,nxs:nxe,nys:nye,5)+w(nzs:nz,nxs:nxe,nys:nye,10)

bx=bx*sqrt(mu)*1.0e4
by=by*sqrt(mu)*1.0e4
bz=bz*sqrt(mu)*1.0e4

vx(0:nz-nzs,0:nxe-nxs,0:nye-nys)=w(nzs:nz,nxs:nxe,nys:nye,2)
vy(0:nz-nzs,0:nxe-nxs,0:nye-nys)=w(nzs:nz,nxs:nxe,nys:nye,3)
vz(0:nz-nzs,0:nxe-nxs,0:nye-nys)=w(nzs:nz,nxs:nxe,nys:nye,1)

rho(0:nz-nzs,0:nxe-nxs,0:nye-nys)=w(nzs:nz,nxs:nxe,nys:nye,0)+w(nzs:nz,nxs:nxe,nys:nye,9)
;bx(nzs:nz,nxs:nxe,nys:nye)=w(nzs:nz,nxs:nxe,nys:nye,11)
;by(nzs:nz,nxs:nxe,nys:nye)=w(nzs:nz,nxs:nxe,nys:nye,12)
;bz(nzs:nz,nxs:nxe,nys:nye)=w(nzs:nz,nxs:nxe,nys:nye,10)

;************************** E N D   F I L E *******************************

endif


   
   WSet, dataWin


ll=5

loadct, ll


!p.multi=[0,1,0,0,1]

scale=1.0d6

z=dblarr(nz-nzs+1)
x=dblarr(nxe-nxs+1)
y=dblarr(nye-nys+1)


;****************end data from file 

;--------------------
; Make a vector of 16 points, A[i] = 2pi/16: 
A = FINDGEN(17) *  (!PI*2/16.) 
; Define the symbol to be a unit circle with 16 points,  
; and set the filled flag: 
USERSYM, 0.5*COS(A), 0.5*SIN(A), /fill
;--------------------



for i=0,nz-nzs do begin
   if aa eq 1 then begin
      z(i)=xx(nzs+i,0,0,0)/scale
   endif else begin
      z(i)=(nzs+i-nz/2.0)
   endelse 
endfor

for j=0,nxe-nxs do begin
   if aa eq 1 then begin
      x(j)=xx(0,nxs+j,0,1)/scale
   endif else begin
      x(j)=((nxs+j-nxe)/2.0)
   endelse 
endfor

for k=0,nye-nys do begin
   if aa eq 1 then begin
      y(k)=xx(0,0,nys+k,2)/scale
   endif else begin
      y(k)=((nys+k-nye)/2.0)
   endelse  
endfor




sc=0.9d0

mi_x=min(x) ;*sc
ma_x=max(x) ; *sc

mi_y=min(y) ;*sc
ma_y=max(y) ;*sc

mi_z=min(z)
ma_z=max(z)

b_total=sqrt(bx(0,*,*)^2.0+by(0,*,*)^2.0+bz(0,*,*)^2.0)

v_total=sqrt(vx(*,*,*)^2.0+vy(*,*,*)^2.0+vz(*,*,*)^2.0)/rho(*,*,*)

tvframe,b_total, /bar


max_b=max(b_total)
min_b=min(b_total)


;******************************************************

nn=4.0

;******************************************************

xyz_color=dblarr(1)
xyz_color(0)=1.0
color_min=40
color_max=256.0


;******************************************************************
dL=2000.d0/scale
;******************************************************************

delta_x=(ma_x-mi_x)/nn
delta_y=(ma_y-mi_y)/nn

;; xx0=0.1
;; yy0=0.15
;; xx1=0.85
;; yy1=0.90
xx0=0.0
yy0=0.0
xx1=1.99
yy1=1.99
SURFACE, DIST(5), /NODATA, /SAVE, XRANGE=[mi_x, ma_x],  $
   YRANGE=[mi_y, ma_y], ZRANGE=[mi_z, ma_z], XSTYLE=1, $
   YSTYLE=1, ZSTYLE=1, CHARSIZE=2.5, xtitle='x [Mm]',ytitle='y [Mm]', $
   ztitle='z [Mm]', POSITION=[xx0, yy0, xx1, yy1, 0.1, 0.95], az=30.0, ax=30.0
   ;BACKGROUND=FSC_Color('ivory')  

   
 levels = 80
loadct,3, NColors=40, Bottom=3

step = (Max(vz(12,*,*)/rho(12,*,*)) - Min(vz(12,*,*)/rho(12,*,*))) / levels
userLevels = IndGen(levels) * step + Min(vz(12,*,*)/rho(12,*,*))

;step = (Max(v_total(12,*,*)) - Min(v_total(12,*,*))) / levels
;userLevels = IndGen(levels) * step + Min(v_total(12,*,*))

;vz(12,*,*)

Contour, vz(12,*,*)/rho(12,*,*), x, y, /Fill, C_Colors=Indgen(levels)+3, Background=1, /T3D,$
   Levels=userLevels, POSITION=[xx0, yy0, xx1, yy1, 0.1, 0.95], Color=black, ZVALUE = 0.1, $
   XSTYLE=1, YSTYLE=1, ZSTYLE=1, /NOCLIP,CHARSIZE=2.5, XRANGE=[mi_x, ma_x],  $
   YRANGE=[mi_y, ma_y], ZRANGE=[mi_z, ma_z]
  
   
loadct,ll

tvlct, r,g,b, /get

SURFACE, DIST(5), /NODATA, /SAVE, XRANGE=[mi_x, ma_x], /noerase, $
   YRANGE=[mi_y, ma_y], ZRANGE=[mi_z, ma_z], XSTYLE=1, $
   YSTYLE=1, ZSTYLE=1, CHARSIZE=2.5, xtitle='x [Mm]',ytitle='y [Mm]', $
   ztitle='z [Mm]', POSITION=[xx0, yy0, xx1, yy1, 0.1, 0.95], az=30.0, ax=30.0
   ;BACKGROUND=FSC_Color('ivory')  



;AXIS, XAXIS=1, XRANGE=[mi_x, ma_x], /T3D, charsize=2.0, /NODATA
;AXIS, YAXIS=1, XRANGE=[mi_y, ma_y], /T3D, charsize=2.0, /NODATA
ii=0


Rr=(ma_x-mi_x)/2.d0

aas=-0.02d0 ;for zoom
;aas=-0.5d0 ;for zoom





;***************** begin start positions ***************************
if ij eq 0 then begin

x0=ma_x
y0=ma_y
z0=1.0d0*ma_z

x0_s=x0
y0_s=y0
z0_s=z0

y0_pred=y0_s

while y0_s ge mi_y do begin

while (x0_s ge mi_x) do begin

if Rr^2.d0-((y0_s-mi_y-(ma_y-mi_y)/2.d0)^2.d0+(x0_s-mi_x-(ma_x-mi_x)/2.d0)^2.d0) ge aas then begin
  if ij eq 0 then begin
   x0n[ij]=x0_s
   y0n[ij]=y0_s
   z0n[ij]=z0
  endif else begin
   x0n=[x0n,x0_s]
   y0n=[y0n,y0_s]
   z0n=[z0n,z0]
  endelse
  ij=ij+1  
endif

  x0_s=x0_s-delta_x
  x0=x0_s
  y0=y0_pred
  z0=1.d0*ma_z    
endwhile
  y0_s=y0_s-delta_y
  y0=y0_s
  y0_pred=y0
  
  x0=ma_x
  x0_s=ma_x
  z0=1.d0*ma_z      
endwhile

endif

;***************** end start positions ***************************


;***************** begin start NEW positions ***************************
 for ii=0, n_elements(x0n)-1 do begin 
     xc=interpol(dindgen(nxe-nxs+1),x,x0n[ii])
     yc=interpol(dindgen(nye-nys+1),y,y0n[ii])
     zc=interpol(dindgen(nz-nzs+1),z,z0n[ii])

     Vvx=interpolate(vx/rho, zc,xc,yc)
     Vvy=interpolate(vy/rho, zc,xc,yc)
     Vvz=interpolate(vz/rho, zc,xc,yc)
     
     x0n[ii]=x0n[ii]+Vvx*dt/scale
     y0n[ii]=y0n[ii]+Vvy*dt/scale
     z0n[ii]=z0n[ii]+Vvz*dt/scale   
     print, '*', x0n[ii], y0n[ii], z0n[ii], Vvx, VVy, Vvz, ii       

endfor
;stop
;***************** END start NEW positions ***************************

xm=dblarr(1)
ym=dblarr(1)
zm=dblarr(1)


xm(0)=x0n(0)
ym(0)=y0n(0)
zm(0)=z0n(0)

for ii=0, n_elements(x0n)-1 do begin

x0=x0n[ii]
y0=y0n[ii]
z0=z0n[ii]

while ((x0 ge cxmi*mi_x) and (x0 le cxma*ma_x) and $
       (y0 ge cymi*mi_y) and (y0 le cyma*ma_y) and $
       (z0 ge czmi*mi_z) and (z0 le czma*ma_z)) do begin 


     xc=interpol(dindgen(nxe-nxs+1),x,x0)
     yc=interpol(dindgen(nye-nys+1),y,y0)
     zc=interpol(dindgen(nz-nzs+1),z,z0)

;stop
  Bbx=interpolate(bx, zc, xc, yc)
  Bby=interpolate(by, zc, xc, yc)
  Bbz=interpolate(bz, zc, xc, yc)  

   Bb=sqrt(Bbx^2.0+Bby^2.0+Bbz^2.0)

   cl=color_max-(max_b-Bb)/(max_b-min_b)*(color_max-color_min)
   
     
   dx=-dL*Bbx/Bb
   dy=-dL*Bby/Bb
   dz=-dL*Bbz/Bb    


   x0=x0+dx
   y0=y0+dy
   z0=z0+dz



   xm=[xm,x0]
   ym=[ym,y0]
   zm=[zm,z0]
   
   xyz_color=[xyz_color,cl]
   
;    PLOTS, xm,ym,zm, psym=3, /T3D 
endwhile


endfor

;cl=color_min+ijj*1.d0/2.0
cl=color_min+ijj*1.d0/6.0

   
for i=0,n_elements(x0n)-1 do begin

if ijjp eq 0.d0 then begin
   xmt(0)=x0n[i]
   ymt(0)=y0n[i]
   zmt(0)=z0n[i]
endif else begin
   xmt=[xmt,x0n[i]]
   ymt=[ymt,y0n[i]]
   zmt=[zmt,z0n[i]]
endelse

if ijjp eq 0.d0 then begin
  if i eq 0 then xyz_color_point(0)=cl
endif else begin
  if i eq 0 then xyz_color_point=[xyz_color_point,cl]
endelse

ijjp=ijjp+1

endfor


ijj=ijj+1.d0


PLOTS, [ma_x,ma_x], [ma_y,ma_y], [mi_z, ma_z], /T3D
PLOTS, xm,ym,zm, color=xyz_color, psym=8, /T3D 

PLOTS, [mi_x,mi_x], [mi_y,mi_y], [mi_z, ma_z], /T3D
PLOTS, [ma_x,ma_x], [mi_y,mi_y], [mi_z, ma_z], /T3D

surf=dblarr(nxe-nxs+1,nye-nys+1)

nzz=nz
surf(*,*)=sqrt(bx(nzz-1,*,*)^2.d0+by(nzz-1,*,*)^2.d0+bz(nzz-1,*,*)^2.d0)



CONTOUR, surf, x, y, /T3D, /noerase, zvalue=0.95,/follow, $
         XRANGE=[mi_x, ma_x], YRANGE=[mi_y, ma_y], $
	 CHARSIZE=2.0,  $        ;levels=[15.1,15.4,15.8,16.2,16.6,17.0,17.5], $ ; levels for 0.95
	 POSITION=[xx0, yy0, xx1, yy1, 0.95, 0.95], $
	 C_charsize=2.d0, xstyle=1, ystyle=1

;--------------------
; Make a vector of 16 points, A[i] = 2pi/16: 
A = FINDGEN(17) *  (!PI*2/16.) 
; Define the symbol to be a unit circle with 16 points,  
; and set the filled flag: 
USERSYM, 1.5*COS(A), 1.5*SIN(A), /fill
;--------------------


;****************** begin recalculating and plots top points ******************

np=10 ; kolichesto tochek for trajectory

;np=21 ; kolichesto tochek for trajectory

p=intarr(np)

; tochki for trajectory

;; p[0]= 7
;; p[1]=10
;; p[2]=12

for i=0, np-1 do p[i]=i

nnp=n_elements(xmt)/n_elements(x0n)
npp=n_elements(x0n)
print, '#############', nnp

for  i=0, np-1 do begin ; kolichestvo tochek

xxc=XYZ_coord(xmt,nnp,npp,p[i])
yyc=XYZ_coord(ymt,nnp,npp,p[i])
zzc=XYZ_coord(zmt,nnp,npp,p[i])

PLOTS, xxc,yyc,zzc, color=xyz_color_point, psym=8, /T3D

endfor

;****************** end recalculating and plots top points ******************




     
ss='time ='+strTrim(string(FORMAT='(6F10.2)', time),2)+' it ='+strTrim(string(it),1)
 xyouts,50,20, ss, /device, color=200	


indexs=strtrim(ni,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   


;bar 1 = 'source'  2 = 'field'
for bar=1,2 do begin

case bar of
   1: begin   ;*******'source' 
	      tek_color
	      xyouts, 525, 160, 'Vz [m/s]', orientation=90, /device, charsize=1.5

   
              loadct,3 

              px = [465, 525]  
              py = [40, 280] 
              mina = min(vz(12,*,*)/rho(12,*,*))
              maxa = max(vz(12,*,*)/rho(12,*,*))
	     end 
 
   2: begin  ;*******'field' 
	      tek_color
	      xyouts, 525, 370, 'B [Gauss]', orientation=90, /device, charsize=1.5

   
             tvlct, r,g,b
           
             px = [465, 525]     ;Position of frame in device units
             py = [290,530]
             mina=min_b   
             maxa=max_b
            end 
   
endcase   
;******************* begin bar source *******************

sx = px[1]-px[0]                ;Size of frame in device units
sy = py[1]-py[0]


   bx    = fix(px[0]+sx*1.2)
   by    = fix(py[0])
   bsx   = fix(sx*0.2)
   bsy   = fix(sy)

   barpos= [bx,by,bx+bsx,by+bsy]



     mm=max([abs(mina),abs(maxa)])     
        barim=findgen(1,bsy)/(bsy-1)*(maxa-mina)+mina
        barim=rebin(barim,bsx,bsy,/sample)
        if keyword_set(noscale) then begin
           tv, 0>barim<mcol ,bx,by,/device
        endif else begin

	     bb=barim
       if bar eq 1 then begin 
        for i=0,bsx-1 do begin
	 for j=0, bsy-1 do begin
	  if barim[i,j] lt 0.d0 then begin
	    bb[i,j]=127.d0*(barim[i,j]+abs(mina))/abs(mina)
	  endif else begin
   	    bb[i,j]=128.d0+127.d0*barim[i,j]/abs(maxa)
	  endelse
         endfor
	 endfor
	endif else begin 

        for i=0,bsx-1 do begin
	 for j=0, bsy-1 do begin

  	    bb[i,j]=color_min+(color_max-color_min)*barim[i,j]/abs(maxa)
	  
	 endfor
	endfor  
        endelse

        tv, bb,bx,by, /device
        endelse
	
;loadct,3

BRANGE=float([mina,maxa])


   plot,[0,1],BRANGE,/nodata,/noerase,pos=barpos,/device,xsty=5,ysty=5
   plots,[bx+bsx,bx,bx,bx+bsx],[by,by,by+bsy,by+bsy],/device
   axis,yaxis=1,bx+bsx,/device,yrange=BRANGE,ystyle=1   $
     ,  ytitle=BTITLE $
     ,  TICKLEN=-.15,charsize=1.0,YTICKS=BTICKS,YMINOR=BMINOR ;, YTICKFORMAT=YT

endfor ; bar
;******************* end bar source *******************

image_p = TVRD_24()

write_png,traj+indexss+'.png',image_p, red,green, blue

;********************** begin side trajectory *****************************
   WSet, 1
   !p.multi = [0,3,1,0,1]


;--------------------
; Make a vector of 16 points, A[i] = 2pi/16: 
A = FINDGEN(17) *  (!PI*2/16.) 
; Define the symbol to be a unit circle with 16 points,  
; and set the filled flag: 
USERSYM, 0.6*COS(A), 0.6*SIN(A), /fill
;--------------------


   css=2.0
   plot, xxc[0:n_elements(xxc)-1],yyc[0:n_elements(yyc)-1], title='x-y', /nodata, $
         charsize=css, XRANGE=[mi_x, ma_x], YRANGE=[mi_y, ma_y], $
	 xtitle='x', ytitle='y'

for i=0, np-1 do begin ; kolichestvo tochek

   xxc=XYZ_coord(xmt,nnp,npp,p[i])
   yyc=XYZ_coord(ymt,nnp,npp,p[i])
   plots, xxc,yyc, color=xyz_color_point, psym=8
   
endfor

   plot, xxc[0:n_elements(xxc)-1],zzc[0:n_elements(zzc)-1], title='x-z', /nodata, $
         charsize=css, XRANGE=[mi_x, ma_x], YRANGE=[0.8, 1.0], $
	 xtitle='x', ytitle='z'	 
   
for i=0, np-1 do begin; kolichestvo tochek

   xxc=XYZ_coord(xmt,nnp,npp,p[i])
   zzc=XYZ_coord(zmt,nnp,npp,p[i])	 
   plots, xxc,zzc, color=xyz_color_point, psym=8
   
endfor

   plot, yyc[0:n_elements(yyc)-1],zzc[0:n_elements(zzc)-1], title='y-z', /nodata, $
         charsize=css, XRANGE=[mi_y, ma_y], YRANGE=[0.8, 1.0], $
	 xtitle='y', ytitle='z'	 

for i=0, np-1 do begin ; kolichestvo tochek

   yyc=XYZ_coord(ymt,nnp,npp,p[i])
   zzc=XYZ_coord(zmt,nnp,npp,p[i])	 
   plots, yyc,zzc, color=xyz_color_point, psym=8
   
endfor   


ss='time ='+strTrim(string(FORMAT='(6F10.2)', time),2)+' it ='+strTrim(string(it),1)
 xyouts,20,10, ss, /device, color=200	

;********************** end side trajectory *****************************	 



image_p = TVRD_24()

write_png,traj_proj+indexss+'.png',image_p, red,green, blue
   
   Wset,2
;   Window, 2, Title='Vz', xsize=320,ysize=300, $
;           XPOS = 850, YPOS = 750

!P.multi=[0,1,1,0,1]

hh=111

hight='h='+strTrim(string(FORMAT='(6F10.2)',xx(hh,1,1,0)/1.0d6),1)+' [Mm]'



tvframe,w(hh,*,*,1)/(w(hh,*,*,0)+w(hh,*,*,9)),/sample, /bar, title='Vz '+hight,$
        xtitle='x', ytitle='y',charsize=1.0, CT='dPdT'

ss='time ='+strTrim(string(FORMAT='(6F10.2)', time),2)+' it ='+strTrim(string(it),1)
 xyouts,10,6, ss, /device, color=200	


image_p = TVRD_24()
write_png,vz_slice+indexss+'.png',image_p, red,green, blue



ni=ni+1


endwhile

i_files=i_files+1

endwhile 
end
