; INPUTS:
;        data = Two dimensional array of numerical type
; KEYWORD PARAMETERS:
;     XRANGE   : Array with at least 2 elements giving the range for
;                the annotation of the x-axis.
;                Only min(XRANGE) and max(XRANGE) are used.
;     YRANGE   : Same as XRANGE but for y-axis.
;     POSITION : (output) position parameter for the axis in device units,
;                may be used in subsequent calls to PLOT or CONTOUR.
;     /sample  : If set, nearest neighbourhood method is used in
;                resizing the data array.
;                Default : averaging or linear interpolation.
;     /center  : If set, tickmarks are placed at the center of each pixel.
;     /aspect  : If set, the aspect ratio of the data array is preserved
;     /noscale : If set, data is not scaled
;     /bar     : If set, an intensity bar is displayed.
;     BRANGE   : Range for annotation of the intensity bar.
;     BTITLE   : Title of the intensity bar.
;     BTICKS,BMINOR : Control the number of tickmarks for the intensity bar.
; Standard plotting keywords :
;      TITLE,SUBTITLE,XTITLE,YTITLE,TICKLEN,CHARSIZE,XSTYLE,YSTYLE,
;      XTICKS,YTICKS,XMINOR,YMINOR
;                Normally they are just passed through,

pro tvframe , a    $
    , sample=sample , center=center , aspect=aspect , noscale=noscale    $
    , POSITION=POSITION    $
    , XRANGE=XRANGI , YRANGE=YRANGI   $
    , TITLE=TITLE , XTITLE=XTITLE , YTITLE=YTITLE , SUBTITLE=SUBTITLE   $
    , TICKLEN=TICKLEN , CHARSIZE=CHARSIZE   $
    , XTICKS=XTICKS , YTICKS=YTICKS , XMINOR=XMINOR , YMINOR=YMINOR $
    , XSTYLE=XSTYLI , YSTYLE=YSTYLI   $
    , bar=bar , BTITLE=BTITLE , BRANGE=BRANGE  $
    , BTICKS=BTICKS , BMINOR=BMINOR,YTICKFORMAT=YTICKFORMAT,CT=CT

on_error,2
sa=size(a) 
if sa(0) ne 2 then begin
print,'* Non-2D data, reformatting...'
a=reform(a)
endif

sa=size(a)
if sa(0) ne 2 then goto, errout



mina=min(a,max=maxa)
;
; set keyword parameters to default values if not present
if n_elements(XRANGI) eq 0 then XRANGI=[0,sa(1)-1]
if n_elements(YRANGI) eq 0 then YRANGI=[0,sa(2)-1]
if n_elements(   TITLE) eq 0 then    TITLE=''
if n_elements(SUBTITLE) eq 0 then SUBTITLE=''
if n_elements(  XTITLE) eq 0 then   XTITLE=''
if n_elements(  YTITLE) eq 0 then   YTITLE=''
if n_elements(TICKLEN) eq 0 then TICKLEN=-.01
if n_elements(XTICKS) eq 0 then XTICKS=0
if n_elements(XMINOR) eq 0 then XMINOR=0
if n_elements(YTICKS) eq 0 then YTICKS=0
if n_elements(YMINOR) eq 0 then YMINOR=0
if n_elements(CHARSIZE) eq 0 then CHARSIZE=1.0
;XSTYLE=1  &  YSTYLE=1
if n_elements(XSTYLI) eq 1 then XSTYLE=( 1 or XSTYLI ) and 29
if n_elements(YSTYLI) eq 1 then YSTYLE=( 1 or YSTYLI ) and 29
if n_elements(BTITLE) eq 0 then BTITLE=''
if n_elements(BRANGE) eq 0 then BRANGE=float([mina,maxa])
if n_elements(BTICKS) eq 0 then BTICKS=0
if n_elements(BMINOR) eq 0 then BMINOR=0
if n_elements(YTICKFORMAT) eq 0 then YTICKFORMAT=''
if n_elements(CT) eq 0 then CT='VpVd'
;

XRANGE=XRANGI;float(minmax(XRANGI))
YRANGE=YRANGI;reverse(float(minmax(YRANGI)))

;
if keyword_set(center) then begin
    xunit=0.5*(XRANGE(1)-XRANGE(0))/float(sa(1)-1)
    yunit=0.5*(YRANGE(1)-YRANGE(0))/float(sa(2)-1)
    XRANGE(0)=XRANGE(0)-xunit  &  XRANGE(1)=XRANGE(1)+xunit
    YRANGE(0)=YRANGE(0)-yunit  &  YRANGE(1)=YRANGE(1)+yunit
endif else begin
    xunit=(XRANGE(1)-XRANGE(0))/float(sa(1)-1)
    yunit=(YRANGE(1)-YRANGE(0))/float(sa(2)-1)
    XRANGE(1)=XRANGE(1)+xunit
    YRANGE(1)=YRANGE(1)+yunit
endelse
;
XRANGE=XRANGI;float(minmax(XRANGI))
YRANGE=YRANGI;reverse(float(minmax(YRANGI)))
;loadct,4

plot,[0],[0],xrange=XRANGE,yrange=YRANGE,/nodata $
     ,  TITLE=' ',XTITLE=XTITLE,YTITLE=YTITLE    $
     ,  SUBTITLE=SUBTITLE      $
     ,  TICKLEN=TICKLEN,CHARSIZE=CHARSIZE   $
     ,  color=!p.background, xsty=4,ysty=5

case CT of                                                           
 'VpVd'    : mixct1
 'dPdT'    : mixct7
 'mixct3'  : mixct3 
 'mixct4'  : mixct4
 'mixct6'  : mixct6
 'mixct8'  : mixct8  
 'mixct10'  : mixct10   
 'mixct17' : mixct17   
 'mixct16' : mixct16    
 'rho0'    : loadct,3
 'P0'      : loadct,3
 'loadct1' : loadct,1
 'loadct2' : loadct,2 
 'loadct3' : loadct,3
 'loadct4' : loadct,4 
 'loadct5' : loadct,5
 'loadct6' : loadct,6 
 'loadct7' : loadct,7
 'loadct8' : loadct,8 
 'loadct9' : loadct,9
 'loadct10': loadct,10 
 'loadct11': loadct,11
 'loadct12': loadct,12 
 'loadct13': loadct,13
 'loadct14': loadct,14 
 'loadct15': loadct,15
 'loadct16': loadct,16 
 'loadct17': loadct,17
 'loadct18': loadct,18 
 'loadct19': loadct,19
 'loadct20': loadct,20 
 'loadct21': loadct,21
 'loadct22': loadct,22 
 'loadct23': loadct,23
 'loadct24': loadct,24 
endcase   


px = !x.window * !d.x_vsize     ;Position of frame in device units
py = !y.window * !d.y_vsize
sx = px(1)-px(0)                ;Size of frame in device units
sy = py(1)-py(0)
if keyword_set(bar) then sx = sx/1.25
if keyword_set(aspect) then begin
     f = float(sa(1))/sa(2)*sy/sx
     if f ge 1. then sy=sy/f else sx=sx*f
     sx=fix(sx)
endif
POSITION = [px(0),py(0),px(0)+sx,py(0)+sy]
if keyword_set(bar) then begin
   bx    = fix(px(0)+sx*1.04)
   by    = fix(py(0))
   bsx   = fix(sx*0.08)
   bsy   = fix(sy)
   barpos= [bx,by,bx+bsx,by+bsy]
endif
;
mcol=( !D.N_COLORS - 1) > 0
;

;nnx=256
;nny=256

nnx=1500
nny=1500

	mm=max([abs(mina),abs(maxa)]) 
   ;  b=a
     b=congrid(a,nnx,nny,/interp)
 
     nx=n_elements(b(*,1))
     ny=n_elements(b(1,*))

        mm=max([abs(min(b)),abs(max(b))])

        bb=b
        for i=0,nx-1 do begin
	 for j=0, ny-1 do begin
	  if b[i,j] lt 0.d0 then begin
	    bb[i,j]=127.d0*(b[i,j]+abs(mina))/abs(mina)
	  endif else begin
   	    bb[i,j]=128.d0+127.d0*b[i,j]/abs(maxa)
	  endelse
	 endfor
	endfor 	
       
       case CT of
       'rho0'    : tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'P0'      : tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct1' : tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct2' : tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct3' : tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct4' : tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct5' : tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct6' : tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct7' : tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct8' : tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct9' : tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct10': tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct11': tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct12': tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct13': tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct14': tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct15': tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct16': tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct17': tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct18': tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct19': tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct20': tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct21': tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct22': tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct23': tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       'loadct24': tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
       else      : tv, bb ,px(0),py(0),xsize=sx,ysize=sy,/device
       endcase 


     mm=max([abs(mina),abs(maxa)])     

        barim=findgen(1,255)/(255-1)*(maxa-mina)+mina

        barim=rebin(barim,255,255)
	
	     bb=barim
     	     
        for i=0,255-1 do begin
	 for j=0, 255-1 do begin
	  if barim[i,j] lt 0.d0 then begin
	    bb[i,j]=127.d0*(barim[i,j]+abs(mina))/abs(mina)
	  endif else begin
   	    bb[i,j]=128.d0+127.d0*barim[i,j]/abs(maxa)
	  endelse
	  
	 endfor
	endfor 

       case CT of
       'rho0'      : tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device
       'P0'        : tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device
       'loadct1' : tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct2' : tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct3' : tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct4' : tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct5' : tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct6' : tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct7' : tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct8' : tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct9' : tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct10': tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct11': tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct12': tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct13': tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct14': tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct15': tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct16': tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device
       'loadct17': tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct18': tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct19': tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct20': tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct21': tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct22': tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct23': tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       'loadct24': tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device       
       
       else        : tv, bb ,bx,by,xsize=bsx,ysize=bsy,/device
       endcase 
	

	
 if (mm gt 10.0) then YT='(I6)'
 if (mm le 10.0) then YT='(F10.1)'
 if (mm le 1.0) then YT='(F10.2)'
 if (mm le 0.1) then YT='(F10.3)'
 if (mm le 0.01) then YT='(F10.4)' 
if (mm le 0.001) then YT='' 

 
;
loadct,3
if keyword_set(bar) then begin
   plot,[0,1],BRANGE,/nodata,/noerase,pos=barpos,/device,xsty=5,ysty=5
   plots,[bx+bsx,bx,bx,bx+bsx],[by,by,by+bsy,by+bsy],/device
   axis,yaxis=1,bx+bsx,/device,yrange=BRANGE,ystyle=2   $
     ,  ytitle=BTITLE $
     ,  TICKLEN=-.15,CHARSIZE=CHARSIZE,YTICKS=BTICKS,YMINOR=BMINOR ;, YTICKFORMAT=YT
endif

plot,[0],[0],xrange=XRANGE,yrange=YRANGE,/nodata,/noerase,xstyle=1,ystyle=1 $
     ,  POSITION=POSITION  ,/device    $
     ,  TITLE=TITLE,XTITLE=XTITLE,YTITLE=YTITLE    $
     ,  SUBTITLE=SUBTITLE      $
     ,  TICKLEN=TICKLEN,CHARSIZE=CHARSIZE, XCharSize=1.2, YCharSize=1.2   $
     ,  XTICKS=XTICKS,XMINOR=XMINOR , YTICKS=YTICKS,YMINOR=YMINOR

return
;
errout: print,' TVFRAME: unable to reformat the data. The data must be 3-dimensional!'
return
end


