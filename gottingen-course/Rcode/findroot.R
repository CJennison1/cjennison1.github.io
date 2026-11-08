findroot=function(fofx,yt,xstart,delx,eps,par)
{
 # Finds the value of x for which fofx(x)=y where fofx is an
 # increasing function of x
 #
 # IN
 #
 #   fofx     the function whose root is to be found
 #   yt       the target value
 #   xstart   starting value for search over x
 #   delx     initial increment in seach over x
 #   eps      tolerance: root will have fofx within eps of yt
 #   par      additional parameters to be passed to fofx
 #
 # OUT
 #
 #   xroot   solution of fofx(x,par)=yt

 x0=xstart
 y0=fofx(x0,par)

 if(y0<yt)
 {
  xlo=x0
  ylo=y0
  xhi=x0+delx
  yhi=fofx(xhi,par)
  #  print(c("A1: xlo,xhi,ylo,yhi=",xlo,xhi,ylo,yhi))
  while(yhi<yt)
  {
   xlo=xhi
   ylo=yhi
   xhi=xlo+delx
   yhi=fofx(xhi,par)
  #  print(c("A2: xlo,xhi,ylo,yhi=",xlo,xhi,ylo,yhi))
  }
 } else
 {
  xhi=x0
  yhi=y0
  xlo=x0-delx
  ylo=fofx(xlo,par)
  while(ylo>yt)
  {
   xhi=xlo
   yhi=ylo
   xlo=xhi-delx
   ylo=fofx(xlo,par)
   # print(c("A3: xlo,xhi,ylo,yhi=",xlo,xhi,ylo,yhi))
  }
 }
  # print(c("A4: xlo,xhi,ylo,yhi=",xlo,xhi,ylo,yhi))

 while((yhi-ylo)>eps)
 {
  xmid=(xlo+xhi)/2
  ymid=fofx(xmid,par)
  if(ymid<yt)
  {
   xlo=xmid
   ylo=ymid
  } else
  {
   xhi=xmid
   yhi=ymid
  }
  # print(c("B: xlo,xhi,ylo,yhi=",xlo,xhi,ylo,yhi))
 }
 xroot=xlo+(xhi-xlo)*(yt-ylo)/(yhi-ylo)
 froot=fofx(xroot,par)
 # print(c("C: xroot,froot=",xroot,froot))

 return(xroot)
}