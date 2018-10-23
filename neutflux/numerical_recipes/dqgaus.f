      SUBROUTINE dqgaus(func,a,b,ss)
      REAL*8 a,b,ss,func
      EXTERNAL func
      INTEGER j
      REAL*8 dx,xm,xr,w(5),x(5)
      SAVE w,x
      DATA w/.2955242247D0,.2692667193D0,.2190863625D0,.1494513491D0,
     *.0666713443D0/
      DATA x/.1488743389D0,.4333953941D0,.6794095682D0,.8650633666D0,
     *.9739065285D0/
      xm=0.5D0*(b+a)
      xr=0.5D0*(b-a)
      ss=0.D0
      do 11 j=1,5
        dx=xr*x(j)
        ss=ss+w(j)*(func(xm+dx)+func(xm-dx))
11    continue
      ss=xr*ss
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Qi.
