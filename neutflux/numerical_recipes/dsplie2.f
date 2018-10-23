      SUBROUTINE dsplie2(x1a,x2a,ya,m,n,y2a)
      INTEGER m,n,NN
      REAL*8 x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=100)
CU    USES spline
      INTEGER j,k
      REAL*8 y2tmp(NN),ytmp(NN)
      do 13 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
11      continue
        call dspline(x2a,ytmp,n,1.e30,1.e30,y2tmp)
        do 12 k=1,n
          y2a(j,k)=y2tmp(k)
12      continue
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Qi.
