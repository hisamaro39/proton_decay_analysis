      SUBROUTINE dsplin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      INTEGER m,n,NN
      REAL*8 x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=100)
CU    USES dspline,dsplint
      INTEGER j,k
      REAL*8 y2tmp(NN),ytmp(NN),yytmp(NN)
      do 12 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
11      continue
        call dsplint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
12    continue
      call dspline(x1a,yytmp,m,1.e30,1.e30,y2tmp)
      call dsplint(x1a,yytmp,y2tmp,m,x1,y)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Qi.
