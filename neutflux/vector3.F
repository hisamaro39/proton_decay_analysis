c $Id: vector3.F,v 1.1 2007-11-10 09:35:29 skrep Exp $
c In this file, there are fandamental routines 
c for calculation of 3-dimension vectors
c
c    S.Hatakeyama, Tohoku Univ.
c
c The routines are follows
c/////////////////////////////////////////////////////
c subroutine v3set(v,x,y,z)   : set elements of vector
c subroutine v3copy(v1,v2)    : copy v1 to v2
c subroutine v3unit(v1,v2)    : v2 is unit vector of v1
c subroutine v3add(v1,v2,v3)  : v3 = v1 + v2
c subroutine v3sub(v1,v2,v3)  : v3 = v1 - v2
c subroutine v3mean(v1,v2,v3) : v3 = (v1 + v2)/2.
c subroutine v3mult(c,v1,v2)  : v2 = c * v1 
c subroutine v3rot(v1,th,phi,v2): v2 is th,phi rotaion of v1
c real function v3times(v1,v2): inner product v1 * v2
c real function v3abs(v)      : abs v
c real function v3cos(v1,v2)  : cosine of v1 and v2
c real function v3len(v1,v2)  : length of v1 and v2
c/////////////////////////////////////////////////////
c vector type is REAL for Fortran
c but routine begin with 'd' is double precision(i.e. real*8)
c---------------------------------------------------
      subroutine v3set(v,x,y,z)
c
c input:  v,x,y,z
c output: v
c---------------------------------------------------
      real v(3),x,y,z
c
      v(1) = x
      v(2) = y
      v(3) = z
c
      end
c---------------------------------------------------
      subroutine v3copy(v1,v2)
c
c input:  v1
c output: v2
c---------------------------------------------------
      real v1(3),v2(3)
c
      v2(1) = v1(1)
      v2(2) = v1(2)
      v2(3) = v1(3)
c
      end
c---------------------------------------------------
      subroutine v3add(v1,v2,v3)
c
c input:  v1, v2
c output: v3
c---------------------------------------------------
      real v1(3),v2(3),v3(3)
c
      v3(1) = v1(1) + v2(1)
      v3(2) = v1(2) + v2(2)
      v3(3) = v1(3) + v2(3)
c
      end
c---------------------------------------------------
      subroutine v3sub(v1,v2,v3)
c
c input:  v1, v2
c output: v3
c---------------------------------------------------
      real v1(3),v2(3),v3(3)
c
      v3(1) = v1(1) - v2(1)
      v3(2) = v1(2) - v2(2)
      v3(3) = v1(3) - v2(3)
c
      end
c---------------------------------------------------
      subroutine v3mean(v1,v2,v3)
c
c input:  v1, v2
c output: v3
c---------------------------------------------------
      real v1(3),v2(3),v3(3)
c
      v3(1) = ( v1(1) + v2(1) )/2.
      v3(2) = ( v1(2) + v2(2) )/2.
      v3(3) = ( v1(3) + v2(3) )/2.
c
      end
c---------------------------------------------------
      subroutine v3mult(c,v1,v2)
c
c input:  c, v1
c output: v2
c---------------------------------------------------
      real c,v1(3),v2(3)
c
      v2(1) = c * v1(1)
      v2(2) = c * v1(2)
      v2(3) = c * v1(3)
c
      end
c---------------------------------------------------
      real function v3times(v1,v2)
c
c input:  v1, v2
c output: v3times
c---------------------------------------------------
      real v1(3),v2(3)
c
      v3times = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
c
      return
      end
c---------------------------------------------------
      real function v3abs(v)
c
c input:  v
c output: v3abs
c---------------------------------------------------
      real v(3)
c
      if(.not.( v(1) .eq. 0. .and. v(2) .eq. 0. .and.
     &      v(3) .eq. 0. )) goto 23000
         v3abs = 0.
         return
23000 continue
      v3abs = sqrt( v(1)*v(1) + v(2)*v(2) + v(3)*v(3) )
c
      return
      end
c---------------------------------------------------
      subroutine v3unit(v1,v2)
c
c input:  v1, v2
c output: v3unit
c---------------------------------------------------
      real v1(3),v2(3)
      real abs_v1
      external v3abs
c
      abs_v1 = v3abs(v1)
c
      if(.not.( abs_v1 .eq. 0. ))goto 23002
         v2(1) = 0.
         v2(2) = 0.
         v2(3) = 0.
         return
23002 continue
      v2(1) = v1(1)/abs_v1
      v2(2) = v1(2)/abs_v1
      v2(3) = v1(3)/abs_v1
c
      end
c---------------------------------------------------
      real function v3cos(v1,v2)
c
c input:  v1, v2
c output: v3cos
c---------------------------------------------------
      real v1(3),v2(3)
      external v3times,v3abs
      real abs_v1,abs_v2
c
      abs_v1 = v3abs(v1)
      abs_v2 = v3abs(v2)
      if(.not.( abs_v1 * abs_v2 .eq. 0. ))goto 23004
         v3cos = 0.
         return
23004 continue
      v3cos = v3times(v1,v2) / (abs_v1 * abs_v2)
c
      return
      end
c---------------------------------------------------
      real function v3len(v1,v2)
c
c input:  v1, v2
c output: v3len
c---------------------------------------------------
      real v1(3),v2(3)
      real v3(3)
      external v3abs
c
      call v3sub(v1,v2,v3)
      v3len = v3abs(v3)
c
      return
      end

c---------------------------------------------------
      subroutine v3rot(v1,th,phi,v2)
c
c input:  v1,th,phi
c output: v2
c---------------------------------------------------
      real v1(3),v2(3)
      real th,phi
c
      v2(1) = cos(th)*cos(phi)*v1(1)
     &      - cos(th)*sin(phi)*v1(2)
     &      + sin(th)*cos(phi)*v1(3)

      v2(2) = cos(th)*sin(phi)*v1(1)
     &      + cos(th)*cos(phi)*v1(2)
     &      + sin(th)*sin(phi)*v1(3)

      v2(3) = - sin(th)*v1(1)
     &        - sin(th)*v1(2)
     &        + cos(th)*v1(3)
c
      end

c===================================================
c under here, double precision routine
c===================================================
c---------------------------------------------------
      subroutine dv3set(v,x,y,z)
c
c input:  v,x,y,z
c output: v
c---------------------------------------------------
      real*8 v(3),x,y,z
c
      v(1) = x
      v(2) = y
      v(3) = z
c
      end
c---------------------------------------------------
      subroutine dv3copy(v1,v2)
c
c input:  v1
c output: v2
c---------------------------------------------------
      real*8 v1(3),v2(3)
c
      v2(1) = v1(1)
      v2(2) = v1(2)
      v2(3) = v1(3)
c
      end
c---------------------------------------------------
      subroutine dv3add(v1,v2,v3)
c
c input:  v1, v2
c output: v3
c---------------------------------------------------
      real*8 v1(3),v2(3),v3(3)
c
      v3(1) = v1(1) + v2(1)
      v3(2) = v1(2) + v2(2)
      v3(3) = v1(3) + v2(3)
c
      end
c---------------------------------------------------
      subroutine dv3sub(v1,v2,v3)
c
c input:  v1, v2
c output: v3
c---------------------------------------------------
      real*8 v1(3),v2(3),v3(3)
c
      v3(1) = v1(1) - v2(1)
      v3(2) = v1(2) - v2(2)
      v3(3) = v1(3) - v2(3)
c
      end
c---------------------------------------------------
      subroutine dv3mean(v1,v2,v3)
c
c input:  v1, v2
c output: v3
c---------------------------------------------------
      real*8 v1(3),v2(3),v3(3)
c
      v3(1) = ( v1(1) + v2(1) )/2.
      v3(2) = ( v1(2) + v2(2) )/2.
      v3(3) = ( v1(3) + v2(3) )/2.
c
      end
c---------------------------------------------------
      subroutine dv3mult(c,v1,v2)
c
c input:  c, v1
c output: v2
c---------------------------------------------------
      real*8 c,v1(3),v2(3)
c
      v2(1) = c * v1(1)
      v2(2) = c * v1(2)
      v2(3) = c * v1(3)
c
      end
c---------------------------------------------------
      real*8 function dv3times(v1,v2)
c
c input:  v1, v2
c output: v3times
c---------------------------------------------------
      real*8 v1(3),v2(3)
c
      dv3times = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
c
      return
      end
c---------------------------------------------------
      real*8 function dv3abs(v)
c
c input:  v
c output: v3abs
c---------------------------------------------------
      real*8 v(3)
c
      if(.not.( v(1) .eq. 0. .and. v(2) .eq. 0. .and.
     &      v(3) .eq. 0. )) goto 23000
         dv3abs = 0.
         return
23000 continue
      dv3abs = sqrt( v(1)*v(1) + v(2)*v(2) + v(3)*v(3) )
c
      return
      end
c---------------------------------------------------
      subroutine dv3unit(v1,v2)
c
c input:  v1, v2
c output: v3unit
c---------------------------------------------------
      real*8 v1(3),v2(3)
      real*8 abs_v1
      external dv3abs
      real*8 dv3abs
c
      abs_v1 = dv3abs(v1)
c
      if(.not.( abs_v1 .eq. 0. ))goto 23002
         v2(1) = 0.
         v2(2) = 0.
         v2(3) = 0.
         return
23002 continue
      v2(1) = v1(1)/abs_v1
      v2(2) = v1(2)/abs_v1
      v2(3) = v1(3)/abs_v1
c
      end
c---------------------------------------------------
      real*8 function dv3cos(v1,v2)
c
c input:  v1, v2
c output: v3cos
c---------------------------------------------------
      real*8 v1(3),v2(3)
      external dv3times,dv3abs
      real*8 dv3times,dv3abs
      real*8 abs_v1,abs_v2
c
      abs_v1 = dv3abs(v1)
      abs_v2 = dv3abs(v2)
      if(.not.( abs_v1 * abs_v2 .eq. 0. ))goto 23004
         v3cos = 0.
         return
23004 continue
      dv3cos = dv3times(v1,v2) / (abs_v1 * abs_v2)
c
      return
      end
c---------------------------------------------------
      real*8 function dv3len(v1,v2)
c
c input:  v1, v2
c output: v3len
c---------------------------------------------------
      real*8 v1(3),v2(3)
      real*8 v3(3)
      external dv3abs
      real*8 dv3abs
c
      call dv3sub(v1,v2,v3)
      dv3len = dv3abs(v3)
c
      return
      end

c---------------------------------------------------
      subroutine dv3rot(v1,th,phi,v2)
c
c input:  v1,th,phi
c output: v2
c---------------------------------------------------
      real*8 v1(3),v2(3)
      real*8 th,phi
c
      v2(1) = cos(th)*cos(phi)*v1(1)
     &      - cos(th)*sin(phi)*v1(2)
     &      + sin(th)*cos(phi)*v1(3)

      v2(2) = cos(th)*sin(phi)*v1(1)
     &      + cos(th)*cos(phi)*v1(2)
     &      + sin(th)*sin(phi)*v1(3)

      v2(3) = - sin(th)*v1(1)
     &        - sin(th)*v1(2)
     &        + cos(th)*v1(3)
c
      end
