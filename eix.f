c*********************************************************************72
c
cc EIX computes the exponential integral Ei(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    10 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision EI, the function value.
c
     subroutine eix ( x, ei )
	  implicit none

      double precision ei
      double precision ga
      integer k
      double precision r
      double precision x

      if ( x .eq. 0.0D+00 ) then

        ei = -1.0D+300

      else if ( x .le. 40.0D+00 ) then

        ei = 1.0D+00
        r = 1.0D+00
        do k = 1, 100
          r = r * k * x / ( k + 1.0D+00 )**2
          ei = ei + r
          if ( abs ( r / ei ) .le. 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        ga = 0.5772156649015328D+00
        ei = ga + log ( x ) + x * ei

      else

        ei = 1.0D+00
        r = 1.0D+00
        do k = 1, 20
          r = r * k / x
          ei = ei + r
        end do
        ei = exp ( x ) / x * ei

      end if

      
      end
      
	  
	  
	  
