       implicit none

      double precision ei
      double precision ga
      integer k
      double precision r
      double precision x
	x=-5.0

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

      write *, 'result is',ei

      end
      
