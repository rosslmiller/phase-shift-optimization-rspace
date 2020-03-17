
       implicit none

       real ( kind = 8 ) a0
       complex ( kind = 8 ) ce1
       complex ( kind = 8 ) cr
       complex ( kind = 8 ) ct
       complex ( kind = 8 ) ct0
       real ( kind = 8 ) el
       integer ( kind = 4 ) k
       real ( kind = 8 ) pi
       real ( kind = 8 ) x
       complex ( kind = 8 ) z
          z=5.0
       pi = 3.141592653589793D+00
       el = 0.5772156649015328D+00
       x = real ( z, kind = 8 )
       a0 = abs ( z )

       if ( a0 == 0.0D+00 ) then
       ce1 = cmplx ( 1.0D+300, 0.0D+00, kind = 8 )
       else if ( a0 <= 10.0D+00.or.(x<0.0D+00 .and. a0 < 20.0))then
       ce1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
       cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
       do k = 1, 150
       cr = - cr * k * z / ( k + 1.0D+00 )**2
       ce1 = ce1 + cr
       if ( abs ( cr ) <= abs ( ce1 ) * 1.0D-15 ) then
        exit
       end if
       end do

       ce1 = - el - log ( z ) + z * ce1

       else

       ct0 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
       do k = 120, 1, -1
        ct0 = k / ( 1.0D+00 + k / ( z + ct0 ) )
       end do
       ct = 1.0D+00 / ( z + ct0 )

       ce1 = exp ( - z ) * ct
       if ( x <= 0.0D+00 .and. imag ( z ) == 0.0D+00 ) then
       ce1 = ce1 - pi * cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
       end if

       end if

        write *, 'result is',-ce1
       end