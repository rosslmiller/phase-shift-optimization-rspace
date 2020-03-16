subroutine gauss(ax,bx,n,x,w)
        implicit double precision (a-h,o-z)
        parameter(pi=3.1415926535897932d0)
        dimension x(n),w(n)

        eps = 1.d-15
        do 2 i=1,n
        z=dcos(pi*(i-.25d0)/(n+.5d0))
    1   z1=z-alegf(n,0,z,0)/alegf(n,1,z,0)*dsqrt((1.d0+z)*(1.d0-z))
        if (dabs(z1-z).gt.eps) then
        z=z1
        goto 1
        end if
        x(i) = z1 
    2   w(i) = 2.d0/(alegf(n,1,z1,0)**2)
c
c        at this point, x and w are the gauss points and weights 
c        for the interval (-1,+1); x is descending.
c
c        now, re-scale gauss points and weights for interval (ax,bx)
c        such that x is ascending.
c
      alpha=0.5d0*(ax+bx)
      beta=0.5d0*(bx-ax)
      do 3 j=1,n
      x(j)=alpha-beta*x(j)
    3 w(j)=beta*w(j)
c
        return
        end
