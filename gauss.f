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
cern	  c315	    version    05/03/68 alegf	     94 	       c
      double precision function alegf(l,m,y,norm)
c     norm = 0 , unnormalised legendre functions
c     norm = 1 , normalised legendre functions
      implicit double precision (a-h,o-z)
      ma = iabs(m)
      if (l - ma) 1,2,2
    1 alegf = 0.d0
      return
    2 if (m) 3,4,4
    3 ib = l - ma + 1
      ie = l + ma
      prod = 1.d0
      do 5 i = ib,ie
      a = i
    5 prod = prod*a
      fact = 1.0d0/prod
      go to 6
    4 fact = 1.d0
    6 mm = ma
      if (ma - 1) 7,7,8
    7 p0 = 1.0d0
      if (ma) 9,9,10
    9 p1 = y
      go to 11
   10 p0 = 0.0d0
      p1 = dsqrt(1.00 - y**2)
   11 if (l - 1) 12,13,14
   12 alegf = p0
      go to 23
   13 alegf = p1 * fact
      go to 23
   14 pnmn1 = p0
      pn = p1
      am = ma
      k = l - 1
      do 15 n = 1,k
      an = n
      pnpl1 = (1.d0/(an-am+1.d0))*((2.d0*an+1.d0)*y*pn-(an+am)*pnmn1)
      pnmn1 = pn
   15 pn = pnpl1
   16 alegf = pnpl1*fact
      go to 23
    8 z2 = (1.d0-y**2)
      am = ma
      ham = am/2.d0
      zham = z2**ham
      if (ma.eq.(l-1)) go to 17
      nb = l+1
      ne = 2*l
      prod = 1.d0
      do 18 ni = nb,ne
      ai = ni
   18 prod = prod*ai
      denom = 2.d0**l
      dnn = prod/denom
      if (ma.eq.l) go to 19
   17 ne = 2*l-1
      prod = 1.d0
      do 20 ni = l,ne
      ai = ni
   20 prod = prod*ai
      denom = 2.d0**(l-1)
      dnm1n = y*prod/denom
      if (ma.eq.(l-1)) go to 21
      me = l-1-ma
      do 22 mn = 1,me
      an = l
      am = l-1-mn
      dnm2n=(1.d0/((an-am)*(an+am+1.d0)))*
     *      (2.d0*(am+1.d0)*y*dnm1n-z2*dnn)
      dnn = dnm1n
   22 dnm1n = dnm2n
      alegf = dnm2n*zham*fact
      go to 23
   19 alegf = dnn*zham*fact
      go to 23
   21 alegf = dnm1n*zham*fact
   23 if (norm.eq.1) go to 24
      return
   24 b = l
      if (m) 25,26,25
   25 ms = -m/mm
      jb = l - mm + 1
      je = l + mm
      prod = 1.d0
      do 27 i = jb,je
      a = i
   27 prod = prod*a
      factor = 0.5d0*(2.d0*b+1.d0)*(prod**ms)
      go to 28
   26 factor = 0.5d0*(2.d0*b+1.d0)
   28 factor = dsqrt(factor)
      alegf = alegf*factor
      return
      end
