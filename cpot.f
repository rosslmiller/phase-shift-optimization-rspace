c        call pot
c
      implicit real*8 (a-h,o-z)
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
c
c        arguments and values
c
      common /cpot/   v(6),xmev,ymev
      common /cpts/   q(97),cca,n1,ix,iy
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
      character*4 label
c        this has been the end of the arguments and values of the
c        potential subroutine
c
c
      dimension vv(6)
      integer hel
      data lsj/'lsj '/,hel/'hel '/
      logical switch
      data switch/.false./
c
c
c
c
      open (unit=5,file='dpot.d')
      open (unit=6,file='pot.d')
      open (unit=7,file='pchpot.d')
c
c
c
c        set parameters
c
c
c
      kread=5
      kwrite=6
c        the following parameter is actually not needed
      kpunch=7
c
c
c        the following statement implies that the potential
c        will be given in the lsj-formalism
      heform=.false.
      name=lsj
      if (heform) name=hel
c        the following three statements mean, that for each j the
c        potential will be given for the states singlet, triplet
c        and the four coupled cases
      sing=.true.
      trip=.true.
      coup=.true.
c
c
      q(1)=10.d0
      q(2)=100.d0
      q(3)=1000.d0
c
      n=3
c
      n1=n+1
      q(n1)=q0mev
c
c
c
c
      do 595 j1=1,7
      j=j1-1
      write (kwrite,10000)
10000 format (/)
c
      do 495 ix=1,n
      xmev=q(ix)
c
c
      do 395 iy=1,n
c**** do 395 iy=ix,ix
      ymev=q(iy)
c
c
c
c
      call chrqq
c
c
c
      if (switch) go to 350
      switch=.true.
      write (kwrite,10001) name
10001 format 
     1 (//'    j      x         y          ',7x,a4,'- s t a t e s'/
     2 1h ,78(1h-)//)
c
c
  350 continue
c
c
c
      write (kwrite,10000)
c
      write (kwrite,10002) j,xmev,ymev
10002 format (' ',i5,2d15.6)
c
c
      do 385 iv=1,6
      vv(iv)=v(iv)
  385 continue
c
  395 write (kwrite,10004) vv
10004 format (1h ,6d13.6)
c
c
  495 continue
c
c
  595 continue
c
c
      stop
      end
