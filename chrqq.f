         subroutine chrqq
c
c           this is the 1-interval version that has been
c           determined to be the best on 8/29/02.
c
c
c        chiral nn potential in r-space space (= chrnn)
c        Bessel-transformed into momentum space (= chrqq).
c
c
      implicit real*8 (a-h,o-z)
c
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
c        arguments and values of this subroutine
c
      common /cpot/   v(6),xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep,label
c
c
c        this has been the end of the common-blocks containing
c        the arguments and values of this subroutine.
c
c        specifications for these two common blocks
c
      logical heform,sing,trip,coup,endep
c
c
c        `classical' r-space common blocks
      common /readw/kwrit1,kpunc1,kread1,kread2
      common /nquant/omega(5,23),iq(4),mlauf,index
      common /facen/ccf,e,x,xl,n
      common /cpotr1/vh(3),pquad(2)
c        this was the end of the `classical' r-space common blocks
c        note: the classical name for /cpotr1/ was /cpot/.
c
c        additional r-space common block
      common /cpotr2/vr(5),pqr(5),pqpr(5)
c
c        the following common block is the internal one of chrnn;
c        wn (= nucl. mass) is needed here.
      common /cchr/c(20,300),wn,tlamb,ga,fpi,
     1            cb1,cb2,cb3,cb4,cd12,cd3,cd5,cd145,
     2            ic(20,300),mgg(300),ime,
     3            indt(300)
c
c
      dimension vv(6),pq(6)
      dimension vrr(5,900),pqrr(5,900)
      dimension xg(900),wg(900)
      dimension xr(900),wr(900),xrr(900),xll(900)
      dimension xrrw(900),xrx(900),xry(900)
      dimension sbjx(900),sbpx(900),sbmx(900)
      dimension sbjy(900),sbpy(900),sbmy(900)
      dimension vl(4),adminv(4,4),ldminv(4),mdminv(4)
      character*4 name(3),nname(15)
c
      logical indind/.false./
      data jj/-1/
      data hbarc/197.32698d0/
      data pi/3.141592653589793d0/
      save
c
c
c
c
      if (indind) go to 30
      indind=.true.
c
      endep=.false.
c
c
      kwrit1=kwrite
      kread1=kread
      index=1
      mlauf=1
c
c
c        read in parameters for Bessel transformation
c
      write (kwrite,10000)
      read  (kread ,10001) name,nname
      write (kwrite,10002) name,nname
      read  (kread ,10003) name,nn
      write (kwrite,10004) name,nn
      read  (kread ,10005) name,ax,dx,ex
      write (kwrite,10006) name,ax,dx,ex
      read  (kread ,10005) name,qcut
      write (kwrite,10006) name,qcut
10000 format (//' chrqq: Bessel transformation of chiral r-space pot'
     1       /' ',50(1h-)/
     2        ' parameters used for Bessel transformation:'
     3       /' ',41(1h-))
10001 format (2a4,a2,15a4)
10002 format (' ',2a4,a2,15a4)
10003 format (2a4,a2,20i3)
10004 format (' ',2a4,a2,20i3)
10005 format (2a4,a2,6f10.4)
10006 format (' ',2a4,a2,6f10.4)
c
c
c        prepare integration points and weights
c        for Bessel transformation
c
      call gauss(0.d0,1.d0,nn,xg,wg)
c
      pih=pi/2.d0
      do 25 i=1,nn
      xx=pih*xg(i)
      xll(i) = dx*dtan(xx) + ax
      if (xll(i).gt.ex) go to 27
      dc=1.d0/dcos(xx)
      wr(i) = dx*pih*wg(i)*dc*dc/hbarc
      xr(i) = xll(i)/hbarc
      xrr(i) = xr(i)*xr(i)
   25 xrrw(i) = xrr(i)*wr(i)*2.d0/pi
   27 nnx=i-1
c
c
c
c
   30 if (j.eq.jj) go to 50
      jj=j
c
      jp=j+1
      jm=j-1
c
      iq(4)=j
      call omegam
c
c
c        call potential
c        and store potential matrix elements
c
c
      do 35 i=1,nnx
      xl=xll(i)
c
      call chrnn
c
      do 35 iv=1,5
      vrr(iv,i)=vr(iv)
   35 pqrr(iv,i)=pqr(iv)*2.d0/wn
c
c
c     prepare transformation to helicity states in case
c     it will be requested
c
      if (j.eq.0) go to 50
      aj=dble(j)
      aj1=dble(j+1)
      a2j1=dble(2*j+1)
      aaj6=dsqrt(aj*aj1)
c
c        coefficient matrix for the translations into lsj formalism
c
      adminv(1,1)=aj1
      adminv(1,2)=aj
      adminv(1,3)=-aaj6
      adminv(1,4)=-aaj6
      adminv(2,1)=aj
      adminv(2,2)=aj1
      adminv(2,3)=aaj6
      adminv(2,4)=aaj6
      adminv(3,1)=aaj6
      adminv(3,2)=-aaj6
      adminv(3,3)=aj1
      adminv(3,4)=-aj
      adminv(4,1)=aaj6
      adminv(4,2)=-aaj6
      adminv(4,3)=-aj
      adminv(4,4)=aj1
c
c       inversion
c
      call dminv (adminv,4,deter,ldminv,mdminv)
c
c
c
   50 continue
      do 55 iv=1,6
   55 v(iv)=0.d0
c
c
c         cut off potential for large momenta
c         if you wish to do so
c
      if (qcut.ne.0.d0) then
      if (xmev.gt.qcut.or.ymev.gt.qcut) return
      end if
c
c
c
c        prepare arguments for Bessel functions
c
      do 115 i=1,nnx
      xrx(i)=xr(i)*xmev
  115 xry(i)=xr(i)*ymev
c
c        call Bessel functions
c
      do 125 i=1,nnx
      sbjx(i)=sbess(j,xrx(i))
      sbjy(i)=sbess(j,xry(i))
  125 continue
      do 135 i=1,nnx
      sbpx(i)=sbess(jp,xrx(i))
      sbpy(i)=sbess(jp,xry(i))
  135 continue
c
      if (jm.ge.0) then
      do 145 i=1,nnx
      sbmx(i)=sbess(jm,xrx(i))
      sbmy(i)=sbess(jm,xry(i))
  145 continue
      else
      do 155 i=1,nnx
      sbmx(i)=0.d0
      sbmy(i)=0.d0
  155 continue
      end if
c
c
c        integrate
c        ---------
c
      do 205 iv=1,6
      vv(iv)=0.d0
  205 pq(iv)=0.d0
c
      do 215 i=1,nnx
c
c
      vv(1)=vv(1)+xrrw(i)*sbjx(i)*vrr(1,i)*sbjy(i)
      vv(2)=vv(2)+xrrw(i)*sbjx(i)*vrr(2,i)*sbjy(i)
      vv(3)=vv(3)+xrrw(i)*sbpx(i)*vrr(4,i)*sbpy(i)
      vv(4)=vv(4)+xrrw(i)*sbmx(i)*vrr(3,i)*sbmy(i)
      vv(5)=vv(5)+xrrw(i)*sbpx(i)*vrr(5,i)*sbmy(i)
      vv(6)=vv(6)+xrrw(i)*sbmx(i)*vrr(5,i)*sbpy(i)
c
c
      pq(1)=pq(1)+xrrw(i)*sbjx(i)*pqrr(1,i)*sbjy(i)
      pq(2)=pq(2)+xrrw(i)*sbjx(i)*pqrr(2,i)*sbjy(i)
      pq(3)=pq(3)+xrrw(i)*sbpx(i)*pqrr(4,i)*sbpy(i)
      pq(4)=pq(4)+xrrw(i)*sbmx(i)*pqrr(3,i)*sbmy(i)
      pq(5)=pq(5)+xrrw(i)*sbpx(i)*pqrr(5,i)*sbmy(i)
      pq(6)=pq(6)+xrrw(i)*sbmx(i)*pqrr(5,i)*sbpy(i)
c
  215 continue
c
c
      facpq=(xmev*xmev+ymev*ymev)*0.5d0
c
      do 255 iv=1,6
  255 v(iv)=vv(iv)+facpq*pq(iv)
c
c
c
c
      if (j.eq.0.or..not.heform) go to 8900
c
c
c         translation into (combinations of) helicity states
c         if requested
c
c
      do 8505 i=1,4
 8505 vl(i)=v(i+2)
c
      do 8520 ii=1,4
      iii=ii+2
      v(iii)=0.d0
c
      do 8515 i=1,4
 8515 v(iii)=v(iii)+adminv(ii,i)*vl(i)
 8520 v(iii)=v(iii)*a2j1
c
c
c
c
 8900 continue
      return
      end
      function sbess (l,z)
c**** sphaerische besselfunctionen ****
c**** es muss sein l groesser gleich null und kleiner gleich
c**** me minus zwei ****
c**** c von k ist die doppelfakultaet von 2*k plus 1, wobei k
c**** gleich mue plus l ist ****
c**** cc von mue ist das produkt 2 ** mue mal mue fakultaet ****
c**** ccc von mue ist der koeffizient des mue-ten gliedes ****
c**** die genauigkeit ist nur kritisch im uebergangsbereich ****
c**** sie betraegt dort bis l gleich 18   14 stellen,
c**** bis l gleich 23   13 stellen
c**** bis l gleich 30 12 stellen ****
c**** der output dieses programms wurde verglichen mit  tables of sphe-
c**** rical besselfunctions mathem. tables projekt  national bureau of
c**** standards, columbia university press, die dortige tabelle ist
c**** 10-stellig, in diesen 10 stllen herrscht uebereinstimmung ****
c
c**** moeglicherweise enthaelt der programmteil einen fehler,
c**** der benutzt wird, wenn die dimension ueberschritten wird.
c**** mit me=102 wuerde das jedoch nur auswirkungen haben
c**** fuer j_l mit l groesser 100.
c
c
      implicit real*8 (a-h,o-z)
      data ke/1/,mcce/1/,ll/-1/
      data me/102/
      real*8 c(102),cc(102),ccc(102)
      data c/3.d0,101*0.d0/,cc/2.d0,101*0.d0/
      data tol/1.d-16/
c
c
c
      if (z.ne.0.d0) go to 2
c**** der fall z gleich null ****
      if (l.eq.0) go to 1
      sbess=0.d0
      return
    1 sbess=1.d0
      return
c
    2 if (l.eq.0) go to 20
c
c**** berechnung von 2*l plus 1 doppelfakultaet ****
      if (l.le.ke) go to 4
      ka=ke+1
      ke=l
      do 3 k=ka,ke
    3 c(k)=c(k-1)*dfloat(2*k+1)
c
c**** entscheidung welches rechenverfahren gewaehlt wird ****
    4 if (z.ge.c(l )**(1.d0/dfloat(l))) go to 20
c
c
c**** berechnung mit reihe ****
c
c
c**** nulltes glied und toleranz ****
      zq=z**2
      zzq=-zq
      sbess=1.d0/c(l)
      tolrel=tol*sbess
c
      if (l.eq.ll) go to 5
      ll=l
      ma=2
      go to 7
c
c**** rechnung, wenn koeffizienten schon vorhanden ****
c
    5 do 6 mue=1,mme
      s=zzq*ccc(mue)
      sbess=sbess+s
      if (dabs(s).lt.tolrel) go to 13
    6 zzq=zzq*(-zq)
      if ((l+mme).eq.me) go to 15
      ma=mme+1
      go to 10
c
c**** rechnung, wenn koeffizienten zugleich berechnet werden ****
c
    7 k=l+1
      if (k.le.ke) go to 9
      ka=ke+1
      ke=k
      do 8 kk=ka,ke
    8 c(kk)=c(kk-1)*dfloat(2*kk+1)
c
c**** erstes glied ****
    9 ccc(1)=1.d0/(c(k)*cc(1))
      sbess=sbess+zzq*ccc(1)
      zzq=zzq*(-zq)
c
c**** die weiteren glieder ****
   10 do 11 mue=ma,me
      if (mue.gt.mcce) cc(mue)=cc(mue-1)*2.d0*dfloat(mue)
      k=l+mue
      if (k.gt.ke) c(k)=c(k-1)*dfloat(2*k+1)
      ccc(mue)=1.d0/(c(k)*cc(mue))
      s=zzq*ccc(mue)
      sbess=sbess+s
      if (dabs(s).lt.tolrel) go to 12
      if (k.eq.me) go to 14
   11 zzq=zzq*(-zq)
c**** festhalten der indices und vorbereitung der fortsetzung ****
      ke=k
      mcce=me
      mme=me
      amue=dfloat(me)
      aa=dfloat(2*me+1)
      go to 16
c
c**** festhalten der indices und sprung zum schluss ****
   12 if (mue.gt.mcce) mcce=mue
      if (k.gt.ke) ke=k
      mme=mue
      go to 13
c
c**** festhalten der indices und vorbereitung der fortsetzung ****
   14 ke=me
      if (mue.gt.mcce) mcce=mue
      mme=mue
      amue=dfloat(mue)
      aa=dfloat(2*me+1)
      go to 16
c
c**** fortsetzung der rechnung bei ueberschreiten der dimension ****
C**** THIS PART OF THE CODE MAY BE WRONG **************
   15 amue=dfloat(mme)
      aa=dfloat(2*me+1)
c
   16 amue=amue+1.d0
      aa=aa+2.d0
      s=s*(-zq/(2.d0*aa*amue))
      sbess=sbess+s
      if (dabs(s).gt.tolrel) go to 16
c
c**** schlussrechnung ****
   13 sbess=sbess*z**l
      return
c
c
c
c**** rekursive berechnung ****
c
c
   20 zinv=1.d0/z
      sa=dsin(z)*zinv
      if (l.ne.0) go to 21
      sbess=sa
      return
   21 sbess=(sa-dcos(z))*zinv
      if (l.eq.1) return
      do 22 i=2,l
      ss=dfloat(2*i-1)*zinv*sbess-sa
      sa=sbess
   22 sbess=ss
c
      return
      end
c********************************************************************
c name:    dminv
c        programmbibliothek rhrz bonn        28/11/78       dminv
c                                            fortran iv     ibm 370/168
c
c purpose:
c
c invert a matrix
c
c usage:   call dminv (a,n,d,l,m)
c
c parameters:
c
c a:       input matrix, destroyed in computation and replaced by
c          resultant inverse.
c          double precision required.
c
c n:       order of matrix a
c
c d:       resultant determinant
c          double precision required.
c
c l:       work vector of length n
c
c m:       work vector of length n
c
c remarks: matrix a must be a general matrix
c
c method:
c
c the standard gauss-jordan method is used. the determinant
c is also calculated. a determinant of zero indicates that
c the matrix is singular.
c
c programs required:
c          none
c
c author:  ibm, ssp iii
c
c**********************************************************************
      subroutine dminv (a,n,d,l,m)
      implicit real*8 (a-h,o-z)
      dimension a(1),l(1),m(1)

c
c
c        search for largest element
c
      d=1.d0
      nk=-n
      do 80 k=1,n
      nk=nk+n
      l(k)=k
      m(k)=k
      kk=nk+k
      biga=a(kk)
      do 20 j=k,n
      iz=n*(j-1)
      do 20 i=k,n
      ij=iz+i
   10 if (dabs(biga)-dabs(a(ij)))  15,20,20
   15 biga=a(ij)
      l(k)=i
      m(k)=j
   20 continue
c
c        interchange rows
c
      j=l(k)
      if(j-k) 35,35,25
   25 ki=k-n
      do 30 i=1,n
      ki=ki+n
      hold=-a(ki)
      ji=ki-k+j
      a(ki)=a(ji)
   30 a(ji) =hold
c
c        interchange columns
c
   35 i=m(k)
      if(i-k) 45,45,38
   38 jp=n*(i-1)
      do 40 j=1,n
      jk=nk+j
      ji=jp+j
      hold=-a(jk)
      a(jk)=a(ji)
   40 a(ji) =hold
c
c        divide column by minus pivot (value of pivot element is
c        contained in biga)
c
   45 if(biga) 48,46,48
   46 d=0.d0
      return
   48 do 55 i=1,n
      if(i-k) 50,55,50
   50 ik=nk+i
      a(ik)=a(ik)/(-biga)
   55 continue
c
c        reduce matrix
c
      do 65 i=1,n
      ik=nk+i
      hold=a(ik)
      ij=i-n
      do 65 j=1,n
      ij=ij+n
      if(i-k) 60,65,60
   60 if(j-k) 62,65,62
   62 kj=ij-i+k
      a(ij)=hold*a(kj)+a(ij)
   65 continue
c
c        divide row by pivot
c
      kj=k-n
      do 75 j=1,n
      kj=kj+n
      if(j-k) 70,75,70
   70 a(kj)=a(kj)/biga
   75 continue
c
c        product of pivots
c
      d=d*biga
c
c        replace pivot by reciprocal
c
      a(kk)=1.d0/biga
   80 continue
c
c        final row and column interchange
c
      k=n
  100 k=(k-1)
      if(k) 150,150,105
  105 i=l(k)
      if(i-k) 120,120,108
  108 jq=n*(k-1)
      jr=n*(i-1)
      do 110 j=1,n
      jk=jq+j
      hold=a(jk)
      ji=jr+j
      a(jk)=-a(ji)
  110 a(ji) =hold
  120 j=m(k)
      if(j-k) 100,100,125
  125 ki=k-n
      do 130 i=1,n
      ki=ki+n
      hold=a(ki)
      ji=ki-k+j
      a(ki)=-a(ji)
  130 a(ji) =hold
      go to 100
  150 return
      end
