c**** this package is consistently in double precision ******
c**** July 20, 1993
c
c**** names of new IMSL routines have been inserted on 11/5/93.
c
c
      subroutine deuter (pot,s,us,usf,a,b,aa,z,zz,l,m)
c
c
c
c         program computes the deuteron properties for a given momentum
c         space nn-potential
c
c
c         author:  r. machleidt
c                  institut fuer theoretische kern-physik der
c                  universitaet bonn
c                  nussallee 14-16
c                  d - 5300  bonn
c                  w. germany
c
c
c         pot is the name of the potential subroutine
c         dimension:
c         the dimension of s,us,usf,z and zz is n;
c         dimension a(2n*2n),aa(2*n*(n+1));
c         b, z and l and m have the dimension 2*n.
c
c
c         input-data:
c         n - number of gausspoints
c         iprop=1 n0n-relativistic propagator
c              =2 relativistic propagator
c         ipotfc=0: no special factor
c         ipotfc.ne.0: potential gets an additional factor m/e
c         iwavfc=0: no special factor
c         iwavfc.ne.0: wavefunctions get an additional factor
c                      dsqrt(e/m)
c         c - factor for transformation of gausspoints
c         wn - mass of nucleon
c         ee: experimental binding energy of the deuteron
c         deltae: energy variation in first iteration
c         tole: tolerance for final result of binding energy
c         ipunch=0 - no punching of wave functions
c         ipunch.ne.0 - punching of wave functions
c
c
      implicit real*8 (a-h,o-z)
      external pot
c
c        common blocks
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
c        arguments of the potential subroutine pot being called in this
c        program
c
      common /cpot/   v(6),xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep,label
      common /cpoted/ q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,
     1                znrl,zrel,smev,noced
      common /cpts/   q(97),c,n,ix,iy
c        specifications for these arguments
      logical heform,sing,trip,coup,endep
      logical noced
c
c
c        further specifications
c
c
      dimension s(1),us(1),usf(1),a(1),b(1),aa(1),z(1),zz(1),l(1),m(1)
      dimension r(6),ur(6),wr(6)
      data pih/1.570796326794897d0/
      data uf/197.327d0/
      real*4 eps
      data eps/1.0e-20/
      character*4 name(3),nname(17)
      character*4 charn,blanks
      data charn/'n   '/,blanks/'    '/
      logical inded
      logical indpch,indpts
      logical switch
      data inded/.false./
      data indpch/.false./,indpts/.false./
      data switch/.false./
      logical inddeu
c
c
c
c
c
10000 format (2a4,a2,20i3)
10001 format (1h ,2a4,a2,20i3)
10002 format (2a4,a2,6f10.3)
10003 format (1h ,2a4,a2,6f10.3)
10004 format (1h //' input-parameters for deuter'/1h ,27(1h-))
10005 format (//' transf gauss pts and wghts for c =',f8.2,
     1', cut =',f8.2,' and n =',i3/' points')
10006 format (7x,4f15.4)
10007 format (' weights')
10008 format (17a4)
10009 format (1h ,17a4)
10010 format (2a4,a2,3f20.10)
10011 format (1h ,2a4,a2,3f20.10)
10012 format (2a4,a2,10i6)
10013 format (1h ,2a4,a2,10i6)
10014 format (2a4,a2,6f10.5)
10015 format (1h ,2a4,a2,6f10.5)
10020 format (1h //' deuteron properties:'/1h ,20(1h-)/1h ,28x,'binding
     1energy (mev)',15x,'determinant'//30x,f14.8,d32.6)
10021 format (30x,f14.8,d32.6)
10022 format (1h /' final result',f31.8,d32.6)
10023 format (1h ///' error in deuter: matrix-inversion-error index =',
     1 i4///)
10024 format (1h /' control of omi. row:   sum',d19.6/' largest number o
     1f omi. row',d19.6)
10025 format (1h //34x,'wavefunction'//
     1 7x,'q(mev)',13x,'psi 0',18x,'psi 2',15x,'rho'/)
10026 format (1h ,d15.6,2d23.13,d15.6)
10027 format (1h /' deuteron d-state probability',f20.13)
10028 format (1h //' deuteron s-state probability',f20.13)
10029 format (1h /' sum',25x,f20.13)
10041 format (4d20.13)
10040 format ('deuteron wave functions in q-space for: ',a4,
     1' c =',f6.0,' n =',i3)
10030 format (1h /' quadrupolmoment of deuteron ',f15.8)
10050 format (1h /' asymptotic s- and d/s-state ',2f15.8,'   at r =',
     1 f6.2)
c
c
c
c
      write (kwrite,10004)
      read  (kread ,10008) nname
      write (kwrite,10009) nname
      read  (kread ,10000) name,n,nc
      write (kwrite,10001) name,n,nc
      read  (kread ,10000) name,iprop
      write (kwrite,10001) name,iprop
      read  (kread ,10000) name,ipotfc
      write (kwrite,10001) name,ipotfc
      read  (kread ,10000) name,iwavfc
      write (kwrite,10001) name,iwavfc
      read  (kread ,10000) name,iendep
      write (kwrite,10001) name,iendep
      read  (kread ,10002) name,c,cut
      write (kwrite,10003) name,c,cut
      read  (kread ,10014) name,wn1,wn2
      write (kwrite,10015) name,wn1,wn2
      read  (kread ,10010) name,ee
      write (kwrite,10011) name,ee
      read  (kread ,10010) name,deltae
      write (kwrite,10011) name,deltae
      read  (kread ,10010) name,tole
      write (kwrite,10011) name,tole
      read  (kread ,10000) name,ipoint
      write (kwrite,10001) name,ipoint
      read  (kread ,10000) name,ipunch
      write (kwrite,10001) name,ipunch
      read  (kread ,10000) name,nr
      write (kwrite,10001) name,nr
      read  (kread ,10002) name,(r(i),i=1,nr)
      write (kwrite,10003) name,(r(i),i=1,nr)
      read  (kread ,10002) name,h1,h2,qlin
      write (kwrite,10003) name,h1,h2,qlin
      read  (kread ,10012) name,iiee
      write (kwrite,10013) name,iiee
c
c
c
c
c        preparation of constants
c
c
      if (iendep.ne.0) inded=.true.
      if (ipoint.ne.0) indpts=.true.
      if (ipunch.ne.0) indpch=.true.
c
c
c        calculate the reduced mass for the two nucleon masses given
c        and define the nucleon mass to be used in this code as
c        twice the reduced mass
c
      wn=2.d0*wn1*wn2/(wn1+wn2)
c
c
      wnq=wn*wn
      wnr=1.d0/wn
      wn2=2.d0*wn
      dwn2=1.d0/wn2
c
      wurz2=dsqrt(2.d0)
      wurz6=dsqrt(6.d0)
      wurz6r=1.d0/wurz6
      dpi4=1.d0/(8.d0*pih)
c
      f=0.25d0*dsqrt(2.d0)
      ff=uf*uf*0.1d0*dsqrt(2.d0)
c
      j=1
      heform=.true.
      sing=.false.
      trip=.false.
      coup=.true.
      qfmev=0.d0
      pmev=0.d0
      uanspp=0.d0
      wsnspp=wn
      ucnspp=0.d0
      udnspp=0.d0
c
c
c        get gauss points and weights
c
      nn=n
      if (cut.eq.0.d0) then
      acut=1.d0
      else
      acut=cut
      endif
      call gset(0.d0,acut,nn,q,s)
      n=nn
      if (nc.eq.0) nc=n
      nq=n*n
      nq2=2*nq
      nm1=n-1
      n2=2*n
      nq2pn=nq2+n
      n2p1=n2+1
      n2m1=n2-1
      nq2mn=nq2-n
c
c        transform gauss points und weights
c
      do 209 i=1,n
      if (cut.eq.0.d0) then
      xx=pih*q(i)
      q(i)=c*dtan(xx)
      s(i)=pih*c/dcos(xx)**2*s(i)
      end if
      qq=q(i)*q(i)
      us(i)=qq*s(i)
      usf(i)=us(i)
      go to (201,202),iprop
  201 z(i)=qq*wnr
      go to 209
  202 z(i)=2.d0*dsqrt(wnq+qq)-wn2
  209 continue
      if (.not.indpts) go to 301
      write (kwrite,10005) c,cut,n
      write (kwrite,10006) (q(i),i=1,n)
      write (kwrite,10007)
      write (kwrite,10006) (s(i),i=1,n)
      go to 301
c
c
c
c
c        iterations for energy
c
c
  300 da=dd
      if (.not.inded.or..not.endep) go to 400
c
  301 q0qmev=-ee*wn
      znrl=-ee
      zrel=wn2-ee
      smev=wn2-ee
      noced=.false.
c
      zrelh=zrel*0.5d0
c
c
c        get potential matrix
c
c
      iii=0
      do 331 ix=1,n
      xmev=q(ix)
c
c
      do 331 iy=ix,n
      iaa=iii*4
      ymev=q(iy)
c
c
      call pot
c
c
      do 321 iv=1,4
      ivv=iv+2
  321 aa(iv+iaa)=v(ivv)
c
  331 iii=iii+1
c
c
c        build up matrix for determinant
c
c
  400 do 415 i=1,n
      zz(i)=1.d0/(z(i)+ee)
      if (ipotfc.eq.0) go to 415
      eq=dsqrt(wnq+q(i)*q(i))
      fax=dsqrt((eq+zrelh)*dwn2)
      zz(i)=zz(i)*fax
      usf(i)=us(i)*fax
  415 continue
c
      iii=0
      do 431 i=1,n
      in=(i-1)*n2
      do 431 ii=i,n
c
      iaa=iii*4
      iaa1=iaa+1
      iaa2=iaa+2
      iaa3=iaa+3
      iaa4=iaa+4
      iin=(ii-1)*n2
      il=iin+i
      iil=in+ii
      u=zz(i)*usf(ii)
      uu=zz(ii)*usf(i)
c
c
      a(il)=aa(iaa1)*u
      a(iil)=aa(iaa1)*uu
      a(il+nq2pn)=aa(iaa2)*u
      a(iil+nq2pn)=aa(iaa2)*uu
      a(il+nq2)=aa(iaa3)*u
      a(iil+n)=aa(iaa3)*uu
      a(il+n)=aa(iaa4)*u
      a(iil+nq2)=aa(iaa4)*uu
c
c
  431 iii=iii+1
c
c
      do 441 i=1,n2
      i1=i+(i-1)*n2
  441 a(i1)=a(i1)+1.d0
c
c
c        compute determinant
c
c
      call dminv(a,n2,dd,l,m)
c
c
c        write energy and determinant
c        and guess energy for next iteration
c
      if (switch) go to 550
      switch=.true.
      write (kwrite,10020) ee,dd
      ea=ee
      ee=ee-dsign(deltae,dd)
      go to 300
c
c
  550 deltae=(ee-ea)/(dd-da)*dd
      if (dabs(deltae).lt.tole) go to 590
      write (kwrite,10021) ee,dd
      ea=ee
      ee=ee-deltae
      go to 300
c
c
c        write final result for energy
c
  590 write (kwrite,10022) ee,dd
c
c
c        compute wave functions
c
c
c        build up matrix to be inverted
c
c
      iii=0
      do 649 i=1,n
      ins=(i-1)*n2m1
      do 649 ii=i,n
c
      iaa=iii*4
      iaa1=iaa+1
      iaa2=iaa+2
      iaa3=iaa+3
      iaa4=iaa+4
      iins=(ii-1)*n2m1
      ils=iins+i
      iils=ins+ii
      u=zz(i)*usf(ii)
      uu=zz(ii)*usf(i)
c
      a(ils)=aa(iaa1)*u
      a(iils)=aa(iaa1)*uu
c
      if (i.eq.nc) go to 620
      a(ils+n)=aa(iaa4)*u
      a(iils+nq2mn)=aa(iaa4)*uu
c
  620 if (ii.eq.nc) go to 625
      a(ils+nq2mn)=aa(iaa3)*u
      a(iils+n)=aa(iaa3)*uu
      go to 630
  625 b(i)=-aa(iaa3)*u
      z(i)=aa(iaa3)*uu
c
  630 if (ii.eq.nc) go to 635
      if (i.eq.nc) go to 649
      a(ils+nq2)=aa(iaa2)*u
      a(iils+nq2)=aa(iaa2)*uu
      go to 649
  635 if (i.eq.nc) go to 640
      b(n+i)=-aa(iaa2)*u
      z(n+i)=aa(iaa2)*uu
      go to 649
  640 znc=aa(iaa2)*uu+1.d0
c
  649 iii=iii+1
c
c
      do 659 i=1,n2m1
      i1=i+(i-1)*n2m1
  659 a(i1)=a(i1)+1.d0
c
c
c        solve system of linear equations
c
c
      call dgelg(b,a,n2m1,1,eps,ier)
c
      if (ier.ne.0) write (kwrite,10023) ier
c
c
c        prepare check of omitted line
c
      nd=n-nc
      if (nd.eq.0) go to 687
      do 685 i=1,nd
      i1=n2p1-i
      z(i1)=z(i1-1)
  685 b(i1)=b(i1-1)
  687 b(n+nc)=1.d0
      z(n+nc)=znc
c
c
c        check of omitted line
c
      akont=0.d0
      grsum=0.d0
c
      do 719 i=1,n2
c
      zb=z(i)*b(i)
      akont=akont+zb
      if(grsum-dabs(zb)) 718,719,719
  718 grsum=dabs(zb)
  719 continue
c
      write (kwrite,10024) akont,grsum
c
c
c        normalize wave function
c
c
      wnormq=0.d0
      do 725 i=1,n
      ni=n+i
      if (iwavfc.eq.0) go to 724
      eq=dsqrt(q(i)*q(i)+wnq)
      rz=wn2/(eq+zrelh)
      rz=dsqrt(rz)
      b(i)=b(i)*rz
      b(ni)=b(ni)*rz
  724 if (i.eq.n) go to 726
  725 wnormq=wnormq+0.5d0*us(i)*(b(i)*b(i)+b(ni)*b(ni))
  726 wnorma=1.d0/dsqrt(wnormq)
c
c
      do 729 i=1,n2
  729 b(i)=b(i)*wnorma
c
c
c        write wave function
c
c
      write (kwrite,10025)
c
      do 739 i=1,n
      ni=n+i
      z(i)=(b(i)+wurz2*b(ni))*wurz6r
      zz(i)=(-wurz2*b(i)+b(ni))*wurz6r
  739 usf(i)=(z(i)*z(i)+zz(i)*zz(i))*dpi4
c
c        sign of wave functions
      fz=1.d0
      fzz=1.d0
      if (z(3).lt.0.d0) fz=-1.d0
      if (zz(3).lt.0.d0) fzz=-1.d0
      do 745 i=1,n
      z(i)=fz*z(i)
      zz(i)=fzz*zz(i)
  745 write (kwrite,10026) q(i),z(i),zz(i),usf(i)
      z(n)=0.d0
      zz(n)=0.d0
c
c
c        punch wave functions, if requested
c
c
      if (.not.indpch) go to 750
      write (kpunch,10040) label,c,n
      name(1)=charn
      name(2)=blanks
      name(3)=blanks
      write (kpunch,10000) name,n
      write (kpunch,10041) (q(i),s(i),z(i),zz(i),i=1,n)
c
c
c        compute d-state probability
c
c
  750 pd=0.d0
      sint=0.d0
      do 759 i=1,nm1
      sint=sint+us(i)*z(i)*z(i)
  759 pd=pd+us(i)*zz(i)*zz(i)
c
      sum=sint+pd
      write (kwrite,10028) sint
      write (kwrite,10027) pd
      write (kwrite,10029) sum
c
c
c        compute quadrupolmoment
c
c
c         derive wave functions
      q12=q(1)-q(2)
      us(1)=(z(1)-z(2))/q12
      usf(1)=(zz(1)-zz(2))/q12
      q12=q(nm1)-q(n)
      us(n)=(z(nm1)-z(n))/q12
      usf(n)=(zz(nm1)-zz(n))/q12
      del01=us(1)
      del21=usf(1)
      do 763 i=2,nm1
      i1=i+1
      q12=q(i)-q(i1)
      del02=(z(i)-z(i1))/q12
      del22=(zz(i)-zz(i1))/q12
      us(i)=0.5d0*(del01+del02)
      usf(i)=0.5d0*(del21+del22)
      del01=del02
  763 del21=del22
c
c        integrate
      qdp=0.d0
      do 765 i=1,n
      qq=q(i)*q(i)
  765 qdp=qdp+s(i)*(qq*us(i)*usf(i)+3.d0*q(i)*zz(i)*us(i)
     1+f*(qq*usf(i)*usf(i)+6.d0*zz(i)*zz(i)))
c
      qdp=-ff*qdp
c
c        write quadrupolmoment
      write (kwrite,10030) qdp
c
c
c
c
c        asymptotic d/s
c
      if (nr.eq.0) go to 1000
      inddeu=.true.
      call deubtr (z,zz,r,ur,wr,h1,h2,qlin,iiee,nr,inddeu)
      alpha=dsqrt(ee*wn)/uf
      do 815 i=1,nr
      alphar=alpha*r(i)
      as=ur(i)*dexp(alphar)
      eta=wr(i)/(ur(i)*(1.d0+3.d0/alphar+3.d0/(alphar*alphar)))
  815 write (kwrite,10050) as,eta,r(i)
c
c
 1000 stop
      end
      subroutine deubtr (uq,wq,r,ur,wr,h1,h2,qlin,iiee,nr,inddeu)
c
c        besseltransformation of deuteron wave functions
c
c
      implicit real *8 (a-h,o-z)
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
      common /cpts/   q(97),cx,n,ix,iy
c
c
      dimension uq(1),wq(1),r(1),ur(1),wr(1)
      dimension qknot(194),bscoeu(97)
      dimension bscoew(97)
      data pih/1.570796326794897d0/
      data uf/197.327d0/
      character*4 name(3),nname(17)
      character*4 charnr,blanks
      data charnr/'nr  '/,blanks/'    '/
      logical inddeu
      logical indint
      logical indwrt
      data indwrt/.false./
c
c
c
10000 format (1h //' input-parameters for deubtr'/1h ,27(1h-))
10001 format (17a4)
10002 format (1h ,17a4)
10003 format (2a4,a2,20i3)
10004 format (1h ,2a4,a2,20i3)
10005 format (2a4,a2,6f10.2)
10006 format (1h ,2a4,a2,6f10.2)
10007 format (2a4,a2,10i6)
10008 format (1h ,2a4,a2,10i6)
10009 format (7f10.4)
10010 format (1h ,7f10.4)
10020 format (4d20.13)
10021 format (3d23.13)
10030 format ('deuteron wave functions in r-space for: ',7a4)
10040 format (1h ///' r-space wave functions'//)
10100 format (' iteration of integrals'/' l',3x,'r',7x,'ii',12x,
     1 'qint',16x,'psi(qint)',i2,15x,'phi(r)',i2/)
10101 format (i2,f7.3,i7,3d22.13)
c
c
c
c
      iex=0
      if (inddeu) go to 50
c
      write (kwrite,10000)
      read  (kread ,10001) nname
      write (kwrite,10002) nname
      read  (kread ,10003) name,nr
      write (kwrite,10004) name,nr
      read  (kread ,10009) (r(i),i=1,nr)
      write (kwrite,10010) (r(i),i=1,nr)
c        linear interpolation starts in the next interval after
c        qlin; it is advisable to have linear interpolation for
c        the last two intervalls (an intervall is the range between
c        two gausspoints of the q-integration)'
c        further recommended values for high accuracy (in brackets
c        a sufficient accuracy is given):
c        h1 = 72. (48.)
c        h2 = 5.  (10.)
c        iiee = 6000  (3000)
      read  (kread ,10005) name,h1,h2,qlin
      write (kwrite,10006) name,h1,h2,qlin
      read  (kread ,10007) name,iiee
      write (kwrite,10008) name,iiee
      read  (kread ,10003) name,iwrite
      write (kwrite,10004) name,iwrite
      read  (kread ,10003) name,ipunch
      write (kwrite,10004) name,ipunch
c
      if (iwrite.ne.0) indwrt=.true.
c
c
c        read in q-space deuteron wave function from data set
c
      kda1=kda(1)
      read  (kda1  ,10001) nname
      write (kwrite,10002) nname
      read  (kda1  ,10003) name,n
      write (kwrite,10004) name,n
      read  (kda1  ,10020) (q(i),a,uq(i),wq(i),i=1,n)
      write (kwrite,10021) (q(i),uq(i),wq(i),i=1,n)
c
c
c
c
   50 pi2=4.d0*pih
      anh1=uf*pi2/h2
      pihr=1.d0/pih
      wpihr=dsqrt(pihr)
      wuf=dsqrt(uf)
      wufr=1.d0/wuf
      c2d3=2.d0/3.d0
      c4d3=2.d0*c2d3
      nnh1=h1
      qn=q(n)
      qnm1=q(n-1)
      qnm2=q(n-2)
      if (qnm1.lt.qlin) qnm1=qn
      if (qnm2.lt.qlin) qnm2=qn
c
c
c
c
c        l-loop
c
c
c
c
      korder=5
      call dbsnak (n,q,korder,qknot)
c
c
      do 595 k=1,3,2
c
      l=k-1
c
      if (indwrt)
     1write (6,10100) l,l
c
c
c        prepare interpolation
      if (l.eq.2) go to 120
      call dbsint (n,q,uq,korder,qknot,bscoeu)
      psinm1=uq(n-1)
      psinm2=uq(n-2)
      if (qn.eq.qnm1) go to 115
      ss1=(uq(n)-uq(n-1))/(qn-qnm1)
  115 if (qnm2.eq.qnm1) go to 200
      ss2=(uq(n-1)-uq(n-2))/(qnm1-qnm2)
      go to 200
  120 continue
      call dbsint (n,q,wq,korder,qknot,bscoew)
      psinm1=wq(n-1)
      psinm2=wq(n-2)
      if (qn.eq.qnm1) go to 125
      ss1=(wq(n)-wq(n-1))/(qn-qnm1)
  125 if (qnm2.eq.qnm1) go to 200
      ss2=(wq(n-1)-wq(n-2))/(qnm1-qnm2)
c
c
c
c
c        loop of arguments r
c
c
c
c
  200 do 585 i=1,nr
c
c
      if (indwrt) write (kwrite,10002)
c
c
      phil=0.d0
      if (r(i).eq.0.d0) go to 550
c
c
c
c        nh1 is the number of integration pionts for an interval of
c        length 2pi of the argument of the bessel-function
      nh1=anh1/r(i)
      nh1=(nh1/2)*2
      if (nh1.lt.nnh1) nh1=nnh1
      h=pi2/dfloat(nh1)
      iie=(iiee/nh1+1)*nh1
      iie=iie+nh1/4
      iiex=iie-10*nh1
      ufdrwf=uf/r(i)
c
c
      hint=0.d0
      indint=.false.
c
c
c
c
c        integration loop of besseltransformation
c
c
c
c
      do 495 ii=1,iie
c
      hint=hint+h
      qint=hint*ufdrwf
      qintq=qint*qint
c
      if (indint) go to 306
      gew=c4d3
      indint=.true.
      go to 307
  306 gew=c2d3
      indint=.false.
  307 continue
c
c
c
c
c        interpolate
c
c
      if (qint.ge.qnm1) go to 315
      if (qint.lt.qnm2) go to 320
c
c
c        interpolate linearly for the second but last interval
c
      psil=psinm2+ss2*(qint-qnm2)
      go to 340
c
c
c        interpolate linearly for the last interval
c
  315 if (qint.ge.qn) go to 317
      psil=psinm1+ss1*(qint-qnm1)
      go to 340
  317 psil=0.d0
      go to 340
c
c
c        interpolate by spline
c
  320 if (l.eq.2) go to 330
      psil=dbsval(qint,korder,qknot,n,bscoeu)
      go to 340
  330 continue
      psil=dbsval(qint,korder,qknot,n,bscoew)
c
c
c
c
c        integrate
c
c
  340 phil=phil+qintq*sbess(l,hint)*psil*gew
c
c
c
c
      if (ii.lt.iiex) go to 495
      if (mod(ii,nh1).eq.0.and.indwrt)
     1write (6,10101) l,r(i),ii,qint,psil,phil
  495 continue
c        this has been the end of the integration loop of the
c        bessel-transformation
c
c
c
  550 if (l.eq.2) go to 560
      ur(i)=phil*h*wpihr*wufr
      go to 585
  560 wr(i)=phil*h*wpihr*wufr
c
  585 continue
c        this has been the end of the loop of arguments r
c
c
c
  595 continue
c        this has been the end of the l-loop
c
c
c
      if (ipunch.eq.0) go to 650
      write (kpunch,10030) (nname(i),i=11,17)
      name(1)=charnr
      name(2)=blanks
      name(3)=blanks
      write (kpunch,10003) name,nr
      write (kpunch,10021) (r(i),ur(i),wr(i),i=1,nr)
c
  650 if (inddeu) go to 1000
      write (kwrite,10040)
      write (kwrite,10021) (r(i),ur(i),wr(i),i=1,nr)
c
c
 1000 return
      end
      subroutine deucff
      implicit real*8 (a-h,o-z)
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
      common /ccff/ eb,wn,anorm,indren
      logical indren
      data pi/3.1415926535897932d0/
      integer name(3),nname(17)
      logical index
      data index/.false./
c
10000 format (1h //' input-parameters for deucff'/1h ,27(1h-))
10001 format (17a4)
10002 format (1h ,17a4)
10003 format (2a4,a2,20i3)
10004 format (1h ,2a4,a2,20i3)
10005 format (2a4,a2,6f10.4)
10006 format (1h ,2a4,a2,6f10.4)
10007 format (2a4,a2,3f20.10)
10008 format (1h ,2a4,a2,3f20.10)
10100 format (' '///'   q(fm-1)         cff',14x,'norm',
     1 10x,'pd',10x,'quad',10x,'rms'//)
10101 format (1h ,f10.3,d16.6,f18.13,3f12.8)
      write (kwrite,10000)
      read  (kread ,10001) nname
      write (kwrite,10002) nname
      read  (kread ,10003) name,nh
      write (kwrite,10004) name,nh
      read  (kread ,10005) name,re
      write (kwrite,10006) name,re
      read  (kread ,10005) name,dq
      write (kwrite,10006) name,dq
      read  (kread ,10005) name,qe
      write (kwrite,10006) name,qe
      read  (kread ,10007) name,eb
      write (kwrite,10008) name,eb
      read  (kread ,10005) name,wn
      write (kwrite,10006) name,wn
      read  (kread ,10003) name,iren
      write (kwrite,10004) name,iren
      indren=.false.
      cq1=dsqrt(2.d0)
      cq2=cq1/4.d0
      cq1=cq1/10.d0
      q=0.d0
      iqe=qe/dq
      do 200 iq=1,iqe
      q=q+dq
      qq=q*q
      qx=q
      if (q.lt.1.d0) qx=1.d0
      hr=2.d0*pi/(qx*dfloat(nh))
      nh2=2*nh
      ie=re/hr
      ie=ie/nh2
      ie=(ie+1)*nh2
      r=0.d0
      cff=0.d0
      fakt=2.d0
      if (index) go to 90
      anorm=0.d0
      quad=0.d0
      pd=0.d0
      rms=0.d0
   90 do 100 i=1,ie
      if(fakt.ne.4.d0) go to 111
      fakt=2.d0
      go to 112
  111 fakt=4.d00
  112 continue
      r=r+hr
      call dwave(r,u,w)
c
c
c**** prepare integrand ****
c
      ww=w*w
      ai=u*u+ww
      aif=ai*fakt
      xi=q*r*0.5d0
c
c
c**** integrals ****
c
      cff=cff+aif*dsin(xi)/xi
c
      if (index) go to 100
      anorm=anorm+aif
      pd=pd+fakt*ww
      rr=r*r
      quad=quad+fakt*rr*(u*w-cq2*ww)
      rms=rms+rr*aif
c
  100 continue
c
      hr3=hr/3.d0
      cff=cff*hr3
c
      if (index) go to 250
      index=.true.
      anorm=anorm*hr3
      pd=pd*hr3
      quad=quad*hr3*cq1
      rms=dsqrt(rms*hr3)*0.5d0
c
      write (6,10100)
      write (6,10101) q,cff,anorm,pd,quad,rms
      go to 200
c
c
  250 write (6,10101) q,qq,cff
c
c
  200 continue
c
c
c
c
c        renormalize wave functions if requested
c
      if (iren.eq.0) go to 1000
      indren=.true.
      call dwave (r,u,w)
c
c
 1000 stop
      end
      subroutine dwave (r,u,w)
c
c**** deuteron wave functions in r-space ****
c
c
      implicit real*8 (a-h,o-z)
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
      common /ccff/ eb,wn,anorm,indren
      logical indren
c
      real*8 ra(100),ua(100),wa(100)
      dimension rknot(200),bscoeu(100)
      dimension bscoew(100)
      data uf/197.327d0/
      character*4 name(3),nname(17)
      character*4 reno,charnr,blanks
      data reno/'reno'/,charnr/'nr  '/,blanks/'    '/
      logical index
      data index/.false./
c
c
c**** statement functions ****
c
      f3(b,x)=1.d0+3.d0/(b*x)+3.d0/(b*b*x*x)
c
c
c
c
10001 format (17a4)
10002 format (1h ,17a4)
10003 format (2a4,a2,20i3)
10004 format (1h ,2a4,a2,20i3)
10020 format (3d23.13)
10021 format (1h ,3d23.13)
c
c
c
c
      if (index) go to 100
      index=.true.
c
c
c        read in r-space deuteron wave function from data set
c
      kda1=kda(1)
      read  (kda1  ,10001) nname
      write (kwrite,10002) nname
      read  (kda1  ,10003) name,nr
      write (kwrite,10004) name,nr
      read  (kda1  ,10020) (ra(i),ua(i),wa(i),i=1,nr)
      write (kwrite,10021) (ra(i),ua(i),wa(i),i=1,nr)
c
c
c
c
c        compute asymptotic normalizations
c
c
      na=nr-5
      alpha=dsqrt(eb*wn)/uf
      write (6,10100)
10100 format ('1'///5x,'r(fm)',11x,'as',16x,'eta'//)
      do 10 i=na,nr
      alphar=alpha*ra(i)
      as=ua(i)*dexp(alphar)
      eta=wa(i)/(ua(i)*f3(alpha,ra(i)))
   10 write (6,10101) ra(i),as,eta
10101 format (1h ,f10.3,2f18.8)
c
c
c        prepare interpolation
c
      korder=5
      call dbsnak (nr,ra,korder,rknot)
      call dbsint (nr,ra,ua,korder,rknot,bscoeu)
      call dbsint (nr,ra,wa,korder,rknot,bscoew)
      ran=ra(nr)
c
c
  100 if (indren) go to 200
c
c
c**** interpolate wave functions for r ****
c
c
      if (r.gt.ran) go to 151
c
      u=dbsval(r,korder,rknot,nr,bscoeu)
      w=dbsval(r,korder,rknot,nr,bscoew)
c
      go to 152
c
c
  151 eex=1.d0/dexp(alpha*r)
      u=as*eex
      w=eta*f3(alpha,r)*u
c
c
  152 return
c
c
c
c
c        renormalize wave functions
c
  200 dnorm=1.d0/dsqrt(anorm)
      do 215 i=1,nr
      ua(i)=ua(i)*dnorm
  215 wa(i)=wa(i)*dnorm
      nname(17)=reno
      write (kpunch,10001) nname
      name(1)=charnr
      name(2)=blanks
      name(3)=blanks
      write (kpunch,10003) name,nr
      write (kpunch,10020) (ra(i),ua(i),wa(i),i=1,nr)
c
c
c
      stop
c
c
      end
      subroutine dgelg(r,a,m,n,eps,ier)
c
c
      implicit real*8 (a-h,o-z)
      dimension a(1),r(1)
      real*4 eps
c
c
c
c
      if(m)23,23,1
c
c     search for greatest element in matrix a
    1 ier=0
      piv=0.d0
      mm=m*m
      nm=n*m
      do 3 l=1,mm
      tb=dabs(a(l))
      if(tb-piv)3,3,2
    2 piv=tb
      i=l
    3 continue
      tol=eps*piv
c     a(i) is pivot element. piv contains the absolute value of a(i).
c
c
c     start elimination loop
      lst=1
      do 17 k=1,m
c
c     test on singularity
      if(piv)23,23,4
    4 if(ier)7,5,7
    5 if(piv-tol)6,6,7
    6 ier=k-1
    7 pivi=1.d0/a(i)
      j=(i-1)/m
      i=i-j*m-k
      j=j+1-k
c     i+k is row-index, j+k column-index of pivot element
c
c     pivot row reduction and row interchange in right hand side r
      do 8 l=k,nm,m
      ll=l+i
      tb=pivi*r(ll)
      r(ll)=r(l)
    8 r(l)=tb
c
c     is elimination terminated
      if(k-m)9,18,18
c
c     column interchange in matrix a
    9 lend=lst+m-k
      if(j)12,12,10
   10 ii=j*m
      do 11 l=lst,lend
      tb=a(l)
      ll=l+ii
      a(l)=a(ll)
   11 a(ll)=tb
c
c     row interchange and pivot row reduction in matrix a
   12 do 13 l=lst,mm,m
      ll=l+i
      tb=pivi*a(ll)
      a(ll)=a(l)
   13 a(l)=tb
c
c     save column interchange information
      a(lst)=j
c
c     element reduction and next pivot search
      piv=0.d0
      lst=lst+1
      j=0
      do 16 ii=lst,lend
      pivi=-a(ii)
      ist=ii+m
      j=j+1
      do 15 l=ist,mm,m
      ll=l-j
      a(l)=a(l)+pivi*a(ll)
      tb=dabs(a(l))
      if(tb-piv)15,15,14
   14 piv=tb
      i=l
   15 continue
      do 16 l=k,nm,m
      ll=l+j
   16 r(ll)=r(ll)+pivi*r(l)
   17 lst=lst+m
c     end of elimination loop
c
c
c     back substitution and back interchange
   18 if(m-1)23,22,19
   19 ist=mm+m
      lst=m+1
      do 21 i=2,m
      ii=lst-i
      ist=ist-lst
      l=ist-m
      l=a(l)+.5d0
      do 21 j=ii,nm,m
      tb=r(j)
      ll=j
      do 20 k=ist,mm,m
      ll=ll+1
   20 tb=tb-a(k)*r(ll)
      k=j+l
      r(j)=r(k)
   21 r(k)=tb
   22 return
c
c
c     error return
   23 ier=-1
      return
      end
      subroutine gset(ax,bx,n,z,w)
c
c
c        this code has been obtained from the CERN computer library
c        in the year of the lord 1972.
c
c
      implicit real*8 (a-h,o-z)
c
c     n-point gauss zeros and weights for the interval (ax,bx) are
c           stored in  arrays z and w respectively.
c
      dimension     a(273),x(273),ktab(96)
      dimension z(2),w(2)
c
c-----table of initial subscripts for n=2(1)16(4)96
      data ktab(2)/1/
      data ktab(3)/2/
      data ktab(4)/4/
      data ktab(5)/6/
      data ktab(6)/9/
      data ktab(7)/12/
      data ktab(8)/16/
      data ktab(9)/20/
      data ktab(10)/25/
      data ktab(11)/30/
      data ktab(12)/36/
      data ktab(13)/42/
      data ktab(14)/49/
      data ktab(15)/56/
      data ktab(16)/64/
      data ktab(20)/72/
      data ktab(24)/82/
      data ktab(28)/82/
      data ktab(32)/94/
      data ktab(36)/94/
      data ktab(40)/110/
      data ktab(44)/110/
      data ktab(48)/130/
      data ktab(52)/130/
      data ktab(56)/130/
      data ktab(60)/130/
      data ktab(64)/154/
      data ktab(68)/154/
      data ktab(72)/154/
      data ktab(76)/154/
      data ktab(80)/186/
      data ktab(84)/186/
      data ktab(88)/186/
      data ktab(92)/186/
      data ktab(96)/226/
c
c-----table of abscissae (x) and weights (a) for interval (-1,+1).
c
c**** n=2
      data x(1)/0.577350269189626  d0/, a(1)/1.000000000000000  d0/
c**** n=3
      data x(2)/0.774596669241483  d0/, a(2)/0.555555555555556  d0/
      data x(3)/0.000000000000000  d0/, a(3)/0.888888888888889  d0/
c**** n=4
      data x(4)/0.861136311594053  d0/, a(4)/0.347854845137454  d0/
      data x(5)/0.339981043584856  d0/, a(5)/0.652145154862546  d0/
c**** n=5
      data x(6)/0.906179845938664  d0/, a(6)/0.236926885056189  d0/
      data x(7)/0.538469310105683  d0/, a(7)/0.478628670499366  d0/
      data x(8)/0.000000000000000  d0/, a(8)/0.568888888888889  d0/
c**** n=6
      data x(9)/0.932469514203152  d0/, a(9)/0.171324492379170  d0/
      data x(10)/0.661209386466265 d0/, a(10)/0.360761573048139 d0/
      data x(11)/0.238619186083197 d0/, a(11)/0.467913934572691 d0/
c**** n=7
      data x(12)/0.949107912342759 d0/, a(12)/0.129484966168870 d0/
      data x(13)/0.741531185599394 d0/, a(13)/0.279705391489277 d0/
      data x(14)/0.405845151377397 d0/, a(14)/0.381830050505119 d0/
      data x(15)/0.000000000000000 d0/, a(15)/0.417959183673469 d0/
c**** n=8
      data x(16)/0.960289856497536 d0/, a(16)/0.101228536290376 d0/
      data x(17)/0.796666477413627 d0/, a(17)/0.222381034453374 d0/
      data x(18)/0.525532409916329 d0/, a(18)/0.313706645877887 d0/
      data x(19)/0.183434642495650 d0/, a(19)/0.362683783378362 d0/
c**** n=9
      data x(20)/0.968160239507626 d0/, a(20)/0.081274388361574 d0/
      data x(21)/0.836031107326636 d0/, a(21)/0.180648160694857 d0/
      data x(22)/0.613371432700590 d0/, a(22)/0.260610696402935 d0/
      data x(23)/0.324253423403809 d0/, a(23)/0.312347077040003 d0/
      data x(24)/0.000000000000000 d0/, a(24)/0.330239355001260 d0/
c**** n=10
      data x(25)/0.973906528517172 d0/, a(25)/0.066671344308688 d0/
      data x(26)/0.865063366688985 d0/, a(26)/0.149451349150581 d0/
      data x(27)/0.679409568299024 d0/, a(27)/0.219086362515982 d0/
      data x(28)/0.433395394129247 d0/, a(28)/0.269266719309996 d0/
      data x(29)/0.148874338981631 d0/, a(29)/0.295524224714753 d0/
c**** n=11
      data x(30)/0.978228658146057 d0/, a(30)/0.055668567116174 d0/
      data x(31)/0.887062599768095 d0/, a(31)/0.125580369464905 d0/
      data x(32)/0.730152005574049 d0/, a(32)/0.186290210927734 d0/
      data x(33)/0.519096129206812 d0/, a(33)/0.233193764591990 d0/
      data x(34)/0.269543155952345 d0/, a(34)/0.262804544510247 d0/
      data x(35)/0.000000000000000 d0/, a(35)/0.272925086777901 d0/
c**** n=12
      data x(36)/0.981560634246719 d0/, a(36)/0.047175336386512 d0/
      data x(37)/0.904117256370475 d0/, a(37)/0.106939325995318 d0/
      data x(38)/0.769902674194305 d0/, a(38)/0.160078328543346 d0/
      data x(39)/0.587317954286617 d0/, a(39)/0.203167426723066 d0/
      data x(40)/0.367831498998180 d0/, a(40)/0.233492536538355 d0/
      data x(41)/0.125233408511469 d0/, a(41)/0.249147045813403 d0/
c**** n=13
      data x(42)/0.984183054718588 d0/, a(42)/0.040484004765316 d0/
      data x(43)/0.917598399222978 d0/, a(43)/0.092121499837728 d0/
      data x(44)/0.801578090733310 d0/, a(44)/0.138873510219787 d0/
      data x(45)/0.642349339440340 d0/, a(45)/0.178145980761946 d0/
      data x(46)/0.448492751036447 d0/, a(46)/0.207816047536889 d0/
      data x(47)/0.230458315955135 d0/, a(47)/0.226283180262897 d0/
      data x(48)/0.000000000000000 d0/, a(48)/0.232551553230874 d0/
c**** n=14
      data x(49)/0.986283808696812 d0/, a(49)/0.035119460331752 d0/
      data x(50)/0.928434883663574 d0/, a(50)/0.080158087159760 d0/
      data x(51)/0.827201315069765 d0/, a(51)/0.121518570687903 d0/
      data x(52)/0.687292904811685 d0/, a(52)/0.157203167158194 d0/
      data x(53)/0.515248636358154 d0/, a(53)/0.185538397477938 d0/
      data x(54)/0.319112368927890 d0/, a(54)/0.205198463721296 d0/
      data x(55)/0.108054948707344 d0/, a(55)/0.215263853463158 d0/
c**** n=15
      data x(56)/0.987992518020485 d0/, a(56)/0.030753241996117 d0/
      data x(57)/0.937273392400706 d0/, a(57)/0.070366047488108 d0/
      data x(58)/0.848206583410427 d0/, a(58)/0.107159220467172 d0/
      data x(59)/0.724417731360170 d0/, a(59)/0.139570677926154 d0/
      data x(60)/0.570972172608539 d0/, a(60)/0.166269205816994 d0/
      data x(61)/0.394151347077563 d0/, a(61)/0.186161000015562 d0/
      data x(62)/0.201194093997435 d0/, a(62)/0.198431485327111 d0/
      data x(63)/0.000000000000000 d0/, a(63)/0.202578241925561 d0/
c**** n=16
      data x(64)/0.989400934991650 d0/, a(64)/0.027152459411754 d0/
      data x(65)/0.944575023073233 d0/, a(65)/0.062253523938648 d0/
      data x(66)/0.865631202387832 d0/, a(66)/0.095158511682493 d0/
      data x(67)/0.755404408355003 d0/, a(67)/0.124628971255534 d0/
      data x(68)/0.617876244402644 d0/, a(68)/0.149595988816577 d0/
      data x(69)/0.458016777657227 d0/, a(69)/0.169156519395003 d0/
      data x(70)/0.281603550779259 d0/, a(70)/0.182603415044924 d0/
      data x(71)/0.095012509837637 d0/, a(71)/0.189450610455069 d0/
c**** n=20
      data x(72)/0.993128599185094 d0/, a(72)/0.017614007139152 d0/
      data x(73)/0.963971927277913 d0/, a(73)/0.040601429800386 d0/
      data x(74)/0.912234428251325 d0/, a(74)/0.062672048334109 d0/
      data x(75)/0.839116971822218 d0/, a(75)/0.083276741576704 d0/
      data x(76)/0.746331906460150 d0/, a(76)/0.101930119817240 d0/
      data x(77)/0.636053680726515 d0/, a(77)/0.118194531961518 d0/
      data x(78)/0.510867001950827 d0/, a(78)/0.131688638449176 d0/
      data x(79)/0.373706088715419 d0/, a(79)/0.142096109318382 d0/
      data x(80)/0.227785851141645 d0/, a(80)/0.149172986472603 d0/
      data x(81)/0.076526521133497 d0/, a(81)/0.152753387130725 d0/
c**** n=24
      data x(82)/0.995187219997021 d0/, a(82)/0.012341229799987 d0/
      data x(83)/0.974728555971309 d0/, a(83)/0.028531388628933 d0/
      data x(84)/0.938274552002732 d0/, a(84)/0.044277438817419 d0/
      data x(85)/0.886415527004401 d0/, a(85)/0.059298584915436 d0/
      data x(86)/0.820001985973902 d0/, a(86)/0.073346481411080 d0/
      data x(87)/0.740124191578554 d0/, a(87)/0.086190161531953 d0/
      data x(88)/0.648093651936975 d0/, a(88)/0.097618652104113 d0/
      data x(89)/0.545421471388839 d0/, a(89)/0.107444270115965 d0/
      data x(90)/0.433793507626045 d0/, a(90)/0.115505668053725 d0/
      data x(91)/0.315042679696163 d0/, a(91)/0.121670472927803 d0/
      data x(92)/0.191118867473616 d0/, a(92)/0.125837456346828 d0/
      data x(93)/0.064056892862605 d0/, a(93)/0.127938195346752 d0/
c**** n=32
      data x(94)/0.997263861849481 d0/, a(94)/0.007018610009470 d0/
      data x(95)/0.985611511545268 d0/, a(95)/0.016274394730905 d0/
      data x(96)/0.964762255587506 d0/, a(96)/0.025392065309262 d0/
      data x(97)/0.934906075937739 d0/, a(97)/0.034273862913021 d0/
      data x(98)/0.896321155766052 d0/, a(98)/0.042835898022226 d0/
      data x(99)/0.849367613732569 d0/, a(99)/0.050998059262376 d0/
      data x(100)/0.794483795967942d0/, a(100)/0.058684093478535d0/
      data x(101)/0.732182118740289d0/, a(101)/0.065822222776361d0/
      data x(102)/0.663044266930215d0/, a(102)/0.072345794108848d0/
      data x(103)/0.587715757240762d0/, a(103)/0.078193895787070d0/
      data x(104)/0.506899908932229d0/, a(104)/0.083311924226946d0/
      data x(105)/0.421351276130635d0/, a(105)/0.087652093004403d0/
      data x(106)/0.331868602282127d0/, a(106)/0.091173878695763d0/
      data x(107)/0.239287362252137d0/, a(107)/0.093844399080804d0/
      data x(108)/0.144471961582796d0/, a(108)/0.095638720079274d0/
      data x(109)/0.048307665687738d0/, a(109)/0.096540088514727d0/
c**** n=40
      data x(110)/0.998237709710559d0/, a(110)/0.004521277098533d0/
      data x(111)/0.990726238699457d0/, a(111)/0.010498284531152d0/
      data x(112)/0.977259949983774d0/, a(112)/0.016421058381907d0/
      data x(113)/0.957916819213791d0/, a(113)/0.022245849194166d0/
      data x(114)/0.932812808278676d0/, a(114)/0.027937006980023d0/
      data x(115)/0.902098806968874d0/, a(115)/0.033460195282547d0/
      data x(116)/0.865959503212259d0/, a(116)/0.038782167974472d0/
      data x(117)/0.824612230833311d0/, a(117)/0.043870908185673d0/
      data x(118)/0.778305651426519d0/, a(118)/0.048695807635072d0/
      data x(119)/0.727318255189927d0/, a(119)/0.053227846983936d0/
      data x(120)/0.671956684614179d0/, a(120)/0.057439769099391d0/
      data x(121)/0.612553889667980d0/, a(121)/0.061306242492928d0/
      data x(122)/0.549467125095128d0/, a(122)/0.064804013456601d0/
      data x(123)/0.483075801686178d0/, a(123)/0.067912045815233d0/
      data x(124)/0.413779204371605d0/, a(124)/0.070611647391286d0/
      data x(125)/0.341994090825758d0/, a(125)/0.072886582395804d0/
      data x(126)/0.268152185007253d0/, a(126)/0.074723169057968d0/
      data x(127)/0.192697580701371d0/, a(127)/0.076110361900626d0/
      data x(128)/0.116084070675255d0/, a(128)/0.077039818164247d0/
      data x(129)/0.038772417506050d0/, a(129)/0.077505947978424d0/
c**** n=48
      data x(130)/0.998771007252426d0/, a(130)/0.003153346052305d0/
      data x(131)/0.993530172266350d0/, a(131)/0.007327553901276d0/
      data x(132)/0.984124583722826d0/, a(132)/0.011477234579234d0/
      data x(133)/0.970591592546247d0/, a(133)/0.015579315722943d0/
      data x(134)/0.952987703160430d0/, a(134)/0.019616160457355d0/
      data x(135)/0.931386690706554d0/, a(135)/0.023570760839324d0/
      data x(136)/0.905879136715569d0/, a(136)/0.027426509708356d0/
      data x(137)/0.876572020274247d0/, a(137)/0.031167227832798d0/
      data x(138)/0.843588261624393d0/, a(138)/0.034777222564770d0/
      data x(139)/0.807066204029442d0/, a(139)/0.038241351065830d0/
      data x(140)/0.767159032515740d0/, a(140)/0.041545082943464d0/
      data x(141)/0.724034130923814d0/, a(141)/0.044674560856694d0/
      data x(142)/0.677872379632663d0/, a(142)/0.047616658492490d0/
      data x(143)/0.628867396776513d0/, a(143)/0.050359035553854d0/
      data x(144)/0.577224726083972d0/, a(144)/0.052890189485193d0/
      data x(145)/0.523160974722233d0/, a(145)/0.055199503699984d0/
      data x(146)/0.466902904750958d0/, a(146)/0.057277292100403d0/
      data x(147)/0.408686481990716d0/, a(147)/0.059114839698395d0/
      data x(148)/0.348755886292160d0/, a(148)/0.060704439165893d0/
      data x(149)/0.287362487355455d0/, a(149)/0.062039423159892d0/
      data x(150)/0.224763790394689d0/, a(150)/0.063114192286254d0/
      data x(151)/0.161222356068891d0/, a(151)/0.063924238584648d0/
      data x(152)/0.097004699209462d0/, a(152)/0.064466164435950d0/
      data x(153)/0.032380170962869d0/, a(153)/0.064737696812683d0/
c**** n=64
      data x(154)/0.999305041735772d0/, a(154)/0.001783280721696d0/
      data x(155)/0.996340116771955d0/, a(155)/0.004147033260562d0/
      data x(156)/0.991013371476744d0/, a(156)/0.006504457968978d0/
      data x(157)/0.983336253884625d0/, a(157)/0.008846759826363d0/
      data x(158)/0.973326827789910d0/, a(158)/0.011168139460131d0/
      data x(159)/0.961008799652053d0/, a(159)/0.013463047896718d0/
      data x(160)/0.946411374858402d0/, a(160)/0.015726030476024d0/
      data x(161)/0.929569172131939d0/, a(161)/0.017951715775697d0/
      data x(162)/0.910522137078502d0/, a(162)/0.020134823153530d0/
      data x(163)/0.889315445995114d0/, a(163)/0.022270173808383d0/
      data x(164)/0.865999398154092d0/, a(164)/0.024352702568710d0/
      data x(165)/0.840629296252580d0/, a(165)/0.026377469715054d0/
      data x(166)/0.813265315122797d0/, a(166)/0.028339672614259d0/
      data x(167)/0.783972358943341d0/, a(167)/0.030234657072402d0/
      data x(168)/0.752819907260531d0/, a(168)/0.032057928354851d0/
      data x(169)/0.719881850171610d0/, a(169)/0.033805161837141d0/
      data x(170)/0.685236313054233d0/, a(170)/0.035472213256882d0/
      data x(171)/0.648965471254657d0/, a(171)/0.037055128540240d0/
      data x(172)/0.611155355172393d0/, a(172)/0.038550153178615d0/
      data x(173)/0.571895646202634d0/, a(173)/0.039953741132720d0/
      data x(174)/0.531279464019894d0/, a(174)/0.041262563242623d0/
      data x(175)/0.489403145707052d0/, a(175)/0.042473515123653d0/
      data x(176)/0.446366017253464d0/, a(176)/0.043583724529323d0/
      data x(177)/0.402270157963991d0/, a(177)/0.044590558163756d0/
      data x(178)/0.357220158337668d0/, a(178)/0.045491627927418d0/
      data x(179)/0.311322871990210d0/, a(179)/0.046284796581314d0/
      data x(180)/0.264687162208767d0/, a(180)/0.046968182816210d0/
      data x(181)/0.217423643740007d0/, a(181)/0.047540165714830d0/
      data x(182)/0.169644420423992d0/, a(182)/0.047999388596458d0/
      data x(183)/0.121462819296120d0/, a(183)/0.048344762234802d0/
      data x(184)/0.072993121787799d0/, a(184)/0.048575467441503d0/
      data x(185)/0.024350292663424d0/, a(185)/0.048690957009139d0/
c**** n=80
      data x(186)/0.999553822651630d0/, a(186)/0.001144950003186d0/
      data x(187)/0.997649864398237d0/, a(187)/0.002663533589512d0/
      data x(188)/0.994227540965688d0/, a(188)/0.004180313124694d0/
      data x(189)/0.989291302499755d0/, a(189)/0.005690922451403d0/
      data x(190)/0.982848572738629d0/, a(190)/0.007192904768117d0/
      data x(191)/0.974909140585727d0/, a(191)/0.008683945269260d0/
      data x(192)/0.965485089043799d0/, a(192)/0.010161766041103d0/
      data x(193)/0.954590766343634d0/, a(193)/0.011624114120797d0/
      data x(194)/0.942242761309872d0/, a(194)/0.013068761592401d0/
      data x(195)/0.928459877172445d0/, a(195)/0.014493508040509d0/
      data x(196)/0.913263102571757d0/, a(196)/0.015896183583725d0/
      data x(197)/0.896675579438770d0/, a(197)/0.017274652056269d0/
      data x(198)/0.878722567678213d0/, a(198)/0.018626814208299d0/
      data x(199)/0.859431406663111d0/, a(199)/0.019950610878141d0/
      data x(200)/0.838831473580255d0/, a(200)/0.021244026115782d0/
      data x(201)/0.816954138681463d0/, a(201)/0.022505090246332d0/
      data x(202)/0.793832717504605d0/, a(202)/0.023731882865930d0/
      data x(203)/0.769502420135041d0/, a(203)/0.024922535764115d0/
      data x(204)/0.744000297583597d0/, a(204)/0.026075235767565d0/
      data x(205)/0.717365185362099d0/, a(205)/0.027188227500486d0/
      data x(206)/0.689637644342027d0/, a(206)/0.028259816057276d0/
      data x(207)/0.660859898986119d0/, a(207)/0.029288369583267d0/
      data x(208)/0.631075773046871d0/, a(208)/0.030272321759557d0/
      data x(209)/0.600330622829751d0/, a(209)/0.031210174188114d0/
      data x(210)/0.568671268122709d0/, a(210)/0.032100498673487d0/
      data x(211)/0.536145920897131d0/, a(211)/0.032941939397645d0/
      data x(212)/0.502804111888784d0/, a(212)/0.033733214984611d0/
      data x(213)/0.468696615170544d0/, a(213)/0.034473120451753d0/
      data x(214)/0.433875370831756d0/, a(214)/0.035160529044747d0/
      data x(215)/0.398393405881969d0/, a(215)/0.035794393953416d0/
      data x(216)/0.362304753499487d0/, a(216)/0.036373749905835d0/
      data x(217)/0.325664370747701d0/, a(217)/0.036897714638276d0/
      data x(218)/0.288528054884511d0/, a(218)/0.037365490238730d0/
      data x(219)/0.250952358392272d0/, a(219)/0.037776364362001d0/
      data x(220)/0.212994502857666d0/, a(220)/0.038129711314477d0/
      data x(221)/0.174712291832646d0/, a(221)/0.038424993006959d0/
      data x(222)/0.136164022809143d0/, a(222)/0.038661759774076d0/
      data x(223)/0.097408398441584d0/, a(223)/0.038839651059051d0/
      data x(224)/0.058504437152420d0/, a(224)/0.038958395962769d0/
      data x(225)/0.019511383256793d0/, a(225)/0.039017813656306d0/
c**** n=96
      data x(226)/0.999689503883230d0/, a(226)/0.000796792065552d0/
      data x(227)/0.998364375863181d0/, a(227)/0.001853960788946d0/
      data x(228)/0.995981842987209d0/, a(228)/0.002910731817934d0/
      data x(229)/0.992543900323762d0/, a(229)/0.003964554338444d0/
      data x(230)/0.988054126329623d0/, a(230)/0.005014202742927d0/
      data x(231)/0.982517263563014d0/, a(231)/0.006058545504235d0/
      data x(232)/0.975939174585136d0/, a(232)/0.007096470791153d0/
      data x(233)/0.968326828463264d0/, a(233)/0.008126876925698d0/
      data x(234)/0.959688291448742d0/, a(234)/0.009148671230783d0/
      data x(235)/0.950032717784437d0/, a(235)/0.010160770535008d0/
      data x(236)/0.939370339752755d0/, a(236)/0.011162102099838d0/
      data x(237)/0.927712456722308d0/, a(237)/0.012151604671088d0/
      data x(238)/0.915071423120898d0/, a(238)/0.013128229566961d0/
      data x(239)/0.901460635315852d0/, a(239)/0.014090941772314d0/
      data x(240)/0.886894517402420d0/, a(240)/0.015038721026994d0/
      data x(241)/0.871388505909296d0/, a(241)/0.015970562902562d0/
      data x(242)/0.854959033434601d0/, a(242)/0.016885479864245d0/
      data x(243)/0.837623511228187d0/, a(243)/0.017782502316045d0/
      data x(244)/0.819400310737931d0/, a(244)/0.018660679627411d0/
      data x(245)/0.800308744139140d0/, a(245)/0.019519081140145d0/
      data x(246)/0.780369043867433d0/, a(246)/0.020356797154333d0/
      data x(247)/0.759602341176647d0/, a(247)/0.021172939892191d0/
      data x(248)/0.738030643744400d0/, a(248)/0.021966644438744d0/
      data x(249)/0.715676812348967d0/, a(249)/0.022737069658329d0/
      data x(250)/0.692564536642171d0/, a(250)/0.023483399085926d0/
      data x(251)/0.668718310043916d0/, a(251)/0.024204841792364d0/
      data x(252)/0.644163403784967d0/, a(252)/0.024900633222483d0/
      data x(253)/0.618925840125468d0/, a(253)/0.025570036005349d0/
      data x(254)/0.593032364777572d0/, a(254)/0.026212340735672d0/
      data x(255)/0.566510418561397d0/, a(255)/0.026826866725591d0/
      data x(256)/0.539388108324357d0/, a(256)/0.027412962726029d0/
      data x(257)/0.511694177154667d0/, a(257)/0.027970007616848d0/
      data x(258)/0.483457973920596d0/, a(258)/0.028497411065085d0/
      data x(259)/0.454709422167743d0/, a(259)/0.028994614150555d0/
      data x(260)/0.425478988407300d0/, a(260)/0.029461089958167d0/
      data x(261)/0.395797649828908d0/, a(261)/0.029896344136328d0/
      data x(262)/0.365696861472313d0/, a(262)/0.030299915420827d0/
      data x(263)/0.335208522892625d0/, a(263)/0.030671376123669d0/
      data x(264)/0.304364944354496d0/, a(264)/0.031010332586313d0/
      data x(265)/0.273198812591049d0/, a(265)/0.031316425596861d0/
      data x(266)/0.241743156163840d0/, a(266)/0.031589330770727d0/
      data x(267)/0.210031310460567d0/, a(267)/0.031828758894411d0/
      data x(268)/0.178096882367618d0/, a(268)/0.032034456231992d0/
      data x(269)/0.145973714654896d0/, a(269)/0.032206204794030d0/
      data x(270)/0.113695850110665d0/, a(270)/0.032343822568575d0/
      data x(271)/0.081297495464425d0/, a(271)/0.032447163714064d0/
      data x(272)/0.048812985136049d0/, a(272)/0.032516118713868d0/
      data x(273)/0.016276744849602d0/, a(273)/0.032550614492363d0/
c
c
c-----test n
      alpha=0.5d0*(ax+bx)
      beta=0.5d0*(bx-ax)
      if( n.lt.1 .or. n.gt.96 ) go to 100
      if(n.ne.1) go to 1
      z(1)=alpha
      w(1)=bx-ax
      return
c
    1 if (n.le.16) go to 3
      if (n.gt.24) go to 4
      n=4*(n/4)
      go to 3
    4 if (n.gt.48) go to 5
      n=8*(n/8)
      go to 3
    5 n=16*(n/16)
c
c----- set k equal to initial subscript and store results
    3 k=ktab(n)
      m=n/2
      do 2 j=1,m
      jtab=k-1+j
      wtemp=beta*a(jtab)
      delta=beta*x(jtab)
      z(j)=alpha-delta
      w(j)=wtemp
      jp=n+1-j
      z(jp)=alpha+delta
      w(jp)=wtemp
    2 continue
      if((n-m-m).eq.0) return
      z(m+1)=alpha
      jmid=k+m
      w(m+1)=beta*a(jmid)
      return
c
  100 zn=n
      write(6,200) zn
  200 format(/////' error in gset. n has the non-permissible value',
     1e11.3/' execution terminated.')
      stop
      end
