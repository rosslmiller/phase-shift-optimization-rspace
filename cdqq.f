c        call deuter
c     program cdqq
      implicit real*8 (a-h,o-z)
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
      dimension s(96),us(96),usf(96),a(36864),b(192),aa(18624)
      dimension z(192),zz(96),l(192),m(192)
      external chrqq
      open (unit=5,file='ddqq.d')
      open (unit=6,file='deuqq.d')
      open (unit=7,file='dwqqq.d')
      kread=5
      kwrite=6
      kpunch=7
      call deuter (chrqq,s,us,usf,a,b,aa,z,zz,l,m)
      end
