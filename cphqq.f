c     program cphqq
c
      implicit real*8 (a-h,o-z)
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c**** dimension vv(6*n),s(n),u(n),a(2*(n+1)*2*(n+1)),b(4*(n+1)),
c**** aa(6*n/2*(n+1)),qq(n),eq(n) ****
c**** real*8 vv(576),s(96),u(96),a(37636),b(388),aa(27936),
c****1 qq(96),eq(96)
      real*8 vv(6000),s(1000),u(1000),a(4008004),b(4004),
     1 aa(3003000),qq(1000),eq(1000)
c
      external chrqq
c
      open (unit=5,file='dphqq.d')
      open (unit=6,file='phqq.d')
      open (unit=7,file='pchqq.d')
      open (unit=8,file='sm99_np.1000')
      open (unit=11,file='stoks1.d')
      open (unit=12,file='stoks2.d')
      open (unit=13,file='phpperr.d')
c
      kread=5
      kwrite=6
      kpunch=7
      kda(8)=8
      kda(1)=11
      kda(2)=12
      kda(3)=13
c
      call phases (chrqq,vv,s,u,a,b,aa,qq,eq)
c
      end
