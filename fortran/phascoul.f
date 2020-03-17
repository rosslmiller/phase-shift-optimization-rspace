c
c*********** this version of phascoul uses the 
c*********** subroutine gauss for the gauss points and weights.
c*********** 7/31/07
c
c*********** watch it: common block /cpts/ is changed!!!
c
c
         subroutine phases (pot,vv,s,u,a,b,aa,qq,eq)                    
c                                                                       
c                                                                       
c        the file of this code is called  phascoul.f                    
c        and contains besides phases all subroutines needed             
c                                                                       
c                                                                       
c        routine computes phase-shifts of nucleon-nucleon scattering    
c        for a given nn-potential pot                                   
c        this version of the code includes the coulomb distorsion option
c                                                                       
c                                                                       
c        author: r. machleidt                                           
c                institut fuer theoretische kernphysik bonn             
c                nussallee 14-16                                        
c                d-5300 bonn, w. germany                                
c                                                                       
c                                                                       
      implicit real*8 (a-h,o-z)                                         
      external pot                                                      
      common /crdwrt/ kread,kwrite,kpunch,kda(9)                        
c                                                                       
c        common block for coulomb subroutine                            
      common /ccoul/ wn,rcoul,ncoul,iprop,indcou,indrel,indnnb          
c                                                                       
c        common block for bonn full model interface subroutine          
      common /cnnex/ iiienn                                             
c                                                                       
c        arguments of the potential subroutine pot being called in this 
c        program                                                        
c                                                                       
      common /cpot/   v(6),xmev,ymev                                    
      common /cstate/ j,heform,sssing,tttrip,cccoup,endep,label         
      common /cpoted/ q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,    
     1                znrl,zrel,smev,noced                              
c**** common /cpts/   q(97),c,n1,ix,iy                                  
      common /cpts/   q(1001),c,n1,ix,iy                                  
c                                                                       
c        specifications for these common blocks                         
      logical indcou,indrel,indnnb                                      
      logical heform,sssing,tttrip,cccoup,endep                         
      logical noced                                                     
c                                                                       
c        further specifications                                         
c                                                                       
      real*8 vv(1),s(1),u(1),a(1),b(1),aa(1),qq(1),eq(1)                
c**** dimension aan(27936)                                              
c**** dimension aac(27936)                                              
      dimension aan(3003000)                                            
      dimension aac(3003000)                                            
      real*8 delta(5),ddelta(5)                                         
      data delta/5*0.d0/                                                
      real*8 rb(2)                                                      
      real*8 r(6)                                                       
      real*8 rb1(3),rb2(3)                                              
      data rb1/3*0.d0/,rb2/3*0.d0/                                      
      data uf/197.32697d0/                                                
      data pih/1.570796326794897d0/                                     
      real*4 ops                                                        
      data ops/1.e-15/                                                  
      real*4 chitot                                                     
      real*4 expch(3)                                                   
c                                                                       
c        the following dimension for at most 80 elabs                   
      data memax/81/                                                    
      character*4 llps(80,5)                                            
      real*8 elab(81),q0(80),q0q(80),eq0(80),arel(80)                   
      real*4 deld(4,80,5),delr(4,80,5),chi(80,5)                        
c                                                                       
c        the following dimension for at most j up to 19                 
      integer nj(20)                                                    
c                                                                       
c                                                                       
c        the following dimensions are for coulomb corrections           
c                                                                       
      dimension dss(2,2),fc0(2,2),gc0(2,2),fb0(2,2),gb0(2,2),hm3(2,2)   
      dimension fc(2,2),gc(2,2),fd(2,2),gd(2,2),hm1(2,2),hm2(2,2)       
      dimension ale(3,3),ble(3,3),cle(3,3)                              
c                                                                       
c                                                                       
      character*4 name(3),nname(15)                                     
      data nalt/-1/                                                     
c                                                                       
      character*4 ichar                                                 
      character*1 state(2),state3                                       
      character*1 multi(4)                                              
      data multi/'1',3*'3'/                                             
      integer ldel(4)                                                   
      data ldel/0,0,-1,1/                                               
      character*1 spd(50)                                               
      data spd/'s','p','d','f','g','h','i','k','l','m','n',             
     1'o','q','r','t','u','v','w','x','y','z','a','b','c','e',25*' '/   
      character*4 charmr                                                
      character*4 charm                                                 
      character*4 charlb                                                
      data charmr/'merr'/                                               
      data charm/'m   '/                                                
      data charlb/'labl'/                                               
      character*1 chars,chart                                           
      data chars/'s'/,chart/'t'/                                        
      character*1 lab1,lab2                                             
      character*4 blanks                                                
      character*4 nnnp                                                  
      character*4 nnpp                                                  
      data blanks/'    '/                                               
      data nnnp/'n-p '/                                                 
      data nnpp/'p-p '/                                                 
      character*1 blank                                                 
      character*1 eps                                                   
      character*1 eeps                                                  
      data blank/' '/                                                   
      data eps/'e'/                                                     
      data eeps/'e'/                                                    
      logical sing,trip,coup                                            
      logical ssing,ttrip,ccoup                                         
      logical indbrn                                                    
      logical indrma                                                    
      logical indrsw                                                    
      logical indqua                                                    
      logical indwrt                                                    
      logical indpts                                                    
      logical indshp                                                    
      logical indpch                                                    
      logical inddeg                                                    
      logical inderr                                                    
      logical inderg                                                    
      logical indp1                                                     
      logical indpl                                                     
      data indbrn/.false./                                              
      data indrma/.false./                                              
      data indrsw/.false./                                              
      data indqua/.false./                                              
      data indwrt/.false./                                              
      data indpts/.false./                                              
      data indshp/.false./                                              
      data indpch/.false./                                              
      data inddeg/.false./                                              
      data inderr/.false./                                              
      data inderg/.false./                                              
      data indp1/.false./                                               
      data indpl/.false./                                               
      logical indj                                                      
      logical indelb                                                    
c                                                                       
c                                                                       
c        datum and parameter for coulomb                                
c                                                                       
      data alfinv/137.0359895d0/                                        
      data ifail/0/                                                     
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c        arrays used for applying the error matrix                      
c        -----------------------------------------                      
c                                                                       
      dimension elaber(10),npherr(10),errma(24,24,10)                   
      dimension pherr(24,10),diff(24)                                   
      character*3 namer1(24)                                            
      character*5 namer2(24,10)                                         
      dimension delj(5,5,80)                                            
      data delj/2000*0.d0/                                              
      dimension phdel(24,10)                                            
      data phdel/240*0.d0/                                              
c                                                                       
c        for pp-only error matrix stoks1.d                              
      dimension iverr1(24),jerr1(24)                                    
      data iverr1/0,1,3,2,3,1,5,4,2,3,1,5,4,11*0/                       
      data jerr1 /0,0,0,1,2,2,2,2,3,4,4,4,4,11*0/                       
c                                                                       
c        for pp+np error matrix stoks2.d                                
      dimension iverr2(24,10),jerr2(24,10)                              
      data iverr2                                                       
     1           /0,2,3,1,3,5,4,2,1,3,5,4,2,0,0,0,0,0,0,0,0,0,0,0,      
     2            0,2,3,1,3,5,4,2,1,3,5,4,2,0,0,0,0,0,0,0,0,0,0,0,      
     3            0,2,3,2,3,3,5,4,2,1,3,5,4,2,0,0,0,0,0,0,0,0,0,0,      
     4            0,2,3,2,3,1,1,3,5,4,2,1,3,5,4,2,0,0,0,0,0,0,0,0,      
     5            0,2,3,2,3,5,1,1,3,5,4,2,1,3,5,4,2,0,0,0,0,0,0,0,      
     6            0,2,3,2,1,3,5,1,1,3,5,4,2,1,3,5,4,2,0,0,0,0,0,0,      
     7            0,2,3,2,1,3,5,1,1,3,5,4,2,1,3,5,4,2,0,0,0,0,0,0,      
     8            0,2,3,2,1,3,5,4,2,1,1,3,5,4,2,1,3,5,4,2,0,0,0,0,      
     9            0,2,3,2,1,3,5,4,2,1,3,5,4,1,1,3,5,4,2,1,3,5,4,2,      
     *            0,2,3,2,1,3,5,4,2,1,3,5,4,1,1,3,5,4,2,1,3,5,4,2/      
      data jerr2                                                        
     1           /0,0,0,0,1,1,1,1,2,3,3,3,3,4,0,0,0,0,0,0,0,0,0,0,      
     2            0,0,0,0,1,1,1,1,2,3,3,3,3,4,0,0,0,0,0,0,0,0,0,0,      
     3            0,0,0,1,2,1,1,1,1,2,3,3,3,3,4,0,0,0,0,0,0,0,0,0,      
     4            0,0,0,1,2,0,1,1,1,1,2,3,3,3,3,4,0,0,0,0,0,0,0,0,      
     5            0,0,0,1,2,2,0,1,1,1,1,2,3,3,3,3,4,0,0,0,0,0,0,0,      
     6            0,0,0,1,2,2,2,0,1,1,1,1,2,3,3,3,3,4,0,0,0,0,0,0,      
     7            0,0,0,1,2,2,2,0,1,1,1,1,2,3,3,3,3,4,0,0,0,0,0,0,      
     8            0,0,0,1,2,2,2,2,3,0,1,1,1,1,2,3,3,3,3,4,0,0,0,0,      
     9            0,0,0,1,2,2,2,2,3,4,4,4,4,0,1,1,1,1,2,3,3,3,3,4,      
     *            0,0,0,1,2,2,2,2,3,4,4,4,4,0,1,1,1,1,2,3,3,3,3,4/      
c                                                                       
      logical indema                                                    
      data indema/.false./                                              
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c        data from phase-shift analyses                                 
c        ------------------------------                                 
c                                                                       
c                                                                       
c        number of phase-shift analyses                                 
      data nps/7/                                                       
      integer mps(7)                                                    
      logical indps(7)                                                  
      data indps/7*.false./                                             
c        names of phase-shift analyses                                  
      character*4 labps(7)                                              
      data labps /'    ','vs35','bg90','ni90','n93p','n93n','ni95'/     
c                                                                       
c      '    ': this analysis is read in                                 
c                                                                       
c        vs35: arndt analysis vs35 for n-p, said program 5/27/91        
c                                                                       
c        bg90: n-p, d.v. bugg, phys. rev. c 41(1990) 2708               
c                                                                       
c        ni90: nijmegen p-p single-energy analysis,                     
c              bergervoet et al., phys.rev. c 41 (1990), 1435           
c                                                                       
c        n93p: nijmegen p-p multi-energy analysis,                      
c              stoks et al., phys. rev. c48, 792 (1993);                
c                    1s0, 1. - 50. mev, corrected acc. to               
c                    stoks et al., phys. rev. c49, 2950 (1994),         
c                    and stoks, priv comm. 11/6/94.                     
c                                                                       
c        n93n: nijmegen n-p multi-energy analysis,                      
c               stoks et al., phys.rev. c 48, 792 (1993)                
c                                                                       
c        ni95: nijmegen n-p single-energy analysis,                     
c              stoks privat communication of 4/4/95                     
c                                                                       
c                                                                       
c        number of data                                                 
      integer mpsa(5,10,7)                                              
c                                                                       
c        number of '    ' data                                          
      data   mpsa/50*0,                                                 
c                                                                       
c        number of vs35 data                                            
c        j = 0                                                          
     1            7,0,7,0,0,                                            
c        j = 1                                                          
     1            7,7,7,5,6,                                            
c        j = 2                                                          
     1            4,5,7,4,6,                                            
c        j = 3                                                          
     1            1,2,4,3,4,                                            
c        j = 4                   j ge 5                                 
     1            0,3,1,0,1,     25*0,                                  
c                                                                       
c        number of bg90 data                                            
c        j = 0                                                          
     1            3,0,3,0,0,                                            
c        j = 1                                                          
     1            3,3,3,3,3,                                            
c        j = 2                                                          
     1            3,3,3,3,3,                                            
c        j = 3                                                          
     1            3,3,3,3,3,                                            
c        j = 4                                                          
     1            3,3,3,3,3,                                            
c        j = 5                                                          
     1            0,3,3,0,0,                                            
c        j = 6                   j ge 7                                 
     1            0,0,3,0,0,     15*0,                                  
c                                                                       
c        number of ni90 data                                            
c        j = 0                                                          
     1            10,0,8,0,0,                                           
c        j = 1                                                          
     1            0,8,0,0,0,                                            
c        j = 2                                                          
     1            8,0,8,4,6,                                            
c        j = 3                                                          
     1            0,3,0,0,0,                                            
c        j = 4                   j ge 5                                 
     1            2,0,2,2,2,     25*0,                                  
c                                                                       
c        number of n93p data                                            
c        j = 0                                                          
     1            11,0,11,0,0,                                          
c        j = 1                                                          
     1            0,11,0,0,0,                                           
c        j = 2                                                          
     1            11,0,11,11,11,                                        
c        j = 3                                                          
     1            0,11,0,0,0,                                           
c        j = 4                   j ge 5                                 
     1            11,0,11,11,11,     25*0,                              
c                                                                       
c        number of n93n data                                            
c        j = 0                                                          
     1            11,0,11,0,0,                                          
c        j = 1                                                          
     1            11,11,11,11,11,                                       
c        j = 2                                                          
     1            11,11,11,11,11,                                       
c        j = 3                                                          
     1            11,11,11,11,11,                                       
c        j = 4                   j ge 5                                 
     1            11,11,11,11,11,     25*0,                             
c                                                                       
c        number of ni95 data                                            
c        j = 0                                                          
     1            10,0,0,0,0,                                           
c        j = 1                                                          
     1            8,0,9,6,6,                                            
c        j = 2                                                          
     1            0,5,0,0,0,                                            
c        j = 3                                                          
     1            3,0,3,2,3,                                            
c        j = 4                   j ge 5                                 
     1            0,2,0,0,0,          25*0/                             
c                                                                       
c                                                                       
c                                                                       
      real*4 psa(3,50,5,10,7)                                           
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c        phase-shift analysis vl35, ips=1                               
c        --------------------------------                               
c        arndt  analysis vl35 for p-p .   4/09/91                       
c                                                                       
c        note: this analysis is overwritten by the one                  
c              that is read in                                          
c                                                                       
c                                                                       
c        j = 0                                                          
c        1s0                                                            
      data psa(1,1,1,1,1)/ 0.15000e+02/                                 
      data psa(2,1,1,1,1)/ 0.52810e+02/, psa(3,1,1,1,1)/ 0.55000e+00/   
      data psa(1,2,1,1,1)/ 0.25000e+02/                                 
      data psa(2,2,1,1,1)/ 0.48420e+02/, psa(3,2,1,1,1)/ 0.23000e+00/   
      data psa(1,3,1,1,1)/ 0.50000e+02/                                 
      data psa(2,3,1,1,1)/ 0.39040e+02/, psa(3,3,1,1,1)/ 0.09000e+00/   
      data psa(1,4,1,1,1)/ 0.10000e+03/                                 
      data psa(2,4,1,1,1)/ 0.24620e+02/, psa(3,4,1,1,1)/ 0.55000e+00/   
      data psa(1,5,1,1,1)/ 0.15000e+03/                                 
      data psa(2,5,1,1,1)/ 0.14480e+02/, psa(3,5,1,1,1)/ 0.34000e+00/   
      data psa(1,6,1,1,1)/ 0.20000e+03/                                 
      data psa(2,6,1,1,1)/ 0.74100e+01/, psa(3,6,1,1,1)/ 0.38000e+00/   
      data psa(1,7,1,1,1)/ 0.30000e+03/                                 
      data psa(2,7,1,1,1)/-0.71100e+01/, psa(3,7,1,1,1)/ 0.32000e+00/   
c        3p0                                                            
      data psa(1,1,3,1,1)/ 0.15000e+02/                                 
      data psa(2,1,3,1,1)/ 0.54600e+01/, psa(3,1,3,1,1)/ 0.14800e+01/   
      data psa(1,2,3,1,1)/ 0.25000e+02/                                 
      data psa(2,2,3,1,1)/ 0.86700e+01/, psa(3,2,3,1,1)/ 0.42000e+00/   
      data psa(1,3,3,1,1)/ 0.50000e+02/                                 
      data psa(2,3,3,1,1)/ 0.11600e+02/, psa(3,3,3,1,1)/ 0.12000e+00/   
      data psa(1,4,3,1,1)/ 0.10000e+03/                                 
      data psa(2,4,3,1,1)/ 0.11660e+02/, psa(3,4,3,1,1)/ 0.16600e+01/   
      data psa(1,5,3,1,1)/ 0.15000e+03/                                 
      data psa(2,5,3,1,1)/ 0.54100e+01/, psa(3,5,3,1,1)/ 0.28000e+00/   
      data psa(1,6,3,1,1)/ 0.20000e+03/                                 
      data psa(2,6,3,1,1)/ 0.38000e+00/, psa(3,6,3,1,1)/ 0.33000e+00/   
      data psa(1,7,3,1,1)/ 0.30000e+03/                                 
      data psa(2,7,3,1,1)/-0.97600e+01/, psa(3,7,3,1,1)/ 0.40000e+00/   
c        j = 1                                                          
c        3p1                                                            
      data psa(1,1,2,2,1)/ 0.15000e+02/                                 
      data psa(2,1,2,2,1)/-0.31200e+01/, psa(3,1,2,2,1)/ 0.77000e+00/   
      data psa(1,2,2,2,1)/ 0.25000e+02/                                 
      data psa(2,2,2,2,1)/-0.48500e+01/, psa(3,2,2,2,1)/ 0.20000e+00/   
      data psa(1,3,2,2,1)/ 0.50000e+02/                                 
      data psa(2,3,2,2,1)/-0.82300e+01/, psa(3,3,2,2,1)/ 0.05000e+00/   
      data psa(1,4,2,2,1)/ 0.10000e+03/                                 
      data psa(2,4,2,2,1)/-0.13540e+02/, psa(3,4,2,2,1)/ 0.33000e+00/   
      data psa(1,5,2,2,1)/ 0.15000e+03/                                 
      data psa(2,5,2,2,1)/-0.17580e+02/, psa(3,5,2,2,1)/ 0.08000e+00/   
      data psa(1,6,2,2,1)/ 0.20000e+03/                                 
      data psa(2,6,2,2,1)/-0.21290e+02/, psa(3,6,2,2,1)/ 0.17000e+00/   
      data psa(1,7,2,2,1)/ 0.30000e+03/                                 
      data psa(2,7,2,2,1)/-0.28280e+02/, psa(3,7,2,2,1)/ 0.24000e+00/   
c        j = 2                                                          
c        1d2                                                            
      data psa(1,1,1,3,1)/ 0.10000e+03/                                 
      data psa(2,1,1,3,1)/ 0.37000e+01/, psa(3,1,1,3,1)/ 0.15000e+00/   
      data psa(1,2,1,3,1)/ 0.15000e+03/                                 
      data psa(2,2,1,3,1)/ 0.52700e+01/, psa(3,2,1,3,1)/ 0.08000e+00/   
      data psa(1,3,1,3,1)/ 0.20000e+03/                                 
      data psa(2,3,1,3,1)/ 0.67600e+01/, psa(3,3,1,3,1)/ 0.11000e+00/   
      data psa(1,4,1,3,1)/ 0.30000e+03/                                 
      data psa(2,4,1,3,1)/ 0.94600e+01/, psa(3,4,1,3,1)/ 0.11000e+00/   
c        3p2                                                            
      data psa(1,1,3,3,1)/ 0.15000e+02/                                 
      data psa(2,1,3,3,1)/ 0.11500e+01/, psa(3,1,3,3,1)/ 0.14000e+00/   
      data psa(1,2,3,3,1)/ 0.25000e+02/                                 
      data psa(2,2,3,3,1)/ 0.23700e+01/, psa(3,2,3,3,1)/ 0.05000e+00/   
      data psa(1,3,3,3,1)/ 0.50000e+02/                                 
      data psa(2,3,3,3,1)/ 0.57680e+01/, psa(3,3,3,3,1)/ 0.02000e+00/   
      data psa(1,4,3,3,1)/ 0.10000e+03/                                 
      data psa(2,4,3,3,1)/ 0.10440e+02/, psa(3,4,3,3,1)/ 0.24000e+00/   
      data psa(1,5,3,3,1)/ 0.15000e+03/                                 
      data psa(2,5,3,3,1)/ 0.14070e+02/, psa(3,5,3,3,1)/ 0.05000e+00/   
      data psa(1,6,3,3,1)/ 0.20000e+03/                                 
      data psa(2,6,3,3,1)/ 0.15620e+02/, psa(3,6,3,3,1)/ 0.10000e+00/   
      data psa(1,7,3,3,1)/ 0.30000e+03/                                 
      data psa(2,7,3,3,1)/ 0.17330e+02/, psa(3,7,3,3,1)/ 0.14000e+00/   
c        3f2                                                            
      data psa(1,1,4,3,1)/ 0.10000e+03/                                 
      data psa(2,1,4,3,1)/ 0.90000e+00/, psa(3,1,4,3,1)/ 0.14000e+00/   
      data psa(1,2,4,3,1)/ 0.15000e+03/                                 
      data psa(2,2,4,3,1)/ 0.12200e+01/, psa(3,2,4,3,1)/ 0.05000e+00/   
      data psa(1,3,4,3,1)/ 0.20000e+03/                                 
      data psa(2,3,4,3,1)/ 0.11200e+01/, psa(3,3,4,3,1)/ 0.10000e+00/   
      data psa(1,4,4,3,1)/ 0.30000e+03/                                 
      data psa(2,4,4,3,1)/ 0.72000e+00/, psa(3,4,4,3,1)/ 0.15000e+00/   
c        e2                                                             
      data psa(1,1,5,3,1)/ 0.25000e+02/                                 
      data psa(2,1,5,3,1)/-0.73000e+00/, psa(3,1,5,3,1)/ 0.04000e+00/   
      data psa(1,2,5,3,1)/ 0.50000e+02/                                 
      data psa(2,2,5,3,1)/-0.15800e+01/, psa(3,2,5,3,1)/ 0.01000e+00/   
      data psa(1,3,5,3,1)/ 0.10000e+03/                                 
      data psa(2,3,5,3,1)/-0.27300e+01/, psa(3,3,5,3,1)/ 0.15000e+00/   
      data psa(1,4,5,3,1)/ 0.15000e+03/                                 
      data psa(2,4,5,3,1)/-0.28000e+01/, psa(3,4,5,3,1)/ 0.05000e+00/   
      data psa(1,5,5,3,1)/ 0.20000e+03/                                 
      data psa(2,5,5,3,1)/-0.28600e+01/, psa(3,5,5,3,1)/ 0.09000e+00/   
      data psa(1,6,5,3,1)/ 0.30000e+03/                                 
      data psa(2,6,5,3,1)/-0.24200e+01/, psa(3,6,5,3,1)/ 0.09000e+00/   
c        j = 3                                                          
c        3f3                                                            
      data psa(1,1,2,4,1)/ 0.20000e+03/                                 
      data psa(2,1,2,4,1)/-0.24600e+01/, psa(3,1,2,4,1)/ 0.11000e+00/   
      data psa(1,2,2,4,1)/ 0.30000e+03/                                 
      data psa(2,2,2,4,1)/-0.27700e+01/, psa(3,2,2,4,1)/ 0.14000e+00/   
c        j = 4                                                          
c        1g4                                                            
c        no data                                                        
c        3f4                                                            
      data psa(1,1,3,5,1)/ 0.30000e+03/                                 
      data psa(2,1,3,5,1)/ 0.27300e+01/, psa(3,1,3,5,1)/ 0.09000e+00/   
c        3h4                                                            
c        no data                                                        
c        e4                                                             
      data psa(1,1,5,5,1)/ 0.30000e+03/                                 
      data psa(2,1,5,5,1)/-0.14400e+01/, psa(3,1,5,5,1)/ 0.50000e-01/   
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c        phase-shift analysis vs35, ips=2                               
c        --------------------------------                               
c        vs35: arndt analysis vs35 for n-p, said program 5/27/91        
c                                                                       
c                                                                       
c                                                                       
c        j = 0                                                          
c        1s0                                                            
      data psa(1,1,1,1,2)/ 0.15000e+02/                                 
      data psa(2,1,1,1,2)/ 0.55980e+02/, psa(3,1,1,1,2)/ 0.55000e+00/   
      data psa(1,2,1,1,2)/ 0.25000e+02/                                 
      data psa(2,2,1,1,2)/ 0.50450e+02/, psa(3,2,1,1,2)/ 0.23000e+00/   
      data psa(1,3,1,1,2)/ 0.50000e+02/                                 
      data psa(2,3,1,1,2)/ 0.40260e+02/, psa(3,3,1,1,2)/ 0.09000e+00/   
      data psa(1,4,1,1,2)/ 0.10000e+03/                                 
      data psa(2,4,1,1,2)/ 0.25330e+02/, psa(3,4,1,1,2)/ 0.55000e+00/   
      data psa(1,5,1,1,2)/ 0.15000e+03/                                 
      data psa(2,5,1,1,2)/ 0.14870e+02/, psa(3,5,1,1,2)/ 0.34000e+00/   
      data psa(1,6,1,1,2)/ 0.20000e+03/                                 
      data psa(2,6,1,1,2)/ 0.75000e+01/, psa(3,6,1,1,2)/ 0.38000e+00/   
      data psa(1,7,1,1,2)/ 0.30000e+03/                                 
      data psa(2,7,1,1,2)/-0.74900e+01/, psa(3,7,1,1,2)/ 0.32000e+00/   
c        3p0                                                            
      data psa(1,1,3,1,2)/ 0.15000e+02/                                 
      data psa(2,1,3,1,2)/ 0.57400e+01/, psa(3,1,3,1,2)/ 0.14800e+01/   
      data psa(1,2,3,1,2)/ 0.25000e+02/                                 
      data psa(2,2,3,1,2)/ 0.89600e+01/, psa(3,2,3,1,2)/ 0.42000e+00/   
      data psa(1,3,3,1,2)/ 0.50000e+02/                                 
      data psa(2,3,3,1,2)/ 0.11860e+02/, psa(3,3,3,1,2)/ 0.12000e+00/   
      data psa(1,4,3,1,2)/ 0.10000e+03/                                 
      data psa(2,4,3,1,2)/ 0.11780e+02/, psa(3,4,3,1,2)/ 0.16600e+01/   
      data psa(1,5,3,1,2)/ 0.15000e+03/                                 
      data psa(2,5,3,1,2)/ 0.53900e+01/, psa(3,5,3,1,2)/ 0.28000e+00/   
      data psa(1,6,3,1,2)/ 0.20000e+03/                                 
      data psa(2,6,3,1,2)/ 0.22000e+00/, psa(3,6,3,1,2)/ 0.33000e+00/   
      data psa(1,7,3,1,2)/ 0.30000e+03/                                 
      data psa(2,7,3,1,2)/-0.10130e+02/, psa(3,7,3,1,2)/ 0.40000e+00/   
c        j = 1                                                          
c        1p1                                                            
      data psa(1,1,1,2,2)/ 0.15000e+02/                                 
      data psa(2,1,1,2,2)/-0.38800e+01/, psa(3,1,1,2,2)/ 0.10900e+01/   
      data psa(1,2,1,2,2)/ 0.25000e+02/                                 
      data psa(2,2,1,2,2)/-0.65100e+01/, psa(3,2,1,2,2)/ 0.18000e+00/   
      data psa(1,3,1,2,2)/ 0.50000e+02/                                 
      data psa(2,3,1,2,2)/-0.86400e+01/, psa(3,3,1,2,2)/ 0.20000e+00/   
      data psa(1,4,1,2,2)/ 0.10000e+03/                                 
      data psa(2,4,1,2,2)/-0.12350e+02/, psa(3,4,1,2,2)/ 0.16700e+01/   
      data psa(1,5,1,2,2)/ 0.15000e+03/                                 
      data psa(2,5,1,2,2)/-0.16300e+02/, psa(3,5,1,2,2)/ 0.11100e+01/   
      data psa(1,6,1,2,2)/ 0.20000e+03/                                 
      data psa(2,6,1,2,2)/-0.22690e+02/, psa(3,6,1,2,2)/ 0.52000e+00/   
      data psa(1,7,1,2,2)/ 0.30000e+03/                                 
      data psa(2,7,1,2,2)/-0.28730e+02/, psa(3,7,1,2,2)/ 0.58000e+00/   
c        3p1                                                            
      data psa(1,1,2,2,2)/ 0.15000e+02/                                 
      data psa(2,1,2,2,2)/-0.32800e+01/, psa(3,1,2,2,2)/ 0.77000e+00/   
      data psa(1,2,2,2,2)/ 0.25000e+02/                                 
      data psa(2,2,2,2,2)/-0.50500e+01/, psa(3,2,2,2,2)/ 0.20000e+00/   
      data psa(1,3,2,2,2)/ 0.50000e+02/                                 
      data psa(2,3,2,2,2)/-0.85000e+01/, psa(3,3,2,2,2)/ 0.05000e+00/   
      data psa(1,4,2,2,2)/ 0.10000e+03/                                 
      data psa(2,4,2,2,2)/-0.13890e+02/, psa(3,4,2,2,2)/ 0.33000e+00/   
      data psa(1,5,2,2,2)/ 0.15000e+03/                                 
      data psa(2,5,2,2,2)/-0.17990e+02/, psa(3,5,2,2,2)/ 0.08000e+00/   
      data psa(1,6,2,2,2)/ 0.20000e+03/                                 
      data psa(2,6,2,2,2)/-0.21740e+02/, psa(3,6,2,2,2)/ 0.17000e+00/   
      data psa(1,7,2,2,2)/ 0.30000e+03/                                 
      data psa(2,7,2,2,2)/-0.28780e+02/, psa(3,7,2,2,2)/ 0.24000e+00/   
c        3s1                                                            
      data psa(1,1,3,2,2)/ 0.15000e+02/                                 
      data psa(2,1,3,2,2)/ 0.91170e+02/, psa(3,1,3,2,2)/ 0.28600e+01/   
      data psa(1,2,3,2,2)/ 0.25000e+02/                                 
      data psa(2,2,3,2,2)/ 0.80780e+02/, psa(3,2,3,2,2)/ 0.55000e+00/   
      data psa(1,3,3,2,2)/ 0.50000e+02/                                 
      data psa(2,3,3,2,2)/ 0.62120e+02/, psa(3,3,3,2,2)/ 0.32000e+00/   
      data psa(1,4,3,2,2)/ 0.10000e+03/                                 
      data psa(2,4,3,2,2)/ 0.42750e+02/, psa(3,4,3,2,2)/ 0.69000e+00/   
      data psa(1,5,3,2,2)/ 0.15000e+03/                                 
      data psa(2,5,3,2,2)/ 0.27860e+02/, psa(3,5,3,2,2)/ 0.47000e+00/   
      data psa(1,6,3,2,2)/ 0.20000e+03/                                 
      data psa(2,6,3,2,2)/ 0.20530e+02/, psa(3,6,3,2,2)/ 0.35000e+00/   
      data psa(1,7,3,2,2)/ 0.30000e+03/                                 
      data psa(2,7,3,2,2)/ 0.55900e+01/, psa(3,7,3,2,2)/ 0.60000e+00/   
c        3d1                                                            
      data psa(1,1,4,2,2)/ 0.15000e+02/                                 
      data psa(2,1,4,2,2)/-0.12000e+01/, psa(3,1,4,2,2)/ 0.04000e+00/   
      data psa(1,2,4,2,2)/ 0.10000e+03/                                 
      data psa(2,2,4,2,2)/-0.13110e+02/, psa(3,2,4,2,2)/ 0.52000e+00/   
      data psa(1,3,4,2,2)/ 0.15000e+03/                                 
      data psa(2,3,4,2,2)/-0.15340e+02/, psa(3,3,4,2,2)/ 0.38000e+00/   
      data psa(1,4,4,2,2)/ 0.20000e+03/                                 
      data psa(2,4,4,2,2)/-0.18960e+02/, psa(3,4,4,2,2)/ 0.29000e+00/   
      data psa(1,5,4,2,2)/ 0.30000e+03/                                 
      data psa(2,5,4,2,2)/-0.24910e+02/, psa(3,5,4,2,2)/ 0.42000e+00/   
c        e1                                                             
      data psa(1,1,5,2,2)/ 0.25000e+02/                                 
      data psa(2,1,5,2,2)/ 0.20500e+01/, psa(3,1,5,2,2)/ 0.31000e+00/   
      data psa(1,2,5,2,2)/ 0.50000e+02/                                 
      data psa(2,2,5,2,2)/ 0.28700e+01/, psa(3,2,5,2,2)/ 0.21000e+00/   
      data psa(1,3,5,2,2)/ 0.10000e+03/                                 
      data psa(2,3,5,2,2)/-0.04000e+00/, psa(3,3,5,2,2)/ 0.87000e+00/   
      data psa(1,4,5,2,2)/ 0.15000e+03/                                 
      data psa(2,4,5,2,2)/ 0.41100e+01/, psa(3,4,5,2,2)/ 0.50000e+00/   
      data psa(1,5,5,2,2)/ 0.20000e+03/                                 
      data psa(2,5,5,2,2)/ 0.42400e+01/, psa(3,5,5,2,2)/ 0.29000e+00/   
      data psa(1,6,5,2,2)/ 0.30000e+03/                                 
      data psa(2,6,5,2,2)/ 0.58700e+01/, psa(3,6,5,2,2)/ 0.28000e+00/   
c        j = 2                                                          
c        1d2                                                            
      data psa(1,1,1,3,2)/ 0.10000e+03/                                 
      data psa(2,1,1,3,2)/ 0.38500e+01/, psa(3,1,1,3,2)/ 0.15000e+00/   
      data psa(1,2,1,3,2)/ 0.15000e+03/                                 
      data psa(2,2,1,3,2)/ 0.54700e+01/, psa(3,2,1,3,2)/ 0.08000e+00/   
      data psa(1,3,1,3,2)/ 0.20000e+03/                                 
      data psa(2,3,1,3,2)/ 0.70200e+01/, psa(3,3,1,3,2)/ 0.11000e+00/   
      data psa(1,4,1,3,2)/ 0.30000e+03/                                 
      data psa(2,4,1,3,2)/ 0.97900e+01/, psa(3,4,1,3,2)/ 0.11000e+00/   
c        3d2                                                            
      data psa(1,1,2,3,2)/ 0.50000e+02/                                 
      data psa(2,1,2,3,2)/ 0.92300e+01/, psa(3,1,2,3,2)/ 0.11000e+00/   
      data psa(1,2,2,3,2)/ 0.10000e+03/                                 
      data psa(2,2,2,3,2)/ 0.17770e+02/, psa(3,2,2,3,2)/ 0.97000e+00/   
      data psa(1,3,2,3,2)/ 0.15000e+03/                                 
      data psa(2,3,2,3,2)/ 0.24730e+02/, psa(3,3,2,3,2)/ 0.56000e+00/   
      data psa(1,4,2,3,2)/ 0.20000e+03/                                 
      data psa(2,4,2,3,2)/ 0.24110e+02/, psa(3,4,2,3,2)/ 0.27000e+00/   
      data psa(1,5,2,3,2)/ 0.30000e+03/                                 
      data psa(2,5,2,3,2)/ 0.23930e+02/, psa(3,5,2,3,2)/ 0.28000e+00/   
c        3p2                                                            
      data psa(1,1,3,3,2)/ 0.15000e+02/                                 
      data psa(2,1,3,3,2)/ 0.12600e+01/, psa(3,1,3,3,2)/ 0.14000e+00/   
      data psa(1,2,3,3,2)/ 0.25000e+02/                                 
      data psa(2,2,3,3,2)/ 0.25400e+01/, psa(3,2,3,3,2)/ 0.05000e+00/   
      data psa(1,3,3,3,2)/ 0.50000e+02/                                 
      data psa(2,3,3,3,2)/ 0.60400e+01/, psa(3,3,3,3,2)/ 0.02000e+00/   
      data psa(1,4,3,3,2)/ 0.10000e+03/                                 
      data psa(2,4,3,3,2)/ 0.10850e+02/, psa(3,4,3,3,2)/ 0.24000e+00/   
      data psa(1,5,3,3,2)/ 0.15000e+03/                                 
      data psa(2,5,3,3,2)/ 0.14560e+02/, psa(3,5,3,3,2)/ 0.05000e+00/   
      data psa(1,6,3,3,2)/ 0.20000e+03/                                 
      data psa(2,6,3,3,2)/ 0.16170e+02/, psa(3,6,3,3,2)/ 0.10000e+00/   
      data psa(1,7,3,3,2)/ 0.30000e+03/                                 
      data psa(2,7,3,3,2)/ 0.17930e+02/, psa(3,7,3,3,2)/ 0.14000e+00/   
c        3f2                                                            
      data psa(1,1,4,3,2)/ 0.10000e+03/                                 
      data psa(2,1,4,3,2)/ 0.86000e+00/, psa(3,1,4,3,2)/ 0.14000e+00/   
      data psa(1,2,4,3,2)/ 0.15000e+03/                                 
      data psa(2,2,4,3,2)/ 0.11800e+01/, psa(3,2,4,3,2)/ 0.05000e+00/   
      data psa(1,3,4,3,2)/ 0.20000e+03/                                 
      data psa(2,3,4,3,2)/ 0.10800e+01/, psa(3,3,4,3,2)/ 0.10000e+00/   
      data psa(1,4,4,3,2)/ 0.30000e+03/                                 
      data psa(2,4,4,3,2)/ 0.68000e+00/, psa(3,4,4,3,2)/ 0.15000e+00/   
c        e2                                                             
      data psa(1,1,5,3,2)/ 0.25000e+02/                                 
      data psa(2,1,5,3,2)/-0.73000e+00/, psa(3,1,5,3,2)/ 0.04000e+00/   
      data psa(1,2,5,3,2)/ 0.50000e+02/                                 
      data psa(2,2,5,3,2)/-0.15900e+01/, psa(3,2,5,3,2)/ 0.01000e+00/   
      data psa(1,3,5,3,2)/ 0.10000e+03/                                 
      data psa(2,3,5,3,2)/-0.27700e+01/, psa(3,3,5,3,2)/ 0.15000e+00/   
      data psa(1,4,5,3,2)/ 0.15000e+03/                                 
      data psa(2,4,5,3,2)/-0.28700e+01/, psa(3,4,5,3,2)/ 0.05000e+00/   
      data psa(1,5,5,3,2)/ 0.20000e+03/                                 
      data psa(2,5,5,3,2)/-0.29400e+01/, psa(3,5,5,3,2)/ 0.09000e+00/   
      data psa(1,6,5,3,2)/ 0.30000e+03/                                 
      data psa(2,6,5,3,2)/-0.25300e+01/, psa(3,6,5,3,2)/ 0.09000e+00/   
c        j = 3                                                          
c        1f3                                                            
      data psa(1,1,1,4,2)/ 0.30000e+03/                                 
      data psa(2,1,1,4,2)/-0.59000e+01/, psa(3,1,1,4,2)/ 0.14000e+00/   
c        3f3                                                            
      data psa(1,1,2,4,2)/ 0.20000e+03/                                 
      data psa(2,1,2,4,2)/-0.24400e+01/, psa(3,1,2,4,2)/ 0.11000e+00/   
      data psa(1,2,2,4,2)/ 0.30000e+03/                                 
      data psa(2,2,2,4,2)/-0.27800e+01/, psa(3,2,2,4,2)/ 0.14000e+00/   
c        3d3                                                            
      data psa(1,1,3,4,2)/ 0.10000e+03/                                 
      data psa(2,1,3,4,2)/ 0.11800e+01/, psa(3,1,3,4,2)/ 0.46000e+00/   
      data psa(1,2,3,4,2)/ 0.15000e+03/                                 
      data psa(2,2,3,4,2)/ 0.65000e+00/, psa(3,2,3,4,2)/ 0.47000e+00/   
      data psa(1,3,3,4,2)/ 0.20000e+03/                                 
      data psa(2,3,3,4,2)/ 0.40600e+01/, psa(3,3,3,4,2)/ 0.23000e+00/   
      data psa(1,4,3,4,2)/ 0.30000e+03/                                 
      data psa(2,4,3,4,2)/ 0.42300e+01/, psa(3,4,3,4,2)/ 0.28000e+00/   
c        3g3                                                            
      data psa(1,1,4,4,2)/ 0.15000e+03/                                 
      data psa(2,1,4,4,2)/-0.29900e+01/, psa(3,1,4,4,2)/ 0.35000e+00/   
      data psa(1,2,4,4,2)/ 0.20000e+03/                                 
      data psa(2,2,4,4,2)/-0.31100e+01/, psa(3,2,4,4,2)/ 0.20000e+00/   
      data psa(1,3,4,4,2)/ 0.30000e+03/                                 
      data psa(2,3,4,4,2)/-0.50400e+01/, psa(3,3,4,4,2)/ 0.26000e+00/   
c        e3                                                             
      data psa(1,1,5,4,2)/ 0.10000e+03/                                 
      data psa(2,1,5,4,2)/ 0.28300e+01/, psa(3,1,5,4,2)/ 0.29000e+00/   
      data psa(1,2,5,4,2)/ 0.15000e+03/                                 
      data psa(2,2,5,4,2)/ 0.31200e+01/, psa(3,2,5,4,2)/ 0.37000e+00/   
      data psa(1,3,5,4,2)/ 0.20000e+03/                                 
      data psa(2,3,5,4,2)/ 0.56500e+01/, psa(3,3,5,4,2)/ 0.13000e+00/   
      data psa(1,4,5,4,2)/ 0.30000e+03/                                 
      data psa(2,4,5,4,2)/ 0.65400e+01/, psa(3,4,5,4,2)/ 0.13000e+00/   
c        j = 4                                                          
c        1g4                                                            
c      no data                                                          
c        3g4                                                            
      data psa(1,1,2,5,2)/ 0.15000e+03/                                 
      data psa(2,1,2,5,2)/ 0.54500e+01/, psa(3,1,2,5,2)/ 0.41000e+00/   
      data psa(1,2,2,5,2)/ 0.20000e+03/                                 
      data psa(2,2,2,5,2)/ 0.54700e+01/, psa(3,2,2,5,2)/ 0.20000e+00/   
      data psa(1,3,2,5,2)/ 0.30000e+03/                                 
      data psa(2,3,2,5,2)/ 0.74200e+01/, psa(3,3,2,5,2)/ 0.20000e+00/   
c        3f4                                                            
      data psa(1,1,3,5,2)/ 0.30000e+03/                                 
      data psa(2,1,3,5,2)/ 0.27300e+01/, psa(3,1,3,5,2)/ 0.09000e+00/   
c        e4                                                             
      data psa(1,1,5,5,2)/ 0.30000e+03/                                 
      data psa(2,1,5,5,2)/-0.14300e+01/, psa(3,1,5,5,2)/ 0.05000e+00/   
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c        phase-shift analysis bg90, ips=3                               
c        --------------------------------                               
c        bg90: n-p, d.v. bugg, phys. rev. c 41(1990) 2708               
c                                                                       
c                                                                       
c        j = 0                                                          
c        1s0 np                                                         
      data psa(1,1,1,1,3)/ 0.14200e+03/                                 
      data psa(2,1,1,1,3)/ 0.16060e+02/, psa(3,1,1,1,3)/ 0.56000e+00/   
      data psa(1,2,1,1,3)/ 0.21000e+03/                                 
      data psa(2,2,1,1,3)/ 0.37700e+01/, psa(3,2,1,1,3)/ 0.31000e+00/   
      data psa(1,3,1,1,3)/ 0.32500e+03/                                 
      data psa(2,3,1,1,3)/-0.95400e+01/, psa(3,3,1,1,3)/ 0.24000e+00/   
c        3p0                                                            
      data psa(1,1,3,1,3)/ 0.14200e+03/                                 
      data psa(2,1,3,1,3)/ 0.56700e+01/, psa(3,1,3,1,3)/ 0.45000e+00/   
      data psa(1,2,3,1,3)/ 0.21000e+03/                                 
      data psa(2,2,3,1,3)/-0.21800e+01/, psa(3,2,3,1,3)/ 0.31000e+00/   
      data psa(1,3,3,1,3)/ 0.32500e+03/                                 
      data psa(2,3,3,1,3)/-0.13840e+02/, psa(3,3,3,1,3)/ 0.25000e+00/   
c        j = 1                                                          
c        1p1                                                            
      data psa(1,1,1,2,3)/ 0.14200e+03/                                 
      data psa(2,1,1,2,3)/-0.15620e+02/, psa(3,1,1,2,3)/ 0.81000e+00/   
      data psa(1,2,1,2,3)/ 0.21000e+03/                                 
      data psa(2,2,1,2,3)/-0.25180e+02/, psa(3,2,1,2,3)/ 0.43000e+00/   
      data psa(1,3,1,2,3)/ 0.32500e+03/                                 
      data psa(2,3,1,2,3)/-0.32350e+02/, psa(3,3,1,2,3)/ 0.42000e+00/   
c        3p1                                                            
      data psa(1,1,2,2,3)/ 0.14200e+03/                                 
      data psa(2,1,2,2,3)/-0.17260e+02/, psa(3,1,2,2,3)/ 0.14000e+00/   
      data psa(1,2,2,2,3)/ 0.21000e+03/                                 
      data psa(2,2,2,2,3)/-0.21710e+02/, psa(3,2,2,2,3)/ 0.15000e+00/   
      data psa(1,3,2,2,3)/ 0.32500e+03/                                 
      data psa(2,3,2,2,3)/-0.30490e+02/, psa(3,3,2,2,3)/ 0.16000e+00/   
c        3s1                                                            
      data psa(1,1,3,2,3)/ 0.14200e+03/                                 
      data psa(2,1,3,2,3)/ 0.31920e+02/, psa(3,1,3,2,3)/ 0.67000e+00/   
      data psa(1,2,3,2,3)/ 0.21000e+03/                                 
      data psa(2,2,3,2,3)/ 0.18740e+02/, psa(3,2,3,2,3)/ 0.22000e+00/   
      data psa(1,3,3,2,3)/ 0.32500e+03/                                 
      data psa(2,3,3,2,3)/ 0.23400e+01/, psa(3,3,3,2,3)/ 0.30000e+00/   
c        3d1                                                            
      data psa(1,1,4,2,3)/ 0.14200e+03/                                 
      data psa(2,1,4,2,3)/-0.15830e+02/, psa(3,1,4,2,3)/ 0.82000e+00/   
      data psa(1,2,4,2,3)/ 0.21000e+03/                                 
      data psa(2,2,4,2,3)/-0.19260e+02/, psa(3,2,4,2,3)/ 0.19000e+00/   
      data psa(1,3,4,2,3)/ 0.32500e+03/                                 
      data psa(2,3,4,2,3)/-0.25360e+02/, psa(3,3,4,2,3)/ 0.17000e+00/   
c        e1                                                             
      data psa(1,1,5,2,3)/ 0.14200e+03/                                 
      data psa(2,1,5,2,3)/ 0.22900e+01/, psa(3,1,5,2,3)/ 0.83000e+00/   
      data psa(1,2,5,2,3)/ 0.21000e+03/                                 
      data psa(2,2,5,2,3)/ 0.39200e+01/, psa(3,2,5,2,3)/ 0.19000e+00/   
      data psa(1,3,5,2,3)/ 0.32500e+03/                                 
      data psa(2,3,5,2,3)/ 0.45800e+01/, psa(3,3,5,2,3)/ 0.20000e+00/   
c        j = 2                                                          
c        1d2                                                            
      data psa(1,1,1,3,3)/ 0.14200e+03/                                 
      data psa(2,1,1,3,3)/ 0.52300e+01/, psa(3,1,1,3,3)/ 0.18000e+00/   
      data psa(1,2,1,3,3)/ 0.21000e+03/                                 
      data psa(2,2,1,3,3)/ 0.76100e+01/, psa(3,2,1,3,3)/ 0.10000e+00/   
      data psa(1,3,1,3,3)/ 0.32500e+03/                                 
      data psa(2,3,1,3,3)/ 0.98500e+01/, psa(3,3,1,3,3)/ 0.09000e+00/   
c        3d2                                                            
      data psa(1,1,2,3,3)/ 0.14200e+03/                                 
      data psa(2,1,2,3,3)/ 0.22700e+02/, psa(3,1,2,3,3)/ 0.66000e+00/   
      data psa(1,2,2,3,3)/ 0.21000e+03/                                 
      data psa(2,2,2,3,3)/ 0.24500e+02/, psa(3,2,2,3,3)/ 0.22000e+00/   
      data psa(1,3,2,3,3)/ 0.32500e+03/                                 
      data psa(2,3,2,3,3)/ 0.24230e+02/, psa(3,3,2,3,3)/ 0.18000e+00/   
c        3p2                                                            
      data psa(1,1,3,3,3)/ 0.14200e+03/                                 
      data psa(2,1,3,3,3)/ 0.13950e+02/, psa(3,1,3,3,3)/ 0.10000e+00/   
      data psa(1,2,3,3,3)/ 0.21000e+03/                                 
      data psa(2,2,3,3,3)/ 0.16040e+02/, psa(3,2,3,3,3)/ 0.13000e+00/   
      data psa(1,3,3,3,3)/ 0.32500e+03/                                 
      data psa(2,3,3,3,3)/ 0.17030e+02/, psa(3,3,3,3,3)/ 0.11000e+00/   
c        3f2                                                            
      data psa(1,1,4,3,3)/ 0.14200e+03/                                 
      data psa(2,1,4,3,3)/ 0.92000e+00/, psa(3,1,4,3,3)/ 0.21000e+00/   
      data psa(1,2,4,3,3)/ 0.21000e+03/                                 
      data psa(2,2,4,3,3)/ 0.13800e+01/, psa(3,2,4,3,3)/ 0.15000e+00/   
      data psa(1,3,4,3,3)/ 0.32500e+03/                                 
      data psa(2,3,4,3,3)/ 0.77000e+00/, psa(3,3,4,3,3)/ 0.11000e+00/   
c        e2                                                             
      data psa(1,1,5,3,3)/ 0.14200e+03/                                 
      data psa(2,1,5,3,3)/-0.28000e+01/, psa(3,1,5,3,3)/ 0.70000e-01/   
      data psa(1,2,5,3,3)/ 0.21000e+03/                                 
      data psa(2,2,5,3,3)/-0.25600e+01/, psa(3,2,5,3,3)/ 0.08000e+00/   
      data psa(1,3,5,3,3)/ 0.32500e+03/                                 
      data psa(2,3,5,3,3)/-0.22900e+01/, psa(3,3,5,3,3)/ 0.09000e+00/   
c        j = 3                                                          
c        1f3                                                            
      data psa(1,1,1,4,3)/ 0.14200e+03/                                 
      data psa(2,1,1,4,3)/-0.31100e+01/, psa(3,1,1,4,3)/ 0.15000e+00/   
      data psa(1,2,1,4,3)/ 0.21000e+03/                                 
      data psa(2,2,1,4,3)/-0.39700e+01/, psa(3,2,1,4,3)/ 0.23000e+00/   
      data psa(1,3,1,4,3)/ 0.32500e+03/                                 
      data psa(2,3,1,4,3)/-0.58400e+01/, psa(3,3,1,4,3)/ 0.12000e+00/   
c        3f3                                                            
      data psa(1,1,2,4,3)/ 0.14200e+03/                                 
      data psa(2,1,2,4,3)/-0.19900e+01/, psa(3,1,2,4,3)/ 0.17000e+00/   
      data psa(1,2,2,4,3)/ 0.21000e+03/                                 
      data psa(2,2,2,4,3)/-0.29900e+01/, psa(3,2,2,4,3)/ 0.11000e+00/   
      data psa(1,3,2,4,3)/ 0.32500e+03/                                 
      data psa(2,3,2,4,3)/-0.28200e+01/, psa(3,3,2,4,3)/ 0.09000e+00/   
c        3d3                                                            
      data psa(1,1,3,4,3)/ 0.14200e+03/                                 
      data psa(2,1,3,4,3)/ 0.16700e+01/, psa(3,1,3,4,3)/ 0.41000e+00/   
      data psa(1,2,3,4,3)/ 0.21000e+03/                                 
      data psa(2,2,3,4,3)/ 0.32200e+01/, psa(3,2,3,4,3)/ 0.20000e+00/   
      data psa(1,3,3,4,3)/ 0.32500e+03/                                 
      data psa(2,3,3,4,3)/ 0.33600e+01/, psa(3,3,3,4,3)/ 0.20000e+00/   
c        3g3                                                            
      data psa(1,1,4,4,3)/ 0.14200e+03/                                 
      data psa(2,1,4,4,3)/-0.16200e+01/, psa(3,1,4,4,3)/ 0.30000e+00/   
      data psa(1,2,4,4,3)/ 0.21000e+03/                                 
      data psa(2,2,4,4,3)/-0.35000e+01/, psa(3,2,4,4,3)/ 0.23000e+00/   
      data psa(1,3,4,4,3)/ 0.32500e+03/                                 
      data psa(2,3,4,4,3)/-0.39700e+01/, psa(3,3,4,4,3)/ 0.17000e+00/   
c        e3                                                             
      data psa(1,1,5,4,3)/ 0.14200e+03/                                 
      data psa(2,1,5,4,3)/ 0.46100e+01/, psa(3,1,5,4,3)/ 0.24000e+00/   
      data psa(1,2,5,4,3)/ 0.21000e+03/                                 
      data psa(2,2,5,4,3)/ 0.60200e+01/, psa(3,2,5,4,3)/ 0.08000e+00/   
      data psa(1,3,5,4,3)/ 0.32500e+03/                                 
      data psa(2,3,5,4,3)/ 0.73700e+01/, psa(3,3,5,4,3)/ 0.08000e+00/   
c        j = 4                                                          
c        1g4                                                            
      data psa(1,1,1,5,3)/ 0.14200e+03/                                 
      data psa(2,1,1,5,3)/ 0.72000e+00/, psa(3,1,1,5,3)/ 0.70000e-01/   
      data psa(1,2,1,5,3)/ 0.21000e+03/                                 
      data psa(2,2,1,5,3)/ 0.98000e+00/, psa(3,2,1,5,3)/ 0.06000e+00/   
      data psa(1,3,1,5,3)/ 0.32500e+03/                                 
      data psa(2,3,1,5,3)/ 0.15600e+01/, psa(3,3,1,5,3)/ 0.06000e+00/   
c        3g4                                                            
      data psa(1,1,2,5,3)/ 0.14200e+03/                                 
      data psa(2,1,2,5,3)/ 0.39500e+01/, psa(3,1,2,5,3)/ 0.58000e+00/   
      data psa(1,2,2,5,3)/ 0.21000e+03/                                 
      data psa(2,2,2,5,3)/ 0.54400e+01/, psa(3,2,2,5,3)/ 0.19000e+00/   
      data psa(1,3,2,5,3)/ 0.32500e+03/                                 
      data psa(2,3,2,5,3)/ 0.71700e+01/, psa(3,3,2,5,3)/ 0.15000e+00/   
c        3f4                                                            
      data psa(1,1,3,5,3)/ 0.14200e+03/                                 
      data psa(2,1,3,5,3)/ 0.85000e+00/, psa(3,1,3,5,3)/ 0.10000e+00/   
      data psa(1,2,3,5,3)/ 0.21000e+03/                                 
      data psa(2,2,3,5,3)/ 0.17900e+01/, psa(3,2,3,5,3)/ 0.10000e+00/   
      data psa(1,3,3,5,3)/ 0.32500e+03/                                 
      data psa(2,3,3,5,3)/ 0.29600e+01/, psa(3,3,3,5,3)/ 0.07000e+00/   
c        3h4                                                            
      data psa(1,1,4,5,3)/ 0.14200e+03/                                 
      data psa(2,1,4,5,3)/ 0.24000e+00/, psa(3,1,4,5,3)/ 0.03000e+00/   
      data psa(1,2,4,5,3)/ 0.21000e+03/                                 
      data psa(2,2,4,5,3)/ 0.33000e+00/, psa(3,2,4,5,3)/ 0.05000e+00/   
      data psa(1,3,4,5,3)/ 0.32500e+03/                                 
      data psa(2,3,4,5,3)/ 0.57000e+00/, psa(3,3,4,5,3)/ 0.05000e+00/   
c        e4                                                             
      data psa(1,1,5,5,3)/ 0.14200e+03/                                 
      data psa(2,1,5,5,3)/-0.83000e+00/, psa(3,1,5,5,3)/ 0.30000e-01/   
      data psa(1,2,5,5,3)/ 0.21000e+03/                                 
      data psa(2,2,5,5,3)/-0.11000e+01/, psa(3,2,5,5,3)/ 0.30000e-01/   
      data psa(1,3,5,5,3)/ 0.32500e+03/                                 
      data psa(2,3,5,5,3)/-0.14700e+01/, psa(3,3,5,5,3)/ 0.04000e+00/   
c        j = 5                                                          
c        3h5                                                            
      data psa(1,1,2,6,3)/ 0.14200e+03/                                 
      data psa(2,1,2,6,3)/-0.50000e+00/, psa(3,1,2,6,3)/ 0.03000e+00/   
      data psa(1,2,2,6,3)/ 0.21000e+03/                                 
      data psa(2,2,2,6,3)/-0.75000e+00/, psa(3,2,2,6,3)/ 0.04000e+00/   
      data psa(1,3,2,6,3)/ 0.32500e+03/                                 
      data psa(2,3,2,6,3)/-0.10400e+01/, psa(3,3,2,6,3)/ 0.05000e+00/   
c        3g5                                                            
      data psa(1,1,3,6,3)/ 0.14200e+03/                                 
      data psa(2,1,3,6,3)/-0.13000e+00/, psa(3,1,3,6,3)/ 0.22000e+00/   
      data psa(1,2,3,6,3)/ 0.21000e+03/                                 
      data psa(2,2,3,6,3)/ 0.24000e+00/, psa(3,2,3,6,3)/ 0.17000e+00/   
      data psa(1,3,3,6,3)/ 0.32500e+03/                                 
      data psa(2,3,3,6,3)/-0.54000e+00/, psa(3,3,3,6,3)/ 0.12000e+00/   
c        j = 6                                                          
c        3h6                                                            
      data psa(1,1,3,7,3)/ 0.14200e+03/                                 
      data psa(2,1,3,7,3)/ 0.13000e+00/, psa(3,1,3,7,3)/ 0.03000e+00/   
      data psa(1,2,3,7,3)/ 0.21000e+03/                                 
      data psa(2,2,3,7,3)/ 0.24000e+00/, psa(3,2,3,7,3)/ 0.04000e+00/   
      data psa(1,3,3,7,3)/ 0.32500e+03/                                 
      data psa(2,3,3,7,3)/ 0.56000e+00/, psa(3,3,3,7,3)/ 0.05000e+00/   
c                                                                       
c                                                                       
c *****  phase-shift analysis  ni90,ips=4                               
c        --------------------------------                               
c         ni90: nijmegen p-p single-energy analysis,                    
c               bergervoet et al., phys.rev. c 41 (1990), 1435          
c                                                                       
c                                                                       
c                                                                       
c        j = 0                                                          
c        1s0 pp                                                         
      data psa(1,1,1,1,4)/ 0.25000e+02/                                 
      data psa(2,1,1,1,4)/ 0.48836e+02/, psa(3,1,1,1,4)/ 0.13320e+00/   
      data psa(1,2,1,1,4)/ 0.50000e+02/                                 
      data psa(2,2,1,1,4)/ 0.39086e+02/, psa(3,2,1,1,4)/ 0.10340e+00/   
      data psa(1,3,1,1,4)/ 0.10000e+03/                                 
      data psa(2,3,1,1,4)/ 0.24722e+02/, psa(3,3,1,1,4)/ 0.47900e+00/   
      data psa(1,4,1,1,4)/ 0.15000e+03/                                 
      data psa(2,4,1,1,4)/ 0.14306e+02/, psa(3,4,1,1,4)/ 0.37400e+00/   
      data psa(1,5,1,1,4)/ 0.21500e+03/                                 
      data psa(2,5,1,1,4)/ 0.42050e+01/, psa(3,5,1,1,4)/ 0.31100e+00/   
      data psa(1,6,1,1,4)/ 0.32000e+03/                                 
      data psa(2,6,1,1,4)/-0.81420e+01/, psa(3,6,1,1,4)/ 0.34900e+00/   
      data psa(1,7,1,1,4)/ 0.10000e+02/                                 
      data psa(2,7,1,1,4)/ 0.55089e+02/, psa(3,7,1,1,4)/ 0.06810e+00/   
      data psa(1,8,1,1,4)/ 0.10000e+01/                                 
      data psa(2,8,1,1,4)/ 0.32681e+02/, psa(3,8,1,1,4)/ 0.00940e+00/   
      data psa(1,9,1,1,4)/ 0.50000e+01/                                 
      data psa(2,9,1,1,4)/ 0.54546e+02/, psa(3,9,1,1,4)/ 0.08710e+00/   
      data psa(1,10,1,1,4)/ 0.38254e+00/                                
      data psa(2,10,1,1,4)/ 0.14509e+02/, psa(3,10,1,1,4)/ 0.00250e+00/ 
c        3p0                                                            
      data psa(1,1,3,1,4)/ 0.25000e+02/                                 
      data psa(2,1,3,1,4)/ 0.84968e+01/, psa(3,1,3,1,4)/ 0.60570e+00/   
      data psa(1,2,3,1,4)/ 0.50000e+02/                                 
      data psa(2,2,3,1,4)/ 0.11254e+02/, psa(3,2,3,1,4)/ 0.30310e+00/   
      data psa(1,3,3,1,4)/ 0.10000e+03/                                 
      data psa(2,3,3,1,4)/ 0.10690e+02/, psa(3,3,3,1,4)/ 0.11970e+01/   
      data psa(1,4,3,1,4)/ 0.15000e+03/                                 
      data psa(2,4,3,1,4)/ 0.49660e+01/, psa(3,4,3,1,4)/ 0.27500e+00/   
      data psa(1,5,3,1,4)/ 0.21500e+03/                                 
      data psa(2,5,3,1,4)/-0.18100e+01/, psa(3,5,3,1,4)/ 0.39500e+00/   
      data psa(1,6,3,1,4)/ 0.32000e+03/                                 
      data psa(2,6,3,1,4)/-0.12495e+02/, psa(3,6,3,1,4)/ 0.46800e+00/   
      data psa(1,7,3,1,4)/ 0.10000e+02/                                 
      data psa(2,7,3,1,4)/ 0.35426e+01/, psa(3,7,3,1,4)/ 0.07290e+00/   
      data psa(1,8,3,1,4)/ 0.50000e+01/                                 
      data psa(2,8,3,1,4)/ 0.16144e+01/, psa(3,8,3,1,4)/ 0.09270e+00/   
c        j = 1                                                          
c        3p1                                                            
      data psa(1,1,2,2,4)/ 0.25000e+02/                                 
      data psa(2,1,2,2,4)/-0.47272e+01/, psa(3,1,2,2,4)/ 0.19640e+00/   
      data psa(1,2,2,2,4)/ 0.50000e+02/                                 
      data psa(2,2,2,2,4)/-0.82992e+01/, psa(3,2,2,2,4)/ 0.04800e+00/   
      data psa(1,3,2,2,4)/ 0.10000e+03/                                 
      data psa(2,3,2,2,4)/-0.13429e+02/, psa(3,3,2,2,4)/ 0.23600e+00/   
      data psa(1,4,2,2,4)/ 0.15000e+03/                                 
      data psa(2,4,2,2,4)/-0.17626e+02/, psa(3,4,2,2,4)/ 0.09400e+00/   
      data psa(1,5,2,2,4)/ 0.21500e+03/                                 
      data psa(2,5,2,2,4)/-0.22152e+02/, psa(3,5,2,2,4)/ 0.20400e+00/   
      data psa(1,6,2,2,4)/ 0.32000e+03/                                 
      data psa(2,6,2,2,4)/-0.29292e+02/, psa(3,6,2,2,4)/ 0.29200e+00/   
      data psa(1,7,2,2,4)/ 0.10000e+02/                                 
      data psa(2,7,2,2,4)/-0.21048e+01/, psa(3,7,2,2,4)/ 0.02660e+00/   
      data psa(1,8,2,2,4)/ 0.50000e+01/                                 
      data psa(2,8,2,2,4)/-0.90720e+00/, psa(3,8,2,2,4)/ 0.02820e+00/   
c        j = 2                                                          
c        1d2                                                            
      data psa(1,1,1,3,4)/ 0.25000e+02/                                 
      data psa(2,1,1,3,4)/ 0.76650e+00/, psa(3,1,1,3,4)/ 0.92200e-01/   
      data psa(1,2,1,3,4)/ 0.50000e+02/                                 
      data psa(2,2,1,3,4)/ 0.17153e+01/, psa(3,2,1,3,4)/ 0.01410e+00/   
      data psa(1,3,1,3,4)/ 0.10000e+03/                                 
      data psa(2,3,1,3,4)/ 0.36300e+01/, psa(3,3,1,3,4)/ 0.10000e+00/   
      data psa(1,4,1,3,4)/ 0.15000e+03/                                 
      data psa(2,4,1,3,4)/ 0.56910e+01/, psa(3,4,1,3,4)/ 0.07300e+00/   
      data psa(1,5,1,3,4)/ 0.21500e+03/                                 
      data psa(2,5,1,3,4)/ 0.77660e+01/, psa(3,5,1,3,4)/ 0.11500e+00/   
      data psa(1,6,1,3,4)/ 0.32000e+03/                                 
      data psa(2,6,1,3,4)/ 0.98710e+01/, psa(3,6,1,3,4)/ 0.11500e+00/   
      data psa(1,7,1,3,4)/ 0.10000e+02/                                 
      data psa(2,7,1,3,4)/ 0.17660e+00/, psa(3,7,1,3,4)/ 0.01080e+00/   
      data psa(1,8,1,3,4)/ 0.50000e+01/                                 
      data psa(2,8,1,3,4)/ 0.03950e+00/, psa(3,8,1,3,4)/ 0.00960e+00/   
c        3p2                                                            
      data psa(1,1,3,3,4)/ 0.25000e+02/                                 
      data psa(2,1,3,3,4)/ 0.23712e+01/, psa(3,1,3,3,4)/ 0.14250e+00/   
      data psa(1,2,3,3,4)/ 0.50000e+02/                                 
      data psa(2,2,3,3,4)/ 0.58834e+01/, psa(3,2,3,3,4)/ 0.06290e+00/   
      data psa(1,3,3,3,4)/ 0.10000e+03/                                 
      data psa(2,3,3,3,4)/ 0.10736e+02/, psa(3,3,3,3,4)/ 0.20700e+00/   
      data psa(1,4,3,3,4)/ 0.15000e+03/                                 
      data psa(2,4,3,3,4)/ 0.13943e+02/, psa(3,4,3,3,4)/ 0.06900e+00/   
      data psa(1,5,3,3,4)/ 0.21500e+03/                                 
      data psa(2,5,3,3,4)/ 0.15846e+02/, psa(3,5,3,3,4)/ 0.12000e+00/   
      data psa(1,6,3,3,4)/ 0.32000e+03/                                 
      data psa(2,6,3,3,4)/ 0.17508e+02/, psa(3,6,3,3,4)/ 0.13300e+00/   
      data psa(1,7,3,3,4)/ 0.10000e+02/                                 
      data psa(2,7,3,3,4)/ 0.66010e+00/, psa(3,7,3,3,4)/ 0.01720e+00/   
      data psa(1,8,3,3,4)/ 0.50000e+01/                                 
      data psa(2,8,3,3,4)/ 0.22720e+00/, psa(3,8,3,3,4)/ 0.01770e+00/   
c        3f2                                                            
      data psa(1,1,4,3,4)/ 0.15000e+03/                                 
      data psa(2,1,4,3,4)/ 0.10100e+01/, psa(3,1,4,3,4)/ 0.08900e+00/   
      data psa(1,2,4,3,4)/ 0.21500e+03/                                 
      data psa(2,2,4,3,4)/ 0.11830e+01/, psa(3,2,4,3,4)/ 0.16800e+00/   
      data psa(1,3,4,3,4)/ 0.32000e+03/                                 
      data psa(2,3,4,3,4)/ 0.13210e+01/, psa(3,3,4,3,4)/ 0.16000e+00/   
      data psa(1,4,4,3,4)/ 0.10000e+03/                                 
      data psa(2,4,4,3,4)/ 0.10830e+01/, psa(3,4,4,3,4)/ 0.11100e+00/   
c        e2                                                             
      data psa(1,1,5,3,4)/ 0.50000e+02/                                 
      data psa(2,1,5,3,4)/-0.17186e+01/, psa(3,1,5,3,4)/ 0.02210e+00/   
      data psa(1,2,5,3,4)/ 0.10000e+03/                                 
      data psa(2,2,5,3,4)/-0.25620e+01/, psa(3,2,5,3,4)/ 0.10500e+00/   
      data psa(1,3,5,3,4)/ 0.15000e+03/                                 
      data psa(2,3,5,3,4)/-0.28540e+01/, psa(3,3,5,3,4)/ 0.45000e-01/   
      data psa(1,4,5,3,4)/ 0.21500e+03/                                 
      data psa(2,4,5,3,4)/-0.25500e+01/, psa(3,4,5,3,4)/ 0.08500e+00/   
      data psa(1,5,5,3,4)/ 0.32000e+03/                                 
      data psa(2,5,5,3,4)/-0.24100e+01/, psa(3,5,5,3,4)/ 0.09400e+00/   
      data psa(1,6,5,3,4)/ 0.25000e+02/                                 
      data psa(2,6,5,3,4)/-0.95820e+00/, psa(3,6,5,3,4)/ 0.13030e+00/   
c        j = 3                                                          
c        3f3                                                            
      data psa(1,1,2,4,4)/ 0.15000e+03/                                 
      data psa(2,1,2,4,4)/-0.18810e+01/, psa(3,1,2,4,4)/ 0.11800e+00/   
      data psa(1,2,2,4,4)/ 0.21500e+03/                                 
      data psa(2,2,2,4,4)/-0.24550e+01/, psa(3,2,2,4,4)/ 0.10900e+00/   
      data psa(1,3,2,4,4)/ 0.32000e+03/                                 
      data psa(2,3,2,4,4)/-0.29190e+01/, psa(3,3,2,4,4)/ 0.11300e+00/   
c        j = 4                                                          
c        1g4                                                            
      data psa(1,1,1,5,4)/ 0.21500e+03/                                 
      data psa(2,1,1,5,4)/ 0.98600e+00/, psa(3,1,1,5,4)/ 0.04800e+00/   
      data psa(1,2,1,5,4)/ 0.32000e+03/                                 
      data psa(2,2,1,5,4)/ 0.15650e+01/, psa(3,2,1,5,4)/ 0.06500e+00/   
c        3f4                                                            
      data psa(1,1,3,5,4)/ 0.21500e+03/                                 
      data psa(2,1,3,5,4)/ 0.18590e+01/, psa(3,1,3,5,4)/ 0.10200e+00/   
      data psa(1,2,3,5,4)/ 0.32000e+03/                                 
      data psa(2,2,3,5,4)/ 0.30560e+01/, psa(3,2,3,5,4)/ 0.09100e+00/   
c        3h4                                                            
      data psa(1,1,4,5,4)/ 0.21500e+03/                                 
      data psa(2,1,4,5,4)/ 0.29200e+00/, psa(3,1,4,5,4)/ 0.06900e+00/   
      data psa(1,2,4,5,4)/ 0.32000e+03/                                 
      data psa(2,2,4,5,4)/ 0.67400e+00/, psa(3,2,4,5,4)/ 0.07900e+00/   
c        e4                                                             
      data psa(1,1,5,5,4)/ 0.21500e+03/                                 
      data psa(2,1,5,5,4)/-0.11620e+01/, psa(3,1,5,5,4)/ 0.05800e+00/   
      data psa(1,2,5,5,4)/ 0.32000e+03/                                 
      data psa(2,2,5,5,4)/-0.15570e+01/, psa(3,2,5,5,4)/ 0.05100e+00/   
c                                                                       
c                                                                       
c *****  phase-shift analysis  n93p,ips=5                               
c        --------------------------------                               
c         n93p: nijmegen p-p multi-energy analysis,                     
c               stoks et al., phys. rev. c48, 792 (1993);               
c                    1s0, 1. - 50. mev, corrected acc. to               
c                    stoks et al., phys. rev. c49, 2950 (1994),         
c                    and stoks, priv comm. 11/6/94.                     
c                                                                       
c                                                                       
c                                                                       
c        j = 0                                                          
c        1s0 pp                                                         
      data psa(1,1,1,1,5)/ 0.25000e+02/                                 
      data psa(2,1,1,1,5)/ 0.48664e+02/, psa(3,1,1,1,5)/ 0.03900e+00/   
      data psa(1,2,1,1,5)/ 0.50000e+02/                                 
      data psa(2,2,1,1,5)/ 0.38920e+02/, psa(3,2,1,1,5)/ 0.04900e+00/   
      data psa(1,3,1,1,5)/ 0.10000e+03/                                 
      data psa(2,3,1,1,5)/ 0.24980e+02/, psa(3,3,1,1,5)/ 0.08000e+00/   
      data psa(1,4,1,1,5)/ 0.15000e+03/                                 
      data psa(2,4,1,1,5)/ 0.14770e+02/, psa(3,4,1,1,5)/ 0.13000e+00/   
      data psa(1,5,1,1,5)/ 0.20000e+03/                                 
      data psa(2,5,1,1,5)/ 0.65700e+01/, psa(3,5,1,1,5)/ 0.16000e+00/   
      data psa(1,6,1,1,5)/ 0.30000e+03/                                 
      data psa(2,6,1,1,5)/-0.61400e+01/, psa(3,6,1,1,5)/ 0.25000e+00/   
      data psa(1,7,1,1,5)/ 0.10000e+02/                                 
      data psa(2,7,1,1,5)/ 0.55220e+02/, psa(3,7,1,1,5)/ 0.02500e+00/   
      data psa(1,8,1,1,5)/ 0.10000e+01/                                 
      data psa(2,8,1,1,5)/ 0.32770e+02/, psa(3,8,1,1,5)/ 0.00500e+00/   
      data psa(1,9,1,1,5)/ 0.50000e+01/                                 
      data psa(2,9,1,1,5)/ 0.54850e+02/, psa(3,9,1,1,5)/ 0.01700e+00/   
      data psa(1,10,1,1,5)/ 0.25000e+03/                                
      data psa(2,10,1,1,5)/-0.30000e+00/, psa(3,10,1,1,5)/ 0.18000e+00/ 
      data psa(1,11,1,1,5)/ 0.35000e+03/                                
      data psa(2,11,1,1,5)/-0.11110e+02/, psa(3,11,1,1,5)/ 0.46000e+00/ 
c        3p0                                                            
      data psa(1,1,3,1,5)/ 0.25000e+02/                                 
      data psa(2,1,3,1,5)/ 0.85750e+01/, psa(3,1,3,1,5)/ 0.05300e+00/   
      data psa(1,2,3,1,5)/ 0.50000e+02/                                 
      data psa(2,2,3,1,5)/ 0.11470e+02/, psa(3,2,3,1,5)/ 0.09000e+00/   
      data psa(1,3,3,1,5)/ 0.10000e+03/                                 
      data psa(2,3,3,1,5)/ 0.94500e+01/, psa(3,3,3,1,5)/ 0.11000e+00/   
      data psa(1,4,3,1,5)/ 0.15000e+03/                                 
      data psa(2,4,3,1,5)/ 0.47400e+01/, psa(3,4,3,1,5)/ 0.14000e+00/   
      data psa(1,5,3,1,5)/ 0.20000e+03/                                 
      data psa(2,5,3,1,5)/-0.37000e+00/, psa(3,5,3,1,5)/ 0.17000e+00/   
      data psa(1,6,3,1,5)/ 0.30000e+03/                                 
      data psa(2,6,3,1,5)/-0.10390e+02/, psa(3,6,3,1,5)/ 0.33000e+00/   
      data psa(1,7,3,1,5)/ 0.10000e+02/                                 
      data psa(2,7,3,1,5)/ 0.37290e+01/, psa(3,7,3,1,5)/ 0.01700e+00/   
      data psa(1,8,3,1,5)/ 0.50000e+01/                                 
      data psa(2,8,3,1,5)/ 0.15820e+01/, psa(3,8,3,1,5)/ 0.00600e+00/   
      data psa(1,9,3,1,5)/ 0.10000e+01/                                 
      data psa(2,9,3,1,5)/ 0.13400e+00/, psa(3,9,3,1,5)/ 0.00050e+00/   
      data psa(1,10,3,1,5)/ 0.25000e+03/                                
      data psa(2,10,3,1,5)/-0.54300e+01/, psa(3,10,3,1,5)/ 0.21000e+00/ 
      data psa(1,11,3,1,5)/ 0.35000e+03/                                
      data psa(2,11,3,1,5)/-0.15300e+02/, psa(3,11,3,1,5)/ 0.57000e+00/ 
c        j = 1                                                          
c        3p1                                                            
      data psa(1,1,2,2,5)/ 0.25000e+02/                                 
      data psa(2,1,2,2,5)/-0.49320e+01/, psa(3,1,2,2,5)/ 0.00800e+00/   
      data psa(1,2,2,2,5)/ 0.50000e+02/                                 
      data psa(2,2,2,2,5)/-0.83170e+01/, psa(3,2,2,2,5)/ 0.01700e+00/   
      data psa(1,3,2,2,5)/ 0.10000e+03/                                 
      data psa(2,3,2,2,5)/-0.13258e+02/, psa(3,3,2,2,5)/ 0.03200e+00/   
      data psa(1,4,2,2,5)/ 0.15000e+03/                                 
      data psa(2,4,2,2,5)/-0.17434e+02/, psa(3,4,2,2,5)/ 0.04500e+00/   
      data psa(1,5,2,2,5)/ 0.20000e+03/                                 
      data psa(2,5,2,2,5)/-0.21250e+02/, psa(3,5,2,2,5)/ 0.07000e+00/   
      data psa(1,6,2,2,5)/ 0.30000e+03/                                 
      data psa(2,6,2,2,5)/-0.27990e+02/, psa(3,6,2,2,5)/ 0.19000e+00/   
      data psa(1,7,2,2,5)/ 0.10000e+02/                                 
      data psa(2,7,2,2,5)/-0.20600e+01/, psa(3,7,2,2,5)/ 0.00200e+00/   
      data psa(1,8,2,2,5)/ 0.50000e+01/                                 
      data psa(2,8,2,2,5)/-0.90200e+00/, psa(3,8,2,2,5)/ 0.00100e+00/   
      data psa(1,9,2,2,5)/ 0.10000e+01/                                 
      data psa(2,9,2,2,5)/-0.08100e+00/, psa(3,9,2,2,5)/ 0.00050e+00/   
      data psa(1,10,2,2,5)/ 0.25000e+03/                                
      data psa(2,10,2,2,5)/-0.24770e+02/, psa(3,10,2,2,5)/ 0.12000e+00/ 
      data psa(1,11,2,2,5)/ 0.35000e+03/                                
      data psa(2,11,2,2,5)/-0.30890e+02/, psa(3,11,2,2,5)/ 0.27000e+00/ 
c        j = 2                                                          
c        1d2                                                            
      data psa(1,1,1,3,5)/ 0.25000e+02/                                 
      data psa(2,1,1,3,5)/ 0.69600e+00/, psa(3,1,1,3,5)/ 0.00100e+00/   
      data psa(1,2,1,3,5)/ 0.50000e+02/                                 
      data psa(2,2,1,3,5)/ 0.17110e+01/, psa(3,2,1,3,5)/ 0.00400e+00/   
      data psa(1,3,1,3,5)/ 0.10000e+03/                                 
      data psa(2,3,1,3,5)/ 0.37900e+01/, psa(3,3,1,3,5)/ 0.01800e+00/   
      data psa(1,4,1,3,5)/ 0.15000e+03/                                 
      data psa(2,4,1,3,5)/ 0.56060e+01/, psa(3,4,1,3,5)/ 0.03300e+00/   
      data psa(1,5,1,3,5)/ 0.20000e+03/                                 
      data psa(2,5,1,3,5)/ 0.70580e+01/, psa(3,5,1,3,5)/ 0.04500e+00/   
      data psa(1,6,1,3,5)/ 0.30000e+03/                                 
      data psa(2,6,1,3,5)/ 0.94200e+01/, psa(3,6,1,3,5)/ 0.08000e+00/   
      data psa(1,7,1,3,5)/ 0.10000e+02/                                 
      data psa(2,7,1,3,5)/ 0.16500e+00/, psa(3,7,1,3,5)/ 0.00050e+00/   
      data psa(1,8,1,3,5)/ 0.50000e+01/                                 
      data psa(2,8,1,3,5)/ 0.04300e+00/, psa(3,8,1,3,5)/ 0.00050e+00/   
      data psa(1,9,1,3,5)/ 0.10000e+01/                                 
      data psa(2,9,1,3,5)/ 0.00100e+00/, psa(3,9,1,3,5)/ 0.00050e+00/   
      data psa(1,10,1,3,5)/ 0.25000e+03/                                
      data psa(2,10,1,3,5)/ 0.82700e+01/, psa(3,10,1,3,5)/ 0.06000e+00/ 
      data psa(1,11,1,3,5)/ 0.35000e+03/                                
      data psa(2,11,1,3,5)/ 0.10690e+02/, psa(3,11,1,3,5)/ 0.14000e+00/ 
c        3p2                                                            
      data psa(1,1,3,3,5)/ 0.25000e+02/                                 
      data psa(2,1,3,3,5)/ 0.24910e+01/, psa(3,1,3,3,5)/ 0.00800e+00/   
      data psa(1,2,3,3,5)/ 0.50000e+02/                                 
      data psa(2,2,3,3,5)/ 0.58550e+01/, psa(3,2,3,3,5)/ 0.01600e+00/   
      data psa(1,3,3,3,5)/ 0.10000e+03/                                 
      data psa(2,3,3,3,5)/ 0.11013e+02/, psa(3,3,3,3,5)/ 0.02500e+00/   
      data psa(1,4,3,3,5)/ 0.15000e+03/                                 
      data psa(2,4,3,3,5)/ 0.13982e+02/, psa(3,4,3,3,5)/ 0.03900e+00/   
      data psa(1,5,3,3,5)/ 0.20000e+03/                                 
      data psa(2,5,3,3,5)/ 0.15630e+02/, psa(3,5,3,3,5)/ 0.05200e+00/   
      data psa(1,6,3,3,5)/ 0.30000e+03/                                 
      data psa(2,6,3,3,5)/ 0.17170e+02/, psa(3,6,3,3,5)/ 0.10000e+00/   
      data psa(1,7,3,3,5)/ 0.10000e+02/                                 
      data psa(2,7,3,3,5)/ 0.65110e+00/, psa(3,7,3,3,5)/ 0.00200e+00/   
      data psa(1,8,3,3,5)/ 0.50000e+01/                                 
      data psa(2,8,3,3,5)/ 0.21400e+00/, psa(3,8,3,3,5)/ 0.00100e+00/   
      data psa(1,9,3,3,5)/ 0.10000e+01/                                 
      data psa(2,9,3,3,5)/ 0.01400e+00/, psa(3,9,3,3,5)/ 0.00050e+00/   
      data psa(1,10,3,3,5)/ 0.25000e+03/                                
      data psa(2,10,3,3,5)/ 0.16590e+02/, psa(3,10,3,3,5)/ 0.07000e+00/ 
      data psa(1,11,3,3,5)/ 0.35000e+03/                                
      data psa(2,11,3,3,5)/ 0.17540e+02/, psa(3,11,3,3,5)/ 0.15000e+00/ 
c        3f2                                                            
      data psa(1,1,4,3,5)/ 0.25000e+02/                                 
      data psa(2,1,4,3,5)/ 0.10500e+00/, psa(3,1,4,3,5)/ 0.00050e+00/   
      data psa(1,2,4,3,5)/ 0.50000e+02/                                 
      data psa(2,2,4,3,5)/ 0.33800e+00/, psa(3,2,4,3,5)/ 0.00050e+00/   
      data psa(1,3,4,3,5)/ 0.10000e+03/                                 
      data psa(2,3,4,3,5)/ 0.81700e+00/, psa(3,3,4,3,5)/ 0.00400e+00/   
      data psa(1,4,4,3,5)/ 0.15000e+03/                                 
      data psa(2,4,4,3,5)/ 0.11970e+01/, psa(3,4,4,3,5)/ 0.01400e+00/   
      data psa(1,5,4,3,5)/ 0.20000e+03/                                 
      data psa(2,5,4,3,5)/ 0.14240e+01/, psa(3,5,4,3,5)/ 0.03400e+00/   
      data psa(1,6,4,3,5)/ 0.30000e+03/                                 
      data psa(2,6,4,3,5)/ 0.13400e+01/, psa(3,6,4,3,5)/ 0.11000e+00/   
      data psa(1,7,4,3,5)/ 0.10000e+02/                                 
      data psa(2,7,4,3,5)/ 0.01300e+00/, psa(3,7,4,3,5)/ 0.00050e+00/   
      data psa(1,8,4,3,5)/ 0.50000e+01/                                 
      data psa(2,8,4,3,5)/ 0.00200e+00/, psa(3,8,4,3,5)/ 0.00050e+00/   
      data psa(1,9,4,3,5)/ 0.10000e+01/                                 
      data psa(2,9,4,3,5)/ 0.00000e+00/, psa(3,9,4,3,5)/ 0.00050e+00/   
      data psa(1,10,4,3,5)/ 0.25000e+03/                                
      data psa(2,10,4,3,5)/ 0.14700e+01/, psa(3,10,4,3,5)/ 0.06000e+00/ 
      data psa(1,11,4,3,5)/ 0.35000e+03/                                
      data psa(2,11,4,3,5)/ 0.10400e+01/, psa(3,11,4,3,5)/ 0.16000e+00/ 
c        e2                                                             
      data psa(1,1,5,3,5)/ 0.25000e+02/                                 
      data psa(2,1,5,3,5)/-0.81000e+00/, psa(3,1,5,3,5)/ 0.00100e+00/   
      data psa(1,2,5,3,5)/ 0.50000e+02/                                 
      data psa(2,2,5,3,5)/-0.17120e+01/, psa(3,2,5,3,5)/ 0.00400e+00/   
      data psa(1,3,5,3,5)/ 0.10000e+03/                                 
      data psa(2,3,5,3,5)/-0.26590e+01/, psa(3,3,5,3,5)/ 0.01700e+00/   
      data psa(1,4,5,3,5)/ 0.15000e+03/                                 
      data psa(2,4,5,3,5)/-0.28730e+01/, psa(3,4,5,3,5)/ 0.02900e+00/   
      data psa(1,5,5,3,5)/ 0.20000e+03/                                 
      data psa(2,5,5,3,5)/-0.27590e+01/, psa(3,5,5,3,5)/ 0.03700e+00/   
      data psa(1,6,5,3,5)/ 0.30000e+03/                                 
      data psa(2,6,5,3,5)/-0.23400e+01/, psa(3,6,5,3,5)/ 0.09000e+00/   
      data psa(1,7,5,3,5)/ 0.10000e+02/                                 
      data psa(2,7,5,3,5)/-0.20000e+00/, psa(3,7,5,3,5)/ 0.00050e+00/   
      data psa(1,8,5,3,5)/ 0.50000e+01/                                 
      data psa(2,8,5,3,5)/-0.05200e+00/, psa(3,8,5,3,5)/ 0.00050e+00/   
      data psa(1,9,5,3,5)/ 0.10000e+01/                                 
      data psa(2,9,5,3,5)/-0.00100e+00/, psa(3,9,5,3,5)/ 0.00050e+00/   
      data psa(1,10,5,3,5)/ 0.25000e+03/                                
      data psa(2,10,5,3,5)/-0.25420e+01/, psa(3,10,5,3,5)/ 0.04600e+00/ 
      data psa(1,11,5,3,5)/ 0.35000e+03/                                
      data psa(2,11,5,3,5)/-0.22100e+01/, psa(3,11,5,3,5)/ 0.11000e+00/ 
c        j = 3                                                          
c        3f3                                                            
      data psa(1,1,2,4,5)/ 0.25000e+02/                                 
      data psa(2,1,2,4,5)/-0.23100e+00/, psa(3,1,2,4,5)/ 0.00050e+00/   
      data psa(1,2,2,4,5)/ 0.50000e+02/                                 
      data psa(2,2,2,4,5)/-0.69000e+00/, psa(3,2,2,4,5)/ 0.00050e+00/   
      data psa(1,3,2,4,5)/ 0.10000e+03/                                 
      data psa(2,3,2,4,5)/-0.15170e+01/, psa(3,3,2,4,5)/ 0.00300e+00/   
      data psa(1,4,2,4,5)/ 0.15000e+03/                                 
      data psa(2,4,2,4,5)/-0.21000e+01/, psa(3,4,2,4,5)/ 0.01000e+00/   
      data psa(1,5,2,4,5)/ 0.20000e+03/                                 
      data psa(2,5,2,4,5)/-0.24870e+01/, psa(3,5,2,4,5)/ 0.02500e+00/   
      data psa(1,6,2,4,5)/ 0.30000e+03/                                 
      data psa(2,6,2,4,5)/-0.28400e+01/, psa(3,6,2,4,5)/ 0.11000e+00/   
      data psa(1,7,2,4,5)/ 0.10000e+02/                                 
      data psa(2,7,2,4,5)/-0.03200e+00/, psa(3,7,2,4,5)/ 0.00050e+00/   
      data psa(1,8,2,4,5)/ 0.50000e+01/                                 
      data psa(2,8,2,4,5)/-0.00500e+00/, psa(3,8,2,4,5)/ 0.00050e+00/   
      data psa(1,9,2,4,5)/ 0.10000e+01/                                 
      data psa(2,9,2,4,5)/-0.00000e+00/, psa(3,9,2,4,5)/ 0.00050e+00/   
      data psa(1,10,2,4,5)/ 0.25000e+03/                                
      data psa(2,10,2,4,5)/-0.27240e+01/, psa(3,10,2,4,5)/ 0.04900e+00/ 
      data psa(1,11,2,4,5)/ 0.35000e+03/                                
      data psa(2,11,2,4,5)/-0.28700e+01/, psa(3,11,2,4,5)/ 0.13000e+00/ 
c        j = 4                                                          
c        1g4                                                            
      data psa(1,1,1,5,5)/ 0.25000e+02/                                 
      data psa(2,1,1,5,5)/ 0.04000e+00/, psa(3,1,1,5,5)/ 0.00050e+00/   
      data psa(1,2,1,5,5)/ 0.50000e+02/                                 
      data psa(2,2,1,5,5)/ 0.15200e+00/, psa(3,2,1,5,5)/ 0.00050e+00/   
      data psa(1,3,1,5,5)/ 0.10000e+03/                                 
      data psa(2,3,1,5,5)/ 0.41800e+00/, psa(3,3,1,5,5)/ 0.00100e+00/   
      data psa(1,4,1,5,5)/ 0.15000e+03/                                 
      data psa(2,4,1,5,5)/ 0.70000e+00/, psa(3,4,1,5,5)/ 0.00300e+00/   
      data psa(1,5,1,5,5)/ 0.20000e+03/                                 
      data psa(2,5,1,5,5)/ 0.99300e+00/, psa(3,5,1,5,5)/ 0.01000e+00/   
      data psa(1,6,1,5,5)/ 0.30000e+03/                                 
      data psa(2,6,1,5,5)/ 0.15030e+01/, psa(3,6,1,5,5)/ 0.04800e+00/   
      data psa(1,7,1,5,5)/ 0.10000e+02/                                 
      data psa(2,7,1,5,5)/ 0.00300e+00/, psa(3,7,1,5,5)/ 0.00050e+00/   
      data psa(1,8,1,5,5)/ 0.50000e+01/                                 
      data psa(2,8,1,5,5)/ 0.00000e+00/, psa(3,8,1,5,5)/ 0.00050e+00/   
      data psa(1,9,1,5,5)/ 0.10000e+01/                                 
      data psa(2,9,1,5,5)/ 0.00000e+00/, psa(3,9,1,5,5)/ 0.00050e+00/   
      data psa(1,10,1,5,5)/ 0.25000e+03/                                
      data psa(2,10,1,5,5)/ 0.12720e+01/, psa(3,10,1,5,5)/ 0.02400e+00/ 
      data psa(1,11,1,5,5)/ 0.35000e+03/                                
      data psa(2,11,1,5,5)/ 0.16400e+01/, psa(3,11,1,5,5)/ 0.08000e+00/ 
c        3f4                                                            
      data psa(1,1,3,5,5)/ 0.25000e+02/                                 
      data psa(2,1,3,5,5)/ 0.02000e+00/, psa(3,1,3,5,5)/ 0.00050e+00/   
      data psa(1,2,3,5,5)/ 0.50000e+02/                                 
      data psa(2,2,3,5,5)/ 0.10800e+00/, psa(3,2,3,5,5)/ 0.00100e+00/   
      data psa(1,3,3,5,5)/ 0.10000e+03/                                 
      data psa(2,3,3,5,5)/ 0.47800e+00/, psa(3,3,3,5,5)/ 0.00700e+00/   
      data psa(1,4,3,5,5)/ 0.15000e+03/                                 
      data psa(2,4,3,5,5)/ 0.10320e+01/, psa(3,4,3,5,5)/ 0.02200e+00/   
      data psa(1,5,3,5,5)/ 0.20000e+03/                                 
      data psa(2,5,3,5,5)/ 0.16780e+01/, psa(3,5,3,5,5)/ 0.03900e+00/   
      data psa(1,6,3,5,5)/ 0.30000e+03/                                 
      data psa(2,6,3,5,5)/ 0.28900e+01/, psa(3,6,3,5,5)/ 0.06000e+00/   
      data psa(1,7,3,5,5)/ 0.10000e+02/                                 
      data psa(2,7,3,5,5)/ 0.00100e+00/, psa(3,7,3,5,5)/ 0.00050e+00/   
      data psa(1,8,3,5,5)/ 0.50000e+01/                                 
      data psa(2,8,3,5,5)/ 0.00000e+00/, psa(3,8,3,5,5)/ 0.00050e+00/   
      data psa(1,9,3,5,5)/ 0.10000e+01/                                 
      data psa(2,9,3,5,5)/ 0.00000e+00/, psa(3,9,3,5,5)/ 0.00050e+00/   
      data psa(1,10,3,5,5)/ 0.25000e+03/                                
      data psa(2,10,3,5,5)/ 0.23250e+01/, psa(3,10,3,5,5)/ 0.05100e+00/ 
      data psa(1,11,3,5,5)/ 0.35000e+03/                                
      data psa(2,11,3,5,5)/ 0.33000e+01/, psa(3,11,3,5,5)/ 0.11000e+00/ 
c        3h4                                                            
      data psa(1,1,4,5,5)/ 0.25000e+02/                                 
      data psa(2,1,4,5,5)/ 0.00400e+00/, psa(3,1,4,5,5)/ 0.00050e+00/   
      data psa(1,2,4,5,5)/ 0.50000e+02/                                 
      data psa(2,2,4,5,5)/ 0.02600e+00/, psa(3,2,4,5,5)/ 0.00050e+00/   
      data psa(1,3,4,5,5)/ 0.10000e+03/                                 
      data psa(2,3,4,5,5)/ 0.10800e+00/, psa(3,3,4,5,5)/ 0.00050e+00/   
      data psa(1,4,4,5,5)/ 0.15000e+03/                                 
      data psa(2,4,4,5,5)/ 0.21100e+00/, psa(3,4,4,5,5)/ 0.00050e+00/   
      data psa(1,5,4,5,5)/ 0.20000e+03/                                 
      data psa(2,5,4,5,5)/ 0.32100e+00/, psa(3,5,4,5,5)/ 0.00050e+00/   
      data psa(1,6,4,5,5)/ 0.30000e+03/                                 
      data psa(2,6,4,5,5)/ 0.52600e+00/, psa(3,6,4,5,5)/ 0.00500e+00/   
      data psa(1,7,4,5,5)/ 0.10000e+02/                                 
      data psa(2,7,4,5,5)/ 0.00000e+00/, psa(3,7,4,5,5)/ 0.00050e+00/   
      data psa(1,8,4,5,5)/ 0.50000e+01/                                 
      data psa(2,8,4,5,5)/ 0.00000e+00/, psa(3,8,4,5,5)/ 0.00050e+00/   
      data psa(1,9,4,5,5)/ 0.10000e+01/                                 
      data psa(2,9,4,5,5)/ 0.00000e+00/, psa(3,9,4,5,5)/ 0.00050e+00/   
      data psa(1,10,4,5,5)/ 0.25000e+03/                                
      data psa(2,10,4,5,5)/ 0.42800e+00/, psa(3,10,4,5,5)/ 0.00050e+00/ 
      data psa(1,11,4,5,5)/ 0.35000e+03/                                
      data psa(2,11,4,5,5)/ 0.60800e+00/, psa(3,11,4,5,5)/ 0.00050e+00/ 
c        e4                                                             
      data psa(1,1,5,5,5)/ 0.25000e+02/                                 
      data psa(2,1,5,5,5)/-0.04900e+00/, psa(3,1,5,5,5)/ 0.00050e+00/   
      data psa(1,2,5,5,5)/ 0.50000e+02/                                 
      data psa(2,2,5,5,5)/-0.19500e+00/, psa(3,2,5,5,5)/ 0.00050e+00/   
      data psa(1,3,5,5,5)/ 0.10000e+03/                                 
      data psa(2,3,5,5,5)/-0.53900e+00/, psa(3,3,5,5,5)/ 0.00050e+00/   
      data psa(1,4,5,5,5)/ 0.15000e+03/                                 
      data psa(2,4,5,5,5)/-0.84900e+00/, psa(3,4,5,5,5)/ 0.00050e+00/   
      data psa(1,5,5,5,5)/ 0.20000e+03/                                 
      data psa(2,5,5,5,5)/-0.11080e+01/, psa(3,5,5,5,5)/ 0.00050e+00/   
      data psa(1,6,5,5,5)/ 0.30000e+03/                                 
      data psa(2,6,5,5,5)/-0.14700e+01/, psa(3,6,5,5,5)/ 0.00050e+00/   
      data psa(1,7,5,5,5)/ 0.10000e+02/                                 
      data psa(2,7,5,5,5)/-0.00400e+00/, psa(3,7,5,5,5)/ 0.00050e+00/   
      data psa(1,8,5,5,5)/ 0.50000e+01/                                 
      data psa(2,8,5,5,5)/-0.00000e+00/, psa(3,8,5,5,5)/ 0.00050e+00/   
      data psa(1,9,5,5,5)/ 0.10000e+01/                                 
      data psa(2,9,5,5,5)/-0.00000e+00/, psa(3,9,5,5,5)/ 0.00050e+00/   
      data psa(1,10,5,5,5)/ 0.25000e+03/                                
      data psa(2,10,5,5,5)/-0.13140e+01/, psa(3,10,5,5,5)/ 0.00050e+00/ 
      data psa(1,11,5,5,5)/ 0.35000e+03/                                
      data psa(2,11,5,5,5)/-0.15880e+01/, psa(3,11,5,5,5)/ 0.00100e+00/ 
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c *****  phase-shift analysis  n93n, ips=6                              
c        --------------------------------                               
c         n93n: nijmegen n-p multi-energy analysis,                     
c               stoks et al., phys.rev. c 48, 792 (1993)                
c                                                                       
c                                                                       
c                                                                       
c        j = 0                                                          
c        1s0 np                                                         
      data psa(1,1,1,1,6)/ 0.25000e+02/                                 
      data psa(2,1,1,1,6)/ 0.50900e+02/, psa(3,1,1,1,6)/ 0.19000e+00/   
      data psa(1,2,1,1,6)/ 0.50000e+02/                                 
      data psa(2,2,1,1,6)/ 0.40540e+02/, psa(3,2,1,1,6)/ 0.28000e+00/   
      data psa(1,3,1,1,6)/ 0.10000e+03/                                 
      data psa(2,3,1,1,6)/ 0.26780e+02/, psa(3,3,1,1,6)/ 0.38000e+00/   
      data psa(1,4,1,1,6)/ 0.15000e+03/                                 
      data psa(2,4,1,1,6)/ 0.16940e+02/, psa(3,4,1,1,6)/ 0.41000e+00/   
      data psa(1,5,1,1,6)/ 0.20000e+03/                                 
      data psa(2,5,1,1,6)/ 0.89400e+01/, psa(3,5,1,1,6)/ 0.39000e+00/   
      data psa(1,6,1,1,6)/ 0.30000e+03/                                 
      data psa(2,6,1,1,6)/-0.44600e+01/, psa(3,6,1,1,6)/ 0.43000e+00/   
      data psa(1,7,1,1,6)/ 0.10000e+02/                                 
      data psa(2,7,1,1,6)/ 0.59960e+02/, psa(3,7,1,1,6)/ 0.11000e+00/   
      data psa(1,8,1,1,6)/ 0.10000e+01/                                 
      data psa(2,8,1,1,6)/ 0.62068e+02/, psa(3,8,1,1,6)/ 0.03000e+00/   
      data psa(1,9,1,1,6)/ 0.50000e+01/                                 
      data psa(2,9,1,1,6)/ 0.63630e+02/, psa(3,9,1,1,6)/ 0.08000e+00/   
      data psa(1,10,1,1,6)/ 0.25000e+03/                                
      data psa(2,10,1,1,6)/ 0.19600e+01/, psa(3,10,1,1,6)/ 0.37000e+00/ 
      data psa(1,11,1,1,6)/ 0.35000e+03/                                
      data psa(2,11,1,1,6)/-0.10590e+02/, psa(3,11,1,1,6)/ 0.62000e+00/ 
c        3p0                                                            
      data psa(1,1,3,1,6)/ 0.25000e+02/                                 
      data psa(2,1,3,1,6)/ 0.81300e+01/, psa(3,1,3,1,6)/ 0.05300e+00/   
      data psa(1,2,3,1,6)/ 0.50000e+02/                                 
      data psa(2,2,3,1,6)/ 0.10700e+02/, psa(3,2,3,1,6)/ 0.09000e+00/   
      data psa(1,3,3,1,6)/ 0.10000e+03/                                 
      data psa(2,3,3,1,6)/ 0.84600e+01/, psa(3,3,3,1,6)/ 0.11000e+00/   
      data psa(1,4,3,1,6)/ 0.15000e+03/                                 
      data psa(2,4,3,1,6)/ 0.36900e+01/, psa(3,4,3,1,6)/ 0.14000e+00/   
      data psa(1,5,3,1,6)/ 0.20000e+03/                                 
      data psa(2,5,3,1,6)/-0.14400e+01/, psa(3,5,3,1,6)/ 0.17000e+00/   
      data psa(1,6,3,1,6)/ 0.30000e+03/                                 
      data psa(2,6,3,1,6)/-0.11470e+02/, psa(3,6,3,1,6)/ 0.33000e+00/   
      data psa(1,7,3,1,6)/ 0.10000e+02/                                 
      data psa(2,7,3,1,6)/ 0.36500e+01/, psa(3,7,3,1,6)/ 0.01700e+00/   
      data psa(1,8,3,1,6)/ 0.50000e+01/                                 
      data psa(2,8,3,1,6)/ 0.16300e+01/, psa(3,8,3,1,6)/ 0.00600e+00/   
      data psa(1,9,3,1,6)/ 0.10000e+01/                                 
      data psa(2,9,3,1,6)/ 0.18000e+00/, psa(3,9,3,1,6)/ 0.0050e+00/    
      data psa(1,10,3,1,6)/ 0.25000e+03/                                
      data psa(2,10,3,1,6)/-0.65100e+01/, psa(3,10,3,1,6)/ 0.21000e+00/ 
      data psa(1,11,3,1,6)/ 0.35000e+03/                                
      data psa(2,11,3,1,6)/-0.16390e+02/, psa(3,11,3,1,6)/ 0.57000e+00/ 
c        j = 1                                                          
c        1p1                                                            
      data psa(1,1,1,2,6)/ 0.25000e+02/                                 
      data psa(2,1,1,2,6)/-0.63110e+01/, psa(3,1,1,2,6)/ 0.03900e+00/   
      data psa(1,2,1,2,6)/ 0.50000e+02/                                 
      data psa(2,2,1,2,6)/-0.96700e+01/, psa(3,2,1,2,6)/ 0.08000e+00/   
      data psa(1,3,1,2,6)/ 0.10000e+03/                                 
      data psa(2,3,1,2,6)/-0.14520e+02/, psa(3,3,1,2,6)/ 0.14000e+00/   
      data psa(1,4,1,2,6)/ 0.15000e+03/                                 
      data psa(2,4,1,2,6)/-0.18650e+02/, psa(3,4,1,2,6)/ 0.16000e+00/   
      data psa(1,5,1,2,6)/ 0.20000e+03/                                 
      data psa(2,5,1,2,6)/-0.22180e+02/, psa(3,5,1,2,6)/ 0.18000e+00/   
      data psa(1,6,1,2,6)/ 0.30000e+03/                                 
      data psa(2,6,1,2,6)/-0.27580e+02/, psa(3,6,1,2,6)/ 0.22000e+00/   
      data psa(1,7,1,2,6)/ 0.10000e+02/                                 
      data psa(2,7,1,2,6)/-0.30390e+01/, psa(3,7,1,2,6)/ 0.01200e+00/   
      data psa(1,8,1,2,6)/ 0.50000e+01/                                 
      data psa(2,8,1,2,6)/-0.14870e+01/, psa(3,8,1,2,6)/ 0.00400e+00/   
      data psa(1,9,1,2,6)/ 0.10000e+01/                                 
      data psa(2,9,1,2,6)/-0.18700e+00/, psa(3,9,1,2,6)/ 0.00050e+00/   
      data psa(1,10,1,2,6)/ 0.25000e+03/                                
      data psa(2,10,1,2,6)/-0.25130e+02/, psa(3,10,1,2,6)/ 0.20000e+00/ 
      data psa(1,11,1,2,6)/ 0.35000e+03/                                
      data psa(2,11,1,2,6)/-0.29660e+02/, psa(3,11,1,2,6)/ 0.33000e+00/ 
c        3p1                                                            
      data psa(1,1,2,2,6)/ 0.25000e+02/                                 
      data psa(2,1,2,2,6)/-0.48800e+01/, psa(3,1,2,2,6)/ 0.00800e+00/   
      data psa(1,2,2,2,6)/ 0.50000e+02/                                 
      data psa(2,2,2,2,6)/-0.82500e+01/, psa(3,2,2,2,6)/ 0.01700e+00/   
      data psa(1,3,2,2,6)/ 0.10000e+03/                                 
      data psa(2,3,2,2,6)/-0.13240e+02/, psa(3,3,2,2,6)/ 0.03200e+00/   
      data psa(1,4,2,2,6)/ 0.15000e+03/                                 
      data psa(2,4,2,2,6)/-0.17460e+02/, psa(3,4,2,2,6)/ 0.04500e+00/   
      data psa(1,5,2,2,6)/ 0.20000e+03/                                 
      data psa(2,5,2,2,6)/-0.21300e+02/, psa(3,5,2,2,6)/ 0.07000e+00/   
      data psa(1,6,2,2,6)/ 0.30000e+03/                                 
      data psa(2,6,2,2,6)/-0.28070e+02/, psa(3,6,2,2,6)/ 0.19000e+00/   
      data psa(1,7,2,2,6)/ 0.10000e+02/                                 
      data psa(2,7,2,2,6)/-0.20600e+01/, psa(3,7,2,2,6)/ 0.00500e+00/   
      data psa(1,8,2,2,6)/ 0.50000e+01/                                 
      data psa(2,8,2,2,6)/-0.94000e+00/, psa(3,8,2,2,6)/ 0.00500e+00/   
      data psa(1,9,2,2,6)/ 0.10000e+01/                                 
      data psa(2,9,2,2,6)/-0.11000e+00/, psa(3,9,2,2,6)/ 0.0050e+00/    
      data psa(1,10,2,2,6)/ 0.25000e+03/                                
      data psa(2,10,2,2,6)/-0.24840e+02/, psa(3,10,2,2,6)/ 0.12000e+00/ 
      data psa(1,11,2,2,6)/ 0.35000e+03/                                
      data psa(2,11,2,2,6)/-0.30970e+02/, psa(3,11,2,2,6)/ 0.27000e+00/ 
c        3s1                                                            
      data psa(1,1,3,2,6)/ 0.25000e+02/                                 
      data psa(2,1,3,2,6)/ 0.80630e+02/, psa(3,1,3,2,6)/ 0.07000e+00/   
      data psa(1,2,3,2,6)/ 0.50000e+02/                                 
      data psa(2,2,3,2,6)/ 0.62770e+02/, psa(3,2,3,2,6)/ 0.10000e+00/   
      data psa(1,3,3,2,6)/ 0.10000e+03/                                 
      data psa(2,3,3,2,6)/ 0.43230e+02/, psa(3,3,3,2,6)/ 0.14000e+00/   
      data psa(1,4,3,2,6)/ 0.15000e+03/                                 
      data psa(2,4,3,2,6)/ 0.30720e+02/, psa(3,4,3,2,6)/ 0.14000e+00/   
      data psa(1,5,3,2,6)/ 0.20000e+03/                                 
      data psa(2,5,3,2,6)/ 0.21220e+02/, psa(3,5,3,2,6)/ 0.15000e+00/   
      data psa(1,6,3,2,6)/ 0.30000e+03/                                 
      data psa(2,6,3,2,6)/ 0.66000e+01/, psa(3,6,3,2,6)/ 0.23000e+00/   
      data psa(1,7,3,2,6)/ 0.10000e+02/                                 
      data psa(2,7,3,2,6)/ 1.02611e+02/, psa(3,7,3,2,6)/ 0.03500e+00/   
      data psa(1,8,3,2,6)/ 0.50000e+01/                                 
      data psa(2,8,3,2,6)/ 1.18178e+02/, psa(3,8,3,2,6)/ 0.02100e+00/   
      data psa(1,9,3,2,6)/ 0.10000e+01/                                 
      data psa(2,9,3,2,6)/ 1.47747e+02/, psa(3,9,3,2,6)/ 0.01000e+00/   
      data psa(1,10,3,2,6)/ 0.25000e+03/                                
      data psa(2,10,3,2,6)/ 0.13390e+02/, psa(3,10,3,2,6)/ 0.17000e+00/ 
      data psa(1,11,3,2,6)/ 0.35000e+03/                                
      data psa(2,11,3,2,6)/ 0.50200e+00/, psa(3,11,3,2,6)/ 0.32000e+00/ 
c        3d1                                                            
      data psa(1,1,4,2,6)/ 0.25000e+02/                                 
      data psa(2,1,4,2,6)/-0.27990e+01/, psa(3,1,4,2,6)/ 0.00600e+00/   
      data psa(1,2,4,2,6)/ 0.50000e+02/                                 
      data psa(2,2,4,2,6)/-0.64330e+01/, psa(3,2,4,2,6)/ 0.01700e+00/   
      data psa(1,3,4,2,6)/ 0.10000e+03/                                 
      data psa(2,3,4,2,6)/-0.12230e+02/, psa(3,3,4,2,6)/ 0.05000e+00/   
      data psa(1,4,4,2,6)/ 0.15000e+03/                                 
      data psa(2,4,4,2,6)/-0.16480e+02/, psa(3,4,4,2,6)/ 0.09000e+00/   
      data psa(1,5,4,2,6)/ 0.20000e+03/                                 
      data psa(2,5,4,2,6)/-0.19710e+02/, psa(3,5,4,2,6)/ 0.14000e+00/   
      data psa(1,6,4,2,6)/ 0.30000e+03/                                 
      data psa(2,6,4,2,6)/-0.24140e+02/, psa(3,6,4,2,6)/ 0.18000e+00/   
      data psa(1,7,4,2,6)/ 0.10000e+02/                                 
      data psa(2,7,4,2,6)/-0.67700e+00/, psa(3,7,4,2,6)/ 0.00100e+00/   
      data psa(1,8,4,2,6)/ 0.50000e+01/                                 
      data psa(2,8,4,2,6)/-0.18300e+00/, psa(3,8,4,2,6)/ 0.00050e+00/   
      data psa(1,9,4,2,6)/ 0.10000e+01/                                 
      data psa(2,9,4,2,6)/-0.00500e+00/, psa(3,9,4,2,6)/ 0.00050e+00/   
      data psa(1,10,4,2,6)/ 0.25000e+03/                                
      data psa(2,10,4,2,6)/-0.22210e+02/, psa(3,10,4,2,6)/ 0.16000e+00/ 
      data psa(1,11,4,2,6)/ 0.35000e+03/                                
      data psa(2,11,4,2,6)/-0.25570e+02/, psa(3,11,4,2,6)/ 0.30000e+00/ 
c        e1                                                             
      data psa(1,1,5,2,6)/ 0.25000e+02/                                 
      data psa(2,1,5,2,6)/ 0.17930e+01/, psa(3,1,5,2,6)/ 0.02500e+00/   
      data psa(1,2,5,2,6)/ 0.50000e+02/                                 
      data psa(2,2,5,2,6)/ 0.21090e+01/, psa(3,2,5,2,6)/ 0.04800e+00/   
      data psa(1,3,5,2,6)/ 0.10000e+03/                                 
      data psa(2,3,5,2,6)/ 0.24200e+01/, psa(3,3,5,2,6)/ 0.09000e+00/   
      data psa(1,4,5,2,6)/ 0.15000e+03/                                 
      data psa(2,4,5,2,6)/ 0.27500e+01/, psa(3,4,5,2,6)/ 0.11000e+00/   
      data psa(1,5,5,2,6)/ 0.20000e+03/                                 
      data psa(2,5,5,2,6)/ 0.31300e+01/, psa(3,5,5,2,6)/ 0.12000e+00/   
      data psa(1,6,5,2,6)/ 0.30000e+03/                                 
      data psa(2,6,5,2,6)/ 0.40300e+01/, psa(3,6,5,2,6)/ 0.17000e+00/   
      data psa(1,7,5,2,6)/ 0.10000e+02/                                 
      data psa(2,7,5,2,6)/ 0.11590e+01/, psa(3,7,5,2,6)/ 0.01000e+00/   
      data psa(1,8,5,2,6)/ 0.50000e+01/                                 
      data psa(2,8,5,2,6)/ 0.67200e+00/, psa(3,8,5,2,6)/ 0.00400e+00/   
      data psa(1,9,5,2,6)/ 0.10000e+01/                                 
      data psa(2,9,5,2,6)/ 0.10500e+00/, psa(3,9,5,2,6)/ 0.00100e+00/   
      data psa(1,10,5,2,6)/ 0.25000e+03/                                
      data psa(2,10,5,2,6)/ 0.35600e+01/, psa(3,10,5,2,6)/ 0.13000e+00/ 
      data psa(1,11,5,2,6)/ 0.35000e+03/                                
      data psa(2,11,5,2,6)/ 0.45700e+01/, psa(3,11,5,2,6)/ 0.25000e+00/ 
c        j = 2                                                          
c        1d2                                                            
      data psa(1,1,1,3,6)/ 0.25000e+02/                                 
      data psa(2,1,1,3,6)/ 0.68000e+00/, psa(3,1,1,3,6)/ 0.00500e+00/   
      data psa(1,2,1,3,6)/ 0.50000e+02/                                 
      data psa(2,2,1,3,6)/ 0.17300e+01/, psa(3,2,1,3,6)/ 0.00500e+00/   
      data psa(1,3,1,3,6)/ 0.10000e+03/                                 
      data psa(2,3,1,3,6)/ 0.39000e+01/, psa(3,3,1,3,6)/ 0.01800e+00/   
      data psa(1,4,1,3,6)/ 0.15000e+03/                                 
      data psa(2,4,1,3,6)/ 0.57900e+01/, psa(3,4,1,3,6)/ 0.03300e+00/   
      data psa(1,5,1,3,6)/ 0.20000e+03/                                 
      data psa(2,5,1,3,6)/ 0.72900e+01/, psa(3,5,1,3,6)/ 0.04500e+00/   
      data psa(1,6,1,3,6)/ 0.30000e+03/                                 
      data psa(2,6,1,3,6)/ 0.96900e+01/, psa(3,6,1,3,6)/ 0.08000e+00/   
      data psa(1,7,1,3,6)/ 0.10000e+02/                                 
      data psa(2,7,1,3,6)/ 0.16000e+00/, psa(3,7,1,3,6)/ 0.0050e+00/    
      data psa(1,8,1,3,6)/ 0.50000e+01/                                 
      data psa(2,8,1,3,6)/ 0.04000e+00/, psa(3,8,1,3,6)/ 0.0050e+00/    
      data psa(1,9,1,3,6)/ 0.10000e+01/                                 
      data psa(2,9,1,3,6)/ 0.00000e+00/, psa(3,9,1,3,6)/ 0.0050e+00/    
      data psa(1,10,1,3,6)/ 0.25000e+03/                                
      data psa(2,10,1,3,6)/ 0.85300e+01/, psa(3,10,1,3,6)/ 0.06000e+00/ 
      data psa(1,11,1,3,6)/ 0.35000e+03/                                
      data psa(2,11,1,3,6)/ 0.10960e+02/, psa(3,11,1,3,6)/ 0.14000e+00/ 
c        3d2                                                            
      data psa(1,1,2,3,6)/ 0.25000e+02/                                 
      data psa(2,1,2,3,6)/ 0.37080e+01/, psa(3,1,2,3,6)/ 0.00300e+00/   
      data psa(1,2,2,3,6)/ 0.50000e+02/                                 
      data psa(2,2,2,3,6)/ 0.89660e+01/, psa(3,2,2,3,6)/ 0.01700e+00/   
      data psa(1,3,2,3,6)/ 0.10000e+03/                                 
      data psa(2,3,2,3,6)/ 0.17280e+02/, psa(3,3,2,3,6)/ 0.06000e+00/   
      data psa(1,4,2,3,6)/ 0.15000e+03/                                 
      data psa(2,4,2,3,6)/ 0.22130e+02/, psa(3,4,2,3,6)/ 0.10000e+00/   
      data psa(1,5,2,3,6)/ 0.20000e+03/                                 
      data psa(2,5,2,3,6)/ 0.24510e+02/, psa(3,5,2,3,6)/ 0.11000e+00/   
      data psa(1,6,2,3,6)/ 0.30000e+03/                                 
      data psa(2,6,2,3,6)/ 0.25450e+02/, psa(3,6,2,3,6)/ 0.12000e+00/   
      data psa(1,7,2,3,6)/ 0.10000e+02/                                 
      data psa(2,7,2,3,6)/ 0.84600e+00/, psa(3,7,2,3,6)/ 0.00050e+00/   
      data psa(1,8,2,3,6)/ 0.50000e+01/                                 
      data psa(2,8,2,3,6)/ 0.22200e+00/, psa(3,8,2,3,6)/ 0.00050e+00/   
      data psa(1,9,2,3,6)/ 0.10000e+01/                                 
      data psa(2,9,2,3,6)/ 0.00600e+00/, psa(3,9,2,3,6)/ 0.00050e+00/   
      data psa(1,10,2,3,6)/ 0.25000e+03/                                
      data psa(2,10,2,3,6)/ 0.25400e+02/, psa(3,10,2,3,6)/ 0.11000e+00/ 
      data psa(1,11,2,3,6)/ 0.35000e+03/                                
      data psa(2,11,2,3,6)/ 0.25080e+02/, psa(3,11,2,3,6)/ 0.19000e+00/ 
c        3p2                                                            
      data psa(1,1,3,3,6)/ 0.25000e+02/                                 
      data psa(2,1,3,3,6)/ 0.25600e+01/, psa(3,1,3,3,6)/ 0.00800e+00/   
      data psa(1,2,3,3,6)/ 0.50000e+02/                                 
      data psa(2,2,3,3,6)/ 0.58900e+01/, psa(3,2,3,3,6)/ 0.01600e+00/   
      data psa(1,3,3,3,6)/ 0.10000e+03/                                 
      data psa(2,3,3,3,6)/ 0.10940e+02/, psa(3,3,3,3,6)/ 0.02500e+00/   
      data psa(1,4,3,3,6)/ 0.15000e+03/                                 
      data psa(2,4,3,3,6)/ 0.13840e+02/, psa(3,4,3,3,6)/ 0.03900e+00/   
      data psa(1,5,3,3,6)/ 0.20000e+03/                                 
      data psa(2,5,3,3,6)/ 0.15460e+02/, psa(3,5,3,3,6)/ 0.05200e+00/   
      data psa(1,6,3,3,6)/ 0.30000e+03/                                 
      data psa(2,6,3,3,6)/ 0.16950e+02/, psa(3,6,3,3,6)/ 0.10000e+00/   
      data psa(1,7,3,3,6)/ 0.10000e+02/                                 
      data psa(2,7,3,3,6)/ 0.71000e+00/, psa(3,7,3,3,6)/ 0.00500e+00/   
      data psa(1,8,3,3,6)/ 0.50000e+01/                                 
      data psa(2,8,3,3,6)/ 0.25000e+00/, psa(3,8,3,3,6)/ 0.00500e+00/   
      data psa(1,9,3,3,6)/ 0.10000e+01/                                 
      data psa(2,9,3,3,6)/ 0.02000e+00/, psa(3,9,3,3,6)/ 0.0050e+00/    
      data psa(1,10,3,3,6)/ 0.25000e+03/                                
      data psa(2,10,3,3,6)/ 0.16390e+02/, psa(3,10,3,3,6)/ 0.07000e+00/ 
      data psa(1,11,3,3,6)/ 0.35000e+03/                                
      data psa(2,11,3,3,6)/ 0.17310e+02/, psa(3,11,3,3,6)/ 0.15000e+00/ 
c        3f2                                                            
      data psa(1,1,4,3,6)/ 0.25000e+02/                                 
      data psa(2,1,4,3,6)/ 0.09000e+00/, psa(3,1,4,3,6)/ 0.0050e+00/    
      data psa(1,2,4,3,6)/ 0.50000e+02/                                 
      data psa(2,2,4,3,6)/ 0.30000e+00/, psa(3,2,4,3,6)/ 0.0050e+00/    
      data psa(1,3,4,3,6)/ 0.10000e+03/                                 
      data psa(2,3,4,3,6)/ 0.76000e+00/, psa(3,3,4,3,6)/ 0.00500e+00/   
      data psa(1,4,4,3,6)/ 0.15000e+03/                                 
      data psa(2,4,4,3,6)/ 0.11200e+01/, psa(3,4,4,3,6)/ 0.01400e+00/   
      data psa(1,5,4,3,6)/ 0.20000e+03/                                 
      data psa(2,5,4,3,6)/ 0.13300e+01/, psa(3,5,4,3,6)/ 0.03400e+00/   
      data psa(1,6,4,3,6)/ 0.30000e+03/                                 
      data psa(2,6,4,3,6)/ 0.11900e+01/, psa(3,6,4,3,6)/ 0.11000e+00/   
      data psa(1,7,4,3,6)/ 0.10000e+02/                                 
      data psa(2,7,4,3,6)/ 0.01000e+00/, psa(3,7,4,3,6)/ 0.0050e+00/    
      data psa(1,8,4,3,6)/ 0.50000e+01/                                 
      data psa(2,8,4,3,6)/ 0.00000e+00/, psa(3,8,4,3,6)/ 0.0050e+00/    
      data psa(1,9,4,3,6)/ 0.10000e+01/                                 
      data psa(2,9,4,3,6)/ 0.00000e+00/, psa(3,9,4,3,6)/ 0.0050e+00/    
      data psa(1,10,4,3,6)/ 0.25000e+03/                                
      data psa(2,10,4,3,6)/ 0.13500e+01/, psa(3,10,4,3,6)/ 0.06000e+00/ 
      data psa(1,11,4,3,6)/ 0.35000e+03/                                
      data psa(2,11,4,3,6)/ 0.87000e+00/, psa(3,11,4,3,6)/ 0.16000e+00/ 
c        e2                                                             
      data psa(1,1,5,3,6)/ 0.25000e+02/                                 
      data psa(2,1,5,3,6)/-0.76000e+00/, psa(3,1,5,3,6)/ 0.00500e+00/   
      data psa(1,2,5,3,6)/ 0.50000e+02/                                 
      data psa(2,2,5,3,6)/-0.16300e+01/, psa(3,2,5,3,6)/ 0.00500e+00/   
      data psa(1,3,5,3,6)/ 0.10000e+03/                                 
      data psa(2,3,5,3,6)/-0.25800e+01/, psa(3,3,5,3,6)/ 0.01700e+00/   
      data psa(1,4,5,3,6)/ 0.15000e+03/                                 
      data psa(2,4,5,3,6)/-0.28000e+01/, psa(3,4,5,3,6)/ 0.02900e+00/   
      data psa(1,5,5,3,6)/ 0.20000e+03/                                 
      data psa(2,5,5,3,6)/-0.27000e+01/, psa(3,5,5,3,6)/ 0.03700e+00/   
      data psa(1,6,5,3,6)/ 0.30000e+03/                                 
      data psa(2,6,5,3,6)/-0.23000e+01/, psa(3,6,5,3,6)/ 0.09000e+00/   
      data psa(1,7,5,3,6)/ 0.10000e+02/                                 
      data psa(2,7,5,3,6)/-0.18000e+00/, psa(3,7,5,3,6)/ 0.0050e+00/    
      data psa(1,8,5,3,6)/ 0.50000e+01/                                 
      data psa(2,8,5,3,6)/-0.05000e+00/, psa(3,8,5,3,6)/ 0.0050e+00/    
      data psa(1,9,5,3,6)/ 0.10000e+01/                                 
      data psa(2,9,5,3,6)/-0.00000e+00/, psa(3,9,5,3,6)/ 0.0050e+00/    
      data psa(1,10,5,3,6)/ 0.25000e+03/                                
      data psa(2,10,5,3,6)/-0.24900e+01/, psa(3,10,5,3,6)/ 0.04600e+00/ 
      data psa(1,11,5,3,6)/ 0.35000e+03/                                
      data psa(2,11,5,3,6)/-0.21800e+01/, psa(3,11,5,3,6)/ 0.11000e+00/ 
c        j = 3                                                          
c        1f3                                                            
      data psa(1,1,1,4,6)/ 0.25000e+02/                                 
      data psa(2,1,1,4,6)/-0.41500e+00/, psa(3,1,1,4,6)/ 0.00050e+00/   
      data psa(1,2,1,4,6)/ 0.50000e+02/                                 
      data psa(2,2,1,4,6)/-0.11010e+01/, psa(3,2,1,4,6)/ 0.00050e+00/   
      data psa(1,3,1,4,6)/ 0.10000e+03/                                 
      data psa(2,3,1,4,6)/-0.20890e+01/, psa(3,3,1,4,6)/ 0.00300e+00/   
      data psa(1,4,1,4,6)/ 0.15000e+03/                                 
      data psa(2,4,1,4,6)/-0.27020e+01/, psa(3,4,1,4,6)/ 0.01100e+00/   
      data psa(1,5,1,4,6)/ 0.20000e+03/                                 
      data psa(2,5,1,4,6)/-0.32350e+01/, psa(3,5,1,4,6)/ 0.02700e+00/   
      data psa(1,6,1,4,6)/ 0.30000e+03/                                 
      data psa(2,6,1,4,6)/-0.47200e+01/, psa(3,6,1,4,6)/ 0.08000e+00/   
      data psa(1,7,1,4,6)/ 0.10000e+02/                                 
      data psa(2,7,1,4,6)/-0.06600e+00/, psa(3,7,1,4,6)/ 0.00050e+00/   
      data psa(1,8,1,4,6)/ 0.50000e+01/                                 
      data psa(2,8,1,4,6)/-0.01100e+00/, psa(3,8,1,4,6)/ 0.00050e+00/   
      data psa(1,9,1,4,6)/ 0.10000e+01/                                 
      data psa(2,9,1,4,6)/-0.00000e+00/, psa(3,9,1,4,6)/ 0.00050e+00/   
      data psa(1,10,1,4,6)/ 0.25000e+03/                                
      data psa(2,10,1,4,6)/-0.38800e+01/, psa(3,10,1,4,6)/ 0.05000e+00/ 
      data psa(1,11,1,4,6)/ 0.35000e+03/                                
      data psa(2,11,1,4,6)/-0.58000e+01/, psa(3,11,1,4,6)/ 0.12000e+00/ 
c        3f3                                                            
      data psa(1,1,2,4,6)/ 0.25000e+02/                                 
      data psa(2,1,2,4,6)/-0.20000e+00/, psa(3,1,2,4,6)/ 0.0050e+00/    
      data psa(1,2,2,4,6)/ 0.50000e+02/                                 
      data psa(2,2,2,4,6)/-0.62000e+00/, psa(3,2,2,4,6)/ 0.0050e+00/    
      data psa(1,3,2,4,6)/ 0.10000e+03/                                 
      data psa(2,3,2,4,6)/-0.14100e+01/, psa(3,3,2,4,6)/ 0.00500e+00/   
      data psa(1,4,2,4,6)/ 0.15000e+03/                                 
      data psa(2,4,2,4,6)/-0.19800e+01/, psa(3,4,2,4,6)/ 0.01000e+00/   
      data psa(1,5,2,4,6)/ 0.20000e+03/                                 
      data psa(2,5,2,4,6)/-0.23600e+01/, psa(3,5,2,4,6)/ 0.02500e+00/   
      data psa(1,6,2,4,6)/ 0.30000e+03/                                 
      data psa(2,6,2,4,6)/-0.27300e+01/, psa(3,6,2,4,6)/ 0.11000e+00/   
      data psa(1,7,2,4,6)/ 0.10000e+02/                                 
      data psa(2,7,2,4,6)/-0.03000e+00/, psa(3,7,2,4,6)/ 0.0050e+00/    
      data psa(1,8,2,4,6)/ 0.50000e+01/                                 
      data psa(2,8,2,4,6)/-0.00000e+00/, psa(3,8,2,4,6)/ 0.0050e+00/    
      data psa(1,9,2,4,6)/ 0.10000e+01/                                 
      data psa(2,9,2,4,6)/-0.00000e+00/, psa(3,9,2,4,6)/ 0.0050e+00/    
      data psa(1,10,2,4,6)/ 0.25000e+03/                                
      data psa(2,10,2,4,6)/-0.26000e+01/, psa(3,10,2,4,6)/ 0.04900e+00/ 
      data psa(1,11,2,4,6)/ 0.35000e+03/                                
      data psa(2,11,2,4,6)/-0.27600e+01/, psa(3,11,2,4,6)/ 0.13000e+00/ 
c        3d3                                                            
      data psa(1,1,3,4,6)/ 0.25000e+02/                                 
      data psa(2,1,3,4,6)/ 0.04800e+00/, psa(3,1,3,4,6)/ 0.00050e+00/   
      data psa(1,2,3,4,6)/ 0.50000e+02/                                 
      data psa(2,2,3,4,6)/ 0.32400e+00/, psa(3,2,3,4,6)/ 0.00300e+00/   
      data psa(1,3,3,4,6)/ 0.10000e+03/                                 
      data psa(2,3,3,4,6)/ 0.14570e+01/, psa(3,3,3,4,6)/ 0.01400e+00/   
      data psa(1,4,3,4,6)/ 0.15000e+03/                                 
      data psa(2,4,3,4,6)/ 0.27360e+01/, psa(3,4,3,4,6)/ 0.03300e+00/   
      data psa(1,5,3,4,6)/ 0.20000e+03/                                 
      data psa(2,5,3,4,6)/ 0.37400e+01/, psa(3,5,3,4,6)/ 0.06000e+00/   
      data psa(1,6,3,4,6)/ 0.30000e+03/                                 
      data psa(2,6,3,4,6)/ 0.46200e+01/, psa(3,6,3,4,6)/ 0.12000e+00/   
      data psa(1,7,3,4,6)/ 0.10000e+02/                                 
      data psa(2,7,3,4,6)/ 0.00600e+00/, psa(3,7,3,4,6)/ 0.00050e+00/   
      data psa(1,8,3,4,6)/ 0.50000e+01/                                 
      data psa(2,8,3,4,6)/ 0.00200e+00/, psa(3,8,3,4,6)/ 0.00050e+00/   
      data psa(1,9,3,4,6)/ 0.10000e+01/                                 
      data psa(2,9,3,4,6)/ 0.00000e+00/, psa(3,9,3,4,6)/ 0.00050e+00/   
      data psa(1,10,3,4,6)/ 0.25000e+03/                                
      data psa(2,10,3,4,6)/ 0.43700e+01/, psa(3,10,3,4,6)/ 0.09000e+00/ 
      data psa(1,11,3,4,6)/ 0.35000e+03/                                
      data psa(2,11,3,4,6)/ 0.46000e+01/, psa(3,11,3,4,6)/ 0.16000e+00/ 
c        3g3                                                            
      data psa(1,1,4,4,6)/ 0.25000e+02/                                 
      data psa(2,1,4,4,6)/-0.05300e+00/, psa(3,1,4,4,6)/ 0.00050e+00/   
      data psa(1,2,4,4,6)/ 0.50000e+02/                                 
      data psa(2,2,4,4,6)/-0.25800e+00/, psa(3,2,4,4,6)/ 0.00050e+00/   
      data psa(1,3,4,4,6)/ 0.10000e+03/                                 
      data psa(2,3,4,4,6)/-0.93400e+00/, psa(3,3,4,4,6)/ 0.00050e+00/   
      data psa(1,4,4,4,6)/ 0.15000e+03/                                 
      data psa(2,4,4,4,6)/-0.17370e+01/, psa(3,4,4,4,6)/ 0.00050e+00/   
      data psa(1,5,4,4,6)/ 0.20000e+03/                                 
      data psa(2,5,4,4,6)/-0.25340e+01/, psa(3,5,4,4,6)/ 0.00100e+00/   
      data psa(1,6,4,4,6)/ 0.30000e+03/                                 
      data psa(2,6,4,4,6)/-0.39020e+01/, psa(3,6,4,4,6)/ 0.00400e+00/   
      data psa(1,7,4,4,6)/ 0.10000e+02/                                 
      data psa(2,7,4,4,6)/-0.00400e+00/, psa(3,7,4,4,6)/ 0.00050e+00/   
      data psa(1,8,4,4,6)/ 0.50000e+01/                                 
      data psa(2,8,4,4,6)/-0.00000e+00/, psa(3,8,4,4,6)/ 0.00050e+00/   
      data psa(1,9,4,4,6)/ 0.10000e+01/                                 
      data psa(2,9,4,4,6)/-0.00000e+00/, psa(3,9,4,4,6)/ 0.00050e+00/   
      data psa(1,10,4,4,6)/ 0.25000e+03/                                
      data psa(2,10,4,4,6)/-0.32650e+01/, psa(3,10,4,4,6)/ 0.00100e+00/ 
      data psa(1,11,4,4,6)/ 0.35000e+03/                                
      data psa(2,11,4,4,6)/-0.44400e+01/, psa(3,11,4,4,6)/ 0.00700e+00/ 
c        e3                                                             
      data psa(1,1,5,4,6)/ 0.25000e+02/                                 
      data psa(2,1,5,4,6)/ 0.54900e+00/, psa(3,1,5,4,6)/ 0.00050e+00/   
      data psa(1,2,5,4,6)/ 0.50000e+02/                                 
      data psa(2,2,5,4,6)/ 0.16000e+01/, psa(3,2,5,4,6)/ 0.00050e+00/   
      data psa(1,3,5,4,6)/ 0.10000e+03/                                 
      data psa(2,3,5,4,6)/ 0.34690e+01/, psa(3,3,5,4,6)/ 0.00200e+00/   
      data psa(1,4,5,4,6)/ 0.15000e+03/                                 
      data psa(2,4,5,4,6)/ 0.48040e+01/, psa(3,4,5,4,6)/ 0.00800e+00/   
      data psa(1,5,5,4,6)/ 0.20000e+03/                                 
      data psa(2,5,5,4,6)/ 0.57230e+01/, psa(3,5,5,4,6)/ 0.02000e+00/   
      data psa(1,6,5,4,6)/ 0.30000e+03/                                 
      data psa(2,6,5,4,6)/ 0.68000e+01/, psa(3,6,5,4,6)/ 0.06000e+00/   
      data psa(1,7,5,4,6)/ 0.10000e+02/                                 
      data psa(2,7,5,4,6)/ 0.08100e+00/, psa(3,7,5,4,6)/ 0.00050e+00/   
      data psa(1,8,5,4,6)/ 0.50000e+01/                                 
      data psa(2,8,5,4,6)/ 0.01300e+00/, psa(3,8,5,4,6)/ 0.00050e+00/   
      data psa(1,9,5,4,6)/ 0.10000e+01/                                 
      data psa(2,9,5,4,6)/ 0.00000e+00/, psa(3,9,5,4,6)/ 0.00050e+00/   
      data psa(1,10,5,4,6)/ 0.25000e+03/                                
      data psa(2,10,5,4,6)/ 0.63570e+01/, psa(3,10,5,4,6)/ 0.03700e+00/ 
      data psa(1,11,5,4,6)/ 0.35000e+03/                                
      data psa(2,11,5,4,6)/ 0.71300e+01/, psa(3,11,5,4,6)/ 0.09000e+00/ 
c        j = 4                                                          
c        1g4                                                            
      data psa(1,1,1,5,6)/ 0.25000e+02/                                 
      data psa(2,1,1,5,6)/ 0.03000e+00/, psa(3,1,1,5,6)/ 0.0050e+00/    
      data psa(1,2,1,5,6)/ 0.50000e+02/                                 
      data psa(2,2,1,5,6)/ 0.13000e+00/, psa(3,2,1,5,6)/ 0.0050e+00/    
      data psa(1,3,1,5,6)/ 0.10000e+03/                                 
      data psa(2,3,1,5,6)/ 0.39000e+00/, psa(3,3,1,5,6)/ 0.00500e+00/   
      data psa(1,4,1,5,6)/ 0.15000e+03/                                 
      data psa(2,4,1,5,6)/ 0.68000e+00/, psa(3,4,1,5,6)/ 0.00500e+00/   
      data psa(1,5,1,5,6)/ 0.20000e+03/                                 
      data psa(2,5,1,5,6)/ 0.98000e+00/, psa(3,5,1,5,6)/ 0.01000e+00/   
      data psa(1,6,1,5,6)/ 0.30000e+03/                                 
      data psa(2,6,1,5,6)/ 0.15200e+01/, psa(3,6,1,5,6)/ 0.04800e+00/   
      data psa(1,7,1,5,6)/ 0.10000e+02/                                 
      data psa(2,7,1,5,6)/ 0.00000e+00/, psa(3,7,1,5,6)/ 0.0050e+00/    
      data psa(1,8,1,5,6)/ 0.50000e+01/                                 
      data psa(2,8,1,5,6)/ 0.00000e+00/, psa(3,8,1,5,6)/ 0.0050e+00/    
      data psa(1,9,1,5,6)/ 0.10000e+01/                                 
      data psa(2,9,1,5,6)/ 0.00000e+00/, psa(3,9,1,5,6)/ 0.0050e+00/    
      data psa(1,10,1,5,6)/ 0.25000e+03/                                
      data psa(2,10,1,5,6)/ 0.12800e+01/, psa(3,10,1,5,6)/ 0.02400e+00/ 
      data psa(1,11,1,5,6)/ 0.35000e+03/                                
      data psa(2,11,1,5,6)/ 0.16700e+01/, psa(3,11,1,5,6)/ 0.08000e+00/ 
c        3g4                                                            
      data psa(1,1,2,5,6)/ 0.25000e+02/                                 
      data psa(2,1,2,5,6)/ 0.16900e+00/, psa(3,1,2,5,6)/ 0.00050e+00/   
      data psa(1,2,2,5,6)/ 0.50000e+02/                                 
      data psa(2,2,2,5,6)/ 0.71600e+00/, psa(3,2,2,5,6)/ 0.00050e+00/   
      data psa(1,3,2,5,6)/ 0.10000e+03/                                 
      data psa(2,3,2,5,6)/ 0.21540e+01/, psa(3,3,2,5,6)/ 0.00050e+00/   
      data psa(1,4,2,5,6)/ 0.15000e+03/                                 
      data psa(2,4,2,5,6)/ 0.36180e+01/, psa(3,4,2,5,6)/ 0.00050e+00/   
      data psa(1,5,2,5,6)/ 0.20000e+03/                                 
      data psa(2,5,2,5,6)/ 0.49870e+01/, psa(3,5,2,5,6)/ 0.00050e+00/   
      data psa(1,6,2,5,6)/ 0.30000e+03/                                 
      data psa(2,6,2,5,6)/ 0.73370e+01/, psa(3,6,2,5,6)/ 0.00050e+00/   
      data psa(1,7,2,5,6)/ 0.10000e+02/                                 
      data psa(2,7,2,5,6)/ 0.01400e+00/, psa(3,7,2,5,6)/ 0.00050e+00/   
      data psa(1,8,2,5,6)/ 0.50000e+01/                                 
      data psa(2,8,2,5,6)/ 0.00100e+00/, psa(3,8,2,5,6)/ 0.00050e+00/   
      data psa(1,9,2,5,6)/ 0.10000e+01/                                 
      data psa(2,9,2,5,6)/ 0.00000e+00/, psa(3,9,2,5,6)/ 0.00050e+00/   
      data psa(1,10,2,5,6)/ 0.25000e+03/                                
      data psa(2,10,2,5,6)/ 0.62320e+01/, psa(3,10,2,5,6)/ 0.00050e+00/ 
      data psa(1,11,2,5,6)/ 0.35000e+03/                                
      data psa(2,11,2,5,6)/ 0.82940e+01/, psa(3,11,2,5,6)/ 0.00050e+00/ 
c        3f4                                                            
      data psa(1,1,3,5,6)/ 0.25000e+02/                                 
      data psa(2,1,3,5,6)/ 0.02000e+00/, psa(3,1,3,5,6)/ 0.0050e+00/    
      data psa(1,2,3,5,6)/ 0.50000e+02/                                 
      data psa(2,2,3,5,6)/ 0.10000e+00/, psa(3,2,3,5,6)/ 0.0050e+00/    
      data psa(1,3,3,5,6)/ 0.10000e+03/                                 
      data psa(2,3,3,5,6)/ 0.45000e+00/, psa(3,3,3,5,6)/ 0.00700e+00/   
      data psa(1,4,3,5,6)/ 0.15000e+03/                                 
      data psa(2,4,3,5,6)/ 0.99000e+00/, psa(3,4,3,5,6)/ 0.02200e+00/   
      data psa(1,5,3,5,6)/ 0.20000e+03/                                 
      data psa(2,5,3,5,6)/ 0.16300e+01/, psa(3,5,3,5,6)/ 0.03900e+00/   
      data psa(1,6,3,5,6)/ 0.30000e+03/                                 
      data psa(2,6,3,5,6)/ 0.28100e+01/, psa(3,6,3,5,6)/ 0.06000e+00/   
      data psa(1,7,3,5,6)/ 0.10000e+02/                                 
      data psa(2,7,3,5,6)/ 0.00000e+00/, psa(3,7,3,5,6)/ 0.0050e+00/    
      data psa(1,8,3,5,6)/ 0.50000e+01/                                 
      data psa(2,8,3,5,6)/ 0.00000e+00/, psa(3,8,3,5,6)/ 0.0050e+00/    
      data psa(1,9,3,5,6)/ 0.10000e+01/                                 
      data psa(2,9,3,5,6)/ 0.00000e+00/, psa(3,9,3,5,6)/ 0.0050e+00/    
      data psa(1,10,3,5,6)/ 0.25000e+03/                                
      data psa(2,10,3,5,6)/ 0.22600e+01/, psa(3,10,3,5,6)/ 0.05100e+00/ 
      data psa(1,11,3,5,6)/ 0.35000e+03/                                
      data psa(2,11,3,5,6)/ 0.32100e+01/, psa(3,11,3,5,6)/ 0.11000e+00/ 
c        3h4                                                            
      data psa(1,1,4,5,6)/ 0.25000e+02/                                 
      data psa(2,1,4,5,6)/ 0.00000e+00/, psa(3,1,4,5,6)/ 0.0050e+00/    
      data psa(1,2,4,5,6)/ 0.50000e+02/                                 
      data psa(2,2,4,5,6)/ 0.02000e+00/, psa(3,2,4,5,6)/ 0.0050e+00/    
      data psa(1,3,4,5,6)/ 0.10000e+03/                                 
      data psa(2,3,4,5,6)/ 0.09000e+00/, psa(3,3,4,5,6)/ 0.0050e+00/    
      data psa(1,4,4,5,6)/ 0.15000e+03/                                 
      data psa(2,4,4,5,6)/ 0.19000e+00/, psa(3,4,4,5,6)/ 0.0050e+00/    
      data psa(1,5,4,5,6)/ 0.20000e+03/                                 
      data psa(2,5,4,5,6)/ 0.29000e+00/, psa(3,5,4,5,6)/ 0.0050e+00/    
      data psa(1,6,4,5,6)/ 0.30000e+03/                                 
      data psa(2,6,4,5,6)/ 0.48000e+00/, psa(3,6,4,5,6)/ 0.0050e+00/    
      data psa(1,7,4,5,6)/ 0.10000e+02/                                 
      data psa(2,7,4,5,6)/ 0.00000e+00/, psa(3,7,4,5,6)/ 0.0050e+00/    
      data psa(1,8,4,5,6)/ 0.50000e+01/                                 
      data psa(2,8,4,5,6)/ 0.00000e+00/, psa(3,8,4,5,6)/ 0.0050e+00/    
      data psa(1,9,4,5,6)/ 0.10000e+01/                                 
      data psa(2,9,4,5,6)/ 0.00000e+00/, psa(3,9,4,5,6)/ 0.0050e+00/    
      data psa(1,10,4,5,6)/ 0.25000e+03/                                
      data psa(2,10,4,5,6)/ 0.39000e+00/, psa(3,10,4,5,6)/ 0.0050e+00/  
      data psa(1,11,4,5,6)/ 0.35000e+03/                                
      data psa(2,11,4,5,6)/ 0.56000e+00/, psa(3,11,4,5,6)/ 0.0050e+00/  
c        e4                                                             
      data psa(1,1,5,5,6)/ 0.25000e+02/                                 
      data psa(2,1,5,5,6)/-0.04000e+00/, psa(3,1,5,5,6)/ 0.0050e+00/    
      data psa(1,2,5,5,6)/ 0.50000e+02/                                 
      data psa(2,2,5,5,6)/-0.17000e+00/, psa(3,2,5,5,6)/ 0.0050e+00/    
      data psa(1,3,5,5,6)/ 0.10000e+03/                                 
      data psa(2,3,5,5,6)/-0.49000e+00/, psa(3,3,5,5,6)/ 0.0050e+00/    
      data psa(1,4,5,5,6)/ 0.15000e+03/                                 
      data psa(2,4,5,5,6)/-0.79000e+00/, psa(3,4,5,5,6)/ 0.0050e+00/    
      data psa(1,5,5,5,6)/ 0.20000e+03/                                 
      data psa(2,5,5,5,6)/-0.10500e+01/, psa(3,5,5,5,6)/ 0.0050e+00/    
      data psa(1,6,5,5,6)/ 0.30000e+03/                                 
      data psa(2,6,5,5,6)/-0.14200e+01/, psa(3,6,5,5,6)/ 0.0050e+00/    
      data psa(1,7,5,5,6)/ 0.10000e+02/                                 
      data psa(2,7,5,5,6)/-0.00000e+00/, psa(3,7,5,5,6)/ 0.0050e+00/    
      data psa(1,8,5,5,6)/ 0.50000e+01/                                 
      data psa(2,8,5,5,6)/-0.00000e+00/, psa(3,8,5,5,6)/ 0.0050e+00/    
      data psa(1,9,5,5,6)/ 0.10000e+01/                                 
      data psa(2,9,5,5,6)/-0.00000e+00/, psa(3,9,5,5,6)/ 0.0050e+00/    
      data psa(1,10,5,5,6)/ 0.25000e+03/                                
      data psa(2,10,5,5,6)/-0.12600e+01/, psa(3,10,5,5,6)/ 0.0050e+00/  
      data psa(1,11,5,5,6)/ 0.35000e+03/                                
      data psa(2,11,5,5,6)/-0.15400e+01/, psa(3,11,5,5,6)/ 0.0050e+00/  
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c *****  phase-shift analysis  ni95, ips=7                              
c        --------------------------------                               
c        ni95: nijmegen n-p single-energy analysis,                     
c              stoks privat communication of 4/4/95                     
c                                                                       
c                                                                       
c                                                                       
c        j = 0                                                          
c        1s0 np                                                         
      data psa(1,1,1,1,7)/ 0.25000e+02/                                 
      data psa(2,1,1,1,7)/ 0.50698e+02/, psa(3,1,1,1,7)/ 0.63750e+00/   
      data psa(1,2,1,1,7)/ 0.50000e+02/                                 
      data psa(2,2,1,1,7)/ 0.41361e+02/, psa(3,2,1,1,7)/ 0.90010e+00/   
      data psa(1,3,1,1,7)/ 0.10000e+03/                                 
      data psa(2,3,1,1,7)/ 0.30682e+02/, psa(3,3,1,1,7)/ 3.09080e+00/   
      data psa(1,4,1,1,7)/ 0.15000e+03/                                 
      data psa(2,4,1,1,7)/ 0.17839e+02/, psa(3,4,1,1,7)/ 1.33250e+00/   
      data psa(1,5,1,1,7)/ 0.21500e+03/                                 
      data psa(2,5,1,1,7)/ 0.74792e+01/, psa(3,5,1,1,7)/ 0.93860e+00/   
      data psa(1,6,1,1,7)/ 0.32000e+03/                                 
      data psa(2,6,1,1,7)/-0.72577e+01/, psa(3,6,1,1,7)/ 0.57180e+00/   
      data psa(1,7,1,1,7)/ 0.10000e+02/                                 
      data psa(2,7,1,1,7)/ 0.62482e+02/, psa(3,7,1,1,7)/ 1.26390e+00/   
      data psa(1,8,1,1,7)/ 0.10000e+01/                                 
      data psa(2,8,1,1,7)/ 0.62119e+02/, psa(3,8,1,1,7)/ 0.79510e+00/   
      data psa(1,9,1,1,7)/ 0.50000e+01/                                 
      data psa(2,9,1,1,7)/ 0.63485e+02/, psa(3,9,1,1,7)/ 0.15390e+00/   
      data psa(1,10,1,1,7)/ 0.38254e+00/                                
      data psa(2,10,1,1,7)/ 0.54575e+02/, psa(3,10,1,1,7)/ 0.00590e+00/ 
c        j = 1                                                          
c        1p1                                                            
      data psa(1,1,1,2,7)/ 0.25000e+02/                                 
      data psa(2,1,1,2,7)/-0.66855e+01/, psa(3,1,1,2,7)/ 0.12990e+00/   
      data psa(1,2,1,2,7)/ 0.50000e+02/                                 
      data psa(2,2,1,2,7)/-0.93598e+01/, psa(3,2,1,2,7)/ 0.16850e+00/   
      data psa(1,3,1,2,7)/ 0.10000e+03/                                 
      data psa(2,3,1,2,7)/-0.10599e+02/, psa(3,3,1,2,7)/ 1.72220e+00/   
      data psa(1,4,1,2,7)/ 0.15000e+03/                                 
      data psa(2,4,1,2,7)/-0.17940e+02/, psa(3,4,1,2,7)/ 0.56340e+00/   
      data psa(1,5,1,2,7)/ 0.21500e+03/                                 
      data psa(2,5,1,2,7)/-0.22695e+02/, psa(3,5,1,2,7)/ 0.47010e+00/   
      data psa(1,6,1,2,7)/ 0.32000e+03/                                 
      data psa(2,6,1,2,7)/-0.28314e+02/, psa(3,6,1,2,7)/ 0.39230e+00/   
      data psa(1,7,1,2,7)/ 0.10000e+02/                                 
      data psa(2,7,1,2,7)/-0.30280e+01/, psa(3,7,1,2,7)/ 0.27840e+00/   
      data psa(1,8,1,2,7)/ 0.50000e+01/                                 
      data psa(2,8,1,2,7)/-0.16774e+01/, psa(3,8,1,2,7)/ 1.40680e+00/   
c        3s1                                                            
      data psa(1,1,3,2,7)/ 0.25000e+02/                                 
      data psa(2,1,3,2,7)/ 0.79950e+02/, psa(3,1,3,2,7)/ 0.59910e+00/   
      data psa(1,2,3,2,7)/ 0.50000e+02/                                 
      data psa(2,2,3,2,7)/ 0.62217e+02/, psa(3,2,3,2,7)/ 0.36340e+00/   
      data psa(1,3,3,2,7)/ 0.10000e+03/                                 
      data psa(2,3,3,2,7)/ 0.41615e+02/, psa(3,3,3,2,7)/ 1.07070e+00/   
      data psa(1,4,3,2,7)/ 0.15000e+03/                                 
      data psa(2,4,3,2,7)/ 0.30119e+02/, psa(3,4,3,2,7)/ 0.34830e+00/   
      data psa(1,5,3,2,7)/ 0.21500e+03/                                 
      data psa(2,5,3,2,7)/ 0.18742e+02/, psa(3,5,3,2,7)/ 0.28820e+00/   
      data psa(1,6,3,2,7)/ 0.32000e+03/                                 
      data psa(2,6,3,2,7)/ 0.43715e+01/, psa(3,6,3,2,7)/ 0.38220e+00/   
      data psa(1,7,3,2,7)/ 0.10000e+02/                                 
      data psa(2,7,3,2,7)/-0.75187e+02/, psa(3,7,3,2,7)/ 0.69910e+00/   
      data psa(1,8,3,2,7)/ 0.10000e+01/                                 
      data psa(2,8,3,2,7)/-0.32326e+02/, psa(3,8,3,2,7)/ 0.23120e+00/   
      data psa(1,9,3,2,7)/ 0.38254e+00/                                 
      data psa(2,9,3,2,7)/-0.20614e+02/, psa(3,9,3,2,7)/ 0.00380e+00/   
c        3d1                                                            
      data psa(1,1,4,2,7)/ 0.25000e+02/                                 
      data psa(2,1,4,2,7)/-0.28200e+01/, psa(3,1,4,2,7)/ 0.01960e+00/   
      data psa(1,2,4,2,7)/ 0.50000e+02/                                 
      data psa(2,2,4,2,7)/-0.64909e+01/, psa(3,2,4,2,7)/ 0.06590e+00/   
      data psa(1,3,4,2,7)/ 0.10000e+03/                                 
      data psa(2,3,4,2,7)/-0.12065e+02/, psa(3,3,4,2,7)/ 0.30860e+00/   
      data psa(1,4,4,2,7)/ 0.15000e+03/                                 
      data psa(2,4,4,2,7)/-0.16072e+02/, psa(3,4,4,2,7)/ 0.28920e+00/   
      data psa(1,5,4,2,7)/ 0.21500e+03/                                 
      data psa(2,5,4,2,7)/-0.20590e+02/, psa(3,5,4,2,7)/ 0.25580e+00/   
      data psa(1,6,4,2,7)/ 0.32000e+03/                                 
      data psa(2,6,4,2,7)/-0.24833e+02/, psa(3,6,4,2,7)/ 0.27510e+00/   
c        e1                                                             
      data psa(1,1,5,2,7)/ 0.25000e+02/                                 
      data psa(2,1,5,2,7)/ 0.21379e+01/, psa(3,1,5,2,7)/ 0.33160e+00/   
      data psa(1,2,5,2,7)/ 0.50000e+02/                                 
      data psa(2,2,5,2,7)/ 0.25689e+01/, psa(3,2,5,2,7)/ 0.36210e+00/   
      data psa(1,3,5,2,7)/ 0.10000e+03/                                 
      data psa(2,3,5,2,7)/ 0.32701e+01/, psa(3,3,5,2,7)/ 1.48150e+00/   
      data psa(1,4,5,2,7)/ 0.15000e+03/                                 
      data psa(2,4,5,2,7)/ 0.24684e+01/, psa(3,4,5,2,7)/ 0.46960e+00/   
      data psa(1,5,5,2,7)/ 0.21500e+03/                                 
      data psa(2,5,5,2,7)/ 0.27004e+01/, psa(3,5,5,2,7)/ 0.32720e+00/   
      data psa(1,6,5,2,7)/ 0.32000e+03/                                 
      data psa(2,6,5,2,7)/ 0.46504e+01/, psa(3,6,5,2,7)/ 0.25290e+00/   
c        j = 2                                                          
c        3d2                                                            
      data psa(1,1,2,3,7)/ 0.50000e+02/                                 
      data psa(2,1,2,3,7)/ 0.91369e+01/, psa(3,1,2,3,7)/ 0.13840e+00/   
      data psa(1,2,2,3,7)/ 0.10000e+03/                                 
      data psa(2,2,2,3,7)/ 0.18641e+02/, psa(3,2,2,3,7)/ 0.51090e+00/   
      data psa(1,3,2,3,7)/ 0.15000e+03/                                 
      data psa(2,3,2,3,7)/ 0.22828e+02/, psa(3,3,2,3,7)/ 0.24550e+00/   
      data psa(1,4,2,3,7)/ 0.21500e+03/                                 
      data psa(2,4,2,3,7)/ 0.25218e+02/, psa(3,4,2,3,7)/ 0.23080e+00/   
      data psa(1,5,2,3,7)/ 0.32000e+03/                                 
      data psa(2,5,2,3,7)/ 0.25229e+02/, psa(3,5,2,3,7)/ 0.24540e+00/   
c        j = 3                                                          
c        1f3                                                            
      data psa(1,1,1,4,7)/ 0.15000e+03/                                 
      data psa(2,1,1,4,7)/-0.26700e+01/, psa(3,1,1,4,7)/ 0.16670e+00/   
      data psa(1,2,1,4,7)/ 0.21500e+03/                                 
      data psa(2,2,1,4,7)/-0.33818e+01/, psa(3,2,1,4,7)/ 0.11410e+00/   
      data psa(1,3,1,4,7)/ 0.32000e+03/                                 
      data psa(2,3,1,4,7)/-0.51727e+01/, psa(3,3,1,4,7)/ 0.10970e+00/   
c        3d3                                                            
      data psa(1,1,3,4,7)/ 0.15000e+03/                                 
      data psa(2,1,3,4,7)/ 0.21549e+01/, psa(3,1,3,4,7)/ 0.29950e+00/   
      data psa(1,2,3,4,7)/ 0.21500e+03/                                 
      data psa(2,2,3,4,7)/ 0.39382e+01/, psa(3,2,3,4,7)/ 0.21730e+00/   
      data psa(1,3,3,4,7)/ 0.32000e+03/                                 
      data psa(2,3,3,4,7)/ 0.48425e+01/, psa(3,3,3,4,7)/ 0.21000e+00/   
c        3g3                                                            
      data psa(1,1,4,4,7)/ 0.21500e+03/                                 
      data psa(2,1,4,4,7)/-0.25256e+01/, psa(3,1,4,4,7)/ 0.23030e+00/   
      data psa(1,2,4,4,7)/ 0.32000e+03/                                 
      data psa(2,2,4,4,7)/-0.44874e+01/, psa(3,2,4,4,7)/ 0.19770e+00/   
c        e3                                                             
      data psa(1,1,5,4,7)/ 0.15000e+03/                                 
      data psa(2,1,5,4,7)/ 0.42086e+01/, psa(3,1,5,4,7)/ 0.16420e+00/   
      data psa(1,2,5,4,7)/ 0.21500e+03/                                 
      data psa(2,2,5,4,7)/ 0.59509e+01/, psa(3,2,5,4,7)/ 0.11560e+00/   
      data psa(1,3,5,4,7)/ 0.32000e+03/                                 
      data psa(2,3,5,4,7)/ 0.69315e+01/, psa(3,3,5,4,7)/ 0.09300e+00/   
c        j = 4                                                          
c        3g4                                                            
      data psa(1,1,2,5,7)/ 0.21500e+03/                                 
      data psa(2,1,2,5,7)/ 0.54862e+01/, psa(3,1,2,5,7)/ 0.11180e+00/   
      data psa(1,2,2,5,7)/ 0.32000e+03/                                 
      data psa(2,2,2,5,7)/ 0.77934e+01/, psa(3,2,2,5,7)/ 0.12590e+00/   
c                                                                       
c                                                                       
c                                                                       
c                                                                       
10000 format (2a4,a2,20i3)                                              
10001 format (1h ,2a4,a2,20i3)                                          
10002 format (2a4,a2,6f10.5)                                            
10003 format (1h ,2a4,a2,6f10.5)                                        
10004 format (' input-parameters for phases'/1h ,27(1h-))               
10005 format (//' transf gauss pts and wghts for c =',f8.2,             
     1', cut =',f8.2,' and n =',i3/' points')                           
10006 format (7x,4f15.4)                                                
10007 format (' weights')                                               
10008 format (2a4,a2,15a4)                                              
10009 format (1h ,2a4,a2,15a4)                                          
10010 format (1h ,f8.2,f16.5,2x,3f9.3,1x,a4,f16.2,f21.7,6x,3f8.4)       
10011 format (/1h ,60x,a1,1x,a1,i2,3x,a4/)                              
10012 format (/1h ,60x,'1 s 0',3x,'icoul =',i2/)                        
10013 format (///' error in phases. matrix inversion error index =',    
     1 i5///)                                                           
10014 format (1h ,2a4,a2,6f10.5)                                        
10015 format (1h ,49x,'total chi**2',f13.2/50x,'chi**2/datum',f13.2/    
     1 50x,'for',i5,' p. of data')                                      
10016 format (1h /' low energy parameters    a',a1,' =',                
     1 f10.4,'    r',a1,' =',f10.4)                                     
10021 format ('state',5x,a1,a1,i1)                                      
10022 format ('state',5x,a1,i1)                                         
10024 format (a4,1x,2a1,i1,2x,3d20.7)                                   
10025 format (a4,1x,2a1,i1,i2,d20.7,2d20.7)                             
10026 format (80x)                                                      
10031 format ('nx',8x,' 17')                                            
10032 format ('namex',5x,'lab. energy (mev)')                           
10033 format ('ny',8x,' 17')                                            
10034 format ('namey',5x,'phase shift (rad)')                           
10035 format ('ny',8x,' 22')                                            
10036 format ('namey',5x,'mixing parameter (rad)')                      
10037 format ('namey',5x,'phase shift (deg)')                           
10038 format ('namey',5x,'mixing parameter (deg)')                      
10040 format ('phase-shifts for    ',15a4)                              
10041 format (f7.3,i3,5d14.6)                                           
10050 format (/' elab(mev)',3x,3(2a1,i2,8x),' e',i2,8x,2a1,i2)          
10051 format (1h ,f8.2,5e14.6)                                          
10052 format (/' elab(mev)',3x,'1s0',46x,'3p0')                         
10053 format (///' p h a s e - s h i f t s (radians)'/1h ,33(1h-))      
10054 format (///' p h a s e - s h i f t s (degrees)'/1h ,33(1h-))      
10055 format (f7.2,5f12.5)                                              
10061 format (//2x,a1,1x,a1,i2,3x,a4/)                                  
10062 format (/1h ,' 1 s 0',3x,'icoul =',i2/)                           
10072 format (2a4,a2,6f10.2)                                            
10073 format (1h ,2a4,a2,6f10.2)                                        
10100 format (4d20.13)                                                  
10110 format (a4,2x,2a1,'-',a1,i1,3d20.12)                              
10111 format (i3)                                                       
10112 format (f10.4)                                                    
10113 format (4d20.12)                                                  
10114 format ()                                                         
10115 format (/' blatt-biedenharn conventions used')                    
c                                                                       
10200 format (f9.2,8x,f9.2,f7.2,f9.2,f7.2)                              
c                                                                       
11000 format (i3)                                                       
11001 format (f10.5,i3)                                                 
11002 format (21x,a3,2x,f20.16)                                         
11022 format (21x,a5,2x,f20.16)                                         
11003 format (6d13.6)                                                   
11005 format (15x,5f10.6)                                               
11006 format (//' warning in phascoul: the range given for j'/          
     1 ' is insufficient to apply the error matrix. the chi**2'/        
     2 ' in regard to the error matrix is not calculated.'/             
     3 ' execution continud.'//)                                        
11010 format (///' ch**2 using the nijmegen pp-only error matrix'/      
     1           ' ---------------------------------------------'/)     
11020 format (///' ch**2 using the nijmegen pp+np error matrix'/        
     1           ' -------------------------------------------'/)       
11011 format (/' single-energy bin',f10.5,' mev'/                       
     1 ' partial contributions:')                                       
11012 format (' ',a3,3f15.5)                                            
11032 format (' ',a5,3f15.5)                                            
11013 format ('  chi**2 for',f10.5,' mev s. e. bin:',3f13.5)            
11014 format (/' total chi**2 using the error matrix:',3f13.5/          
     1         ' ------------------------------------')                 
11024 format (/' total chi**2 for the first six bins:',26x,f13.5/       
     1         ' ------------------------------------')                 
11100 format (5d16.9)                                                   
c                                                                       
c                                                                       
c                                                                       
      go to 8888
c                                                                       
c        read in phase shift analysis #1, ips=1                         
c        --------------------------------------                         
c                                                                       
c                                                                       
      kda8=kda(8)                                                       
      read (kda8,10000) labps(1)                                        
c**** write(kwrite,10000) labps(1)                                      
      read (kda8,10000) name,nskip                                      
      j1=0                                                              
   21 j1=j1+1                                                           
      iv=0                                                              
   22 iv=iv+1                                                           
c                                                                       
      if (j1.eq.1) then                                                 
         if (iv.eq.2) iv=3                                              
         if (iv.eq.4) go to 21                                          
      end if                                                            
c                                                                       
      k=0                                                               
c                                                                       
c                                                                       
c        skip line                                                      
c                                                                       
      do i=1,nskip                                                      
      read (kda8,10008) name,nname                                      
c**** write(kwrite,10008) name,nname                                    
      enddo                                                             
c                                                                       
c                                                                       
   23 k=k+1                                                             
c                                                                       
c                                                                       
      read (kda8,10200) psa(1,k,iv,j1,1),                               
     2                  psa(2,k,iv,j1,1),                               
     3                  psa(3,k,iv,j1,1)                                
c                                                                       
c**** write(kwrite,10200) psa(1,k,iv,j1,1),                             
c    2                  psa(2,k,iv,j1,1),                               
c    3                  psa(3,k,iv,j1,1)                                
c                                                                       
c                                                                       
c                                                                       
      if (psa(1,k,iv,j1,1).ne.0.e0) go to 23                            
      mpsa(iv,j1,1)=k-1                                                 
c                                                                       
      if (iv.lt.5) go to 22                                             
      if (j1.lt.7) go to 21                                             
c                                                                       
c        this has been the end of reading ips=1                         
c        --------------------------------------                         
c
 8888 continue
c                                                                       
c                                                                       
c        read and write input parameters for this program               
c                                                                       
c                                                                       
      write (kwrite,10004)                                              
      read  (kread, 10008) name,nname                                   
      write (kwrite,10009) name,nname                                   
c                                                                       
      read  (kread ,10000) name,jb,je                                   
      write (kwrite,10001) name,jb,je                                   
      read  (kread ,10000) name,jborn                                   
      write (kwrite,10001) name,jborn                                   
      jb1=jb+1                                                          
      je1=je+1                                                          
      jee1=je1                                                          
      if (jee1.gt.20) jee1=20                                           
      read  (kread ,10000) name,(nj(j1),j1=1,jee1)                      
      write (kwrite,10001) name,(nj(j1),j1=1,jee1)                      
      read  (kread ,10000) name,ising,itrip,icoup                       
      write (kwrite,10001) name,ising,itrip,icoup                       
c        iprop=1: non-relativistic propagator in scattering equation;   
c        iprop=2: relativistic propagator.                              
      read  (kread ,10000) name,iprop                                   
      write (kwrite,10001) name,iprop                                   
c        iphrel=1: non-relativistic phase-relation;                     
c        iphrel=2: relativistic phase-relation.                         
      read  (kread ,10000) name,iphrel                                  
      write (kwrite,10001) name,iphrel                                  
c        ipotfc is presently a dead parameter                           
      read  (kread ,10000) name,ipotfc                                  
      write (kwrite,10001) name,ipotfc                                  
      read  (kread ,10072) name,c,cut                                   
      write (kwrite,10073) name,c,cut                                   
c        wn1 is projectile mass, wn2 is target mass;                    
c        i.e., in neutron-proton scattering, wn1 is the neutron mass.   
      read  (kread ,10002) name,wn1,wn2                                 
      write (kwrite,10003) name,wn1,wn2                                 
      write (kwrite,10003)                                              
      do 1 k=1,memax                                                    
      read  (kread ,10002) name,elab(k)                                 
      write (kwrite,10014) name,elab(k)                                 
      if (elab(k).eq.0.d0) go to 2                                      
    1 continue                                                          
    2 melab=k-1                                                         
c        mps(ips).ne.0: phase-shift analysis ips is to be considered    
      read  (kread ,10000) name,mps                                     
      write (kwrite,10001) name,mps                                     
c        icou=0: no coulomb corrections                                 
c        icou=1: coulomb corrections                                    
c        irel=1: relativistic coulomb potential (breit, pr 99, 1581     
c                (1955); austen + deswart, prl 50, 2039 (1983);         
c                bergervoet et al., prc 38, 15 (1988)).                 
      read  (kread ,10000) name,icou,irel                               
      write (kwrite,10001) name,icou,irel                               
      read  (kread ,10000) name,ncoul                                   
      write (kwrite,10001) name,ncoul                                   
      read  (kread ,10002) name,rcoul                                   
      write (kwrite,10003) name,rcoul                                   
c        irma=1: the r-matrix is punched as vector in lsj-states;       
c        irma=2: the r-matrix is punched in short form,                 
c                for bremsstrahlungs-calculations;                      
c        iqua=1: a quadratic r-matrix is evaluated and punched          
c                       in lsj-states                                   
c        note: the matrix b has to be as large as the matrix a          
c              if iqua=1.                                               
      read  (kread ,10000) name,irma,iqua                               
      write (kwrite,10001) name,irma,iqua                               
c        ipoint=0: the transformed gauss-points and -weights are        
c                  not printed;                                         
c        ipoint.ne.0: ... are printed.                                  
      read  (kread ,10000) name,ipoint                                  
      write (kwrite,10001) name,ipoint                                  
c        iwrite=1: the phase-shifts are printed in a compact way        
c                  (in degrees or radians, no chi-square);              
c        iwrite=2: the phase-shifts are printed in degrees and with     
c                  chi-square;                                          
c        iwrite=3: the phase-shifts are printed in degrees and in       
c                  radians and with chi-square.                         
      read  (kread ,10000) name,iwrite                                  
      write (kwrite,10001) name,iwrite                                  
c                                                                       
c        ierrma=0: the error matrix is not used to calculate the chi**2;
c        ierrma=1: the pp-only error matrix is used to calculate the chi
c                  a file phpperr.d is produced;                        
c        ierrma=2: the pp+np error matrix is  used to calculate the chi*
c                  the file phpperr.d must exist.                       
c                                                                       
c        ipunch=0: no punching of phase-shifts;                         
c        ipunch=1: punching for computation of observables              
c                  (5 states in one data-card);                         
c        ipunch=2 - 4: punching for the plot-routine:                   
c        ipunch=2: short punch of predicted phase shifts.               
c                  (i.e. no extra cards for state and labels of axes);  
c        ipunch=3: long punch with experimental data, but without       
c                  predicted phase shifts.                              
c        ipunch=4: long punch without experimental data and without     
c                  predicted phase shifts.                              
      read  (kread ,10000) name,ierrma,ipunch                           
      write (kwrite,10001) name,ierrma,ipunch                           
c        ideg=0: phase-shifts are printed in radians, stapp convention; 
c        ideg.ne.0: phase-shifts are printed in degrees:                
c                   ideg=1: blatt-bidenharn convention,                 
c                   ideg=2: stapp convention.                           
      read  (kread ,10000) name,ideg                                    
      write (kwrite,10001) name,ideg                                    
c                                                                       
c                                                                       
c        read-in error matrix                                           
c        --------------------                                           
c                                                                       
      if (ierrma.eq.1.or.ierrma.eq.2) indema=.true.                     
      if (indema) then                                                  
c                                                                       
c        check if indema makes sense                                    
      if (jb.ne.0.or.je.lt.4) then                                      
      write (kwrite,11006)                                              
      indema=.false.                                                    
      go to 150                                                         
      end if                                                            
c                                                                       
      kdaa=kda(ierrma)                                                  
      read (kdaa,10008) name,nname                                      
c**** write(kwrite,10008) name,nname                                    
      read (kdaa,11005) (phdel(2,kk),kk=1,5)                            
c**** write(kwrite,11005) (phdel(2,kk),kk=1,5)                          
      read (kdaa,11000) melerr                                          
c**** write(kwrite,11000) melerr                                        
      do 125 kk=1,melerr                                                
      read (kdaa,10008) name,nname                                      
c**** write(kwrite,10008) name,nname                                    
      read (kdaa,11001) elaber(kk),npherr(kk)                           
c**** write(kwrite,11001) elaber(kk),npherr(kk)                         
      npherk=npherr(kk)                                                 
      do 105 i=1,npherk                                                 
      go to (101,102),ierrma                                            
  101 read (kdaa,11002) namer1(i),pherr(i,kk)                           
c**** write(kwrite,11002) namer1(i),pherr(i,kk)                         
      go to 105                                                         
  102 read (kdaa,11022) namer2(i,kk),pherr(i,kk)                        
c**** write(kwrite,11022) namer2(i,kk),pherr(i,kk)                      
  105 continue                                                          
      do 115 i=1,npherk                                                 
      read (kdaa,11003) (errma(i,ii,kk),ii=1,i)                         
c**** write(kwrite,11003) (errma(i,ii,kk),ii=1,i)                       
  115 continue                                                          
  125 continue                                                          
c                                                                       
c                                                                       
      kda3=kda(3)                                                       
      if (ierrma.eq.2) then                                             
      read (kda3,11100) delj                                            
      end if                                                            
c                                                                       
c                                                                       
      end if                                                            
c                                                                       
c                                                                       
  150 continue                                                          
c                                                                       
c                                                                       
c        preparation of constants                                       
c                                                                       
c                                                                       
      sssing=.true.                                                     
      tttrip=.true.                                                     
      cccoup=.true.                                                     
      sing=.false.                                                      
      trip=.false.                                                      
      coup=.false.                                                      
      if (ising.ne.0) sing=.true.                                       
      if (itrip.ne.0) trip=.true.                                       
      if (icoup.ne.0) coup=.true.                                       
      ssing=sing                                                        
      ttrip=trip                                                        
      ccoup=coup                                                        
      heform=.true.                                                     
      if (irma.ne.0) indrma=.true.                                      
      if (iqua.ne.0) indqua=.true.                                      
      if (ipoint.ne.0) indpts=.true.                                    
      indcou=.false.                                                    
      if (icou.ne.0) indcou=.true.                                      
      indrel=.false.                                                    
      if (irel.ne.0) indrel=.true.                                      
      if (ipunch.eq.1) indp1=.true.                                     
      if (ipunch.ge.2) indpch=.true.                                    
      if (ipunch.eq.2) indshp=.true.                                    
      if (ipunch.eq.3) inderr=.true.                                    
      if (ideg.ne.0) inddeg=.true.                                      
      do 4 ips=1,nps                                                    
      if (mps(ips).ne.0) indps(ips)=.true.                              
    4 continue                                                          
c                                                                       
      if(jb.gt.1.or.elab(1).gt.5.d0.or.elab(2).gt.5.d0.or.melab.lt.3)   
     1 go to 3                                                          
      inderg=.true.                                                     
    3 lppg=1000                                                         
      line=lppg                                                         
c                                                                       
c                                                                       
c        calculate reduced mass from the two given nucleon masses       
c        and use twice reduced mass as nucleon mass in this code        
c                                                                       
      wn=2.d0*wn1*wn2/(wn1+wn2)                                         
c                                                                       
      wnh=wn*0.5d0                                                      
      wnq=wn*wn                                                         
      wp=wn*pih                                                         
      chitot=0.d0                                                       
      ichi=0                                                            
      rd=90.d0/pih                                                      
      dr=1.d0/rd                                                        
c                                                                       
c                                                                       
      indnnb=.false.                                                    
      gamma0= wn/2.d0/alfinv                                            
      gamma=0.d0                                                        
      if(indcou) then                                                   
         write(kwrite,*) '    '                                         
         write(kwrite,*) ' coulomb potential in t=1 partial waves '     
         write(kwrite,*) ' rcoul = ',rcoul,' ncoul = ',ncoul            
      endif                                                             
      do 6 i1=1,2                                                       
      do 6 i2=1,2                                                       
      fc( i1,i2)=0.d0                                                   
      fc0(i1,i2)=0.d0                                                   
      gc (i1,i2)=0.d0                                                   
      gc0(i1,i2)=0.d0                                                   
      fd (i1,i2)=0.d0                                                   
      fb0(i1,i2)=0.d0                                                   
      gd (i1,i2)=0.d0                                                   
 6    gb0(i1,i2)=0.d0                                                   
c                                                                       
c                                                                       
      uanspp=0.d0                                                       
      wsnspp=wn                                                         
      ucnspp=0.d0                                                       
      udnspp=0.d0                                                       
      qfmev=0.d0                                                        
      pmev=0.d0                                                         
c                                                                       
c                                                                       
c        prepare energies                                               
c                                                                       
      do 10 k=1,melab                                                   
      q0q(k)=elab(k)*wn2**2*(elab(k)+2.d0*wn1)                          
     1       /((wn1+wn2)**2+2.d0*wn2*elab(k))                           
      q0(k)=dsqrt(q0q(k))                                               
      eq0(k)=dsqrt(q0q(k)+wnq)                                          
   10 arel(k)=(wnq+2.*q0q(k))/(wn*eq0(k))                               
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c        loop of total angular momentum j                               
c        --------------------------------                               
c        --------------------------------                               
c                                                                       
c                                                                       
c                                                                       
c                                                                       
      do 2000 j1=jb1,je1                                                
c                                                                       
c                                                                       
      indj=.false.                                                      
      j2=j1+1                                                           
      j=j1-1                                                            
      aj=dfloat(j)                                                      
      aj1=dfloat(j+1)                                                   
      a2j1=dfloat(2*j+1)                                                
      d2j1=1.d0/a2j1                                                    
      arjj1=dsqrt(aj*aj1)                                               
      aaj=arjj1                                                         
c                                                                       
c                                                                       
      if (indcou) then                                                  
          if (mod(j,2).eq.0) then                                       
             sing=ssing                                                 
             trip=.false.                                               
             coup=ccoup                                                 
             delta(2)=0.                                                
             rb1(2)=0.                                                  
             rb2(2)=0.                                                  
          else                                                          
             sing=.false.                                               
             trip=ttrip                                                 
             coup=.false.                                               
             delta(1)=0.                                                
             delta(3)=0.                                                
             delta(4)=0.                                                
             delta(5)=0.                                                
             rb1(1)=0.                                                  
             rb2(1)=0.                                                  
             rb1(3)=0.                                                  
             rb2(3)=0.                                                  
          endif                                                         
      endif                                                             
c                                                                       
c                                                                       
      if (j.ge.jborn) indbrn=.true.                                     
c                                                                       
c                                                                       
c                                                                       
c                                                                       
      if (j1.gt.20) go to 300                                           
c        number of gausspoints for this j                               
      if (nj(j1).eq.0) nj(j1)=nj(j)                                     
      if (nj(j1).eq.0) nj(j1)=nj(1)                                     
      n=nj(j1)                                                          
      if (n.eq.nalt) go to 300                                          
c                                                                       
c                                                                       
c        get gauss points and weights                                   
c                                                                       
      if (cut.eq.0.d0) then                                             
      acut=1.d0                                                         
      else                                                              
      acut=cut                                                          
      endif                                                             
c**** call gset (0.d0,acut,n,u,s)                                       
      call gauss(0.d0,acut,n,u,s)                                       
c                                                                       
      nalt=n                                                            
      n1=n+1                                                            
      n2=2*n1                                                           
      n3=3*n1                                                           
      n4=4*n1                                                           
      nx=n1                                                             
      if (indqua) nx=nx*nx                                              
      nx2=2*nx                                                          
      nx2mn=nx2-n1                                                      
      nx2pn=nx2+n1                                                      
      nx4=4*nx                                                          
      nx4mn=nx4-n1                                                      
c                                                                       
       iiienn=n1*(n/2+1)*melab                                          
c                                                                       
c        transform gauss points and weights                             
c                                                                       
      do 201 i=1,n                                                      
c                                                                       
      if (cut.eq.0.d0) then                                             
      xx=pih*u(i)                                                       
c        transformed gauss point                                        
      q(i)=dtan(xx)*c                                                   
c        transformed gauss weight                                       
      dc=1.d0/dcos(xx)                                                  
      s(i)=pih*c*dc*dc*s(i)                                             
      else                                                              
      q(i)=u(i)                                                         
      end if                                                            
c                                                                       
      qq(i)=q(i)*q(i)                                                   
      eq(i)=dsqrt(qq(i)+wnq)                                            
c                                                                       
  201 continue                                                          
c                                                                       
c                                                                       
      if (.not.indpts) go to 250                                        
c                                                                       
c        write gauss points and weights                                 
c                                                                       
      if (j.eq.jb) go to 204                                            
      if (line.lt.lppg-5-n/4) go to 202                                 
      call headln (iwrite)                                              
      line=0                                                            
  202 line=line+5+n/4                                                   
  204 write (kwrite,10005) c,cut,n                                      
      write (kwrite,10006) (q(i),i=1,n)                                 
c                                                                       
c                                                                       
c        punch gausspoints                                              
c                                                                       
c                                                                       
      if (irma.eq.0) go to 206                                          
c                                                                       
      write (kpunch,10111) n                                            
      write (kpunch,10113) (q(i),i=1,n)                                 
c                                                                       
c                                                                       
  206 if (j.eq.jb) go to 205                                            
      if (line.lt.lppg-2-n/4) go to 203                                 
      call headln (iwrite)                                              
      line=0                                                            
  203 line=line+2+n/4                                                   
  205 write (kwrite,10007)                                              
      write (kwrite,10006) (s(i),i=1,n)                                 
c                                                                       
c                                                                       
  250 continue                                                          
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c        loop of elabs                                                  
c        -------------                                                  
c        -------------                                                  
c                                                                       
c                                                                       
c                                                                       
c                                                                       
  300 do 1000 k=1,melab                                                 
c                                                                       
c                                                                       
c        define starting energy                                         
c                                                                       
      noced=.false.                                                     
      zrel=2.d0*eq0(k)                                                  
      znrl=q0q(k)/wn                                                    
      smev=2.d0*eq0(k)                                                  
      q0qmev=q0q(k)                                                     
      q(n1)=q0(k)                                                       
c                                                                       
c                                                                       
      if (indrma) write (kpunch,10112) elab(k)                          
c                                                                       
c                                                                       
      if (indbrn.and..not.indqua) go to 500                             
c                                                                       
c                                                                       
c        check if right potential matrix does already exist             
      if (indj.and..not.endep) go to 450                                
      indj=.true.                                                       
c                                                                       
c                                                                       
c        compute potential matrix                                       
c        ------------------------                                       
c                                                                       
c                                                                       
      iii=0                                                             
      do 405 ix=1,n                                                     
c                                                                       
c                                                                       
      xmev=q(ix)                                                        
c                                                                       
c                                                                       
      do 405 iy=ix,n                                                    
c                                                                       
c                                                                       
      ymev=q(iy)                                                        
c                                                                       
c                                                                       
      call pot                                                          
c                                                                       
c                                                                       
      iaa=iii*6                                                         
      iii=iii+1                                                         
      do 405 iv=1,6                                                     
  405 aa(iv+iaa)=v(iv)                                                  
c                                                                       
c                                                                       
      if (indcou)  then                                                 
      iii=0                                                             
      do 415 ix=1,n                                                     
      xmev=q(ix)                                                        
      do 415 iy=ix,n                                                    
      ymev=q(iy)                                                        
c                                                                       
      call coul(1)                                                      
c                                                                       
      iaa=iii*6                                                         
      iii=iii+1                                                         
      do 415 iv=1,6                                                     
      aac(iv+iaa)=v(iv)                                                 
  415 aan(iv+iaa)=aa(iv+iaa)                                            
      end if                                                            
c                                                                       
c                                                                       
  450 if (.not.indcou) go to 500                                        
      facrel=1.d0                                                       
      if (indrel) facrel=(wnq+2.d0*q0q(k))/(wn*eq0(k))                  
      iii=0                                                             
      do 455 ix=1,n                                                     
      do 455 iy=ix,n                                                    
      iaa=iii*6                                                         
      iii=iii+1                                                         
      do 455 iv=1,6                                                     
  455 aa(iv+iaa)=aan(iv+iaa)+aac(iv+iaa)*facrel                         
c                                                                       
c                                                                       
c                                                                       
c        compute potential vector                                       
c        ------------------------                                       
c                                                                       
c                                                                       
  500 if (indbrn.and..not.indrma) go to 510                             
c                                                                       
      xmev=q0(k)                                                        
      ix=n1                                                             
c                                                                       
c                                                                       
      do 501 iy=1,n                                                     
c                                                                       
c                                                                       
      ymev=q(iy)                                                        
c                                                                       
c                                                                       
      call pot                                                          
c        exchange v5 and v6, as q0 and q were exchanged                 
      v5=v(5)                                                           
      v(5)=v(6)                                                         
      v(6)=v5                                                           
c                                                                       
c                                                                       
      if (indcou)  call coul(0)                                         
c                                                                       
      do 501 iv=1,6                                                     
      ivv=iy+(iv-1)*n                                                   
  501 vv(ivv)=v(iv)                                                     
c                                                                       
c                                                                       
c        compute potential element                                      
c        -------------------------                                      
c                                                                       
c                                                                       
  510 xmev=q0(k)                                                        
      ymev=q0(k)                                                        
      ix=n1                                                             
      iy=n1                                                             
c                                                                       
c                                                                       
      call pot                                                          
c                                                                       
c                                                                       
      if (indcou)  call coul(0)                                         
c                                                                       
c                                                                       
c                                                                       
c        compute factor for the phase relation requested                
c                                                                       
      go to (601,602,603),iphrel   
  601 wpq0=-wp*q0(k)                  
      go to 610                     
  602 wpq0=-pih*eq0(k)*q0(k)       
      go to 610                
  603 wpq0=-pih*q0(k)*wnq/eq0(k)       
  610 continue
c                                                                       
c                                                                       
      if (indbrn.and..not.indqua) go to 700                             
c                                                                       
c                                                                       
c        compute propagator                                             
c        ------------------                                             
c                                                                       
c                                                                       
      uq0=0.d0                                                          
      do 620 i=1,n                                                      
      sdq=s(i)/(qq(i)-q0q(k))                                           
c                                                                       
c        respect typ of propagator requested                            
c                                                                       
      go to (621,622),iprop                                             
  621 u(i)=sdq*qq(i)*wn                                                 
      go to 620                                                         
  622 u(i)=sdq*qq(i)*(eq(i)+eq0(k))*0.5d0                               
c                                                                       
  620 uq0=uq0+sdq                                                       
c                                                                       
c                                                                       
      uq0=uq0+dlog(dabs((cut+q0(k))/(cut-q0(k))))/(2.d0*q0(k))          
c                                                                       
c                                                                       
      go to (631,632),iprop                                             
  631 uq0=-uq0*q0q(k)*wn                                                
      go to 700                                                         
  632 uq0=-uq0*q0q(k)*eq0(k)                                            
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c        build up matrix to be inverted                                 
c        ------------------------------                                 
c                                                                       
c                                                                       
  700 ni=0                                                              
      nii=0                                                             
      nv=n1                                                             
      mv=1                                                              
      if (indqua) mv=mv*n1                                              
      ib=0                                                              
      eins=1.d0                                                         
c                                                                       
c                                                                       
      if (.not.sing) go to 720                                          
      iv=1                                                              
      go to 770                                                         
c                                                                       
c                                                                       
  720 if (.not.trip.or.j.eq.0) go to 730                                
      iv=2                                                              
      go to 770                                                         
c                                                                       
c                                                                       
  730 if (.not.coup) go to 900                                          
      iv=3                                                              
      if (j.eq.0) go to 770                                             
      nv=n2                                                             
      mv=2                                                              
      if (indqua) mv=mv*n1                                              
      go to 770                                                         
c                                                                       
  740 if (j.eq.0) go to 800                                             
      iv=4                                                              
      ib=n3                                                             
      ni=n1                                                             
      nii=n1                                                            
      go to 770                                                         
c                                                                       
  750 iv=5                                                              
      ivi=6                                                             
      ib=n2                                                             
      ni=0                                                              
      nii=n1                                                            
      eins=0.d0                                                         
      go to 770                                                         
c                                                                       
  760 iv=6                                                              
      ivi=5                                                             
      ib=n1                                                             
      ni=n1                                                             
      nii=0                                                             
c                                                                       
c                                                                       
c                                                                       
c                                                                       
  770 iii=0                                                             
      if (iv.le.4) ivi=iv                                               
      igg=(iv-1)*n                                                      
      i1=(nii+n)*nv                                                     
      i2=(nii-1)*nv                                                     
c                                                                       
c                                                                       
      if (indbrn.and..not.indrma) go to 785                             
c                                                                       
c                                                                       
      do 780 i=1,n                                                      
      i3=i*nv                                                           
      i4=ni+i                                                           
c                                                                       
c                                                                       
      do 781 ii=i,n                                                     
      iaa=iii*6                                                         
      iii=iii+1                                                         
      i5=i2+i3+ni+ii                                                    
      i6=ivi+iaa                                                        
      i7=i2+i4+ii*nv                                                    
      i8=iv+iaa                                                         
      if (i.eq.ii) go to 782                                            
c                                                                       
c        matrix a                                                       
      a(i7)=aa(i8)*u(ii)                                                
      a(i5)=aa(i6)*u(i)                                                 
      if (.not.indqua) go to 781                                        
c                                                                       
c        matrix b                                                       
      b(i7)=aa(i8)                                                      
      b(i5)=aa(i6)                                                      
      go to 781                                                         
c        diagonal element                                               
  782 a(i7)=aa(i8)*u(i)+eins                                            
      if (.not.indqua) go to 781                                        
      b(i7)=aa(i8)                                                      
  781 continue                                                          
c                                                                       
c        last column                                                    
      i9=i1+i4                                                          
      i10=i+igg                                                         
      a(i9)=vv(i10)*uq0                                                 
c        last row                                                       
      i11=i2+i3+ni+n1                                                   
      ivv=i+(ivi-1)*n                                                   
      a(i11)=vv(ivv)*u(i)                                               
      if (.not.indqua) go to 783                                        
      b(i9)=vv(i10)                                                     
      b(i11)=vv(ivv)                                                    
      go to 780                                                         
c                                                                       
c        vector b                                                       
  783 b(ib+i)=vv(i+igg)                                                 
c                                                                       
  780 continue                                                          
c                                                                       
c                                                                       
c        last element                                                   
      i12=i1+ni+n1                                                      
      a(i12)=v(iv)*uq0+eins                                             
      if (.not.indqua) go to 785                                        
      b(i12)=v(iv)                                                      
      go to 790                                                         
  785 b(ib+n1)=v(iv)                                                    
c                                                                       
c                                                                       
c                                                                       
c                                                                       
  790 go to (800,800,740,750,760,800),iv                                
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c        invert matrix                                                  
c        -------------                                                  
c                                                                       
c                                                                       
c                                                                       
c                                                                       
  800 if (indbrn) go to 801                                             
      call dgelg (b,a,nv,mv,ops,ier)                                    
c                                                                       
c                                                                       
      if (ier.lt.0) write(kwrite,10013) ier                             
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c        compute phase shifts                                           
c        --------------------                                           
c                                                                       
c                                                                       
c                                                                       
c                                                                       
  801 if (iv.gt.2.and.j.ne.0) go to 820                                 
c                                                                       
c                                                                       
c        uncoupled cases                                                
c                                                                       
      r0=b(nx)*wpq0                                                     
c                                                                       
      ivj=iabs(iv-j)                                                    
c                                                                       
c    match asymptotic coulomb wave functions (for t=1 states)           
c                                                                       
      if(indcou.and.mod(ivj,2).eq.1) then                               
        rho=q0(k)*rcoul/uf                                              
        gamma=gamma0/q0(k)                                              
      if (indrel) gamma=gamma*arel(k)                                   
        xln=dfloat(j)                                                   
        if(iv.gt.2) xln=dfloat(j+1)                                     
        call klein(rho,0.d0,xln,fp,gp,fdp,gdp,ifail)                    
        call klein(rho,gamma,xln,fp1,gp1,fdp1,gdp1,ifail)               
        xx=r0*(gp*fdp1-fp1*gdp)+(fp*fdp1-fp1*fdp)                       
        r0=xx/(r0*(gp1*gdp-gp*gdp1)+(gp1*fdp-fp*gdp1))                  
      endif                                                             
c                                                                       
        delta(iv)=datan(r0)                                             
c                                                                       
      if (.not.indrma) go to 809                                        
c                                                                       
c                                                                       
c                                                                       
c        punch r-matrix                                                 
c                                                                       
c                                                                       
      if (irma.ne.2) go to 803                                          
c                                                                       
c                                                                       
c        in short form for bremsstrahlungs-calculations                 
c                                                                       
c                                                                       
      write (kpunch,10113) (b(i),i=1,n1)                                
      go to 809                                                         
c                                                                       
c                                                                       
  803 state(1)=multi(iv)                                                
      state(2)=spd(j1+ldel(iv))                                         
      if (j.eq.0.and.iv.eq.3) state(2)=spd(2)                           
      state3=state(2)                                                   
      if (indqua) write (kpunch,10110) label,state,state3,j             
c        write on-shell element                                         
      if (indqua) write (kpunch,10113) b(nx)                            
      if (indqua) write (kpunch,10114)                                  
      do 805 i=1,n1                                                     
      if (indqua) go to 804                                             
c                                                                       
c        punch r-matrix as vector                                       
      write (kpunch,10110) label,state,state3,j,q(i),q(n1),b(i)         
      go to 805                                                         
c                                                                       
c        punch quadratic r-matrix                                       
  804 i1=(i-1)*n1                                                       
      write (kpunch,10113) (b(i1+ii),ii=1,n1)                           
      write (kpunch,10114)                                              
  805 continue                                                          
c                                                                       
c                                                                       
c    calculation of low energy parameter (if requested)                 
c                                                                       
 809  if(inderg.and.mv.eq.1.and.j.le.1.and.k.le.3) then                 
        im=mod(k,3)+1                                                   
        if(indcou) then                                                 
          cf=4.d0*pih*gamma/(dexp(4.d0*pih*gamma)-1.d0)                 
          hgamma=hfunc(gamma)                                           
        else                                                            
          cf=1.d0                                                       
          hgamma=0.d0                                                   
        endif                                                           
        ale(im,iv)=q0(k)*q0(k)/uf/uf/2.d0                               
        if(j.eq.1.or.iv.ne.1) then                                      
          ble(im,iv)=(q0(k)/uf)**2*(1.d0+gamma**2)                      
     1              *(cf*q0(k)/r0+2.d0*q0(k)*gamma*hgamma)/uf           
        else                                                            
          ble(im,iv)=(cf*q0(k)/r0+2.d0*q0(k)*gamma*hgamma)/uf           
        endif                                                           
        cle(im,iv)=(q0(k)/uf)**4                                        
        if(im.eq.1) then                                                
          call equ(ale(1,iv),ble(1,iv),cle(1,iv),rb2(iv),xx)            
          rb1(iv)=-1.d0/xx                                              
        endif                                                           
      endif                                                             
c                                                                       
c                                                                       
c                                                                       
      go to (720,730,900),iv                                            
c                                                                       
c                                                                       
c                                                                       
c        coupled cases                                                  
c                                                                       
  820 continue                                                          
c                                                                       
c    conversion to lsj - basis                                          
c                                                                       
c**** if(heform) then                                                   
         r(3)=b(nx2mn)                                                  
         r(4)=b(nx4)                                                    
         r(5)=b(nx4mn)                                                  
         r(6)=b(nx2)                                                    
         r34=(r(3)-r(4))*aaj                                            
         r56=(r(5)+r(6))*aaj                                            
         r2=(aj1*r(3)+aj*r(4)-r56)*d2j1                                 
         r0=(aj*r(3)+aj1*r(4)+r56)*d2j1                                 
         r1=(r34+aj1*r(5)-aj*r(6))*d2j1                                 
c**** else                                                              
c****    r0=b(nx4)                                                      
c****    r1=(b(nx2)+b(nx2))/2.d0                                        
c****    r2=b(nx2mn)                                                    
c**** endif                                                             
c                                                                       
c    matches to asymptotic coulomb wave functions (for t=1 states)      
c                                                                       
      if(indcou.and.mod(j,2).eq.0) then                                 
        rho=q0(k)*rcoul/uf                                              
        gamma=gamma0/q0(k)                                              
      if (indrel) gamma=gamma*arel(k)                                   
        dss(1,1)=wpq0*r0                                                
        dss(1,2)=wpq0*r1                                                
        dss(2,1)=dss(1,2)                                               
        dss(2,2)=wpq0*r2                                                
        xln=dfloat(j-1)                                                 
        call klein(rho,0.d0,xln,fc0(1,1),gc0(1,1),fb0(1,1),gb0(1,1)     
     1  ,ifail)                                                         
        call klein(rho,gamma,xln,fc(1,1),gc(1,1),fd(1,1),gd(1,1),ifail) 
        xln=dfloat(j+1)                                                 
        call klein(rho,0.d0,xln,fc0(2,2),gc0(2,2),fb0(2,2),gb0(2,2)     
     1  ,ifail)                                                         
        call klein(rho,gamma,xln,fc(2,2),gc(2,2),fd(2,2),gd(2,2),ifail) 
        call mult(gc0,dss,hm1)                                          
        call mult(gb0,dss,hm2)                                          
        do 1515 ij=1,2                                                  
        do 1515 jj=1,2                                                  
        hm1(ij,jj)=fc0(ij,jj)+hm1(ij,jj)                                
1515    hm2(ij,jj)=fb0(ij,jj)+hm2(ij,jj)                                
        call inv(hm2,hm3)                                               
        call mult(hm1,hm3,hm2)                                          
        call mult(hm2,fd,hm1)                                           
        call mult(hm2,gd,hm3)                                           
        do 1516 ij=1,2                                                  
        do 1516 jj=1,2                                                  
        hm1(ij,jj)=hm1(ij,jj)-fc(ij,jj)                                 
1516    hm3(ij,jj)=gc(ij,jj)-hm3(ij,jj)                                 
        call inv(hm3,hm2)                                               
        call mult(hm2,hm1,dss)                                          
        r0=dss(1,1)/wpq0                                                
        r1=dss(1,2)/wpq0                                                
        r2=dss(2,2)/wpq0                                                
      endif                                                             
c                                                                       
c                                                                       
      rr=2.d0*r1/(r0-r2)                                                
      delta(5)=datan(rr)/2.d0                                           
      rr=(r0-r2)*dsqrt(1.d0+rr*rr)                                      
      if (inderg.and.j.eq.1.and.k.le.2)                                 
     1 rb(k)=-2.d0*q0(k)/(wpq0*(r0+r2+rr))                              
      delta(3)=datan((r0+r2+rr)*wpq0*0.5d0)                             
      delta(4)=datan((r0+r2-rr)*wpq0*0.5d0)                             
      if (ideg.eq.1) go to 830                                          
c        so far the delta(..) have been the blatt-biedenharn phase-     
c        shifts, transform the delta(..) now into bar-phase-shifts      
c        according to stapp et al. phys.rev. 105 (1957) 302             
      pp=delta(3)+delta(4)                                              
      pm=delta(3)-delta(4)                                              
      d52=2.d0*delta(5)                                                 
      pm1=dtan(pm)*dcos(d52)                                            
      pm1=datan(pm1)                                                    
      pm2=dtan(d52)*dsin(pm1)                                           
      delta(5)=0.5d0*datan(pm2)                                         
      delta(3)=0.5d0*(pp+pm1)                                           
      delta(4)=0.5d0*(pp-pm1)                                           
      if (j.eq.1.and.elab(k).lt.50..and.delta(3).lt.0.) then            
      delta(3)=delta(3)+2.*pih                                          
      delta(5)=-delta(5)                                                
      end if                                                            
  830 continue                                                          
c                                                                       
      if (.not.indrma) go to 900                                        
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c        transform the r-matrix into lsj-form                           
c                                                                       
      na=1                                                              
      if (indqua) na=n1                                                 
      do 835 i=1,na                                                     
      i1=(i-1)*n2                                                       
      do 835 ii=i,n1                                                    
      i3=i1+ii                                                          
      i4=nx2pn+i3                                                       
      i5=nx2+i3                                                         
      i6=n1+i3                                                          
      r(3)=b(i3)                                                        
      r(4)=b(i4)                                                        
      r(5)=b(i5)                                                        
      r(6)=b(i6)                                                        
      r34=(r(3)-r(4))*aaj                                               
      r56=(r(5)+r(6))*aaj                                               
      b(i6)=(aj1*r(3)+aj*r(4)-r56)*d2j1                                 
      b(i3)=(aj*r(3)+aj1*r(4)+r56)*d2j1                                 
      b(i5)=(r34+aj1*r(5)-aj*r(6))*d2j1                                 
      b(i4)=(r34-aj*r(5)+aj1*r(6))*d2j1                                 
  835 continue                                                          
c                                                                       
c        punch r-matrix                                                 
      ivx=0                                                             
      do 840 iv1=3,4                                                    
      do 840 iv2=3,4                                                    
      ivx=ivx+1                                                         
      state(1)=multi(iv1)                                               
      state(2)=spd(j1+ldel(iv1))                                        
      state3=spd(j1+ldel(iv2))                                          
      go to (831,832,833,834),ivx                                       
  831 ny=0                                                              
      go to 836                                                         
  832 ny=nx2pn                                                          
      go to 836                                                         
  833 ny=nx2                                                            
      go to 836                                                         
  834 ny=n1                                                             
  836 if (irma.ne.2) go to 837                                          
c                                                                       
c        punch r-matrix in short form,                                  
c        for bremsstrahlungs-calculations                               
c                                                                       
      ny1=ny+1                                                          
      nyn1=ny+n1                                                        
      write (kpunch,10113) (b(i),i=ny1,nyn1)                            
      go to 840                                                         
c                                                                       
c                                                                       
  837 if (indqua) write (kpunch,10110) label,state,state3,j             
c        write on-shell element                                         
      if (indqua) write (kpunch,10113) b(n*n2+ny+n1)                    
      if (indqua) write (kpunch,10114)                                  
      do 839 i=1,n1                                                     
      if (indqua) go to 838                                             
c                                                                       
c        punch r-matrix as vector                                       
      ii1=ny+i                                                          
      write (kpunch,10110) label,state,state3,j,q(i),q(n1),b(ii1)       
      go to 839                                                         
c                                                                       
c        punch quadratic r-matrix                                       
  838 i1=(i-1)*n2+ny                                                    
      write (kpunch,10113) (b(i1+ii),ii=1,n1)                           
      write (kpunch,10114)                                              
  839 continue                                                          
  840 continue                                                          
c                                                                       
c                                                                       
c                                                                       
c                                                                       
  900 continue                                                          
c                                                                       
c                                                                       
c        store phase shifts for application to error matrix             
c                                                                       
      if (.not.indema) go to 960                                        
      if (j1.gt.5) go to 960                                            
      go to (951,952),ierrma                                            
  951 do 955 iv=1,5                                                     
  955 delj(iv,j1,k)=delta(iv)*rd                                        
      go to 960                                                         
  952 if (mod(j,2).eq.0) then                                           
      if (j.eq.0) then                                                  
      delj(2,j1,k)=delj(1,j1,k)                                         
      delj(1,j1,k)=delta(1)*rd                                          
      else                                                              
      delj(2,j1,k)=delta(2)*rd                                          
      end if                                                            
      else                                                              
      delj(1,j1,k)=delta(1)*rd                                          
      delj(3,j1,k)=delta(3)*rd                                          
      delj(4,j1,k)=delta(4)*rd                                          
      delj(5,j1,k)=delta(5)*rd                                          
c**** if (j.eq.1.and.elab(k).lt.25.d0.and.delj(3,j1,k).lt.0.d0) then    
c**** delj(3,j1,k)=delj(3,j1,k)+180.d0                                  
c**** delj(5,j1,k)=-delj(5,j1,k)                                        
c**** end if                                                            
      end if                                                            
c                                                                       
c                                                                       
  960 if (.not.indp1) go to 930                                         
c                                                                       
c                                                                       
c        punch phase-shifts in short form (e.g. for computation of      
c        observables)                                                   
c                                                                       
      if (indpl) go to 931                                              
      write (kpunch,10040) label                                        
      indpl=.true.                                                      
  931 write (kpunch,10041) elab(k),j,delta                              
c                                                                       
c                                                                       
c                                                                       
c                                                                       
  930 if (iwrite.gt.1) go to 920                                        
c                                                                       
c                                                                       
c        write phase-shifts in short form                               
c                                                                       
      if (indwrt) go to 921                                             
      indwrt=.true.                                                     
      if (inddeg) go to 942                                             
      write (kwrite,10053)                                              
      go to 943                                                         
  942 write (kwrite,10054)                                              
  943 continue                                                          
c                                                                       
  921 if (k.ne.1) go to 923                                             
      if (ideg.eq.1) write (kwrite,10115)                               
      if (j.ne.0) go to 922                                             
      write (kwrite,10052)                                              
      go to 923                                                         
  922 write (kwrite,10050) multi(1),spd(j1),j,                          
     1                     multi(2),spd(j1),j,                          
     2                     multi(2),spd(j),j,                           
     3                     j,                                           
     4                     multi(2),spd(j2),j                           
  923 if (j.ne.0) go to 924                                             
      del3p0=delta(3)                                                   
      delta(3)=0.                                                       
      delta(4)=del3p0                                                   
  924 if (.not.inddeg) go to 926                                        
      do 925 iv=1,5                                                     
  925 ddelta(iv)=delta(iv)*rd                                           
      write (kwrite,10055) elab(k),(ddelta(iv),iv=1,3),ddelta(5),       
     1 ddelta(4)                                                        
      go to 1000                                                        
  926 write (kwrite,10051) elab(k),(delta(iv),iv=1,3),delta(5),delta(4) 
      go to 1000                                                        
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c        store phase shifts and compute chi square                      
c        -----------------------------------------                      
c                                                                       
c                                                                       
  920 do 903 iv=1,5                                                     
      if (.not.sing.and.iv.eq.1) go to 903                              
      if (.not.trip.and.iv.eq.2) go to 903                              
      if (.not.coup.and.iv.ge.3) go to 903                              
      if (j.eq.0.and.iv.ge.4) go to 903                                 
      if (j.eq.0.and.iv.eq.2) go to 903                                 
c                                                                       
c                                                                       
      delr(1,k,iv)=delta(iv)                                            
c                                                                       
      delta(iv)=delta(iv)*rd                                            
c                                                                       
      deld(1,k,iv)=delta(iv)                                            
c                                                                       
c                                                                       
      indelb=.false.                                                    
      do 907 ips=1,nps                                                  
      if (.not.indps(ips)) go to 904                                    
      kke=mpsa(iv,j1,ips)                                               
      if (kke.eq.0) go to 904                                           
c                                                                       
c                                                                       
c        look for identical experimental energy                         
c                                                                       
      do 901 kk=1,kke                                                   
      if (elab(k).eq.psa(1,kk,iv,j1,ips)) go to 902                     
  901 continue                                                          
c                                                                       
c        an experimental value does not exist                           
c                                                                       
  904 if (indelb) go to 907                                             
      llps(k,iv)=blanks                                                 
      ichid=0                                                           
      chi(k,iv)=0.                                                      
      delr(2,k,iv)=0.                                                   
      delr(3,k,iv)=0.                                                   
      delr(4,k,iv)=0.                                                   
      deld(2,k,iv)=0.                                                   
      deld(3,k,iv)=0.                                                   
      deld(4,k,iv)=0.                                                   
      go to 907                                                         
c                                                                       
c        experimental value exists                                      
c                                                                       
c        store label of phase-shift analysis                            
  902 indelb=.true.                                                     
      llps(k,iv)=labps(ips)                                             
      ichid=1                                                           
c        compute chi square                                             
      chi(k,iv)=((delta(iv)-psa(2,kk,iv,j1,ips))/psa(3,kk,iv,j1,ips     
     1 ))**2                                                            
c                                                                       
c        compute and store error bar                                    
c                                                                       
  906 deld(2,k,iv)=psa(2,kk,iv,j1,ips)                                  
      deld(3,k,iv)=psa(2,kk,iv,j1,ips)+psa(3,kk,iv,j1,ips)              
      deld(4,k,iv)=psa(2,kk,iv,j1,ips)-psa(3,kk,iv,j1,ips)              
c                                                                       
      delr(2,k,iv)=deld(2,k,iv)*dr                                      
      delr(3,k,iv)=deld(3,k,iv)*dr                                      
      delr(4,k,iv)=deld(4,k,iv)*dr                                      
c                                                                       
c                                                                       
  907 continue                                                          
      ichi=ichi+ichid                                                   
      chitot=chitot+chi(k,iv)                                           
c                                                                       
c                                                                       
  903 continue                                                          
c                                                                       
c                                                                       
c                                                                       
 1000 continue                                                          
c        this has been the end of the elab loop                         
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c        low energy parameters                                          
c        ---------------------                                          
c                                                                       
c                                                                       
      if (.not.inderg) go to 1100                                       
      if (j.gt.1) go to 1100                                            
      if (j.eq.0) then                                                  
c                                                                       
c   for 1s0 low energy parameters have been already calculated above    
c                                                                       
        lab1=blanks                                                     
        lab2=blanks                                                     
      else                                                              
          if(coup) then                                                 
        rb2(3)=(-4./wn*(rb(2)-rb(1))/(elab(2)-elab(1)))*uf              
        rb1(3)=rb(1)+wn/4.*elab(1)*rb2(3)/uf                            
        rb1(3)=1./rb1(3)*uf                                             
        lab1=blanks                                                     
        lab2=blanks                                                     
                endif                                                   
      endif                                                             
      if (iwrite.gt.1) go to 1120                                       
      do 1005 iv=1,3                                                    
 1005 write (kwrite,10016) lab1,rb1(iv),lab2,rb2(iv)                    
c                                                                       
c                                                                       
c                                                                       
c        write phase shifts                                             
c        ------------------                                             
c                                                                       
c                                                                       
 1100 if (iwrite.le.1) go to 2000                                       
 1120 do 1101 iv=1,5                                                    
      if (.not.sing.and.iv.eq.1) go to 1101                             
      if (.not.trip.and.iv.eq.2) go to 1101                             
      if (.not.coup.and.iv.ge.3) go to 1101                             
      if (j.eq.0.and.iv.ge.4) go to 1101                                
      if (j.eq.0.and.iv.eq.2) go to 1101                                
c                                                                       
c                                                                       
c        check line number and write state                              
c                                                                       
      if (line.lt.lppg-4) go to 1103                                    
      call headln (iwrite)                                              
      line=0                                                            
 1103 line=line+4                                                       
c                                                                       
      if (ideg.eq.1) write (kwrite,10115)                               
c                                                                       
      ichar=blanks                                                      
      if (indcou) ichar=nnpp                                            
      if (iv.ne.5) go to 1105                                           
      state(1)=blank                                                    
      state(2)=eeps                                                     
      go to  1112                                                       
 1105 state(1)=multi(iv)                                                
      ispd=j1+ldel(iv)                                                  
      if (j.ne.0) go to 1108                                            
      if (iv.ge.3) ispd=2                                               
 1108 state(2)=spd(ispd)                                                
 1112 go to (1141,1141,1142),iwrite                                     
 1141 write (kwrite,10061) state,j,ichar                                
      go to 1143                                                        
 1142 write (kwrite,10011) state,j,ichar                                
 1143 continue                                                          
c                                                                       
c                                                                       
      do 1107 k=1,melab                                                 
c                                                                       
c        check line number and write one line of phase shifts           
c                                                                       
      if (line.lt.lppg) go to 1104                                      
      call headln (iwrite)                                              
      line=0                                                            
 1104 line=line+1                                                       
      go to (1131,1131,1132),iwrite                                     
 1131 write (kwrite,10010) elab(k),(deld(i,k,iv),i=1,4),llps(k,iv),     
     1 chi(k,iv)                                                        
      go to 1107                                                        
 1132 write (kwrite,10010) elab(k),(deld(i,k,iv),i=1,4),llps(k,iv),     
     1 chi(k,iv),(delr(ii,k,iv),ii=1,4)                                 
 1107 continue                                                          
c                                                                       
c                                                                       
c        write low energy parameters                                    
c        ---------------------------                                    
c                                                                       
      if (.not.inderg) go to 1101                                       
      if (j.ne.0) go to 1102                                            
      lab1=blanks                                                       
      lab2=blanks                                                       
      go to 1111                                                        
 1102 if (j.ge.2.or.iv.gt.3) go to 1101                                 
      lab1=blanks                                                       
      lab2=blanks                                                       
 1111 if (line.lt.lppg-1) go to 1110                                    
      call headln (iwrite)                                              
      line=0                                                            
 1110 line=line+2                                                       
      write (kwrite,10016) lab1,rb1(iv),lab2,rb2(iv)                    
c                                                                       
 1101 continue                                                          
c                                                                       
c                                                                       
c        write chi**2                                                   
c        ------------                                                   
c                                                                       
c                                                                       
 1300 if (ichi.eq.0) go to 1302                                         
      chipd=chitot/dfloat(ichi)                                         
      go to 1303                                                        
 1302 chipd=0.d0                                                        
 1303 if (line.lt.lppg-3) go to 1301                                    
      call headln (iwrite)                                              
      line=0                                                            
 1301 line=line+4                                                       
      write (kwrite,10015) chitot,chipd,ichi                            
c                                                                       
c                                                                       
      if (.not.indpch) go to 2000                                       
c                                                                       
c        punch                                                          
c        -----                                                          
c                                                                       
      do 1201 iv=1,5                                                    
      if (.not.sing.and.iv.eq.1) go to 1201                             
      if (.not.trip.and.iv.eq.2) go to 1201                             
      if (.not.coup.and.iv.ge.3) go to 1201                             
      if (j.eq.0.and.iv.ge.4) go to 1201                                
      if (j.eq.0.and.iv.eq.2) go to 1201                                
      name(2)=blanks                                                    
      name(3)=blanks                                                    
c                                                                       
c                                                                       
c        define state                                                   
c                                                                       
      if (iv.ne.5) go to 1204                                           
c                                                                       
c        case of epsilon                                                
c                                                                       
      state(1)=blank                                                    
      state(2)=eeps                                                     
c        branch if short punch                                          
      if (indshp) go to 1202                                            
c        punch epsilon                                                  
      write (kpunch,10022) eps,j                                        
      go to 1215                                                        
c                                                                       
c        case of partial wave state                                     
c                                                                       
 1204 state(1)=multi(iv)                                                
      ispd=j1+ldel(iv)                                                  
      if (j.ne.0) go to 1211                                            
      if (iv.eq.3) go to 1210                                           
      state(1)=multi(1)                                                 
      go to 1211                                                        
 1210 ispd=2                                                            
 1211 state(2)=spd(ispd)                                                
c                                                                       
c                                                                       
c        branch if short punch                                          
c                                                                       
      if (indshp) go to 1202                                            
c                                                                       
c                                                                       
c        punch partial wave state                                       
c                                                                       
      write (kpunch,10021) state,j                                      
c                                                                       
c        punch plot information for axes                                
 1215 write (kpunch,10026)                                              
      write (kpunch,10031)                                              
      write (kpunch,10032)                                              
      if (iv.eq.5) go to 1213                                           
      write (kpunch,10033)                                              
      if (inddeg) go to 1225                                            
      write (kpunch,10034)                                              
      go to 1205                                                        
 1225 write (kpunch,10037)                                              
      go to 1205                                                        
 1213 write (kpunch,10035)                                              
      if (inddeg) go to 1227                                            
      write (kpunch,10036)                                              
      go to 1205                                                        
 1227 write (kpunch,10038)                                              
 1205 write (kpunch,10026)                                              
c                                                                       
c        punch experimental error bars                                  
c                                                                       
      m1=0                                                              
      if (.not.inderr) go to 1233                                       
      do 1232 ips=1,nps                                                 
      if (.not.indps(ips)) go to 1232                                   
      mm1=mpsa(iv,j1,ips)                                               
      if (mm1.eq.0) go to 1232                                          
      m1=m1+mm1                                                         
 1232 continue                                                          
 1233 m2=m1                                                             
      name(1)=charmr                                                    
      write (kpunch,10000) name,m1,m2                                   
      if (m1.eq.0) go to 1203                                           
c                                                                       
      do 1239 ips=1,nps                                                 
      if (.not.indps(ips)) go to 1239                                   
      mm1=mpsa(iv,j1,ips)                                               
      if (mm1.eq.0) go to 1239                                          
c                                                                       
c        parameter for symbol in the middel of error bar                
      ierr=ips                                                          
      do 1235 ke=1,mm1                                                  
      expch(1)=psa(1,ke,iv,j1,ips)                                      
      expch(2)=psa(2,ke,iv,j1,ips)                                      
      expch(3)=psa(3,ke,iv,j1,ips)                                      
      if (inddeg) go to 1234                                            
      expch(2)=expch(2)*dr                                              
      expch(3)=expch(3)*dr                                              
 1234 write (kpunch,10025) labps(ips),state,j,ierr,expch                
 1235 continue                                                          
 1239 continue                                                          
c                                                                       
 1203 continue                                                          
      if (inderr) go to 1201                                            
      if (ipunch.eq.4) go to 1201                                       
c                                                                       
c        punch phase shifts                                             
c                                                                       
 1202 name(1)=charlb                                                    
      write (kpunch,10008) name,label                                   
      m1=melab                                                          
      m2=m1                                                             
      name(1)=charm                                                     
      write (kpunch,10000) name,m1,m2                                   
      do 1249 k=1,melab                                                 
      if (inddeg) go to 1208                                            
      expch(1)=delr(1,k,iv)                                             
      go to 1207                                                        
 1208 expch(1)=deld(1,k,iv)                                             
 1207 write (kpunch,10024) label,state,j,elab(k),expch(1)               
 1249 continue                                                          
c                                                                       
 1201 continue                                                          
c                                                                       
c                                                                       
c                                                                       
 2000 continue                                                          
c        this has been the end of the j loop                            
c                                                                       
c                                                                       
c                                                                       
c        calculate chi**2 using the error matrix                        
c        ---------------------------------------                        
c                                                                       
c                                                                       
      if (indema) then                                                  
c                                                                       
      go to (2001,2002),ierrma                                          
 2001 write (kwrite,11010)                                              
      go to 2003                                                        
 2002 write (kwrite,11020)                                              
 2003 continue                                                          
c                                                                       
      chisee=0.d0                                                       
      chidd=0.d0                                                        
      chierr=0.d0                                                       
      chier6=0.d0                                                       
      do 2095 k=1,melab                                                 
      do 2015 kk=1,melerr                                               
      if (elab(k).eq.elaber(kk)) go to 2020                             
 2015 continue                                                          
      go to 2095                                                        
 2020 write (kwrite,11011) elaber(kk)                                   
c                                                                       
c                                                                       
c        calculate 3p central force combinations for kk=1 or 2          
c        and store it where usually 3p0 is stored                       
c                                                                       
      if (kk.le.2) then                                                 
      delj(3,1,k)                                                       
     1      =(delj(3,1,k)                                               
     2  +3.d0*delj(2,2,k)                                               
     3  +5.d0*delj(3,3,k))/9.d0                                         
      end if                                                            
c                                                                       
c                                                                       
      chikk=0.d0                                                        
      npherk=npherr(kk)                                                 
      do 2035 i=1,npherk                                                
      if (i.eq.1) then                                                  
      diff(i)=1.d0                                                      
      else                                                              
c                                                                       
c                                                                       
c        difference between model-phase shift and analysis              
c        -------------------------------------------------              
c                                                                       
      go to (2021,2022),ierrma                                          
c                                                                       
c        use pp-only error matrix, stoks1.d                             
c                                                                       
 2021 diff(i)=pherr(i,kk)+phdel(i,kk)                                   
     1       -delj(iverr1(i),jerr1(i)+1,k)                              
      go to 2023                                                        
c                                                                       
c        use pp+np error matrix, stoks2.d                               
c                                                                       
 2022 diff(i)=pherr(i,kk)+phdel(i,kk)                                   
     1       -delj(iverr2(i,kk),jerr2(i,kk)+1,k)                        
 2023 continue                                                          
c                                                                       
c                                                                       
      end if                                                            
      chisub=0.d0                                                       
      do 2025 ii=1,i                                                    
      if (i.eq.ii) then                                                 
      fachi=1.d0                                                        
      else                                                              
      fachi=2.d0                                                        
      end if                                                            
c                                                                       
c        sum up error matrix                                            
c                                                                       
      chidel=fachi*errma(i,ii,kk)*diff(i)*diff(ii)                      
c                                                                       
 2025 chisub=chisub+chidel                                              
      deldel=chisub-chidel                                              
      go to (2031,2032),ierrma                                          
 2031 write (kwrite,11012) namer1(i),chidel,deldel,chisub               
      go to 2033                                                        
 2032 write (kwrite,11032) namer2(i,kk),chidel,deldel,chisub            
 2033 continue                                                          
      if (i.eq.1) then                                                  
      chise=chidel                                                      
      end if                                                            
 2035 chikk=chikk+chisub                                                
      deldel=chikk-chise                                                
      write (kwrite,11013) elaber(kk),chise,deldel,chikk                
      chierr=chierr+chikk                                               
      if (kk.le.6) chier6=chier6+chikk                                  
      chisee=chisee+chise                                               
      chidd=chidd+deldel                                                
 2095 continue                                                          
      write (kwrite,11014) chisee,chidd,chierr                          
      write (kwrite,11024) chier6                                       
c                                                                       
c                                                                       
      if (ierrma.eq.1) then                                             
      write (kda3,11100) delj                                           
      end if                                                            
c                                                                       
c                                                                       
      end if                                                            
c                                                                       
c                                                                       
c                                                                       
c                                                                       
      stop                                                              
      end                                                               
      subroutine headln (iwrite)                                        
c                                                                       
c        writes headline for phase-shift program                        
c                                                                       
      common /crdwrt/ kread,kwrite,kpunch,kda(9)                        
c                                                                       
c                                                                       
10002 format (//' p h a s e - s h i f t s'/1x,23(1h-)/                  
     1 20x,'d e g r e e s'/                                             
     2 ' elab(mev)    theoretical  experimental  upper   lower          
     3     chi**2'/1h ,74(1h-))                                         
10003 format (//52x,'p h a s e - s h i f t s'/52x,23(1h-)/              
     1 20x,'d e g r e e s',58x,'r a d i a n s'/                         
     2 ' elab(mev)    theoretical  experimental  upper   lower          
     3     chi**2',10x,'theoretical   experimental  upper   lower'/     
     4 1h ,125(1h-))                                                    
c                                                                       
c                                                                       
      go to (1,2,3),iwrite                                              
    2 write (kwrite,10002)                                              
      go to 1                                                           
    3 write (kwrite,10003)                                              
c                                                                       
c                                                                       
    1 return                                                            
      end                                                               
c                                                                       
c   haidenbauer program !                                               
c   subroutines concerned with coulomb potential                        
c                                                                       
c                                                                       
c   partial wave projection of truncated coulomb potential              
c                                                                       
      subroutine coul(iopt)                                             
c                                                                       
      implicit real*8 (a-h,o-z)                                         
c                                                                       
c                                                                       
c                                                                       
c                                                                       
      common /crdwrt/ kread,kwrite,kpunch,kda(9)                        
c                                                                       
c        arguments of the potential subroutine pot being called in this 
c        program                                                        
c                                                                       
      common /cpot/  v(6),xmev,ymev                                     
      common /ccoul/ wn,rc,nc,iprop,indcou,indrel,indnnb                
c                                                                       
      common /cstate/ j,heform,sing,trip,coup,endep,label               
c**** common /cpts/   q(97),c,n1,ix,iy                                  
      common /cpts/   q(1001),c,n1,ix,iy                                  
c                                                                       
      logical indcou,indrel,indnnb                                      
      logical heform,sing,trip,coup,endep                               
c                                                                       
c     further specifications                                            
c                                                                       
      dimension rr(1000),wr(1000),rbesx(1000),rbesy(1000),vv(1000)                 
      logical index                                                     
c                                                                       
c     general constants                                                 
c                                                                       
      data pih/1.570796326794897d0/                                     
      data uf/197.327d0/                                                
c                                                                       
c                                                                       
      data index/.false./                                               
      data jj/-1/                                                       
c                                                                       
c                                                                       
c                                                                       
c     return if no coulomb correction is to be applied                  
c                                                                       
      if(.not.indcou) return                                            
c                                                                       
c                                                                       
      if(index) goto 100                                                
c                                                                       
      index=.true.                                                      
      alfa= 1.d0/137.036d0                                              
      ffcoul=alfa/(pih*uf*uf)                                           
      wnq=wn*wn                                                         
c                                                                       
      exysq=1.d0                                                        
c                                                                       
c**** call gset (0.d0,rc,nc,rr,wr)                                      
      call gauss(0.d0,rc,nc,rr,wr)                                      
c                                                                       
 100  continue                                                          
c                                                                       
c                                                                       
      q0q=q(n1)**2                                                      
      eq0=dsqrt(wnq+q0q)                                                
      facoul=ffcoul                                                     
      if (indrel.and.iopt.eq.0) facoul=facoul*(wnq+2.d0*q0q)/(wn*eq0)   
      if(indnnb) facoul=-facoul                                         
c                                                                       
c                                                                       
c                                                                       
      if (j.ne.jj) then                                                 
         jj=j                                                           
         if(mod(j,2).eq.0.or.indnnb) then                               
            lna=j-1                                                     
            lne=j+1                                                     
            if(lna.lt.0) lna=0                                          
         else                                                           
            lna=j                                                       
            lne=lna                                                     
         endif                                                          
      endif                                                             
c                                                                       
      do 1050 iv=1,6                                                    
 1050 vv(iv)=0.d0                                                       
c                                                                       
c                                                                       
c                                                                       
      do 5 ln=lna,lne                                                   
      if(lna.eq.lne) then                                               
         iv=2                                                           
      else                                                              
         if(ln.lt.j) then                                               
            iv=4                                                        
         else if(ln.eq.j) then                                          
            iv=1                                                        
         else                                                           
            iv=3                                                        
         endif                                                          
      endif                                                             
c                                                                       
c                                                                       
        qx=xmev/uf                                                      
        call bessel(qx,rr,rbesx,nc,ln)                                  
        qy=ymev/uf                                                      
        call bessel(qy,rr,rbesy,nc,ln)                                  
          do 15 i=1,nc                                                  
           vv(iv)=vv(iv)+rbesx(i)*rbesy(i)*rr(i)*wr(i)                  
 15       continue                                                      
          vv(iv)=vv(iv)*facoul                                          
 5    continue                                                          
c                                                                       
c                                                                       
      if (heform) then                                                  
c                                                                       
c        transformation into helicity - basis (if requested)            
c                                                                       
         aj=dfloat(j)                                                   
         aj1=dfloat(j+1)                                                
         d2j1=1.d0/dfloat(2*j+1)                                        
         arjj1=dsqrt(aj*aj1)                                            
c                                                                       
        v3=vv(3)                                                        
        v4=vv(4)                                                        
        v34=-arjj1*(v3-v4)                                              
        vv(3)=d2j1*(aj1*v3+aj*v4)                                       
        vv(4)=d2j1*(aj*v3+aj1*v4)                                       
        vv(5)=d2j1*v34                                                  
        vv(6)=d2j1*v34                                                  
c                                                                       
      endif                                                             
c                                                                       
c**** if(iprop.eq.2) then                                               
c****    ex=dsqrt(wnq+xmev*xmev)                                        
c****    ey=dsqrt(wnq+ymev*ymev)                                        
c****    exysq=2.d0*wn/dsqrt(ex+eq0)/dsqrt(ey+eq0)                      
c        exysq=     wn/dsqrt(ex)/dsqrt(ey)                              
c        exysq=     wn/eq0                                              
c**** endif                                                             
c                                                                       
c                                                                       
c     add matrix elements                                               
c                                                                       
      do 1000 iv=1,6                                                    
      if (iopt.eq.0) then                                               
      v(iv)=v(iv)+vv(iv)*exysq                                          
      else                                                              
      v(iv)=vv(iv)*exysq                                                
      end if                                                            
 1000 continue                                                          
c                                                                       
      return                                                            
      end                                                               
c                                                                       
c                                                                       
c   calculation of bessel functions                                     
c                                                                       
      subroutine bessel(sk,rr,rbes,nc,ln)                               
c                                                                       
c**** this routine has been changed on 3/23/92                          
c**** by building sbess in.                                             
c                                                                       
      implicit real*8(a-h,o-z)                                          
      dimension rr(nc),rbes(nc)                                         
      common /crdwrt/ kread,kwrite,kpunch,kda(9)                        
c                                                                       
c                                                                       
         do 11 i=1,nc                                                   
         z=sk*rr(i)                                                     
         rbes(i)=sbess(ln,z)                                            
 11      continue                                                       
c                                                                       
      return                                                            
      end                                                               
c                                                                       
c                                                                       
c    calculation of coulomb wave functions f(eta,rho) and g(eta,rho)    
c                                                                       
c    xx:   rho         eta1:      eta         xlambd:  float(l)         
c                                                                       
      subroutine klein(xx,eta1,xlamda,fc,gc,fcp,gcp,ifail)              
      implicit real*8 (a-h,o-z)                                         
            common /steed/ paccq,nfp,npq,iexp                           
            logical xlturn                                              
c  ***                                                                  
c  ***  real coulomb functions for x .gt. 0.,eta real,xl=lamda.gt.-1.   
c  ***  program in comp phys commun 1981              a r barnett       
c  ***  suitable accur = 10**-16 if 56-bit mantissa, 10**-14 if 48-bit  
c  ***  common block is for information only.  not required in code.    
c  ***                                                                  
      data zero,one,two,ten2,abort /0.0d0,1.0d0,2.0d0,1.0d2,2.0d4/      
c                                                                       
c                                                                       
c                                                                       
                 accur = 1.0d-16                                        
      ifail = 0                                                         
      iexp  = 1                                                         
      nfp   = 0                                                         
      npq   = 0                                                         
      eta   = eta1                                                      
      acc   = accur                                                     
      acc4  = acc*ten2**2                                               
      acch  = dsqrt(acc)                                                
      gjwkb = zero                                                      
      paccq = one                                                       
         if(dabs(xx) .le. acch)              go to 100                  
      x     =    xx                                                     
             if( xx  .gt. zero)              go to 1                    
             if(eta  .eq. zero)              go to 100                  
      x     =   -xx                                                     
      eta   =   -eta                                                    
    1    if(xlamda .le. -one) xl = -xlamda - one                        
         if(xlamda .gt. -one) xl =  xlamda                              
      pk    = xl + one                                                  
      e2ll1 = eta*eta + xl*pk                                           
      xlturn= x*(x - two*eta) .lt. xl*pk                                
c                                                                       
c  ***   evaluate cf1  =  f   =  fprime(xl,x,eta)/f(xl,x,eta)           
c                                                                       
      xi  = one/x                                                       
      fcl = one                                                         
      px  = pk + abort                                                  
    2 ek  = eta/pk                                                      
      f   = (ek + pk*xi)*fcl + (fcl - one)*xi                           
      pk1 =  pk + one                                                   
c  ***  test ensures b1 .ne. zero for negative eta: fixup is exact.     
        if(dabs(eta*x + pk*pk1) .gt. acc)    go to 3                    
             fcl = (one + ek*ek)/(one + (eta/pk1)**2)                   
             pk  =  two + pk                                            
      go to 2                                                           
    3 d   =  one/((pk + pk1)*(xi + ek/pk1))                             
      df  = -fcl*(one + ek*ek)*d                                        
            if(fcl .ne. one ) fcl = -one                                
            if(d   .lt. zero) fcl = -fcl                                
      f   = f  + df                                                     
c                                                                       
c  ***  begin cf1 loop on pk = k = lamda + 1                            
c                                                                       
      p     = one                                                       
    4 pk    = pk1                                                       
        pk1 = pk1 + one                                                 
        ek  = eta / pk                                                  
        tk  = (pk + pk1)*(xi + ek/pk1)                                  
        d   = tk - d*(one + ek*ek)                                      
              if(dabs(d) .gt. acch)          go to 5                    
              write (6,1000) d,df,acch,pk,ek,eta,x                      
              p = p  +   one                                            
              if( p .gt. two)                go to 110                  
    5 d     = one/d                                                     
              if (d .lt. zero) fcl = -fcl                               
        df  = df*(d*tk - one)                                           
        f   = f  + df                                                   
              if(pk .gt. px)                 go to 110                  
      if(dabs(df) .ge. dabs(f)*acc)          go to 4                    
              nfp = pk - xl - 1                                         
      if(xlturn) call jwkb(x,eta,xl,fjwkb,gjwkb,iexp)                   
      if(iexp .gt. 1 .or. gjwkb .gt. one/(acch*ten2))  go to 7          
c                                                                       
c  ***   evaluate cf2 = p + i.q  again using steed's algorithm          
c                                                                       
      ta = two*abort                                                    
      p  = zero                                                         
      q  = one - eta*xi                                                 
      pk = zero                                                         
      ar = -e2ll1                                                       
      ai = eta                                                          
      br = two*(x - eta)                                                
      bi = two                                                          
      wi = two*eta                                                      
      dr =  br/(br*br + bi*bi)                                          
      di = -bi/(br*br + bi*bi)                                          
      dp = -xi*(ar*di + ai*dr)                                          
      dq =  xi*(ar*dr - ai*di)                                          
    6 p     = p  + dp                                                   
         q  = q  + dq                                                   
         pk = pk + two                                                  
         ar = ar + pk                                                   
         ai = ai + wi                                                   
         bi = bi + two                                                  
         d  = ar*dr - ai*di + br                                        
         di = ai*dr + ar*di + bi                                        
         c  = one/(d*d + di*di)                                         
         dr =  c*d                                                      
         di = -c*di                                                     
         a  = br*dr - bi*di - one                                       
         b  = bi*dr + br*di                                             
         c  = dp*a  - dq*b                                              
         dq = dp*b  + dq*a                                              
         dp = c                                                         
         if(pk .gt. ta)                      go to 120                  
      if(dabs(dp)+dabs(dq).ge.(dabs(p)+dabs(q))*acc)   go to 6          
                     npq   = pk/two                                     
                     paccq = acc/(two*dmin1(dabs(q),one))               
                     if(dabs(p) .gt. dabs(q)) paccq = paccq*dabs(p)     
c                                                                       
c *** solve for fp,g,gp and normalise f  at lamda = xl                  
c                                                                       
            if(q .le. acc4*dabs(p))          go to 130                  
      gam = (f - p)/q                                                   
      w   = fcl/dsqrt((f - p)*gam + q)                                  
            go to 8                                                     
    7 w   = fjwkb                                                       
      gc  = gjwkb                                                       
      gcp = f*gc - one/w                                                
            go to 9                                                     
    8 gc  = w*gam                                                       
      gcp = w*(p*gam - q)                                               
    9 fc  = w                                                           
      fcp = w*f                                                         
      return                                                            
 1000 format(/' cf1 accuracy loss: d,df,acch,k,eta/k,eta,x ='/,1pd11.2) 
c                                                                       
c  ***   error messages                                                 
c                                                                       
  100 ifail = -1                                                        
      write(6,2000) xx,acch,eta,xl                                      
 2000 format(' for xx = ',1pd10.3,' try small-x asymptotic solutions',3x
     *,'sqrt(acc),eta,xl = ',d8.1,0pf10.2,f9.2)                         
      return                                                            
  110 ifail =  1                                                        
      write (6,2010) abort,eta,x,xl,f,df,pk,px,acc                      
 2010 format(/' cf1 has failed to converge after ',f10.0,' iterations', 
     *'  eta,x,xl = ',3f10.3,/,' f,df,pk,px,accur  = ',1p,5d13.4//)     
      return                                                            
  120 ifail =  2                                                        
      write (6,2020) abort,eta,x,xl,p,q,dp,dq,acc                       
 2020 format(/' cf2 has failed to converge after ',f10.0,' iterations', 
     *'  eta,x,xl = ',3f10.3,/,' p,q,dp,dq,accur = ',1p,5d12.3//)       
      return                                                            
  130 ifail =  3                                                        
      write (6,2030) p,q,acc, eta,x,xl                                  
 2030 format(2x,' final q.le./p/*acc*10**4,  p,q,acc = ',1p,3d12.3,     
     *'  eta,x,xl = ',0p,3f10.3/)                                       
      return                                                            
      end                                                               
      subroutine jwkb(xx,eta1,xl,fjwkb,gjwkb,iexp)                      
      real*8    dzero,xx,eta1,xl,fjwkb,gjwkb                            
      common   /print/ iprnt                                            
c *** computes jwkb approximations to coulomb functions    for xl.ge. 0 
c *** as modified by biedenharn et al. phys rev 97 (1955) 542-554       
c *** calls dmax1,sqrt,alog,exp,atan2,float,int   a.r.barnett feb 1981  
      data    zero,half,one,six,ten/ 0.0e0, 0.5e0, 1.0e0, 6.0e0, 1.0e1/ 
      data   dzero, rl35, aloge  /0.0d0, 35.0e0, 0.43429 45 e0/         
c                                                                       
c                                                                       
c                                                                       
c                                                                       
      x    = xx                                                         
      eta  = eta1                                                       
      gh2  = x*(eta + eta - x)                                          
      xll1 = dmax1(xl*xl + xl,dzero)                                    
      if(gh2 + xll1 .le. zero) return                                   
      hll = xll1 + six/rl35                                             
      hl  = sqrt(hll)                                                   
      sl  = eta/hl + hl/x                                               
      rl2 = one + eta*eta/hll                                           
      gh  = sqrt(gh2 + hll)/x                                           
      phi = x*gh - half*( hl*log((gh + sl)**2/rl2) - log(gh) )          
          if(eta .ne. zero) phi = phi - eta* atan2(x*gh,x - eta)        
      phi10 = phi * aloge                                               
      iexp  =   int(-phi10)                                             
      if(iexp .gt. 70) gjwkb = ten**(-phi10 -  float(iexp))             
      if(iexp .le. 70) gjwkb = exp(-phi)                                
      if(iexp .le. 70) iexp  = 0                                        
      fjwkb = half/(gh*gjwkb)                                           
      if(iprnt .eq. 1) write (6,1000) iexp,fjwkb,gjwkb,eta,x,xl         
 1000 format(/,' jwkb calcs  iexp =',i8,1pd13.4,' = f ',2x,d13.4,' = g',
     *6x,' eta,x,xl = ',0pf10.3,2f9.2)                                  
      return                                                            
      end                                                               
c                                                                       
c                                                                       
      subroutine mult(a,b,c)                                            
      implicit real*8(a-h,o-z)                                          
      dimension a(2,2),b(2,2),c(2,2)                                    
      c(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)                                
      c(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)                                
      c(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)                                
      c(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)                                
      return                                                            
      end                                                               
      subroutine inv(a,b)                                               
      implicit real*8(a-h,o-z)                                          
      dimension a(2,2),b(2,2)                                           
      det=a(1,1)*a(2,2)-a(1,2)*a(2,1)                                   
      b(1,1)=a(2,2)/det                                                 
      b(1,2)=-a(1,2)/det                                                
      b(2,1)=-a(2,1)/det                                                
      b(2,2)=a(1,1)/det                                                 
      return                                                            
      end                                                               
c                                                                       
c     solves equation   y=r+x*k**2+p*k**4                               
c                                                                       
      subroutine equ(a,b,c,r,x)                                         
      implicit real*8(a-h,o-z)                                          
      dimension a(3),b(3),c(3)                                          
      a1=a(1)-a(2)                                                      
      a2=a(1)-a(3)                                                      
      c1=c(1)-c(2)                                                      
      c2=c(1)-c(3)                                                      
      b1=b(1)-b(2)                                                      
      b2=b(1)-b(3)                                                      
      zz=(b1*a2-b2*a1)/(c1*a2-c2*a1)                                    
      r=(b1-zz*c1)/a1                                                   
      x=b(1)-zz*c(1)-r*a(1)                                             
      return                                                            
      end                                                               
c                                                                       
c  calculation of function for low energy expansion                     
c                                                                       
      function hfunc(ga)                                                
      implicit real*8(a-h,o-z)                                          
      data ce/0.577215664901533d0/                                      
      hm=0.d0                                                           
      fhm=0.d0                                                          
1     continue                                                          
      hm=hm+1.d0                                                        
      fad=ga*ga/(hm*hm+ga*ga)                                           
      fad=fad/hm                                                        
      fhm=fhm+fad                                                       
      if(dabs(fad/fhm).gt.1.d-12) goto 1                                
      hfunc=-ce-dlog(ga)+fhm                                            
      return                                                            
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
