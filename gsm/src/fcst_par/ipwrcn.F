      function ipwrcn(cn)                                                       
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    ipwrcn      get expected order of magnitude of a field.        
c   prgmmr: mark iredell     org: w/nmc23    date: 91-03-15                     
c                                                                               
c abstract: integer function to return the expected                             
c           order of magnitude of a field given its name.                       
c                                                                               
c program history log:                                                          
c   91-03-15  mark iredell                                                      
c                                                                               
c usage:    i = ipwrcn (cn)                                                     
c   input argument list:                                                        
c     cn       - character*8 field name                                         
c                                                                               
c   output argument list:                                                       
c     ipwrcn   - expected order of magnitude                                    
c                                                                               
c   subprograms called:                                                         
c     (isc8)   - character*8 string search                                      
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray ymp.                                                         
c                                                                               
c$$$                                                                            
      parameter(nl=56)                                                          
      character*8 cn                                                            
      character*8 cnl(nl)                                                       
      integer ipl(0:nl)                                                         
      save cnl,ipl                                                              
      data ipl(0)/0/                                                            
      data cnl( 1)/'u       '/, ipl( 1)/0/                                      
      data cnl( 2)/'v       '/, ipl( 2)/0/                                      
      data cnl( 3)/'tv      '/, ipl( 3)/0/                                      
      data cnl( 4)/'q       '/, ipl( 4)/-3/                                     
      data cnl( 5)/'vot**2  '/, ipl( 5)/-12/                                    
      data cnl( 6)/'div**2  '/, ipl( 6)/-12/                                    
      data cnl( 7)/'omega   '/, ipl( 7)/-5/                                     
      data cnl( 8)/'t       '/, ipl( 8)/0/                                      
      data cnl( 9)/'rh      '/, ipl( 9)/0/                                      
      data cnl(10)/'ke      '/, ipl(10)/0/                                      
      data cnl(11)/'pd      '/, ipl(11)/-1/                                     
      data cnl(12)/'dtconv  '/, ipl(12)/-5/                                     
      data cnl(13)/'dtlarg  '/, ipl(13)/-5/                                     
      data cnl(14)/'dtshal  '/, ipl(14)/-5/                                     
      data cnl(15)/'dtvrdf  '/, ipl(15)/-5/                                     
      data cnl(16)/'dqconv  '/, ipl(16)/-8/                                     
      data cnl(17)/'dqshal  '/, ipl(17)/-8/                                     
      data cnl(18)/'dqvrdf  '/, ipl(18)/-8/                                     
      data cnl(19)/'duvrdf  '/, ipl(19)/-5/                                     
      data cnl(20)/'dvvrdf  '/, ipl(20)/-5/                                     
      data cnl(21)/'dthsw   '/, ipl(21)/-5/                                     
      data cnl(22)/'dthlw   '/, ipl(22)/-5/                                     
cdg3..   inserted multi-layered cloud here, rather than at the end,             
cdg3..     since it seems to be just before the single lyr data...kac           
      data cnl(23)/'cloud   '/, ipl(23)/0/                                      
      data cnl(24)/'cvcld   '/, ipl(24)/0/                                      
      data cnl(25)/'rain    '/, ipl(25)/-3/                                     
      data cnl(26)/'rainc   '/, ipl(26)/-3/                                     
      data cnl(27)/'dtsfc   '/, ipl(27)/0/                                      
      data cnl(28)/'dqsfc   '/, ipl(28)/0/                                      
      data cnl(29)/'dusfc   '/, ipl(29)/-3/                                     
      data cnl(30)/'dvsfc   '/, ipl(30)/-3/                                     
      data cnl(31)/'rcov    '/, ipl(31)/0/                                      
      data cnl(32)/'rcovc   '/, ipl(32)/0/                                      
      data cnl(33)/'tskin   '/, ipl(33)/0/                                      
      data cnl(34)/'wetness '/, ipl(34)/0/                                      
      data cnl(35)/'snow    '/, ipl(35)/1/                                      
      data cnl(36)/'tg1     '/, ipl(36)/0/                                      
      data cnl(37)/'tg2     '/, ipl(37)/0/                                      
      data cnl(38)/'tg3     '/, ipl(38)/0/                                      
      data cnl(39)/'sfcsw   '/, ipl(39)/0/                                      
      data cnl(40)/'sfclw   '/, ipl(40)/0/                                      
      data cnl(41)/'zorl    '/, ipl(41)/0/                                      
      data cnl(42)/'slmsk   '/, ipl(42)/0/                                      
      data cnl(43)/'ps      '/, ipl(43)/-1/                                     
      data cnl(44)/'tvs     '/, ipl(44)/0/                                      
      data cnl(45)/'qs      '/, ipl(45)/-3/                                     
      data cnl(46)/'ts      '/, ipl(46)/0/                                      
      data cnl(47)/'rhs     '/, ipl(47)/0/                                      
      data cnl(48)/'dugwd   '/, ipl(48)/-3/                                     
      data cnl(49)/'dvgwd   '/, ipl(49)/-3/                                     
      data cnl(50)/'ua      '/, ipl(50)/10/                                     
      data cnl(51)/'uap     '/, ipl(51)/10/                                     
      data cnl(52)/'duasfc  '/, ipl(52)/3/                                      
      data cnl(53)/'duagwd  '/, ipl(53)/3/                                      
      data cnl(54)/'duamtn  '/, ipl(54)/3/                                      
      data cnl(55)/'ep      '/, ipl(55)/0/                                      
      data cnl(56)/'cldwork '/, ipl(56)/0/                                      
c                                                                               
      n=isc8(nl,cnl,cn)                                                         
      ipwrcn=ipl(n)                                                             
      return                                                                    
      end                                                                       
