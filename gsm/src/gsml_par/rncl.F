      subroutine rncl(len, pl, rno, clf)                                        
c                                                                               
c                                                                               
c                                                                               
c*********************************************************************          
c*********************** aries   model *******************************          
c********************   subroutine  rncl  ****************************          
c********************** 01 february 1991 *****************************          
c*********************************************************************          
c                                                                               
c     parameter (p5=50.0,  p8=900.0, pt8=0.8, pt2=0.2)                          
c     parameter (p5=200.0,  p8=800.0, pt8=0.4, pt2=0.6)                         
c     parameter (p5=400.0,  p8=800.0, pt8=0.5, pt2=0.5)                         
c     parameter (p5=200.0,  p8=800.0, pt8=0.8, pt2=0.2)                         
      parameter (p5=400.0,  p8=800.0, pt8=0.8, pt2=0.2)                         
c     parameter (p5=500.0,  p8=800.0, pt8=0.5, pt2=0.5)                         
      parameter (pfac=pt2/(p8-p5))                                              
c                                                                               
      parameter (p4=400.0,    p6=950.0)                                         
      parameter (crtm=5.0e-7, crtmi=1.0/crtm)                                   
      parameter (cfac=1.0/(p6-p4))                                              
c                                                                               
c                                                                               
      dimension pl(len),  rno(len), clf(len)                                    
c                                                                               
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c                                                                               
      do 10 i=1,len                                                             
      rno(i) = 1.0                                                              
      clf(i) = 1.0                                                              
c                                                                               
      if (pl(i) .ge. p5 .and. pl(i) .le. p8) then                               
          rno(i) = (p8-pl(i))*pfac + pt8                                        
      elseif (pl(i) .gt. p8 ) then                                              
          rno(i) = pt8                                                          
      endif                                                                     
      if (pl(i) .ge. p4 .and. pl(i) .le. p6) then                               
         clf(i) = (p6-pl(i))*cfac                                               
      elseif (pl(i) .gt. p6 ) then                                              
         clf(i) = 0.0                                                           
      endif                                                                     
   10 continue                                                                  
c                                                                               
      return                                                                    
      end                                                                       
