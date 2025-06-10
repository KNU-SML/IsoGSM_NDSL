      subroutine gau2ll(gauin,imxin,jmxin,                                      
     1                  xlonw,xlatn,dxlon,dxlat,regout,imxout,jmxout,           
     2                  undef,luptr)                                            
c                                                                               
      save                                                                      
c                                                                               
c  interpolation from lat/lon grid to other lat/lon grid                        
c                                                                               
      dimension gauin (imxin,jmxin)                                             
c                                                                               
      dimension regout(imxout,jmxout)                                           
      dimension gaul(500),regl(500)                                             
      dimension iindx1(1000)                                                    
      dimension iindx2(1000)                                                    
      dimension jindx1(500)                                                     
      dimension jindx2(500)                                                     
      dimension ddx(1000)                                                       
      dimension ddy(500)                                                        
c                                                                               
      data ifp/0/                                                               
c                                                                               
      if(ifp.ne.0) go to 111                                                    
      ifp=1                                                                     
c                                                                               
      write(luptr,*) 'imxin=',imxin,' jmxin=',jmxin                             
      write(luptr,*) 'xlatn=',xlatn,' xlonw=',xlonw,                            
     1               ' dxlat=',dxlat,' dxlon=',dxlon                            
      write(luptr,*) 'imxout=',imxout,' jmxout=',jmxout                         
c                                                                               
      call gaulat(gaul,jmxin)                                                   
c                                                                               
      do 20 j=1,jmxout                                                          
      regl(j)=xlatn-float(j-1)*dxlat                                            
   20 continue                                                                  
c                                                                               
      dxin =360./float(imxin)                                                   
c                                                                               
      do 30 i=1,imxout                                                          
      alamd=xlonw+float(i-1)*dxlon                                              
      if(alamd.lt.0.) alamd=360.+alamd                                          
      i1=alamd/dxin+1.001                                                       
      if(i1.gt.imxin) i1=1                                                      
      iindx1(i)=i1                                                              
      i2=i1+1                                                                   
      if(i2.gt.imxin) i2=1                                                      
      iindx2(i)=i2                                                              
      ddx(i)=(alamd-float(i1-1)*dxin)/dxin                                      
   30 continue                                                                  
c                                                                               
      do 40 j=1,jmxout                                                          
      aphi=regl(j)                                                              
      do 50 jj=1,jmxin                                                          
      if(aphi.lt.gaul(jj)) go to 50                                             
      j2=jj                                                                     
      go to 42                                                                  
   50 continue                                                                  
      j2=jmxin                                                                  
   42 continue                                                                  
      if(j2.gt.2) go to 43                                                      
      j1=1                                                                      
      j2=2                                                                      
      go to 44                                                                  
   43 continue                                                                  
      if(j2.le.jmxin) go to 45                                                  
      j1=jmxin-1                                                                
      j2=jmxin                                                                  
      go to 44                                                                  
   45 continue                                                                  
      j1=j2-1                                                                   
   44 continue                                                                  
      jindx1(j)=j1                                                              
      jindx2(j)=j2                                                              
      ddy(j)=(aphi-gaul(j1))/(gaul(j2)-gaul(j1))                                
   40 continue                                                                  
c                                                                               
  111 continue                                                                  
c                                                                               
c     write(luptr,*) 'iindx1'                                                   
c     write(luptr,*) (iindx1(n),n=1,imxout)                                     
c     write(luptr,*) 'iindx2'                                                   
c     write(luptr,*) (iindx2(n),n=1,imxout)                                     
c     write(luptr,*) 'jindx1'                                                   
c     write(luptr,*) (jindx1(n),n=1,jmxout)                                     
c     write(luptr,*) 'jindx2'                                                   
c     write(luptr,*) (jindx2(n),n=1,jmxout)                                     
c     write(luptr,*) 'ddy'                                                      
c     write(luptr,*) (ddy(n),n=1,jmxout)                                        
c     write(luptr,*) 'ddx'                                                      
c     write(luptr,*) (ddx(n),n=1,jmxout)                                        
c                                                                               
      do 60 j=1,jmxout                                                          
      y=ddy(j)                                                                  
      j1=jindx1(j)                                                              
      j2=jindx2(j)                                                              
      do 60 i=1,imxout                                                          
      x=ddx(i)                                                                  
      i1=iindx1(i)                                                              
      i2=iindx2(i)                                                              
      if(gauin(i1,j1).eq.undef.or.gauin(i2,j1).eq.undef.or.                     
     1   gauin(i1,j2).eq.undef.or.gauin(i2,j2).eq.undef) then                   
        regout(i,j)=undef                                                       
      else                                                                      
        regout(i,j)=(1.-x)*(1.-y)*gauin(i1,j1)+(1.-y)*x*gauin(i2,j1)+           
     1           (1.-x)*y*gauin(i1,j2)+x*y*gauin(i2,j2)                         
      endif                                                                     
   60 continue                                                                  
c                                                                               
      sum1=0.                                                                   
      sum2=0.                                                                   
      do 70 i=1,imxin                                                           
      sum1=sum1+gauin(i,1)                                                      
      sum2=sum2+gauin(i,jmxin)                                                  
   70 continue                                                                  
      sum1=sum1/float(imxin)                                                    
      sum2=sum2/float(imxin)                                                    
      do i=1,imxin                                                              
        if(gauin(i,1).eq.undef) sum1=undef                                      
        if(gauin(i,jmxin).eq.undef) sum2=undef                                  
      enddo                                                                     
c                                                                               
      do 80 i=1,imxout                                                          
      if(abs(regl(1)).eq.90.) then                                              
        regout(i,     1)=sum1                                                   
      endif                                                                     
      if(abs(regl(jmxout)).eq.90.) then                                         
        regout(i,jmxout)=sum2                                                   
      endif                                                                     
   80 continue                                                                  
c                                                                               
      return                                                                    
      end                                                                       
