      subroutine gpvs                                                           
      save                                                                      
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: gpvs         compute saturation vapor pressure table              
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: compute saturation vapor pressure table as a function of            
c   temperature for function fpvs. the water model assumes a perfect gas        
c   constant specific heats for gas and liquid, and neglects                    
c   the volume of the liquid. the ice option is no longer included.             
c   the model does account for the variation of the latent heat                 
c   of condensation with temperature. the clausius-clapeyron equation           
c   is integrated from the triple point to get the formula                      
c       pvs=psatk*(tr**xa)*exp(xb*(1.-tr))                                      
c   where tr is ttp/t and other values are physical constants                   
c   defined in parameter statements in the code.                                
c   the current implementation computes a table with a length                   
c   of 1501 for temperatures ranging from 180. to 330. kelvin.                  
c                                                                               
c usage:  call gpvs                                                             
c                                                                               
c common blocks:                                                                
c   compvs   - scaling parameters and table for function fpvs.                  
c                                                                               
c                                                                               
c$$$                                                                            
      parameter(nx=1501)                                                        
      dimension tbpvs(nx)                                                       
      common/compvs/ c1xpvs,c2xpvs,anxpvs,tbpvs                                 
c-----------------------------------------------------------------------        
      cp = 1004.6                                                               
      rd = 287.05                                                               
      rv = 461.5                                                                
      ttp = 273.16                                                              
      hvap = 2.5e6                                                              
      psatk = 0.61078                                                           
      cliq = 4185.5                                                             
      cvap = 1846.                                                              
      dldt = cvap - cliq                                                        
      xa = -dldt / rv                                                           
      xb = xa + hvap / (rv*ttp)                                                 
      xmin=180.0                                                                
      xmax=330.0                                                                
      xinc=(xmax-xmin)/(nx-1)                                                   
      c1xpvs=1.-xmin/xinc                                                       
      c2xpvs=1./xinc                                                            
      anxpvs=nx-0.01                                                            
      do 10 jx=1,nx                                                             
        x=xmin+(jx-1)*xinc                                                      
        t=x                                                                     
        tr=ttp/t                                                                
        tbpvs(jx)=psatk*(tr**xa)*exp(xb*(1.-tr))                                
10    continue                                                                  
      print 90,xmin,xmax,xinc                                                   
90    format(' gpvs. saturation vapor pressure table created'/                  
     &       '       for temperature range and increment of ',3f8.3)            
      return                                                                    
      end                                                                       
