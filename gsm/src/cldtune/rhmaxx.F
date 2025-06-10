#define LEVS levs
#define LONF lonf
      subroutine rhmaxx(ps,rh,t,vvel,si,sl,vvcldin,rhmdl,xlatrd,                
     1                  nundg,nok,nvvel)                                        
c                                                                               
#include <paramodel.h>
      save                                                                      
c...                                                                            
      parameter(idim=lonf_,kdim=levs_,mcld=3)                                   
      parameter(kdimp=kdim+1,levm1=kdim-1,levm2=kdim-2)                         
c...                                                                            
c...     part of cloud parameterization code to get max model rh                
c         in the 'mcld' mrf cloud domains...                                    
c30sep93    use new yh/kac cloud code..                                         
c...                                                                            
c...   input - ps(cb),rh(frac),t(deg k)                                         
c             - si,sl=sigma at levels,layers.                                   
c       output- rh max(-fraction-for l,m,h cloud)                               
      dimension ps(idim),si(kdimp),sl(kdim)                                     
      dimension rh(idim,kdim),t(idim,kdim),rhmdl(idim,mcld)                     
      dimension vvel(idim,kdim)                                                 
      dimension ptopc(mcld,2),ptop(mcld),cc(9),p(9)                             
c====>    workspace......                                                       
      dimension prslv(idim,kdimp),delp(idim,kdim)                               
      dimension theta(idim,kdim),dthdp(idim,kdim),                              
     2          prsly(idim,kdim)                                                
      dimension kbot(idim,mcld),ktop(idim,mcld)                                 
      dimension ivv(mcld),vvcld(2,mcld)                                         
c                                                                               
      dimension nundg(mcld),nok(mcld),nvvel(mcld)                               
c                                                                               
      logical usemax                                                            
      data usemax/.true./                                                       
c                                                                               
c...      jm slingo tables for convective camt and precip rate,based on         
c       cc=a+b*ln(p)--we ll simply linearly interpolate for now(sep87)          
c==>    done in kuo (tables will be changed in the future....)                  
      data cc/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/                               
      data p /.14,.31,.7,1.6,3.4,7.7,17.0,38.0,85.0 /                           
c--     top of each cld type(l,m,h) domain (cb) - for psfc=100. model           
c         tops vary by latitude 1st mcld valus=0-45 deg, 2nd set=polar          
c...note.. these differ from af hgts, cause we lowered ptopc(2,  )              
c..          values to try to catch more af high cloud                          
      data ptopc/64.2,35.0,15.0, 75.0,50.0,15.0/                                
      data ivv/1,0,0/                                                           
c     data vvcld/.0004,-.0006,.0005,-.0005,0.,-.00005/                          
      data vvcld/.005,0.0,0.0,0.0,0.0,0.0/                                      
      data kfir/0/                                                              
c                                                                               
      if(vvcldin.gt.-9.99.and.vvcld(1,1).ne.vvcldin) then                       
        vvcld(1,1)=vvcldin                                                      
        print *,'vvcld(1,1) modified to ',vvcld(1,1)                            
      elseif(vvcldin.eq.-9.99) then                                             
        print *,'vvcld(1,1)',vvcld(1,1)                                         
        vvcldin=vvcld(1,1)                                                      
      endif                                                                     
c===>    begin here ................................................            
      rocp = 287.05e0 / 1005.0e0                                                
c... default llyr .. for cloud domain calculation                               
      kl = kdimp                                                                
c... topmost noncloud layer will be the one at or above lowest                  
c         0.1 of the atmosphere..                                               
      xthk = 0.                                                                 
      do 202 k=1,kdim                                                           
       xthk = xthk + si(k) - si(k+1)                                            
       if (xthk.lt.0.1) go to 202                                               
       kl = k                                                                   
       go to 204                                                                
  202 continue                                                                  
  204 llyr = kl                                                                 
      if (kfir.le.0) print 208,llyr                                             
  208 format(1h ,'...llyr=',i4)                                                 
c--    compute model lyr at top of each stratiform domain                       
      ls = llyr+2                                                               
      if (ls.lt.kdim) go to 206                                                 
      print 207                                                                 
  207 format(1h ,' ************** bad llyr,so no clouds computed **')           
      return                                                                    
  206 continue                                                                  
c....   coming in have rh in all layers, forced to be between 0.and 1.          
c...  find top pressure for each cloud domain, varies with latitude             
      pi = 3.14159265                                                           
      fac = amax1(0.0, 4.0*xlatrd/pi-1.0)                                       
      do 86 k=1,mcld                                                            
       ptop(k) = ptopc(k,1) + (ptopc(k,2)-ptopc(k,1)) * fac                     
   86 continue                                                                  
c...   in order to do atmosphere above lyr 2, uncomment next 2 lines..          
c.      i.e. if we want low domain defined down to layer 2                      
call  llyr = 2                                                                  
call  ls = llyr + 2                                                             
c..       set the array which will contain max rh to missing                    
      do 5 k=1,mcld                                                             
        do 5 i=1,idim                                                           
          rhmdl(i,k) = -1.                                                      
    5 continue                                                                  
      do 6 k=1,mcld                                                             
        do 6 i=1,idim                                                           
          ktop(i,k) = 0                                                         
          kbot(i,k) = 0                                                         
    6 continue                                                                  
c--    get pressure at all model levels(not layers)                             
      do 10 k=1,kdimp                                                           
        do 10 i=1,idim                                                          
          prslv(i,k) = ps(i) * si(k)                                            
   10 continue                                                                  
c--    get pressure thick for all model lyrs                                    
      do 12 k=1,kdim                                                            
        do 12 i=1,idim                                                          
          delp(i,k) = prslv(i,k)-prslv(i,k+1)                                   
   12 continue                                                                  
c                                                                               
      do 30 k=1,mcld                                                            
        do 25 i=1,idim                                                          
c--   if underground cld then skip (ktop=kbot=0)                                
          if (ptop(k).ge.prslv(i,llyr+1)) go to 25                              
          do 15 kkk=ls,kdimp                                                    
            if (prslv(i,kkk).gt.ptop(k)) go to 15                               
c--    model layer..                                                            
            ken = kkk - 1                                                       
            go to 20                                                            
   15     continue                                                              
          print 111,i,k                                                         
          go to 25                                                              
  111 format (1h ,'  cannot find top layer for cloud at i=',i5,' type=',        
     1        i5)                                                               
   20     ktop(i,k) = ken                                                       
          if (ken.ge.levm2) ktop(i,k) = levm2 - 1                               
   25   continue                                                                
   30 continue                                                                  
c--    set the bottom of each domain..                                          
      do 40 i=1,idim                                                            
c--    find first nonzero top                                                   
        do 35 k=1,mcld                                                          
          if (ktop(i,k).eq.0) go to 35                                          
          kt = k                                                                
          go to 36                                                              
   35   continue                                                                
c--    no room for any cld type                                                 
        go to 40                                                                
   36   continue                                                                
        kbot(i,kt) = llyr + 1                                                   
        if (kt.eq.mcld) go to 40                                                
        kt1 = kt + 1                                                            
        do 37 k=kt1,mcld                                                        
          kbot(i,k) = ktop(i,k-1) + 1                                           
   37   continue                                                                
   40 continue                                                                  
c---     compute potential temp and layer pressure(mb)                          
      do 65 i=1,idim                                                            
        do 65 k=1,kdim                                                          
          prsly(i,k) = ps(i) * sl(k) * 10. e0                                   
          exnr = (prsly(i,k)/1000. e0) ** rocp                                  
          theta(i,k) = t(i,k) / exnr                                            
   65 continue                                                                  
c===>    compute dthdp,then look in hi cloud regime for stratosphere            
c      by assuming trop above dthdp lt-.25(obtained from looking at             
c      one case)--don t allow clouds here.....                                  
      do 70 i=1,idim                                                            
        do 70 k=1,levm1                                                         
          dthdp(i,k)=(theta(i,k+1)-theta(i,k))/(prsly(i,k+1)-prsly(i,k))        
   70 continue                                                                  
      do 73 i=1,idim                                                            
        if (ktop(i,mcld).eq.0) go to 73                                         
        kt = ktop(i,mcld)                                                       
        kb = kbot(i,mcld)                                                       
        kzero = 0                                                               
        do 71 kkk=kb,kt                                                         
          if (dthdp(i,kkk).gt.-.25 e0) go to 71                                 
          kzero = kkk                                                           
          go to 74                                                              
   71   continue                                                                
   74   if (kzero.lt.kb.or.kzero.gt.kt) go to 73                                
        do 72 kkk=kzero,kt                                                      
          rh  (i,kkk) = 0. e0                                                   
   72   continue                                                                
   73 continue                                                                  
      if (kfir.le.0) then                                                       
       do 1993 kcl=1,mcld                                                       
        print 1994,kcl,ktop(1,kcl),kbot(1,kcl)                                  
 1993  continue                                                                 
 1994  format (1h ,'..first point..cldtype=',i2,'..toplyr,botlyr=',2i4)         
      end if                                                                    
c--    refine the low cloud estimates for frontal(vvel) and boundary            
c        layer(lapse rate inversion) processes. j slingo finds them             
c        to be (generally) mutually exclusive....                               
c--    vvel is vertical motion in sigma layers(convert to mb/sec,               
c        which=100 pa/sec-the units in julia s paper)                           
      do 1100 kc=1,mcld                                                         
        if (ivv(kc).le.0) go to 1100                                            
        do 60 i=1,idim                                                          
        kt = ktop(i,kc)                                                         
        kb = kbot(i,kc)                                                         
        if (kt.eq.0) go to 60                                                   
          do 55 k=kb,kt                                                         
            omeg = 10. e0 *vvel(i,k)                                            
c....    if we will vv filter the clds, set rh to missing so its not            
c.         used in the frequency distribution....                               
            if(omeg.ge.vvcld(1,kc)) rh(i,k) = -1.                               
cc          vv= 0.0 e0                                                          
cc          if(omeg.lt.vvcld(1,kc)) vv=(vvcld(1,kc)-omeg)/                      
cc   1                                 (vvcld(1,kc)-vvcld(2,kc))                
cc          if(vv.gt.1.0 e0) vv = 1.0 e0                                        
cc          camt(i,k) = camt(i,k) * vv                                          
   55     continue                                                              
   60   continue                                                                
 1100 continue                                                                  
c--  now get maximum rh  amt layer for each type--l,m,h                         
c--   filter out clds which will be removed by vv                               
c                                                                               
c  rhmdl=-1 ... underground                                                     
c  rhmdl=-2 ... vertical motion filtered                                        
c                                                                               
      do 95 i=1,idim                                                            
        do 94 kcl=1,mcld                                                        
c                                                                               
          kt=ktop(i,kcl)                                                        
          kb=kbot(i,kcl)                                                        
          if(kt.eq.0) then                                                      
            rhmdl(i,kcl)=-1.                                                    
            nundg(kcl)=nundg(kcl)+1                                             
            go to 94                                                            
          endif                                                                 
ctune                                                                           
          if(usemax) then                                                       
            rrhmx = -0.000001                                                   
            krhlyr = 0                                                          
            do 1090 k=kb,kt                                                     
              if(rh(i,k).le.rrhmx) go to 1090                                   
              krhlyr = k                                                        
              rrhmx = rh(i,k)                                                   
 1090       continue                                                            
            if (krhlyr.gt.0) then                                               
              rhmdl(i,kcl) = rh(i,krhlyr)                                       
              nok(kcl)=nok(kcl)+1                                               
            else                                                                
              rhmdl(i,kcl)=-2.                                                  
              nvvel(kcl)=nvvel(kcl)+1                                           
            end if                                                              
c                                                                               
c  weighted mean                                                                
c                                                                               
          else                                                                  
            rhav=0.                                                             
            wei=0.                                                              
            imiss=1                                                             
            do k=kb,kt                                                          
              if(rh(i,k).ge.0.0) then                                           
                rhav=rhav+rh(i,k)*(si(k)-si(k+1))                               
                wei=wei+(si(k)-si(k+1))                                         
                imiss=imiss*0                                                   
              elseif(rh(i,k).lt.-0.5) then                                      
                imiss=imiss*1                                                   
              endif                                                             
            enddo                                                               
            if(imiss.eq.0) then                                                 
              if(wei.le.0.) then                                                
                print *,'sum of weight in the vertical is zero'                 
                call abort                                                      
              endif                                                             
              rhmdl(i,kcl)=rhav/wei                                             
              nok(kcl)=nok(kcl)+1                                               
            else                                                                
              rhmdl(i,kcl)=-2.                                                  
              nvvel(kcl)=nvvel(kcl)+1                                           
            endif                                                               
          endif                                                                 
c                                                                               
   94   continue                                                                
   95 continue                                                                  
ctune                                                                           
      kfir = 10                                                                 
      return                                                                    
      end                                                                       
