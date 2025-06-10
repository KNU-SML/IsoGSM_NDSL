#define LEVS levs
      subroutine co2int(itape,t15a,t15b,t22,ratio,ir,nmethd)                    
cnov89                                                                          
c     *********************************************************                 
c       changes to data read  and format see co222     ***                      
c          ..... k.campana   march 1988,october 1988                            
c       changes to pass itape,and if ir=4,read 1 co2 rec..kac nov89             
c     *********************************************************                 
c       co2int interpolates carbon dioxide transmission functions               
c  from the 109 level grid,for which the transmission functions                 
c  have been pre-calculated, to the grid structure specified by the             
c  user.                                                                        
c                                                                               
c        method:                                                                
c                                                                               
c      co2int is employable for two purposes: 1) to obtain transmis-            
c  sivities between any 2 of an array of user-defined pressures; and            
c  2) to obtain layer-mean transmissivities between any 2 of an array           
c  of user-defined pressure layers.to clarify these two purposes,see            
c  the diagram and discussion below.                                            
c      co2int may be used to execute only one purpose at one time.              
c                                                                               
c     let p be an array of user-defined pressures                               
c     and pd be user-defined pressure layers.                                   
c                                                                               
c       - - - - - - - - -   pd(i-1) ---                                         
c                                     ^                                         
c       -----------------   p(i)      ^  pressure layer i  (plm(i))             
c                                     ^                                         
c       - - - - - - - - -   pd(i)  ---                                          
c                                     ^                                         
c       -----------------   p(i+1)    ^  pressure layer i+1 (plm(i+1))          
c                                     ^                                         
c       - - - - - - - - -   pd(i+1)---                                          
c            ...                          (the notation used is                 
c            ...                          consistent with the code)             
c            ...                                                                
c      - - - - - - - - -    pd(j-1)                                             
c                                                                               
c      -----------------    p(j)                                                
c                                                                               
c      - - - - - - - - -    pd(j)                                               
c                                                                               
c      purpose 1:   the transmissivity between specific pressures               
c      p(i) and p(j) ,tau(p(i),p(j))  is computed by this program.              
c      in this mode,there is no reference to layer pressures pd                 
c      (pd,plm are not inputted).                                               
c                                                                               
c      purpose 2:   the layer-mean transmissivity between a layer-              
c      mean pressure plm(j) and pressure layer i is given by                    
c         taulm(plm(i),plm(j)). it is computed by the integral                  
c                                                                               
c                           pd(i)                                               
c                           ----                                                
c             1             ^                                                   
c        -------------  *   ^   tau ( p',plm(j) )  dp'                          
c        pd(i)-pd(i-1)      ^                                                   
c                        ----                                                   
c                        pd(i-1)                                                
c                                                                               
c           the layer-mean pressure plm(i) is specified by the user.            
c        for many purposes,plm will be chosen to be the average                 
c        pressure in the layer-ie,plm(i)=0.5*(pd(i-1)+pd(i)).                   
c           for layer-mean transmissivities,the user thus inputs                
c        a pressure array (pd) defining the pressure layers and an              
c        array (plm) defining the layer-mean pressures.the calculation          
c        does not depend on the p array used for purpose 1 (p is not            
c        inputted).                                                             
c                                                                               
c            the following paragraphs depict the utilization of this            
c       code when used to compute transmissivities between specific             
c       pressures. later paragraphs describe additional features needed         
c       for layer-mean transmissivities.                                        
c                                                                               
c          for a given co2 mixing ratio and standard temperature                
c      profile,a table of transmission functions for a fixed grid               
c     of atmospheric pressures has been pre-calculated.                         
c      the standard temperature profile is computed from the us                 
c     standard atmosphere (1977) table.additionally, the                        
c     same transmission functions have been pre-calculated for a                
c     temperature profile increased and decreased (at all levels)               
c     by 25 degrees.                                                            
c         this program reads in the prespecified transmission functions         
c     and a user-supplied pressure grid (p(i)) and calculates trans-            
c     mission functions ,tau(p(i),p(j)), for all p(i) s and p(j) s.             
c     a logarithmic interpolation scheme is used.                               
c         this method is repeated for the three temperature profiles            
c     given above .therefore outputs from the program are three tables          
c     of transmission functions for the user-supplied pressure grid.            
c     the existence of the three tables permits subsequent interpo-             
c     lation to a user-supplied temperature profile using the method            
c     described in the reference.see limitations section if the                 
c     user desires to obtain only 1 table of transmissivities.                  
c                                                                               
c     modifications for layer-mean transmissivities:                            
c          the pressures inputted are the layer-mean pressures,pd,              
c     and the layer-mean pressures ,plm. a series of transmissivities           
c     (tau(p'',plm(j)) are computed and the integral given in the               
c     discussion of purpose 2 is computed.for plm(i) not equal to               
c     plm(j) simpson s rule is used with 5 points. if plm(i)=plm(j)             
c     (the -nearby layer- case) a 49-point quadrature is used for               
c     greater accuracy.the output is in taulm(plm(i),plm(j)).                   
c        note:                                                                  
c     taulm is not a symmetrical matrix. for the array element                  
c     taulm(plm(i),plm(j)),the inner(first,most rapidly varying)                
c     dimension is the varying layer-mean pressure,plm(i);the outer             
c     (second) dimension is the fixed layer-mean pressure plm(j).               
c     thus the element taulm(2,3) is the transmission function between          
c     the fixed pressure plm(3)  and the pressure layer having an averag        
c     pressure of plm(2).                                                       
c         also note that no quadrature is performed over the layer              
c     between the smallest nonzero pressure and zero pressure;                  
c     taulm is taulm(0,plm(j)) in this case,and taulm(0,0)=1.                   
c                                                                               
c                                                                               
c             reference:                                                        
c         s.b.fels and m.d.schwarzkopf,-an efficient accurate                   
c     algorithm for calculating co2 15 um band cooling rates-,journal           
c     of geophysical research,vol.86,no. c2, pp.1205-1232,1981.                 
c        modifications to the algorithm have been made by the authors;          
c     contact s.b.f.or m.d.s. for further details.a note to j.g.r.              
c     is planned to document these changes.                                     
c                                                                               
c            author:    m.daniel schwarzkopf                                    
c                                                                               
c            date:      14 july 1983                                            
c                                                                               
c            address:                                                           
c                                                                               
c                      g.f.d.l.                                                 
c                      p.o.box 308                                              
c                      princeton,n.j.08540                                      
c                      u.s.a.                                                   
c            telephone:  (609) 452-6521                                         
c                                                                               
c            information on tape: this source is the first file                 
c        on this tape.the six files that follow are co2 trans-                  
c        missivities for the 500-850 cm-1 interval for co2                      
c        concentrations of 330 ppmv (1x) ,660 ppmv (2x), and                    
c        1320 ppmv (4x). the files are arranged as follows:                     
c          file 2   1x,consolidated using b(250) weighting fctn.                
c          file 3   1x,consolidated with no weighting fctn.                     
c          file 4   2x,consolidated using b(250) weighting fctn.                
c          file 5   2x,consolidated with no weighting fctn.                     
c          file 6   4x,consolidated using b(250) weighting fctn.                
c          file 7   4x,consolidated with no weighting fctn.                     
c            files 2,4,6 are recommended for use in obtaining                   
c        transmission functions for use in heating rate                         
c        computations;they correspond to the transmissivities                   
c        discussed in the 1980 paper.files 3,5,7 are provided                   
c        to facilitate comparison with observation and with other               
c        calculations.                                                          
c                                                                               
c            program language: fortran 1977,including parameter                 
c        and program statements.the program is written on a                     
c        cyber 170-730.see the section on limitations for                       
c        adaptations to other machines.                                         
c                                                                               
c          input units,formats and format statement nos:                        
c                                                                               
c   unit no    variables       format      statement no.    type                
c      5        p (purpose 1)  (5e16.9)        201         cards                
c      5        pd (purpose 2) (5e16.9)        201         cards                
c      5        plm(purpose 2) (5e16.9)        201         cards                
c      5        nmethd         (i3)            202         cards                
c      20       transa         (4f20.14)       102          tape                
cnov89                                                                          
c      itape    transa         (4f20.14)       102          tape                
cnov89                                                                          
c                                                                               
c         output units,formats and format statement nos:                        
c                                                                               
c   unit no    variables       format     statement no.                         
c      6         trnfct        (1x,8f15.8)     301         print                
c      22        trnfct        (4f20.14)       102          tape                
c                                                                               
c            parameter inputs:                                                  
c     a) nlevls    : nlevls is an (integer) parameter denoting                  
c        the number of nonzero pressure levels for purpose 1                    
c        or the number of nonzero layer pressures needed to                     
c        specify the pressure layers(purpose 2) in the output                   
c        grid. for example,in purpose 1,if p=0,100,1000,nlevls=2.               
c        if,in purpose 2,pd=0,100,500,1000,the number of nonzero                
c        pressure layers=2,so nlevls=2                                          
c           in the code as written,nlevls=40; the user should                   
c        change this value to a user-specified value.                           
c     b) nlp1,nlp2 : integer parameters defined as: nlp1=nlevls+1;              
c        nlp2=nlevls+2.                                                         
c           see limitations for code modifications if parameter                 
c        statements are not allowed on your machine.                            
c                                                                               
c            inputs:                                                            
c                                                                               
c     a) transa    : the 109x109 grid of transmission functions                 
c            transa is a  double precision real array.                          
c                                                                               
c           transa  is read from file 20. this file contains 3                  
c     records,as follows:                                                       
c        1)   transa, standard temperature profile                              
c        3)   transa, standard temperatures + 25 deg                            
c        5)   transa, standard temperatures - 25 deg                            
c                                                                               
c     b)   nmethd: an integer whose value is either 1 (if co2int is             
c       to be used for purpose 1) or 2 (if co2int is to be used for             
c       purpose 2).                                                             
c                                                                               
c     c)     p,pd,plm :                                                         
c          p is a real array (length nlp1) specifying the pressure              
c       grid at which transmission functions are to be computed for             
c       purpose 1.the dimension  of p is  in millibars.the                      
c       following limitations will be explained more                            
c       in the section on limitations: p(1) must be zero; p(nlp1),the           
c       largest pressure, must not exceed 1165 millibars.                       
c         pd is a real array (length nlp2) specifying the pressure              
c       layers for which layer-averaged transmission functions are              
c       to be computed.the dimension of pd is millibars.the limitations         
c       for pd are the same as for p,and are given in the section on            
c       limitations.                                                            
c         plm is a real array (length nlp2) specifying the layer-mean           
c       pressures. the dimension of plm is millibars. the limitations           
c       for plm are the same as for p,and are given in the section on           
c       limitations.pd is read in before plm.                                   
c                                                                               
c          note: again,we note that the user will input either p (for           
c       purpose 1) or pd and plm(for purpose 2) but not both.                   
c                                                                               
c                                                                               
c                                                                               
c                                                                               
c           limitations:                                                        
c     1)       p(1)=0.,pd(1)=0.,plm(1)=0. the top pressure level                
c       must be zero,or the top pressure layer must be bounded by zero.         
c       the top layer-mean pressure (plm(1)) must be zero; no                   
c       quadrature is done on the top pressure layer.even if one is             
c       not interested in the transmission function between 0 and p(j),         
c       one must include such a level.                                          
c     2)      pd(nlp2)=p(nlp1) is less than or equal to 1165 mb.                
c       extrapolation to higher pressures is not possible.                      
c     3)      if program is not permitted on your compiler,                     
c       simply delete the line.                                                 
c     4)      if parameter is not permitted,do the following:                   
c            1) delete all parameter statements in co2int                       
c            2) at the point where nmethod is read in,add:                      
c                read (5,202) nlevls                                            
c                nlp1=nlevls+1                                                  
c                nlp2=nlevls+2                                                  
c            3) change dimension and/or common statements defining              
c              arrays trns,delta,p,pd,trnfct,ps,pds,plm in co2int.              
c              the numerical value of (nlevls+1) should be inserted             
c              in dimension or common statements for trns,delta,                
c              p,trnfct,ps,plm; the numerical value of (nlevls+2)               
c              in dimension or common statements for pd,pds.                    
c      5)    parameter (nlevls=40) and the other parameter                      
c       statements are written in cdc fortran; on other machines the            
c       same statement may be written differently,for example as                
c       parameter   nlevls=40                                                   
c      6) -double precision- is used instead of -real*8- ,due to                
c       requirements of cdc fortan.                                             
c      7) the statement -do 400 kkk=1,3- controls the number of                 
c       transmissivity output matrices porduced by the program.to               
c       produce 1 output matrix,delete this statement.                          
c                                                                               
c     output:                                                                   
c         a) trnfct is an (nlp1,nlp1) real array of the transmission            
c     functions appropriate to your array. it is to be saved on file 22.        
c     the procedure for saving may be modified; as given here,the               
c     output is in card image form with a format of (4f20.14).                  
c                                                                               
c         b)  printed  output is a listing of trnfct on unit 6, in              
c     the format (1x,8f15.8) (format statement 301). the user may               
c     modify or eliminate this at will.                                         
c                                                                               
c      ************   function interpolater routine  *****************          
c                                                                               
c                                                                               
c     ******   the following parameter gives the number of     *******          
c     ******           data levels in the model                *******          
c     ****************************************************************          
#include <paramodel.h>
      parameter          (nlevls=levs_)                                         
c     ****************************************************************          
      parameter          (nlp1=nlevls+1,nlp2=nlevls+2)                          
      common/input/p1,p2,trnslo,ia,ja,n                                         
      common/press/pa(109)                                                      
      common/tran/ transa(109,109)                                              
      common / output / trns(nlp1,nlp1)                                         
      common/inputp/p(nlp1),pd(nlp2)                                            
      dimension ps(nlp1),pds(nlp2),plm(nlp1)                                    
      dimension nrtab(3)                                                        
      dimension t15a(nlp2,2),t15b(nlp1)                                         
      dimension t22(nlp1,nlp1,3)                                                
      data nrtab/1,2,4/                                                         
c***********************************                                            
c   the following are the input formats                                         
100   format (4f20.14)                                                          
201   format (5e16.9)                                                           
202   format (i3)                                                               
co222   203   format (f12.6,i2)                                                 
203   format (f12.6)                                                            
c    the following are the output formats                                       
102   format (4f20.14)                                                          
301   format (1x,8f15.8)                                                        
c                                                                               
ccc   rewind 15                                                                 
ccc   rewind 20                                                                 
cnov89                                                                          
      rewind itape                                                              
cnov89                                                                          
ccc   rewind 22                                                                 
c                                                                               
c     calculation of pa -the -table- of 109 grid pressures                      
c     note-this code must not be changed by the user^^^^^^^^^                   
      pa(1)=0.                                                                  
      fact15=10.**(1./15.)                                                      
      fact30=10.**(1./30.)                                                      
      pa(2)=1.0e-3                                                              
      do 231 i=2,76                                                             
      pa(i+1)=pa(i)*fact15                                                      
231   continue                                                                  
      do 232 i=77,108                                                           
      pa(i+1)=pa(i)*fact30                                                      
232   continue                                                                  
c                                                                               
      n=25                                                                      
      nlv=nlevls                                                                
      nlp1v=nlp1                                                                
      nlp2v=nlp2                                                                
c     read in the co2 mixing ratio(in units of 330 ppmv),and an index           
c     giving the frequency range of the lbl data                                
co222    read (5,203) ratio,ir                                                  
ccc         ir = 1                                                              
ccc         read (5,203) ratio                                                  
co222   ***********************************                                     
c***values for ir*****                                                          
c          ir=1     consol. lbl trans. =490-850                                 
c          ir=2     consol. lbl trans. =490-670                                 
c          ir=3     consol. lbl trans. =670-850                                 
c          ir=4     consol. lbl trans. =2270-2380                               
c*** ir must be 1,2,3 or 4 for the pgm. to work                                 
c     also read in the method no.(1 or 2)                                       
ccc         read (5,202) nmethd                                                 
coct92if (ratio.eq.1.0) go to 621                                               
coct92if (ratio.ge.2.0) go to 622                                               
coct92if (ratio.eq.4.0) go to 623                                               
coct92if (ratio.gt.1.0.and.ratio.lt.2.0) go to 624                              
coct92if (ratio.gt.2.0.and.ratio.lt.4.0) go to 625                              
coct92stop 8746                                                                 
cnov89  621   itap1=20                                                          
621   itap1=itape                                                               
cnov89                                                                          
      ntap=1                                                                    
      go to 630                                                                 
622   itap1=21                                                                  
      ntap=1                                                                    
      go to 630                                                                 
623   itap1=22                                                                  
      ntap=1                                                                    
      go to 630                                                                 
cnov89   624   itap1=20                                                         
624   itap1=itape                                                               
cnov89                                                                          
      ntap=2                                                                    
      ratstd=2.0                                                                
      ratsm=1.0                                                                 
      go to 630                                                                 
625   itap1=21                                                                  
      ntap=2                                                                    
      ratstd=4.0                                                                
      ratsm=2.0                                                                 
630   continue                                                                  
      if (nmethd.eq.2) go to 502                                                
c   *****cards for purpose 1(nmethd=1)                                          
ccc         read (15,201) (p(i),i=1,nlp1)                                       
      do 300 i=1,nlp1                                                           
        p(i)=t15b(i)                                                            
  300 continue                                                                  
      do 801 i=1,nlp1                                                           
      ps(i)=p(i)                                                                
801   continue                                                                  
      go to 503                                                                 
502   continue                                                                  
c  *****cards for purpose 2(nmethd=2)                                           
ccc         read (15,201) (pd(i),i=1,nlp2)                                      
ccc         read (15,201) (plm(i),i=1,nlp1)                                     
      do 303 i=1,nlp2                                                           
        pd(i)=t15a(i,1)                                                         
  303 continue                                                                  
      do 302 i=1,nlp1                                                           
        plm(i)=t15a(i,2)                                                        
  302 continue                                                                  
      do 802 i=1,nlp1                                                           
      pds(i)=pd(i+1)                                                            
      ps(i)=plm(i)                                                              
802   continue                                                                  
c                                                                               
503   continue                                                                  
c  *****do loop controlling number of output matrices                           
cnov89                                                                          
cnov89    do 400 kkk=1,3                                                        
      icloop = 3                                                                
      if (ir.eq.4) icloop = 1                                                   
      do 400 kkk=1,icloop                                                       
cnov89                                                                          
      if (ntap.eq.2) call rctrns(itap1,ratstd,ratsm,ratio)                      
c  **********************                                                       
      if (nmethd.eq.2) go to 505                                                
c   *****cards for purpose 1(nmethd=1)                                          
      do 803 i=1,nlp1                                                           
      p(i)=ps(i)                                                                
803   continue                                                                  
      go to 506                                                                 
505   continue                                                                  
c  *****cards for purpose 2(nmethd=2)                                           
      do 804 i=1,nlp1                                                           
      pd(i)=pds(i)                                                              
      p(i)=ps(i)                                                                
804   continue                                                                  
c                                                                               
506   continue                                                                  
      ia=108                                                                    
      iap=ia+1                                                                  
cnov89   if (ntap.eq.1) read (20,100) ((transa(i,j),i=1,109),j=1,109)           
      if (ntap.eq.1) read (itape,100) ((transa(i,j),i=1,109),j=1,109)           
cnov89                                                                          
      do 4 i=1,iap                                                              
      transa(i,i)=1.0                                                           
    4 continue                                                                  
      call coeint(ratio,ir)                                                     
      do 805 i=1,nlp1                                                           
      do 805 j=1,nlp1                                                           
      trns(j,i)=1.00                                                            
805   continue                                                                  
      do 10 i=1,nlp1                                                            
      do 20 j=1,i                                                               
      if (i.eq.j) go to 20                                                      
      p1=p(j)                                                                   
      p2=p(i)                                                                   
      call sintr2                                                               
      trns(j,i)=trnslo                                                          
20    continue                                                                  
10    continue                                                                  
      do 47 i=1,nlp1                                                            
      do 47 j=i,nlp1                                                            
      trns(j,i)=trns(i,j)                                                       
47    continue                                                                  
c  *****this is the end of purpose 1 calculations                               
      if (nmethd.eq.1) go to 2872                                               
c                                                                               
      do 51 j=1,nlp1                                                            
      do 52 i=2,nlp1                                                            
      ia=i                                                                      
      ja=j                                                                      
      n=25                                                                      
      if (i.ne.j) n=3                                                           
      call quadsr(nlv,nlp1v,nlp2v,p,pd,trns)                                    
52    continue                                                                  
51    continue                                                                  
c  *****this is the end of purpose 2 calculations                               
2872  continue                                                                  
c                                                                               
      write (6,301) ((trns(i,j),i=1,nlp1),j=1,nlp1)                             
ccc         write (22,102) ((trns(i,j),i=1,nlp1),j=1,nlp1)                      
      do 304 j=1,nlp1                                                           
       do 304 i=1,nlp1                                                          
        t22(i,j,kkk) = trns(i,j)                                                
  304 continue                                                                  
400   continue                                                                  
      return                                                                    
      end                                                                       
