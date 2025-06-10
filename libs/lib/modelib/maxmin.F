#define NPES npes
      subroutine maxmin(a,len,k,k1,k2,ch)                                       
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:  maxmin                                                           
c   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06          
c                                                                               
c abstract:  do print maximum and minimum of a given array.                     
c                                                                               
c program history log:                                                          
c                                                                               
c usage:    call  maxmin(a,len,k,k1,k2,ch)                                      
c   input argument list:                                                        
c     a 	- array for computing max and min (len,k)                              
c     len	- the first dimension of a                                            
c     k 	- the second dimension of a                                            
c     k1	- lower limit of second dimension to print                             
c     k2	- upper limit to print                                                 
c     ch 	- charcter string to print                                            
c                 fpl and fml                                                   
c                                                                               
c   output argument list:                                                       
c                                                                               
c   input files: none                                                           
c                                                                               
c   output files:                                                               
c     standard output                                                           
c                                                                               
c   subprograms called:                                                         
c     intrinsic functions: amax1 amin1                                          
c                                                                               
c   remark: none                                                                
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray ymp.                                                         
c                                                                               
c$$$ 
      dimension a(len,k) 
      character ch*(*)
c                    
      do j=k1,k2
      aamax = a(1,j)
      aamin = a(1,j)
      do m=1,len
      aamax = max( aamax, a(m,j) )
      aamin = min( aamin, a(m,j) )
      enddo
      if (k1.ne.k2.or.k1.gt.1) then
      print 100,ch,': max=',aamax,' min=',aamin,' at k=',j
  100 format(A12,A6,E12.4,A5,E12.5,A6,I3)
      else
      print 101,ch,': max=',aamax,' min=',aamin
  101 format(A12,A6,E12.4,A5,E12.5)
      endif
      enddo
      return
      end
