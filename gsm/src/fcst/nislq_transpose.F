#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVH levh
#define LEVP1 levp1
#define LEVS levs
#define JCAP jcap
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define NCLDG ncldg
#define NPES npes
#define NROW nrow
#define NTRAC ntrac
#define NTOTAL ntotal
#ifdef MP
#define LEVSS levsp_
#else
#define LEVSS levs_
#endif
#define G g
#define PI pi


#include <define.h>
#ifdef NISLQ
           subroutine nislq_transpose_we2ns(a,b,levv,nsize)
c-------------------------------------------------------------------------------
c
c mpi transport from full dimension of west-east to full dimension of 
c north-south with latitude shuffl.
c
c-------------------------------------------------------------------------------
#include <paramodel.h>
#include <constant.h>
#include <comnislq.h>
#ifdef MP
#include <commpi.h>
#endif /* MP end */
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
        integer ::  levv,nsize    
        integer ::   i,j,k,n,mn,jj,lat1,lat2,ierr
        real   , dimension(lonfull,levv,latpart)     ::  a
        real   , dimension(latfull,levv,lonpart)     ::  b
#ifdef MP
        real(_mpi_real_),dimension(2*levv*mysllatlen*lonhalf) :: works
        real(_mpi_real_),dimension(2*levv*mysllonlen*lathalf) :: workr
#else
        real   , dimension(2*levv*mysllatlen*lonhalf)  ::  works
        real   , dimension(2*levv*mysllonlen*lathalf)  ::  workr
#endif
        integer, dimension(nsize)  ::  lensend,lenrecv         
        integer, dimension(nsize)  ::  locsend,locrecv
c
           mn=0
           do n = 1,nsize 
             locsend(n)=mn
             do j = 1,mysllatlen
               do k = 1,levv
                 do i = sllonstr(n),sllonstr(n)+sllonlen(n)-1
                   mn = mn + 1
                   works(mn)=a(i        ,k,j)
                 enddo
                 do i = sllonstr(n),sllonstr(n)+sllonlen(n)-1
                   mn = mn + 1
                   works(mn)=a(i+lonhalf,k,j)
                 enddo
               enddo
             enddo
             lensend(n)=mn-locsend(n)
           enddo
c
           mn=0
           do n = 1,nsize
             locrecv(n)=mn
             lenrecv(n)=levv*sllatlen(n)*mysllonlen*2
             mn=mn+lenrecv(n)
           enddo
c
#ifdef MP
         call mpi_barrier (comm_column,ierr)
         call mpi_alltoallv(works,lensend,locsend,MPI_REAL
     &,    workr,lenrecv,locrecv,MPI_REAL,comm_column,ierr)
        
#define WORKR workr
#else
#define WORKR works
#endif /* MP end */
           mn=0
           jj=0
           do n = 1,nsize
             do j = 1,sllatlen(n)
               jj = jj + 1
               lat1=truej(jj)
               lat2=latfull+1-lat1
               do k = 1,levv
                 do i = 1,mysllonlen
                   mn = mn + 1
                   b(lat1,k,i) = WORKR(mn)
                 enddo
                 do i = 1,mysllonlen
                   mn = mn + 1
                   b(lat2,k,i) = WORKR(mn)
                 enddo
                enddo
             enddo
           enddo
         
c
#undef WORKR
           return
           end subroutine nislq_transpose_we2ns
#endif /* NISLQ end */
c-------------------------------------------------------------------------------
c
c-------------------------------------------------------------------------------
#ifdef NISLQ
           subroutine nislq_transpose_ns2we(a,b,levv,nsize)
c-------------------------------------------------------------------------------
c
c mpi transport from full dimension of west-east to full dimension of 
c north-south with latitude shuffl.
c
c-------------------------------------------------------------------------------
#include <paramodel.h>
#ifdef MP
#include <commpi.h>
#endif /* MP end */
#include <comnislq.h>
c-------------------------------------------------------------------------------
           integer     ::  levv,nsize                  
           integer     ::  i,j,k,n,mn,jj,lat1,lat2,ierr
           real   , dimension(latfull,levv,lonpart)     ::  a
           real   , dimension(lonfull,levv,latpart)     ::  b
#ifdef MP
           real(_mpi_real_),dimension(2*levv*mysllonlen*lathalf)::works
           real(_mpi_real_),dimension(2*levv*mysllatlen*lonhalf)::workr
#else
           real   , dimension(2*levv*mysllonlen*lathalf)  ::  works
           real   , dimension(2*levv*mysllatlen*lonhalf)  ::  workr
#endif
           integer, dimension(nsize)    ::  lensend,lenrecv          
           integer, dimension(nsize)    ::  locsend,locrecv
c
           mn=0
           jj=0
           do n = 1,nsize
             locsend(n)=mn
             do j = 1,sllatlen(n)
               jj=jj+1
               lat1=truej(jj)
               lat2=latfull+1-lat1
               do k = 1,levv
                 do i = 1,mysllonlen
                   mn = mn + 1
                   works(mn) = a(lat1,k,i)
                 enddo
                 do i = 1,mysllonlen
                   mn = mn + 1
                   works(mn) = a(lat2,k,i)
                 enddo
               enddo
             enddo
             lensend(n)=mn-locsend(n)
           enddo
c
           mn=0
           do n = 1,nsize
             locrecv(n)=mn
             lenrecv(n)=levv*mysllatlen*sllonlen(n)*2
             mn=mn+lenrecv(n)
           enddo
c
#ifdef MP
           call mpi_barrier (comm_column,ierr)
           call mpi_alltoallv(works,lensend,locsend,MPI_REAL
     &,               workr,lenrecv,locrecv,MPI_REAL,comm_column,ierr)
#define WORKR workr
#else
#define WORKR works
#endif /* MP end */
c
           mn=0
           do n = 1,nsize 
             do j = 1,mysllatlen
               do k = 1,levv
                 do i = sllonstr(n),sllonstr(n)+sllonlen(n)-1
                   mn = mn + 1
                   b(i        ,k,j) = WORKR(mn)
                 enddo
                 do i = sllonstr(n),sllonstr(n)+sllonlen(n)-1
                   mn = mn + 1
                   b(i+lonhalf,k,j) = WORKR(mn)
                 enddo
               enddo
             enddo
           enddo
#undef WORKR
c
           return
           end subroutine nislq_transpose_ns2we
c-------------------------------------------------------------------------------
#endif /* NISLQ end */
