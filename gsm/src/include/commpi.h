c
      include 'mpif.h'
      integer npes,ncol,nrow,mype,master,msgtag
     *       ,myrow,mycol,comm_row,comm_column
     *       ,levstr,levlen,lonstr,lonlen,latstr,latlen
     *       ,lwvstr,lwvlen,lntstr,lntlen,lnpstr,lnplen
     *       ,lwvdef,latdef
#ifdef COUPLE
     *       ,mpi_comm_gsm
#endif
      common /commpi/ npes,ncol,nrow,mype,master,msgtag
     *               ,myrow,mycol,comm_row,comm_column
     *               ,levstr(0:npes_-1),levlen(0:npes_-1)
     *               ,lonstr(0:npes_-1),lonlen(0:npes_-1)
     *               ,latstr(0:npes_-1),latlen(0:npes_-1)
     *               ,lwvstr(0:npes_-1),lwvlen(0:npes_-1)
     *               ,lntstr(0:npes_-1),lntlen(0:npes_-1)
     *               ,lnpstr(0:npes_-1),lnplen(0:npes_-1)
     *               ,lwvdef(jcap_+1)  ,latdef(latg_/2)
#ifdef COUPLE
     *               ,mpi_comm_gsm
#endif
