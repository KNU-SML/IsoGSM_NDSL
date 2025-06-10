c
      parameter (nae=7,ndm=5,ncm1=6,ncm2=4,ncm=ncm1+ncm2,ncf=3)
      common/swaer99i/ idm(l,nae)
      common/swaer99/  haer(ndm,nae),dz(l),hh(l),hz(lp1), taur(l,nbd)
     1,              tauray(12),sigref(ndm,nae)
      common/opac01/ 
     1     ext0(ncm1,nbd),sca0(ncm1,nbd),ssa0(ncm1,nbd),estr(nbd)
     2,    aext(ncf,ncm2,nbd),bsca(ncf,ncm2,nbd),asf0(ncm1,nbd)
     3,    cssa(ncf,ncm2,nbd),dasf(ncf,ncm2,nbd),abpw(ncm2,nbd)
