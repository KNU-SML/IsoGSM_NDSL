LIB = gsml_par.a

OBJS = accdia.o acritn.o aeros8.o aeros.o albaer.o albsno.o alcalc.o ambmsv.o ambmsv_hyb.o b10exps.o b10kdis.o bartri.o bdpkap.o bnmc.o canres.o cdate.o cfcexps.o cfckdis.o ch4exps.o ch4kdis.o cldadj.o cldjms.o cldjmsn.o cldprp.o cldprpi.o cloud.o cmpind.o cnvcldn.o co2exps.o co2kdis.o column.o comexps.o comkdis.o compjd.o conexps.o conrad.o coszmn.o crhtab.o csnow.o dbaeros8.o dbaeros.o dcyc2.o dellnp.o devap.o dzuvle.o epsilo.o epslon.o eqm_deg.o esatx.o fcstim.o fl2i.o flxco2.o freq.o frh2o.o frkin.o ftilonf.o funcdf.o funckt.o fwindx.o gaeros8.o gaeros.o gbphys.o gcljms.o getdia.o getrh.o gfdlrd.o gftlonf.o ggozrm.o glats.o gncpus.o gozrin.o gpln2i.o gradini.o gradpre.o grddf.o grdkt.o grlwsw.o grphys.o grrad1.o grrad2.o gwdps.o gwdrag.o h2oexps.o h2okdis.o hconst.o hrt.o hrtice.o hstep.o hydro.o iminv.o inddia.o insur2.o irrad1.o ismax.o ismin.o isorsv.o isrcheq.o isrchflt.o ktsoil.o liftix.o lrgscl.o lwrchou.o moninp.o mstadb.o mtntrq.o n2oexps.o n2okdis.o nopac.o o3int.o o3intn.o o3phys.o o3sbuv.o omega.o omegas.o osutm.o oz2d.o ozon2d.o penman.o pln2i.o poly.o progt2.o prtime.o putdia.o rascnp_1d.o rascnp.o rascnv.o rasv.o rdlwsw.o redprm.o rmsgt.o rncl.o rosr12.o row1ns.o seaicetm.o setes.o sethyb.o seto3p.o setras.o setrasv1_1d.o setsig.o sfcdiag.o sfcdif.o sfcdrv.o sfc_noah.o sfc_ocean.o sfc_osu.o sfc_seaice.o sflx.o shalcv.o shflx.o sig2mw.o sig2p.o sig2pt.o sig2tp.o sincoslat.o smflx.o snfrac.o snksrc.o snopac.o snow_new.o snowpack.o snowz0.o solar.o solir.o soliri.o soluv.o soluvi.o spcoef.o spfmax.o spline.o srt.o sstep.o subpgb.o sums2i.o sumtop.o sundry.o sunprm.o swflux.o swr95.o table.o tablup.o tbnd.o tdfcnd.o thsat.o tmpavg.o transo.o transp.o tridi2.o tridin.o twlt.o wdfcnd.o wridia.o wriflx.o wryte.o zerdia.o zerflx.o 

SRCS = ${OBJS:o=F}

MACHINE = linux
SHELL = /bin/sh
CP = cp
RM = rm
CD = cd
CHMOD = chmod

CC = cc
CPP = gcc -E
F77 = mpif90
FORT_FLAGS = -r8 -O2 -convert big_endian -I/usr/local/mpi/intel18/mvapich2-2.2/include
LOAD_FLAGS = -r8 -O2 -convert big_endian -L/opt/gm/lib
INCLUDE_DIR = -I../.. -I../include -I/usr/local/mpi/intel18/mvapich2-2.2/include
LIB_DIR = /rhome/yeong/test/svrtest/libs/lib
INSTALL_DIR = /rhome/yeong/test/svrtest/gsm/bin
JCAP = 62
LEVS = 28
NPES = 36

.DEFAULT :
	ln -fs ../gsml/$@ $@

${LIB} : ${SRCS} ${OBJS}
	${AR} -r $@ ${OBJS}
	${CP} $@ ${INSTALL_DIR}
	${CP} $@ ${INSTALL_DIR}/gsml_t${JCAP}k${LEVS}_n${NPES}.a ; \
		${RM} -f ${INSTALL_DIR}/gsml_par.a ; \
	ln -fs ${INSTALL_DIR}/gsml_t${JCAP}k${LEVS}_n${NPES}.a ${INSTALL_DIR}/gsml_par.a

${OBJS} : ${SRCS}

.SUFFIXES : .o .f .F

.F.f :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f ; rm $*.i

.F.o :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f ; rm $*.i
	${F77} ${FORT_FLAGS} -c $*.f
.f.o :
	${F77} ${FORT_FLAGS} -c $*.f

clean :
	- ${RM} -f *.f *.o .mod core *.a

dist_clean :
	- ${RM} -f *.f *.o *.mod *.a core Makefile
	- ${RM} -rf rii_files
