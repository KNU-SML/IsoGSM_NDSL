EXECUTABLE = fcst_par.x

OBJS = basefnam.o cloud2.o cloud3.o cloud5.o cloud6.o cof2grd.o dampux.o deldfsp.o deldif.o dfini.o dfkt.o filtr1.o filtr2.o fl2ip.o flpflm.o gau2ll.o getcon.o getij.o gfidi_hyb.o gfidiu.o gloopa.o gloopb.o gloopp.o gloopr.o greduceg.o gridiff.o gsicdf.o gsmend.o gsm.o gsmini.o gsmstep.o impadj.o ipwrcn.o isc8.o kenput.o mltdia.o mltflx.o nislq.o nislq_mass_advect.o nislq_mono_advect.o nislq_transpose.o progtm99.o rdsig.o reduceintp.o reducepick.o rowsep.o sicdif.o sph_fft_trans.o uvsums.o uvtodz.o wri_cldrad.o wrisig.o wriznl.o zertmm.o 

SRCS = ${OBJS:.o=.F}

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
INSTALLDLIB = /rhome/yeong/test/svrtest/gsm/bin
INSTALL_LIB = /rhome/yeong/test/svrtest/gsm/bin
EXTRA_LIBS = 
EXTRA_LIBS_FCST = 
JCAP = 62
LEVS = 28
NPES = 36

LIBS = ${INSTALLDLIB}/gsml_par.a ${INSTALLDLIB}/sfcl_par.a ${INSTALLDLIB}/mpi.a ${INSTALLDLIB}/share.a ${LIB_DIR}/modelib.a ${LIB_DIR}/w3lib.a ${EXTRA_LIBS_FCST}

.DEFAULT :
	ln -fs ../fcst/$@ $@

${EXECUTABLE} : ${OBJS} ${LIBS}
	if [ linux = mac ] ; then \
		ranlib -s ${LIBS} ; \
	fi
	${F77} ${LOAD_FLAGS} -o $@ ${OBJS} ${LIBS} ${EXTRA_LIBS}
	${CP} $@ ${INSTALL_DIR}/fcst_t${JCAP}k${LEVS}_n${NPES}.x ; ${RM} -f ${INSTALL_DIR}/fcst_par.x ; \
ln -fs ${INSTALL_DIR}/fcst_t${JCAP}k${LEVS}_n${NPES}.x ${INSTALL_DIR}/fcst_par.x

${OBJS} : ${SRCS}

depend :
	${LIB_DIR}/../etc/mkdep

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
	- ${RM} -f *.f *.o .mod ${EXECUTABLE} core

dist_clean :
	- ${RM} -f *.f *.o *.mod *.x core Makefile
	- ${RM} -rf rii_files
