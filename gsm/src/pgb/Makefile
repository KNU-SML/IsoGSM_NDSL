EXECUTABLE = pgb.x
OBJS = bsslz1.o c2z.o dz2uv.o elat.o gau2ll.o gaulat.o getrh.o gradq.o gspc.o gtbits.o hydro.o iminv.o liftix.o llspc.o maxfac.o mpfdef.o nntprt.o omega.o panaly.o pgb1.o pgb.o pleg.o polext.o psynth.o psynthv.o rdsgh.o rdss.o rowsep.o sig2mw.o sig2p.o sig2tp.o spcoef.o spfmax.o spline.o sundry.o sunprm.o trss.o uv2dz.o wryte.o 

MACHINE = linux
SHELL = /bin/sh
CP = cp
RM = rm
CD = cd
CHMOD = chmod

CC = cc
CPP = gcc -E
F77 = ifort
FORT_FLAGS = -r8 -O2 -convert big_endian -shared-intel -mcmodel=medium
LOAD_FLAGS = -r8 -O2 -convert big_endian
INCLUDE_DIR = -I../.. -I../include -I/usr/local/mpi/intel18/mvapich2-2.2/include
LIB_DIR = /rhome/yeong/test/svrtest/libs/lib
INSTALL_DIR = /rhome/yeong/test/svrtest/gsm/bin
INSTALLDLIB = /rhome/yeong/test/svrtest/gsm/bin
EXTRA_LIBS = 
EXTRA_LIBS_PGB = 

LIBS = ${INSTALLDLIB}/share.a ${LIB_DIR}/w3lib.a ${LIB_DIR}/modelib.a \
		${EXTRA_LIBS_PGB}

#EXP_PARENT = @EXP_PARENT@
#.DEFAULT :
#	ln -fs ${EXP_PARENT}/src/pgb/$@ $@

${EXECUTABLE} : ${OBJS} ${LIBS}
	if [ linux = mac ] ; then \
		ranlib -s ${LIBS} ; \
	fi
	${F77} ${LOAD_FLAGS} -o $@ ${OBJS} ${LIBS} ${EXTRA_LIBS}
	${CP} $@ ${INSTALL_DIR}

depend :
	${LIB_DIR}/../etc/mkdep

.SUFFIXES : .o .f .F

.F.f :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
#	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f ; rm $*.i

.F.o :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
#	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f ; rm $*.i
	${F77} ${FORT_FLAGS} -c $*.f

.f.o :
	${F77} ${FORT_FLAGS} -c $*.f

clean :
	- ${RM} -f *.f *.o ${EXECUTABLE} core

dist_clean :
	- ${RM} -f *.f *.o *.x core Makefile
	- ${RM} -rf rii_files
