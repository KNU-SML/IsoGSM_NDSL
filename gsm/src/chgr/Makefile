EXECUTABLE = chgr.x
OBJS = bsslz1.o changers.o chgres.o cmpind.o dztouv.o epslon.o exmois.o exmoisn.o extrap.o fftidim.o fl22a.o gaulat.o gl2gl.o glats.o gozrmt.o greadi.o gwrite.o hbartr.o maklab.o maxmn.o msum2a.o newhyb.o newps.o newsig.o pln2t.o poly.o psum2a.o rmsgti.o rmsgto.o setakbk.o setsig.o sg2sg.o shufls.o shuflv.o sums2a.o sumv2a.o symasa.o transp.o trispl.o valts.o 

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
EXTRA_LIBS = 
EXTRA_LIBS_CHGR = 
IJCAP = 62
ILEVS = 28

LIBS = ${LIB_DIR}/w3lib.a ${LIB_DIR}/modelib.a ${EXTRA_LIBS_CHGR}

#EXP_PARENT = @EXP_PARENT@
#.DEFAULT :
#	ln -fs ${EXP_PARENT}/src/chgr/$@ $@

${EXECUTABLE} : ${OBJS} ${LIBS}
	if [ linux = mac ] ; then \
		ranlib -s ${LIBS} ; \
	fi
	${F77} ${LOAD_FLAGS} -o $@ ${OBJS} ${LIBS} ${EXTRA_LIBS}
	${CP} $@ ${INSTALL_DIR}/chgr.x
	${CP} $@ ${INSTALL_DIR}/chgr_t${IJCAP}k${ILEVS}.x ; ${RM} -f ${INSTALL_DIR}/chgr.x ; \
	ln -fs ${INSTALL_DIR}/chgr_t${IJCAP}k${ILEVS}.x ${INSTALL_DIR}/chgr.x

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
	- ${RM} -f *.f *.o ${EXECUTABLE} core ${INSTALL_DIR}/${EXECUTABLE}

dist_clean :
	- ${RM} -f *.f *.o *.x core Makefile
	- ${RM} -rf rii_files
