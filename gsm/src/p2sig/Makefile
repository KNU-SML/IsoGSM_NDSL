EXECUTABLE = p2sig.x
OBJS = cmpind.o exmois2.o fftidim.o fixrd2.o fl22a.o gozrmt.o hbartr.o msum2a.o newps2.o newps3.o newsig2.o p2sig.o pln2t.o prsgrbrd.o psum2a.o rmsgto.o row1ns2.o rowsep2.o setgrd.o sg2sg.o sg2sgrh.o shalfi.o shalfo.o symasa.o todxdy.o 

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

LIBS = ${INSTALL_DIR}/sfcl.a ${LIB_DIR}/w3lib.a ${LIB_DIR}/modelib.a ${INSTALL_DIR}/share.a ${EXTRA_LIBS_CHGR}

#EXP_PARENT = @EXP_PARENT@
#.DEFAULT :
#	ln -fs ${EXP_PARENT}/src/chgr/$@ $@

${EXECUTABLE} : ${OBJS} ${LIBS}
	if [ linux = mac ] ; then \
		ranlib -s ${LIBS} ; \
	fi
	${F77} ${LOAD_FLAGS} -o $@ ${OBJS} ${LIBS} ${EXTRA_LIBS}
	${CP} $@ ${INSTALL_DIR}/p2sig.x

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
	- ${RM} -f *.f *.o ${EXECUTABLE} core ${INSTALL_DIR}/${EXECUTABLE}

dist_clean :
	- ${RM} -f *.f *.o *.x core Makefile
	- ${RM} -rf rii_files
