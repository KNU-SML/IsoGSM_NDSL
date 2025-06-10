EXECUTABLE = mtn.x
OBJS = equallat.o gausslat.o legendre.o makemt.o makeoa.o mtn.o nntprt.o sphert.o tersub.o 

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
EXTRA_LIBS_MTN = 

LIBS = ${LIB_DIR}/modelib.a ${EXTRA_LIBS_MTN}

#EXP_PARENT = @EXP_PARENT@
#.DEFAULT :
#	ln -fs ${EXP_PARENT}/src/mtn/$@ $@

${INSTALL_DIR}/.mtn : ${EXECUTABLE}
	 ./mtn.sh || exit 8

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
	- ${RM} -f *.f *.o *.x core ${INSTALL_DIR}/.mtn
	- ${RM} -f top*.asc

dist_clean :
	- ${RM} -f *.f *.o *.x core Makefile ${INSTALL_DIR}/.mtn
	- ${RM} -f top*.asc
	- ${RM} -rf rii_files
