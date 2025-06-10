EXECUTABLE = cnvaer.x
OBJS = cnvaer.o glats.o poly.o 

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

LIBS = 

#EXP_PARENT = @EXP_PARENT@
#.DEFAULT :
#	ln -s ${EXP_PARENT}/src/cnvaer/$@ $@

${EXECUTABLE} : ${OBJS} ${LIBS}
	${F77} ${LOAD_FLAGS} -o $@ ${OBJS} ${LIBS}
	${CP} $@ ${INSTALL_DIR}

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
	- ${RM} -f *.f *.o *.x core

dist_clean :
	- ${RM} -f *.f *.o *.x core Makefile
	- ${RM} -rf rii_files
