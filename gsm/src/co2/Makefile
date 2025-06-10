EXECUTABLE = co2.x
OBJS = antemp.o co2in1.o co2ins.o co2int.o co2o3.o co2ptz.o coeint.o mixrat.o newsig.o ozone.o path.o qintrp.o quadsr.o rctrns.o sigp.o sintr2.o 

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
#	ln -fs ${EXP_PARENT}/src/co2/$@ $@

${INSTALL_DIR}/.co2 : ${EXECUTABLE}
	./co2.sh || exit 8

${EXECUTABLE} : ${OBJS} ${LIBS}
	${F77} ${LOAD_FLAGS} -o $@ ${OBJS} ${LIBS}
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
	- ${RM} -f *.f *.o *.x core ${INSTALL_DIR}/.co2

dist_clean :
	- ${RM} -f *.f *.o *.x core Makefile ${INSTALL_DIR}/.co2
	- ${RM} -rf rii_files
