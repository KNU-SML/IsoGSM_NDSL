LIB = share.a
OBJS = cmpind.o count.o dztouv.o elat.o epslon.o fnam.o fpkap.o fpkapx.o fpvs0.o fpvs.o fpvsq0.o fpvsq.o fpvsx0.o fpvsx.o ftdp.o ftdpq.o ftdpx.o ftdpxg.o fthe.o ftheq.o fthex.o ftlcl.o ftma.o ftmaq.o ftmax.o ftmaxg.o gozrmt.o gpvs.o gribit2.o gribit.o gtdp.o gthe.o gtma.o idsdef.o isrcheqx.o isrchfltx.o isrchnex.o psynth.o psynthv.o qcsli.o row1ns.o rowsep.o 

MACHINE = linux
SHELL = /bin/sh
AR = ar
CP = cp
RM = rm
CD = cd
CHMOD = chmod

CC = cc
CPP = gcc -E
F77 = ifort
FORT_FLAGS = -r8 -O2 -convert big_endian -shared-intel -mcmodel=medium
LOAD_FLAGS = -r8 -O2 -convert big_endian -shared-intel -mcmodel=medium
INCLUDE_DIR = -I../.. -I../include -I/usr/local/mpi/intel18/mvapich2-2.2/include
LIB_DIR = /rhome/yeong/test/svrtest/libs/lib
INSTALLDLIB = /rhome/yeong/test/svrtest/gsm/bin

#EXP_PARENT = @EXP_PARENT@
#.DEFAULT :
#	ln -fs ${EXP_PARENT}/src/share/$@ $@

${LIB} : ${OBJS}
	${AR} -r $@ ${OBJS}
	${CP} $@ ${INSTALLDLIB}

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
	- ${RM} -f *.f *.o *.x *.a core

dist_clean :
	- ${RM} -f *.f *.o *.x *.a core Makefile
	- ${RM} -rf rii_files
