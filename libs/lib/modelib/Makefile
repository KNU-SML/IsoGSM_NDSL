LIB = modelib.a
OBJS = bsslz1.o filtcof.o gaulat.o gl2gl.o glats.o gspc.o iminv.o \
incdte.o maxmin.o xmaxmin.o parzen.o pder.o pleg.o poly.o sg2sg.o trispl.o valts.o \
tridi2.o idsdef.o wryte.o gtbits.o rmaxmin.o sphertpt.o sfcfld.o bdrysmth.o fft99m.o

COBJS = streamio.o

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
LOAD_FLAGS = -r8 -O2 -convert big_endian
INCLUDE_DIR = -I.. -I../..
LIB_DIR = /rhome/yeong/test/svrtest/libs/lib
EXRTA_LIBS = 

${LIB} : ${OBJS} ${COBJS}
	${AR} -r $@ $?
	${CP} $@ ${LIB_DIR} 
	${CHMOD} 755 ${LIB}

.SUFFIXES : .o .f .F .c

.F.f :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i

.F.o :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i
	${F77} ${FORT_FLAGS} -c $*.f

.f.o :
	${F77} ${FORT_FLAGS} -c $*.f

.c.o :
	${CC} ${INCLUDE_DIR} -c -o $*.o $*.c

clean :
	- ${RM} -rf *.f *.o *.a core rii_files

dist_clean :
	- ${CD} ${LIB_DIR}; \
	${RM} -f ${LIB}
