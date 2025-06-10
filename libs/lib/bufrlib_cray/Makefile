LIB = bufr.a
OBJS = abt.o adn30.o bfrini.o chekstab.o chrtrn.o closbf.o closmg.o conwin.o \
copybf.o copymg.o copysb.o cpbfdx.o datebf.o drstpl.o dxinit.o dxmini.o \
elemdx.o getwin.o idn30.o ifxy.o inctab.o invcon.o invwin.o ipkm.o ireadmg.o \
irev.o iupb.o iupm.o jstify.o lmsg.o lstrpc.o lstrps.o makestab.o msgini.o \
msgupd.o msgwrt.o mvb.o nemock.o nemtab.o nemtba.o nemtbb.o nemtbd.o nenuck.o \
newwin.o numbck.o numtab.o nvnwin.o nxtwin.o openbf.o openmb.o openmg.o padit.o \
parseq.o parusr.o parutg.o pkb.o pkc.o pktdd.o rcstpl.o rdbfdx.o rdmgsb.o \
rdtree.o rdusdx.o readdx.o readmg.o readsb.o rewnbf.o rjust.o seqsdx.o status.o \
strcln.o string.o strnum.o strsuc.o tabent.o tabsub.o ufbcnt.o ufbcpy.o ufbdmp.o \
ufbevn.o ufbint.o ufbinx.o ufbpks.o ufbqcd.o ufbqcp.o ufbrew.o ufbrw.o ufbtab.o \
upb.o upbb.o upc.o uptdd.o usrtpl.o vals.o wrdlen.o writdx.o writsb.o wrtree.o \
wtstat.o

SHELL = /bin/sh
AR = ar
CP = cp
RM = rm
CD = cd
CHMOD = chmod

CC = cc
CPP = gcc -E
F77 = ifort
FORT_FLAGS = 
LOAD_FLAGS = 
LIB_DIR = /rhome/yeong/test/svrtest/libs/lib
INCLUDE_DIR = -I.. -I../..

${LIB} : ${OBJS}
	${AR} -r $@ $?
	${CP} $@ ${LIB_DIR} ; \
	${CD} ${LIB_DIR} ; \
	${CHMOD} 755 ${LIB}

.SUFFIXES : .o .f .F

.F.f :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i

.F.o :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i
	${F77} ${FORT_FLAGS} -c $*.f

.f.o :
	${F77} ${FORT_FLAGS} -c $*.f

clean :
	- ${RM} -f *.f *.o *.x core

dist_clean :
	- ${RM} -f *.f *.o *.x core Makefile

# DO NOT DELETE
