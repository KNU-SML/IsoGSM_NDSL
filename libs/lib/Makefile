MACHINE = linux
SHELL = /bin/sh
RM = rm
CD = cd
MKDIR = mkdir
MAKE = make

LIB_DIR = /rhome/yeong/test/svrtest/libs/lib

LIBS = w3lib_intel/w3lib.a ncaru/libncaru.a modelib/modelib.a bufrlib_linux/bufr.a

all : ${LIBS}

w3lib_intel/w3lib.a :
	${CD} w3lib_intel ; ${MAKE} || exit 8

bufrlib_linux/bufr.a :
	${CD} bufrlib_linux ; ${MAKE} || exit 8

ncaru/libncaru.a :
	${CD} ncaru ; ${MAKE} install || exit 8

modelib/modelib.a :
	${CD} modelib ; ${MAKE} || exit 8

clean :
	- ${RM} *.a
	${CD} w3lib_intel ; ${MAKE} clean
	${CD} modelib ; ${MAKE} clean
	${CD} ncaru ; ${MAKE} clean
	if [ -s bufrlib_linux/Makefile ] ; then \
	  ${CD} bufrlib_linux ; ${MAKE} clean ; \
	fi
