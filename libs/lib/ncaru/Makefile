
# Generic makefile for maintaining the ncaru lib and it's utilities
#
#  PWM - SCD Consulting 16Feb94
#
#
# Copyright 1994 University Corporation for Atmospheric Research (UCAR).
# 
# Permission to use, copy, and modify this software and its documentation 
# for any non-commercial purpose is hereby granted without fee, provided 
# that the above copyright notice appear in all copies and that both that 
# copyright notice and this permission notice appear in supporting 
# documentation. UCAR makes no representations about the suitability of 
# this software for any purpose.   It is provided "as is" without express
# or implied warranty.

SHELL = /bin/sh
RM = rm
CD = cd
MAKE = make
ECHO = echo

LIB_DIR = /rhome/yeong/test/svrtest/libs/lib
ETC_DIR = /rhome/yeong/test/svrtest/libs/etc

#SUBDIRS = lib utils
SUBDIRS = lib

MACHINE = linux

all:
	if [ linux != cray -a linux != t90 -a linux != t3e ] ; then \
		for i in $(SUBDIRS) ; \
		do \
				${CD} $$i ;\
				${ECHO} "Making all in ./$$i"; \
				${MAKE} || exit 8 ;\
				${ECHO} "done..."; \
				${CD} .. ;\
		done ;\
	fi

clean:
	@for i in $(SUBDIRS) ; \
	do \
	    ${CD} $$i ;\
	    ${ECHO} "Making clean in ./$$i"; \
	    ${MAKE} clean ;\
	    ${ECHO} "done..."; \
	    ${CD} .. ;\
	done

depend:
	@for i in $(SUBDIRS) ; \
	do \
	    ${CD} $$i ;\
	    ${ECHO} "Making depend in ./$$i"; \
	    ${MAKE} depend || exit 8 ;\
	    ${ECHO} "done..."; \
	    ${CD} .. ;\
	done
	
install:
	if [ linux != cray -a linux != t90 -a linux != t3e ] ; then \
		for i in $(SUBDIRS) ; \
		do \
				${CD} $$i ;\
				${ECHO} "Installing in ./$$i"; \
				${MAKE} install LIBDIR=${LIB_DIR} INSTDIR=${ETC_DIR} || exit 8 ;\
				${ECHO} "done..."; \
				${CD} .. ;\
		done ;\
	fi
