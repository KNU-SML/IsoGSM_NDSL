#
SHELL = /bin/sh
RM = rm
CD = cd
MAKE = make
INSTALL_DIR = /rhome/yeong/test/svrtest/gsm/bin

GSM_PROGS = share co2 mtn chgr cnvaer cnvalb albaer sfcl pgb sgb sfc0 p2sig sigtool gsml fcst
GSM_MP_PROGS = mpi sfcl_par gsml_par fcst_par
RSM_PROGS = share rsml sfcl rsml gsml rmtn rinpg rpgb co2 rsfc0 p2sig mtn cnvaer cnvalb albaer rfcst
RSM_MP_PROGS = rmpi sfcl_par rsml_par gsml_par rpgb_par rfcst_par
SCM_PROGS = share prescm co2 mtn sfcl pgb cnvaer cnvalb albaer gsml fcst
SCM_MPI_PROGS = mpi
GDAS_PROGS = share co2 mtn chgr cnvaer cnvalb albaer sfcl pgb sgb sfc0 p2sig sigtool gsml fcst cldtune combrtrh3 precipadj acqc cqc oiqc prevents ssi combbufr combobog edbuoy eddrib ednosat edrejx edsatw edswnd edupat ednosat fixpaobs grepsat postvents
GDAS_MPI_PROGS = gsml_par fcst_par mpi

gsm_single : 
	for prog in ${GSM_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done

gsm_thread : 
	for prog in ${GSM_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done

gsm_mpi : 
	for prog in ${GSM_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done
	for prog in ${GSM_MP_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done

gsm_hybrid : 
	for prog in ${GSM_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done
	for prog in ${GSM_MP_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done

rsm_single : 
	for prog in ${RSM_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done

rsm_thread : 
	for prog in ${RSM_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done

rsm_mpi : 
	for prog in ${RSM_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done
	for prog in ${RSM_MP_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done

rsm_hybrid : 
	for prog in ${RSM_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done
	for prog in ${RSM_MP_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done

scm :
	for prog in ${SCM_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} ; \
			${CD} .. ; \
		fi ; \
	done

gdas : 
	for prog in ${GDAS_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done
	for prog in ${GDAS_MP_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done

clean_gsm :  
	for prog in ${GSM_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} clean ; \
			${CD} .. ; \
		fi ; \
	done

gsm_mpi_clean :  
	for prog in ${GSM_MP_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} clean ; \
			${CD} .. ; \
		fi ; \
	done

rsm_mpi_clean :  
	for prog in ${RSM_MP_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} clean ; \
			${CD} .. ; \
		fi ; \
	done

clean_rsm :  
	for prog in ${RSM_PROGS} ${RSM_MP_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} clean ; \
			${CD} .. ; \
		fi ; \
	done

clean_scm :
	for prog in ${SCM_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} clean ; \
			${CD} .. ; \
		fi ; \
	done
clean_gdas :  
	for prog in ${GDAS_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} clean ; \
			${CD} .. ; \
		fi ; \
	done
