#!/bin/sh
#
# @ job_type = serial
# @ output = albaer.out.$(jobid).out
# @ error = albaer.out.$(jobid).err
# @ queue

set -x
FTNID=fort.
CON_DIR=/rhome/yeong/test/svrtest/libs/con
INSTALL_DIR=/rhome/yeong/test/svrtest/gsm/bin
BIN_DIR=/rhome/yeong/test/svrtest/gsm/bin
MODEL_DIR=/rhome/yeong/test/svrtest/gsm
#
jcap=62
levs=28
lonf=192
latg=94
echo $jcap $levs $lonf $latg
#
cd $MODEL_DIR/src/cnvalb
#
# cnvalb
#
PROG=cnvalb
#
# === NAMELIST VARIABLES ===
#  imax,jmax    : dimensions for output data (default is model resolution)
#  igau         : =1 output data on gaussing grid
#                 =0 output data on regular lat/lon grid
#  kssn         : N.H. season =1-4 (win,spr,sum,fal), =5 annual mean.
#                 =0 prodece win,spr,sum,and fal all at once
#  clab         : label for output title (a4)
#  ichk         : =1 print at each check point
#                 =0 no check print
#  iout         : =1 output formatted data, scale 1E3
#                 =0 output unformatted data in fraction
#  outdr        : output data file name and path (A48)
#-T40
#  imax=128, jmax=64, igau=1, kssn=0, clab='T040', ichk=1, iout=1
#-T62
#  imax=192, jmax=94, igau=1, kssn=0, clab='T062', ichk=1, iout=1
#-T126
#  imax=384, jmax=190, igau=1, kssn=0, clab='T126', ichk=1, iout=1
#-1DEG
#  imax=360, jmax=180, igau=0, kssn=0, clab='1DEG', ichk=1, iout=1
#
cat > $PROG.parm << EOF
 &NAMEIN
  CLAB="T$jcap",
 &END
EOF
rm -f ${FTNID}11 2>/dev/null
ln -fs $CON_DIR/mdlv1e.asc ${FTNID}11
$BIN_DIR/$PROG.x <$PROG.parm 1>$PROG.out 2>&1
if [ $? -ne 0 ] ; then
	echo error in $PROG.x
	exit 8 
fi
rm -f ${FTNID}11
#
# cnvaer
#
PROG=cnvaer
# === NAMELIST VARIABLES ===
#  imax,jmax    : dimensions for output data
#  (LONR,LATR) fcs   T40:128x64; T62:192x94; T126:384x190 
#              rad       128x64      128x62       256x126
#  igau         : =1 output data on gaussing grid
#                 =0 output data on regular lat/lon grid
#  issn         : =1-4 (win,spr,sum,fal), =5 annual mean.
#  clab         : label for output title (a4)
#  ichk         : =1 print at each check point
#                 =0 no check print
#  iout         : =1 output formatted data, scale 1E3
#                 =0 output unformatted data in fraction
#-T126
#  imax=384, jmax=190, igau=1, issn=5, clab='T126', ichk=1, iout=1
#-1DEG
#  imax=360, jmax=180, igau=0, issn=5, clab='1DEG', ichk=1, iout=1
cat > $PROG.parm << EOF
 &NAMEIN
   CLAB="T$JCAP",
 &END
EOF
rm -f ${FTNID}11 ${FTNID}51 2>/dev/null
ln -fs $CON_DIR/aerdefb.asc ${FTNID}11
ln -fs aer_out.asc                      ${FTNID}51
$BIN_DIR/$PROG.x <$PROG.parm 1>$PROG.out 2>&1
if [ $? -ne 0 ] ; then
	echo error in $PROG.x
	exit 8 
fi
rm -f ${FTNID}11 ${FTNID}51
#
# albaer
#
PROG=albaer
#
# === NAMELIST VARIABLES ===
#  imax,jmax    : dimensions for output data
#  (LONR,LATR) fcs   T40:128x64; T62:192x94; T126:384x190 
#              rad       128x64      128x62       256x126
#  infl         : =1 input data is formatted integer, scaled 1E3
#                 =0 input data is unformatted in fraction
#  ialb,iaer    : num of input albedo/aerosol data sets
#                 monthly=12, seasonal=4, annual mean=1
#-T62
#  imax=192, jmax=94, infl=1, ialb=4, iaer=1
#-T126
#  imax=384, jmax=190, infl=1, ialb=4, iaer=1
#-1DEG
#  imax=360, jmax=180, infl=1, ialb=4, iaer=1
cat > $PROG.parm << EOF
 &NAMEIN
 &END
EOF
rm -f ${FTNID}11 ${FTNID}12 ${FTNID}13 ${FTNID}14 ${FTNID}21 ${FTNID}51 2>/dev/null
ln -fs alb_out.asc.win ${FTNID}11
ln -fs alb_out.asc.spr ${FTNID}12
ln -fs alb_out.asc.sum ${FTNID}13
ln -fs alb_out.asc.fal ${FTNID}14
ln -fs aer_out.asc     ${FTNID}21
ln -fs $INSTALL_DIR/albaer.snl ${FTNID}51
$BIN_DIR/$PROG.x <$PROG.parm 1>$PROG.out 2>&1
if [ $? -ne 0 ] ; then
	echo error in $PROG.x
	exit 8 
fi
rm alb_out.asc.win alb_out.asc.spr alb_out.asc.sum alb_out.asc.fal aer_out.asc
rm -f ${FTNID}11 ${FTNID}12 ${FTNID}13 ${FTNID}14 ${FTNID}21 ${FTNID}51
echo "albaer.snl created on `date`" >$INSTALL_DIR/.albaer
