#!/bin/csh -x
#
#@$-q p32
#@$-r ch060ow
#@$-eo
#@$-lT 2000:00:00
#@$-me
#@$-mu s21437@gen
#@$
setenv F_PROGINF=YES
setenv F_FTRACE=YES
setenv F_ERRCNT=100
#
 setenv F_FF10 /home/short7/s21437/POIFLW/CH060W/AVERAGE/ch060w_uave1.dat
 setenv F_FF11 /home/short7/s21437/POIFLW/CH060W/AVERAGE/ch060w_uave2.dat
 setenv F_FF12 /home/short7/s21437/POIFLW/CH060W/AVERAGE/ch060w_tave1.dat
 setenv F_FF13 /home/short7/s21437/POIFLW/CH060W/AVERAGE/ch060w_tave2.dat
 setenv F_FF14 $HOME/POIFLW/CH060W/TEXTDATA/ch060w_utxt.dat
 setenv F_FF17 $HOME/POIFLW/CH060W/TEXTDATA/ch060w_ttxt.dat
 setenv F_FF18 $HOME/POIFLW/CH060W/TEXTDATA/ch060w_tbtxt.dat
 setenv F_FF19 $HOME/POIFLW/CH060W/TEXTDATA/ch060w_tvtxt.dat
 setenv F_FF20 /home/short7/s21437/SEED/CH060W/ch060w_uAsd1.dat
 setenv F_FF21 /home/short7/s21437/SEED/CH060W/ch060w_uBsd1.dat
 setenv F_FF22 /home/short7/s21437/SEED/CH060W/ch060w_vAsd1.dat
 setenv F_FF23 /home/short7/s21437/SEED/CH060W/ch060w_vBsd1.dat
 setenv F_FF24 /home/short7/s21437/SEED/CH060W/ch060w_wAsd1.dat
 setenv F_FF25 /home/short7/s21437/SEED/CH060W/ch060w_wBsd1.dat
 setenv F_FF26 /home/short7/s21437/SEED/CH060W/ch060w_uAsd2.dat
 setenv F_FF27 /home/short7/s21437/SEED/CH060W/ch060w_uBsd2.dat
 setenv F_FF28 /home/short7/s21437/SEED/CH060W/ch060w_vAsd2.dat
 setenv F_FF29 /home/short7/s21437/SEED/CH060W/ch060w_vBsd2.dat
 setenv F_FF30 /home/short7/s21437/SEED/CH060W/ch060w_wAsd2.dat
 setenv F_FF31 /home/short7/s21437/SEED/CH060W/ch060w_wBsd2.dat
 setenv F_FF32 /home/short7/s21437/SEED/CH060W/ch060w_t1Asd1.dat
 setenv F_FF33 /home/short7/s21437/SEED/CH060W/ch060w_t1Bsd1.dat
 setenv F_FF34 /home/short7/s21437/SEED/CH060W/ch060w_t2Asd1.dat
 setenv F_FF35 /home/short7/s21437/SEED/CH060W/ch060w_t2Bsd1.dat
 setenv F_FF36 /home/short7/s21437/SEED/CH060W/ch060w_t1Asd2.dat
 setenv F_FF37 /home/short7/s21437/SEED/CH060W/ch060w_t1Bsd2.dat
 setenv F_FF38 /home/short7/s21437/SEED/CH060W/ch060w_t2Asd2.dat
 setenv F_FF39 /home/short7/s21437/SEED/CH060W/ch060w_t2Bsd2.dat
 setenv F_FF41 $HOME/POIFLW/CH060W/TEXTDATA/ch060w_utmsrs.dat
 setenv F_FF42 $HOME/POIFLW/CH060W/TEXTDATA/ch060w_ttmsrs.dat
 setenv F_FF51 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_u_A.dat
 setenv F_FF52 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_v_A.dat
 setenv F_FF53 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_w_A.dat
 setenv F_FF54 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_u_B.dat
 setenv F_FF55 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_v_B.dat
 setenv F_FF56 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_w_B.dat
 setenv F_FF57 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_u_C.dat
 setenv F_FF58 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_v_C.dat
 setenv F_FF59 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_w_C.dat
 setenv F_FF60 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_u_D.dat
 setenv F_FF61 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_v_D.dat
 setenv F_FF62 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_w_D.dat
 setenv F_FF71 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_t1A.dat
 setenv F_FF72 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_t2A.dat
 setenv F_FF73 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_t1B.dat
 setenv F_FF74 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_t2B.dat
 setenv F_FF75 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_t1C.dat
 setenv F_FF76 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_t2C.dat
 setenv F_FF77 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_t1D.dat
 setenv F_FF78 /home/short7/s21437/POIFLW/CH060W/TIME/tmsrs_t2D.dat
#
 cd /fshome/home01/s21437/POIFLW/CH060W/SOLV/
 timex ./ch060ow.exe
#
