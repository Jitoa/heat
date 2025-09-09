C=======================================================================
C/                                              2004/09/30  VER 0.10   /
C/                                              2005/12/29  VER 2.00   /
C/                    TURBULENT CHANNEL FLOW                           /
C/                       FOR SX-7 PARALLEL                             /
C/                         WITH OPEN-MP                                /
C/                                                                     /
C/                         <RE_TAU= 60>                                /
C/                      <TEMPERATURE FIELD>                            /
C/                                                                     /
C=======================================================================
C/                                                                     /
C/   PROGRAMER :  HIROYUKI ABE ( SCIENCE UNVERSITY OF TOKYO )          /
C/                TAKAHIRO TSUKAHARA ( TOKYO UNVERSITY OF SCIENCE )    /
C/                                                                     /
C/   FRACTIONAL-STEP METHOD ON STAGGERED GRID (V-V:NON-UNIFORM MESH)   /
C/      4TH-ORDER CENTRAL SCHEME (X,Z DIR.) FOR CONVECTIVE TERM        /
C/      2ND-ORDER CENTRAL SCHEME (Y DIR.CONSISTENT) FOR CONVECTIVE TERM/
C/      4TH-ORDER CENTRAL SCHEME (X,Z DIR.) FOR VISCOUS TERM           /
C/      2ND-ORDER CENTRAL SCHEME FOR THE OTHER TERMS                   /
C/      2ND-ORDER ADAMS-BASHFORTH METHOD FOR TIME ADVANCEMENT          /
C/                                                                     /
C/   TIME ADVANCEMENT                                                  /
C/      2ND-ORDER CRANK-NICOLSON METHOD FOR VISCOUS TERM (Y-DIRECTION) /
C/      2ND-ORDER ADAMS-BASHFORTH METHOD FOR THE OTHERS                /
C/                                                                     /
C/   OTHER METHODS                                                     /
C/      PSEUDOSPECTRAL METHOD FOR POISSON SOLVER                       /
C/                                                                     /
C/         SEED-DATA OPTIONS                                           /
C/          (VEL) REVISED ON 2003.01.08 BY T.TSUKAHARA                 /
C/             UNIT 20-25 => SEEDIN  : BINARY                          /
C/             UNIT 26-31 => SEEDOUT : BINARY                          /
C/          (TEM) REVISED ON 2003.01.08 BY T.TSUKAHARA                 /
C/             UNIT 32-35 => TSEEDIN  : BINARY                         /
C/             UNIT 36-39 => TSEEDOUT : BINARY                         /
C/                                                                     /
C/         TIME-AVERAGE DATA OPTIONS                                   /
C/          (VEL)                                                      /
C/             UNIT 10 => AVEIN  : BINARY                              /
C/             UNIT 11 => AVEOUT : BINARY                              /
C/          (TEM)                                                      /
C/             UNIT 12 => TAVEIN      : BINARY                         /
C/             UNIT 13 => TAVEOUT     : BINARY                         /
C/                                                                     /
C/         TEXT AVERAGE DATA OPTIONS                                   /
C/          (VEL)                                                      /
C/             UNIT 14 => OUTPUT : TEXT                                /
C/                + SHEAR STRESS (UU,VV,WW,UV,UW,VW)                   /
C/                + TRIPLE CORR. (UUU,UUV,VVV,WWV,UVV,UWW)             /
C/                + BUDGET (K,UU,VV,WW,UV) - CONSISTENT FORM           /
C/                + ENERGY SPECTRA                                     /
C/                + TWO-POINT CORRELATION                              /
C/                + SKEWNESS & FLATNESS FACTOR                         /
C/                + RMS VORTICITY FLUCTUATION                          /
C/                + VEL.GRAD.TENSOR II,III-INVARIANT                   /
C/          (TEM)                                                      /
C/             UNIT 17 => TAVE_OUTPUT : TEXT                           /
C/             UNIT 18 => TBUD_OUTPUT : TEXT                           /
C/             UNIT 19 => TVEL_BACK   : TEXT                           /
C/                + MEAN TEMPERATURE , RMS T.FLUCTUATION               /
C/                + TOTAL HEAT FLUX                                    /
C/                + TEMPERATURE VARIANCE (KT)                          /
C/                + HEAT FLUX (UT,VT,WT)                               /
C/                + BUDGET (KT,UT,VT,WT) - CONSISTENT FORM             /
C/                + ENERGY SPECTRA                                     /
C/                + TWO-POINT CORRELATION                              /
C/                + SKEWNESS & FLATNESS FACTOR                         /
C/                                                                     /
C/         TEXT TIME-SERIES DATA OPTIONS                               /
C/          (VEL) REVISED ON 2005.01.12 BY T.TSUKAHARA                 /
C/             UNIT 41 => AVE : TEXT                                   /
C/          (TEM) REVISED ON 2005.01.12 BY T.TSUKAHARA                 /
C/             UNIT 42 => STCS : TEXT                                  /
C=======================================================================
C                                                                       
      PROGRAM CHANNEL                                                   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
      INCLUDE './INCFILE/INCDATA'                                         
      INCLUDE './INCFILE/INCMESH'                                         
C                                                                       
      CHARACTER*32 CODE                                                 
C                                                                       
      COMMON /ANUM/ NAVE                                                
      COMMON /ATIM/ AVETIM                                              
      COMMON /JAVE/ JA,JB,JC,JD,JE,JF,JH,JI,JJ,JK,JL                    
C                                                                       
      COMMON /TDMA/ AA(LY),BB(LY),CC(LY),RK(0:LX,0:KG)                  
      COMMON /TDMAXZ/ AAXZ(LY),BBXZ(LY),CCXZ(LY)                        
      COMMON /TDMAY/ AAY(LY),BBY(LY),CCY(LY)                            
C                                                                       
C----LOCAL VARIABLES-------------------------------------------         
C                                                                       
      DIMENSION U(-4:IG+4,-2:JG+2,-4:KG+4,3)                            
      DIMENSION UT(-4:IG+4,-2:JG+2,-4:KG+4,3)                           
      DIMENSION P(-3:IG+3,-2:JG+2,-3:KG+3)                              
      DIMENSION PHI(-3:IG+3,-2:JG+2,-3:KG+3)                            
      DIMENSION AU(-1:JG+2),AV(-1:JG+2),AW(-1:JG+2),AP(-1:JG+2)         
      DIMENSION AAU(JG),AAV(JG),AAW(JG),AAP(JG)                         
     &         ,AURMS(JG),AVRMS(JG),AWRMS(JG),APRMS(JG)                 
     &         ,ATSST(JG),AREST(JG)                                     
      DIMENSION AU11(JG),AU12(JG),AU13(JG)                              
     &         ,AU22(JG),AU23(JG),AU33(JG)                              
      DIMENSION AAK(JG)                                                 
     &         ,APK(JG),AP11(JG),AP12(JG)                               
     &         ,APIK(JG),API11(JG),API22(JG),API33(JG),API12(JG)        
     &         ,ADK(JG),AD11(JG),AD22(JG),AD33(JG),AD12(JG)             
     &         ,AEK(JG),AE11(JG),AE22(JG),AE33(JG),AE12(JG)             
     &         ,ATAK(JG),ATA11(JG),ATA22(JG),ATA33(JG),ATA12(JG)        
     &         ,APDK(JG),APD22(JG),APD12(JG)                            
     &         ,APS11(JG),APS22(JG),APS33(JG),APS12(JG)                 
      DIMENSION ASTRM(IG,JG,5),ASPAN(KG,JG,5)                           
      DIMENSION ATPC1(IG,JG,4),ATPC3(KG,JG,4)                           
      DIMENSION OMG(JG,3)                                               
      DIMENSION SKEWUV(JG),FLATUV(JG)                                   
C                                                                       
      DIMENSION AUTRI(JG,8),AUQUAD(JG,10),AQAUV(JG,4)                   
      DIMENSION ANQAUV(JG,4)                                            
C                                                                       
      DIMENSION AINV2(JG),AINV3(JG)                                     
      DIMENSION AAINV2(JG),AAINV3(JG)                                   
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /VELO/ UG(-4:IG+4,-2:JG+2,-4:KG+4,3)                       
      COMMON /TEMP/ UTG(-4:IG+4,-2:JG+2,-4:KG+4,3)                      
      COMMON /OLDF/ FFX(0:IG,0:JG,0:KG),FFY(0:IG,0:JG,0:KG)             
     &             ,FFZ(0:IG,0:JG,0:KG)                                 
      COMMON /OLDCN/ CNX_OLD(0:IG,0:JG,0:KG),CNY_OLD(0:IG,0:JG,0:KG)    
     &              ,CNZ_OLD(0:IG,0:JG,0:KG)                            
      COMMON /SCAL/ PG(-3:IG+3,-2:JG+2,-3:KG+3)                         
      COMMON /PHIP/ PHIG(-3:IG+3,-2:JG+2,-3:KG+3)                       
      COMMON /MEAN/ AUG(-1:JG+2),AVG(-1:JG+2)                           
     &             ,AWG(-1:JG+2),APG(-1:JG+2)                           
      COMMON /FFT3/ WA(IG+1,KG,JG,2)                                    
      COMMON /UTAU/ AUTUB,AUTUT                                         
      COMMON /UTAU2/ AUTUB2,AUTUT2                                      
      COMMON /AVE1/ AAUG(JG),AAVG(JG),AAWG(JG),AAPG(JG)                 
      COMMON /AVE2/ AURMSG(JG),AVRMSG(JG),AWRMSG(JG),APRMSG(JG)         
      COMMON /AVE3/ ATSSTG(JG),ARESTG(JG)                               
      COMMON /AVE4/ AU11G(JG),AU12G(JG),AU13G(JG)                       
     &             ,AU22G(JG),AU23G(JG),AU33G(JG)                       
      COMMON /SFUV/ SKEWUVG(JG),FLATUVG(JG)                             
      COMMON /BUDG/ AAKG(JG)                                            
     &         ,APKG(JG),AP11G(JG),AP12G(JG)                            
     &         ,APIKG(JG),API11G(JG),API22G(JG),API33G(JG),API12G(JG)   
     &         ,ADKG(JG),AD11G(JG),AD22G(JG),AD33G(JG),AD12G(JG)        
     &         ,AEKG(JG),AE11G(JG),AE22G(JG),AE33G(JG),AE12G(JG)        
     &         ,ATAKG(JG),ATA11G(JG),ATA22G(JG),ATA33G(JG),ATA12G(JG)   
     &         ,APDKG(JG),APD22G(JG),APD12G(JG)                         
     &         ,APS11G(JG),APS22G(JG),APS33G(JG),APS12G(JG)             
      COMMON /POWR/ ASTRMG(IG,JG,5),ASPANG(KG,JG,5)                     
      COMMON /CORR/ ATPC1G(IG,JG,4),ATPC3G(KG,JG,4)                     
      COMMON /VORT/ OMGG(JG,3)                                          
C                                                                       
      COMMON /UHOA/ AUTRIG(JG,8),AUQUADG(JG,10),AQAUVG(JG,4)            
      COMMON /NQAUV/ ANQAUVG(JG,4)                                      
C                                                                       
      COMMON /AINV/ AINV2G(JG),AINV3G(JG)                               
      COMMON /AAIN/ AAINV2G(JG),AAINV3G(JG)                             
C                                                                       
C----COMBINATION------------------------------------------------        
C                                                                       
      EQUIVALENCE (UG,U),(UTG,UT),(PG,P),(PHIG,PHI)                     
      EQUIVALENCE (AUG,AU),(AVG,AV),(AWG,AW),(APG,AP)                   
      EQUIVALENCE (AAUG,AAU),(AAVG,AAV),(AAWG,AAW),(AAPG,AAP)           
     &           ,(AURMSG,AURMS),(AVRMSG,AVRMS)                         
     &           ,(AWRMSG,AWRMS),(APRMSG,APRMS)                         
     &           ,(ATSSTG,ATSST),(ARESTG,AREST)                         
      EQUIVALENCE (AU11G,AU11),(AU12G,AU12),(AU13G,AU13)                
     &           ,(AU22G,AU22),(AU23G,AU23),(AU33G,AU33)                
      EQUIVALENCE (AAKG,AAK)                                            
     &           ,(APKG,APK),(AP11G,AP11),(AP12G,AP12)                  
     &           ,(APIKG,APIK),(API11G,API11),(API22G,API22)            
     &           ,(API33G,API33),(API12G,API12)                         
     &           ,(ADKG,ADK),(AD11G,AD11),(AD22G,AD22)                  
     &           ,(AD33G,AD33),(AD12G,AD12)                             
     &           ,(AEKG,AEK),(AE11G,AE11),(AE22G,AE22)                  
     &           ,(AE33G,AE33),(AE12G,AE12)                             
     &           ,(ATAKG,ATAK),(ATA11G,ATA11),(ATA22G,ATA22)            
     &           ,(ATA33G,ATA33),(ATA12G,ATA12)                         
     &           ,(APDKG,APDK),(APD22G,APD22),(APD12G,APD12)            
     &           ,(APS11G,APS11),(APS22G,APS22)                         
     &           ,(APS33G,APS33),(APS12G,APS12)                         
      EQUIVALENCE (ASTRMG,ASTRM),(ASPANG,ASPAN)                         
      EQUIVALENCE (ATPC1G,ATPC1),(ATPC3G,ATPC3)                         
      EQUIVALENCE (OMGG,OMG)                                            
      EQUIVALENCE (SKEWUVG,SKEWUV),(FLATUVG,FLATUV)                     
C                                                                       
      EQUIVALENCE (AUTRIG,AUTRI),(AUQUADG,AUQUAD),(AQAUVG,AQAUV)        
      EQUIVALENCE (ANQAUVG,ANQAUV)                                      
C                                                                       
      EQUIVALENCE (AINV2G,AINV2),(AINV3G,AINV3)                         
      EQUIVALENCE (AAINV2G,AAINV2),(AAINV3G,AAINV3)                     
C                                                                       
      INCLUDE './INCFILE/TMAIN_CHOMP'                                     
C                                                                       
C=================================================================      
C                                                                       
      CALL ETIME(WT1)                                                   
C                                                                       
      NN=0                                                              
C---------------------------------- VEROCITY SEED TIME STEP COUNTER     
      NUM=0                                                             
      TIME=0.0D0                                                        
C---------------------------------------- AVERAGE DATA STEP COUNTER     
      NAVE=0                                                            
      AVETIM=0.0D0                                                      
C-------------------------------------- TEMP SEED TIME STEP COUNTER     
      NCALC=0                                                           
      CALTIM=0.0D0                                                      
C----------------------------- TEMP.FIELD AVERAGE DATA STEP COUNTER     
      NSTCS=0                                                           
      STCSTM=0.0D0                                                      
C                                                                       
C------------------------------------- GRID-NUMBER OF AVERAGE POINT     
C     Re_tau=060 (1024*96*512)                                          
      JA=1   !Y+= 0.16 FIRST POINT FROM THE WALL                        
      JB=4   !Y+= 1.18                                                  
      JC=12  !Y+= 5.07                                                  
      JD=19  !Y+=10.35                                                  
      JE=24  !Y+=15.54                                                  
      JF=28  !Y+=20.67                                                  
      JH=31  !Y+=25.14                                                  
      JI=34  !Y+=30.14 MID-HEIGHT                                       
      JJ=39  !Y+=39.56                                                  
      JK=42  !Y+=45.75                                                  
      JL=48  !Y+=58.88 CHANNEL CENTER                                   
C                                                                       
      IERR=0                                                            
C                                                                       
C------------------------------------ SEEDIN INPUT                      
C      CALL SEEDIN(IERR)                                                 
C      CALL TSEEDIN                                                      
C      CALL AVEIN                                                       
C      CALL TAVEIN                                                      
      CALL ZERO                                                        
      CALL CLR                                                         
      NFIN=0                                                           
      IF (IERR.NE.0) GOTO 1                                             
C------------------------------------ SET PARAMETER                     
      RE=120.0D0                                                       
      XL=6.4D0                                                         
      ZL=3.2D0                                                         
      DX=XL/DBLE(IG)                                                   
      DZ=ZL/DBLE(KG)                                                   
C                                                                       
      DT=0.00020D0                                                     
      NDRAW=0                                                           
         NC=10000                                                      
         NG=1000                                                       
C                                                                       
      PR1=0.71D0                                                        
      PR2=0.71D0                                                        
C                                                                       
      CODE='CH060_WIDE_4TH_1024X96X512'                                 
      CALL ETIME(WT2)                                                   
C------------------------------------                                   
      CALL UBC(DIV)                                                     
      CALL FS0                                                          
C                                                                       
      CALL COEF_DEM(PR1,PR2,PPR1,PPR2,PPE1,PPE2,COE1BD,COE2BD)          
C                                                                       
C-- FIRST STEP                                                          
        NN=NN+1                                                         
        NUM=NUM+1                                                       
        NCALC=NCALC+1                                                   
        TIME=TIME+DT                                                    
        AVETIM=AVETIM+DT                                                
        CALTIM=CALTIM+DT                                                
        STCSTM=STCSTM+DT                                                
C                                                                       
        CALL ERG_1STP(PR1,PR2)                                          
        CALL CNTDMAT                                                    
C                                                                       
        CALL CNTDMA                                                     
        CALL FS2                                                        
        CALL UBC(DIV)                                                   
        IF (MOD(NN,10).EQ.0) THEN                                       
         CALL TMSRSU                                                    
         CALL TMSRST                                                    
        ENDIF                                                           
        WRITE(6,2000)                                                   
        WRITE(6,1000)                                                   
     & NN,NUM,NAVE,NCALC,NSTCS,TIME,AVETIM,CALTIM,STCSTM,DIV            
        IF(NDRAW.EQ.0) THEN                                             
          NSS=2                                                         
        ELSE                                                            
          NSS=1                                                         
        ENDIF                                                           
C                                                                       
C-- NDRAW                                                               
C                                                                       
        IF(NDRAW.LT.2) GOTO 15                                          
        DO 10 NST=2,NDRAW                                               
          NN=NN+1                                                       
          NUM=NUM+1                                                     
          TIME=TIME+DT                                                  
          CALL FS1                                                      
          CALL CNTDMA                                                   
          CALL FS2                                                      
          CALL UBC(DIV)                                                 
        WRITE(6,1000)                                                   
     & NN,NUM,NAVE,NCALC,NSTCS,TIME,AVETIM,CALTIM,STCSTM,DIV            
   10   CONTINUE                                                        
C                                                                       
C-- TIME AVERAGE                                                        
C                                                                       
   15   CONTINUE                                                        
        IF(NC.LT.NSS) GOTO 1                                            
        DO 20 NST=NSS,NG                                                
          NN=NN+1                                                       
          NUM=NUM+1                                                     
          NCALC=NCALC+1                                                 
          TIME=TIME+DT                                                  
          AVETIM=AVETIM+DT                                              
          CALTIM=CALTIM+DT                                              
          STCSTM=STCSTM+DT                                              
C                                                                       
          CALL ERG(PR1,PR2)                                             
          CALL CNTDMAT                                                  
         IF (MOD(NN,500).EQ.0) THEN                                     
         NSTCS=NSTCS+1                                                  
          CALL STCS(PPR1,PPR2,NN)                                       
          CALL BUDGEK1(PPR1,COE1BD)                                    
          CALL BUDGEK2(PPR2,COE2BD)                                    
          CALL TPOWSPE                                                 
          CALL TCORREL                                                 
         ENDIF                                                          
C                                                                       
          CALL FS1                                                      
          CALL CNTDMA                                                   
          CALL FS2                                                      
          CALL UBC(DIV)                                                 
          IF (MOD(NN,10).EQ.0) THEN                                     
           CALL TMSRSU                                                  
           CALL TMSRST                                                  
          ENDIF                                                         
         IF (MOD(NN,500).EQ.0) THEN                                     
         NAVE=NAVE+1                                                    
          CALL AVE(NN)                                                  
          CALL BUDGET                                                  
          CALL POWSPE                                                  
          CALL CORREL                                                  
          CALL VRTCTY                                                  
         ENDIF                                                          
   20   CONTINUE                                                        
        WRITE(6,1000)                                                   
     & NN,NUM,NAVE,NCALC,NSTCS,TIME,AVETIM,CALTIM,STCSTM,DIV            
C                                                                       
        IF(NC.LT.NG) GOTO 1                                             
        NF=NC/NG-1                                                      
        DO 21 NSU=1,NF                                                  
        DO 22 NST=1,NG                                                  
          NN=NN+1                                                       
          NUM=NUM+1                                                     
          NCALC=NCALC+1                                                 
          TIME=TIME+DT                                                  
          AVETIM=AVETIM+DT                                              
          CALTIM=CALTIM+DT                                              
          STCSTM=STCSTM+DT                                              
C                                                                       
          CALL ERG(PR1,PR2)                                             
          CALL CNTDMAT                                                  
         IF (MOD(NN,500).EQ.0) THEN                                     
         NSTCS=NSTCS+1                                                  
          CALL STCS(PPR1,PPR2,NN)                                       
          CALL BUDGEK1(PPR1,COE1BD)                                    
          CALL BUDGEK2(PPR2,COE2BD)                                    
          CALL TPOWSPE                                                 
          CALL TCORREL                                                 
         ENDIF                                                          
C                                                                       
          CALL FS1                                                      
          CALL CNTDMA                                                   
          CALL FS2                                                      
          CALL UBC(DIV)                                                 
          IF (MOD(NN,10).EQ.0) THEN                                     
           CALL TMSRSU                                                  
           CALL TMSRST                                                  
          ENDIF                                                         
         IF (MOD(NN,500).EQ.0) THEN                                     
         NAVE=NAVE+1                                                    
          CALL AVE(NN)                                                  
          CALL BUDGET                                                  
          CALL POWSPE                                                  
          CALL CORREL                                                  
          CALL VRTCTY                                                  
         ENDIF                                                          
   22   CONTINUE                                                        
        WRITE(6,1000)                                                   
     & NN,NUM,NAVE,NCALC,NSTCS,TIME,AVETIM,CALTIM,STCSTM,DIV            
   21   CONTINUE                                                        
C------------------------------------                                   
      CALL SEEDOUT(CODE)                                                
      CALL TSEEDOUT                                                     
C------------------------------------ PARALLEL FINISH                   
      CALL ETIME(WT3)                                                   
C                                                                       
      NFIN=NFIN+NC                                                      
      CALL AVEOUT                                                       
      CALL TAVEOUT                                                      
      CALL OUTPUT(CODE)                                                 
      CALL TAVE_OUTPUT(PPE1,PPE2)                                       
      CALL TBUD_OUTPUT                                                  
      CALL TVEL_BACK                                                    
C                                                                       
      CALL ETIME(WT4)                                                   
        WRITE(6,*) ' TIME', ' ', '=',(WT4-WT1), '(sec)'                 
C                                                                       
 1000 FORMAT(I6,I8,I7,I8,I7,4F10.4,E10.3)                               
 2000 FORMAT('    NN     NUM   NAVE   NCALC  NSTCS      TIME    AVETIM  
     &  STCSTM    CALTIM       DIV')                                    
C                                                                       
    1 CONTINUE                                                          
      STOP                                                              
      END                                                               
C                                                                       
C====================================================================   
C     ZERO (VELOCITY FIELD)                                             
C     PROGRAMED BY H.ABE                                                
      SUBROUTINE ZERO                                                   
C     1998.01.01                                                        
C====================================================================   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
C                                                                       
      COMMON /ANUM/ NAVE                                                
      COMMON /ATIM/ AVETIM                                              
C                                                                       
C----LOCAL VARIABLES-------------------------------------------         
C                                                                       
      DIMENSION AAU(JG),AAV(JG),AAW(JG),AAP(JG)                         
     &         ,AURMS(JG),AVRMS(JG),AWRMS(JG),APRMS(JG)                 
     &         ,ATSST(JG),AREST(JG)                                     
      DIMENSION AAK(JG)                                                 
     &         ,APK(JG),AP11(JG),AP12(JG)                               
     &         ,APIK(JG),API11(JG),API22(JG),API33(JG),API12(JG)        
     &         ,ADK(JG),AD11(JG),AD22(JG),AD33(JG),AD12(JG)             
     &         ,AEK(JG),AE11(JG),AE22(JG),AE33(JG),AE12(JG)             
     &         ,ATAK(JG),ATA11(JG),ATA22(JG),ATA33(JG),ATA12(JG)        
     &         ,APDK(JG),APD22(JG),APD12(JG)                            
     &         ,APS11(JG),APS22(JG),APS33(JG),APS12(JG)                 
      DIMENSION AU11(JG),AU12(JG),AU13(JG)                              
     &         ,AU22(JG),AU23(JG),AU33(JG)                              
      DIMENSION ASTRM(IG,JG,5),ASPAN(KG,JG,5)                           
      DIMENSION ATPC1(IG,JG,4),ATPC3(KG,JG,4)                           
      DIMENSION OMG(JG,3)                                               
      DIMENSION SKEWUV(JG),FLATUV(JG)                                   
C                                                                       
      DIMENSION AUTRI(JG,8),AUQUAD(JG,10),AQAUV(JG,4)                   
      DIMENSION ANQAUV(JG,4)                                            
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /UTAU/ AUTUB,AUTUT                                         
      COMMON /UTAU2/ AUTUB2,AUTUT2                                      
      COMMON /AVE1/ AAUG(JG),AAVG(JG),AAWG(JG),AAPG(JG)                 
      COMMON /AVE2/ AURMSG(JG),AVRMSG(JG),AWRMSG(JG),APRMSG(JG)         
      COMMON /AVE3/ ATSSTG(JG),ARESTG(JG)                               
      COMMON /AVE4/ AU11G(JG),AU12G(JG),AU13G(JG)                       
     &             ,AU22G(JG),AU23G(JG),AU33G(JG)                       
      COMMON /BUDG/ AAKG(JG)                                            
     &         ,APKG(JG),AP11G(JG),AP12G(JG)                            
     &         ,APIKG(JG),API11G(JG),API22G(JG),API33G(JG),API12G(JG)   
     &         ,ADKG(JG),AD11G(JG),AD22G(JG),AD33G(JG),AD12G(JG)        
     &         ,AEKG(JG),AE11G(JG),AE22G(JG),AE33G(JG),AE12G(JG)        
     &         ,ATAKG(JG),ATA11G(JG),ATA22G(JG),ATA33G(JG),ATA12G(JG)   
     &         ,APDKG(JG),APD22G(JG),APD12G(JG)                         
     &         ,APS11G(JG),APS22G(JG),APS33G(JG),APS12G(JG)             
      COMMON /POWR/ ASTRMG(IG,JG,5),ASPANG(KG,JG,5)                     
      COMMON /CORR/ ATPC1G(IG,JG,4),ATPC3G(KG,JG,4)                     
      COMMON /VORT/ OMGG(JG,3)                                          
      COMMON /SFUV/ SKEWUVG(JG),FLATUVG(JG)                             
C                                                                       
      COMMON /UHOA/ AUTRIG(JG,8),AUQUADG(JG,10),AQAUVG(JG,4)            
      COMMON /NQAUV/ ANQAUVG(JG,4)                                      
C----COMBINATION------------------------------------------------        
C                                                                       
      EQUIVALENCE (AAUG,AAU),(AAVG,AAV),(AAWG,AAW),(AAPG,AAP)           
     &           ,(AURMSG,AURMS),(AVRMSG,AVRMS)                         
     &           ,(AWRMSG,AWRMS),(APRMSG,APRMS)                         
     &           ,(ATSSTG,ATSST),(ARESTG,AREST)                         
      EQUIVALENCE (AU11G,AU11),(AU12G,AU12),(AU13G,AU13)                
     &           ,(AU22G,AU22),(AU23G,AU23),(AU33G,AU33)                
      EQUIVALENCE (AAKG,AAK)                                            
     &           ,(APKG,APK),(AP11G,AP11),(AP12G,AP12)                  
     &           ,(APIKG,APIK),(API11G,API11),(API22G,API22)            
     &           ,(API33G,API33),(API12G,API12)                         
     &           ,(ADKG,ADK),(AD11G,AD11),(AD22G,AD22)                  
     &           ,(AD33G,AD33),(AD12G,AD12)                             
     &           ,(AEKG,AEK),(AE11G,AE11),(AE22G,AE22)                  
     &           ,(AE33G,AE33),(AE12G,AE12)                             
     &           ,(ATAKG,ATAK),(ATA11G,ATA11),(ATA22G,ATA22)            
     &           ,(ATA33G,ATA33),(ATA12G,ATA12)                         
     &           ,(APDKG,APDK),(APD22G,APD22),(APD12G,APD12)            
     &           ,(APS11G,APS11),(APS22G,APS22)                         
     &           ,(APS33G,APS33),(APS12G,APS12)                         
      EQUIVALENCE (ASTRMG,ASTRM),(ASPANG,ASPAN)                         
      EQUIVALENCE (ATPC1G,ATPC1),(ATPC3G,ATPC3)                         
      EQUIVALENCE (OMGG,OMG)                                            
      EQUIVALENCE (SKEWUVG,SKEWUV),(FLATUVG,FLATUV)                     
C                                                                       
      EQUIVALENCE (AUTRIG,AUTRI),(AUQUADG,AUQUAD),(AQAUVG,AQAUV)        
      EQUIVALENCE (ANQAUVG,ANQAUV)                                      
C                                                                       
C=================================================================      
C---------------------------------------- AVERAGE DATA STEP COUNTER     
      NAVE=0                                                            
      AVETIM=0.0D0                                                      
C                                                                       
      AUTUT=0.0D0                                                       
      AUTUB=0.0D0                                                       
      AUTUT2=0.0D0                                                      
      AUTUB2=0.0D0                                                      
C                                                                       
      DO 10 J=1,JG                                                      
           AAU(J)=0.0D0                                                 
           AAV(J)=0.0D0                                                 
           AAW(J)=0.0D0                                                 
           AAP(J)=0.0D0                                                 
         AURMS(J)=0.0D0                                                 
         AVRMS(J)=0.0D0                                                 
         AWRMS(J)=0.0D0                                                 
         APRMS(J)=0.0D0                                                 
         ATSST(J)=0.0D0                                                 
         AREST(J)=0.0D0                                                 
         AU11(J)=0.0D0                                                  
         AU12(J)=0.0D0                                                  
         AU13(J)=0.0D0                                                  
         AU22(J)=0.0D0                                                  
         AU23(J)=0.0D0                                                  
         AU33(J)=0.0D0                                                  
C                                                                       
           AAK(J)=0.0D0                                                 
           APK(J)=0.0D0                                                 
          AP11(J)=0.0D0                                                 
          AP12(J)=0.0D0                                                 
           AEK(J)=0.0D0                                                 
          AE11(J)=0.0D0                                                 
          AE22(J)=0.0D0                                                 
          AE33(J)=0.0D0                                                 
          AE12(J)=0.0D0                                                 
           ADK(J)=0.0D0                                                 
          AD11(J)=0.0D0                                                 
          AD22(J)=0.0D0                                                 
          AD33(J)=0.0D0                                                 
          AD12(J)=0.0D0                                                 
          APIK(J)=0.0D0                                                 
         API11(J)=0.0D0                                                 
         API22(J)=0.0D0                                                 
         API33(J)=0.0D0                                                 
         API12(J)=0.0D0                                                 
          ATAK(J)=0.0D0                                                 
         ATA11(J)=0.0D0                                                 
         ATA22(J)=0.0D0                                                 
         ATA33(J)=0.0D0                                                 
          APDK(J)=0.0D0                                                 
         APD22(J)=0.0D0                                                 
         APD12(J)=0.0D0                                                 
         APS11(J)=0.0D0                                                 
         APS22(J)=0.0D0                                                 
         APS33(J)=0.0D0                                                 
         APS12(J)=0.0D0                                                 
         ATA12(J)=0.0D0                                                 
C                                                                       
        OMG(J,1)=0.0D0                                                  
        OMG(J,2)=0.0D0                                                  
        OMG(J,3)=0.0D0                                                  
C                                                                       
        SKEWUV(J)=0.0D0                                                 
        FLATUV(J)=0.0D0                                                 
C                                                                       
        AUTRI(J,1)=0.0D0                                                
        AUTRI(J,2)=0.0D0                                                
        AUTRI(J,3)=0.0D0                                                
        AUTRI(J,4)=0.0D0                                                
        AUTRI(J,5)=0.0D0                                                
        AUTRI(J,6)=0.0D0                                                
        AUTRI(J,7)=0.0D0                                                
        AUTRI(J,8)=0.0D0                                                
C                                                                       
        AUQUAD(J,1)=0.0D0                                               
        AUQUAD(J,2)=0.0D0                                               
        AUQUAD(J,3)=0.0D0                                               
        AUQUAD(J,4)=0.0D0                                               
        AUQUAD(J,5)=0.0D0                                               
        AUQUAD(J,6)=0.0D0                                               
        AUQUAD(J,7)=0.0D0                                               
        AUQUAD(J,8)=0.0D0                                               
        AUQUAD(J,9)=0.0D0                                               
        AUQUAD(J,10)=0.0D0                                              
C                                                                       
        AQAUV(J,1)=0.0D0                                                
        AQAUV(J,2)=0.0D0                                                
        AQAUV(J,3)=0.0D0                                                
        AQAUV(J,4)=0.0D0                                                
C                                                                       
   10 CONTINUE                                                          
C                                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
      DO 30 J=1,JG                                                      
      DO 30 I=1,IG                                                      
        ASTRM(I,J,1)=0.0D0                                              
        ASTRM(I,J,2)=0.0D0                                              
        ASTRM(I,J,3)=0.0D0                                              
        ASTRM(I,J,4)=0.0D0                                              
        ASTRM(I,J,5)=0.0D0                                              
        ATPC1(I,J,1)=0.0D0                                              
        ATPC1(I,J,2)=0.0D0                                              
        ATPC1(I,J,3)=0.0D0                                              
        ATPC1(I,J,4)=0.0D0                                              
   30 CONTINUE                                                          
C                                                                       
!$OMP PARALLEL DO PRIVATE(K)                                            
      DO 32 J=1,JG                                                      
      DO 32 K=1,KG                                                      
        ASPAN(K,J,1)=0.0D0                                              
        ASPAN(K,J,2)=0.0D0                                              
        ASPAN(K,J,3)=0.0D0                                              
        ASPAN(K,J,4)=0.0D0                                              
        ASPAN(K,J,5)=0.0D0                                              
        ATPC3(K,J,1)=0.0D0                                              
        ATPC3(K,J,2)=0.0D0                                              
        ATPC3(K,J,3)=0.0D0                                              
        ATPC3(K,J,4)=0.0D0                                              
   32 CONTINUE                                                          
C                                                                       
      DO 70 J=1,JG                                                      
         ANQAUV(J,1)=0.0D0                                              
         ANQAUV(J,2)=0.0D0                                              
         ANQAUV(J,3)=0.0D0                                              
         ANQAUV(J,4)=0.0D0                                              
   70 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     CALC. TEMPORARY VELOCIYT AT FIRST STEP (VELOCITY FIELD)           
C     PROGRAMED BY H.ABE                                                
      SUBROUTINE FS0                                                    
C     1998.01.01                                                        
C     2000.08.05 REVISED BY H.ABE                                       
C====================================================================   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
C                                                                       
C----LOCAL VARIABLES-------------------------------------------         
C                                                                       
      INCLUDE './INCFILE/INCDATA'                                         
      INCLUDE './INCFILE/INCMESH'                                         
      COMMON /FD2C/ CEN2_F(0:JG,2),CEN2_M(0:JG,2),CEN2_B(0:JG,2)        
      COMMON /TDMA/ AA(LY),BB(LY),CC(LY),RK(0:LX,0:KG)                  
      COMMON /TDMAXZ/ AAXZ(LY),BBXZ(LY),CCXZ(LY)                        
      COMMON /TDMAY/ AAY(LY),BBY(LY),CCY(LY)                            
      DIMENSION U(-4:IG+4,-2:JG+2,-4:KG+4,3)                            
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /VELO/ UG(-4:IG+4,-2:JG+2,-4:KG+4,3)                       
      COMMON /OLDF/ FFX(0:IG,0:JG,0:KG),FFY(0:IG,0:JG,0:KG)             
     &             ,FFZ(0:IG,0:JG,0:KG)                                 
      COMMON /OLDCN/ CNX_OLD(0:IG,0:JG,0:KG),CNY_OLD(0:IG,0:JG,0:KG)    
     &              ,CNZ_OLD(0:IG,0:JG,0:KG)                            
C                                                                       
C----COMBINATION------------------------------------------------        
C                                                                       
      EQUIVALENCE (UG,U)                                                
C                                                                       
C=================================================================      
C                                                                       
C-- PRE-CALC.                                                           
      PDX=1.0D0/DX                                                      
      PDZ=1.0D0/DZ                                                      
      PDX2=1.0D0/(DX**2)                                                
      PDZ2=1.0D0/(DZ**2)                                                
      PRE=1.0D0/RE                                                      
C                                                                       
      DO 10 J=0,JG+1                                                    
        DR2(J)=DR(J)+DR(J+1)                                            
        DRR2(J)=DRR(J)+DRR(J+1)                                         
        PDR(J)=1.0D0/DR(J)                                              
        PDRR(J)=1.0D0/DRR(J)                                            
        PDR2(J)=1.0D0/(DR(J)+DR(J+1))                                   
        PDR3(J)=1.0D0/(DR(J)*DR(J+1)*(DR(J)+DR(J+1)))                   
        PDRR3(J)=1.0D0/(DRR(J)*DRR(J+1)*(DRR(J)+DRR(J+1)))              
   10 CONTINUE                                                          
C                                                                       
      DO 20 J=0,JG                                                      
       CEN2_F(J,1)=2.0D0/( DR(J+1)*( DR(J)+DR(J+1) ) )                  
       CEN2_F(J,2)=2.0D0/( DRR(J+1)*( DRR(J)+DRR(J+1) ) )               
       CEN2_M(J,1)=-2.0D0/( DR(J)*DR(J+1) )                             
       CEN2_M(J,2)=-2.0D0/( DRR(J)*DRR(J+1) )                           
       CEN2_B(J,1)=2.0D0/( DR(J)*( DR(J)+DR(J+1) ) )                    
       CEN2_B(J,2)=2.0D0/( DRR(J)*( DRR(J)+DRR(J+1) ) )                 
   20 CONTINUE                                                          
C                                                                       
      PAI=4.0D0*DATAN(1.0D0)                                            
C                                                                       
      RL1=1.0D0/288.D0*(PDX**2)                                         
      RL2=2.0D0*PAI/DBLE(IG)                                            
      RL3=4.0D0*PAI/DBLE(IG)                                            
      RL4=6.0D0*PAI/DBLE(IG)                                            
C                                                                       
      RN1=1.0D0/288.D0*(PDZ**2)                                         
      RN2=2.0D0*PAI/DBLE(KG)                                            
      RN3=4.0D0*PAI/DBLE(KG)                                            
      RN4=6.0D0*PAI/DBLE(KG)                                            
C                                                                       
C-------------------------------------------------------------          
C                                                                       
!$OMP PARALLEL DO PRIVATE(L)                                            
      DO 30 N=0,KG-1                                                    
      DO 30 L=0,IG-1                                                    
        RK(L,N)=RL1*(730.0D0-783.0D0*DCOS(DBLE(L)*RL2)                  
     &          +54.0D0*DCOS(DBLE(L)*RL3)-1.0D0*DCOS(DBLE(L)*RL4))      
     &          +RN1*(730.0D0-783.0D0*DCOS(DBLE(N)*RN2)                 
     &          +54.0D0*DCOS(DBLE(N)*RN3)-1.0D0*DCOS(DBLE(N)*RN4))      
   30 CONTINUE                                                          
      RK(0,0)=0.0D0                                                     
C                                                                       
      DO 40 J=1,JG                                                      
        AA(J)=-1.0D0/DRR(J  )/DR(J)                                     
        BB(J)=(DRR(J)+DRR(J+1))/DRR(J)/DRR(J+1)/DR(J)                   
        CC(J)=-1.0D0/DRR(J+1)/DR(J)                                     
C        AA(J)=-2.0D0/DRR(J  )/(DRR(J)+DRR(J+1))                        
C        BB(J)=+2.0D0/DRR(J  )/DRR(J+1)                                 
C        CC(J)=-2.0D0/DRR(J+1)/(DRR(J)+DRR(J+1))                        
        CCXZ(J)=2.0D0/DRR(J  )/(DRR(J)+DRR(J+1))                        
        BBXZ(J)=2.0D0/DRR(J  )/DRR(J+1)+2.0D0*RE/DT                     
        AAXZ(J)=2.0D0/DRR(J+1)/(DRR(J)+DRR(J+1))                        
        CCY(J)=2.0D0/DR(J  )/(DR(J)+DR(J+1))                            
        BBY(J)=2.0D0/DR(J  )/DR(J+1)+2.0D0*RE/DT                        
        AAY(J)=2.0D0/DR(J+1)/(DR(J)+DR(J+1))                            
   40 CONTINUE                                                          
C                                                                       
C                                                                       
C-------------------------------------------------------------          
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 50 J=1,JG                                                      
      DO 50 K=1,KG                                                      
      DO 50 I=1,IG                                                      
       FFX(I,J,K)=0.0D0                                                 
       FFY(I,J,K)=0.0D0                                                 
       FFZ(I,J,K)=0.0D0                                                 
  50  CONTINUE                                                          
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,UC1,UC2,UC3,CONX,WC1,WC2,WC3,CONZ,UVISX   
!$OMP&                   ,WVISX,VVISX,UVISZ,WVISZ,VVISZ)                
      DO 60 J=1,JG                                                      
      DO 60 K=1,KG                                                      
      DO 60 I=1,IG                                                      
C--- STREAMWISE DIRECTION                                               
         UC1=                                                           
     &      (9.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I,J,K,1)+U(I-1,J,K,1))*0.5D0)-           
     &      (1.0D0/8.0D0)*((U(I+1,J,K,1)+U(I-2,J,K,1))*0.5D0))*         
     &      (U(I,J,K,1)-U(I-1,J,K,1))*PDX                               
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I,J,K,1)+U(I+1,J,K,1))*0.5D0)-           
     &      (1.0D0/8.0D0)*((U(I+2,J,K,1)+U(I-1,J,K,1))*0.5D0))*         
     &      (U(I+1,J,K,1)-U(I,J,K,1))*PDX                               
     &      )*0.5D0                                                     
     &      -                                                           
     &      (1.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I-1,J,K,1)+U(I-2,J,K,1))*0.5D0)-         
     &      (1.0D0/8.0D0)*((U(I,J,K,1)+U(I-3,J,K,1))*0.5D0))*           
     &      (U(I,J,K,1)-U(I-3,J,K,1))/3.0D0*PDX                         
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I+1,J,K,1)+U(I+2,J,K,1))*0.5D0)-         
     &      (1.0D0/8.0D0)*((U(I,J,K,1)+U(I+3,J,K,1))*0.5D0))*           
     &      (U(I+3,J,K,1)-U(I,J,K,1))/3.0D0*PDX                         
     &      )*0.5D0                                                     
        UC2=0.25                                                        
     &      *((U(I+1,J  ,K,2)+U(I,J  ,K,2))                             
     &       *(U(I  ,J+1,K,1)-U(I,J  ,K,1))*PDRR(J+1)                   
     &       +(U(I+1,J-1,K,2)+U(I,J-1,K,2))                             
     &       *(U(I  ,J  ,K,1)-U(I,J-1,K,1))*PDRR(J))                    
         UC3=                                                           
     &      (9.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I+1,J,K,3)+U(I,J,K,3))*0.5D0)-           
     &      (1.0D0/8.0D0)*((U(I+2,J,K,3)+U(I-1,J,K,3))*0.5D0))*         
     &      (U(I,J,K+1,1)-U(I,J,K,1))*PDZ                               
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I+1,J,K-1,3)+U(I,J,K-1,3))*0.5D0)-       
     &      (1.0D0/8.0D0)*((U(I+2,J,K-1,3)+U(I-1,J,K-1,3))*0.5D0))*     
     &      (U(I,J,K,1)-U(I,J,K-1,1))*PDZ                               
     &      )*0.5D0                                                     
     &      -                                                           
     &      (1.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I+1,J,K+1,3)+U(I,J,K+1,3))*0.5D0)-       
     &      (1.0D0/8.0D0)*((U(I+2,J,K+1,3)+U(I-1,J,K+1,3))*0.5D0))*     
     &      (U(I,J,K+3,1)-U(I,J,K,1))/3.0D0*PDZ                         
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I+1,J,K-2,3)+U(I,J,K-2,3))*0.5D0)-       
     &      (1.0D0/8.0D0)*((U(I+2,J,K-2,3)+U(I-1,J,K-2,3))*0.5D0))*     
     &      (U(I,J,K,1)-U(I,J,K-3,1))/3.0D0*PDZ                         
     &      )*0.5D0                                                     
        CONX=UC1+UC2+UC3                                                
C                                                                       
        UVISX=(-U(I+2,J,K,1)+16.0D0*U(I+1,J,K,1)-30.0D0*U(I,J,K,1)      
     &          +16.0D0*U(I-1,J,K,1)-U(I-2,J,K,1))/12.0D0               
     &          *PDX2*PRE                                               
C                                                                       
        WVISX=(-U(I,J,K+2,1)+16.0D0*U(I,J,K+1,1)-30.0D0*U(I,J,K,1)      
     &          +16.0D0*U(I,J,K-1,1)-U(I,J,K-2,1))/12.0D0               
     &          *PDZ2*PRE                                               
C                                                                       
        VVISX=(DRR(J+1)*U(I,J-1,K,1)-DRR2(J)*U(I,J,K,1)                 
     &          +DRR(J)*U(I,J+1,K,1))*2.0*PDRR3(J)                      
C                                                                       
        FFX(I,J,K)=2.0D0*RE*(CONX-UVISX-WVISX)-VVISX                    
     &             -2.0D0*RE/DT*U(I,J,K,1)-2.0D0*RE*2.0D0               
C                                                                       
        CNX_OLD(I,J,K)=CONX-UVISX-WVISX                                 
C                                                                       
C--- SPANWISE DIRECTION                                                 
C                                                                       
         WC1=                                                           
     &      (9.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I,J,K+1,1)+U(I,J,K,1))*0.5D0)-           
     &      (1.0D0/8.0D0)*((U(I,J,K+2,1)+U(I,J,K-1,1))*0.5D0))*         
     &      (U(I+1,J,K,3)-U(I,J,K,3))*PDX                               
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I-1,J,K+1,1)+U(I-1,J,K,1))*0.5D0)-       
     &      (1.0D0/8.0D0)*((U(I-1,J,K+2,1)+U(I-1,J,K-1,1))*0.5D0))*     
     &      (U(I,J,K,3)-U(I-1,J,K,3))*PDX                               
     &      )*0.5D0                                                     
     &      -                                                           
     &      (1.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I+1,J,K+1,1)+U(I+1,J,K,1))*0.5D0)-       
     &      (1.0D0/8.0D0)*((U(I+1,J,K+2,1)+U(I+1,J,K-1,1))*0.5D0))*     
     &      (U(I+3,J,K,3)-U(I,J,K,3))/3.0D0*PDX                         
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I-2,J,K+1,1)+U(I-2,J,K,1))*0.5D0)-       
     &      (1.0D0/8.0D0)*((U(I-2,J,K+2,1)+U(I-2,J,K-1,1))*0.5D0))*     
     &      (U(I,J,K,3)-U(I-3,J,K,3))/3.0D0*PDX                         
     &      )*0.5D0                                                     
        WC2=0.25                                                        
     &      *((U(I,J  ,K+1,2)+U(I,J  ,K,2))                             
     &       *(U(I,J+1,K  ,3)-U(I,J  ,K,3))*PDRR(J+1)                   
     &       +(U(I,J-1,K+1,2)+U(I,J-1,K,2))                             
     &       *(U(I,J  ,K  ,3)-U(I,J-1,K,3))*PDRR(J))                    
         WC3=                                                           
     &      (9.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I,J,K+1,3)+U(I,J,K,3))*0.5D0)-           
     &      (1.0D0/8.0D0)*((U(I,J,K+2,3)+U(I,J,K-1,3))*0.5D0))*         
     &      (U(I,J,K+1,3)-U(I,J,K,3))*PDZ                               
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I,J,K,3)+U(I,J,K-1,3))*0.5D0)-           
     &      (1.0D0/8.0D0)*((U(I,J,K+1,3)+U(I,J,K-2,3))*0.5D0))*         
     &      (U(I,J,K,3)-U(I,J,K-1,3))*PDZ                               
     &      )*0.5D0                                                     
     &      -                                                           
     &      (1.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I,J,K+2,3)+U(I,J,K+1,3))*0.5D0)-         
     &      (1.0D0/8.0D0)*((U(I,J,K+3,3)+U(I,J,K,3))*0.5D0))*           
     &      (U(I,J,K+3,3)-U(I,J,K,3))/3.0D0*PDZ                         
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I,J,K-1,3)+U(I,J,K-2,3))*0.5D0)-         
     &      (1.0D0/8.0D0)*((U(I,J,K,3)+U(I,J,K-3,3))*0.5D0))*           
     &      (U(I,J,K,3)-U(I,J,K-3,3))/3.0D0*PDZ                         
     &      )*0.5D0                                                     
        CONZ=WC1+WC2+WC3                                                
C                                                                       
        UVISZ=(-U(I+2,J,K,3)+16.0D0*U(I+1,J,K,3)-30.0D0*U(I,J,K,3)      
     &          +16.0D0*U(I-1,J,K,3)-U(I-2,J,K,3))/12.0D0               
     &          *PDX2*PRE                                               
C                                                                       
        WVISZ=(-U(I,J,K+2,3)+16.0D0*U(I,J,K+1,3)-30.0D0*U(I,J,K,3)      
     &          +16.0D0*U(I,J,K-1,3)-U(I,J,K-2,3))/12.0D0               
     &          *PDZ2*PRE                                               
C                                                                       
        VVISZ=(DRR(J+1)*U(I,J-1,K,3)-DRR2(J)*U(I,J,K,3)                 
     &            +DRR(J)*U(I,J+1,K,3))*2.0*PDRR3(J)                    
C                                                                       
        FFZ(I,J,K)=2.0D0*RE*(CONZ-UVISZ-WVISZ)-VVISZ                    
     &             -2.0D0*RE/DT*U(I,J,K,3)                              
C                                                                       
        CNZ_OLD(I,J,K)=CONZ-UVISZ-WVISZ                                 
C                                                                       
   60 CONTINUE                                                          
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,VC1,VC2,VC3,CONY,UVISY,WVISY,VVISY)       
      DO 70 J=1,JG-1                                                    
      DO 70 K=1,KG                                                      
      DO 70 I=1,IG                                                      
C--- WALL TO WALL DIRECTION                                             
         VC1=                                                           
     &      (9.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*(DR(J)*U(I,J+1,K,1)+DR(J+1)*U(I,J,K,1))      
     &                     *PDR2(J)-                                    
     &       (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I,J+2,K,1)          
     &                      +(2.0D0*DR(J+1)+DR(J+2))*U(I,J-1,K,1))      
     &                     /(DR(J+2)+2.0D0*(DR(J+1)+DR(J))+DR(J-1))     
     &       )                                                          
     &       *(U(I+1,J,K,2)-U(I,J,K,2))*PDX                             
     &      +                                                           
     &      ((9.0D0/8.0D0)*(DR(J)*U(I-1,J+1,K,1)+DR(J+1)*U(I-1,J,K,1))  
     &                     *PDR2(J)-                                    
     &       (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I-1,J+2,K,1)        
     &                      +(2.0D0*DR(J+1)+DR(J+2))*U(I-1,J-1,K,1))    
     &                     /(DR(J+2)+2.0D0*(DR(J+1)+DR(J))+DR(J-1))     
     &       )                                                          
     &       *(U(I,J,K,2)-U(I-1,J,K,2))*PDX                             
     &      )*0.5D0                                                     
     &      -                                                           
     &      (1.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*(DR(J)*U(I+1,J+1,K,1)+DR(J-1)*U(I+1,J,K,1))  
     &                     *PDR2(J)-                                    
     &      (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I+1,J+2,K,1)         
     &                     +(2.0D0*DR(J+1)+DR(J+2))*U(I+1,J-1,K,1))     
     &       )                                                          
     &      *(U(I+3,J,K,2)-U(I,J,K,2))/3.0D0*PDX                        
     &      +                                                           
     &      ((9.0D0/8.0D0)*(DR(J)*U(I-2,J+1,K,1)+DR(J-1)*U(I-2,J,K,1))  
     &                     *PDR2(J)-                                    
     &      (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I-2,J+2,K,1)         
     &                     +(2.0D0*DR(J+1)+DR(J+2))*U(I-2,J-1,K,1))     
     &       )                                                          
     &      *(U(I,J,K,2)-U(I-3,J,K,2))/3.0D0*PDX                        
     &      )*0.5D0                                                     
        VC2=0.5*PDR2(J)                                                 
     &      *((U(I,J+1,K,2)+U(I,J,K,2))                                 
     &       *(U(I,J+1,K,2)-U(I,J,K,2))*PDR(J+1)*DR(J)                  
     &       +(U(I,J,K,2)+U(I,J-1,K,2))                                 
     &       *(U(I,J,K,2)-U(I,J-1,K,2))*PDR(J)*DR(J+1))                 
         VC3=                                                           
     &      (9.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*(DR(J)*U(I,J+1,K,3)+DR(J+1)*U(I,J,K,3))      
     &                     *PDR2(J)-                                    
     &       (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I,J+2,K,3)          
     &                      +(2.0D0*DR(J+1)+DR(J+2))*U(I,J-1,K,3))      
     &                     /(DR(J+2)+2.0D0*(DR(J+1)+DR(J))+DR(J-1))     
     &       )                                                          
     &       *(U(I,J,K+1,2)-U(I,J,K,2))*PDZ                             
     &      +                                                           
     &      ((9.0D0/8.0D0)*(DR(J)*U(I,J+1,K-1,3)+DR(J+1)*U(I,J,K-1,3))  
     &                     *PDR2(J)-                                    
     &       (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I,J+2,K-1,3)        
     &                      +(2.0D0*DR(J+1)+DR(J+2))*U(I,J-1,K-1,3))    
     &                     /(DR(J+2)+2.0D0*(DR(J+1)+DR(J))+DR(J-1))     
     &       )                                                          
     &       *(U(I,J,K,2)-U(I,J,K-1,2))*PDZ                             
     &      )*0.5D0                                                     
     &      -                                                           
     &      (1.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*(DR(J)*U(I,J+1,K+1,3)+DR(J-1)*U(I,J,K+1,3))  
     &                     *PDR2(J)-                                    
     &      (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I,J+2,K+1,3)         
     &                     +(2.0D0*DR(J+1)+DR(J+2))*U(I,J-1,K+1,3))     
     &       )                                                          
     &      *(U(I,J,K+3,2)-U(I,J,K,2))/3.0D0*PDZ                        
     &      +                                                           
     &      ((9.0D0/8.0D0)*(DR(J)*U(I,J+1,K-2,3)+DR(J-1)*U(I,J,K-2,3))  
     &                     *PDR2(J)-                                    
     &      (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I,J+2,K-2,3)         
     &                     +(2.0D0*DR(J+1)+DR(J+2))*U(I,J-1,K-2,3))     
     &       )                                                          
     &      *(U(I,J,K,2)-U(I,J,K-3,2))/3.0D0*PDZ                        
     &      )*0.5D0                                                     
        CONY=VC1+VC2+VC3                                                
C                                                                       
        UVISY=(-U(I+2,J,K,2)+16.0D0*U(I+1,J,K,2)-30.0D0*U(I,J,K,2)      
     &          +16.0D0*U(I-1,J,K,2)-U(I-2,J,K,2))/12.0D0               
     &          *PDX2*PRE                                               
C                                                                       
        WVISY=(-U(I,J,K+2,2)+16.0D0*U(I,J,K+1,2)-30.0D0*U(I,J,K,2)      
     &          +16.0D0*U(I,J,K-1,2)-U(I,J,K-2,2))/12.0D0               
     &          *PDZ2*PRE                                               
C                                                                       
        VVISY=(DR(J+1)*U(I,J-1,K,2)-DR2(J)*U(I,J,K,2)                   
     &          +DR(J)*U(I,J+1,K,2))*2.0*PDR3(J)                        
C                                                                       
        FFY(I,J,K)=2.0D0*RE*(CONY-UVISY-WVISY)-VVISY                    
     &             -2.0D0*RE/DT*U(I,J,K,2)                              
C                                                                       
        CNY_OLD(I,J,K)=CONY-UVISY-WVISY                                 
C                                                                       
   70 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     CALC. TEMPORARY VELOCIYT (VELOCITY FIELD)                         
      SUBROUTINE FS1                                                    
C     1998.01.01                                                        
C====================================================================   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
C                                                                       
C----LOCAL VARIABLES-------------------------------------------         
C                                                                       
      INCLUDE './INCFILE/INCDATA'                                         
      INCLUDE './INCFILE/INCMESH'                                         
      DIMENSION U(-4:IG+4,-2:JG+2,-4:KG+4,3)                            
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /VELO/ UG(-4:IG+4,-2:JG+2,-4:KG+4,3)                       
      COMMON /OLDF/ FFX(0:IG,0:JG,0:KG),FFY(0:IG,0:JG,0:KG)             
     &             ,FFZ(0:IG,0:JG,0:KG)                                 
      COMMON /OLDCN/ CNX_OLD(0:IG,0:JG,0:KG),CNY_OLD(0:IG,0:JG,0:KG)    
     &              ,CNZ_OLD(0:IG,0:JG,0:KG)                            
C                                                                       
C----COMBINATION------------------------------------------------        
C                                                                       
      EQUIVALENCE (UG,U)                                                
C                                                                       
C=================================================================      
C                                                                       
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,UC1,UC2,UC3,CONX,UVISX,WVISX,VVISX        
!$OMP&                   ,WC1,WC2,WC3,CONZ,UVISZ,WVISZ,VVISZ)           
      DO 10 J=1,JG                                                      
      DO 10 K=1,KG                                                      
      DO 10 I=1,IG                                                      
C--- STREAMWISE DIRECTION                                               
         UC1=                                                           
     &      (9.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I,J,K,1)+U(I-1,J,K,1))*0.5D0)-           
     &      (1.0D0/8.0D0)*((U(I+1,J,K,1)+U(I-2,J,K,1))*0.5D0))*         
     &      (U(I,J,K,1)-U(I-1,J,K,1))*PDX                               
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I,J,K,1)+U(I+1,J,K,1))*0.5D0)-           
     &      (1.0D0/8.0D0)*((U(I+2,J,K,1)+U(I-1,J,K,1))*0.5D0))*         
     &      (U(I+1,J,K,1)-U(I,J,K,1))*PDX                               
     &      )*0.5D0                                                     
     &      -                                                           
     &      (1.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I-1,J,K,1)+U(I-2,J,K,1))*0.5D0)-         
     &      (1.0D0/8.0D0)*((U(I,J,K,1)+U(I-3,J,K,1))*0.5D0))*           
     &      (U(I,J,K,1)-U(I-3,J,K,1))/3.0D0*PDX                         
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I+1,J,K,1)+U(I+2,J,K,1))*0.5D0)-         
     &      (1.0D0/8.0D0)*((U(I,J,K,1)+U(I+3,J,K,1))*0.5D0))*           
     &      (U(I+3,J,K,1)-U(I,J,K,1))/3.0D0*PDX                         
     &      )*0.5D0                                                     
        UC2=0.25                                                        
     &      *((U(I+1,J  ,K,2)+U(I,J  ,K,2))                             
     &       *(U(I  ,J+1,K,1)-U(I,J  ,K,1))*PDRR(J+1)                   
     &       +(U(I+1,J-1,K,2)+U(I,J-1,K,2))                             
     &       *(U(I  ,J  ,K,1)-U(I,J-1,K,1))*PDRR(J))                    
         UC3=                                                           
     &      (9.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I+1,J,K,3)+U(I,J,K,3))*0.5D0)-           
     &      (1.0D0/8.0D0)*((U(I+2,J,K,3)+U(I-1,J,K,3))*0.5D0))*         
     &      (U(I,J,K+1,1)-U(I,J,K,1))*PDZ                               
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I+1,J,K-1,3)+U(I,J,K-1,3))*0.5D0)-       
     &      (1.0D0/8.0D0)*((U(I+2,J,K-1,3)+U(I-1,J,K-1,3))*0.5D0))*     
     &      (U(I,J,K,1)-U(I,J,K-1,1))*PDZ                               
     &      )*0.5D0                                                     
     &      -                                                           
     &      (1.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I+1,J,K+1,3)+U(I,J,K+1,3))*0.5D0)-       
     &      (1.0D0/8.0D0)*((U(I+2,J,K+1,3)+U(I-1,J,K+1,3))*0.5D0))*     
     &      (U(I,J,K+3,1)-U(I,J,K,1))/3.0D0*PDZ                         
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I+1,J,K-2,3)+U(I,J,K-2,3))*0.5D0)-       
     &      (1.0D0/8.0D0)*((U(I+2,J,K-2,3)+U(I-1,J,K-2,3))*0.5D0))*     
     &      (U(I,J,K,1)-U(I,J,K-3,1))/3.0D0*PDZ                         
     &      )*0.5D0                                                     
        CONX=UC1+UC2+UC3                                                
C                                                                       
        UVISX=(-U(I+2,J,K,1)+16.0D0*U(I+1,J,K,1)-30.0D0*U(I,J,K,1)      
     &          +16.0D0*U(I-1,J,K,1)-U(I-2,J,K,1))/12.0D0               
     &          *PDX2*PRE                                               
C                                                                       
        WVISX=(-U(I,J,K+2,1)+16.0D0*U(I,J,K+1,1)-30.0D0*U(I,J,K,1)      
     &          +16.0D0*U(I,J,K-1,1)-U(I,J,K-2,1))/12.0D0               
     &          *PDZ2*PRE                                               
C                                                                       
        VVISX=(DRR(J+1)*U(I,J-1,K,1)-DRR2(J)*U(I,J,K,1)                 
     &          +DRR(J)*U(I,J+1,K,1))*2.0*PDRR3(J)                      
C                                                                       
        FFX(I,J,K)=RE*(3.0D0*(CONX-UVISX-WVISX)-CNX_OLD(I,J,K))         
     &             -VVISX-2.0D0*RE/DT*U(I,J,K,1)-2.0D0*RE*2.0D0         
C                                                                       
        CNX_OLD(I,J,K)=CONX-UVISX-WVISX                                 
C                                                                       
C--- SPANWISE DIRECTION                                                 
C                                                                       
         WC1=                                                           
     &      (9.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I,J,K+1,1)+U(I,J,K,1))*0.5D0)-           
     &      (1.0D0/8.0D0)*((U(I,J,K+2,1)+U(I,J,K-1,1))*0.5D0))*         
     &      (U(I+1,J,K,3)-U(I,J,K,3))*PDX                               
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I-1,J,K+1,1)+U(I-1,J,K,1))*0.5D0)-       
     &      (1.0D0/8.0D0)*((U(I-1,J,K+2,1)+U(I-1,J,K-1,1))*0.5D0))*     
     &      (U(I,J,K,3)-U(I-1,J,K,3))*PDX                               
     &      )*0.5D0                                                     
     &      -                                                           
     &      (1.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I+1,J,K+1,1)+U(I+1,J,K,1))*0.5D0)-       
     &      (1.0D0/8.0D0)*((U(I+1,J,K+2,1)+U(I+1,J,K-1,1))*0.5D0))*     
     &      (U(I+3,J,K,3)-U(I,J,K,3))/3.0D0*PDX                         
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I-2,J,K+1,1)+U(I-2,J,K,1))*0.5D0)-       
     &      (1.0D0/8.0D0)*((U(I-2,J,K+2,1)+U(I-2,J,K-1,1))*0.5D0))*     
     &      (U(I,J,K,3)-U(I-3,J,K,3))/3.0D0*PDX                         
     &      )*0.5D0                                                     
        WC2=0.25                                                        
     &      *((U(I,J  ,K+1,2)+U(I,J  ,K,2))                             
     &       *(U(I,J+1,K  ,3)-U(I,J  ,K,3))*PDRR(J+1)                   
     &       +(U(I,J-1,K+1,2)+U(I,J-1,K,2))                             
     &       *(U(I,J  ,K  ,3)-U(I,J-1,K,3))*PDRR(J))                    
         WC3=                                                           
     &      (9.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I,J,K+1,3)+U(I,J,K,3))*0.5D0)-           
     &      (1.0D0/8.0D0)*((U(I,J,K+2,3)+U(I,J,K-1,3))*0.5D0))*         
     &      (U(I,J,K+1,3)-U(I,J,K,3))*PDZ                               
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I,J,K,3)+U(I,J,K-1,3))*0.5D0)-           
     &      (1.0D0/8.0D0)*((U(I,J,K+1,3)+U(I,J,K-2,3))*0.5D0))*         
     &      (U(I,J,K,3)-U(I,J,K-1,3))*PDZ                               
     &      )*0.5D0                                                     
     &      -                                                           
     &      (1.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*((U(I,J,K+2,3)+U(I,J,K+1,3))*0.5D0)-         
     &      (1.0D0/8.0D0)*((U(I,J,K+3,3)+U(I,J,K,3))*0.5D0))*           
     &      (U(I,J,K+3,3)-U(I,J,K,3))/3.0D0*PDZ                         
     &      +                                                           
     &      ((9.0D0/8.0D0)*((U(I,J,K-1,3)+U(I,J,K-2,3))*0.5D0)-         
     &      (1.0D0/8.0D0)*((U(I,J,K,3)+U(I,J,K-3,3))*0.5D0))*           
     &      (U(I,J,K,3)-U(I,J,K-3,3))/3.0D0*PDZ                         
     &      )*0.5D0                                                     
        CONZ=WC1+WC2+WC3                                                
C                                                                       
        UVISZ=(-U(I+2,J,K,3)+16.0D0*U(I+1,J,K,3)-30.0D0*U(I,J,K,3)      
     &          +16.0D0*U(I-1,J,K,3)-U(I-2,J,K,3))/12.0D0               
     &          *PDX2*PRE                                               
C                                                                       
        WVISZ=(-U(I,J,K+2,3)+16.0D0*U(I,J,K+1,3)-30.0D0*U(I,J,K,3)      
     &          +16.0D0*U(I,J,K-1,3)-U(I,J,K-2,3))/12.0D0               
     &          *PDZ2*PRE                                               
C                                                                       
        VVISZ=(DRR(J+1)*U(I,J-1,K,3)-DRR2(J)*U(I,J,K,3)                 
     &          +DRR(J)*U(I,J+1,K,3))*2.0*PDRR3(J)                      
C                                                                       
        FFZ(I,J,K)=RE*(3.0D0*(CONZ-UVISZ-WVISZ)-CNZ_OLD(I,J,K))         
     &             -VVISZ-2.0D0*RE/DT*U(I,J,K,3)                        
C                                                                       
        CNZ_OLD(I,J,K)=CONZ-UVISZ-WVISZ                                 
C                                                                       
   10 CONTINUE                                                          
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,VC1,VC2,VC3,CONY,UVISY,WVISY,VVISY)       
      DO 20 J=1,JG-1                                                    
      DO 20 K=1,KG                                                      
      DO 20 I=1,IG                                                      
C--- WALL TO WALL DIRECTION                                             
         VC1=                                                           
     &      (9.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*(DR(J)*U(I,J+1,K,1)+DR(J+1)*U(I,J,K,1))      
     &                     *PDR2(J)-                                    
     &       (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I,J+2,K,1)          
     &                      +(2.0D0*DR(J+1)+DR(J+2))*U(I,J-1,K,1))      
     &                     /(DR(J+2)+2.0D0*(DR(J+1)+DR(J))+DR(J-1))     
     &       )                                                          
     &       *(U(I+1,J,K,2)-U(I,J,K,2))*PDX                             
     &      +                                                           
     &      ((9.0D0/8.0D0)*(DR(J)*U(I-1,J+1,K,1)+DR(J+1)*U(I-1,J,K,1))  
     &                     *PDR2(J)-                                    
     &       (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I-1,J+2,K,1)        
     &                      +(2.0D0*DR(J+1)+DR(J+2))*U(I-1,J-1,K,1))    
     &                     /(DR(J+2)+2.0D0*(DR(J+1)+DR(J))+DR(J-1))     
     &       )                                                          
     &       *(U(I,J,K,2)-U(I-1,J,K,2))*PDX                             
     &      )*0.5D0                                                     
     &      -                                                           
     &      (1.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*(DR(J)*U(I+1,J+1,K,1)+DR(J-1)*U(I+1,J,K,1))  
     &                     *PDR2(J)-                                    
     &      (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I+1,J+2,K,1)         
     &                     +(2.0D0*DR(J+1)+DR(J+2))*U(I+1,J-1,K,1))     
     &       )                                                          
     &      *(U(I+3,J,K,2)-U(I,J,K,2))/3.0D0*PDX                        
     &      +                                                           
     &      ((9.0D0/8.0D0)*(DR(J)*U(I-2,J+1,K,1)+DR(J-1)*U(I-2,J,K,1))  
     &                     *PDR2(J)-                                    
     &      (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I-2,J+2,K,1)         
     &                     +(2.0D0*DR(J+1)+DR(J+2))*U(I-2,J-1,K,1))     
     &       )                                                          
     &      *(U(I,J,K,2)-U(I-3,J,K,2))/3.0D0*PDX                        
     &      )*0.5D0                                                     
        VC2=0.5*PDR2(J)                                                 
     &      *((U(I,J+1,K,2)+U(I,J,K,2))                                 
     &       *(U(I,J+1,K,2)-U(I,J,K,2))*PDR(J+1)*DR(J)                  
     &       +(U(I,J,K,2)+U(I,J-1,K,2))                                 
     &       *(U(I,J,K,2)-U(I,J-1,K,2))*PDR(J)*DR(J+1))                 
         VC3=                                                           
     &      (9.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*(DR(J)*U(I,J+1,K,3)+DR(J+1)*U(I,J,K,3))      
     &                     *PDR2(J)-                                    
     &       (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I,J+2,K,3)          
     &                      +(2.0D0*DR(J+1)+DR(J+2))*U(I,J-1,K,3))      
     &                     /(DR(J+2)+2.0D0*(DR(J+1)+DR(J))+DR(J-1))     
     &       )                                                          
     &       *(U(I,J,K+1,2)-U(I,J,K,2))*PDZ                             
     &      +                                                           
     &      ((9.0D0/8.0D0)*(DR(J)*U(I,J+1,K-1,3)+DR(J+1)*U(I,J,K-1,3))  
     &                     *PDR2(J)-                                    
     &       (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I,J+2,K-1,3)        
     &                      +(2.0D0*DR(J+1)+DR(J+2))*U(I,J-1,K-1,3))    
     &                     /(DR(J+2)+2.0D0*(DR(J+1)+DR(J))+DR(J-1))     
     &       )                                                          
     &       *(U(I,J,K,2)-U(I,J,K-1,2))*PDZ                             
     &      )*0.5D0                                                     
     &      -                                                           
     &      (1.0D0/8.0D0)*(                                             
     &      ((9.0D0/8.0D0)*(DR(J)*U(I,J+1,K+1,3)+DR(J-1)*U(I,J,K+1,3))  
     &                     *PDR2(J)-                                    
     &      (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I,J+2,K+1,3)         
     &                     +(2.0D0*DR(J+1)+DR(J+2))*U(I,J-1,K+1,3))     
     &       )                                                          
     &      *(U(I,J,K+3,2)-U(I,J,K,2))/3.0D0*PDZ                        
     &      +                                                           
     &      ((9.0D0/8.0D0)*(DR(J)*U(I,J+1,K-2,3)+DR(J-1)*U(I,J,K-2,3))  
     &                     *PDR2(J)-                                    
     &      (1.0D0/8.0D0)*((2.0D0*DR(J)+DR(J-1))*U(I,J+2,K-2,3)         
     &                     +(2.0D0*DR(J+1)+DR(J+2))*U(I,J-1,K-2,3))     
     &       )                                                          
     &      *(U(I,J,K,2)-U(I,J,K-3,2))/3.0D0*PDZ                        
     &      )*0.5D0                                                     
        CONY=VC1+VC2+VC3                                                
C                                                                       
        UVISY=(-U(I+2,J,K,2)+16.0D0*U(I+1,J,K,2)-30.0D0*U(I,J,K,2)      
     &          +16.0D0*U(I-1,J,K,2)-U(I-2,J,K,2))/12.0D0               
     &          *PDX2*PRE                                               
C                                                                       
        WVISY=(-U(I,J,K+2,2)+16.0D0*U(I,J,K+1,2)-30.0D0*U(I,J,K,2)      
     &          +16.0D0*U(I,J,K-1,2)-U(I,J,K-2,2))/12.0D0               
     &          *PDZ2*PRE                                               
C                                                                       
        VVISY=(DR(J+1)*U(I,J-1,K,2)-DR2(J)*U(I,J,K,2)                   
     &          +DR(J)*U(I,J+1,K,2))*2.0*PDR3(J)                        
C                                                                       
        FFY(I,J,K)=RE*(3.0D0*(CONY-UVISY-WVISY)-CNY_OLD(I,J,K))         
     &             -VVISY-2.0D0*RE/DT*U(I,J,K,2)                        
C                                                                       
        CNY_OLD(I,J,K)=CONY-UVISY-WVISY                                 
C                                                                       
   20 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     CALC. TEMPORARY VELOCITY WITH THE USE OF TDMA (Y-DIRECTION)       
C     1997.10.31 H.ABE                                                  
      SUBROUTINE CNTDMA                                                 
C     1998.01.01                                                        
C====================================================================   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
      INCLUDE './INCFILE/INCDATA'                                         
      INCLUDE './INCFILE/INCMESH'                                         
C                                                                       
C----LOCAL VARIABLES-------------------------------------------         
C                                                                       
      COMMON /TDMAXZ/ AAXZ(LY),BBXZ(LY),CCXZ(LY)                        
      COMMON /TDMAY/ AAY(LY),BBY(LY),CCY(LY)                            
      DIMENSION UT(-4:IG+4,-2:JG+2,-4:KG+4,3)                           
      DIMENSION B(0:IG,JG),D(0:IG,JG,2),E(0:IG,JG)                      
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /TEMP/ UTG(-4:IG+4,-2:JG+2,-4:KG+4,3)                      
      COMMON /OLDF/ FFX(0:IG,0:JG,0:KG),FFY(0:IG,0:JG,0:KG)             
     &             ,FFZ(0:IG,0:JG,0:KG)                                 
C                                                                       
C----COMBINATION------------------------------------------------        
C                                                                       
      EQUIVALENCE (UTG,UT)                                              
C==============================================================         
C----- X,Z-DIRECTION                                                    
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,J,B,D)                                      
      DO 110 K=1,KG                                                     
       DO 120 J=1,JG                                                    
       DO 120 I=1,IG                                                    
        B(I,J)=BBXZ(J)                                                  
        D(I,J,1)=FFX(I,J,K)                                             
        D(I,J,2)=FFZ(I,J,K)                                             
 120   CONTINUE                                                         
C                                                                       
       DO 130 I=1,IG                                                    
        B(I,1)=BBXZ(1)+CCXZ(1)                                          
        B(I,JG)=BBXZ(JG)+AAXZ(JG)                                       
 130   CONTINUE                                                         
C                                                                       
       CALL TDMAUW(AAXZ,B,CCXZ,D)                                       
C                                                                       
       DO 140 J=1,JG                                                    
       DO 140 I=1,IG                                                    
        UT(I,J,K,1)=D(I,J,1)                                            
        UT(I,J,K,3)=D(I,J,2)                                            
 140   CONTINUE                                                         
 110  CONTINUE                                                          
C                                                                       
C----- Y-DIRECTION                                                      
C                                                                       
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,J,B,E)                                      
      DO 300 K=1,KG                                                     
       DO 310 J=1,JG                                                    
       DO 310 I=1,IG                                                    
        B(I,J)=BBY(J)                                                   
        E(I,J)=FFY(I,J,K)                                               
 310   CONTINUE                                                         
C                                                                       
       CALL TDMAV(AAY,B,CCY,E)                                          
C                                                                       
       DO 320 J=1,JG                                                    
       DO 320 I=1,IG                                                    
        UT(I,J,K,2)=E(I,J)                                              
 320   CONTINUE                                                         
 300  CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     CONTENTS : TDMA (Y-DIRECTION) SUBROUTINE                          
C     1998.01.01 H.ABE                                                  
C     CALC. TEMPORARY VELOCITY WITH THE USE OF TDMA (Y-DIRECTION)       
C     1997.12.22 H.ABE                                                  
      SUBROUTINE TDMAUW(A,B,C,D)                                        
C==============================================================         
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
C                                                                       
      DIMENSION A(JG+1),B(0:IG,JG),C(JG+1),D(0:IG,JG,2),E(0:IG,JG)      
C                                                                       
      DO 200 I=1,IG                                                     
       E(I,1)=A(1)/B(I,1)                                               
       D(I,1,1)=-D(I,1,1)/B(I,1)                                        
       D(I,1,2)=-D(I,1,2)/B(I,1)                                        
  200 CONTINUE                                                          
C                                                                       
      DO 210 J=2,JG                                                     
      DO 210 I=1,IG                                                     
       B(I,J)=B(I,J)-C(J)*E(I,J-1)                                      
       E(I,J)=A(J)/B(I,J)                                               
       D(I,J,1)=(-D(I,J,1)+C(J)*D(I,J-1,1))/B(I,J)                      
       D(I,J,2)=(-D(I,J,2)+C(J)*D(I,J-1,2))/B(I,J)                      
  210 CONTINUE                                                          
C                                                                       
      DO 220 J=JG-1,1,-1                                                
      DO 220 I=1,IG                                                     
       D(I,J,1)=D(I,J,1)+D(I,J+1,1)*E(I,J)                              
       D(I,J,2)=D(I,J,2)+D(I,J+1,2)*E(I,J)                              
  220 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C==============================================================         
C     CALC. TEMPORARY VELOCITY WITH THE USE OF TDMA (Y-DIRECTION)       
C     1997.12.22 H.ABE                                                  
      SUBROUTINE TDMAV(A,B,C,D)                                         
C==============================================================         
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
C                                                                       
      DIMENSION A(JG+1),B(0:IG,JG),C(JG+1),D(0:IG,JG),E(0:IG,JG)        
C                                                                       
      DO 400 I=1,IG                                                     
       E(I,1)=A(1)/B(I,1)                                               
       D(I,1)=-D(I,1)/B(I,1)                                            
  400 CONTINUE                                                          
C                                                                       
      DO 410 J=2,JG                                                     
      DO 410 I=1,IG                                                     
       B(I,J)=B(I,J)-C(J)*E(I,J-1)                                      
       E(I,J)=A(J)/B(I,J)                                               
       D(I,J)=(-D(I,J)+C(J)*D(I,J-1))/B(I,J)                            
  410 CONTINUE                                                          
C                                                                       
      DO 420 I=1,IG                                                     
       D(I,JG)=0.0D0                                                    
  420 CONTINUE                                                          
C                                                                       
      DO 430 J=JG-1,1,-1                                                
      DO 430 I=1,IG                                                     
        D(I,J)=D(I,J)+D(I,J+1)*E(I,J)                                   
  430 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     POISSON EQUATION FOR NWT                                          
C     SUBROUTINE (FS2,FFT3DJ,STDMA)                                     
C     1998.01.01                                                        
C=================================================== POISSON EQ. FOR NWT
C     FOR Crank-Nicolson (y-direction)                                  
C     1998.01.03                                                        
C     Programed by H.Abe                                                
      SUBROUTINE FS2                                                    
C====================================================================   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
C                                                                       
C----LOCAL VARIABLES-------------------------------------------         
C                                                                       
      INCLUDE './INCFILE/INCDATA'                                         
      INCLUDE './INCFILE/INCMESH'                                         
      COMMON /FD2C/ CEN2_F(0:JG,2),CEN2_M(0:JG,2),CEN2_B(0:JG,2)        
      DIMENSION CPR(0:IG,JG,0:KG),CQR(0:IG,JG,0:KG)                     
     &         ,CPI(0:IG,JG,0:KG),CQI(0:IG,JG,0:KG)                     
      DIMENSION B(0:IG,JG),D(0:IG,JG,2)                                 
      DIMENSION U(-4:IG+4,-2:JG+2,-4:KG+4,3)                            
      DIMENSION UT(-4:IG+4,-2:JG+2,-4:KG+4,3)                           
      DIMENSION P(-3:IG+3,-2:JG+2,-3:KG+3)                              
      DIMENSION PHI(-3:IG+3,-2:JG+2,-3:KG+3)                            
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /VELO/ UG(-4:IG+4,-2:JG+2,-4:KG+4,3)                       
      COMMON /TEMP/ UTG(-4:IG+4,-2:JG+2,-4:KG+4,3)                      
      COMMON /SCAL/ PG(-3:IG+3,-2:JG+2,-3:KG+3)                         
      COMMON /PHIP/ PHIG(-3:IG+3,-2:JG+2,-3:KG+3)                       
      COMMON /TDMA/ AA(LY),BB(LY),CC(LY),RK(0:LX,0:KG)                  
      COMMON /FFT3/ WA(IG+1,KG,JG,2)                                    
C                                                                       
C----COMBINATION------------------------------------------------        
C                                                                       
      EQUIVALENCE (UG,U),(UTG,UT),(PG,P)                                
      EQUIVALENCE (PHIG,PHI)                                            
C                                                                       
C=================================================================      
C                                                                       
C                                                                       
C---- W-W DIRECTION B.C. NON-SLIP                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
      DO 10 K=1,KG                                                      
      DO 10 I=1,IG                                                      
        UT(I,0 ,K,2)=0.0D0                                              
        UT(I,JG,K,2)=0.0D0                                              
   10 CONTINUE                                                          
C----------------------------------------------------                   
C---- STREAM DIRECTION B.C. PERIODIC                                    
!$OMP PARALLEL DO PRIVATE(K)                                            
      DO 11 J=-2,JG+2                                                   
      DO 11 K=-4,KG+4                                                   
        UT(-4,J,K,1)=UT(IG-4,J,K,1)                                     
        UT(-3,J,K,1)=UT(IG-3,J,K,1)                                     
        UT(-2,J,K,1)=UT(IG-2,J,K,1)                                     
        UT(-1,J,K,1)=UT(IG-1,J,K,1)                                     
        UT(0,J,K,1)=UT(IG,J,K,1)                                        
        UT(IG+1,J,K,1)=UT(1,J,K,1)                                      
        UT(IG+2,J,K,1)=UT(2,J,K,1)                                      
        UT(IG+3,J,K,1)=UT(3,J,K,1)                                      
        UT(IG+4,J,K,1)=UT(4,J,K,1)                                      
        UT(-4,J,K,2)=UT(IG-4,J,K,2)                                     
        UT(-3,J,K,2)=UT(IG-3,J,K,2)                                     
        UT(-2,J,K,2)=UT(IG-2,J,K,2)                                     
        UT(-1,J,K,2)=UT(IG-1,J,K,2)                                     
        UT( 0,J,K,2)=UT(IG,J,K,2)                                       
        UT(IG+1,J,K,2)=UT(1,J,K,2)                                      
        UT(IG+2,J,K,2)=UT(2,J,K,2)                                      
        UT(IG+3,J,K,2)=UT(3,J,K,2)                                      
        UT(IG+4,J,K,2)=UT(4,J,K,2)                                      
        UT(-4,J,K,3)=UT(IG-4,J,K,3)                                     
        UT(-3,J,K,3)=UT(IG-3,J,K,3)                                     
        UT(-2,J,K,3)=UT(IG-2,J,K,3)                                     
        UT(-1,J,K,3)=UT(IG-1,J,K,3)                                     
        UT( 0,J,K,3)=UT(IG,J,K,3)                                       
        UT(IG+1,J,K,3)=UT(1,J,K,3)                                      
        UT(IG+2,J,K,3)=UT(2,J,K,3)                                      
        UT(IG+3,J,K,3)=UT(3,J,K,3)                                      
        UT(IG+4,J,K,3)=UT(4,J,K,3)                                      
   11 CONTINUE                                                          
C---- SPAN DIRECTION B.C. PERIODIC                                      
!$OMP PARALLEL DO PRIVATE(I)                                            
      DO 12 J=-2,JG+2                                                   
      DO 12 I=-4,IG+4                                                   
        UT(I,J,-4,1)=UT(I,J,KG-4,1)                                     
        UT(I,J,-3,1)=UT(I,J,KG-3,1)                                     
        UT(I,J,-2,1)=UT(I,J,KG-2,1)                                     
        UT(I,J,-1,1)=UT(I,J,KG-1,1)                                     
        UT(I,J,0,1)=UT(I,J,KG,1)                                        
        UT(I,J,KG+1,1)=UT(I,J,1,1)                                      
        UT(I,J,KG+2,1)=UT(I,J,2,1)                                      
        UT(I,J,KG+3,1)=UT(I,J,3,1)                                      
        UT(I,J,KG+4,1)=UT(I,J,4,1)                                      
        UT(I,J,-4,2)=UT(I,J,KG-4,2)                                     
        UT(I,J,-3,2)=UT(I,J,KG-3,2)                                     
        UT(I,J,-2,2)=UT(I,J,KG-2,2)                                     
        UT(I,J,-1,2)=UT(I,J,KG-1,2)                                     
        UT(I,J,0,2)=UT(I,J,KG,2)                                        
        UT(I,J,KG+1,2)=UT(I,J,1,2)                                      
        UT(I,J,KG+2,2)=UT(I,J,2,2)                                      
        UT(I,J,KG+3,2)=UT(I,J,3,2)                                      
        UT(I,J,KG+4,2)=UT(I,J,4,2)                                      
        UT(I,J,-4,3)=UT(I,J,KG-4,3)                                     
        UT(I,J,-3,3)=UT(I,J,KG-3,3)                                     
        UT(I,J,-2,3)=UT(I,J,KG-2,3)                                     
        UT(I,J,-1,3)=UT(I,J,KG-1,3)                                     
        UT(I,J,0,3)=UT(I,J,KG,3)                                        
        UT(I,J,KG+1,3)=UT(I,J,1,3)                                      
        UT(I,J,KG+2,3)=UT(I,J,2,3)                                      
        UT(I,J,KG+3,3)=UT(I,J,3,3)                                      
        UT(I,J,KG+4,3)=UT(I,J,4,3)                                      
   12 CONTINUE                                                          
C                                                                       
C                                                                       
C------------------------------------ SCHUMANN METHOD                   
C                                                                       
      PDT=1.0/DT                                                        
      W=1.0/DBLE(IG*KG)                                                 
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 110 J=1,JG                                                     
      DO 110 K=1,KG                                                     
      DO 110 I=1,IG                                                     
         WA(I,K,J,1)=((9.0D0/8.0D0)*(UT(I,J,K,1)-UT(I-1,J,K,1))*PDX     
     &                -(1.0D0/8.0D0)*(UT(I+1,J,K,1)-UT(I-2,J,K,1))      
     &                                /3.0D0*PDX                        
     &              +(UT(I,J,K,2)-UT(I,J-1,K,2))*PDR(J)                 
     &              +(9.0D0/8.0D0)*(UT(I,J,K,3)-UT(I,J,K-1,3))*PDZ      
     &                -(1.0D0/8.0D0)*(UT(I,J,K+1,3)-UT(I,J,K-2,3))      
     &                                /3.0D0*PDZ                        
     &                )*PDT                                             
C        WA(I,K,J,1)=((UT(I,J,K,1)-UT(I-1,J,K,1))*PDX                   
C     &              +(UT(I,J,K,2)-UT(I,J-1,K,2))*PDR(J)                
C     &              +(UT(I,J,K,3)-UT(I,J,K-1,3))*PDZ)*PDT              
        WA(I,K,J,2)=0.0D0                                               
  110 CONTINUE                                                          
C                                                                       
C                                                                       
      ISN=1                                                             
      CALL FFT3DJ(ISN)                                                  
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 112 J=1,JG                                                     
      DO 112 K=1,KG                                                     
      DO 112 I=1,IG                                                     
        WA(I,K,J,1)=-WA(I,K,J,1)*W                                      
        WA(I,K,J,2)=-WA(I,K,J,2)*W                                      
  112 CONTINUE                                                          
C--------------------------------                                       
!$OMP PARALLEL DO PRIVATE(I,J)                                          
        DO 113 K=1,KG                                                   
        DO 113 J=1,JG                                                   
        DO 113 I=1,IG                                                   
          CQR(I,J,K)=WA(I,K,J,1)                                        
          CQI(I,J,K)=WA(I,K,J,2)                                        
  113   CONTINUE                                                        
C                                                                       
C--------------------------------                                       
C-- L=1 AND N=1                                                         
      CPR(1,1,1)=0.0D0                                                  
      CPI(1,1,1)=0.0D0                                                  
      CPR(1,2,1)=CQR(1,1,1)/CC(1)                                       
      CPI(1,2,1)=CQI(1,1,1)/CC(1)                                       
      DO 131 J=3,JG                                                     
        CPR(1,J,1)=(-BB(J-1)*CPR(1,J-1,1)-AA(J-1)*CPR(1,J-2,1)          
     &             +CQR(1,J-1,1))/CC(J-1)                               
        CPI(1,J,1)=(-BB(J-1)*CPI(1,J-1,1)-AA(J-1)*CPI(1,J-2,1)          
     &             +CQI(1,J-1,1))/CC(J-1)                               
  131 CONTINUE                                                          
C                                                                       
C-- L=2,IG AND N=1                                                      
      N=1                                                               
C     IS=2                                                              
!$OMP PARALLEL DO PRIVATE(L)                                            
      DO 240 J=1,JG                                                     
        DO 240 L=2,IG                                                   
          B(L,J)=BB(J)+RK(L-1,N-1)                                      
          D(L,J,1)=CQR(L,J,N)                                           
          D(L,J,2)=CQI(L,J,N)                                           
  240 CONTINUE                                                          
      DO 250 L=2,IG                                                     
          B(L,1)=B(L,1)+AA(1)                                           
         B(L,JG)=B(L,JG)+CC(JG)                                         
  250 CONTINUE                                                          
C                                                                       
C                                                                       
      CALL STDMA(2,AA,B,CC,D)                                           
!$OMP PARALLEL DO PRIVATE(L)                                            
      DO 230 J=1,JG                                                     
        DO 230 L=2,IG                                                   
          CPR(L,J,N)=D(L,J,1)                                           
          CPI(L,J,N)=D(L,J,2)                                           
  230 CONTINUE                                                          
C                                                                       
C                                                                       
C-- L=1,IG AND N=2,KG                                                   
C     IS=1                                                              
!$OMP PARALLEL DO PRIVATE(J,L,B,D)                                      
      DO 320 N=2,KG                                                     
      DO 340 J=1,JG                                                     
        DO 340 L=1,IG                                                   
          B(L,J)=BB(J)+RK(L-1,N-1)                                      
          D(L,J,1)=CQR(L,J,N)                                           
          D(L,J,2)=CQI(L,J,N)                                           
  340 CONTINUE                                                          
      DO 350 L=1,IG                                                     
          B(L,1)=B(L,1)+AA(1)                                           
         B(L,JG)=B(L,JG)+CC(JG)                                         
  350 CONTINUE                                                          
      CALL STDMA(1,AA,B,CC,D)                                           
      DO 330 J=1,JG                                                     
        DO 330 L=1,IG                                                   
          CPR(L,J,N)=D(L,J,1)                                           
          CPI(L,J,N)=D(L,J,2)                                           
  330 CONTINUE                                                          
  320 CONTINUE                                                          
C                                                                       
C--------------------------------                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
        DO 171 J=1,JG                                                   
        DO 171 K=1,KG                                                   
        DO 171 I=1,IG                                                   
          WA(I,K,J,1)=CPR(I,J,K)                                        
          WA(I,K,J,2)=CPI(I,J,K)                                        
  171   CONTINUE                                                        
C                                                                       
C--------------------------------                                       
      ISN=-1                                                            
      CALL FFT3DJ(ISN)                                                  
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 172 J=1,JG                                                     
      DO 172 K=1,KG                                                     
      DO 172 I=1,IG                                                     
        PHI(I,J,K)=WA(I,K,J,1)                                          
  172 CONTINUE                                                          
C                                                                       
C                                                                       
C----------------------------- PHI B.C.                                 
C-- NEUMANN (V-V)                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
        DO 600 K=-3,KG+3                                                
        DO 600 I=-3,IG+3                                                
          PHI(I,-2,K)=PHI(I,3,K)                                        
          PHI(I,-1,K)=PHI(I,2,K)                                        
          PHI(I,0,K)=PHI(I,1,K)                                         
          PHI(I,JG+1,K)=PHI(I,JG,K)                                     
          PHI(I,JG+2,K)=PHI(I,JG-1,K)                                   
  600   CONTINUE                                                        
C----------------------------------------------------                   
C-- PERIODIC (SPANWISE)                                                 
!$OMP PARALLEL DO PRIVATE(I)                                            
      DO 41 J=1,JG                                                      
      DO 41 I=-3,IG+3                                                   
        PHI(I,J,KG+3)=PHI(I,J,3)                                        
        PHI(I,J,KG+2)=PHI(I,J,2)                                        
        PHI(I,J,KG+1)=PHI(I,J,1)                                        
        PHI(I,J,0)=PHI(I,J,KG)                                          
        PHI(I,J,-1)=PHI(I,J,KG-1)                                       
        PHI(I,J,-2)=PHI(I,J,KG-2)                                       
        PHI(I,J,-3)=PHI(I,J,KG-3)                                       
   41 CONTINUE                                                          
C-- PERIODIC (STREAMWISE)                                               
!$OMP PARALLEL DO PRIVATE(K)                                            
      DO 42 J=1,JG                                                      
      DO 42 K=-3,KG+3                                                   
        PHI(IG+3,J,K)=PHI(3,J,K)                                        
        PHI(IG+2,J,K)=PHI(2,J,K)                                        
        PHI(IG+1,J,K)=PHI(1,J,K)                                        
        PHI(0,J,K)=PHI(IG,J,K)                                          
        PHI(-1,J,K)=PHI(IG-1,J,K)                                       
        PHI(-2,J,K)=PHI(IG-2,J,K)                                       
        PHI(-3,J,K)=PHI(IG-3,J,K)                                       
   42 CONTINUE                                                          
C                                                                       
C                                                                       
C-------------------------------------- VELOCITY  NEW <- OLD            
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 50 J=1,JG                                                      
      DO 50 K=1,KG                                                      
      DO 50 I=1,IG                                                      
        U(I,J,K,1)=UT(I,J,K,1)-DT*(                                     
     &            (9.0D0/8.0D0)*(PHI(I+1,J,K)-PHI(I,J,K))*PDX           
     &           -(1.0D0/8.0D0)*(PHI(I+2,J,K)-PHI(I-1,J,K))/3.0D0*PDX   
     &                              )                                   
        U(I,J,K,2)=UT(I,J,K,2)-DT*(PHI(I,J+1,K)-PHI(I,J,K))*PDRR(J+1)   
        U(I,J,K,3)=UT(I,J,K,3)-DT*(                                     
     &            (9.0D0/8.0D0)*(PHI(I,J,K+1)-PHI(I,J,K))*PDZ           
     &           -(1.0D0/8.0D0)*(PHI(I,J,K+2)-PHI(I,J,K-1))/3.0D0*PDZ   
     &                              )                                   
   50 CONTINUE                                                          
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 51 J=1,JG                                                      
      DO 51 K=1,KG                                                      
      DO 51 I=1,IG                                                      
       P(I,J,K)=PHI(I,J,K)-0.5D0*DT*PRE*(                               
     &           CEN2_F(J,2)*PHI(I,J+1,K)+CEN2_M(J,2)*PHI(I,J,K)        
     &          +CEN2_B(J,2)*PHI(I,J-1,K))                              
   51 CONTINUE                                                          
C                                                                       
C                                                                       
C----------------------------- PRESSURE B.C.                            
C-- NEUMANN (V-V)                                                       
C                                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
      DO 620 K=-3,KG+3                                                  
      DO 620 I=-3,IG+3                                                  
        P(I,-2,K)=P(I,3,K)                                              
        P(I,-1,K)=P(I,2,K)                                              
        P(I,0,K)=P(I,1,K)                                               
        P(I,JG+1,K)=P(I,JG,K)                                           
        P(I,JG+2,K)=P(I,JG-1,K)                                         
  620 CONTINUE                                                          
C----------------------------------------------------                   
C-- PERIODIC (SPANWISE)                                                 
!$OMP PARALLEL DO PRIVATE(I)                                            
      DO 54 J=1,JG                                                      
      DO 54 I=-3,IG+3                                                   
        P(I,J,KG+3)=P(I,J,3)                                            
        P(I,J,KG+2)=P(I,J,2)                                            
        P(I,J,KG+1)=P(I,J,1)                                            
        P(I,J,0)=P(I,J,KG)                                              
        P(I,J,-1)=P(I,J,KG-1)                                           
        P(I,J,-2)=P(I,J,KG-2)                                           
        P(I,J,-3)=P(I,J,KG-3)                                           
   54 CONTINUE                                                          
C-- PERIODIC (STREAMWISE)                                               
!$OMP PARALLEL DO PRIVATE(K)                                            
      DO 55 J=1,JG                                                      
      DO 55 K=-3,KG+3                                                   
        P(IG+3,J,K)=P(3,J,K)                                            
        P(IG+2,J,K)=P(2,J,K)                                            
        P(IG+1,J,K)=P(1,J,K)                                            
        P(0,J,K)=P(IG,J,K)                                              
        P(-1,J,K)=P(IG-1,J,K)                                           
        P(-2,J,K)=P(IG-2,J,K)                                           
        P(-3,J,K)=P(IG-3,J,K)                                           
   55 CONTINUE                                                          
C                                                                       
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     FAST FOURIE TRANSFORMATION (VELOCITY FIELD)                       
      SUBROUTINE FFT3DJ(ISN)                                            
C     1998.01.03                                                        
C====================================================================   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
      PARAMETER (MX=LIG,MY=LKG)                                         
      PARAMETER (NX=IG,NY=KG,NR=JG)                                     
      PARAMETER (NX1=IG+1,NY1=KG+1)                                     
      PARAMETER (NP=IG,MP=LIG)                                          
C                                                                       
C----LOCAL VARIABLES-------------------------------------------         
C                                                                       
      DIMENSION WA1(NX1,NY,2),WA2(NX1,NY,2)                             
      DIMENSION WB1(NY1,NX,2),WB2(NY1,NX,2)                             
      DIMENSION WNP2(2,NP,2)                                            
      DIMENSION N2K(0:MP)                                               
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /FFT3/ WA(IG+1,KG,JG,2)                                    
C                                                                       
C=================================================================      
C                                                                       
      CALL FINIT( WNP2, N2K, NX, NY, NP, MP, ISN )                      
      IA = 1                                                            
      IB = 2                                                            
C                                                                       
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,J,WA1,WA2,WB1,WB2) FIRSTPRIVATE(IA,IB)      
      DO 100 IR = 1, NR                                                 
C                                                                       
        DO 10 I = 1, NY                                                 
            WA1(NX1,I,IA)  =  0.0D0                                     
            WA2(NX1,I,IA)  =  0.0D0                                     
   10 CONTINUE                                                          
        DO 20 I = 1, NY                                                 
          DO 20 J = 1, NX                                               
            WA1(J,I,IA)  =  WA(J,I,IR,1)                                
            WA2(J,I,IA)  =  WA(J,I,IR,2)                                
   20 CONTINUE                                                          
C                                                                       
      CALL FSTEP( WA1, WA2, WNP2(1,1,2)                                 
     1            , N2K, NX1, NY, MY, NP, MP, IA, IB )                  
C                                                                       
          DO 30 I = 1, NX                                               
            WB1(NY1,I,IA)  =  0.0D0                                     
            WB2(NY1,I,IA)  =  0.0D0                                     
   30 CONTINUE                                                          
        DO 40 I = 1, NY                                                 
          DO 40 J = 1, NX                                               
            WB1(I,J,IA)  =  WA1(J,I,IA)                                 
            WB2(I,J,IA)  =  WA2(J,I,IA)                                 
   40 CONTINUE                                                          
C                                                                       
      CALL FSTEP( WB1, WB2, WNP2(1,1,1)                                 
     1            , N2K, NY1, NX, MX, NP, MP, IA, IB )                  
C                                                                       
        DO 60 I = 1, NY                                                 
          DO 60 J = 1, NX                                               
            WA(J,I,IR,1)  =  WB1(I,J,IA)                                
            WA(J,I,IR,2)  =  WB2(I,J,IA)                                
   60 CONTINUE                                                          
  100 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C--------------------------------------------------------------         
      SUBROUTINE FINIT( WNP2, N2K, NX, NY, NP, MP, ISN )                
C--------------------------------------------------------------         
      INCLUDE './INCFILE/INCINIT96'                                       
      DIMENSION  WNP2( 2, NP, 2)                                        
      DIMENSION  N2K( 0 : MP )                                          
      PARAMETER( PI02 =  3.1415926535897933D0 * 2.0D0 )                 
C                                                                       
      DO  100  K = 0, MP                                                
        N2K( K ) = 2 ** K                                               
  100 CONTINUE                                                          
       IF( ISN.LT.0 ) THEN                                              
C        SK = PI02 / NX                                                 
         SK = PI02 / DBLE(NX)                                           
         DO 200 L = 1,2                                                 
          IF( L.GE.2 )THEN                                              
C           SK = PI02 / NY                                              
            SK = PI02 / DBLE(NY)                                        
           END IF                                                       
          DO 160 I = 1, NP                                              
            WORKS      =  SK * DBLE( I -1 )                             
            WNP2(1,I,L)   = DCOS( WORKS )                               
            WNP2(2,I,L)   = DSIN( WORKS )                               
  160 CONTINUE                                                          
  200 CONTINUE                                                          
        GO TO 500                                                       
       END IF                                                           
       IF(ISN.GT.0) THEN                                                
C          SK =  PI02 / NX                                              
           SK = PI02 / DBLE(NX)                                         
           DO 300 L = 1, 2                                              
            IF( L.GE.2 )THEN                                            
C             SK = PI02 / NY                                            
              SK = PI02 / DBLE(NY)                                      
             END IF                                                     
             DO 280 I = 1, NP                                           
              WORKS      =    SK * DBLE( I -1 )                         
              WNP2(1,I,L)   =   DCOS( WORKS )                           
              WNP2(2,I,L)   =  -DSIN( WORKS )                           
  280   CONTINUE                                                        
  300   CONTINUE                                                        
         GO TO 500                                                      
           END IF                                                       
         IF( ISN.EQ.0 ) THEN                                            
            WRITE(6,*) ' *****  ISN .EQ. 0, THEN JOB STOP  ******* '    
              STOP                                                      
        END IF                                                          
  500   CONTINUE                                                        
      RETURN                                                            
      END                                                               
C                                                                       
C  --------------------------------------------------------------       
      SUBROUTINE FSTEP( WA1, WA2, WNP, N2K, NX1, NY, MY, NP, MP         
     1                   , IA, IB )                                     
C  --------------------------------------------------------------       
      INCLUDE './INCFILE/INCINIT96'                                       
      DIMENSION  WA1( NX1, NY, 2 )                                      
      DIMENSION  WA2( NX1, NY, 2 )                                      
      DIMENSION  WNP( 2, NP)                                            
      DIMENSION  N2K( 0 : MP )                                          
        NX = NX1 - 1                                                    
        MM = MY                                                         
        NW2  =  NY / 2                                                  
C                                                                       
      IASAVE  =  IA                                                     
      IBSAVE  =  IB                                                     
C                                                                       
        DO 200 K = 1, MM-1                                              
          IIMAX  =  N2K( MM - K )                                       
          N2K1   =  N2K(  K - 1 )                                       
          IDP    =  N2K( MM - K + 1 )                                   
          IP01   =  0                                                   
          IP02   =  0                                                   
          DO 150 IL  = 1, N2K1                                          
            DO 100 I = 1, IIMAX                                         
              I1    =  I    + IP01                                      
              I2    =  I    + IP02                                      
              I3    =  I1   + NW2                                       
              I4    =  I2   + IIMAX                                     
              INP   =  N2K1 * (I-1) + 1                                 
              WNP1  =  WNP(1,INP)                                       
              WNP2  =  WNP(2,INP)                                       
              DO 80 J = 1, NX                                           
                WAI2              =  WA1(J,I2,IA)                       
                WAI4              =  WA1(J,I4,IA)                       
                WBI2              =  WA2(J,I2,IA)                       
                WBI4              =  WA2(J,I4,IA)                       
                WA1(J,I1,IB)  =  WAI2 + WAI4                            
                WA2(J,I1,IB)  =  WBI2 + WBI4                            
                WREAL1            =  WAI2 - WAI4                        
                WIMAG1            =  WBI2 - WBI4                        
                WA1(J,I3,IB)  =  WNP1 * WREAL1 - WNP2 * WIMAG1          
                WA2(J,I3,IB)  =  WNP1 * WIMAG1 + WNP2 * WREAL1          
   80         CONTINUE                                                  
  100       CONTINUE                                                    
            IP01  =  IP01 + IIMAX                                       
            IP02  =  IP02 + IDP                                         
  150     CONTINUE                                                      
          IWORK  =  IA                                                  
          IA     =  IB                                                  
          IB     =  IWORK                                               
  200 CONTINUE                                                          
C                                                                       
      IWORK  =  IA                                                      
      IA     =  IB                                                      
      IB     =  IWORK                                                   
      IIMAX              =  N2K( MM - 1 )                               
        DO 400 I = 1, IIMAX                                             
          I1  =  2 * I - 1                                              
          I2  =  I  + IIMAX                                             
          I3  =  I1 + 1                                                 
          DO 400 J = 1, NX                                              
            WAI1              =  WA1(J,I1,IB)                           
            WAI3              =  WA1(J,I3,IB)                           
            WBI1              =  WA2(J,I1,IB)                           
            WBI3              =  WA2(J,I3,IB)                           
            WA1(J, I,IA)  =  WAI1 + WAI3                                
            WA1(J,I2,IA)  =  WAI1 - WAI3                                
            WA2(J, I,IA)  =  WBI1 + WBI3                                
            WA2(J,I2,IA)  =  WBI1 - WBI3                                
  400 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C====================================================================   
C*********************************************** TDMA OF 2D-DATA FOR NWT
      SUBROUTINE STDMA(IS,A,B,C,D)                                      
      INCLUDE './INCFILE/INCINIT96'                                       
      DIMENSION A(JG+1),B(0:IG,JG),C(JG+1),D(0:IG,JG,2),E(0:IG,JG)      
C                                                                       
      DO 100 I=IS,IG                                                    
        E(I,1)=C(1)/B(I,1)                                              
        D(I,1,1)=D(I,1,1)/B(I,1)                                        
        D(I,1,2)=D(I,1,2)/B(I,1)                                        
  100 CONTINUE                                                          
C                                                                       
      DO 210 J=2,JG                                                     
        DO 110 I=IS,IG                                                  
          B(I,J)=B(I,J)-A(J)*E(I,J-1)                                   
          E(I,J)=C(J)/B(I,J)                                            
          D(I,J,1)=(D(I,J,1)-A(J)*D(I,J-1,1))/B(I,J)                    
          D(I,J,2)=(D(I,J,2)-A(J)*D(I,J-1,2))/B(I,J)                    
  110   CONTINUE                                                        
  210 CONTINUE                                                          
C                                                                       
      DO 220 J=JG-1,1,-1                                                
        DO 120 I=IS,IG                                                  
          D(I,J,1)=D(I,J,1)-D(I,J+1,1)*E(I,J)                           
          D(I,J,2)=D(I,J,2)-D(I,J+1,2)*E(I,J)                           
  120   CONTINUE                                                        
  220 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     BOUNDARY CONDITION FOR VELOCITY FIELD                             
      SUBROUTINE UBC(DIV)                                               
C     1998.01.01                                                        
C====================================================================   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
      INCLUDE './INCFILE/INCDATA'                                         
      INCLUDE './INCFILE/INCMESH'                                         
C                                                                       
C----LOCAL VARIABLES-------------------------------------------         
C                                                                       
      DIMENSION U(-4:IG+4,-2:JG+2,-4:KG+4,3)                            
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /VELO/ UG(-4:IG+4,-2:JG+2,-4:KG+4,3)                       
C                                                                       
C----COMBINATION------------------------------------------------        
C                                                                       
      EQUIVALENCE (UG,U)                                                
C                                                                       
C==============================================================         
C                                                                       
C                                                                       
C--- X-PERIODIC                                                         
!$OMP PARALLEL DO PRIVATE(K)                                            
      DO 10 J=-2,JG+2                                                   
      DO 10 K=-4,KG+4                                                   
        U(  -4,J,K,1)=U(IG-4,J,K,1)                                     
        U(  -3,J,K,1)=U(IG-3,J,K,1)                                     
        U(  -2,J,K,1)=U(IG-2,J,K,1)                                     
        U(  -1,J,K,1)=U(IG-1,J,K,1)                                     
        U(   0,J,K,1)=U(IG  ,J,K,1)                                     
        U(IG+1,J,K,1)=U(   1,J,K,1)                                     
        U(IG+2,J,K,1)=U(   2,J,K,1)                                     
        U(IG+3,J,K,1)=U(   3,J,K,1)                                     
        U(IG+4,J,K,1)=U(   4,J,K,1)                                     
        U(  -4,J,K,2)=U(IG-4,J,K,2)                                     
        U(  -3,J,K,2)=U(IG-3,J,K,2)                                     
        U(  -2,J,K,2)=U(IG-2,J,K,2)                                     
        U(  -1,J,K,2)=U(IG-1,J,K,2)                                     
        U(   0,J,K,2)=U(IG  ,J,K,2)                                     
        U(IG+1,J,K,2)=U(   1,J,K,2)                                     
        U(IG+2,J,K,2)=U(   2,J,K,2)                                     
        U(IG+3,J,K,2)=U(   3,J,K,2)                                     
        U(IG+4,J,K,2)=U(   4,J,K,2)                                     
        U(  -4,J,K,3)=U(IG-4,J,K,3)                                     
        U(  -3,J,K,3)=U(IG-3,J,K,3)                                     
        U(  -2,J,K,3)=U(IG-2,J,K,3)                                     
        U(  -1,J,K,3)=U(IG-1,J,K,3)                                     
        U(   0,J,K,3)=U(IG  ,J,K,3)                                     
        U(IG+1,J,K,3)=U(   1,J,K,3)                                     
        U(IG+2,J,K,3)=U(   2,J,K,3)                                     
        U(IG+3,J,K,3)=U(   3,J,K,3)                                     
        U(IG+4,J,K,3)=U(   4,J,K,3)                                     
   10 CONTINUE                                                          
C                                                                       
C--- Z-PERIODIC                                                         
!$OMP PARALLEL DO PRIVATE(I)                                            
      DO 20 J=-2,JG+2                                                   
      DO 20 I=-4,IG+4                                                   
        U(I,J,  -4,1)=U(I,J,KG-4,1)                                     
        U(I,J,  -3,1)=U(I,J,KG-3,1)                                     
        U(I,J,  -2,1)=U(I,J,KG-2,1)                                     
        U(I,J,  -1,1)=U(I,J,KG-1,1)                                     
        U(I,J,   0,1)=U(I,J,KG  ,1)                                     
        U(I,J,KG+1,1)=U(I,J,   1,1)                                     
        U(I,J,KG+2,1)=U(I,J,   2,1)                                     
        U(I,J,KG+3,1)=U(I,J,   3,1)                                     
        U(I,J,KG+4,1)=U(I,J,   4,1)                                     
        U(I,J,  -4,2)=U(I,J,KG-4,2)                                     
        U(I,J,  -3,2)=U(I,J,KG-3,2)                                     
        U(I,J,  -2,2)=U(I,J,KG-2,2)                                     
        U(I,J,  -1,2)=U(I,J,KG-1,2)                                     
        U(I,J,   0,2)=U(I,J,KG  ,2)                                     
        U(I,J,KG+1,2)=U(I,J,   1,2)                                     
        U(I,J,KG+2,2)=U(I,J,   2,2)                                     
        U(I,J,KG+3,2)=U(I,J,   3,2)                                     
        U(I,J,KG+4,2)=U(I,J,   4,2)                                     
        U(I,J,  -4,3)=U(I,J,KG-4,3)                                     
        U(I,J,  -3,3)=U(I,J,KG-3,3)                                     
        U(I,J,  -2,3)=U(I,J,KG-2,3)                                     
        U(I,J,  -1,3)=U(I,J,KG-1,3)                                     
        U(I,J,   0,3)=U(I,J,KG  ,3)                                     
        U(I,J,KG+1,3)=U(I,J,   1,3)                                     
        U(I,J,KG+2,3)=U(I,J,   2,3)                                     
        U(I,J,KG+3,3)=U(I,J,   3,3)                                     
        U(I,J,KG+4,3)=U(I,J,   4,3)                                     
   20 CONTINUE                                                          
C                                                                       
C--- Y-NONSLIP                                                          
!$OMP PARALLEL DO PRIVATE(I)                                            
      DO 600 K=-4,KG+4                                                  
      DO 600 I=-4,IG+4                                                  
        U(I,  -1,K,1)=-U(I,   2,K,1)                                    
        U(I,  -1,K,2)= U(I,   1,K,2)                                    
        U(I,  -1,K,3)=-U(I,   2,K,3)                                    
        U(I,   0,K,1)=-U(I,   1,K,1)                                    
        U(I,   0,K,2)= 0.0D0                                            
        U(I,   0,K,3)=-U(I,   1,K,3)                                    
        U(I,JG  ,K,2)= 0.0D0                                            
        U(I,JG+1,K,1)=-U(I,JG  ,K,1)                                    
        U(I,JG+1,K,2)= U(I,JG-1,K,2)                                    
        U(I,JG+1,K,3)=-U(I,JG  ,K,3)                                    
        U(I,JG+2,K,1)=-U(I,JG-1,K,1)                                    
        U(I,JG+2,K,3)=-U(I,JG-1,K,3)                                    
  600 CONTINUE                                                          
C                                                                       
C------------------------------------ CONTINUITY CHECKER                
      DIV=0.0D0                                                         
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,QQQQ,QQQ) REDUCTION(MAX:DIV)              
      DO 40 J=1,JG                                                      
      DO 40 K=1,KG                                                      
      DO 40 I=1,IG                                                      
        QQQQ=(9.0D0/8.0D0)*(U(I,J,K,1)-U(I-1,J,K,1))*PDX                
     &                -(1.0D0/8.0D0)*(U(I+1,J,K,1)-U(I-2,J,K,1))        
     &                                /3.0D0*PDX                        
     &       +(U(I,J,K,2)-U(I,J-1,K,2))*PDR(J)                          
     &       +(9.0D0/8.0D0)*(U(I,J,K,3)-U(I,J,K-1,3))*PDZ               
     &                -(1.0D0/8.0D0)*(U(I,J,K+1,3)-U(I,J,K-2,3))        
     &                                /3.0D0*PDZ                        
        QQQ=DABS(QQQQ)                                                  
        DIV=DMAX1(DIV,QQQ)                                              
   40 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     TIME AVERAGE (VELOCITY FIELD)                                     
C     PROGRAMED BY H.ABE                                                
      SUBROUTINE AVE(NN)                                                
C     1998.01.01                                                        
C     REVISED ON 2001.03.04  BY H.ABE (CON_AVE:UV, PDF)                 
C     REVISED ON 2000.10.27  BY H.ABE (CON_AVE:UV)                      
C     REVISED ON 2000.10.16  BY H.ABE (QUADRUPLE CORR.)                 
C     REVISED ON 2000.10.16  BY H.ABE (QUADRANT ANALYSIS: UV)           
C     REVISED BY H.ABE (2000.09.29) => ADD: SKEWUV & FLATUV             
C     REVISED BY H.ABE (2000.09.24) => UV,UW,VW                         
C====================================================================   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
      INCLUDE './INCFILE/INCDATA'                                         
      INCLUDE './INCFILE/INCMESH'                                         
C                                                                       
      CHARACTER NUMM*3                                                  
      COMMON /ANUM/ NAVE                                                
C                                                                       
C----LOCAL VARIABLES-------------------------------------------         
C                                                                       
      DIMENSION U(-4:IG+4,-2:JG+2,-4:KG+4,3)                            
      DIMENSION UT(-4:IG+4,-2:JG+2,-4:KG+4,3)                           
C                                                                       
      DIMENSION P(-3:IG+3,-2:JG+2,-3:KG+3)                              
      DIMENSION AU(-1:JG+2),AV(-1:JG+2),AW(-1:JG+2),AP(-1:JG+2)         
      DIMENSION AAU(JG),AAV(JG),AAW(JG),AAP(JG)                         
     &         ,AURMS(JG),AVRMS(JG),AWRMS(JG),APRMS(JG)                 
     &         ,ATSST(JG),AREST(JG)                                     
      DIMENSION URMS(JG),VRMS(JG),WRMS(JG),PRMS(JG)                     
      DIMENSION TSST(JG),REST(JG)                                       
      DIMENSION U11(0:LY),U12(0:LY),U13(0:LY)                           
     &         ,U22(0:LY),U23(0:LY),U33(0:LY)                           
      DIMENSION AU11(JG),AU12(JG),AU13(JG)                              
     &         ,AU22(JG),AU23(JG),AU33(JG)                              
      DIMENSION U111(0:LY),U112(0:LY),U222(0:LY)                        
      DIMENSION U332(0:LY),U122(0:LY),U133(0:LY)                        
      DIMENSION U333(0:LY),PPP(0:LY)                                    
      DIMENSION SKEWUV(JG),FLATUV(JG)                                   
      DIMENSION SUV(JG,2),FUV(JG,2)                                     
C                                                                       
      DIMENSION U1111(0:LY),U1112(0:LY),U2222(0:LY)                     
      DIMENSION U1332(0:LY),U1222(0:LY),U1122(0:LY)                     
      DIMENSION U1133(0:LY),U2233(0:LY)                                 
      DIMENSION U3333(0:LY),PPPP(0:LY)                                  
      DIMENSION AUTRI(JG,8),AUQUAD(JG,10),AQAUV(JG,4)                   
      DIMENSION NQAUV(JG,4)                                             
      DIMENSION ANQAUV(JG,4)                                            
C                                                                       
      DIMENSION DINV2(0:IG,0:JG,0:KG)                                   
      DIMENSION AINV2(JG),AINV3(JG)                                     
      DIMENSION AAINV2(JG),AAINV3(JG)                                   
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /VELO/ UG(-4:IG+4,-2:JG+2,-4:KG+4,3)                       
      COMMON /TEMP/ UTG(-4:IG+4,-2:JG+2,-4:KG+4,3)                      
      COMMON /SCAL/ PG(-3:IG+3,-2:JG+2,-3:KG+3)                         
      COMMON /MEAN/ AUG(-1:JG+2),AVG(-1:JG+2)                           
     &             ,AWG(-1:JG+2),APG(-1:JG+2)                           
      COMMON /UTAU/ AUTUB,AUTUT                                         
      COMMON /UTAU2/ AUTUB2,AUTUT2                                      
      COMMON /AVE1/ AAUG(JG),AAVG(JG),AAWG(JG),AAPG(JG)                 
      COMMON /AVE2/ AURMSG(JG),AVRMSG(JG),AWRMSG(JG),APRMSG(JG)         
      COMMON /AVE3/ ATSSTG(JG),ARESTG(JG)                               
      COMMON /AVE4/ AU11G(JG),AU12G(JG),AU13G(JG)                       
     &             ,AU22G(JG),AU23G(JG),AU33G(JG)                       
      COMMON /SFUV/ SKEWUVG(JG),FLATUVG(JG)                             
C                                                                       
      COMMON /UHOA/ AUTRIG(JG,8),AUQUADG(JG,10),AQAUVG(JG,4)            
      COMMON /NQAUV/ ANQAUVG(JG,4)                                      
C                                                                       
      COMMON /AINV/ AINV2G(JG),AINV3G(JG)                               
      COMMON /AAIN/ AAINV2G(JG),AAINV3G(JG)                             
C                                                                       
C----COMBINATION------------------------------------------------        
C                                                                       
      EQUIVALENCE (UG,U),(UTG,UT),(PG,P)                                
      EQUIVALENCE (AUG,AU),(AVG,AV),(AWG,AW),(APG,AP)                   
      EQUIVALENCE (AAUG,AAU),(AAVG,AAV),(AAWG,AAW),(AAPG,AAP)           
     &           ,(AURMSG,AURMS),(AVRMSG,AVRMS)                         
     &           ,(AWRMSG,AWRMS),(APRMSG,APRMS)                         
     &           ,(ATSSTG,ATSST),(ARESTG,AREST)                         
      EQUIVALENCE (AU11G,AU11),(AU12G,AU12),(AU13G,AU13)                
     &           ,(AU22G,AU22),(AU23G,AU23),(AU33G,AU33)                
      EQUIVALENCE (SKEWUVG,SKEWUV),(FLATUVG,FLATUV)                     
C                                                                       
      EQUIVALENCE (AUTRIG,AUTRI),(AUQUADG,AUQUAD),(AQAUVG,AQAUV)        
      EQUIVALENCE (ANQAUVG,ANQAUV)                                      
C                                                                       
      EQUIVALENCE (AINV2G,AINV2),(AINV3G,AINV3)                         
      EQUIVALENCE (AAINV2G,AAINV2),(AAINV3G,AAINV3)                     
C                                                                       
C=================================================================      
C                                                                       
C                                                                       
      DO 10 J=1,JG                                                      
        AU(J)=0.0D0                                                     
        AV(J)=0.0D0                                                     
        AW(J)=0.0D0                                                     
        AP(J)=0.0D0                                                     
   10 CONTINUE                                                          
C                                                                       
      DO 11 J=1,JG                                                      
        URMS(J)=0.0D0                                                   
        VRMS(J)=0.0D0                                                   
        WRMS(J)=0.0D0                                                   
        PRMS(J)=0.0D0                                                   
        REST(J)=0.0D0                                                   
        TSST(J)=0.0D0                                                   
        SUV(J,1)=0.0D0                                                  
        SUV(J,2)=0.0D0                                                  
        FUV(J,1)=0.0D0                                                  
        FUV(J,2)=0.0D0                                                  
   11 CONTINUE                                                          
C                                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
        DO 12 K=1,KG                                                    
        DO 12 I=1,IG                                                    
          U(I,0,K,2)=0.0D0                                              
   12   CONTINUE                                                        
!$OMP PARALLEL DO PRIVATE(I)                                            
        DO 13 K=1,KG                                                    
        DO 13 I=1,IG                                                    
          U(I,JG,K,2)=0.0D0                                             
   13   CONTINUE                                                        
C----------------------------------------------------                   
C                                                                       
      DO 14 J=0,JG+1                                                    
        U11(J)=0.0D0                                                    
        U12(J)=0.0D0                                                    
        U13(J)=0.0D0                                                    
        U22(J)=0.0D0                                                    
        U23(J)=0.0D0                                                    
        U33(J)=0.0D0                                                    
        U111(J)=0.0D0                                                   
        U112(J)=0.0D0                                                   
        U222(J)=0.0D0                                                   
        U332(J)=0.0D0                                                   
        U122(J)=0.0D0                                                   
        U133(J)=0.0D0                                                   
        U333(J)=0.0D0                                                   
        PPP(J)=0.0D0                                                    
C                                                                       
        U1111(J)=0.0D0                                                  
        U1112(J)=0.0D0                                                  
        U2222(J)=0.0D0                                                  
        U1332(J)=0.0D0                                                  
        U1222(J)=0.0D0                                                  
        U1122(J)=0.0D0                                                  
        U1133(J)=0.0D0                                                  
        U2233(J)=0.0D0                                                  
        U3333(J)=0.0D0                                                  
        PPPP(J)=0.0D0                                                   
C                                                                       
        AINV2(J)=0.0D0                                                  
        AINV3(J)=0.0D0                                                  
   14 CONTINUE                                                          
C                                                                       
      DO 15 J=1,JG                                                      
        NQAUV(J,1)=0                                                    
        NQAUV(J,2)=0                                                    
        NQAUV(J,3)=0                                                    
        NQAUV(J,4)=0                                                    
   15 CONTINUE                                                          
C                                                                       
      DOOP=1.0D0/DBLE(IG*KG)                                            
C-------------------------------------- PLANE AVE.                      
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 20 J=1,JG                                                      
      DO 20 K=1,KG                                                      
      DO 20 I=1,IG                                                      
        AU(J)=U(I,J,K,1)+AU(J)                                          
        AV(J)=U(I,J,K,2)+AV(J)                                          
        AW(J)=U(I,J,K,3)+AW(J)                                          
        AP(J)=P(I,J,K)+AP(J)                                            
   20 CONTINUE                                                          
C                                                                       
      DO 21 J=1,JG                                                      
        AU(J)=AU(J)*DOOP                                                
        AV(J)=AV(J)*DOOP                                                
        AW(J)=AW(J)*DOOP                                                
        AP(J)=AP(J)*DOOP                                                
   21 CONTINUE                                                          
C                                                                       
!$OMP PARALLEL SECTIONS                                                 
!$OMP SECTION                                                           
      AU(  -1)=-AU(2)                                                   
      AV(  -1)= AV(1)                                                   
      AW(  -1)=-AW(2)                                                   
!$OMP SECTION                                                           
      AU(   0)=-AU(1)                                                   
      AV(   0)= 0.0D0                                                   
      AW(   0)=-AW(1)                                                   
      AP(   0)= AP(1)                                                   
C                                                                       
!$OMP SECTION                                                           
      AV(  JG)= 0.0D0                                                   
      AU(JG+1)=-AU(JG)                                                  
      AV(JG+1)= AV(JG-1)                                                
      AW(JG+1)=-AW(JG)                                                  
      AP(JG+1)= AP(JG)                                                  
!$OMP SECTION                                                           
      AU(JG+2)=-AU(JG-1)                                                
      AW(JG+2)=-AW(JG-1)                                                
!$OMP END PARALLEL SECTIONS                                             
C                                                                       
C----------------------------------------------------                   
      UTAUT=DSQRT(((4.0D0*DRR(JG)**2+4.0D0*DRR(JG)*DRR(JG+1)+           
     &      DRR(JG+1)**2)                                               
     &      *AUG(JG)-DRR(JG+1)**2*AUG(JG-1))                            
     &      /(DRR(JG+1)*(DRR(JG+1)+2.0D0*DRR(JG))*DRR(JG))/RE)          
      UTAUB=DSQRT(((DRR(1)**2+4.0D0*DRR(1)*DRR(2)+4.0D0                 
     &    *DRR(2)**2)*AUG(1)                                            
     &    -DRR(1)**2*AUG(2))/(DRR(1)*(DRR(1)+2.0D0*DRR(2))*DRR(2))/RE)  
      UTAUT2=DSQRT(2.0*AUG(JG)/DRR(JG+1)/RE)                            
      UTAUB2=DSQRT(2.0*AUG(1 )/DRR(   1)/RE)                            
C                                                                       
C========================= UT ZERO-CLEAR                                
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 22 J=-2,JG+2                                                   
      DO 22 K=-4,KG+4                                                   
      DO 22 I=-4,IG+4                                                   
        UT(I,J,K,1)=0.0D0                                               
        UT(I,J,K,2)=0.0D0                                               
        UT(I,J,K,3)=0.0D0                                               
   22 CONTINUE                                                          
C                                                                       
C========================= PRODUCING UV                                 
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 23 J=1,JG                                                      
      DO 23 K=1,KG                                                      
      DO 23 I=1,IG                                                      
        UT(I,J,K,1)=((U(I,J+1,K,1)-AU(J+1))*DR(J)/(DR(J)+DR(J+1))       
     &              +(U(I,J,K,1)-AU(J))*DR(J+1)/(DR(J)+DR(J+1)))        
     &         *0.5*((U(I,J,K,2)-AV(J))+(U(I+1,J,K,2)-AV(J)))           
   23 CONTINUE                                                          
C                                                                       
C-------------------------------------- PLANE AVE.                      
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 30 J=1,JG                                                      
      DO 30 K=1,KG                                                      
      DO 30 I=1,IG                                                      
        URMS(J)=URMS(J)+(U(I,J,K,1)-AU(J))**2                           
        VRMS(J)=VRMS(J)+(U(I,J,K,2)-AV(J))**2                           
        WRMS(J)=WRMS(J)+(U(I,J,K,3)-AW(J))**2                           
        PRMS(J)=PRMS(J)+(P(I,J,K)-AP(J))**2                             
C                                                                       
        REST(J)=REST(J)+UT(I,J,K,1)                                     
C                                                                       
        SUV(J,1)=SUV(J,1)+UT(I,J,K,1)**3                                
        SUV(J,2)=SUV(J,2)+UT(I,J,K,1)**2                                
        FUV(J,1)=FUV(J,1)+UT(I,J,K,1)**4                                
        FUV(J,2)=FUV(J,2)+UT(I,J,K,1)**2                                
C                                                                       
C========================= UV CONVENTIONAL AVE.                         
       IF(UT(I,J,K,1).GE.0) THEN                                        
        GO TO 1000                                                      
       ELSE                                                             
        GO TO 1100                                                      
       ENDIF                                                            
 1000   IF((U(I,J,K,2)-AV(J)).GE.0) THEN                                
         AQAUV(J,1)=AQAUV(J,1)+UT(I,J,K,1)                              
         NQAUV(J,1)=NQAUV(J,1)+1                                        
        ELSE                                                            
         AQAUV(J,3)=AQAUV(J,3)+UT(I,J,K,1)                              
         NQAUV(J,3)=NQAUV(J,3)+1                                        
        ENDIF                                                           
        GO TO 1200                                                      
 1100   IF((U(I,J,K,2)-AV(J)).GE.0) THEN                                
         AQAUV(J,2)=AQAUV(J,2)+UT(I,J,K,1)                              
         NQAUV(J,2)=NQAUV(J,2)+1                                        
        ELSE                                                            
         AQAUV(J,4)=AQAUV(J,4)+UT(I,J,K,1)                              
         NQAUV(J,4)=NQAUV(J,4)+1                                        
        ENDIF                                                           
        GO TO 1200                                                      
 1200  CONTINUE                                                         
C                                                                       
   30 CONTINUE                                                          
C                                                                       
C========================= UT ZERO-CLEAR                                
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 24 J=-2,JG+2                                                   
      DO 24 K=-4,KG+4                                                   
      DO 24 I=-4,IG+4                                                   
        UT(I,J,K,1)=0.0D0                                               
        UT(I,J,K,2)=0.0D0                                               
        UT(I,J,K,3)=0.0D0                                               
   24 CONTINUE                                                          
C                                                                       
C-------------------------------------- U',V',W'                        
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 40 J=-1,JG+2                                                   
      DO 40 K=-1,KG+2                                                   
      DO 40 I=-1,IG+2                                                   
        UT(I,J,K,1)=U(I,J,K,1)-AU(J)                                    
        UT(I,J,K,2)=U(I,J,K,2)-AV(J)                                    
        UT(I,J,K,3)=U(I,J,K,3)-AW(J)                                    
   40 CONTINUE                                                          
C                                                                       
C                                                                       
      DO 31 J=1,JG                                                      
        URMS(J)=URMS(J)*DOOP                                            
        VRMS(J)=VRMS(J)*DOOP                                            
        WRMS(J)=WRMS(J)*DOOP                                            
        PRMS(J)=PRMS(J)*DOOP                                            
        REST(J)=REST(J)*DOOP                                            
        TSST(J)=REST(J)-(AU(J+1)-AU(J))/DRR(J+1)/RE                     
   31 CONTINUE                                                          
C                                                                       
      DO 32 J=1,JG-1                                                    
        SUV(J,1)=SUV(J,1)*DOOP                                          
        SUV(J,2)=DSQRT((SUV(J,2)*DOOP)**3)                              
        FUV(J,1)=FUV(J,1)*DOOP                                          
        FUV(J,2)=(FUV(J,2)*DOOP)**2                                     
        SKEWUV(J)=SKEWUV(J)+SUV(J,1)/SUV(J,2)                           
        FLATUV(J)=FLATUV(J)+FUV(J,1)/FUV(J,2)                           
   32 CONTINUE                                                          
C                                                                       
      DO 33 J=1,JG                                                      
         ANQAUV(J,1)=ANQAUV(J,1)+NQAUV(J,1)*DOOP                        
         ANQAUV(J,2)=ANQAUV(J,2)+NQAUV(J,2)*DOOP                        
         ANQAUV(J,3)=ANQAUV(J,3)+NQAUV(J,3)*DOOP                        
         ANQAUV(J,4)=ANQAUV(J,4)+NQAUV(J,4)*DOOP                        
   33 CONTINUE                                                          
C--------------------------------------                                 
!$OMP PARALLEL DO PRIVATE(I,K,VGT11,VGT12,VGT13,VGT21,VGT22,VGT23       
!$OMP&                   ,VGT31,VGT32,VGT33)                            
      DO 51 J= 1, JG                                                    
      DO 51 K= 1, KG                                                    
      DO 51 I= 1, IG                                                    
        VGT11=(UT(I,J,K,1)-UT(I-1,J,K,1))*PDX                           
        VGT12=0.25D0*((UT(I  ,J+1,K,1)-UT(I  ,J,K,1))*PDRR(J+1)         
     &               +(UT(I  ,J,K,1)-UT(I  ,J-1,K,1))*PDRR(J)           
     &               +(UT(I-1,J+1,K,1)-UT(I-1,J,K,1))*PDRR(J+1)         
     &               +(UT(I-1,J,K,1)-UT(I-1,J-1,K,1))*PDRR(J) )         
        VGT13=0.25D0*PDZ*(UT(I,J,K+1,1)+UT(I-1,J,K+1,1)                 
     &                   -UT(I,J,K-1,1)-UT(I-1,J,K-1,1))                
        VGT21=0.25D0*PDX*(UT(I+1,J,K,2)+UT(I+1,J-1,K,2)                 
     &                   -UT(I-1,J,K,2)-UT(I-1,J-1,K,2))                
        VGT22=(UT(I,J,K,2)-UT(I,J-1,K,2))*PDR(J)                        
        VGT23=0.25D0*PDZ*(UT(I,J,K+1,2)+UT(I,J-1,K+1,2)                 
     &                   -UT(I,J,K-1,2)-UT(I,J-1,K-1,2))                
        VGT31=0.25D0*PDX*(UT(I+1,J,K,3)+UT(I+1,J,K-1,3)                 
     &                   -UT(I-1,J,K,3)-UT(I-1,J,K-1,3))                
        VGT32=0.25D0*((UT(I,J+1,K  ,3)-UT(I,J,K  ,3))*PDRR(J+1)         
     &               +(UT(I,J,K  ,3)-UT(I,J-1,K  ,3))*PDRR(J)           
     &               +(UT(I,J+1,K-1,3)-UT(I,J,K-1,3))*PDRR(J+1)         
     &               +(UT(I,J,K-1,3)-UT(I,J-1,K-1,3))*PDRR(J) )         
        VGT33=(UT(I,J,K,3)-UT(I,J,K-1,3))*PDZ                           
C                                                                       
        AINV2(J)=AINV2(J)+((VGT11**2)+(VGT22**2)+(VGT33**2)             
     &              +2.0D0*(VGT12*VGT21+VGT23*VGT32+VGT31*VGT13))       
        AINV3(J)=AINV3(J)+((VGT11**3)+(VGT22**3)+(VGT33**3)             
     &              +3.0D0*( VGT11*VGT12*VGT21+VGT11*VGT13*VGT31        
     &                      +VGT12*VGT22*VGT21+VGT13*VGT33*VGT31        
     &                      +VGT22*VGT23*VGT32+VGT23*VGT33*VGT32)       
     &              +6.0D0*VGT12*VGT23*VGT31)                           
C                                                                       
        DINV2(I,J,K)=PRE*PRE*((VGT11**2)+(VGT22**2)+(VGT33**2)          
     &              +2.0D0*(VGT12*VGT21+VGT23*VGT32+VGT31*VGT13))       
   51 CONTINUE                                                          
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 41 J=0,JG+1                                                    
      DO 41 K=1,KG                                                      
      DO 41 I=1,IG                                                      
       U11(J)=U11(J)+U(I,J,K,1)*U(I,J,K,1)-AU(J)*AU(J)                  
       U22(J)=U22(J)+U(I,J,K,2)*U(I,J,K,2)-AV(J)*AV(J)                  
       U33(J)=U33(J)+U(I,J,K,3)*U(I,J,K,3)-AW(J)*AW(J)                  
C                                                                       
       U12(J)=U12(J)+(UT(I,J,K,1)+UT(I-1,J,K,1))*0.5D0                  
     &              *(UT(I,J,K,2)+UT(I,J-1,K,2))*0.5D0                  
       U13(J)=U13(J)+(UT(I,J,K,1)+UT(I-1,J,K,1))*0.5D0                  
     &              *(UT(I,J,K,3)+UT(I,J,K-1,3))*0.5D0                  
       U23(J)=U23(J)+(UT(I,J,K,2)+UT(I,J-1,K,2))*0.5D0                  
     &              *(UT(I,J,K,3)+UT(I,J,K-1,3))*0.5D0                  
C============================ TRIPLE CORR.                              
       U111(J)=U111(J)+UT(I,J,K,1)**3                                   
       U112(J)=U112(J)+((UT(I,J,K,1)+UT(I-1,J,K,1))*0.5D0)**2           
     &                *(UT(I,J,K,2)+UT(I,J-1,K,2))*0.5D0                
       U222(J)=U222(J)+UT(I,J,K,2)**3                                   
       U332(J)=U332(J)+((UT(I,J,K,3)+UT(I,J,K-1,3))*0.5D0)**2           
     &                *(UT(I,J,K,2)+UT(I,J-1,K,2))*0.5D0                
       U122(J)=U122(J)+(UT(I,J,K,1)+UT(I-1,J,K,1))*0.5D0                
     &                *((UT(I,J,K,2)+UT(I,J-1,K,2))*0.5D0)**2           
       U133(J)=U133(J)+(UT(I,J,K,1)+UT(I-1,J,K,1))*0.5D0                
     &                *((UT(I,J,K,3)+UT(I,J,K-1,3))*0.5D0)**2           
       U333(J)=U333(J)+((UT(I,J,K,3)+UT(I,J,K-1,3))*0.5D0)**3           
       PPP(J)=PPP(J)+(P(I,J,K)-AP(J))**3                                
C============================ QUADRUPLE CORR.                           
       U1111(J)=U1111(J)+UT(I,J,K,1)**4                                 
       U1112(J)=U1112(J)+((UT(I,J,K,1)+UT(I-1,J,K,1))*0.5D0)**3         
     &                *(UT(I,J,K,2)+UT(I,J-1,K,2))*0.5D0                
       U2222(J)=U2222(J)+UT(I,J,K,2)**4                                 
       U1332(J)=U1332(J)+((UT(I,J,K,1)+UT(I-1,J,K,1))*0.5D0)            
     &                *((UT(I,J,K,3)+UT(I,J,K-1,3))*0.5D0)**2           
     &                *(UT(I,J,K,2)+UT(I,J-1,K,2))*0.5D0                
       U1222(J)=U1222(J)+(UT(I,J,K,1)+UT(I-1,J,K,1))*0.5D0              
     &                *((UT(I,J,K,2)+UT(I,J-1,K,2))*0.5D0)**3           
       U1122(J)=U1122(J)+((UT(I,J,K,1)+UT(I-1,J,K,1))*0.5D0)**2         
     &                *((UT(I,J,K,2)+UT(I,J-1,K,2))*0.5D0)**2           
       U1133(J)=U1133(J)+((UT(I,J,K,1)+UT(I-1,J,K,1))*0.5D0)**2         
     &                *((UT(I,J,K,3)+UT(I,J,K-1,3))*0.5D0)**2           
       U2233(J)=U2233(J)+((UT(I,J,K,2)+UT(I,J-1,K,2))*0.5D0)**2         
     &                *((UT(I,J,K,3)+UT(I,J,K-1,3))*0.5D0)**2           
       U3333(J)=U3333(J)+((UT(I,J,K,3)+UT(I,J,K-1,3))*0.5D0)**4         
       PPPP(J)=PPPP(J)+(P(I,J,K)-AP(J))**4                              
   41 CONTINUE                                                          
C                                                                       
      DV3=1.0D0/3.0D0                                                   
      DO 42 J=0,JG+1                                                    
        AINV2(J)=AINV2(J)*PRE*PRE    *0.5D0                             
        AINV3(J)=AINV3(J)*PRE*PRE*PRE*DV3                               
   42 CONTINUE                                                          
C                                                                       
      DO 43 J=0,JG+1                                                    
        U11(J)=U11(J)*DOOP                                              
        U12(J)=U12(J)*DOOP                                              
        U13(J)=U13(J)*DOOP                                              
        U22(J)=U22(J)*DOOP                                              
        U23(J)=U23(J)*DOOP                                              
        U33(J)=U33(J)*DOOP                                              
        U111(J)=U111(J)*DOOP                                            
        U112(J)=U112(J)*DOOP                                            
        U222(J)=U222(J)*DOOP                                            
        U332(J)=U332(J)*DOOP                                            
        U122(J)=U122(J)*DOOP                                            
        U133(J)=U133(J)*DOOP                                            
        U333(J)=U333(J)*DOOP                                            
        PPP(J)=PPP(J)*DOOP                                              
C                                                                       
        U1111(J)=U1111(J)*DOOP                                          
        U1112(J)=U1112(J)*DOOP                                          
        U2222(J)=U2222(J)*DOOP                                          
        U1332(J)=U1332(J)*DOOP                                          
        U1222(J)=U1222(J)*DOOP                                          
        U1122(J)=U1122(J)*DOOP                                          
        U1133(J)=U1133(J)*DOOP                                          
        U2233(J)=U2233(J)*DOOP                                          
        U3333(J)=U3333(J)*DOOP                                          
        PPPP(J)=PPPP(J)*DOOP                                            
C                                                                       
        AINV2(J)=AINV2(J)*DOOP                                          
        AINV3(J)=AINV3(J)*DOOP                                          
   43 CONTINUE                                                          
C                                                                       
C------------------------------------ TIME AVE.                         
      AUTUB=AUTUB+UTAUB                                                 
      AUTUT=AUTUT+UTAUT                                                 
      AUTUB2=AUTUB2+UTAUB2                                              
      AUTUT2=AUTUT2+UTAUT2                                              
C                                                                       
      DO 80 J=1,JG                                                      
        AAU(J)=AU(J)+AAU(J)                                             
        AAV(J)=AV(J)+AAV(J)                                             
        AAW(J)=AW(J)+AAW(J)                                             
        AAP(J)=AP(J)+AAP(J)                                             
        AURMS(J)=URMS(J)+AURMS(J)                                       
        AVRMS(J)=VRMS(J)+AVRMS(J)                                       
        AWRMS(J)=WRMS(J)+AWRMS(J)                                       
        APRMS(J)=PRMS(J)+APRMS(J)                                       
        ATSST(J)=TSST(J)+ATSST(J)                                       
        AREST(J)=REST(J)+AREST(J)                                       
C                                                                       
        AU11(J)=AU11(J)+U11(J)                                          
        AU12(J)=AU12(J)+U12(J)                                          
        AU13(J)=AU13(J)+U13(J)                                          
        AU22(J)=AU22(J)+U22(J)                                          
        AU23(J)=AU23(J)+U23(J)                                          
        AU33(J)=AU33(J)+U33(J)                                          
        AUTRI(J,1)=AUTRI(J,1)+U111(J)                                   
        AUTRI(J,2)=AUTRI(J,2)+U112(J)                                   
        AUTRI(J,3)=AUTRI(J,3)+U222(J)                                   
        AUTRI(J,4)=AUTRI(J,4)+U332(J)                                   
        AUTRI(J,5)=AUTRI(J,5)+U122(J)                                   
        AUTRI(J,6)=AUTRI(J,6)+U133(J)                                   
        AUTRI(J,7)=AUTRI(J,7)+U333(J)                                   
        AUTRI(J,8)=AUTRI(J,8)+PPP(J)                                    
C                                                                       
        AUQUAD(J,1)=AUQUAD(J,1)+U1111(J)                                
        AUQUAD(J,2)=AUQUAD(J,2)+U1112(J)                                
        AUQUAD(J,3)=AUQUAD(J,3)+U2222(J)                                
        AUQUAD(J,4)=AUQUAD(J,4)+U1332(J)                                
        AUQUAD(J,5)=AUQUAD(J,5)+U1222(J)                                
        AUQUAD(J,6)=AUQUAD(J,6)+U1122(J)                                
        AUQUAD(J,7)=AUQUAD(J,7)+U1133(J)                                
        AUQUAD(J,8)=AUQUAD(J,8)+U2233(J)                                
        AUQUAD(J,9)=AUQUAD(J,9)+U3333(J)                                
        AUQUAD(J,10)=AUQUAD(J,10)+PPPP(J)                               
C                                                                       
        AAINV2(J)=AAINV2(J)+AINV2(J)                                    
        AAINV3(J)=AAINV3(J)+AINV3(J)                                    
   80 CONTINUE                                                          
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      IF(MOD(NN,1000).EQ.0) THEN                                        
       PAVE=1.0/DBLE(NAVE)                                              
       TAUTOP=AUTUT*PAVE                                                
       TAUBOT=AUTUB*PAVE                                                
       AAUGCL=AAUG(JG/2)*PAVE                                           
C                                                                       
       UBULK=0.0D0                                                      
       ABULK=0.0D0                                                      
C                                                                       
      DO 100 J=1,JG                                                     
        UBULK=UBULK+DR(J)* AUG(J)                                       
        ABULK=ABULK+DR(J)*AAUG(J)*PAVE                                  
  100 CONTINUE                                                          
C                                                                       
      OPEN(41,ACCESS='APPEND')                                          
      WRITE(41,'(I9, E13.6, 6F8.4)')                                    
     & NUM, TIME, TAUTOP, TAUBOT, UBULK, ABULK, AUG(JG/2), AAUGCL       
      ENDIF                                                             
      CLOSE(41)                                                         
C                                                                       
C--------------                                                         
C                                                                       
C      NNN=NUM/500                                                      
C      NONE=MOD(NNN,10)                                                 
C      NTEN=MOD(NNN,100)-MOD(NNN,10)                                    
C      NHUN=MOD(NNN,1000)-MOD(NNN,100)                                  
C      NUMM=CHAR(48+NHUN/100)//CHAR(48+NTEN/10)//CHAR(48+NONE)          
CC                                                                      
C      OPEN(91, FILE='/fshome/home01/s21437/POIFLW/CH064W/MOV/ch064w_vel
C     &//NUMM//'.dat'                                                   
C     &, STATUS='UNKNOWN')                                              
C       DO 910 K=1,KG,4                                                 
C       DO 910 J=1,JG/2,2                                               
C       DO 910 I=1,IG,4                                                 
C        WRITE(91,'(2E14.6)') UT(I,J,K,1),DINV2(I,J,K)*1000             
C  910  CONTINUE                                                        
C      CLOSE(91)                                                        
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     BUDGET : CONSISTENT (VELOCITY FIELD)                              
C     PROGRAMED BY H.ABE                                                
      SUBROUTINE BUDGET                                                 
C     1997.01.02                                                        
C====================================================================   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
C                                                                       
C----LOCAL VARIABLES-------------------------------------------         
C                                                                       
      INCLUDE './INCFILE/INCDATA'                                         
      INCLUDE './INCFILE/INCMESH'                                         
      DIMENSION U(-4:IG+4,-2:JG+2,-4:KG+4,3)                            
      DIMENSION UT(-4:IG+4,-2:JG+2,-4:KG+4,3)                           
      DIMENSION P(-3:IG+3,-2:JG+2,-3:KG+3)                              
      DIMENSION AU(-1:JG+2),AV(-1:JG+2),AW(-1:JG+2),AP(-1:JG+2)         
      DIMENSION AK(0:JG)                                                
     & ,PK(0:JG),P11(0:JG),P12(0:JG)                                    
     & ,PIK(0:JG),PI11(0:JG),PI22(0:JG),PI33(0:JG),PI12(0:JG)           
     & ,DK(0:JG),D11(0:JG),D22(0:JG),D33(0:JG),D12(0:JG)                
     & ,EK(0:JG),E11(0:JG),E22(0:JG),E33(0:JG),E12(0:JG)                
     & ,TAK(0:JG),TA11(0:JG),TA22(0:JG),TA33(0:JG),TA12(0:JG)           
     & ,PDK(0:JG),PD22(0:JG),PD12(0:JG)                                 
     & ,PS11(0:JG),PS22(0:JG),PS33(0:JG),PS12(0:JG)                     
      DIMENSION U11(0:JG),U22(0:JG),U33(0:JG)                           
      DIMENSION PI22G(0:JG),D22G(0:JG),E22G(0:JG),TA22G(0:JG)           
C                                                                       
      DIMENSION AAK(JG)                                                 
     &         ,APK(JG),AP11(JG),AP12(JG)                               
     &         ,APIK(JG),API11(JG),API22(JG),API33(JG),API12(JG)        
     &         ,ADK(JG),AD11(JG),AD22(JG),AD33(JG),AD12(JG)             
     &         ,AEK(JG),AE11(JG),AE22(JG),AE33(JG),AE12(JG)             
     &         ,ATAK(JG),ATA11(JG),ATA22(JG),ATA33(JG),ATA12(JG)        
     &         ,APDK(JG),APD22(JG),APD12(JG)                            
     &         ,APS11(JG),APS22(JG),APS33(JG),APS12(JG)                 
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /VELO/ UG(-4:IG+4,-2:JG+2,-4:KG+4,3)                       
      COMMON /TEMP/ UTG(-4:IG+4,-2:JG+2,-4:KG+4,3)                      
      COMMON /SCAL/ PG(-3:IG+3,-2:JG+2,-3:KG+3)                         
      COMMON /MEAN/ AUG(-1:JG+2),AVG(-1:JG+2)                           
     &             ,AWG(-1:JG+2),APG(-1:JG+2)                           
      COMMON /BUDG/ AAKG(JG)                                            
     &         ,APKG(JG),AP11G(JG),AP12G(JG)                            
     &         ,APIKG(JG),API11G(JG),API22G(JG),API33G(JG),API12G(JG)   
     &         ,ADKG(JG),AD11G(JG),AD22G(JG),AD33G(JG),AD12G(JG)        
     &         ,AEKG(JG),AE11G(JG),AE22G(JG),AE33G(JG),AE12G(JG)        
     &         ,ATAKG(JG),ATA11G(JG),ATA22G(JG),ATA33G(JG),ATA12G(JG)   
     &         ,APDKG(JG),APD22G(JG),APD12G(JG)                         
     &         ,APS11G(JG),APS22G(JG),APS33G(JG),APS12G(JG)             
C                                                                       
C----COMBINATION------------------------------------------------        
C                                                                       
      EQUIVALENCE (UG,U),(UTG,UT),(PG,P)                                
      EQUIVALENCE (AUG,AU),(AVG,AV),(AWG,AW),(APG,AP)                   
      EQUIVALENCE (PI22G,PI22),(D22G,D22),(E22G,E22),(TA22G,TA22)       
      EQUIVALENCE (AAKG,AAK)                                            
     &           ,(APKG,APK),(AP11G,AP11),(AP12G,AP12)                  
     &           ,(APIKG,APIK),(API11G,API11),(API22G,API22)            
     &           ,(API33G,API33),(API12G,API12)                         
     &           ,(ADKG,ADK),(AD11G,AD11),(AD22G,AD22)                  
     &           ,(AD33G,AD33),(AD12G,AD12)                             
     &           ,(AEKG,AEK),(AE11G,AE11),(AE22G,AE22)                  
     &           ,(AE33G,AE33),(AE12G,AE12)                             
     &           ,(ATAKG,ATAK),(ATA11G,ATA11),(ATA22G,ATA22)            
     &           ,(ATA33G,ATA33),(ATA12G,ATA12)                         
     &           ,(APDKG,APDK),(APD22G,APD22),(APD12G,APD12)            
     &           ,(APS11G,APS11),(APS22G,APS22)                         
     &           ,(APS33G,APS33),(APS12G,APS12)                         
C                                                                       
C=================================================================      
C                                                                       
      DO 10 J=0,JG                                                      
          AK(J)=0.0D0                                                   
          PK(J)=0.0D0                                                   
         P11(J)=0.0D0                                                   
         P12(J)=0.0D0                                                   
          EK(J)=0.0D0                                                   
         E11(J)=0.0D0                                                   
         E22(J)=0.0D0                                                   
         E33(J)=0.0D0                                                   
         E12(J)=0.0D0                                                   
          DK(J)=0.0D0                                                   
         D11(J)=0.0D0                                                   
         D22(J)=0.0D0                                                   
         D33(J)=0.0D0                                                   
         D12(J)=0.0D0                                                   
         PIK(J)=0.0D0                                                   
        PI11(J)=0.0D0                                                   
        PI22(J)=0.0D0                                                   
        PI33(J)=0.0D0                                                   
        PI12(J)=0.0D0                                                   
         TAK(J)=0.0D0                                                   
        TA11(J)=0.0D0                                                   
        TA22(J)=0.0D0                                                   
        TA33(J)=0.0D0                                                   
        TA12(J)=0.0D0                                                   
         PDK(J)=0.0D0                                                   
        PD22(J)=0.0D0                                                   
        PD12(J)=0.0D0                                                   
        PS11(J)=0.0D0                                                   
        PS22(J)=0.0D0                                                   
        PS33(J)=0.0D0                                                   
        PS12(J)=0.0D0                                                   
         U11(J)=0.0D0                                                   
         U33(J)=0.0D0                                                   
   10 CONTINUE                                                          
C                                                                       
      DO 11 J=0,JG                                                      
         U22(J)=0.0D0                                                   
   11 CONTINUE                                                          
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 109 J=1,JG                                                     
      DO 109 K=1,KG                                                     
      DO 109 I=1,IG                                                     
        U22(J)=U22(J)+UT(I,J,K,2)**2                                    
  109 CONTINUE                                                          
C                                                                       
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,UC2E,UC2W,TA111,TA112,TA113               
!$OMP& ,TA221,TA222,TA223,TA331,TA332,TA333,TA121,TA122,TA123           
!$OMP& ,D111,D112,D113,D221,D222,D223,D331,D332,D333,D12E               
!$OMP& ,D12W,D12T,D12B,D12N,D12S,D12P,D121,D122,D123                    
!$OMP& ,E121,E122,E123)                                                 
      DO 110 J=1,JG                                                     
      DO 110 K=1,KG                                                     
      DO 110 I=1,IG                                                     
C-------------------------------------- K                               
        U11(J)=U11(J)+UT(I,J,K,1)**2                                    
        U33(J)=U33(J)+UT(I,J,K,3)**2                                    
C-------------------------------------- PRODUCTION                      
        P11(J)=P11(J)                                                   
     &        +0.5*UT(I,J,K,1)                                          
     &         *((UT(I+1,J,K,2)+UT(I,J,K,2))                            
     &         *(AU(J+1)-AU(J))*PDRR(J+1)                               
     &         +(UT(I+1,J-1,K,2)+UT(I,J-1,K,2))                         
     &         *(AU(J)-AU(J-1))*PDRR(J))                                
        UC2E=0.25*((UT(I+1,J,K,2)+UT(I,J,K,2))                          
     &             *(AU(J+1)-AU(J))*PDRR(J+1)                           
     &            +(UT(I+1,J-1,K,2)+UT(I,J-1,K,2))                      
     &             *(AU(J)-AU(J-1))*PDRR(J))                            
        UC2W=0.25*((UT(I,J,K,2)+UT(I-1,J,K,2))                          
     &             *(AU(J+1)-AU(J))*PDRR(J+1)                           
     &            +(UT(I,J-1,K,2)+UT(I-1,J-1,K,2))                      
     &             *(AU(J)-AU(J-1))*PDRR(J))                            
        P12(J)=P12(J)+0.25*(UC2E+UC2W)*(UT(I,J,K,2)+UT(I,J-1,K,2))      
C------------------------------------- TURBULENT DIFFUSION (ANALYTICAL) 
        TA111=((0.5*(UT(I+1,J,K,1)+UT(I,J,K,1)))**3                     
     &        -(0.5*(UT(I-1,J,K,1)+UT(I,J,K,1)))**3)*PDX                
        TA112=((((DR(J)*UT(I,J+1,K,1)+DR(J+1)*UT(I,J,K,1))*PDR2(J))**2) 
     &        *(0.5*(UT(I,J,K,2)+UT(I+1,J,K,2)))                        
     &       -(((DR(J-1)*UT(I,J,K,1)+DR(J)*UT(I,J-1,K,1))*PDR2(J-1))**2)
     &        *(0.5*(UT(I,J-1,K,2)+UT(I+1,J-1,K,2))))*PDR(J)            
        TA113=(((0.5*(UT(I,J,K+1,1)+UT(I,J,K,1)))**2)                   
     &        *(0.5*(UT(I,J,K,3)+UT(I+1,J,K,3)))                        
     &       -((0.5*(UT(I,J,K,1)+UT(I,J,K-1,1)))**2)                    
     &        *(0.5*(UT(I,J,K-1,3)+UT(I+1,J,K-1,3))))*PDZ               
        TA221=(((0.5*(UT(I+1,J,K,2)+UT(I,J,K,2)))**2)                   
     &        *((DR(J)*UT(I,J+1,K,1)+DR(J+1)*UT(I,J,K,1))*PDR2(J))      
     &       -((0.5*(UT(I,J,K,2)+UT(I-1,J,K,2)))**2)                    
     &        *((DR(J)*UT(I-1,J+1,K,1)+DR(J+1)*UT(I-1,J,K,1))*PDR2(J))  
     &        )*PDX                                                     
        TA222=(((0.5*(UT(I,J+1,K,2)+UT(I,J,K,2)))**3)                   
     &        -((0.5*(UT(I,J-1,K,2)+UT(I,J,K,2)))**3))*PDRR(J+1)        
        TA223=(((0.5*(UT(I,J,K+1,2)+UT(I,J,K,2)))**2)                   
     &        *((DR(J)*UT(I,J+1,K,3)+DR(J+1)*UT(I,J,K,3))*PDR2(J))      
     &       -((0.5*(UT(I,J,K,2)+UT(I,J,K-1,2)))**2)                    
     &        *((DR(J)*UT(I,J+1,K-1,3)+DR(J+1)*UT(I,J,K-1,3))*PDR2(J))  
     &        )*PDZ                                                     
        TA331=(((0.5*(UT(I+1,J,K,3)+UT(I,J,K,3)))**2)                   
     &        *(0.5*(UT(I,J,K,1)+UT(I,J,K+1,1)))                        
     &       -((0.5*(UT(I-1,J,K,3)+UT(I,J,K,3)))**2)                    
     &        *(0.5*(UT(I-1,J,K,1)+UT(I-1,J,K+1,1))))*PDX               
        TA332=((((DR(J)*UT(I,J+1,K,3)+DR(J+1)*UT(I,J,K,3))*PDR2(J))**2) 
     &        *(0.5*(UT(I,J,K+1,2)+UT(I,J,K,2)))                        
     &       -(((DR(J-1)*UT(I,J,K,3)+DR(J)*UT(I,J-1,K,3))*PDR2(J-1))**2)
     &        *(0.5*(UT(I,J-1,K+1,2)+UT(I,J-1,K,2))))*PDR(J)            
        TA333=((0.5*(UT(I,J,K+1,3)+UT(I,J,K,3)))**3                     
     &        -(0.5*(UT(I,J,K-1,3)+UT(I,J,K,3)))**3)*PDZ                
        TA121=0.25*PDX                                                  
     &      *((UT(I,J,K,1)**2)                                          
     &       *(UT(I+1,J,K,2)+UT(I,J,K,2)+UT(I+1,J-1,K,2)+UT(I,J-1,K,2)) 
     &      -(UT(I-1,J,K,1)**2)                                         
     &       *(UT(I,J,K,2)+UT(I-1,J,K,2)+UT(I,J-1,K,2)+UT(I-1,J-1,K,2)))
        TA122=(0.5*PDR2(J)*(DR(J)*(UT(I,J+1,K,1)+UT(I-1,J+1,K,1))       
     &                  +DR(J+1)*(UT(I,J,K,1)+UT(I-1,J,K,1)))           
     &                 *(UT(I,J,K,2)**2)                                
     &        -0.5*PDR2(J-1)*(DR(J-1)*(UT(I,J,K,1)+UT(I-1,J,K,1))       
     &                  +DR(J)*(UT(I,J-1,K,1)+UT(I-1,J-1,K,1)))         
     &                 *(UT(I,J-1,K,2)**2))*PDR(J)                      
        TA123=0.0625*PDZ                                                
     &       *(UT(I,J,K,3)                                              
     &         *(UT(I,J,K+1,1)+UT(I-1,J,K+1,1)                          
     &          +UT(I,J,K,1)+UT(I-1,J,K,1))                             
     &         *(UT(I,J,K+1,2)+UT(I,J-1,K+1,2)                          
     &          +UT(I,J,K,2)+UT(I,J-1,K,2))                             
     &        -UT(I,J,K-1,3)                                            
     &         *(UT(I,J,K,1)+UT(I-1,J,K,1)                              
     &          +UT(I,J,K-1,1)+UT(I-1,J,K-1,1))                         
     &         *(UT(I,J,K,2)+UT(I,J-1,K,2)                              
     &          +UT(I,J,K-1,2)+UT(I,J-1,K-1,2)))                        
        TA11(J)=TA11(J)+TA111+TA112+TA113                               
        TA22(J)=TA22(J)+TA221+TA222+TA223                               
        TA33(J)=TA33(J)+TA331+TA332+TA333                               
        TA12(J)=TA12(J)+TA121+TA122+TA123                               
C-------------------------------------- VELOCITY-PRESSURE GRAD.         
C      p => p' by H.ABE (1998.01.10)                                    
C                                                                       
        PI11(J)=PI11(J)+UT(I,J,K,1)*((P(I+1,J,K)-AP(J))                 
     &                              -(P(I,J,K)-AP(J)))*PDX              
        PI22(J)=PI22(J)+UT(I,J,K,2)*((P(I,J+1,K)-AP(J+1))               
     &                              -(P(I,J,K)-AP(J)))*PDRR(J+1)        
        PI33(J)=PI33(J)+UT(I,J,K,3)*(P(I,J,K+1)-P(I,J,K))*PDZ           
        PI12(J)=PI12(J)                                                 
     &         +0.25*(UT(I,J,K,1)+UT(I-1,J,K,1))                        
     &         *(((P(I,J+1,K)-AP(J+1))-(P(I,J,K)-AP(J)))*PDRR(J+1)      
     &          +((P(I,J,K)-AP(J))-(P(I,J-1,K)-AP(J-1)))*PDRR(J))       
     &         +0.25*(UT(I,J,K,2)+UT(I,J-1,K,2))                        
     &         *(((P(I+1,J,K)-AP(J))-(P(I,J,K)-AP(J)))*PDX              
     &          +((P(I,J,K)-AP(J))-(P(I-1,J,K)-AP(J)))*PDX)             
C                                                                       
C-------------------------------------- VISCOUS DIFFUSION               
        D111=PDX2*((UT(I+1,J,K,1)**2)-2.0*(UT(I,J,K,1)**2)              
     &            +(UT(I-1,J,K,1)**2))                                  
        D112=2.0*PDRR3(J)*(DRR(J+1)*(UT(I,J-1,K,1)**2)                  
     &                     -DRR2(J)*(UT(I,J,K,1)**2)                    
     &                      +DRR(J)*(UT(I,J+1,K,1)**2))                 
        D113=PDZ2*((UT(I,J,K+1,1)**2)-2.0*(UT(I,J,K,1)**2)              
     &            +(UT(I,J,K-1,1)**2))                                  
        D221=PDX2*((UT(I+1,J,K,2)**2)-2.0*(UT(I,J,K,2)**2)              
     &            +(UT(I-1,J,K,2)**2))                                  
        D222=2.0*PDR3(J)*(DR(J+1)*(UT(I,J-1,K,2)**2)                    
     &                    -DR2(J)*(UT(I,J,K,2)**2)                      
     &                     +DR(J)*(UT(I,J+1,K,2)**2))                   
        D223=PDZ2*((UT(I,J,K+1,2)**2)-2.0*(UT(I,J,K,2)**2)              
     &            +(UT(I,J,K-1,2)**2))                                  
        D331=PDX2*((UT(I+1,J,K,3)**2)-2.0*(UT(I,J,K,3)**2)              
     &            +(UT(I-1,J,K,3)**2))                                  
        D332=2.0*PDRR3(J)*(DRR(J+1)*(UT(I,J-1,K,3)**2)                  
     &                     -DRR2(J)*(UT(I,J,K,3)**2)                    
     &                      +DRR(J)*(UT(I,J+1,K,3)**2))                 
        D333=PDZ2*((UT(I,J,K+1,3)**2)-2.0*(UT(I,J,K,3)**2)              
     &            +(UT(I,J,K-1,3)**2))                                  
        D11(J)=D11(J)+D111+D112+D113                                    
        D22(J)=D22(J)+D221+D222+D223                                    
        D33(J)=D33(J)+D331+D332+D333                                    
        D12E=0.25*(UT(I+1,J,K,1)+UT(I,J,K,1))                           
     &           *(UT(I+1,J,K,2)+UT(I+1,J-1,K,2))                       
        D12W=0.25*(UT(I-1,J,K,1)+UT(I-2,J,K,1))                         
     &           *(UT(I-1,J,K,2)+UT(I-1,J-1,K,2))                       
        D12T=0.25*(UT(I,J,K+1,1)+UT(I-1,J,K+1,1))                       
     &           *(UT(I,J,K+1,2)+UT(I,J-1,K+1,2))                       
        D12B=0.25*(UT(I,J,K-1,1)+UT(I-1,J,K-1,1))                       
     &           *(UT(I,J,K-1,2)+UT(I,J-1,K-1,2))                       
        D12N=0.25*(UT(I,J+1,K,1)+UT(I-1,J+1,K,1))                       
     &           *(UT(I,J+1,K,2)+UT(I,J,K,2))                           
        D12S=0.25*(UT(I,J-1,K,1)+UT(I-1,J-1,K,1))                       
     &           *(UT(I,J-1,K,2)+UT(I,J-2,K,2))                         
        D12P=0.25*(UT(I,J,K,1)+UT(I-1,J,K,1))                           
     &           *(UT(I,J,K,2)+UT(I,J-1,K,2))                           
        D121=PDX2*(D12E-2.0*D12P+D12W)                                  
        D122=2.0*PDRR3(J)*(DRR(J+1)*D12S-DRR2(J)*D12P+DRR(J)*D12N)      
        D123=PDZ2*(D12T-2.0*D12P+D12B)                                  
        D12(J)=D12(J)+D121+D122+D123                                    
C-------------------------------------- DISSIPATION                     
        E11(J)=E11(J)                                                   
     &        +((UT(I+1,J,K,1)-UT(I,J,K,1))*PDX)**2                     
     &        +((UT(I,J,K,1)-UT(I-1,J,K,1))*PDX)**2                     
     &        +((UT(I,J+1,K,1)-UT(I,J,K,1))*PDRR(J+1))**2               
     &        +((UT(I,J,K,1)-UT(I,J-1,K,1))*PDRR(J))**2                 
     &        +((UT(I,J,K+1,1)-UT(I,J,K,1))*PDZ)**2                     
     &        +((UT(I,J,K,1)-UT(I,J,K-1,1))*PDZ)**2                     
        E22(J)=E22(J)                                                   
     &        +((UT(I+1,J,K,2)-UT(I,J,K,2))*PDX)**2                     
     &        +((UT(I,J,K,2)-UT(I-1,J,K,2))*PDX)**2                     
     &        +((UT(I,J,K+1,2)-UT(I,J,K,2))*PDZ)**2                     
     &        +((UT(I,J,K,2)-UT(I,J,K-1,2))*PDZ)**2                     
     &        +2.0*(DR(J)*PDR2(J)                                       
     &              *((UT(I,J+1,K,2)-UT(I,J,K,2))*PDR(J+1))**2          
     &             +DR(J+1)*PDR2(J)                                     
     &              *((UT(I,J,K,2)-UT(I,J-1,K,2))*PDR(J))**2)           
        E33(J)=E33(J)                                                   
     &        +((UT(I+1,J,K,3)-UT(I,J,K,3))*PDX)**2                     
     &        +((UT(I,J,K,3)-UT(I-1,J,K,3))*PDX)**2                     
     &        +((UT(I,J+1,K,3)-UT(I,J,K,3))*PDRR(J+1))**2               
     &        +((UT(I,J,K,3)-UT(I,J-1,K,3))*PDRR(J))**2                 
     &        +((UT(I,J,K+1,3)-UT(I,J,K,3))*PDZ)**2                     
     &        +((UT(I,J,K,3)-UT(I,J,K-1,3))*PDZ)**2                     
        E121=0.25*PDX2*((UT(I+1,J,K,1)-UT(I-1,J,K,1))                   
     &                 *(UT(I+1,J,K,2)+UT(I+1,J-1,K,2)                  
     &                  -UT(I,J,K,2)-UT(I,J-1,K,2))                     
     &                 +(UT(I,J,K,1)-UT(I-2,J,K,1))                     
     &                 *(UT(I,J,K,2)+UT(I,J-1,K,2)                      
     &                  -UT(I-1,J,K,2)-UT(I-1,J-1,K,2)))                
        E122=(UT(I,J+1,K,1)+UT(I-1,J+1,K,1)                             
     &       -UT(I,J,K,1)-UT(I-1,J,K,1))*PDRR(J+1)*0.5                  
     &      *(UT(I,J+1,K,2)-UT(I,J-1,K,2))*PDR2(J)                      
     &      +(UT(I,J,K,1)+UT(I-1,J,K,1)                                 
     &       -UT(I,J-1,K,1)-UT(I-1,J-1,K,1))*PDRR(J)*0.5                
     &      *(UT(I,J,K,2)-UT(I,J-2,K,2))*PDR2(J-1)                      
        E123=0.25*PDZ2*((UT(I,J,K+1,1)+UT(I-1,J,K+1,1)                  
     &                  -UT(I,J,K,1)-UT(I-1,J,K,1))                     
     &                 *(UT(I,J,K+1,2)+UT(I,J-1,K+1,2)                  
     &                  -UT(I,J,K,2)-UT(I,J-1,K,2))                     
     &                 +(UT(I,J,K,1)+UT(I-1,J,K,1)                      
     &                  -UT(I,J,K-1,1)-UT(I-1,J,K-1,1))                 
     &                 *(UT(I,J,K,2)+UT(I,J-1,K,2)                      
     &                  -UT(I,J,K-1,2)-UT(I,J-1,K-1,2)))                
        E12(J)=E12(J)+E121+E122+E123                                    
C                                                                       
C-------------------------------------- PRESSURE DIFFUSION              
C                                                                       
C                                   (1998.07.27, K. HONMA)              
        PD22(J)=PD22(J)                                                 
     &         +0.5D0*(                                                 
     &               (UT(I,J+1,K,2)+UT(I,J,K,2))*(P(I,J+1,K)-AP(J+1))   
     &              -(UT(I,J,K,2)+UT(I,J-1,K,2))*(P(I,J,K)-AP(J))       
     &                                                    )*PDRR(J+1)   
        PD12(J)=PD12(J)                                                 
     &   +(                                                             
     &       0.5D0*( DR(J)*(UT(I,J+1,K,1)+UT(I-1,J+1,K,1))              
     &            +DR(J+1)*(UT(I,J,K,1)+UT(I-1,J,K,1))                  
     &              )*PDR2(J)                                           
     &       *( DR(J)*(P(I,J+1,K)-AP(J+1))+DR(J+1)*(P(I,J,K)-AP(J))     
     &           )*PDR2(J)                                              
     &    -                                                             
     &       0.5D0*( DR(J-1)*(UT(I,J,K,1)+UT(I-1,J,K,1))                
     &            +DR(J)*(UT(I,J-1,K,1)+UT(I-1,J-1,K,1))                
     &              )*PDR2(J-1)                                         
     &       *( DR(J-1)*(P(I,J,K)-AP(J))+DR(J)*(P(I,J-1,K)-AP(J-1))     
     &           )*PDR2(J-1)                                            
     &      )*PDR(J)                                                    
     &   +(                                                             
     &    0.25D0*(UT(I,J,K,2)+UT(I+1,J,K,2)                             
     &           +UT(I,J-1,K,2)+UT(I+1,J-1,K,2)                         
     &           )                                                      
     &    *0.5D0*( (P(I+1,J,K)-AP(J))+(P(I,J,K)-AP(J))                  
     &           )                                                      
     &   -0.25D0*(UT(I,J,K,2)+UT(I-1,J,K,2)                             
     &           +UT(I,J-1,K,2)+UT(I-1,J-1,K,2)                         
     &           )                                                      
     &    *0.5D0*( (P(I,J,K)-AP(J))+(P(I-1,J,K)-AP(J))                  
     &           )                                                      
     &      )*PDX                                                       
C                                                                       
C----------------------------------------- PRESSURE STRAIN              
C                                                                       
C                                   (1998.07.27, K. HONMA)              
       PS11(J)=PS11(J)                                                  
     &     +0.5D0*(                                                     
     &             (P(I+1,J,K)-AP(J))*(UT(I+1,J,K,1)-UT(I,J,K,1))*PDX   
     &            +(P(I,J,K)-AP(J))*(UT(I,J,K,1)-UT(I-1,J,K,1))*PDX     
     &                                                              )   
       PS22(J)=PS22(J)                                                  
     &   +0.5D0*(                                                       
     &         (P(I,J+1,K)-AP(J+1))*(UT(I,J+1,K,2)-UT(I,J,K,2))*PDR(J+1)
     &        +(P(I,J,K)-AP(J))*(UT(I,J,K,2)-UT(I,J-1,K,2))*PDR(J)      
     &                                                                ) 
       PS33(J)=PS33(J)                                                  
     &        +0.5D0*(                                                  
     &              (P(I,J,K+1)-AP(J))*(UT(I,J,K+1,3)-UT(I,J,K,3))*PDZ  
     &             +(P(I,J,K)-AP(J))*(UT(I,J,K,3)-UT(I,J,K-1,3))*PDZ    
     &                                                                ) 
       PS12(J)=PS12(J)                                                  
     &   +0.5D0*(  PDR2(J)*                                             
     &         ( DR(J)*(P(I,J+1,K)-AP(J+1))+DR(J+1)*(P(I,J,K)-AP(J)) )  
     &         *0.5D0*( (UT(I,J+1,K,1)+UT(I-1,J+1,K,1))                 
     &               -(UT(I,J,K,1)+UT(I-1,J,K,1)) )*PDRR(J+1)           
     &          +PDR2(J-1)*                                             
     &         ( DR(J-1)*(P(I,J,K)-AP(J))+DR(J)*(P(I,J-1,K)-AP(J-1)) )  
     &         *0.5D0*( (UT(I,J,K,1)+UT(I-1,J,K,1))                     
     &               -(UT(I,J-1,K,1)+UT(I-1,J-1,K,1)) )*PDRR(J)        )
     &   +0.5D0*( 0.5D0*( (P(I+1,J,K)-AP(J))+(P(I,J,K)-AP(J)) )         
     &         *0.5D0*( (UT(I+1,J-1,K,2)+UT(I+1,J,K,2))                 
     &               -(UT(I,J-1,K,2)+UT(I,J,K,2)) )*PDX                 
     &         +0.5D0*( (P(I,J,K)-AP(J))+(P(I-1,J,K)-AP(J)) )           
     &         *0.5D0*( (UT(I,J,K,2)+UT(I,J-1,K,2))                     
     &               -(UT(I-1,J,K,2)+UT(I-1,J-1,K,2)) )*PDX           ) 
  110 CONTINUE                                                          
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 120 K=1,KG                                                     
      DO 120 I=1,IG                                                     
        D22(0)=D22(0)+2.0*(UT(I,1,K,2)*PDR(1))**2                       
  120 CONTINUE                                                          
      E22(0)=D22(0)                                                     
C                                                                       
C                                                                       
      DOOP=1.0D0/DBLE(IG*KG)                                            
C                                                                       
      DO 200 J=0,JG                                                     
        P11(J)=-1.0*P11(J)*DOOP                                         
        P12(J)=-1.0*P12(J)*DOOP                                         
        TA11(J)=-1.0*TA11(J)*DOOP                                       
        TA22(J)=-1.0*TA22(J)*DOOP                                       
        TA33(J)=-1.0*TA33(J)*DOOP                                       
        TA12(J)=-1.0*TA12(J)*DOOP                                       
        PI11(J)=-2.0*PI11(J)*DOOP                                       
        PI22(J)=-2.0*PI22(J)*DOOP                                       
        PI33(J)=-2.0*PI33(J)*DOOP                                       
        PI12(J)=-1.0*PI12(J)*DOOP                                       
        D11(J)=D11(J)*DOOP                                              
        D22(J)=D22(J)*DOOP                                              
        D33(J)=D33(J)*DOOP                                              
        D12(J)=D12(J)*DOOP                                              
        E11(J)=E11(J)*DOOP                                              
        E22(J)=E22(J)*DOOP                                              
        E33(J)=E33(J)*DOOP                                              
        E12(J)=E12(J)*DOOP                                              
        PDK(J)=PDK(J)*DOOP                                              
        PD22(J)=-2.0D0*PD22(J)*DOOP                                     
        PD12(J)=-1.0D0*PD12(J)*DOOP                                     
        PS11(J)=2.0D0*PS11(J)*DOOP                                      
        PS22(J)=2.0D0*PS22(J)*DOOP                                      
        PS33(J)=2.0D0*PS33(J)*DOOP                                      
        PS12(J)=1.0D0*PS12(J)*DOOP                                      
  200 CONTINUE                                                          
C                                                                       
      DO 210 J=1,JG                                                     
        AK(J)=0.5*(U11(J)+0.5*(U22(J)+U22(J-1))+U33(J))*DOOP            
        PK(J)=P11(J)*0.5                                                
        TAK(J)=0.5*(TA11(J)+0.5*(TA22(J)+TA22(J-1))+TA33(J))            
        PIK(J)=0.5*(PI11(J)+0.5*(PI22(J)+PI22(J-1))+PI33(J))            
        DK(J)=0.5*(D11(J)+0.5*(D22(J)+D22(J-1))+D33(J))                 
        EK(J)=0.5*(E11(J)+0.5*(E22(J)+E22(J-1))+E33(J))                 
        PDK(J)=0.5D0*(PD22(J)+PD22(J-1))                                
  210 CONTINUE                                                          
C                                                                       
C------------------------------------------------------ TIME AVE.       
      DO 300 J=1,JG                                                     
        AAK(J)=AAK(J)+AK(J)                                             
        APK(J)=APK(J)+PK(J)                                             
        AP11(J)=AP11(J)+P11(J)                                          
        AP12(J)=AP12(J)+P12(J)                                          
        AEK(J)=AEK(J)+EK(J)                                             
        AE11(J)=AE11(J)+E11(J)                                          
        AE22(J)=AE22(J)+E22(J)                                          
        AE33(J)=AE33(J)+E33(J)                                          
        AE12(J)=AE12(J)+E12(J)                                          
        ADK(J)=ADK(J)+DK(J)                                             
        AD11(J)=AD11(J)+D11(J)                                          
        AD22(J)=AD22(J)+D22(J)                                          
        AD33(J)=AD33(J)+D33(J)                                          
        AD12(J)=AD12(J)+D12(J)                                          
        APIK(J)=APIK(J)+PIK(J)                                          
        API11(J)=API11(J)+PI11(J)                                       
        API22(J)=API22(J)+PI22(J)                                       
        API33(J)=API33(J)+PI33(J)                                       
        API12(J)=API12(J)+PI12(J)                                       
        ATAK(J)=ATAK(J)+TAK(J)                                          
        ATA11(J)=ATA11(J)+TA11(J)                                       
        ATA22(J)=ATA22(J)+TA22(J)                                       
        ATA33(J)=ATA33(J)+TA33(J)                                       
        ATA12(J)=ATA12(J)+TA12(J)                                       
        APDK(J)=APDK(J)+PDK(J)                                          
        APD22(J)=APD22(J)+PD22(J)                                       
        APD12(J)=APD12(J)+PD12(J)                                       
        APS11(J)=APS11(J)+PS11(J)                                       
        APS22(J)=APS22(J)+PS22(J)                                       
        APS33(J)=APS33(J)+PS33(J)                                       
        APS12(J)=APS12(J)+PS12(J)                                       
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     ENERGY SPECTRA (VELOCITY FIELD)                                   
C     PROGRAMED BY H.ABE                                                
      SUBROUTINE POWSPE                                                 
C     REVISED ON 2004.01.12 BY T.TSUKAHARA (ARRANGED FOR SX-7)          
C     REVISED ON 2003.11.08 BY Y.SEKI (SPECTRA OF UV)                   
C     REVISED ON 2001.03.03 BY H.ABE (JN=1)                             
C     REVISED ON 2001.01.08 BY H.ABE (ADD: P, UV)                       
C     2000.08.26 REVISED BY H.ABE (SPREAD REGION -> SPREAD DO)          
C====================================================================   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
C                                                                       
      COMMON /JAVE/ JA,JB,JC,JD,JE,JF,JH,JI,JJ,JK,JL                    
C                                                                       
C----LOCAL VARIABLES-------------------------------------------         
C                                                                       
      INCLUDE './INCFILE/INCMESH'                                         
      DIMENSION U(-4:IG+4,-2:JG+2,-4:KG+4,3)                            
      DIMENSION UT(-4:IG+4,-2:JG+2,-4:KG+4,3)                           
      DIMENSION P(-3:IG+3,-2:JG+2,-3:KG+3)                              
      DIMENSION AU(-1:JG+2),AV(-1:JG+2),AW(-1:JG+2),AP(-1:JG+2)         
      DIMENSION ASTRM(IG,JG,5),ASPAN(KG,JG,5)                           
      DIMENSION STRM(IG,JG,5),SPAN(KG,JG,5)                             
      DIMENSION CR1(LX,LZ),CI1(LX,LZ)                                   
     &         ,CR2(LX,LZ),CI2(LX,LZ)                                   
     &         ,CR3(LX,LZ),CI3(LX,LZ)                                   
     &         ,CR4(LX,LZ),CI4(LX,LZ)                                   
      DIMENSION TRIGS(2*(IG+KG)),WK(2*LX*LZ),IFAX(40)                   
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /VELO/ UG(-4:IG+4,-2:JG+2,-4:KG+4,3)                       
      COMMON /TEMP/ UTG(-4:IG+4,-2:JG+2,-4:KG+4,3)                      
      COMMON /SCAL/ PG(-3:IG+3,-2:JG+2,-3:KG+3)                         
      COMMON /MEAN/ AUG(-1:JG+2),AVG(-1:JG+2)                           
     &             ,AWG(-1:JG+2),APG(-1:JG+2)                           
      COMMON /POWR/ ASTRMG(IG,JG,5),ASPANG(KG,JG,5)                     
C                                                                       
C----COMBINATION------------------------------------------------        
C                                                                       
      EQUIVALENCE (UG,U),(UTG,UT),(PG,P)                                
      EQUIVALENCE (AUG,AU),(AVG,AV),(AWG,AW),(APG,AP)                   
      EQUIVALENCE (ASTRMG,ASTRM),(ASPANG,ASPAN)                         
C                                                                       
C=================================================================      
C                                                                       
C                                                                       
      PW=1.0/DBLE(IG*KG)                                                
C                                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
      DO 10 J=1,JG                                                      
      DO 10 I=1,IG                                                      
        STRM(I,J,1)=0.0D0                                               
        STRM(I,J,2)=0.0D0                                               
        STRM(I,J,3)=0.0D0                                               
        STRM(I,J,4)=0.0D0                                               
        STRM(I,J,5)=0.0D0                                               
   10 CONTINUE                                                          
C                                                                       
!$OMP PARALLEL DO PRIVATE(K)                                            
      DO 11 J=1,JG                                                      
      DO 11 K=1,KG                                                      
        SPAN(K,J,1)=0.0D0                                               
        SPAN(K,J,2)=0.0D0                                               
        SPAN(K,J,3)=0.0D0                                               
        SPAN(K,J,4)=0.0D0                                               
        SPAN(K,J,5)=0.0D0                                               
   11 CONTINUE                                                          
C-------------------------------------- U',V',W'                        
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 12 J=-1,JG+2                                                   
      DO 12 K=-1,KG+2                                                   
      DO 12 I=-1,IG+2                                                   
        UT(I,J,K,1)=U(I,J,K,1)-AU(J)                                    
        UT(I,J,K,2)=(U(I,J,K,2)-AV(J)                                   
     &              +U(I,J-1,K,2)-AV(J-1))*0.5D0                        
        UT(I,J,K,3)=U(I,J,K,3)-AW(J)                                    
   12 CONTINUE                                                          
C                                                                       
C--- Y:JA                                                               
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 500 J=JA,JA                                                    
       IF (J.EQ.JA) THEN                                                
        DO 20 K=1,KG                                                    
        DO 20 I=1,IG                                                    
          CR1(I,K)=UT(I,JA,K,1)                                         
          CR2(I,K)=UT(I,JA,K,2)                                         
          CR3(I,K)=UT(I,JA,K,3)                                         
          CR4(I,K)=P(I,JA,K)-AP(JA)                                     
          CI1(I,K)=0.0D0                                                
          CI2(I,K)=0.0D0                                                
          CI3(I,K)=0.0D0                                                
          CI4(I,K)=0.0D0                                                
   20   CONTINUE                                                        
       ENDIF                                                            
  500 CONTINUE                                                          
C                                                                       
      ISW=1                                                             
C      CALL DFC2FB(IG,KG,CR1,CI1,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR2,CI2,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR3,CI3,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR4,CI4,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,WA1,WA2,WA3,WA4,WA5)                      
      DO 501 J=JA,JA                                                    
       IF (J.EQ.JA) THEN                                                
        DO 21 K=1,KG                                                    
        DO 21 I=1,IG                                                    
         WA1=(CR1(I,K)*PW)**2+(CI1(I,K)*PW)**2                          
         WA2=(CR2(I,K)*PW)**2+(CI2(I,K)*PW)**2                          
         WA3=(CR3(I,K)*PW)**2+(CI3(I,K)*PW)**2                          
         WA4=(CR4(I,K)*PW)**2+(CI4(I,K)*PW)**2                          
         WA5=CR1(I,K)*PW*CR2(I,K)*PW+CI1(I,K)*PW*CI2(I,K)*PW            
C                                                                       
         STRM(I,JA,1)=STRM(I,JA,1)+WA1                                  
         STRM(I,JA,2)=STRM(I,JA,2)+WA2                                  
         STRM(I,JA,3)=STRM(I,JA,3)+WA3                                  
         STRM(I,JA,4)=STRM(I,JA,4)+WA4                                  
         STRM(I,JA,5)=STRM(I,JA,5)+WA5                                  
C                                                                       
         SPAN(K,JA,1)=SPAN(K,JA,1)+WA1                                  
         SPAN(K,JA,2)=SPAN(K,JA,2)+WA2                                  
         SPAN(K,JA,3)=SPAN(K,JA,3)+WA3                                  
         SPAN(K,JA,4)=SPAN(K,JA,4)+WA4                                  
         SPAN(K,JA,5)=SPAN(K,JA,5)+WA5                                  
   21   CONTINUE                                                        
       ENDIF                                                            
  501 CONTINUE                                                          
C                                                                       
C--- Y:JB                                                               
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 502 J=JB,JB                                                    
       IF (J.EQ.JB) THEN                                                
        DO 22 K=1,KG                                                    
        DO 22 I=1,IG                                                    
          CR1(I,K)=UT(I,JB,K,1)                                         
          CR2(I,K)=UT(I,JB,K,2)                                         
          CR3(I,K)=UT(I,JB,K,3)                                         
          CR4(I,K)=P(I,JB,K)-AP(JB)                                     
          CI1(I,K)=0.0D0                                                
          CI2(I,K)=0.0D0                                                
          CI3(I,K)=0.0D0                                                
          CI4(I,K)=0.0D0                                                
   22   CONTINUE                                                        
       ENDIF                                                            
  502 CONTINUE                                                          
C                                                                       
      ISW=1                                                             
C      CALL DFC2BF(IG,KG,CR1,CI1,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR2,CI2,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR3,CI3,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR4,CI4,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,WA1,WA2,WA3,WA4,WA5)                      
      DO 503 J=JB,JB                                                    
       IF (J.EQ.JB) THEN                                                
        DO 23 K=1,KG                                                    
        DO 23 I=1,IG                                                    
         WA1=(CR1(I,K)*PW)**2+(CI1(I,K)*PW)**2                          
         WA2=(CR2(I,K)*PW)**2+(CI2(I,K)*PW)**2                          
         WA3=(CR3(I,K)*PW)**2+(CI3(I,K)*PW)**2                          
         WA4=(CR4(I,K)*PW)**2+(CI4(I,K)*PW)**2                          
         WA5=CR1(I,K)*PW*CR2(I,K)*PW+CI1(I,K)*PW*CI2(I,K)*PW            
C                                                                       
         STRM(I,JB,1)=STRM(I,JB,1)+WA1                                  
         STRM(I,JB,2)=STRM(I,JB,2)+WA2                                  
         STRM(I,JB,3)=STRM(I,JB,3)+WA3                                  
         STRM(I,JB,4)=STRM(I,JB,4)+WA4                                  
         STRM(I,JB,5)=STRM(I,JB,5)+WA5                                  
C                                                                       
         SPAN(K,JB,1)=SPAN(K,JB,1)+WA1                                  
         SPAN(K,JB,2)=SPAN(K,JB,2)+WA2                                  
         SPAN(K,JB,3)=SPAN(K,JB,3)+WA3                                  
         SPAN(K,JB,4)=SPAN(K,JB,4)+WA4                                  
         SPAN(K,JB,5)=SPAN(K,JB,5)+WA5                                  
   23   CONTINUE                                                        
       ENDIF                                                            
  503 CONTINUE                                                          
C                                                                       
C--- Y:JC                                                               
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 504 J=JC,JC                                                    
       IF (J.EQ.JC) THEN                                                
        DO 24 K=1,KG                                                    
        DO 24 I=1,IG                                                    
          CR1(I,K)=UT(I,JC,K,1)                                         
          CR2(I,K)=UT(I,JC,K,2)                                         
          CR3(I,K)=UT(I,JC,K,3)                                         
          CR4(I,K)=P(I,JC,K)-AP(JC)                                     
          CI1(I,K)=0.0D0                                                
          CI2(I,K)=0.0D0                                                
          CI3(I,K)=0.0D0                                                
          CI4(I,K)=0.0D0                                                
   24   CONTINUE                                                        
       ENDIF                                                            
  504 CONTINUE                                                          
C                                                                       
      ISW=1                                                             
C      CALL DFC2BF(IG,KG,CR1,CI1,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR2,CI2,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR3,CI3,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR4,CI4,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,WA1,WA2,WA3,WA4,WA5)                      
      DO 505 J=JC,JC                                                    
       IF (J.EQ.JC) THEN                                                
        DO 25 K=1,KG                                                    
        DO 25 I=1,IG                                                    
         WA1=(CR1(I,K)*PW)**2+(CI1(I,K)*PW)**2                          
         WA2=(CR2(I,K)*PW)**2+(CI2(I,K)*PW)**2                          
         WA3=(CR3(I,K)*PW)**2+(CI3(I,K)*PW)**2                          
         WA4=(CR4(I,K)*PW)**2+(CI4(I,K)*PW)**2                          
         WA5=CR1(I,K)*PW*CR2(I,K)*PW+CI1(I,K)*PW*CI2(I,K)*PW            
C                                                                       
         STRM(I,JC,1)=STRM(I,JC,1)+WA1                                  
         STRM(I,JC,2)=STRM(I,JC,2)+WA2                                  
         STRM(I,JC,3)=STRM(I,JC,3)+WA3                                  
         STRM(I,JC,4)=STRM(I,JC,4)+WA4                                  
         STRM(I,JC,5)=STRM(I,JC,5)+WA5                                  
C                                                                       
         SPAN(K,JC,1)=SPAN(K,JC,1)+WA1                                  
         SPAN(K,JC,2)=SPAN(K,JC,2)+WA2                                  
         SPAN(K,JC,3)=SPAN(K,JC,3)+WA3                                  
         SPAN(K,JC,4)=SPAN(K,JC,4)+WA4                                  
         SPAN(K,JC,5)=SPAN(K,JC,5)+WA5                                  
   25   CONTINUE                                                        
       ENDIF                                                            
  505 CONTINUE                                                          
C                                                                       
C--- Y:JD                                                               
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 506 J=JD,JD                                                    
       IF (J.EQ.JD) THEN                                                
        DO 26 K=1,KG                                                    
        DO 26 I=1,IG                                                    
          CR1(I,K)=UT(I,JD,K,1)                                         
          CR2(I,K)=UT(I,JD,K,2)                                         
          CR3(I,K)=UT(I,JD,K,3)                                         
          CR4(I,K)=P(I,JD,K)-AP(JD)                                     
          CI1(I,K)=0.0D0                                                
          CI2(I,K)=0.0D0                                                
          CI3(I,K)=0.0D0                                                
          CI4(I,K)=0.0D0                                                
   26   CONTINUE                                                        
       ENDIF                                                            
  506 CONTINUE                                                          
C                                                                       
      ISW=1                                                             
C      CALL DFC2BF(IG,KG,CR1,CI1,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR2,CI2,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR3,CI3,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR4,CI4,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,WA1,WA2,WA3,WA4,WA5)                      
      DO 507 J=JD,JD                                                    
       IF (J.EQ.JD) THEN                                                
        DO 27 K=1,KG                                                    
        DO 27 I=1,IG                                                    
         WA1=(CR1(I,K)*PW)**2+(CI1(I,K)*PW)**2                          
         WA2=(CR2(I,K)*PW)**2+(CI2(I,K)*PW)**2                          
         WA3=(CR3(I,K)*PW)**2+(CI3(I,K)*PW)**2                          
         WA4=(CR4(I,K)*PW)**2+(CI4(I,K)*PW)**2                          
         WA5=CR1(I,K)*PW*CR2(I,K)*PW+CI1(I,K)*PW*CI2(I,K)*PW            
C                                                                       
         STRM(I,JD,1)=STRM(I,JD,1)+WA1                                  
         STRM(I,JD,2)=STRM(I,JD,2)+WA2                                  
         STRM(I,JD,3)=STRM(I,JD,3)+WA3                                  
         STRM(I,JD,4)=STRM(I,JD,4)+WA4                                  
         STRM(I,JD,5)=STRM(I,JD,5)+WA5                                  
C                                                                       
         SPAN(K,JD,1)=SPAN(K,JD,1)+WA1                                  
         SPAN(K,JD,2)=SPAN(K,JD,2)+WA2                                  
         SPAN(K,JD,3)=SPAN(K,JD,3)+WA3                                  
         SPAN(K,JD,4)=SPAN(K,JD,4)+WA4                                  
         SPAN(K,JD,5)=SPAN(K,JD,5)+WA5                                  
   27   CONTINUE                                                        
       ENDIF                                                            
  507 CONTINUE                                                          
C                                                                       
C--- Y:JE                                                               
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 508 J=JE,JE                                                    
       IF (J.EQ.JE) THEN                                                
        DO 28 K=1,KG                                                    
        DO 28 I=1,IG                                                    
          CR1(I,K)=UT(I,JE,K,1)                                         
          CR2(I,K)=UT(I,JE,K,2)                                         
          CR3(I,K)=UT(I,JE,K,3)                                         
          CR4(I,K)=P(I,JE,K)-AP(JE)                                     
          CI1(I,K)=0.0D0                                                
          CI2(I,K)=0.0D0                                                
          CI3(I,K)=0.0D0                                                
          CI4(I,K)=0.0D0                                                
   28   CONTINUE                                                        
       ENDIF                                                            
  508 CONTINUE                                                          
C                                                                       
      ISW=1                                                             
C      CALL DFC2BF(IG,KG,CR1,CI1,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR2,CI2,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR3,CI3,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR4,CI4,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,WA1,WA2,WA3,WA4,WA5)                      
      DO 509 J=JE,JE                                                    
       IF (J.EQ.JE) THEN                                                
        DO 29 K=1,KG                                                    
        DO 29 I=1,IG                                                    
         WA1=(CR1(I,K)*PW)**2+(CI1(I,K)*PW)**2                          
         WA2=(CR2(I,K)*PW)**2+(CI2(I,K)*PW)**2                          
         WA3=(CR3(I,K)*PW)**2+(CI3(I,K)*PW)**2                          
         WA4=(CR4(I,K)*PW)**2+(CI4(I,K)*PW)**2                          
         WA5=CR1(I,K)*PW*CR2(I,K)*PW+CI1(I,K)*PW*CI2(I,K)*PW            
C                                                                       
         STRM(I,JE,1)=STRM(I,JE,1)+WA1                                  
         STRM(I,JE,2)=STRM(I,JE,2)+WA2                                  
         STRM(I,JE,3)=STRM(I,JE,3)+WA3                                  
         STRM(I,JE,4)=STRM(I,JE,4)+WA4                                  
         STRM(I,JE,5)=STRM(I,JE,5)+WA5                                  
C                                                                       
         SPAN(K,JE,1)=SPAN(K,JE,1)+WA1                                  
         SPAN(K,JE,2)=SPAN(K,JE,2)+WA2                                  
         SPAN(K,JE,3)=SPAN(K,JE,3)+WA3                                  
         SPAN(K,JE,4)=SPAN(K,JE,4)+WA4                                  
         SPAN(K,JE,5)=SPAN(K,JE,5)+WA5                                  
   29   CONTINUE                                                        
       ENDIF                                                            
  509 CONTINUE                                                          
C                                                                       
C--- Y:JF                                                               
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 510 J=JF,JF                                                    
       IF (J.EQ.JF) THEN                                                
        DO 30 K=1,KG                                                    
        DO 30 I=1,IG                                                    
          CR1(I,K)=UT(I,JF,K,1)                                         
          CR2(I,K)=UT(I,JF,K,2)                                         
          CR3(I,K)=UT(I,JF,K,3)                                         
          CR4(I,K)=P(I,JF,K)-AP(JF)                                     
          CI1(I,K)=0.0D0                                                
          CI2(I,K)=0.0D0                                                
          CI3(I,K)=0.0D0                                                
          CI4(I,K)=0.0D0                                                
   30   CONTINUE                                                        
       ENDIF                                                            
  510 CONTINUE                                                          
C                                                                       
      ISW=1                                                             
C      CALL DFC2BF(IG,KG,CR1,CI1,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR2,CI2,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR3,CI3,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR4,CI4,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,WA1,WA2,WA3,WA4,WA5)                      
      DO 511 J=JF,JF                                                    
       IF (J.EQ.JF) THEN                                                
        DO 31 K=1,KG                                                    
        DO 31 I=1,IG                                                    
         WA1=(CR1(I,K)*PW)**2+(CI1(I,K)*PW)**2                          
         WA2=(CR2(I,K)*PW)**2+(CI2(I,K)*PW)**2                          
         WA3=(CR3(I,K)*PW)**2+(CI3(I,K)*PW)**2                          
         WA4=(CR4(I,K)*PW)**2+(CI4(I,K)*PW)**2                          
         WA5=CR1(I,K)*PW*CR2(I,K)*PW+CI1(I,K)*PW*CI2(I,K)*PW            
C                                                                       
         STRM(I,JF,1)=STRM(I,JF,1)+WA1                                  
         STRM(I,JF,2)=STRM(I,JF,2)+WA2                                  
         STRM(I,JF,3)=STRM(I,JF,3)+WA3                                  
         STRM(I,JF,4)=STRM(I,JF,4)+WA4                                  
         STRM(I,JF,5)=STRM(I,JF,5)+WA5                                  
C                                                                       
         SPAN(K,JF,1)=SPAN(K,JF,1)+WA1                                  
         SPAN(K,JF,2)=SPAN(K,JF,2)+WA2                                  
         SPAN(K,JF,3)=SPAN(K,JF,3)+WA3                                  
         SPAN(K,JF,4)=SPAN(K,JF,4)+WA4                                  
         SPAN(K,JF,5)=SPAN(K,JF,5)+WA5                                  
   31   CONTINUE                                                        
       ENDIF                                                            
  511 CONTINUE                                                          
C                                                                       
C--- Y:JH                                                               
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 512 J=JH,JH                                                    
       IF (J.EQ.JH) THEN                                                
        DO 32 K=1,KG                                                    
        DO 32 I=1,IG                                                    
          CR1(I,K)=UT(I,JH,K,1)                                         
          CR2(I,K)=UT(I,JH,K,2)                                         
          CR3(I,K)=UT(I,JH,K,3)                                         
          CR4(I,K)=P(I,JH,K)-AP(JH)                                     
          CI1(I,K)=0.0D0                                                
          CI2(I,K)=0.0D0                                                
          CI3(I,K)=0.0D0                                                
          CI4(I,K)=0.0D0                                                
   32   CONTINUE                                                        
       ENDIF                                                            
  512 CONTINUE                                                          
C                                                                       
      ISW=1                                                             
C      CALL DFC2BF(IG,KG,CR1,CI1,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR2,CI2,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR3,CI3,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR4,CI4,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,WA1,WA2,WA3,WA4,WA5)                      
      DO 513 J=JH,JH                                                    
       IF (J.EQ.JH) THEN                                                
        DO 33 K=1,KG                                                    
        DO 33 I=1,IG                                                    
         WA1=(CR1(I,K)*PW)**2+(CI1(I,K)*PW)**2                          
         WA2=(CR2(I,K)*PW)**2+(CI2(I,K)*PW)**2                          
         WA3=(CR3(I,K)*PW)**2+(CI3(I,K)*PW)**2                          
         WA4=(CR4(I,K)*PW)**2+(CI4(I,K)*PW)**2                          
         WA5=CR1(I,K)*PW*CR2(I,K)*PW+CI1(I,K)*PW*CI2(I,K)*PW            
C                                                                       
         STRM(I,JH,1)=STRM(I,JH,1)+WA1                                  
         STRM(I,JH,2)=STRM(I,JH,2)+WA2                                  
         STRM(I,JH,3)=STRM(I,JH,3)+WA3                                  
         STRM(I,JH,4)=STRM(I,JH,4)+WA4                                  
         STRM(I,JH,5)=STRM(I,JH,5)+WA5                                  
C                                                                       
         SPAN(K,JH,1)=SPAN(K,JH,1)+WA1                                  
         SPAN(K,JH,2)=SPAN(K,JH,2)+WA2                                  
         SPAN(K,JH,3)=SPAN(K,JH,3)+WA3                                  
         SPAN(K,JH,4)=SPAN(K,JH,4)+WA4                                  
         SPAN(K,JH,5)=SPAN(K,JH,5)+WA5                                  
   33   CONTINUE                                                        
       ENDIF                                                            
  513 CONTINUE                                                          
C                                                                       
C--- Y:JI                                                               
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 514 J=JI,JI                                                    
       IF (J.EQ.JI) THEN                                                
        DO 34 K=1,KG                                                    
        DO 34 I=1,IG                                                    
          CR1(I,K)=UT(I,JI,K,1)                                         
          CR2(I,K)=UT(I,JI,K,2)                                         
          CR3(I,K)=UT(I,JI,K,3)                                         
          CR4(I,K)=P(I,JI,K)-AP(JI)                                     
          CI1(I,K)=0.0D0                                                
          CI2(I,K)=0.0D0                                                
          CI3(I,K)=0.0D0                                                
          CI4(I,K)=0.0D0                                                
   34   CONTINUE                                                        
       ENDIF                                                            
  514 CONTINUE                                                          
C                                                                       
      ISW=1                                                             
C      CALL DFC2BF(IG,KG,CR1,CI1,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR2,CI2,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR3,CI3,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR4,CI4,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,WA1,WA2,WA3,WA4,WA5)                      
      DO 515 J=JI,JI                                                    
       IF (J.EQ.JI) THEN                                                
        DO 35 K=1,KG                                                    
        DO 35 I=1,IG                                                    
         WA1=(CR1(I,K)*PW)**2+(CI1(I,K)*PW)**2                          
         WA2=(CR2(I,K)*PW)**2+(CI2(I,K)*PW)**2                          
         WA3=(CR3(I,K)*PW)**2+(CI3(I,K)*PW)**2                          
         WA4=(CR4(I,K)*PW)**2+(CI4(I,K)*PW)**2                          
         WA5=CR1(I,K)*PW*CR2(I,K)*PW+CI1(I,K)*PW*CI2(I,K)*PW            
C                                                                       
         STRM(I,JI,1)=STRM(I,JI,1)+WA1                                  
         STRM(I,JI,2)=STRM(I,JI,2)+WA2                                  
         STRM(I,JI,3)=STRM(I,JI,3)+WA3                                  
         STRM(I,JI,4)=STRM(I,JI,4)+WA4                                  
         STRM(I,JI,5)=STRM(I,JI,5)+WA5                                  
C                                                                       
         SPAN(K,JI,1)=SPAN(K,JI,1)+WA1                                  
         SPAN(K,JI,2)=SPAN(K,JI,2)+WA2                                  
         SPAN(K,JI,3)=SPAN(K,JI,3)+WA3                                  
         SPAN(K,JI,4)=SPAN(K,JI,4)+WA4                                  
         SPAN(K,JI,5)=SPAN(K,JI,5)+WA5                                  
   35   CONTINUE                                                        
       ENDIF                                                            
  515 CONTINUE                                                          
C                                                                       
C--- Y:JJ                                                               
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 516 J=JJ,JJ                                                    
       IF (J.EQ.JJ) THEN                                                
        DO 36 K=1,KG                                                    
        DO 36 I=1,IG                                                    
          CR1(I,K)=UT(I,JJ,K,1)                                         
          CR2(I,K)=UT(I,JJ,K,2)                                         
          CR3(I,K)=UT(I,JJ,K,3)                                         
          CR4(I,K)=P(I,JJ,K)-AP(JJ)                                     
          CI1(I,K)=0.0D0                                                
          CI2(I,K)=0.0D0                                                
          CI3(I,K)=0.0D0                                                
          CI4(I,K)=0.0D0                                                
   36   CONTINUE                                                        
       ENDIF                                                            
  516 CONTINUE                                                          
C                                                                       
      ISW=1                                                             
C      CALL DFC2BF(IG,KG,CR1,CI1,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR2,CI2,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR3,CI3,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR4,CI4,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,WA1,WA2,WA3,WA4,WA5)                      
      DO 517 J=JJ,JJ                                                    
       IF (J.EQ.JJ) THEN                                                
        DO 37 K=1,KG                                                    
        DO 37 I=1,IG                                                    
         WA1=(CR1(I,K)*PW)**2+(CI1(I,K)*PW)**2                          
         WA2=(CR2(I,K)*PW)**2+(CI2(I,K)*PW)**2                          
         WA3=(CR3(I,K)*PW)**2+(CI3(I,K)*PW)**2                          
         WA4=(CR4(I,K)*PW)**2+(CI4(I,K)*PW)**2                          
         WA5=CR1(I,K)*PW*CR2(I,K)*PW+CI1(I,K)*PW*CI2(I,K)*PW            
C                                                                       
         STRM(I,JJ,1)=STRM(I,JJ,1)+WA1                                  
         STRM(I,JJ,2)=STRM(I,JJ,2)+WA2                                  
         STRM(I,JJ,3)=STRM(I,JJ,3)+WA3                                  
         STRM(I,JJ,4)=STRM(I,JJ,4)+WA4                                  
         STRM(I,JJ,5)=STRM(I,JJ,5)+WA5                                  
C                                                                       
         SPAN(K,JJ,1)=SPAN(K,JJ,1)+WA1                                  
         SPAN(K,JJ,2)=SPAN(K,JJ,2)+WA2                                  
         SPAN(K,JJ,3)=SPAN(K,JJ,3)+WA3                                  
         SPAN(K,JJ,4)=SPAN(K,JJ,4)+WA4                                  
         SPAN(K,JJ,5)=SPAN(K,JJ,5)+WA5                                  
   37   CONTINUE                                                        
       ENDIF                                                            
  517 CONTINUE                                                          
C                                                                       
C--- Y:JK                                                               
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 518 J=JK,JK                                                    
       IF (J.EQ.JK) THEN                                                
        DO 38 K=1,KG                                                    
        DO 38 I=1,IG                                                    
          CR1(I,K)=UT(I,JK,K,1)                                         
          CR2(I,K)=UT(I,JK,K,2)                                         
          CR3(I,K)=UT(I,JK,K,3)                                         
          CR4(I,K)=P(I,JK,K)-AP(JK)                                     
          CI1(I,K)=0.0D0                                                
          CI2(I,K)=0.0D0                                                
          CI3(I,K)=0.0D0                                                
          CI4(I,K)=0.0D0                                                
   38   CONTINUE                                                        
       ENDIF                                                            
  518 CONTINUE                                                          
C                                                                       
      ISW=1                                                             
C      CALL DFC2BF(IG,KG,CR1,CI1,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR2,CI2,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR3,CI3,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR4,CI4,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,WA1,WA2,WA3,WA4,WA5)                      
      DO 519 J=JK,JK                                                    
       IF (J.EQ.JK) THEN                                                
        DO 39 K=1,KG                                                    
        DO 39 I=1,IG                                                    
         WA1=(CR1(I,K)*PW)**2+(CI1(I,K)*PW)**2                          
         WA2=(CR2(I,K)*PW)**2+(CI2(I,K)*PW)**2                          
         WA3=(CR3(I,K)*PW)**2+(CI3(I,K)*PW)**2                          
         WA4=(CR4(I,K)*PW)**2+(CI4(I,K)*PW)**2                          
         WA5=CR1(I,K)*PW*CR2(I,K)*PW+CI1(I,K)*PW*CI2(I,K)*PW            
C                                                                       
         STRM(I,JK,1)=STRM(I,JK,1)+WA1                                  
         STRM(I,JK,2)=STRM(I,JK,2)+WA2                                  
         STRM(I,JK,3)=STRM(I,JK,3)+WA3                                  
         STRM(I,JK,4)=STRM(I,JK,4)+WA4                                  
         STRM(I,JK,5)=STRM(I,JK,5)+WA5                                  
C                                                                       
         SPAN(K,JK,1)=SPAN(K,JK,1)+WA1                                  
         SPAN(K,JK,2)=SPAN(K,JK,2)+WA2                                  
         SPAN(K,JK,3)=SPAN(K,JK,3)+WA3                                  
         SPAN(K,JK,4)=SPAN(K,JK,4)+WA4                                  
         SPAN(K,JK,5)=SPAN(K,JK,5)+WA5                                  
   39   CONTINUE                                                        
       ENDIF                                                            
  519 CONTINUE                                                          
C                                                                       
C--- Y:JL                                                               
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 520 J=JL,JL                                                    
       IF (J.EQ.JL) THEN                                                
        DO 40 K=1,KG                                                    
        DO 40 I=1,IG                                                    
          CR1(I,K)=UT(I,JL,K,1)                                         
          CR2(I,K)=UT(I,JL,K,2)                                         
          CR3(I,K)=UT(I,JL,K,3)                                         
          CR4(I,K)=P(I,JL,K)-AP(JL)                                     
          CI1(I,K)=0.0D0                                                
          CI2(I,K)=0.0D0                                                
          CI3(I,K)=0.0D0                                                
          CI4(I,K)=0.0D0                                                
   40   CONTINUE                                                        
       ENDIF                                                            
  520 CONTINUE                                                          
C                                                                       
      ISW=1                                                             
C      CALL DFC2BF(IG,KG,CR1,CI1,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR2,CI2,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR3,CI3,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C      CALL DFC2BF(IG,KG,CR4,CI4,LX,LZ,ISW,IFAX,TRIGS,WK,IERR)           
      IF(IERR.NE.0) STOP                                                
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K,WA1,WA2,WA3,WA4,WA5)                      
      DO 521 J=JL,JL                                                    
       IF (J.EQ.JL) THEN                                                
        DO 41 K=1,KG                                                    
        DO 41 I=1,IG                                                    
         WA1=(CR1(I,K)*PW)**2+(CI1(I,K)*PW)**2                          
         WA2=(CR2(I,K)*PW)**2+(CI2(I,K)*PW)**2                          
         WA3=(CR3(I,K)*PW)**2+(CI3(I,K)*PW)**2                          
         WA4=(CR4(I,K)*PW)**2+(CI4(I,K)*PW)**2                          
         WA5=CR1(I,K)*PW*CR2(I,K)*PW+CI1(I,K)*PW*CI2(I,K)*PW            
C                                                                       
         STRM(I,JL,1)=STRM(I,JL,1)+WA1                                  
         STRM(I,JL,2)=STRM(I,JL,2)+WA2                                  
         STRM(I,JL,3)=STRM(I,JL,3)+WA3                                  
         STRM(I,JL,4)=STRM(I,JL,4)+WA4                                  
         STRM(I,JL,5)=STRM(I,JL,5)+WA5                                  
C                                                                       
         SPAN(K,JL,1)=SPAN(K,JL,1)+WA1                                  
         SPAN(K,JL,2)=SPAN(K,JL,2)+WA2                                  
         SPAN(K,JL,3)=SPAN(K,JL,3)+WA3                                  
         SPAN(K,JL,4)=SPAN(K,JL,4)+WA4                                  
         SPAN(K,JL,5)=SPAN(K,JL,5)+WA5                                  
   41   CONTINUE                                                        
       ENDIF                                                            
  521 CONTINUE                                                          
C                                                                       
C-------------------------------------- TIME AVE.                       
C====================================== STREAMWISE                      
C                                                                       
        DO 100 I=1,IG                                                   
          ASTRM(I,JA,1)=ASTRM(I,JA,1)+STRM(I,JA,1)                      
          ASTRM(I,JA,2)=ASTRM(I,JA,2)+STRM(I,JA,2)                      
          ASTRM(I,JA,3)=ASTRM(I,JA,3)+STRM(I,JA,3)                      
          ASTRM(I,JA,4)=ASTRM(I,JA,4)+STRM(I,JA,4)                      
          ASTRM(I,JA,5)=ASTRM(I,JA,5)+STRM(I,JA,5)                      
  100   CONTINUE                                                        
C                                                                       
        DO 101 I=1,IG                                                   
          ASTRM(I,JB,1)=ASTRM(I,JB,1)+STRM(I,JB,1)                      
          ASTRM(I,JB,2)=ASTRM(I,JB,2)+STRM(I,JB,2)                      
          ASTRM(I,JB,3)=ASTRM(I,JB,3)+STRM(I,JB,3)                      
          ASTRM(I,JB,4)=ASTRM(I,JB,4)+STRM(I,JB,4)                      
          ASTRM(I,JB,5)=ASTRM(I,JB,5)+STRM(I,JB,5)                      
  101   CONTINUE                                                        
C                                                                       
        DO 102 I=1,IG                                                   
          ASTRM(I,JC,1)=ASTRM(I,JC,1)+STRM(I,JC,1)                      
          ASTRM(I,JC,2)=ASTRM(I,JC,2)+STRM(I,JC,2)                      
          ASTRM(I,JC,3)=ASTRM(I,JC,3)+STRM(I,JC,3)                      
          ASTRM(I,JC,4)=ASTRM(I,JC,4)+STRM(I,JC,4)                      
          ASTRM(I,JC,5)=ASTRM(I,JC,5)+STRM(I,JC,5)                      
  102   CONTINUE                                                        
C                                                                       
        DO 103 I=1,IG                                                   
          ASTRM(I,JD,1)=ASTRM(I,JD,1)+STRM(I,JD,1)                      
          ASTRM(I,JD,2)=ASTRM(I,JD,2)+STRM(I,JD,2)                      
          ASTRM(I,JD,3)=ASTRM(I,JD,3)+STRM(I,JD,3)                      
          ASTRM(I,JD,4)=ASTRM(I,JD,4)+STRM(I,JD,4)                      
          ASTRM(I,JD,5)=ASTRM(I,JD,5)+STRM(I,JD,5)                      
  103   CONTINUE                                                        
C                                                                       
        DO 104 I=1,IG                                                   
          ASTRM(I,JE,1)=ASTRM(I,JE,1)+STRM(I,JE,1)                      
          ASTRM(I,JE,2)=ASTRM(I,JE,2)+STRM(I,JE,2)                      
          ASTRM(I,JE,3)=ASTRM(I,JE,3)+STRM(I,JE,3)                      
          ASTRM(I,JE,4)=ASTRM(I,JE,4)+STRM(I,JE,4)                      
          ASTRM(I,JE,5)=ASTRM(I,JE,5)+STRM(I,JE,5)                      
  104   CONTINUE                                                        
C                                                                       
        DO 105 I=1,IG                                                   
          ASTRM(I,JF,1)=ASTRM(I,JF,1)+STRM(I,JF,1)                      
          ASTRM(I,JF,2)=ASTRM(I,JF,2)+STRM(I,JF,2)                      
          ASTRM(I,JF,3)=ASTRM(I,JF,3)+STRM(I,JF,3)                      
          ASTRM(I,JF,4)=ASTRM(I,JF,4)+STRM(I,JF,4)                      
          ASTRM(I,JF,5)=ASTRM(I,JF,5)+STRM(I,JF,5)                      
  105   CONTINUE                                                        
C                                                                       
        DO 106 I=1,IG                                                   
          ASTRM(I,JH,1)=ASTRM(I,JH,1)+STRM(I,JH,1)                      
          ASTRM(I,JH,2)=ASTRM(I,JH,2)+STRM(I,JH,2)                      
          ASTRM(I,JH,3)=ASTRM(I,JH,3)+STRM(I,JH,3)                      
          ASTRM(I,JH,4)=ASTRM(I,JH,4)+STRM(I,JH,4)                      
          ASTRM(I,JH,5)=ASTRM(I,JH,5)+STRM(I,JH,5)                      
  106   CONTINUE                                                        
C                                                                       
        DO 107 I=1,IG                                                   
          ASTRM(I,JI,1)=ASTRM(I,JI,1)+STRM(I,JI,1)                      
          ASTRM(I,JI,2)=ASTRM(I,JI,2)+STRM(I,JI,2)                      
          ASTRM(I,JI,3)=ASTRM(I,JI,3)+STRM(I,JI,3)                      
          ASTRM(I,JI,4)=ASTRM(I,JI,4)+STRM(I,JI,4)                      
          ASTRM(I,JI,5)=ASTRM(I,JI,5)+STRM(I,JI,5)                      
  107   CONTINUE                                                        
C                                                                       
        DO 108 I=1,IG                                                   
          ASTRM(I,JJ,1)=ASTRM(I,JJ,1)+STRM(I,JJ,1)                      
          ASTRM(I,JJ,2)=ASTRM(I,JJ,2)+STRM(I,JJ,2)                      
          ASTRM(I,JJ,3)=ASTRM(I,JJ,3)+STRM(I,JJ,3)                      
          ASTRM(I,JJ,4)=ASTRM(I,JJ,4)+STRM(I,JJ,4)                      
          ASTRM(I,JJ,5)=ASTRM(I,JJ,5)+STRM(I,JJ,5)                      
  108   CONTINUE                                                        
C                                                                       
        DO 109 I=1,IG                                                   
          ASTRM(I,JK,1)=ASTRM(I,JK,1)+STRM(I,JK,1)                      
          ASTRM(I,JK,2)=ASTRM(I,JK,2)+STRM(I,JK,2)                      
          ASTRM(I,JK,3)=ASTRM(I,JK,3)+STRM(I,JK,3)                      
          ASTRM(I,JK,4)=ASTRM(I,JK,4)+STRM(I,JK,4)                      
          ASTRM(I,JK,5)=ASTRM(I,JK,5)+STRM(I,JK,5)                      
  109   CONTINUE                                                        
C                                                                       
        DO 110 I=1,IG                                                   
          ASTRM(I,JL,1)=ASTRM(I,JL,1)+STRM(I,JL,1)                      
          ASTRM(I,JL,2)=ASTRM(I,JL,2)+STRM(I,JL,2)                      
          ASTRM(I,JL,3)=ASTRM(I,JL,3)+STRM(I,JL,3)                      
          ASTRM(I,JL,4)=ASTRM(I,JL,4)+STRM(I,JL,4)                      
          ASTRM(I,JL,5)=ASTRM(I,JL,5)+STRM(I,JL,5)                      
  110   CONTINUE                                                        
C                                                                       
C====================================== SPANWISE                        
C                                                                       
        DO 200 K=1,KG                                                   
          ASPAN(K,JA,1)=ASPAN(K,JA,1)+SPAN(K,JA,1)                      
          ASPAN(K,JA,2)=ASPAN(K,JA,2)+SPAN(K,JA,2)                      
          ASPAN(K,JA,3)=ASPAN(K,JA,3)+SPAN(K,JA,3)                      
          ASPAN(K,JA,4)=ASPAN(K,JA,4)+SPAN(K,JA,4)                      
          ASPAN(K,JA,5)=ASPAN(K,JA,5)+SPAN(K,JA,5)                      
  200   CONTINUE                                                        
C                                                                       
        DO 201 K=1,KG                                                   
          ASPAN(K,JB,1)=ASPAN(K,JB,1)+SPAN(K,JB,1)                      
          ASPAN(K,JB,2)=ASPAN(K,JB,2)+SPAN(K,JB,2)                      
          ASPAN(K,JB,3)=ASPAN(K,JB,3)+SPAN(K,JB,3)                      
          ASPAN(K,JB,4)=ASPAN(K,JB,4)+SPAN(K,JB,4)                      
          ASPAN(K,JB,5)=ASPAN(K,JB,5)+SPAN(K,JB,5)                      
  201   CONTINUE                                                        
C                                                                       
        DO 202 K=1,KG                                                   
          ASPAN(K,JC,1)=ASPAN(K,JC,1)+SPAN(K,JC,1)                      
          ASPAN(K,JC,2)=ASPAN(K,JC,2)+SPAN(K,JC,2)                      
          ASPAN(K,JC,3)=ASPAN(K,JC,3)+SPAN(K,JC,3)                      
          ASPAN(K,JC,4)=ASPAN(K,JC,4)+SPAN(K,JC,4)                      
          ASPAN(K,JC,5)=ASPAN(K,JC,5)+SPAN(K,JC,5)                      
  202   CONTINUE                                                        
C                                                                       
        DO 203 K=1,KG                                                   
          ASPAN(K,JD,1)=ASPAN(K,JD,1)+SPAN(K,JD,1)                      
          ASPAN(K,JD,2)=ASPAN(K,JD,2)+SPAN(K,JD,2)                      
          ASPAN(K,JD,3)=ASPAN(K,JD,3)+SPAN(K,JD,3)                      
          ASPAN(K,JD,4)=ASPAN(K,JD,4)+SPAN(K,JD,4)                      
          ASPAN(K,JD,5)=ASPAN(K,JD,5)+SPAN(K,JD,5)                      
  203   CONTINUE                                                        
C                                                                       
        DO 204 K=1,KG                                                   
          ASPAN(K,JE,1)=ASPAN(K,JE,1)+SPAN(K,JE,1)                      
          ASPAN(K,JE,2)=ASPAN(K,JE,2)+SPAN(K,JE,2)                      
          ASPAN(K,JE,3)=ASPAN(K,JE,3)+SPAN(K,JE,3)                      
          ASPAN(K,JE,4)=ASPAN(K,JE,4)+SPAN(K,JE,4)                      
          ASPAN(K,JE,5)=ASPAN(K,JE,5)+SPAN(K,JE,5)                      
  204   CONTINUE                                                        
C                                                                       
        DO 205 K=1,KG                                                   
          ASPAN(K,JF,1)=ASPAN(K,JF,1)+SPAN(K,JF,1)                      
          ASPAN(K,JF,2)=ASPAN(K,JF,2)+SPAN(K,JF,2)                      
          ASPAN(K,JF,3)=ASPAN(K,JF,3)+SPAN(K,JF,3)                      
          ASPAN(K,JF,4)=ASPAN(K,JF,4)+SPAN(K,JF,4)                      
          ASPAN(K,JF,5)=ASPAN(K,JF,5)+SPAN(K,JF,5)                      
  205   CONTINUE                                                        
C                                                                       
        DO 206 K=1,KG                                                   
          ASPAN(K,JH,1)=ASPAN(K,JH,1)+SPAN(K,JH,1)                      
          ASPAN(K,JH,2)=ASPAN(K,JH,2)+SPAN(K,JH,2)                      
          ASPAN(K,JH,3)=ASPAN(K,JH,3)+SPAN(K,JH,3)                      
          ASPAN(K,JH,4)=ASPAN(K,JH,4)+SPAN(K,JH,4)                      
          ASPAN(K,JH,5)=ASPAN(K,JH,5)+SPAN(K,JH,5)                      
  206   CONTINUE                                                        
C                                                                       
        DO 207 K=1,KG                                                   
          ASPAN(K,JI,1)=ASPAN(K,JI,1)+SPAN(K,JI,1)                      
          ASPAN(K,JI,2)=ASPAN(K,JI,2)+SPAN(K,JI,2)                      
          ASPAN(K,JI,3)=ASPAN(K,JI,3)+SPAN(K,JI,3)                      
          ASPAN(K,JI,4)=ASPAN(K,JI,4)+SPAN(K,JI,4)                      
          ASPAN(K,JI,5)=ASPAN(K,JI,5)+SPAN(K,JI,5)                      
  207   CONTINUE                                                        
C                                                                       
        DO 208 K=1,KG                                                   
          ASPAN(K,JJ,1)=ASPAN(K,JJ,1)+SPAN(K,JJ,1)                      
          ASPAN(K,JJ,2)=ASPAN(K,JJ,2)+SPAN(K,JJ,2)                      
          ASPAN(K,JJ,3)=ASPAN(K,JJ,3)+SPAN(K,JJ,3)                      
          ASPAN(K,JJ,4)=ASPAN(K,JJ,4)+SPAN(K,JJ,4)                      
          ASPAN(K,JJ,5)=ASPAN(K,JJ,5)+SPAN(K,JJ,5)                      
  208   CONTINUE                                                        
C                                                                       
        DO 209 K=1,KG                                                   
          ASPAN(K,JK,1)=ASPAN(K,JK,1)+SPAN(K,JK,1)                      
          ASPAN(K,JK,2)=ASPAN(K,JK,2)+SPAN(K,JK,2)                      
          ASPAN(K,JK,3)=ASPAN(K,JK,3)+SPAN(K,JK,3)                      
          ASPAN(K,JK,4)=ASPAN(K,JK,4)+SPAN(K,JK,4)                      
          ASPAN(K,JK,5)=ASPAN(K,JK,5)+SPAN(K,JK,5)                      
  209   CONTINUE                                                        
C                                                                       
        DO 210 K=1,KG                                                   
          ASPAN(K,JL,1)=ASPAN(K,JL,1)+SPAN(K,JL,1)                      
          ASPAN(K,JL,2)=ASPAN(K,JL,2)+SPAN(K,JL,2)                      
          ASPAN(K,JL,3)=ASPAN(K,JL,3)+SPAN(K,JL,3)                      
          ASPAN(K,JL,4)=ASPAN(K,JL,4)+SPAN(K,JL,4)                      
          ASPAN(K,JL,5)=ASPAN(K,JL,5)+SPAN(K,JL,5)                      
  210   CONTINUE                                                        
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     TWO-POINT CORRELATION (VELOCITY FIELD)                            
C     PROGRAMED BY H.ABE                                                
      SUBROUTINE CORREL                                                 
C     REVISED ON 2001.03.03 BY H.ABE (JN=1)                             
C     REVISED ON 2001.01.08 BY H.ABE (JA TO JM)                         
C     2000.08.26 REVISED BY H.ABE (SPREAD REGION -> SPREAD DO)          
C====================================================================   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
C                                                                       
      COMMON /JAVE/ JA,JB,JC,JD,JE,JF,JH,JI,JJ,JK,JL                    
C                                                                       
C----LOCAL VARIABLES-------------------------------------------         
C                                                                       
      DIMENSION U(-4:IG+4,-2:JG+2,-4:KG+4,3)                            
      DIMENSION P(-3:IG+3,-2:JG+2,-3:KG+3)                              
      DIMENSION UT(-4:IG+4,-2:JG+2,-4:KG+4,3)                           
      DIMENSION AU(-1:JG+2),AV(-1:JG+2),AW(-1:JG+2),AP(-1:JG+2)         
      DIMENSION ATPC1(IG,JG,4),ATPC3(KG,JG,4)                           
      DIMENSION R1A1(IG,2),R2A1(IG,2),R3A1(IG,2)                        
     &         ,R1B1(IG,2),R2B1(IG,2),R3B1(IG,2)                        
     &         ,R1C1(IG,2),R2C1(IG,2),R3C1(IG,2)                        
     &         ,R1D1(IG,2),R2D1(IG,2),R3D1(IG,2)                        
     &         ,R1E1(IG,2),R2E1(IG,2),R3E1(IG,2)                        
     &         ,R1F1(IG,2),R2F1(IG,2),R3F1(IG,2)                        
     &         ,R1H1(IG,2),R2H1(IG,2),R3H1(IG,2)                        
     &         ,R1I1(IG,2),R2I1(IG,2),R3I1(IG,2)                        
     &         ,R1J1(IG,2),R2J1(IG,2),R3J1(IG,2)                        
     &         ,R1K1(IG,2),R2K1(IG,2),R3K1(IG,2)                        
     &         ,R1L1(IG,2),R2L1(IG,2),R3L1(IG,2)                        
     &         ,R1A3(KG,2),R2A3(KG,2),R3A3(KG,2)                        
     &         ,R1B3(KG,2),R2B3(KG,2),R3B3(KG,2)                        
     &         ,R1C3(KG,2),R2C3(KG,2),R3C3(KG,2)                        
     &         ,R1D3(KG,2),R2D3(KG,2),R3D3(KG,2)                        
     &         ,R1E3(KG,2),R2E3(KG,2),R3E3(KG,2)                        
     &         ,R1F3(KG,2),R2F3(KG,2),R3F3(KG,2)                        
     &         ,R1H3(KG,2),R2H3(KG,2),R3H3(KG,2)                        
     &         ,R1I3(KG,2),R2I3(KG,2),R3I3(KG,2)                        
     &         ,R1J3(KG,2),R2J3(KG,2),R3J3(KG,2)                        
     &         ,R1K3(KG,2),R2K3(KG,2),R3K3(KG,2)                        
     &         ,R1L3(KG,2),R2L3(KG,2),R3L3(KG,2)                        
C                                                                       
      DIMENSION R4A1(IG,2),R4B1(IG,2),R4C1(IG,2)                        
     &         ,R4D1(IG,2),R4E1(IG,2),R4F1(IG,2)                        
     &         ,R4H1(IG,2),R4I1(IG,2),R4J1(IG,2)                        
     &         ,R4K1(IG,2),R4L1(IG,2)                                   
     &         ,R4A3(KG,2),R4B3(KG,2),R4C3(KG,2)                        
     &         ,R4D3(KG,2),R4E3(KG,2),R4F3(KG,2)                        
     &         ,R4H3(KG,2),R4I3(KG,2),R4J3(KG,2)                        
     &         ,R4K3(KG,2),R4L3(KG,2)                                   
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /VELO/ UG(-4:IG+4,-2:JG+2,-4:KG+4,3)                       
      COMMON /SCAL/ PG(-3:IG+3,-2:JG+2,-3:KG+3)                         
      COMMON /TEMP/ UTG(-4:IG+4,-2:JG+2,-4:KG+4,3)                      
      COMMON /MEAN/ AUG(-1:JG+2),AVG(-1:JG+2)                           
     &             ,AWG(-1:JG+2),APG(-1:JG+2)                           
      COMMON /CORR/ ATPC1G(IG,JG,4),ATPC3G(KG,JG,4)                     
C                                                                       
C----COMBINATION------------------------------------------------        
C                                                                       
      EQUIVALENCE (UG,U),(UTG,UT),(PG,P)                                
      EQUIVALENCE (AUG,AU),(AVG,AV),(AWG,AW),(APG,AP)                   
      EQUIVALENCE (ATPC1G,ATPC1),(ATPC3G,ATPC3)                         
C                                                                       
C=================================================================      
C                                                                       
C                                                                       
C=================================== ZERO-CLEAR                         
C----------------------------------- STREAMWISE                         
C                                                                       
        DO 10 I=1,IG                                                    
          R1A1(I,1)=0.0D0                                               
          R1A1(I,2)=0.0D0                                               
          R2A1(I,1)=0.0D0                                               
          R2A1(I,2)=0.0D0                                               
          R3A1(I,1)=0.0D0                                               
          R3A1(I,2)=0.0D0                                               
          R4A1(I,1)=0.0D0                                               
          R4A1(I,2)=0.0D0                                               
   10   CONTINUE                                                        
C                                                                       
        DO 11 I=1,IG                                                    
          R1B1(I,1)=0.0D0                                               
          R1B1(I,2)=0.0D0                                               
          R2B1(I,1)=0.0D0                                               
          R2B1(I,2)=0.0D0                                               
          R3B1(I,1)=0.0D0                                               
          R3B1(I,2)=0.0D0                                               
          R4B1(I,1)=0.0D0                                               
          R4B1(I,2)=0.0D0                                               
   11   CONTINUE                                                        
C                                                                       
        DO 12 I=1,IG                                                    
          R1C1(I,1)=0.0D0                                               
          R1C1(I,2)=0.0D0                                               
          R2C1(I,1)=0.0D0                                               
          R2C1(I,2)=0.0D0                                               
          R3C1(I,1)=0.0D0                                               
          R3C1(I,2)=0.0D0                                               
          R4C1(I,1)=0.0D0                                               
          R4C1(I,2)=0.0D0                                               
   12   CONTINUE                                                        
C                                                                       
        DO 13 I=1,IG                                                    
          R1D1(I,1)=0.0D0                                               
          R1D1(I,2)=0.0D0                                               
          R2D1(I,1)=0.0D0                                               
          R2D1(I,2)=0.0D0                                               
          R3D1(I,1)=0.0D0                                               
          R3D1(I,2)=0.0D0                                               
          R4D1(I,1)=0.0D0                                               
          R4D1(I,2)=0.0D0                                               
   13   CONTINUE                                                        
C                                                                       
        DO 14 I=1,IG                                                    
          R1E1(I,1)=0.0D0                                               
          R1E1(I,2)=0.0D0                                               
          R2E1(I,1)=0.0D0                                               
          R2E1(I,2)=0.0D0                                               
          R3E1(I,1)=0.0D0                                               
          R3E1(I,2)=0.0D0                                               
          R4E1(I,1)=0.0D0                                               
          R4E1(I,2)=0.0D0                                               
   14   CONTINUE                                                        
C                                                                       
        DO 15 I=1,IG                                                    
          R1F1(I,1)=0.0D0                                               
          R1F1(I,2)=0.0D0                                               
          R2F1(I,1)=0.0D0                                               
          R2F1(I,2)=0.0D0                                               
          R3F1(I,1)=0.0D0                                               
          R3F1(I,2)=0.0D0                                               
          R4F1(I,1)=0.0D0                                               
          R4F1(I,2)=0.0D0                                               
   15   CONTINUE                                                        
C                                                                       
        DO 16 I=1,IG                                                    
          R1H1(I,1)=0.0D0                                               
          R1H1(I,2)=0.0D0                                               
          R2H1(I,1)=0.0D0                                               
          R2H1(I,2)=0.0D0                                               
          R3H1(I,1)=0.0D0                                               
          R3H1(I,2)=0.0D0                                               
          R4H1(I,1)=0.0D0                                               
          R4H1(I,2)=0.0D0                                               
   16   CONTINUE                                                        
C                                                                       
        DO 17 I=1,IG                                                    
          R1I1(I,1)=0.0D0                                               
          R1I1(I,2)=0.0D0                                               
          R2I1(I,1)=0.0D0                                               
          R2I1(I,2)=0.0D0                                               
          R3I1(I,1)=0.0D0                                               
          R3I1(I,2)=0.0D0                                               
          R4I1(I,1)=0.0D0                                               
          R4I1(I,2)=0.0D0                                               
   17   CONTINUE                                                        
C                                                                       
        DO 18 I=1,IG                                                    
          R1J1(I,1)=0.0D0                                               
          R1J1(I,2)=0.0D0                                               
          R2J1(I,1)=0.0D0                                               
          R2J1(I,2)=0.0D0                                               
          R3J1(I,1)=0.0D0                                               
          R3J1(I,2)=0.0D0                                               
          R4J1(I,1)=0.0D0                                               
          R4J1(I,2)=0.0D0                                               
   18   CONTINUE                                                        
C                                                                       
        DO 19 I=1,IG                                                    
          R1K1(I,1)=0.0D0                                               
          R1K1(I,2)=0.0D0                                               
          R2K1(I,1)=0.0D0                                               
          R2K1(I,2)=0.0D0                                               
          R3K1(I,1)=0.0D0                                               
          R3K1(I,2)=0.0D0                                               
          R4K1(I,1)=0.0D0                                               
          R4K1(I,2)=0.0D0                                               
   19   CONTINUE                                                        
C                                                                       
        DO 20 I=1,IG                                                    
          R1L1(I,1)=0.0D0                                               
          R1L1(I,2)=0.0D0                                               
          R2L1(I,1)=0.0D0                                               
          R2L1(I,2)=0.0D0                                               
          R3L1(I,1)=0.0D0                                               
          R3L1(I,2)=0.0D0                                               
          R4L1(I,1)=0.0D0                                               
          R4L1(I,2)=0.0D0                                               
   20   CONTINUE                                                        
C                                                                       
C----------------------------------- SPANWISE                           
C                                                                       
        DO 30 K=1,KG                                                    
          R1A3(K,1)=0.0D0                                               
          R1A3(K,2)=0.0D0                                               
          R2A3(K,1)=0.0D0                                               
          R2A3(K,2)=0.0D0                                               
          R3A3(K,1)=0.0D0                                               
          R3A3(K,2)=0.0D0                                               
          R4A3(K,1)=0.0D0                                               
          R4A3(K,2)=0.0D0                                               
   30   CONTINUE                                                        
C                                                                       
        DO 31 K=1,KG                                                    
          R1B3(K,1)=0.0D0                                               
          R1B3(K,2)=0.0D0                                               
          R2B3(K,1)=0.0D0                                               
          R2B3(K,2)=0.0D0                                               
          R3B3(K,1)=0.0D0                                               
          R3B3(K,2)=0.0D0                                               
          R4B3(K,1)=0.0D0                                               
          R4B3(K,2)=0.0D0                                               
   31   CONTINUE                                                        
C                                                                       
        DO 32 K=1,KG                                                    
          R1C3(K,1)=0.0D0                                               
          R1C3(K,2)=0.0D0                                               
          R2C3(K,1)=0.0D0                                               
          R2C3(K,2)=0.0D0                                               
          R3C3(K,1)=0.0D0                                               
          R3C3(K,2)=0.0D0                                               
          R4C3(K,1)=0.0D0                                               
          R4C3(K,2)=0.0D0                                               
   32   CONTINUE                                                        
C                                                                       
        DO 33 K=1,KG                                                    
          R1D3(K,1)=0.0D0                                               
          R1D3(K,2)=0.0D0                                               
          R2D3(K,1)=0.0D0                                               
          R2D3(K,2)=0.0D0                                               
          R3D3(K,1)=0.0D0                                               
          R3D3(K,2)=0.0D0                                               
          R4D3(K,1)=0.0D0                                               
          R4D3(K,2)=0.0D0                                               
   33   CONTINUE                                                        
C                                                                       
        DO 34 K=1,KG                                                    
          R1E3(K,1)=0.0D0                                               
          R1E3(K,2)=0.0D0                                               
          R2E3(K,1)=0.0D0                                               
          R2E3(K,2)=0.0D0                                               
          R3E3(K,1)=0.0D0                                               
          R3E3(K,2)=0.0D0                                               
          R4E3(K,1)=0.0D0                                               
          R4E3(K,2)=0.0D0                                               
   34   CONTINUE                                                        
C                                                                       
        DO 35 K=1,KG                                                    
          R1F3(K,1)=0.0D0                                               
          R1F3(K,2)=0.0D0                                               
          R2F3(K,1)=0.0D0                                               
          R2F3(K,2)=0.0D0                                               
          R3F3(K,1)=0.0D0                                               
          R3F3(K,2)=0.0D0                                               
          R4F3(K,1)=0.0D0                                               
          R4F3(K,2)=0.0D0                                               
   35   CONTINUE                                                        
C                                                                       
        DO 36 K=1,KG                                                    
          R1H3(K,1)=0.0D0                                               
          R1H3(K,2)=0.0D0                                               
          R2H3(K,1)=0.0D0                                               
          R2H3(K,2)=0.0D0                                               
          R3H3(K,1)=0.0D0                                               
          R3H3(K,2)=0.0D0                                               
          R4H3(K,1)=0.0D0                                               
          R4H3(K,2)=0.0D0                                               
   36   CONTINUE                                                        
C                                                                       
        DO 37 K=1,KG                                                    
          R1I3(K,1)=0.0D0                                               
          R1I3(K,2)=0.0D0                                               
          R2I3(K,1)=0.0D0                                               
          R2I3(K,2)=0.0D0                                               
          R3I3(K,1)=0.0D0                                               
          R3I3(K,2)=0.0D0                                               
          R4I3(K,1)=0.0D0                                               
          R4I3(K,2)=0.0D0                                               
   37   CONTINUE                                                        
C                                                                       
        DO 38 K=1,KG                                                    
          R1J3(K,1)=0.0D0                                               
          R1J3(K,2)=0.0D0                                               
          R2J3(K,1)=0.0D0                                               
          R2J3(K,2)=0.0D0                                               
          R3J3(K,1)=0.0D0                                               
          R3J3(K,2)=0.0D0                                               
          R4J3(K,1)=0.0D0                                               
          R4J3(K,2)=0.0D0                                               
   38   CONTINUE                                                        
C                                                                       
        DO 39 K=1,KG                                                    
          R1K3(K,1)=0.0D0                                               
          R1K3(K,2)=0.0D0                                               
          R2K3(K,1)=0.0D0                                               
          R2K3(K,2)=0.0D0                                               
          R3K3(K,1)=0.0D0                                               
          R3K3(K,2)=0.0D0                                               
          R4K3(K,1)=0.0D0                                               
          R4K3(K,2)=0.0D0                                               
   39   CONTINUE                                                        
C                                                                       
        DO 40 K=1,KG                                                    
          R1L3(K,1)=0.0D0                                               
          R1L3(K,2)=0.0D0                                               
          R2L3(K,1)=0.0D0                                               
          R2L3(K,2)=0.0D0                                               
          R3L3(K,1)=0.0D0                                               
          R3L3(K,2)=0.0D0                                               
          R4L3(K,1)=0.0D0                                               
          R4L3(K,2)=0.0D0                                               
   40   CONTINUE                                                        
C                                                                       
C----------------------------------- STREAMWISE                         
C                                                                       
!$OMP PARALLEL DO PRIVATE(K)                                            
        DO 50 I=1,IG                                                    
        DO 50 K=1,KG                                                    
         R1A1(I,1)=R1A1(I,1)+UT(1,JA,K,1)*UT(I,JA,K,1)                  
         R1A1(I,2)=R1A1(I,2)+UT(1,JA,K,1)**2                            
         R2A1(I,1)=R2A1(I,1)+(UT(1,JA,K,2)+UT(1,JA-1,K,2))              
     &                     *(UT(I,JA,K,2)+UT(I,JA-1,K,2))*0.25D0        
         R2A1(I,2)=R2A1(I,2)+((UT(1,JA,K,2)+UT(1,JA-1,K,2))*0.5D0)**2   
         R3A1(I,1)=R3A1(I,1)+UT(1,JA,K,3)*UT(I,JA,K,3)                  
         R3A1(I,2)=R3A1(I,2)+UT(1,JA,K,3)**2                            
         R4A1(I,1)=R4A1(I,1)+(P(1,JA,K)-AP(JA))*(P(I,JA,K)-AP(JA))      
         R4A1(I,2)=R4A1(I,2)+(P(1,JA,K)-AP(JA))**2                      
   50   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(K)                                            
        DO 51 I=1,IG                                                    
        DO 51 K=1,KG                                                    
         R1B1(I,1)=R1B1(I,1)+UT(1,JB,K,1)*UT(I,JB,K,1)                  
         R1B1(I,2)=R1B1(I,2)+UT(1,JB,K,1)**2                            
         R2B1(I,1)=R2B1(I,1)+(UT(1,JB,K,2)+UT(1,JB-1,K,2))              
     &                       *(UT(I,JB,K,2)+UT(I,JB-1,K,2))*0.25D0      
         R2B1(I,2)=R2B1(I,2)+((UT(1,JB,K,2)+UT(1,JB-1,K,2))*0.5D0)**2   
         R3B1(I,1)=R3B1(I,1)+UT(1,JB,K,3)*UT(I,JB,K,3)                  
         R3B1(I,2)=R3B1(I,2)+UT(1,JB,K,3)**2                            
         R4B1(I,1)=R4B1(I,1)+(P(1,JB,K)-AP(JB))*(P(I,JB,K)-AP(JB))      
         R4B1(I,2)=R4B1(I,2)+(P(1,JB,K)-AP(JB))**2                      
   51   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(K)                                            
        DO 52 I=1,IG                                                    
        DO 52 K=1,KG                                                    
         R1C1(I,1)=R1C1(I,1)+UT(1,JC,K,1)*UT(I,JC,K,1)                  
         R1C1(I,2)=R1C1(I,2)+UT(1,JC,K,1)**2                            
         R2C1(I,1)=R2C1(I,1)+(UT(1,JC,K,2)+UT(1,JC-1,K,2))              
     &                       *(UT(I,JC,K,2)+UT(I,JC-1,K,2))*0.25D0      
         R2C1(I,2)=R2C1(I,2)+((UT(1,JC,K,2)+UT(1,JC-1,K,2))*0.5D0)**2   
         R3C1(I,1)=R3C1(I,1)+UT(1,JC,K,3)*UT(I,JC,K,3)                  
         R3C1(I,2)=R3C1(I,2)+UT(1,JC,K,3)**2                            
         R4C1(I,1)=R4C1(I,1)+(P(1,JC,K)-AP(JC))*(P(I,JC,K)-AP(JC))      
         R4C1(I,2)=R4C1(I,2)+(P(1,JC,K)-AP(JC))**2                      
   52   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(K)                                            
        DO 53 I=1,IG                                                    
        DO 53 K=1,KG                                                    
         R1D1(I,1)=R1D1(I,1)+UT(1,JD,K,1)*UT(I,JD,K,1)                  
         R1D1(I,2)=R1D1(I,2)+UT(1,JD,K,1)**2                            
         R2D1(I,1)=R2D1(I,1)+(UT(1,JD,K,2)+UT(1,JD-1,K,2))              
     &                       *(UT(I,JD,K,2)+UT(I,JD-1,K,2))*0.25D0      
         R2D1(I,2)=R2D1(I,2)+((UT(1,JD,K,2)+UT(1,JD-1,K,2))*0.5D0)**2   
         R3D1(I,1)=R3D1(I,1)+UT(1,JD,K,3)*UT(I,JD,K,3)                  
         R3D1(I,2)=R3D1(I,2)+UT(1,JD,K,3)**2                            
         R4D1(I,1)=R4D1(I,1)+(P(1,JD,K)-AP(JD))*(P(I,JD,K)-AP(JD))      
         R4D1(I,2)=R4D1(I,2)+(P(1,JD,K)-AP(JD))**2                      
   53   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(K)                                            
        DO 54 I=1,IG                                                    
        DO 54 K=1,KG                                                    
         R1E1(I,1)=R1E1(I,1)+UT(1,JE,K,1)*UT(I,JE,K,1)                  
         R1E1(I,2)=R1E1(I,2)+UT(1,JE,K,1)**2                            
         R2E1(I,1)=R2E1(I,1)+(UT(1,JE,K,2)+UT(1,JE-1,K,2))              
     &                       *(UT(I,JE,K,2)+UT(I,JE-1,K,2))*0.25D0      
         R2E1(I,2)=R2E1(I,2)+((UT(1,JE,K,2)+UT(1,JE-1,K,2))*0.5D0)**2   
         R3E1(I,1)=R3E1(I,1)+UT(1,JE,K,3)*UT(I,JE,K,3)                  
         R3E1(I,2)=R3E1(I,2)+UT(1,JE,K,3)**2                            
         R4E1(I,1)=R4E1(I,1)+(P(1,JE,K)-AP(JE))*(P(I,JE,K)-AP(JE))      
         R4E1(I,2)=R4E1(I,2)+(P(1,JE,K)-AP(JE))**2                      
   54   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(K)                                            
        DO 55 I=1,IG                                                    
        DO 55 K=1,KG                                                    
         R1F1(I,1)=R1F1(I,1)+UT(1,JF,K,1)*UT(I,JF,K,1)                  
         R1F1(I,2)=R1F1(I,2)+UT(1,JF,K,1)**2                            
         R2F1(I,1)=R2F1(I,1)+(UT(1,JF,K,2)+UT(1,JF-1,K,2))              
     &                       *(UT(I,JF,K,2)+UT(I,JF-1,K,2))*0.25D0      
         R2F1(I,2)=R2F1(I,2)+((UT(1,JF,K,2)+UT(1,JF-1,K,2))*0.5D0)**2   
         R3F1(I,1)=R3F1(I,1)+UT(1,JF,K,3)*UT(I,JF,K,3)                  
         R3F1(I,2)=R3F1(I,2)+UT(1,JF,K,3)**2                            
         R4F1(I,1)=R4F1(I,1)+(P(1,JF,K)-AP(JF))*(P(I,JF,K)-AP(JF))      
         R4F1(I,2)=R4F1(I,2)+(P(1,JF,K)-AP(JF))**2                      
   55   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(K)                                            
        DO 56 I=1,IG                                                    
        DO 56 K=1,KG                                                    
         R1H1(I,1)=R1H1(I,1)+UT(1,JH,K,1)*UT(I,JH,K,1)                  
         R1H1(I,2)=R1H1(I,2)+UT(1,JH,K,1)**2                            
         R2H1(I,1)=R2H1(I,1)+(UT(1,JH,K,2)+UT(1,JH-1,K,2))              
     &                       *(UT(I,JH,K,2)+UT(I,JH-1,K,2))*0.25D0      
         R2H1(I,2)=R2H1(I,2)+((UT(1,JH,K,2)+UT(1,JH-1,K,2))*0.5D0)**2   
         R3H1(I,1)=R3H1(I,1)+UT(1,JH,K,3)*UT(I,JH,K,3)                  
         R3H1(I,2)=R3H1(I,2)+UT(1,JH,K,3)**2                            
         R4H1(I,1)=R4H1(I,1)+(P(1,JH,K)-AP(JH))*(P(I,JH,K)-AP(JH))      
         R4H1(I,2)=R4H1(I,2)+(P(1,JH,K)-AP(JH))**2                      
   56   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(K)                                            
        DO 57 I=1,IG                                                    
        DO 57 K=1,KG                                                    
         R1I1(I,1)=R1I1(I,1)+UT(1,JI,K,1)*UT(I,JI,K,1)                  
         R1I1(I,2)=R1I1(I,2)+UT(1,JI,K,1)**2                            
         R2I1(I,1)=R2I1(I,1)+(UT(1,JI,K,2)+UT(1,JI-1,K,2))              
     &                       *(UT(I,JI,K,2)+UT(I,JI-1,K,2))*0.25D0      
         R2I1(I,2)=R2I1(I,2)+((UT(1,JI,K,2)+UT(1,JI-1,K,2))*0.5D0)**2   
         R3I1(I,1)=R3I1(I,1)+UT(1,JI,K,3)*UT(I,JI,K,3)                  
         R3I1(I,2)=R3I1(I,2)+UT(1,JI,K,3)**2                            
         R4I1(I,1)=R4I1(I,1)+(P(1,JI,K)-AP(JI))*(P(I,JI,K)-AP(JI))      
         R4I1(I,2)=R4I1(I,2)+(P(1,JI,K)-AP(JI))**2                      
   57   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(K)                                            
        DO 58 I=1,IG                                                    
        DO 58 K=1,KG                                                    
         R1J1(I,1)=R1J1(I,1)+UT(1,JJ,K,1)*UT(I,JJ,K,1)                  
         R1J1(I,2)=R1J1(I,2)+UT(1,JJ,K,1)**2                            
         R2J1(I,1)=R2J1(I,1)+(UT(1,JJ,K,2)+UT(1,JJ-1,K,2))              
     &                       *(UT(I,JJ,K,2)+UT(I,JJ-1,K,2))*0.25D0      
         R2J1(I,2)=R2J1(I,2)+((UT(1,JJ,K,2)+UT(1,JJ-1,K,2))*0.5D0)**2   
         R3J1(I,1)=R3J1(I,1)+UT(1,JJ,K,3)*UT(I,JJ,K,3)                  
         R3J1(I,2)=R3J1(I,2)+UT(1,JJ,K,3)**2                            
         R4J1(I,1)=R4J1(I,1)+(P(1,JJ,K)-AP(JJ))*(P(I,JJ,K)-AP(JJ))      
         R4J1(I,2)=R4J1(I,2)+(P(1,JJ,K)-AP(JJ))**2                      
   58   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(K)                                            
        DO 59 I=1,IG                                                    
        DO 59 K=1,KG                                                    
         R1K1(I,1)=R1K1(I,1)+UT(1,JK,K,1)*UT(I,JK,K,1)                  
         R1K1(I,2)=R1K1(I,2)+UT(1,JK,K,1)**2                            
         R2K1(I,1)=R2K1(I,1)+(UT(1,JK,K,2)+UT(1,JK-1,K,2))              
     &                       *(UT(I,JK,K,2)+UT(I,JK-1,K,2))*0.25D0      
         R2K1(I,2)=R2K1(I,2)+((UT(1,JK,K,2)+UT(1,JK-1,K,2))*0.5D0)**2   
         R3K1(I,1)=R3K1(I,1)+UT(1,JK,K,3)*UT(I,JK,K,3)                  
         R3K1(I,2)=R3K1(I,2)+UT(1,JK,K,3)**2                            
         R4K1(I,1)=R4K1(I,1)+(P(1,JK,K)-AP(JK))*(P(I,JK,K)-AP(JK))      
         R4K1(I,2)=R4K1(I,2)+(P(1,JK,K)-AP(JK))**2                      
   59   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(K)                                            
        DO 60 I=1,IG                                                    
        DO 60 K=1,KG                                                    
         R1L1(I,1)=R1L1(I,1)+UT(1,JL,K,1)*UT(I,JL,K,1)                  
         R1L1(I,2)=R1L1(I,2)+UT(1,JL,K,1)**2                            
         R2L1(I,1)=R2L1(I,1)+(UT(1,JL,K,2)+UT(1,JL-1,K,2))              
     &                       *(UT(I,JL,K,2)+UT(I,JL-1,K,2))*0.25D0      
         R2L1(I,2)=R2L1(I,2)+((UT(1,JL,K,2)+UT(1,JL-1,K,2))*0.5D0)**2   
         R3L1(I,1)=R3L1(I,1)+UT(1,JL,K,3)*UT(I,JL,K,3)                  
         R3L1(I,2)=R3L1(I,2)+UT(1,JL,K,3)**2                            
         R4L1(I,1)=R4L1(I,1)+(P(1,JL,K)-AP(JL))*(P(I,JL,K)-AP(JL))      
         R4L1(I,2)=R4L1(I,2)+(P(1,JL,K)-AP(JL))**2                      
   60   CONTINUE                                                        
C                                                                       
C----------------------------------- SPANWISE                           
C                                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
        DO 70 K=1,KG                                                    
        DO 70 I=1,IG                                                    
         R1A3(K,1)=R1A3(K,1)+UT(I,JA,1,1)*UT(I,JA,K,1)                  
         R1A3(K,2)=R1A3(K,2)+UT(I,JA,1,1)**2                            
         R2A3(K,1)=R2A3(K,1)+(UT(I,JA,1,2)+UT(I,JA-1,1,2))              
     &                       *(UT(I,JA,K,2)+UT(I,JA-1,K,2))*0.25D0      
         R2A3(K,2)=R2A3(K,2)+((UT(I,JA,1,2)+UT(I,JA-1,1,2))*0.5D0)**2   
         R3A3(K,1)=R3A3(K,1)+UT(I,JA,1,3)*UT(I,JA,K,3)                  
         R3A3(K,2)=R3A3(K,2)+UT(I,JA,1,3)**2                            
         R4A3(K,1)=R4A3(K,1)+(P(I,JA,1)-AP(JA))*(P(I,JA,K)-AP(JA))      
         R4A3(K,2)=R4A3(K,2)+(P(I,JA,1)-AP(JA))**2                      
   70   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
        DO 71 K=1,KG                                                    
        DO 71 I=1,IG                                                    
         R1B3(K,1)=R1B3(K,1)+UT(I,JB,1,1)*UT(I,JB,K,1)                  
         R1B3(K,2)=R1B3(K,2)+UT(I,JB,1,1)**2                            
         R2B3(K,1)=R2B3(K,1)+(UT(I,JB,1,2)+UT(I,JB-1,1,2))              
     &                       *(UT(I,JB,K,2)+UT(I,JB-1,K,2))*0.25D0      
         R2B3(K,2)=R2B3(K,2)+((UT(I,JB,1,2)+UT(I,JB-1,1,2))*0.5D0)**2   
         R3B3(K,1)=R3B3(K,1)+UT(I,JB,1,3)*UT(I,JB,K,3)                  
         R3B3(K,2)=R3B3(K,2)+UT(I,JB,1,3)**2                            
         R4B3(K,1)=R4B3(K,1)+(P(I,JB,1)-AP(JB))*(P(I,JB,K)-AP(JB))      
         R4B3(K,2)=R4B3(K,2)+(P(I,JB,1)-AP(JB))**2                      
   71   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
        DO 72 K=1,KG                                                    
        DO 72 I=1,IG                                                    
         R1C3(K,1)=R1C3(K,1)+UT(I,JC,1,1)*UT(I,JC,K,1)                  
         R1C3(K,2)=R1C3(K,2)+UT(I,JC,1,1)**2                            
         R2C3(K,1)=R2C3(K,1)+(UT(I,JC,1,2)+UT(I,JC-1,1,2))              
     &                       *(UT(I,JC,K,2)+UT(I,JC-1,K,2))*0.25D0      
         R2C3(K,2)=R2C3(K,2)+((UT(I,JC,1,2)+UT(I,JC-1,1,2))*0.5D0)**2   
         R3C3(K,1)=R3C3(K,1)+UT(I,JC,1,3)*UT(I,JC,K,3)                  
         R3C3(K,2)=R3C3(K,2)+UT(I,JC,1,3)**2                            
         R4C3(K,1)=R4C3(K,1)+(P(I,JC,1)-AP(JC))*(P(I,JC,K)-AP(JC))      
         R4C3(K,2)=R4C3(K,2)+(P(I,JC,1)-AP(JC))**2                      
   72   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
        DO 73 K=1,KG                                                    
        DO 73 I=1,IG                                                    
         R1D3(K,1)=R1D3(K,1)+UT(I,JD,1,1)*UT(I,JD,K,1)                  
         R1D3(K,2)=R1D3(K,2)+UT(I,JD,1,1)**2                            
         R2D3(K,1)=R2D3(K,1)+(UT(I,JD,1,2)+UT(I,JD-1,1,2))              
     &                       *(UT(I,JD,K,2)+UT(I,JD-1,K,2))*0.25D0      
         R2D3(K,2)=R2D3(K,2)+((UT(I,JD,1,2)+UT(I,JD-1,1,2))*0.5D0)**2   
         R3D3(K,1)=R3D3(K,1)+UT(I,JD,1,3)*UT(I,JD,K,3)                  
         R3D3(K,2)=R3D3(K,2)+UT(I,JD,1,3)**2                            
         R4D3(K,1)=R4D3(K,1)+(P(I,JD,1)-AP(JD))*(P(I,JD,K)-AP(JD))      
         R4D3(K,2)=R4D3(K,2)+(P(I,JD,1)-AP(JD))**2                      
   73   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
        DO 74 K=1,KG                                                    
        DO 74 I=1,IG                                                    
         R1E3(K,1)=R1E3(K,1)+UT(I,JE,1,1)*UT(I,JE,K,1)                  
         R1E3(K,2)=R1E3(K,2)+UT(I,JE,1,1)**2                            
         R2E3(K,1)=R2E3(K,1)+(UT(I,JE,1,2)+UT(I,JE-1,1,2))              
     &                       *(UT(I,JE,K,2)+UT(I,JE-1,K,2))*0.25D0      
         R2E3(K,2)=R2E3(K,2)+((UT(I,JE,1,2)+UT(I,JE-1,1,2))*0.5D0)**2   
         R3E3(K,1)=R3E3(K,1)+UT(I,JE,1,3)*UT(I,JE,K,3)                  
         R3E3(K,2)=R3E3(K,2)+UT(I,JE,1,3)**2                            
         R4E3(K,1)=R4E3(K,1)+(P(I,JE,1)-AP(JE))*(P(I,JE,K)-AP(JE))      
         R4E3(K,2)=R4E3(K,2)+(P(I,JE,1)-AP(JE))**2                      
   74   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
        DO 75 K=1,KG                                                    
        DO 75 I=1,IG                                                    
         R1F3(K,1)=R1F3(K,1)+UT(I,JF,1,1)*UT(I,JF,K,1)                  
         R1F3(K,2)=R1F3(K,2)+UT(I,JF,1,1)**2                            
         R2F3(K,1)=R2F3(K,1)+(UT(I,JF,1,2)+UT(I,JF-1,1,2))              
     &                       *(UT(I,JF,K,2)+UT(I,JF-1,K,2))*0.25D0      
         R2F3(K,2)=R2F3(K,2)+((UT(I,JF,1,2)+UT(I,JF-1,1,2))*0.5D0)**2   
         R3F3(K,1)=R3F3(K,1)+UT(I,JF,1,3)*UT(I,JF,K,3)                  
         R3F3(K,2)=R3F3(K,2)+UT(I,JF,1,3)**2                            
         R4F3(K,1)=R4F3(K,1)+(P(I,JF,1)-AP(JF))*(P(I,JF,K)-AP(JF))      
         R4F3(K,2)=R4F3(K,2)+(P(I,JF,1)-AP(JF))**2                      
   75   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
        DO 76 K=1,KG                                                    
        DO 76 I=1,IG                                                    
         R1H3(K,1)=R1H3(K,1)+UT(I,JH,1,1)*UT(I,JH,K,1)                  
         R1H3(K,2)=R1H3(K,2)+UT(I,JH,1,1)**2                            
         R2H3(K,1)=R2H3(K,1)+(UT(I,JH,1,2)+UT(I,JH-1,1,2))              
     &                       *(UT(I,JH,K,2)+UT(I,JH-1,K,2))*0.25D0      
         R2H3(K,2)=R2H3(K,2)+((UT(I,JH,1,2)+UT(I,JH-1,1,2))*0.5D0)**2   
         R3H3(K,1)=R3H3(K,1)+UT(I,JH,1,3)*UT(I,JH,K,3)                  
         R3H3(K,2)=R3H3(K,2)+UT(I,JH,1,3)**2                            
         R4H3(K,1)=R4H3(K,1)+(P(I,JH,1)-AP(JH))*(P(I,JH,K)-AP(JH))      
         R4H3(K,2)=R4H3(K,2)+(P(I,JH,1)-AP(JH))**2                      
   76   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
        DO 77 K=1,KG                                                    
        DO 77 I=1,IG                                                    
         R1I3(K,1)=R1I3(K,1)+UT(I,JI,1,1)*UT(I,JI,K,1)                  
         R1I3(K,2)=R1I3(K,2)+UT(I,JI,1,1)**2                            
         R2I3(K,1)=R2I3(K,1)+(UT(I,JI,1,2)+UT(I,JI-1,1,2))              
     &                       *(UT(I,JI,K,2)+UT(I,JI-1,K,2))*0.25D0      
         R2I3(K,2)=R2I3(K,2)+((UT(I,JI,1,2)+UT(I,JI-1,1,2))*0.5D0)**2   
         R3I3(K,1)=R3I3(K,1)+UT(I,JI,1,3)*UT(I,JI,K,3)                  
         R3I3(K,2)=R3I3(K,2)+UT(I,JI,1,3)**2                            
         R4I3(K,1)=R4I3(K,1)+(P(I,JI,1)-AP(JI))*(P(I,JI,K)-AP(JI))      
         R4I3(K,2)=R4I3(K,2)+(P(I,JI,1)-AP(JI))**2                      
   77   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
        DO 78 K=1,KG                                                    
        DO 78 I=1,IG                                                    
         R1J3(K,1)=R1J3(K,1)+UT(I,JJ,1,1)*UT(I,JJ,K,1)                  
         R1J3(K,2)=R1J3(K,2)+UT(I,JJ,1,1)**2                            
         R2J3(K,1)=R2J3(K,1)+(UT(I,JJ,1,2)+UT(I,JJ-1,1,2))              
     &                       *(UT(I,JJ,K,2)+UT(I,JJ-1,K,2))*0.25D0      
         R2J3(K,2)=R2J3(K,2)+((UT(I,JJ,1,2)+UT(I,JJ-1,1,2))*0.5D0)**2   
         R3J3(K,1)=R3J3(K,1)+UT(I,JJ,1,3)*UT(I,JJ,K,3)                  
         R3J3(K,2)=R3J3(K,2)+UT(I,JJ,1,3)**2                            
         R4J3(K,1)=R4J3(K,1)+(P(I,JJ,1)-AP(JJ))*(P(I,JJ,K)-AP(JJ))      
         R4J3(K,2)=R4J3(K,2)+(P(I,JJ,1)-AP(JJ))**2                      
   78   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
        DO 79 K=1,KG                                                    
        DO 79 I=1,IG                                                    
         R1K3(K,1)=R1K3(K,1)+UT(I,JK,1,1)*UT(I,JK,K,1)                  
         R1K3(K,2)=R1K3(K,2)+UT(I,JK,1,1)**2                            
         R2K3(K,1)=R2K3(K,1)+(UT(I,JK,1,2)+UT(I,JK-1,1,2))              
     &                       *(UT(I,JK,K,2)+UT(I,JK-1,K,2))*0.25D0      
         R2K3(K,2)=R2K3(K,2)+((UT(I,JK,1,2)+UT(I,JK-1,1,2))*0.5D0)**2   
         R3K3(K,1)=R3K3(K,1)+UT(I,JK,1,3)*UT(I,JK,K,3)                  
         R3K3(K,2)=R3K3(K,2)+UT(I,JK,1,3)**2                            
         R4K3(K,1)=R4K3(K,1)+(P(I,JK,1)-AP(JK))*(P(I,JK,K)-AP(JK))      
         R4K3(K,2)=R4K3(K,2)+(P(I,JK,1)-AP(JK))**2                      
   79   CONTINUE                                                        
C                                                                       
!$OMP PARALLEL DO PRIVATE(I)                                            
        DO 80 K=1,KG                                                    
        DO 80 I=1,IG                                                    
         R1L3(K,1)=R1L3(K,1)+UT(I,JL,1,1)*UT(I,JL,K,1)                  
         R1L3(K,2)=R1L3(K,2)+UT(I,JL,1,1)**2                            
         R2L3(K,1)=R2L3(K,1)+(UT(I,JL,1,2)+UT(I,JL-1,1,2))              
     &                       *(UT(I,JL,K,2)+UT(I,JL-1,K,2))*0.25D0      
         R2L3(K,2)=R2L3(K,2)+((UT(I,JL,1,2)+UT(I,JL-1,1,2))*0.5D0)**2   
         R3L3(K,1)=R3L3(K,1)+UT(I,JL,1,3)*UT(I,JL,K,3)                  
         R3L3(K,2)=R3L3(K,2)+UT(I,JL,1,3)**2                            
         R4L3(K,1)=R4L3(K,1)+(P(I,JL,1)-AP(JL))*(P(I,JL,K)-AP(JL))      
         R4L3(K,2)=R4L3(K,2)+(P(I,JL,1)-AP(JL))**2                      
   80   CONTINUE                                                        
C                                                                       
C----------------------------------- STREAMWISE                         
C                                                                       
        DO 90 I=1,IG                                                    
         ATPC1(I,JA,1)=ATPC1(I,JA,1)+R1A1(I,1)/R1A1(I,2)                
         ATPC1(I,JA,2)=ATPC1(I,JA,2)+R2A1(I,1)/R2A1(I,2)                
         ATPC1(I,JA,3)=ATPC1(I,JA,3)+R3A1(I,1)/R3A1(I,2)                
         ATPC1(I,JA,4)=ATPC1(I,JA,4)+R4A1(I,1)/R4A1(I,2)                
   90   CONTINUE                                                        
C                                                                       
        DO 91 I=1,IG                                                    
         ATPC1(I,JB,1)=ATPC1(I,JB,1)+R1B1(I,1)/R1B1(I,2)                
         ATPC1(I,JB,2)=ATPC1(I,JB,2)+R2B1(I,1)/R2B1(I,2)                
         ATPC1(I,JB,3)=ATPC1(I,JB,3)+R3B1(I,1)/R3B1(I,2)                
         ATPC1(I,JB,4)=ATPC1(I,JB,4)+R4B1(I,1)/R4B1(I,2)                
   91   CONTINUE                                                        
C                                                                       
        DO 92 I=1,IG                                                    
         ATPC1(I,JC,1)=ATPC1(I,JC,1)+R1C1(I,1)/R1C1(I,2)                
         ATPC1(I,JC,2)=ATPC1(I,JC,2)+R2C1(I,1)/R2C1(I,2)                
         ATPC1(I,JC,3)=ATPC1(I,JC,3)+R3C1(I,1)/R3C1(I,2)                
         ATPC1(I,JC,4)=ATPC1(I,JC,4)+R4C1(I,1)/R4C1(I,2)                
   92   CONTINUE                                                        
C                                                                       
        DO 93 I=1,IG                                                    
         ATPC1(I,JD,1)=ATPC1(I,JD,1)+R1D1(I,1)/R1D1(I,2)                
         ATPC1(I,JD,2)=ATPC1(I,JD,2)+R2D1(I,1)/R2D1(I,2)                
         ATPC1(I,JD,3)=ATPC1(I,JD,3)+R3D1(I,1)/R3D1(I,2)                
         ATPC1(I,JD,4)=ATPC1(I,JD,4)+R4D1(I,1)/R4D1(I,2)                
   93   CONTINUE                                                        
C                                                                       
        DO 94 I=1,IG                                                    
         ATPC1(I,JE,1)=ATPC1(I,JE,1)+R1E1(I,1)/R1E1(I,2)                
         ATPC1(I,JE,2)=ATPC1(I,JE,2)+R2E1(I,1)/R2E1(I,2)                
         ATPC1(I,JE,3)=ATPC1(I,JE,3)+R3E1(I,1)/R3E1(I,2)                
         ATPC1(I,JE,4)=ATPC1(I,JE,4)+R4E1(I,1)/R4E1(I,2)                
   94   CONTINUE                                                        
C                                                                       
        DO 95 I=1,IG                                                    
         ATPC1(I,JF,1)=ATPC1(I,JF,1)+R1F1(I,1)/R1F1(I,2)                
         ATPC1(I,JF,2)=ATPC1(I,JF,2)+R2F1(I,1)/R2F1(I,2)                
         ATPC1(I,JF,3)=ATPC1(I,JF,3)+R3F1(I,1)/R3F1(I,2)                
         ATPC1(I,JF,4)=ATPC1(I,JF,4)+R4F1(I,1)/R4F1(I,2)                
   95   CONTINUE                                                        
C                                                                       
        DO 96 I=1,IG                                                    
         ATPC1(I,JH,1)=ATPC1(I,JH,1)+R1H1(I,1)/R1H1(I,2)                
         ATPC1(I,JH,2)=ATPC1(I,JH,2)+R2H1(I,1)/R2H1(I,2)                
         ATPC1(I,JH,3)=ATPC1(I,JH,3)+R3H1(I,1)/R3H1(I,2)                
         ATPC1(I,JH,4)=ATPC1(I,JH,4)+R4H1(I,1)/R4H1(I,2)                
   96   CONTINUE                                                        
C                                                                       
        DO 97 I=1,IG                                                    
         ATPC1(I,JI,1)=ATPC1(I,JI,1)+R1I1(I,1)/R1I1(I,2)                
         ATPC1(I,JI,2)=ATPC1(I,JI,2)+R2I1(I,1)/R2I1(I,2)                
         ATPC1(I,JI,3)=ATPC1(I,JI,3)+R3I1(I,1)/R3I1(I,2)                
         ATPC1(I,JI,4)=ATPC1(I,JI,4)+R4I1(I,1)/R4I1(I,2)                
   97   CONTINUE                                                        
C                                                                       
        DO 98 I=1,IG                                                    
         ATPC1(I,JJ,1)=ATPC1(I,JJ,1)+R1J1(I,1)/R1J1(I,2)                
         ATPC1(I,JJ,2)=ATPC1(I,JJ,2)+R2J1(I,1)/R2J1(I,2)                
         ATPC1(I,JJ,3)=ATPC1(I,JJ,3)+R3J1(I,1)/R3J1(I,2)                
         ATPC1(I,JJ,4)=ATPC1(I,JJ,4)+R4J1(I,1)/R4J1(I,2)                
   98   CONTINUE                                                        
C                                                                       
        DO 99 I=1,IG                                                    
         ATPC1(I,JK,1)=ATPC1(I,JK,1)+R1K1(I,1)/R1K1(I,2)                
         ATPC1(I,JK,2)=ATPC1(I,JK,2)+R2K1(I,1)/R2K1(I,2)                
         ATPC1(I,JK,3)=ATPC1(I,JK,3)+R3K1(I,1)/R3K1(I,2)                
         ATPC1(I,JK,4)=ATPC1(I,JK,4)+R4K1(I,1)/R4K1(I,2)                
   99   CONTINUE                                                        
C                                                                       
        DO 100 I=1,IG                                                   
         ATPC1(I,JL,1)=ATPC1(I,JL,1)+R1L1(I,1)/R1L1(I,2)                
         ATPC1(I,JL,2)=ATPC1(I,JL,2)+R2L1(I,1)/R2L1(I,2)                
         ATPC1(I,JL,3)=ATPC1(I,JL,3)+R3L1(I,1)/R3L1(I,2)                
         ATPC1(I,JL,4)=ATPC1(I,JL,4)+R4L1(I,1)/R4L1(I,2)                
  100   CONTINUE                                                        
C                                                                       
C----------------------------------- SPANWISE                           
C                                                                       
        DO 110 K=1,KG                                                   
         ATPC3(K,JA,1)=ATPC3(K,JA,1)+R1A3(K,1)/R1A3(K,2)                
         ATPC3(K,JA,2)=ATPC3(K,JA,2)+R2A3(K,1)/R2A3(K,2)                
         ATPC3(K,JA,3)=ATPC3(K,JA,3)+R3A3(K,1)/R3A3(K,2)                
         ATPC3(K,JA,4)=ATPC3(K,JA,4)+R4A3(K,1)/R4A3(K,2)                
  110   CONTINUE                                                        
C                                                                       
        DO 111 K=1,KG                                                   
         ATPC3(K,JB,1)=ATPC3(K,JB,1)+R1B3(K,1)/R1B3(K,2)                
         ATPC3(K,JB,2)=ATPC3(K,JB,2)+R2B3(K,1)/R2B3(K,2)                
         ATPC3(K,JB,3)=ATPC3(K,JB,3)+R3B3(K,1)/R3B3(K,2)                
         ATPC3(K,JB,4)=ATPC3(K,JB,4)+R4B3(K,1)/R4B3(K,2)                
  111   CONTINUE                                                        
C                                                                       
        DO 112 K=1,KG                                                   
         ATPC3(K,JC,1)=ATPC3(K,JC,1)+R1C3(K,1)/R1C3(K,2)                
         ATPC3(K,JC,2)=ATPC3(K,JC,2)+R2C3(K,1)/R2C3(K,2)                
         ATPC3(K,JC,3)=ATPC3(K,JC,3)+R3C3(K,1)/R3C3(K,2)                
         ATPC3(K,JC,4)=ATPC3(K,JC,4)+R4C3(K,1)/R4C3(K,2)                
  112   CONTINUE                                                        
C                                                                       
        DO 113 K=1,KG                                                   
         ATPC3(K,JD,1)=ATPC3(K,JD,1)+R1D3(K,1)/R1D3(K,2)                
         ATPC3(K,JD,2)=ATPC3(K,JD,2)+R2D3(K,1)/R2D3(K,2)                
         ATPC3(K,JD,3)=ATPC3(K,JD,3)+R3D3(K,1)/R3D3(K,2)                
         ATPC3(K,JD,4)=ATPC3(K,JD,4)+R4D3(K,1)/R4D3(K,2)                
  113   CONTINUE                                                        
C                                                                       
        DO 114 K=1,KG                                                   
         ATPC3(K,JE,1)=ATPC3(K,JE,1)+R1E3(K,1)/R1E3(K,2)                
         ATPC3(K,JE,2)=ATPC3(K,JE,2)+R2E3(K,1)/R2E3(K,2)                
         ATPC3(K,JE,3)=ATPC3(K,JE,3)+R3E3(K,1)/R3E3(K,2)                
         ATPC3(K,JE,4)=ATPC3(K,JE,4)+R4E3(K,1)/R4E3(K,2)                
  114   CONTINUE                                                        
C                                                                       
        DO 115 K=1,KG                                                   
         ATPC3(K,JF,1)=ATPC3(K,JF,1)+R1F3(K,1)/R1F3(K,2)                
         ATPC3(K,JF,2)=ATPC3(K,JF,2)+R2F3(K,1)/R2F3(K,2)                
         ATPC3(K,JF,3)=ATPC3(K,JF,3)+R3F3(K,1)/R3F3(K,2)                
         ATPC3(K,JF,4)=ATPC3(K,JF,4)+R4F3(K,1)/R4F3(K,2)                
  115   CONTINUE                                                        
C                                                                       
        DO 116 K=1,KG                                                   
         ATPC3(K,JH,1)=ATPC3(K,JH,1)+R1H3(K,1)/R1H3(K,2)                
         ATPC3(K,JH,2)=ATPC3(K,JH,2)+R2H3(K,1)/R2H3(K,2)                
         ATPC3(K,JH,3)=ATPC3(K,JH,3)+R3H3(K,1)/R3H3(K,2)                
         ATPC3(K,JH,4)=ATPC3(K,JH,4)+R4H3(K,1)/R4H3(K,2)                
  116   CONTINUE                                                        
C                                                                       
        DO 117 K=1,KG                                                   
         ATPC3(K,JI,1)=ATPC3(K,JI,1)+R1I3(K,1)/R1I3(K,2)                
         ATPC3(K,JI,2)=ATPC3(K,JI,2)+R2I3(K,1)/R2I3(K,2)                
         ATPC3(K,JI,3)=ATPC3(K,JI,3)+R3I3(K,1)/R3I3(K,2)                
         ATPC3(K,JI,4)=ATPC3(K,JI,4)+R4I3(K,1)/R4I3(K,2)                
  117   CONTINUE                                                        
C                                                                       
        DO 118 K=1,KG                                                   
         ATPC3(K,JJ,1)=ATPC3(K,JJ,1)+R1J3(K,1)/R1J3(K,2)                
         ATPC3(K,JJ,2)=ATPC3(K,JJ,2)+R2J3(K,1)/R2J3(K,2)                
         ATPC3(K,JJ,3)=ATPC3(K,JJ,3)+R3J3(K,1)/R3J3(K,2)                
         ATPC3(K,JJ,4)=ATPC3(K,JJ,4)+R4J3(K,1)/R4J3(K,2)                
  118   CONTINUE                                                        
C                                                                       
        DO 119 K=1,KG                                                   
         ATPC3(K,JK,1)=ATPC3(K,JK,1)+R1K3(K,1)/R1K3(K,2)                
         ATPC3(K,JK,2)=ATPC3(K,JK,2)+R2K3(K,1)/R2K3(K,2)                
         ATPC3(K,JK,3)=ATPC3(K,JK,3)+R3K3(K,1)/R3K3(K,2)                
         ATPC3(K,JK,4)=ATPC3(K,JK,4)+R4K3(K,1)/R4K3(K,2)                
  119   CONTINUE                                                        
C                                                                       
        DO 120 K=1,KG                                                   
         ATPC3(K,JL,1)=ATPC3(K,JL,1)+R1L3(K,1)/R1L3(K,2)                
         ATPC3(K,JL,2)=ATPC3(K,JL,2)+R2L3(K,1)/R2L3(K,2)                
         ATPC3(K,JL,3)=ATPC3(K,JL,3)+R3L3(K,1)/R3L3(K,2)                
         ATPC3(K,JL,4)=ATPC3(K,JL,4)+R4L3(K,1)/R4L3(K,2)                
  120   CONTINUE                                                        
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     VORTICITY (VELOCITY FIELD)                                        
C     PROGRAMED BY H.ABE                                                
      SUBROUTINE VRTCTY                                                 
C     1997.12.29                                                        
C     REVISED BY H.ABE ON 2000.09.24  (OMGU => UT)                      
C====================================================================   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
      INCLUDE './INCFILE/INCMESH'                                         
C                                                                       
C----LOCAL VARIABLES-------------------------------------------         
C                                                                       
      COMMON /DAT1/ RE,DT,DX,DZ                                         
C                                                                       
      DIMENSION U(-4:IG+4,-2:JG+2,-4:KG+4,3)                            
      DIMENSION P(-3:IG+3,-2:JG+2,-3:KG+3)                              
      DIMENSION UT(-4:IG+4,-2:JG+2,-4:KG+4,3)                           
      DIMENSION OMG(JG,3)                                               
      DIMENSION OMG2(JG,3),OMG3(JG,3)                                   
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /VELO/ UG(-4:IG+4,-2:JG+2,-4:KG+4,3)                       
      COMMON /SCAL/ PG(-3:IG+3,-2:JG+2,-3:KG+3)                         
      COMMON /TEMP/ UTG(-4:IG+4,-2:JG+2,-4:KG+4,3)                      
      COMMON /VORT/ OMGG(JG,3)                                          
C                                                                       
C----COMBINATION------------------------------------------------        
C                                                                       
      EQUIVALENCE (UG,U),(PG,P)                                         
      EQUIVALENCE (UTG,UT)                                              
      EQUIVALENCE (OMGG,OMG)                                            
C                                                                       
C=================================================================      
C                                                                       
C                                                                       
C======================================= VORT. FLUC. & AVE.             
C                                                                       
      DO 10 J=1,JG                                                      
        OMG2(J,1)=0.0D0                                                 
        OMG2(J,2)=0.0D0                                                 
        OMG2(J,3)=0.0D0                                                 
        OMG3(J,1)=0.0D0                                                 
        OMG3(J,2)=0.0D0                                                 
        OMG3(J,3)=0.0D0                                                 
   10 CONTINUE                                                          
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 50 J=-2,JG+2                                                   
      DO 50 K=-4,KG+4                                                   
      DO 50 I=-4,IG+4                                                   
        UT(I,J,K,1)=0.0D0                                               
        UT(I,J,K,2)=0.0D0                                               
        UT(I,J,K,3)=0.0D0                                               
   50 CONTINUE                                                          
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 11 J=1,JG                                                      
      DO 11 K=1,KG                                                      
      DO 11 I=1,IG                                                      
        UT(I,J,K,1)=(U(I,J,K+1,2)-U(I,J,K,2))*PDZ                       
     &               -(U(I,J+1,K,3)-U(I,J,K,3))*PDRR(J+1)               
        UT(I,J,K,2)=(U(I,J,K+1,1)-U(I,J,K,1))*PDZ                       
     &               -(U(I+1,J,K,3)-U(I,J,K,3))*PDX                     
        UT(I,J,K,3)=(U(I+1,J,K,2)-U(I,J,K,2))*PDX                       
     &               -(U(I,J+1,K,1)-U(I,J,K,1))*PDRR(J+1)               
   11 CONTINUE                                                          
C                                                                       
      DOOP=1.0/DBLE(IG*KG)                                              
C                                                                       
C-------------------------------------- PLANE AVE.                      
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 20 J=1,JG                                                      
      DO 20 K=1,KG                                                      
      DO 20 I=1,IG                                                      
        OMG2(J,1)=OMG2(J,1)+UT(I,J,K,1)                                 
        OMG2(J,2)=OMG2(J,2)+UT(I,J,K,2)                                 
        OMG2(J,3)=OMG2(J,3)+UT(I,J,K,3)                                 
   20 CONTINUE                                                          
C                                                                       
      DO 21 J=1,JG                                                      
        OMG2(J,1)=OMG2(J,1)*DOOP                                        
        OMG2(J,2)=OMG2(J,2)*DOOP                                        
        OMG2(J,3)=OMG2(J,3)*DOOP                                        
   21 CONTINUE                                                          
C                                                                       
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 30 J=1,JG                                                      
      DO 30 K=1,KG                                                      
      DO 30 I=1,IG                                                      
        OMG3(J,1)=OMG3(J,1)+(UT(I,J,K,1)-OMG2(J,1))**2                  
        OMG3(J,2)=OMG3(J,2)+(UT(I,J,K,2)-OMG2(J,2))**2                  
        OMG3(J,3)=OMG3(J,3)+(UT(I,J,K,3)-OMG2(J,3))**2                  
   30 CONTINUE                                                          
C                                                                       
      DO 31 J=1,JG                                                      
        OMG3(J,1)=OMG3(J,1)*DOOP                                        
        OMG3(J,2)=OMG3(J,2)*DOOP                                        
        OMG3(J,3)=OMG3(J,3)*DOOP                                        
   31 CONTINUE                                                          
C                                                                       
C-------------------------------------- TIME AVE.                       
C                                                                       
      DO 40 J=1,JG                                                      
        OMG(J,1)=OMG(J,1)+OMG3(J,1)                                     
        OMG(J,2)=OMG(J,2)+OMG3(J,2)                                     
        OMG(J,3)=OMG(J,3)+OMG3(J,3)                                     
   40 CONTINUE                                                          
C                                                                       
C-------------------------------------- ZERO-CLR                        
!$OMP PARALLEL DO PRIVATE(I,K)                                          
      DO 60 J=-2,JG+2                                                   
      DO 60 K=-4,KG+4                                                   
      DO 60 I=-4,IG+4                                                   
        UT(I,J,K,1)=0.0D0                                               
        UT(I,J,K,2)=0.0D0                                               
        UT(I,J,K,3)=0.0D0                                               
   60 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     INPUT TIME-AVERAGE DATA (VELOCITY FIELD)                          
C     PROGRAMED BY H.ABE                                                
C====================================================================   
C     INPUT TIME-AVERAGE DATA (VELOCITY FIELD)                          
C     PROGRAMED BY H.ABE                                                
      SUBROUTINE AVEIN                                                  
C     REVISED ON 2001.03.03 BY H.ABE                                    
C     1998.01.03                                                        
C====================================================================   
C                                                                       
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
      INCLUDE './INCFILE/INCDATA'                                         
C      INCLUDE './INCFILE/INCMESH'                                        
C                                                                       
      COMMON /ANUM/ NAVE                                                
      COMMON /ATIM/ AVETIM                                              
C                                                                       
      COMMON /JAVE/ JA,JB,JC,JD,JE,JF,JH,JI,JJ,JK,JL                    
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /UTAU/ AUTUB,AUTUT                                         
      COMMON /UTAU2/ AUTUB2,AUTUT2                                      
      COMMON /AVE1/ AAUG(JG),AAVG(JG),AAWG(JG),AAPG(JG)                 
      COMMON /AVE2/ AURMSG(JG),AVRMSG(JG),AWRMSG(JG),APRMSG(JG)         
      COMMON /AVE3/ ATSSTG(JG),ARESTG(JG)                               
      COMMON /AVE4/ AU11G(JG),AU12G(JG),AU13G(JG)                       
     &             ,AU22G(JG),AU23G(JG),AU33G(JG)                       
      COMMON /SFUV/ SKEWUVG(JG),FLATUVG(JG)                             
      COMMON /BUDG/ AAKG(JG)                                            
     &         ,APKG(JG),AP11G(JG),AP12G(JG)                            
     &         ,APIKG(JG),API11G(JG),API22G(JG),API33G(JG),API12G(JG)   
     &         ,ADKG(JG),AD11G(JG),AD22G(JG),AD33G(JG),AD12G(JG)        
     &         ,AEKG(JG),AE11G(JG),AE22G(JG),AE33G(JG),AE12G(JG)        
     &         ,ATAKG(JG),ATA11G(JG),ATA22G(JG),ATA33G(JG),ATA12G(JG)   
     &         ,APDKG(JG),APD22G(JG),APD12G(JG)                         
     &         ,APS11G(JG),APS22G(JG),APS33G(JG),APS12G(JG)             
      COMMON /POWR/ ASTRMG(IG,JG,5),ASPANG(KG,JG,5)                     
      COMMON /CORR/ ATPC1G(IG,JG,4),ATPC3G(KG,JG,4)                     
      COMMON /VORT/ OMGG(JG,3)                                          
C                                                                       
      COMMON /UHOA/ AUTRIG(JG,8),AUQUADG(JG,10),AQAUVG(JG,4)            
      COMMON /NQAUV/ ANQAUVG(JG,4)                                      
C                                                                       
      COMMON /AAIN/ AAINV2G(JG),AAINV3G(JG)                             
C                                                                       
C=================================================================      
C                                                                       
C      OPEN(10)                                                         
      READ(10) NFIN                                                     
      READ(10) NAVE,AVETIM                                              
      READ(10) AUTUB,AUTUT                                              
      READ(10) AUTUB2,AUTUT2                                            
      READ(10) AAPG,AAWG,AAUG                                           
      READ(10) APRMSG,AWRMSG,AURMSG                                     
      READ(10) ATSSTG,ARESTG,AAVG,AVRMSG                                
      READ(10) AU11G,AU22G,AU33G                                        
      READ(10) AU12G,AU13G                                              
      READ(10) AU23G,AAKG                                               
      READ(10) APKG,ATAKG,APIKG                                         
      READ(10) ADKG,AEKG,APDKG                                          
      READ(10) AP11G,ATA11G,API11G                                      
      READ(10) AD11G,AE11G,APS11G                                       
      READ(10) ATA22G,API22G                                            
      READ(10) AD22G,AE22G,APD22G,APS22G                                
      READ(10) ATA33G,API33G                                            
      READ(10) AD33G,AE33G,APS33G                                       
      READ(10) AP12G,ATA12G,API12G                                      
      READ(10) AD12G,AE12G,APD12G,APS12G                                
      READ(10) ASTRMG                                                   
      READ(10) ASPANG                                                   
      READ(10) ATPC1G                                                   
      READ(10) ATPC3G                                                   
      READ(10) OMGG                                                     
      READ(10) SKEWUVG,FLATUVG                                          
      READ(10) AUTRIG,AUQUADG,AQAUVG                                    
      READ(10) ANQAUVG                                                  
      READ(10) AAINV2G, AAINV3G                                         
C      CLOSE(10)                                                        
C                                                                       
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     OUTPUT TIME-AVERAGE DATA (VELOCITY FIELD)                         
C     PROGRAMED BY H.ABE                                                
      SUBROUTINE AVEOUT                                                 
C     REVISED ON 2001.03.03 BY H.ABE                                    
C     1998.01.03                                                        
C====================================================================   
C                                                                       
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
      INCLUDE './INCFILE/INCDATA'                                         
C      INCLUDE './INCFILE/INCMESH'                                        
C                                                                       
      COMMON /ANUM/ NAVE                                                
      COMMON /ATIM/ AVETIM                                              
C                                                                       
      COMMON /JAVE/ JA,JB,JC,JD,JE,JF,JH,JI,JJ,JK,JL                    
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /UTAU/ AUTUB,AUTUT                                         
      COMMON /UTAU2/ AUTUB2,AUTUT2                                      
      COMMON /AVE1/ AAUG(JG),AAVG(JG),AAWG(JG),AAPG(JG)                 
      COMMON /AVE2/ AURMSG(JG),AVRMSG(JG),AWRMSG(JG),APRMSG(JG)         
      COMMON /AVE3/ ATSSTG(JG),ARESTG(JG)                               
      COMMON /AVE4/ AU11G(JG),AU12G(JG),AU13G(JG)                       
     &             ,AU22G(JG),AU23G(JG),AU33G(JG)                       
      COMMON /SFUV/ SKEWUVG(JG),FLATUVG(JG)                             
      COMMON /BUDG/ AAKG(JG)                                            
     &         ,APKG(JG),AP11G(JG),AP12G(JG)                            
     &         ,APIKG(JG),API11G(JG),API22G(JG),API33G(JG),API12G(JG)   
     &         ,ADKG(JG),AD11G(JG),AD22G(JG),AD33G(JG),AD12G(JG)        
     &         ,AEKG(JG),AE11G(JG),AE22G(JG),AE33G(JG),AE12G(JG)        
     &         ,ATAKG(JG),ATA11G(JG),ATA22G(JG),ATA33G(JG),ATA12G(JG)   
     &         ,APDKG(JG),APD22G(JG),APD12G(JG)                         
     &         ,APS11G(JG),APS22G(JG),APS33G(JG),APS12G(JG)             
      COMMON /POWR/ ASTRMG(IG,JG,5),ASPANG(KG,JG,5)                     
      COMMON /CORR/ ATPC1G(IG,JG,4),ATPC3G(KG,JG,4)                     
      COMMON /VORT/ OMGG(JG,3)                                          
C                                                                       
      COMMON /UHOA/ AUTRIG(JG,8),AUQUADG(JG,10),AQAUVG(JG,4)            
      COMMON /NQAUV/ ANQAUVG(JG,4)                                      
C                                                                       
      COMMON /AAIN/ AAINV2G(JG),AAINV3G(JG)                             
C                                                                       
C=================================================================      
C                                                                       
C                                                                       
C      OPEN(11)                                                         
      WRITE(11) NFIN                                                    
      WRITE(11) NAVE,AVETIM                                             
      WRITE(11) AUTUB,AUTUT                                             
      WRITE(11) AUTUB2,AUTUT2                                           
      WRITE(11) AAPG,AAWG,AAUG                                          
      WRITE(11) APRMSG,AWRMSG,AURMSG                                    
      WRITE(11) ATSSTG,ARESTG,AAVG,AVRMSG                               
      WRITE(11) AU11G,AU22G,AU33G                                       
      WRITE(11) AU12G,AU13G                                             
      WRITE(11) AU23G,AAKG                                              
      WRITE(11) APKG,ATAKG,APIKG                                        
      WRITE(11) ADKG,AEKG,APDKG                                         
      WRITE(11) AP11G,ATA11G,API11G                                     
      WRITE(11) AD11G,AE11G,APS11G                                      
      WRITE(11) ATA22G,API22G                                           
      WRITE(11) AD22G,AE22G,APD22G,APS22G                               
      WRITE(11) ATA33G,API33G                                           
      WRITE(11) AD33G,AE33G,APS33G                                      
      WRITE(11) AP12G,ATA12G,API12G                                     
      WRITE(11) AD12G,AE12G,APD12G,APS12G                               
      WRITE(11) ASTRMG                                                  
      WRITE(11) ASPANG                                                  
      WRITE(11) ATPC1G                                                  
      WRITE(11) ATPC3G                                                  
      WRITE(11) OMGG                                                    
      WRITE(11) SKEWUVG,FLATUVG                                         
      WRITE(11) AUTRIG,AUQUADG,AQAUVG                                   
      WRITE(11) ANQAUVG                                                 
      WRITE(11) AAINV2G, AAINV3G                                        
C      CLOSE(11)                                                        
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     OUTPUT TIME-AVERAGE DATA (VELOCITY FIELD)                         
C     PROGRAMED BY H.ABE                                                
      SUBROUTINE OUTPUT(CODE)                                           
C     1998.01.01                                                        
C     REVISED ON 2004.12.20 T.TSUKAHARA                                 
C     REVISED ON 2001.03.04 H.ABE                                       
C====================================================================   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
      INCLUDE './INCFILE/INCDATA'                                         
      INCLUDE './INCFILE/INCMESH'                                         
C                                                                       
      CHARACTER*32 CODE                                                 
      CHARACTER*71 CPRM(102),CMTD(3)                                    
      CHARACTER*52 CMSP(2),CMCR(2)                                      
C                                                                       
      DIMENSION Y(JG),YY(JG)                                            
      DIMENSION RESK(JG),RES11(JG),RES22(JG),RES33(JG),RES12(JG)        
      DIMENSION SKEW(JG,4),FLAT(JG,4)                                   
C                                                                       
      COMMON /ANUM/ NAVE                                                
      COMMON /ATIM/ AVETIM                                              
      COMMON /JAVE/ JA,JB,JC,JD,JE,JF,JH,JI,JJ,JK,JL                    
C                                                                       
C----GLOBAL VARIABLES------------------------------------------         
C                                                                       
      COMMON /MEAN/ AUG(-1:JG+2),AVG(-1:JG+2)                           
     &             ,AWG(-1:JG+2),APG(-1:JG+2)                           
      COMMON /UTAU/ AUTUB,AUTUT                                         
      COMMON /UTAU2/ AUTUB2,AUTUT2                                      
      COMMON /AVE1/ AAUG(JG),AAVG(JG),AAWG(JG),AAPG(JG)                 
      COMMON /AVE2/ AURMSG(JG),AVRMSG(JG),AWRMSG(JG),APRMSG(JG)         
      COMMON /AVE3/ ATSSTG(JG),ARESTG(JG)                               
      COMMON /AVE4/ AU11G(JG),AU12G(JG),AU13G(JG)                       
     &             ,AU22G(JG),AU23G(JG),AU33G(JG)                       
      COMMON /SFUV/ SKEWUVG(JG),FLATUVG(JG)                             
      COMMON /BUDG/ AAKG(JG)                                            
     &         ,APKG(JG),AP11G(JG),AP12G(JG)                            
     &         ,APIKG(JG),API11G(JG),API22G(JG),API33G(JG),API12G(JG)   
     &         ,ADKG(JG),AD11G(JG),AD22G(JG),AD33G(JG),AD12G(JG)        
     &         ,AEKG(JG),AE11G(JG),AE22G(JG),AE33G(JG),AE12G(JG)        
     &         ,ATAKG(JG),ATA11G(JG),ATA22G(JG),ATA33G(JG),ATA12G(JG)   
     &         ,APDKG(JG),APD22G(JG),APD12G(JG)                         
     &         ,APS11G(JG),APS22G(JG),APS33G(JG),APS12G(JG)             
      COMMON /POWR/ ASTRMG(IG,JG,5),ASPANG(KG,JG,5)                     
      COMMON /CORR/ ATPC1G(IG,JG,4),ATPC3G(KG,JG,4)                     
      COMMON /VORT/ OMGG(JG,3)                                          
C                                                                       
      COMMON /UHOA/ AUTRIG(JG,8),AUQUADG(JG,10),AQAUVG(JG,4)            
      COMMON /NQAUV/ ANQAUVG(JG,4)                                      
C                                                                       
      COMMON /AAIN/ AAINV2G(JG),AAINV3G(JG)                             
C                                                                       
C=================================================================      
C                                                                       
C                                                                       
      PAVE=1.0/DBLE(NAVE)                                               
      PRE=1.0/RE                                                        
      PRERE=1.0/(RE**2)                                                 
      DOOP=1.0D0/DBLE(IG*KG)                                            
C                                                                       
      AUTUT=AUTUT*PAVE                                                  
      AUTUB=AUTUB*PAVE                                                  
      AUTUT2=AUTUT2*PAVE                                                
      AUTUB2=AUTUB2*PAVE                                                
      PVISC=(AUTUT+AUTUB)*0.5D0*RE                                      
C                                                                       
      UBULK=0.0D0                                                       
      ABULK=0.0D0                                                       
C                                                                       
      DO 100 J=1,JG                                                     
          AAUG(J)=AAUG(J)*PAVE                                          
          AAVG(J)=AAVG(J)*PAVE                                          
          AAWG(J)=AAWG(J)*PAVE                                          
          AAPG(J)=AAPG(J)*PAVE                                          
        AURMSG(J)=DSQRT(AURMSG(J)*PAVE)                                 
        AVRMSG(J)=DSQRT(AVRMSG(J)*PAVE)                                 
        AWRMSG(J)=DSQRT(AWRMSG(J)*PAVE)                                 
        APRMSG(J)=DSQRT(APRMSG(J)*PAVE)                                 
        ATSSTG(J)=ATSSTG(J)*PAVE                                        
        ARESTG(J)=ARESTG(J)*PAVE                                        
        AU11G(J)=AU11G(J)*PAVE                                          
        AU12G(J)=AU12G(J)*PAVE                                          
        AU13G(J)=AU13G(J)*PAVE                                          
        AU22G(J)=AU22G(J)*PAVE                                          
        AU23G(J)=AU23G(J)*PAVE                                          
        AU33G(J)=AU33G(J)*PAVE                                          
C                                                                       
        UBULK=UBULK+DR(J)* AUG(J)                                       
        ABULK=ABULK+DR(J)*AAUG(J)                                       
C                                                                       
        AUTRIG(J,1)=AUTRIG(J,1)*PAVE                                    
        AUTRIG(J,2)=AUTRIG(J,2)*PAVE                                    
        AUTRIG(J,3)=AUTRIG(J,3)*PAVE                                    
        AUTRIG(J,4)=AUTRIG(J,4)*PAVE                                    
        AUTRIG(J,5)=AUTRIG(J,5)*PAVE                                    
        AUTRIG(J,6)=AUTRIG(J,6)*PAVE                                    
        AUTRIG(J,7)=AUTRIG(J,7)*PAVE                                    
        AUTRIG(J,8)=AUTRIG(J,8)*PAVE                                    
C                                                                       
        AUQUADG(J,1)=AUQUADG(J,1)*PAVE                                  
        AUQUADG(J,2)=AUQUADG(J,2)*PAVE                                  
        AUQUADG(J,3)=AUQUADG(J,3)*PAVE                                  
        AUQUADG(J,4)=AUQUADG(J,4)*PAVE                                  
        AUQUADG(J,5)=AUQUADG(J,5)*PAVE                                  
        AUQUADG(J,6)=AUQUADG(J,6)*PAVE                                  
        AUQUADG(J,7)=AUQUADG(J,7)*PAVE                                  
        AUQUADG(J,8)=AUQUADG(J,8)*PAVE                                  
        AUQUADG(J,9)=AUQUADG(J,9)*PAVE                                  
        AUQUADG(J,10)=AUQUADG(J,10)*PAVE                                
C                                                                       
        AQAUVG(J,1)=AQAUVG(J,1)*DOOP*PAVE                               
        AQAUVG(J,2)=AQAUVG(J,2)*DOOP*PAVE                               
        AQAUVG(J,3)=AQAUVG(J,3)*DOOP*PAVE                               
        AQAUVG(J,4)=AQAUVG(J,4)*DOOP*PAVE                               
C                                                                       
        ANQAUVG(J,1)=ANQAUVG(J,1)*PAVE                                  
        ANQAUVG(J,2)=ANQAUVG(J,2)*PAVE                                  
        ANQAUVG(J,3)=ANQAUVG(J,3)*PAVE                                  
        ANQAUVG(J,4)=ANQAUVG(J,4)*PAVE                                  
C                                                                       
        AAKG(J)=AAKG(J)*PAVE                                            
        APKG(J)=APKG(J)*PAVE*PRE                                        
        AP11G(J)=AP11G(J)*PAVE*PRE                                      
        AP12G(J)=AP12G(J)*PAVE*PRE                                      
        ATAKG(J)=ATAKG(J)*PAVE*PRE                                      
        ATA11G(J)=ATA11G(J)*PAVE*PRE                                    
        ATA22G(J)=ATA22G(J)*PAVE*PRE                                    
        ATA33G(J)=ATA33G(J)*PAVE*PRE                                    
        ATA12G(J)=ATA12G(J)*PAVE*PRE                                    
        APIKG(J)=APIKG(J)*PAVE*PRE                                      
        API11G(J)=API11G(J)*PAVE*PRE                                    
        API22G(J)=API22G(J)*PAVE*PRE                                    
        API33G(J)=API33G(J)*PAVE*PRE                                    
        API12G(J)=API12G(J)*PAVE*PRE                                    
        ADKG(J)=ADKG(J)*PAVE*PRERE                                      
        AD11G(J)=AD11G(J)*PAVE*PRERE                                    
        AD22G(J)=AD22G(J)*PAVE*PRERE                                    
        AD33G(J)=AD33G(J)*PAVE*PRERE                                    
        AD12G(J)=AD12G(J)*PAVE*PRERE                                    
        AEKG(J)=AEKG(J)*PAVE*PRERE                                      
        AE11G(J)=AE11G(J)*PAVE*PRERE                                    
        AE22G(J)=AE22G(J)*PAVE*PRERE                                    
        AE33G(J)=AE33G(J)*PAVE*PRERE                                    
        AE12G(J)=AE12G(J)*PAVE*PRERE                                    
        RESK(J)=APKG(J)+ATAKG(J)+APIKG(J)+ADKG(J)-AEKG(J)               
        RES11(J)=AP11G(J)+ATA11G(J)+API11G(J)+AD11G(J)-AE11G(J)         
        RES22(J)=ATA22G(J)+API22G(J)+AD22G(J)-AE22G(J)                  
        RES33(J)=ATA33G(J)+API33G(J)+AD33G(J)-AE33G(J)                  
        RES12(J)=AP12G(J)+ATA12G(J)+API12G(J)+AD12G(J)-AE12G(J)         
C                                                                       
        APDKG(J)=APDKG(J)*PAVE*PRE                                      
        APD22G(J)=APD22G(J)*PAVE*PRE                                    
        APD12G(J)=APD12G(J)*PAVE*PRE                                    
        APS11G(J)=APS11G(J)*PAVE*PRE                                    
        APS22G(J)=APS22G(J)*PAVE*PRE                                    
        APS33G(J)=APS33G(J)*PAVE*PRE                                    
        APS12G(J)=APS12G(J)*PAVE*PRE                                    
C                                                                       
        OMGG(J,1)=DSQRT(OMGG(J,1)*PAVE)*PRE                             
        OMGG(J,2)=DSQRT(OMGG(J,2)*PAVE)*PRE                             
        OMGG(J,3)=DSQRT(OMGG(J,3)*PAVE)*PRE                             
C                                                                       
        ANQAUVG(J,1)=ANQAUVG(J,1)*PAVE                                  
        ANQAUVG(J,2)=ANQAUVG(J,2)*PAVE                                  
        ANQAUVG(J,3)=ANQAUVG(J,3)*PAVE                                  
        ANQAUVG(J,4)=ANQAUVG(J,4)*PAVE                                  
C                                                                       
        AAINV2G(J)=AAINV2G(J)*PAVE                                      
        AAINV3G(J)=AAINV3G(J)*PAVE                                      
C                                                                       
C========> FOR NEW_MESH                                                 
        Y(J)=R(J)*RE                                                    
        YY(J)=RR(J)*RE                                                  
  100 CONTINUE                                                          
C                                                                       
      DO 101 I=1,IG                                                     
        ASTRMG(I,JA,1)=ASTRMG(I,JA,1)*PAVE                              
        ASTRMG(I,JA,2)=ASTRMG(I,JA,2)*PAVE                              
        ASTRMG(I,JA,3)=ASTRMG(I,JA,3)*PAVE                              
        ASTRMG(I,JA,4)=ASTRMG(I,JA,4)*PAVE                              
        ASTRMG(I,JA,5)=ASTRMG(I,JA,5)*PAVE                              
C                                                                       
        ASTRMG(I,JB,1)=ASTRMG(I,JB,1)*PAVE                              
        ASTRMG(I,JB,2)=ASTRMG(I,JB,2)*PAVE                              
        ASTRMG(I,JB,3)=ASTRMG(I,JB,3)*PAVE                              
        ASTRMG(I,JB,4)=ASTRMG(I,JB,4)*PAVE                              
        ASTRMG(I,JB,5)=ASTRMG(I,JB,5)*PAVE                              
C                                                                       
        ASTRMG(I,JC,1)=ASTRMG(I,JC,1)*PAVE                              
        ASTRMG(I,JC,2)=ASTRMG(I,JC,2)*PAVE                              
        ASTRMG(I,JC,3)=ASTRMG(I,JC,3)*PAVE                              
        ASTRMG(I,JC,4)=ASTRMG(I,JC,4)*PAVE                              
        ASTRMG(I,JC,5)=ASTRMG(I,JC,5)*PAVE                              
C                                                                       
        ASTRMG(I,JD,1)=ASTRMG(I,JD,1)*PAVE                              
        ASTRMG(I,JD,2)=ASTRMG(I,JD,2)*PAVE                              
        ASTRMG(I,JD,3)=ASTRMG(I,JD,3)*PAVE                              
        ASTRMG(I,JD,4)=ASTRMG(I,JD,4)*PAVE                              
        ASTRMG(I,JD,5)=ASTRMG(I,JD,5)*PAVE                              
C                                                                       
        ASTRMG(I,JE,1)=ASTRMG(I,JE,1)*PAVE                              
        ASTRMG(I,JE,2)=ASTRMG(I,JE,2)*PAVE                              
        ASTRMG(I,JE,3)=ASTRMG(I,JE,3)*PAVE                              
        ASTRMG(I,JE,4)=ASTRMG(I,JE,4)*PAVE                              
        ASTRMG(I,JE,5)=ASTRMG(I,JE,5)*PAVE                              
C                                                                       
        ASTRMG(I,JF,1)=ASTRMG(I,JF,1)*PAVE                              
        ASTRMG(I,JF,2)=ASTRMG(I,JF,2)*PAVE                              
        ASTRMG(I,JF,3)=ASTRMG(I,JF,3)*PAVE                              
        ASTRMG(I,JF,4)=ASTRMG(I,JF,4)*PAVE                              
        ASTRMG(I,JF,5)=ASTRMG(I,JF,5)*PAVE                              
C                                                                       
        ASTRMG(I,JH,1)=ASTRMG(I,JH,1)*PAVE                              
        ASTRMG(I,JH,2)=ASTRMG(I,JH,2)*PAVE                              
        ASTRMG(I,JH,3)=ASTRMG(I,JH,3)*PAVE                              
        ASTRMG(I,JH,4)=ASTRMG(I,JH,4)*PAVE                              
        ASTRMG(I,JH,5)=ASTRMG(I,JH,5)*PAVE                              
C                                                                       
        ASTRMG(I,JI,1)=ASTRMG(I,JI,1)*PAVE                              
        ASTRMG(I,JI,2)=ASTRMG(I,JI,2)*PAVE                              
        ASTRMG(I,JI,3)=ASTRMG(I,JI,3)*PAVE                              
        ASTRMG(I,JI,4)=ASTRMG(I,JI,4)*PAVE                              
        ASTRMG(I,JI,5)=ASTRMG(I,JI,5)*PAVE                              
C                                                                       
        ASTRMG(I,JJ,1)=ASTRMG(I,JJ,1)*PAVE                              
        ASTRMG(I,JJ,2)=ASTRMG(I,JJ,2)*PAVE                              
        ASTRMG(I,JJ,3)=ASTRMG(I,JJ,3)*PAVE                              
        ASTRMG(I,JJ,4)=ASTRMG(I,JJ,4)*PAVE                              
        ASTRMG(I,JJ,5)=ASTRMG(I,JJ,5)*PAVE                              
C                                                                       
        ASTRMG(I,JK,1)=ASTRMG(I,JK,1)*PAVE                              
        ASTRMG(I,JK,2)=ASTRMG(I,JK,2)*PAVE                              
        ASTRMG(I,JK,3)=ASTRMG(I,JK,3)*PAVE                              
        ASTRMG(I,JK,4)=ASTRMG(I,JK,4)*PAVE                              
        ASTRMG(I,JK,5)=ASTRMG(I,JK,5)*PAVE                              
C                                                                       
        ASTRMG(I,JL,1)=ASTRMG(I,JL,1)*PAVE                              
        ASTRMG(I,JL,2)=ASTRMG(I,JL,2)*PAVE                              
        ASTRMG(I,JL,3)=ASTRMG(I,JL,3)*PAVE                              
        ASTRMG(I,JL,4)=ASTRMG(I,JL,4)*PAVE                              
        ASTRMG(I,JL,5)=ASTRMG(I,JL,5)*PAVE                              
C                                                                       
        ATPC1G(I,JA,1)=ATPC1G(I,JA,1)*PAVE                              
        ATPC1G(I,JB,1)=ATPC1G(I,JB,1)*PAVE                              
        ATPC1G(I,JC,1)=ATPC1G(I,JC,1)*PAVE                              
        ATPC1G(I,JD,1)=ATPC1G(I,JD,1)*PAVE                              
        ATPC1G(I,JE,1)=ATPC1G(I,JE,1)*PAVE                              
        ATPC1G(I,JF,1)=ATPC1G(I,JF,1)*PAVE                              
        ATPC1G(I,JH,1)=ATPC1G(I,JH,1)*PAVE                              
        ATPC1G(I,JI,1)=ATPC1G(I,JI,1)*PAVE                              
        ATPC1G(I,JJ,1)=ATPC1G(I,JJ,1)*PAVE                              
        ATPC1G(I,JK,1)=ATPC1G(I,JK,1)*PAVE                              
        ATPC1G(I,JL,1)=ATPC1G(I,JL,1)*PAVE                              
C                                                                       
        ATPC1G(I,JA,2)=ATPC1G(I,JA,2)*PAVE                              
        ATPC1G(I,JB,2)=ATPC1G(I,JB,2)*PAVE                              
        ATPC1G(I,JC,2)=ATPC1G(I,JC,2)*PAVE                              
        ATPC1G(I,JD,2)=ATPC1G(I,JD,2)*PAVE                              
        ATPC1G(I,JE,2)=ATPC1G(I,JE,2)*PAVE                              
        ATPC1G(I,JF,2)=ATPC1G(I,JF,2)*PAVE                              
        ATPC1G(I,JH,2)=ATPC1G(I,JH,2)*PAVE                              
        ATPC1G(I,JI,2)=ATPC1G(I,JI,2)*PAVE                              
        ATPC1G(I,JJ,2)=ATPC1G(I,JJ,2)*PAVE                              
        ATPC1G(I,JK,2)=ATPC1G(I,JK,2)*PAVE                              
        ATPC1G(I,JL,2)=ATPC1G(I,JL,2)*PAVE                              
C                                                                       
        ATPC1G(I,JA,3)=ATPC1G(I,JA,3)*PAVE                              
        ATPC1G(I,JB,3)=ATPC1G(I,JB,3)*PAVE                              
        ATPC1G(I,JC,3)=ATPC1G(I,JC,3)*PAVE                              
        ATPC1G(I,JD,3)=ATPC1G(I,JD,3)*PAVE                              
        ATPC1G(I,JE,3)=ATPC1G(I,JE,3)*PAVE                              
        ATPC1G(I,JF,3)=ATPC1G(I,JF,3)*PAVE                              
        ATPC1G(I,JH,3)=ATPC1G(I,JH,3)*PAVE                              
        ATPC1G(I,JI,3)=ATPC1G(I,JI,3)*PAVE                              
        ATPC1G(I,JJ,3)=ATPC1G(I,JJ,3)*PAVE                              
        ATPC1G(I,JK,3)=ATPC1G(I,JK,3)*PAVE                              
        ATPC1G(I,JL,3)=ATPC1G(I,JL,3)*PAVE                              
C                                                                       
        ATPC1G(I,JA,4)=ATPC1G(I,JA,4)*PAVE                              
        ATPC1G(I,JB,4)=ATPC1G(I,JB,4)*PAVE                              
        ATPC1G(I,JC,4)=ATPC1G(I,JC,4)*PAVE                              
        ATPC1G(I,JD,4)=ATPC1G(I,JD,4)*PAVE                              
        ATPC1G(I,JE,4)=ATPC1G(I,JE,4)*PAVE                              
        ATPC1G(I,JF,4)=ATPC1G(I,JF,4)*PAVE                              
        ATPC1G(I,JH,4)=ATPC1G(I,JH,4)*PAVE                              
        ATPC1G(I,JI,4)=ATPC1G(I,JI,4)*PAVE                              
        ATPC1G(I,JJ,4)=ATPC1G(I,JJ,4)*PAVE                              
        ATPC1G(I,JK,4)=ATPC1G(I,JK,4)*PAVE                              
        ATPC1G(I,JL,4)=ATPC1G(I,JL,4)*PAVE                              
  101 CONTINUE                                                          
C                                                                       
      DO 102 K=1,KG                                                     
        ASPANG(K,JA,1)=ASPANG(K,JA,1)*PAVE                              
        ASPANG(K,JA,2)=ASPANG(K,JA,2)*PAVE                              
        ASPANG(K,JA,3)=ASPANG(K,JA,3)*PAVE                              
        ASPANG(K,JA,4)=ASPANG(K,JA,4)*PAVE                              
        ASPANG(K,JA,5)=ASPANG(K,JA,5)*PAVE                              
C                                                                       
        ASPANG(K,JB,1)=ASPANG(K,JB,1)*PAVE                              
        ASPANG(K,JB,2)=ASPANG(K,JB,2)*PAVE                              
        ASPANG(K,JB,3)=ASPANG(K,JB,3)*PAVE                              
        ASPANG(K,JB,4)=ASPANG(K,JB,4)*PAVE                              
        ASPANG(K,JB,5)=ASPANG(K,JB,5)*PAVE                              
C                                                                       
        ASPANG(K,JC,1)=ASPANG(K,JC,1)*PAVE                              
        ASPANG(K,JC,2)=ASPANG(K,JC,2)*PAVE                              
        ASPANG(K,JC,3)=ASPANG(K,JC,3)*PAVE                              
        ASPANG(K,JC,4)=ASPANG(K,JC,4)*PAVE                              
        ASPANG(K,JC,5)=ASPANG(K,JC,5)*PAVE                              
C                                                                       
        ASPANG(K,JD,1)=ASPANG(K,JD,1)*PAVE                              
        ASPANG(K,JD,2)=ASPANG(K,JD,2)*PAVE                              
        ASPANG(K,JD,3)=ASPANG(K,JD,3)*PAVE                              
        ASPANG(K,JD,4)=ASPANG(K,JD,4)*PAVE                              
        ASPANG(K,JD,5)=ASPANG(K,JD,5)*PAVE                              
C                                                                       
        ASPANG(K,JE,1)=ASPANG(K,JE,1)*PAVE                              
        ASPANG(K,JE,2)=ASPANG(K,JE,2)*PAVE                              
        ASPANG(K,JE,3)=ASPANG(K,JE,3)*PAVE                              
        ASPANG(K,JE,4)=ASPANG(K,JE,4)*PAVE                              
        ASPANG(K,JE,5)=ASPANG(K,JE,5)*PAVE                              
C                                                                       
        ASPANG(K,JF,1)=ASPANG(K,JF,1)*PAVE                              
        ASPANG(K,JF,2)=ASPANG(K,JF,2)*PAVE                              
        ASPANG(K,JF,3)=ASPANG(K,JF,3)*PAVE                              
        ASPANG(K,JF,4)=ASPANG(K,JF,4)*PAVE                              
        ASPANG(K,JF,5)=ASPANG(K,JF,5)*PAVE                              
C                                                                       
        ASPANG(K,JH,1)=ASPANG(K,JH,1)*PAVE                              
        ASPANG(K,JH,2)=ASPANG(K,JH,2)*PAVE                              
        ASPANG(K,JH,3)=ASPANG(K,JH,3)*PAVE                              
        ASPANG(K,JH,4)=ASPANG(K,JH,4)*PAVE                              
        ASPANG(K,JH,5)=ASPANG(K,JH,5)*PAVE                              
C                                                                       
        ASPANG(K,JI,1)=ASPANG(K,JI,1)*PAVE                              
        ASPANG(K,JI,2)=ASPANG(K,JI,2)*PAVE                              
        ASPANG(K,JI,3)=ASPANG(K,JI,3)*PAVE                              
        ASPANG(K,JI,4)=ASPANG(K,JI,4)*PAVE                              
        ASPANG(K,JI,5)=ASPANG(K,JI,5)*PAVE                              
C                                                                       
        ASPANG(K,JJ,1)=ASPANG(K,JJ,1)*PAVE                              
        ASPANG(K,JJ,2)=ASPANG(K,JJ,2)*PAVE                              
        ASPANG(K,JJ,3)=ASPANG(K,JJ,3)*PAVE                              
        ASPANG(K,JJ,4)=ASPANG(K,JJ,4)*PAVE                              
        ASPANG(K,JJ,5)=ASPANG(K,JJ,5)*PAVE                              
C                                                                       
        ASPANG(K,JK,1)=ASPANG(K,JK,1)*PAVE                              
        ASPANG(K,JK,2)=ASPANG(K,JK,2)*PAVE                              
        ASPANG(K,JK,3)=ASPANG(K,JK,3)*PAVE                              
        ASPANG(K,JK,4)=ASPANG(K,JK,4)*PAVE                              
        ASPANG(K,JK,5)=ASPANG(K,JK,5)*PAVE                              
C                                                                       
        ASPANG(K,JL,1)=ASPANG(K,JL,1)*PAVE                              
        ASPANG(K,JL,2)=ASPANG(K,JL,2)*PAVE                              
        ASPANG(K,JL,3)=ASPANG(K,JL,3)*PAVE                              
        ASPANG(K,JL,4)=ASPANG(K,JL,4)*PAVE                              
        ASPANG(K,JL,5)=ASPANG(K,JL,5)*PAVE                              
C                                                                       
        ATPC3G(K,JA,1)=ATPC3G(K,JA,1)*PAVE                              
        ATPC3G(K,JB,1)=ATPC3G(K,JB,1)*PAVE                              
        ATPC3G(K,JC,1)=ATPC3G(K,JC,1)*PAVE                              
        ATPC3G(K,JD,1)=ATPC3G(K,JD,1)*PAVE                              
        ATPC3G(K,JE,1)=ATPC3G(K,JE,1)*PAVE                              
        ATPC3G(K,JF,1)=ATPC3G(K,JF,1)*PAVE                              
        ATPC3G(K,JH,1)=ATPC3G(K,JH,1)*PAVE                              
        ATPC3G(K,JI,1)=ATPC3G(K,JI,1)*PAVE                              
        ATPC3G(K,JJ,1)=ATPC3G(K,JJ,1)*PAVE                              
        ATPC3G(K,JK,1)=ATPC3G(K,JK,1)*PAVE                              
        ATPC3G(K,JL,1)=ATPC3G(K,JL,1)*PAVE                              
C                                                                       
        ATPC3G(K,JA,2)=ATPC3G(K,JA,2)*PAVE                              
        ATPC3G(K,JB,2)=ATPC3G(K,JB,2)*PAVE                              
        ATPC3G(K,JC,2)=ATPC3G(K,JC,2)*PAVE                              
        ATPC3G(K,JD,2)=ATPC3G(K,JD,2)*PAVE                              
        ATPC3G(K,JE,2)=ATPC3G(K,JE,2)*PAVE                              
        ATPC3G(K,JF,2)=ATPC3G(K,JF,2)*PAVE                              
        ATPC3G(K,JH,2)=ATPC3G(K,JH,2)*PAVE                              
        ATPC3G(K,JI,2)=ATPC3G(K,JI,2)*PAVE                              
        ATPC3G(K,JJ,2)=ATPC3G(K,JJ,2)*PAVE                              
        ATPC3G(K,JK,2)=ATPC3G(K,JK,2)*PAVE                              
        ATPC3G(K,JL,2)=ATPC3G(K,JL,2)*PAVE                              
C                                                                       
        ATPC3G(K,JA,3)=ATPC3G(K,JA,3)*PAVE                              
        ATPC3G(K,JB,3)=ATPC3G(K,JB,3)*PAVE                              
        ATPC3G(K,JC,3)=ATPC3G(K,JC,3)*PAVE                              
        ATPC3G(K,JD,3)=ATPC3G(K,JD,3)*PAVE                              
        ATPC3G(K,JE,3)=ATPC3G(K,JE,3)*PAVE                              
        ATPC3G(K,JF,3)=ATPC3G(K,JF,3)*PAVE                              
        ATPC3G(K,JH,3)=ATPC3G(K,JH,3)*PAVE                              
        ATPC3G(K,JI,3)=ATPC3G(K,JI,3)*PAVE                              
        ATPC3G(K,JJ,3)=ATPC3G(K,JJ,3)*PAVE                              
        ATPC3G(K,JK,3)=ATPC3G(K,JK,3)*PAVE                              
        ATPC3G(K,JL,3)=ATPC3G(K,JL,3)*PAVE                              
C                                                                       
        ATPC3G(K,JA,4)=ATPC3G(K,JA,4)*PAVE                              
        ATPC3G(K,JB,4)=ATPC3G(K,JB,4)*PAVE                              
        ATPC3G(K,JC,4)=ATPC3G(K,JC,4)*PAVE                              
        ATPC3G(K,JD,4)=ATPC3G(K,JD,4)*PAVE                              
        ATPC3G(K,JE,4)=ATPC3G(K,JE,4)*PAVE                              
        ATPC3G(K,JF,4)=ATPC3G(K,JF,4)*PAVE                              
        ATPC3G(K,JH,4)=ATPC3G(K,JH,4)*PAVE                              
        ATPC3G(K,JI,4)=ATPC3G(K,JI,4)*PAVE                              
        ATPC3G(K,JJ,4)=ATPC3G(K,JJ,4)*PAVE                              
        ATPC3G(K,JK,4)=ATPC3G(K,JK,4)*PAVE                              
        ATPC3G(K,JL,4)=ATPC3G(K,JL,4)*PAVE                              
  102 CONTINUE                                                          
C                                                                       
C---------------SKEWNESS & FLATNESS FACTOR                              
C                                                                       
      DO 103 J=1,JG-1                                                   
        SKEW(J,1)=AUTRIG(J,1)/AURMSG(J)**3                              
        SKEW(J,2)=AUTRIG(J,3)/AVRMSG(J)**3                              
        SKEW(J,3)=AUTRIG(J,7)/AWRMSG(J)**3                              
        SKEW(J,4)=AUTRIG(J,8)/APRMSG(J)**3                              
        FLAT(J,1)=AUQUADG(J,1)/AURMSG(J)**4                             
        FLAT(J,2)=AUQUADG(J,3)/AVRMSG(J)**4                             
        FLAT(J,3)=AUQUADG(J,9)/AWRMSG(J)**4                             
        FLAT(J,4)=AUQUADG(J,10)/APRMSG(J)**4                            
  103 CONTINUE                                                          
C                                                                       
      DO 104 J=JG,JG                                                    
        SKEW(J,1)=AUTRIG(J,1)/AURMSG(J)**3                              
        SKEW(J,3)=AUTRIG(J,7)/AWRMSG(J)**3                              
        SKEW(J,4)=AUTRIG(J,8)/APRMSG(J)**3                              
        FLAT(J,1)=AUQUADG(J,1)/AURMSG(J)**4                             
        FLAT(J,3)=AUQUADG(J,9)/AWRMSG(J)**4                             
        FLAT(J,4)=AUQUADG(J,10)/APRMSG(J)**4                            
  104 CONTINUE                                                          
C                                                                       
      OO=0.0D0                                                          
C                                                                       
      CMTD( 1)='#======================================================'
      CMTD( 2)='# DNS OF TURBULENT CHANNEL FLOW : '//CODE               
      CMTD( 3)='# TURBULENT STATISTICS IN VELOCITY FIELD'               
      CPRM( 1)='# STRM W-W SPAN  TIMESTEP  REYNOLDS  TI'                
      CPRM( 2)='# X-LENGTH  Z-LENGTH  DRAWNUM FININUM'                  
      CPRM( 3)='#  ALP       R1        R2        DX        DZ'          
      CPRM( 4)='#  UTAUB     UTAUT'                                     
      CPRM( 5)='#J R DR RR DRR'                                         
      CPRM( 6)='#YPLUS P-AVE W-AVE U-AVE'                               
      CPRM( 7)='#YPLUS P-RMS W-RMS U-RMS'                               
      CPRM( 8)='#YPLUS TOTAL-STRESS RE-STRESS'                          
      CPRM( 9)='#YPLUS UU VV WW '                                       
      CPRM(10)='#YPLUS UV UW '                                          
      CPRM(11)='#YPLUS VW TURBULENT-ENERGY'                             
      CPRM(12)='#BUDGET-K   YPLUS  PRO  T-DIF  V-P  P-DIF'              
      CPRM(13)='#BUDGET-K   YPLUS  V-DIF  DISS  RES'                    
      CPRM(14)='#BUDGET-UU  YPLUS  PRO  T-DIF  V-P'                     
      CPRM(15)='#BUDGET-UU  YPLUS  P-S  V-DIF  DISS  RES'               
      CPRM(16)='#BUDGET-VV  YPLUS  PRO  T-DIF  V-P  P-DIF'              
      CPRM(17)='#BUDGET-VV  YPLUS  P-S  V-DIF  DISS  RES'               
      CPRM(18)='#BUDGET-WW  YPLUS  PRO  T-DIF  V-P'                     
      CPRM(19)='#BUDGET-WW  YPLUS  P-S  V-DIF  DISS  RES'               
      CPRM(20)='#BUDGET-UV  YPLUS  PRO  T-DIF  V-P  P-DIF'              
      CPRM(21)='#BUDGET-UV  YPLUS  P-S  V-DIF  DISS  RES'               
C                                                                       
      CMSP(1)='# I STRM-1 STRM-2 STRM-3 STRM-P STRM-UV ENE.SPE. Y+='    
      CMSP(2)='# K SPAN-1 SPAN-2 SPAN-3 SPAN-P SPAN-UV ENE.SPE. Y+='    
      CMCR(1)='# XPLUS STRM-1 STRM-2 STRM-3 STRM-P  T.P.C. Y+='         
      CMCR(2)='# ZPLUS SPAN-1 SPAN-2 SPAN-3 SPAN-P  T.P.C. Y+='         
C                                                                       
      CPRM(22)='#YPLUS SKEW-1 SKEW-2 SKEW-3 SKEW-P  SKEWNESS FACTOR'    
      CPRM(23)='#YPLUS FLAT-1 FLAT-2 FLAT-3 FLAT-P  FLATNESS FACTOR'    
      CPRM(24)='#YPLUS OMG-1 OMG-2 OMG-3  RMS VORTICITY FLUCTUATION'    
      CPRM(25)='#YPLUS V-AVE V-RMS VV SKEW-2 FLAT-2'                    
      CPRM(26)='#YPLUS U111 U112 U222'                                  
      CPRM(27)='#YPLUS U332 U122 U133'                                  
      CPRM(28)='#YPLUS SKEWUV FLATUV'                                   
      CPRM(29)='#  UTAUB2    UTAUT2'                                    
      CPRM(30)='#CALC.NUM.     CALC.TIME      AVE.NUM.     AVE. TIME'   
      CPRM(31)='#YPLUS SKEW-UV FLAT-UV SKEWNESS & FLATNESS FACTOR'      
      CPRM(32)='#YPLUS UUUU UUUV UWWV UVVV'                             
      CPRM(33)='#YPLUS UUVV UUWW VVWW'                                  
      CPRM(34)='#YPLUS VV VVV VVVV'                                     
      CPRM(35)='#YPLUS UV-1 UV-2 UV-3 UV-4 QUADRANT ANALYSIS'           
      CPRM(36)='#YPLUS UV-1 UV-2 UV-3 UV-4 QA-PDF'                      
      CPRM(37)='#YPLUS INVARIANT-II -III VEL.GRAD.TENSOR'               
C                                                                       
C      OPEN(14)                                                         
      WRITE(14,1010) CMTD(1)                                            
      WRITE(14,1010) CMTD(2)                                            
      WRITE(14,1010) CMTD(1)                                            
      WRITE(14,1010) CPRM(30)                                           
      WRITE(14,1000) NUM,TIME,NAVE,AVETIM                               
      WRITE(14,1010) CPRM(1)                                            
      WRITE(14,1020) IG,JG,KG,DT,RE,TI                                  
      WRITE(14,1010) CPRM(2)                                            
      WRITE(14,1030) XL,ZL,NDRAW,NFIN                                   
      WRITE(14,1010) CPRM(3)                                            
      WRITE(14,1040) ALP,R1,R2,DX,DZ                                    
      WRITE(14,1010) CPRM(4)                                            
      WRITE(14,1050) AUTUB,AUTUT                                        
      WRITE(14,1010) CPRM(29)                                           
      WRITE(14,1050) AUTUB2,AUTUT2                                      
      WRITE(14,1010) CMTD(3)                                            
 1000 FORMAT(I10,E18.10,I10,E18.10)                                     
 1010 FORMAT(A71)                                                       
 1020 FORMAT(3I5,F10.7,F10.1,F10.5)                                     
 1030 FORMAT(2F10.7,I5,I10)                                             
 1040 FORMAT(5F10.7)                                                    
 1050 FORMAT(2F10.7)                                                    
C                                                                       
      WRITE(14,1010) CPRM(5)                                            
      DO 10 J=0,JG+1                                                    
        WRITE(14,1060) J,R(J),DR(J),RR(J),DRR(J)                        
   10 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(6)                                            
      DO 11 J=1,JG                                                      
        WRITE(14,1070) YY(J),AAPG(J),AAWG(J),AAUG(J)                    
   11 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(7)                                            
      DO 12 J=1,JG                                                      
        WRITE(14,1070) YY(J),APRMSG(J),AWRMSG(J),AURMSG(J)              
   12 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(8)                                            
      DO 13 J=1,JG                                                      
        WRITE(14,1090) Y(J),ATSSTG(J),ARESTG(J)                         
   13 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(9)                                            
      DO 14 J=1,JG                                                      
        WRITE(14,1070) YY(J),AU11G(J),AU22G(J),AU33G(J)                 
   14 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(10)                                           
      DO 15 J=1,JG                                                      
        WRITE(14,1090) YY(J),AU12G(J),AU13G(J)                          
   15 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(11)                                           
      DO 16 J=1,JG                                                      
        WRITE(14,1090) YY(J),AU23G(J),AAKG(J)                           
   16 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(12)                                           
      DO 17 J=1,JG                                                      
        WRITE(14,1080) YY(J),APKG(J),ATAKG(J),APIKG(J),APDKG(J)         
   17 CONTINUE                                                          
      WRITE(14,1010) CPRM(13)                                           
      DO 18 J=1,JG                                                      
        WRITE(14,1070) YY(J),ADKG(J),AEKG(J),RESK(J)                    
   18 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(14)                                           
      DO 19 J=1,JG                                                      
        WRITE(14,1070) YY(J),AP11G(J),ATA11G(J),API11G(J)               
   19 CONTINUE                                                          
      WRITE(14,1010) CPRM(15)                                           
      DO 20 J=1,JG                                                      
        WRITE(14,1080) YY(J),API11G(J),AD11G(J),AE11G(J),RES11(J)       
   20 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(16)                                           
      DO 21 J=1,JG                                                      
        WRITE(14,1080) Y(J),OO,ATA22G(J),API22G(J),APD22G(J)            
   21 CONTINUE                                                          
      WRITE(14,1010) CPRM(17)                                           
      DO 22 J=1,JG                                                      
        WRITE(14,1080) Y(J),APS22G(J),AD22G(J),AE22G(J),RES22(J)        
   22 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(18)                                           
      DO 23 J=1,JG                                                      
        WRITE(14,1070) YY(J),OO,ATA33G(J),API33G(J)                     
   23 CONTINUE                                                          
      WRITE(14,1010) CPRM(19)                                           
      DO 24 J=1,JG                                                      
        WRITE(14,1080) YY(J),APS33G(J),AD33G(J),AE33G(J),RES33(J)       
   24 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(20)                                           
      DO 25 J=1,JG                                                      
        WRITE(14,1080) YY(J),AP12G(J),ATA12G(J),API12G(J),APD12G(J)     
   25 CONTINUE                                                          
      WRITE(14,1010) CPRM(21)                                           
      DO 26 J=1,JG                                                      
        WRITE(14,1080) YY(J),APS12G(J),AD12G(J),AE12G(J),RES12(J)       
   26 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(1),RR(JA)*PVISC                       
      DO 27 I=2,IG/2+1                                                  
        WRITE(14,1100) I-1,ASTRMG(I,JA,1),ASTRMG(I,JA,2)                
     &                ,ASTRMG(I,JA,3),ASTRMG(I,JA,4),ASTRMG(I,JA,5)     
   27 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(1),RR(JB)*PVISC                       
      DO 28 I=2,IG/2+1                                                  
        WRITE(14,1100) I-1,ASTRMG(I,JB,1),ASTRMG(I,JB,2)                
     &                ,ASTRMG(I,JB,3),ASTRMG(I,JB,4),ASTRMG(I,JB,5)     
   28 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(1),RR(JC)*PVISC                       
      DO 29 I=2,IG/2+1                                                  
        WRITE(14,1100) I-1,ASTRMG(I,JC,1),ASTRMG(I,JC,2)                
     &                ,ASTRMG(I,JC,3),ASTRMG(I,JC,4),ASTRMG(I,JC,5)     
   29 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(1),RR(JD)*PVISC                       
      DO 30 I=2,IG/2+1                                                  
        WRITE(14,1100) I-1,ASTRMG(I,JD,1),ASTRMG(I,JD,2)                
     &                ,ASTRMG(I,JD,3),ASTRMG(I,JD,4),ASTRMG(I,JD,5)     
   30 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(1),RR(JE)*PVISC                       
      DO 31 I=2,IG/2+1                                                  
        WRITE(14,1100) I-1,ASTRMG(I,JE,1),ASTRMG(I,JE,2)                
     &                ,ASTRMG(I,JE,3),ASTRMG(I,JE,4),ASTRMG(I,JE,5)     
   31 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(1),RR(JF)*PVISC                       
      DO 32 I=2,IG/2+1                                                  
        WRITE(14,1100) I-1,ASTRMG(I,JF,1),ASTRMG(I,JF,2)                
     &                ,ASTRMG(I,JF,3),ASTRMG(I,JF,4),ASTRMG(I,JF,5)     
   32 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(1),RR(JH)*PVISC                       
      DO 33 I=2,IG/2+1                                                  
        WRITE(14,1100) I-1,ASTRMG(I,JH,1),ASTRMG(I,JH,2)                
     &                ,ASTRMG(I,JH,3),ASTRMG(I,JH,4),ASTRMG(I,JH,5)     
   33 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(1),RR(JI)*PVISC                       
      DO 34 I=2,IG/2+1                                                  
        WRITE(14,1100) I-1,ASTRMG(I,JI,1),ASTRMG(I,JI,2)                
     &                ,ASTRMG(I,JI,3),ASTRMG(I,JI,4),ASTRMG(I,JI,5)     
   34 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(1),RR(JJ)*PVISC                       
      DO 35 I=2,IG/2+1                                                  
        WRITE(14,1100) I-1,ASTRMG(I,JJ,1),ASTRMG(I,JJ,2)                
     &                ,ASTRMG(I,JJ,3),ASTRMG(I,JJ,4),ASTRMG(I,JJ,5)     
   35 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(1),RR(JK)*PVISC                       
      DO 36 I=2,IG/2+1                                                  
        WRITE(14,1100) I-1,ASTRMG(I,JK,1),ASTRMG(I,JK,2)                
     &                ,ASTRMG(I,JK,3),ASTRMG(I,JK,4),ASTRMG(I,JK,5)     
   36 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(1),RR(JL)*PVISC                       
      DO 37 I=2,IG/2+1                                                  
        WRITE(14,1100) I-1,ASTRMG(I,JL,1),ASTRMG(I,JL,2)                
     &                ,ASTRMG(I,JL,3),ASTRMG(I,JL,4),ASTRMG(I,JL,5)     
   37 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(2),RR(JA)*PVISC                       
      DO 43 K=2,KG/2+1                                                  
        WRITE(14,1100) K-1,ASPANG(K,JA,1),ASPANG(K,JA,2)                
     &                ,ASPANG(K,JA,3),ASPANG(K,JA,4),ASPANG(K,JA,5)     
   43 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(2),RR(JB)*PVISC                       
      DO 44 K=2,KG/2+1                                                  
        WRITE(14,1100) K-1,ASPANG(K,JB,1),ASPANG(K,JB,2)                
     &                ,ASPANG(K,JB,3),ASPANG(K,JB,4),ASPANG(K,JB,5)     
   44 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(2),RR(JC)*PVISC                       
      DO 45 K=2,KG/2+1                                                  
        WRITE(14,1100) K-1,ASPANG(K,JC,1),ASPANG(K,JC,2)                
     &                ,ASPANG(K,JC,3),ASPANG(K,JC,4),ASPANG(K,JC,5)     
   45 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(2),RR(JD)*PVISC                       
      DO 46 K=2,KG/2+1                                                  
        WRITE(14,1100) K-1,ASPANG(K,JD,1),ASPANG(K,JD,2)                
     &                ,ASPANG(K,JD,3),ASPANG(K,JD,4),ASPANG(K,JD,5)     
   46 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(2),RR(JE)*PVISC                       
      DO 47 K=2,KG/2+1                                                  
        WRITE(14,1100) K-1,ASPANG(K,JE,1),ASPANG(K,JE,2)                
     &                ,ASPANG(K,JE,3),ASPANG(K,JE,4),ASPANG(K,JE,5)     
   47 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(2),RR(JF)*PVISC                       
      DO 48 K=2,KG/2+1                                                  
        WRITE(14,1100) K-1,ASPANG(K,JF,1),ASPANG(K,JF,2)                
     &                ,ASPANG(K,JF,3),ASPANG(K,JF,4),ASPANG(K,JF,5)     
   48 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(2),RR(JH)*PVISC                       
      DO 49 K=2,KG/2+1                                                  
        WRITE(14,1100) K-1,ASPANG(K,JH,1),ASPANG(K,JH,2)                
     &                ,ASPANG(K,JH,3),ASPANG(K,JH,4),ASPANG(K,JH,5)     
   49 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(2),RR(JI)*PVISC                       
      DO 50 K=2,KG/2+1                                                  
        WRITE(14,1100) K-1,ASPANG(K,JI,1),ASPANG(K,JI,2)                
     &                ,ASPANG(K,JI,3),ASPANG(K,JI,4),ASPANG(K,JI,5)     
   50 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(2),RR(JJ)*PVISC                       
      DO 51 K=2,KG/2+1                                                  
        WRITE(14,1100) K-1,ASPANG(K,JJ,1),ASPANG(K,JJ,2)                
     &                ,ASPANG(K,JJ,3),ASPANG(K,JJ,4),ASPANG(K,JJ,5)     
   51 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(2),RR(JK)*PVISC                       
      DO 52 K=2,KG/2+1                                                  
        WRITE(14,1100) K-1,ASPANG(K,JK,1),ASPANG(K,JK,2)                
     &                ,ASPANG(K,JK,3),ASPANG(K,JK,4),ASPANG(K,JK,5)     
   52 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A52,F8.3)') CMSP(2),RR(JL)*PVISC                       
      DO 53 K=2,KG/2+1                                                  
        WRITE(14,1100) K-1,ASPANG(K,JL,1),ASPANG(K,JL,2)                
     &                ,ASPANG(K,JL,3),ASPANG(K,JL,4),ASPANG(K,JL,5)     
   53 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(1),RR(JA)*PVISC                       
      DO 59 I=1,IG/2+1                                                  
        WRITE(14,1080) DBLE(I-1)*DX*RE,ATPC1G(I,JA,1),ATPC1G(I,JA,2)    
     &                ,ATPC1G(I,JA,3),ATPC1G(I,JA,4)                    
   59 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(1),RR(JB)*PVISC                       
      DO 60 I=1,IG/2+1                                                  
        WRITE(14,1080) DBLE(I-1)*DX*RE,ATPC1G(I,JB,1),ATPC1G(I,JB,2)    
     &                ,ATPC1G(I,JB,3),ATPC1G(I,JB,4)                    
   60 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(1),RR(JC)*PVISC                       
      DO 61 I=1,IG/2+1                                                  
        WRITE(14,1080) DBLE(I-1)*DX*RE,ATPC1G(I,JC,1),ATPC1G(I,JC,2)    
     &                ,ATPC1G(I,JC,3),ATPC1G(I,JC,4)                    
   61 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(1),RR(JD)*PVISC                       
      DO 62 I=1,IG/2+1                                                  
        WRITE(14,1080) DBLE(I-1)*DX*RE,ATPC1G(I,JD,1),ATPC1G(I,JD,2)    
     &                ,ATPC1G(I,JD,3),ATPC1G(I,JD,4)                    
   62 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(1),RR(JE)*PVISC                       
      DO 63 I=1,IG/2+1                                                  
        WRITE(14,1080) DBLE(I-1)*DX*RE,ATPC1G(I,JE,1),ATPC1G(I,JE,2)    
     &                ,ATPC1G(I,JE,3),ATPC1G(I,JE,4)                    
   63 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(1),RR(JF)*PVISC                       
      DO 64 I=1,IG/2+1                                                  
        WRITE(14,1080) DBLE(I-1)*DX*RE,ATPC1G(I,JF,1),ATPC1G(I,JF,2)    
     &                ,ATPC1G(I,JF,3),ATPC1G(I,JF,4)                    
   64 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(1),RR(JH)*PVISC                       
      DO 65 I=1,IG/2+1                                                  
        WRITE(14,1080) DBLE(I-1)*DX*RE,ATPC1G(I,JH,1),ATPC1G(I,JH,2)    
     &                ,ATPC1G(I,JH,3),ATPC1G(I,JH,4)                    
   65 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(1),RR(JI)*PVISC                       
      DO 66 I=1,IG/2+1                                                  
        WRITE(14,1080) DBLE(I-1)*DX*RE,ATPC1G(I,JI,1),ATPC1G(I,JI,2)    
     &                ,ATPC1G(I,JI,3),ATPC1G(I,JI,4)                    
   66 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(1),RR(JJ)*PVISC                       
      DO 67 I=1,IG/2+1                                                  
        WRITE(14,1080) DBLE(I-1)*DX*RE,ATPC1G(I,JJ,1),ATPC1G(I,JJ,2)    
     &                ,ATPC1G(I,JJ,3),ATPC1G(I,JJ,4)                    
   67 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(1),RR(JK)*PVISC                       
      DO 68 I=1,IG/2+1                                                  
        WRITE(14,1080) DBLE(I-1)*DX*RE,ATPC1G(I,JK,1),ATPC1G(I,JK,2)    
     &                ,ATPC1G(I,JK,3),ATPC1G(I,JK,4)                    
   68 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(1),RR(JL)*PVISC                       
      DO 69 I=1,IG/2+1                                                  
        WRITE(14,1080) DBLE(I-1)*DX*RE,ATPC1G(I,JL,1),ATPC1G(I,JL,2)    
     &                ,ATPC1G(I,JL,3),ATPC1G(I,JL,4)                    
   69 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(2),RR(JA)*PVISC                       
      DO 75 K=1,KG/2+1                                                  
        WRITE(14,1080) DBLE(K-1)*DZ*RE,ATPC3G(K,JA,1),ATPC3G(K,JA,2)    
     &                ,ATPC3G(K,JA,3),ATPC3G(K,JA,4)                    
   75 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(2),RR(JB)*PVISC                       
      DO 76 K=1,KG/2+1                                                  
        WRITE(14,1080) DBLE(K-1)*DZ*RE,ATPC3G(K,JB,1),ATPC3G(K,JB,2)    
     &                ,ATPC3G(K,JB,3),ATPC3G(K,JB,4)                    
   76 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(2),RR(JC)*PVISC                       
      DO 77 K=1,KG/2+1                                                  
        WRITE(14,1080) DBLE(K-1)*DZ*RE,ATPC3G(K,JC,1),ATPC3G(K,JC,2)    
     &                ,ATPC3G(K,JC,3),ATPC3G(K,JC,4)                    
   77 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(2),RR(JD)*PVISC                       
      DO 78 K=1,KG/2+1                                                  
        WRITE(14,1080) DBLE(K-1)*DZ*RE,ATPC3G(K,JD,1),ATPC3G(K,JD,2)    
     &                ,ATPC3G(K,JD,3),ATPC3G(K,JD,4)                    
   78 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(2),RR(JE)*PVISC                       
      DO 79 K=1,KG/2+1                                                  
        WRITE(14,1080) DBLE(K-1)*DZ*RE,ATPC3G(K,JE,1),ATPC3G(K,JE,2)    
     &                ,ATPC3G(K,JE,3),ATPC3G(K,JE,4)                    
   79 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(2),RR(JF)*PVISC                       
      DO 80 K=1,KG/2+1                                                  
        WRITE(14,1080) DBLE(K-1)*DZ*RE,ATPC3G(K,JF,1),ATPC3G(K,JF,2)    
     &                ,ATPC3G(K,JF,3),ATPC3G(K,JF,4)                    
   80 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(2),RR(JH)*PVISC                       
      DO 81 K=1,KG/2+1                                                  
        WRITE(14,1080) DBLE(K-1)*DZ*RE,ATPC3G(K,JH,1),ATPC3G(K,JH,2)    
     &                ,ATPC3G(K,JH,3),ATPC3G(K,JH,4)                    
   81 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(2),RR(JI)*PVISC                       
      DO 82 K=1,KG/2+1                                                  
        WRITE(14,1080) DBLE(K-1)*DZ*RE,ATPC3G(K,JI,1),ATPC3G(K,JI,2)    
     &                ,ATPC3G(K,JI,3),ATPC3G(K,JI,4)                    
   82 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(2),RR(JJ)*PVISC                       
      DO 83 K=1,KG/2+1                                                  
        WRITE(14,1080) DBLE(K-1)*DZ*RE,ATPC3G(K,JJ,1),ATPC3G(K,JJ,2)    
     &                ,ATPC3G(K,JJ,3),ATPC3G(K,JJ,4)                    
   83 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(2),RR(JK)*PVISC                       
      DO 84 K=1,KG/2+1                                                  
        WRITE(14,1080) DBLE(K-1)*DZ*RE,ATPC3G(K,JK,1),ATPC3G(K,JK,2)    
     &                ,ATPC3G(K,JK,3),ATPC3G(K,JK,4)                    
   84 CONTINUE                                                          
C                                                                       
      WRITE(14,'(A47,F8.3)') CMCR(2),RR(JL)*PVISC                       
      DO 85 K=1,KG/2+1                                                  
        WRITE(14,1080) DBLE(K-1)*DZ*RE,ATPC3G(K,JL,1),ATPC3G(K,JL,2)    
     &                ,ATPC3G(K,JL,3),ATPC3G(K,JL,4)                    
   85 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(22)                                           
      DO 91 J=1,JG                                                      
        WRITE(14,1080) YY(J),SKEW(J,1),SKEW(J,2),SKEW(J,3)              
     &                ,SKEW(J,4)                                        
   91 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(23)                                           
      DO 92 J=1,JG                                                      
        WRITE(14,1080) YY(J),FLAT(J,1),FLAT(J,2),FLAT(J,3)              
     &                ,FLAT(J,4)                                        
   92 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(24)                                           
      DO 703 J=1,JG                                                     
        WRITE(14,1070) YY(J),OMGG(J,1),OMGG(J,2),OMGG(J,3)              
  703 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(25)                                           
      DO 94 J=1,JG                                                      
        WRITE(14,1110) Y(J),AAVG(J),AVRMSG(J),AU22G(J),SKEW(J,2)        
     &                ,FLAT(J,2)                                        
   94 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(26)                                           
      DO 95 J=1,JG                                                      
        WRITE(14,1070) YY(J),AUTRIG(J,1),AUTRIG(J,2),AUTRIG(J,3)        
   95 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(27)                                           
      DO 96 J=1,JG                                                      
        WRITE(14,1070) YY(J),AUTRIG(J,4),AUTRIG(J,5),AUTRIG(J,6)        
   96 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(28)                                           
      DO 97 J=1,JG                                                      
        WRITE(14,1090) Y(J),SKEWUVG(J),FLATUVG(J)                       
   97 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(32)                                           
      DO 98 J=1,JG                                                      
        WRITE(14,1080) YY(J),AUQUADG(J,1),AUQUADG(J,2)                  
     &                ,AUQUADG(J,4),AUQUADG(J,5)                        
   98 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(33)                                           
      DO 99 J=1,JG                                                      
        WRITE(14,1070) YY(J),AUQUADG(J,6),AUQUADG(J,7)                  
     &                ,AUQUADG(J,8)                                     
   99 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(34)                                           
      DO 200 J=1,JG                                                     
        WRITE(14,1070) Y(J),AU22G(J),AUTRIG(J,3),AUQUADG(J,3)           
  200 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(35)                                           
      DO 201 J=1,JG                                                     
        WRITE(14,1080) Y(J),AQAUVG(J,1),AQAUVG(J,2)                     
     &                ,AQAUVG(J,3),AQAUVG(J,4)                          
  201 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(36)                                           
      DO 202 J=1,JG                                                     
        WRITE(14,1080) Y(J),ANQAUVG(J,1),ANQAUVG(J,2)                   
     &                ,ANQAUVG(J,3),ANQAUVG(J,4)                        
  202 CONTINUE                                                          
C                                                                       
      WRITE(14,1010) CPRM(37)                                           
      DO 203 J=1,JG                                                     
        WRITE(14,1090) YY(J),AAINV2G(J),AAINV3G(J)                      
  203 CONTINUE                                                          
C                                                                       
 1060 FORMAT(I3,4E14.6)                                                 
 1070 FORMAT(4E14.6)                                                    
 1080 FORMAT(5E14.6)                                                    
 1090 FORMAT(3E14.6)                                                    
 1100 FORMAT(I3,5E14.6)                                                 
 1110 FORMAT(6E14.6)                                                    
C                                                                       
C      CLOSE(14)                                                        
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C     SEED DATA (VELOCITY FIELD) : INPUT & OUTPUT                       
C====================================================================   
C                                                                       
      SUBROUTINE SEEDIN(IERR)                                           
C                                                                       
C====================================================================   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
      INCLUDE './INCFILE/INCDATA'                                         
      INCLUDE './INCFILE/INCMESH'                                         
C                                                                       
      CHARACTER*71 CMTD                                                 
C                                                                       
      DIMENSION U(-4:IG+4,-2:JG+2,-4:KG+4,3)                            
      COMMON /VELO/ UG(-4:IG+4,-2:JG+2,-4:KG+4,3)                       
C                                                                       
      DIMENSION UASDIN(-4:IG/2,-2:JG+2,-4:KG+4)                         
      DIMENSION UBSDIN( IG/2+4,-2:JG+2,-4:KG+4)                         
      DIMENSION VASDIN(-4:IG/2,-2:JG+2,-4:KG+4)                         
      DIMENSION VBSDIN( IG/2+4,-2:JG+2,-4:KG+4)                         
      DIMENSION WASDIN(-4:IG/2,-2:JG+2,-4:KG+4)                         
      DIMENSION WBSDIN( IG/2+4,-2:JG+2,-4:KG+4)                         
      DIMENSION UASDING(-4:IG/2,-2:JG+2,-4:KG+4)                        
      DIMENSION UBSDING( IG/2+4,-2:JG+2,-4:KG+4)                        
      DIMENSION VASDING(-4:IG/2,-2:JG+2,-4:KG+4)                        
      DIMENSION VBSDING( IG/2+4,-2:JG+2,-4:KG+4)                        
      DIMENSION WASDING(-4:IG/2,-2:JG+2,-4:KG+4)                        
      DIMENSION WBSDING( IG/2+4,-2:JG+2,-4:KG+4)                        
C                                                                       
C----COMBINATION------------------------------------------------        
C                                                                       
      EQUIVALENCE (UG,U)                                                
      EQUIVALENCE (UASDING,UASDIN),(UBSDING,UBSDIN)                     
     &           ,(VASDING,VASDIN),(VBSDING,VBSDIN)                     
     &           ,(WASDING,WASDIN),(WBSDING,WBSDIN)                     
C                                                                       
C=====================================================================  
C------------------------------- DATAIN                                 
C      OPEN(20)                                                         
      READ(20) NUM,TIME                                                 
      READ(20) LKGD,JGD,LIGD,DT,RE,TI                                   
      READ(20) XL,ZL,NDRAW,NFIN                                         
      READ(20) DZ,DX,R1,R2,ALP                                          
 1010 FORMAT(A71)                                                       
C                                                                       
      IERR=0                                                            
      IF ((LIGD.NE.LIG).OR.(JGD.NE.JG).OR.(LKGD.NE.LKG)) THEN           
        IERR=1000                                                       
        CMTD=' MESH DATA ERROR'                                         
        WRITE(6,1010) CMTD                                              
        GOTO 1                                                          
      ENDIF                                                             
C                                                                       
      READ(20) R,DR,RR,DRR                                              
      READ(20) UASDING                                                  
      READ(21) UBSDING                                                  
      READ(22) VASDING                                                  
      READ(23) VBSDING                                                  
      READ(24) WASDING                                                  
      READ(25) WBSDING                                                  
C      CLOSE(20)                                                        
C                                                                       
      DO 10 J= -2, JG+2                                                 
      DO 10 I= -4, IG/2                                                 
      DO 10 K= -4, KG+4                                                 
        U(I,J,K,1)=UASDIN(I,J,K)                                        
        U(I,J,K,2)=VASDIN(I,J,K)                                        
        U(I,J,K,3)=WASDIN(I,J,K)                                        
   10 CONTINUE                                                          
C                                                                       
      DO 11 J= -2, JG+2                                                 
      DO 11 I= IG/2+1, IG+4                                             
      DO 11 K= -4, KG+4                                                 
        U(I,J,K,1)=UBSDIN(I-IG/2,J,K)                                   
        U(I,J,K,2)=VBSDIN(I-IG/2,J,K)                                   
        U(I,J,K,3)=WBSDIN(I-IG/2,J,K)                                   
   11 CONTINUE                                                          
C                                                                       
    1 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C                                                                       
      SUBROUTINE SEEDOUT(CODE)                                          
C                                                                       
C====================================================================   
C                                                                       
      INCLUDE './INCFILE/INCINIT96'                                       
      INCLUDE './INCFILE/INCDATA'                                         
      INCLUDE './INCFILE/INCMESH'                                         
C                                                                       
      DIMENSION U(-4:IG+4,-2:JG+2,-4:KG+4,3)                            
      COMMON /VELO/ UG(-4:IG+4,-2:JG+2,-4:KG+4,3)                       
C                                                                       
      DIMENSION UASDOUT(-4:IG/2,-2:JG+2,-4:KG+4)                        
      DIMENSION UBSDOUT( IG/2+4,-2:JG+2,-4:KG+4)                        
      DIMENSION VASDOUT(-4:IG/2,-2:JG+2,-4:KG+4)                        
      DIMENSION VBSDOUT( IG/2+4,-2:JG+2,-4:KG+4)                        
      DIMENSION WASDOUT(-4:IG/2,-2:JG+2,-4:KG+4)                        
      DIMENSION WBSDOUT( IG/2+4,-2:JG+2,-4:KG+4)                        
      DIMENSION UASDOUTG(-4:IG/2,-2:JG+2,-4:KG+4)                       
      DIMENSION UBSDOUTG( IG/2+4,-2:JG+2,-4:KG+4)                       
      DIMENSION VASDOUTG(-4:IG/2,-2:JG+2,-4:KG+4)                       
      DIMENSION VBSDOUTG( IG/2+4,-2:JG+2,-4:KG+4)                       
      DIMENSION WASDOUTG(-4:IG/2,-2:JG+2,-4:KG+4)                       
      DIMENSION WBSDOUTG( IG/2+4,-2:JG+2,-4:KG+4)                       
C                                                                       
C----COMBINATION------------------------------------------------        
C                                                                       
      EQUIVALENCE (UG,U)                                                
      EQUIVALENCE (UASDOUTG,UASDOUT),(UBSDOUTG,UBSDOUT)                 
     &           ,(VASDOUTG,VASDOUT),(VBSDOUTG,VBSDOUT)                 
     &           ,(WASDOUTG,WASDOUT),(WBSDOUTG,WBSDOUT)                 
C                                                                       
C=====================================================================  
C                                                                       
      DO 10 J= -2, JG+2                                                 
      DO 10 I= -4, IG/2                                                 
      DO 10 K= -4, KG+4                                                 
        UASDOUT(I,J,K)=U(I,J,K,1)                                       
        VASDOUT(I,J,K)=U(I,J,K,2)                                       
        WASDOUT(I,J,K)=U(I,J,K,3)                                       
   10 CONTINUE                                                          
C                                                                       
      DO 11 J= -2, JG+2                                                 
      DO 11 I= IG/2+1, IG+4                                             
      DO 11 K= -4, KG+4                                                 
        UBSDOUT(I-IG/2,J,K)=U(I,J,K,1)                                  
        VBSDOUT(I-IG/2,J,K)=U(I,J,K,2)                                  
        WBSDOUT(I-IG/2,J,K)=U(I,J,K,3)                                  
   11 CONTINUE                                                          
C                                                                       
C------------------------------- DATAOUT                                
C      OPEN(26)                                                         
      WRITE(26) NUM,TIME                                                
      WRITE(26) LKG,JG,LIG,DT,RE,TI                                     
      WRITE(26) XL,ZL,NDRAW,NFIN                                        
      WRITE(26) DZ,DX,R1,R2,ALP                                         
      WRITE(26) R,DR,RR,DRR                                             
      WRITE(26) UASDOUTG                                                
      WRITE(27) UBSDOUTG                                                
      WRITE(28) VASDOUTG                                                
      WRITE(29) VBSDOUTG                                                
      WRITE(30) WASDOUTG                                                
      WRITE(31) WBSDOUTG                                                
C                                                                       
C      CLOSE(26)                                                        
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C====================================================================   
C
      INCLUDE './INCFILE/TFIELD_CHOMP'
      INCLUDE './INCFILE/INCTMSRS'