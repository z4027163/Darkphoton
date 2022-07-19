C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C *****************************************************************
C ************* SUBROUTINE FOR THE SUSY COUPLINGS *****************
C *****************************************************************
      SUBROUTINE SUSYCP(TGBET)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      DOUBLE PRECISION LA1,LA2,LA3,LA4,LA5,LA6,LA7,LA3T
      COMPLEX*16 YF0
      DIMENSION MST(2),MSB(2),MSL(2),MSU(2),MSD(2),MSE(2),MSN(2),
     .          GCEN(2,2),GCTB(2,2),GLEE(2,2),GLTT(2,2),GLBB(2,2),
     .          GHEE(2,2),GHTT(2,2),GHBB(2,2)
      COMMON/FLAG/IPOLE
      COMMON/MODEL/IMODEL
      COMMON/MASSES/AMS,AMC,AMB,AMT
      COMMON/ALS/XLAMBDA,AMCA,AMBA,AMTA,N0A
      COMMON/HMASS/AMSM,AMA,AML,AMH,AMCH,AMAR
      COMMON/CHIMASS/AMCHI
      COMMON/HSELF/LA1,LA2,LA3,LA4,LA5,LA6,LA7
      COMMON/BREAK/AMSQ,AMUR,AMDR,AU,AD,AMU,AM2
      COMMON/BREAKGLU/AMGLU
      COMMON/PARAM/GF,AMW,AMZ
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      COMMON/SFER1ST/MQL1,MUR1,MDR1,MEL1,MER1
      COMMON/GLUINO/AMGLUINO,AMSB1,AMSB2,STHB,CTHB,
     .              XLBB(2,2),XHBB(2,2),XABB(2,2),
     .              AMST1,AMST2,STHT,CTHT,
     .              XLTT(2,2),XHTT(2,2),XATT(2,2)
      AMC0 = AMAC
      AMB0 = AMAB
      AMT0 = AMAT
      AMAC = AMC
      AMAB = AMB
      AMAT = AMT
      PI=4*DATAN(1D0)
      V=1.D0/DSQRT(DSQRT(2.D0)*GF)
      BET=DATAN(TGBET)
      SB = DSIN(BET)
      CB = DCOS(BET)
      AMAR = AMA
C  ============ HEAVIEST CHARGINO MASS NEEDED FOR SUBH ========== 
      AMCHI2=AM2**2+AMU**2+2.D0*AMW**2+DSQRT((AM2**2-AMU**2)**2
     .      +4.D0*AMW**4*DCOS(2.D0*BET)**2+4.D0*AMW**2*
     .      (AM2**2+AMU**2+2.D0*AMU*AM2*DSIN(2.D0*BET) ) ) 
      AMCHI=DSQRT(0.5D0*AMCHI2)
C ===============================================================
C ========== RUNNING MASSES
      IF(IMODEL.EQ.1)THEN
       CALL SUBH1(AMA,TGBET,AMSQ,AMUR,AMDR,AMT,AU,AD,AMU,AMCHI,
     .            AMLR,AMHR,AMCH,SA,CA,TANBA,AMGLU)
      ELSEIF(IMODEL.EQ.2)THEN
       CALL SUBH2(AMA,TGBET,AMSQ,AMUR,AMT,AU,AD,AMU,
     .            AMLR,AMHR,AMCH,SA,CA,TANBA)
      ELSEIF(IMODEL.EQ.3)THEN
C--Use Carena et al. for everything not included in FeynHiggs....
       B=DATAN(TGBET)
       MQL1 = AMSQ
       MUR1 = AMSQ
       MDR1 = AMSQ
       MEL1 = AMSQ
       MER1 = AMSQ
       AMEL = AMSQ
       AMER = AMSQ
       AL = 0
c      TSC = (AMSQ+AMUR+AMDR)/3
c      BSC = (AMSQ+AMUR+AMDR)/3
       TSC = AMT
       BSC = AMT
       CALL SFERMION(TSC,BSC,AMSQ,AMUR,AMDR,AMEL,AMER,AL,AU,AD,AMU,
     .                    MST,MSB,MSL,MSU,MSD,MSE,MSN, 
     .                    GLEE,GLTT,GLBB,GHEE,GHTT,GHBB,
     .                    GAEE,GATT,GABB,GCEN,GCTB)
       CALL SUBH2(AMA,TGBET,AMSQ,AMUR,AMT,AU,AD,AMU,
     .            AMLR,AMHR,AMCH,SA,CA,TANBA)
       CALL FEYNHIGGS(AMA,TGBET,AMT,AMST1,AMST2,STHT,AMSB1,
     .                AMSB2,STHB,AMU,AMGLU,AM2,AMLR,AMHR,SA,CA)
      ENDIF
      LA3T=LA3+LA4+LA5
      AMA2=AMAR**2
      AML2=AMLR**2
      AMH2=AMHR**2
      AMP2=AMCH**2
C ========== HIGGS COUPLINGS 
      SBMA = SB*CA-CB*SA
      CBMA = CB*CA+SB*SA
      SBPA = SB*CA+CB*SA
      CBPA = CB*CA-SB*SA
      S2A = 2*SA*CA
      C2A = CA**2-SA**2
      S2B = 2*SB*CB
      C2B = CB**2-SB**2
      GLZZ = 1/V/2*AML2*SBMA
      GHZZ = 1/V/2*AMH2*CBMA
      GLWW = 2*GLZZ
      GHWW = 2*GHZZ
      GLAZ = 1/V*(AML2-AMA2)*CBMA
      GHAZ = -1/V*(AMH2-AMA2)*SBMA
      GLPW = -1/V*(AMP2-AML2)*CBMA
      GLMW = GLPW
      GHPW = 1/V*(AMP2-AMH2)*SBMA
      GHMW = GHPW
      GAPW = 1/V*(AMP2-AMA2)
      GAMW = -GAPW
      GHHH = V/2*(LA1*CA**3*CB + LA2*SA**3*SB + LA3T*SA*CA*SBPA
     .     + LA6*CA**2*(3*SA*CB+CA*SB) + LA7*SA**2*(3*CA*SB+SA*CB))
      GLLL = -V/2*(LA1*SA**3*CB - LA2*CA**3*SB + LA3T*SA*CA*CBPA
     .     - LA6*SA**2*(3*CA*CB-SA*SB) + LA7*CA**2*(3*SA*SB-CA*CB))
      GLHH = -3*V/2*(LA1*CA**2*CB*SA - LA2*SA**2*SB*CA
     .     + LA3T*(SA**3*CB-CA**3*SB+2*SBMA/3)
     .     - LA6*CA*(CB*C2A-SA*SBPA) - LA7*SA*(C2A*SB+CA*SBPA))
      GHLL = 3*V/2*(LA1*SA**2*CB*CA + LA2*CA**2*SB*SA
     .     + LA3T*(SA**3*SB+CA**3*CB-2*CBMA/3)
     .     - LA6*SA*(CB*C2A+CA*CBPA) + LA7*CA*(C2A*SB+SA*CBPA))
      GLAA = -V/2*(LA1*SB**2*CB*SA - LA2*CB**2*SB*CA
     .     - LA3T*(SB**3*CA-CB**3*SA) + 2*LA5*SBMA
     .     - LA6*SB*(CB*SBPA+SA*C2B) - LA7*CB*(C2B*CA-SB*SBPA))
      GHAA = V/2*(LA1*SB**2*CB*CA + LA2*CB**2*SB*SA
     .     + LA3T*(SB**3*SA+CB**3*CA) - 2*LA5*CBMA
     .     - LA6*SB*(CB*CBPA+CA*C2B) + LA7*CB*(SB*CBPA+SA*C2B))
      GLPM = 2*GLAA + V*(LA5 - LA4)*SBMA
      GHPM = 2*GHAA + V*(LA5 - LA4)*CBMA
      GLZZ = 2*GLZZ
      GHZZ = 2*GHZZ
      GLLL = 6*GLLL
      GHHH = 6*GHHH
      GLHH = 2*GLHH
      GHLL = 2*GHLL
      GLAA = 2*GLAA
      GHAA = 2*GHAA
      XNORM = AMZ**2/V
      GLLL = GLLL/XNORM
      GHLL = GHLL/XNORM
      GLHH = GLHH/XNORM
      GHHH = GHHH/XNORM
      GHAA = GHAA/XNORM
      GLAA = GLAA/XNORM
      GLPM = GLPM/XNORM
      GHPM = GHPM/XNORM
      GAT=1.D0/TGBET
      GAB=TGBET
      GLT=CA/SB
      GLB=-SA/CB
      GHT=SA/SB
      GHB=CA/CB
      GZAL=-CBMA
      GZAH=SBMA
      GLVV=SBMA
      GHVV=CBMA
      B=BET
      A=DATAN(SA/CA)
      IF(CA.LT.0D0)THEN
       IF(SA.LT.0D0)THEN
        A = A-PI
       ELSE
        A = A+PI
       ENDIF
      ENDIF
C ===============================================================
C ========== POLE MASSES 
      IF(IMODEL.EQ.1)THEN
      IF(IPOLE.EQ.1) THEN 
       MT=RUNM(AMT,6)
       MB=RUNM(AMT,5)
       SW2=1.D0-AMW**2/AMZ**2
C===== STOP MASSES
       MSTL2=AMSQ**2+(0.5D0-2.D0/3.D0*SW2)*AMZ**2*DCOS(2.D0*B)
       MSTR2=AMUR**2+2.D0/3.D0*SW2*AMZ**2*DCOS(2.D0*B)
       MLRT=AU-AMU/TGBET
       DELT=(MSTL2-MSTR2)**2+4*MT**2*MLRT**2
       MST12=MT**2+0.5D0*(MSTL2+MSTR2-DSQRT(DELT))
       MST22=MT**2+0.5D0*(MSTL2+MSTR2+DSQRT(DELT))
        IF(MST12.LT.0.D0)GOTO 111
       MST(1)=DSQRT(MST12)
       MST(2)=DSQRT(MST22)
       IF(MSTL2.EQ.MSTR2) THEN
        THET = PI/4
       ELSE
        THET=0.5D0*DATAN(2.D0*MT*MLRT / (MSTL2-MSTR2) )
        IF(MSTL2.GT.MSTR2) THET = THET + PI/2
       ENDIF
       CST= DCOS(THET)
       SST= DSIN(THET)
C===== SBOTTOM MASSES
       MSBL2=AMSQ**2+(-0.5D0+1.D0/3.D0*SW2)*AMZ**2*DCOS(2.D0*B)
       MSBR2=AMDR**2-1.D0/3.D0*SW2*AMZ**2*DCOS(2.D0*B)
       MLRB=AD-AMU*TGBET
       DELB=(MSBL2-MSBR2)**2+4*MB**2*MLRB**2
       MSB12=MB**2+0.5D0*(MSBL2+MSBR2-DSQRT(DELB))
       MSB22=MB**2+0.5D0*(MSBL2+MSBR2+DSQRT(DELB))
        IF(MSB12.LT.0.D0)GOTO 111
       MSB(1)=DSQRT(MSB12)
       MSB(2)=DSQRT(MSB22)
       IF(MSBL2.EQ.MSBR2) THEN
        THEB = PI/4
       ELSE
        THEB=0.5D0*DATAN(2.D0*MB*MLRB / (MSBL2-MSBR2) )
        IF(MSBL2.GT.MSBR2) THEB = THEB + PI/2
       ENDIF
       CSB= DCOS(THEB)
       SSB= DSIN(THEB)
C===== LIGHT HIGGS COUPLINGS 
       GLTT(1,1)=-SBPA*(0.5D0*CST**2-2.D0/3.D0*SW2*DCOS(2*THET) )
     .     +MT**2/AMZ**2*GLT + MT*SST*CST/AMZ**2*(AU*GLT+AMU*GHT)
       GLTT(2,2)=-SBPA*(0.5D0*SST**2+2.D0/3.D0*SW2*DCOS(2*THET) )
     .     +MT**2/AMZ**2*GLT - MT*SST*CST/AMZ**2*(AU*GLT+AMU*GHT)
       GLTT(1,2)=-2*SBPA*SST*CST*(2.D0/3.D0*SW2-0.25D0)
     .     + MT*DCOS(2*THET)/2.D0/AMZ**2*(AU*GLT+AMU*GHT)
       GLTT(2,1)=-2*SBPA*SST*CST*(2.D0/3.D0*SW2-0.25D0)
     .     + MT*DCOS(2*THET)/2.D0/AMZ**2*(AU*GLT+AMU*GHT)
       GLBB(1,1)=-SBPA*(-0.5D0*CSB**2+1.D0/3.D0*SW2*DCOS(2*THEB))
     .     +MB**2/AMZ**2*GLB + MB*SSB*CSB/AMZ**2*(AD*GLB-AMU*GHB)
       GLBB(2,2)=-SBPA*(-0.5D0*SSB**2-1.D0/3.D0*SW2*DCOS(2*THEB))
     .     +MB**2/AMZ**2*GLB - MB*SSB*CSB/AMZ**2*(AD*GLB-AMU*GHB)
       GLBB(1,2)=-2*SBPA*SSB*CSB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/AMZ**2*(AD*GLB-AMU*GHB)
       GLBB(2,1)=-2*SBPA*SSB*CSB*(-1.D0/3.D0*SW2+0.25D0)
     .     + MB*DCOS(2*THEB)/2.D0/AMZ**2*(AD*GLB-AMU*GHB)
C===== HEAVY HIGGS COUPLINGS 
       GHTT(1,1)=CBPA*(0.5D0*CST**2-2.D0/3.D0*SW2*DCOS(2*THET))
     .     +MT**2/AMZ**2*GHT + MT*SST*CST/AMZ**2*(AU*GHT-AMU*GLT)
       GHTT(2,2)=CBPA*(0.5D0*SST**2+2.D0/3.D0*SW2*DCOS(2*THET))
     .     +MT**2/AMZ**2*GHT - MT*SST*CST/AMZ**2*(AU*GHT-AMU*GLT)
       GHTT(1,2)=2*CBPA*SST*CST*(2.D0/3.D0*SW2-0.25D0)
     .     +MT*DCOS(2*THET)/2.D0/AMZ**2*(AU*GHT-AMU*GLT)
       GHTT(2,1)=2*CBPA*SST*CST*(2.D0/3.D0*SW2-0.25D0)
     .     + MT*DCOS(2*THET)/2.D0/AMZ**2*(AU*GHT-AMU*GLT)
       GHBB(1,1)=CBPA*(-0.5D0*CSB**2+1.D0/3.D0*SW2*DCOS(2*THEB))
     .     +MB**2/AMZ**2*GHB + MB*SSB*CSB/AMZ**2*(AD*GHB+AMU*GLB)
       GHBB(2,2)=CBPA*(-0.5D0*SSB**2-1.D0/3.D0*SW2*DCOS(2*THEB))
     .     + MB**2/AMZ**2*GHB - MB*SSB*CSB/AMZ**2*(AD*GHB+AMU*GLB)
       GHBB(1,2)=2*CBPA*SSB*CSB*(-1.D0/3.D0*SW2+0.25D0)
     .     + MB*DCOS(2*THEB)/2.D0/AMZ**2*(AD*GHB+AMU*GLB)
       GHBB(2,1)=2*CBPA*SSB*CSB*(-1.D0/3.D0*SW2+0.25D0)
     .     + MB*DCOS(2*THEB)/2.D0/AMZ**2*(AD*GHB+AMU*GLB)
C===== PSEUDOSCALAR HIGGS COUPLINGS 
       GATT=MT/2.D0/AMZ**2*(AMU+AU*GAT) 
       GABB=MB/2.D0/AMZ**2*(AMU+AD*GAB) 
C======= LOOP CORRECTIONS  
       XDLT=GF/(2.D0*DSQRT(2.D0)*PI**2)*GLT**2*(-2.D0*MT**2+0.5D0*AML2)
     .     *DREAL(YF0(MT,MT,AML2))
     .     *3*MT**2
       XDLB=GF/(2.D0*DSQRT(2.D0)*PI**2)*GLB**2*(-2.D0*MB**2+0.5D0*AML2)
     .     *DREAL(YF0(MB,MB,AML2))
     .     *3*MB**2
C--BUG IN CARENA ET AL. FIXED
     .     +GF/(2.D0*DSQRT(2.D0)*PI**2)*GLB**2*(0.5D0*AML2)
     .     *DLOG(MB**2/MT**2)
     .     *3*MB**2
       XDHT=GF/(2.D0*DSQRT(2.D0)*PI**2)*GHT**2*(-2.D0*MT**2+0.5D0*AMH2)
     .     *DREAL(YF0(MT,MT,AMH2))
     .     *3*MT**2
       XDHB=GF/(2.D0*DSQRT(2.D0)*PI**2)*GHB**2*(-2.D0*MB**2+0.5D0*AMH2)
     .     *DREAL(YF0(MB,MB,AMH2))
     .     *3*MB**2
C--BUG IN CARENA ET AL. FIXED
     .     +GF/(2.D0*DSQRT(2.D0)*PI**2)*GHB**2*(0.5D0*AMH2)
     .     *DLOG(MB**2/MT**2)
     .     *3*MB**2
       XDAT=GF/(2.D0*DSQRT(2.D0)*PI**2)*GAT**2*(-0.5D0*AMA2)
     .     *DREAL(YF0(MT,MT,AMA2))
     .     *3*MT**2
       XDAB=GF/(2.D0*DSQRT(2.D0)*PI**2)*GAB**2*(-0.5D0*AMA2)
     .     *DREAL(YF0(MB,MB,AMA2))
     .     *3*MB**2
C--BUG IN CARENA ET AL. FIXED
     .     +GF/(2.D0*DSQRT(2.D0)*PI**2)*GAB**2*(-0.5D0*AMA2)
     .     *DLOG(MB**2/MT**2)
     .     *3*MB**2
       XDLST=0.D0
       XDLSB=0.D0
       XDHST=0.D0
       XDHSB=0.D0
         DO 311 I=1,2
         DO 311 J=1,2
       XDLST=XDLST+GF/(2.D0*DSQRT(2.D0)*PI**2)*GLTT(I,J)**2*
     .       DREAL(YF0(MST(I),MST(J),AML2))
     .     *3*AMZ**4
       XDLSB=XDLSB+GF/(2.D0*DSQRT(2.D0)*PI**2)*GLBB(I,J)**2*
     .       DREAL(YF0(MSB(I),MSB(J),AML2))
     .    *3*AMZ**4
       XDHST=XDHST+GF/(2.D0*DSQRT(2.D0)*PI**2)*GHTT(I,J)**2*
     .       DREAL(YF0(MST(I),MST(J),AMH2))
     .     *3*AMZ**4
       XDHSB=XDHSB+GF/(2.D0*DSQRT(2.D0)*PI**2)*GHBB(I,J)**2*
     .       DREAL(YF0(MSB(I),MSB(J),AMH2))
     .     *3*AMZ**4
311    CONTINUE
       XDAST=GF/(1.D0*DSQRT(2.D0)*PI**2)*GATT**2*
     .       DREAL(YF0(MST(1),MST(2),AMA2))
     .     *3*AMZ**4
       XDASB=GF/(1.D0*DSQRT(2.D0)*PI**2)*GABB**2*
     .       DREAL(YF0(MSB(1),MSB(2),AMA2))
     .     *3*AMZ**4
      
       AML=DSQRT(AML2+XDLT+XDLB+XDLST+XDLSB)
       AMH=DSQRT(AMH2+XDHT+XDHB+XDHST+XDHSB)  
       AMA=DSQRT(AMA2+XDAT+XDAB+XDAST+XDASB)  
      ELSE
       AML=AMLR
       AMH=AMHR     
       AMA=AMAR     
      ENDIF 
      ELSE
       AML=AMLR
       AMH=AMHR
       AMA=AMAR
      ENDIF
      AMAC = AMC0
      AMAB = AMB0
      AMAT = AMT0
      RETURN
111   STOP
      END

C ===================== THE FUNCTION F0 ===============
      COMPLEX*16 FUNCTION YF0(M1,M2,QSQ)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 CD,CR,CQ2,IEPS,CBET,CXX
      M1SQ = M1*M1
      M2SQ = M2*M2
      AQSQ = DABS(QSQ)
      IEPS = DCMPLX(1.D0,1.D-12)
      CQ2 = QSQ*IEPS
      CD = (M1SQ-M2SQ)/CQ2
      CR = CDSQRT((1+CD)**2 - 4*M1SQ/CQ2)
      IF(QSQ.EQ.0.D0) THEN
       YF0 = 0.D0
      ELSE
       IF(M1.EQ.M2) THEN
        YF0 = -2.D0 + CR*CDLOG(-(1+CR)/(1-CR))
       ELSE
        CBET = CDSQRT(1-4*M1*M2/(CQ2 - (M1-M2)**2))
        CXX = (CBET-1)/(CBET+1)
        YF0 = -1 + ((QSQ+M2SQ-M1SQ)/2/QSQ - M2SQ/(M2SQ-M1SQ))
     .                                           *DLOG(M2SQ/M1SQ)
     .     - (QSQ-(M1-M2)**2)/QSQ*CBET*CDLOG(CXX)
       ENDIF
      ENDIF
      RETURN
      END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     THIS PROGRAM COMPUTES THE RENORMALIZATION GROUP IMPROVED
C     VALUES OF HIGGS MASSES AND COUPLINGS IN THE MSSM.
C
C     INPUT: MA,TANB = TAN(BETA),MQ,MUR,MDR,MTOP,AU,AD,MU,MCHI
C
C     ALL MASSES IN GEV UNITS. MA IS THE CP-ODD HIGGS MASS,
C     MTOP IS THE PHYSICAL TOP MASS, MQ AND MUR/MDR ARE THE SOFT
C     SUPERSYMMETRY BREAKING MASS PARAMETERS OF LEFT HANDED
C     AND RIGHT HANDED STOPS RESPECTIVELY, AU AND AD ARE THE
C     STOP AND SBOTTOM TRILINEAR SOFT BREAKING TERMS,
C     RESPECTIVELY,  AND MU IS THE SUPERSYMMETRIC
C     HIGGS MASS PARAMETER. WE USE THE  CONVENTIONS FROM
C     THE PHYSICS REPORT OF HABER AND KANE: LEFT RIGHT
C     STOP MIXING TERM PROPORTIONAL TO (AU - MU/TANB).
C     MCHI IS THE HEAVIEST CHARGINO MASS. 
C     WE USE AS INPUT TANB DEFINED AT THE SCALE MTOP.

C     OUTPUT: MH,HM,MCH, SA = SIN(ALPHA), CA= COS(ALPHA), TANBA
C     WHERE MHP AND HPM ARE THE LIGHTEST AND HEAVIEST CP-EVEN
C     HIGGS MASSES, MHCH IS THE CHARGED HIGGS MASS AND
C     ALPHA IS THE HIGGS MIXING ANGLE.
C     TANBA IS THE ANGLE TANB AT THE CP-ODD HIGGS MASS SCALE.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Program based on the work by M. Carena, M. Quiros
c       and C.E.M. Wagner, "Effective potential methods and
c       the Higgs mass spectrum in the MSSM", Nucl. Phys.
c       B461 (1996) 407. 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SUBH1(MA,TANB,MQ,MUR,MD,MTOP,AU,AD,MU,MCHI0,
     *                 MHP,HMP,MCH,SA,CA,TANBA,MGLU)

      IMPLICIT REAL*8(A-H,L,M,O-Z)
      DIMENSION VH(2,2),M2(2,2),M2P(2,2)
      COMMON/PARAM/GF,AMW,AMZ
      COMMON/HSELF/LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,LAMBDA5,
     .             LAMBDA6,LAMBDA7

      MCHI = MCHI0
      TANBA = TANB
      TANBT = TANB
      
      PI = 4*DATAN(1D0)
      MZ = AMZ
      MW = AMW
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = ALPHAS(AMZ,2)
      ALPHA3  = ALPHAS(MTOP,2)
      MB      = RUNM(MTOP,5)
      RMTOP   = RUNM(MTOP,6)

      TQ = LOG((MQ**2+MTOP**2)/MTOP**2)
      TU = LOG((MUR**2 + MTOP**2)/MTOP**2)
      TD = LOG((MD**2 + MTOP**2)/MTOP**2)
      SINB = TANB/DSQRT(1.D0 + TANB**2)
      COSB = SINB/TANB

      IF(MA.GT.MTOP)
     *       TANBA = TANB*(1.D0-3.D0/32.D0/PI**2*
     *       (RMTOP**2/V**2/SINB**2-MB**2/V**2/COSB**2)*
     *       DLOG(MA**2/MTOP**2))
      IF(MA.LT.MTOP.OR.MA.EQ.MTOP) TANBT = TANBA

      SINB = TANBT/DSQRT(1.D0 + TANBT**2)
      COSB = 1.D0/DSQRT(1.D0 + TANBT**2)
      COS2B = (TANBT**2 - 1.D0)/(TANBT**2 + 1.D0)
      G1 = DSQRT(ALPHA1*4.D0*PI)
      G2 = DSQRT(ALPHA2*4.D0*PI)
      G3 = DSQRT(ALPHA3*4.D0*PI)
      HU = RMTOP/V/SINB
      HD =  MB/V/COSB
C

      IF(MQ.GT.MUR) TP = TQ - TU
      IF(MQ.LT.MUR.OR.MQ.EQ.MUR) TP = TU - TQ
      IF(MQ.GT.MUR) TDP = TU
      IF(MQ.LT.MUR.OR.MQ.EQ.MUR) TDP = TQ
      IF(MQ.GT.MD) TPD = TQ - TD
      IF(MQ.LT.MD.OR.MQ.EQ.MD) TPD = TD - TQ
      IF(MQ.GT.MD) TDPD = TD
      IF(MQ.LT.MD.OR.MQ.EQ.MD) TDPD = TQ

      IF(MQ.GT.MD) DLAMBDA1 = 6./96./PI**2*G1**2*HD**2*TPD
      IF(MQ.LT.MD.OR.MQ.EQ.MD) DLAMBDA1 = 3./32./PI**2*
     * HD**2*(G1**2/3.+G2**2)*TPD

      IF(MQ.GT.MUR) DLAMBDA2 =12./96./PI**2*G1**2*HU**2*TP
      IF(MQ.LT.MUR.OR.MQ.EQ.MUR) DLAMBDA2 = 3./32./PI**2*
     * HU**2*(-G1**2/3.+G2**2)*TP

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  dlambdap1 and dlambdap2 are the new log corrections due to
c  the presence of the gluino mass. They are in general very small,
c  and only present if there is a hierarchy of masses between the
c  two stops.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        dlambdap2 = 0
        tglu = log(mglu**2/mtop**2)

        if(mglu.lt.mur.or.mglu.lt.mq) then
        if(mq.gt.mur.and.mglu.gt.mur) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tq**2-tglu**2)
        endif

        if(mq.gt.mur.and.mglu.lt.mur) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tq**2-tu**2)
        endif

        if(mq.gt.mur.and.mglu.eq.mur) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tq**2-tu**2)
        endif

        if(mur.gt.mq.and.mglu.gt.mq) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tu**2-tglu**2)
        endif

        if(mur.gt.mq.and.mglu.lt.mq) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tu**2-tq**2)
        endif

        if(mur.gt.mq.and.mglu.eq.mq) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tu**2-tq**2)
        endif
        endif

      DLAMBDA3 = 0.
      DLAMBDA4 = 0.

      IF(MQ.GT.MD) DLAMBDA3 = -1./32./PI**2*G1**2*HD**2*TPD
      IF(MQ.LT.MD.OR.MQ.EQ.MD) DLAMBDA3 = 3./64./PI**2*HD**2*
     *(G2**2-G1**2/3.)*TPD
      
      IF(MQ.GT.MUR) DLAMBDA3 = DLAMBDA3 - 
     *1./16./PI**2*G1**2*HU**2*TP
      IF(MQ.LT.MUR.OR.MQ.EQ.MUR) DLAMBDA3 = DLAMBDA3 + 
     * 3./64./PI**2*HU**2*(G2**2+G1**2/3.)*TP

      IF(MQ.LT.MUR) DLAMBDA4 = -3./32./PI**2*G2**2*HU**2*TP
      IF(MQ.LT.MD) DLAMBDA4 = DLAMBDA4 - 3./32./PI**2*G2**2*
     *                        HD**2*TPD
C
      LAMBDA1 = ((G1**2 + G2**2)/4.)*
     *(1.-3.*HD**2*(TPD + TDPD)/8./PI**2)
     *+(3.*HD**4./16./PI**2) *TPD*(1.   
     *+ (3.*HD**2/2. + HU**2/2.       
     *- 8.*G3**2) * (TPD + 2.*TDPD)/16./PI**2) 
     *+(3.*HD**4./8./PI**2) *TDPD*(1.  + (3.*HD**2/2. + HU**2/2.       
     *- 8.*G3**2) * TDPD/16./PI**2) + DLAMBDA1 
C
      LAMBDA2 = ((G1**2 + G2**2)/4.)*(1.-3.*HU**2*
     *(TP + TDP)/8./PI**2)
     *+(3.*HU**4./16./PI**2) *TP*(1.   
     *+ (3.*HU**2/2. + HD**2/2.       
     *- 8.*G3**2) * (TP + 2.*TDP)/16./PI**2) 
     *+(3.*HU**4./8./PI**2) *TDP*(1. + (3.*HU**2/2. + HD**2/2.       
     *- 8.*G3**2) * TDP/16./PI**2) + DLAMBDA2  + DLAMBDAP2
C
      LAMBDA3 = ((G2**2 - G1**2)/4.)*(1.-3.*
     *(HU**2)*(TP + TDP)/16./PI**2 -3.*
     *(HD**2)*(TPD + TDPD)/16./PI**2) +DLAMBDA3 
C
      LAMBDA4 = (- G2**2/2.)*(1.
     *-3.*(HU**2)*(TP + TDP)/16./PI**2
     *-3.*(HD**2)*(TPD + TDPD)/16./PI**2) +DLAMBDA4
C     
	LAMBDA5 = 0.
	LAMBDA6 = 0.
	LAMBDA7 = 0.


C
C     THIS IS THE CONTRIBUTION FROM LIGHT CHARGINOS/NEUTRALINOS
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  	 MSSUSY=DSQRT(0.5D0*(MQ**2+MUR**2)+MTOP**2)
	IF(MCHI.GT.MSSUSY)GOTO 3790
	IF(MCHI.LT.MTOP) MCHI=MTOP
	TCHAR=LOG(MSSUSY**2/MCHI**2)
	DELTAL12=(9./64./PI**2*G2**4+5./192./PI**2*G1**4)*TCHAR
	DELTAL3P4=(3./64./PI**2*G2**4+7./192./PI**2*G1**4
     *       +4./32/PI**2*G1**2*G2**2)*TCHAR
	DELTAM112=2.*DELTAL12*V**2*COSB**2
	DELTAM222=2.*DELTAL12*V**2*SINB**2
	DELTAM122=2.*DELTAL3P4*V**2*SINB*COSB
C--EXTENSION OF CARENA ET AL.: TRAFO MASS MATRIX -> LAMBDA_I
        DLAM1 = DELTAM112/2.D0/V**2/COSB**2
        DLAM2 = DELTAM222/2.D0/V**2/SINB**2
        DLAM3 = DELTAM122/2.D0/V**2/SINB/COSB
     .        *(G1**2-G2**2)/(G1**2+G2**2)
        DLAM4 = DELTAM122/2.D0/V**2/SINB/COSB
     .        *(2*G2**2)/(G1**2+G2**2)
        LAMBDA1 = LAMBDA1+DLAM1
        LAMBDA2 = LAMBDA2+DLAM2
        LAMBDA3 = LAMBDA3+DLAM3
        LAMBDA4 = LAMBDA4+DLAM4
C--END OF EXTENSION
 3790	CONTINUE
CCCCCCCCCCCCCCC    END OF CHARGINOS AND NEUTRALINOS  CCCCCCCCCCCC 


C--EXTENSION OF CARENA ET AL.: TRAFO MASS MATRIX -> LAMBDA_I
      CALL GFUN(MA,TANBA,MQ,MUR,MD,MTOP,AU,AD,MU,MGLU,
     *                 DLAM1,DLAM2,DLAM3,DLAM4,DLAM5,DLAM6,DLAM7)

      LAMBDA1 = LAMBDA1+DLAM1
      LAMBDA2 = LAMBDA2+DLAM2
      LAMBDA3 = LAMBDA3+DLAM3
      LAMBDA4 = LAMBDA4+DLAM4
      LAMBDA5 = LAMBDA5+DLAM5
      LAMBDA6 = LAMBDA6+DLAM6
      LAMBDA7 = LAMBDA7+DLAM7
      
      M2(1,1) = 2.*V**2*(LAMBDA1*COSB**2+2.*LAMBDA6*
     *COSB*SINB + LAMBDA5*SINB**2) + MA**2*SINB**2
      M2(2,2) = 2.*V**2*(LAMBDA5*COSB**2+2.*LAMBDA7*
     *COSB*SINB + LAMBDA2*SINB**2) + MA**2*COSB**2
      M2(1,2) = 2.*V**2*(LAMBDA6*COSB**2+(LAMBDA3+LAMBDA4)*
     *COSB*SINB + LAMBDA7*SINB**2) - MA**2*SINB*COSB
      M2(2,1) = M2(1,2)

      M2P(1,1) = M2(1,1)
      M2P(2,2) = M2(2,2)
      M2P(1,2) = M2(1,2)
      M2P(2,1) = M2(2,1)

C--END OF EXTENSION

      TRM2P  = M2P(1,1) + M2P(2,2)
      DETM2P = M2P(1,1)*M2P(2,2) - M2P(1,2)*M2P(2,1)

      MH2P = (TRM2P - DSQRT(TRM2P**2 - 4.D0* DETM2P))/2.D0
      HM2P = (TRM2P + DSQRT(TRM2P**2 - 4.D0* DETM2P))/2.D0
C !!!!!!!!!!!!!!!!!!!
      MCH2=MA**2+(LAMBDA5-LAMBDA4)*V**2
C !!!!!!!!!!!!!!!!!!!
      MCH=DSQRT(MCH2)
      HMP = DSQRT(HM2P) 
      IF(MH2P.LT.0.)GOTO 5555
      MHP = DSQRT(MH2P) 
C
      SIN2ALPHA = 2.*M2P(1,2)/DSQRT(TRM2P**2-4.D0*DETM2P)
      COS2ALPHA = (M2P(1,1)-M2P(2,2))/DSQRT(TRM2P**2-4.D0*DETM2P)
      IF(COS2ALPHA.GT.0.) ALPHA = DASIN(SIN2ALPHA)/2.D0
      IF(COS2ALPHA.LT.0.) ALPHA = -PI/2.D0-DASIN(SIN2ALPHA)/2.D0
      SA = DSIN(ALPHA)
      CA = DCOS(ALPHA)  
      SQBMA = (SINB*CA - COSB*SA)**2

5555  RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCC NON DEGENERATE STOP/SBOTTOM EFFECTS CCCCCCCCC
C
        SUBROUTINE GFUN(MA,TANB,MQ,MUR,MD,MTOP,AT,AB,MU,MGLU,
     *                     DLAM1,DLAM2,DLAM3,DLAM4,DLAM5,DLAM6,DLAM7)
        IMPLICIT REAL*8 (A-H,L,M,O-Z)
        DIMENSION VH(2,2),VH1(2,2),VH2(2,2),
     *            VH3T(2,2),VH3B(2,2),AL(2,2)
        COMMON/PARAM/GF,AMW,AMZ
        G(X,Y) = 2.D0 - (X+Y)/(X-Y)*DLOG(X/Y)

        IF(DABS(MU).LT.0.000001) MU = 0.000001
        MQ2   = MQ**2
        MUR2  = MUR**2
        MD2   = MD**2
        TANBA = TANB
        SINBA = TANBA/DSQRT(TANBA**2+1.D0)
        COSBA = SINBA/TANBA        
        SINB = TANB/DSQRT(TANB**2+1.D0)
        COSB = SINB/TANB

      MB = RUNM(MTOP,5)
      PI = 4*DATAN(1D0)
      MZ = AMZ
      MW = AMW
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = ALPHAS(AMZ,2)
      ALPHA3  = ALPHAS(MTOP,2)

      G1 = DSQRT(ALPHA1*4.*PI)
      G2 = DSQRT(ALPHA2*4.*PI)
      G3 = DSQRT(ALPHA3*4.*PI)
      
        IF(MQ.GT.MUR) MST = MQ
        IF(MUR.GT.MQ.OR.MUR.EQ.MQ) MST = MUR
        MSUSYT = DSQRT(MST**2  + MTOP**2)

	IF(MQ.GT.MD) MSB = MQ
	IF(MD.GT.MQ.OR.MD.EQ.MQ) MSB = MD
	MSUSYB = DSQRT(MSB**2 + MB**2)

	TT = LOG(MSUSYT**2/MTOP**2)
	TB = LOG(MSUSYB**2/MTOP**2)

        RMTOP   = RUNM(MTOP,6)

        HT = RMTOP/V/SINB
        HTST = RMTOP/V
        HB =  MB/V/COSB
        G32 = ALPHA3*4.*PI

        BT2 = -(8.*G32 - 9.*HT**2/2. - HB**2/2.)/(4.*PI)**2
	BB2 = -(8.*G32 - 9.*HB**2/2. - HT**2/2.)/(4.*PI)**2
        AL2 = 3./8./PI**2*HT**2
        BT2ST = -(8.*G32 - 9.*HTST**2/2.)/(4.*PI)**2
        ALST = 3./8./PI**2*HTST**2
        AL1 = 3./8./PI**2*HB**2

        AL(1,1) = AL1
        AL(1,2) = (AL2+AL1)/2.
        AL(2,1) = (AL2+AL1)/2.
        AL(2,2) = AL2

	IF(MA.GT.MTOP) THEN
        VI = V*(1. + 3./32./PI**2*HTST**2*LOG(MTOP**2/MA**2))
        H1I = VI*COSBA
        H2I = VI*SINBA
        H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYT**2))**.25
        H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYT**2))**.25
        H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYB**2))**.25
        H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYB**2))**.25
	ELSE
	VI =  V
	H1I = VI*COSB
	H2I = VI*SINB
        H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYT**2))**.25
        H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYT**2))**.25
        H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYB**2))**.25
        H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYB**2))**.25
	END IF

        TANBST = H2T/H1T
        SINBT = TANBST/(1.+TANBST**2)**.5
        COSBT = SINBT/TANBST

        TANBSB = H2B/H1B
        SINBB = TANBSB/(1.+TANBSB**2)**.5
        COSBB = SINBB/TANBSB

      CALL DELMB(MA,TANB,MQ,MUR,MD,AT,AB,MU,MGLU,
     .           MTOP,DELTAMT,DELTAMB,STOP12,STOP22,SBOT12,SBOT22)

        IF(STOP22.LT.0.) GOTO 4237
        IF(SBOT22.LT.0.) GOTO 4237

        STOP1 = STOP12**.5
        STOP2 = STOP22**.5
        SBOT1 = SBOT12**.5
        SBOT2 = SBOT22**.5

        mtop4 = rmtop**4.*(1.+2.*bt2*tt- al2*tt - 4.*deltamt)
c     * /(1.+deltamt)**4.
        mbot4 = mb**4.*(1.+2.*bb2*tb - al1*tb)
     * /(1.+deltamb)**4.
        MTOP2 = DSQRT(MTOP4)
        MBOT2 = DSQRT(MBOT4)
        mb = mb/(1+deltamb)

        VH1(1,1) = 1./TANBST
        VH1(2,1) = -1.
        VH1(1,2) = -1.
        VH1(2,2) = TANBST
        VH2(1,1) = TANBST
        VH2(1,2) = -1.
        VH2(2,1) = -1.
        VH2(2,2) = 1./TANBST

C CCCCCCCCCCCCCCCCCCCCCCCCCCC  D-terms CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	STW=SW

	F1T=(MQ2-MUR2)/(STOP12-STOP22)*(.5-4./3.*STW)*
     *         LOG(STOP1/STOP2)
     *        +(.5-2./3.*STW)*LOG(STOP1*STOP2/(MQ2+MTOP2))
     *        + 2./3.*STW*LOG(STOP1*STOP2/(MUR2+MTOP2))

	F1B=(MQ2-MD2)/(SBOT12-SBOT22)*(-.5+2./3.*STW)*
     *        LOG(SBOT1/SBOT2)
     *        +(-.5+1./3.*STW)*LOG(SBOT1*SBOT2/(MQ2+MBOT2))
     *        - 1./3.*STW*LOG(SBOT1*SBOT2/(MD2+MBOT2))

	F2T=1/(STOP12-STOP22)*
     *         (-.5*LOG(STOP12/STOP22)
     *        +(4./3.*STW-.5)*(MQ2-MUR2)/(STOP12-STOP22)*
     *         G(STOP12,STOP22))

	F2B=1/(SBOT12-SBOT22)*
     *         (.5*LOG(SBOT12/SBOT22)
     *        +(-2./3.*STW+.5)*(MQ2-MD2)/(SBOT12-SBOT22)*
     *        G(SBOT12,SBOT22))

C*************************************************************
C
C--EXTENSION OF CARENA ET AL.: TRAFO MASS MATRIX -> LAMBDA_I
C
C TRAFOS APPROXIMATE -> EXACT:
C
C (i)  1/M_{SUSY}^2 -> LOG(M1^2/M2^2) / (M1^2-M2^2)
C
C (ii) 1/M_{SUSY}^4 -> -6 G(M1^2,M2^2) / (M1^2-M2^2)^2
C
C Then use results of Phys. Lett. B355 (1995) 209 in order to
C obtain the results for lambda_1 - lambda_7 according to
C Nucl. Phys. B461 (1996) 407. Perform a full evolution from
C M_SUSY -> m_t for lambdas (anomalous dimensions, v_i).
C
C - ht^2*hb^2 terms neglected in lambda_3,4 (according to
C   Nucl. Phys. B461 (1996) 407)
C
C*************************************************************

        DLAM1T = MTOP4/(SINBT**4)*(MU**2/(STOP1**2
     *    -STOP2**2))**2*G(STOP12,STOP22)
     *  - MZ**2*MTOP2*MU**2/TANBST**2*F2T/COSBT**2

        DLAM1B = MBOT4/(COSBB**4)*(LOG(SBOT1**2*SBOT2**2/
     *    (MQ2+MBOT2)/(MD2+MBOT2))
     *    + 2*AB**2/(SBOT1**2-SBOT2**2)*LOG(SBOT1**2/SBOT2**2))
     *  + MBOT4/(COSBB**4)*(AB**2/
     *    (SBOT1**2-SBOT2**2))**2*G(SBOT12,SBOT22)
     *  + MZ**2*(2*MBOT2*F1B-MBOT2*AB**2*F2B)/COSBB**2

        DLAM2T = MTOP4/(SINBT**4)*(LOG(STOP1**2*STOP2**2/
     *    (MQ2+MTOP2)/(MUR2+MTOP2))
     *  + 2*AT**2/(STOP1**2-STOP2**2)*LOG(STOP1**2/STOP2**2))
     *  + MTOP4/(SINBT**4)*(AT**2/
     *    (STOP1**2-STOP2**2))**2*G(STOP12,STOP22)
     *  + MZ**2*(-2*MTOP2*F1T+MTOP2*AT**2*F2T)/SINBT**2
 
        DLAM2B = MBOT4/(COSBB**4)*MU**4/(SBOT1**2
     *    -SBOT2**2)**2*G(SBOT12,SBOT22)
     *    + MZ**2*MBOT2*MU**2*TANBSB**2*F2B/SINBB**2
 
        DLAM3T = MTOP4/(SINBT**4)*
     *    MU**2/(STOP1**2-STOP2**2)*(LOG(STOP1**2/STOP2**2)/2.D0
     *  + AT**2/(STOP1**2-STOP2**2)*G(STOP12,STOP22))
     *  + MZ**2*(MTOP2/TANBST*F1T-MTOP2*(AT**2-MU**2)/TANBST/2.*F2T)
     *    /SINBT/COSBT/2
c    *  + MTOP2*MBOT2/(SINBT**2*COSBB**2)*(
c    *    LOG(STOP1**2*STOP2**2/(MQ2+MTOP2)/(MUR2+MTOP2))
c    *  + LOG(SBOT1**2*SBOT2**2/(MQ2+MBOT2)/(MD2+MBOT2))
c    *  + ((AT+AB)**2/2-MU**2)*(
c    *      1.D0/(STOP1**2-SBOT1**2)*LOG(STOP1**2/SBOT1**2)
c    *    + 1.D0/(STOP2**2-SBOT2**2)*LOG(STOP2**2/SBOT2**2))
c    *  - (MU**2-AT*AB)**2*(
c    *    - 1.D0/(STOP1**2-SBOT1**2)**2*G(STOP12,SBOT12)
c    *    - 1.D0/(STOP2**2-SBOT2**2)**2*G(STOP22,SBOT22)))

        DLAM3B = MBOT4/(COSBB**4)*MU**2/(SBOT1**2-SBOT2**2)*(
     *    LOG(SBOT1**2/SBOT2**2)/2.D0
     *  + AB**2/(SBOT1**2-SBOT2**2)*G(SBOT12,SBOT22))
     *  + MZ**2*(-MBOT2*TANBSB*F1B+MBOT2*(AB**2-MU**2)*TANBSB/2.*F2B)
     *    /SINBB/COSBB/2

        DLAM4T = MTOP4/(SINBT**4)*
     *    MU**2/(STOP1**2-STOP2**2)*(LOG(STOP1**2/STOP2**2)/2.D0
     *  + AT**2/(STOP1**2-STOP2**2)*G(STOP12,STOP22))
     *  + MZ**2*(MTOP2/TANBST*F1T-MTOP2*(AT**2-MU**2)/TANBST/2.*F2T)
     *    /SINBT/COSBT/2
c    *  - MTOP2*MBOT2/(SINBT**2*COSBB**2)*(
c    *    LOG(STOP1**2*STOP2**2/(MQ2+MTOP2)/(MUR2+MTOP2))
c    *  + LOG(SBOT1**2*SBOT2**2/(MQ2+MBOT2)/(MD2+MBOT2))
c    *  + ((AT+AB)**2/2-MU**2)*(
c    *      1.D0/(STOP1**2-SBOT1**2)*LOG(STOP1**2/SBOT1**2)
c    *    + 1.D0/(STOP2**2-SBOT2**2)*LOG(STOP2**2/SBOT2**2))
c    *  - (MU**2-AT*AB)**2*(
c    *    - 1.D0/(STOP1**2-SBOT1**2)**2*G(STOP12,SBOT12)
c    *    - 1.D0/(STOP2**2-SBOT2**2)**2*G(STOP22,SBOT22)))

        DLAM4B = MBOT4/(COSBB**4)*MU**2/(SBOT1**2-SBOT2**2)*(
     *    LOG(SBOT1**2/SBOT2**2)/2.D0
     *  + AB**2/(SBOT1**2-SBOT2**2)*G(SBOT12,SBOT22))
     *  + MZ**2*(-MBOT2*TANBSB*F1B+MBOT2*(AB**2-MU**2)*TANBSB/2.*F2B)
     *    /SINBB/COSBB/2

        DLAM5T = MTOP4/(SINBT**4)*
     *    (MU**2*AT**2)/(STOP1**2-STOP2**2)**2*G(STOP12,STOP22)

        DLAM5B = MBOT4/(COSBB**4)*
     *    (MU**2*AB**2)/(SBOT1**2-SBOT2**2)**2*G(SBOT12,SBOT22)

        DLAM6T = MTOP4/(SINBT**4)*
     *    (-MU**3*AT)/(STOP1**2-STOP2**2)**2*G(STOP12,STOP22)
     *  + MZ**2*MTOP2*MU*AT/TANBST*F2T/(2*SINBT*COSBT)

        DLAM6B = MBOT4/(COSBB**4)*MU*AB*
     *    (-1.D0/(SBOT1**2-SBOT2**2)*LOG(SBOT1**2/SBOT2**2)
     *    -AB**2/(SBOT1**2-SBOT2**2)**2*G(SBOT12,SBOT22))
     *  - MZ**2*(-MBOT2*AB*MU*TANBSB*F2B)/(2*SINBB*COSBB)

        DLAM7T = MTOP4/(SINBT**4)*MU*AT*
     *    (-1.D0/(STOP1**2-STOP2**2)*LOG(STOP1**2/STOP2**2)
     *    -AT**2/(STOP1**2-STOP2**2)**2*G(STOP12,STOP22))
     *  - MZ**2*MTOP2*AT*MU/TANBST*F2T/(2*SINBT*COSBT)

        DLAM7B = MBOT4/(COSBB**4)*
     *    (-MU**3*AB)/(SBOT1**2-SBOT2**2)**2*G(SBOT12,SBOT22)
     *    - MZ**2*MBOT2*MU*AB*TANBSB*F2B/(2*SINBB*COSBB)

       TQ = LOG((MQ2 + MTOP2)/MTOP2)
       TU = LOG((MUR2+MTOP2)/MTOP2)
       TQD = LOG((MQ2 + MB**2)/MB**2)
       TD = LOG((MD2+MB**2)/MB**2)

        FACT = 3.D0/(16.D0*PI**2*(H1T**2+H2T**2)**2)
        FACB = 3.D0/(16.D0*PI**2*(H1B**2+H2B**2)**2)

        DLAM1 = FACT*DLAM1T*(1.-AL1*TT) + FACB*DLAM1B*(1.-AL1*TB)

        DLAM2 = FACT*DLAM2T*(1.-AL2*TT) + FACB*DLAM2B*(1.-AL2*TB)

        DLAM3 = FACT*DLAM3T*(1.-(AL1+AL2)/2*TT)
     *        + FACB*DLAM3B*(1.-(AL1+AL2)/2*TB)

        DLAM4 = FACT*DLAM4T*(1.-(AL1+AL2)/2*TT)
     *        + FACB*DLAM4B*(1.-(AL1+AL2)/2*TB)

        DLAM5 = FACT*DLAM5T*(1.-(AL1+AL2)/2*TT)
     *        + FACB*DLAM5B*(1.-(AL1+AL2)/2*TB)

        DLAM6 = FACT*DLAM6T*(1.-(3*AL1+AL2)/4*TT)
     *        + FACB*DLAM6B*(1.-(3*AL1+AL2)/4*TB)

        DLAM7 = FACT*DLAM7T*(1.-(AL1+3*AL2)/4*TT)
     *        + FACB*DLAM7B*(1.-(AL1+3*AL2)/4*TB)

        FACTOR = 1.D0
        DLAM1 = DLAM1 * FACTOR
        DLAM2 = DLAM2 * FACTOR
        DLAM3 = DLAM3 * FACTOR
        DLAM4 = DLAM4 * FACTOR
        DLAM5 = DLAM5 * FACTOR
        DLAM6 = DLAM6 * FACTOR
        DLAM7 = DLAM7 * FACTOR

C--END OF EXTENSION

        GOTO 4236
 4237   CONTINUE

        DLAM1 = -1.D+15
        DLAM2 = -1.D+15
        DLAM3 = -1.D+15
        DLAM4 = -1.D+15
        DLAM5 = -1.D+15
        DLAM6 = -1.D+15
        DLAM7 = -1.D+15

4236    RETURN
        END

      FUNCTION T(X,Y,Z)
      implicit real*8(a-h,l,m,o-z)
      if(x.eq.y) x = x - 0.00001
      if(x.eq.z) x = x - 0.00002
      if(y.eq.z) y = y - 0.00003
c       write(*,*) 'xyz',x,y,z
      T = (X**2*Y**2*log(X**2/Y**2) + X**2*Z**2*log(Z**2/X**2)
     * + Y**2*Z**2*log(Y**2/Z**2))/((X**2-Y**2)*(Y**2-Z**2)*(X**2-Z**2))
      return
      end

      SUBROUTINE DELMB(MA,TANB,MQ,MUR,MD,AT,AB,MU,MGLU,
     .           MTOP,DELTAMT,DELTAMB,STOP12,STOP22,SBOT12,SBOT22)
        IMPLICIT REAL*8 (A-H,L,M,O-Z)
        COMMON/PARAM/GF,AMW,AMZ

        IF(DABS(MU).LT.0.000001) MU = 0.000001
        MQ2   = MQ**2
        MUR2  = MUR**2
        MD2   = MD**2
        TANBA = TANB
        SINBA = TANBA/DSQRT(TANBA**2+1.D0)
        COSBA = SINBA/TANBA        
        SINB = TANB/DSQRT(TANB**2+1.D0)
        COSB = SINB/TANB

      RMTOP = RUNM(MTOP,6)
      MB = RUNM(MTOP,5)
      PI = 4*DATAN(1D0)
      MZ = AMZ
      MW = AMW
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = ALPHAS(AMZ,2)
      ALPHA3  = ALPHAS(MTOP,2)

      G1 = DSQRT(ALPHA1*4.*PI)
      G2 = DSQRT(ALPHA2*4.*PI)
      G3 = DSQRT(ALPHA3*4.*PI)
      
        IF(MQ.GT.MUR) MST = MQ
        IF(MUR.GT.MQ.OR.MUR.EQ.MQ) MST = MUR
        MSUSYT = DSQRT(MST**2  + MTOP**2)

	IF(MQ.GT.MD) MSB = MQ
	IF(MD.GT.MQ.OR.MD.EQ.MQ) MSB = MD
	MSUSYB = DSQRT(MSB**2 + MB**2)

	TT = LOG(MSUSYT**2/MTOP**2)
	TB = LOG(MSUSYB**2/MTOP**2)

        HT = RMTOP/V/SINB
        HTST = RMTOP/V
        HB =  MB/V/COSB
        G32 = ALPHA3*4.*PI

        BT2 = -(8.*G32 - 9.*HT**2/2. - HB**2/2.)/(4.*PI)**2
	BB2 = -(8.*G32 - 9.*HB**2/2. - HT**2/2.)/(4.*PI)**2
        AL2 = 3./8./PI**2*HT**2
        BT2ST = -(8.*G32 - 9.*HTST**2/2.)/(4.*PI)**2
        ALST = 3./8./PI**2*HTST**2
        AL1 = 3./8./PI**2*HB**2

        IF(MA.GT.MTOP) THEN
        VI = V*(1. + 3./32./PI**2*HTST**2*LOG(MTOP**2/MA**2))
        H1I = VI*COSBA
        H2I = VI*SINBA
        H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYT**2))**.25
        H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYT**2))**.25
        H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYB**2))**.25
        H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYB**2))**.25
        ELSE
        VI =  V
        H1I = VI*COSB
        H2I = VI*SINB
        H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYT**2))**.25
        H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYT**2))**.25
        H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYB**2))**.25
        H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYB**2))**.25
        END IF

        TANBST = H2T/H1T
        SINBT = TANBST/(1.+TANBST**2)**.5
        COSBT = SINBT/TANBST

        TANBSB = H2B/H1B
        SINBB = TANBSB/(1.+TANBSB**2)**.5
        COSBB = SINBB/TANBSB

        deltamt = 0
        deltamb = 0

        mtop4 = rmtop**4.*(1.+2.*bt2*tt- al2*tt - 4.*deltamt)
c     * /(1.+deltamt)**4.
        mbot4 = mb**4.*(1.+2.*bb2*tb - al1*tb)
     * /(1.+deltamb)**4.
        MTOP2 = DSQRT(MTOP4)
	MBOT2 = DSQRT(MBOT4)

        STOP12 = (MQ2 + MUR2)*.5 + MTOP2 
     *   +1./8.*(G2**2+G1**2)*(H1T**2-H2T**2)
     *   +(((G2**2-5.*G1**2/3.)/4.*(H1T**2-H2T**2) +
     *   MQ2 - MUR2)**2*0.25 + MTOP2*(AT-MU/TANBST)**2)**.5

        STOP22 = (MQ2 + MUR2)*.5 + MTOP2 
     *  +1./8.*(G2**2+G1**2)*(H1T**2-H2T**2) 
     *   - (((G2**2-5.*G1**2/3.)/4.*(H1T**2-H2T**2) +
     *  MQ2 - MUR2)**2*0.25 
     *  + MTOP2*(AT-MU/TANBST)**2)**.5

        IF(STOP22.LT.0.) GOTO 4237

        SBOT12 = (MQ2 + MD2)*.5  
     *   - 1./8.*(G2**2+G1**2)*(H1B**2-H2B**2)
     *  + (((G1**2/3.-G2**2)/4.*(H1B**2-H2B**2) +
     *  MQ2 - MD2)**2*0.25 + MBOT2*(AB-MU*TANBSB)**2)**.5

        SBOT22 = (MQ2 + MD2)*.5  
     *   - 1./8.*(G2**2+G1**2)*(H1B**2-H2B**2)
     *   - (((G1**2/3.-G2**2)/4.*(H1B**2-H2B**2) +
     *   MQ2 - MD2)**2*0.25 + MBOT2*(AB-MU*TANBSB)**2)**.5

        IF(SBOT22.LT.0.) GOTO 4237

        STOP1 = STOP12**.5
        STOP2 = STOP22**.5
        SBOT1 = SBOT12**.5
        SBOT2 = SBOT22**.5

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Here is the definition of deltamb and deltamt, which
c     are the vertex corrections to the bottom and top quark
c     mass, keeping the dominant QCD and top Yukawa coupling
c     induced corrections.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        deltamb = -2*alpha3/3./pi*mglu*(ab-mu*tanb)*
     *  T(sbot1,sbot2,mglu)
     *  + ht**2/(4.*pi)**2*(at-mu/tanb)*mu*tanb*
     *  T(stop1,stop2,mu)


        deltamt = -2.*alpha3/3./pi*(at-mu/tanb)*mglu*
     *  T(stop1,stop2,mglu)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Here the new values of the top and bottom quark masses at
c   the scale MS are defined, to be used in the effective
c   potential approximation. They are just the old ones, but
c   including the finite corrections deltamt and deltamb.
c   The deltamb corrections can become large and are resummed
c   to all orders, as suggested in the two recent works by M. Carena,
c   S. Mrenna and C.E.M. Wagner, as well as in the work by M. Carena,
c   D. Garcia, U. Nierste and C.E.M. Wagner, to appear. The top
c   quark mass corrections are small and are kept in the perturbative
c   formulation. The function T(X,Y,Z) is necessary for the calculation.
c   the entries are masses and NOT their squares !
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        mtop4 = rmtop**4.*(1.+2.*bt2*tt- al2*tt - 4.*deltamt)
c     * /(1.+deltamt)**4.
        mbot4 = mb**4.*(1.+2.*bb2*tb - al1*tb)
     * /(1.+deltamb)**4.
        MTOP2 = DSQRT(MTOP4)
	MBOT2 = DSQRT(MBOT4)

        STOP12 = (MQ2 + MUR2)*.5 + MTOP2 
     *   +1./8.*(G2**2+G1**2)*(H1T**2-H2T**2)
     *   +(((G2**2-5.*G1**2/3.)/4.*(H1T**2-H2T**2) +
     *   MQ2 - MUR2)**2*0.25 + MTOP2*(AT-MU/TANBST)**2)**.5

        STOP22 = (MQ2 + MUR2)*.5 + MTOP2 
     *  +1./8.*(G2**2+G1**2)*(H1T**2-H2T**2) 
     *   - (((G2**2-5.*G1**2/3.)/4.*(H1T**2-H2T**2) +
     *  MQ2 - MUR2)**2*0.25 
     *  + MTOP2*(AT-MU/TANBST)**2)**.5

        IF(STOP22.LT.0.) GOTO 4237

        SBOT12 = (MQ2 + MD2)*.5  
     *   - 1./8.*(G2**2+G1**2)*(H1B**2-H2B**2)
     *  + (((G1**2/3.-G2**2)/4.*(H1B**2-H2B**2) +
     *  MQ2 - MD2)**2*0.25 + MBOT2*(AB-MU*TANBSB)**2)**.5

        SBOT22 = (MQ2 + MD2)*.5  
     *   - 1./8.*(G2**2+G1**2)*(H1B**2-H2B**2)
     *   - (((G1**2/3.-G2**2)/4.*(H1B**2-H2B**2) +
     *   MQ2 - MD2)**2*0.25 + MBOT2*(AB-MU*TANBSB)**2)**.5

4237    RETURN
        END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       End of program from M. Carena, M. Quiros and C.E.M. Wagner.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE AMHAMA(ICASE,MH,TANB)
C--CALCULATION OF PSEUDOSCALAR HIGGS MASS FROM HIGGS MASS MH
C--ICASE=0: MH=PSEUDOSCALAR MASS
C--ICASE=1: MH=LIGHT SCALAR MASS
C--ICASE=2: MH=HEAVY SCALAR MASS
C--ICASE=3: MH=CHARGED HIGGS MASS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8(A-H,L,M,O-Z)
      DIMENSION VH(2,2),M2(2,2),M2P(2,2)
      COMMON/HMASS/AMSM,AMA,AML,AMH,AMCH,AMAR
      COMMON/MODEL/IMODEL
      MA0 = 0.D0
      MA1 = 1.D8
      IF(IMODEL.EQ.3)MA0 = 1.D0
      IF(IMODEL.EQ.3)MA1 = 1.D5
      IF(ICASE.EQ.0)THEN
       MA = MH
      ELSE
       IF(ICASE.EQ.1)THEN
        AMA = MA0
        CALL SUSYCP(TANB)
        AMH0 = AML
        AMA = MA1
        CALL SUSYCP(TANB)
        AMH1 = AML
        IF(MH.LT.AMH0)THEN
         WRITE(6,*)'HIGGS MASS OUT OF RANGE. TAKING LOWER LIMIT'
         MA = MA0
         AMA = MA
         CALL SUSYCP(TANB)
         RETURN
        ENDIF
        IF(MH.GT.AMH1)THEN
         WRITE(6,*)'HIGGS MASS OUT OF RANGE. TAKING UPPER LIMIT'
         MA = MA1
         AMA = MA
         CALL SUSYCP(TANB)
         RETURN
        ENDIF
       ENDIF
       IF(ICASE.EQ.2)THEN
        AMA = MA0
        CALL SUSYCP(TANB)
        AMH0 = AMH
        IF(MH.LT.AMH0)THEN
         WRITE(6,*)'HIGGS MASS OUT OF RANGE. TAKING LOWER LIMIT'
         MA = MA0
         AMA = MA
         CALL SUSYCP(TANB)
         RETURN
        ENDIF
       ENDIF
       DEL0 = 1.D-4
c      MA0 = 0.D0
c      MA1 = 1.D4
1      MA = (MA0+MA1)/2
C      CALL SUBH(MA,TANB,MQ,MUR,MD,MTOP,AU,AD,MU,MCHI0,
C    *                 MHP,HMP,MCH,SA,CA,TANBA)
       AMA = MA
       MX0 = MX
       CALL SUSYCP(TANB)
       IF(ICASE.EQ.1)THEN
        MX = AML
       ELSEIF(ICASE.EQ.2)THEN
        MX = AMH
       ELSEIF(ICASE.EQ.3)THEN
        MX = AMCH
       ENDIF
       DEL = DABS(MA1 - MA0)/MA
       IF(MX.EQ.MH) DEL = 0
       IF(DEL.GT.DEL0) THEN
        IF(MX.GT.MH) MA1 = MA
        IF(MX.LT.MH) MA0 = MA
        GOTO 1
       ENDIF
       IF(DEL.NE.0)THEN
       FAC = 1
       MAXI = DINT(FAC*MA+0.5D0)/FAC
C      CALL SUBH(MAXI,TANB,MQ,MUR,MD,MTOP,AU,AD,MU,MCHI0,
C    *                 MHP,HMP,MCH,SA,CA,TANBA)
       AMA = MAXI
       CALL SUSYCP(TANB)
       IF(ICASE.EQ.1)THEN
        MX = AML
       ELSEIF(ICASE.EQ.2)THEN
        MX = AMH
       ELSEIF(ICASE.EQ.3)THEN
        MX = AMCH
       ENDIF
       IF(MX.EQ.MH)THEN
        MA = MAXI
       ELSE
        DEL0 = 1.D-8
2       MA = (MA0+MA1)/2
C       CALL SUBH(MA,TANB,MQ,MUR,MD,MTOP,AU,AD,MU,MCHI0,
C    *                  MHP,HMP,MCH,SA,CA,TANBA)
        AMA = MA
        MX0 = MX
        CALL SUSYCP(TANB)
        IF(ICASE.EQ.1)THEN
         MX = AML
        ELSEIF(ICASE.EQ.2)THEN
         MX = AMH
        ELSEIF(ICASE.EQ.3)THEN
         MX = AMCH
        ENDIF
        DEL = DABS(MA1 - MA0)/MA
        IF(MX.EQ.MH) DEL = 0
        IF(DEL.GT.DEL0) THEN
         IF(MX.GT.MH) MA1 = MA
         IF(MX.LT.MH) MA0 = MA
         GOTO 2
        ENDIF
       ENDIF
      ENDIF
      ENDIF
      AMA = MA
      CALL SUSYCP(TANB)
      RETURN
      END

      DOUBLE PRECISION FUNCTION RUNM(Q,NF)
C--RUNNING QUARK MASSES: Q = SCALE, NF = 3,..,6 QUARK TYPE
      PARAMETER (NN=6)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZETA3 = 1.202056903159594D0)
      DIMENSION AM(NN),YMSB(NN)
      COMMON/ALS/XLAMBDA,AMCA,AMBA,AMTA,N0A
      COMMON/MASSES/AMS,AMC,AMB,AMT
      COMMON/STRANGE/AMSB
      SAVE ISTRANGE
      B0(NF)=(33.D0-2.D0*NF)/12D0
      B1(NF) = (102D0-38D0/3D0*NF)/16D0
      B2(NF) = (2857D0/2D0-5033D0/18D0*NF+325D0/54D0*NF**2)/64D0
      G0(NF) = 1D0
      G1(NF) = (202D0/3D0-20D0/9D0*NF)/16D0
      G2(NF) = (1249D0-(2216D0/27D0+160D0/3D0*ZETA3)*NF
     .       - 140D0/81D0*NF**2)/64D0
      C1(NF) = G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2
      C2(NF) = ((G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2)**2
     .       + G2(NF)/B0(NF) + B1(NF)**2*G0(NF)/B0(NF)**3
     .       - B1(NF)*G1(NF)/B0(NF)**2 - B2(NF)*G0(NF)/B0(NF)**2)/2D0
      TRAN(X,XK)=1D0+4D0/3D0*ALPHAS(X,2)/PI+XK*(ALPHAS(X,2)/PI)**2
      CQ(X,NF)=(2D0*B0(NF)*X)**(G0(NF)/B0(NF))
     .            *(1D0+C1(NF)*X+C2(NF)*X**2)
      DATA ISTRANGE/0/
      PI=4D0*DATAN(1D0)
      ACC = 1.D-8
      AM(1) = 0
      AM(2) = 0
C--------------------------------------------
      IMSBAR = 0
      NNLO = 0
      IF(IMSBAR.EQ.1)THEN
       IF(ISTRANGE.EQ.0)THEN
C--STRANGE POLE MASS FROM MSBAR-MASS AT 1 GEV
        AMSD = XLAMBDA
        AMSU = 1.D8
123     AMS  = (AMSU+AMSD)/2
        AM(3) = AMS
        XMSB = AMS/CQ(ALPHAS(AMS,2)/PI,3)
     .            *CQ(ALPHAS(1.D0,2)/PI,3)/TRAN(AMS,0D0)
        DD = (XMSB-AMSB)/AMSB
        IF(DABS(DD).GE.ACC)THEN
         IF(DD.LE.0.D0)THEN
          AMSD = AM(3)
         ELSE
          AMSU = AM(3)
         ENDIF
         GOTO 123
        ENDIF
        ISTRANGE=1
       ENDIF
       AM(3) = AMSB
      ELSE
       AMS=AMSB
       AM(3) = AMS
      ENDIF
C--------------------------------------------
      AM(3) = AMSB
      AM(4) = AMC
      AM(5) = AMB
      AM(6) = AMT
      XK = 16.11D0
      DO 1 I=1,NF-1
       XK = XK - 1.04D0*(1.D0-AM(I)/AM(NF))
1     CONTINUE
      IF(NF.GE.4)THEN
       XMSB = AM(NF)/TRAN(AM(NF),0D0)
       XMHAT = XMSB/CQ(ALPHAS(AM(NF),2)/PI,NF)
      ELSE
       XMSB = 0
       XMHAT = 0
      ENDIF
      YMSB(3) = AMSB
      IF(NF.EQ.3)THEN
       YMSB(4) = YMSB(3)*CQ(ALPHAS(AM(4),2)/PI,3)/
     .                   CQ(ALPHAS(1.D0,2)/PI,3)
       YMSB(5) = YMSB(4)*CQ(ALPHAS(AM(5),2)/PI,4)/
     .                   CQ(ALPHAS(AM(4),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .                   CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.4)THEN
       YMSB(4) = XMSB
       YMSB(5) = YMSB(4)*CQ(ALPHAS(AM(5),2)/PI,4)/
     .                   CQ(ALPHAS(AM(4),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .                   CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.5)THEN
       YMSB(5) = XMSB
       YMSB(4) = YMSB(5)*CQ(ALPHAS(AM(4),2)/PI,4)/
     .                   CQ(ALPHAS(AM(5),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .                   CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.6)THEN
       YMSB(6) = XMSB
       YMSB(5) = YMSB(6)*CQ(ALPHAS(AM(5),2)/PI,5)/
     .                   CQ(ALPHAS(AM(6),2)/PI,5)
       YMSB(4) = YMSB(5)*CQ(ALPHAS(AM(4),2)/PI,4)/
     .                   CQ(ALPHAS(AM(5),2)/PI,4)
      ENDIF
      IF(Q.LT.AMC)THEN
       N0=3
       Q0 = 1.D0
      ELSEIF(Q.LE.AMB)THEN
       N0=4
       Q0 = AMC
      ELSEIF(Q.LE.AMT)THEN
       N0=5
       Q0 = AMB
      ELSE
       N0=6
       Q0 = AMT
      ENDIF
      IF(NNLO.EQ.1.AND.NF.GT.3)THEN
       XKFAC = TRAN(AM(NF),0D0)/TRAN(AM(NF),XK)
      ELSE
       XKFAC = 1D0
      ENDIF
      RUNM = YMSB(N0)*CQ(ALPHAS(Q,2)/PI,N0)/
     .               CQ(ALPHAS(Q0,2)/PI,N0)
     .       * XKFAC
      RETURN
      END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     THIS PROGRAM COMPUTES THE RENORMALIZATION GROUP IMPROVED
C     VALUES OF HIGGS MASSES AND COUPLINGS IN THE MSSM.
C
C     INPUT: MA,TANB = TAN(BETA),MQ,MUR,MTOP,AU,AD,MU.
C
C     ALL MASSES IN GEV UNITS. MA IS THE CP-ODD HIGGS MASS,
C     MTOP IS THE PHYSICAL TOP MASS, MQ AND MUR ARE THE SOFT
C     SUPERSYMMETRY BREAKING MASS PARAMETERS OF LEFT HANDED
C     AND RIGHT HANDED STOPS RESPECTIVELY, AU AND AD ARE THE
C     STOP AND SBOTTOM TRILINEAR SOFT BREAKING TERMS,
C     RESPECTIVELY,  AND MU IS THE SUPERSYMMETRIC
C     HIGGS MASS PARAMETER. WE USE THE  CONVENTIONS FROM
C     THE PHYSICS REPORT OF HABER AND KANE: LEFT RIGHT
C     STOP MIXING TERM PROPORTIONAL TO (AU - MU/TANB).
C
C     WE USE AS INPUT TANB DEFINED AT THE SCALE MTOP.
C
C     OUTPUT: MH,HM,MHCH, SA = SIN(ALPHA), CA= COS(ALPHA), TANBA
C
C     WHERE MH AND HM ARE THE LIGHTEST AND HEAVIEST CP-EVEN
C     HIGGS MASSES, MHCH IS THE CHARGED HIGGS MASS AND
C     ALPHA IS THE HIGGS MIXING ANGLE.
C
C     TANBA IS THE ANGLE TANB AT THE CP-ODD HIGGS MASS SCALE.
C
C     RANGE OF VALIDITY:
C
C    (STOP1**2 - STOP2**2)/(STOP2**2 + STOP1**2) < 0.5
C    (SBOT1**2 - SBOT2**2)/(SBOT2**2 + SBOT2**2) < 0.5
C
C     WHERE STOP1, STOP2, SBOT1 AND SBOT2 ARE THE STOP AND
C     ARE THE SBOTTOM  MASS EIGENVALUES, RESPECTIVELY. THIS
C     RANGE AUTOMATICALLY EXCLUDES THE EXISTENCE OF TACHYONS.
C
C
C     FOR THE CHARGED HIGGS MASS COMPUTATION, THE METHOD IS
C     VALID IF
C
C     2 * |MB * AD* TANB|  < M_SUSY**2,  2 * |MTOP * AU| < M_SUSY**2
C
C     2 * |MB * MU * TANB| < M_SUSY**2,  2 * |MTOP * MU| < M_SUSY**2
C
C     WHERE M_SUSY**2 IS THE AVERAGE OF THE SQUARED STOP MASS
C     EIGENVALUES, M_SUSY**2 = (STOP1**2 + STOP2**2)/2. THE SBOTTOM
C     MASSES HAVE BEEN ASSUMED TO BE OF ORDER OF THE STOP ONES.
C
C     M_SUSY**2 = (MQ**2 + MUR**2)*0.5 + MTOP**2
C
C     PROGRAM BASED ON THE WORK BY M. CARENA, J.R. ESPINOSA,
C     M. QUIROS AND C.E.M. WAGNER, PHYS. LETT. B355 (1995) 209
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SUBH2(MA,TANB,MQ,MUR,MTOP,AU,AD,MU,MH,HM,
     * MHCH,SA,CA,TANBA)
      IMPLICIT REAL*8(A-H,L,M,O-Z)
      COMMON/PARAM/GF,AMW,AMZ
      COMMON/HSELF/LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,LAMBDA5,
     .             LAMBDA6,LAMBDA7

C     MZ = 91.18
C     ALPHA1 = 0.0101
C     ALPHA2 = 0.0337
C     ALPHA3Z = 0.12
C     V = 174.1
C     PI = 3.14159
      TANBA = TANB
      TANBT = TANB

C     MBOTTOM(MTOP) = 3. GEV
C     MB = 3.
C     ALPHA3 = ALPHA3Z/(1. +(11. - 10./3.)/4./PI*ALPHA3Z*
C    *LOG(MTOP**2/MZ**2))

C     RMTOP= RUNNING TOP QUARK MASS
C     RMTOP = MTOP/(1.+4.*ALPHA3/3./PI)
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      MB = RUNM(MTOP,5)
      PI = 4*DATAN(1D0)
      MZ = AMZ
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = ALPHAS(AMZ,2)
      ALPHA3  = ALPHAS(MTOP,2)
      RMTOP   = RUNM(MTOP,6)
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C      RMTOP=MTOP
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      MS = ((MQ**2 + MUR**2)/2. + MTOP**2)**.5
      T = LOG(MS**2/MTOP**2)
      SINB = TANB/((1. + TANB**2)**.5)
      COSB = SINB/TANB
C      IF(MA.LE.MTOP) TANBA = TANBT
      IF(MA.GT.MTOP)
     *TANBA = TANBT*(1.-3./32./PI**2*
     *(RMTOP**2/V**2/SINB**2-MB**2/V**2/COSB**2)*
     *LOG(MA**2/MTOP**2))

      SINBT = TANBT/((1. + TANBT**2)**.5)
      COSBT = 1./((1. + TANBT**2)**.5)
      COS2BT = (TANBT**2 - 1.)/(TANBT**2 + 1.)
      G1 = (ALPHA1*4.*PI)**.5
      G2 = (ALPHA2*4.*PI)**.5
      G3 = (ALPHA3*4.*PI)**.5
      HU = RMTOP/V/SINBT
      HD =  MB/V/COSBT

C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C      G3=0
C      HU=0
C      HD=0
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      XAU = (2.*AU**2/MS**2)*(1. - AU**2/12./MS**2)
      XAD = (2.*AD**2/MS**2)*(1. - AD**2/12./MS**2)
      AUD = (-6.*MU**2/MS**2 - ( MU**2- AD*AU)**2/MS**4.
     *+ 3.*(AU + AD)**2/MS**2)/6.
      LAMBDA1 = ((G1**2 + G2**2)/4.)*(1.-3.*HD**2*T/8./PI**2)
     *+(3.*HD**4/8./PI**2) * (T + XAD/2. + (3.*HD**2/2. + HU**2/2.
     *- 8.*G3**2) * (XAD*T + T**2)/16./PI**2)
     *-(3.*HU**4* MU**4/96./PI**2/MS**4) * (1+ (9.*HU**2 -5.* HD**2
     *-  16.*G3**2) *T/16./PI**2)
      LAMBDA2 = ((G1**2 + G2**2)/4.)*(1.-3.*HU**2*T/8./PI**2)
     *+(3.*HU**4/8./PI**2) * (T + XAU/2. + (3.*HU**2/2. + HD**2/2.
     *- 8.*G3**2) * (XAU*T + T**2)/16./PI**2)
     *-(3.*HD**4* MU**4/96./PI**2/MS**4) * (1+ (9.*HD**2 -5.* HU**2
     *-  16.*G3**2) *T/16./PI**2)
      LAMBDA3 = ((G2**2 - G1**2)/4.)*(1.-3.*
     *(HU**2 + HD**2)*T/16./PI**2)
     *+(6.*HU**2*HD**2/16./PI**2) * (T + AUD/2. + (HU**2 + HD**2
     *- 8.*G3**2) * (AUD*T + T**2)/16./PI**2)
     *+(3.*HU**4/96./PI**2) * (3.*MU**2/MS**2 - MU**2*AU**2/
     *MS**4)* (1.+ (6.*HU**2 -2.* HD**2/2.
     *-  16.*G3**2) *T/16./PI**2)
     *+(3.*HD**4/96./PI**2) * (3.*MU**2/MS**2 - MU**2*AD**2/
     *MS**4)*(1.+ (6.*HD**2 -2.* HU**2
     *-  16.*G3**2) *T/16./PI**2)
      LAMBDA4 = (- G2**2/2.)*(1.-3.*(HU**2 + HD**2)*T/16./PI**2)
     *-(6.*HU**2*HD**2/16./PI**2) * (T + AUD/2. + (HU**2 + HD**2
     *- 8.*G3**2) * (AUD*T + T**2)/16./PI**2)
     *+(3.*HU**4/96./PI**2) * (3.*MU**2/MS**2 - MU**2*AU**2/
     *MS**4)*
     *(1+ (6.*HU**2 -2.* HD**2
     *-  16.*G3**2) *T/16./PI**2)
     *+(3.*HD**4/96./PI**2) * (3.*MU**2/MS**2 - MU**2*AD**2/
     *MS**4)*
     *(1+ (6.*HD**2 -2.* HU**2/2.
     *-  16.*G3**2) *T/16./PI**2)
      LAMBDA5 = -(3.*HU**4* MU**2*AU**2/96./PI**2/MS**4) *
     * (1- (2.*HD**2 -6.* HU**2 + 16.*G3**2) *T/16./PI**2)
     *-(3.*HD**4* MU**2*AD**2/96./PI**2/MS**4) *
     * (1- (2.*HU**2 -6.* HD**2 + 16.*G3**2) *T/16./PI**2)
      LAMBDA6 = (3.*HU**4* MU**3*AU/96./PI**2/MS**4) *
     * (1- (7.*HD**2/2. -15.* HU**2/2. + 16.*G3**2) *T/16./PI**2)
     *+(3.*HD**4* MU *(AD**3/MS**3 - 6.*AD/MS )/96./PI**2/MS) *
     * (1- (HU**2/2. -9.* HD**2/2. + 16.*G3**2) *T/16./PI**2)
      LAMBDA7 = (3.*HD**4* MU**3*AD/96./PI**2/MS**4) *
     * (1- (7.*HU**2/2. -15.* HD**2/2. + 16.*G3**2) *T/16./PI**2)
     *+(3.*HU**4* MU *(AU**3/MS**3 - 6.*AU/MS )/96./PI**2/MS) *
     * (1- (HD**2/2. -9.* HU**2/2. + 16.*G3**2) *T/16./PI**2)
      TRM2 = MA**2 + 2.*V**2* (LAMBDA1* COSBT**2 +
     *2.* LAMBDA6*SINBT*COSBT
     *+ LAMBDA5*SINBT**2 + LAMBDA2* SINBT**2 + 2.* LAMBDA7*SINBT*COSBT
     *+ LAMBDA5*COSBT**2)
      DETM2 = 4.*V**4*(-(SINBT*COSBT*(LAMBDA3 + LAMBDA4) +
     *LAMBDA6*COSBT**2
     *+ LAMBDA7* SINBT**2)**2 + (LAMBDA1* COSBT**2 +
     *2.* LAMBDA6* COSBT*SINBT
     *+ LAMBDA5*SINBT**2)*(LAMBDA2* SINBT**2 +2.* LAMBDA7* COSBT*SINBT
     *+ LAMBDA5*COSBT**2)) + MA**2*2.*V**2 *
     *((LAMBDA1* COSBT**2 +2.*
     *LAMBDA6* COSBT*SINBT + LAMBDA5*SINBT**2)*COSBT**2 +
     *(LAMBDA2* SINBT**2 +2.* LAMBDA7* COSBT*SINBT + LAMBDA5*COSBT**2)
     **SINBT**2
     * +2.*SINBT*COSBT* (SINBT*COSBT*(LAMBDA3
     * + LAMBDA4) + LAMBDA6*COSBT**2
     *+ LAMBDA7* SINBT**2))

      MH2 = (TRM2 - (TRM2**2 - 4.* DETM2)**.5)/2.
      HM2 = (TRM2 + (TRM2**2 - 4.* DETM2)**.5)/2.
      HM = HM2**.5
      MH = MH2**.5
      MHCH2 = MA**2 + (LAMBDA5 - LAMBDA4)* V**2
      MHCH = MHCH2**.5
      MHCH = MHCH2**.5

      SINALPHA = SQRT(((TRM2**2 - 4.* DETM2)**.5) -
     * ((2.*V**2*(LAMBDA1* COSBT**2 + 2.*
     *LAMBDA6* COSBT*SINBT
     *+ LAMBDA5*SINBT**2) + MA**2*SINBT**2)
     *- (2.*V**2*(LAMBDA2* SINBT**2 +2.* LAMBDA7* COSBT*SINBT
     *+ LAMBDA5*COSBT**2) + MA**2*COSBT**2)))/
     *SQRT(((TRM2**2 - 4.* DETM2)**.5))/2.**.5

      COSALPHA = (2.*(2.*V**2*(SINBT*COSBT*(LAMBDA3 + LAMBDA4) +
     *LAMBDA6*COSBT**2 + LAMBDA7* SINBT**2) -
     *MA**2*SINBT*COSBT))/2.**.5/
     *SQRT(((TRM2**2 - 4.* DETM2)**.5)*
     *(((TRM2**2 - 4.* DETM2)**.5) -
     * ((2.*V**2*(LAMBDA1* COSBT**2 + 2.*
     *LAMBDA6* COSBT*SINBT
     *+ LAMBDA5*SINBT**2) + MA**2*SINBT**2)
     *- (2.*V**2*(LAMBDA2* SINBT**2 +2.* LAMBDA7* COSBT*SINBT
     *+ LAMBDA5*COSBT**2) + MA**2*COSBT**2))))

      SA = -SINALPHA
      CA = -COSALPHA

 2242 RETURN
      END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C     ****************************************************************
C	  CHARGINO AND NEUTRALINO MASS MATRICES AND COUPLINGS
C     ****************************************************************
      SUBROUTINE GAUGINO(MU,M2,B,A,MC,MN,XMN,AC1,AC2,AC3,AN1,AN2
     .                 ,AN3,ACNL,ACNR,AGDL,AGDA,AGDH,AGDC)            
      IMPLICIT REAL*8(A-H,K-Z)
      COMPLEX*16 CXA,CXB,CXC,CXD,CX1,CX2,CX3
      DIMENSION MC(2),MN(4),XMN(4),Z(4,4),ZX(4,4),U(2,2),V(2,2),
     .          QQ(4,4),SS(4,4),S(2,2),Q(2,2),AC1(2,2),AC2(2,2),
     .          AC3(2,2),AN1(4,4),AN2(4,4),AN3(4,4),ACNL(2,4),
     .          ACNR(2,4),IORD(4),IREM(2)
      DIMENSION X(2,2)
      DIMENSION YMN(4),YZ(4,4),XMC(2),BU(2),BV(2)
      DIMENSION AGDL(4),AGDA(4),AGDH(4),AGDC(2)
c     COMMON/PARAM/GF,ALPH,AMTAU,AMMUON,MZ,MW
      COMMON/PARAM/GF,MW,MZ
      CW=MW/MZ
      SW=DSQRT(1-CW**2)
      PI=4.D0*DATAN(1.D0)
      SB=DSIN(B)
      CB=DCOS(B)
      TW=SW/CW
      M1=5.D0/3.D0*TW**2*M2
C     ************  NEUTRALINO MASSES AND MATRIX ELEMENTS ***********
      EPS=-1.D-10
      XC2=(M1*M2-MZ**2-MU**2)-3.D0/8.D0*(M1+M2)**2
      XC3=-1.D0/8.D0*(M1+M2)**3+1.D0/2.D0*(M1+M2)*(M1*M2-MZ**2
     .    -MU**2)+(M1+M2)*MU**2+(M1*CW**2+M2*SW**2)*MZ**2
     .    -MU*MZ**2*DSIN(2.D0*B)
      XC4=+(M1*CW**2+M2*SW**2)*MU*MZ**2*DSIN(2.D0*B)-M1*M2*MU**2
     .    +1.D0/4.D0*(M1+M2)*( (M1+M2)*MU**2+(M1*CW**2+M2*SW**2)
     .    *MZ**2-MU*MZ**2*DSIN(2.D0*B) )+1.D0/16.D0*(M1+M2)**2*
     .    (M1*M2-MZ**2-MU**2)-3.D0/256.D0*(M1+M2)**4
      XS=-XC3**2-2.D0/27.D0*XC2**3+8.D0/3.D0*XC2*XC4
      XU=-1.D0/3.D0*XC2**2-4.D0*XC4
      CXD=(-4*XU**3-27*XS**2)*DCMPLX(1.D0,EPS)
      CXC=1.D0/2.D0*(-XS+DCMPLX(0.D0,1.D0)*CDSQRT(CXD/27.D0))
      CXA=DREAL(CXC**(1.D0/3.D0))*DCMPLX(1.D0,-EPS)
      CXB=8.D0*CXA-8.D0/3.D0*XC2*DCMPLX(1.D0,-EPS)
C     *********** MASSES AND COUPLINGS:
      X0=(M1+M2)/4.D0
      CX1= CXA/2.D0-XC2/6.D0*DCMPLX(1.D0,-EPS)
      CX2=-CXA/2.D0-XC2/3.D0*DCMPLX(1.D0,-EPS)
      CX3=XC3*DCMPLX(1.D0,-EPS)/CDSQRT(CXB)
      XMN(1)=X0-CDABS(CDSQRT(CX1))+CDABS(CDSQRT(CX2+CX3))
      XMN(2)=X0+CDABS(CDSQRT(CX1))-CDABS(CDSQRT(CX2-CX3))
      XMN(3)=X0-CDABS(CDSQRT(CX1))-CDABS(CDSQRT(CX2+CX3))
      XMN(4)=X0+CDABS(CDSQRT(CX1))+CDABS(CDSQRT(CX2-CX3))
      DO 10 I=1,4
      MN(I)=DABS(XMN(I))
      YMN(I)=XMN(I)
      ZX(I,2)=-CW/SW*(M1-XMN(I))/(M2-XMN(I))
      ZX(I,3)=(MU*(M2-XMN(I))*(M1-XMN(I))-MZ**2*SB*CB*((M1-M2)*CW**2
     .       +M2-XMN(I)))/MZ/(M2-XMN(I))/SW/(MU*CB+XMN(I)*SB)
      ZX(I,4)=(-XMN(I)*(M2-XMN(I))*(M1-XMN(I))-MZ**2*CB*CB*((M1-M2)
     .       *CW**2+M2-XMN(I)))/MZ/(M2-XMN(I))/SW/(MU*CB+XMN(I)*SB)
      ZX(I,1)=1.D0/DSQRT(1.D0+ZX(I,2)**2+ZX(I,3)**2+ZX(I,4)**2) 
      YZ(I,1)=ZX(I,1)
      YZ(I,2)=ZX(I,2)*ZX(I,1)
      YZ(I,3)=ZX(I,3)*ZX(I,1)
      YZ(I,4)=ZX(I,4)*ZX(I,1)
 10   CONTINUE
C     *************  ORDERING THE DISORDER ******************
      XX0 = DMIN1(MN(1),MN(2),MN(3),MN(4))
      XX1 = DMAX1(MN(1),MN(2),MN(3),MN(4))
      IDUMMY = 1
      DO I = 1,4
       IF(MN(I).EQ.XX0)THEN
        IORD(1) = I
       ELSEIF(MN(I).EQ.XX1)THEN
        IORD(4) = I
       ELSE
        IREM(IDUMMY) = I
        IDUMMY = IDUMMY+1
       ENDIF
      ENDDO
      IF(MN(IREM(1)).LE.MN(IREM(2)))THEN
       IORD(2) = IREM(1)
       IORD(3) = IREM(2)
      ELSE
       IORD(2) = IREM(2)
       IORD(3) = IREM(1)
      ENDIF
C 
      DO 98 J=1,4
      I=IORD(J)
      XMN(J)=YMN(I)
      MN(J) =DABS(YMN(I))
        DO I1=1,4
        Z(J,I1)=YZ(I,I1)
        ENDDO
 98   CONTINUE
C     ************  NEUTRALINO COUPLINGS TO HIGGS BOSONS ***********
	DO 11 I=1,4
	DO 11 J=1,4
	QQ(I,J)=1.D0/2.D0*(Z(I,3)*(Z(J,2)-TW*Z(J,1))+Z(J,3)*
     .		(Z(I,2)-TW*Z(I,1)))
	SS(I,J)=1.D0/2.D0*(Z(I,4)*(Z(J,2)-TW*Z(J,1))+Z(J,4)*
     .		(Z(I,2)-TW*Z(I,1)))
 11	CONTINUE
	DO 21 I=1,4
	DO 21 J=1,4
	AN1(I,J)= QQ(I,J)*DCOS(A)-SS(I,J)*DSIN(A)
	AN2(I,J)=-QQ(I,J)*DSIN(A)-SS(I,J)*DCOS(A)
	AN3(I,J)= QQ(I,J)*DSIN(B)-SS(I,J)*DCOS(B)
 21	CONTINUE

C       ************* CHARGINO MASSES AND MATRIX ELEMENTS ***********
	DELTA=DABS(B-.25*PI)
	DDD=MU*DCOS(B)+M2*DSIN(B)
	CCC=MU*DSIN(B)+M2*DCOS(B)
	IF(DELTA.LT.0.01D0) THEN
	PHIM=PI/4.D0-.5D0*DATAN((M2-MU)/(2.D0*MW))
	PHIP=PHIM
	ELSE IF	(DABS(CCC).LT.1.D-5) THEN
	PHIM=0.D0
	PHIP=DATAN(DSQRT(2.D0)*MW*DSIN(B)/(M2+1.D-5))
	ELSE IF	(DABS(DDD).LT.1.D-5) THEN
	PHIP=0.D0
	PHIM=DATAN(DSQRT(2.D0)*MW*DCOS(B)/(M2+1.D-5))
	ELSE
	RAD=DSQRT((M2**2-MU**2)**2+4.D0*MW**4*DCOS(2.D0*B)**2
     +	+4.D0*MW**2*(M2**2+MU**2+2.D0*M2*MU*DSIN(2.D0*B)))
	PHIP=DATAN((RAD-(M2**2-MU**2+2.D0*MW**2*DCOS(2.D0*B)))
     +	/(2.D0*DSQRT(2.D0)*MW*(MU*DCOS(B)+M2*DSIN(B))))
	PHIM=DATAN((RAD-(M2**2-MU**2-2.D0*MW**2*DCOS(2.D0*B)))
     +	/(2.D0*DSQRT(2.D0)*MW*(MU*DSIN(B)+M2*DCOS(B))))
	ENDIF
	CP=DCOS(PHIP)
	SP=DSIN(PHIP)
	CM=DCOS(PHIM)
	SM=DSIN(PHIM)
C  MY CONVENTION
	U(2,2)=CM
	U(2,1)=-SM
	U(1,2)=SM
	U(1,1)=CM
	V(1,1)=CP
	V(1,2)=SP
	V(2,1)=-SP
	V(2,2)=CP
        X(1,1)=M2
        X(1,2)=DSQRT(2.D0)*MW*DSIN(B)
        X(2,1)=DSQRT(2.D0)*MW*DCOS(B)
        X(2,2)=MU
 555    CONTINUE
       XMC(1)=(U(1,1)*X(1,1)+U(1,2)*X(2,1))*V(1,1)
     .       +(U(1,1)*X(1,2)+U(1,2)*X(2,2))*V(1,2)
       XMC(2)=(U(2,1)*X(1,1)+U(2,2)*X(2,1))*V(2,1)
     .       +(U(2,1)*X(1,2)+U(2,2)*X(2,2))*V(2,2)
        IF(XMC(1).LT.0.D0) THEN
	V(1,1)=-CP
	V(1,2)=-SP
	V(2,1)=-SP
	V(2,2)=CP
        GOTO 555
        ENDIF
        IF(XMC(2).LT.0.D0) THEN
	V(1,1)=CP
	V(1,2)=SP
	V(2,1)=SP
	V(2,2)=-CP
        GOTO 555
        ENDIF
        IF(XMC(1).GT.XMC(2)) THEN
        MTEMP=XMC(1)
        XMC(1)=XMC(2)
        XMC(2)=MTEMP
        DO J=1,2
        BU(J)=U(1,J)
        U(1,J)=U(2,J)
        U(2,J)=BU(J)
        BV(J)=V(1,J)
        V(1,J)=V(2,J)
        V(2,J)=BV(J)
        ENDDO
        ENDIF        
        MC(1)=DABS(XMC(1))
        MC(2)=DABS(XMC(2))

C     ************  CHARGINO COUPLINGS TO HIGGS BOSONS ***********
	DO 12 I=1,2
	DO 12 J=1,2
	Q(I,J)=DSQRT(1.D0/2.D0)*U(J,2)*V(I,1)
	S(I,J)=DSQRT(1.D0/2.D0)*U(J,1)*V(I,2)
 12	CONTINUE
	DO 22 I=1,2
	DO 22 J=1,2	
	AC1(I,J)= Q(I,J)*DCOS(A)+S(I,J)*DSIN(A)
	AC2(I,J)=-Q(I,J)*DSIN(A)+S(I,J)*DCOS(A)
	AC3(I,J)= Q(I,J)*DSIN(B)+S(I,J)*DCOS(B)
 22	CONTINUE
C     **** CHARGINO-NEUTRALINO COUPLINGS TO CHARGED HIGGS BOSONS 
	DO 13 I=1,2
	DO 13 J=1,4
        ACNL(I,J)=DCOS(B)*(Z(J,4)*V(I,1)+(Z(J,2)+Z(J,1)*TW)
     .       *V(I,2)/DSQRT(2.D0)) 
        ACNR(I,J)=DSIN(B)*(Z(J,3)*U(I,1)-(Z(J,2)+Z(J,1)*TW)
     .       *U(I,2)/DSQRT(2.D0)) 
 13     CONTINUE

C   ************* HIGGS--NEUTRALINO--GOLDSTINO COUPLINGS
      DO 51 I=1,4
      AGDL(I)=Z(I,3)*DSIN(A)-Z(I,4)*DCOS(A)
      AGDH(I)=Z(I,3)*DCOS(A)+Z(I,4)*DSIN(A)
      AGDA(I)=Z(I,3)*DSIN(B)+Z(I,4)*DCOS(B)
 51   CONTINUE
C
C   ************* CHARGED HIGGS--CHARGINO--GOLDSTINO COUPLINGS
      AGDC(1)=DSQRT( V(1,2)**2*DCOS(B)**2+ U(1,2)**2*DSIN(B)**2 )
      AGDC(2)=DSQRT( V(2,2)**2*DCOS(B)**2+ U(2,2)**2*DSIN(B)**2 )

       RETURN
       END

C   ****************************************************************
C     SUBROUTINE FOR SFERMION MASSES, MIXING AND COUPLINGS 
C   ****************************************************************

       SUBROUTINE SFERMION(TSC,BSC,MQL,MUR,MDR,MEL,MER,AL,AT,AB,MU,
     .                    MST,MSB,MSL,MSU,MSD,MSE,MSN, 
     .                    GLEE,GLTT,GLBB,GHEE,GHTT,GHBB,
     .                    GAEE,GATT,GABB,GCEN,GCTB)

      IMPLICIT REAL*8(A-H,K-Z)
      DIMENSION MST(2),MSB(2),MSL(2),MSU(2),MSD(2),MSE(2),MSN(2),
     .          GCEN(2,2),GCTB(2,2),GLEE(2,2),GLTT(2,2),GLBB(2,2),
     .          GHEE(2,2),GHTT(2,2),GHBB(2,2)
      COMMON/MASSES/AMS,AMC,AMB,AMT
c     COMMON/PARAM/GF,ALPH,AMTAU,AMMUON,MZ,MW
      COMMON/PARAM/GF,MW,MZ
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,GHHH,GLLL,GHLL,
     .            GLHH,GHAA,GLAA,GLVV,GHVV,GLPM,GHPM,B,A
      COMMON/SFER1ST/MQL1,MUR1,MDR1,MEL1,MER1
      COMMON/BREAKGLU/AMGLU
      COMMON/GLUINO/AMGLUINO,XMSB1,XMSB2,STHB,CTHB,
     .              XLBB(2,2),XHBB(2,2),XABB(2,2),
     .              XMST1,XMST2,STHT,CTHT,
     .              XLTT(2,2),XHTT(2,2),XATT(2,2)
C
      PI = 4*DATAN(1.D0)
      SW2=1.D0-MW**2/MZ**2
      TB=DTAN(B)
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     MB = RUNM(AMT,5)
c     MT = RUNM(AMT,6)
c     ML = AMTAU
c     MT = RUNM(TSC,6)
c     MB = RUNM(BSC,5)
      MT = AMT
      MB = AMB
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ML = 1.7771D0
C FIRST TWO GENERATIONS:  NO MIXING INCLUDED 
C UP SQUARKS: 
      MSTL2=MQL1**2+(0.5D0-2.D0/3.D0*SW2)*MZ**2*DCOS(2.D0*B)
      MSTR2=MUR1**2+2.D0/3.D0*SW2*MZ**2*DCOS(2.D0*B) 
      MSU(1)=DSQRT(MSTL2)
      MSU(2)=DSQRT(MSTR2)
C DOWN SQUARKS
      MSBL2=MQL1**2+(-0.5D0+1.D0/3.D0*SW2)*MZ**2*DCOS(2.D0*B)
      MSBR2=MDR1**2-1.D0/3.D0*SW2*MZ**2*DCOS(2.D0*B) 
      MSD(1)=DSQRT(MSBL2)
      MSD(2)=DSQRT(MSBR2)
C SLEPTONS
      MSEL2=MEL1**2+(-0.5D0+SW2)*MZ**2*DCOS(2.D0*B)
      MSER2=MER1**2- SW2*MZ**2*DCOS(2.D0*B) 
      MSNL2=MEL1**2+0.5D0*MZ**2*DCOS(2.D0*B)
      MSE(1)=DSQRT(MSEL2)
      MSE(2)=DSQRT(MSER2)
      MSN(1)=DSQRT(MSNL2)
      MSN(2)=1.D+15

C NOW THE THIRD GENERATION
C
C STOP MASSES/MIXING
C
      MSTL2=MQL**2+(0.5D0-2.D0/3.D0*SW2)*MZ**2*DCOS(2.D0*B)
      MSTR2=MUR**2+2.D0/3.D0*SW2*MZ**2*DCOS(2.D0*B) 
      MLRT=AT-MU/TB
      DELT=(MSTL2-MSTR2)**2+4*MT**2*MLRT**2
      MST12=MT**2+0.5D0*(MSTL2+MSTR2-DSQRT(DELT))
      MST22=MT**2+0.5D0*(MSTL2+MSTR2+DSQRT(DELT))
        IF(MST12.LT.0.D0)THEN 
      PRINT *, 'MSTOP**2 is negative!!!!'
      GOTO 111 
      ELSE 
      MST(1)=DSQRT(MST12)
      MST(2)=DSQRT(MST22)
      IF(MSTL2.EQ.MSTR2) THEN
       THET = PI/4
      ELSE
       THET=0.5D0*DATAN(2.D0*MT*MLRT / (MSTL2-MSTR2) )
       IF(MSTL2.GT.MSTR2) THET = THET + PI/2
      ENDIF
        ENDIF 
      CT= DCOS(THET)
      ST= DSIN(THET) 
C
C SBOTTOM MASSES/MIXING
C
      MSBL2=MQL**2+(-0.5D0+1.D0/3.D0*SW2)*MZ**2*DCOS(2.D0*B)
      MSBR2=MDR**2-1.D0/3.D0*SW2*MZ**2*DCOS(2.D0*B) 
      MLRB=AB-MU*TB
      DELB=(MSBL2-MSBR2)**2+4*MB**2*MLRB**2
      MSB12=MB**2+0.5D0*(MSBL2+MSBR2-DSQRT(DELB))
      MSB22=MB**2+0.5D0*(MSBL2+MSBR2+DSQRT(DELB))
        IF(MSB12.LT.0.D0)THEN
      PRINT *, 'MSBOT**2 is negative!!!!'
      GOTO 111
        ELSE
      MSB(1)=DSQRT(MSB12)
      MSB(2)=DSQRT(MSB22)
      IF(MSBL2.EQ.MSBR2) THEN
       THEB = PI/4
      ELSE
       THEB=0.5D0*DATAN(2.D0*MB*MLRB / (MSBL2-MSBR2) )
       IF(MSBL2.GT.MSBR2) THEB = THEB + PI/2
      ENDIF
        ENDIF  
      CB= DCOS(THEB)
      SB= DSIN(THEB) 
c
c     write(6,*)'s_t, c_t: ',st,ct
c     write(6,*)'s_b, c_b: ',sb,cb

C
C  STAU MASSES/MIXING
C
      MSEL2=MEL**2+(-0.5D0+SW2)*MZ**2*DCOS(2.D0*B)
      MSER2=MER**2- SW2*MZ**2*DCOS(2.D0*B) 
      MSNL2=MEL**2+0.5D0*MZ**2*DCOS(2.D0*B)
      MLRE=AL-MU*TB
      DELE=(MSEL2-MSER2)**2+4*ML**2*MLRE**2
      MSE12=ML**2+0.5D0*(MSEL2+MSER2-DSQRT(DELE))
      MSE22=ML**2+0.5D0*(MSEL2+MSER2+DSQRT(DELE))
        IF(MSE12.LT.0.D0)THEN
      PRINT *, 'MSTAU**2 is negative!!!!'
      GOTO 111
        ELSE
      MSL(1)=DSQRT(MSE12)
      MSL(2)=DSQRT(MSE22)
      IF(MSEL2.EQ.MSER2) THEN
       THEL = PI/4
      ELSE
       THEL=0.5D0*DATAN(2.D0*ML*MLRE / (MSEL2-MSER2) )
       IF(MSEL2.GT.MSER2) THEL = THEL + PI/2
      ENDIF
        ENDIF  
      CL= DCOS(THEL)
      SL= DSIN(THEL) 
C
C LIGHT CP--EVEN HIGGS COUPLINGS TO STOPS
C 
      GLTT(1,1)=-DSIN(B+A)*(0.5D0*CT**2-2.D0/3.D0*SW2*DCOS(2*THET)) 
     .    + MT**2/MZ**2*GLT + MT*ST*CT/MZ**2*(AT*GLT+MU*GHT)
      GLTT(2,2)=-DSIN(B+A)*(0.5D0*ST**2+2.D0/3.D0*SW2*DCOS(2*THET))
     .    + MT**2/MZ**2*GLT - MT*ST*CT/MZ**2*(AT*GLT+MU*GHT)
      GLTT(1,2)=-2*DSIN(B+A)*ST*CT*(2.D0/3.D0*SW2-0.25D0)
     .    + MT*DCOS(2*THET)/2.D0/MZ**2*(AT*GLT+MU*GHT) 
      GLTT(2,1)=-2*DSIN(B+A)*ST*CT*(2.D0/3.D0*SW2-0.25D0)
     .    + MT*DCOS(2*THET)/2.D0/MZ**2*(AT*GLT+MU*GHT) 
C
C LIGHT CP--EVEN HIGGS COUPLINGS TO SBOTTOMS
C
      GLBB(1,1)=-DSIN(B+A)*(-0.5D0*CB**2+1.D0/3.D0*SW2*DCOS(2*THEB)) 
     .    + MB**2/MZ**2*GLB + MB*SB*CB/MZ**2*(AB*GLB-MU*GHB)
      GLBB(2,2)=-DSIN(B+A)*(-0.5D0*SB**2-1.D0/3.D0*SW2*DCOS(2*THEB)) 
     .    + MB**2/MZ**2*GLB - MB*SB*CB/MZ**2*(AB*GLB-MU*GHB)
      GLBB(1,2)=-2*DSIN(B+A)*SB*CB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/MZ**2*(AB*GLB-MU*GHB) 
      GLBB(2,1)=-2*DSIN(B+A)*SB*CB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/MZ**2*(AB*GLB-MU*GHB) 

C
C LIGHT CP--EVEN HIGGS COUPLINGS TO STAU'S 
C
      GLEE(1,1)=-DSIN(B+A)*(-0.5D0*CL**2+SW2*DCOS(2*THEL)) 
     .    + ML**2/MZ**2*GLB + ML*SL*CL/MZ**2*(AL*GLB-MU*GHB)
      GLEE(2,2)=-DSIN(B+A)*(-0.5D0*SL**2-SW2*DCOS(2*THEL)) 
     .    + ML**2/MZ**2*GLB - ML*SL*CL/MZ**2*(AL*GLB-MU*GHB)
      GLEE(1,2)=-2*DSIN(B+A)*SL*CL*(-SW2+0.25D0)
     .    + ML*DCOS(2*THEL)/2.D0/MZ**2*(AL*GLB-MU*GHB) 
      GLEE(2,1)=-2*DSIN(B+A)*SL*CL*(-SW2+0.25D0)
     .    + ML*DCOS(2*THEL)/2.D0/MZ**2*(AL*GLB-MU*GHB) 
C
C HEAVY CP--EVEN HIGGS COUPLINGS TO STOPS
C
      GHTT(1,1)=DCOS(B+A)*(0.5D0*CT**2-2.D0/3.D0*SW2*DCOS(2*THET)) 
     .    + MT**2/MZ**2*GHT + MT*ST*CT/MZ**2*(AT*GHT-MU*GLT)
      GHTT(2,2)= DCOS(B+A)*(0.5D0*ST**2+2.D0/3.D0*SW2*DCOS(2*THET))
     .    + MT**2/MZ**2*GHT - MT*ST*CT/MZ**2*(AT*GHT-MU*GLT)
      GHTT(1,2)=2*DCOS(B+A)*ST*CT*(2.D0/3.D0*SW2-0.25D0)
     .    + MT*DCOS(2*THET)/2.D0/MZ**2*(AT*GHT-MU*GLT) 
      GHTT(2,1)=2*DCOS(B+A)*ST*CT*(2.D0/3.D0*SW2-0.25D0)
     .    + MT*DCOS(2*THET)/2.D0/MZ**2*(AT*GHT-MU*GLT) 
C
C HEAVY CP--EVEN HIGGS COUPLINGS TO SBOTTOMS
C
      GHBB(1,1)= DCOS(B+A)*(-0.5D0*CB**2+1.D0/3.D0*SW2*DCOS(2*THEB)) 
     .    + MB**2/MZ**2*GHB + MB*SB*CB/MZ**2*(AB*GHB+MU*GLB)
      GHBB(2,2)= DCOS(B+A)*(-0.5D0*SB**2-1.D0/3.D0*SW2*DCOS(2*THEB)) 
     .    + MB**2/MZ**2*GHB - MB*SB*CB/MZ**2*(AB*GHB+MU*GLB)
      GHBB(1,2)=2*DCOS(B+A)*SB*CB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/MZ**2*(AB*GHB+MU*GLB) 
      GHBB(2,1)=2*DCOS(B+A)*SB*CB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/MZ**2*(AB*GHB+MU*GLB) 
C
C HEAVY CP--EVEN HIGGS COUPLINGS TO STAU'S 
C
      GHEE(1,1)= DCOS(B+A)*(-0.5D0*CL**2+SW2*DCOS(2*THEL)) 
     .    + ML**2/MZ**2*GHB + ML*SL*CL/MZ**2*(AL*GHB+MU*GLB)
      GHEE(2,2)= DCOS(B+A)*(-0.5D0*SL**2-SW2*DCOS(2*THEL)) 
     .    + ML**2/MZ**2*GHB - ML*SL*CL/MZ**2*(AL*GHB+MU*GLB)
      GHEE(1,2)=2*DCOS(B+A)*SL*CL*(-SW2+0.25D0)
     .    + ML*DCOS(2*THEL)/2.D0/MZ**2*(AL*GHB+MU*GLB) 
      GHEE(2,1)=2*DCOS(B+A)*SL*CL*(-SW2+0.25D0)
     .    + ML*DCOS(2*THEL)/2.D0/MZ**2*(AL*GHB+MU*GLB) 

C
C PSEUDOSCALAR COUPLINGS 
C
      GATT=MT/2.D0/MZ**2*(MU+AT*GAT) 
      GABB=MB/2.D0/MZ**2*(MU+AB*GAB) 
      GAEE=ML/2.D0/MZ**2*(MU+AL*GAB) 
C
C CHARGED HIGGS COUPLINGS STOPS/SBOTTOMS 
C
      CLL=(MW**2*DSIN(2*B)-MT**2*GAT-MB**2*GAB)/DSQRT(2.D0)/MW**2
      CRR=-MT*MB*(GAT+GAB)/DSQRT(2.D0)/MW**2
      CLR=-MB*(MU+AB*GAB)/DSQRT(2.D0)/MW**2
      CRL=-MT*(MU+AT*GAT)/DSQRT(2.D0)/MW**2
      GCTB(1,1)=+CT*CB*CLL+ST*SB*CRR+CT*SB*CLR+ST*CB*CRL
      GCTB(1,2)=-CT*SB*CLL+ST*CB*CRR+CT*CB*CLR-ST*SB*CRL
      GCTB(2,1)=-ST*CB*CLL+CT*SB*CRR-ST*SB*CLR+CT*CB*CRL
      GCTB(2,2)=+ST*SB*CLL+CT*CB*CRR-ST*CB*CLR-CT*SB*CRL

C
C CHARGED HIGGS COUPLINGS TAU'S AND NEUTRINOS 
C
      CLL=(MW**2*DSIN(2*B)-ML**2*GAB)/DSQRT(2.D0)/MW**2
      CLR=-ML*(MU+AL*GAB)/DSQRT(2.D0)/MW**2
      GCEN(1,1)=CL*CLL+SL*CLR
      GCEN(1,2)=-SL*CLL+CL*CLR
      GCEN(2,1)=0.D0
      GCEN(2,2)=0.D0 

C--FILL COMMON BLOCK GLUINO FOR SUSY-QCD CORRECTIONS TO
C  HIGGS -> BB, SQUARKS
      XMST1 = MST(1)
      XMST2 = MST(2)
      XMSB1 = MSB(1)
      XMSB2 = MSB(2)
      STHT = ST
      CTHT = CT
      STHB = SB
      CTHB = CB
       DO I=1,2
        DO J=1,2
         XLBB(I,J) = GLBB(I,J)
         XHBB(I,J) = GHBB(I,J)
         XABB(I,J) = 0
         XLTT(I,J) = GLTT(I,J)
         XHTT(I,J) = GHTT(I,J)
         XATT(I,J) = 0
        ENDDO
       ENDDO
       XABB(1,2) = GABB
       XABB(2,1) = -GABB
       XATT(1,2) = GATT
       XATT(2,1) = -GATT
       AMGLUINO = AMGLU

c     write(6,*)mst(1),mst(2),gltt(1,1)*mz**2,gltt(2,2)*mz**2,
c    .         gltt(1,2)*mz**2,glt,ght,st,ct,dtan(b),dtan(a),mz,mt,at,mu
      RETURN 
111   STOP
      END 

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE BOTSUSY(GLB,GHB,GAB,XGLB,XGHB,XGAB,ASH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ICASE = 0
      CALL DMBAPP(ICASE,DGLB,DGHB,DGAB,ASH)
      XGLB = GLB*(1+DGLB)
      XGHB = GHB*(1+DGHB)
      XGAB = GAB*(1+DGAB)
      RETURN
      END

      SUBROUTINE DMBAPP(ICASE,DGLB,DGHB,DGAB,ASH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/PARAM/GF,AMW,AMZ
      COMMON/MASSES/AMS,AMC,AMB,AMT
      COMMON/HMASS/AMSM,AMA,AMHL,AMHH,AMCH,AMAR
      COMMON/GLUINO/AMG,AMSB1,AMSB2,STH,CTH,
     .              GLBB(2,2),GHBB(2,2),GABB(2,2),
     .              AMST1,AMST2,STHT,CTHT,
     .              GLTT(2,2),GHTT(2,2),GATT(2,2)
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      COMMON/BREAK/AMSQ,AMUR,AMDR,AU,AD,AMU,AM2
      PI = 4*DATAN(1.D0)
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      TANB = DTAN(B)
      TANA = DTAN(A)
      SB = TANB/DSQRT(1+TANB**2)
      AT = AU
      AB = AD
      RMTOP   = RUNM(AMT,6)
      HT = RMTOP/V/SB
      STOP1 = AMST1
      STOP2 = AMST2
      SBOT1 = AMSB1
      SBOT2 = AMSB2

      IF(ICASE.EQ.0)THEN
       DELTAMB = 2*ASH/3/PI*AMG*AMU*TANB*T(SBOT1,SBOT2,AMG)
     *         /(1-2*ASH/3/PI*AMG*AB*T(SBOT1,SBOT2,AMG))
     * + HT**2/(4*PI)**2*(AT-AMU/TANB)*AMU*TANB*
     * T(STOP1,STOP2,AMU)
       DGLB = -DELTAMB/(1+DELTAMB)*(1+1/TANA/TANB)
       DGHB = -DELTAMB/(1+DELTAMB)*(1-TANA/TANB)
       DGAB = -DELTAMB/(1+DELTAMB)*(1+1/TANB**2)
      ELSE
       DELTAMB = 2*ASH/3/PI*AMG*AMU*TANB*T(SBOT1,SBOT2,AMG)
       DGLB = -DELTAMB*(1+1/TANA/TANB)
       DGHB = -DELTAMB*(1-TANA/TANB)
       DGAB = -DELTAMB*(1+1/TANB**2)
      ENDIF
      RETURN
      END
