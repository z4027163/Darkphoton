C==============================================================
C
C                    ****************
C                    * Version 4.50 *
C                    ****************
C
C This is the program HIGLU, which calculates
C the total Higgs production cross section via gluon fusion
C including the next-to-leading QCD corrections at hadron
C colliders as well as the gluonic decay mode of neutral Higgs
C bosons at the same level of accuracy. It needs an input file
C defined as unit 95, which name has to be defined in the first
C OPEN statement of this program. The output is directed to
C unit 97. The program is based on the results of the following
C publication and references therein:
C
C M. Spira, A. Djouadi, D. Graudenz and P.M. Zerwas,
C Nucl. Phys. B453 (1995) 17.
C
C The program has been written by M. Spira and is documented
C in the publication
C
C M. Spira, Report DESY T-95-05 (October 1995), hep-ph/9510347
C
C In case of questions, problems or suggestions send an e-mail
C to michael.spira@psi.ch or spira@desy.de, please.
C
C M. Spira, January 13, 2016.
C
C ================ It uses as input parameters in higlu.in:
C (the rest of the input parameters is defined in hdecay.in - see comment
C card at beginning of hdecay.f.)
C
C  PROCESS: = 0: gg --> H
C           = 1:  H --> gg
C
C      ELW: = 0: no electroweak corrections
C           = 1: NLO electroweak corrections for SM Higgs
C           = 2: approximate mixed electroweak/QCD corrections for SM Higgs
C
C  DISTRIB: = 0: total cross section
C           = 1: d^2 sigma/dp_T/dy (LO)
C           = 2: d sigma/dp_T (LO)
C           = 3: sigma (p_Tmax > p_T > p_Tmin) (LO)
C
C   P_T:    transverse momentum [GeV] (DISTRIB = 1,2)
C     Y:    rapidity (DISTRIB = 1)
C   PTMIN:  minimal transverse momentum (DISTRIB = 3)
C   PTMAX:  maximal transverse momentum (DISTRIB = 3)
C
C COLLIDER: = 0: pp collider
C           = 1: ppbar collider
C
C   ENERGY: = 0: hadronic energy [TeV]
C
C       G_B: scaling factor of bottom Yukawa coupling
C       G_T: scaling factor of top Yukawa coupling
C       G_C: scaling factor of charm Yukawa coupling
C     G_SB1: scaling factor of Higgs coupling to sbottom_1-sbottom_1
C     G_SB2: scaling factor of Higgs coupling to sbottom_2-sbottom_2
C     G_ST1: scaling factor of Higgs coupling to stop_1-stop_1
C     G_ST2: scaling factor of Higgs coupling to stop_2-stop_2
C       C_G: Wilson coefficient of dim6 operator for a point-like Hgg coupling
C C_G-SCALE: Input scale of C_G
C
C This point-like coupling is defined in terms of the Lagrangian
C
C   L = alpha_s/pi c_g G^{a\mu\nu}G^a_{\mu\nu} H/v
C
C     SILH: = 0: non-linear HEFT
C           = 1: linear SILH-HEFT (omitting squares and products of dim6
C                operators)
C
C      G_BP: scaling factor of SM4 Higgs Yukawa coupling coupling to b'
C      G_TP: scaling factor of SM4 Higgs Yukawa coupling coupling to t'
C
C  TOPMASS: = 0: use top pole mass inside loops
C           = 1: use top MSbar mass inside loops
C
C      MU_1: scale contribution of top MSbar mass in units of M_Higgs
C      MU_2: constant scale contribution of top MSbar mass
C            The total scale is MU_1*M_Higgs + MU_2
C
C  BOTMASS: = 0: use bottom pole mass inside loops
C           = 1: use bottom MSbar mass inside loops
C
C      MU_1: scale contribution of bottom MSbar mass in units of M_Higgs
C      MU_2: constant scale contribution of bottom MSbar mass
C            The total scale is MU_1*M_Higgs + MU_2
C
C  CHMMASS: = 0: use charm pole mass inside loops
C           = 1: use charm MSbar mass inside loops
C
C      MU_1: scale contribution of charm MSbar mass in units of M_Higgs
C      MU_2: constant scale contribution of charm MSbar mass
C            The total scale is MU_1*M_Higgs + MU_2
C
C     TYPE: = 1: calculate cross section for the heavy scalar
C           = 2: calculate cross section for the pseudoscalar
C           = 3: calculate cross section for the light scalar
C                (default choice for the SM Higgs is TYPE = 1)
C
C INDIVIDU: = 0: Higgs masses are determined from M_A=M_HIGGS
C           = 1: Higgs masses are fixed by the input M_HIGGS
C
C   N_HIGGS: Higgs mass [GeV]
C
C      MU_1: renormalization scale contribution in units of M_Higgs
C      MU_2: constant renormalization scale contribution
C            The total scale is MU_1*M_Higgs + MU_2
C       Q_1: factorization scale contribution in units of M_Higgs
C       Q_2: constant factorization scale contribution
C            The total scale is Q_1*M_Higgs + Q_2
C
C     LOOP  = 1: (consistent) LO calculation
C           = 2: (consistent) NLO calculation
C           = 3: (consistent) NNLO calculation
C
C    CHOICE = 1: define alpha_s via alpha_s(M_Z)
C           = 2: define alpha_s via Lambda_NF
C
C     N_EXT: number of external light flavours included in H -> gg
C
C    ABSERR: absolute error anticipated for the Vegas integration
C    POINTS: number of points used for the Vegas integration
C     ITMAX: number of iterations used for the Vegas integration
C     PRINT: style of output from Vegas
C            =  0: no output of individual iterations
C            =  1: full output of individual iterations
C            = 10: output of individual iterations as table
C
C    SCHEME = 0: MSbar factorization scheme for PDFs
C           = 1: DIS factorization scheme for PDFs
C
C     STFUN = 0: use LHAPDF (up to version 5.9.1)
C           = 1: use GRV98 PDFs
C           = 2: use CTEQ6 PDFs
C
C   PDFNAME: name of LHAPDFset for LO, NLO and NNLO separately
C      NSET: set of LHAPDFset for LO, NLO and NNLO separately
C
C    ALPHAS: value of used alpha_s(M_Z) at LO, NLO and NNLO separately
C            (for CHOICE = 1)
C
C        NF: number of flavours of following values of lambda_NF
C    LAMBDA: values of lambda_NF [GeV] for LO, NLO and NNLO separately
C            (for CHOICE = 2)
C
C==============================================================
      PROGRAM HIGLU
      PARAMETER (NIN=95, NOUT=97, NIN0=98)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      CHARACTER*5 NAME,FNAME
      CHARACTER*12 OPNAME,IPNAME
      CHARACTER*1 HIGNAME(3)
      CHARACTER*100 PDFNAME1,PDFNAME2,PDFNAME3,PDFNAME,PATHNAME
      INTEGER NHIG(3)
      DIMENSION HIGH(24),HIGA(30)
      DOUBLE PRECISION YLTT(2,2),YLBB(2,2),YHTT(2,2),YHBB(2,2)
      DOUBLE PRECISION MST(2),MSB(2),MSL(2),MSU(2),MSD(2),MSE(2),MSN(2),
     .          MSN1(2),
     .          GCEN(2,2),GCTB(2,2),GLEE(2,2),GLTT(2,2),GLBB(2,2),
     .          GHEE(2,2),GHTT(2,2),GHBB(2,2)
      DOUBLE PRECISION ZMST(2),ZMSB(2),ZHTT(2,2),ZHBB(2,2)
      COMMON/PROCESS/IPROC
      COMMON/INTPOL/KORD
C     COMMON/STFU/ISET
      COMMON/PDFLIB/IPDFLIB,NGROUP,NSET,ISETERR
      COMMON/LUGG/VTAU,VSC
      COMMON/MASS/AMH,AMQ,S
      COMMON/SUSY/AMB1,AMC1,FACT,FACB,FACC,FACCG,FACCTG,ISILH
      COMMON/SUSYT/SCALCG
      COMMON/SUSY0/AMST1,AMST2,AMSB1,AMSB2,FACST1,FACST2,FACSB1,FACSB2
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/ETACUT/ETA0
      COMMON/GRENZ/IGRENZ
      COMMON/GRENZSQ/IGRENZSQ
      COMMON/PARTS/IPRTGG,IPRTGQ,ICHGG,ICHGQ
      COMMON/INTRO/INTEG
      COMMON/DCADRE/DO1,UP1,AERR1,RERR1
      COMMON/FACSC/ISCHEME
      COMMON/SMOOTH/DLT1,DLT2
      COMMON/ALS/XLAMBDA,AMAC,AMAB,AMAT,N0
      COMMON/ALSN3LO/N3LO
      COMMON/MASSES/AMS,AMC0,AMB0,AMT0
      COMMON/STRANGE/AMSB
      COMMON/PARSC/XKAPM,XKAPQ
      COMMON/FACSHIFT/RHOFAC
      COMMON/HIGGS/IHIGGS
      COMMON/FILE/NAME
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/SUSYP/GF,ALPHA,SW2,TGBET,AMTQ,AMZ,AMSQ
      COMMON/BREAK/AMSQ0,AMUR,AMDR,AU,AD,AMU,AM2
      COMMON/BREAKGLU/AMGLU
      COMMON/RESULT/VRES,VERR,DUM1,DUM2
      COMMON/QLIM/QMAX,QMIN
      COMMON/COLLI/ICOLL
      COMMON/BORN/CF0,CF0T,CF0B,CF0C,CF0G
      COMMON/BORN00/CF0G00
      COMMON/BORN4/CF0T4,CF0B4
      COMMON/CONST/ZETA2,ZETA3
      COMMON/PTY/PT,ETA,NLOOP
      COMMON/PTCUT/PTMIN,PTMAX
      COMMON/MODEL/IMODEL,ISUSY,I2HDM
      COMMON/SM4/AMT4,AMB4,FACT4,FACB4,ISM4,IGGELW
      COMMON/HMASS/AMSM,AMA,AMHL,AMHH,AMCH,AMAR
      COMMON/PARAM0/GF0,AMW0,AMZ0
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      COMMON/FLAGZ/IZG1,IN3LO
      COMMON/SLHA_vals_HDEC/islhai,islhao
      COMMON/HMASS_HDEC/AMSM_X,AMA_X,AML_X,AMH_X,AMCH_X,AMAR_X
      EXTERNAL DLUGG,DLUGQ
      EXTERNAL DGG,DGQ,DQQ
      EXTERNAL DVGG,DVGQ,DVQQ
      EXTERNAL EGG,EGQ,EQQ
      EXTERNAL EVGG,EVGQ,EVQQ
      EXTERNAL D9VGG
      EXTERNAL FPTYGG,FPTYGQ,FPTYQQ
      EXTERNAL FPTGG,FPTGQ,FPTQQ
      EXTERNAL TPTGG,TPTGQ,TPTQQ
      EXTERNAL D2VGG,D2VGQ,D2VQQB,D2VQQP,D2VQQ
c     RAT1(A1,A0,NF)=((33-2*NF)/12.D0+(153-19*NF)/24.D0*A1)
c    .              /((33-2*NF)/12.D0+(153-19*NF)/24.D0*A0)
c     RAT2(A1,A0,NF)=((33-2*NF)/12.D0+(153-19*NF)/24.D0*A1
c    .               +(2857-5033/9.D0*NF+325/27.D0*NF**2)/128.D0*A1**2)
c    .              /((33-2*NF)/12.D0+(153-19*NF)/24.D0*A0
c    .               +(2857-5033/9.D0*NF+325/27.D0*NF**2)/128.D0*A0**2)
      RAT1(A1,A0,NF)=1
      RAT2(A1,A0,NF)=1

      DATA HIGH/ 50.D0,100.D0,150.D0,200.D0,250.D0,300.D0,320.D0,
     .          340.D0,350.D0,360.D0,380.D0,400.D0,450.D0,500.D0,
     .          550.D0,600.D0,650.D0,700.D0,750.D0,800.D0,850.D0,
     .          900.D0,950.D0,1000.D0/
      DATA HIGA/ 50.D0,100.D0,150.D0,200.D0,250.D0,300.D0,320.D0,
     .          340.D0,348.D0,349.0D0,349.5D0,350.D0,350.5D0,351.D0,
     .          352.D0,360.D0,380.D0,400.D0,450.D0,500.D0,
     .          550.D0,600.D0,650.D0,700.D0,750.D0,800.D0,850.D0,
     .          900.D0,950.D0,1000.D0/

C--VERSION OF INPUT FILE: 0 = PPH  1 = HIGLU
      IVERSION = 1
 
C--OUTPUT-FILE
      OPNAME='higlu.out'
      OPEN(NOUT,FILE=OPNAME)

C--INITIALIZE RANDOM-GENERATOR
      CALL RSTART(12,34,56,78)
 
      PI=4.D0*DATAN(1.D0)
C--FACTOR GEV**(-2) --> PB
      GEVPB=389379338.D0
C--FERMI-CONSTANT
C     GF=1.16637D-5
C--ZETA(4)
      ZETA4 = PI**4/90.D0

C--INTIALIZE HDECAY COMMON BLOCKS
      CALL HDECINI(200.D0)
 
C--READ INPUT-FILE
      IF(IVERSION.EQ.0)THEN
C--INPUT-FILE
      OPEN(NIN,FILE='pph.in')
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)IPROC
      READ(NIN,100)IELW
      READ(NIN,100)ISQCD
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)NCOLL
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)W
c     READ(NIN,*)
c     READ(NIN,*)
c     READ(NIN,*)
c     READ(NIN,*)
c     READ(NIN,100)ISUSY
c     READ(NIN,100)ISM4
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)INDMH
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)FACT0
      READ(NIN,101)FACB0
      READ(NIN,101)FACC0
      READ(NIN,101)FACSB10
      READ(NIN,101)FACSB20
      READ(NIN,101)FACST10
      READ(NIN,101)FACST20
      READ(NIN,101)FACCG0
      READ(NIN,101)FACCTG0
      READ(NIN,101)SCALCG
      READ(NIN,100)ISILH0
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)FACB4
      READ(NIN,101)FACT4
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)ITOPSCHEME
      READ(NIN,101)TOPMU1
      READ(NIN,101)TOPMU2
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)IBOTSCHEME
      READ(NIN,101)BOTMU1
      READ(NIN,101)BOTMU2
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)ICHMSCHEME
      READ(NIN,101)DCHMMU1
      READ(NIN,101)DCHMMU2
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)NHIGGS
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)AMHBEG
      READ(NIN,101)AMHEND
      READ(NIN,100)MULTMH
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)ICHANNEL
      READ(NIN,100)ISOFT2L
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)ISCHEME
      READ(NIN,100)IPDFLIB
      READ(NIN,103)PATHNAME
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,103)PDFNAME1
      READ(NIN,100)NSET1
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,103)PDFNAME2
      READ(NIN,100)NSET2
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,103)PDFNAME3
      READ(NIN,100)NSET3
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)IALPS
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)ALSMZ1
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)ALSMZ2
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)ALSMZ3
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)N0
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)XLAMBDA1
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)XLAMBDA2
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)XLAMBDA3
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)IGRENZ0
      READ(NIN,100)IGRENZSQ0
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)NLOOP00
      READ(NIN,100)NFULL
      READ(NIN,100)NRESUM0
      READ(NIN,100)MCOLL0
      READ(NIN,100)IZG1
      READ(NIN,100)IN3LO
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)NFEXT
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)NRHO
      READ(NIN,101)RHOMBEG
      READ(NIN,101)RHOQBEG
      READ(NIN,101)RHOEND
      READ(NIN,100)MULTRHO
      READ(NIN,101)RHOFAC
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)INTEG
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)AERR
      READ(NIN,101)RERR
      READ(NIN,101)FACV
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)VAERR
      READ(NIN,100)IVPNT
      READ(NIN,100)IVITM
      READ(NIN,100)IVPRN
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)EPST
      READ(NIN,101)EPSV
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)REPS
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)DLT1
      READ(NIN,101)DLT2
      READ(NIN,100)KORD
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)NBER
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)IDIFF
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)IPTY
      READ(NIN,101)PTMIN
      READ(NIN,101)PTMAX
      READ(NIN,101)PTSTEP
      READ(NIN,101)ETA
      AMC = AMC0
      AMB = AMB0
      AMT = AMT0
      IF(FACCG0.NE.0.D0.OR.FACCTG.NE.0.D0) NRESUM0 = 0

      if(islhai.ne.0)then
       INDMH = 0
       AMHBEG = AMA
       AMHEND = AMA
      endif

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ETA0 = ETA
      PT = PTMIN
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     READ(55,*)NSET2
cc    NSET1 = NSET2
c     NSET3 = NSET2
c     FACB0 = 0
c     FACC0 = 0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      LOFULL = NFULL
      FACSQ = 1
      IF(FACSB1.EQ.0.D0.AND.FACSB2.EQ.0.D0.AND.
     .   FACST1.EQ.0.D0.AND.FACST2.EQ.0.D0) THEN
       FACSQ = 0
      ENDIF
      FACSM4 = 1
      IF(ISM4.EQ.0)THEN
       FACSM4 = 0
      ELSEIF(FACT4.EQ.0.D0.AND.FACB4.EQ.0.D0) THEN
       FACSM4 = 0
      ENDIF
      FACM = 1
      IF(FACSQ.EQ.0.D0.AND.FACSM4.EQ.0.D0) FACM = 0
      IGRZ = IGRENZ0 * IGRENZSQ0
      ILIMIT = 1
      IF(FACM.EQ.0.D0.AND.IGRENZ0.EQ.0) ILIMIT = 0
      IF(FACM.NE.0.D0.AND.IGRZ.EQ.0) ILIMIT = 0
c     IF(NLOOP00.LE.2.OR.ILIMIT.NE.0.OR.
      IF(NLOOP00.LE.2.OR.(FACT0.EQ.0.D0.AND.FACCG0.EQ.0.D0
     .   .AND.FACM.EQ.0.D0)) NFULL = 0
      IF(FACT0.EQ.0.D0.AND.FACCG0.EQ.0.D0.AND.FACM.EQ.0.D0
     .   .AND.NLOOP00.GT.2) NLOOP00 = 2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      N3LO0 = 0
      N3LO = 0
      IF(NLOOP00.GT.3)THEN
       NLOOP00 = 3
       N3LO0 = 1
      ENDIF
      IF(NRESUM0.NE.0)THEN
       N3LO0 = 1
      ENDIF
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      PDFNAME = PDFNAME2
      NSET    = NSET2
      ALSMZ   = ALSMZ2
      XLAMBDA = XLAMBDA2
      NLOOP   = NLOOP00
      N3LO = N3LO0
      NRESUM1 = NRESUM0
      NRESUM2 = NRESUM0
      MCOLL = MCOLL0
      IF(NLOOP00.EQ.1) PDFNAME = PDFNAME1
      IF(NLOOP00.EQ.1) NSET    = NSET1
      IF(NLOOP00.EQ.1) ALSMZ   = ALSMZ1
      IF(NLOOP00.EQ.1) XLAMBDA = XLAMBDA1
      IF(NFULL.EQ.0)THEN
       NLOOP = NLOOP00
       IF(NLOOP00.GE.3) PDFNAME = PDFNAME3
       IF(NLOOP00.GE.3) NSET    = NSET3
       IF(NLOOP00.GE.3) ALSMZ   = ALSMZ3
       IF(NLOOP00.GE.3) XLAMBDA = XLAMBDA3
      ELSE
       NLOOP = 2
       N3LO = 0
      ENDIF
      IGRENZ   = IGRENZ0
      IGRENZSQ = IGRENZSQ0
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      IF(ISM4.EQ.1)THEN
       ISUSY = 0
       ISQCD = 0
      ENDIF
      SQCD_TOPR = 0
      SQCD_TOPI = 0
      SQCD_BOTR = 0
      SQCD_BOTI = 0
      DSQCD_TOPR = 0
      DSQCD_TOPI = 0
      DSQCD_BOTR = 0
      DSQCD_BOTI = 0
      IF(ISQCD.EQ.1)THEN
       OPEN(NIN0,FILE='susyqcd.in')
       READ(NIN0,*)
       READ(NIN0,*)SQCD_TOPR,DSQCD_TOPR
       READ(NIN0,*)SQCD_TOPI,DSQCD_TOPI
       READ(NIN0,*)
       READ(NIN0,*)SQCD_BOTR,DSQCD_BOTR
       READ(NIN0,*)SQCD_BOTI,DSQCD_BOTI
       CLOSE(NIN0)
      ELSEIF(ISQCD.EQ.2)THEN
       OPEN(NIN0,FILE='babis.in')
       READ(NIN0,*)ZMST(1),ZMST(2),ZMSB(1),ZMSB(2)
       READ(NIN0,*)ZHT,ZHB
       READ(NIN0,*)ZHTT(1,1),ZHTT(1,2)
       READ(NIN0,*)ZHTT(2,1),ZHTT(2,2)
       READ(NIN0,*)ZHBB(1,1),ZHBB(1,2)
       READ(NIN0,*)ZHBB(2,1),ZHBB(2,2)
       READ(NIN0,*)
       READ(NIN0,*)
       READ(NIN0,*)SQCD_TOPR,DSQCD_TOPR
       READ(NIN0,*)SQCD_TOPI,DSQCD_TOPI
       READ(NIN0,*)SQCDA_TOPR
       READ(NIN0,*)SQCDA_TOPI
       READ(NIN0,*)
       READ(NIN0,*)SQCD_BOTR,DSQCD_BOTR
       READ(NIN0,*)SQCD_BOTI,DSQCD_BOTI
       READ(NIN0,*)SQCDA_BOTR
       READ(NIN0,*)SQCDA_BOTI
       READ(NIN0,*)
       READ(NIN0,*)ZMH
       CLOSE(NIN0)
      ENDIF
      IF(NCOLL.EQ.1)ICOLL=-1
      IF(NCOLL.EQ.0)ICOLL=1
      IMODEL = ISUSY
      AMS = AMSB
      AMB1 = AMB
      AMB00 = AMB1
      AMC1 = AMC
      AMC00 = AMC1
      NHIG(1) = INT(NHIGGS/100)
      NHIG(2) = INT((NHIGGS-100*NHIG(1))/10)
      NHIG(3) = INT(NHIGGS-100*NHIG(1)-10*NHIG(2))
      SW2 = 1.D0-AMW0**2/AMZ0**2
      GF0 = GF
      AMW = AMW0
      AMZ = AMZ0
      ALPHA=1.D0/128.9D0
      AMSQ0=AMSQ
      AMQBEG = AMT
      AMQEND = AMT
      MULTMQ = 1
      IF(NRHO.EQ.2)THEN
        RHOBEG=RHOQBEG
      ELSE
        RHOBEG=RHOMBEG
      ENDIF
C--DEFINE NUMBER OF MASSLESS FLAVORS
      NFD=NFEXT

      ELSE

C--INPUT-FILE
      OPEN(NIN,FILE='higlu.in')

C--READ INPUT-FILE
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)IPROC0
      READ(NIN,100)IELW
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)IDIST0
      READ(NIN,101)PT0
      READ(NIN,101)Y0
      READ(NIN,101)PTMIN0
      READ(NIN,101)PTMAX0
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)NCOLL
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)W
c     READ(NIN,*)
c     READ(NIN,*)
c     READ(NIN,*)
c     READ(NIN,*)
c     READ(NIN,100)ISUSY
c     READ(NIN,100)ISM4
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)FACT0
      READ(NIN,101)FACB0
      READ(NIN,101)FACC0
      READ(NIN,101)FACSB10
      READ(NIN,101)FACSB20
      READ(NIN,101)FACST10
      READ(NIN,101)FACST20
      READ(NIN,101)FACCG0
c     READ(NIN,101)FACCTG0
      FACCTG = 0
      READ(NIN,101)SCALCG
      READ(NIN,100)ISILH0
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)FACB4
      READ(NIN,101)FACT4
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)ITOPSCHEME
      READ(NIN,101)TOPMU1
      READ(NIN,101)TOPMU2
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)IBOTSCHEME
      READ(NIN,101)BOTMU1
      READ(NIN,101)BOTMU2
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)ICHMSCHEME
      READ(NIN,101)DCHMMU1
      READ(NIN,101)DCHMMU2
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)IHIGGS0
      READ(NIN,100)INDMH
      READ(NIN,101)AMH
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)S1REN
      READ(NIN,101)S2REN
      READ(NIN,101)S1FAC
      READ(NIN,101)S2FAC
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)NLOOP00
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)IALPS
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)NFEXT
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)VAERR
      READ(NIN,100)IVPNT
      READ(NIN,100)IVITM
      READ(NIN,100)IVPRN
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)ISCHEME
      READ(NIN,100)IPDFLIB
      READ(NIN,103)PATHNAME
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,103)PDFNAME1
      READ(NIN,100)NSET1
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,103)PDFNAME2
      READ(NIN,100)NSET2
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,103)PDFNAME3
      READ(NIN,100)NSET3
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)ALSMZ1
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)ALSMZ2
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)ALSMZ3
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,100)N0
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)XLAMBDA1
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)XLAMBDA2
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,101)XLAMBDA3
      AMC = AMC0
      AMB = AMB0
      AMT = AMT0

C--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C--CORRECTIONS BELOW THIS LINE AT YOUR OWN RISK!!!!!!!
C--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C--ENERGY IN GEV
      W = 1000.D0*W

C--NUMBER OF FLAVORS FOR LAMBDA
      IF(N0.GT.6.OR.N0.LT.3)THEN
       WRITE(NOUT,*)'N_F NOT POSSIBLE. TAKING DEFAULT N_F = 4'
       N0 = 4
      ENDIF

C--ELECTROWEAK CORRECTIONS ONLY FOR SM
      IF(IELW.NE.0)THEN
       IF(ISUSY.NE.0.OR.FACT0.NE.FACB0.OR.
     .     DABS(FACSB10)+DABS(FACSB20)+DABS(FACST10)
     .    +DABS(FACST20).NE.0.D0)THEN
        WRITE(NOUT,*)'ELW. CORRECTIONS NOT AVAILABLE. TAKING ELW = 0'
        IELW = 0
       ENDIF
      ENDIF
      IF(ISM4.EQ.1.AND.ISUSY.NE.0)THEN
      WRITE(NOUT,*)'4TH SUSY-GENERATION NOT AVAILABLE. TAKING MODEL = 0'
       ISUSY = 0
      ENDIF

C--TRANSLATE FLAG FOR PROCESS
      IF(IPROC0.EQ.0)THEN
       IPROC=1
      ELSEIF(IPROC0.EQ.1) THEN
       IPROC=0
C--NUMBER OF EXTERNAL LIGHT FLAVORS FOR HIGGS --> GG
       IF(NFEXT.GT.5.OR.NFEXT.LT.3)THEN
        WRITE(NOUT,*)'N_FEXT NOT POSSIBLE. TAKING DEFAULT N_FEXT = 5'
        NFEXT = 5
       ENDIF
      ELSE
       WRITE(NOUT,*)'PROCESS NOT POSSIBLE. TAKING DEFAULT PROCESS = 0'
       IPROC=1
      ENDIF

C--ADJUST YUKAWA COUPLINGS OF SM
C     IF(ISUSY.EQ.0)THEN
C      FACT=1.D0
C      FACB=1.D0
C      FACC=1.D0
C      FACCG=0.D0
C      IHIGGS0=1
C      ISILH=1
C     ENDIF

      NHIG(1)=0
      NHIG(2)=0
      NHIG(3)=0
C--TRANSLATE FLAG FOR HIGGS TYPE
      IF(IHIGGS0.EQ.1) THEN
       IHIGGS=3
      ELSEIF(IHIGGS0.EQ.2) THEN
       IHIGGS=1
      ELSEIF(IHIGGS0.EQ.3) THEN
c      IF(ISUSY.NE.0.AND.ISUSY.NE.3)THEN
       IF(ISUSY.NE.0.OR.I2HDM.NE.0)THEN
        IHIGGS=2
       ELSE
        WRITE(NOUT,*)'HIGGS TYPE NOT POSSIBLE. TAKING DEFAULT TYPE = 1'
        IHIGGS=3
       ENDIF
      ELSE
       WRITE(NOUT,*)'HIGGS TYPE NOT POSSIBLE. TAKING DEFAULT TYPE = 1'
       IHIGGS=3
      ENDIF
      NHIG(IHIGGS)=1

C--TRANSLATE FLAG FOR COLLIDER
      IF(NCOLL.EQ.1) THEN
       ICOLL=-1
      ELSEIF(NCOLL.EQ.0)THEN
       ICOLL=1
      ELSE
       WRITE(NOUT,*)'COLLIDER NOT POSSIBLE. TAKING DEFAULT COLLDER = 0'
       ICOLL=1
      ENDIF

C--DEFINE NUMBER OF MASSLESS FLAVORS
      NFD=NFEXT

C--FLAG FOR PRINTING VEGAS ITERATIONS
      IF(IVPRN.NE.0.AND.IVPRN.NE.1.AND.IVPRN.NE.10)THEN
       WRITE(NOUT,*)'PRINT NOT POSSIBLE. TAKING DEFAULT PRINT = 0'
       IVPRN = 0
      ENDIF

C--INITIALIZATION OF FLAGS AND VARIABLES
      ISQCD = 0
      ICHANNEL=1111
      ISOFT2L=0
      IGRENZ0=0
      IGRENZSQ0=0
      AMQBEG=AMT
      AMQEND=AMT
      MULTMQ=1
      AMHBEG=AMH
      AMHEND=AMH
      MULTMH=1
      RHOFAC=1.D0
      NRHO=1
      IMODEL = ISUSY
      AMS = AMSB
      AMB1 = AMB
      AMB00 = AMB1
      AMC1 = AMC
      AMC00 = AMC1
      SW2 = 1.D0-AMW0**2/AMZ0**2
      GF0 = GF
      AMW = AMW0
      AMZ = AMZ0
      ALPHA=1.D0/128.90D0
      AMSQ0=AMSQ
      EPST=1.D-8
      EPSV=1.D-8
      REPS=1.D-8
      KORD=2
      NBER=18
      INTEG=3
      IDIFF=0
      IF(IDIST0.EQ.0)THEN
       IPTY=0
      ELSEIF(IDIST0.EQ.1)THEN
       IPTY=3
      ELSEIF(IDIST0.EQ.2)THEN
       IPTY=2
      ELSE
       IPTY=1
      ENDIF
C     PTMIN = PT0
C     PTMAX = PT0
      PTMIN = PTMIN0
      PTMAX = PTMAX0
      PTSTEP = 100.D0
      PT = PT0
      ETA = Y0

      if(islhai.ne.0)then
       INDMH = 0
       AMHBEG = AMA
       AMHEND = AMA
      endif

      IZG1 = 0
      IN3LO = 0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(FACCG0.NE.0.D0.OR.FACCTG.NE.0.D0) NRESUM0 = 0
      NRESUM0 = 0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      NFULL = 1
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     IGRENZ0=1
c     NFULL = 0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      LOFULL = NFULL
      FACSQ = 1
      IF(FACSB1.EQ.0.D0.AND.FACSB2.EQ.0.D0.AND.
     .   FACST1.EQ.0.D0.AND.FACST2.EQ.0.D0) THEN
       FACSQ = 0
      ENDIF
      FACSM4 = 1
      IF(ISM4.EQ.0)THEN
       FACSM4 = 0
      ELSEIF(FACT4.EQ.0.D0.AND.FACB4.EQ.0.D0) THEN
       FACSM4 = 0
      ENDIF
      FACM = 1
      IF(FACSQ.EQ.0.D0.AND.FACSM4.EQ.0.D0) FACM = 0
      IGRZ = IGRENZ0 * IGRENZSQ0
      ILIMIT = 1
      IF(FACM.EQ.0.D0.AND.IGRENZ0.EQ.0) ILIMIT = 0
      IF(FACM.NE.0.D0.AND.IGRZ.EQ.0) ILIMIT = 0
c     IF(NLOOP00.LE.2.OR.ILIMIT.NE.0.OR.
      IF(NLOOP00.LE.2.OR.(FACT0.EQ.0.D0.AND.FACCG0.EQ.0.D0
     .   .AND.FACM.EQ.0.D0)) NFULL = 0
      IF(FACT0.EQ.0.D0.AND.FACCG0.EQ.0.D0.AND.FACM.EQ.0.D0
     .   .AND.NLOOP00.GT.2) NLOOP00 = 2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      N3LO0 = 0
      N3LO = 0
      IF(NLOOP00.GT.3)THEN
       NLOOP00 = 3
       N3LO0 = 1
      ENDIF
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      PDFNAME = PDFNAME2
      NSET = NSET2
      ALSMZ   = ALSMZ2
      XLAMBDA = XLAMBDA2
      NLOOP   = NLOOP00
      N3LO = N3LO0
      IF(NLOOP00.EQ.1) PDFNAME = PDFNAME1
      IF(NLOOP00.EQ.1) NSET    = NSET1
      IF(NLOOP00.EQ.1) ALSMZ   = ALSMZ1
      IF(NLOOP00.EQ.1) XLAMBDA = XLAMBDA1
      IF(NFULL.EQ.0)THEN
       NLOOP = NLOOP00
       IF(NLOOP00.GE.3) PDFNAME = PDFNAME3
       IF(NLOOP00.GE.3) NSET    = NSET3
       IF(NLOOP00.GE.3) ALSMZ   = ALSMZ3
       IF(NLOOP00.GE.3) XLAMBDA = XLAMBDA3
      ELSE
       NLOOP = 2
       N3LO = 0
      ENDIF
      IGRENZ   = IGRENZ0
      IGRENZSQ = IGRENZSQ0
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C     IF(IPROC.EQ.1)THEN
        QMIN=0.D0
        QMAX=1.D20
C     ENDIF

      ENDIF
C--END READING

C--INITIALIZE PDFSET
C     CALL PDFSET(PATHNAME,PDFNAME)

      FACT   = FACT0
      FACB   = FACB0
      FACC   = FACC0
      FACCG  = FACCG0
      V = 1/DSQRT(DSQRT(2.D0)*GF)
      FACCTG = FACCTG0 * AMT**2/V**2
      ISILH  = ISILH0
      IF(IHIGGS.EQ.1) FACCG = 0
      IF(IHIGGS.EQ.1) FACCTG = 0
      IF(IHIGGS.EQ.1) ISILH = 0
      IF(ISUSY.EQ.0) THEN
       FACST1 = 0
       FACST2 = 0
       FACSB1 = 0
       FACSB2 = 0
      ELSE
       FACST1 = FACST10
       FACST2 = FACST20
       FACSB1 = FACSB10
       FACSB2 = FACSB20
      ENDIF

      IVDPT=1
      IF(IPROC.EQ.1)THEN
       IF(IGRENZ.EQ.1)THEN
        IVDIM=2
        IVDPT=1
       ELSE
        IVDIM=3
        IVDPT=1
       ENDIF
      ELSE
       IVDIM=2
      ENDIF

      QMIN=0.D0
      QMAX=1.D50
      IF(IPROC.EQ.1)THEN
       IF(INTEG.NE.3)THEN
C--INITIALIZE GRID FOR LUMINOSITIES
        CALL LUMINI
       ELSE
        QMIN=0.D0
        QMAX=1.D50
       ENDIF
      ENDIF

C--INITIALIZE COEFFICIENTS FOR POLYLOGARITHMS
      CALL BERNINI(NBER)

      HIGNAME(1)='A'
      HIGNAME(2)='h'
      HIGNAME(3)='H'
 
C--TOTAL ENERGY
      S=W**2

      ICHVIRT=INT(ICHANNEL/1000)
      ICHGG=INT(ICHANNEL/100)-10*ICHVIRT
      ICHGQ=INT(ICHANNEL/10)-100*ICHVIRT-10*ICHGG
      ICHQQ=ICHANNEL-1000*ICHVIRT-100*ICHGG-10*ICHGQ
      IPRTGG=0
      IPRTGQ=0
 
      IF(IVERSION.EQ.0)THEN
      IF(IPROC.EQ.1)THEN
       IF(IPDFLIB.EQ.1) WRITE(NOUT,1201)ISCHEME
       IF(IPDFLIB.EQ.2) WRITE(NOUT,1202)ISCHEME
c      IF(IPDFLIB.EQ.0) WRITE(NOUT,120)PDFNAME(7:100),NSET,ISCHEME
       IF(IPDFLIB.EQ.0)THEN
        IF(NFULL.EQ.0)THEN
         WRITE(NOUT,120)PDFNAME(1:100),NSET,ISCHEME
        ELSE
         WRITE(NOUT,1205)PDFNAME1(1:100),NSET1,ISCHEME
         WRITE(NOUT,1206)PDFNAME2(1:100),NSET2,ISCHEME
         WRITE(NOUT,1207)PDFNAME3(1:100),NSET3,ISCHEME
        ENDIF
       ENDIF
       WRITE(NOUT,*)
      ENDIF
      WRITE(NOUT,113)EPST,EPSV
      IF(INTEG.EQ.1)THEN
        WRITE(NOUT,*)
        WRITE(NOUT,161)
        WRITE(NOUT,160)
        WRITE(NOUT,114)AERR,RERR
      ELSE
        WRITE(NOUT,*)
        WRITE(NOUT,162)
        WRITE(NOUT,160)
        WRITE(NOUT,165)IVPNT,IVITM
      ENDIF
      WRITE(NOUT,*)
      WRITE(NOUT,115)IGRENZ
      WRITE(NOUT,180)
      WRITE(NOUT,181)RHOFAC
      WRITE(NOUT,*)
      WRITE(NOUT,117)DLT1,DLT2
      WRITE(NOUT,*)
      WRITE(NOUT,126)REPS,NLOOP00+N3LO0
      WRITE(NOUT,*)
      ELSE
        WRITE(NOUT,*)
        WRITE(NOUT,5162)
        WRITE(NOUT,5160)
        WRITE(NOUT,5166)VAERR
        WRITE(NOUT,5165)IVPNT,IVITM
        WRITE(NOUT,*)
      ENDIF
 
      IF(MULTMQ.GT.1)THEN
        DLTMQ=(AMQEND-AMQBEG)/DFLOAT(MULTMQ-1)
      ELSE
        DLTMQ=0.D0
      ENDIF
 
      IF(MULTMH.GT.1)THEN
        DLTMH=(AMHEND-AMHBEG)/DFLOAT(MULTMH-1)
      ELSE
        DLTMH=0.D0
      ENDIF
 
      DO 9999 IQ=1,MULTMQ
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     DO 9999 ITGB=1,4
c      IF(ITGB.EQ.1)TGBET = 5
c      IF(ITGB.EQ.2)TGBET = 10
c      IF(ITGB.EQ.3)TGBET = 20
c      IF(ITGB.EQ.4)TGBET = 30
c      AU = DSQRT(6.D0)*AMSQ + AMU/TGBET
c      AD = AU
c      write(6,*)AU
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
C--TOP MASS
        AMQ=AMQBEG+DLTMQ*(IQ-1)
        AMTQ=AMQ
        AMT=AMQ
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C       AMQ = AMQ*1.D2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C--MASSES FOR ALPHA_S (TOP QUARK DECOUPLED)
        AMAC=AMC0
        AMAB=AMB0
        AMAT=AMT0
c       ACC=1.D-8
        ACC=1.D-10
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c  LO alpha_s
        XLAMBDA=XITLA(1,ALSMZ,ACC)
        N0=5
        CALL ALSINI(ACC)
        AMAT=AMT0*1.D8
        ALS_LO = ALPHAS(125.D0,1)
c NLO alpha_s
        XLAMBDA=XITLA(2,ALSMZ,ACC)
        N0=5
        CALL ALSINI(ACC)
        AMAT=AMT0*1.D8
        ALS_NLO = ALPHAS(125.D0,2)
c       write(6,*)XLAMBDA,ALS_NLO
c NNLO alpha_s
        XLAMBDA=XITLA(3,ALSMZ3,ACC)
        N0=5
        CALL ALSINI(ACC)
        AMAT=AMT0*1.D8
        ALS_NNLO = ALPHAS(125.D0,3)
c       write(6,*)XLAMBDA,ALS_NNLO
c N3LO alpha_s
        XLAMBDA=XITLA(4,ALSMZ3,ACC)
        N0=5
        CALL ALSINI(ACC)
        AMAT=AMT0*1.D8
        ALS_N3LO = ALPHAS(125.D0,4)
c       ALS_N3LO = ALPHAS(AMZ,4)
c       write(6,*)XLAMBDA,ALS_N3LO
c       write(6,*)ALS_NNLO,ALS_N3LO
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        IF(IALPS.EQ.1) THEN
         XLAMBDA=XITLA(NLOOP,ALSMZ,ACC)
         N0=5
        ENDIF
        AMAT=AMT0
        CALL ALSINI(ACC)
        AMAT=AMT0*1.D8
c     write(6,*)AMT0,AMAT
c     write(6,*)'alpha_s: ',alphas(amz,nloop),amz,nloop
c     write(6,*)'alpha_s: ',alphas(amh/2,nloop),amh/2,nloop
c     write(6,*)'alpha_s: ',alphas(amh,nloop),amh,nloop
        ALS0MZ = ALPHAS(AMZ,NLOOP)
c       ALSMZ = ALS0MZ
      IF(IVERSION.EQ.0)THEN
      IF(IPROC.EQ.1)THEN
        IF(NFULL.EQ.0)THEN
         WRITE(NOUT,118)W,N0,XLAMBDA
         WRITE(NOUT,*)
         WRITE(NOUT,1460)ALS0MZ
        ELSE
         WRITE(NOUT,11801)W,N0,XLAMBDA1
         WRITE(NOUT,11802)N0,XLAMBDA2,N0,XLAMBDA3
         WRITE(NOUT,*)
         WRITE(NOUT,1461)ALSMZ1
         WRITE(NOUT,1462)ALSMZ2
         WRITE(NOUT,1463)ALSMZ3
        ENDIF
      ELSE
        WRITE(NOUT,146)NFD,N0,XLAMBDA
        WRITE(NOUT,*)
      ENDIF
        WRITE(NOUT,*)
        WRITE(NOUT,119)AMB,AMQ,AMC
        IF(ISM4.NE.0)WRITE(NOUT,1191)AMB4,AMT4
        IF(FACCG.NE.0.D0.OR.FACCTG.NE.0.D0)THEN
         WRITE(NOUT,*)
         IF(ISILH.EQ.0)THEN
          WRITE(NOUT,1189)
         ELSE
          WRITE(NOUT,1188)
         ENDIF
        ENDIF
        WRITE(NOUT,*)
      ELSE
      IF(IPROC.EQ.1)THEN
        WRITE(NOUT,1184)
        WRITE(NOUT,1185)
        WRITE(NOUT,*)
        IF(ICOLL.EQ.1)THEN
         WRITE(NOUT,1180)
         WRITE(NOUT,1182)
        ELSE
         WRITE(NOUT,1181)
         WRITE(NOUT,1183)
        ENDIF
        WRITE(NOUT,5118)W/1000.D0
        WRITE(NOUT,*)
      ELSE
        WRITE(NOUT,1186)
        WRITE(NOUT,1187)
        WRITE(NOUT,*)
        WRITE(NOUT,5146)NFEXT
        WRITE(NOUT,*)
      ENDIF
        IF(NLOOP00.GE.1) WRITE(NOUT,1260)N0,XLAMBDA1,ALSMZ1
        IF(NLOOP00.GE.2) WRITE(NOUT,1261)N0,XLAMBDA2,ALSMZ2
        IF(NLOOP00.GE.3) WRITE(NOUT,1262)N0,XLAMBDA3,ALSMZ3
        WRITE(NOUT,*)
        WRITE(NOUT,5119)AMQ,AMB
        WRITE(NOUT,51190)AMC
        IF(ISM4.NE.0)WRITE(NOUT,51191)AMB4,AMT4
        IF(FACCG.NE.0.D0.OR.FACCTG.NE.0.D0)THEN
         WRITE(NOUT,*)
         IF(ISILH.EQ.0)THEN
          WRITE(NOUT,1189)
         ELSE
          WRITE(NOUT,1188)
         ENDIF
        ENDIF
        WRITE(NOUT,*)
        IF(ISUSY.EQ.1)THEN
         WRITE(NOUT,51291)TGBET
         WRITE(NOUT,51290)
         WRITE(NOUT,*)
        ELSEIF(ISUSY.EQ.2)THEN
         WRITE(NOUT,51292)TGBET
         WRITE(NOUT,51290)
         WRITE(NOUT,*)
        ELSEIF(ISUSY.EQ.3)THEN
         WRITE(NOUT,51293)TGBET
         WRITE(NOUT,51290)
         WRITE(NOUT,*)
        ENDIF
      ENDIF
C--DECAY WIDTH USING ALPHA_S WITH HEAVY QUARKS DECOUPLED
       IF(IPROC.EQ.0)THEN
C--USING ALPHA_S WITH HEAVY QUARKS DECOUPLED
        IF(NFD.EQ.5)THEN
         AMAT = AMT0*1.D8
        ELSEIF(NFD.EQ.4)THEN
         AMAT = AMT0*1.D8
         AMAB = AMB0*1.D8
        ELSEIF(NFD.EQ.3)THEN
         AMAT = AMT0*1.D8
         AMAB = AMB0*1.D8
         AMAC = AMC0*1.D8
        ENDIF
       ENDIF
 
        IF(ISUSY.NE.0)THEN
          IHIGBEG=1
          IHIGEND=3
        ELSE
          IHIGBEG=1
          IHIGEND=3
        ENDIF
 
      DO 9998 IGLOBAL=IHIGBEG,IHIGEND
       IHIGGS=IGLOBAL
       IFLAG = NHIG(IGLOBAL)
       IF(IFLAG.NE.0) THEN
        WRITE(NOUT,130)HIGNAME(IHIGGS)
        WRITE(NOUT,112)
 
      DO 9997 IH=1,MULTMH
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     DO 9997 IH=1,7
c      IF(IH.EQ.1)AMH=125
c      IF(IH.EQ.2)AMH=125.5D0
c      IF(IH.EQ.3)AMH=126
c      IF(IH.EQ.4)AMH=150
c      IF(IH.EQ.5)AMH=200
c      IF(IH.EQ.6)AMH=250
c      IF(IH.EQ.7)AMH=300
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
C--HIGGS MASS
      AMH=AMHBEG+DLTMH*(IH-1)
      AMSM = AMH
      AMA = AMH

      NRESUM1 = NRESUM0
      NRESUM2 = NRESUM0
      FACB = FACB0
      FACC = FACC0

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     IF(ISM4.NE.0) THEN
c      AMT4 = AMB4 + 50*(1+DLOG(AMH/115)/5)
c      WRITE(NOUT,1191)AMB4,AMT4
c     ENDIF
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      amb = amb00
      amc = amc00

      NLOOP0 = NLOOP
      IF(NLOOP0.EQ.1)THEN
       ALSMZ0 = ALSMZ
       XLAMBDA0 = XLAMBDA
       N00=N0
       NLOOP = 3
       ALSMZ = 0.118D0
       XLAMBDA=XITLA(NLOOP,ALSMZ,ACC)
c      XLAMBDA=0.226D0
       N0=5
       AMAT=AMT0
       CALL ALSINI(ACC)
       AMAT=AMT0*1.D8
c     write(6,*)AMT0,AMAT
c     write(6,*)'alpha_s: ',alphas(amz,nloop),amz,nloop
c     write(6,*)'alpha_s: ',alphas(amh/2,nloop),amh/2,nloop
c     write(6,*)'alpha_s: ',alphas(amh,nloop),amh,nloop
      ENDIF

C--CALCULATE PSEUDOSCALAR HIGGS MASS FROM AMH
      IF(ISUSY.NE.0)THEN
       IF(INDMH.EQ.1)THEN
        CALL AMHAMA_HDEC(IHIGGS-1,AMH,TGBET)
        AMA = AMA_X
       ELSE
        AMA=AMH
        CALL HDECTRAFO1
        CALL SUSYCP_HDEC(TGBET)
        CALL HDECTRAFO2
        IF(IHIGGS.EQ.1)THEN
         AMH = AMA
        ELSEIF(IHIGGS.EQ.2)THEN
         AMH = AMHL
        ELSEIF(IHIGGS.EQ.3)THEN
         AMH = AMHH
        ENDIF
       ENDIF
      ENDIF
C--INTIALIZE HDECAY COMMON BLOCKS
      CALL HDECINI(AMA)
      AMAT=AMT0*1.D8

      FACT = FACT0
      FACB = FACB0
      FACC = FACC0
      FACCG= FACCG0
      V = 1/DSQRT(DSQRT(2.D0)*GF)
      FACCTG = FACCTG0 * AMT**2/V**2
      ISILH= ISILH0
      IF(IHIGGS.EQ.1) FACCG = 0
      IF(IHIGGS.EQ.1) FACCTG = 0
      IF(IHIGGS.EQ.1) ISILH = 0
      IF(I2HDM.NE.0)THEN
       IF(IHIGGS.EQ.1)THEN
        FACT = FACT0*GAT
        FACB = FACB0*GAB
        FACC = FACC0*GAT
        FACELW = 0
       ELSEIF(IHIGGS.EQ.2)THEN
        FACT = FACT0*GLT
        FACB = FACB0*GLB
        FACC = FACC0*GLT
        FACELW = GLVV
       ELSE
        FACT = FACT0*GHT
        FACB = FACB0*GHB
        FACC = FACC0*GHT
        FACELW = GHVV
       ENDIF
      ENDIF
      if(itopscheme.eq.1)then
       topscale = topmu1*amh + topmu2
       rmt = runm_hdec(topscale,6)
c      fact = fact0 * rmt/amt
       amt = rmt
       amtq=amt
       amq=amt
       write(nout,6000)topscale,rmt
      endif
      if(ibotscheme.eq.1)then
       botscale = botmu1*amh + botmu2
       rmb = runm_hdec(botscale,5)
c      facb = facb0 * rmb/amb
       amb = rmb
       write(nout,6001)botscale,rmb
       amb1 = amb
       amb0 = amb
       amb00= amb
      endif
      if(ichmscheme.eq.1)then
       dchmscale = dchmmu1*amh + dchmmu2
       rmc = runm_hdec(dchmscale,4)
c      facc = facc0 * rmc/amc
       amc = rmc
       write(nout,6002)dchmscale,rmc
       amc1 = amc
       amc0 = amc
       amc00= amc
      endif

C--CALCULATE SUSY-COUPLINGS
      IF(ISUSY.NE.0)THEN
       CALL HDECTRAFO1
       CALL SUSYCP_HDEC(TGBET)
       CALL HDECTRAFO2
       AMEL = 1000.D0
       AMER = 1000.D0
       AL = 0
       TSC = (AMSQ+AMUR+AMDR)/3
       BSC = (AMSQ+AMUR+AMDR)/3
       CALL SFERMION_HDEC(TSC,BSC,AMSQ,AMUR,AMDR,AMEL,AMER,AL,AU,AD,AMU,
     .                   MST,MSB,MSL,MSU,MSD,MSE,MSN,MSN1,
     .                   GLEE,GLTT,GLBB,GHEE,GHTT,GHBB,
     .                   GAEE,GATT,GABB,GCEN,GCTB)
       AMST1 = MST(1)
       AMST2 = MST(2)
       AMSB1 = MSB(1)
       AMSB2 = MSB(2)
       TSC = (AMSQ+AMUR+AMDR)/3
       BSC = (AMSQ+AMUR+AMDR)/3
       CALL SFERMION_HDEC(TSC,BSC,AMSQ,AMUR,AMDR,AMEL,AMER,AL,AU,AD,AMU,
     .                    MST,MSB,MSL,MSU,MSD,MSE,MSN,MSN1,
     .                    GLEE,GLTT,GLBB,GHEE,GHTT,GHBB,
     .                    GAEE,GATT,GABB,GCEN,GCTB)
       do i=1,2
        do j=1,2
         yltt(i,j)=gltt(i,j)*amz**2
         ylbb(i,j)=glbb(i,j)*amz**2
         yhtt(i,j)=ghtt(i,j)*amz**2
         yhbb(i,j)=ghbb(i,j)*amz**2
        enddo
       enddo
c      write(6,*)'m_hh=  ',amhh
c      write(6,*)'g_ht:  ',ght
c      write(6,*)'g_hb:  ',ghb
c      write(6,*)'m_st:  ',amst1,amst2
c      write(6,*)'m_sb:  ',amsb1,amsb2
c      write(6,*)'gl_st: ',yltt
c      write(6,*)'gl_sb: ',ylbb
c      write(6,*)'gh_st: ',yhtt
c      write(6,*)'gh_sb: ',yhbb
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      write(6,*)amhl,glt,glb
c      write(6,*)amhh,ght,ghb
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       QSUSY = (AMSB1+AMSB2+AMGLU)/3
       CALL DMBAPP_HDEC(1,DGLB,DGHB,DGAB,QSUSY,1)
       DMBFAC = ALPHAS_HDEC(QSUSY,3)/PI
       CALL BOTSUSY_HDEC(GLB,GHB,GAB,XGLB,XGHB,XGAB,QSUSY,1)
c      GLB = XGLB
c      GHB = XGHB
c      GAB = XGAB
       if(ibotscheme.eq.2) then
        susyscalt = (mst(1)+mst(2))/2
        ratt = runm_hdec(susyscalt,6)/amq
        glt = glt*ratt
        ght = ght*ratt
        gat = gat*ratt
        susyscalb = (msb(1)+msb(2))/2
        ratb = runm_hdec(susyscalb,5)/amb
        glb = xglb*ratb
        ghb = xghb*ratb
        gab = xgab*ratb
c       write(6,*)'stop: ',susyscalt,runm_hdec(susyscalt,6),amq
c       write(6,*)'sbot: ',susyscalb,runm_hdec(susyscalb,5),amb
       endif
       if(isqcd.eq.2)then
c--Babis' conventions
        amhl = zmh
        amhh = zmh
        ama  = zmh
        glt = zht
        glb = zhb
        ght = zht
        ghb = zhb
        gat = zht
        gab = zhb
        mst(1) = zmst(1)
        mst(2) = zmst(2)
        msb(1) = zmsb(1)
        msb(2) = zmsb(2)
        amst1 = mst(1)
        amst2 = mst(2)
        amsb1 = msb(1)
        amsb2 = msb(2)
        do i=1,2
         do j=1,2
          gltt(i,j) = zhtt(i,j)/amz**2
          ghtt(i,j) = zhtt(i,j)/amz**2
c         gatt(i,j) = zhtt(i,j)/amz**2
          glbb(i,j) = zhbb(i,j)/amz**2
          ghbb(i,j) = zhbb(i,j)/amz**2
c         gabb(i,j) = zhbb(i,j)/amz**2
         enddo
        enddo
c       write(6,*)'Babis: ',mst(1),mst(2),msb(1),msb(2)
c       write(6,*)zhtt(1,1),zhtt(1,2),zhtt(2,2)
c       write(6,*)zhbb(1,1),zhbb(1,2),zhbb(2,2)
       endif
       IF(IHIGGS.EQ.1)THEN
        AMH = AMA
        FACT = FACT0*GAT
        FACB = FACB0*GAB
        FACC = FACC0*GAT
        DGDMB = DGAB/DMBFAC
       ELSEIF(IHIGGS.EQ.2)THEN
        AMH = AMHL
        FACT = FACT0*GLT
        FACB = FACB0*GLB
        FACC = FACC0*GLT
        FACST1 = FACST10*GLTT(1,1)*AMZ**2/MST(1)**2
        FACST2 = FACST20*GLTT(2,2)*AMZ**2/MST(2)**2
        FACSB1 = FACSB10*GLBB(1,1)*AMZ**2/MSB(1)**2
        FACSB2 = FACSB20*GLBB(2,2)*AMZ**2/MSB(2)**2
        DGDMB = DGLB/DMBFAC
       ELSEIF(IHIGGS.EQ.3)THEN
        AMH = AMHH
        FACT = FACT0*GHT
        FACB = FACB0*GHB
        FACC = FACC0*GHT
        FACST1 = FACST10*GHTT(1,1)*AMZ**2/MST(1)**2
        FACST2 = FACST20*GHTT(2,2)*AMZ**2/MST(2)**2
        FACSB1 = FACSB10*GHBB(1,1)*AMZ**2/MSB(1)**2
        FACSB2 = FACSB20*GHBB(2,2)*AMZ**2/MSB(2)**2
        DGDMB = DGHB/DMBFAC
       ENDIF
      ELSE
        TGBET=1.D0
c       FACT = 1
c       FACB = 1
c       FACC = 1
c       FACCG = 0
c       FACST1 = 0
c       FACST2 = 0
c       FACSB1 = 0
c       FACSB2 = 0
c       ISILH = 0
        AMST1 = 1000.D0
        AMST2 = 1000.D0
        AMSB1 = 1000.D0
        AMSB2 = 1000.D0
      ENDIF
      FACT00 = FACT
      FACB00 = FACB
      FACC00 = FACC
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c       FACST1 = 0
c       FACST2 = 0
c       FACSB1 = 0
c       FACSB2 = 0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      IF(NLOOP0.EQ.1)THEN
       ALSMZ = ALSMZ0
       XLAMBDA = XLAMBDA0
       NLOOP = NLOOP0
       N0=N00
       AMAT=AMT0
       CALL ALSINI(ACC)
       AMAT = AMT0*1.D8
      ENDIF

      AMAT = AMT0
c     write(6,*)AMT0,AMAT
c     write(6,*)'alpha_s: ',alphas(amz,nloop),amz,nloop
c     write(6,*)'alpha_s: ',alphas(amh/2,nloop),amh/2,nloop
c     write(6,*)'alpha_s: ',alphas(amh,nloop),amh,nloop
c     write(6,*)'alpha_s: ',alphas(amh*100,nloop),amh*100,nloop
c     write(6,*)'alpha_s: ',alphas(amh/100,nloop),amh/100,nloop
      AMAT = AMT0*1.D8
c     write(6,*)AMT0,AMAT
c     write(6,*)'alpha_s: ',alphas(amz,nloop),amz,nloop
c     write(6,*)'alpha_s: ',alphas(amh/2,nloop),amh/2,nloop
c     write(6,*)'alpha_s: ',alphas(amh,nloop),amh,nloop
c     write(6,*)'alpha_s: ',alphas(amh*100,nloop),amh*100,nloop
c     write(6,*)'alpha_s: ',alphas(amh/100,nloop),amh/100,nloop
c     write(6,*)'alpha_s: ',als_nlo,2
c     write(6,*)'alpha_s: ',als_lo,1

      IF(IVERSION.NE.0)THEN
C--SCALES
       SCREN = S1REN*AMH + S2REN
       SCFAC = S1FAC*AMH + S2FAC
       RHOMBEG=SCREN/AMH
       RHOQBEG=SCFAC/AMH
       RHOEND=RHOMBEG
       MULTRHO=1
       IF(NRHO.EQ.2)THEN
         RHOBEG=RHOQBEG
       ELSE
         RHOBEG=RHOMBEG
       ENDIF
      ENDIF
      IF(MULTRHO.GT.1)THEN
        DLTRHO=(RHOEND-RHOBEG)/DFLOAT(MULTRHO-1)
      ELSE
        DLTRHO=0.D0
      ENDIF
 
      DO 9997 IR=1,MULTRHO
C--SCALE-FACTOR IN UNITS OF HIGGS MASS
      XKAP=RHOBEG+DLTRHO*(IR-1)
      IF(NRHO.EQ.1)THEN
C--SCALE-FACTOR FOR RENORMALIZATION SCALE
        XKAPM=XKAP
C--SCALE-FACTOR FOR FACTORIZATION SCALE
        XKAPQ=RHOQBEG 
      ELSEIF(NRHO.EQ.2)THEN
        XKAPM=RHOMBEG 
        XKAPQ=XKAP
      ELSE
        XKAPM=XKAP
        XKAPQ=XKAP
      ENDIF
 
      NRESUM1 = NRESUM0
      NRESUM2 = NRESUM0
      FACB = FACB00
      FACC = FACC00

      TH=AMH**2/S
C--CALCULATE SCALES: QM = RENORMALIZATION-SCALE
C--                  QQ = FACTORIZATION-SCALE
      CALL SCALES(QM,QQ)
 
C--DEFINE LIMITS OF INTEGRATIONS
      IF(IPROC.EQ.1)THEN
C--INTEGRATION OVER TAU
       DO=TH
      ELSE
       DO=EPST
      ENDIF
C--INTEGRATION OVER SCATTERING ANGLE
      UP=1.D0-EPST
      DO1=EPSV
      UP1=1.D0-EPSV
      AERR1=AERR*FACV
      RERR1=RERR*FACV
 
C--CALCULATE LOWEST ORDER FORM FACTOR
      IF(IHIGGS.NE.1)THEN
C--SCALAR HIGGS
C--TOTAL FORM FACTOR
        CF0=FACT*CFBORN(AMH,AMQ)+FACB*CFBORN(AMH,AMB)
     .     +FACC*CFBORN(AMH,AMC)+12*FACCG
     .     +FACCTG*CFBORNT(AMH,AMQ)
     .     +FACST1*CSBORN(AMH,AMST1)+FACST2*CSBORN(AMH,AMST2)
     .     +FACSB1*CSBORN(AMH,AMSB1)+FACSB2*CSBORN(AMH,AMSB2)
        CF0S = CF0
        write(6,*)'R_LO = ',CF0/FACT
C--TOP FORM FACTOR
        CF0T=FACT*CFBORN(AMH,AMQ)
        CF0ST1=FACST1*CSBORN(AMH,AMST1)
        CF0ST2=FACST2*CSBORN(AMH,AMST2)
c       write(6,*)'ratio = ',AMH,AMQ,cdabs(CFBORN(AMH,AMQ))**2
C--BOTTOM FORM FACTOR
        CF0B=FACB*CFBORN(AMH,AMB)
        CF0SB1=FACSB1*CSBORN(AMH,AMSB1)
        CF0SB2=FACSB2*CSBORN(AMH,AMSB2)
C--CHARM FORM FACTOR
        CF0C=FACC*CFBORN(AMH,AMC)
C--DIM6 FORM FACTOR
        CF0G00=12*FACCG
        CF0G0=12*FACCG+FACCTG*CFBORNT(AMH,AMQ)
        CF0G=CF0G0
        CF0TP=(FACT-1)*CFBORN(AMH,AMQ)
        CF0BP=(FACB-1)*CFBORN(AMH,AMB)
        CF0CP=(FACC-1)*CFBORN(AMH,AMC)
        CF0P=(FACT-1)*CFBORN(AMH,AMQ)+(FACB-1)*CFBORN(AMH,AMB)
     .      +(FACC-1)*CFBORN(AMH,AMC)
        if(igrenz.eq.1) CF0P=(FACT-1)*CFBORN(AMH,AMQ)
        ratnnlo = cdabs(cf0t+cf0g)**2/cdabs(cf0t+cf0b+cf0c+cf0g)**2
        if(isilh.ne.0)then
         ratnnlo = (cdabs(cf0t+cf0g)**2-cdabs(cf0p+cf0g)**2)
     .           / (cdabs(cf0t+cf0b+cf0c+cf0g)**2-cdabs(cf0p+cf0g)**2)
         if(igrenz.eq.1) ratnnlo = 1
        endif
        if(ism4.ne.0)then
         cf0 = cf0 + fact4*cfborn(amh,amt4) + facb4*cfborn(amh,amb4)
         cf0t4=fact4*cfborn(amh,amt4)
         cf0b4=facb4*cfborn(amh,amb4)
         ratnnlo = cdabs(cf0t+cf0t4+cf0b4+cf0g)**2
     .           / cdabs(cf0t+cf0b+cf0c+cf0t4+cf0b4+cf0g)**2
        if(isilh.ne.0)
     .   ratnnlo=(cdabs(cf0t+cf0t4+cf0b4+cf0g)**2-cdabs(cf0p+cf0g)**2)
     .          / (cdabs(cf0t+cf0b+cf0c+cf0t4+cf0b4+cf0g)**2
     .             -cdabs(cf0p+cf0g)**2)
        endif
        ctrun = 0
        cbrun = 0
        ccrun = 0
        if(itopscheme.eq.1)then
         dshift = dlog(topscale**2/amt**2) + 4/3.d0
         ctrun=-cdct(amh,amt)*dshift
c        write(6,*)'als: ',alphas(91.187d0,nloop)
c        write(6,*)'top: ',cdct(amh,amt),dshift,alphas(qm,nloop)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c        hh = 1.d-4
c        cdum = (cfborn(amh,amt*(1+hh))-cfborn(amh,amt*(1-hh)))/2/hh
c    .        / cfborn(amh,amt)
c        write(6,*)'top: ',(cdct(amh,amt)/2/cdum)
c        cdum = (cfborn(amh,amb*(1+hh))-cfborn(amh,amb*(1-hh)))/2/hh
c    .        / cfborn(amh,amb)
c        write(6,*)'bot: ',(cdct(amh,amb)/2/cdum)
c        cdum = (cfborn(amh,amc*(1+hh))-cfborn(amh,amc*(1-hh)))/2/hh
c    .        / cfborn(amh,amc)
c        write(6,*)'chm: ',(cdct(amh,amc)/2/cdum)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        endif
        if(ibotscheme.eq.1)then
         dshift = dlog(botscale**2/amb**2) + 4/3.d0
         cbrun=-cdct(amh,amb)*dshift
c        write(6,*)'bot: ',cdct(amh,amb),dshift,alphas(qm,nloop)
c        write(6,*)'dum: ',-cdct(amh,amb)*dshift
c        write(6,*)'dum: ',-cdct(amh,amb)*dshift*alphas(qm,nloop)/pi
c        write(6,*)'run: ',cbrun
c        write(6,*)'run: ',cbrun*alphas(qm,nloop)/pi
c        write(6,*)'cf0t: ',cf0t
c        write(6,*)'cf0b: ',cf0b
c        write(6,*)'cf0:  ',cf0
        elseif(ibotscheme.eq.2)then
         susyscalt = (mst(1)+mst(2))/2
         susyscalb = (msb(1)+msb(2))/2
         susysct = (mst(1)+mst(2)+dabs(amglu))/3
         susyscb = (msb(1)+msb(2)+dabs(amglu))/3
         ctrun = dlog(susyscalt**2/amq**2) + 4/3.d0
         cbrun = dlog(susyscalb**2/amb**2) + 4/3.d0
         ccrun = dlog(susyscalt**2/amc**2) + 4/3.d0
c        write(6,*)'stop: ',susyscalt,ctrun,alphas_hdec(susysct,3)
c        write(6,*)'sbot: ',susyscalb,cbrun,alphas_hdec(susyscb,3)
c        write(6,*)
        else
        endif
        if(ichmscheme.eq.1)then
         dshift = dlog(dchmscale**2/amc**2) + 4/3.d0
         ccrun=-cdct(amh,amc)*dshift
c        write(6,*)'chm: ',cdct(amh,amc),dshift,alphas(qm,nloop)
        endif
      ELSE
C--PSEUDOSCALAR HIGGS
        CF0=1.5D0*(FACT*CFBORNA(AMH,AMQ)+FACB*CFBORNA(AMH,AMB)
     .            +FACC*CFBORNA(AMH,AMC))
        CF0S = CF0
        CF0G0=0
        CF0G=CF0G0
C--TOP FORM FACTOR
        CF0T=1.5D0*FACT*CFBORNA(AMH,AMQ)
C--BOTTOM FORM FACTOR
        CF0B=1.5D0*FACB*CFBORNA(AMH,AMB)
C--CHARM FORM FACTOR
        CF0C=1.5D0*FACC*CFBORNA(AMH,AMC)
        CF0ST1 = 0
        CF0ST2 = 0
        CF0SB1 = 0
        CF0SB2 = 0
        ratnnlo = cdabs(cf0t)**2/cdabs(cf0t+cf0b+cf0c)**2
        if(ism4.ne.0)then
         cf0 = cf0 + fact4*cfborna(amh,amt4) + facb4*cfborna(amh,amb4)
         cf0t4=fact4*cfborna(amh,amt4)
         cf0b4=facb4*cfborna(amh,amb4)
         ratnnlo = cdabs(cf0t+cf0t4+cf0b4)**2
     .           / cdabs(cf0t+cf0b+cf0c+cf0t4+cf0b4)**2
        endif
        ctrun = 0
        cbrun = 0
        ccrun = 0
        if(itopscheme.eq.1)then
         dshift = dlog(topscale**2/amt**2) + 4/3.d0
         ctrun=-cdcta(amh,amt)*dshift
        endif
        if(ibotscheme.eq.1)then
         dshift = dlog(botscale**2/amb**2) + 4/3.d0
         cbrun=-cdcta(amh,amb)*dshift
        elseif(ibotscheme.eq.2)then
         susyscalt = (mst(1)+mst(2))/2
         susyscalb = (msb(1)+msb(2))/2
         ctrun = dlog(susyscalt**2/amq**2) + 4/3.d0
         cbrun = dlog(susyscalb**2/amb**2) + 4/3.d0
         ccrun = dlog(susyscalt**2/amc**2) + 4/3.d0
c        write(6,*)'stop: ',susyscalt,ctrun,alphas(susyscalt,2)
c        write(6,*)'sbot: ',susyscalb,cbrun,alphas(susyscalb,2)
c        write(6,*)
        endif
        if(ichmscheme.eq.1)then
         dshift = dlog(dchmscale**2/amc**2) + 4/3.d0
         ccrun=-cdcta(amh,amc)*dshift
        endif
      ENDIF

      IF(IVERSION.EQ.0)THEN
      IF(IPROC.EQ.1)THEN
       WRITE(NOUT,121)AMH,XKAPM,XKAPQ
      ELSE
       WRITE(NOUT,1210)AMH,XKAPM
      ENDIF
      IF(ISUSY.NE.0)THEN
        WRITE(NOUT,1294)AMA
        WRITE(NOUT,129)FACB,FACT,TGBET
        WRITE(NOUT,1295)FACC,FACCG
        WRITE(NOUT,1296)FACCTG0
        IF(IHIGGS.NE.1.AND.
     .  DABS(FACSB1)+DABS(FACSB2)+DABS(FACST1)+DABS(FACST2).NE.0.D0)THEN
         WRITE(NOUT,1290)MST(1),MST(2)
         WRITE(NOUT,1291)FACST1,FACST2
         WRITE(NOUT,1292)MSB(1),MSB(2)
         WRITE(NOUT,1293)FACSB1,FACSB2
        ENDIF
      ELSE
        WRITE(NOUT,132)FACB,FACT,FACC
        WRITE(NOUT,1320)FACCG,FACCTG0
        IF(ISM4.NE.0) WRITE(NOUT,1321)FACB4,FACT4
        IF(IHIGGS.NE.1.AND.
     .  DABS(FACSB1)+DABS(FACSB2)+DABS(FACST1)+DABS(FACST2).NE.0.D0)THEN
         WRITE(NOUT,1290)AMST1,AMST2
         WRITE(NOUT,1291)FACST1,FACST2
         WRITE(NOUT,1292)AMSB1,AMSB2
         WRITE(NOUT,1293)FACSB1,FACSB2
        ENDIF
      ENDIF
      ELSE
      WRITE(NOUT,5121)HIGNAME(IHIGGS),AMH
      IF(ISUSY.NE.0.AND.IHIGGS.NE.1)THEN
       WRITE(NOUT,1219)AMA
      ENDIF
      IF(IPROC.EQ.1)THEN
       IF(IPTY.EQ.0)THEN
        WRITE(NOUT,1211)SCREN,SCFAC
       ELSE
        AMHT = DSQRT(AMH**2+PT**2)
        DUMR = S1REN*AMHT + S2REN
        DUMF = S1FAC*AMHT + S2FAC
        WRITE(NOUT,1211)DUMR,DUMF
       ENDIF
      ELSE
       WRITE(NOUT,1212)SCREN
      ENDIF
      WRITE(NOUT,5132)HIGNAME(IHIGGS),FACB,HIGNAME(IHIGGS),FACT
      WRITE(NOUT,51320)HIGNAME(IHIGGS),FACC,HIGNAME(IHIGGS),FACCG
      WRITE(NOUT,51322)HIGNAME(IHIGGS),FACCTG0
      IF(ISM4.NE.0)
     .WRITE(NOUT,51321)HIGNAME(IHIGGS),FACB4,HIGNAME(IHIGGS),FACT4
      IF(IHIGGS.NE.1.AND.
     . DABS(FACSB1)+DABS(FACSB2)+DABS(FACST1)+DABS(FACST2).NE.0.D0)THEN
       WRITE(NOUT,5133)HIGNAME(IHIGGS),FACST1,HIGNAME(IHIGGS),FACST2
       WRITE(NOUT,5134)HIGNAME(IHIGGS),FACSB1,HIGNAME(IHIGGS),FACSB2
       WRITE(NOUT,11290)AMST1,AMST2
       WRITE(NOUT,11292)AMSB1,AMSB2
      ENDIF
      ENDIF
      IF(IPROC.EQ.1)THEN
C--LOWEST ORDER CROSS SECTION
        SIG0=GEVPB*GF/288.D0/DSQRT(2.D0)/PI*CDABS(CF0)**2
        IF(ISILH.NE.0)
     .   SIG0=GEVPB*GF/288.D0/DSQRT(2.D0)/PI*(CDABS(CF0)**2
     .                                       -CDABS(CF0P+CF0G)**2)
        IF(INTEG.EQ.3)THEN
          VTAU=TH
          VSC=QQ
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          IF(NLOOP00.NE.1.AND.LOFULL.NE.0)THEN
           NLOOP = 1
           N3LO = 0
           ALSMZ = ALSMZ1
           XLAMBDA = XLAMBDA1
           PDFNAME = PDFNAME1
           NSET = NSET1
           IF(IALPS.EQ.1) THEN
            XLAMBDA=XITLA(NLOOP,ALSMZ,ACC)
            N0=5
           ENDIF
           AMAT=AMT0
           CALL ALSINI(ACC)
           AMAT=AMT0*1.D8
           CALL PDFSET(PATHNAME,PDFNAME)
           CALL VEGAS(DLUGG,VAERR,1,IVPNT,IVITM,IVPRN,1)
           DLUM0=TH*VRES*ALPHAS(QM,NLOOP)**2
           SIGBORN00=SIG0*DLUM0
           DSIGB00=SIG0*DLUM0/VRES*VERR
           NLOOP = 2
           N3LO = 0
           ALSMZ = ALSMZ2
           XLAMBDA = XLAMBDA2
           PDFNAME = PDFNAME2
           NSET = NSET2
           IF(IALPS.EQ.1) THEN
            XLAMBDA=XITLA(NLOOP,ALSMZ,ACC)
            N0=5
           ENDIF
           AMAT=AMT0
           CALL ALSINI(ACC)
           AMAT=AMT0*1.D8
           CALL PDFSET(PATHNAME,PDFNAME)
           CALL VEGAS(DLUGG,VAERR,1,IVPNT,IVITM,IVPRN,1)
           DLUM0=TH*VRES*ALPHAS(QM,NLOOP)**2
           SIGBORN=SIG0*DLUM0
           DSIGB=SIG0*DLUM0/VRES*VERR
          ELSE
           CALL PDFSET(PATHNAME,PDFNAME)
           CALL VEGAS(DLUGG,VAERR,1,IVPNT,IVITM,IVPRN,1)
           DLUM0=TH*VRES*ALPHAS(QM,NLOOP)**2
           SIGBORN=SIG0*DLUM0
           DSIGB=SIG0*DLUM0/VRES*VERR
           SIGBORN00 = SIGBORN
           DSIGB00   = DSIGB
          ENDIF
          ratnlo = sigborn/sigborn00
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         IF(IVERSION.EQ.0)THEN
          WRITE(NOUT,1700)SIGBORN00,DABS(DSIGB00)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c         WRITE(6,*)CDABS(CF0)**2,CDABS(CF0)**2-CDABS(CF0G)**2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         ELSE
          WRITE(NOUT,5170)SIGBORN00,DABS(DSIGB00)
         ENDIF
        ELSE
          DLUM0=TH*DLUMGG(TH,QQ)*ALPHAS(QM,NLOOP)**2
          SIGBORN=SIG0*DLUM0
          DSIGB=0.D0
          WRITE(NOUT,122)SIGBORN
        ENDIF
C--COUPLING OF RADIATIVE CORRECTIONS
        DCPL0=ALPHAS(QM,NLOOP)**3/PI*TH
      ELSE
C--LOWEST ORDER DECAY WIDTH
        GAM0=1.D3*GF/36.D0/DSQRT(2.D0)*ALPHAS(QM,NLOOP)**2/PI**3
     .       *AMH**3*CDABS(CF0)**2
        IF(ISILH.NE.0)
     .   GAM0=1.D3*GF/36.D0/DSQRT(2.D0)*ALPHAS(QM,NLOOP)**2/PI**3
     .       *AMH**3*(CDABS(CF0)**2-CDABS(CF0P+CF0G)**2)
         IF(IVERSION.EQ.0)THEN
          WRITE(NOUT,140)GAM0
         ELSE
          WRITE(NOUT,5140)GAM0
         ENDIF
      ENDIF
C--LOGARITHMS OF SCALES
      SLG0=DLOG(QM**2/AMH**2)
      SLG1=DLOG(QQ**2/AMH**2)
      SLGT=DLOG(QM**2/AMT0**2)

c!!!!!!!
c     SLG1= 0
c     SLG0= 0
c!!!!!!!
 
C***************************************************************
C         SOFT GLUON RESUMMATION
C***************************************************************

      XK9GG=0.D0
      DK9GG=0.D0
      IF(ISOFT2L.NE.0)THEN
C--2-LOOP SOFT GLUON EXPANSION
       INTEG0=INTEG
       INTEG=3
       XK9GG=0.D0
       DK9GG=0.D0
       NC=3
       B0=(11*NC-2*DNF(QM))/12.D0
C--VIRTUAL CORRECTIONS
C--LIMIT OF HEAVY TOP MASS
C--SCALAR HIGGS
        XC9VIRT=NC*(NC/6.D0*PI**4-12*NC*ZETA4-2909/216.D0*B0)
     .         +NC*(SLG0**2*(-2*NC*ZETA2)
     .             +SLG0*(2*B0*ZETA2 - 8*NC*ZETA3))
       IF(IHIGGS.EQ.1)THEN
C--PSEUDOSCALAR HIGGS
        XC9VIRT=XC9VIRT
     .         + NC**2/9.D0*6.D0*(6.D0/2+PI**2)
     .         + 6*NC*B0/3.D0*SLG0
       ENDIF
       IF(IPROC.EQ.1)THEN
C--GLUON FUSION
        XK9VIRT=(XC9VIRT)*(ALPHAS(QM,NLOOP)/PI)**2
        WRITE(NOUT,500)XK9VIRT,0.D0
        IF(ISOFT2L.EQ.2)THEN
C--SCALAR HIGGS
         DEL = 203.D0/12.D0
         PK9VT=(NC**2/9.D0*DEL*(DEL/2+PI**2)
     .        +NC*B0*2909.D0/216.D0
     .        +NC*(SLG0**2*(11.D0/12.D0*B0+121.D0/72.D0*NC)
     .            +SLG0*(203.D0/36.D0*B0+11.D0/3.D0*NC*ZETA2
     .                  +2233.D0/216.D0*NC)))
     .        *(ALPHAS(QM,NLOOP)/PI)**2
       IF(IHIGGS.EQ.1)THEN
C--PSEUDOSCALAR HIGGS
         PK9VT=PK9VT
     .        +NC*(203.D0/18.D0*NC + 11.D0/3.D0*NC*SLG0)
     .        *(ALPHAS(QM,NLOOP)/PI)**2
       ENDIF
         WRITE(NOUT,501)PK9VT
        ENDIF
C--REAL CORRECTIONS
        ICALL1=0
        ICALL2=0
        IFAIL=0
        IF66=0
        I9DIM=2
        IPRTGG=0
C--GG --> HG
         CALL VEGAS(D9VGG,VAERR,I9DIM,IVPNT,IVITM,IVPRN,1)
         IER=-1
         XK9GG=VRES
         DK9GG=VERR/DLUM0
C--REST-TERMS FROM PLUS-DISTRIBUTIONS
         P0 =  DLOG(1.D0-TH)
         P1 = -DLOG(1.D0-TH)**2/2
         P2 =  DLOG(1.D0-TH)**3/3
         P3 = -DLOG(1.D0-TH)**4/4
         REST = 8*NC**2*P3 - 4*NC*B0*P2 - 4/3.D0*NC**2*PI**2*P1
     .        + 16*NC**2*ZETA3*P0
     .        + NC*(SLG0**2*(4*NC*P1 - B0*P0)
     .        - SLG0*(12*NC*P2 - 4*B0*P1 - 4*NC*ZETA2*P0))
         IF(IHIGGS.EQ.1)REST=REST+6*4/3.D0*NC**2*(P1-P0/2*SLG0)
         XK9GG=XK9GG+REST*DLUM0*(ALPHAS(QM,NLOOP)/PI)**2
         XK9GG=XK9GG/DLUM0
         DK9GG=XK9GG*DSQRT((DK9GG/XK9GG)**2+(DSIGB/SIGBORN)**2)
         IGG=IER
         WRITE(NOUT,502)XK9GG,DK9GG,IGG,ICALL1,ICALL2,IFAIL,IF66
         IF(ISOFT2L.EQ.2)THEN
          XT9GG=XK9GG
          DT9GG=DK9GG
          DO 511 I=1,3
           IF(I.EQ.3.AND.SLG0.EQ.0.D0)THEN
            XI9GG = 0
            DI9GG = 0
           ELSE
           ICALL1=0
           IPRTGG=I
           CALL VEGAS(D9VGG,VAERR,IVDPT,IVPNT,IVITM,IVPRN,1)
           IER=-1
           IF(VRES.EQ.0) THEN
            XI9GG=0
            DI9GG=0
           ELSE
            XI9GG=VRES
            DI9GG=VERR/DLUM0
            XI9GG=XI9GG/DLUM0
            DI9GG=XI9GG*DSQRT((DI9GG/XI9GG)**2+(DSIGB/SIGBORN)**2)
           ENDIF
C--REST-TERMS FROM PLUS-DISTRIBUTIONS
           IF(I.EQ.1)THEN
             REST = 8*NC**2*P3 - 4*NC*B0*P2 - 4/3.D0*NC**2*PI**2*P1
     .            + 16*NC**2*ZETA3*P0
     .            + NC*(SLG0**2*(4*NC*P1 - B0*P0)
     .            - SLG0*(12*NC*P2 - 4*B0*P1 - 4*NC*ZETA2*P0))
             IF(IHIGGS.EQ.1)REST=REST+6*4/3.D0*NC**2*(P1-P0/2*SLG0)
             XI9GG=XI9GG+REST*(ALPHAS(QM,NLOOP)/PI)**2
           ELSEIF(I.EQ.2)THEN
C--REST-TERMS FROM PLUS-DISTRIBUTIONS
             DEL = 203.D0/12.D0
             REST=DEL*4/3.D0*NC**2*P1
     .           + NC*(SLG0**2*(-11.D0/3.D0*NC*P0)
     .             - SLG0*(-22.D0/3.D0*NC*P1 + 203.D0/18.D0*NC*P0))
             XI9GG=XI9GG+REST*(ALPHAS(QM,NLOOP)/PI)**2
           ELSEIF(I.EQ.3)THEN
C--REST-TERMS FROM PLUS-DISTRIBUTIONS
             PGG1 = 3.D0/2*(- 20/9.D0*P0)
     .            + 9*(67/9.D0-2*ZETA2)*P0
c!!!!!!!!!!!!!!!!!!!!!!!
c            PGG1 = 0
c!!!!!!!!!!!!!!!!!!!!!!!
             REST = SLG0*(- PGG1)/2
             XI9GG=XI9GG+REST*(ALPHAS(QM,NLOOP)/PI)**2
C--DELTA-TERMS FROM AP-SPLITTING FUNCTION AND ALPHA_S
             TF = DNF(QM)/2
             B1 = (153 - 19*DNF(QM))/24.D0
             DELTA = 9*(8/3.D0 + 3*ZETA3) - 4/3.D0*TF - 4/3.D0*NC*TF
c!!!!!!!!!!!!!!!!!!!!!!!
c            B1 = 0
c            DELTA = 0
c!!!!!!!!!!!!!!!!!!!!!!!
             DELTA = (2*B1-DELTA/2)*SLG0*(ALPHAS(QM,NLOOP)/PI)**2
             XI9GG=XI9GG+DELTA
           ENDIF
           ENDIF
           IGG=IER
           IF(I.EQ.1)THEN 
             WRITE(NOUT,503)XI9GG,DI9GG,IGG,ICALL1
           ELSEIF(I.EQ.2)THEN 
             WRITE(NOUT,504)XI9GG,DI9GG,IGG,ICALL1
           ELSEIF(I.EQ.3)THEN 
             WRITE(NOUT,505)XI9GG,DI9GG,IGG,ICALL1
           ENDIF
511       CONTINUE
          IPRTGG=0
         ENDIF
         XK9GG = XK9VIRT + XK9GG
         WRITE(NOUT,506)XK9GG,DK9GG
       ENDIF
       INTEG=INTEG0
      ENDIF
C***************************************************************
C      END OF SOFT GLUON RESUMMATION
C***************************************************************

c     write(6,*)AMT0,AMAT
c     write(6,*)'alpha_s = ',ALPHAS(AMZ,NLOOP),AMZ,NLOOP
c     write(6,*)'alpha_s = ',ALPHAS(QM,NLOOP),QM,NLOOP
 
      XKVIRT=0.D0
      XKGG=0.D0
      XKGQ=0.D0
      XKQQ=0.D0
      DKVIRT=0.D0
      DKGG=0.D0
      DKGQ=0.D0
      DKQQ=0.D0
      IPRTGG=0
      IPRTGQ=0
      IF(IPTY.EQ.0)THEN
C*****************************************
C--TOTAL CROSS SECTIONS
C*****************************************************************
C
C IPART: = 0: SIG_TBC^NLO WITH NLO PDF+ALPHA_S
C        = 1: SIG_T(INFTY)^NLO WITH NLO PDF+ALPHA_S
C        = 2: SIG_T(INFTY)^NNLO WITH NNLO PDF+ALPHA_S
C
C*****************************************************************
      IF(NLOOP.NE.1.AND.ICHANNEL.NE.0)THEN
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       IPART = 0
c--NLL running of point-like Higgs coupling to gluons
       CF0G = CF0G0
     .      + CF0G00*(RAT1(ALPHAS(QM,NLOOP),ALPHAS(SCALCG,NLOOP),5)-1)
       CF0X = CF0T+CF0G
       CF0 = CF0S - CF0G0 + CF0G
       CF0XP = CF0TP+CF0G
       SIG0=GEVPB*GF/288.D0/DSQRT(2.D0)/PI*CDABS(CF0)**2
       IF(ISILH.NE.0)
     .  SIG0=GEVPB*GF/288.D0/DSQRT(2.D0)/PI*(CDABS(CF0)**2
     .                                      -CDABS(CF0P+CF0G)**2)
       IF(ISM4.NE.0) CF0X = CF0T+CF0T4+CF0B4+CF0G
c      RATX = CDABS(CF0X)**2/CDABS(CF0)**2
c      IF(ISILH.NE.0) RATX = (CDABS(CF0X)**2-CDABS(CF0P+CF0G)**2)
c    .                     / (CDABS(CF0)**2-CDABS(CF0P+CF0G0)**2)
c----------------------------------------------------------------
9991   CONTINUE
       IF(IPART.EQ.0)THEN
        IGRENZ = IGRENZ0
        IGRENZSQ = IGRENZSQ0
        N3LO = 0
c--NLL running of point-like Higgs coupling to gluons
        CF0G = CF0G0
     .       + CF0G00*(RAT1(ALPHAS(QM,NLOOP),ALPHAS(SCALCG,NLOOP),5)-1)
c       WRITE(6,*)'rat = ',RAT1(ALPHAS(QM,NLOOP),ALPHAS(SCALCG,NLOOP),5)
        CF0X = CF0T+CF0G
        CF0 = CF0S - CF0G0 + CF0G
        IF(ISM4.NE.0) CF0X = CF0T+CF0T4+CF0B4+CF0G
        SIG0=GEVPB*GF/288.D0/DSQRT(2.D0)/PI*CDABS(CF0)**2
        IF(ISILH.NE.0)
     .   SIG0=GEVPB*GF/288.D0/DSQRT(2.D0)/PI
     .                            *(CDABS(CF0)**2-CDABS(CF0P+CF0G)**2)
        DFF = DSIGB/SIGBORN
        SIGBORN=SIG0*DLUM0
        DSIGB=DFF * SIGBORN
        ratnlo = sigborn/sigborn00
       ENDIF
       IF(IPART.EQ.1)THEN
        IGRENZ = 1
        IGRENZSQ = 1
        N3LO = 0
        FACB = 0
        FACC = 0
c--NLL running of point-like Higgs coupling to gluons
        CF0G = CF0G0
     .       + CF0G00*(RAT1(ALPHAS(QM,NLOOP),ALPHAS(SCALCG,NLOOP),5)-1)
        CF0X = CF0T+CF0G
        IF(ISM4.NE.0) CF0X = CF0T+CF0T4+CF0B4+CF0G
        SIG0=GEVPB*GF/288.D0/DSQRT(2.D0)/PI*CDABS(CF0X)**2
        IF(ISILH.NE.0)
     .   SIG0=GEVPB*GF/288.D0/DSQRT(2.D0)/PI*(CDABS(CF0X)**2
     .                                       -CDABS(CF0TP+CF0G)**2)
        DFF = DSIGB/SIGBORN
        SIGBORN=SIG0*DLUM0
        DSIGB=DFF * SIGBORN
        CF0 = CF0X
        CF0B = 0
        CF0C = 0
       ELSEIF(IPART.EQ.2)THEN
        IGRENZ = 1
        IGRENZSQ = 1
        FACB = 0
        FACC = 0
        NLOOP = 3
        N3LO = N3LO0
        ALSMZ = ALSMZ3
        XLAMBDA = XLAMBDA3
        IF(IALPS.EQ.1) THEN
         XLAMBDA=XITLA(NLOOP,ALSMZ,ACC)
         N0=5
        ENDIF
        AMAT=AMT0
        CALL ALSINI(ACC)
        AMAT=AMT0*1.D8
        PDFNAME = PDFNAME3
        NSET = NSET3
        CALL PDFSET(PATHNAME,PDFNAME)
        CALL VEGAS(DLUGG,VAERR,1,IVPNT,IVITM,IVPRN,1)
        DLUM0=TH*VRES*ALPHAS(QM,NLOOP)**2
c--NNLL running of point-like Higgs coupling to gluons
        CF0G = CF0G0
     .       + CF0G00*(RAT2(ALPHAS(QM,NLOOP),ALPHAS(SCALCG,NLOOP),5)-1)
        CF0X = CF0T+CF0G
        if(ism4.ne.0)then
         cf0t4=fact4*cfborn(amh,amt4)
         cf0b4=facb4*cfborn(amh,amb4)
         cf0x = cf0t+cf0t4+cf0b4+cf0g
        endif
c----------------------------------------------------------------
        SIG0=GEVPB*GF/288.D0/DSQRT(2.D0)/PI*CDABS(CF0X)**2
        IF(ISILH.NE.0)
     .   SIG0=GEVPB*GF/288.D0/DSQRT(2.D0)/PI*(CDABS(CF0X)**2
     .                                       -CDABS(CF0TP+CF0G)**2)
        SIGBORN=SIG0*DLUM0
        DSIGB=SIG0*DLUM0/VRES*VERR
        CF0 = CF0X
        CF0B = 0
        CF0C = 0
        ratnnlo = sigborn/sigborn00
       ENDIF
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       IF(ICHVIRT.NE.0)THEN
C--VIRTUAL CORRECTIONS
        IF(IGRENZ.EQ.1)THEN
C--LIMIT OF HEAVY TOP MASS AND SMALL BOTTOM MASS
c        CDLR=CDLOG(-AMH**2/AMB**2/DCMPLX(1.D0,-REPS))
c        CVISO=CF0B*(5.D0/36.D0*CDLR**2-4.D0/3.D0*CDLR)
c        CDLR=CDLOG(-AMH**2/AMC**2/DCMPLX(1.D0,-REPS))
c        CVISO=CVISO+CF0C*(5.D0/36.D0*CDLR**2-4.D0/3.D0*CDLR)
         IF(IHIGGS.NE.1)THEN
C--SCALAR HIGGS
c         XCON=5.D0
c         CVISO=CVISO+(CF0B+CF0C)*XCON
c         CVISO=CF0T*11.D0/2.D0+CVISO
          CVISO=CF0T*11.D0/2.D0
          if(ism4.ne.0)then
           cviso=(cf0t4+cf0b4)*11.d0/2.d0+cviso
          endif
          XCVIRT=DREAL(CVISO/CF0)
          IF(ISILH.NE.0)THEN
           XCVIRT = 11/2.D0*(1 - DREAL(CF0G)/(CDABS(FACT+CF0G)**2
     .                                       -CDABS(FACT-1+CF0G)**2))
          ENDIF
         ELSE
C--PSEUDOSCALAR HIGGS
          XCON=3.8D0
          CVISO=CVISO+(CF0B+CF0C)*XCON
c         CVISO=CF0T*6.D0+CVISO
          CVISO=CF0T*6.D0
          if(ism4.ne.0)then
           cviso=(cf0t4+cf0b4)*6.d0+cviso
          endif
          XCVIRT=DREAL(CVISO/CF0)
         ENDIF
        ELSE
         RHOT=AMH**2/AMQ**2
         RHOB=AMH**2/AMB**2
         RHOC=AMH**2/AMC**2
         RHOST1=AMH**2/AMST1**2
         RHOST2=AMH**2/AMST2**2
         RHOSB1=AMH**2/AMSB1**2
         RHOSB2=AMH**2/AMSB2**2
         IF(IGRENZSQ.NE.0.D0)THEN
          CKSQ = 9.D0
          XCVIRT=DREAL(((CKOF(RHOT)+CTRUN)*CF0T+(CKOF(RHOB)+CBRUN)*CF0B
     .                 +(CKOF(RHOC)+CCRUN)*CF0C
     .                 +CKSQ*(CF0ST1+CF0ST2+CF0SB1+CF0SB2))/CF0)
          IF(ISILH.NE.0)THEN
           RATIO = CDABS(CF0)**2/(CDABS(CF0)**2-CDABS(CF0P+CF0G)**2)
           RAT0  = RATIO*CDABS(CF0P+CF0G)**2/CDABS(CF0)**2
           XCVIRT=RATIO*DREAL(((CKOF(RHOT)+CTRUN)*CF0T
     .                        +(CKOF(RHOB)+CBRUN)*CF0B
     .                        +(CKOF(RHOC)+CCRUN)*CF0C)/CF0)
           IF(RAT0.NE.0.D0) XCVIRT = XCVIRT
     .          -RAT0*DREAL(((CKOF(RHOT)+CTRUN)*(CF0TP)
     .                      +(CKOF(RHOB)+CBRUN)*(CF0BP)
     .                      +(CKOF(RHOC)+CCRUN)*(CF0CP))/(CF0P+CF0G))
          ENDIF
         ELSE
          XCVIRT=DREAL(((CKOF(RHOT)+CTRUN)*CF0T+(CKOF(RHOB)+CBRUN)*CF0B
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c         XCVIRT=DREAL(((CTRUN)*CF0T+(CBRUN)*CF0B
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     .                 +(CKOF(RHOC)+CCRUN)*CF0C
     .                 +CKOFS(RHOST1)*CF0ST1+CKOFS(RHOST2)*CF0ST2
     .                 +CKOFS(RHOSB1)*CF0SB1+CKOFS(RHOSB2)*CF0SB2)/CF0)
          IF(ISILH.NE.0)THEN
           RATIO = CDABS(CF0)**2/(CDABS(CF0)**2-CDABS(CF0P+CF0G)**2)
           RAT0  = RATIO*CDABS(CF0P+CF0G)**2/CDABS(CF0)**2
           XCVIRT=RATIO*DREAL(((CKOF(RHOT)+CTRUN)*CF0T
     .                        +(CKOF(RHOB)+CBRUN)*CF0B
     .                        +(CKOF(RHOC)+CCRUN)*CF0C)/CF0)
           IF(RAT0.NE.0.D0) XCVIRT = XCVIRT
     .          -RAT0*DREAL(((CKOF(RHOT)+CTRUN)*(CF0TP)
     .                      +(CKOF(RHOB)+CBRUN)*(CF0BP)
     .                      +(CKOF(RHOC)+CCRUN)*(CF0CP))/(CF0P+CF0G))
          ENDIF
         ENDIF
         if(ism4.ne.0)then
          rhot4=amh**2/amt4**2
          rhob4=amh**2/amb4**2
          xcvirt=dreal(((ckof(rhot)+ctrun)*cf0t+(ckof(rhob)+cbrun)*cf0b
     .                 +(ckof(rhoc)+ccrun)*cf0c
     .          + ckof(rhot4)*cf0t4+ckof(rhob4)*cf0b4)/cf0)
         endif
        ENDIF
C################  BABIS  ##############################################
        XCSQCD=0
        FACSQCD = 1
        IF(ISQCD.NE.0)THEN
          CSQCD_TOP = DCMPLX(SQCD_TOPR,SQCD_TOPI)*FACSQCD
          CSQCD_BOT = DCMPLX(SQCD_BOTR,SQCD_BOTI)*FACSQCD
          CDSQCD_TOP = DCMPLX(DSQCD_TOPR,DSQCD_TOPI)*FACSQCD
          CDSQCD_BOT = DCMPLX(DSQCD_BOTR,DSQCD_BOTI)*FACSQCD
          XCSQCD = 2*DREAL((CSQCD_TOP*CF0T+CSQCD_BOT*CF0B)/CF0)
c         write(6,*)'SUSYQCD: ',2*CSQCD_TOP*CF0T/CF0*ALPHAS(QM,NLOOP)/PI
c    .                         ,2*CSQCD_BOT*CF0B/CF0*ALPHAS(QM,NLOOP)/PI
c    .                         ,CF0T,CF0B,CF0,ALPHAS(QM,NLOOP)/PI
          CSQCDA_TOP = DCMPLX(SQCDA_TOPR,SQCDA_TOPI)*FACSQCD
          CSQCDA_BOT = DCMPLX(SQCDA_BOTR,SQCDA_BOTI)*FACSQCD
          XCSQCDAT = 2*DREAL((CSQCDA_TOP*CF0T)/CF0)
          XCSQCDAB = 2*DREAL((CSQCDA_BOT*CF0B)/CF0)
          DCSQCD = 2*DSQRT(CDABS(CDSQCD_TOP*CF0T/CF0)**2
     .                    +CDABS(CDSQCD_BOT*CF0B/CF0)**2)
          XCSQCDT = 2*DREAL((CSQCD_TOP*CF0T)/CF0)
          DCSQCDT = 2*DSQRT(CDABS(CDSQCD_TOP*CF0T/CF0)**2)
          XCDMB = 2*DREAL(DGDMB*CF0B/CF0)
        ENDIF
C#######################################################################
        IF(IPROC.EQ.1)THEN
C--GLUON FUSION
         XKVIRT=(XCVIRT+PI**2+(33.D0-2.D0*DNF(QM))/6.D0*SLG0)
     .           *ALPHAS(QM,NLOOP)/PI
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c         CF00T=FACT
c         CF00ST1=FACST1/4
c         CF00ST2=FACST2/4
c         CF00=CF00T+CF00ST1+CF00ST2
c         XCVIRTB=DREAL((11/2.D0*CF00T+9*(CF00ST1+CF00ST2))/CF00)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C################  BABIS  ##############################################
          XCVIRTB=XCVIRT
          XCVIRT0=DREAL(((CKOF(RHOT)+CTRUN)*CF0T/FACT
     .           +(CKOF(RHOB)+CBRUN)*CF0B/FACB
     .           +(CKOF(RHOC)+CCRUN)*CF0C/FACC)
     .           / (CF0T/FACT+CF0B/FACB+CF0C/FACC+CF0G))
          CF00=CF0T+CF0ST1+CF0ST2
          XCVIRTT=DREAL(((CKOF(RHOT)+CTRUN)*CF0T
     .               +CKOFS(RHOST1)*CF0ST1+CKOFS(RHOST2)*CF0ST2)/CF00)
          XCSQCDT0= 2*DREAL((CSQCD_TOP*CF0T)/CF00)
          DCSQCDT0= 2*DSQRT(CDABS(CDSQCD_TOP*CF0T/CF00)**2)
          XCSQCDAT0= 2*DREAL((CSQCDA_TOP*CF0T)/CF00)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         XBABIS=(XCVIRTB+(33.D0-2.D0*DNF(QM))/6.D0*SLG0)
     .           *ALPHAS(QM,NLOOP)/PI
         XBABIST=(XCVIRTT+(33.D0-2.D0*DNF(QM))/6.D0*SLG0)
     .           *ALPHAS(QM,NLOOP)/PI
         XBABISSM=(XCVIRT0+(33.D0-2.D0*DNF(QM))/6.D0*SLG0)
     .           *ALPHAS(QM,NLOOP)/PI
C#######################################################################
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         XKVIRT = XKVIRT*ratnlo + ratnlo-1
         IF(IVERSION.EQ.0)THEN
          IF(IPART.EQ.0) WRITE(NOUT,131)XKVIRT
         ELSE
          SIGVIRT=XKVIRT*SIGBORN00
          DSIGV=XKVIRT*DSIGB00
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          IF(IPART.EQ.0) WRITE(NOUT,5131)SIGVIRT,DABS(DSIGV)
         ENDIF
C################  BABIS  ##############################################
         IF(ISQCD.NE.0)THEN
          XKSQCD=XCSQCD*ALPHAS(QM,NLOOP)/PI
          DKSQCD=DCSQCD*ALPHAS(QM,NLOOP)/PI
          XKSQCDT=XCSQCDT*ALPHAS(QM,NLOOP)/PI
          DKSQCDT=DCSQCDT*ALPHAS(QM,NLOOP)/PI
          XKSQCDT0=XCSQCDT0*ALPHAS(QM,NLOOP)/PI
          DKSQCDT0=DCSQCDT0*ALPHAS(QM,NLOOP)/PI
          XKSQCDAT=XCSQCDAT*ALPHAS(QM,NLOOP)/PI
          XKSQCDAT0=XCSQCDAT0*ALPHAS(QM,NLOOP)/PI
          XKSQCDAB=XCSQCDAB*ALPHAS(QM,NLOOP)/PI
          XKSQCDA=XKSQCDAT+XKSQCDAB
          XKDMB =XCDMB *ALPHAS(QM,NLOOP)/PI
          IF(IVERSION.EQ.0)THEN
           IF(ISQCD.EQ.2) THEN
            XBABIS0= 1 + XBABIS
            XBABIS0T= 1 + XBABIST
            DBABIS0= 0
            XBABIS = 1 + XBABIS + XKSQCD
            DBABIS = DSQRT(DBABIS0**2 + DKSQCD**2)
            XBABIST = 1 + XBABIST + XKSQCDT0
            DBABIST = DSQRT(DBABIS0**2 + DKSQCDT0**2)
            XBABISA = XBABIS0 + XKSQCDA
            XBABISAT = XBABIS0T + XKSQCDAT0
            XBABISAB = XBABIS0 + XKSQCDT + XKSQCDAB
            XMATSM0 = SIG0 * AMH**4*ALPHAS(QM,NLOOP)**2/PI/GEVPB
c    .              * CDABS(CF0T+CF0B+CF0C)**2/CDABS(CF0)**2
            XMATMSSM0 = SIG0 * AMH**4*ALPHAS(QM,NLOOP)**2/PI/GEVPB
            XMATVSM = XMATSM0 * XBABISSM
            XMATV  = XMATMSSM0 * XBABIS
            XMATVA = XMATMSSM0 * XBABISA
           ELSE
            XBABIS = 0
            DBABIS = 0
           ENDIF
           WRITE(NOUT,1319)XKSQCD,DKSQCD
c          WRITE(NOUT,1317)XBABIS,DBABIS
c          WRITE(NOUT,13173)XBABIST,DBABIST
c          WRITE(NOUT,13170)XBABISA
c          WRITE(NOUT,13171)XBABISAT
c          WRITE(NOUT,13172)XBABISAB
c          WRITE(NOUT,1318)XKDMB
           WRITE(6,*)'M_H = ',AMH,' GeV'
           WRITE(6,170)SIGBORN,DABS(DSIGB)
           WRITE(6,1319)XKSQCD,DKSQCD
c          WRITE(6,1317)XBABIS,DBABIS
c          WRITE(6,13173)XBABIST
c          WRITE(6,13170)XBABISA
c          WRITE(6,13171)XBABISAT
c          WRITE(6,13172)XBABISAB
           WRITE(6,*)'M_T         = ',AMQ
           WRITE(6,*)'M_B         = ',AMB
           WRITE(6,*)'K_SM        = ',(1+XBABISSM)*(ALS_NLO/ALS_LO)**2
           WRITE(6,*)'M^2_SM_LO   = ',XMATSM0
     .               * (ALS_LO/ALS_NLO)**2
           WRITE(6,*)'M^2_SM_V    = ',XMATVSM
           WRITE(6,*)'M^2_MSSM_LO = ',XMATMSSM0
     .               * (ALS_LO/ALS_NLO)**2
           WRITE(6,*)'M^2_MSSM_V  = ',XMATV
           WRITE(6,*)'K_MSSM      = ',XBABIS*(ALS_NLO/ALS_LO)**2,
     .               ' +- ',DBABIS*(ALS_NLO/ALS_LO)**2
           WRITE(6,*)'K_MSSM_TOP  = ',XBABIST*(ALS_NLO/ALS_LO)**2,
     .               ' +- ',DBABIST*(ALS_NLO/ALS_LO)**2
           WRITE(6,*)'K_MSSM_APP  = ',XBABISA*(ALS_NLO/ALS_LO)**2
           WRITE(6,*)'K_MSSM_APPT = ',XBABISAT*(ALS_NLO/ALS_LO)**2
           WRITE(6,*)'K_MSSM_APPB = ',XBABISAB*(ALS_NLO/ALS_LO)**2
           WRITE(6,*)'K_NOSUSY    = ',XBABIS0*(ALS_NLO/ALS_LO)**2
c          write(6,*)CSQCD_TOP,CSQCD_BOT
c          write(6,*)CF0,CF0T,CF0B,CF0ST1,CF0ST2,CF0SB1,CF0SB2
c          write(6,*) 2*DREAL((CSQCD_TOP*CF0T)/CF0)*ALPHAS(QM,NLOOP)/PI,
c    .                2*DREAL((CSQCD_BOT*CF0B)/CF0)*ALPHAS(QM,NLOOP)/PI
c          write(6,*) ALPHAS(QM,NLOOP)
c          WRITE(6,1318)XKDMB
          ENDIF
         ENDIF
C#######################################################################
         IF(ICHVIRT.EQ.2)THEN
         PKVT=(PI**2+(33.D0-2.D0*DNF(QM))/6.D0*SLG0)*ALPHAS(QM,NLOOP)/PI
          WRITE(NOUT,200)PKVT
          PKVT=(203.D0/12.D0+33.D0/6.D0*SLG0)*ALPHAS(QM,NLOOP)/PI
          WRITE(NOUT,2001)PKVT
          PKVT=XCVIRT*ALPHAS(QM,NLOOP)/PI
          WRITE(NOUT,201)PKVT
         ENDIF
        ELSE
C--GLUONIC HIGGS DECAYS
         IF(IDIFF.EQ.1)THEN
          IF(IHIGGS.NE.1)THEN
           DUMMY=XCVIRT-5.5D0
          ELSE
           DUMMY=XCVIRT-6.D0
          ENDIF
          DD=0.D0
          WRITE(NOUT,150)DUMMY,DD
         ENDIF
         XKVIRT=XCVIRT+PI**2+(33.D0-2.D0*DNF(QM))/6.D0*SLG0
         IF(IVERSION.EQ.0)THEN
          WRITE(NOUT,141)XKVIRT
         ELSE
          WRITE(NOUT,5141)XKVIRT
         ENDIF
        ENDIF
       ENDIF
C--REAL CORRECTIONS
       IF(ICHGG.NE.0)THEN
        ICALL1=0
        ICALL2=0
        IFAIL=0
        IF66=0
        IF(IPROC.EQ.1)THEN
C--GG --> HG
         IF(INTEG.EQ.1)THEN
          XKGG=DCADR1(DGG,DO,UP,AERR,RERR,ERR,IER)
          DKGG=ERR/DLUM0
         ELSE
          CALL VEGAS(DVGG,VAERR,IVDIM,IVPNT,IVITM,IVPRN,1)
          IER=-1
          XKGG=VRES
          DKGG=VERR/DLUM0
         ENDIF
C--REST-TERMS FROM PLUS-DISTRIBUTIONS
         REST=6.D0*(-DLOG(1.D0-TH)*SLG1+DLOG(1.D0-TH)**2)
         XKGG=XKGG+REST*DLUM0*ALPHAS(QM,NLOOP)/PI
         XKGG=XKGG/DLUM0
         DKGG=XKGG*DSQRT((DKGG/XKGG)**2+(DSIGB/SIGBORN)**2)
C--DELTA-TERMS FROM P_GG SPLITTING-FUNCTION
         DELTA=-(33.D0-2.D0*DNF(QQ))/6.D0*SLG1*ALPHAS(QM,NLOOP)/PI
         XKGG=XKGG+DELTA
         IGG=IER
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         XKGG = XKGG * ratnlo
         DKGG = DKGG * ratnlo
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         IF(IVERSION.EQ.0)THEN
          IF(IPART.EQ.0)
     .     WRITE(NOUT,123)XKGG,DKGG,IGG,ICALL1,ICALL2,IFAIL,IF66
         ELSE
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          SIGGG=XKGG*SIGBORN00
          DSIGGG=SIGGG*DSQRT((DKGG/XKGG)**2-(DSIGB00/SIGBORN00)**2)
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          IF(IPART.EQ.0) WRITE(NOUT,5123)SIGGG,DABS(DSIGGG)
         ENDIF
         IF(ICHGG.EQ.2)THEN
          XTGG=XKGG
          DTGG=DKGG
          DO 301 I=1,3
           ICALL1=0
           IPRTGG=I
           IF(INTEG.EQ.1)THEN
            XIGG=DCADR1(DGG,DO,UP,AERR,RERR,ERR,IER)
            DXIGG=ERR/DLUM0
           ELSE
            CALL VEGAS(DVGG,VAERR,IVDPT,IVPNT,IVITM,IVPRN,1)
            IER=-1
            XIGG=VRES
            DXIGG=VERR/DLUM0
           ENDIF
           XIGG=XIGG/DLUM0
           DXIGG=XIGG*DSQRT((DXIGG/XIGG)**2+(DSIGB/SIGBORN)**2)
C--REST-TERMS FROM PLUS-DISTRIBUTIONS
           IF(I.EQ.1)THEN
             REST=6.D0*(-DLOG(1.D0-TH)*SLG1) + 6.D0*DLOG(1.D0-TH)**2
             XIGG=XIGG+REST*ALPHAS(QM,NLOOP)/PI
C--DELTA-TERMS FROM P_GG SPLITTING-FUNCTION
             DELTA=-(33.D0-2.D0*DNF(QQ))/6.D0*SLG1*ALPHAS(QM,NLOOP)/PI
             XIGG=XIGG+DELTA
           ENDIF
           IGG=IER
           IF(I.EQ.1)THEN 
             WRITE(NOUT,202)XIGG,DXIGG,IGG,ICALL1
           ELSEIF(I.EQ.2)THEN 
             WRITE(NOUT,203)XIGG,DXIGG,IGG,ICALL1
           ELSEIF(I.EQ.3)THEN 
             WRITE(NOUT,2030)XIGG,DXIGG,IGG,ICALL1
           ENDIF
           XTGG=XTGG-XIGG
           DTGG=DSQRT(DTGG**2+DXIGG**2)
301       CONTINUE
          WRITE(NOUT,205)XTGG,DTGG
          IPRTGG=0
         ENDIF
         xrgg = 0
         drgg = 0
         if(nresum1.ne.0.and.igrenz.eq.0)then
          loop0 = 2
          xsig0=gevpb*gf/288.d0/dsqrt(2.d0)/pi*cdabs(cf0t)**2
     .         * alphas(qm,loop0)**2
          massive = 1
c         call soft(mcoll,massive,loop0,xsig0,itopscheme,vaerr,ivpnt,
c    .              ivitm,ivprn,xcorr1,dcorr1)
          massive = 0
c         call soft(mcoll,massive,loop0,xsig0,itopscheme,vaerr,ivpnt,
c    .              ivitm,ivprn,xcorr2,dcorr2)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c         xcorr2 = 0
c         dcorr2 = 0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          xcorr = xcorr1 - xcorr2
          dcorr = dsqrt(dcorr1**2+dcorr2**2)
          xrgg = xcorr * ratnlo/sigborn
          drgg = dcorr * ratnlo/sigborn
c         xkgg = xkgg + xrgg
c         dkgg = dsqrt(dkgg**2 + drgg**2)
          nresum1 = 0
         endif
        ELSE
C-- H --> GGG
         IF(IGRENZ.EQ.1)THEN
          XKGG=0.D0
          DKGG=0.D0
          IGG=0
          ICALL1=0
          ICALL2=0
          IFAIL=0
          IF66=0
         ELSE
          IF(INTEG.EQ.1)THEN
           XKGG=DCADR1(EGG,DO,UP,AERR,RERR,ERR,IER)
           DKGG=ERR
          ELSE
           CALL VEGAS(EVGG,VAERR,IVDIM,IVPNT,IVITM,IVPRN,1)
           IER=-1
           XKGG=VRES
           DKGG=VERR
          ENDIF
          IF(IDIFF.EQ.1)THEN
           WRITE(NOUT,150)XKGG,DKGG
          ENDIF
          IGG=IER
         ENDIF
C--ADD LIMIT OF HEAVY TOP MASS
         XKGG=XKGG+73.D0/4.D0-PI**2
         IF(IVERSION.EQ.0)THEN
          WRITE(NOUT,142)XKGG,DKGG,IGG,ICALL1,ICALL2,IFAIL,IF66
         ELSE
          WRITE(NOUT,5142)XKGG,DABS(DKGG)
         ENDIF
        ENDIF
       ENDIF
 
       IF(ICHGQ.NE.0)THEN
        ICALL1=0
        ICALL2=0
        IFAIL=0
        IF66=0
        IF(IPROC.EQ.1)THEN
C-- GQ --> HQ
         IF(INTEG.EQ.1)THEN
          XKGQ=DCADR1(DGQ,DO,UP,AERR,RERR,ERR,IER)
          DKGQ=ERR
         ELSE
           IF(INTEG.EQ.3)THEN
            VTAU=TH
            VSC=QQ
            CALL VEGAS(DLUGQ,VAERR,1,IVPNT,IVITM,IVPRN,1)
            DLUMQ=VRES
            DDLUMQ=VERR
           ELSE
            DLUMQ=DLUMGQ(TH,QQ)
            DDLUMQ=0.D0
           ENDIF
           ICALL1=0
           ICALL2=0
           CALL VEGAS(DVGQ,VAERR,IVDIM,IVPNT,IVITM,IVPRN,1)
           IER=-1
           XKGQ=VRES
           DKGQ=VERR
         ENDIF
         IF(ISCHEME.EQ.1)THEN
C--REST-TERMS FROM PLUS-DISTRIBUTIONS
          REST=-2.D0/3.D0*(-DLOG(1.D0-TH)**2
     .       +(3.D0-TH-TH**2/2.D0)*DLOG(1.D0-TH)
     .       +9.D0/2.D0*TH+5.D0/4.D0*TH**2)
          XKGQ=XKGQ+DCPL0*REST*DLUMQ
          DKGQ=DSQRT(DKGQ**2+(DCPL0*REST*DDLUMQ)**2)
         ENDIF
         XKGQ=XKGQ/DLUM0
         DKGQ=DKGQ/DLUM0
         DKGQ=XKGQ*DSQRT((DKGQ/XKGQ)**2+(DSIGB/SIGBORN)**2)
         IGQ=IER
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         XKGQ = XKGQ * ratnlo
         DKGQ = DKGQ * ratnlo
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         IF(IVERSION.EQ.0)THEN
          IF(IPART.EQ.0)
     .     WRITE(NOUT,124)XKGQ,DKGQ,IGQ,ICALL1,ICALL2,IFAIL,IF66
         ELSE
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          SIGGQ=XKGQ*SIGBORN00
          DSIGGQ=SIGGQ*DSQRT((DKGQ/XKGQ)**2-(DSIGB00/SIGBORN00)**2)
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          IF(IPART.EQ.0) WRITE(NOUT,5124)SIGGQ,DABS(DSIGGQ)
         ENDIF
         IF(ICHGQ.EQ.2)THEN
          XTGQ=XKGQ
          DTGQ=DKGQ
          DO 302 I=1,2
           ICALL1=0
           IPRTGQ=I
           IF(INTEG.EQ.1)THEN
            XIGQ=DCADR1(DGQ,DO,UP,AERR,RERR,ERR,IER)
            DIGQ=ERR
           ELSE
            CALL VEGAS(DVGQ,VAERR,IVDPT,IVPNT,IVITM,IVPRN,1)
            IER=-1
            XIGQ=VRES
            DIGQ=VERR
           ENDIF
           IF(I.EQ.1.AND.ISCHEME.EQ.1)THEN
C--REST-TERMS FROM PLUS-DISTRIBUTIONS
            REST=-2.D0/3.D0*(-DLOG(1.D0-TH)**2
     .         +(3.D0-TH-TH**2/2.D0)*DLOG(1.D0-TH)
     .         +9.D0/2.D0*TH+5.D0/4.D0*TH**2)
            XIGQ=XIGQ+DCPL0*REST*DLUMQ
            DIGQ=DSQRT(DIGQ**2+(DCPL0*REST*DDLUMQ)**2)
           ENDIF
           XIGQ=XIGQ/DLUM0
           DIGQ=DIGQ/DLUM0
           DIGQ=XIGQ*DSQRT((DIGQ/XIGQ)**2+(DSIGB/SIGBORN)**2)
           IGQ=IER
           IF(I.EQ.1)THEN 
             WRITE(NOUT,206)XIGQ,DIGQ,IGQ,ICALL1
           ELSEIF(I.EQ.2)THEN 
             WRITE(NOUT,207)XIGQ,DIGQ,IGQ,ICALL1
           ENDIF
           XTGQ=XTGQ-XIGQ
           DTGQ=DSQRT(DTGQ**2+DIGQ**2)
302       CONTINUE
          WRITE(NOUT,208)XTGQ,DTGQ
          IPRTGQ=0
         ENDIF
        ELSE
C-- H --> GQQBAR
         IF(IGRENZ.EQ.1)THEN
          XKGQ=0.D0
          DKGQ=0.D0
          IGQ=0
          ICALL1=0
          ICALL2=0
          IFAIL=0
          IF66=0
         ELSE
          IF(INTEG.EQ.1)THEN
           XKGQ=DCADR1(EGQ,DO,UP,AERR,RERR,ERR,IER)
           DKGQ=ERR
          ELSE
           CALL VEGAS(EVGQ,VAERR,IVDIM,IVPNT,IVITM,IVPRN,1)
           IER=-1
           XKGQ=VRES
           DKGQ=VERR
          ENDIF
          IF(IDIFF.EQ.1)THEN
           WRITE(NOUT,150)XKGQ,DKGQ
          ENDIF
          IGQ=IER
         ENDIF
C--ADD LIMIT OF HEAVY TOP MASS
         XKGQ=XKGQ-7.D0/6.D0
         IF(IVERSION.EQ.0)THEN
          WRITE(NOUT,143)XKGQ,DKGQ,IGQ,ICALL1,ICALL2,IFAIL,IF66
         ELSE
          WRITE(NOUT,5143)XKGQ,DABS(DKGQ)
         ENDIF
        ENDIF
       ENDIF
 
       IF(ICHQQ.NE.0)THEN
        ICALL1=0
        ICALL2=0
        IFAIL=0
        IF66=0
        IF(IPROC.EQ.1)THEN
C-- QQBAR --> HG
         IF(INTEG.EQ.1)THEN
          XKQQ=DCADR1(DQQ,DO,UP,AERR,RERR,ERR,IER)
          DKQQ=ERR/DLUM0
         ELSE
          CALL VEGAS(DVQQ,VAERR,IVDIM,IVPNT,IVITM,IVPRN,1)
          IER=-1
          XKQQ=VRES
          DKQQ=VERR/DLUM0
         ENDIF
         XKQQ=XKQQ/DLUM0
         IQQ=IER
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         XKQQ = XKQQ * ratnlo
         DKQQ = DKQQ * ratnlo
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         IF(IVERSION.EQ.0)THEN
          IF(IPART.EQ.0)
     .     WRITE(NOUT,125)XKQQ,DKQQ,IQQ,ICALL1,ICALL2,IFAIL,IF66
         ELSE
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          SIGQQ=XKQQ*SIGBORN00
          DSIGQQ=SIGQQ*DSQRT((DKQQ/XKQQ)**2-(DSIGB00/SIGBORN00)**2)
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          IF(IPART.EQ.0) WRITE(NOUT,5125)SIGQQ,DABS(DSIGQQ)
         ENDIF
        ELSE
C-- H --> GTTBAR
         IF(IGRENZ.EQ.1.OR.NFD.LE.5)THEN
          XKQQ=0.D0
          DKQQ=0.D0
          IQQ=0
          ICALL1=0
          ICALL2=0
          IFAIL=0
          IF66=0
         ELSE
         IF(INTEG.EQ.1)THEN
           XKQQ=DCADR1(EQQ,DO,UP,AERR,RERR,ERR,IER)
           DKQQ=ERR
         ELSE
           IF(AMH.GT.2.D0*AMQ)THEN
            CALL VEGAS(EVQQ,VAERR,IVDIM,IVPNT,IVITM,IVPRN,1)
            XKQQ=VRES
            DKQQ=VERR
           ENDIF
           IER=-1
         ENDIF
          IQQ=IER
         ENDIF
         IF(IVERSION.EQ.0)THEN
          WRITE(NOUT,148)XKQQ,DKQQ,IQQ,ICALL1,ICALL2,IFAIL,IF66
         ELSE
c         WRITE(NOUT,5148)XKQQ,DABS(DKQQ)
         ENDIF
        ENDIF
        if(iversion.eq.0)then
         if(ipart.eq.0)
     .    write(nout,1230)xrgg,drgg,igg,icall1,icall2,ifail,if66
        else
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         if(xrgg.ne.0.d0)then
          siggg=xrgg*sigborn00
          dsiggg=siggg*dsqrt((drgg/xrgg)**2-(dsigb00/sigborn00)**2)
          if(ipart.eq.0) write(nout,5123)siggg,dabs(dsiggg)
         else
          siggg=0
          dsiggg=0
         endif
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        endif
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        IF(NFULL.NE.0)THEN
         IF(IPART.EQ.0)THEN
          SIGBORN0 = SIGBORN
          DSIGB0   = DSIGB
c-------------------------------------------------------------------------
          SIGVIRT0 = XKVIRT*SIGBORN
          DSIGV0   = XKVIRT*DSIGB
          SIGGG0   = XKGG*SIGBORN
          DSIGGG0  = SIGGG0*DSQRT((DKGG/XKGG)**2-(DSIGB/SIGBORN)**2)
          SIGGQ0   = XKGQ*SIGBORN
          DSIGGQ0  = SIGGQ0*DSQRT((DKGQ/XKGQ)**2-(DSIGB/SIGBORN)**2)
          SIGQQ0   = XKQQ*SIGBORN
          DSIGQQ0  = SIGQQ0*DSQRT((DKQQ/XKQQ)**2-(DSIGB/SIGBORN)**2)
          IF(XRGG.NE.0.D0)THEN
           SIGR0    = XRGG*SIGBORN
           DSIGR0   = SIGR0*DSQRT((DRGG/XRGG)**2-(DSIGB/SIGBORN)**2)
          ELSE
           SIGR0    = 0
           DSIGR0   = 0
          ENDIF
c-------------------------------------------------------------------------
c         SIGVIRT0 = XKVIRT*SIGBORN00
c         DSIGV0   = XKVIRT*DSIGB00
c         write(6,*)'0: ',sigvirt0
c         SIGGG0   = XKGG*SIGBORN00
c         DSIGGG0  = SIGGG0*DSQRT((DKGG/XKGG)**2-(DSIGB00/SIGBORN00)**2)
c         SIGGQ0   = XKGQ*SIGBORN00
c         DSIGGQ0  = SIGGQ0*DSQRT((DKGQ/XKGQ)**2-(DSIGB00/SIGBORN00)**2)
c         SIGQQ0   = XKQQ*SIGBORN00
c         DSIGQQ0  = SIGQQ0*DSQRT((DKQQ/XKQQ)**2-(DSIGB00/SIGBORN00)**2)
c         IF(XRGG.NE.0.D0)THEN
c          SIGR0    = XRGG*SIGBORN00
c          DSIGR0   = SIGR0*DSQRT((DRGG/XRGG)**2-(DSIGB00/SIGBORN00)**2)
c         ELSE
c          SIGR0    = 0
c          DSIGR0   = 0
c         ENDIF
c-------------------------------------------------------------------------
         ENDIF
         IF(IPART.EQ.1)THEN
          SIGBORN1 = SIGBORN
          DSIGB1   = DSIGB
c-------------------------------------------------------------------------
          SIGVIRT1 = XKVIRT*SIGBORN
          DSIGV1   = XKVIRT*DSIGB
          SIGGG1   = XKGG*SIGBORN
          DSIGGG1  = SIGGG1*DSQRT((DKGG/XKGG)**2-(DSIGB/SIGBORN)**2)
          SIGGQ1   = XKGQ*SIGBORN
          DSIGGQ1  = SIGGQ1*DSQRT((DKGQ/XKGQ)**2-(DSIGB/SIGBORN)**2)
          SIGQQ1   = XKQQ*SIGBORN
          DSIGQQ1  = SIGQQ1*DSQRT((DKQQ/XKQQ)**2-(DSIGB/SIGBORN)**2)
          IF(XRGG.NE.0.D0)THEN
           SIGR1    = XRGG*SIGBORN
           DSIGR1   = SIGR1*DSQRT((DRGG/XRGG)**2-(DSIGB/SIGBORN)**2)
          ELSE
           SIGR1    = 0
           DSIGR1   = 0
          ENDIF
c-------------------------------------------------------------------------
c         SIGVIRT1 = XKVIRT*SIGBORN00
c         DSIGV1   = XKVIRT*DSIGB00
c         write(6,*)'1: ',sigvirt1,sigvirt1*sigborn0/sigborn00,
c    .                    sigvirt1*sigborn1/sigborn00
c         SIGGG1   = XKGG*SIGBORN00
c         DSIGGG1  = SIGGG1*DSQRT((DKGG/XKGG)**2-(DSIGB00/SIGBORN00)**2)
c         SIGGQ1   = XKGQ*SIGBORN00
c         DSIGGQ1  = SIGGQ1*DSQRT((DKGQ/XKGQ)**2-(DSIGB00/SIGBORN00)**2)
c         SIGQQ1   = XKQQ*SIGBORN00
c         DSIGQQ1  = SIGQQ1*DSQRT((DKQQ/XKQQ)**2-(DSIGB00/SIGBORN00)**2)
c         IF(XRGG.NE.0.D0)THEN
c          SIGR1    = XRGG*SIGBORN00
c          DSIGR1   = SIGR1*DSQRT((DRGG/XRGG)**2-(DSIGB00/SIGBORN00)**2)
c         ELSE
c          SIGR1    = 0
c          DSIGR1   = 0
c         ENDIF
c-------------------------------------------------------------------------
         ENDIF
         IF(IPART.EQ.2)THEN
          SIGBORN2 = SIGBORN
          DSIGB2   = DSIGB
c-------------------------------------------------------------------------
          SIGVIRT2 = XKVIRT*SIGBORN
          DSIGV2   = XKVIRT*DSIGB
          SIGGG2   = XKGG*SIGBORN
          DSIGGG2  = SIGGG2*DSQRT((DKGG/XKGG)**2-(DSIGB/SIGBORN)**2)
          SIGGQ2   = XKGQ*SIGBORN
          DSIGGQ2  = SIGGQ2*DSQRT((DKGQ/XKGQ)**2-(DSIGB/SIGBORN)**2)
          SIGQQ2   = XKQQ*SIGBORN
          DSIGQQ2  = SIGQQ2*DSQRT((DKQQ/XKQQ)**2-(DSIGB/SIGBORN)**2)
          IF(XRGG.NE.0.D0)THEN
           SIGR2    = XRGG*SIGBORN
           DSIGR2   = SIGR2*DSQRT((DRGG/XRGG)**2-(DSIGB/SIGBORN)**2)
          ELSE
           SIGR2    = 0
           DSIGR2   = 0
          ENDIF
c-------------------------------------------------------------------------
c         SIGVIRT2 = XKVIRT*SIGBORN00
c         DSIGV2   = XKVIRT*DSIGB00
c         write(6,*)'2: ',sigvirt2
c         SIGGG2   = XKGG*SIGBORN00
c         DSIGGG2  = SIGGG2*DSQRT((DKGG/XKGG)**2-(DSIGB00/SIGBORN00)**2)
c         SIGGQ2   = XKGQ*SIGBORN00
c         DSIGGQ2  = SIGGQ2*DSQRT((DKGQ/XKGQ)**2-(DSIGB00/SIGBORN00)**2)
c         SIGQQ2   = XKQQ*SIGBORN00
c         DSIGQQ2  = SIGQQ2*DSQRT((DKQQ/XKQQ)**2-(DSIGB00/SIGBORN00)**2)
c         IF(XRGG.NE.0.D0)THEN
c          SIGR2    = XRGG*SIGBORN00
c          DSIGR2   = SIGR2*DSQRT((DRGG/XRGG)**2-(DSIGB00/SIGBORN00)**2)
c         ELSE
c          SIGR2    = 0
c          DSIGR2   = 0
c         ENDIF
c-------------------------------------------------------------------------
          xkggd = xkgg/ratnlo
          xkgqd = xkgq/ratnlo
          xkqqd = xkqq/ratnlo
          xrggd = xrgg/ratnlo
         ENDIF
         IF(ILIMIT.EQ.1.AND.IPART.EQ.0)THEN
          SIGBORN1 = SIGBORN0
          DSIGB1   = DSIGB0
          SIGVIRT1 = SIGVIRT0
          DSIGV1   = DSIGV0
          SIGGG1   = SIGGG0
          DSIGGG1  = DSIGGG0
          SIGGQ1   = SIGGQ0
          DSIGGQ1  = DSIGGQ0
          SIGQQ1   = SIGQQ0
          DSIGQQ1  = DSIGQQ0
          IPART = IPART+2
         ELSE
          IPART = IPART+1
         ENDIF
         IF(IPART.LE.2) GOTO 9991

         SIGVIRT  = SIGVIRT0
         DSIGV    = DSIGVIRT0
         SIGGG    = SIGGG0
         DSIGGG   = DSIGGG0
         SIGGQ    = SIGGQ0
         DSIGGQ   = DSIGGQ0
         SIGQQ    = SIGQQ0
         DSIGQQ   = DSIGQQ0
         XKVIRT   = SIGVIRT0/SIGBORN0
         XKGG     = SIGGG0/SIGBORN0
         DKGG     = DSIGGG0/SIGBORN0
         XKGQ     = SIGGQ0/SIGBORN0
         DKGQ     = DSIGGQ0/SIGBORN0
         XKQQ     = SIGQQ0/SIGBORN0
         DKQQ     = DSIGQQ0/SIGBORN0
         XRGG     = SIGR0/SIGBORN0
         DRGG     = DSIGR0/SIGBORN0

         XCORVIRT = (SIGBORN2 + SIGVIRT2 - SIGBORN1 - SIGVIRT1)/SIGBORN0
         DCORVIRT = DSQRT(DSIGV2**2 + DSIGV1**2)/SIGBORN0
         XCORGG   = (SIGGG2 - SIGGG1)/SIGBORN0
         DCORGG   = DSQRT(DSIGGG2**2 + DSIGGG1**2)/SIGBORN0
         XCORGQ   = (SIGGQ2 - SIGGQ1)/SIGBORN0
         DCORGQ   = DSQRT(DSIGGQ2**2 + DSIGGQ1**2)/SIGBORN0
         XCORQQB  = (SIGQQ2 - SIGQQ1)/SIGBORN0
         DCORQQB  = DSQRT(DSIGQQ2**2 + DSIGQQ1**2)/SIGBORN0
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c        XCORVIRT = 0
c        DCORVIRT = 0
c        XCORGG   = 0
c        DCORGG   = 0
c        XCORGQ   = 0
c        DCORGQ   = 0
c        XCORQQB  = 0
c        DCORQQB  = 0
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

         SIGBORN = SIGBORN0
         DSIGB   = DSIGB0
        ELSE
         XCORVIRT = 0
         DCORVIRT = 0
         XCORGG   = 0
         DCORGG   = 0
         XCORGQ   = 0
         DCORGQ   = 0
         XCORQQB  = 0
         DCORQQB  = 0
         SIGBORN0 = SIGBORN
         DSIGB0   = DSIGB
         SIGBORN1 = SIGBORN
         DSIGB1   = DSIGB
         SIGBORN2 = SIGBORN
         DSIGB2   = DSIGB
        ENDIF
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        ENDIF
C***************************************************************
C         NNLO CORRECTIONS
C***************************************************************

      IGRENZ = 1
      FACB = 0
      FACC = 0

      xk2gg=0.d0
      dk2gg=0.d0
      xk2gq=0.d0
      dk2gq=0.d0
      xk2qq=0.d0
      dk2qq=0.d0
      xk2qqp=0.d0
      dk2qqp=0.d0
      xk2qqb=0.d0
      dk2qqb=0.d0
      if(nloop.eq.3)then
C--NNLO QCD CORRECTIONS
       dlgfac = -slg1
       dlgt = slgt-slg0+slg1
       dlgt0= slgt-slg0
       nf = 5
       integ0=integ
       integ=3
       xk2gg=0.d0
       dk2gg=0.d0
       xk2gq=0.d0
       dk2gq=0.d0
       xk2qq=0.d0
       dk2qq=0.d0
       xk2qqp=0.d0
       dk2qqp=0.d0
       xk2qqb=0.d0
       dk2qqb=0.d0
C--VIRTUAL CORRECTIONS
C--LIMIT OF HEAVY TOP MASS
C--SCALAR HIGGS
      xc2virt=(-(12960*dlgfac**2*zeta2-720*dlgfac*nf*zeta2-1320*dlgfac*
     . nf+11880*dlgfac*zeta2-61560*dlgfac*zeta3+9720*dlgfac-480*dlgt*
     . nf-1710*dlgt+1200*nf*zeta2-600*nf*zeta3+5945*nf+324*zeta2**2-
     . 47880*zeta2+29700*zeta3-56995))/720
c--separation of renormalization and factorization scales
      xc2virt = xc2virt + (11/2.d0+pi**2+(33-2*nf)/6.d0*slg1)
     .                  * (33-2*nf)/4.d0*(slg0-slg1)
     .        + (153-19*nf)/12.d0*(slg0-slg1)
     .        + (33-2*nf)**2/48.d0*(slg0-slg1)**2
c      write(6,*)'point-like: ',xc2virt
       if(ihiggs.ne.1)then
       if(ism4.ne.0)then
        if(fact+fact4+facb4.ne.0.d0)then
         xc2virt = xc2virt -
     .       (77/288.d0*(fact+fact4+facb4-1)
     .     + (2/3.d0*5+19/8.d0)/(fact+fact4+facb4)
     .     * (facb4*dlog(amb4**2/amt0**2)+fact4*dlog(amt4**2/amt0**2)))
c    .            (77/288.d0*2 + (2/3.d0*5+19/8.d0)/3
c    .          * (dlog(amb4**2/amt0**2)+dlog(amt4**2/amt0**2)))
        else
         xc2virt = xc2virt - 77/288.d0*(fact+fact4+facb4-1)
        endif
       endif
       xc2virt0= xc2virt
c--matching of point-like coupling to effective Lagrangian
       del1 = -2*dreal(cf0g*(cf0t+cf0g))*(11*pi**2/4+2777/288.d0
     .                              +19*dlgt/16+nf*(dlgt/3-67/96.d0)
     .                              +11/4.d0*dlgfac*(33-2*nf)/12.d0)
c      del2 = -cf0g**2*(11*pi**2/2+1933/72.d0+19*dlgt/8
c    .                 +nf*(2*dlgt/3-67/48.d0)
c    .                 +11/2.d0*dlgfac*(33-2*nf)/12.d0)
       del2 = cf0g*(2*cf0t+3*cf0g)*(11/4.d0)**2
       xc2virt = xc2virt+(del1+del2)/cdabs(cf0t+cf0g)**2
       if(isilh.ne.0)then
        dc2 = 2777/288.d0+19*dlgt/16+nf*(dlgt/3-67/96.d0)
        dc1 = 11/4.d0
        dd1 = 2*(dc2+dc1*(pi**2+dlgfac*(33-2*nf)/12.d0))
        xc2virt = xc2virt0 - dreal(dd1*cf0g*(cf0t-cf0tp))
     .          / (cdabs(cf0t+cf0g)**2-cdabs(cf0tp+cf0g)**2)
       endif
       endif
       if(ihiggs.eq.1)then
C--PSEUDOSCALAR HIGGS
        xc2virt=xc2virt + nf*(dlgt0/3-21/16.D0)
     .                  + 3*zeta2+1939/144.D0-19*dlgt0/8
     .                  - dlgfac/2 * (33-2*nf)/12.d0
c      write(6,*)'point-like: ',xc2virt
C last term from alpha_s and gluon PDFs
       endif
       if(iproc.eq.1)then
C--GLUON FUSION
        xk2virt=(xc2virt)*(alphas(qm,nloop)/pi)**2
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        xk2virt=xk2virt*ratnnlo
        xk2virt = xk2virt + xcorvirt
        dk2virt = dcorvirt
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if(iversion.eq.0)then
         write(nout,500)xk2virt,dk2virt
        else
         sig2virt=xk2virt*sigborn00
         dsig2virt=sig2virt*dsqrt((dsigb00/sigborn00)**2
     .                           +(dk2virt/xk2virt)**2)
         write(nout,5500)sig2virt,dabs(dsig2virt)
        endif
C--REAL CORRECTIONS
        icall1=0
        icall2=0
        ifail=0
        if66=0
        i2dim=2
        iprtgg=0
C--GG --> HG(G)
        call vegas(d2vgg,vaerr,i2dim,ivpnt,ivitm,ivprn,1)
        ier=-1
        xk2gg=vres
        dk2gg=verr/dlum0
C--REST-TERMS FROM PLUS-DISTRIBUTIONS
        P0 =  DLOG(1.D0-TH)
        P1 = -DLOG(1.D0-TH)**2/2
        P2 =  DLOG(1.D0-TH)**3/3
        P3 = -DLOG(1.D0-TH)**4/4
        rgg0=(18*dlgfac**2*nf-297*dlgfac**2-60*dlgfac*nf-1620*dlgfac*
     .   zeta2+2394*dlgfac-72*nf*zeta2+56*nf+1188*zeta2+6318*zeta3-1212
     .   )/36
        rgg1=(108*dlgfac**2+6*dlgfac*nf-99*dlgfac-10*nf-270*zeta2+399)/
     .   3
        rgg2=2*nf-33+108*dlgfac
        rgg3=72
        rest = rgg3*p3 + rgg2*p2 + rgg1*p1 + rgg0*p0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if(ihiggs.eq.1)rest=rest+6*p1+6*p0*dlgfac/2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        xk2gg=xk2gg+rest*dlum0*(alphas(qm,nloop)/pi)**2
        xk2gg=xk2gg/dlum0
c--separation of renormalization and factorization scales
        xk2gg = xk2gg + xkggd*(33-2*nf)/4.d0*(slg0-slg1)
     .                       *alphas(qm,nloop)/pi
c--matching of point-like gluon coupling
        delta = -11/2.d0*dreal(cf0g*(cf0t+cf0g))/cdabs(cf0t+cf0g)**2
     .        * alphas(qm,nloop)/pi * xkggd
        xk2gg0= xk2gg
        xk2gg = xk2gg + delta
        if(isilh.ne.0)then
         delta = -11/2.d0*dreal(cf0g*(cf0t-cf0tp))
     .         / (cdabs(cf0t+cf0g)**2-cdabs(cf0tp+cf0g)**2)
     .         * alphas(qm,nloop)/pi * xkggd
         xk2gg = xk2gg0 + delta
        endif
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        xk2gg=xk2gg*ratnnlo
        dk2gg=dk2gg*ratnnlo
c       dk2gg=xk2gg*dsqrt((dk2gg/xk2gg)**2-(dsigb/sigborn)**2)
        xk2gg = xk2gg + xcorgg
        dk2gg=dsqrt(dk2gg**2+dcorgg**2)
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        igg=ier
        if(iversion.eq.0)then
         write(nout,502)xk2gg,dk2gg,igg,icall1,icall2,ifail,if66
        else
         sig2gg=xk2gg*sigborn00
         dsig2gg=sig2gg*dsqrt((dk2gg/xk2gg)**2-(dsigb00/sigborn00)**2)
         write(nout,5502)sig2gg,dabs(dsig2gg)
        endif
C--GQ --> HQ(G)
        icall1=0
        icall2=0
        ifail=0
        if66=0
        i2dim=2
        iprtgg=0
        call vegas(d2vgq,vaerr,i2dim,ivpnt,ivitm,ivprn,1)
        ier=-1
        xk2gq=vres/dlum0
        dk2gq=verr/dlum0
c--separation of renormalization and factorization scales
        xk2gq = xk2gq +xkgqd*(33-2*nf)/4.d0*(slg0-slg1)
     .                      *alphas(qm,nloop)/pi
c--matching of point-like gluon coupling
        delta = -11/2.d0*dreal(cf0g*(cf0t+cf0g))/cdabs(cf0t+cf0g)**2
     .        * alphas(qm,nloop)/pi * xkgqd
        xk2gq0= xk2gq
        xk2gq = xk2gq + delta
        if(isilh.ne.0)then
         delta = -11/2.d0*dreal(cf0g*(cf0t-cf0tp))
     .         / (cdabs(cf0t+cf0g)**2-cdabs(cf0tp+cf0g)**2)
     .         * alphas(qm,nloop)/pi * xkgqd
         xk2gq = xk2gq0 + delta
        endif
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        xk2gq=xk2gq*ratnnlo
        dk2gq=dk2gq*ratnnlo
        xk2gq = xk2gq + xcorgq
        dk2gq=dsqrt(dk2gq**2+dcorgq**2)
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        igq=ier
        if(iversion.eq.0)then
         write(nout,5021)xk2gq,dk2gq,igq,icall1,icall2,ifail,if66
        else
         sig2gq=xk2gq*sigborn00
         dsig2gq=sig2gq*dsqrt((dk2gq/xk2gq)**2-(dsigb00/sigborn00)**2)
         write(nout,55021)sig2gq,dabs(dsig2gq)
        endif
C--QQB --> HG(G)
        icall1=0
        icall2=0
        ifail=0
        if66=0
        i2dim=2
        iprtgg=0
        call vegas(d2vqqb,vaerr,i2dim,ivpnt,ivitm,ivprn,1)
        ier=-1
        xk2qqb=vres/dlum0
        dk2qqb=verr/dlum0
c--separation of renormalization and factorization scales
        xk2qqb =xk2qqb+xkqqd*(33-2*nf)/4.d0*(slg0-slg1)
     .                      *alphas(qm,nloop)/pi
c--matching of point-like gluon coupling
        delta = -11/2.d0*dreal(cf0g*(cf0t+cf0g))/cdabs(cf0t+cf0g)**2
     .        * alphas(qm,nloop)/pi * xkqqd
        xk2qqb0= xk2qqb
        xk2qqb = xk2qqb + delta
        if(isilh.ne.0)then
         delta = -11/2.d0*dreal(cf0g*(cf0t-cf0tp))
     .         / (cdabs(cf0t+cf0g)**2-cdabs(cf0tp+cf0g)**2)
     .         * alphas(qm,nloop)/pi * xkqqd
         xk2qqb = xk2qqb0 + delta
        endif
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        xk2qqb=xk2qqb*ratnnlo
        dk2qqb=dk2qqb*ratnnlo
        xk2qqb = xk2qqb + xcorqqb
        dk2qqb =dsqrt(dk2qqb**2+dcorqqb**2)
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        iqqb=ier
        if(iversion.eq.0)then
         write(nout,5022)xk2qqb,dk2qqb,iqqb,icall1,icall2,ifail,if66
        else
         sig2qqb=xk2qqb*sigborn00
         dsig2qqb=sig2qqb*dsqrt((dk2qqb/xk2qqb)**2
     .                         -(dsigb00/sigborn00)**2)
         write(nout,55022)sig2qqb,dabs(dsig2qqb)
        endif
C--QQP --> HQQP
        icall1=0
        icall2=0
        ifail=0
        if66=0
        i2dim=2
        iprtgg=0
        call vegas(d2vqqp,vaerr,i2dim,ivpnt,ivitm,ivprn,1)
        ier=-1
        xk2qqp=vres/dlum0
        dk2qqp=verr/dlum0
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        xk2qqp=xk2qqp*ratnnlo
        dk2qqp=dk2qqp*ratnnlo
        iqqp=ier
        if(iversion.eq.0)then
         write(nout,5024)xk2qqp,dk2qqp,iqqp,icall1,icall2,ifail,if66
        else
         sig2qqp=xk2qqp*sigborn00
         dsig2qqp=sig2qqp*dsqrt((dk2qqp/xk2qqp)**2
     .                         -(dsigb00/sigborn00)**2)
         write(nout,55024)sig2qqp,dabs(dsig2qqp)
        endif
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C--QQ --> HQQ
        icall1=0
        icall2=0
        ifail=0
        if66=0
        i2dim=2
        iprtgg=0
        call vegas(d2vqq,vaerr,i2dim,ivpnt,ivitm,ivprn,1)
        ier=-1
        xk2qq=vres/dlum0
        dk2qq=verr/dlum0
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        xk2qq=xk2qq*ratnnlo
        dk2qq=dk2qq*ratnnlo
        iqqb=ier
        if(iversion.eq.0)then
         write(nout,5023)xk2qq,dk2qq,iqq,icall1,icall2,ifail,if66
        else
         sig2qq=xk2qq*sigborn00
         dsig2qq=sig2qq*dsqrt((dk2qq/xk2qq)**2-(dsigb00/sigborn00)**2)
         write(nout,55023)sig2qq,dabs(dsig2qq)
        endif
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        xkres = 0
        dkres = 0
        if(nresum2.ne.0)then
         massive = 0
         loop0 = 4
         xsig0=gevpb*gf/288.d0/dsqrt(2.d0)/pi*cdabs(cf0t)**2
     .        * alphas(qm,loop0)**2
c        call soft(mcoll,massive,loop0,xsig0,itopscheme,vaerr,ivpnt,
c    .             ivitm,ivprn,xkres,dkres)
c        xkres=xkres*ratnnlo/sigborn
c        dkres=dkres*ratnnlo/sigborn
         xkres=xkres/sigborn00
         dkres=dkres/sigborn00
         nresum2 = 0
         if(iversion.eq.0)then
          write(nout,5020)xkres,dkres,igg,icall1,icall2,ifail,if66
         else
          sigres=xkres*sigborn00
          dsigres=sigres*dsqrt((dkres/xkres)**2-(dsigb00/sigborn00)**2)
          write(nout,55020)sigres,dabs(dsigres)
         endif
        endif

        xk1tot=xkgg+xkgq+xkqq+xkvirt
        dk1tot=dsqrt(dkgg**2+dkgq**2+dkqq**2)
        xk2tot = xk2virt+xk2gg+xk2gq+xk2qqb+xk2qqp+xk2qq
        dk2tot = xk2tot*dsqrt((dk2gg**2 + dk2gq**2 + dk2qqb**2
     .         + dk2qqp**2+dk2qq**2)/xk2tot**2
     .         - (dsigb/sigborn)**2)
        if(iversion.eq.0)then
         write(nout,507)xk1tot,dk1tot
         write(nout,5070)xrgg,drgg
         write(nout,506)xk2tot,dk2tot
         write(nout,5060)xkres,dkres
        endif
       endif
       integ=integ0
      endif
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         WRITE(6,*)
         WRITE(6,*)'Born: '
         WRITE(6,*)SIGBORN00,SIGBORN1,SIGBORN2
c        WRITE(6,*)'rat = ',RAT2(ALPHAS(QM,3),ALPHAS(SCALCG,3),5)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C***************************************************************
C      END OF NNLO CORRECTIONS
C***************************************************************

C--ELW CORRECTIONS
        XKELW = 0
        IF(ISUSY.EQ.0)THEN
         IF(ISM4.EQ.0) THEN
          DMIX = 0
          IF(IELW.EQ.2)THEN
           CW2 = AMW**2/AMZ**2
           SW2 = 1-CW2
           XLAEW = 3*GF*AMW**2/8/DSQRT(2.D0)/PI**2
     .           * (2/CW2*(5/4.D0-7*SW2/3+22*SW2**2/9)+4)
           DMIX = -19*XLAEW/6*ALPHAS(QM,NLOOP)/PI
          ENDIF
          XKELW = GLGL_ELW(AMT0,AMH) + DMIX
         ELSE
          XKELW = GLGL_ELW4(IGGELW,AMH)
c         write(6,*)'elw.: ',AMH,XKELW*100
         ENDIF
        ENDIF
        IF(IELW.EQ.0) XKELW = 0

        IF(I2HDM.NE.0) XKELW = FACELW*XKELW

        IF(IPROC.EQ.1)THEN
         IF(ISUSY.EQ.0.AND.IELW.NE.0)THEN
          IF(IVERSION.EQ.0)THEN
           WRITE(NOUT,1310)XKELW
           WRITE(NOUT,1311)XKELW-DMIX
           WRITE(NOUT,1312)DMIX
          ELSE
           SIGELW = 1+XKELW
           WRITE(NOUT,51310)SIGELW
          ENDIF
         ENDIF
C--TOTAL CORRECTION
         XK1=1.D0+XKGG+XKGQ+XKQQ+XKVIRT
         XK2=1.D0+XKGG+XKGQ+XKQQ+XKVIRT + XK2TOT
         XK1R=XK1+XRGG
         XK2R=XK2+XRGG+XKRES
         XKTOT=1.D0+XKGG+XKGQ+XKQQ+XKVIRT + XK2TOT + XRGG + XKRES
         XKTOT=XKTOT*(1+XKELW)
C--ERROR
         DK1=DSQRT(DKGG**2+DKGQ**2+DKQQ**2)
         DK2=DSQRT(DKGG**2+DKGQ**2+DKQQ**2 + DK2TOT**2)
         DK1R=DSQRT(DK1**2+DRGG**2)
         DK2R=DSQRT(DK2**2+DRGG**2+DKRES**2)
         DKTOT=DSQRT(DKGG**2+DKGQ**2+DKQQ**2+DK2TOT**2
     .              +DRGG**2+DKRES**2)*(1+XKELW)
C--TOTAL CROSS SECTION
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c        FACN3LO = 1
c        IF(N3LO.NE.0)THEN
c         XK2 = 1 + (XK2-1)*ALPHAS(QM,4)/ALPHAS(QM,3)
c         XKTOT = XK2*(1+XKELW)
c         DK2 = DK2*ALPHAS(QM,4)/ALPHAS(QM,3)
c         DKTOT = DKTOT*ALPHAS(QM,4)/ALPHAS(QM,3)
c         FACN3LO = (ALPHAS(QM,4)/ALPHAS(QM,3))**2
c        ENDIF
c        SIG1=SIGBORN00*XK1
c        SIG2=SIGBORN00*XK2 * FACN3LO
c        SIGTOT=SIGBORN00*XKTOT * FACN3LO
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         SIG1=SIGBORN00*XK1
         SIG2=SIGBORN00*XK2
         SIG1R=SIGBORN00*XK1R
         SIG2R=SIGBORN00*XK2R
         SIGTOT=SIGBORN00*XKTOT
C--ERROR
         DSIG1=DSQRT(DSIGB00**2+SIGBORN00**2*DK1**2)
         DSIG2=DSQRT(DSIGB00**2+SIGBORN00**2*DK2**2)
         DSIG1R=DSQRT(DSIGB00**2+SIGBORN00**2*DK1R**2)
         DSIG2R=DSQRT(DSIGB00**2+SIGBORN00**2*DK2R**2)
         DSIGTOT=DSQRT(DSIGB00**2+SIGBORN00**2*DKTOT**2)
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         IF(IVERSION.EQ.0)THEN
          WRITE(NOUT,1270)XKTOT/(1+XKELW),DKTOT/(1+XKELW)
          WRITE(NOUT,127)XKTOT,DKTOT
          WRITE(6,127)XKTOT,DKTOT
          WRITE(NOUT,1700)SIGBORN00,DABS(DSIGB00)
          WRITE(NOUT,1281)SIG1,DSIG1
          IF(NRESUM0.NE.0)WRITE(NOUT,12810)SIG1R,DSIG1R
          IF(NLOOP00.GT.2)WRITE(NOUT,1282)SIG2,DSIG2
          IF(NLOOP00.GT.2.AND.NRESUM0.NE.0)WRITE(NOUT,12820)SIG2R,DSIG2R
          WRITE(NOUT,128)SIGTOT,DSIGTOT
         ELSE
          IF(NLOOP00.EQ.2)THEN
           WRITE(NOUT,5127)SIGBORN00,DSIGB00
           WRITE(NOUT,5128)SIGTOT,DABS(DSIGTOT)
          ELSE
           WRITE(NOUT,5127)SIGBORN00,DSIGB00
           WRITE(NOUT,5128)SIG1,DSIG1
           IF(NLOOP00.GT.2) WRITE(NOUT,5129)SIG2,DSIG2
           WRITE(NOUT,5130)SIGTOT,DABS(DSIGTOT)
          ENDIF
         ENDIF
        ELSE
C--TOTAL COEFFICIENT
         XKTOT=XKGG+XKGQ*NFD+XKQQ+XKVIRT
C--ERROR
         DKTOT=DSQRT(DKGG**2+DKGQ**2*NFD**2+DKQQ**2)
         IF(IDIFF.EQ.1)THEN
          IF(IHIGGS.EQ.1)THEN
           DUMMY=XKTOT-97.D0/4.D0+7.D0/6.D0*NFD
          ELSE
           DUMMY=XKTOT-95.D0/4.D0+7.D0/6.D0*NFD
          ENDIF
          WRITE(NOUT,150)DUMMY,DKTOT
         ENDIF
         IF(IVERSION.EQ.0)THEN
          WRITE(NOUT,144)XKTOT,DKTOT
         ELSE
          WRITE(NOUT,5144)XKTOT,DABS(DKTOT)
         ENDIF
C--ELW CORRECTIONS
         IF(ISUSY.EQ.0)THEN
          IF(IVERSION.EQ.0)THEN
           WRITE(NOUT,1410)XKELW
          ELSE
           WRITE(NOUT,51410)XKELW
          ENDIF
         ENDIF
C--TOTAL CORRECTION
         XKTOT=1.D0+XKTOT*ALPHAS(QM,NLOOP)/PI
C--ERROR
         DKTOT=DKTOT*ALPHAS(QM,NLOOP)/PI
C--TOTAL DECAY WIDTH
         GAMTOT=GAM0*XKTOT*(1+XKELW)
C--ERROR
         DGAMTOT=GAM0*DKTOT
         IF(IVERSION.EQ.0)THEN
          WRITE(NOUT,147)XKTOT,DKTOT
          WRITE(NOUT,1470)XKTOT*(1+XKELW),DKTOT*(1+XKELW)
          WRITE(NOUT,145)GAMTOT,DGAMTOT
         ELSE
          WRITE(NOUT,5145)GAMTOT,DABS(DGAMTOT)
         ENDIF
       ENDIF
      ENDIF
      ELSEIF(IPTY.EQ.1)THEN
C*****************************************
C--SIGMA(PTMAX > PT > PTMIN)
C*****************************************

C--CHECK FOR PT
      write(6,*)'SIGMA (PTMAX > PT > PTMIN)'
      ETA = -9.999D30
      PTMAX0 = DSQRT((S+AMH**2)**2 - 4*AMH**2*S)/2.D0/DSQRT(S)
      IF(PTMIN.LT.0.D0)THEN
       WRITE(6,*)'PTMIN OUT OF RANGE!! TAKING PTMIN = 0...'
       PTMIN = 0
      ENDIF
      IF(PTMAX.GT.PTMAX0)THEN
       WRITE(6,*)'PTMAX OUT OF RANGE!! TAKING PTMAX = ',PTMAX0,' GEV...'
       PTMAX = PTMAX0
      ENDIF
      WRITE(NOUT,3210)PTMIN,PTMAX

      IF(ICHANNEL.NE.0)THEN
C--REAL CORRECTIONS
       IF(ICHGG.NE.0)THEN
        ICALL1=0
        ICALL2=0
        IFAIL=0
        IF66=0
        IVDIM = 3
        IF(IPROC.EQ.1)THEN
C--GG --> HG
         CALL VEGAS(TPTGG,VAERR,IVDIM,IVPNT,IVITM,IVPRN,1)
         IER=-1
         XKGG=VRES*SIG0
         DKGG=VERR*SIG0
         IGG=IER
         WRITE(NOUT,3230)XKGG,DKGG
        ENDIF
       ENDIF
 
       IF(ICHGQ.NE.0)THEN
        ICALL1=0
        ICALL2=0
        IFAIL=0
        IF66=0
        IF(IPROC.EQ.1)THEN
C-- GQ --> HQ
         ICALL1=0
         ICALL2=0
         CALL VEGAS(TPTGQ,VAERR,IVDIM,IVPNT,IVITM,IVPRN,1)
         IER=-1
         XKGQ=VRES*SIG0
         DKGQ=VERR*SIG0
         IGQ=IER
         WRITE(NOUT,3240)XKGQ,DKGQ
        ENDIF
       ENDIF
 
       IF(ICHQQ.NE.0)THEN
        ICALL1=0
        ICALL2=0
        IFAIL=0
        IF66=0
        IF(IPROC.EQ.1)THEN
C-- QQBAR --> HG
         CALL VEGAS(TPTQQ,VAERR,IVDIM,IVPNT,IVITM,IVPRN,1)
         IER=-1
         XKQQ=VRES
         XKQQ=VRES*SIG0
         DKQQ=VERR*SIG0
         IQQ=IER
         WRITE(NOUT,3250)XKQQ,DKQQ
C--TOTAL CORRECTION
         XKTOT=XKGG+XKGQ+XKQQ
C--ERROR
         DKTOT=DSQRT(DKGG**2+DKGQ**2+DKQQ**2)
         WRITE(NOUT,3280)XKTOT,DKTOT
        ENDIF
       ENDIF
      ENDIF
      ELSEIF(IPTY.EQ.2)THEN
C*****************************************
C--DSIGMA/DPT
C*****************************************

C     DO 9990 XPT = PTMIN,PTMAX,PTSTEP
      IMAX = (PTMAX-PTMIN)/PTSTEP
      DO 9990 IJK = 0,IMAX
       XPT = PTMIN + PTSTEP*IJK
       PT = XPT

C--CHECK FOR PT
      write(6,*)'DSIGMA/DPT'
      ETA = -9.999D30
      PTMAX = DSQRT((S+AMH**2)**2 - 4*AMH**2*S)/2.D0/DSQRT(S)
      IF(PT.GE.PTMAX.OR.PT.LE.0.D0)THEN
       WRITE(6,*)'PT OUT OF RANGE!!!'
       STOP
      ENDIF
      WRITE(NOUT,321)PT
      AMHT = DSQRT(AMH**2+PT**2)
      DUMR = XKAPM*AMHT
      WRITE(6,*)'Scale = ',DUMR,' GeV'

      IF(ICHANNEL.NE.0)THEN
C--REAL CORRECTIONS
       IF(ICHGG.NE.0)THEN
        ICALL1=0
        ICALL2=0
        IFAIL=0
        IF66=0
        IVDIM = 2
        IF(IPROC.EQ.1)THEN
C--GG --> HG
         CALL VEGAS(FPTGG,VAERR,IVDIM,IVPNT,IVITM,IVPRN,1)
         IER=-1
         XKGG=VRES*SIG0
         DKGG=VERR*SIG0
         IGG=IER
         WRITE(NOUT,323)XKGG,DKGG
        ENDIF
       ENDIF
 
       IF(ICHGQ.NE.0)THEN
        ICALL1=0
        ICALL2=0
        IFAIL=0
        IF66=0
        IF(IPROC.EQ.1)THEN
C-- GQ --> HQ
         ICALL1=0
         ICALL2=0
         CALL VEGAS(FPTGQ,VAERR,IVDIM,IVPNT,IVITM,IVPRN,1)
         IER=-1
         XKGQ=VRES*SIG0
         DKGQ=VERR*SIG0
         IGQ=IER
         WRITE(NOUT,324)XKGQ,DKGQ
        ENDIF
       ENDIF
 
       IF(ICHQQ.NE.0)THEN
        ICALL1=0
        ICALL2=0
        IFAIL=0
        IF66=0
        IF(IPROC.EQ.1)THEN
C-- QQBAR --> HG
         CALL VEGAS(FPTQQ,VAERR,IVDIM,IVPNT,IVITM,IVPRN,1)
         IER=-1
         XKQQ=VRES
         XKQQ=VRES*SIG0
         DKQQ=VERR*SIG0
         IQQ=IER
         WRITE(NOUT,325)XKQQ,DKQQ
C--TOTAL CORRECTION
         XKTOT=XKGG+XKGQ+XKQQ
C--ERROR
         DKTOT=DSQRT(DKGG**2+DKGQ**2+DKQQ**2)
         WRITE(NOUT,328)XKTOT,DKTOT
        ENDIF
       ENDIF
      ENDIF
9990  CONTINUE
      ELSE
C*****************************************
C--D2SIGMA/DPT/DY
C*****************************************

      write(6,*)'D2SIGMA/DPT/DY'
C--CHECK FOR PT,ETA
      XCSH = (DEXP(ETA)+DEXP(-ETA))/2.D0
      PTMAX = DSQRT((S+AMH**2)**2 - 4*AMH**2*S*XCSH**2)
     .               /2.D0/DSQRT(S)/XCSH
      Y0 = (S+AMH**2)/2.D0/DSQRT(S*(AMH**2+PT**2))
      YMAX = DLOG(Y0 + DSQRT(Y0**2 -1.D0))
      IF(PT.GE.PTMAX.OR.PT.LE.0.D0)THEN
       WRITE(6,*)'PT OUT OF RANGE!!!'
       write(6,*)pt,ptmax,eta,ymax
       STOP
      ENDIF
      IF(DABS(ETA).GE.YMAX)THEN
       WRITE(6,*)'Y OUT OF RANGE!!!'
       STOP
      ENDIF
      WRITE(NOUT,320)PT,ETA

      IF(ICHANNEL.NE.0)THEN
C--REAL CORRECTIONS
       IF(ICHGG.NE.0)THEN
        ICALL1=0
        ICALL2=0
        IFAIL=0
        IF66=0
        IVDIM = 1
        IF(IPROC.EQ.1)THEN
C--GG --> HG
         CALL VEGAS(FPTYGG,VAERR,IVDIM,IVPNT,IVITM,IVPRN,1)
         IER=-1
         XKGG=VRES*SIG0
         DKGG=VERR*SIG0
         IGG=IER
         WRITE(NOUT,323)XKGG,DKGG
        ENDIF
       ENDIF
 
       IF(ICHGQ.NE.0)THEN
        ICALL1=0
        ICALL2=0
        IFAIL=0
        IF66=0
        IF(IPROC.EQ.1)THEN
C-- GQ --> HQ
         ICALL1=0
         ICALL2=0
         CALL VEGAS(FPTYGQ,VAERR,IVDIM,IVPNT,IVITM,IVPRN,1)
         IER=-1
         XKGQ=VRES*SIG0
         DKGQ=VERR*SIG0
         IGQ=IER
         WRITE(NOUT,324)XKGQ,DKGQ
        ENDIF
       ENDIF
 
       IF(ICHQQ.NE.0)THEN
        ICALL1=0
        ICALL2=0
        IFAIL=0
        IF66=0
        IF(IPROC.EQ.1)THEN
C-- QQBAR --> HG
         CALL VEGAS(FPTYQQ,VAERR,IVDIM,IVPNT,IVITM,IVPRN,1)
         IER=-1
         XKQQ=VRES
         XKQQ=VRES*SIG0
         DKQQ=VERR*SIG0
         IQQ=IER
         WRITE(NOUT,325)XKQQ,DKQQ
C--TOTAL CORRECTION
         XKTOT=XKGG+XKGQ+XKQQ
C--ERROR
         DKTOT=DSQRT(DKGG**2+DKGQ**2+DKQQ**2)
         WRITE(NOUT,328)XKTOT,DKTOT
        ENDIF
       ENDIF
      ENDIF
      ENDIF
 
      WRITE(NOUT,*)
 
100   FORMAT(10X,I20)
101   FORMAT(10X,G20.10)
102   FORMAT(11X,A5)
103   FORMAT(11X,A100)
110   FORMAT(1X,'MH',10X,'SIGBORN',5X,'KGG',9X,'KGQ',9X,'KQQ')
111   FORMAT(1X,'RHO',9X,'SIGBORN',5X,'KGG',9X,'KGQ',9X,'KQQ')
112   FORMAT(1X,60('='))
113   FORMAT(1X,'EPSTAU = ',G15.6,5X,'EPSV   = ',G15.6)
114   FORMAT(1X,'ABSERR = ',G15.6,5X,'RELERR = ',G15.6)
115   FORMAT(1X,'LIMIT  = ',I5,15X,'SCALE  =  RHO * MH')
180   FORMAT(9X,'STRUCTURE-FUNCTIONS: SCALE  =  RHO * RHOFAC * MH')
181   FORMAT(1X,'RHOFAC = ',G15.6)
117   FORMAT(1X,'DLTALS = ',G15.6,5X,'DLTNF  = ',G15.6)
118   FORMAT(1X,'ENERGY = ',G15.6,3X,'LAMBDA_',I1,' = ',G15.6)
11801 FORMAT(1X,'ENERGY       = ',G15.6,5X,'LO-LAMBDA_',I1,' = ',G15.6)
11802 FORMAT(1X,'NLO-LAMBDA_',I1,' = ',G15.6,3X,
     .         'NNLO-LAMBDA_',I1,' = ',G15.6)
119   FORMAT(1X,'MB     = ',G15.6,5X,'MT     = ',G15.6,1X,
     .          'MC     = ',G15.6)
1191  FORMAT(1X,'MBP    = ',G15.6,5X,'MTP    = ',G15.6)
C120  FORMAT(1X,'STFUN  = ',I2,3X,I2,13X,'SCHEME = ',I1)
120   FORMAT(1X,'STFUN  = ',A40,'SET = ',I2,5X,'SCHEME = ',I1)
1205  FORMAT(1X,'  LO-STFUN  = ',A40,'SET = ',I2,5X,'SCHEME = ',I1)
1206  FORMAT(1X,' NLO-STFUN  = ',A40,'SET = ',I2,5X,'SCHEME = ',I1)
1207  FORMAT(1X,'NNLO-STFUN  = ',A40,'SET = ',I2,5X,'SCHEME = ',I1)
1201  FORMAT(1X,'STFUN  = GRV(HO)',13X,'SCHEME = ',I1)
1202  FORMAT(1X,'STFUN  = CTEQ6  ',13X,'SCHEME = ',I1)
121   FORMAT(1X,'MH     = ',G15.6,5X,'RHOM   = ',G15.6,1X,
     .          'RHOQ   = ',G15.6)
1210  FORMAT(1X,'MH     = ',G15.6,5X,'RHOM   = ',G15.6)
122   FORMAT(1X,'SIGBORN  =  ',G15.6,2X,'PB')
170   FORMAT(1X,'SIGBORN  =  ',G15.6,1X,'+- ',G15.6,2X,'PB')
1700  FORMAT(1X,'SIGLO  = ',G15.6,1X,'+- ',G15.6,1X,'PB')
123   FORMAT(1X,'KGG    = ',G15.6,1X,'+- ',G15.6,2X,I3,
     .       2X,I5,2X,I7,2X,I5,2X,I5)
1230  FORMAT(1X,'KMRES  = ',G15.6,1X,'+- ',G15.6,2X,I3,
     .       2X,I5,2X,I7,2X,I5,2X,I5)
124   FORMAT(1X,'KGQ    = ',G15.6,1X,'+- ',G15.6,2X,I3,
     .       2X,I5,2X,I7,2X,I5,2X,I5)
125   FORMAT(1X,'KQQB   = ',G15.6,1X,'+- ',G15.6,2X,I3,
     .       2X,I5,2X,I7,2X,I5,2X,I5)
126   FORMAT(1X,'RHOIM  = ',G15.6,5X,'LOOP   = ',I1)
1270  FORMAT(1X,'KQCD   = ',G15.6,1X,'+- ',G15.6)
127   FORMAT(1X,'KTOT   = ',G15.6,1X,'+- ',G15.6)
128   FORMAT(1X,'SIGTOT = ',G15.6,1X,'+- ',G15.6,1X,'PB')
1281  FORMAT(1X,'SIGNLO = ',G15.6,1X,'+- ',G15.6,1X,'PB')
1282  FORMAT(1X,'SIGNNLO= ',G15.6,1X,'+- ',G15.6,1X,'PB')
12810 FORMAT(1X,'SIG1RES= ',G15.6,1X,'+- ',G15.6,1X,'PB')
12820 FORMAT(1X,'SIG2RES= ',G15.6,1X,'+- ',G15.6,1X,'PB')
129   FORMAT(1X,'COUPB  = ',G15.6,5X,'COUPT  = ',G15.6,2X,
     .      'TG(BETA) = ',G9.3)
1295  FORMAT(1X,'COUPC  = ',G15.6,5X,'COUPG  = ',G15.6)
1296  FORMAT(1X,'COUPCTG= ',G15.6)
1290  FORMAT(1X,'MSTOP1 = ',G15.6,5X,'MSTOP2 = ',G15.6)
1291  FORMAT(1X,'COUPST1= ',G15.6,5X,'COUPST2= ',G15.6)
1292  FORMAT(1X,'MSBOT1 = ',G15.6,5X,'MSBOT2 = ',G15.6)
1293  FORMAT(1X,'COUPSB1= ',G15.6,5X,'COUPSB2= ',G15.6)
1294  FORMAT(1X,'MA     = ',G15.6)
130   FORMAT(1X,'HIGGS  = ',A1)
131   FORMAT(1X,'KVIRT  = ',G15.6)
1319  FORMAT(1X,'KSQCD  = ',G15.6,1X,'+- ',G15.6)
1317  FORMAT(1X,'BABIS  = ',G15.6,1X,'+- ',G15.6)
13170 FORMAT(1X,'APP    = ',G15.6)
13171 FORMAT(1X,'TOP_APP= ',G15.6)
13172 FORMAT(1X,'BOT_APP= ',G15.6)
13173 FORMAT(1X,'TOP    = ',G15.6)
1318  FORMAT(1X,'K_DMB  = ',G15.6)
1310  FORMAT(1X,'K_ELW  = ',G15.6)
1311  FORMAT(1X,'K_ELW1 = ',G15.6)
1312  FORMAT(1X,'K_ELW2 = ',G15.6)
132   FORMAT(1X,'COUPB  = ',G15.6,5X,'COUPT  = ',G15.6,1X,
     .          'COUPC  = ',G15.6)
1320  FORMAT(1X,'COUPG  = ',G15.6,5X,'COUPCTG= ',G15.6)
1321  FORMAT(1X,'COUPBP = ',G15.6,5X,'COUPTP = ',G15.6)
140   FORMAT(1X,'GAMBORN  =  ',G15.6,2X,'MEV')
141   FORMAT(1X,'CVIRT  = ',G15.6)
1410  FORMAT(1X,'C_ELW  = ',G15.6)
142   FORMAT(1X,'CGGG   = ',G15.6,1X,'+- ',G15.6,2X,I3,
     .       2X,I5,2X,I7,2X,I5,2X,I5)
143   FORMAT(1X,'CGQQ   = ',G15.6,1X,'+- ',G15.6,2X,I3,
     .       2X,I5,2X,I7,2X,I5,2X,I5)
148   FORMAT(1X,'CGTT   = ',G15.6,1X,'+- ',G15.6,2X,I3,
     .       2X,I5,2X,I7,2X,I5,2X,I5)
144   FORMAT(1X,'CTOT   = ',G15.6,1X,'+- ',G15.6)
145   FORMAT(1X,'GAMTOT = ',G15.6,1X,'+- ',G15.6,1X,'MEV')
146   FORMAT(1X,'NF     = ',I15,3X,'LAMBDA_',I1,' = ',G15.6)
1460  FORMAT(1X,'ALPHA_S(M_Z) = ',G15.6)
1461  FORMAT(1X,'  LO-ALPHA_S(M_Z) = ',G15.6)
1462  FORMAT(1X,' NLO-ALPHA_S(M_Z) = ',G15.6)
1463  FORMAT(1X,'NNLO-ALPHA_S(M_Z) = ',G15.6)
147   FORMAT(1X,'KQCD   = ',G15.6,1X,'+- ',G15.6)
1470  FORMAT(1X,'KTOT   = ',G15.6,1X,'+- ',G15.6)
150   FORMAT(1X,'DIFF   = ',G15.6,1X,'+- ',G15.6)
160   FORMAT(1X,'=======')
161   FORMAT(1X,'DCADRE:')
162   FORMAT(1X,'VEGAS: (ERROR-CODE = -1)')
165   FORMAT(1X,'POINTS = ',I10,10X,'ITERATIONS = ',I5)
200   FORMAT(1X,'  PI2  = ',G15.6)
2001  FORMAT(1X,'DALPHA = ',G15.6)
201   FORMAT(1X,'  C0   = ',G15.6)
202   FORMAT(1X,'ALP/BET= ',G15.6,1X,'+- ',G15.6,2X,I3,2X,I5)
203   FORMAT(1X,'DGAMMA = ',G15.6,1X,'+- ',G15.6,2X,I3,2X,I5)
2030  FORMAT(1X,'DDELTA = ',G15.6,1X,'+- ',G15.6,2X,I3,2X,I5)
205   FORMAT(1X,'  REST = ',G15.6,1X,'+- ',G15.6)
206   FORMAT(1X,'  PGQS = ',G15.6,1X,'+- ',G15.6,2X,I3,2X,I5)
207   FORMAT(1X,'  PGQL = ',G15.6,1X,'+- ',G15.6,2X,I3,2X,I5)
208   FORMAT(1X,'  D0GQ = ',G15.6,1X,'+- ',G15.6)

500   FORMAT(1X,'K2VIRT = ',G15.6,1X,'+- ',G15.6)
501   FORMAT(1X,' DALPHA = ',G15.6)
5020  FORMAT(1X,'KGGRES = ',G15.6,1X,'+- ',G15.6,2X,I3,
     .       2X,I5,2X,I7,2X,I5,2X,I5)
502   FORMAT(1X,'K2GG   = ',G15.6,1X,'+- ',G15.6,2X,I3,
     .       2X,I5,2X,I7,2X,I5,2X,I5)
5021  FORMAT(1X,'K2GQ   = ',G15.6,1X,'+- ',G15.6,2X,I3,
     .       2X,I5,2X,I7,2X,I5,2X,I5)
5022  FORMAT(1X,'K2QQB  = ',G15.6,1X,'+- ',G15.6,2X,I3,
     .       2X,I5,2X,I7,2X,I5,2X,I5)
5023  FORMAT(1X,'K2QQ   = ',G15.6,1X,'+- ',G15.6,2X,I3,
     .       2X,I5,2X,I7,2X,I5,2X,I5)
5024  FORMAT(1X,'K2QQP  = ',G15.6,1X,'+- ',G15.6,2X,I3,
     .       2X,I5,2X,I7,2X,I5,2X,I5)
5500  FORMAT(1X,'SIG2_VIRT = ',G15.6,1X,'+- ',G15.6,2X,'PB')
55020 FORMAT(1X,'SIG_GGRES = ',G15.6,1X,'+- ',G15.6,2X,'PB')
5502  FORMAT(1X,'SIG2_GG   = ',G15.6,1X,'+- ',G15.6,2X,'PB')
55021 FORMAT(1X,'SIG2_GQ   = ',G15.6,1X,'+- ',G15.6,2X,'PB')
55022 FORMAT(1X,'SIG2_QQB  = ',G15.6,1X,'+- ',G15.6,2X,'PB')
55023 FORMAT(1X,'SIG2_QQ   = ',G15.6,1X,'+- ',G15.6,2X,'PB')
55024 FORMAT(1X,'SIG2_QQP  = ',G15.6,1X,'+- ',G15.6,2X,'PB')
503   FORMAT(1X,'  BETA  = ',G15.6,1X,'+- ',G15.6,2X,I3,2X,I5)
504   FORMAT(1X,' DALPHA = ',G15.6,1X,'+- ',G15.6,2X,I3,2X,I5)
505   FORMAT(1X,' DDELTA = ',G15.6,1X,'+- ',G15.6,2X,I3,2X,I5)
506   FORMAT(1X,'K(2)   = ',G15.6,1X,'+- ',G15.6)
507   FORMAT(1X,'K(1)   = ',G15.6,1X,'+- ',G15.6)
5060  FORMAT(1X,'K(2)RES= ',G15.6,1X,'+- ',G15.6)
5070  FORMAT(1X,'K(1)RES= ',G15.6,1X,'+- ',G15.6)

320   FORMAT(1X,'PT     = ',G15.6,' GEV     Y  = ',G15.6)
321   FORMAT(1X,'PT     = ',G15.6,' GEV')
3210  FORMAT(1X,'PTMIN  = ',G15.6,' GEV',3X,'PTMAX  = ',G15.6,' GEV')
323   FORMAT(1X,'SGG    = ',G15.6,1X,'+- ',G15.6,1X,'PB/GEV')
324   FORMAT(1X,'SGQ    = ',G15.6,1X,'+- ',G15.6,1X,'PB/GEV')
325   FORMAT(1X,'SQQB   = ',G15.6,1X,'+- ',G15.6,1X,'PB/GEV')
328   FORMAT(1X,'SIGTOT = ',G15.6,1X,'+- ',G15.6,1X,'PB/GEV')
3230  FORMAT(1X,'SGG    = ',G15.6,1X,'+- ',G15.6,1X,'PB')
3240  FORMAT(1X,'SGQ    = ',G15.6,1X,'+- ',G15.6,1X,'PB')
3250  FORMAT(1X,'SQQB   = ',G15.6,1X,'+- ',G15.6,1X,'PB')
3280  FORMAT(1X,'SIGTOT = ',G15.6,1X,'+- ',G15.6,1X,'PB')
 
5162  FORMAT(1X,'VEGAS:')
5160  FORMAT(1X,'======')
5165  FORMAT(1X,'POINTS    = ',I10,10X,'ITERATIONS   = ',I5)
5166  FORMAT(1X,'ABSERR    = ',G15.6)
1188  FORMAT(1X,'SILH LAGRANGIAN')
1189  FORMAT(1X,'NON-LINEAR LAGRANGIAN')
1184  FORMAT(1X,'GLUON FUSION: GG --> HIGGS')
1185  FORMAT(1X,'==========================')
1186  FORMAT(1X,'HIGGS --> GG')
1187  FORMAT(1X,'============')
1180  FORMAT(1X,'P P COLLIDER')
1182  FORMAT(1X,'============')
1181  FORMAT(1X,'P PBAR COLLIDER')
1183  FORMAT(1X,'===============')
5118   FORMAT(1X,'ENERGY        = ',G15.6,' TEV')
5146  FORMAT(1X,'NF_EXT    = ',I5)
1260  FORMAT(1X,'  LO-LAMBDA_',I1,' = ',G15.6,' GEV',2X,
     .          '  LO-ALPHA_S (M_Z) = ',G15.6)
1261  FORMAT(1X,' NLO-LAMBDA_',I1,' = ',G15.6,' GEV',2X,
     .          ' NLO-ALPHA_S (M_Z) = ',G15.6)
1262  FORMAT(1X,'NNLO-LAMBDA_',I1,' = ',G15.6,' GEV',2X,
     .          'NNLO-ALPHA_S (M_Z) = ',G15.6)

1211  FORMAT(1X,'REN-SCALE = ',G15.6,'  GEV',3X,'FAC-SCALE = ',
     .       G15.6,'  GEV')
1212  FORMAT(1X,'REN-SCALE = ',G15.6,' GEV')

5119  FORMAT(1X,'T-MASS        = ',G15.6,' GEV',5X,'B-MASS    = ',
     .       G15.6,' GEV')
51190 FORMAT(1X,'C-MASS        = ',G15.6,' GEV')
51191 FORMAT(1X,'TP-MASS       = ',G15.6,' GEV',5X,'BP-MASS   = ',
     .       G15.6,' GEV')
5121  FORMAT(1X,'M_',A1,'       = ',G15.6,'  GEV')
1219  FORMAT(1X,'M_A       = ',G15.6,'  GEV')
51291 FORMAT(1X,'MSSM (subhpole):  TG(BETA) = ',G9.3)
51292 FORMAT(1X,'MSSM (subh):  TG(BETA) = ',G9.3)
51293 FORMAT(1X,'MSSM (FeynHiggsFast):  TG(BETA) = ',G9.3)
51290 FORMAT(1X,'====')
5132  FORMAT(1X,'G^',A1,'_B     = ',G15.6,2X,'G^',A1,'_T     = ',G15.6)
51320 FORMAT(1X,'G^',A1,'_C     = ',G15.6,2X,'C^',A1,'_G     = ',G15.6)
51322 FORMAT(1X,'G^',A1,'_CTG   = ',G15.6)
51321 FORMAT(1X,'G^',A1,'_BP    = ',G15.6,2X,'G^',A1,'_TP    = ',G15.6)
5133  FORMAT(1X,'G^',A1,'_ST1   = ',G15.6,2X,'G^',A1,'_ST2   = ',G15.6)
5134  FORMAT(1X,'G^',A1,'_SB1   = ',G15.6,2X,'G^',A1,'_SB2   = ',G15.6)
11290 FORMAT(1X,'MSTOP1    = ',G11.6,' GEV',2X,
     .          'MSTOP2    = ',G11.6,' GEV')
11292 FORMAT(1X,'MSBOT1    = ',G11.6,' GEV',2X,
     .          'MSBOT2    = ',G11.6,' GEV')
5170  FORMAT(1X,'SIG0      = ',G15.6,1X,'+- ',G15.6,2X,'PB')
5131  FORMAT(1X,'SIG1_VIRT = ',G15.6,1X,'+- ',G15.6,2X,'PB')
51310 FORMAT(1X,'FAC_ELW   = ',G15.6)
5123  FORMAT(1X,'SIG1_GG   = ',G15.6,1X,'+- ',G15.6,2X,'PB')
5124  FORMAT(1X,'SIG1_GQ   = ',G15.6,1X,'+- ',G15.6,2X,'PB')
5125  FORMAT(1X,'SIG1_QQB  = ',G15.6,1X,'+- ',G15.6,2X,'PB')
5127  FORMAT(1X,'SIG_LO    = ',G15.6,1X,'+- ',G15.6,2X,'PB')
5128  FORMAT(1X,'SIG_NLO   = ',G15.6,1X,'+- ',G15.6,2X,'PB')
5129  FORMAT(1X,'SIG_NNLO  = ',G15.6,1X,'+- ',G15.6,2X,'PB')
5130  FORMAT(1X,'SIG_TOT   = ',G15.6,1X,'+- ',G15.6,2X,'PB')
5140  FORMAT(1X,'GAM_LO    = ',G15.6,2X,'MEV')
5141  FORMAT(1X,'E_VIRT    = ',G15.6)
51410 FORMAT(1X,'E_ELW     = ',G15.6)
5142  FORMAT(1X,'E_GGG     = ',G15.6,1X,'+- ',G15.6)
5143  FORMAT(1X,'E_GQQB    = ',G15.6,1X,'+- ',G15.6)
5148  FORMAT(1X,'E_GTTB    = ',G15.6,1X,'+- ',G15.6)
5144  FORMAT(1X,'E_TOT     = ',G15.6,1X,'+- ',G15.6)
5145  FORMAT(1X,'GAM_NLO   = ',G15.6,1X,'+- ',G15.6,2X,'MEV')

6000  FORMAT(1X,'T-SCALE   = ',G15.6,'GEV',3X,'MSBAR T-MASS = ',
     .       G15.6,' GEV')
6001  FORMAT(1X,'B-SCALE   = ',G15.6,'GEV',3X,'MSBAR B-MASS = ',
     .       G15.6,' GEV')
6002  FORMAT(1X,'C-SCALE   = ',G15.6,'GEV',3X,'MSBAR C-MASS = ',
     .       G15.6,' GEV')

9997  CONTINUE
      ENDIF
9998  CONTINUE
9999  CONTINUE
 
      CLOSE(NIN)
      CLOSE(NOUT)
      STOP
      END
 
      DOUBLE PRECISION FUNCTION D9HGG(TAUH)
C--INTEGRAND FOR GG --> HG; 2-LOOP SOFT GLUON TERM
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMPLEX*16 LI2
      COMMON/MASS/AMH,AMQ,S
      COMMON/FACSC/ISCHEME
      COMMON/PARSC/XKAPM,XKAPQ
      COMMON/SUSY/AMB,AMC,FACT,FACB,FACC,FACCG,FACCTG,ISILH
      COMMON/CONST/ZETA2,ZETA3
      COMMON/PARTS/IPRTGG,IPRTGQ,ICHGG,ICHGQ
      COMMON/HIGGS/IHIGGS
      COMMON/PTY/AA,BB,NLOOP
      SP(X) = DREAL(LI2(DCMPLX(X,0.D0)))
      PI=4.D0*DATAN(1.D0)
      TH=AMH**2/S
      TAU=TH/TAUH
C--TOTAL SUBENERGY
      SHAT=AMH**2/TAUH
C--CALL SCALES
      CALL SCALES(QM,QQ)
C--PLUS-DISTRIBUTIONS
      P0=1.D0/(1.D0-TAUH)
      P1=DLOG(1.D0-TAUH)/(1.D0-TAUH)
      P2=DLOG(1.D0-TAUH)**2/(1.D0-TAUH)
      P3=DLOG(1.D0-TAUH)**3/(1.D0-TAUH)
      O0=0.D0
C     O0=1.D0
      O1=DLOG(1.D0-TAUH)
      O2=DLOG(1.D0-TAUH)**2
      O3=DLOG(1.D0-TAUH)**3
C--LOGARITHMS OF FACTORIZATION-SCALE
      SLG=DLOG(QQ**2/SHAT)
      SLG0=DLOG(QQ**2/AMH**2)

c!!!!!!!
c     SLG = DLOG(AMH**2/SHAT)
c     SLG0= 0
c!!!!!!!

C--COUPLINGS
      DCPL=ALPHAS(QM,NLOOP)**4/PI**2*TH/TAUH**2
      DCPL0=ALPHAS(QM,NLOOP)**4/PI**2*TH
      NC=3
      B0=(11*NC-2*DNF(QM))/12.D0
      IF(IPRTGG.EQ.0)THEN
       D9HGG=(DCPL*DLUMGG(TAU,QQ)-DCPL0*DLUMGG(TH,QQ))
     .       *(8*NC**2*P3 - 4*NC*B0*P2 - 4/3.D0*NC**2*PI**2*P1
     .         + 16*NC**2*ZETA3*P0
     .      + NC*(SLG0**2*(4*NC*P1 - B0*P0)
     .      - SLG0*(12*NC*P2 - 4*B0*P1 - 4*NC*ZETA2*P0)))
     .       -2.D0*DCPL*DLUMGG(TAU,QQ)
     .       *(8*NC**2*O3 - 4*NC*(B0+NC)*O2 - 4/3.D0*NC**2*PI**2*O1
     .      - NC*(SLG0**2*(- 4*NC*O1 + (NC+B0)*O0)
     .      - SLG0*(- 12*NC*O2 + (4*B0+4*NC)*O1 + 4*NC*ZETA2*O0)))
       IF(IHIGGS.EQ.1)D9HGG = D9HGG
     .     + (DCPL*DLUMGG(TAU,QQ)-DCPL0*DLUMGG(TH,QQ))
     .       *6*4/3.D0*NC**2*(P1 - P0/2*SLG0)
     .     + DCPL*DLUMGG(TAU,QQ)
     .       *6*4/3.D0*NC**2*(- 2*O1 + O0*SLG0)
      ELSEIF(IPRTGG.EQ.1)THEN
       D9HGG=(DCPL*DLUMGG(TAU,QQ)-DCPL0*DLUMGG(TH,QQ))
     .       *(8*NC**2*P3 - 4*NC*B0*P2 - 4/3.D0*NC**2*PI**2*P1
     .         + 16*NC**2*ZETA3*P0
     .      + NC*(SLG0**2*(4*NC*P1 - B0*P0)
     .      - SLG0*(12*NC*P2 - 4*B0*P1 - 4*NC*ZETA2*P0)))
       IF(IHIGGS.EQ.1)D9HGG = D9HGG +
     .       (DCPL*DLUMGG(TAU,QQ)-DCPL0*DLUMGG(TH,QQ))
     .       *6*4/3.D0*NC**2*(P1 - P0/2*SLG0)
      ELSEIF(IPRTGG.EQ.2)THEN
       DEL=203/12.D0
       D9HGG=(DCPL*DLUMGG(TAU,QQ)-DCPL0*DLUMGG(TH,QQ))
     .       *(DEL*4/3.D0*NC**2*P1
     .      + NC*(SLG0**2*(-11.D0/3.D0*NC*P0)
     .        - SLG0*(-22.D0/3.D0*NC*P1 + 203.D0/18.D0*NC*P0)))
      ELSEIF(IPRTGG.EQ.3)THEN
       Z = TAUH

       DIFF = 2.D0*DCPL*DLUMGG(TAU,QQ)*(
     .        NC*(SLG0**2*(NC+B0) - SLG0*4*NC*ZETA2))
       IF(IHIGGS.EQ.1)DIFF = DIFF 
     .     + DCPL*DLUMGG(TAU,QQ)*6*4/3.D0*NC**2*SLG0

      RES1 = (24*ZETA2*NC*(Z**4 - 2*Z**3 + 3*Z**2 - 4*Z + 2) + 12*DLOG
     .(Z)**2*NC*( - Z**4 + 4*Z**3 - 3*Z**2 - 1) + 24*DLOG(Z)*
     .DLOG( - Z + 1)*NC*(3*Z**4 - 10*Z**3 + 9*Z**2 - 2*Z + 3)
     . + 6*DLOG(Z)*( - 2*Z**4*B0 - 22*Z**4*NC + 4*Z**3*B0 + 46
     .*Z**3*NC - 6*Z**2*B0 - 35*Z**2*NC + 4*Z*B0 + 22*Z*NC
     . - 2*B0 - 11*NC) + 72*DLOG( - Z + 1)**2*NC*( - Z**4
     . + 2*Z**3 - 3*Z**2 + 4*Z - 2) + 6*DLOG( - Z + 1)*(4*Z
     .**4*B0 + 33*Z**4*NC - 8*Z**3*B0 - 68*Z**3*NC + 12*Z**
     .2*B0 + 70*Z**2*NC - 16*Z*B0 - 76*Z*NC + 8*B0 + 41*NC)
     . + 96*SP( - Z + 1)*Z*NC*( - Z**2 + 1) - 11*Z**4*B0 - 170
     .*Z**4*NC + 44*Z**3*B0 + 284*Z**3*NC - 66*Z**2*B0 - 228*Z
     .**2*NC + 44*Z*B0 + 284*Z*NC - 11*B0 - 170*NC)/(6*(Z - 1)
     .)
 
      RES2 = (6*DLOG(Z)*NC*(Z**4 - 4*Z**3 + 3*Z**2 + 1) + 12*DLOG( - Z
     . + 1)*NC*( - Z**4 + 2*Z**3 - 3*Z**2 + 4*Z - 2) + 3*Z**4*
     .B0 + 22*Z**4*NC - 6*Z**3*B0 - 40*Z**3*NC + 9*Z**2*B0 +
     .36*Z**2*NC - 12*Z*B0 - 46*Z*NC + 6*B0 + 28*NC)/(3*(Z - 1
     .))
 
       TF = DNF(QM)/2
       PGG = 1/(1-Z) + 1/Z - 2 + Z*(1-Z)
       PGGM = 1/(1+Z) - 1/Z - 2 - Z*(1+Z)
       SS2 = SP(Z/(1+Z)) - SP(1/(1+Z))
     .     - (DLOG(1+Z)**2 - DLOG(Z/(1+Z))**2)/2
       PGG1 = 4.D0/3*TF*(-16 + 8*Z + 20.D0/3*Z**2 + 4/3.D0/Z
     .                  - (6+10*Z)*DLOG(Z) - 2*(1+Z)*DLOG(Z)**2)
     .      + 3.D0*TF*(2 - 2*Z + 26/9.D0*(Z**2-1/Z)
     .                - 4/3.D0*(1+Z)*DLOG(Z) - 20/9.D0*PGG)
     .      + 9*(27/2.D0*(1-Z) + 67/9.D0*(Z**2-1/Z)
     .          - (25/3.D0-11/3.D0*Z+44/3.D0*Z**2)*DLOG(Z)
     .          + 4*(1+Z)*DLOG(Z)**2 + 2*PGGM*SS2
     .          + (67/9.D0-4*DLOG(Z)*DLOG(1-Z)+DLOG(Z)**2-2*ZETA2)
     .                                                         *PGG)
       SGG1 = 3.D0*TF*(- 20/9.D0*P0)
     .      + 9*(67/9.D0-2*ZETA2)*P0
c!!!!!!!!!!!!!!!!!!!!!!!
c      PGG1=0
c      SGG1=0
c!!!!!!!!!!!!!!!!!!!!!!!
c      RES1=0
c!!!!!!!!!!!!!!!!!!!!!!!

       RES = NC*SLG0*(SLG0*RES2 - RES1) + SLG0*(- Z*PGG1)/2
       SING = SLG0*(- SGG1)/2
       D9HGG = DCPL*DLUMGG(TAU,QQ) * RES
     .       - DCPL0*DLUMGG(TH,QQ) * SING
     .       + DIFF
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION D9VGG(X)
C--GG --> HG: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(2)
      COMMON/MASS/AMH,AMQ,S
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/GRENZ/IGRENZ
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/INTRO/INTEG
      COMMON/PARTDN/XVAR
      COMMON/PARTS/IPRTGG,IPRTGQ,ICHGG,ICHGQ
      XVAR=X(2)
      ICALL1=ICALL1+1
      ONE=1.D0-EPST
      TH=AMH**2/S
      TAUH=TH+(ONE-TH)*X(1)
      V2GG=D9HGG(TAUH)
      D9VGG=V2GG*(ONE-TH)
      RETURN
      END
 
c***************************************************************
c         NNLO CORRECTIONS
c***************************************************************
      double precision function d2hgg(tauh)
C--INTEGRAND FOR GG --> HG(G): 2-LOOP TERM
      implicit double precision (a-b,d-h,o-z), complex*16 (c)
      complex*16 li2,li3,s12
      common/mass/amh,amq,s
      common/facsc/ischeme
      common/parsc/xkapm,xkapq
      common/susy/amb,amc,fact,facb,facc,faccg,facctg,isilh
      common/const/zeta2,zeta3
      common/parts/iprtgg,iprtgq,ichgg,ichgq
      common/higgs/ihiggs
      common/pty/aa,bb,nloop
      sli2(x) = dreal(li2(dcmplx(x,0.d0)))
      sli3(x) = dreal(li3(dcmplx(x,0.d0)))
      ss12(x) = dreal(s12(dcmplx(x,0.d0)))
      pi=4.d0*datan(1.d0)
      th=amh**2/s
      tau=th/tauh
C--TOTAL SUBENERGY
      shat=amh**2/tauh
C--CALL SCALES
      call scales(qm,qq)
C--PLUS-DISTRIBUTIONS
      p0=1.d0/(1.d0-tauh)
      p1=dlog(1.d0-tauh)/(1.d0-tauh)
      p2=dlog(1.d0-tauh)**2/(1.d0-tauh)
      p3=dlog(1.d0-tauh)**3/(1.d0-tauh)
C--LOGARITHMS OF FACTORIZATION-SCALE
      slg=dlog(qq**2/shat)
      slg0=dlog(qq**2/amh**2)
      nf = 5
      z = tauh
      dlgfac = -slg0

C--COUPLINGS
      dcpl=alphas(qm,nloop)**4/pi**2*th/tauh**2
      dcpl0=alphas(qm,nloop)**4/pi**2*th
      rgg0=(18*dlgfac**2*nf-297*dlgfac**2-60*dlgfac*nf-1620*dlgfac*
     . zeta2+2394*dlgfac-72*nf*zeta2+56*nf+1188*zeta2+6318*zeta3-1212
     . )/36
      rgg1=(108*dlgfac**2+6*dlgfac*nf-99*dlgfac-10*nf-270*zeta2+399)/
     . 3
      rgg2=2*nf-33+108*dlgfac
      rgg3=72
      d2hgg=(dcpl*dlumgg(tau,qq)-dcpl0*dlumgg(th,qq))
     .     *(rgg3*p3 + rgg2*p2 + rgg1*p1 + rgg0*p0)
      if(ihiggs.eq.1)d2hgg = d2hgg
     .               +(dcpl*dlumgg(tau,qq)-dcpl0*dlumgg(th,qq))*6*p1
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     .               +(dcpl*dlumgg(tau,qq)-dcpl0*dlumgg(th,qq))*6*p0
     .               *dlgfac/2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ans5=-10692*dlgfac**2*z**2-5832*dlgfac**2*z+85536*dlgfac**2-
     . 9420*dlgfac*nf*z**4-8088*dlgfac*nf*z**3-432*dlgfac*nf*z**2+
     . 3768*dlgfac*nf*z+5532*dlgfac*nf-69984*dlgfac*z**4*zeta2+392850
     . *dlgfac*z**4-23328*dlgfac*z**3*zeta2+104976*dlgfac*z**3-93312*
     . dlgfac*z**2*zeta2+86184*dlgfac*z**2-139968*dlgfac*z*zeta2+
     . 67392*dlgfac*z-11664*dlgfac*zeta2-306666*dlgfac-4896*nf*z**4*
     . zeta2+17774*nf*z**4-4032*nf*z**3*zeta2+13973*nf*z**3-2592*nf*z
     . **2*zeta2-819*nf*z**2-1152*nf*z*zeta2-9941*nf*z+2304*nf*zeta2-
     . 12923*nf+342144*z**4*zeta2+256608*z**4*zeta3-610659*z**4+95256
     . *z**3*zeta2+23328*z**3*zeta3+40230*z**3+34020*z**2*zeta2+
     . 209952*z**2*zeta3+33156*z**2-15552*z*zeta2+443232*z*zeta3-
     . 127494*z-296460*zeta2-5832*zeta3+490239
      ans4=46656*sli3(-z)*z**4+46656*sli3(-z)*z**3-17496*sli3(-z)*z**2-
     . 11664*sli3(-z)*z-5832*sli3(-z)+46656*sli3((-(z-1))/(z+1))*z**4+
     . 93312*sli3((-(z-1))/(z+1))*z**3+139968*sli3((-(z-1))/(z+1))*z**2
     . +93312*sli3((-(z-1))/(z+1))*z+46656*sli3((-(z-1))/(z+1))-46656*
     . sli3((z-1)/(z+1))*z**4-93312*sli3((z-1)/(z+1))*z**3-139968*sli3((
     . z-1)/(z+1))*z**2-93312*sli3((z-1)/(z+1))*z-46656*sli3((z-1)/(z+1
     . ))+93312*log(-(z-1))**3*z**4+93312*log(-(z-1))**3*z**2+186624*
     . log(-(z-1))**3*z+23328*log(z+1)**2*log(z)*z**4+46656*log(z+1)
     . **2*log(z)*z**3+17496*log(z+1)**2*log(z)*z**2+11664*log(z+1)**
     . 2*log(z)*z+5832*log(z+1)**2*log(z)+46656*ss12(-z)*z**4+93312*
     .ss12(-z)*z**3+34992*ss12(-z)*z**2+23328*ss12(-z)*z+11664*ss12(-z)+
     . 1224*dlgfac**2*nf*z**4+1008*dlgfac**2*nf*z**3+648*dlgfac**2*nf
     . *z**2+288*dlgfac**2*nf*z-576*dlgfac**2*nf-96228*dlgfac**2*z**4
     . -15552*dlgfac**2*z**3+ans5
      ans6=(z-1)
      ans3=ans4*ans6
      ans2=-ans3
      ans10=12*((27*(432*z**4*zeta2-2559*z**4+144*z**3*zeta2-401*z**3
     . +576*z**2*zeta2-529*z**2+864*z*zeta2-660*z+72*zeta2+2027-144*(
     . z**2-z+2)*(z+1)*dlgfac**2*z+144*(z**2+z+1)**2*log(z+1)*log(z))
     . +(1570*z**3-249*z**2+294*z-922)*(z+1)*nf-12*(34*nf*z**3-6*nf*z
     . **2+24*nf*z-16*nf-2970*z**3+3132*z**2-3429*z+2673)*(z+1)*
     . dlgfac)*(z-1)+12*(24*dlgfac*nf*z**3-24*dlgfac*nf*z+1134*dlgfac
     . *z**4-3564*dlgfac*z**3+3402*dlgfac*z**2-972*dlgfac*z+1134*
     . dlgfac+50*nf*z**4-122*nf*z**3+42*nf*z**2+14*nf*z+34*nf-5049*z
     . **4+9855*z**3-8424*z**2+5697*z-2376)*(z+1)*log(z)-36*(27*(9*z
     . **5-19*z**4-3*z**3+19*z**2+3*z+7)+8*(z+1)**2*(z-1)*nf*z)*log(z
     . )**2)*log(-(z-1))
      ans9=-6*(27*(756*z**5*zeta2-4318*z**5-1764*z**4*zeta2+2637*z**4
     . -432*z**3*zeta2+508*z**3+1728*z**2*zeta2-836*z**2+432*z*zeta2+
     . 3278*z+576*zeta2-2333)+(2384*z**4+576*z**3*zeta2-3041*z**3+333
     . *z**2-576*z*zeta2-598*z+1282)*(z+1)*nf-144*(nf*z**3-nf*z+27*z
     . **4-108*z**3+81*z**2+27)*(z+1)*dlgfac**2-12*(50*nf*z**4-122*nf
     . *z**3+42*nf*z**2+14*nf*z+34*nf-5049*z**4+9531*z**3-7533*z**2+
     . 4833*z-2079)*(z+1)*dlgfac)*log(z)+1944*(12*(2*log(-(z-1))+
     . dlgfac)*(z**2+z+1)**2+(11*z**3+21*z**2-12*z-14)*(z+1)-6*(2*z**
     . 3+3*z**2+1)*(2*z+1)*log(z+1)+3*(2*z**4-15*z**2-10*z-5)*log(z))
     . *(z-1)*sli2(-z)+ans10
      ans8=36*((96*(2*log(-(z-1))+dlgfac)*(nf-54)*(z+1)*(z-1)*z-(4*nf
     . *z**4+227*nf*z**3+21*nf*z**2-302*nf*z+68*nf+3861*z**4-7560*z**
     . 3+567*z**2+10638*z-7803))*(z+1)+6*(27*(z**5-11*z**4+12*z**2+12
     . )-8*(z+1)**2*(z-1)*nf*z)*log(z))*sli2(-(z-1))+72*(2*(27*(110*z
     . **3-116*z**2+127*z-99-36*(z**2-z+2)*dlgfac*z)-2*(17*z**3-3*z**
     . 2+12*z-8)*nf)*(z-1)+3*(16*nf*z**3-16*nf*z+837*z**4-2538*z**3+
     . 2511*z**2-810*z+837)*log(z))*(z+1)*log(-(z-1))**2+ans9
      ans7=72*(27*(9*z**5-17*z**4-7*z**3+18*z**2+7*z+4)+10*(z+1)**2*(
     . z-1)*nf*z)*log(z)**3+108*(54*(11*z**5+67*z**4+73*z**3-69*z**2-
     . 75*z-5)-(65*z**2+62*z+2)*(z+1)*(z-1)*nf)*sli3(-(z-1))-972*(3*((
     . 4*z**4+8*z**3+27*z**2+18*z+9)*log(z)**2+2*(2*z**3+3*z**2+1)*(2
     . *z+1)*zeta2)-2*((11*z**3+21*z**2-12*z-14)*(z+1)+12*(z**2+z+1)
     . **2*dlgfac)*log(z))*(z-1)*log(z+1)-18*(24*(27*(3*z**5-8*z**4-3
     . *z**3+8*z**2+3*z+2)+4*(z+1)**2*(z-1)*nf*z)*dlgfac+(132*nf*z**4
     . -351*nf*z**3+117*nf*z**2+52*nf*z+68*nf-10098*z**4+22329*z**3-
     . 18873*z**2+10503*z-4158)*(z+1))*log(z)**2+ans8
      ans1=108*(108*(7*z**5-39*z**4-20*z**3+38*z**2+18*z+14)+(31*z**2
     . +34*z-2)*(z+1)*(z-1)*nf)*ss12(-(z-1))+ans2+ans7
      rggf=ans1/(1296*(z+1)*(z-1))
      d2hgg=d2hgg + dcpl*dlumgg(tau,qq)*rggf
      if(ihiggs.eq.1)then
       rgga = nf*(2*z/3*log(z)**2+z*log(z)-(10*z-1)*(z-1)/6)
     .      - 9*z*log(z)**2-6*z*(z**2-z+2)*log(1-z)
     .      + 3*(2*z**4-4*z**3+13*z**2+z-10)/2*log(z)/(z-1)
     .      + (z-1)*(11*z**2+35*z-154)/4
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     .      - dlgfac/2*6*z*(2-z*(1-z))
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       d2hgg=d2hgg + dcpl*dlumgg(tau,qq)*rgga
      endif
      return
      end
 
      double precision function d2vgg(x)
C--GG --> HG(G): INTEGRAND FOR VEGAS-INTEGRATION
      implicit double precision (a-b,d-h,o-z), complex*16 (c)
      dimension x(2)
      common/mass/amh,amq,s
      common/cut/epst,epsv,reps
      common/grenz/igrenz
      common/calls/icall1,icall2,ifail,if66
      common/intro/integ
      common/partdn/xvar
      common/parts/iprtgg,iprtgq,ichgg,ichgq
      xvar=x(2)
      icall1=icall1+1
      one=1.d0-epst
      th=amh**2/s
      tauh=th+(one-th)*x(1)
      v2gg=d2hgg(tauh)
      d2vgg=v2gg*(one-th)
      return
      end

      double precision function d2hgq(tauh)
C--INTEGRAND FOR GQ --> HQ(G): 2-LOOP TERM
      implicit double precision (a-b,d-h,o-z), complex*16 (c)
      complex*16 li2,li3,s12
      common/mass/amh,amq,s
      common/facsc/ischeme
      common/parsc/xkapm,xkapq
      common/susy/amb,amc,fact,facb,facc,faccg,facctg,isilh
      common/const/zeta2,zeta3
      common/parts/iprtgg,iprtgq,ichgg,ichgq
      common/higgs/ihiggs
      common/pty/aa,bb,nloop
      sli2(x) = dreal(li2(dcmplx(x,0.d0)))
      sli3(x) = dreal(li3(dcmplx(x,0.d0)))
      ss12(x) = dreal(s12(dcmplx(x,0.d0)))
      pi=4.d0*datan(1.d0)
      th=amh**2/s
      tau=th/tauh
C--TOTAL SUBENERGY
      shat=amh**2/tauh
C--CALL SCALES
      call scales(qm,qq)
C--LOGARITHMS OF FACTORIZATION-SCALE
      slg=dlog(qq**2/shat)
      slg0=dlog(qq**2/amh**2)
      nf = 5
      z = tauh
      dlgfac = -slg0

C--COUPLINGS
      dcpl=alphas(qm,nloop)**4/pi**2*th/tauh**2
      ans5=77112*ss12(-(z-1))*z**2+110160*ss12(-(z-1))*z+73008*ss12(-(z-
     . 1))-216*dlgfac**2*nf*z**2+432*dlgfac**2*nf*z-432*dlgfac**2*nf-
     . 3888*dlgfac**2*z**3+864*dlgfac**2*z**2-31320*dlgfac**2*z+37260
     . *dlgfac**2+1368*dlgfac*nf*z**2-2736*dlgfac*nf*z+2088*dlgfac*nf
     . +4704*dlgfac*z**3+3456*dlgfac*z**2*zeta2+13284*dlgfac*z**2-
     . 14688*dlgfac*z*zeta2+109944*dlgfac*z+6912*dlgfac*zeta2-147588*
     . dlgfac-2148*nf*z**2+5016*nf*z-3180*nf+15696*z**3*zeta2-6148*z
     . **3+15228*z**2*zeta2-36504*z**2*zeta3-7179*z**2+84672*z*zeta2+
     . 61344*z*zeta3-197580*z-109692*zeta2-73008*zeta3+210115-54*(2*(
     . 186*dlgfac+nf)*(z**2-2*z+2)+288*z**3+111*z**2+2278*z-2592-2*(
     . 553*z**2+190*z+642)*log(z))*log(-(z-1))**2+72*(54*(2*log(z)-
     . dlgfac-2*log(-(z-1)))*(z**2+2*z+2)+4*z**3+33*z**2+222*z+166)*
     . sli2(-z)
      ans4=2808*log(-(z-1))*nf+9408*log(-(z-1))*z**3+6912*log(-(z-1))
     . *z**2*zeta2+31752*log(-(z-1))*z**2-29376*log(-(z-1))*z*zeta2+
     . 206064*log(-(z-1))*z+13824*log(-(z-1))*zeta2-284052*log(-(z-1)
     . )+5220*log(z)**3*z**2+5112*log(z)**3*z+3888*log(z)**3-14472*
     . log(z)**2*dlgfac*z**2-13824*log(z)**2*dlgfac*z-11664*log(z)**2
     . *dlgfac-216*log(z)**2*nf*z**2+432*log(z)**2*nf*z-432*log(z)**2
     . *nf-8136*log(z)**2*z**3+8424*log(z)**2*z**2-52704*log(z)**2*z+
     . 30132*log(z)**2+12096*log(z)*dlgfac**2*z**2+10800*log(z)*
     . dlgfac**2*z+11664*log(z)*dlgfac**2+432*log(z)*dlgfac*nf*z**2-
     . 864*log(z)*dlgfac*nf*z+864*log(z)*dlgfac*nf+14400*log(z)*
     . dlgfac*z**3-9504*log(z)*dlgfac*z**2+74304*log(z)*dlgfac*z-
     . 60264*log(z)*dlgfac-1368*log(z)*nf*z**2+2736*log(z)*nf*z-2088*
     . log(z)*nf-14784*log(z)*z**3-44064*log(z)*z**2*zeta2-35730*log(
     . z)*z**2-44064*log(z)*z*zeta2-172800*log(z)*z-38016*log(z)*
     . zeta2+146508*log(z)+ans5
      ans3=-7776*log(-(z-1))*log(z+1)*log(z)*z**2-15552*log(-(z-1))*
     . log(z+1)*log(z)*z-15552*log(-(z-1))*log(z+1)*log(z)-33696*log(
     . -(z-1))*log(z)**2*z**2-18144*log(-(z-1))*log(z)**2*z-32832*log
     . (-(z-1))*log(z)**2+57888*log(-(z-1))*log(z)*dlgfac*z**2+24192*
     . log(-(z-1))*log(z)*dlgfac*z+65664*log(-(z-1))*log(z)*dlgfac+
     . 432*log(-(z-1))*log(z)*nf*z**2-864*log(-(z-1))*log(z)*nf*z+864
     . *log(-(z-1))*log(z)*nf+28800*log(-(z-1))*log(z)*z**3-12528*log
     . (-(z-1))*log(z)*z**2+146880*log(-(z-1))*log(z)*z-119880*log(-(
     . z-1))*log(z)-6696*log(-(z-1))*dlgfac**2*z**2+13392*log(-(z-1))
     . *dlgfac**2*z-13392*log(-(z-1))*dlgfac**2-432*log(-(z-1))*
     . dlgfac*nf*z**2+864*log(-(z-1))*dlgfac*nf*z-864*log(-(z-1))*
     . dlgfac*nf-15552*log(-(z-1))*dlgfac*z**3-432*log(-(z-1))*dlgfac
     . *z**2-133056*log(-(z-1))*dlgfac*z+148392*log(-(z-1))*dlgfac+
     . 1944*log(-(z-1))*nf*z**2-3456*log(-(z-1))*nf*z+ans4
      ans2=72*(4*z**3+33*z**2+222*z+166-54*(z**2+2*z+2)*dlgfac+81*(z
     . **2+2*z+2)*log(z))*log(z+1)*log(z)+67824*sli2(-(z-1))*log(-(z-1
     . ))*z**2+144288*sli2(-(z-1))*log(-(z-1))*z+35424*sli2(-(z-1))*log
     . (-(z-1))+16416*sli2(-(z-1))*log(z)*z**2+13824*sli2(-(z-1))*log(z
     . )*z+26784*sli2(-(z-1))*log(z)+34992*sli2(-(z-1))*dlgfac*z**2+
     . 69984*sli2(-(z-1))*dlgfac*z+19872*sli2(-(z-1))*dlgfac+5328*sli2(-
     . (z-1))*z**3-25596*sli2(-(z-1))*z**2-31536*sli2(-(z-1))*z+93780*
     . sli2(-(z-1))-73656*sli3(-(z-1))*z**2-163728*sli3(-(z-1))*z-47088*
     .sli3(-(z-1))-3888*sli3(-z)*z**2-7776*sli3(-z)*z-7776*sli3(-z)+7776
     . *sli3((-(z-1))/(z+1))*z**2+15552*sli3((-(z-1))/(z+1))*z+15552*
     . sli3((-(z-1))/(z+1))-7776*sli3((z-1)/(z+1))*z**2-15552*sli3((z-1)
     . /(z+1))*z-15552*sli3((z-1)/(z+1))-13212*log(-(z-1))**3*z**2+
     . 26424*log(-(z-1))**3*z-26424*log(-(z-1))**3+ans3
      ans1=-ans2
      rgqf=ans1/1944
      d2hgq=dcpl*dlumgq(tau,qq)*rgqf
      if(ihiggs.eq.1)then
       rgqa = 2*(-2*z+2+z**2)/3*log(1-z)-28*z/9*log(z)**2
     .      + (22/3.d0+10*z-z**2/3)*log(z)+17*z**2/6-191*z/9+337/18.d0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     .      + dlgfac/2*2*(2-2*z+z**2)/3
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       d2hgq=d2hgq + dcpl*dlumgq(tau,qq)*rgqa
      endif
      return
      end
 
      double precision function d2vgq(x)
C--GQ --> HQ(G): INTEGRAND FOR VEGAS-INTEGRATION
      implicit double precision (a-b,d-h,o-z), complex*16 (c)
      dimension x(2)
      common/mass/amh,amq,s
      common/cut/epst,epsv,reps
      common/grenz/igrenz
      common/calls/icall1,icall2,ifail,if66
      common/intro/integ
      common/partdn/xvar
      common/parts/iprtgg,iprtgq,ichgg,ichgq
      xvar=x(2)
      icall1=icall1+1
      one=1.d0-epst
      th=amh**2/s
      tauh=th+(one-th)*x(1)
      v2gq=d2hgq(tauh)
      d2vgq=v2gq*(one-th)
      return
      end

      double precision function d2hqqb(tauh)
C--INTEGRAND FOR QQB --> HG(G): 2-LOOP TERM
      implicit double precision (a-b,d-h,o-z), complex*16 (c)
      complex*16 li2,li3,s12
      common/mass/amh,amq,s
      common/facsc/ischeme
      common/parsc/xkapm,xkapq
      common/susy/amb,amc,fact,facb,facc,faccg,facctg,isilh
      common/const/zeta2,zeta3
      common/parts/iprtgg,iprtgq,ichgg,ichgq
      common/higgs/ihiggs
      common/pty/aa,bb,nloop
      sli2(x) = dreal(li2(dcmplx(x,0.d0)))
      sli3(x) = dreal(li3(dcmplx(x,0.d0)))
      ss12(x) = dreal(s12(dcmplx(x,0.d0)))
      pi=4.d0*datan(1.d0)
      th=amh**2/s
      tau=th/tauh
C--TOTAL SUBENERGY
      shat=amh**2/tauh
C--CALL SCALES
      call scales(qm,qq)
C--LOGARITHMS OF FACTORIZATION-SCALE
      slg=dlog(qq**2/shat)
      slg0=dlog(qq**2/amh**2)
      nf = 5
      z = tauh
      dlgfac = -slg0

C--COUPLINGS
      dcpl=alphas(qm,nloop)**4/pi**2*th/tauh**2
      ans4=-6*((384*z**2-967*z-75-8*(z-1)**2*nf)*(z-1)+18*(z+2)**2*
     . log(z)**2-8*(8*z**2-25*z-19)*(z-1)*dlgfac+2*(64*z**3-189*z**2+
     . 72*z+76-18*(z+2)**2*dlgfac)*log(z))*log(-(z-1))+3*(768*z**3-
     . 288*z*zeta2+1448*z-288*zeta2+135-(72*zeta2+2111)*z**2+18*(z+2)
     . **2*dlgfac**2-8*(4*z**3-12*z**2+15*z-3)*nf-2*(64*z**3-141*z**2
     . +24*z+108)*dlgfac)*log(z)
      ans3=2*(3*(3*(log(z)**3+24*ss12(-(z-1))-24*sli3(-(z-1)))*(z+2)**2
     . -(208*z**3*zeta2-783*z**3-582*z**2*zeta2-12*z**2*zeta3+2708*z
     . **2+44*zeta2-24*zeta3+373)-18*(z+3)*(z-1)*dlgfac**2+6*(log(z+1
     . )**2*log(z)+2*ss12(-z)+3*sli3(-z))*(z**2+2*z+2))-2*(41*z**2-88*z
     . +23)*(z-1)*nf+18*(4*zeta3+383-64*zeta2)*z+12*(2*(13*z**2-35*z-
     . 14)*(z-1)+9*(z+2)**2*log(z))*log(-(z-1))**2)-3*(396*z**2-1055*
     . z-63-24*(z-1)**2*nf)*(z-1)*dlgfac+6*(2*(44*z**3-81*z**2+39*z+
     . 27)-9*(z+2)**2*dlgfac)*log(z)**2-18*((3*log(z)**2-2*zeta2)*(z
     . **2+2*z+2)+4*(6*z**2+z+2)*log(z)*z)*log(z+1)+36*((2*log(z+1)-3
     . *log(z))*(z**2+2*z+2)-2*(6*z**2+z+2)*z)*sli2(-z)+12*(9*(log(z)+
     . 2*dlgfac+4*log(-(z-1)))*(z+2)**2+10*z**3-87*z**2+42*z+12)*sli2(
     . -(z-1))+ans4
      ans2=2*ans3
      ans1=-ans2
      rqqbf=ans1/243
      d2hqqb=dcpl*dlumqqb(tau,qq)*rqqbf
      if(ihiggs.eq.1)then
       rqqba = nf*(-32*z*log(z)/27+16*(z**2-1)/27)+32*z/27*log(z)**2
     .       + 32*(3+8*z)/27*log(z)-16*(z-1)/27*(z**2+10*z+11)
       d2hqqb=d2hqqb + dcpl*dlumqqb(tau,qq)*rqqba
      endif
      return
      end

      double precision function d2vqqb(x)
C--QQB --> HG(G): INTEGRAND FOR VEGAS-INTEGRATION
      implicit double precision (a-b,d-h,o-z), complex*16 (c)
      dimension x(2)
      common/mass/amh,amq,s
      common/cut/epst,epsv,reps
      common/grenz/igrenz
      common/calls/icall1,icall2,ifail,if66
      common/intro/integ
      common/partdn/xvar
      common/parts/iprtgg,iprtgq,ichgg,ichgq
      xvar=x(2)
      icall1=icall1+1
      one=1.d0-epst
      th=amh**2/s
      tauh=th+(one-th)*x(1)
      v2qqb=d2hqqb(tauh)
      d2vqqb=v2qqb*(one-th)
      return
      end
 
      double precision function d2hqqp(tauh)
C--INTEGRAND FOR QQB --> HG(G): 2-LOOP TERM
      implicit double precision (a-b,d-h,o-z), complex*16 (c)
      complex*16 li2,li3,s12
      common/mass/amh,amq,s
      common/facsc/ischeme
      common/parsc/xkapm,xkapq
      common/susy/amb,amc,fact,facb,facc,faccg,facctg,isilh
      common/const/zeta2,zeta3
      common/parts/iprtgg,iprtgq,ichgg,ichgq
      common/higgs/ihiggs
      common/pty/aa,bb,nloop
      sli2(x) = dreal(li2(dcmplx(x,0.d0)))
      sli3(x) = dreal(li3(dcmplx(x,0.d0)))
      ss12(x) = dreal(s12(dcmplx(x,0.d0)))
      pi=4.d0*datan(1.d0)
      th=amh**2/s
      tau=th/tauh
C--TOTAL SUBENERGY
      shat=amh**2/tauh
C--CALL SCALES
      call scales(qm,qq)
C--LOGARITHMS OF FACTORIZATION-SCALE
      slg=dlog(qq**2/shat)
      slg0=dlog(qq**2/amh**2)
      nf = 5
      z = tauh
      dlgfac = -slg0

C--COUPLINGS
      dcpl=alphas(qm,nloop)**4/pi**2*th/tauh**2
      rqqpf=(-2*(2*(log(z)**3+24*ss12(-(z-1))-24*sli3(-(z-1)))*(z+2)**2
     . +3*(16*z*zeta2-11*z+48*zeta2-105)*(z-1)-12*(z+3)*(z-1)*dlgfac
     . **2+9*(5*z+17)*(z-1)*dlgfac-24*(2*(z+3)*(z-1)-(z+2)**2*log(z))
     . *log(-(z-1))**2-6*(2*(z**2+4*z-3)+(z+2)**2*dlgfac)*log(z)**2+
     . 12*((log(z)+2*dlgfac+4*log(-(z-1)))*(z+2)**2-(z**2+4*z-6))*sli2
     . (-(z-1))+6*(3*(5*z+17)*(z-1)-2*(z+2)**2*log(z)**2-8*(z+3)*(z-1
     . )*dlgfac+2*(5*z**2+8*z-12+2*(z+2)**2*dlgfac)*log(z))*log(-(z-1
     . ))-3*(32*z*zeta2+44*z+32*zeta2-59+(8*zeta2+29)*z**2-2*(z+2)**2
     . *dlgfac**2-2*(5*z**2+8*z-12)*dlgfac)*log(z)))/27
      d2hqqp=dcpl*dlumqqp(tau,qq)*rqqpf
      if(ihiggs.eq.1)then
       rqqpa =-16*z/9*log(z)**2+(16*z/3+32/9.d0)*log(z)+8*(z-1)/9*(z-11)
       d2hqqp=d2hqqp + dcpl*dlumqqp(tau,qq)*rqqpa
      endif
      return
      end
 
      double precision function d2vqqp(x)
C--QQP --> HQQP: INTEGRAND FOR VEGAS-INTEGRATION
      implicit double precision (a-b,d-h,o-z), complex*16 (c)
      dimension x(2)
      common/mass/amh,amq,s
      common/cut/epst,epsv,reps
      common/grenz/igrenz
      common/calls/icall1,icall2,ifail,if66
      common/intro/integ
      common/partdn/xvar
      common/parts/iprtgg,iprtgq,ichgg,ichgq
      xvar=x(2)
      icall1=icall1+1
      one=1.d0-epst
      th=amh**2/s
      tauh=th+(one-th)*x(1)
      v2qqp=d2hqqp(tauh)
      d2vqqp=v2qqp*(one-th)
      return
      end

      double precision function d2hqq(tauh)
C--INTEGRAND FOR QQB --> HG(G): 2-LOOP TERM
      implicit double precision (a-b,d-h,o-z), complex*16 (c)
      complex*16 li2,li3,s12
      common/mass/amh,amq,s
      common/facsc/ischeme
      common/parsc/xkapm,xkapq
      common/susy/amb,amc,fact,facb,facc,faccg,facctg,isilh
      common/const/zeta2,zeta3
      common/parts/iprtgg,iprtgq,ichgg,ichgq
      common/higgs/ihiggs
      common/pty/aa,bb,nloop
      sli2(x) = dreal(li2(dcmplx(x,0.d0)))
      sli3(x) = dreal(li3(dcmplx(x,0.d0)))
      ss12(x) = dreal(s12(dcmplx(x,0.d0)))
      pi=4.d0*datan(1.d0)
      th=amh**2/s
      tau=th/tauh
C--TOTAL SUBENERGY
      shat=amh**2/tauh
C--CALL SCALES
      call scales(qm,qq)
C--LOGARITHMS OF FACTORIZATION-SCALE
      slg=dlog(qq**2/shat)
      slg0=dlog(qq**2/amh**2)
      nf = 5
      z = tauh
      dlgfac = -slg0

C--COUPLINGS
      dcpl=alphas(qm,nloop)**4/pi**2*th/tauh**2
      rqqf=(2*(2*((24*sli3(-(z-1))-log(z)**3)*(z+2)**2-(24*z*zeta2-27*
     . z+72*zeta2-160)*(z-1)+6*(z+3)*(z-1)*dlgfac**2)-9*(5*z+17)*(z-1
     . )*dlgfac-4*(13*z**2+46*z+50)*ss12(-(z-1))+24*(2*(z+3)*(z-1)-(z+
     . 2)**2*log(z))*log(-(z-1))**2+2*(9*z**2+26*z-18+3*(z+2)**2*
     . dlgfac)*log(z)**2+(24*z**2*zeta2+69*z**2+96*z*zeta2+124*z+96*
     . zeta2-177-6*(z+2)**2*dlgfac**2-6*(5*z**2+8*z-12)*dlgfac)*log(z
     . )-12*((log(z)+2*dlgfac+4*log(-(z-1)))*(z+2)**2-(z**2+4*z-6))*
     . sli2(-(z-1))-6*(3*(5*z+17)*(z-1)-2*(z+2)**2*log(z)**2-8*(z+3)*(
     . z-1)*dlgfac+2*(5*z**2+8*z-12+2*(z+2)**2*dlgfac)*log(z))*log(-(
     . z-1))))/27
      d2hqq=dcpl*dlumqq(tau,qq)*rqqf
      if(ihiggs.eq.1)then
       rqqa = -64*z/27*log(z)**2+(176*z/27+32/9.d0)*log(z)
     .      + 8*(z-1)/27*(3*z-37)
       d2hqq=d2hqq + dcpl*dlumqq(tau,qq)*rqqa
      endif
      return
      end
 
      double precision function d2vqq(x)
C--QQ --> HQQ: INTEGRAND FOR VEGAS-INTEGRATION
      implicit double precision (a-b,d-h,o-z), complex*16 (c)
      dimension x(2)
      common/mass/amh,amq,s
      common/cut/epst,epsv,reps
      common/grenz/igrenz
      common/calls/icall1,icall2,ifail,if66
      common/intro/integ
      common/partdn/xvar
      common/parts/iprtgg,iprtgq,ichgg,ichgq
      xvar=x(2)
      icall1=icall1+1
      one=1.d0-epst
      th=amh**2/s
      tauh=th+(one-th)*x(1)
      v2qq=d2hqq(tauh)
      d2vqq=v2qq*(one-th)
      return
      end

c***************************************************************
c         END OF NNLO CORRECTIONS
c***************************************************************
 
      DOUBLE PRECISION FUNCTION DHGG(TAUH,XX)
C--INTEGRAND FOR GG --> HG; XX = MASS DEPENDENT PART
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/MASS/AMH,AMQ,S
      COMMON/FACSC/ISCHEME
      COMMON/PARSC/XKAPM,XKAPQ
      COMMON/SUSY/AMB,AMC,FACT,FACB,FACC,FACCG,FACCTG,ISILH
      COMMON/PARTS/IPRTGG,IPRTGQ,ICHGG,ICHGQ
      COMMON/PTY/AA,BB,NLOOP
      PI=4.D0*DATAN(1.D0)
      TH=AMH**2/S
      TAU=TH/TAUH
C--TOTAL SUBENERGY
      SHAT=AMH**2/TAUH
C--CALL SCALES
      CALL SCALES(QM,QQ)
C--ALTARELLI-PARISI KERNEL P_GG WITHOUT PLUS-DISTRIBUTION
C--AND DELTA-TERM
      PGGN=6.D0*(1.D0/TAUH-2.D0+TAUH*(1.D0-TAUH))
C--PLUS-DISTRIBUTION OF P_GG
      PGGP=6.D0/(1.D0-TAUH)
C--LOGARITHMS OF FACTORIZATION-SCALE
      SLG=DLOG(QQ**2/SHAT)
      SLG0=DLOG(QQ**2/AMH**2)
C--COUPLINGS
      DCPL=ALPHAS(QM,NLOOP)**3/PI*TH/TAUH**2
      DCPL0=ALPHAS(QM,NLOOP)**3/PI*TH
      IF(ISCHEME.EQ.1)THEN
C--SHIFT FOR DIS-FACTORIZATION SCHEME
        FGG=-DNF(QQ)
     .       *((TAUH**2+(1.D0-TAUH)**2)*DLOG((1.D0-TAUH)/TAUH)
     .       +8.D0*TAUH*(1.D0-TAUH)-1.D0)
      ELSE
C--SHIFT FOR MSBAR-FACTORIZATION SCHEME
        FGG=0.D0
      ENDIF
      IF(IPRTGG.EQ.0)THEN
       DHGG=DCPL*DLUMGG(TAU,QQ)*(-SLG*TAUH*PGGN-TAUH*FGG
     .     -12.D0*TAUH*(2.D0-TAUH*(1.D0-TAUH))*DLOG(1.D0-TAUH)
     .     +XX)
     .     -(DCPL*TAUH*DLUMGG(TAU,QQ)*SLG
     .         -DCPL0*DLUMGG(TH,QQ)*SLG0)*PGGP
     .     +12.D0*(DCPL*DLUMGG(TAU,QQ)-DCPL0*DLUMGG(TH,QQ))
     .             *DLOG(1.D0-TAUH)/(1.D0-TAUH)
      ELSEIF(IPRTGG.EQ.1)THEN
C--SCHEME ALPHA AND BETA
       DHGG=12.D0*(DCPL*DLUMGG(TAU,QQ)-DCPL0*DLUMGG(TH,QQ))
     .             *DLOG(1.D0-TAUH)/(1.D0-TAUH)
     .     -(DCPL*DLUMGG(TAU,QQ)
     .         -DCPL0*DLUMGG(TH,QQ))*SLG0*PGGP
      ELSEIF(IPRTGG.EQ.2)THEN
C--SCHEME GAMMA - BETA
       DHGG=-12.D0*DCPL*DLUMGG(TAU,QQ)
     .     *(2.D0*DLOG(1.D0-TAUH))
C    .     *(2.D0*DLOG(1.D0-TAUH) - SLG0)
      ELSEIF(IPRTGG.EQ.3)THEN
C--SCHEME DELTA - GAMMA
       DHGG=6.D0*DCPL*DLUMGG(TAU,QQ)*SLG0*TAUH*(2-TAUH*(1-TAUH))
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DHGQ(TAUH,XX)
C--INTEGRAND FOR GQ --> HQ; XX = MASS DEPENDENT PART
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/MASS/AMH,AMQ,S
      COMMON/FACSC/ISCHEME
      COMMON/PARSC/XKAPM,XKAPQ
      COMMON/PARTS/IPRTGG,IPRTGQ,ICHGG,ICHGQ
      COMMON/PTY/AA,BB,NLOOP
      PI=4.D0*DATAN(1.D0)
      TH=AMH**2/S
      TAU=TH/TAUH
C--TOTAL SUBENERGY
      SHAT=AMH**2/TAUH
C--CALL SCALES
      CALL SCALES(QM,QQ)
C--ALTARELLI-PARISI KERNEL P_GQ
      PGQTH=4.D0/3.D0*(1.D0+(1.D0-TAUH)**2)
C--LOGARITHM OF FACTORIZATION-SCALE
      SLG=DLOG(QQ**2/SHAT)
C--COUPLINGS
      DCPL=ALPHAS(QM,NLOOP)**3/PI*TH/TAUH**2
      DCPL0=ALPHAS(QM,NLOOP)**3/PI*TH
      IF(ISCHEME.EQ.1)THEN
C--SHIFT FOR DIS-FACTORIZATION SCHEME
        FGQ=-4.D0/3.D0*((1.D0+TAUH**2)/(1.D0-TAUH)
     .      *(DLOG((1.D0-TAUH)/TAUH)-3.D0/4.D0)
     .      +(9.D0+5.D0*TAUH)/4.D0)
      ELSE
C--SHIFT FOR MSBAR-FACTORIZATION SCHEME
        FGQ=0.D0
      ENDIF
      IF(IPRTGQ.EQ.0)THEN
        DHGQ=DCPL*DLUMGQ(TAU,QQ)*((-SLG/2.D0+DLOG(1.D0-TAUH))*PGQTH
     .                        +XX)
     .      -(DCPL*TAUH*DLUMGQ(TAU,QQ)-DCPL0*DLUMGQ(TH,QQ))/2.D0*FGQ
      ELSEIF(IPRTGQ.EQ.1)THEN
        DHGQ=-DCPL*DLUMGQ(TAU,QQ)*SLG/2.D0*PGQTH
     .      -(DCPL*TAUH*DLUMGQ(TAU,QQ)-DCPL0*DLUMGQ(TH,QQ))/2.D0*FGQ
      ELSEIF(IPRTGQ.EQ.2)THEN
        DHGQ=DCPL*DLUMGQ(TAU,QQ)*DLOG(1.D0-TAUH)*PGQTH
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DHQQ(TAUH,XX)
C--INTEGRAND FOR QQBAR --> HG; XX = MASS DEPENDENT PART
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/PARSC/XKAPM,XKAPQ
      COMMON/MASS/AMH,AMQ,S
      COMMON/PTY/AA,BB,NLOOP
      PI=4.D0*DATAN(1.D0)
      TH=AMH**2/S
      TAU=TH/TAUH
C--TOTAL SUBENERGY
      SHAT=AMH**2/TAUH
C--CALL SCALES
      CALL SCALES(QM,QQ)
C--COUPLING
      DCPL=ALPHAS(QM,NLOOP)**3/PI*TH/TAUH**2
      DHQQ=DCPL*DLUMQQB(TAU,QQ)*XX
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DGG(TAUH)
C--GG --> HG: INTEGRAND FOR TAU-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
C--MASS DEPENDENT PART
      HH=D0GG(TAUH)
      DGG=DHGG(TAUH,HH)
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DGQ(TAUH)
C--GQ --> HQ: INTEGRAND FOR TAU-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/MASS/AMH,AMQ,S
C--MASS DEPENDENT PART
      HH=D0GQ(TAUH)
      DGQ=DHGQ(TAUH,HH)
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DQQ(TAUH)
C--QQBAR --> HG: INTEGRAND FOR TAU-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/MASS/AMH,AMQ,S
C--MASS DEPENDENT PART
      HH=D0QQ(TAUH)
      DQQ=DHQQ(TAUH,HH)
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DLUGG(X)
C--INTEGRAND FOR VEGAS-INTEGRATION OF GG-LUMINOSITY
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/PARTDN/XVAR
      COMMON/LUGG/VTAU,VSC
      XVAR=X
      TAU=VTAU
      Q=VSC
      ICALL1=ICALL1+1
      DLUGG=DLUMGG(TAU,Q)
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DLUGQ(X)
C--INTEGRAND FOR VEGAS-INTEGRATION OF GQ-LUMINOSITY
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/PARTDN/XVAR
      COMMON/LUGG/VTAU,VSC
      XVAR=X
      TAU=VTAU
      Q=VSC
      ICALL1=ICALL1+1
      DLUGQ=DLUMGQ(TAU,Q)
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DVGG(X)
C--GG --> HG: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(3)
      COMMON/MASS/AMH,AMQ,S
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/GRENZ/IGRENZ
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/INTRO/INTEG
      COMMON/PARTDN/XVAR
      COMMON/PARTS/IPRTGG,IPRTGQ,ICHGG,ICHGQ
      IF(INTEG.EQ.3) THEN
       IF(IGRENZ.EQ.1.OR.ICHGG.EQ.2) THEN
        XVAR=X(2)
       ELSE
        XVAR=X(3)
       ENDIF
      ENDIF
      ICALL1=ICALL1+1
      ONE=1.D0-EPST
      TH=AMH**2/S
      TAUH=TH+(ONE-TH)*X(1)
      IF(IGRENZ.EQ.1.OR.ICHGG.EQ.2) THEN
        VGG=DGG(TAUH)
        DVGG=VGG*(ONE-TH)
      ELSE
        V=EPSV+(1.D0-2.D0*EPSV)*X(2)
        HH=DIGG(TAUH,V)
        VGG=DHGG(TAUH,HH)
        DVGG=VGG*(1.D0-2.D0*EPSV)*(ONE-TH)
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DVGQ(X)
C--GQ --> HQ: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(3)
      COMMON/MASS/AMH,AMQ,S
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/GRENZ/IGRENZ
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/INTRO/INTEG
      COMMON/PARTDN/XVAR
      COMMON/PARTS/IPRTGG,IPRTGQ,ICHGG,ICHGQ
      IF(INTEG.EQ.3) THEN
       IF(IGRENZ.EQ.1.OR.ICHGQ.EQ.2) THEN
        XVAR=X(2)
       ELSE
        XVAR=X(3)
       ENDIF
      ENDIF
      ICALL1=ICALL1+1
      ONE=1.D0-EPST
      TH=AMH**2/S
      TAUH=TH+(ONE-TH)*X(1)
      IF(IGRENZ.EQ.1.OR.ICHGQ.EQ.2) THEN
        VGQ=DGQ(TAUH)
        DVGQ=VGQ*(ONE-TH)
      ELSE
        V=EPSV+(1.D0-2.D0*EPSV)*X(2)
        HH=DIGQ(TAUH,V)
        VGQ=DHGQ(TAUH,HH)
        DVGQ=VGQ*(1.D0-2.D0*EPSV)*(ONE-TH)
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION DVQQ(X)
C--QQBAR --> HG: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(3)
      COMMON/MASS/AMH,AMQ,S
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/GRENZ/IGRENZ
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/INTRO/INTEG
      COMMON/PARTDN/XVAR
      IF(INTEG.EQ.3) THEN
       IF(IGRENZ.EQ.1) THEN
        XVAR=X(2)
       ELSE
        XVAR=X(3)
       ENDIF
      ENDIF
      ICALL1=ICALL1+1
      ONE=1.D0-EPST
      TH=AMH**2/S
      TAUH=TH+(ONE-TH)*X(1)
      IF(IGRENZ.EQ.1)THEN
        VQQ=DQQ(TAUH)
        DVQQ=VQQ*(ONE-TH)
      ELSE
        V=EPSV+(1.D0-2.D0*EPSV)*X(2)
        HH=DIQQ(TAUH,V)
        VQQ=DHQQ(TAUH,HH)
        DVQQ=VQQ*(1.D0-2.D0*EPSV)*(ONE-TH)
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION D0QQ(TAUH)
C--QQBAR --> HG
C--MASS DEPENDENT PART: INTEGRATION OVER SCATTERING ANGLE
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/GRENZ/IGRENZ
      COMMON/INTGRLP/TAUH1,IFUN
      COMMON/DCADRE/DO,UP,AERR,RERR
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/MASS/AMH,AMQ,S
      COMMON/SUSY/AMB,AMC,FACT,FACB,FACC,FACCG,FACCTG,ISILH
      COMMON/BORN/CF0,CF0T,CF0B,CF0C,CF0G
      COMMON/BORN4/CF0T4,CF0B4
      COMMON/SM4/AMT4,AMB4,FACT4,FACB4,ISM4,IGGELW
      EXTERNAL DFUN
      PI=4.D0*DATAN(1.D0)
      TAUH1=TAUH
      IFUN=3
      ICALL1=ICALL1+1
      IF(IGRENZ.EQ.1)THEN
        RHO=AMH**2/AMB**2
        DLR=DLOG(RHO)
C--LIMIT OF HEAVY TOP MASS AND SMALL BOTTOM MASS
        D0QQ=CDABS((CF0T+CF0G)/CF0)**2*32.D0/27.D0*(1.D0-TAUH)**3
     .      +DREAL(DCONJG(CF0T+CF0G)*CF0B)/CDABS(CF0)**2*
     .       128.D0/27.D0*TAUH*(1.D0-TAUH)**2*DLOG(TAUH)/DLR
     .      +CDABS(CF0B/CF0)**2*
     .       128.D0/27.D0*TAUH**2*(1.D0-TAUH)*DLOG(TAUH)**2/DLR**2
        IF(ISILH.NE.0)
     .  D0QQ=(CDABS(CF0T+CF0G)**2-CF0G**2)/(CDABS(CF0)**2-CF0G**2)
     .      *32.D0/27.D0*(1.D0-TAUH)**3
     .      +DREAL(DCONJG(CF0T+CF0G)*CF0B)/(CDABS(CF0)**2-CF0G**2)*
     .       128.D0/27.D0*TAUH*(1.D0-TAUH)**2*DLOG(TAUH)/DLR
     .      +CDABS(CF0B)**2/(CDABS(CF0)**2-CF0G**2)*
     .       128.D0/27.D0*TAUH**2*(1.D0-TAUH)*DLOG(TAUH)**2/DLR**2
       if(ism4.ne.0)then
        d0qq=cdabs((cf0t+cf0t4+cf0b4+cf0g)/cf0)**2*32.d0/27.d0
     .                                            *(1.d0-tauh)**3
     .      +dreal(dconjg(cf0t+cf0t4+cf0b4+cf0g)*cf0b)/cdabs(cf0)**2*
     .       128.d0/27.d0*tauh*(1.d0-tauh)**2*dlog(tauh)/dlr
     .      +cdabs(cf0b/cf0)**2*
     .       128.d0/27.d0*tauh**2*(1.d0-tauh)*dlog(tauh)**2/dlr**2
        if(isilh.ne.0)
     .   d0qq=(cdabs(cf0t+cf0t4+cf0b4+cf0g)**2-cf0g**2)
     .       /(cdabs(cf0)**2-cf0g**2)*32.d0/27.d0
     .                                            *(1.d0-tauh)**3
     .      +dreal(dconjg(cf0t+cf0t4+cf0b4+cf0g)*cf0b)
     .           /(cdabs(cf0)**2-cf0g**2)*
     .       128.d0/27.d0*tauh*(1.d0-tauh)**2*dlog(tauh)/dlr
     .      +cdabs(cf0b)**2/(cdabs(cf0)**2-cf0g**2)*
     .       128.d0/27.d0*tauh**2*(1.d0-tauh)*dlog(tauh)**2/dlr**2
       endif
      ELSE
        D0QQ=DCADR2(DFUN,DO,UP,AERR,RERR,ERR,IER)
        IF(IER.GT.66)IFAIL=IFAIL+1
        IF(IER.EQ.66)IF66=IF66+1
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION D0GQ(TAUH)
C--GQ --> HQ
C--MASS DEPENDENT PART: INTEGRATION OVER SCATTERING ANGLE
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/GRENZ/IGRENZ
      COMMON/INTGRLP/TAUH1,IFUN
      COMMON/DCADRE/DO,UP,AERR,RERR
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/HIGGS/IHIGGS
      COMMON/MASS/AMH,AMQ,S
      COMMON/SUSY/AMB,AMC,FACT,FACB,FACC,FACCG,FACCTG,ISILH
      COMMON/BORN/CF0,CF0T,CF0B,CF0C,CF0G
      COMMON/BORN4/CF0T4,CF0B4
      COMMON/SM4/AMT4,AMB4,FACT4,FACB4,ISM4,IGGELW
      EXTERNAL DFUN
      PI=4.D0*DATAN(1.D0)
      TAUH1=TAUH
      IFUN=2
      ICALL1=ICALL1+1
      IF(IGRENZ.EQ.1)THEN
        DLR=DLOG(AMH**2/AMB**2)
C--LIMIT OF HEAVY TOP MASS AND SMALL BOTTOM MASS
        D0GQ=CDABS((CF0T+CF0G)/CF0)**2*
     .       (-1.D0+2.D0*TAUH-TAUH**2/3.D0)
     .      +DREAL(DCONJG(CF0T+CF0G)*CF0B)/CDABS(CF0)**2*
     .       2.D0/3.D0*(TAUH**2+(1.D0+(1.D0-TAUH)**2)
     .      *(-2.D0/3.D0*DLR-2.D0*DLOG((1.D0-TAUH)/TAUH)))
     .      +CDABS(CF0B/CF0)**2*
     .       2.D0/3.D0*(TAUH**2+(1.D0+(1.D0-TAUH)**2)
     .      *(-7.D0/15.D0*DLR-DLOG((1.D0-TAUH)/TAUH)))
        IF(ISILH.NE.0)
     .  D0GQ=(CDABS(CF0T+CF0G)**2-CF0G**2)/(CDABS(CF0)**2-CF0G**2)*
     .       (-1.D0+2.D0*TAUH-TAUH**2/3.D0)
     .      +DREAL(DCONJG(CF0T+CF0G)*CF0B)/(CDABS(CF0)**2-CF0G**2)*
     .       2.D0/3.D0*(TAUH**2+(1.D0+(1.D0-TAUH)**2)
     .      *(-2.D0/3.D0*DLR-2.D0*DLOG((1.D0-TAUH)/TAUH)))
     .      +CDABS(CF0B)**2/(CDABS(CF0)**2-CF0G**2)*
     .       2.D0/3.D0*(TAUH**2+(1.D0+(1.D0-TAUH)**2)
     .      *(-7.D0/15.D0*DLR-DLOG((1.D0-TAUH)/TAUH)))
       if(ism4.ne.0)then
        d0gq=cdabs((cf0t+cf0t4+cf0b4+cf0g)/cf0)**2*
     .       (-1.d0+2.d0*tauh-tauh**2/3.d0)
     .      +dreal(dconjg(cf0t+cf0t4+cf0b4+cf0g)*cf0b)/cdabs(cf0)**2*
     .       2.d0/3.d0*(tauh**2+(1.d0+(1.d0-tauh)**2)
     .      *(-2.d0/3.d0*dlr-2.d0*dlog((1.d0-tauh)/tauh)))
     .      +cdabs(cf0b/cf0)**2*
     .       2.d0/3.d0*(tauh**2+(1.d0+(1.d0-tauh)**2)
     .      *(-7.d0/15.d0*dlr-dlog((1.d0-tauh)/tauh)))
        if(isilh.ne.0)
     .   d0gq=(cdabs(cf0t+cf0t4+cf0b4+cf0g)**2-cf0g**2)
     .       /(cdabs(cf0)**2-cf0g**2)*
     .       (-1.d0+2.d0*tauh-tauh**2/3.d0)
     .      +dreal(dconjg(cf0t+cf0t4+cf0b4+cf0g)*cf0b)
     .           /(cdabs(cf0)**2-cf0g**2)*
     .       2.d0/3.d0*(tauh**2+(1.d0+(1.d0-tauh)**2)
     .      *(-2.d0/3.d0*dlr-2.d0*dlog((1.d0-tauh)/tauh)))
     .      +cdabs(cf0b)**2/(cdabs(cf0)**2-cf0g**2)*
     .       2.d0/3.d0*(tauh**2+(1.d0+(1.d0-tauh)**2)
     .      *(-7.d0/15.d0*dlr-dlog((1.d0-tauh)/tauh)))
       endif
      ELSE
        D0GQ=DCADR2(DFUN,DO,UP,AERR,RERR,ERR,IER)
        IF(IER.GT.66)IFAIL=IFAIL+1
        IF(IER.EQ.66)IF66=IF66+1
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION D0GG(TAUH)
C--GG --> HG
C--MASS DEPENDENT PART: INTEGRATION OVER SCATTERING ANGLE
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/GRENZ/IGRENZ
      COMMON/INTGRLP/TAUH1,IFUN
      COMMON/DCADRE/DO,UP,AERR,RERR
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/HIGGS/IHIGGS
      COMMON/MASS/AMH,AMQ,S
      COMMON/SUSY/AMB,AMC,FACT,FACB,FACC,FACCG,FACCTG,ISILH
      COMMON/BORN/CF0,CF0T,CF0B,CF0C,CF0G
      COMMON/BORN4/CF0T4,CF0B4
      COMMON/SM4/AMT4,AMB4,FACT4,FACB4,ISM4,IGGELW
      EXTERNAL DFUN
      PI=4.D0*DATAN(1.D0)
      TAUH1=TAUH
      IFUN=1
      ICALL1=ICALL1+1
      IF(IGRENZ.EQ.1)THEN
        RHO=AMH**2/AMB**2
        DLGR=DLOG(RHO)
        DLG1T=DLOG(1.D0-TAUH)
        DLGDT=DLOG(TAUH)/(1.D0-TAUH)
C--LIMIT OF HEAVY TOP MASS AND SMALL BOTTOM MASS
        D0GG=CDABS((CF0T+CF0G)/CF0)**2*
     .       (-11.D0/2.D0*(1.D0-TAUH)**3)
     .      +DREAL(DCONJG(CF0T+CF0G)*CF0B)/CDABS(CF0)**2*
     .         (2.D0*DLGR*( - TAUH**2 + 2.D0*TAUH - 2.D0)
     .         + 6.D0*DLG1T*( - TAUH**2 + 2.D0*TAUH - 2.D0)
     .         + 12.D0*DLGDT*( - TAUH**4 + TAUH**3
     .                         - 2.D0*TAUH + 1.D0))
     .      +CDABS(CF0B/CF0)**2*
     .         (2.D0*DLGR*( - 5.D0*TAUH**2 + 7.D0*TAUH - 7.D0)
     .         + 30.D0*DLG1T*( - TAUH**2 + TAUH - 1.D0)
     .         + 10.D0*DLGDT*( - 6.D0*TAUH**4 + 5.D0*TAUH**3
     .                         - 2.D0*TAUH**2 - 6.D0*TAUH + 3.D0)
     .         )/5.D0
        IF(ISILH.NE.0)
     .  D0GG=(CDABS(CF0T+CF0G)**2-CF0G**2)/(CDABS(CF0)**2-CF0G**2)*
     .       (-11.D0/2.D0*(1.D0-TAUH)**3)
     .      +DREAL(DCONJG(CF0T+CF0G)*CF0B)/(CDABS(CF0)**2-CF0G**2)*
     .         (2.D0*DLGR*( - TAUH**2 + 2.D0*TAUH - 2.D0)
     .         + 6.D0*DLG1T*( - TAUH**2 + 2.D0*TAUH - 2.D0)
     .         + 12.D0*DLGDT*( - TAUH**4 + TAUH**3
     .                         - 2.D0*TAUH + 1.D0))
     .      +CDABS(CF0B)**2/(CDABS(CF0)**2-CF0G**2)*
     .         (2.D0*DLGR*( - 5.D0*TAUH**2 + 7.D0*TAUH - 7.D0)
     .         + 30.D0*DLG1T*( - TAUH**2 + TAUH - 1.D0)
     .         + 10.D0*DLGDT*( - 6.D0*TAUH**4 + 5.D0*TAUH**3
     .                         - 2.D0*TAUH**2 - 6.D0*TAUH + 3.D0)
     .         )/5.D0
       if(ism4.ne.0)then
        d0gg=cdabs((cf0t+cf0t4+cf0b4+cf0g)/cf0)**2*
     .       (-11.d0/2.d0*(1.d0-tauh)**3)
     .      +dreal(dconjg(cf0t+cf0t4+cf0b4+cf0g)*cf0b)/cdabs(cf0)**2*
     .         (2.d0*dlgr*( - tauh**2 + 2.d0*tauh - 2.d0)
     .         + 6.d0*dlg1t*( - tauh**2 + 2.d0*tauh - 2.d0)
     .         + 12.d0*dlgdt*( - tauh**4 + tauh**3
     .                         - 2.d0*tauh + 1.d0))
     .      +cdabs(cf0b/cf0)**2*
     .         (2.d0*dlgr*( - 5.d0*tauh**2 + 7.d0*tauh - 7.d0)
     .         + 30.d0*dlg1t*( - tauh**2 + tauh - 1.d0)
     .         + 10.d0*dlgdt*( - 6.d0*tauh**4 + 5.d0*tauh**3
     .                         - 2.d0*tauh**2 - 6.d0*tauh + 3.d0)
     .         )/5.d0
        if(isilh.ne.0)
     .  d0gg=(cdabs(cf0t+cf0t4+cf0b4+cf0g)**2-cf0g**2)
     .       /(cdabs(cf0)**2-cf0g**2)*
     .       (-11.d0/2.d0*(1.d0-tauh)**3)
     .      +dreal(dconjg(cf0t+cf0t4+cf0b4+cf0g)*cf0b)
     .           /(cdabs(cf0)**2-cf0g**2)*
     .         (2.d0*dlgr*( - tauh**2 + 2.d0*tauh - 2.d0)
     .         + 6.d0*dlg1t*( - tauh**2 + 2.d0*tauh - 2.d0)
     .         + 12.d0*dlgdt*( - tauh**4 + tauh**3
     .                         - 2.d0*tauh + 1.d0))
     .      +cdabs(cf0b)**2/(cdabs(cf0)**2-cf0g**2)*
     .         (2.d0*dlgr*( - 5.d0*tauh**2 + 7.d0*tauh - 7.d0)
     .         + 30.d0*dlg1t*( - tauh**2 + tauh - 1.d0)
     .         + 10.d0*dlgdt*( - 6.d0*tauh**4 + 5.d0*tauh**3
     .                         - 2.d0*tauh**2 - 6.d0*tauh + 3.d0)
     .         )/5.d0
       endif
      ELSE
        D0GG=DCADR2(DFUN,DO,UP,AERR,RERR,ERR,IER)
        IF(IER.GT.66)IFAIL=IFAIL+1
        IF(IER.EQ.66)IF66=IF66+1
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION EVGG(X)
C--H --> GGG: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(2)
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/INTRO/INTEG
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      ICALL1=ICALL1+1
      X1=EPST+(1.D0-2.D0*EPST)*X(1)
      X2=EPST+(1.D0-2.D0*EPST)*X(2)
      VGG=EIGG(X1,X2)
      EVGG=VGG*(1.D0-2.D0*EPST)**2
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION EVGQ(X)
C--H --> GQQBAR: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(2)
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/INTRO/INTEG
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      ICALL1=ICALL1+1
      X1=EPST+(1.D0-2.D0*EPST)*X(1)
      X2=EPST+(1.D0-2.D0*EPST)*X(2)
      VGQ=EIGQ(X1,X2)
      EVGQ=VGQ*(1.D0-2.D0*EPST)**2
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION EVQQ(X)
C--H --> GTTBAR: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(2)
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/INTRO/INTEG
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      ICALL1=ICALL1+1
      X1=EPST+(1.D0-2.D0*EPST)*X(1)
      X2=EPST+(1.D0-2.D0*EPST)*X(2)
      VQQ=EIQQ(X1,X2)
      EVQQ=VQQ*(1.D0-2.D0*EPST)**2
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION EGG(TAUH)
C--H --> GGG: FIRST INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/INTGRLP/TAUH1,IFUN
      COMMON/DCADRE/DO,UP,AERR,RERR
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/HIGGS/IHIGGS
      EXTERNAL DFUN
      TAUH1=TAUH
      IFUN=1
      ICALL1=ICALL1+1
      EGG=DCADR2(DFUN,DO,UP,AERR,RERR,ERR,IER)
      IF(IER.GT.66)IFAIL=IFAIL+1
      IF(IER.EQ.66)IF66=IF66+1
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION EGQ(TAUH)
C--H --> GQQBAR: FIRST INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/INTGRLP/TAUH1,IFUN
      COMMON/DCADRE/DO,UP,AERR,RERR
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/HIGGS/IHIGGS
      EXTERNAL DFUN
      TAUH1=TAUH
      IFUN=2
      ICALL1=ICALL1+1
      EGQ=DCADR2(DFUN,DO,UP,AERR,RERR,ERR,IER)
      IF(IER.GT.66)IFAIL=IFAIL+1
      IF(IER.EQ.66)IF66=IF66+1
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION EQQ(TAUH)
C--H --> GTTBAR: FIRST INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/INTGRLP/TAUH1,IFUN
      COMMON/DCADRE/DO,UP,AERR,RERR
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/HIGGS/IHIGGS
      EXTERNAL DFUN
      TAUH1=TAUH
      IFUN=3
      ICALL1=ICALL1+1
      EQQ=DCADR2(DFUN,DO,UP,AERR,RERR,ERR,IER)
      IF(IER.GT.66)IFAIL=IFAIL+1
      IF(IER.EQ.66)IF66=IF66+1
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DFUN(V)
C--FUNCTION FOR FIRST INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/PROCESS/IPROC
      COMMON/INTGRLP/TAUH,IFUN
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      ICALL2=ICALL2+1
      IF(IPROC.EQ.1)THEN
       IF(IFUN.EQ.1)THEN
         DFUN=DIGG(TAUH,V)
       ELSEIF(IFUN.EQ.2)THEN
         DFUN=DIGQ(TAUH,V)
       ELSE
         DFUN=DIQQ(TAUH,V)
       ENDIF
      ELSE
       IF(IFUN.EQ.1)THEN
         DFUN=EIGG(TAUH,V)
       ELSEIF(IFUN.EQ.2)THEN
         DFUN=EIGQ(TAUH,V)
       ELSEIF(IFUN.EQ.3)THEN
         DFUN=EIQQ(TAUH,V)
       ENDIF
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DIQQ(TAUH,V)
C--QQBAR --> HG: MASS-DEPENDENT INTEGRAND
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/HIGGS/IHIGGS
      COMMON/MASS/AMH,AMQ,S
      COMMON/SUSY/AMB,AMC,FACT,FACB,FACC,FACCG,FACCTG,ISILH
      COMMON/SUSY0/AMST1,AMST2,AMSB1,AMSB2,FACST1,FACST2,FACSB1,FACSB2
      COMMON/GRENZSQ/IGRENZSQ
      COMMON/SM4/AMT4,AMB4,FACT4,FACB4,ISM4,IGGELW
      COMMON/BORN/CCF0,CCF0T,CCF0B,CCF0C,CF0G
      COMMON/BORN00/CF0G00
C--MANDELSTAM-VARIABLES
      SS=AMH**2/TAUH
      TT=-SS*(1.D0-TAUH)*V
      UU=-SS*(1.D0-TAUH)*(1.D0-V)
      CSQQ=DCMPLX(0.D0,0.D0)
      CSQQP=DCMPLX(0.D0,0.D0)
C--FORM FACTOR
      IF(FACT.NE.0.D0)THEN
        CSQQ=CSQQ+FACT*CIGQ(TT,SS,UU,AMQ)
        CSQQP=CSQQP+(FACT-1)*CIGQ(TT,SS,UU,AMQ)
        IF(FACST1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          CSQQ=CSQQ+FACST1*CLGQ(TT,SS,UU,AMST1)*CSBORN(AMH,AMST1)
     .             +FACST2*CLGQ(TT,SS,UU,AMST2)*CSBORN(AMH,AMST2)
         ELSE
          CSQQ=CSQQ+FACST1*CIGQS(TT,SS,UU,AMST1)
     .             +FACST2*CIGQS(TT,SS,UU,AMST2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACB.NE.0.D0)THEN
        CSQQ=CSQQ+FACB*CIGQ(TT,SS,UU,AMB)
        CSQQP=CSQQP+(FACB-1)*CIGQ(TT,SS,UU,AMB)
        IF(FACSB1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          CSQQ=CSQQ+FACSB1*CLGQ(TT,SS,UU,AMSB1)*CSBORN(AMH,AMSB1)
     .             +FACSB2*CLGQ(TT,SS,UU,AMSB2)*CSBORN(AMH,AMSB2)
         ELSE
          CSQQ=CSQQ+FACSB1*CIGQS(TT,SS,UU,AMSB1)
     .             +FACSB2*CIGQS(TT,SS,UU,AMSB2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACC.NE.0.D0)THEN
        CSQQ=CSQQ+FACC*CIGQ(TT,SS,UU,AMC)
        CSQQP=CSQQP+(FACC-1)*CIGQ(TT,SS,UU,AMC)
      ENDIF
      IF(FACCG.NE.0.D0)THEN
        CSQQ=CSQQ-4*CF0G00/3
        CSQQP=CSQQP-4*CF0G00/3
      ENDIF
      IF(FACCTG.NE.0.D0)THEN
        CSQQ=CSQQ+FACCTG*CIGQT(TT,SS,UU,AMQ)
        CSQQP=CSQQP+FACCTG*CIGQT(TT,SS,UU,AMQ)
      ENDIF
      if(ism4.ne.0)then
       csqq=csqq+fact4*cigq(tt,ss,uu,amt4)+facb4*cigq(tt,ss,uu,amb4)
       csqqp=csqqp+(fact4-1)*cigq(tt,ss,uu,amt4)
     .            +(facb4-1)*cigq(tt,ss,uu,amb4)
      endif
C--CURRENT FACTOR
      HQQ=(UU**2+TT**2)/SS
C--LOWEST ORDER FORM FACTOR
      IF(IHIGGS.NE.1)THEN
C--SCALAR HIGGS
        CF0=FACT*CFBORN(AMH,AMQ)+FACB*CFBORN(AMH,AMB)
     .     +FACC*CFBORN(AMH,AMC)+CF0G
     .     +FACST1*CSBORN(AMH,AMST1)+FACST2*CSBORN(AMH,AMST2)
     .     +FACSB1*CSBORN(AMH,AMSB1)+FACSB2*CSBORN(AMH,AMSB2)
        CF0P=(FACT-1)*CFBORN(AMH,AMQ)+(FACB-1)*CFBORN(AMH,AMB)
     .      +(FACC-1)*CFBORN(AMH,AMC)+CF0G
        if(ism4.ne.0)then
         cf0=fact*cfborn(amh,amq)+facb*cfborn(amh,amb)
     .      +facc*cfborn(amh,amc)+cf0g
     .      +fact4*cfborn(amh,amt4)+facb4*cfborn(amh,amb4)
         cf0p=(fact-1)*cfborn(amh,amq)+(facb-1)*cfborn(amh,amb)
     .      +(facc-1)*cfborn(amh,amc)+cf0g
     .      +(fact4-1)*cfborn(amh,amt4)+(facb4-1)*cfborn(amh,amb4)
        endif
      ELSE
C--PSEUDOSCALAR HIGGS
        CSQQP = 0
        CF0P = 0
        CF0=FACT*CFBORNA(AMH,AMQ)+FACB*CFBORNA(AMH,AMB)
     .     +FACC*CFBORNA(AMH,AMC)
        if(ism4.ne.0)then
         cf0=fact*cfborna(amh,amq)+facb*cfborna(amh,amb)
     .      +facc*cfborna(amh,amc)
     .      +fact4*cfborna(amh,amt4)+facb4*cfborna(amh,amb4)
        endif
      ENDIF
      DIQQ=CDABS(CSQQ)**2*HQQ/CDABS(CF0)**2
     .     *(1.D0-TAUH)/SS
      IF(ISILH.NE.0) DIQQ=(CDABS(CSQQ)**2-CDABS(CSQQP)**2)*HQQ
     .                   / (CDABS(CF0)**2-CDABS(CF0P)**2)
     .                   *(1.D0-TAUH)/SS
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DIGQ(TAUH,V)
C--GQ --> HQ: MASS-DEPENDENT INTEGRAND
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/HIGGS/IHIGGS
      COMMON/MASS/AMH,AMQ,S
      COMMON/SUSY/AMB,AMC,FACT,FACB,FACC,FACCG,FACCTG,ISILH
      COMMON/SUSY0/AMST1,AMST2,AMSB1,AMSB2,FACST1,FACST2,FACSB1,FACSB2
      COMMON/GRENZSQ/IGRENZSQ
      COMMON/SM4/AMT4,AMB4,FACT4,FACB4,ISM4,IGGELW
      COMMON/BORN/CCF0,CCF0T,CCF0B,CCF0C,CF0G
      COMMON/BORN00/CF0G00
C--MANDELSTAM-VARIABLES
      SS=AMH**2/TAUH
      TT=-SS*(1.D0-TAUH)*V
      UU=-SS*(1.D0-TAUH)*(1.D0-V)
      CSGQ=DCMPLX(0.D0,0.D0)
      CSGQP=DCMPLX(0.D0,0.D0)
C--FORM FACTOR
      IF(FACT.NE.0.D0)THEN
        CSGQ=CSGQ+FACT*CIGQ(SS,TT,UU,AMQ)
        CSGQP=CSGQP+(FACT-1)*CIGQ(SS,TT,UU,AMQ)
        IF(FACST1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          CSGQ=CSGQ+FACST1*CLGQ(SS,TT,UU,AMST1)*CSBORN(AMH,AMST1)
     .             +FACST2*CLGQ(SS,TT,UU,AMST2)*CSBORN(AMH,AMST2)
         ELSE
          CSGQ=CSGQ+FACST1*CIGQS(SS,TT,UU,AMST1)
     .             +FACST2*CIGQS(SS,TT,UU,AMST2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACB.NE.0.D0)THEN
        CSGQ=CSGQ+FACB*CIGQ(SS,TT,UU,AMB)
        CSGQP=CSGQP+(FACB-1)*CIGQ(SS,TT,UU,AMB)
        IF(FACSB1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          CSGQ=CSGQ+FACSB1*CLGQ(SS,TT,UU,AMSB1)*CSBORN(AMH,AMSB1)
     .             +FACSB2*CLGQ(SS,TT,UU,AMSB2)*CSBORN(AMH,AMSB2)
         ELSE
          CSGQ=CSGQ+FACSB1*CIGQS(SS,TT,UU,AMSB1)
     .             +FACSB2*CIGQS(SS,TT,UU,AMSB2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACC.NE.0.D0)THEN
        CSGQ=CSGQ+FACC*CIGQ(SS,TT,UU,AMC)
        CSGQP=CSGQP+(FACC-1)*CIGQ(SS,TT,UU,AMC)
      ENDIF
      IF(FACCG.NE.0.D0)THEN
        CSGQ=CSGQ-4*CF0G00/3
        CSGQP=CSGQP-4*CF0G00/3
      ENDIF
      IF(FACCTG.NE.0.D0)THEN
        CSGQ=CSGQ+FACCTG*CIGQT(SS,TT,UU,AMQ)
        CSGQP=CSGQP+FACCTG*CIGQT(SS,TT,UU,AMQ)
      ENDIF
      if(ism4.ne.0)then
       csgq=csgq+fact4*cigq(ss,tt,uu,amt4)+facb4*cigq(ss,tt,uu,amb4)
       csgqp=csgqp+(fact4-1)*cigq(ss,tt,uu,amt4)
     .            +(facb4-1)*cigq(ss,tt,uu,amb4)
      endif
C--CURRENT FACTOR
      HGQ=-(3.D0*(UU**2+SS**2))/(8.D0*TT)
C--LOWEST ORDER FORM FACTOR
      IF(IHIGGS.NE.1)THEN
C--SCALAR HIGGS
        CF0=FACT*CFBORN(AMH,AMQ)+FACB*CFBORN(AMH,AMB)
     .     +FACC*CFBORN(AMH,AMC)+CF0G
     .     +FACST1*CSBORN(AMH,AMST1)+FACST2*CSBORN(AMH,AMST2)
     .     +FACSB1*CSBORN(AMH,AMSB1)+FACSB2*CSBORN(AMH,AMSB2)
        CF0P=(FACT-1)*CFBORN(AMH,AMQ)+(FACB-1)*CFBORN(AMH,AMB)
     .      +(FACC-1)*CFBORN(AMH,AMC)+CF0G
        if(ism4.ne.0)then
         cf0=fact*cfborn(amh,amq)+facb*cfborn(amh,amb)
     .      +facc*cfborn(amh,amc)+cf0g
     .      +fact4*cfborn(amh,amt4)+facb4*cfborn(amh,amb4)
         cf0p=(fact-1)*cfborn(amh,amq)+(facb-1)*cfborn(amh,amb)
     .      +(facc-1)*cfborn(amh,amc)+cf0g
     .      +(fact4-1)*cfborn(amh,amt4)+(facb4-1)*cfborn(amh,amb4)
        endif
      ELSE
C--PSEUDOSCALAR HIGGS
        CSGQP = 0
        CF0P = 0
        CF0=FACT*CFBORNA(AMH,AMQ)+FACB*CFBORNA(AMH,AMB)
     .     +FACC*CFBORNA(AMH,AMC)
        if(ism4.ne.0)then
         cf0=fact*cfborna(amh,amq)+facb*cfborna(amh,amb)
     .      +facc*cfborna(amh,amc)
     .      +fact4*cfborna(amh,amt4)+facb4*cfborna(amh,amb4)
        endif
      ENDIF
      DADD=2.D0/3.D0*TAUH**2
      DIGQ=CDABS(CSGQ)**2*HGQ/CDABS(CF0)**2
     .     *(1.D0-TAUH)/SS+DADD
     .     -2.D0/3.D0*(1.D0+(1.D0-TAUH)**2)/V
      IF(ISILH.NE.0) DIGQ=(CDABS(CSGQ)**2-CDABS(CSGQP)**2)*HGQ
     .                   / (CDABS(CF0)**2-CDABS(CF0P)**2)
     .     *(1.D0-TAUH)/SS+DADD
     .     -2.D0/3.D0*(1.D0+(1.D0-TAUH)**2)/V
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DIGG(TAUH,V)
C--GG --> HG: MASS-DEPENDENT INTEGRAND
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/HIGGS/IHIGGS
      COMMON/MASS/AMH,AMQ,S
      COMMON/SUSY/AMB,AMC,FACT,FACB,FACC,FACCG,FACCTG,ISILH
      COMMON/SUSY0/AMST1,AMST2,AMSB1,AMSB2,FACST1,FACST2,FACSB1,FACSB2
      COMMON/GRENZSQ/IGRENZSQ
      COMMON/SM4/AMT4,AMB4,FACT4,FACB4,ISM4,IGGELW
      COMMON/BORN/CCF0,CCF0T,CCF0B,CCF0C,CF0G
      COMMON/BORN00/CF0G00
C--MANDELSTAM-VARIABLES
      SS=AMH**2/TAUH
      TT=-SS*(1.D0-TAUH)*V
      UU=-SS*(1.D0-TAUH)*(1.D0-V)
      C5=DCMPLX(0.D0,0.D0)
      C2=DCMPLX(0.D0,0.D0)
      C3=DCMPLX(0.D0,0.D0)
      C4=DCMPLX(0.D0,0.D0)
      C5P=DCMPLX(0.D0,0.D0)
      C2P=DCMPLX(0.D0,0.D0)
      C3P=DCMPLX(0.D0,0.D0)
      C4P=DCMPLX(0.D0,0.D0)
C--FORM FACTOR
      IF(FACT.NE.0.D0)THEN
        C5=C5+FACT*CFGG5(SS,TT,UU,AMQ)
        C2=C2+FACT*CFGG2(SS,TT,UU,AMQ)
        C3=C3+FACT*CFGG2(TT,SS,UU,AMQ)
        C4=C4+FACT*CFGG2(UU,TT,SS,AMQ)
        C5P=C5P+(FACT-1)*CFGG5(SS,TT,UU,AMQ)
        C2P=C2P+(FACT-1)*CFGG2(SS,TT,UU,AMQ)
        C3P=C3P+(FACT-1)*CFGG2(TT,SS,UU,AMQ)
        C4P=C4P+(FACT-1)*CFGG2(UU,TT,SS,AMQ)
        IF(FACST1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          C5=C5+FACST1*CFGL5(SS,TT,UU,AMST1)*CSBORN(AMH,AMST1)
     .         +FACST2*CFGL5(SS,TT,UU,AMST2)*CSBORN(AMH,AMST2)
          C2=C2+FACST1*CFGL2(SS,TT,UU,AMST1)*CSBORN(AMH,AMST1)
     .         +FACST2*CFGL2(SS,TT,UU,AMST2)*CSBORN(AMH,AMST2)
          C3=C3+FACST1*CFGL2(TT,SS,UU,AMST1)*CSBORN(AMH,AMST1)
     .         +FACST2*CFGL2(TT,SS,UU,AMST2)*CSBORN(AMH,AMST2)
          C4=C4+FACST1*CFGL2(UU,TT,SS,AMST1)*CSBORN(AMH,AMST1)
     .         +FACST2*CFGL2(UU,TT,SS,AMST2)*CSBORN(AMH,AMST2)
         ELSE
          C5=C5+FACST1*CFGG5S(SS,TT,UU,AMST1)
     .         +FACST2*CFGG5S(SS,TT,UU,AMST2)
          C2=C2+FACST1*CFGG2S(SS,TT,UU,AMST1)
     .         +FACST2*CFGG2S(SS,TT,UU,AMST2)
          C3=C3+FACST1*CFGG2S(TT,SS,UU,AMST1)
     .         +FACST2*CFGG2S(TT,SS,UU,AMST2)
          C4=C4+FACST1*CFGG2S(UU,TT,SS,AMST1)
     .         +FACST2*CFGG2S(UU,TT,SS,AMST2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACB.NE.0.D0)THEN
        C5=C5+FACB*CFGG5(SS,TT,UU,AMB)
        C2=C2+FACB*CFGG2(SS,TT,UU,AMB)
        C3=C3+FACB*CFGG2(TT,SS,UU,AMB)
        C4=C4+FACB*CFGG2(UU,TT,SS,AMB)
        C5P=C5P+(FACB-1)*CFGG5(SS,TT,UU,AMB)
        C2P=C2P+(FACB-1)*CFGG2(SS,TT,UU,AMB)
        C3P=C3P+(FACB-1)*CFGG2(TT,SS,UU,AMB)
        C4P=C4P+(FACB-1)*CFGG2(UU,TT,SS,AMB)
        IF(FACSB1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          C5=C5+FACSB1*CFGL5(SS,TT,UU,AMSB1)*CSBORN(AMH,AMSB1)
     .         +FACSB2*CFGL5(SS,TT,UU,AMSB2)*CSBORN(AMH,AMSB2)
          C2=C2+FACSB1*CFGL2(SS,TT,UU,AMSB1)*CSBORN(AMH,AMSB1)
     .         +FACSB2*CFGL2(SS,TT,UU,AMSB2)*CSBORN(AMH,AMSB2)
          C3=C3+FACSB1*CFGL2(TT,SS,UU,AMSB1)*CSBORN(AMH,AMSB1)
     .         +FACSB2*CFGL2(TT,SS,UU,AMSB2)*CSBORN(AMH,AMSB2)
          C4=C4+FACSB1*CFGL2(UU,TT,SS,AMSB1)*CSBORN(AMH,AMSB1)
     .         +FACSB2*CFGL2(UU,TT,SS,AMSB2)*CSBORN(AMH,AMSB2)
         ELSE
          C5=C5+FACSB1*CFGG5S(SS,TT,UU,AMSB1)
     .         +FACSB2*CFGG5S(SS,TT,UU,AMSB2)
          C2=C2+FACSB1*CFGG2S(SS,TT,UU,AMSB1)
     .         +FACSB2*CFGG2S(SS,TT,UU,AMSB2)
          C3=C3+FACSB1*CFGG2S(TT,SS,UU,AMSB1)
     .         +FACSB2*CFGG2S(TT,SS,UU,AMSB2)
          C4=C4+FACSB1*CFGG2S(UU,TT,SS,AMSB1)
     .         +FACSB2*CFGG2S(UU,TT,SS,AMSB2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACC.NE.0.D0)THEN
        C5=C5+FACC*CFGG5(SS,TT,UU,AMC)
        C2=C2+FACC*CFGG2(SS,TT,UU,AMC)
        C3=C3+FACC*CFGG2(TT,SS,UU,AMC)
        C4=C4+FACC*CFGG2(UU,TT,SS,AMC)
        C5P=C5P+(FACC-1)*CFGG5(SS,TT,UU,AMC)
        C2P=C2P+(FACC-1)*CFGG2(SS,TT,UU,AMC)
        C3P=C3P+(FACC-1)*CFGG2(TT,SS,UU,AMC)
        C4P=C4P+(FACC-1)*CFGG2(UU,TT,SS,AMC)
      ENDIF
      IF(FACCG.NE.0.D0)THEN
        C5=C5+CF0G00*(SS+TT+UU)**2
        C2=C2+CF0G00*SS**2
        C3=C3+CF0G00*TT**2
        C4=C4+CF0G00*UU**2
        C5P=C5P+CF0G00*(SS+TT+UU)**2
        C2P=C2P+CF0G00*SS**2
        C3P=C3P+CF0G00*TT**2
        C4P=C4P+CF0G00*UU**2
      ENDIF
      IF(FACCTG.NE.0.D0)THEN
        C5=C5+FACCTG*CFGG5T(SS,TT,UU,AMQ)
        C2=C2+FACCTG*CFGG2T(SS,TT,UU,AMQ)
        C3=C3+FACCTG*CFGG2T(TT,SS,UU,AMQ)
        C4=C4+FACCTG*CFGG2T(UU,TT,SS,AMQ)
        C5P=C5P+FACCTG*CFGG5T(SS,TT,UU,AMQ)
        C2P=C2P+FACCTG*CFGG2T(SS,TT,UU,AMQ)
        C3P=C3P+FACCTG*CFGG2T(TT,SS,UU,AMQ)
        C4P=C4P+FACCTG*CFGG2T(UU,TT,SS,AMQ)
      ENDIF
      if(ism4.ne.0)then
       C5=C5+fact4*CFGG5(SS,TT,UU,AMT4)
       C2=C2+fact4*CFGG2(SS,TT,UU,AMT4)
       C3=C3+fact4*CFGG2(TT,SS,UU,AMT4)
       C4=C4+fact4*CFGG2(UU,TT,SS,AMT4)
       C5=C5+facb4*CFGG5(SS,TT,UU,AMB4)
       C2=C2+facb4*CFGG2(SS,TT,UU,AMB4)
       C3=C3+facb4*CFGG2(TT,SS,UU,AMB4)
       C4=C4+facb4*CFGG2(UU,TT,SS,AMB4)
       C5P=C5P+(fact4-1)*CFGG5(SS,TT,UU,AMT4)
       C2P=C2P+(fact4-1)*CFGG2(SS,TT,UU,AMT4)
       C3P=C3P+(fact4-1)*CFGG2(TT,SS,UU,AMT4)
       C4P=C4P+(fact4-1)*CFGG2(UU,TT,SS,AMT4)
       C5P=C5P+(facb4-1)*CFGG5(SS,TT,UU,AMB4)
       C2P=C2P+(facb4-1)*CFGG2(SS,TT,UU,AMB4)
       C3P=C3P+(facb4-1)*CFGG2(TT,SS,UU,AMB4)
       C4P=C4P+(facb4-1)*CFGG2(UU,TT,SS,AMB4)
      endif
      FGG=CDABS(C5)**2+CDABS(C2)**2+CDABS(C3)**2+CDABS(C4)**2
      FGGP=CDABS(C5P)**2+CDABS(C2P)**2+CDABS(C3P)**2+CDABS(C4P)**2
      IF(ISILH.NE.0)FGG = FGG - FGGP
      FGG=FGG/SS**4
C--LOWEST ORDER FORM FACTOR
      IF(IHIGGS.NE.1)THEN
C--SCALAR HIGGS
        CF0=FACT*CFBORN(AMH,AMQ)+FACB*CFBORN(AMH,AMB)
     .     +FACC*CFBORN(AMH,AMC)+CF0G
     .     +FACST1*CSBORN(AMH,AMST1)+FACST2*CSBORN(AMH,AMST2)
     .     +FACSB1*CSBORN(AMH,AMSB1)+FACSB2*CSBORN(AMH,AMSB2)
        CF0P=(FACT-1)*CFBORN(AMH,AMQ)+(FACB-1)*CFBORN(AMH,AMB)
     .      +(FACC-1)*CFBORN(AMH,AMC)+CF0G
        if(ism4.ne.0)then
         cf0=fact*cfborn(amh,amq)+facb*cfborn(amh,amb)
     .      +facc*cfborn(amh,amc)+cf0g
     .      +fact4*cfborn(amh,amt4)+facb4*cfborn(amh,amb4)
         cf0p=(fact-1)*cfborn(amh,amq)+(facb-1)*cfborn(amh,amb)
     .      +(facc-1)*cfborn(amh,amc)+cf0g
     .      +(fact4-1)*cfborn(amh,amt4)+(facb4-1)*cfborn(amh,amb4)
        endif
      ELSE
C--PSEUDOSCALAR HIGGS
        FGGP = 0
        CF0P = 0
        CF0=FACT*CFBORNA(AMH,AMQ)+FACB*CFBORNA(AMH,AMB)
     .     +FACC*CFBORNA(AMH,AMC)
        if(ism4.ne.0)then
         cf0=fact*cfborna(amh,amq)+facb*cfborna(amh,amb)
     .      +facc*cfborna(amh,amc)
     .      +fact4*cfborna(amh,amt4)+facb4*cfborna(amh,amb4)
        endif
      ENDIF
      DIGG=3.D0*(FGG/CDABS(CF0)**2
     .          -1.D0-TAUH**4-(1.D0-TAUH)**4)/(1.D0-TAUH)/V
      IF(ISILH.NE.0) DIGG=3.D0*(FGG/(CDABS(CF0)**2-CDABS(CF0P)**2)
     .          -1.D0-TAUH**4-(1.D0-TAUH)**4)/(1.D0-TAUH)/V
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION EIQQ(X,Y)
C--H --> GTTBAR: MASS-DEPENDENT INTEGRAND
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/HIGGS/IHIGGS
      COMMON/MASS/AMH,AMQ,S
      COMMON/SUSY/AMB,AMC,FACT,FACB,FACC,FACCG,FACCTG,ISILH
      COMMON/SUSY0/AMST1,AMST2,AMSB1,AMSB2,FACST1,FACST2,FACSB1,FACSB2
      COMMON/GRENZSQ/IGRENZSQ
      COMMON/BORN/CCF0,CCF0T,CCF0B,CCF0C,CF0G
      COMMON/BORN00/CF0G00
      IF(AMH.LE.2.D0*AMQ)THEN
       EIQQ=0.D0
       RETURN
      ENDIF
C--MASSIVE THREE-PARTICLE PHASE SPACE
      XMU=2.D0*AMQ/AMH
      X2=1.D0-(1.D0-XMU)*X
      BETA=DSQRT(1.D0-XMU**2/X2**2)
      AA=2.D0*(1.D0-X2)/(2.D0-X2*(1.D0-BETA))
      BB=4.D0*(1.D0-X2)*X2*BETA/((2.D0-X2)**2-X2**2*BETA**2)
      X3=AA+BB*Y
      X1=2.D0-X2-X3
C--MANDELSTAM-VARIABLES: S1=SS-MT**2, U1=UU-MT**2
      TT=AMH**2*(1.D0-X3)
      S1=AMH**2*(1.D0-X2)
      U1=AMH**2*(1.D0-X1)
      CSQQ=DCMPLX(0.D0,0.D0)
      CSQQP=DCMPLX(0.D0,0.D0)
C--FORM FACTOR
      IF(FACT.NE.0.D0)THEN
        CSQQ=CSQQ+FACT*CIGQ(S1,TT,U1,AMQ)
        CSQQP=CSQQP+(FACT-1)*CIGQ(S1,TT,U1,AMQ)
        IF(FACST1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          CSQQ=CSQQ+FACST1*CLGQ(S1,TT,U1,AMST1)*CSBORN(AMH,AMST1)
     .             +FACST2*CLGQ(S1,TT,U1,AMST2)*CSBORN(AMH,AMST2)
         ELSE
          CSQQ=CSQQ+FACST1*CIGQS(S1,TT,U1,AMST1)
     .             +FACST2*CIGQS(S1,TT,U1,AMST2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACB.NE.0.D0)THEN
        CSQQ=CSQQ+FACB*CIGQ(S1,TT,U1,AMB)
        CSQQP=CSQQP+(FACB-1)*CIGQ(S1,TT,U1,AMB)
        IF(FACSB1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          CSQQ=CSQQ+FACSB1*CLGQ(S1,TT,U1,AMSB1)*CSBORN(AMH,AMSB1)
     .             +FACSB2*CLGQ(S1,TT,U1,AMSB2)*CSBORN(AMH,AMSB2)
         ELSE
          CSQQ=CSQQ+FACSB1*CIGQS(S1,TT,U1,AMSB1)
     .             +FACSB2*CIGQS(S1,TT,U1,AMSB2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACC.NE.0.D0)THEN
        CSQQ=CSQQ+FACC*CIGQ(S1,TT,U1,AMC)
        CSQQP=CSQQP+(FACC-1)*CIGQ(S1,TT,U1,AMC)
      ENDIF
      IF(FACCG.NE.0.D0)THEN
        CSQQ=CSQQ-4*CF0G00/3
        CSQQP=CSQQP-4*CF0G00/3
      ENDIF
      IF(FACCTG.NE.0.D0)THEN
        CSQQ=CSQQ+FACCTG*CIGQT(S1,TT,U1,AMQ)
        CSQQP=CSQQP+FACCTG*CIGQT(S1,TT,U1,AMQ)
      ENDIF
C--CURRENT FACTOR
      HQQ=3.D0*(U1**2+S1**2+2.D0*AMQ**2/TT*(S1+U1)**2)/(8.D0*TT)
C--LOWEST ORDER FORM FACTOR
      IF(IHIGGS.NE.1)THEN
C--SCALAR HIGGS
        CF0=FACT*CFBORN(AMH,AMQ)+FACB*CFBORN(AMH,AMB)
     .     +FACC*CFBORN(AMH,AMC)+CF0G
     .     +FACST1*CSBORN(AMH,AMST1)+FACST2*CSBORN(AMH,AMST2)
     .     +FACSB1*CSBORN(AMH,AMSB1)+FACSB2*CSBORN(AMH,AMSB2)
        CF0P=(FACT-1)*CFBORN(AMH,AMQ)+(FACB-1)*CFBORN(AMH,AMB)
     .      +(FACC-1)*CFBORN(AMH,AMC)+CF0G
      ELSE
C--PSEUDOSCALAR HIGGS
        CSQQP=0
        CF0P=0
        CF0=FACT*CFBORNA(AMH,AMQ)+FACB*CFBORNA(AMH,AMB)
     .     +FACC*CFBORNA(AMH,AMC)
      ENDIF
      EIQQ=CDABS(CSQQ)**2*HQQ/CDABS(CF0)**2
      IF(ISILH.NE.0) EIQQ=(CDABS(CSQQ)**2-CDABS(CSQQP)**2)*HQQ
     .                   / (CDABS(CF0)**2-CDABS(CF0P)**2)
      EIQQ=3.D0/4.D0/AMH**2*EIQQ*(1.D0-XMU)*BB
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION EIGQ(X,Y)
C--H --> GQQBAR: MASS-DEPENDENT INTEGRAND
C--DIFFRENCE TO HEAVY TOP LIMIT
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/HIGGS/IHIGGS
      COMMON/MASS/AMH,AMQ,S
      COMMON/SUSY/AMB,AMC,FACT,FACB,FACC,FACCG,FACCTG,ISILH
      COMMON/SUSY0/AMST1,AMST2,AMSB1,AMSB2,FACST1,FACST2,FACSB1,FACSB2
      COMMON/GRENZSQ/IGRENZSQ
      COMMON/BORN/CCF0,CCF0T,CCF0B,CCF0C,CF0G
      COMMON/BORN00/CF0G00
C--MANDELSTAM-VARIABLES:
      SS=AMH**2*(1.D0-X)*(1.D0-Y)
      TT=AMH**2*(1.D0-X)*Y
      UU=AMH**2*X
      CSGQ=DCMPLX(0.D0,0.D0)
      CSGQP=DCMPLX(0.D0,0.D0)
C--FORM FACTOR
      IF(FACT.NE.0.D0)THEN
        CSGQ=CSGQ+FACT*CIGQ(SS,TT,UU,AMQ)
        CSGQP=CSGQP+(FACT-1)*CIGQ(SS,TT,UU,AMQ)
        IF(FACST1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          CSGQ=CSGQ+FACST1*CLGQ(SS,TT,UU,AMST1)*CSBORN(AMH,AMST1)
     .             +FACST2*CLGQ(SS,TT,UU,AMST2)*CSBORN(AMH,AMST2)
         ELSE
          CSGQ=CSGQ+FACST1*CIGQS(SS,TT,UU,AMST1)
     .             +FACST2*CIGQS(SS,TT,UU,AMST2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACB.NE.0.D0)THEN
        CSGQ=CSGQ+FACB*CIGQ(SS,TT,UU,AMB)
        CSGQP=CSGQP+(FACB-1)*CIGQ(SS,TT,UU,AMB)
        IF(FACSB1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          CSGQ=CSGQ+FACSB1*CLGQ(SS,TT,UU,AMSB1)*CSBORN(AMH,AMSB1)
     .             +FACSB2*CLGQ(SS,TT,UU,AMSB2)*CSBORN(AMH,AMSB2)
         ELSE
          CSGQ=CSGQ+FACSB1*CIGQS(SS,TT,UU,AMSB1)
     .             +FACSB2*CIGQS(SS,TT,UU,AMSB2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACC.NE.0.D0)THEN
        CSGQ=CSGQ+FACC*CIGQ(SS,TT,UU,AMC)
        CSGQP=CSGQP+(FACC-1)*CIGQ(SS,TT,UU,AMC)
      ENDIF
      IF(FACCG.NE.0.D0)THEN
        CSGQ=CSGQ-4*CF0G00/3
        CSGQP=CSGQP-4*CF0G00/3
      ENDIF
      IF(FACCTG.NE.0.D0)THEN
        CSGQ=CSGQ+FACCTG*CIGQT(SS,TT,UU,AMQ)
        CSGQP=CSGQP+FACCTG*CIGQT(SS,TT,UU,AMQ)
      ENDIF
C--CURRENT FACTOR
      HGQ=3.D0*(UU**2+SS**2)/(8.D0*TT)
C--LOWEST ORDER FORM FACTOR
      IF(IHIGGS.NE.1)THEN
C--SCALAR HIGGS
        CF0=FACT*CFBORN(AMH,AMQ)+FACB*CFBORN(AMH,AMB)
     .     +FACC*CFBORN(AMH,AMC)+CF0G
     .     +FACST1*CSBORN(AMH,AMST1)+FACST2*CSBORN(AMH,AMST2)
     .     +FACSB1*CSBORN(AMH,AMSB1)+FACSB2*CSBORN(AMH,AMSB2)
        CF0P=(FACT-1)*CFBORN(AMH,AMQ)+(FACB-1)*CFBORN(AMH,AMB)
     .      +(FACC-1)*CFBORN(AMH,AMC)+CF0G
      ELSE
C--PSEUDOSCALAR HIGGS
        CSGQP=0
        CF0P=0
        CF0=FACT*CFBORNA(AMH,AMQ)+FACB*CFBORNA(AMH,AMB)
     .     +FACC*CFBORNA(AMH,AMC)
      ENDIF
      EIGQ=CDABS(CSGQ)**2*HGQ/CDABS(CF0)**2
     .     -2.D0/3.D0*(UU**2+SS**2)/TT
      IF(ISILH.NE.0) EIQQ=(CDABS(CSGQ)**2-CDABS(CSGQP)**2)*HGQ
     .                   / (CDABS(CF0)**2-CDABS(CF0P)**2)
     .                   -2.D0/3.D0*(UU**2+SS**2)/TT
      EIGQ=3.D0/4.D0/AMH**2*EIGQ*(1.D0-X)
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION EIGG(X,Y)
C--H --> GGG: MASS-DEPENDENT INTEGRAND
C--DIFFRENCE TO HEAVY TOP LIMIT
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/HIGGS/IHIGGS
      COMMON/MASS/AMH,AMQ,S
      COMMON/SUSY/AMB,AMC,FACT,FACB,FACC,FACCG,FACCTG,ISILH
      COMMON/SUSY0/AMST1,AMST2,AMSB1,AMSB2,FACST1,FACST2,FACSB1,FACSB2
      COMMON/GRENZSQ/IGRENZSQ
      COMMON/BORN/CCF0,CCF0T,CCF0B,CCF0C,CF0G
      COMMON/BORN00/CF0G00
C--MANDELSTAM-VARIABLES:
      SS=AMH**2*(1.D0-X)*(1.D0-Y)
      TT=AMH**2*(1.D0-X)*Y
      UU=AMH**2*X
      C5=DCMPLX(0.D0,0.D0)
      C2=DCMPLX(0.D0,0.D0)
      C3=DCMPLX(0.D0,0.D0)
      C4=DCMPLX(0.D0,0.D0)
      C5P=DCMPLX(0.D0,0.D0)
      C2P=DCMPLX(0.D0,0.D0)
      C3P=DCMPLX(0.D0,0.D0)
      C4P=DCMPLX(0.D0,0.D0)
C--FORM FACTOR
      IF(FACT.NE.0.D0)THEN
        C5=C5+FACT*CFGG5(SS,TT,UU,AMQ)
        C2=C2+FACT*CFGG2(SS,TT,UU,AMQ)
        C3=C3+FACT*CFGG2(TT,SS,UU,AMQ)
        C4=C4+FACT*CFGG2(UU,TT,SS,AMQ)
        C5P=C5P+(FACT-1)*CFGG5(SS,TT,UU,AMQ)
        C2P=C2P+(FACT-1)*CFGG2(SS,TT,UU,AMQ)
        C3P=C3P+(FACT-1)*CFGG2(TT,SS,UU,AMQ)
        C4P=C4P+(FACT-1)*CFGG2(UU,TT,SS,AMQ)
        IF(FACST1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          C5=C5+FACST1*CFGL5(SS,TT,UU,AMST1)*CSBORN(AMH,AMST1)
     .         +FACST2*CFGL5(SS,TT,UU,AMST2)*CSBORN(AMH,AMST2)
          C2=C2+FACST1*CFGL2(SS,TT,UU,AMST1)*CSBORN(AMH,AMST1)
     .         +FACST2*CFGL2(SS,TT,UU,AMST2)*CSBORN(AMH,AMST2)
          C3=C3+FACST1*CFGL2(TT,SS,UU,AMST1)*CSBORN(AMH,AMST1)
     .         +FACST2*CFGL2(TT,SS,UU,AMST2)*CSBORN(AMH,AMST2)
          C4=C4+FACST1*CFGL2(UU,TT,SS,AMST1)*CSBORN(AMH,AMST1)
     .         +FACST2*CFGL2(UU,TT,SS,AMST2)*CSBORN(AMH,AMST2)
         ELSE
          C5=C5+FACST1*CFGG5S(SS,TT,UU,AMST1)
     .         +FACST2*CFGG5S(SS,TT,UU,AMST2)
          C2=C2+FACST1*CFGG2S(SS,TT,UU,AMST1)
     .         +FACST2*CFGG2S(SS,TT,UU,AMST2)
          C3=C3+FACST1*CFGG2S(TT,SS,UU,AMST1)
     .         +FACST2*CFGG2S(TT,SS,UU,AMST2)
          C4=C4+FACST1*CFGG2S(UU,TT,SS,AMST1)
     .         +FACST2*CFGG2S(UU,TT,SS,AMST2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACB.NE.0.D0)THEN
        C5=C5+FACB*CFGG5(SS,TT,UU,AMB)
        C2=C2+FACB*CFGG2(SS,TT,UU,AMB)
        C3=C3+FACB*CFGG2(TT,SS,UU,AMB)
        C4=C4+FACB*CFGG2(UU,TT,SS,AMB)
        C5P=C5P+(FACB-1)*CFGG5(SS,TT,UU,AMB)
        C2P=C2P+(FACB-1)*CFGG2(SS,TT,UU,AMB)
        C3P=C3P+(FACB-1)*CFGG2(TT,SS,UU,AMB)
        C4P=C4P+(FACB-1)*CFGG2(UU,TT,SS,AMB)
        IF(FACSB1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          C5=C5+FACSB1*CFGL5(SS,TT,UU,AMSB1)*CSBORN(AMH,AMSB1)
     .         +FACSB2*CFGL5(SS,TT,UU,AMSB2)*CSBORN(AMH,AMSB2)
          C2=C2+FACSB1*CFGL2(SS,TT,UU,AMSB1)*CSBORN(AMH,AMSB1)
     .         +FACSB2*CFGL2(SS,TT,UU,AMSB2)*CSBORN(AMH,AMSB2)
          C3=C3+FACSB1*CFGL2(TT,SS,UU,AMSB1)*CSBORN(AMH,AMSB1)
     .         +FACSB2*CFGL2(TT,SS,UU,AMSB2)*CSBORN(AMH,AMSB2)
          C4=C4+FACSB1*CFGL2(UU,TT,SS,AMSB1)*CSBORN(AMH,AMSB1)
     .         +FACSB2*CFGL2(UU,TT,SS,AMSB2)*CSBORN(AMH,AMSB2)
         ELSE
          C5=C5+FACSB1*CFGG5S(SS,TT,UU,AMSB1)
     .         +FACSB2*CFGG5S(SS,TT,UU,AMSB2)
          C2=C2+FACSB1*CFGG2S(SS,TT,UU,AMSB1)
     .         +FACSB2*CFGG2S(SS,TT,UU,AMSB2)
          C3=C3+FACSB1*CFGG2S(TT,SS,UU,AMSB1)
     .         +FACSB2*CFGG2S(TT,SS,UU,AMSB2)
          C4=C4+FACSB1*CFGG2S(UU,TT,SS,AMSB1)
     .         +FACSB2*CFGG2S(UU,TT,SS,AMSB2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACC.NE.0.D0)THEN
        C5=C5+FACC*CFGG5(SS,TT,UU,AMC)
        C2=C2+FACC*CFGG2(SS,TT,UU,AMC)
        C3=C3+FACC*CFGG2(TT,SS,UU,AMC)
        C4=C4+FACC*CFGG2(UU,TT,SS,AMC)
        C5P=C5P+(FACC-1)*CFGG5(SS,TT,UU,AMC)
        C2P=C2P+(FACC-1)*CFGG2(SS,TT,UU,AMC)
        C3P=C3P+(FACC-1)*CFGG2(TT,SS,UU,AMC)
        C4P=C4P+(FACC-1)*CFGG2(UU,TT,SS,AMC)
      ENDIF
      IF(FACCG.NE.0.D0)THEN
        C5=C5+CF0G00*(SS+TT+UU)**2
        C2=C2+CF0G00*SS**2
        C3=C3+CF0G00*TT**2
        C4=C4+CF0G00*UU**2
        C5P=C5P+CF0G00*(SS+TT+UU)**2
        C2P=C2P+CF0G00*SS**2
        C3P=C3P+CF0G00*TT**2
        C4P=C4P+CF0G00*UU**2
      ENDIF
      IF(FACCTG.NE.0.D0)THEN
        C5=C5+FACCTG*CFGG5T(SS,TT,UU,AMQ)
        C2=C2+FACCTG*CFGG2T(SS,TT,UU,AMQ)
        C3=C3+FACCTG*CFGG2T(TT,SS,UU,AMQ)
        C4=C4+FACCTG*CFGG2T(UU,TT,SS,AMQ)
        C5P=C5P+FACCTG*CFGG5T(SS,TT,UU,AMQ)
        C2P=C2P+FACCTG*CFGG2T(SS,TT,UU,AMQ)
        C3P=C3P+FACCTG*CFGG2T(TT,SS,UU,AMQ)
        C4P=C4P+FACCTG*CFGG2T(UU,TT,SS,AMQ)
      ENDIF
      FGG=CDABS(C5)**2+CDABS(C2)**2+CDABS(C3)**2+CDABS(C4)**2
      FGGP=CDABS(C5P)**2+CDABS(C2P)**2+CDABS(C3P)**2+CDABS(C4P)**2
      IF(ISILH.NE.0)FGG = FGG - FGGP
C--LOWEST ORDER FORM FACTOR
      IF(IHIGGS.NE.1)THEN
C--SCALAR HIGGS
        CF0=FACT*CFBORN(AMH,AMQ)+FACB*CFBORN(AMH,AMB)
     .     +FACC*CFBORN(AMH,AMC)+CF0G
     .     +FACST1*CSBORN(AMH,AMST1)+FACST2*CSBORN(AMH,AMST2)
     .     +FACSB1*CSBORN(AMH,AMSB1)+FACSB2*CSBORN(AMH,AMSB2)
        CF0P=(FACT-1)*CFBORN(AMH,AMQ)+(FACB-1)*CFBORN(AMH,AMB)
     .      +(FACC-1)*CFBORN(AMH,AMC)+CF0G
      ELSE
C--PSEUDOSCALAR HIGGS
        FGGP = 0
        CF0P = 0
        CF0=FACT*CFBORNA(AMH,AMQ)+FACB*CFBORNA(AMH,AMB)
     .     +FACC*CFBORNA(AMH,AMC)
      ENDIF
      EIGG=FGG/CDABS(CF0)**2-AMH**8-SS**4-TT**4-UU**4
      IF(ISILH.NE.0) EIGG=(FGG-FGGP)/(CDABS(CF0)**2-CDABS(CF0P)**2)
     .                   -AMH**8-SS**4-TT**4-UU**4
      EIGG=EIGG/AMH**4/TT/UU
      RETURN
      END
 
      COMPLEX*16 FUNCTION CIGQ(SS,TT,UU,AMQ)
C--GQ --> HG: FORM FACTOR
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/MASS/AMH,AMT,S
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/HIGGS/IHIGGS
      SH=SS/AMQ**2
      TH=TT/AMQ**2
      UH=UU/AMQ**2
      RH=SH+TH+UH
 
      IF(IHIGGS.NE.1)THEN
        CIGQ=-(4.D0*CK0(RH)*(UH+SH-4.D0)+4.D0*CK0(TH)*(-UH-SH+4.D0)+
     . 8.D0*CK1(RH)*TH-8.D0*CK1(TH)*TH+8.D0*(UH+SH))/(UH**2+2.D0*UH*
     . SH+SH**2)
      ELSE
        CIGQ=-(8.*CK0(RH)-8.*CK0(TH))/(3.*(UH+SH))
      ENDIF
C     AMH=DSQRT(RH*AMQ**2)
C     CIGQ = CLGQ(SS,TT,UU,AMQ)*CFBORN(AMH,AMQ)
      RETURN
      END
 
      COMPLEX*16 FUNCTION CIGQT(SS,TT,UU,AMQ)
C--GQ --> HG: FORM FACTOR CHROMOMAGNETIC OPERATOR
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION LLR
      COMMON/MASS/AMH,AMT,S
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/HIGGS/IHIGGS
      SH=SS/AMQ**2
      TH=TT/AMQ**2
      UH=UU/AMQ**2
      RH=SH+TH+UH
      SCALE = AMQ
      SCALE = 2*AMQ
C     SCALE = 1.D50*AMQ
      LLR = DLOG(SCALE**2/AMQ**2)
      cigqt=(8*ck0(rh)*(-sh+th-uh)+8*ck0(th)*(sh-th+uh)+4*ck1(rh)*(sh
     . **2+2*sh*uh-th**2+uh**2)+2*ck1(th)*(sh**2+2*sh*uh+2*th**2+uh**
     . 2)+4*(2*llr*sh**2+4*llr*sh*uh+2*llr*uh**2-sh**2-sh*th-2*sh*uh-
     . th*uh-uh**2))/(sh**2+2*sh*uh+uh**2)
      RETURN
      END
 
      COMPLEX*16 FUNCTION CIGQS(SS,TT,UU,AMQ)
C--GQ --> HG: FORM FACTOR SQUARK LOOPS
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/MASS/AMH,AMT,S
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/HIGGS/IHIGGS
      SH=SS/AMQ**2
      TH=TT/AMQ**2
      UH=UU/AMQ**2
      RH=SH+TH+UH
 
      IF(IHIGGS.NE.1)THEN
        CIGQS=(-32*CK0(RH)+32*CK0(TH)+16*CK1(RH)*TH-16*CK1(TH)*TH+16*(SH
     .  +UH))/(SH**2+2*SH*UH+UH**2)
     .  /4
      ELSE
        CIGQS=0
      ENDIF
c     AMH=DSQRT(RH*AMQ**2)
c     CIGQS0 = CLGQ(SS,TT,UU,AMQ)*CSBORN(AMH,AMQ)
c     write(6,*)CIGQS/CIGQS0
      RETURN
      END
 
      COMPLEX*16 FUNCTION CLGQ(SS,TT,UU,AMQ)
C--GQ --> HG: FORM FACTOR FOR MQ -> INFINITY
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/MASS/AMH,AMT,S
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/HIGGS/IHIGGS
      SH=SS/AMQ**2
      TH=TT/AMQ**2
      UH=UU/AMQ**2
      RH=SH+TH+UH
      CLGQ=-4.D0/3.D0
      RETURN
      END
 
      COMPLEX*16 FUNCTION CFGG5(SS,TT,UU,AMQ)
C--GG --> HG: FORM FACTOR C_5
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/HIGGS/IHIGGS
      SH=SS/AMQ**2
      TH=TT/AMQ**2
      UH=UU/AMQ**2
      RH=SH+TH+UH
      IF(IHIGGS.NE.1)THEN
        CFGG5=-(18.D0*CK0(RH)*(-UH-TH-SH+4.D0)+6.D0*CK0(UH)*(UH+TH+
     . SH-4.D0)+6.D0*CK0(TH)*(UH+TH+SH-4.D0)+6.D0*CK0(SH)*(UH+TH+SH
     . -4.D0)+3.D0*CJ(TH,SH,UH)*UH*TH*(UH+TH+SH-4.D0)+3.D0*CJ(SH,UH
     . ,TH)*TH*SH*(UH+TH+SH-4.D0)+3.D0*CJ(SH,TH,UH)*UH*SH*(UH+
     . TH+SH-4.D0)-24.D0*(UH+TH+SH))/2.D0
      ELSE
        CFGG5=-(CJ(TH,SH,UH)*UH*TH+CJ(SH,UH,TH)*TH*SH+CJ(SH,
     . TH,UH)*UH*SH-6.*CK0(RH)+2.*CK0(UH)+2.*CK0(TH)+2.*CK0
     . (SH))*(UH+TH+SH)
      ENDIF
      CFGG5=CFGG5/2.D0*AMQ**4
C     AMH=DSQRT(RH*AMQ**2)
C     CFGG5=CFGL5(SS,TT,UU,AMQ)*CFBORN(AMH,AMQ)
C     write(6,*)cfgg5,CFGL5(SS,TT,UU,AMQ)*CFBORN(AMH,AMQ),
C    .          cfgg5/CFGL5(SS,TT,UU,AMQ)/CFBORN(AMH,AMQ)
      RETURN
      END
 
      COMPLEX*16 FUNCTION CFGG5T(SS,TT,UU,AMQ)
C--GG --> HG: FORM FACTOR C_5 CHROMOMAGNETIC OPERATOR
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION LLR
      COMMON/HIGGS/IHIGGS
      SH=SS/AMQ**2
      TH=TT/AMQ**2
      UH=UU/AMQ**2
      RH=SH+TH+UH
      SCALE = AMQ
      SCALE = 2*AMQ
C     SCALE = 1.D50*AMQ
      LLR = DLOG(SCALE**2/AMQ**2)
      cfgg5t=3*cj(sh,th,uh)*sh*uh*(-sh-2*th-uh)+3*cj(sh,uh,th)*sh*th*(
     . -sh-th-2*uh)+3*cj(th,sh,uh)*th*uh*(-2*sh-th-uh)+18*ck0(rh)*(sh
     . +th+uh)+3*ck0(sh)*(-2*sh+th*uh-2*th-2*uh)+3*ck0(th)*(sh*uh-2*
     . sh-2*th-2*uh)+3*ck0(uh)*(sh*th-2*sh-2*th-2*uh)+6*ck1(rh)*(-sh
     . **2-2*sh*th-2*sh*uh-th**2-2*th*uh-uh**2)-6*ck1(sh)*th*uh-6*ck1
     . (th)*sh*uh-6*ck1(uh)*sh*th+6*(-2*llr*sh**2-4*llr*sh*th-4*llr*
     . sh*uh-2*llr*th**2-4*llr*th*uh-2*llr*uh**2+sh**2+2*sh*th+2*sh*
     . uh+th**2+2*th*uh+uh**2)
      CFGG5T=CFGG5T/2.D0*AMQ**4
      RETURN
      END
 
      COMPLEX*16 FUNCTION CFGG2(SS,TT,UU,AMQ)
C--GG --> HG: FORM FACTOR C_2
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/HIGGS/IHIGGS
      SH=SS/AMQ**2
      TH=TT/AMQ**2
      UH=UU/AMQ**2
      RH=SH+TH+UH
      IF(IHIGGS.NE.1)THEN
        COF00=-(12.D0*(UH*TH-SH**2)*SH)/((UH+SH)*(TH+SH))
        COF0R=(3.D0*(4.D0*UH**3*TH**3+8.D0*UH**3*TH**2*SH+4.D0*UH**3*
     . TH*SH**2+8.D0*UH**2*TH**3*SH+15.D0*UH**2*TH**2*SH**2+4.D0*
     . UH**2*TH**2*SH+8.D0*UH**2*TH*SH**3+8.D0*UH**2*TH*SH**2+
     . UH**2*SH**4-4.D0*UH**2*SH**3+4.D0*UH*TH**3*SH**2+8.D0*UH*
     . TH**2*SH**3+8.D0*UH*TH**2*SH**2+8.D0*UH*TH*SH**4+16.D0*UH*
     . TH*SH**3+4.D0*UH*SH**5-8.D0*UH*SH**4+TH**2*SH**4-4.D0*TH**
     . 2*SH**3+4.D0*TH*SH**5-8.D0*TH*SH**4+3.D0*SH**6-12.D0*SH**5))
     . /((UH+SH)**2*(TH+SH)**2*SH)
        COF0S=-3.D0*(SH-4.D0)
        COF0T=-(3.D0*(4.D0*UH**3*TH+8.D0*UH**2*TH*SH-UH**2*SH**2+4.D0
     . *UH**2*SH+4.D0*UH*TH*SH**2+8.D0*UH*SH**2+SH**4-4.D0*SH**3)
     . )/((UH+SH)**2*SH)
        COF0U=-(3.D0*(4.D0*UH*TH**3+8.D0*UH*TH**2*SH+4.D0*UH*TH*SH**2
     . -TH**2*SH**2+4.D0*TH**2*SH+8.D0*TH*SH**2+SH**4-4.D0*SH**3)
     . )/((TH+SH)**2*SH)
        COF1R=-(12.D0*(UH**2*TH+2.D0*UH**2*SH+UH*TH**2+4.D0*UH*TH*
     . SH+5.D0*UH*SH**2+2.D0*TH**2*SH+5.D0*TH*SH**2+4.D0*SH**3)*UH*
     . TH)/((UH+SH)**2*(TH+SH)**2)
        COF1S=0.D0
        COF1T=(12.D0*(UH+2.D0*SH)*UH*TH)/(UH+SH)**2
        COF1U=(12.D0*(TH+2.D0*SH)*UH*TH)/(TH+SH)**2
        COFJ1=-(3.D0*(SH-4.D0)*UH*SH)/2.D0
        COFJ2=-(3.D0*(SH-4.D0)*TH*SH)/2.D0
        COFJ3=-(3.D0*(4.D0*UH*TH-SH**2+12.D0*SH)*UH*TH)/(2.D0*SH)
      ELSE
        COF00=0.D0
        COF0R=-(2.D0*(UH*TH-UH*SH-TH*SH-3.D0*SH**2)*SH)/((UH+SH)*
     . (TH+SH))
        COF0S=-2.D0*SH
        COF0T=(2.D0*(UH-SH)*SH)/(UH+SH)
        COF0U=(2.D0*(TH-SH)*SH)/(TH+SH)
        COF1R=0.D0
        COF1S=0.D0
        COF1T=0.D0
        COF1U=0.D0
        COFJ1=-UH*SH**2
        COFJ2=-TH*SH**2
        COFJ3=UH*TH*SH
      ENDIF
      CFGG2=CK0(RH)*COF0R+CK0(UH)*COF0U+CK0(TH)*COF0T+CK0(
     . SH)*COF0S+CK1(RH)*COF1R+CK1(UH)*COF1U+CK1(TH)*COF1T+
     . CK1(SH)*COF1S+CJ(TH,SH,UH)*COFJ3+CJ(SH,UH,TH)*COFJ2+
     . CJ(SH,TH,UH)*COFJ1+COF00
      CFGG2=CFGG2/2.D0*AMQ**4
C     AMH=DSQRT(RH*AMQ**2)
C     CFGG2=CFGL2(SS,TT,UU,AMQ)*CFBORN(AMH,AMQ)
      RETURN
      END
 
      COMPLEX*16 FUNCTION CFGG2T(SS,TT,UU,AMQ)
C--GG --> HG: FORM FACTOR C_2 CHROMOMAGNETIC OPERATOR
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION LLR
      COMMON/HIGGS/IHIGGS
      SH=SS/AMQ**2
      TH=TT/AMQ**2
      UH=UU/AMQ**2
      RH=SH+TH+UH
      SCALE = AMQ
      SCALE = 2*AMQ
C     SCALE = 1.D50*AMQ
      LLR = DLOG(SCALE**2/AMQ**2)
      cof00=(6*((sh**2-th*uh)*(sh+th+uh)-2*(sh+th)*(sh+uh)*llr*sh)*sh
     . )/((sh+th)*(sh+uh))
      cof0r=(3*((6*sh**5+sh**4*th*uh+6*sh**4*th+6*sh**4*uh+4*sh**3*th
     . **2*uh-2*sh**3*th**2+4*sh**3*th*uh**2+16*sh**3*th*uh-2*sh**3*
     . uh**2+5*sh**2*th**3*uh-2*sh**2*th**3+12*sh**2*th**2*uh**2+10*
     . sh**2*th**2*uh+5*sh**2*th*uh**3+10*sh**2*th*uh**2-2*sh**2*uh**
     . 3+2*sh*th**4*uh+12*sh*th**3*uh**2+4*sh*th**3*uh+12*sh*th**2*uh
     . **3+6*sh*th**2*uh**2+2*sh*th*uh**4+4*sh*th*uh**3+4*th**4*uh**2
     . +9*th**3*uh**3+2*th**3*uh**2+4*th**2*uh**4+2*th**2*uh**3)*sh+2
     . *(th+uh)*th**3*uh**3))/((sh+th)**2*(sh+uh)**2*sh)
      cof0s=6*(th+uh-sh)
      cof0t=(-3*((2*sh**3+sh**2*th*uh-2*sh**2*th+2*sh**2*uh+2*sh*th**
     . 2*uh+4*sh*th*uh**2+4*sh*th*uh+2*sh*uh**2+4*th**2*uh**2+5*th*uh
     . **3+2*th*uh**2+2*uh**3)*sh+2*(th+uh)*th*uh**3))/((sh+uh)**2*sh
     . )
      cof0u=(-3*((2*sh**3+sh**2*th*uh+2*sh**2*th-2*sh**2*uh+4*sh*th**
     . 2*uh+2*sh*th**2+2*sh*th*uh**2+4*sh*th*uh+5*th**3*uh+2*th**3+4*
     . th**2*uh**2+2*th**2*uh)*sh+2*(th+uh)*th**3*uh))/((sh+th)**2*sh
     . )
      cof1r=(-6*((sh**4+sh**3*th+sh**3*uh+4*sh**2*th*uh+4*sh*th**2*uh
     . +4*sh*th*uh**2+2*th**3*uh+3*th**2*uh**2+2*th*uh**3)*sh+(th+uh)
     . *th**2*uh**2)*(th+uh+sh))/((sh+th)**2*(sh+uh)**2)
      cof1s=0
      cof1t=(6*((sh+2*th+2*uh)*sh+(th+uh)*uh)*th*uh)/(sh+uh)**2
      cof1u=(6*((sh+2*th+2*uh)*sh+(th+uh)*th)*th*uh)/(sh+th)**2
      cofj1=-3*(sh-uh)*sh*uh
      cofj2=-3*(sh-th)*sh*th
      cofj3=(-3*((th*uh+6*th+6*uh)*sh+2*(th+uh)*th*uh)*th*uh)/(2*sh)
      CFGG2T=CK0(RH)*COF0R+CK0(UH)*COF0U+CK0(TH)*COF0T+CK0(
     . SH)*COF0S+CK1(RH)*COF1R+CK1(UH)*COF1U+CK1(TH)*COF1T+
     . CK1(SH)*COF1S+CJ(TH,SH,UH)*COFJ3+CJ(SH,UH,TH)*COFJ2+
     . CJ(SH,TH,UH)*COFJ1+COF00
      CFGG2T=CFGG2T/2.D0*AMQ**4
      RETURN
      END
 
      COMPLEX*16 FUNCTION CFGG5S(SS,TT,UU,AMQ)
C--GG --> HG: FORM FACTOR C_5 SQUARK LOOPS
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/HIGGS/IHIGGS
      SH=SS/AMQ**2
      TH=TT/AMQ**2
      UH=UU/AMQ**2
      RH=SH+TH+UH
      IF(IHIGGS.NE.1)THEN
        CFGG5S=72*CK0(RH)-24*CK0(SH)-24*CK0(TH)-24*CK0(UH)-12*CJ(SH,TH,
     . UH)*SH*UH-12*CJ(SH,UH,TH)*SH*TH-12*CJ(TH,SH,UH)*TH*UH-24*(SH+
     . TH+UH)
      ELSE
        CFGG5S=0
      ENDIF
      CFGG5S=CFGG5S/2.D0*AMQ**4
     .      /4
      RETURN
      END
 
      COMPLEX*16 FUNCTION CFGG2S(SS,TT,UU,AMQ)
C--GG --> HG: FORM FACTOR C_2 SQUARK LOOPS
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/HIGGS/IHIGGS
      SH=SS/AMQ**2
      TH=TT/AMQ**2
      UH=UU/AMQ**2
      RH=SH+TH+UH
      IF(IHIGGS.NE.1)THEN
       COF00=(-24*(SH**2-TH*UH)*SH)/((SH+TH)*(SH+UH))
       COF0R=(24*(3*SH**5-SH**4*TH*UH+2*SH**4*TH+2*SH**4*UH-2*SH**3*TH
     . **2*UH+SH**3*TH**2-2*SH**3*TH*UH**2-4*SH**3*TH*UH+SH**3*UH**2-
     . SH**2*TH**3*UH-4*SH**2*TH**2*UH**2-2*SH**2*TH**2*UH-SH**2*TH*
     . UH**3-2*SH**2*TH*UH**2-2*SH*TH**3*UH**2-2*SH*TH**2*UH**3-SH*TH
     . **2*UH**2-TH**3*UH**3))/((SH+TH)**2*(SH+UH)**2*SH)
       COF0S=(-24)
       COF0T=(-24*(SH**3-SH**2*TH*UH-2*SH**2*UH-2*SH*TH*UH**2-SH*UH**2
     . -TH*UH**3))/((SH+UH)**2*SH)
       COF0U=(-24*(SH**3-SH**2*TH*UH-2*SH**2*TH-2*SH*TH**2*UH-SH*TH**2
     . -TH**3*UH))/((SH+TH)**2*SH)
       COF1R=(24*((4*SH**2+5*SH*TH+5*SH*UH+2*TH**2+4*TH*UH+2*UH**2)*SH
     . +(TH+UH)*TH*UH)*TH*UH)/((SH+TH)**2*(SH+UH)**2)
       COF1S=0
       COF1T=(-24*(2*SH+UH)*TH*UH)/(SH+UH)**2
       COF1U=(-24*(2*SH+TH)*TH*UH)/(SH+TH)**2
       COFJ1=-12*SH*UH
       COFJ2=-12*SH*TH
       COFJ3=(12*(3*SH+TH*UH)*TH*UH)/SH
      ELSE
       COF00=0
       COF0R=0
       COF0S=0
       COF0T=0
       COF0U=0
       COF1R=0
       COF1S=0
       COF1T=0
       COF1U=0
       COFJ1=0
       COFJ2=0
       COFJ3=0
      ENDIF
      CFGG2S=CK0(RH)*COF0R+CK0(UH)*COF0U+CK0(TH)*COF0T+CK0(
     . SH)*COF0S+CK1(RH)*COF1R+CK1(UH)*COF1U+CK1(TH)*COF1T+
     . CK1(SH)*COF1S+CJ(TH,SH,UH)*COFJ3+CJ(SH,UH,TH)*COFJ2+
     . CJ(SH,TH,UH)*COFJ1+COF00
      CFGG2S=CFGG2S/2.D0*AMQ**4
     .      /4
      RETURN
      END
 
      COMPLEX*16 FUNCTION CFGL5(SS,TT,UU,AMQ)
C--GG --> HG: FORM FACTOR C_5 FOR MQ --> INFINITY
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/HIGGS/IHIGGS
      SH=SS/AMQ**2
      TH=TT/AMQ**2
      UH=UU/AMQ**2
      RH=SH+TH+UH
      CFGL5=2.D0*RH**2
      CFGL5=CFGL5/2.D0*AMQ**4
      RETURN
      END
 
      COMPLEX*16 FUNCTION CFGL2(SS,TT,UU,AMQ)
C--GG --> HG: FORM FACTOR C_2 FOR MQ -> INFINITY
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/HIGGS/IHIGGS
      SH=SS/AMQ**2
      TH=TT/AMQ**2
      UH=UU/AMQ**2
      RH=SH+TH+UH
      CFGL2=2*SH**2
      CFGL2=CFGL2/2.D0*AMQ**4
      RETURN
      END
 
      COMPLEX*16 FUNCTION CK0(X)
C--SCALAR INTEGRAL
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/CUT/EPST,EPSV,REPS
      CX=X/DCMPLX(1.D0,-REPS)
      CTAU=4.D0/CX
      CK0=2.D0*CF(CTAU)
      RETURN
      END
 
      COMPLEX*16 FUNCTION CK1(X)
C--SCALAR INTEGRAL
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/CUT/EPST,EPSV,REPS
      CX=X/DCMPLX(1.D0,-REPS)
      CTAU=4.D0/CX
      CK1=2.D0*(1.D0-CG(CTAU))
      RETURN
      END
 
      COMPLEX*16 FUNCTION CSJ(X,Y,Z,V)
C--SCALAR INTEGRAL
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMPLEX*16 LI2
      COMMON/CUT/EPST,EPSV,REPS
      CV=V/DCMPLX(1.D0,-REPS)
      CAP=(1.D0+CDSQRT(1.D0-4.D0/CV))/2.D0
      CAM=1.D0/CV/CAP
      CBP=DCMPLX(1.D0+DSQRT(1.D0+4.D0*Y/X/Z),0.D0)/2.D0
      CBM=-Y/X/Z/CBP
      CSJ=2.D0/X/Z/(CBP-CBM)
     .    *(LI2(CBM/(CBM-CAM))-LI2(CBP/(CBP-CAP))
     .     +LI2(CBM/(CBM-CAP))-LI2(CBP/(CBP-CAM))
     .     +CDLOG(-CBP/CBM)*CDLOG(1.D0+CV*Y/X/Z))
      RETURN
      END
 
      COMPLEX*16 FUNCTION CJ(X,Y,Z)
C--SCALAR BOX-INTEGRAL
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      RH=X+Y+Z
      CJ=CSJ(X,Y,Z,X)+CSJ(X,Y,Z,Z)-CSJ(X,Y,Z,RH)
      RETURN
      END
 
      COMPLEX*16 FUNCTION CFBORN(AMH,AMQ)
C--LOWEST ORDER FORM FACTOR FOR SCALAR HIGGS
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/CUT/EPST,EPSV,REPS
      CTAU=4.D0*AMQ**2/AMH**2*DCMPLX(1.D0,-REPS)
      CFBORN=1.5D0*CTAU*(1.D0+(1.D0-CTAU)*CF(CTAU))
      RETURN
      END
 
      COMPLEX*16 FUNCTION CFBORNT(AMH,AMQ)
C--LOWEST ORDER FORM FACTOR FOR SCALAR HIGGS: CHROMOMAGNETIC OPERATOR
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION LLR
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/SUSYT/SCALCG
      CALL SCALES(QM,QQ)
      SCALE = SCALCG
      SCALE = QM
      SCALE = AMQ
      SCALE = 2*AMQ
      LLR = DLOG(SCALE**2/AMQ**2)
      CTAU=4.D0*AMQ**2/AMH**2*DCMPLX(1.D0,-REPS)
      CFBORNT=3*(CTAU*CF(CTAU)+2*CG(CTAU)-1-2*LLR)
      RETURN
      END
 
      COMPLEX*16 FUNCTION CSBORN(AMH,AMQ)
C--LOWEST ORDER FORM FACTOR FOR SCALAR HIGGS
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/CUT/EPST,EPSV,REPS
      CTAU=4.D0*AMQ**2/AMH**2*DCMPLX(1.D0,-REPS)
      CSBORN=-0.75D0*CTAU*(1.D0-CTAU*CF(CTAU))
      RETURN
      END
 
      COMPLEX*16 FUNCTION CFBORNA(AMH,AMQ)
C--LOWEST ORDER FORM FACTOR FOR PSEUDOSCALAR HIGGS
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/CUT/EPST,EPSV,REPS
      CTAU=4.D0*AMQ**2/AMH**2*DCMPLX(1.D0,-REPS)
      CFBORNA=CTAU*CF(CTAU)
      RETURN
      END

      COMPLEX*16 FUNCTION CDCT(AMH,AMQ)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMPLEX*16 LI2,LI3,S12
      COMMON/CUT/EPST,EPSV,REPS
      CTAU=4.D0*AMQ**2/AMH**2*DCMPLX(1.D0,-REPS)
      CDCT=6.D0*CTAU*((2.D0*CTAU-1.D0)*CF(CTAU)-1.D0-CG(CTAU))
      CDCT=CDCT/CFBORN(AMH,AMQ)
      RETURN
      END

      COMPLEX*16 FUNCTION CDCTA(AMH,AMQ)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMPLEX*16 LI2,LI3,S12
      COMMON/CUT/EPST,EPSV,REPS
      CTAU=4.D0*AMQ**2/AMH**2*DCMPLX(1.D0,-REPS)
      CDCTA=-4.D0*CTAU*(CF(CTAU)-CG(CTAU)/(CTAU-1.D0))
      CDCTA=CDCTA/CFBORNA(AMH,AMQ)
      RETURN
      END
 
      COMPLEX*16 FUNCTION CG(CTAU)
C--SCALAR INTEGRAL
      IMPLICIT COMPLEX*16 (C)
      CETA=CDSQRT(1.D0-CTAU)
      CALP=1.D0+CETA
      CALM=CTAU/CALP
      CG=CDLOG(-CALP/CALM)*CETA/2.D0
      RETURN
      END
 
      COMPLEX*16 FUNCTION CF(CTAU)
C--SCALAR INTEGRAL
      IMPLICIT COMPLEX*16 (C)
      CETA=CDSQRT(1.D0-CTAU)
      CALP=1.D0+CETA
      CALM=CTAU/CALP
      CF=-CDLOG(-CALP/CALM)**2/4.D0
      RETURN
      END
 
      SUBROUTINE SCALES(QM,QQ)
C--CALCULATION OF RENORMALIZATION-SCALE QM AND
C--FACTORIZATION-SCALE QQ
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/MASS/AMH,AMQ,S
      COMMON/PARSC/XKAPM,XKAPQ
      COMMON/QLIM/QMAX,QMIN
      QM=XKAPM*AMH
      QQ=XKAPQ*AMH
C--IF FACTORIZATION-SCALE OUTSIDE THE VALIDITY-RANGE OF THE
C--STRUCTURE-FUNCTIONS: FREEZE THEM
      IF(QQ.GT.QMAX)QQ=QMAX
      IF(QQ.LT.QMIN)QQ=QMIN
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION ALPHAS(Q,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XLB(6)
      COMMON/ALSLAM/XLB1(6),XLB2(6),XLB3(6),XLB4(6)
      COMMON/ALS/XLAMBDA,AMC,AMB,AMT,N0
      COMMON/ALSN3LO/N3LO
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      B2(NF)=27/2.D0*(2857-5033/9.D0*NF+325/27.D0*NF**2)/B0(NF)**3
      B3(NF)= 81*(149753/6.d0+3564*zeta3-(1078361/162.d0+6508*zeta3/27)
     .        *nf+(50065/162.d0+6472*zeta3/81)*nf**2+1093/729.d0*nf**3)
     .      / B0(NF)**4
      ALS1(NF,X)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
      ALS2(NF,X)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
     .          *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB(NF)**2))
     .           /DLOG(X**2/XLB(NF)**2))
      ALS3(NF,X)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
     .          *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB(NF)**2))
     .           /DLOG(X**2/XLB(NF)**2)
     .           +(B1(NF)**2*(DLOG(DLOG(X**2/XLB(NF)**2))**2
     .                      -DLOG(DLOG(X**2/XLB(NF)**2))-1)+B2(NF))
     .           /DLOG(X**2/XLB(NF)**2)**2)
      ALS4(NF,X)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
     .          *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB(NF)**2))
     .           /DLOG(X**2/XLB(NF)**2)
     .           +(B1(NF)**2*(DLOG(DLOG(X**2/XLB(NF)**2))**2
     .                      -DLOG(DLOG(X**2/XLB(NF)**2))-1)+B2(NF))
     .           /DLOG(X**2/XLB(NF)**2)**2
     .           -(B1(NF)**3*(DLOG(DLOG(X**2/XLB(NF)**2))**3
     .                      -5*DLOG(DLOG(X**2/XLB(NF)**2))**2/2
     .                      -2*DLOG(DLOG(X**2/XLB(NF)**2))+1/2.d0)
     .            +3*B1(NF)*B2(NF)*DLOG(DLOG(X**2/XLB(NF)**2))
     .            -B3(NF)/2)/DLOG(X**2/XLB(NF)**2)**3)
      PI=4.D0*DATAN(1.D0)
      ZETA2 = PI**2/6
      ZETA3 = 1.2020569031595942853997381D0
c     write(6,*)'ALS param: ',XLAMBDA,AMC,AMB,AMT,N0
      IF(N.EQ.1)THEN
       DO 1 I=1,6
        XLB(I)=XLB1(I)
1      CONTINUE
      ELSEIF(N.EQ.2)THEN
       DO 2 I=1,6
        XLB(I)=XLB2(I)
2      CONTINUE
      ELSEIF(N.EQ.3)THEN
       DO 3 I=1,6
        XLB(I)=XLB3(I)
3      CONTINUE
      ELSE
       DO 4 I=1,6
        XLB(I)=XLB4(I)
4      CONTINUE
      ENDIF
      IF(Q.LT.AMC)THEN
       NF=3
      ELSEIF(Q.LE.AMB)THEN
       NF=4
      ELSEIF(Q.LE.AMT)THEN
       NF=5
      ELSE
       NF=6
      ENDIF
      IF(N.EQ.1)THEN
        ALPHAS=ALS1(NF,Q)
      ELSEIF(N.EQ.2)THEN
        ALPHAS=ALS2(NF,Q)
      ELSEIF(N.EQ.3)THEN
        ALPHAS=ALS3(NF,Q)
      ELSE
        ALPHAS=ALS4(NF,Q)
      ENDIF
      IF(N3LO.NE.0) ALPHAS=ALS4(NF,Q)
      RETURN
      END

      SUBROUTINE ALSINI(ACC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XLB(6)
      COMMON/ALSLAM/XLB1(6),XLB2(6),XLB3(6),XLB4(6)
      COMMON/ALS/XLAMBDA,AMC,AMB,AMT,N0
      PI=4.D0*DATAN(1.D0)
      XLB1(1)=0D0
      XLB1(2)=0D0
      XLB2(1)=0D0
      XLB2(2)=0D0
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2.D0/25.D0)
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2.D0/23.D0)
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
      ENDIF
      DO 1 I=3,6
       XLB1(I)=XLB(I)
1     CONTINUE
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2.D0/25.D0)
     .             *(2.D0*DLOG(AMC/XLB(3)))**(-107.D0/1875.D0)
       XLB(4)=XITER(AMC,XLB(3),3,XLB(4),4,ACC)
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
     .             *(2.D0*DLOG(AMB/XLB(4)))**(-963.D0/13225.D0)
       XLB(5)=XITER(AMB,XLB(4),4,XLB(5),5,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
     .             *(2.D0*DLOG(AMB/XLB(4)))**(-963.D0/13225.D0)
       XLB(5)=XITER(AMB,XLB(4),4,XLB(5),5,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
     .             *(2.D0*DLOG(AMB/XLB(5)))**(963.D0/14375.D0)
       XLB(4)=XITER(AMB,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2.D0/23.D0)
     .            *(2.D0*DLOG(AMT/XLB(6)))**(321.D0/3703.D0)
       XLB(5)=XITER(AMT,XLB(6),6,XLB(5),5,ACC)
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
     .             *(2.D0*DLOG(AMB/XLB(5)))**(963.D0/14375.D0)
       XLB(4)=XITER(AMB,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
      ENDIF
      DO 2 I=3,6
       XLB2(I)=XLB(I)
2     CONTINUE
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2.D0/25.D0)
     .             *(2.D0*DLOG(AMC/XLB(3)))**(-107.D0/1875.D0)
       XLB(4)=XITER3(AMC,XLB(3),3,XLB(4),4,ACC)
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
     .             *(2.D0*DLOG(AMB/XLB(4)))**(-963.D0/13225.D0)
       XLB(5)=XITER3(AMB,XLB(4),4,XLB(5),5,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER3(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
     .             *(2.D0*DLOG(AMB/XLB(4)))**(-963.D0/13225.D0)
       XLB(5)=XITER3(AMB,XLB(4),4,XLB(5),5,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER3(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER3(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
     .             *(2.D0*DLOG(AMB/XLB(5)))**(963.D0/14375.D0)
       XLB(4)=XITER3(AMB,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER3(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER3(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2.D0/23.D0)
     .            *(2.D0*DLOG(AMT/XLB(6)))**(321.D0/3703.D0)
       XLB(5)=XITER3(AMT,XLB(6),6,XLB(5),5,ACC)
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
     .             *(2.D0*DLOG(AMB/XLB(5)))**(963.D0/14375.D0)
       XLB(4)=XITER3(AMB,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER3(AMC,XLB(4),4,XLB(3),3,ACC)
      ENDIF
      DO 3 I=3,6
       XLB3(I)=XLB(I)
3     CONTINUE
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2.D0/25.D0)
     .             *(2.D0*DLOG(AMC/XLB(3)))**(-107.D0/1875.D0)
       XLB(4)=XITER4(AMC,XLB(3),3,XLB(4),4,ACC)
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
     .             *(2.D0*DLOG(AMB/XLB(4)))**(-963.D0/13225.D0)
       XLB(5)=XITER4(AMB,XLB(4),4,XLB(5),5,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER4(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
     .             *(2.D0*DLOG(AMB/XLB(4)))**(-963.D0/13225.D0)
       XLB(5)=XITER4(AMB,XLB(4),4,XLB(5),5,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER4(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER4(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
     .             *(2.D0*DLOG(AMB/XLB(5)))**(963.D0/14375.D0)
       XLB(4)=XITER4(AMB,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER4(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER4(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2.D0/23.D0)
     .            *(2.D0*DLOG(AMT/XLB(6)))**(321.D0/3703.D0)
       XLB(5)=XITER4(AMT,XLB(6),6,XLB(5),5,ACC)
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
     .             *(2.D0*DLOG(AMB/XLB(5)))**(963.D0/14375.D0)
       XLB(4)=XITER4(AMB,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER4(AMC,XLB(4),4,XLB(3),3,ACC)
      ENDIF
      DO 4 I=3,6
       XLB4(I)=XLB(I)
4     CONTINUE
c     write(6,*)'Lambda =',XLB(4),XLB(5)
      RETURN
      END

      DOUBLE PRECISION FUNCTION XITER(Q,XLB1,NF1,XLB,NF2,ACC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      ALS2(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .              *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .              /DLOG(X**2/XLB**2))
      AA(NF)=12D0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      XIT(A,B,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)))
      PI=4.D0*DATAN(1.D0)
      XLB2=XLB
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB2**2)
      ALP=ALS2(NF1,Q,XLB1)
      A=AA(NF2)/ALP
      B=BB(NF2)*ALP
      XX=XIT(A,B,X)
      XLB2=Q*DEXP(-XX/2.D0)
      Y1=ALS2(NF1,Q,XLB1)
      Y2=ALS2(NF2,Q,XLB2)
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
       XITER=XLB2
      RETURN
      END

      DOUBLE PRECISION FUNCTION XITER3(Q,XLB1,NF1,XLB,NF2,ACC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      B2(NF)=27/2.D0*(2857-5033/9.D0*NF+325/27.D0*NF**2)/B0(NF)**3
      ALS3(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .          *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .           /DLOG(X**2/XLB**2)
     .           +(B1(NF)**2*(DLOG(DLOG(X**2/XLB**2))**2
     .                      -DLOG(DLOG(X**2/XLB**2))-1)+B2(NF))
     .           /DLOG(X**2/XLB**2)**2)
      AA(NF)=12D0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      CC(NF)=B2(NF)/AA(NF)
      XIT(A,B,C,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)
     .          *(1-(A*B*(DLOG(X)**2-DLOG(X)-1)+C/B)/X/DLOG(X))))
      PI=4.D0*DATAN(1.D0)
      XLB2=XLB
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB2**2)
      IF(NF1.LT.NF2)THEN
       DELTA = 7*ALS3(NF1,Q,XLB1)**2/PI**2/24
       ALP=ALS3(NF1,Q,XLB1)*(1+DELTA)
      ELSE
       DELTA = 7*ALS3(NF1,Q,XLB1)**2/PI**2/24
       ALP=ALS3(NF1,Q,XLB1)/(1+DELTA)
      ENDIF
      A=AA(NF2)/ALP
      B=BB(NF2)*ALP
      C=CC(NF2)*ALP
      XX=XIT(A,B,C,X)
      XLB2=Q*DEXP(-XX/2.D0)
      IF(NF1.LT.NF2)THEN
       DELTA = 7*ALS3(NF1,Q,XLB1)**2/PI**2/24
       Y1=ALS3(NF1,Q,XLB1)*(1+DELTA)
       Y2=ALS3(NF2,Q,XLB2)
      ELSE
       DELTA = 7*ALS3(NF1,Q,XLB1)**2/PI**2/24
       Y1=ALS3(NF1,Q,XLB1)/(1+DELTA)
       Y2=ALS3(NF2,Q,XLB2)
      ENDIF
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
       XITER3=XLB2
      RETURN
      END

      DOUBLE PRECISION FUNCTION XITER4(Q,XLB1,NF1,XLB,NF2,ACC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      B2(NF)=27/2.D0*(2857-5033/9.D0*NF+325/27.D0*NF**2)/B0(NF)**3
      B3(NF)= 81*(149753/6.d0+3564*zeta3-(1078361/162.d0+6508*zeta3/27)
     .        *nf+(50065/162.d0+6472*zeta3/81)*nf**2+1093/729.d0*nf**3)
     .      / B0(NF)**4
      ALS4(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .          *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .           /DLOG(X**2/XLB**2)
     .           +(B1(NF)**2*(DLOG(DLOG(X**2/XLB**2))**2
     .                      -DLOG(DLOG(X**2/XLB**2))-1)+B2(NF))
     .           /DLOG(X**2/XLB**2)**2
     .           -(B1(NF)**3*(DLOG(DLOG(X**2/XLB**2))**3
     .                      -5*DLOG(DLOG(X**2/XLB**2))**2/2
     .                      -2*DLOG(DLOG(X**2/XLB**2))+1/2.d0)
     .            +3*B1(NF)*B2(NF)*DLOG(DLOG(X**2/XLB**2))
     .            -B3(NF)/2)/DLOG(X**2/XLB**2)**3)
      AA(NF)=12D0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      CC(NF)=B2(NF)/AA(NF)
      DD(NF)=B3(NF)/AA(NF)
      XIT(A,B,C,D,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)
     .          *(1-(A*B*(DLOG(X)**2-DLOG(X)-1)+C/B)/X/DLOG(X)
     .          +(A**2*B**2*(DLOG(X)**3-5*DLOG(X)**2/2-2*DLOG(X)+1/2.D0)
     .           +3*A*C*DLOG(X)-D/B/2)/X**2/DLOG(X))))
      PI=4.D0*DATAN(1.D0)
      ZETA2 = PI**2/6
      ZETA3 = 1.2020569031595942853997381D0
      XLB2=XLB
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB2**2)
      IF(NF1.LT.NF2)THEN
       DELTA = 7*ALS4(NF1,Q,XLB1)**2/PI**2/24
     .       + (58933/124416.D0+2*ZETA2/3+2*ZETA2/9*DLOG(2.D0)
     .         +80507*ZETA3/27648)*ALS4(NF1,Q,XLB1)**3/PI**3
       ALP=ALS4(NF1,Q,XLB1)*(1+DELTA)
      ELSE
       DELTA = 7*ALS4(NF1,Q,XLB1)**2/PI**2/24
     .       + (58933/124416.D0+2*ZETA2/3+2*ZETA2/9*DLOG(2.D0)
     .         +80507*ZETA3/27648)*ALS4(NF1,Q,XLB1)**3/PI**3
       ALP=ALS4(NF1,Q,XLB1)/(1+DELTA)
      ENDIF
      A=AA(NF2)/ALP
      B=BB(NF2)*ALP
      C=CC(NF2)*ALP
      D=DD(NF2)*ALP
      XX=XIT(A,B,C,D,X)
      XLB2=Q*DEXP(-XX/2.D0)
      IF(NF1.LT.NF2)THEN
       DELTA = 7*ALS4(NF1,Q,XLB1)**2/PI**2/24
     .       + (58933/124416.D0+2*ZETA2/3+2*ZETA2/9*DLOG(2.D0)
     .         +80507*ZETA3/27648)*ALS4(NF1,Q,XLB1)**3/PI**3
       Y1=ALS4(NF1,Q,XLB1)*(1+DELTA)
       Y2=ALS4(NF2,Q,XLB2)
      ELSE
       DELTA = 7*ALS4(NF1,Q,XLB1)**2/PI**2/24
     .       + (58933/124416.D0+2*ZETA2/3+2*ZETA2/9*DLOG(2.D0)
     .         +80507*ZETA3/27648)*ALS4(NF1,Q,XLB1)**3/PI**3
       Y1=ALS4(NF1,Q,XLB1)/(1+DELTA)
       Y2=ALS4(NF2,Q,XLB2)
      ENDIF
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
       XITER4=XLB2
      RETURN
      END

      DOUBLE PRECISION FUNCTION XITLA(NO,ALP,ACC)
C--ITERATION ROUTINE TO DETERMINE IMPROVED LAMBDAS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/SUSYP/GF,ALPHA,SW2,TGBET,AMQ,AMZ,AMS
      COMMON/ALSN3LO/N3LO
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      B2(NF)=27/2.D0*(2857-5033/9.D0*NF+325/27.D0*NF**2)/B0(NF)**3
      B3(NF)= 81*(149753/6.d0+3564*zeta3-(1078361/162.d0+6508*zeta3/27)
     .        *nf+(50065/162.d0+6472*zeta3/81)*nf**2+1093/729.d0*nf**3)
     .      / B0(NF)**4
      ALS2(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .              *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .              /DLOG(X**2/XLB**2))
      ALS3(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .          *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .           /DLOG(X**2/XLB**2)
     .           +(B1(NF)**2*(DLOG(DLOG(X**2/XLB**2))**2
     .                      -DLOG(DLOG(X**2/XLB**2))-1)+B2(NF))
     .           /DLOG(X**2/XLB**2)**2)
      ALS4(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .          *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .           /DLOG(X**2/XLB**2)
     .           +(B1(NF)**2*(DLOG(DLOG(X**2/XLB**2))**2
     .                      -DLOG(DLOG(X**2/XLB**2))-1)+B2(NF))
     .           /DLOG(X**2/XLB**2)**2
     .           -(B1(NF)**3*(DLOG(DLOG(X**2/XLB**2))**3
     .                      -5*DLOG(DLOG(X**2/XLB**2))**2/2
     .                      -2*DLOG(DLOG(X**2/XLB**2))+1/2.d0)
     .            +3*B1(NF)*B2(NF)*DLOG(DLOG(X**2/XLB**2))
     .            -B3(NF)/2)/DLOG(X**2/XLB**2)**3)
      AA(NF)=12D0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      CC(NF)=B2(NF)/AA(NF)
      DD(NF)=B3(NF)/AA(NF)
      XIT(A,B,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)))
      XIT3(A,B,C,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)
     .          *(1-(A*B*(DLOG(X)**2-DLOG(X)-1)+C/B)/X/DLOG(X))))
      XIT4(A,B,C,D,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)
     .          *(1-(A*B*(DLOG(X)**2-DLOG(X)-1)+C/B)/X/DLOG(X)
     .          +(A**2*B**2*(DLOG(X)**3-5*DLOG(X)**2/2-2*DLOG(X)+1/2.D0)
     .           +3*A*C*DLOG(X)-D/B/2)/X**2/DLOG(X))))
      PI=4.D0*DATAN(1.D0)
      ZETA2 = PI**2/6
      ZETA3 = 1.2020569031595942853997381D0
      NF=5
      Q=AMZ
      XLB=Q*DEXP(-AA(NF)/ALP/2.D0)
      IF(NO.EQ.1)GOTO 111
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB**2)
      A=AA(NF)/ALP
      B=BB(NF)*ALP
      C=CC(NF)*ALP
      D=DD(NF)*ALP
      IF(NO.EQ.2)THEN
       XX=XIT(A,B,X)
      ELSEIF(NO.EQ.3)THEN
       XX=XIT3(A,B,C,X)
      ELSE
       XX=XIT4(A,B,C,D,X)
      ENDIF
      IF(N3LO.NE.0) XX=XIT4(A,B,C,D,X)
      XLB=Q*DEXP(-XX/2.D0)
      Y1=ALP
      IF(NO.EQ.2)THEN
       Y2=ALS2(NF,Q,XLB)
      ELSEIF(NO.EQ.3)THEN
       Y2=ALS3(NF,Q,XLB)
      ELSE
       Y2=ALS4(NF,Q,XLB)
      ENDIF
      IF(N3LO.NE.0) Y2=ALS4(NF,Q,XLB)
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
111   XITLA=XLB
      RETURN
      END

      DOUBLE PRECISION FUNCTION DNF(Q)
C--NUMBER OF LIGHT FLAVORS CONTRIBUTING TO THE QCD BETA-FUNCTION
C--Q = SCALE (FOR SMOOTHING AROUND TOP THRESHOLD)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/OUTPUT/NIN,NOUT
      COMMON/SMOOTH/DLT1,DLT
      COMMON/MASS/AMH,AMQ,S
      AMT=AMQ * 1.D8
      ALOW=AMT*(1.D0-DLT)
      AUP =AMT*(1.D0+DLT)
      IF(DLT.EQ.0)THEN
       CORR=1.D0
      ELSE
       XX=(Q-ALOW)/2.D0/DLT/AMT
       CORR=FSMO(XX)
      ENDIF
      IF(Q.LE.ALOW)THEN
        DNF=5.D0
      ELSE
        DNF=5.D0*(1.D0-CORR)+6.D0*CORR
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION FSMO(X)
C--QUADRATIC POLYNOMIAL FOR SMOOTHING
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      IF(X.LE.0.D0)THEN
        FSMO=0.D0
      ELSEIF(X.LE.0.5D0)THEN
        FSMO=2.D0*X**2
      ELSEIF(X.LE.1.D0)THEN
        FSMO=1.D0-2.D0*(1.D0-X)**2
      ELSE
        FSMO=1.D0
      ENDIF
      RETURN
      END
 
      SUBROUTINE LUMINI
C--INITIALIZATION OF LUMINOSITIES BY READING THE GRIDS
      PARAMETER(NN1=10,NN2=500,NN11=NN1+3,NN22=NN2+3)
      PARAMETER(NIN=5)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      CHARACTER*5 NAME,FNAME
      DIMENSION QQ(NN1),TAU(NN2)
      COMMON/LUMPAR/DLGG(NN2,NN1),DLGQ(NN2,NN1),DLQQ(NN2,NN1),
     .              DX(NN2,2),NX(2)
      COMMON/ALS/XLAMBDA,AMAC,AMAB,AMAT,N0
      COMMON/FACSC/ISCHEME
      COMMON/FILE/NAME
      COMMON/QLIM/QMAX,QMIN
 
      OPEN(NIN,FILE=NAME // '.data')
 
C--READ DATA
 
60    READ(NIN,53)FNAME
      IF(FNAME.NE.NAME)GOTO 60
      BACKSPACE(NIN)
 
      READ(NIN,51)XLAMBDA,ISCHEME
      READ(NIN,52)NQ,NTAU
51    FORMAT(25X,G10.4,19X,I10)
52    FORMAT(25X,I10,19X,I10)
53    FORMAT(1X,A5)
 
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
 
C--READ FACTORIZATION-SCALES
      DO 16 I=1,NQ,5
        READ(NIN,50)QQ(I),QQ(I+1),QQ(I+2),QQ(I+3),QQ(I+4)
16    CONTINUE
      QMAX=QQ(NQ)
      QMIN=QQ(1)
 
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
 
C--READ TAU'S
      DO 17 I=1,NTAU,5
        READ(NIN,50)TAU(I),TAU(I+1),TAU(I+2),TAU(I+3),TAU(I+4)
17    CONTINUE
 
      READ(NIN,*)
      READ(NIN,*)
 
      DO 18 I=1,NQ
      READ(NIN,*)
C--READ L^GG'S
      DO 18 J=1,NTAU,5
        READ(NIN,50) DLGG(J,I),DLGG(J+1,I),DLGG(J+2,I),
     .               DLGG(J+3,I),DLGG(J+4,I)
18      CONTINUE
 
      READ(NIN,*)
      READ(NIN,*)
 
      DO 19 I=1,NQ
      READ(NIN,*)
C--READ L^GQ'S
      DO 19 J=1,NTAU,5
        READ(NIN,50) DLGQ(J,I),DLGQ(J+1,I),DLGQ(J+2,I),
     .               DLGQ(J+3,I),DLGQ(J+4,I)
19      CONTINUE
 
      READ(NIN,*)
      READ(NIN,*)
 
      DO 20 I=1,NQ
      READ(NIN,*)
C--READ L^QQBAR'S
      DO 20 J=1,NTAU,5
        READ(NIN,50) DLQQ(J,I),DLQQ(J+1,I),DLQQ(J+2,I),
     .               DLQQ(J+3,I),DLQQ(J+4,I)
20      CONTINUE
 
50    FORMAT(3X,5(D10.4,3X))
 
      DO 777 I=1,NQ
        DX(I,2)=DLOG(QQ(I))
777   CONTINUE
 
      DO 778 I=1,NTAU
        DX(I,1)=DLOG(TAU(I))
778   CONTINUE
 
      NX(1)=NTAU
      NX(2)=NQ
 
      CLOSE(NIN)
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DLUMGG(TAU,Q)
C--GG-LUMINOSITY
      PARAMETER(NN1=10,NN2=500,NN11=NN1+3,NN22=NN2+3)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(2)
      DIMENSION PDF(-6:6)
      COMMON/INTRO/INTEG
      COMMON/PARTDN/XVAR
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/LUMPAR/DLGG(NN2,NN1),DLGQ(NN2,NN1),DLQQ(NN2,NN1),
     .              DX(NN2,2),NX(2)
      COMMON/FACSHIFT/RHOFAC
      QQ=RHOFAC*Q
      IF(INTEG.EQ.3)THEN
        X1=EPST+(1.D0-2.D0*EPST)*XVAR
        X1=-DLOG(TAU)*X1
        X1=DEXP(-X1)
        X2=TAU/X1
        CALL STRUC(X1,QQ,PDF)
        GG=PDF(0)
        CALL STRUC(X2,QQ,PDF)
        DLUMGG=GG*PDF(0)
        DLUMGG=-DLUMGG*DLOG(TAU)/TAU*(1.D0-2.D0*EPST)
      ELSE
        X(1)=DLOG(TAU)
        X(2)=DLOG(QQ)
C--INTERPOLATION OF THE GRID
        DLUMGG = FINT2(X,DX,DLGG,NX)
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DLUMGQ(TAU,Q)
C--GQ-LUMINOSITY
      PARAMETER(NN1=10,NN2=500,NN11=NN1+3,NN22=NN2+3)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(2)
      DIMENSION PDF(-6:6)
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/INTRO/INTEG
      COMMON/PARTDN/XVAR
      COMMON/LUMPAR/DLGG(NN2,NN1),DLGQ(NN2,NN1),DLQQ(NN2,NN1),
     .              DX(NN2,2),NX(2)
      COMMON/FACSHIFT/RHOFAC
      QQ=RHOFAC*Q
      IF(INTEG.EQ.3)THEN
        X1=EPST+(1.D0-2.D0*EPST)*XVAR
        X1=-DLOG(TAU)*X1
        X1=DEXP(-X1)
        X2=TAU/X1
        CALL STRUC(X1,QQ,PDF)
        G1=PDF(0)
        Q1=0.D0
        DO 1 I=1,6
         Q1=Q1+PDF(I)+PDF(-I)
1       CONTINUE
        CALL STRUC(X2,QQ,PDF)
        G2=PDF(0)
        Q2=0.D0
        DO 2 I=1,6
         Q2=Q2+PDF(I)+PDF(-I)
2       CONTINUE
        DLUMGQ=G1*Q2+Q1*G2
        DLUMGQ=-DLUMGQ*DLOG(TAU)/TAU*(1.D0-2.D0*EPST)
      ELSE
        X(1)=DLOG(TAU)
        X(2)=DLOG(QQ)
C--INTERPOLATION OF THE GRID
        DLUMGQ = FINT2(X,DX,DLGQ,NX)
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DLUMQQB(TAU,Q)
C--QQBAR-LUMINOSITY
      PARAMETER(NN1=10,NN2=500,NN11=NN1+3,NN22=NN2+3)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(2)
      DIMENSION PDF1(-6:6),PDF2(-6:6)
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/INTRO/INTEG
      COMMON/PARTDN/XVAR
      COMMON/LUMPAR/DLGG(NN2,NN1),DLGQ(NN2,NN1),DLQQ(NN2,NN1),
     .              DX(NN2,2),NX(2)
      COMMON/COLLI/ICOLL
      COMMON/FACSHIFT/RHOFAC
      QQ=RHOFAC*Q
      IF(INTEG.EQ.3)THEN
        X1=EPST+(1.D0-2.D0*EPST)*XVAR
        X1=-DLOG(TAU)*X1
        X1=DEXP(-X1)
        X2=TAU/X1
        CALL STRUC(X1,QQ,PDF1)
        CALL STRUC(X2,QQ,PDF2)
        DLUMQQB=0.D0
        DO 1 I=1,6
         DLUMQQB=DLUMQQB+PDF1(-I)*PDF2( ICOLL*I)
     .                +PDF1( I)*PDF2(-ICOLL*I)
1       CONTINUE
        DLUMQQB=-DLUMQQB*DLOG(TAU)/TAU*(1.D0-2.D0*EPST)
      ELSE
        X(1)=DLOG(TAU)
        X(2)=DLOG(QQ)
C--INTERPOLATION OF THE GRID
        DLUMQQB = FINT2(X,DX,DLQQ,NX)
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DLUMQQ(TAU,Q)
C--QQ-LUMINOSITY
      PARAMETER(NN1=10,NN2=500,NN11=NN1+3,NN22=NN2+3)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(2)
      DIMENSION PDF1(-6:6),PDF2(-6:6)
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/INTRO/INTEG
      COMMON/PARTDN/XVAR
      COMMON/LUMPAR/DLGG(NN2,NN1),DLGQ(NN2,NN1),DLQQ(NN2,NN1),
     .              DX(NN2,2),NX(2)
      COMMON/COLLI/ICOLL
      COMMON/FACSHIFT/RHOFAC
      QQ=RHOFAC*Q
      IF(INTEG.EQ.3)THEN
        X1=EPST+(1.D0-2.D0*EPST)*XVAR
        X1=-DLOG(TAU)*X1
        X1=DEXP(-X1)
        X2=TAU/X1
        CALL STRUC(X1,QQ,PDF1)
        CALL STRUC(X2,QQ,PDF2)
        DLUMQQ=0.D0
        DO 1 I=1,6
         DLUMQQ=DLUMQQ+PDF1(I)*PDF2( ICOLL*I)
     .                +PDF1(-I)*PDF2(-ICOLL*I)
1       CONTINUE
        DLUMQQ=-DLUMQQ*DLOG(TAU)/TAU*(1.D0-2.D0*EPST)
      ELSE
        X(1)=DLOG(TAU)
        X(2)=DLOG(QQ)
C--INTERPOLATION OF THE GRID
        DLUMQQ = FINT2(X,DX,DLQQ,NX)
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DLUMQQP(TAU,Q)
C--QQP-LUMINOSITY
      PARAMETER(NN1=10,NN2=500,NN11=NN1+3,NN22=NN2+3)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(2)
      DIMENSION PDF1(-6:6),PDF2(-6:6)
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/INTRO/INTEG
      COMMON/PARTDN/XVAR
      COMMON/LUMPAR/DLGG(NN2,NN1),DLGQ(NN2,NN1),DLQQ(NN2,NN1),
     .              DX(NN2,2),NX(2)
      COMMON/COLLI/ICOLL
      COMMON/FACSHIFT/RHOFAC
      QQ=RHOFAC*Q
      IF(INTEG.EQ.3)THEN
        X1=EPST+(1.D0-2.D0*EPST)*XVAR
        X1=-DLOG(TAU)*X1
        X1=DEXP(-X1)
        X2=TAU/X1
        CALL STRUC(X1,QQ,PDF1)
        CALL STRUC(X2,QQ,PDF2)
        DLUMQQP=0.D0
        DO 1 I=1,6
         DO 1 J=1,6
          IF(I.NE.J)THEN
           DLUMQQP=DLUMQQP+PDF1(I)*PDF2( ICOLL*J)+PDF1(I)*PDF2(-ICOLL*J)
     .                   +PDF1(-I)*PDF2(-ICOLL*J)+PDF1(-I)*PDF2(ICOLL*J)
          ENDIF
1       CONTINUE
        DLUMQQP=-DLUMQQP*DLOG(TAU)/TAU*(1.D0-2.D0*EPST)
      ELSE
        X(1)=DLOG(TAU)
        X(2)=DLOG(QQ)
C--INTERPOLATION OF THE GRID
        DLUMQQP = FINT2(X,DX,DLQQ,NX)
      ENDIF
      RETURN
      END
 
      COMPLEX*16 FUNCTION LI2(X)
C--COMPLEX DILOGARITHM (SPENCE-FUNCTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Y,CLI2
      COMMON/CONST/ZETA2,ZETA3
      LI2 = 0.D0
      ZERO=1.D-40
      XR=DREAL(X)
      XI=DIMAG(X)
      R2=XR*XR+XI*XI
      IF(R2.LE.ZERO)THEN
        LI2=X
        RETURN
      ENDIF
      RR=XR/R2
      IF(R2.EQ.1.D0.AND.XI.EQ.0.D0)THEN
        IF(XR.EQ.1.D0)THEN
          LI2=DCMPLX(ZETA2)
        ELSE
          LI2=-DCMPLX(ZETA2/2.D0)
        ENDIF
        RETURN
      ELSEIF(R2.GT.1.D0.AND.RR.GT.0.5D0)THEN
        Y=(X-1.D0)/X
        LI2=CLI2(Y)+ZETA2-CDLOG(X)*CDLOG(1.D0-X)+0.5D0*CDLOG(X)**2
        RETURN
      ELSEIF(R2.GT.1.D0.AND.RR.LE.0.5D0)THEN
        Y=1.D0/X
        LI2=-CLI2(Y)-ZETA2-0.5D0*CDLOG(-X)**2
        RETURN
      ELSEIF(R2.LE.1.D0.AND.XR.GT.0.5D0)THEN
        Y=1.D0-X
        LI2=-CLI2(Y)+ZETA2-CDLOG(X)*CDLOG(1.D0-X)
       RETURN
      ELSEIF(R2.LE.1.D0.AND.XR.LE.0.5D0)THEN
        Y=X
        LI2=CLI2(Y)
        RETURN
      ENDIF
      END
 
      COMPLEX*16 FUNCTION S12(X)
C--COMPLEX TRILOGARITHM S_12
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Y,CS12,CLI3,CLI2
      COMMON/CONST/ZETA2,ZETA3
      S12 = 0.D0
      ONE=1.D0
      HALF=0.5D0
      ZERO=1.D-40
      XR=DREAL(X)
      XI=DIMAG(X)
      R2=XR*XR+XI*XI
      IF(R2.LE.ZERO)THEN
        S12=X**2/4.D0
        RETURN
      ENDIF
      RR=XR/R2
      IF(R2.EQ.1.D0.AND.XI.EQ.0.D0)THEN
        IF(XR.EQ.1.D0)THEN
          S12=DCMPLX(ZETA3)
        ELSE
          S12=DCMPLX(ZETA3/8.D0)
        ENDIF
        RETURN
      ELSEIF(R2.GT.ONE.AND.RR.GT.HALF)THEN
        Y=(X-1.D0)/X
        S12=CLI3(Y)-CS12(Y)-CDLOG(-Y)*CLI2(Y)+ZETA3+CDLOG(1.D0-Y)
     .      *CDLOG(-Y)*CDLOG((Y-1.D0)/Y)/2.D0
     .      -CDLOG(1.D0-Y)**3/6.D0
        RETURN
      ELSEIF(R2.GT.ONE.AND.RR.LE.HALF)THEN
        Y=1.D0/X
        S12=CLI3(Y)-CS12(Y)-CDLOG(-Y)*CLI2(Y)+ZETA3
     .      -CDLOG(-Y)**3/6.D0
        RETURN
      ELSEIF(R2.LE.ONE.AND.XR.GT.HALF)THEN
        Y=1.D0-X
        S12=-CLI3(Y)+CDLOG(Y)*CLI2(Y)+ZETA3+CDLOG(Y)**2
     .      *CDLOG(1.D0-Y)/2.D0
       RETURN
      ELSEIF(R2.LE.ONE.AND.XR.LE.HALF)THEN
        Y=X
        S12=CS12(Y)
        RETURN
      ENDIF
      END
 
      COMPLEX*16 FUNCTION LI3(X)
C--COMPLEX TRILOGARITHM LI_3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Y,CS12,CLI3,CLI2
      COMMON/CONST/ZETA2,ZETA3
      LI3 = 0.D0
      ONE=1.D0
      HALF=0.5D0
      ZERO=1.D-40
      XR=DREAL(X)
      XI=DIMAG(X)
      R2=XR*XR+XI*XI
      IF(R2.LE.ZERO)THEN
        LI3=X
        RETURN
      ENDIF
      RR=XR/R2
      IF(R2.EQ.1.D0.AND.XI.EQ.0.D0)THEN
        IF(XR.EQ.1.D0)THEN
          LI3=DCMPLX(ZETA3)
        ELSE
          LI3=-DCMPLX(ZETA3*3.D0/4.D0)
        ENDIF
        RETURN
      ELSEIF(R2.GT.ONE.AND.RR.GT.HALF)THEN
        Y=(X-1.D0)/X
        LI3=-CS12(Y)-CDLOG(1.D0-Y)*CLI2(Y)+ZETA3-CDLOG(1.D0-Y)**2
     .  *CDLOG(-Y)/2.D0+CDLOG(1.D0-Y)**3/6.D0-ZETA2*CDLOG(1.D0-Y)
        RETURN
      ELSEIF(R2.GT.ONE.AND.RR.LE.HALF)THEN
        Y=1.D0/X
        LI3=CLI3(Y)+ZETA2*CDLOG(-Y)+CDLOG(-Y)**3/6.D0
        RETURN
      ELSEIF(R2.LE.ONE.AND.XR.GT.HALF)THEN
        Y=1.D0-X
        LI3=-CS12(Y)-CDLOG(1.D0-Y)*CLI2(Y)+ZETA3
     .      -CDLOG(Y)*CDLOG(1.D0-Y)**2/2.D0+ZETA2*CDLOG(1.D0-Y)
       RETURN
      ELSEIF(R2.LE.ONE.AND.XR.LE.HALF)THEN
        Y=X
        LI3=CLI3(Y)
        RETURN
      ENDIF
      END
 
      COMPLEX*16 FUNCTION CLI2(X)
C--TAYLOR-EXPANSION FOR COMPLEX DILOGARITHM (SPENCE-FUNCTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Z
      COMMON/BERNOULLI/B2(18),B12(18),B3(18)
      COMMON/POLY/NBER
      N=NBER-1
      Z=-CDLOG(1.D0-X)
      CLI2=B2(NBER)
      DO 111 I=N,1,-1
        CLI2=Z*CLI2+B2(I)
111   CONTINUE
      CLI2=Z**2*CLI2+Z
      RETURN
      END
 
      COMPLEX*16 FUNCTION CS12(X)
C--TAYLOR-EXPANSION FOR COMPLEX TRILOGARITHM S_12
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Z
      COMMON/BERNOULLI/B2(18),B12(18),B3(18)
      COMMON/POLY/NBER
      N=NBER-1
      Z=-CDLOG(1.D0-X)
      CS12=B12(NBER)
      DO 111 I=N,1,-1
        CS12=Z*CS12+B12(I)
111   CONTINUE
      CS12=(Z*CS12+0.25D0)*Z**2
      RETURN
      END
 
      COMPLEX*16 FUNCTION CLI3(X)
C--TAYLOR-EXPANSION FOR COMPLEX TRILOGARITHM LI_3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Z
      COMMON/BERNOULLI/B2(18),B12(18),B3(18)
      COMMON/POLY/NBER
      N=NBER-1
      Z=-CDLOG(1.D0-X)
      CLI3=B3(NBER)
      DO 123 I=N,1,-1
        CLI3=Z*CLI3+B3(I)
123   CONTINUE
      CLI3=CLI3*Z**2+Z
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION FACULT(N)
C--DOUBLE PRECISION VERSION OF FACULTY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      FACULT=1.D0
      IF(N.EQ.0)RETURN
      DO 999 I=1,N
        FACULT=FACULT*DFLOAT(I)
999   CONTINUE
      RETURN
      END
 
      SUBROUTINE BERNINI(N)
C--INITIALIZATION OF COEFFICIENTS FOR POLYLOGARITHMS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(18),PB(19)
      COMMON/BERNOULLI/B2(18),B12(18),B3(18)
      COMMON/CONST/ZETA2,ZETA3
      COMMON/POLY/NBER
 
      NBER=N
      PI=4.D0*DATAN(1.D0)
 
      B(1)=-1.D0/2.D0
      B(2)=1.D0/6.D0
      B(3)=0.D0
      B(4)=-1.D0/30.D0
      B(5)=0.D0
      B(6)=1.D0/42.D0
      B(7)=0.D0
      B(8)=-1.D0/30.D0
      B(9)=0.D0
      B(10)=5.D0/66.D0
      B(11)=0.D0
      B(12)=-691.D0/2730.D0
      B(13)=0.D0
      B(14)=7.D0/6.D0
      B(15)=0.D0
      B(16)=-3617.D0/510.D0
      B(17)=0.D0
      B(18)=43867.D0/798.D0
      ZETA2=PI**2/6.D0
      ZETA3=1.202056903159594D0
 
      DO 995 I=1,18
        B2(I)=B(I)/FACULT(I+1)
        B12(I)=DFLOAT(I+1)/FACULT(I+2)*B(I)/2.D0
        PB(I+1)=B(I)
        B3(I)=0.D0
995   CONTINUE
      PB(1)=1.D0
      DO 996 I=1,18
      DO 996 J=0,I
        B3(I)=B3(I)+PB(J+1)*PB(I-J+1)/FACULT(I-J)/FACULT(J+1)
     .                                            /DFLOAT(I+1)
996   CONTINUE
 
      RETURN
      END
 
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE SUSYCP0(ICASE,AMA,AMH,CPT,CPB)
C--CALCULATION OF SUSY-COUPLINGS FROM PSEUDOSCALAR HIGGS MASS
C--AND TAN(BETA)
C--AMH: WANTED HIGGS MASS
C--CPT: SUSY-COUPLING OF HIGGS TO TOP QUARK
C--CPB: SUSY-COUPLING OF HIGGS TO BOTTOM QUARK
C--ICASE=1: AMH=PSEUDOSCALAR MASS
C--ICASE=2: AMH=LIGHT SCALAR MASS
C--ICASE=3: AMH=HEAVY SCALAR MASS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/SUSYP/GF,ALPHA,SW2,TGBET,AMQ,AMZ,AMS
      PI=4.D0*DATAN(1.D0)
      CW2=1.D0-SW2
      BET=DATAN(TGBET)
      SUSYEPS=3.D0/2.D0*ALPHA/PI/SW2/CW2/DSIN(BET)**2*AMQ**4/AMZ**2
     .             *DLOG(1.D0+AMS**2/AMQ**2)
      DIS=(AMA**2+AMZ**2+SUSYEPS)**2    
     .    -4.D0*AMA**2*AMZ**2*DCOS(2.D0*BET)**2
     .    -4.D0*SUSYEPS*(AMA**2*DSIN(BET)**2+AMZ**2*DCOS(BET)**2)
      IF(DIS.GE.0.D0)THEN
        DELTA=DSQRT(DIS)
        AMHL=DSQRT((AMA**2+AMZ**2+SUSYEPS-DELTA)/2.D0)
        AMHH=DSQRT((AMA**2+AMZ**2+SUSYEPS+DELTA)/2.D0)
        ALP=DATAN(DTAN(2.D0*BET)*(AMA**2+AMZ**2)/(AMA**2-AMZ**2+
     .      SUSYEPS/DCOS(2.D0*BET)))/2.D0
        IF(ALP.GT.0.D0)ALP=ALP-PI/2.D0
        IF(ICASE.EQ.1)THEN
          AMH=AMA
          CPT=1.D0/TGBET
          CPB=TGBET
        ELSEIF(ICASE.EQ.2)THEN
          AMH=AMHL
          CPT=DCOS(ALP)/DSIN(BET)
          CPB=-DSIN(ALP)/DCOS(BET)
        ELSE
          AMH=AMHH
          CPT=DSIN(ALP)/DSIN(BET)
          CPB=DCOS(ALP)/DCOS(BET)
        ENDIF
      ELSE
        AMH=0.D0
        CPT=0.D0
        CPB=0.D0
      ENDIF
      RETURN
      END
 
      SUBROUTINE AMHAMA0(ICASE,AMH,AMA)
C--CALCULATION OF PSEUDOSCALAR HIGGS MASS FROM HIGGS MASS AMH
C--ICASE=1: AMH=PSEUDOSCALAR MASS
C--ICASE=2: AMH=LIGHT SCALAR MASS
C--ICASE=3: AMH=HEAVY SCALAR MASS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/SUSYP/GF,ALPHA,SW2,TGBET,AMQ,AMZ,AMS
      IF(ICASE.NE.1)THEN
       PI=4.D0*DATAN(1.D0)
       CW2=1.D0-SW2
       BET=DATAN(TGBET)
       SUSYEPS=3.D0/2.D0*ALPHA/PI/SW2/CW2/DSIN(BET)**2*AMQ**4/AMZ**2
     .             *DLOG(1.D0+AMS**2/AMQ**2)
       IF(ICASE.EQ.2)THEN
         DIS=AMZ**2*DCOS(2.D0*BET)**2+SUSYEPS*DSIN(BET)**2
         DIS=DIS-AMH**2
       ELSE
         DIS=(AMZ**2+SUSYEPS+DSQRT((AMZ**2+SUSYEPS)**2
     .      -4.D0*SUSYEPS*AMZ**2*DCOS(BET)**2))/2.D0
         DIS=AMH**2-DIS
       ENDIF
       IF(DIS.GE.0.D0)THEN
         AMA=AMH**2*(AMZ**2-AMH**2+SUSYEPS)-SUSYEPS*AMZ**2*DCOS(BET)**2
         AMA=AMA/(AMZ**2*DCOS(2.D0*BET)**2-AMH**2+SUSYEPS*DSIN(BET)**2)
         AMA=DSQRT(AMA)
       ELSEIF(ICASE.EQ.2)THEN
         AMA=1.D8
       ELSE
         AMA=0.D0
       ENDIF
      ELSE
       AMA=AMH
      ENDIF
      RETURN
      END
 
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      BLOCKDATA CVIRT
C--DATA OF VIRTUAL CORRECTIONS TO QUARK LOOPS
      PARAMETER(NN=999)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/RESINT/XX(NN),YHR(NN),YHI(NN),YAR(NN),YAI(NN)

      DATA ( XX(I),I=  1, 40)/
     .   0.000000D+00,   0.100000D+00,   0.200000D+00,   0.300000D+00,
     .   0.400000D+00,   0.500000D+00,   0.600000D+00,   0.800000D+00,
     .   0.900000D+00,   0.100000D+01,   0.110000D+01,   0.120000D+01,
     .   0.130000D+01,   0.150000D+01,   0.160000D+01,   0.170000D+01,
     .   0.180000D+01,   0.190000D+01,   0.200000D+01,   0.210000D+01,
     .   0.220000D+01,   0.230000D+01,   0.240000D+01,   0.250000D+01,
     .   0.260000D+01,   0.270000D+01,   0.280000D+01,   0.300000D+01,
     .   0.310000D+01,   0.320000D+01,   0.330000D+01,   0.340000D+01,
     .   0.350000D+01,   0.360000D+01,   0.370000D+01,   0.380000D+01,
     .   0.390000D+01,   0.391000D+01,   0.392000D+01,   0.393000D+01/

      DATA ( XX(I),I= 41, 80)/
     .   0.394000D+01,   0.395000D+01,   0.396000D+01,   0.397000D+01,
     .   0.398000D+01,   0.399000D+01,   0.399900D+01,   0.399990D+01,
     .   0.399999D+01,   0.400001D+01,   0.400010D+01,   0.400100D+01,
     .   0.401000D+01,   0.402000D+01,   0.403000D+01,   0.404000D+01,
     .   0.405000D+01,   0.406000D+01,   0.407000D+01,   0.408000D+01,
     .   0.409000D+01,   0.410000D+01,   0.420000D+01,   0.430000D+01,
     .   0.440000D+01,   0.450000D+01,   0.460000D+01,   0.470000D+01,
     .   0.480000D+01,   0.490000D+01,   0.500000D+01,   0.510000D+01,
     .   0.520000D+01,   0.530000D+01,   0.540000D+01,   0.550000D+01,
     .   0.560000D+01,   0.570000D+01,   0.580000D+01,   0.590000D+01/

      DATA ( XX(I),I= 81,120)/
     .   0.600000D+01,   0.610000D+01,   0.620000D+01,   0.630000D+01,
     .   0.640000D+01,   0.650000D+01,   0.660000D+01,   0.670000D+01,
     .   0.680000D+01,   0.690000D+01,   0.700000D+01,   0.710000D+01,
     .   0.720000D+01,   0.730000D+01,   0.740000D+01,   0.750000D+01,
     .   0.760000D+01,   0.770000D+01,   0.780000D+01,   0.790000D+01,
     .   0.800000D+01,   0.810000D+01,   0.820000D+01,   0.830000D+01,
     .   0.840000D+01,   0.850000D+01,   0.860000D+01,   0.870000D+01,
     .   0.880000D+01,   0.890000D+01,   0.900000D+01,   0.910000D+01,
     .   0.920000D+01,   0.930000D+01,   0.940000D+01,   0.950000D+01,
     .   0.960000D+01,   0.970000D+01,   0.980000D+01,   0.990000D+01/

      DATA ( XX(I),I=121,160)/
     .   0.100000D+02,   0.105000D+02,   0.110000D+02,   0.115000D+02,
     .   0.120000D+02,   0.125000D+02,   0.130000D+02,   0.135000D+02,
     .   0.140000D+02,   0.145000D+02,   0.150000D+02,   0.155000D+02,
     .   0.160000D+02,   0.165000D+02,   0.170000D+02,   0.175000D+02,
     .   0.180000D+02,   0.185000D+02,   0.190000D+02,   0.195000D+02,
     .   0.200000D+02,   0.210000D+02,   0.220000D+02,   0.230000D+02,
     .   0.240000D+02,   0.250000D+02,   0.260000D+02,   0.270000D+02,
     .   0.280000D+02,   0.290000D+02,   0.300000D+02,   0.310000D+02,
     .   0.320000D+02,   0.330000D+02,   0.340000D+02,   0.350000D+02,
     .   0.360000D+02,   0.370000D+02,   0.380000D+02,   0.390000D+02/

      DATA ( XX(I),I=161,200)/
     .   0.400000D+02,   0.410000D+02,   0.420000D+02,   0.430000D+02,
     .   0.440000D+02,   0.450000D+02,   0.460000D+02,   0.470000D+02,
     .   0.480000D+02,   0.490000D+02,   0.500000D+02,   0.510000D+02,
     .   0.520000D+02,   0.530000D+02,   0.540000D+02,   0.550000D+02,
     .   0.560000D+02,   0.570000D+02,   0.580000D+02,   0.590000D+02,
     .   0.600000D+02,   0.610000D+02,   0.620000D+02,   0.630000D+02,
     .   0.640000D+02,   0.650000D+02,   0.660000D+02,   0.670000D+02,
     .   0.680000D+02,   0.690000D+02,   0.700000D+02,   0.710000D+02,
     .   0.720000D+02,   0.730000D+02,   0.740000D+02,   0.750000D+02,
     .   0.760000D+02,   0.770000D+02,   0.780000D+02,   0.790000D+02/

      DATA ( XX(I),I=201,240)/
     .   0.800000D+02,   0.810000D+02,   0.820000D+02,   0.830000D+02,
     .   0.840000D+02,   0.850000D+02,   0.860000D+02,   0.870000D+02,
     .   0.880000D+02,   0.890000D+02,   0.900000D+02,   0.910000D+02,
     .   0.920000D+02,   0.930000D+02,   0.940000D+02,   0.950000D+02,
     .   0.960000D+02,   0.970000D+02,   0.980000D+02,   0.990000D+02,
     .   0.100000D+03,   0.105000D+03,   0.110000D+03,   0.115000D+03,
     .   0.120000D+03,   0.125000D+03,   0.130000D+03,   0.135000D+03,
     .   0.140000D+03,   0.145000D+03,   0.150000D+03,   0.155000D+03,
     .   0.160000D+03,   0.165000D+03,   0.170000D+03,   0.175000D+03,
     .   0.180000D+03,   0.185000D+03,   0.190000D+03,   0.195000D+03/

      DATA ( XX(I),I=241,280)/
     .   0.200000D+03,   0.210000D+03,   0.220000D+03,   0.230000D+03,
     .   0.240000D+03,   0.250000D+03,   0.260000D+03,   0.270000D+03,
     .   0.280000D+03,   0.290000D+03,   0.300000D+03,   0.310000D+03,
     .   0.320000D+03,   0.330000D+03,   0.340000D+03,   0.350000D+03,
     .   0.360000D+03,   0.370000D+03,   0.380000D+03,   0.390000D+03,
     .   0.400000D+03,   0.410000D+03,   0.420000D+03,   0.430000D+03,
     .   0.440000D+03,   0.450000D+03,   0.460000D+03,   0.470000D+03,
     .   0.480000D+03,   0.490000D+03,   0.500000D+03,   0.550000D+03,
     .   0.600000D+03,   0.650000D+03,   0.700000D+03,   0.750000D+03,
     .   0.800000D+03,   0.850000D+03,   0.900000D+03,   0.950000D+03/

      DATA ( XX(I),I=281,320)/
     .   0.100000D+04,   0.105000D+04,   0.110000D+04,   0.115000D+04,
     .   0.120000D+04,   0.125000D+04,   0.130000D+04,   0.135000D+04,
     .   0.140000D+04,   0.145000D+04,   0.150000D+04,   0.155000D+04,
     .   0.160000D+04,   0.165000D+04,   0.170000D+04,   0.175000D+04,
     .   0.180000D+04,   0.185000D+04,   0.190000D+04,   0.195000D+04,
     .   0.200000D+04,   0.210000D+04,   0.220000D+04,   0.230000D+04,
     .   0.240000D+04,   0.250000D+04,   0.260000D+04,   0.270000D+04,
     .   0.280000D+04,   0.290000D+04,   0.300000D+04,   0.310000D+04,
     .   0.320000D+04,   0.330000D+04,   0.340000D+04,   0.350000D+04,
     .   0.360000D+04,   0.370000D+04,   0.380000D+04,   0.390000D+04/

      DATA ( XX(I),I=321,360)/
     .   0.400000D+04,   0.410000D+04,   0.420000D+04,   0.430000D+04,
     .   0.440000D+04,   0.450000D+04,   0.460000D+04,   0.470000D+04,
     .   0.480000D+04,   0.490000D+04,   0.500000D+04,   0.550000D+04,
     .   0.600000D+04,   0.650000D+04,   0.700000D+04,   0.750000D+04,
     .   0.800000D+04,   0.850000D+04,   0.900000D+04,   0.950000D+04,
     .   0.100000D+05,   0.105000D+05,   0.110000D+05,   0.115000D+05,
     .   0.120000D+05,   0.125000D+05,   0.130000D+05,   0.135000D+05,
     .   0.140000D+05,   0.145000D+05,   0.150000D+05,   0.155000D+05,
     .   0.160000D+05,   0.165000D+05,   0.170000D+05,   0.175000D+05,
     .   0.180000D+05,   0.185000D+05,   0.190000D+05,   0.195000D+05/

      DATA ( XX(I),I=361,400)/
     .   0.200000D+05,   0.210000D+05,   0.220000D+05,   0.230000D+05,
     .   0.240000D+05,   0.250000D+05,   0.260000D+05,   0.270000D+05,
     .   0.280000D+05,   0.290000D+05,   0.300000D+05,   0.310000D+05,
     .   0.320000D+05,   0.330000D+05,   0.340000D+05,   0.350000D+05,
     .   0.360000D+05,   0.370000D+05,   0.380000D+05,   0.390000D+05,
     .   0.400000D+05,   0.410000D+05,   0.420000D+05,   0.430000D+05,
     .   0.440000D+05,   0.450000D+05,   0.460000D+05,   0.470000D+05,
     .   0.480000D+05,   0.490000D+05,   0.500000D+05,   0.550000D+05,
     .   0.600000D+05,   0.650000D+05,   0.700000D+05,   0.750000D+05,
     .   0.800000D+05,   0.850000D+05,   0.900000D+05,   0.950000D+05/

      DATA ( XX(I),I=401,440)/
     .   0.100000D+06,   0.105000D+06,   0.110000D+06,   0.115000D+06,
     .   0.120000D+06,   0.125000D+06,   0.130000D+06,   0.135000D+06,
     .   0.140000D+06,   0.145000D+06,   0.150000D+06,   0.155000D+06,
     .   0.160000D+06,   0.165000D+06,   0.170000D+06,   0.175000D+06,
     .   0.180000D+06,   0.185000D+06,   0.190000D+06,   0.195000D+06,
     .   0.200000D+06,   0.210000D+06,   0.220000D+06,   0.230000D+06,
     .   0.240000D+06,   0.250000D+06,   0.260000D+06,   0.270000D+06,
     .   0.280000D+06,   0.290000D+06,   0.300000D+06,   0.310000D+06,
     .   0.320000D+06,   0.330000D+06,   0.340000D+06,   0.350000D+06,
     .   0.360000D+06,   0.370000D+06,   0.380000D+06,   0.390000D+06/

      DATA ( XX(I),I=441,461)/
     .   0.400000D+06,   0.410000D+06,   0.420000D+06,   0.430000D+06,
     .   0.440000D+06,   0.450000D+06,   0.460000D+06,   0.470000D+06,
     .   0.480000D+06,   0.490000D+06,   0.500000D+06,   0.550000D+06,
     .   0.600000D+06,   0.650000D+06,   0.700000D+06,   0.750000D+06,
     .   0.800000D+06,   0.850000D+06,   0.900000D+06,   0.950000D+06,
     .   0.100000D+07/

C--SCALAR HIGGS

      DATA (YHR(I),I=  1, 40)/
     .   0.153696D+02,   0.153951D+02,   0.154213D+02,   0.154481D+02,
     .   0.154757D+02,   0.155040D+02,   0.155332D+02,   0.155940D+02,
     .   0.156258D+02,   0.156586D+02,   0.156925D+02,   0.157275D+02,
     .   0.157637D+02,   0.158400D+02,   0.158803D+02,   0.159222D+02,
     .   0.159657D+02,   0.160111D+02,   0.160585D+02,   0.161081D+02,
     .   0.161600D+02,   0.162145D+02,   0.162719D+02,   0.163325D+02,
     .   0.163966D+02,   0.164646D+02,   0.165372D+02,   0.166985D+02,
     .   0.167890D+02,   0.168877D+02,   0.169963D+02,   0.171170D+02,
     .   0.172533D+02,   0.174100D+02,   0.175954D+02,   0.178248D+02,
     .   0.181347D+02,   0.181734D+02,   0.182142D+02,   0.182575D+02/

      DATA (YHR(I),I= 41, 80)/
     .   0.183037D+02,   0.183535D+02,   0.184078D+02,   0.184678D+02,
     .   0.185361D+02,   0.186179D+02,   0.187203D+02,   0.187370D+02,
     .   0.187397D+02,   0.187397D+02,   0.187407D+02,   0.187583D+02,
     .   0.188505D+02,   0.189202D+02,   0.189718D+02,   0.190124D+02,
     .   0.190449D+02,   0.190712D+02,   0.190923D+02,   0.191094D+02,
     .   0.191230D+02,   0.191343D+02,   0.191324D+02,   0.190340D+02,
     .   0.188968D+02,   0.187436D+02,   0.185854D+02,   0.184272D+02,
     .   0.182726D+02,   0.181223D+02,   0.179773D+02,   0.178379D+02,
     .   0.177044D+02,   0.175755D+02,   0.174523D+02,   0.173340D+02,
     .   0.172205D+02,   0.171119D+02,   0.170069D+02,   0.169060D+02/

      DATA (YHR(I),I= 81,120)/
     .   0.168092D+02,   0.167155D+02,   0.166296D+02,   0.165386D+02,
     .   0.164551D+02,   0.163742D+02,   0.162951D+02,   0.162194D+02,
     .   0.161458D+02,   0.160745D+02,   0.160057D+02,   0.159387D+02,
     .   0.158736D+02,   0.158106D+02,   0.157493D+02,   0.156897D+02,
     .   0.156310D+02,   0.155753D+02,   0.155202D+02,   0.154678D+02,
     .   0.154148D+02,   0.153639D+02,   0.153141D+02,   0.152655D+02,
     .   0.152188D+02,   0.151723D+02,   0.151270D+02,   0.150829D+02,
     .   0.150399D+02,   0.149983D+02,   0.149571D+02,   0.149162D+02,
     .   0.148769D+02,   0.148380D+02,   0.148001D+02,   0.147630D+02,
     .   0.147268D+02,   0.146917D+02,   0.146564D+02,   0.146222D+02/

      DATA (YHR(I),I=121,160)/
     .   0.145896D+02,   0.144313D+02,   0.142862D+02,   0.141531D+02,
     .   0.140304D+02,   0.139163D+02,   0.138104D+02,   0.137133D+02,
     .   0.136199D+02,   0.135336D+02,   0.134524D+02,   0.133772D+02,
     .   0.133040D+02,   0.132358D+02,   0.131717D+02,   0.131100D+02,
     .   0.130517D+02,   0.129963D+02,   0.129435D+02,   0.128938D+02,
     .   0.128449D+02,   0.127547D+02,   0.126720D+02,   0.125952D+02,
     .   0.125243D+02,   0.124586D+02,   0.123971D+02,   0.123397D+02,
     .   0.122863D+02,   0.122363D+02,   0.121882D+02,   0.121435D+02,
     .   0.121018D+02,   0.120594D+02,   0.120235D+02,   0.119877D+02,
     .   0.119535D+02,   0.119214D+02,   0.118901D+02,   0.118606D+02/

      DATA (YHR(I),I=161,200)/
     .   0.118324D+02,   0.118054D+02,   0.117794D+02,   0.117550D+02,
     .   0.117312D+02,   0.117086D+02,   0.116866D+02,   0.116658D+02,
     .   0.116456D+02,   0.116262D+02,   0.116074D+02,   0.115898D+02,
     .   0.115708D+02,   0.115553D+02,   0.115390D+02,   0.115234D+02,
     .   0.115069D+02,   0.114936D+02,   0.114795D+02,   0.114658D+02,
     .   0.114524D+02,   0.114397D+02,   0.114271D+02,   0.114151D+02,
     .   0.114035D+02,   0.113952D+02,   0.113809D+02,   0.113701D+02,
     .   0.113598D+02,   0.113496D+02,   0.113398D+02,   0.113304D+02,
     .   0.113212D+02,   0.113118D+02,   0.113029D+02,   0.112944D+02,
     .   0.112860D+02,   0.112780D+02,   0.112699D+02,   0.112621D+02/

      DATA (YHR(I),I=201,240)/
     .   0.112545D+02,   0.112473D+02,   0.112401D+02,   0.112333D+02,
     .   0.112262D+02,   0.112198D+02,   0.112129D+02,   0.112068D+02,
     .   0.112006D+02,   0.111946D+02,   0.111880D+02,   0.111824D+02,
     .   0.111771D+02,   0.111715D+02,   0.111657D+02,   0.111610D+02,
     .   0.111555D+02,   0.111508D+02,   0.111457D+02,   0.111411D+02,
     .   0.111356D+02,   0.111142D+02,   0.110943D+02,   0.110760D+02,
     .   0.110602D+02,   0.110450D+02,   0.110321D+02,   0.110202D+02,
     .   0.110092D+02,   0.109992D+02,   0.109884D+02,   0.109818D+02,
     .   0.109735D+02,   0.109674D+02,   0.109598D+02,   0.109557D+02,
     .   0.109502D+02,   0.109457D+02,   0.109417D+02,   0.109376D+02/

      DATA (YHR(I),I=241,280)/
     .   0.109344D+02,   0.109287D+02,   0.109243D+02,   0.109206D+02,
     .   0.109176D+02,   0.109163D+02,   0.109149D+02,   0.109149D+02,
     .   0.109138D+02,   0.109155D+02,   0.109168D+02,   0.109168D+02,
     .   0.109197D+02,   0.109221D+02,   0.109066D+02,   0.109275D+02,
     .   0.109301D+02,   0.109329D+02,   0.109365D+02,   0.109379D+02,
     .   0.109420D+02,   0.109399D+02,   0.109506D+02,   0.109512D+02,
     .   0.109587D+02,   0.109592D+02,   0.109676D+02,   0.109799D+02,
     .   0.109754D+02,   0.109980D+02,   0.109854D+02,   0.110086D+02,
     .   0.110331D+02,   0.110515D+02,   0.110831D+02,   0.111067D+02,
     .   0.111137D+02,   0.111537D+02,   0.111767D+02,   0.111985D+02/

      DATA (YHR(I),I=281,320)/
     .   0.112244D+02,   0.112580D+02,   0.112723D+02,   0.112939D+02,
     .   0.113157D+02,   0.113372D+02,   0.113565D+02,   0.113760D+02,
     .   0.113955D+02,   0.114196D+02,   0.114359D+02,   0.114621D+02,
     .   0.114624D+02,   0.114949D+02,   0.115161D+02,   0.115340D+02,
     .   0.115462D+02,   0.115610D+02,   0.115418D+02,   0.116090D+02,
     .   0.116234D+02,   0.116532D+02,   0.116915D+02,   0.117046D+02,
     .   0.117561D+02,   0.117793D+02,   0.117998D+02,   0.118404D+02,
     .   0.118722D+02,   0.119008D+02,   0.119270D+02,   0.119539D+02,
     .   0.119749D+02,   0.120013D+02,   0.120428D+02,   0.120544D+02,
     .   0.120795D+02,   0.120898D+02,   0.121223D+02,   0.121456D+02/

      DATA (YHR(I),I=321,360)/
     .   0.121634D+02,   0.121898D+02,   0.122172D+02,   0.122386D+02,
     .   0.122635D+02,   0.122877D+02,   0.122954D+02,   0.123050D+02,
     .   0.123233D+02,   0.123454D+02,   0.123801D+02,   0.124654D+02,
     .   0.125601D+02,   0.126519D+02,   0.127199D+02,   0.127902D+02,
     .   0.128598D+02,   0.129767D+02,   0.129954D+02,   0.130571D+02,
     .   0.131151D+02,   0.131736D+02,   0.132275D+02,   0.132820D+02,
     .   0.133335D+02,   0.133835D+02,   0.134319D+02,   0.134789D+02,
     .   0.135246D+02,   0.135691D+02,   0.136124D+02,   0.136546D+02,
     .   0.136957D+02,   0.137359D+02,   0.137751D+02,   0.138134D+02,
     .   0.138510D+02,   0.138876D+02,   0.139235D+02,   0.139586D+02/

      DATA (YHR(I),I=361,400)/
     .   0.139926D+02,   0.140593D+02,   0.141245D+02,   0.141861D+02,
     .   0.142468D+02,   0.143049D+02,   0.143612D+02,   0.144156D+02,
     .   0.144687D+02,   0.145202D+02,   0.145702D+02,   0.146190D+02,
     .   0.146665D+02,   0.147128D+02,   0.147579D+02,   0.148022D+02,
     .   0.148453D+02,   0.148872D+02,   0.149285D+02,   0.149688D+02,
     .   0.150079D+02,   0.150466D+02,   0.150850D+02,   0.151222D+02,
     .   0.151583D+02,   0.151941D+02,   0.152293D+02,   0.152639D+02,
     .   0.152982D+02,   0.153316D+02,   0.153644D+02,   0.155208D+02,
     .   0.156659D+02,   0.158012D+02,   0.159277D+02,   0.160477D+02,
     .   0.161608D+02,   0.162674D+02,   0.163639D+02,   0.164615D+02/

      DATA (YHR(I),I=401,440)/
     .   0.165507D+02,   0.167323D+02,   0.168185D+02,   0.169015D+02,
     .   0.169834D+02,   0.170585D+02,   0.171330D+02,   0.172052D+02,
     .   0.172751D+02,   0.173429D+02,   0.174086D+02,   0.174726D+02,
     .   0.175348D+02,   0.175954D+02,   0.176544D+02,   0.177120D+02,
     .   0.177682D+02,   0.178231D+02,   0.178766D+02,   0.179291D+02,
     .   0.179803D+02,   0.180795D+02,   0.181747D+02,   0.182662D+02,
     .   0.183545D+02,   0.184396D+02,   0.185217D+02,   0.186010D+02,
     .   0.186780D+02,   0.187524D+02,   0.188248D+02,   0.188951D+02,
     .   0.189635D+02,   0.190307D+02,   0.190948D+02,   0.191581D+02,
     .   0.192195D+02,   0.192796D+02,   0.193383D+02,   0.193956D+02/

      DATA (YHR(I),I=441,461)/
     .   0.194517D+02,   0.195066D+02,   0.195604D+02,   0.196129D+02,
     .   0.196644D+02,   0.197149D+02,   0.197643D+02,   0.198131D+02,
     .   0.198606D+02,   0.199076D+02,   0.199535D+02,   0.201718D+02,
     .   0.203734D+02,   0.205605D+02,   0.207353D+02,   0.208990D+02,
     .   0.210539D+02,   0.212002D+02,   0.213389D+02,   0.214705D+02,
     .   0.215521D+02/

      DATA (YHI(I),I=  1, 40)/
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00/

      DATA (YHI(I),I= 41, 80)/
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.103828D-03,   0.956232D-03,   0.101108D-01,
     .   0.939720D-01,   0.184324D+00,   0.270111D+00,   0.351237D+00,
     .   0.429877D+00,   0.506050D+00,   0.578718D+00,   0.649415D+00,
     .   0.714251D+00,   0.796202D+00,   0.134408D+01,   0.177540D+01,
     .   0.211827D+01,   0.239741D+01,   0.262849D+01,   0.282209D+01,
     .   0.298629D+01,   0.312731D+01,   0.324839D+01,   0.335408D+01,
     .   0.344691D+01,   0.352697D+01,   0.359781D+01,   0.366071D+01,
     .   0.371635D+01,   0.376661D+01,   0.380990D+01,   0.384872D+01/

      DATA (YHI(I),I= 81,120)/
     .   0.388367D+01,   0.391480D+01,   0.394363D+01,   0.396735D+01,
     .   0.399014D+01,   0.401029D+01,   0.402624D+01,   0.404140D+01,
     .   0.405531D+01,   0.406688D+01,   0.407779D+01,   0.408696D+01,
     .   0.409483D+01,   0.410147D+01,   0.410726D+01,   0.411202D+01,
     .   0.411582D+01,   0.411879D+01,   0.412116D+01,   0.412387D+01,
     .   0.412368D+01,   0.412350D+01,   0.412368D+01,   0.412299D+01,
     .   0.412206D+01,   0.411998D+01,   0.411814D+01,   0.411566D+01,
     .   0.411291D+01,   0.411029D+01,   0.410560D+01,   0.410255D+01,
     .   0.409883D+01,   0.409429D+01,   0.409001D+01,   0.408554D+01,
     .   0.408035D+01,   0.407563D+01,   0.407019D+01,   0.406465D+01/

      DATA (YHI(I),I=121,160)/
     .   0.405809D+01,   0.403029D+01,   0.399830D+01,   0.396483D+01,
     .   0.392962D+01,   0.389474D+01,   0.385889D+01,   0.382350D+01,
     .   0.378746D+01,   0.375179D+01,   0.371653D+01,   0.368151D+01,
     .   0.364692D+01,   0.361317D+01,   0.357941D+01,   0.354669D+01,
     .   0.351403D+01,   0.348225D+01,   0.345082D+01,   0.341898D+01,
     .   0.338976D+01,   0.333067D+01,   0.327372D+01,   0.321880D+01,
     .   0.316576D+01,   0.311447D+01,   0.306456D+01,   0.301680D+01,
     .   0.296973D+01,   0.292531D+01,   0.288199D+01,   0.283993D+01,
     .   0.279879D+01,   0.277387D+01,   0.272029D+01,   0.268286D+01,
     .   0.264630D+01,   0.261076D+01,   0.257609D+01,   0.254231D+01/

      DATA (YHI(I),I=161,200)/
     .   0.250939D+01,   0.247733D+01,   0.244452D+01,   0.241537D+01,
     .   0.238541D+01,   0.235633D+01,   0.232778D+01,   0.229986D+01,
     .   0.227246D+01,   0.224575D+01,   0.221951D+01,   0.219374D+01,
     .   0.216869D+01,   0.214421D+01,   0.211987D+01,   0.209618D+01,
     .   0.206976D+01,   0.204911D+01,   0.202716D+01,   0.200561D+01,
     .   0.198399D+01,   0.196272D+01,   0.194127D+01,   0.192087D+01,
     .   0.190010D+01,   0.189543D+01,   0.186138D+01,   0.184241D+01,
     .   0.182347D+01,   0.180182D+01,   0.178640D+01,   0.176840D+01,
     .   0.175191D+01,   0.173270D+01,   0.171510D+01,   0.169858D+01,
     .   0.168164D+01,   0.166512D+01,   0.164879D+01,   0.163254D+01/

      DATA (YHI(I),I=201,240)/
     .   0.161641D+01,   0.160101D+01,   0.158550D+01,   0.157001D+01,
     .   0.155500D+01,   0.153987D+01,   0.152506D+01,   0.151075D+01,
     .   0.149651D+01,   0.148231D+01,   0.146604D+01,   0.145443D+01,
     .   0.144050D+01,   0.142725D+01,   0.141303D+01,   0.140051D+01,
     .   0.138766D+01,   0.137455D+01,   0.136181D+01,   0.134906D+01,
     .   0.133671D+01,   0.127592D+01,   0.121827D+01,   0.116199D+01,
     .   0.111162D+01,   0.105920D+01,   0.101241D+01,   0.967089D+00,
     .   0.923058D+00,   0.880397D+00,   0.837579D+00,   0.800124D+00,
     .   0.760050D+00,   0.724849D+00,   0.685016D+00,   0.655055D+00,
     .   0.620992D+00,   0.588618D+00,   0.557638D+00,   0.526851D+00/

      DATA (YHI(I),I=241,280)/
     .   0.497112D+00,   0.439665D+00,   0.385464D+00,   0.333219D+00,
     .   0.283249D+00,   0.236711D+00,   0.190893D+00,   0.148375D+00,
     .   0.107641D+00,   0.674164D-01,   0.282429D-01,  -0.118929D-01,
     .  -0.451581D-01,  -0.800718D-01,  -0.106407D+00,  -0.146132D+00,
     .  -0.177594D+00,  -0.209087D+00,  -0.237939D+00,  -0.268054D+00,
     .  -0.297210D+00,  -0.326472D+00,  -0.349928D+00,  -0.374204D+00,
     .  -0.401362D+00,  -0.427556D+00,  -0.449879D+00,  -0.461232D+00,
     .  -0.497959D+00,  -0.465266D+00,  -0.541772D+00,  -0.645032D+00,
     .  -0.739474D+00,  -0.828596D+00,  -0.905896D+00,  -0.978443D+00,
     .  -0.106542D+01,  -0.111143D+01,  -0.117268D+01,  -0.123118D+01/

      DATA (YHI(I),I=281,320)/
     .  -0.128217D+01,  -0.132113D+01,  -0.137951D+01,  -0.142661D+01,
     .  -0.146928D+01,  -0.151264D+01,  -0.155622D+01,  -0.159162D+01,
     .  -0.162531D+01,  -0.166463D+01,  -0.170018D+01,  -0.173385D+01,
     .  -0.176173D+01,  -0.179582D+01,  -0.182630D+01,  -0.185676D+01,
     .  -0.188124D+01,  -0.190601D+01,  -0.198109D+01,  -0.196626D+01,
     .  -0.199073D+01,  -0.204167D+01,  -0.208667D+01,  -0.213625D+01,
     .  -0.217197D+01,  -0.221764D+01,  -0.223929D+01,  -0.229089D+01,
     .  -0.232265D+01,  -0.236000D+01,  -0.239057D+01,  -0.242162D+01,
     .  -0.245804D+01,  -0.249028D+01,  -0.250018D+01,  -0.254592D+01,
     .  -0.257096D+01,  -0.257866D+01,  -0.263340D+01,  -0.265605D+01/

      DATA (YHI(I),I=321,360)/
     .  -0.267627D+01,  -0.270760D+01,  -0.272511D+01,  -0.274800D+01,
     .  -0.276102D+01,  -0.278071D+01,  -0.281972D+01,  -0.284543D+01,
     .  -0.283352D+01,  -0.287545D+01,  -0.289426D+01,  -0.298833D+01,
     .  -0.307003D+01,  -0.315041D+01,  -0.321780D+01,  -0.329049D+01,
     .  -0.335487D+01,  -0.336054D+01,  -0.345719D+01,  -0.350846D+01,
     .  -0.355815D+01,  -0.360314D+01,  -0.364808D+01,  -0.368901D+01,
     .  -0.372910D+01,  -0.376751D+01,  -0.380438D+01,  -0.383982D+01,
     .  -0.387395D+01,  -0.390684D+01,  -0.393860D+01,  -0.396929D+01,
     .  -0.399898D+01,  -0.402773D+01,  -0.405561D+01,  -0.408269D+01,
     .  -0.410905D+01,  -0.413447D+01,  -0.415930D+01,  -0.418350D+01/

      DATA (YHI(I),I=361,400)/
     .  -0.420826D+01,  -0.425293D+01,  -0.429561D+01,  -0.433733D+01,
     .  -0.437631D+01,  -0.441412D+01,  -0.445041D+01,  -0.448520D+01,
     .  -0.451893D+01,  -0.455133D+01,  -0.458262D+01,  -0.461287D+01,
     .  -0.464213D+01,  -0.467047D+01,  -0.469796D+01,  -0.472493D+01,
     .  -0.475084D+01,  -0.477576D+01,  -0.480023D+01,  -0.482410D+01,
     .  -0.484762D+01,  -0.487028D+01,  -0.489210D+01,  -0.491368D+01,
     .  -0.493501D+01,  -0.495558D+01,  -0.497571D+01,  -0.499540D+01,
     .  -0.501444D+01,  -0.503331D+01,  -0.505178D+01,  -0.513888D+01,
     .  -0.521829D+01,  -0.529127D+01,  -0.535898D+01,  -0.542150D+01,
     .  -0.548015D+01,  -0.553562D+01,  -0.557842D+01,  -0.562670D+01/

      DATA (YHI(I),I=401,440)/
     .  -0.568872D+01,  -0.559614D+01,  -0.563704D+01,  -0.567617D+01,
     .  -0.571222D+01,  -0.574924D+01,  -0.578356D+01,  -0.581671D+01,
     .  -0.584856D+01,  -0.587928D+01,  -0.590900D+01,  -0.593768D+01,
     .  -0.596543D+01,  -0.599231D+01,  -0.601838D+01,  -0.604366D+01,
     .  -0.606824D+01,  -0.609212D+01,  -0.611541D+01,  -0.613797D+01,
     .  -0.616005D+01,  -0.620252D+01,  -0.624308D+01,  -0.628174D+01,
     .  -0.631872D+01,  -0.635416D+01,  -0.638823D+01,  -0.642105D+01,
     .  -0.645264D+01,  -0.648309D+01,  -0.651248D+01,  -0.654093D+01,
     .  -0.656842D+01,  -0.659469D+01,  -0.662087D+01,  -0.664591D+01,
     .  -0.667035D+01,  -0.669401D+01,  -0.671705D+01,  -0.673959D+01/

      DATA (YHI(I),I=441,461)/
     .  -0.676149D+01,  -0.678273D+01,  -0.680348D+01,  -0.682388D+01,
     .  -0.684370D+01,  -0.686313D+01,  -0.688203D+01,  -0.690059D+01,
     .  -0.691884D+01,  -0.693653D+01,  -0.695399D+01,  -0.703614D+01,
     .  -0.711101D+01,  -0.717985D+01,  -0.724357D+01,  -0.730302D+01,
     .  -0.735815D+01,  -0.741018D+01,  -0.745916D+01,  -0.750572D+01,
     .  -0.745858D+01/

C--PSEUDOSCALAR HIGGS

      DATA (YAR(I),I=  1, 40)/
     .   0.158696D+02,   0.159218D+02,   0.159757D+02,   0.160313D+02,
     .   0.160889D+02,   0.161484D+02,   0.162101D+02,   0.163404D+02,
     .   0.164093D+02,   0.164810D+02,   0.165557D+02,   0.166336D+02,
     .   0.167149D+02,   0.168892D+02,   0.169828D+02,   0.170812D+02,
     .   0.171849D+02,   0.172945D+02,   0.174105D+02,   0.175337D+02,
     .   0.176649D+02,   0.178052D+02,   0.179556D+02,   0.181177D+02,
     .   0.182931D+02,   0.184841D+02,   0.186933D+02,   0.191811D+02,
     .   0.194701D+02,   0.197993D+02,   0.201803D+02,   0.206303D+02,
     .   0.211768D+02,   0.218661D+02,   0.227879D+02,   0.241489D+02,
     .   0.266372D+02,   0.270326D+02,   0.274797D+02,   0.279932D+02/

      DATA (YAR(I),I= 41, 80)/
     .   0.285944D+02,   0.293168D+02,   0.302174D+02,   0.314040D+02,
     .   0.331227D+02,   0.361754D+02,   0.471154D+02,   0.588099D+02,
     .   0.708533D+02,   0.707389D+02,   0.585302D+02,   0.464786D+02,
     .   0.348948D+02,   0.316033D+02,   0.297432D+02,   0.284583D+02,
     .   0.274838D+02,   0.267031D+02,   0.260544D+02,   0.255016D+02,
     .   0.250204D+02,   0.246005D+02,   0.219522D+02,   0.205297D+02,
     .   0.195777D+02,   0.188719D+02,   0.183157D+02,   0.178589D+02,
     .   0.174734D+02,   0.171409D+02,   0.168488D+02,   0.165893D+02,
     .   0.163566D+02,   0.161441D+02,   0.159504D+02,   0.157722D+02,
     .   0.156075D+02,   0.154550D+02,   0.153114D+02,   0.151770D+02/

      DATA (YAR(I),I= 81,120)/
     .   0.150511D+02,   0.149318D+02,   0.148239D+02,   0.147125D+02,
     .   0.146119D+02,   0.145157D+02,   0.144222D+02,   0.143341D+02,
     .   0.142498D+02,   0.141686D+02,   0.140915D+02,   0.140168D+02,
     .   0.139450D+02,   0.138759D+02,   0.138093D+02,   0.137450D+02,
     .   0.136821D+02,   0.136227D+02,   0.135645D+02,   0.135097D+02,
     .   0.134537D+02,   0.134006D+02,   0.133491D+02,   0.132990D+02,
     .   0.132512D+02,   0.132035D+02,   0.131574D+02,   0.131127D+02,
     .   0.130693D+02,   0.130276D+02,   0.129858D+02,   0.129452D+02,
     .   0.129061D+02,   0.128673D+02,   0.128298D+02,   0.127932D+02,
     .   0.127575D+02,   0.127231D+02,   0.126885D+02,   0.126549D+02/

      DATA (YAR(I),I=121,160)/
     .   0.126228D+02,   0.124702D+02,   0.123310D+02,   0.122046D+02,
     .   0.120888D+02,   0.119821D+02,   0.118836D+02,   0.117940D+02,
     .   0.117078D+02,   0.116289D+02,   0.115548D+02,   0.114866D+02,
     .   0.114202D+02,   0.113588D+02,   0.113012D+02,   0.112460D+02,
     .   0.111939D+02,   0.111446D+02,   0.110977D+02,   0.110535D+02,
     .   0.110105D+02,   0.109311D+02,   0.108587D+02,   0.107917D+02,
     .   0.107301D+02,   0.106732D+02,   0.106200D+02,   0.105707D+02,
     .   0.105248D+02,   0.104822D+02,   0.104412D+02,   0.104033D+02,
     .   0.103679D+02,   0.103349D+02,   0.103018D+02,   0.102717D+02,
     .   0.102430D+02,   0.102162D+02,   0.101901D+02,   0.101656D+02/

      DATA (YAR(I),I=161,200)/
     .   0.101421D+02,   0.101198D+02,   0.100980D+02,   0.100782D+02,
     .   0.100586D+02,   0.100400D+02,   0.100220D+02,   0.100050D+02,
     .   0.998850D+01,   0.997274D+01,   0.995747D+01,   0.994323D+01,
     .   0.992769D+01,   0.991542D+01,   0.990228D+01,   0.988978D+01,
     .   0.987587D+01,   0.986571D+01,   0.985459D+01,   0.984378D+01,
     .   0.983305D+01,   0.982301D+01,   0.981293D+01,   0.980350D+01,
     .   0.979431D+01,   0.979050D+01,   0.977668D+01,   0.976826D+01,
     .   0.976025D+01,   0.975183D+01,   0.974468D+01,   0.973745D+01,
     .   0.973054D+01,   0.972299D+01,   0.971611D+01,   0.970967D+01,
     .   0.970319D+01,   0.969712D+01,   0.969095D+01,   0.968504D+01/

      DATA (YAR(I),I=201,240)/
     .   0.967920D+01,   0.967386D+01,   0.966839D+01,   0.966325D+01,
     .   0.965788D+01,   0.965311D+01,   0.964786D+01,   0.964339D+01,
     .   0.963882D+01,   0.963439D+01,   0.962914D+01,   0.962534D+01,
     .   0.962147D+01,   0.961741D+01,   0.961292D+01,   0.960977D+01,
     .   0.960566D+01,   0.960235D+01,   0.959864D+01,   0.959531D+01,
     .   0.959126D+01,   0.957614D+01,   0.956225D+01,   0.954937D+01,
     .   0.953906D+01,   0.952858D+01,   0.952050D+01,   0.951315D+01,
     .   0.950643D+01,   0.950040D+01,   0.949334D+01,   0.949056D+01,
     .   0.948552D+01,   0.948295D+01,   0.947822D+01,   0.947754D+01,
     .   0.947505D+01,   0.947339D+01,   0.947215D+01,   0.947070D+01/

      DATA (YAR(I),I=241,280)/
     .   0.947004D+01,   0.946916D+01,   0.946931D+01,   0.946985D+01,
     .   0.947087D+01,   0.947338D+01,   0.947562D+01,   0.947907D+01,
     .   0.948136D+01,   0.948608D+01,   0.949024D+01,   0.949289D+01,
     .   0.949865D+01,   0.950358D+01,   0.949215D+01,   0.951384D+01,
     .   0.951877D+01,   0.952380D+01,   0.952957D+01,   0.953303D+01,
     .   0.953905D+01,   0.953903D+01,   0.955148D+01,   0.955415D+01,
     .   0.956308D+01,   0.956538D+01,   0.957536D+01,   0.958962D+01,
     .   0.958620D+01,   0.961314D+01,   0.959921D+01,   0.962927D+01,
     .   0.965978D+01,   0.968376D+01,   0.972000D+01,   0.974814D+01,
     .   0.975908D+01,   0.980308D+01,   0.982954D+01,   0.985456D+01/

      DATA (YAR(I),I=281,320)/
     .   0.988358D+01,   0.992011D+01,   0.993698D+01,   0.996114D+01,
     .   0.998536D+01,   0.100091D+02,   0.100306D+02,   0.100523D+02,
     .   0.100740D+02,   0.100997D+02,   0.101179D+02,   0.101456D+02,
     .   0.101485D+02,   0.101819D+02,   0.102045D+02,   0.102239D+02,
     .   0.102380D+02,   0.102544D+02,   0.102361D+02,   0.103043D+02,
     .   0.103200D+02,   0.103524D+02,   0.103928D+02,   0.104086D+02,
     .   0.104615D+02,   0.104868D+02,   0.105102D+02,   0.105516D+02,
     .   0.105851D+02,   0.106152D+02,   0.106431D+02,   0.106715D+02,
     .   0.106940D+02,   0.107216D+02,   0.107645D+02,   0.107773D+02,
     .   0.108038D+02,   0.108165D+02,   0.108488D+02,   0.108733D+02/

      DATA (YAR(I),I=321,360)/
     .   0.108925D+02,   0.109196D+02,   0.109480D+02,   0.109705D+02,
     .   0.109964D+02,   0.110216D+02,   0.110301D+02,   0.110409D+02,
     .   0.110612D+02,   0.110834D+02,   0.111183D+02,   0.112078D+02,
     .   0.113058D+02,   0.114004D+02,   0.114719D+02,   0.115449D+02,
     .   0.116169D+02,   0.117359D+02,   0.117572D+02,   0.118210D+02,
     .   0.118809D+02,   0.119413D+02,   0.119969D+02,   0.120531D+02,
     .   0.121062D+02,   0.121576D+02,   0.122075D+02,   0.122559D+02,
     .   0.123029D+02,   0.123486D+02,   0.123931D+02,   0.124365D+02,
     .   0.124788D+02,   0.125200D+02,   0.125603D+02,   0.125997D+02,
     .   0.126382D+02,   0.126758D+02,   0.127126D+02,   0.127487D+02/

      DATA (YAR(I),I=361,400)/
     .   0.127835D+02,   0.128519D+02,   0.129187D+02,   0.129818D+02,
     .   0.130440D+02,   0.131034D+02,   0.131611D+02,   0.132168D+02,
     .   0.132711D+02,   0.133237D+02,   0.133749D+02,   0.134247D+02,
     .   0.134733D+02,   0.135206D+02,   0.135667D+02,   0.136120D+02,
     .   0.136560D+02,   0.136988D+02,   0.137409D+02,   0.137821D+02,
     .   0.138220D+02,   0.138615D+02,   0.139006D+02,   0.139386D+02,
     .   0.139755D+02,   0.140120D+02,   0.140479D+02,   0.140832D+02,
     .   0.141181D+02,   0.141522D+02,   0.141856D+02,   0.143451D+02,
     .   0.144928D+02,   0.146306D+02,   0.147594D+02,   0.148815D+02,
     .   0.149965D+02,   0.151050D+02,   0.152034D+02,   0.153027D+02/

      DATA (YAR(I),I=401,440)/
     .   0.153933D+02,   0.155761D+02,   0.156636D+02,   0.157479D+02,
     .   0.158311D+02,   0.159074D+02,   0.159830D+02,   0.160563D+02,
     .   0.161272D+02,   0.161960D+02,   0.162627D+02,   0.163276D+02,
     .   0.163908D+02,   0.164522D+02,   0.165121D+02,   0.165705D+02,
     .   0.166275D+02,   0.166831D+02,   0.167374D+02,   0.167906D+02,
     .   0.168425D+02,   0.169431D+02,   0.170396D+02,   0.171323D+02,
     .   0.172218D+02,   0.173080D+02,   0.173912D+02,   0.174716D+02,
     .   0.175495D+02,   0.176249D+02,   0.176982D+02,   0.177694D+02,
     .   0.178386D+02,   0.179066D+02,   0.179716D+02,   0.180356D+02,
     .   0.180978D+02,   0.181586D+02,   0.182180D+02,   0.182760D+02/

      DATA (YAR(I),I=441,461)/
     .   0.183327D+02,   0.183883D+02,   0.184427D+02,   0.184958D+02,
     .   0.185480D+02,   0.185990D+02,   0.186490D+02,   0.186984D+02,
     .   0.187464D+02,   0.187940D+02,   0.188404D+02,   0.190612D+02,
     .   0.192650D+02,   0.194542D+02,   0.196309D+02,   0.197963D+02,
     .   0.199529D+02,   0.201007D+02,   0.202408D+02,   0.203737D+02,
     .   0.204582D+02/

      DATA (YAI(I),I=  1, 40)/
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00/

      DATA (YAI(I),I= 41, 80)/
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.165937D+02,   0.163500D+02,   0.157866D+02,
     .   0.144564D+02,   0.138551D+02,   0.134457D+02,   0.131284D+02,
     .   0.128686D+02,   0.126478D+02,   0.124541D+02,   0.122825D+02,
     .   0.121239D+02,   0.120015D+02,   0.110083D+02,   0.104071D+02,
     .   0.997378D+01,   0.963668D+01,   0.936151D+01,   0.912931D+01,
     .   0.892868D+01,   0.875285D+01,   0.859548D+01,   0.845409D+01,
     .   0.832571D+01,   0.820668D+01,   0.809683D+01,   0.799502D+01,
     .   0.789966D+01,   0.781088D+01,   0.772581D+01,   0.764562D+01/

      DATA (YAI(I),I= 81,120)/
     .   0.756969D+01,   0.749755D+01,   0.742783D+01,   0.736254D+01,
     .   0.729977D+01,   0.723938D+01,   0.717960D+01,   0.712279D+01,
     .   0.706858D+01,   0.701550D+01,   0.696474D+01,   0.691534D+01,
     .   0.686729D+01,   0.682052D+01,   0.677522D+01,   0.673111D+01,
     .   0.668834D+01,   0.664607D+01,   0.660531D+01,   0.656609D+01,
     .   0.652629D+01,   0.648770D+01,   0.645091D+01,   0.641455D+01,
     .   0.637888D+01,   0.634364D+01,   0.630966D+01,   0.627602D+01,
     .   0.624308D+01,   0.621104D+01,   0.617807D+01,   0.614776D+01,
     .   0.611735D+01,   0.608704D+01,   0.605769D+01,   0.602882D+01,
     .   0.599987D+01,   0.597188D+01,   0.594413D+01,   0.591675D+01/

      DATA (YAI(I),I=121,160)/
     .   0.588861D+01,   0.576174D+01,   0.564238D+01,   0.553080D+01,
     .   0.542532D+01,   0.532697D+01,   0.523331D+01,   0.514452D+01,
     .   0.506028D+01,   0.497973D+01,   0.490296D+01,   0.482902D+01,
     .   0.475882D+01,   0.469147D+01,   0.462614D+01,   0.456397D+01,
     .   0.450352D+01,   0.444551D+01,   0.438929D+01,   0.433384D+01,
     .   0.428252D+01,   0.418219D+01,   0.408778D+01,   0.399882D+01,
     .   0.391452D+01,   0.383443D+01,   0.375801D+01,   0.368564D+01,
     .   0.361564D+01,   0.354975D+01,   0.348659D+01,   0.342580D+01,
     .   0.336699D+01,   0.332507D+01,   0.325638D+01,   0.320415D+01,
     .   0.315361D+01,   0.310472D+01,   0.305750D+01,   0.301171D+01/

      DATA (YAI(I),I=161,200)/
     .   0.296736D+01,   0.292438D+01,   0.288128D+01,   0.284202D+01,
     .   0.280253D+01,   0.276428D+01,   0.272697D+01,   0.269060D+01,
     .   0.265510D+01,   0.262059D+01,   0.258685D+01,   0.255379D+01,
     .   0.252201D+01,   0.249057D+01,   0.245974D+01,   0.242976D+01,
     .   0.239767D+01,   0.237078D+01,   0.234308D+01,   0.231598D+01,
     .   0.228903D+01,   0.226254D+01,   0.223611D+01,   0.221079D+01,
     .   0.218528D+01,   0.217457D+01,   0.213741D+01,   0.211406D+01,
     .   0.209085D+01,   0.206526D+01,   0.204560D+01,   0.202363D+01,
     .   0.200323D+01,   0.198041D+01,   0.195919D+01,   0.193906D+01,
     .   0.191866D+01,   0.189872D+01,   0.187909D+01,   0.185960D+01/

      DATA (YAI(I),I=201,240)/
     .   0.184034D+01,   0.182179D+01,   0.180327D+01,   0.178480D+01,
     .   0.176692D+01,   0.174894D+01,   0.173140D+01,   0.171433D+01,
     .   0.169740D+01,   0.168058D+01,   0.166198D+01,   0.164770D+01,
     .   0.163128D+01,   0.161560D+01,   0.159913D+01,   0.158416D+01,
     .   0.156907D+01,   0.155369D+01,   0.153877D+01,   0.152384D+01,
     .   0.150946D+01,   0.143855D+01,   0.137167D+01,   0.130702D+01,
     .   0.124860D+01,   0.118908D+01,   0.113536D+01,   0.108362D+01,
     .   0.103362D+01,   0.985383D+00,   0.937598D+00,   0.894882D+00,
     .   0.850328D+00,   0.810513D+00,   0.766864D+00,   0.732489D+00,
     .   0.694695D+00,   0.658687D+00,   0.624212D+00,   0.590172D+00/

      DATA (YAI(I),I=241,280)/
     .   0.557265D+00,   0.493928D+00,   0.434282D+00,   0.377078D+00,
     .   0.322540D+00,   0.271614D+00,   0.221836D+00,   0.175480D+00,
     .   0.131278D+00,   0.876669D-01,   0.453968D-01,   0.260449D-02,
     .  -0.336707D-01,  -0.712215D-01,  -0.990712D-01,  -0.142171D+00,
     .  -0.175867D+00,  -0.209440D+00,  -0.240408D+00,  -0.272259D+00,
     .  -0.303267D+00,  -0.333794D+00,  -0.359626D+00,  -0.385394D+00,
     .  -0.414272D+00,  -0.441650D+00,  -0.465815D+00,  -0.479750D+00,
     .  -0.516504D+00,  -0.489225D+00,  -0.562987D+00,  -0.671934D+00,
     .  -0.771217D+00,  -0.863988D+00,  -0.945448D+00,  -0.102119D+01,
     .  -0.110912D+01,  -0.115945D+01,  -0.122290D+01,  -0.128329D+01/

      DATA (YAI(I),I=281,320)/
     .  -0.133643D+01,  -0.137823D+01,  -0.143730D+01,  -0.148579D+01,
     .  -0.152986D+01,  -0.157442D+01,  -0.161894D+01,  -0.165551D+01,
     .  -0.169033D+01,  -0.173063D+01,  -0.176688D+01,  -0.180168D+01,
     .  -0.182969D+01,  -0.186502D+01,  -0.189631D+01,  -0.192737D+01,
     .  -0.195241D+01,  -0.197778D+01,  -0.204980D+01,  -0.203940D+01,
     .  -0.206432D+01,  -0.211600D+01,  -0.216219D+01,  -0.221162D+01,
     .  -0.224919D+01,  -0.229508D+01,  -0.231773D+01,  -0.236981D+01,
     .  -0.240243D+01,  -0.244023D+01,  -0.247137D+01,  -0.250295D+01,
     .  -0.253941D+01,  -0.257200D+01,  -0.258363D+01,  -0.262860D+01,
     .  -0.265411D+01,  -0.266241D+01,  -0.271666D+01,  -0.273973D+01/

      DATA (YAI(I),I=321,360)/
     .  -0.276024D+01,  -0.279169D+01,  -0.280986D+01,  -0.283299D+01,
     .  -0.284671D+01,  -0.286681D+01,  -0.290494D+01,  -0.293034D+01,
     .  -0.291978D+01,  -0.296114D+01,  -0.298064D+01,  -0.307517D+01,
     .  -0.315780D+01,  -0.323886D+01,  -0.330652D+01,  -0.337925D+01,
     .  -0.344382D+01,  -0.345290D+01,  -0.354709D+01,  -0.359855D+01,
     .  -0.364834D+01,  -0.369356D+01,  -0.373855D+01,  -0.377965D+01,
     .  -0.381984D+01,  -0.385833D+01,  -0.389527D+01,  -0.393078D+01,
     .  -0.396496D+01,  -0.399789D+01,  -0.402968D+01,  -0.406040D+01,
     .  -0.409012D+01,  -0.411890D+01,  -0.414679D+01,  -0.417388D+01,
     .  -0.420024D+01,  -0.422567D+01,  -0.425051D+01,  -0.427470D+01/

      DATA (YAI(I),I=361,400)/
     .  -0.429941D+01,  -0.434408D+01,  -0.438677D+01,  -0.442843D+01,
     .  -0.446740D+01,  -0.450518D+01,  -0.454143D+01,  -0.457618D+01,
     .  -0.460985D+01,  -0.464221D+01,  -0.467345D+01,  -0.470364D+01,
     .  -0.473285D+01,  -0.476113D+01,  -0.478857D+01,  -0.481548D+01,
     .  -0.484133D+01,  -0.486619D+01,  -0.489061D+01,  -0.491441D+01,
     .  -0.493787D+01,  -0.496046D+01,  -0.498224D+01,  -0.500376D+01,
     .  -0.502502D+01,  -0.504554D+01,  -0.506562D+01,  -0.508524D+01,
     .  -0.510424D+01,  -0.512305D+01,  -0.514147D+01,  -0.522829D+01,
     .  -0.530743D+01,  -0.538015D+01,  -0.544760D+01,  -0.550989D+01,
     .  -0.556831D+01,  -0.562354D+01,  -0.566626D+01,  -0.571436D+01/

      DATA (YAI(I),I=401,440)/
     .  -0.577575D+01,  -0.568757D+01,  -0.572828D+01,  -0.576723D+01,
     .  -0.580316D+01,  -0.583996D+01,  -0.587412D+01,  -0.590711D+01,
     .  -0.593880D+01,  -0.596937D+01,  -0.599895D+01,  -0.602748D+01,
     .  -0.605510D+01,  -0.608184D+01,  -0.610778D+01,  -0.613293D+01,
     .  -0.615739D+01,  -0.618116D+01,  -0.620432D+01,  -0.622676D+01,
     .  -0.624873D+01,  -0.629098D+01,  -0.633132D+01,  -0.636978D+01,
     .  -0.640656D+01,  -0.644182D+01,  -0.647570D+01,  -0.650834D+01,
     .  -0.653976D+01,  -0.657004D+01,  -0.659927D+01,  -0.662757D+01,
     .  -0.665490D+01,  -0.668104D+01,  -0.670706D+01,  -0.673197D+01,
     .  -0.675626D+01,  -0.677979D+01,  -0.680271D+01,  -0.682511D+01/

      DATA (YAI(I),I=441,461)/
     .  -0.684689D+01,  -0.686801D+01,  -0.688865D+01,  -0.690893D+01,
     .  -0.692863D+01,  -0.694795D+01,  -0.696675D+01,  -0.698520D+01,
     .  -0.700334D+01,  -0.702094D+01,  -0.703829D+01,  -0.711997D+01,
     .  -0.719440D+01,  -0.726284D+01,  -0.732619D+01,  -0.738528D+01,
     .  -0.744009D+01,  -0.749182D+01,  -0.754051D+01,  -0.758679D+01,
     .  -0.754063D+01/

      END

      BLOCKDATA CVIRTS
C--DATA OF VIRTUAL CORRECTIONS TO SQUARK LOOPS
      PARAMETER(NN=999)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/RESINTS/XX(NN),YHR(NN),YHI(NN),YAR(NN),YAI(NN)

      DATA ( XX(I),I=  1, 40)/
     .   0.000000D+00,   0.100000D+00,   0.200000D+00,   0.300000D+00,
     .   0.400000D+00,   0.500000D+00,   0.600000D+00,   0.800000D+00,
     .   0.900000D+00,   0.100000D+01,   0.110000D+01,   0.120000D+01,
     .   0.130000D+01,   0.150000D+01,   0.160000D+01,   0.170000D+01,
     .   0.180000D+01,   0.190000D+01,   0.200000D+01,   0.210000D+01,
     .   0.220000D+01,   0.230000D+01,   0.240000D+01,   0.250000D+01,
     .   0.260000D+01,   0.270000D+01,   0.280000D+01,   0.300000D+01,
     .   0.310000D+01,   0.320000D+01,   0.330000D+01,   0.340000D+01,
     .   0.350000D+01,   0.360000D+01,   0.370000D+01,   0.380000D+01,
     .   0.390000D+01,   0.391000D+01,   0.392000D+01,   0.393000D+01/

      DATA ( XX(I),I= 41, 80)/
     .   0.394000D+01,   0.395000D+01,   0.396000D+01,   0.397000D+01,
     .   0.398000D+01,   0.399000D+01,   0.399900D+01,   0.399990D+01,
     .   0.399999D+01,   0.400001D+01,   0.400010D+01,   0.400100D+01,
     .   0.401000D+01,   0.402000D+01,   0.403000D+01,   0.404000D+01,
     .   0.405000D+01,   0.406000D+01,   0.407000D+01,   0.408000D+01,
     .   0.409000D+01,   0.410000D+01,   0.420000D+01,   0.430000D+01,
     .   0.440000D+01,   0.450000D+01,   0.460000D+01,   0.470000D+01,
     .   0.480000D+01,   0.490000D+01,   0.500000D+01,   0.510000D+01,
     .   0.520000D+01,   0.530000D+01,   0.540000D+01,   0.550000D+01,
     .   0.560000D+01,   0.570000D+01,   0.580000D+01,   0.590000D+01/

      DATA ( XX(I),I= 81,120)/
     .   0.600000D+01,   0.610000D+01,   0.620000D+01,   0.630000D+01,
     .   0.640000D+01,   0.650000D+01,   0.660000D+01,   0.670000D+01,
     .   0.680000D+01,   0.690000D+01,   0.700000D+01,   0.710000D+01,
     .   0.720000D+01,   0.730000D+01,   0.740000D+01,   0.750000D+01,
     .   0.760000D+01,   0.770000D+01,   0.780000D+01,   0.790000D+01,
     .   0.800000D+01,   0.810000D+01,   0.820000D+01,   0.830000D+01,
     .   0.840000D+01,   0.850000D+01,   0.860000D+01,   0.870000D+01,
     .   0.880000D+01,   0.890000D+01,   0.900000D+01,   0.910000D+01,
     .   0.920000D+01,   0.930000D+01,   0.940000D+01,   0.950000D+01,
     .   0.960000D+01,   0.970000D+01,   0.980000D+01,   0.990000D+01/

      DATA ( XX(I),I=121,160)/
     .   0.100000D+02,   0.105000D+02,   0.110000D+02,   0.115000D+02,
     .   0.120000D+02,   0.125000D+02,   0.130000D+02,   0.135000D+02,
     .   0.140000D+02,   0.145000D+02,   0.150000D+02,   0.155000D+02,
     .   0.160000D+02,   0.165000D+02,   0.170000D+02,   0.175000D+02,
     .   0.180000D+02,   0.185000D+02,   0.190000D+02,   0.195000D+02,
     .   0.200000D+02,   0.210000D+02,   0.220000D+02,   0.230000D+02,
     .   0.240000D+02,   0.250000D+02,   0.260000D+02,   0.270000D+02,
     .   0.280000D+02,   0.290000D+02,   0.300000D+02,   0.310000D+02,
     .   0.320000D+02,   0.330000D+02,   0.340000D+02,   0.350000D+02,
     .   0.360000D+02,   0.370000D+02,   0.380000D+02,   0.390000D+02/

      DATA ( XX(I),I=161,200)/
     .   0.400000D+02,   0.410000D+02,   0.420000D+02,   0.430000D+02,
     .   0.440000D+02,   0.450000D+02,   0.460000D+02,   0.470000D+02,
     .   0.480000D+02,   0.490000D+02,   0.500000D+02,   0.510000D+02,
     .   0.520000D+02,   0.530000D+02,   0.540000D+02,   0.550000D+02,
     .   0.560000D+02,   0.570000D+02,   0.580000D+02,   0.590000D+02,
     .   0.600000D+02,   0.610000D+02,   0.620000D+02,   0.630000D+02,
     .   0.640000D+02,   0.650000D+02,   0.660000D+02,   0.670000D+02,
     .   0.680000D+02,   0.690000D+02,   0.700000D+02,   0.710000D+02,
     .   0.720000D+02,   0.730000D+02,   0.740000D+02,   0.750000D+02,
     .   0.760000D+02,   0.770000D+02,   0.780000D+02,   0.790000D+02/

      DATA ( XX(I),I=201,240)/
     .   0.800000D+02,   0.810000D+02,   0.820000D+02,   0.830000D+02,
     .   0.840000D+02,   0.850000D+02,   0.860000D+02,   0.870000D+02,
     .   0.880000D+02,   0.890000D+02,   0.900000D+02,   0.910000D+02,
     .   0.920000D+02,   0.930000D+02,   0.940000D+02,   0.950000D+02,
     .   0.960000D+02,   0.970000D+02,   0.980000D+02,   0.990000D+02,
     .   0.100000D+03,   0.105000D+03,   0.110000D+03,   0.115000D+03,
     .   0.120000D+03,   0.125000D+03,   0.130000D+03,   0.135000D+03,
     .   0.140000D+03,   0.145000D+03,   0.150000D+03,   0.155000D+03,
     .   0.160000D+03,   0.165000D+03,   0.170000D+03,   0.175000D+03,
     .   0.180000D+03,   0.185000D+03,   0.190000D+03,   0.195000D+03/

      DATA ( XX(I),I=241,280)/
     .   0.200000D+03,   0.210000D+03,   0.220000D+03,   0.230000D+03,
     .   0.240000D+03,   0.250000D+03,   0.260000D+03,   0.270000D+03,
     .   0.280000D+03,   0.290000D+03,   0.300000D+03,   0.310000D+03,
     .   0.320000D+03,   0.330000D+03,   0.340000D+03,   0.350000D+03,
     .   0.360000D+03,   0.370000D+03,   0.380000D+03,   0.390000D+03,
     .   0.400000D+03,   0.410000D+03,   0.420000D+03,   0.430000D+03,
     .   0.440000D+03,   0.450000D+03,   0.460000D+03,   0.470000D+03,
     .   0.480000D+03,   0.490000D+03,   0.500000D+03,   0.550000D+03,
     .   0.600000D+03,   0.650000D+03,   0.700000D+03,   0.750000D+03,
     .   0.800000D+03,   0.850000D+03,   0.900000D+03,   0.950000D+03/

      DATA ( XX(I),I=281,320)/
     .   0.100000D+04,   0.105000D+04,   0.110000D+04,   0.115000D+04,
     .   0.120000D+04,   0.125000D+04,   0.130000D+04,   0.135000D+04,
     .   0.140000D+04,   0.145000D+04,   0.150000D+04,   0.155000D+04,
     .   0.160000D+04,   0.165000D+04,   0.170000D+04,   0.175000D+04,
     .   0.180000D+04,   0.185000D+04,   0.190000D+04,   0.195000D+04,
     .   0.200000D+04,   0.210000D+04,   0.220000D+04,   0.230000D+04,
     .   0.240000D+04,   0.250000D+04,   0.260000D+04,   0.270000D+04,
     .   0.280000D+04,   0.290000D+04,   0.300000D+04,   0.310000D+04,
     .   0.320000D+04,   0.330000D+04,   0.340000D+04,   0.350000D+04,
     .   0.360000D+04,   0.370000D+04,   0.380000D+04,   0.390000D+04/

      DATA ( XX(I),I=321,360)/
     .   0.400000D+04,   0.410000D+04,   0.420000D+04,   0.430000D+04,
     .   0.440000D+04,   0.450000D+04,   0.460000D+04,   0.470000D+04,
     .   0.480000D+04,   0.490000D+04,   0.500000D+04,   0.550000D+04,
     .   0.600000D+04,   0.650000D+04,   0.700000D+04,   0.750000D+04,
     .   0.800000D+04,   0.850000D+04,   0.900000D+04,   0.950000D+04,
     .   0.100000D+05,   0.105000D+05,   0.110000D+05,   0.115000D+05,
     .   0.120000D+05,   0.125000D+05,   0.130000D+05,   0.135000D+05,
     .   0.140000D+05,   0.145000D+05,   0.150000D+05,   0.155000D+05,
     .   0.160000D+05,   0.165000D+05,   0.170000D+05,   0.175000D+05,
     .   0.180000D+05,   0.185000D+05,   0.190000D+05,   0.195000D+05/

      DATA ( XX(I),I=361,400)/
     .   0.200000D+05,   0.210000D+05,   0.220000D+05,   0.230000D+05,
     .   0.240000D+05,   0.250000D+05,   0.260000D+05,   0.270000D+05,
     .   0.280000D+05,   0.290000D+05,   0.300000D+05,   0.310000D+05,
     .   0.320000D+05,   0.330000D+05,   0.340000D+05,   0.350000D+05,
     .   0.360000D+05,   0.370000D+05,   0.380000D+05,   0.390000D+05,
     .   0.400000D+05,   0.410000D+05,   0.420000D+05,   0.430000D+05,
     .   0.440000D+05,   0.450000D+05,   0.460000D+05,   0.470000D+05,
     .   0.480000D+05,   0.490000D+05,   0.500000D+05,   0.550000D+05,
     .   0.600000D+05,   0.650000D+05,   0.700000D+05,   0.750000D+05,
     .   0.800000D+05,   0.850000D+05,   0.900000D+05,   0.950000D+05/

      DATA ( XX(I),I=401,440)/
     .   0.100000D+06,   0.105000D+06,   0.110000D+06,   0.115000D+06,
     .   0.120000D+06,   0.125000D+06,   0.130000D+06,   0.135000D+06,
     .   0.140000D+06,   0.145000D+06,   0.150000D+06,   0.155000D+06,
     .   0.160000D+06,   0.165000D+06,   0.170000D+06,   0.175000D+06,
     .   0.180000D+06,   0.185000D+06,   0.190000D+06,   0.195000D+06,
     .   0.200000D+06,   0.210000D+06,   0.220000D+06,   0.230000D+06,
     .   0.240000D+06,   0.250000D+06,   0.260000D+06,   0.270000D+06,
     .   0.280000D+06,   0.290000D+06,   0.300000D+06,   0.310000D+06,
     .   0.320000D+06,   0.330000D+06,   0.340000D+06,   0.350000D+06,
     .   0.360000D+06,   0.370000D+06,   0.380000D+06,   0.390000D+06/

      DATA ( XX(I),I=441,461)/
     .   0.400000D+06,   0.410000D+06,   0.420000D+06,   0.430000D+06,
     .   0.440000D+06,   0.450000D+06,   0.460000D+06,   0.470000D+06,
     .   0.480000D+06,   0.490000D+06,   0.500000D+06,   0.550000D+06,
     .   0.600000D+06,   0.650000D+06,   0.700000D+06,   0.750000D+06,
     .   0.800000D+06,   0.850000D+06,   0.900000D+06,   0.950000D+06,
     .   0.100000D+07/

C--SCALAR HIGGS

      DATA (YHR(I),I=  1, 40)/
     .   0.188696D+02,   0.189337D+02,   0.189996D+02,   0.190678D+02,
     .   0.191383D+02,   0.192113D+02,   0.192869D+02,   0.194468D+02,
     .   0.195315D+02,   0.196196D+02,   0.197115D+02,   0.198073D+02,
     .   0.199075D+02,   0.201224D+02,   0.202381D+02,   0.203598D+02,
     .   0.204884D+02,   0.206243D+02,   0.207686D+02,   0.209220D+02,
     .   0.210859D+02,   0.212613D+02,   0.214501D+02,   0.216540D+02,
     .   0.218754D+02,   0.221173D+02,   0.223832D+02,   0.230076D+02,
     .   0.233804D+02,   0.238077D+02,   0.243058D+02,   0.248989D+02,
     .   0.256260D+02,   0.265538D+02,   0.278123D+02,   0.297055D+02,
     .   0.332628D+02,   0.338381D+02,   0.344920D+02,   0.352467D+02/

      DATA (YHR(I),I= 41, 80)/
     .   0.361355D+02,   0.372104D+02,   0.385602D+02,   0.403542D+02,
     .   0.429813D+02,   0.477185D+02,   0.652183D+02,   0.844421D+02,
     .   0.104498D+03,   0.104213D+03,   0.837654D+02,   0.637354D+02,
     .   0.448930D+02,   0.396874D+02,   0.367915D+02,   0.348147D+02,
     .   0.333303D+02,   0.321516D+02,   0.311796D+02,   0.303573D+02,
     .   0.296459D+02,   0.290302D+02,   0.252402D+02,   0.232962D+02,
     .   0.220441D+02,   0.211471D+02,   0.204624D+02,   0.199168D+02,
     .   0.194691D+02,   0.190933D+02,   0.187714D+02,   0.184926D+02,
     .   0.182485D+02,   0.180302D+02,   0.178356D+02,   0.176605D+02,
     .   0.175017D+02,   0.173577D+02,   0.172243D+02,   0.171017D+02/

      DATA (YHR(I),I= 81,120)/
     .   0.169888D+02,   0.168839D+02,   0.167896D+02,   0.166953D+02,
     .   0.166109D+02,   0.165314D+02,   0.164546D+02,   0.163835D+02,
     .   0.163166D+02,   0.162528D+02,   0.161931D+02,   0.161360D+02,
     .   0.160816D+02,   0.160298D+02,   0.159806D+02,   0.159335D+02,
     .   0.158881D+02,   0.158453D+02,   0.158040D+02,   0.157656D+02,
     .   0.157262D+02,   0.156894D+02,   0.156544D+02,   0.156205D+02,
     .   0.155883D+02,   0.155563D+02,   0.155259D+02,   0.154965D+02,
     .   0.154682D+02,   0.154412D+02,   0.154138D+02,   0.153883D+02,
     .   0.153636D+02,   0.153391D+02,   0.153157D+02,   0.152930D+02,
     .   0.152707D+02,   0.152496D+02,   0.152284D+02,   0.152079D+02/

      DATA (YHR(I),I=121,160)/
     .   0.151879D+02,   0.150973D+02,   0.150164D+02,   0.149445D+02,
     .   0.148794D+02,   0.148210D+02,   0.147674D+02,   0.147191D+02,
     .   0.146735D+02,   0.146317D+02,   0.145929D+02,   0.145569D+02,
     .   0.145227D+02,   0.144910D+02,   0.144609D+02,   0.144325D+02,
     .   0.144054D+02,   0.143799D+02,   0.143554D+02,   0.143318D+02,
     .   0.143098D+02,   0.142678D+02,   0.142289D+02,   0.141927D+02,
     .   0.141587D+02,   0.141267D+02,   0.140963D+02,   0.140675D+02,
     .   0.140399D+02,   0.140139D+02,   0.139886D+02,   0.139644D+02,
     .   0.139410D+02,   0.139226D+02,   0.138963D+02,   0.138751D+02,
     .   0.138545D+02,   0.138344D+02,   0.138147D+02,   0.137956D+02/

      DATA (YHR(I),I=161,200)/
     .   0.137768D+02,   0.137585D+02,   0.137402D+02,   0.137230D+02,
     .   0.137058D+02,   0.136889D+02,   0.136722D+02,   0.136559D+02,
     .   0.136398D+02,   0.136240D+02,   0.136084D+02,   0.135930D+02,
     .   0.135777D+02,   0.135630D+02,   0.135482D+02,   0.135337D+02,
     .   0.135185D+02,   0.135049D+02,   0.134910D+02,   0.134772D+02,
     .   0.134635D+02,   0.134500D+02,   0.134365D+02,   0.134232D+02,
     .   0.134100D+02,   0.134002D+02,   0.133843D+02,   0.133716D+02,
     .   0.133590D+02,   0.133459D+02,   0.133341D+02,   0.133218D+02,
     .   0.133101D+02,   0.132975D+02,   0.132855D+02,   0.132737D+02,
     .   0.132619D+02,   0.132503D+02,   0.132387D+02,   0.132272D+02/

      DATA (YHR(I),I=201,240)/
     .   0.132158D+02,   0.132045D+02,   0.131933D+02,   0.131821D+02,
     .   0.131710D+02,   0.131601D+02,   0.131491D+02,   0.131383D+02,
     .   0.131276D+02,   0.131169D+02,   0.131059D+02,   0.130957D+02,
     .   0.130853D+02,   0.130749D+02,   0.130644D+02,   0.130543D+02,
     .   0.130441D+02,   0.130340D+02,   0.130239D+02,   0.130139D+02,
     .   0.130039D+02,   0.129550D+02,   0.129075D+02,   0.128611D+02,
     .   0.128164D+02,   0.127722D+02,   0.127296D+02,   0.126880D+02,
     .   0.126471D+02,   0.126072D+02,   0.125678D+02,   0.125299D+02,
     .   0.124922D+02,   0.124557D+02,   0.124192D+02,   0.123844D+02,
     .   0.123497D+02,   0.123157D+02,   0.122823D+02,   0.122494D+02/

      DATA (YHR(I),I=241,280)/
     .   0.122172D+02,   0.121542D+02,   0.120932D+02,   0.120341D+02,
     .   0.119767D+02,   0.119210D+02,   0.118668D+02,   0.118142D+02,
     .   0.117629D+02,   0.117130D+02,   0.116642D+02,   0.116163D+02,
     .   0.115701D+02,   0.115247D+02,   0.114800D+02,   0.114370D+02,
     .   0.113945D+02,   0.113528D+02,   0.113121D+02,   0.112720D+02,
     .   0.112329D+02,   0.111941D+02,   0.111569D+02,   0.111199D+02,
     .   0.110836D+02,   0.110477D+02,   0.110129D+02,   0.109793D+02,
     .   0.109445D+02,   0.109142D+02,   0.108784D+02,   0.107221D+02,
     .   0.105769D+02,   0.104412D+02,   0.103144D+02,   0.101948D+02,
     .   0.100808D+02,   0.997476D+01,   0.987298D+01,   0.977597D+01/

      DATA (YHR(I),I=281,320)/
     .   0.968356D+01,   0.959556D+01,   0.951023D+01,   0.942868D+01,
     .   0.935030D+01,   0.927476D+01,   0.920181D+01,   0.913151D+01,
     .   0.906354D+01,   0.899770D+01,   0.893382D+01,   0.887206D+01,
     .   0.881174D+01,   0.875356D+01,   0.869684D+01,   0.864158D+01,
     .   0.858777D+01,   0.853536D+01,   0.848292D+01,   0.843445D+01,
     .   0.838574D+01,   0.829160D+01,   0.820176D+01,   0.811528D+01,
     .   0.803294D+01,   0.795328D+01,   0.787685D+01,   0.780303D+01,
     .   0.773189D+01,   0.766303D+01,   0.759644D+01,   0.753194D+01,
     .   0.746928D+01,   0.740858D+01,   0.734998D+01,   0.729238D+01,
     .   0.723666D+01,   0.718243D+01,   0.712933D+01,   0.707788D+01/

      DATA (YHR(I),I=321,360)/
     .   0.702760D+01,   0.697845D+01,   0.693068D+01,   0.688385D+01,
     .   0.683817D+01,   0.679343D+01,   0.674932D+01,   0.670635D+01,
     .   0.666455D+01,   0.662324D+01,   0.658300D+01,   0.639226D+01,
     .   0.621794D+01,   0.605732D+01,   0.590835D+01,   0.576950D+01,
     .   0.563954D+01,   0.551788D+01,   0.540222D+01,   0.529316D+01,
     .   0.518963D+01,   0.509115D+01,   0.499719D+01,   0.490741D+01,
     .   0.482142D+01,   0.473892D+01,   0.465964D+01,   0.458334D+01,
     .   0.450981D+01,   0.443885D+01,   0.437029D+01,   0.430398D+01,
     .   0.423976D+01,   0.417752D+01,   0.411713D+01,   0.405849D+01,
     .   0.400150D+01,   0.394607D+01,   0.389212D+01,   0.383957D+01/

      DATA (YHR(I),I=361,400)/
     .   0.378834D+01,   0.368963D+01,   0.359552D+01,   0.350558D+01,
     .   0.341949D+01,   0.333691D+01,   0.325757D+01,   0.318123D+01,
     .   0.310767D+01,   0.303669D+01,   0.296813D+01,   0.290182D+01,
     .   0.283762D+01,   0.277539D+01,   0.271504D+01,   0.265643D+01,
     .   0.259948D+01,   0.254409D+01,   0.249019D+01,   0.243768D+01,
     .   0.238651D+01,   0.233661D+01,   0.228791D+01,   0.224037D+01,
     .   0.219392D+01,   0.214851D+01,   0.210411D+01,   0.206066D+01,
     .   0.201813D+01,   0.197649D+01,   0.193568D+01,   0.174321D+01,
     .   0.156755D+01,   0.140600D+01,   0.125647D+01,   0.111729D+01,
     .   0.987127D+00,   0.864882D+00,   0.749649D+00,   0.640668D+00/

      DATA (YHR(I),I=401,440)/
     .   0.537286D+00,   0.439179D+00,   0.345450D+00,   0.255902D+00,
     .   0.170180D+00,   0.879652D-01,   0.898705D-02,  -0.670007D-01,
     .  -0.140216D+00,  -0.210854D+00,  -0.279090D+00,  -0.345081D+00,
     .  -0.408971D+00,  -0.470889D+00,  -0.530953D+00,  -0.589270D+00,
     .  -0.645940D+00,  -0.701052D+00,  -0.754691D+00,  -0.806932D+00,
     .  -0.857847D+00,  -0.955955D+00,  -0.104949D+01,  -0.113885D+01,
     .  -0.122439D+01,  -0.130644D+01,  -0.138526D+01,  -0.146110D+01,
     .  -0.153417D+01,  -0.160467D+01,  -0.167277D+01,  -0.173864D+01,
     .  -0.180241D+01,  -0.186421D+01,  -0.192417D+01,  -0.198238D+01,
     .  -0.203895D+01,  -0.209396D+01,  -0.214751D+01,  -0.219966D+01/

      DATA (YHR(I),I=441,461)/
     .  -0.225048D+01,  -0.230005D+01,  -0.234842D+01,  -0.239565D+01,
     .  -0.244180D+01,  -0.248690D+01,  -0.253101D+01,  -0.257417D+01,
     .  -0.261642D+01,  -0.265780D+01,  -0.269834D+01,  -0.288957D+01,
     .  -0.306412D+01,  -0.322467D+01,  -0.337329D+01,  -0.351164D+01,
     .  -0.364105D+01,  -0.376259D+01,  -0.387718D+01,  -0.398556D+01,
     .  -0.408837D+01/

      DATA (YHI(I),I=  1, 40)/
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00/

      DATA (YHI(I),I= 41, 80)/
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.000000D+00,   0.000000D+00,   0.000000D+00,
     .   0.000000D+00,   0.277727D+02,   0.271700D+02,   0.257998D+02,
     .   0.227660D+02,   0.214315D+02,   0.205384D+02,   0.198542D+02,
     .   0.192971D+02,   0.188256D+02,   0.184146D+02,   0.180512D+02,
     .   0.177204D+02,   0.174433D+02,   0.153853D+02,   0.141416D+02,
     .   0.132506D+02,   0.125608D+02,   0.120009D+02,   0.115318D+02,
     .   0.111296D+02,   0.107793D+02,   0.104693D+02,   0.101927D+02,
     .   0.994336D+01,   0.971634D+01,   0.950864D+01,   0.931788D+01,
     .   0.914137D+01,   0.897786D+01,   0.882487D+01,   0.868201D+01/

      DATA (YHI(I),I= 81,120)/
     .   0.854794D+01,   0.842218D+01,   0.830063D+01,   0.819100D+01,
     .   0.808463D+01,   0.798376D+01,   0.788746D+01,   0.779573D+01,
     .   0.770869D+01,   0.762524D+01,   0.754548D+01,   0.746919D+01,
     .   0.739591D+01,   0.732549D+01,   0.725786D+01,   0.719285D+01,
     .   0.713069D+01,   0.706976D+01,   0.701163D+01,   0.695512D+01,
     .   0.690097D+01,   0.684821D+01,   0.679774D+01,   0.674870D+01,
     .   0.670074D+01,   0.665467D+01,   0.661012D+01,   0.656665D+01,
     .   0.652444D+01,   0.648330D+01,   0.644302D+01,   0.640492D+01,
     .   0.636727D+01,   0.633065D+01,   0.629511D+01,   0.626045D+01,
     .   0.622639D+01,   0.619319D+01,   0.616119D+01,   0.612979D+01/

      DATA (YHI(I),I=121,160)/
     .   0.609829D+01,   0.595553D+01,   0.582850D+01,   0.571438D+01,
     .   0.561116D+01,   0.551808D+01,   0.543303D+01,   0.535455D+01,
     .   0.528377D+01,   0.521793D+01,   0.515727D+01,   0.510046D+01,
     .   0.504888D+01,   0.500048D+01,   0.495507D+01,   0.491320D+01,
     .   0.487377D+01,   0.483693D+01,   0.480234D+01,   0.476935D+01,
     .   0.473939D+01,   0.468358D+01,   0.463389D+01,   0.458979D+01,
     .   0.455019D+01,   0.451457D+01,   0.448259D+01,   0.445383D+01,
     .   0.442758D+01,   0.440399D+01,   0.438293D+01,   0.436365D+01,
     .   0.434598D+01,   0.433293D+01,   0.431581D+01,   0.430276D+01,
     .   0.429093D+01,   0.428003D+01,   0.427046D+01,   0.426165D+01/

      DATA (YHI(I),I=161,200)/
     .   0.425372D+01,   0.424656D+01,   0.424001D+01,   0.423434D+01,
     .   0.422920D+01,   0.422460D+01,   0.422056D+01,   0.421695D+01,
     .   0.421382D+01,   0.421111D+01,   0.420880D+01,   0.420674D+01,
     .   0.420549D+01,   0.420391D+01,   0.420289D+01,   0.420214D+01,
     .   0.420161D+01,   0.420132D+01,   0.420133D+01,   0.420156D+01,
     .   0.420198D+01,   0.420251D+01,   0.420322D+01,   0.420411D+01,
     .   0.420508D+01,   0.420747D+01,   0.420769D+01,   0.420917D+01,
     .   0.421071D+01,   0.421212D+01,   0.421416D+01,   0.421600D+01,
     .   0.421812D+01,   0.422003D+01,   0.422214D+01,   0.422438D+01,
     .   0.422665D+01,   0.422897D+01,   0.423140D+01,   0.423385D+01/

      DATA (YHI(I),I=201,240)/
     .   0.423636D+01,   0.423892D+01,   0.424154D+01,   0.424416D+01,
     .   0.424691D+01,   0.424959D+01,   0.425243D+01,   0.425522D+01,
     .   0.425808D+01,   0.426094D+01,   0.426375D+01,   0.426684D+01,
     .   0.426972D+01,   0.427273D+01,   0.427573D+01,   0.427870D+01,
     .   0.428182D+01,   0.428479D+01,   0.428789D+01,   0.429093D+01,
     .   0.429414D+01,   0.430963D+01,   0.432538D+01,   0.434116D+01,
     .   0.435707D+01,   0.437273D+01,   0.438840D+01,   0.440393D+01,
     .   0.441927D+01,   0.443441D+01,   0.444943D+01,   0.446409D+01,
     .   0.447856D+01,   0.449286D+01,   0.450679D+01,   0.452072D+01,
     .   0.453429D+01,   0.454762D+01,   0.456076D+01,   0.457366D+01/

      DATA (YHI(I),I=241,280)/
     .   0.458631D+01,   0.461095D+01,   0.463475D+01,   0.465773D+01,
     .   0.467993D+01,   0.470137D+01,   0.472209D+01,   0.474214D+01,
     .   0.476166D+01,   0.478035D+01,   0.479849D+01,   0.481607D+01,
     .   0.483320D+01,   0.484974D+01,   0.486725D+01,   0.488142D+01,
     .   0.489662D+01,   0.491134D+01,   0.492570D+01,   0.493972D+01,
     .   0.495322D+01,   0.496669D+01,   0.497932D+01,   0.499212D+01,
     .   0.500412D+01,   0.501618D+01,   0.502769D+01,   0.503911D+01,
     .   0.505016D+01,   0.506206D+01,   0.507159D+01,   0.512118D+01,
     .   0.516573D+01,   0.520622D+01,   0.524277D+01,   0.527649D+01,
     .   0.530753D+01,   0.533610D+01,   0.536263D+01,   0.538733D+01/

      DATA (YHI(I),I=281,320)/
     .   0.541038D+01,   0.543196D+01,   0.545222D+01,   0.547129D+01,
     .   0.548932D+01,   0.550631D+01,   0.552240D+01,   0.553779D+01,
     .   0.555240D+01,   0.556610D+01,   0.557934D+01,   0.559180D+01,
     .   0.560427D+01,   0.561549D+01,   0.562649D+01,   0.563708D+01,
     .   0.564741D+01,   0.565729D+01,   0.566655D+01,   0.567545D+01,
     .   0.568425D+01,   0.570085D+01,   0.571627D+01,   0.573095D+01,
     .   0.574436D+01,   0.575719D+01,   0.576958D+01,   0.578067D+01,
     .   0.579148D+01,   0.580168D+01,   0.581146D+01,   0.582073D+01,
     .   0.582955D+01,   0.583794D+01,   0.584603D+01,   0.585368D+01,
     .   0.586107D+01,   0.586845D+01,   0.587488D+01,   0.588142D+01/

      DATA (YHI(I),I=321,360)/
     .   0.588778D+01,   0.589372D+01,   0.589955D+01,   0.590516D+01,
     .   0.591063D+01,   0.591585D+01,   0.592085D+01,   0.592578D+01,
     .   0.593078D+01,   0.593517D+01,   0.593951D+01,   0.595985D+01,
     .   0.597744D+01,   0.599283D+01,   0.600662D+01,   0.601885D+01,
     .   0.602986D+01,   0.603990D+01,   0.604897D+01,   0.605729D+01,
     .   0.606493D+01,   0.607199D+01,   0.607853D+01,   0.608461D+01,
     .   0.609027D+01,   0.609557D+01,   0.610053D+01,   0.610519D+01,
     .   0.610958D+01,   0.611373D+01,   0.611764D+01,   0.612135D+01,
     .   0.612487D+01,   0.612821D+01,   0.613140D+01,   0.613443D+01,
     .   0.613732D+01,   0.614009D+01,   0.614273D+01,   0.614526D+01/

      DATA (YHI(I),I=361,400)/
     .   0.614769D+01,   0.615226D+01,   0.615648D+01,   0.616040D+01,
     .   0.616404D+01,   0.616744D+01,   0.617063D+01,   0.617361D+01,
     .   0.617642D+01,   0.617906D+01,   0.618156D+01,   0.618392D+01,
     .   0.618616D+01,   0.618828D+01,   0.619030D+01,   0.619222D+01,
     .   0.619405D+01,   0.619580D+01,   0.619747D+01,   0.619907D+01,
     .   0.620060D+01,   0.620207D+01,   0.620348D+01,   0.620483D+01,
     .   0.620613D+01,   0.620739D+01,   0.620859D+01,   0.620975D+01,
     .   0.621088D+01,   0.621196D+01,   0.621301D+01,   0.621775D+01,
     .   0.622181D+01,   0.622534D+01,   0.622843D+01,   0.623117D+01,
     .   0.623361D+01,   0.623580D+01,   0.623780D+01,   0.623960D+01/

      DATA (YHI(I),I=401,440)/
     .   0.624124D+01,   0.624284D+01,   0.624423D+01,   0.624552D+01,
     .   0.624671D+01,   0.624782D+01,   0.624886D+01,   0.624983D+01,
     .   0.625074D+01,   0.625159D+01,   0.625240D+01,   0.625316D+01,
     .   0.625387D+01,   0.625456D+01,   0.625520D+01,   0.625582D+01,
     .   0.625640D+01,   0.625696D+01,   0.625749D+01,   0.625799D+01,
     .   0.625848D+01,   0.625939D+01,   0.626023D+01,   0.626100D+01,
     .   0.626172D+01,   0.626238D+01,   0.626300D+01,   0.626358D+01,
     .   0.626413D+01,   0.626464D+01,   0.626512D+01,   0.626558D+01,
     .   0.626601D+01,   0.626641D+01,   0.626680D+01,   0.626716D+01,
     .   0.626751D+01,   0.626784D+01,   0.626816D+01,   0.626846D+01/

      DATA (YHI(I),I=441,461)/
     .   0.626875D+01,   0.626903D+01,   0.626929D+01,   0.626954D+01,
     .   0.626979D+01,   0.627002D+01,   0.627025D+01,   0.627046D+01,
     .   0.627067D+01,   0.627087D+01,   0.627107D+01,   0.627194D+01,
     .   0.627269D+01,   0.627333D+01,   0.627389D+01,   0.627439D+01,
     .   0.627483D+01,   0.627522D+01,   0.627557D+01,   0.627590D+01,
     .   0.627620D+01/

      END

      COMPLEX*16 FUNCTION CKOF(RHO)
C--COMPLEX COEFFICIENT OF VIRTUAL CORRECTIONS
C--RHO = (MH/MQ)**2
      PARAMETER(NN=999,N=461,NCUT1=46,NCUT2=54,NCUT3=341)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION XC(NN),YCR(NN),YCI(NN)
      COMMON/RESINT/XX(NN),YHR(NN),YHI(NN),YAR(NN),YAI(NN)
      COMMON/HIGGS/IHIGGS
      COMMON/CUT/EPST,EPSV,REPS
      PI=4.D0*DATAN(1.D0)
      RHO2=XX(NCUT3)
      RHEP=1.D-15
 
C--CALCULATE COMPLEX COEFFICIENT BY INTERPOLATION
      IF(IHIGGS.EQ.1)THEN
C--PSEUDOSCALAR HIGGS
       RHO0=XX(NCUT1)
       RHO1=XX(NCUT2)
       IF(RHO.LE.RHO0.OR.RHO.GE.RHO1)THEN
        IF(RHO.LE.RHO2)THEN
         CKOF=DCMPLX(FINT1(RHO,XX,YAR,N),FINT1(RHO,XX,YAI,N))
        ELSE
         CRHO=RHO2/DCMPLX(1.D0,-RHEP)
         CDLR=CDLOG(-CRHO)
         CLIM=5.D0/36.D0*CDLR**2-4.D0/3.D0*CDLR
         CONST=DCMPLX(YAR(NCUT3),YAI(NCUT3))-CLIM
         CRHO=RHO/DCMPLX(1.D0,-REPS)
         CDLR=CDLOG(-CRHO)
         CLIM=5.D0/36.D0*CDLR**2-4.D0/3.D0*CDLR
         CKOF=CLIM+CONST
        ENDIF
       ELSEIF(RHO.LT.4.D0.AND.RHO.GT.RHO0)THEN
        TAU=4.D0/RHO0
        CONST=DCMPLX(YAR(NCUT1)+16.D0/3.D0*DLOG((TAU-1.D0)/TAU),
     .               YAI(NCUT1))
        CTAU=DCMPLX(4.D0/RHO,0.D0)
        CKOF=-16.D0/3.D0*CDLOG((CTAU-1.D0)/CTAU)+CONST
       ELSEIF(RHO.GT.4.D0.AND.RHO.LT.RHO1)THEN
        TAU=4.D0/RHO1
        BETA=DSQRT(1.D0-TAU)
        CCFF0=DCMPLX(1.D0,4.D0/PI*BETA)
        CLIM=-16.D0/3.D0*DCMPLX(DLOG((1.D0-TAU)/TAU),-PI)/CCFF0
        CONST=DCMPLX(YAR(NCUT2),YAI(NCUT2))-CLIM
        TAU=4.D0/RHO
        BETA=DSQRT(1.D0-TAU)
        CCFF0=DCMPLX(1.D0,4.D0/PI*BETA)
        CLIM=-16.D0/3.D0*DCMPLX(DLOG((1.D0-TAU)/TAU),-PI)/CCFF0
        CKOF=CLIM+CONST
       ELSEIF(RHO.EQ.4.D0)THEN
        CKOF=DCMPLX(1.D5,16.D0/3.D0*PI)
       ENDIF
      ELSE
C--SCALAR HIGGS
       IF(RHO.LE.RHO2)THEN
        CKOF=DCMPLX(FINT1(RHO,XX,YHR,N),FINT1(RHO,XX,YHI,N))
       ELSE
        CRHO=RHO2/DCMPLX(1.D0,-RHEP)
        CDLR=CDLOG(-CRHO)
        CLIM=5.D0/36.D0*CDLR**2-4.D0/3.D0*CDLR
        CONST=DCMPLX(YHR(NCUT3),YHI(NCUT3))-CLIM
        CRHO=RHO/DCMPLX(1.D0,-REPS)
        CDLR=CDLOG(-CRHO)
        CLIM=5.D0/36.D0*CDLR**2-4.D0/3.D0*CDLR
        CKOF=CLIM+CONST
       ENDIF
      ENDIF
C--SUBTRACT UNIVERSAL PI**2-TERM
      CKOF=CKOF-DCMPLX(PI**2,0.D0)
      RETURN
      END
 
      COMPLEX*16 FUNCTION CKOFS(RHO)
C--COMPLEX COEFFICIENT OF VIRTUAL CORRECTIONS
C--RHO = (MH/MQ)**2
c     PARAMETER(NN=999,N=461,NCUT1=46,NCUT2=54,NCUT3=341)
      PARAMETER(NN=999,N=461,NCUT1=46,NCUT2=54,NCUT3=461)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION XC(NN),YCR(NN),YCI(NN)
      COMMON/RESINTS/XX(NN),YHR(NN),YHI(NN),YAR(NN),YAI(NN)
      COMMON/HIGGS/IHIGGS
      COMMON/CUT/EPST,EPSV,REPS
      PI=4.D0*DATAN(1.D0)
      RHO2=XX(NCUT3)
      RHEP=1.D-15
 
C--CALCULATE COMPLEX COEFFICIENT BY INTERPOLATION
      IF(IHIGGS.EQ.1)THEN
C--PSEUDOSCALAR HIGGS
        CKOFS=0
      ELSE
C--SCALAR HIGGS
       IF(RHO.LE.RHO2)THEN
        CKOFS=DCMPLX(FINT1(RHO,XX,YHR,N),FINT1(RHO,XX,YHI,N))
       ELSE
        CRHO=RHO2/DCMPLX(1.D0,-RHEP)
        CDLR=CDLOG(-CRHO)
        CLIM=5.D0/36.D0*CDLR**2-4.D0/3.D0*CDLR
        CONST=DCMPLX(YHR(NCUT3),YHI(NCUT3))-CLIM
        CRHO=RHO/DCMPLX(1.D0,-REPS)
        CDLR=CDLOG(-CRHO)
        CLIM=5.D0/36.D0*CDLR**2-4.D0/3.D0*CDLR
        CKOFS=CLIM+CONST
       ENDIF
      ENDIF
C--SUBTRACT UNIVERSAL PI**2-TERM
      CKOFS=CKOFS-DCMPLX(PI**2,0.D0)
      RETURN
      END
 

C   IMSL ROUTINE NAME   - DCADRE
C
C-----------------------------------------------------------------------
C
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUMERICAL INTEGRATION OF A FUNCTION USING
C                           CAUTIOUS ADAPTIVE ROMBERG EXTRAPOLATION
C
C   USAGE               - FUNCTION DCADRE (F,A,B,AERR,RERR,ERROR,IER)
C
C   ARGUMENTS    DCADRE - ESTIMATE OF THE INTEGRAL OF F(X) FROM A TO B.
C                           (OUTPUT).
C                F      - A SINGLE-ARGUMENT REAL FUNCTION SUBPROGRAM
C                           SUPPLIED BY THE USER. (INPUT)
C                           F MUST BE DECLARED EXTERNAL IN THE
C                           CALLING PROGRAM.
C                A,B    - THE TWO ENDPOINTS OF THE INTERVAL OF
C                           INTEGRATION. (INPUT)
C                AERR   - DESIRED ABSOLUTE ERROR IN THE ANSWER. (INPUT)
C                RERR   - DESIRED RELATIVE ERROR IN THE ANSWER. (INPUT)
C                ERROR  - ESTIMATED BOUND ON THE ABSOLUTE ERROR OF
C                           THE OUTPUT NUMBER, DCADRE. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR(WITH FIX)
C                           IER = 65 IMPLIES THAT ONE OR MORE
C                             SINGULARITIES WERE SUCCESSFULLY HANDLED.
C                           IER = 66 IMPLIES THAT, IN SOME
C                             SUBINTERVAL(S), THE ESTIMATE OF THE
C                             INTEGRAL WAS ACCEPTED MERELY BECAUSE THE
C                             ESTIMATED ERROR WAS SMALL, EVEN THOUGH NO
C                             REGULAR BEHAVIOR WAS RECOGNIZED.
C                         TERMINAL ERROR
C                           IER = 131 INDICATES FAILURE DUE TO
C                             INSUFFICIENT INTERNAL WORKING STORAGE.
C                           IER = 132 INDICATES FAILURE DUE TO
C                             TOO MUCH NOISE IN THE FUNCTION (RELATIVE
C                             TO THE GIVEN ERROR REQUIREMENTS) OR
C                             DUE TO AN ILL-BEHAVED INTEGRAND.
C                           IER = 133 INDICATES THAT RERR IS GREATER
C                             THAN 0.1, OR RERR IS LESS THAN 0.0, OR
C                             RERR IS TOO SMALL FOR THE PRECISION OF
C                             THE MACHINE.
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  DCADRE CAN, IN MANY CASES, HANDLE JUMP
C                DISCONTINUITIES. SEE DOCUMENT REFERENCE FOR FULL
C                DETAILS.
C            2.  THE RELATIVE ERROR PARAMETER RERR MUST BE IN THE
C                INTERVAL (0.0,0.1) INCLUSIVELY. FOR EXAMPLE,
C                RERR = 0.1 INDICATES THAT THE ESTIMATE OF THE
C                INTEGRAL IS TO BE CORRECT TO ONE DIGIT, WHEREAS
C                RERR = .0001 CALLS FOR FOUR DIGITS OF ACCURACY.
C                IF DCADRE DETERMINES THAT THE RELATIVE ACCURACY
C                REQUIREMENTS CANNOT BE SATISFIED, IER IS SET TO
C                133 (RERR SHOULD BE LARGE ENOUGH THAT, WHEN ADDED
C                TO 100.0, THE RESULT IS A NUMBER GREATER THAN
C                100.0).
C            3.  THE ABSOLUTE ERROR PARAMETER, AERR, SHOULD BE NON-
C                NEGATIVE. IN ORDER TO GIVE A REASONABLE VALUE FOR
C                AERR, THE USER MUST KNOW THE APPROXIMATE MAGNITUDE
C                OF THE INTEGRAL BEING COMPUTED. IN MANY CASES IT IS
C                SATISFACTORY TO USE AERR = 0.0. IN THIS CASE, ONLY
C                THE RELATIVE ERROR REQUIREMENT IS SATISFIED IN THE
C                COMPUTATION.
C            4.  EVEN WHEN IER IS NOT EQUAL TO 0, DCADRE RETURNS THE
C                BEST ESTIMATE THAT HAS BEEN COMPUTED.
C                QUOTING FROM THE DOCUMENT REFERENCE- A VERY CAUTIOUS
C                MAN WOULD ACCEPT DCADRE ONLY IF IER IS 0 OR 65. THE
C                MERELY REASONABLE MAN WOULD KEEP THE FAITH EVEN IF
C                IER IS 66. THE ADVENTUROUS MAN IS QUITE OFTEN RIGHT
C                IN ACCEPTING DCADRE EVEN IF IER IS 131 OR 132.
C            5.  DCADRE MAY RETURN WRONG ANSWERS IF F HAS A PERIODIC
C                FACTOR WITH HIGH FREQUENCY AND THE INTERVAL (A,B)
C                CONTAINS AN INTEGRAL NUMBER OF PERIODS. IN THIS CASE
C                THE EASIEST FIX IS TO DIVIDE THE INTERVAL INTO TWO
C                SUBINTERVALS (A,C) AND (C,B) SUCH THAT NEITHER
C                CONTAINS AN INTEGRAL NUMBER OF PERIODS (PICK C AT
C                RANDOM), AND CALL DCADRE TO INTEGRATE OVER EACH
C                SUBINTERVAL.
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DCADR1 (F,A,B,AERR,RERR,ERROR,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      REAL*8             F,A,B,AERR,RERR,ERROR
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IBEGS(30),MAXTS,MAXTBL,MXSTGE,IBEG,II,NNLEFT
      INTEGER            I,N2,III,ISTEP2,IEND,ISTEP,L,LM1,IT,ISTAGE,N
      REAL*8             T(10,10),R(10),AIT(10),DIF(10),RN(4),TS(2049)
      REAL*8             BEGIN(30),FINIS(30),EST(30)
      REAL*8             H2TOL,AITTOL,LENGTH,JUMPTL,ZERO,P1,HALF,ONE
      REAL*8             TWO,FOUR,FOURP5,TEN,HUN,CADRE,AITLOW
      REAL*8             STEPMN,STEPNM,STAGE,CUREST,FNSIZE,HRERR
      REAL*8             PREVER,BEG,FBEG,EDN,FEND,STEP,ASTEP,TABS,HOVN
      REAL*8             FN,SUM,SUMABS,VINT,TABTLM,ERGL,ERGOAL
      REAL*8             ERRA,ERRR,FEXTRP,ERRER,DIFF,SING,FEXTM1
      REAL*8             H2NEXT,SINGNX,SLOPE,FBEG2,ERRET,H2TFEX,FI
      LOGICAL            H2CONV,AITKEN,RIGHT,REGLAR,REGLSV(30)
      EXTERNAL           F
      DATA               AITLOW,H2TOL,AITTOL,JUMPTL,MAXTS,MAXTBL,MXSTGE
     1                   /1.1D0,.15D0,.1D0,.01D0,2049,10,30/
      DATA               RN(1),RN(2),RN(3),RN(4)/
     1                   .7142005D0,.3466282D0,.843751D0,.1263305D0/
      DATA               ZERO,P1,HALF,ONE,TWO,FOUR,FOURP5,TEN,HUN
     1                   /0.0D0,0.1D0,0.5D0,1.0D0,2.0D0,4.0D0,
     2                   4.5D0,10.0D0,100.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      CADRE = ZERO
      ERROR = ZERO
      CUREST = ZERO
      VINT = ZERO
      LENGTH = DABS(B-A)
      IF (LENGTH .EQ. ZERO) GO TO 215
      IF (RERR .GT. P1 .OR. RERR .LT. ZERO) GO TO 210
      HRERR = RERR+HUN
      IF (AERR .EQ. ZERO .AND. HRERR .LE. HUN) GO TO 210
      ERRR = RERR
      ERRA = DABS(AERR)
      STEPMN = LENGTH/(TWO**MXSTGE)
      STEPNM = DMAX1(LENGTH,DABS(A),DABS(B))*TEN
      STAGE = HALF
      ISTAGE = 1
      FNSIZE = ZERO
      PREVER = ZERO
      REGLAR = .FALSE.
C                                  THE GIVEN INTERVAL OF INTEGRATION
C                                    IS THE FIRST INTERVAL CONSIDERED.
      BEG = A
      FBEG = F(BEG)*HALF
      TS(1) = FBEG
      IBEG = 1
      EDN = B
      FEND = F(EDN)*HALF
      TS(2) = FEND
      IEND = 2
    5 RIGHT = .FALSE.
C                                  INVESTIGATION OF A PARTICULAR
C                                    SUBINTERVAL BEGINS AT THIS POINT.
   10 STEP = EDN - BEG
      ASTEP =  DABS(STEP)
      IF (ASTEP .LT. STEPMN) GO TO 205
      HRERR = STEPNM+ASTEP
      IF (HRERR .EQ. STEPNM) GO TO 205
      T(1,1) = FBEG + FEND
      TABS = DABS(FBEG) + DABS(FEND)
      L = 1
      N = 1
      H2CONV = .FALSE.
      AITKEN = .FALSE.
   15 LM1 = L
      L = L + 1
C                                  CALCULATE THE NEXT TRAPEZOID SUM,
C                                    T(L,1), WHICH IS BASED ON *N2* + 1
C                                    EQUISPACED POINTS. HERE,
C                                    N2 = N*2 = 2**(L-1).
      N2 = N+N
      FN = N2
      ISTEP = (IEND - IBEG)/N
      IF (ISTEP .GT. 1) GO TO 25
      II = IEND
      IEND = IEND + N
      IF (IEND .GT. MAXTS) GO TO 200
      HOVN = STEP/FN
      III = IEND
      FI = ONE
      DO 20 I=1,N2,2
         TS(III) = TS(II)
         TS(III-1) = F(EDN - FI * HOVN)
         FI = FI+TWO
         III = III-2
         II = II-1
   20 CONTINUE
      ISTEP = 2
   25 ISTEP2 = IBEG + ISTEP/2
      SUM = ZERO
      SUMABS = ZERO
      DO 30 I=ISTEP2,IEND,ISTEP
         SUM = SUM + TS(I)
         SUMABS = SUMABS + DABS(TS(I))
   30 CONTINUE
      T(L,1) = T(L-1,1)*HALF+SUM/FN
      TABS = TABS*HALF+SUMABS/FN
      N = N2
C                                  GET PRELIMINARY VALUE FOR *VINT*
C                                    FROM LAST TRAPEZOID SUM AND UPDATE
C                                    THE ERROR REQUIREMENT *ERGOAL*
C                                    FOR THIS SUBINTERVAL.
      IT = 1
      VINT = STEP*T(L,1)
      TABTLM = TABS*TEN
      FNSIZE = DMAX1(FNSIZE,DABS(T(L,1)))
      ERGL = ASTEP*FNSIZE*TEN
      ERGOAL = STAGE*DMAX1(ERRA,ERRR*DABS(CUREST+VINT))
C                                  COMPLETE ROW L AND COLUMN L OF *T*
C                                    ARRAY.
      FEXTRP = ONE
      DO 35 I=1,LM1
         FEXTRP = FEXTRP*FOUR
         T(I,L) = T(L,I) - T(L-1,I)
         T(L,I+1) = T(L,I) + T(I,L)/(FEXTRP-ONE)
   35 CONTINUE
      ERRER = ASTEP*DABS(T(1,L))
C                                  PRELIMINARY DECISION PROCEDURE
C                                    IF L = 2 AND T(2,1) = T(1,1),
C                                    GO TO 135 TO FOLLOW UP THE
C                                    IMPRESSION THAT INTERGRAND IS
C                                    STRAIGHT LINE.
      IF (L .GT. 2) GO TO 40
      HRERR = TABS+P1*DABS(T(1,2))
      IF (HRERR .EQ. TABS) GO TO 135
      GO TO 15
C                                  CACULATE NEXT RATIOS FOR
C                                    COLUMNS 1,...,L-2 OF T-TABLE
C                                    RATIO IS SET TO ZERO IF DIFFERENCE
C                                    IN LAST TWO ENTRIES OF COLUMN IS
C                                    ABOUT ZERO
   40 DO 45 I=2,LM1
         DIFF = ZERO
         HRERR = TABTLM+DABS(T(I-1,L))
         IF (HRERR .NE. TABTLM) DIFF = T(I-1,LM1)/T(I-1,L)
         T(I-1,LM1) = DIFF
   45 CONTINUE
      IF (DABS(FOUR-T(1,LM1)) .LE. H2TOL) GO TO 60
      IF (T(1,LM1) .EQ. ZERO) GO TO 55
      IF (DABS(TWO-DABS(T(1,LM1))) .LT. JUMPTL) GO TO 130
      IF (L .EQ. 3) GO TO 15
      H2CONV = .FALSE.
      IF (DABS((T(1,LM1)-T(1,L-2))/T(1,LM1)) .LE. AITTOL) GO TO 75
   50 IF (REGLAR) GO TO 55
      IF (L .EQ. 4) GO TO 15
      HRERR = ERGL+ERRER
   55 IF (ERRER .GT. ERGOAL .AND. HRERR .NE. ERGL) GO TO 175
      GO TO 145
C                                  CAUTIOUS ROMBERG EXTRAPOLATION
   60 IF (H2CONV) GO TO 65
      AITKEN = .FALSE.
      H2CONV = .TRUE.
   65 FEXTRP = FOUR
   70 IT = IT + 1
      VINT = STEP*T(L,IT)
      ERRER = DABS(STEP/(FEXTRP-ONE)*T(IT-1,L))
      IF (ERRER .LE. ERGOAL) GO TO 160
      HRERR = ERGL+ERRER
      IF (HRERR .EQ. ERGL) GO TO 160
      IF (IT .EQ. LM1) GO TO 125
      IF (T(IT,LM1) .EQ. ZERO) GO TO 70
      IF (T(IT,LM1) .LE. FEXTRP) GO TO 125
      IF (DABS(T(IT,LM1)/FOUR-FEXTRP)/FEXTRP .LT. AITTOL)
     1       FEXTRP = FEXTRP*FOUR
      GO TO 70
C                                  INTEGRAND MAY HAVE X**ALPHA TYPE
C                                    SINGULARITY
C                                    RESULTING IN A RATIO OF *SING*  =
C                                    2**(ALPHA + 1)
   75 IF (T(1,LM1) .LT. AITLOW) GO TO 175
      IF (AITKEN) GO TO 80
      H2CONV = .FALSE.
      AITKEN = .TRUE.
   80 FEXTRP = T(L-2,LM1)
      IF (FEXTRP .GT. FOURP5) GO TO 65
      IF (FEXTRP .LT. AITLOW) GO TO 175
      IF (DABS(FEXTRP-T(L-3,LM1))/T(1,LM1) .GT. H2TOL) GO TO 175
      SING = FEXTRP
      FEXTM1 = ONE/(FEXTRP - ONE)
      AIT(1) = ZERO
      DO 85 I=2,L
         AIT(I) = T(I,1) + (T(I,1)-T(I-1,1))*FEXTM1
         R(I) = T(1,I-1)
         DIF(I) = AIT(I) - AIT(I-1)
   85 CONTINUE
      IT = 2
   90 VINT = STEP*AIT(L)
      ERRER = ERRER*FEXTM1
      HRERR = ERGL+ERRER
      IF (ERRER .GT. ERGOAL .AND. HRERR .NE. ERGL) GO TO 95
      IER = MAX0(IER,65)
      GO TO 160
   95 IT = IT + 1
      IF (IT .EQ. LM1) GO TO 125
      IF (IT .GT. 3) GO TO 100
      H2NEXT = FOUR
      SINGNX = SING+SING
  100 IF (H2NEXT .LT. SINGNX) GO TO 105
      FEXTRP = SINGNX
      SINGNX = SINGNX+SINGNX
      GO TO 110
  105 FEXTRP = H2NEXT
      H2NEXT = FOUR*H2NEXT
  110 DO 115 I=IT,LM1
         R(I+1) = ZERO
         HRERR = TABTLM+DABS(DIF(I+1))
         IF (HRERR .NE. TABTLM) R(I+1) = DIF(I)/DIF(I+1)
  115 CONTINUE
      H2TFEX = -H2TOL*FEXTRP
      IF (R(L) - FEXTRP .LT. H2TFEX) GO TO 125
      IF (R(L-1)-FEXTRP .LT. H2TFEX) GO TO 125
      ERRER = ASTEP*DABS(DIF(L))
      FEXTM1 = ONE/(FEXTRP - ONE)
      DO 120 I=IT,L
         AIT(I) = AIT(I) + DIF(I)*FEXTM1
         DIF(I) = AIT(I) - AIT(I-1)
  120 CONTINUE
      GO TO 90
C                                  CURRENT TRAPEZOID SUM AND RESULTING
C                                    EXTRAPOLATED VALUES DID NOT GIVE
C                                    A SMALL ENOUGH *ERRER*.
C                                    NOTE -- HAVING PREVER .LT. ERRER
C                                    IS AN ALMOST CERTAIN SIGN OF
C                                    BEGINNING TROUBLE WITH IN THE FUNC-
C                                    TION VALUES. HENCE, A WATCH FOR,
C                                    AND CONTROL OF, NOISE SHOULD
C                                    BEGIN HERE.
  125 FEXTRP = DMAX1(PREVER/ERRER,AITLOW)
      PREVER = ERRER
      IF (L .LT. 5) GO TO 15
      IF (L-IT .GT. 2 .AND. ISTAGE .LT. MXSTGE) GO TO 170
      ERRET = ERRER/(FEXTRP**(MAXTBL-L))
      HRERR = ERGL+ERRET
      IF (ERRET .GT. ERGOAL .AND. HRERR .NE. ERGL) GO TO 170
      GO TO 15
C                                  INTEGRAND HAS JUMP (SEE NOTES)
  130 HRERR = ERGL+ERRER
      IF (ERRER .GT. ERGOAL .AND. HRERR .NE. ERGL) GO TO 170
C                                  NOTE THAT  2*FN = 2**L
      DIFF = DABS(T(1,L))*(FN+FN)
      GO TO 160
C                                  INTEGRAND IS STRAIGHT LINE
C                                    TEST THIS ASSUMPTION BY COMPARING
C                                    THE VALUE OF THE INTEGRAND AT
C                                    FOUR *RANDOMLY CHOSEN* POINTS WITH
C                                    THE VALUE OF THE STRAIGHT LINE
C                                    INTERPOLATING THE INTEGRAND AT THE
C                                    TWO END POINTS OF THE SUB-INTERVAL.
C                                    IF TEST IS PASSED, ACCEPT *VINT*
  135 SLOPE = (FEND-FBEG)*TWO
      FBEG2 = FBEG+FBEG
      DO 140 I=1,4
         DIFF = DABS(F(BEG+RN(I)*STEP) - FBEG2-RN(I)*SLOPE)
         HRERR = TABTLM+DIFF
         IF(HRERR .NE. TABTLM) GO TO 155
  140 CONTINUE
      GO TO 160
C                                  NOISE MAY BE DOMINANT FEATURE
C                                    ESTIMATE NOISE LEVEL BY COMPARING
C                                    THE VALUE OF THE INTEGRAND AT
C                                    FOUR *RANDOMLY CHOSEN* POINTS WITH
C                                    THE VALUE OF THE STRAIGHT LINE
C                                    INTERPOLATING THE INTEGRAND AT THE
C                                    TWO ENDPOINTS. IF SMALL ENOUGH,
C                                    ACCEPT *VINT*
  145 SLOPE = (FEND-FBEG)*TWO
      FBEG2 = FBEG+FBEG
      I = 1
  150 DIFF = DABS(F(BEG+RN(I)*STEP) - FBEG2-RN(I)*SLOPE)
  155 ERRER = DMAX1(ERRER,ASTEP*DIFF)
      HRERR = ERGL+ERRER
      IF (ERRER .GT. ERGOAL .AND. HRERR .NE. ERGL) GO TO 175
      I = I+1
      IF (I .LE. 4) GO TO 150
      IER = 66
C                                  INTERGRATION OVER CURRENT SUB-
C                                    INTERVAL SUCCESSFUL
C                                    ADD *VINT* TO *DCADRE* AND *ERRER*
C                                    TO *ERROR*, THEN SET UP NEXT SUB-
C                                    INTERVAL, IF ANY.
  160 CADRE = CADRE + VINT
      ERROR = ERROR + ERRER
      IF (RIGHT) GO TO 165
      ISTAGE = ISTAGE - 1
      IF (ISTAGE .EQ. 0) GO TO 220
      REGLAR = REGLSV(ISTAGE)
      BEG = BEGIN(ISTAGE)
      EDN = FINIS(ISTAGE)
      CUREST = CUREST - EST(ISTAGE+1) + VINT
      IEND = IBEG - 1
      FEND = TS(IEND)
      IBEG = IBEGS(ISTAGE)
      GO TO 180
  165 CUREST = CUREST + VINT
      STAGE = STAGE+STAGE
      IEND = IBEG
      IBEG = IBEGS(ISTAGE)
      EDN = BEG
      BEG = BEGIN(ISTAGE)
      FEND = FBEG
      FBEG = TS(IBEG)
      GO TO 5
C                                  INTEGRATION OVER CURRENT SUBINTERVAL
C                                    IS UNSUCCESSFUL. MARK SUBINTERVAL
C                                    FOR FURTHER SUBDIVISION. SET UP
C                                    NEXT SUBINTERVAL.
  170 REGLAR = .TRUE.
  175 IF (ISTAGE .EQ. MXSTGE) GO TO 205
      IF (RIGHT) GO TO 185
      REGLSV(ISTAGE+1) = REGLAR
      BEGIN(ISTAGE) = BEG
      IBEGS(ISTAGE) = IBEG
      STAGE = STAGE*HALF
  180 RIGHT = .TRUE.
      BEG = (BEG+EDN)*HALF
      IBEG = (IBEG+IEND)/2
      TS(IBEG) = TS(IBEG)*HALF
      FBEG = TS(IBEG)
      GO TO 10
  185 NNLEFT = IBEG - IBEGS(ISTAGE)
      IF (IEND+NNLEFT .GE. MAXTS) GO TO 200
      III = IBEGS(ISTAGE)
      II = IEND
      DO 190 I=III,IBEG
         II = II + 1
         TS(II) = TS(I)
  190 CONTINUE
      DO 195 I=IBEG,II
         TS(III) = TS(I)
         III = III + 1
  195 CONTINUE
      IEND = IEND + 1
      IBEG = IEND - NNLEFT
      FEND = FBEG
      FBEG = TS(IBEG)
      FINIS(ISTAGE) = EDN
      EDN = BEG
      BEG = BEGIN(ISTAGE)
      BEGIN(ISTAGE) = EDN
      REGLSV(ISTAGE) = REGLAR
      ISTAGE = ISTAGE + 1
      REGLAR = REGLSV(ISTAGE)
      EST(ISTAGE) = VINT
      CUREST = CUREST + EST(ISTAGE)
      GO TO 5
C                                  FAILURE TO HANDLE GIVEN INTEGRA-
C                                    TION PROBLEM
  200 IER = 131
      GO TO 215
  205 IER = 132
      GO TO 215
  210 IER = 133
  215 CADRE = CUREST + VINT
  220 DCADR1 = CADRE
 9000 CONTINUE
cccccccccccccccccc     IF (IER .NE. 0) CALL UERTST (IER,6HDCADR1)
 9005 RETURN
      END
 
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DCADR2 (F,A,B,AERR,RERR,ERROR,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      REAL*8             F,A,B,AERR,RERR,ERROR
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IBEGS(30),MAXTS,MAXTBL,MXSTGE,IBEG,II,NNLEFT
      INTEGER            I,N2,III,ISTEP2,IEND,ISTEP,L,LM1,IT,ISTAGE,N
      REAL*8             T(10,10),R(10),AIT(10),DIF(10),RN(4),TS(2049)
      REAL*8             BEGIN(30),FINIS(30),EST(30)
      REAL*8             H2TOL,AITTOL,LENGTH,JUMPTL,ZERO,P1,HALF,ONE
      REAL*8             TWO,FOUR,FOURP5,TEN,HUN,CADRE,AITLOW
      REAL*8             STEPMN,STEPNM,STAGE,CUREST,FNSIZE,HRERR
      REAL*8             PREVER,BEG,FBEG,EDN,FEND,STEP,ASTEP,TABS,HOVN
      REAL*8             FN,SUM,SUMABS,VINT,TABTLM,ERGL,ERGOAL
      REAL*8             ERRA,ERRR,FEXTRP,ERRER,DIFF,SING,FEXTM1
      REAL*8             H2NEXT,SINGNX,SLOPE,FBEG2,ERRET,H2TFEX,FI
      LOGICAL            H2CONV,AITKEN,RIGHT,REGLAR,REGLSV(30)
      EXTERNAL           F
      DATA               AITLOW,H2TOL,AITTOL,JUMPTL,MAXTS,MAXTBL,MXSTGE
     1                   /1.1D0,.15D0,.1D0,.01D0,2049,10,30/
      DATA               RN(1),RN(2),RN(3),RN(4)/
     1                   .7142005D0,.3466282D0,.843751D0,.1263305D0/
      DATA               ZERO,P1,HALF,ONE,TWO,FOUR,FOURP5,TEN,HUN
     1                   /0.0D0,0.1D0,0.5D0,1.0D0,2.0D0,4.0D0,
     2                   4.5D0,10.0D0,100.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      CADRE = ZERO
      ERROR = ZERO
      CUREST = ZERO
      VINT = ZERO
      LENGTH = DABS(B-A)
      IF (LENGTH .EQ. ZERO) GO TO 215
      IF (RERR .GT. P1 .OR. RERR .LT. ZERO) GO TO 210
      HRERR = RERR+HUN
      IF (AERR .EQ. ZERO .AND. HRERR .LE. HUN) GO TO 210
      ERRR = RERR
      ERRA = DABS(AERR)
      STEPMN = LENGTH/(TWO**MXSTGE)
      STEPNM = DMAX1(LENGTH,DABS(A),DABS(B))*TEN
      STAGE = HALF
      ISTAGE = 1
      FNSIZE = ZERO
      PREVER = ZERO
      REGLAR = .FALSE.
C                                  THE GIVEN INTERVAL OF INTEGRATION
C                                    IS THE FIRST INTERVAL CONSIDERED.
      BEG = A
      FBEG = F(BEG)*HALF
      TS(1) = FBEG
      IBEG = 1
      EDN = B
      FEND = F(EDN)*HALF
      TS(2) = FEND
      IEND = 2
    5 RIGHT = .FALSE.
C                                  INVESTIGATION OF A PARTICULAR
C                                    SUBINTERVAL BEGINS AT THIS POINT.
   10 STEP = EDN - BEG
      ASTEP =  DABS(STEP)
      IF (ASTEP .LT. STEPMN) GO TO 205
      HRERR = STEPNM+ASTEP
      IF (HRERR .EQ. STEPNM) GO TO 205
      T(1,1) = FBEG + FEND
      TABS = DABS(FBEG) + DABS(FEND)
      L = 1
      N = 1
      H2CONV = .FALSE.
      AITKEN = .FALSE.
   15 LM1 = L
      L = L + 1
C                                  CALCULATE THE NEXT TRAPEZOID SUM,
C                                    T(L,1), WHICH IS BASED ON *N2* + 1
C                                    EQUISPACED POINTS. HERE,
C                                    N2 = N*2 = 2**(L-1).
      N2 = N+N
      FN = N2
      ISTEP = (IEND - IBEG)/N
      IF (ISTEP .GT. 1) GO TO 25
      II = IEND
      IEND = IEND + N
      IF (IEND .GT. MAXTS) GO TO 200
      HOVN = STEP/FN
      III = IEND
      FI = ONE
      DO 20 I=1,N2,2
         TS(III) = TS(II)
         TS(III-1) = F(EDN - FI * HOVN)
         FI = FI+TWO
         III = III-2
         II = II-1
   20 CONTINUE
      ISTEP = 2
   25 ISTEP2 = IBEG + ISTEP/2
      SUM = ZERO
      SUMABS = ZERO
      DO 30 I=ISTEP2,IEND,ISTEP
         SUM = SUM + TS(I)
         SUMABS = SUMABS + DABS(TS(I))
   30 CONTINUE
      T(L,1) = T(L-1,1)*HALF+SUM/FN
      TABS = TABS*HALF+SUMABS/FN
      N = N2
C                                  GET PRELIMINARY VALUE FOR *VINT*
C                                    FROM LAST TRAPEZOID SUM AND UPDATE
C                                    THE ERROR REQUIREMENT *ERGOAL*
C                                    FOR THIS SUBINTERVAL.
      IT = 1
      VINT = STEP*T(L,1)
      TABTLM = TABS*TEN
      FNSIZE = DMAX1(FNSIZE,DABS(T(L,1)))
      ERGL = ASTEP*FNSIZE*TEN
      ERGOAL = STAGE*DMAX1(ERRA,ERRR*DABS(CUREST+VINT))
C                                  COMPLETE ROW L AND COLUMN L OF *T*
C                                    ARRAY.
      FEXTRP = ONE
      DO 35 I=1,LM1
         FEXTRP = FEXTRP*FOUR
         T(I,L) = T(L,I) - T(L-1,I)
         T(L,I+1) = T(L,I) + T(I,L)/(FEXTRP-ONE)
   35 CONTINUE
      ERRER = ASTEP*DABS(T(1,L))
C                                  PRELIMINARY DECISION PROCEDURE
C                                    IF L = 2 AND T(2,1) = T(1,1),
C                                    GO TO 135 TO FOLLOW UP THE
C                                    IMPRESSION THAT INTERGRAND IS
C                                    STRAIGHT LINE.
      IF (L .GT. 2) GO TO 40
      HRERR = TABS+P1*DABS(T(1,2))
      IF (HRERR .EQ. TABS) GO TO 135
      GO TO 15
C                                  CACULATE NEXT RATIOS FOR
C                                    COLUMNS 1,...,L-2 OF T-TABLE
C                                    RATIO IS SET TO ZERO IF DIFFERENCE
C                                    IN LAST TWO ENTRIES OF COLUMN IS
C                                    ABOUT ZERO
   40 DO 45 I=2,LM1
         DIFF = ZERO
         HRERR = TABTLM+DABS(T(I-1,L))
         IF (HRERR .NE. TABTLM) DIFF = T(I-1,LM1)/T(I-1,L)
         T(I-1,LM1) = DIFF
   45 CONTINUE
      IF (DABS(FOUR-T(1,LM1)) .LE. H2TOL) GO TO 60
      IF (T(1,LM1) .EQ. ZERO) GO TO 55
      IF (DABS(TWO-DABS(T(1,LM1))) .LT. JUMPTL) GO TO 130
      IF (L .EQ. 3) GO TO 15
      H2CONV = .FALSE.
      IF (DABS((T(1,LM1)-T(1,L-2))/T(1,LM1)) .LE. AITTOL) GO TO 75
   50 IF (REGLAR) GO TO 55
      IF (L .EQ. 4) GO TO 15
      HRERR = ERGL+ERRER
   55 IF (ERRER .GT. ERGOAL .AND. HRERR .NE. ERGL) GO TO 175
      GO TO 145
C                                  CAUTIOUS ROMBERG EXTRAPOLATION
   60 IF (H2CONV) GO TO 65
      AITKEN = .FALSE.
      H2CONV = .TRUE.
   65 FEXTRP = FOUR
   70 IT = IT + 1
      VINT = STEP*T(L,IT)
      ERRER = DABS(STEP/(FEXTRP-ONE)*T(IT-1,L))
      IF (ERRER .LE. ERGOAL) GO TO 160
      HRERR = ERGL+ERRER
      IF (HRERR .EQ. ERGL) GO TO 160
      IF (IT .EQ. LM1) GO TO 125
      IF (T(IT,LM1) .EQ. ZERO) GO TO 70
      IF (T(IT,LM1) .LE. FEXTRP) GO TO 125
      IF (DABS(T(IT,LM1)/FOUR-FEXTRP)/FEXTRP .LT. AITTOL)
     1       FEXTRP = FEXTRP*FOUR
      GO TO 70
C                                  INTEGRAND MAY HAVE X**ALPHA TYPE
C                                    SINGULARITY
C                                    RESULTING IN A RATIO OF *SING*  =
C                                    2**(ALPHA + 1)
   75 IF (T(1,LM1) .LT. AITLOW) GO TO 175
      IF (AITKEN) GO TO 80
      H2CONV = .FALSE.
      AITKEN = .TRUE.
   80 FEXTRP = T(L-2,LM1)
      IF (FEXTRP .GT. FOURP5) GO TO 65
      IF (FEXTRP .LT. AITLOW) GO TO 175
      IF (DABS(FEXTRP-T(L-3,LM1))/T(1,LM1) .GT. H2TOL) GO TO 175
      SING = FEXTRP
      FEXTM1 = ONE/(FEXTRP - ONE)
      AIT(1) = ZERO
      DO 85 I=2,L
         AIT(I) = T(I,1) + (T(I,1)-T(I-1,1))*FEXTM1
         R(I) = T(1,I-1)
         DIF(I) = AIT(I) - AIT(I-1)
   85 CONTINUE
      IT = 2
   90 VINT = STEP*AIT(L)
      ERRER = ERRER*FEXTM1
      HRERR = ERGL+ERRER
      IF (ERRER .GT. ERGOAL .AND. HRERR .NE. ERGL) GO TO 95
      IER = MAX0(IER,65)
      GO TO 160
   95 IT = IT + 1
      IF (IT .EQ. LM1) GO TO 125
      IF (IT .GT. 3) GO TO 100
      H2NEXT = FOUR
      SINGNX = SING+SING
  100 IF (H2NEXT .LT. SINGNX) GO TO 105
      FEXTRP = SINGNX
      SINGNX = SINGNX+SINGNX
      GO TO 110
  105 FEXTRP = H2NEXT
      H2NEXT = FOUR*H2NEXT
  110 DO 115 I=IT,LM1
         R(I+1) = ZERO
         HRERR = TABTLM+DABS(DIF(I+1))
         IF (HRERR .NE. TABTLM) R(I+1) = DIF(I)/DIF(I+1)
  115 CONTINUE
      H2TFEX = -H2TOL*FEXTRP
      IF (R(L) - FEXTRP .LT. H2TFEX) GO TO 125
      IF (R(L-1)-FEXTRP .LT. H2TFEX) GO TO 125
      ERRER = ASTEP*DABS(DIF(L))
      FEXTM1 = ONE/(FEXTRP - ONE)
      DO 120 I=IT,L
         AIT(I) = AIT(I) + DIF(I)*FEXTM1
         DIF(I) = AIT(I) - AIT(I-1)
  120 CONTINUE
      GO TO 90
C                                  CURRENT TRAPEZOID SUM AND RESULTING
C                                    EXTRAPOLATED VALUES DID NOT GIVE
C                                    A SMALL ENOUGH *ERRER*.
C                                    NOTE -- HAVING PREVER .LT. ERRER
C                                    IS AN ALMOST CERTAIN SIGN OF
C                                    BEGINNING TROUBLE WITH IN THE FUNC-
C                                    TION VALUES. HENCE, A WATCH FOR,
C                                    AND CONTROL OF, NOISE SHOULD
C                                    BEGIN HERE.
  125 FEXTRP = DMAX1(PREVER/ERRER,AITLOW)
      PREVER = ERRER
      IF (L .LT. 5) GO TO 15
      IF (L-IT .GT. 2 .AND. ISTAGE .LT. MXSTGE) GO TO 170
      ERRET = ERRER/(FEXTRP**(MAXTBL-L))
      HRERR = ERGL+ERRET
      IF (ERRET .GT. ERGOAL .AND. HRERR .NE. ERGL) GO TO 170
      GO TO 15
C                                  INTEGRAND HAS JUMP (SEE NOTES)
  130 HRERR = ERGL+ERRER
      IF (ERRER .GT. ERGOAL .AND. HRERR .NE. ERGL) GO TO 170
C                                  NOTE THAT  2*FN = 2**L
      DIFF = DABS(T(1,L))*(FN+FN)
      GO TO 160
C                                  INTEGRAND IS STRAIGHT LINE
C                                    TEST THIS ASSUMPTION BY COMPARING
C                                    THE VALUE OF THE INTEGRAND AT
C                                    FOUR *RANDOMLY CHOSEN* POINTS WITH
C                                    THE VALUE OF THE STRAIGHT LINE
C                                    INTERPOLATING THE INTEGRAND AT THE
C                                    TWO END POINTS OF THE SUB-INTERVAL.
C                                    IF TEST IS PASSED, ACCEPT *VINT*
  135 SLOPE = (FEND-FBEG)*TWO
      FBEG2 = FBEG+FBEG
      DO 140 I=1,4
         DIFF = DABS(F(BEG+RN(I)*STEP) - FBEG2-RN(I)*SLOPE)
         HRERR = TABTLM+DIFF
         IF(HRERR .NE. TABTLM) GO TO 155
  140 CONTINUE
      GO TO 160
C                                  NOISE MAY BE DOMINANT FEATURE
C                                    ESTIMATE NOISE LEVEL BY COMPARING
C                                    THE VALUE OF THE INTEGRAND AT
C                                    FOUR *RANDOMLY CHOSEN* POINTS WITH
C                                    THE VALUE OF THE STRAIGHT LINE
C                                    INTERPOLATING THE INTEGRAND AT THE
C                                    TWO ENDPOINTS. IF SMALL ENOUGH,
C                                    ACCEPT *VINT*
  145 SLOPE = (FEND-FBEG)*TWO
      FBEG2 = FBEG+FBEG
      I = 1
  150 DIFF = DABS(F(BEG+RN(I)*STEP) - FBEG2-RN(I)*SLOPE)
  155 ERRER = DMAX1(ERRER,ASTEP*DIFF)
      HRERR = ERGL+ERRER
      IF (ERRER .GT. ERGOAL .AND. HRERR .NE. ERGL) GO TO 175
      I = I+1
      IF (I .LE. 4) GO TO 150
      IER = 66
C                                  INTERGRATION OVER CURRENT SUB-
C                                    INTERVAL SUCCESSFUL
C                                    ADD *VINT* TO *DCADRE* AND *ERRER*
C                                    TO *ERROR*, THEN SET UP NEXT SUB-
C                                    INTERVAL, IF ANY.
  160 CADRE = CADRE + VINT
      ERROR = ERROR + ERRER
      IF (RIGHT) GO TO 165
      ISTAGE = ISTAGE - 1
      IF (ISTAGE .EQ. 0) GO TO 220
      REGLAR = REGLSV(ISTAGE)
      BEG = BEGIN(ISTAGE)
      EDN = FINIS(ISTAGE)
      CUREST = CUREST - EST(ISTAGE+1) + VINT
      IEND = IBEG - 1
      FEND = TS(IEND)
      IBEG = IBEGS(ISTAGE)
      GO TO 180
  165 CUREST = CUREST + VINT
      STAGE = STAGE+STAGE
      IEND = IBEG
      IBEG = IBEGS(ISTAGE)
      EDN = BEG
      BEG = BEGIN(ISTAGE)
      FEND = FBEG
      FBEG = TS(IBEG)
      GO TO 5
C                                  INTEGRATION OVER CURRENT SUBINTERVAL
C                                    IS UNSUCCESSFUL. MARK SUBINTERVAL
C                                    FOR FURTHER SUBDIVISION. SET UP
C                                    NEXT SUBINTERVAL.
  170 REGLAR = .TRUE.
  175 IF (ISTAGE .EQ. MXSTGE) GO TO 205
      IF (RIGHT) GO TO 185
      REGLSV(ISTAGE+1) = REGLAR
      BEGIN(ISTAGE) = BEG
      IBEGS(ISTAGE) = IBEG
      STAGE = STAGE*HALF
  180 RIGHT = .TRUE.
      BEG = (BEG+EDN)*HALF
      IBEG = (IBEG+IEND)/2
      TS(IBEG) = TS(IBEG)*HALF
      FBEG = TS(IBEG)
      GO TO 10
  185 NNLEFT = IBEG - IBEGS(ISTAGE)
      IF (IEND+NNLEFT .GE. MAXTS) GO TO 200
      III = IBEGS(ISTAGE)
      II = IEND
      DO 190 I=III,IBEG
         II = II + 1
         TS(II) = TS(I)
  190 CONTINUE
      DO 195 I=IBEG,II
         TS(III) = TS(I)
         III = III + 1
  195 CONTINUE
      IEND = IEND + 1
      IBEG = IEND - NNLEFT
      FEND = FBEG
      FBEG = TS(IBEG)
      FINIS(ISTAGE) = EDN
      EDN = BEG
      BEG = BEGIN(ISTAGE)
      BEGIN(ISTAGE) = EDN
      REGLSV(ISTAGE) = REGLAR
      ISTAGE = ISTAGE + 1
      REGLAR = REGLSV(ISTAGE)
      EST(ISTAGE) = VINT
      CUREST = CUREST + EST(ISTAGE)
      GO TO 5
C                                  FAILURE TO HANDLE GIVEN INTEGRA-
C                                    TION PROBLEM
  200 IER = 131
      GO TO 215
  205 IER = 132
      GO TO 215
  210 IER = 133
  215 CADRE = CUREST + VINT
  220 DCADR2 = CADRE
 9000 CONTINUE
cccccccccccccccccc     IF (IER .NE. 0) CALL UERTST (IER,6HDCADR1)
 9005 RETURN
      END
 
      SUBROUTINE VEGAS(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/RESU/RES
      COMMON/BVEG2/NDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
     1,D(50,10),DI(50,10),NXI(50,10)
      DIMENSION XIN(50),R(50),DX(10),IA(10),KG(10),DT(10)
      DIMENSION XL(10),XU(10),QRAN(10),X(10)
      COMMON/RESULT/S1,S2,S3,S4
       EXTERNAL FXN
      DATA XL,XU/10*0.D0,10*1.D0/
      DATA NDMX/50/,ALPH/1.5D0/,ONE/1.D0/,MDS/1/
      IPR=1
      IF(NPRN.GT.0)IPR=0
      NDO=1
      DO 1 J=1,NDIM
1     XI(1,J)=ONE
      ENTRY VEGAS1(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
      NOW=IGRAPH
CS    IF(IGRAPH.GT.0)CALL INPLOT(NOW,F1,W)
      IT=0
      SI=0.D0
      SI2=SI
      SWGT=SI
      SCHI=SI
      SCALLS=SI
      ENTRY VEGAS2(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
      ND=NDMX
      NG=1
      IF(MDS.EQ.0) GO TO 2
      NG=(NCALL*0.5)**(1./NDIM)
      MDS=1
      IF((2*NG-NDMX).LT.0) GO TO 2
      MDS=-1
      NPG=NG/NDMX+1
      ND=NG/NPG
      NG=NPG*ND
2     K=NG**NDIM
      NPG=NCALL/K
      IF(NPG.LT.2)NPG=2
      CALLS=NPG*K
      DXG=ONE/NG
      DV2G=DXG**(2*NDIM)/NPG/NPG/(NPG-ONE)
      XND=ND
      NDM=ND-1
      DXG=DXG*XND
      XJAC=ONE
      DO 3 J=1,NDIM
      DX(J)=XU(J)-XL(J)
3     XJAC=XJAC*DX(J)
      IF(ND.EQ.NDO) GO TO 8
      RC=NDO/XND
      DO 7 J=1,NDIM
      K=0
      XN=0.D0
      DR=XN
      I=K
4     K=K+1
      DR=DR+ONE
      XO=XN
      XN=XI(K,J)
5     IF(RC.GT.DR) GO TO 4
      I=I+1
      DR=DR-RC
      XIN(I)=XN-(XN-XO)*DR
      IF(I.LT.NDM) GO TO 5
      DO 6  I=1,NDM
6     XI(I,J)=XIN(I)
7     XI(ND,J)=ONE
      NDO=ND
      ACC=BCC
      IF(NPRN.NE.0.AND.NPRN.NE.10)PRINT 200,NDIM,CALLS,IT,ITMX
     1,ACC,MDS,ND
8     CONTINUE
      IF(NPRN.EQ.10)PRINT 290,NDIM,CALLS,ITMX,ACC,MDS,ND
      ENTRY VEGAS3(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
9     IT=IT+1
      TI=0.D0
      TSI=TI
CS    IF(IGRAPH.GT.0)CALL REPLOT(NOW,F1,W)
      DO 10 J=1,NDIM
      KG(J)=1
      DO 10 I=1,ND
      NXI(I,J)=0
      D(I,J)=TI
10    DI(I,J)=TI
11    FB=0.D0
      F2B=FB
      K=0
12    K=K+1
      DO 121 J=1,NDIM
121   QRAN(J)=RANDM(0)
      WGT=XJAC
      DO 15 J=1,NDIM
      XN=(KG(J)-QRAN(J))*DXG+ONE
      IA(J)=XN
      IAJ=IA(J)
      IAJ1=IAJ-1
      IF(IAJ.GT.1) GO TO 13
      XO=XI(IAJ,J)
      RC=(XN-IAJ)*XO
      GO TO 14
13    XO=XI(IAJ,J)-XI(IAJ1,J)
      RC=XI(IAJ1,J)+(XN-IAJ)*XO
14    X(J)=XL(J)+RC*DX(J)
15    WGT=WGT*XO*XND
      F=FXN(X)*WGT
      F1=F/CALLS
      W=WGT/CALLS
CS    IF(IGRAPH.GT.0)CALL XPLOT(NOW,F1,W)
      F2=F**2
      FB=FB+F
      F2B=F2B+F2
      DO 16 J=1,NDIM
      IAJ=IA(J)
      NXI(IAJ,J)=NXI(IAJ,J)+1
      DI(IAJ,J)=DI(IAJ,J)+F/CALLS
16    IF(MDS.GE.0)  D(IAJ,J)=D(IAJ,J)+F2
      IF(K.LT.NPG) GO TO 12
      F2B=F2B*NPG
      F2B=SQRT(F2B)
      F2B=(F2B-FB)*(F2B+FB)
      TI=TI+FB
      TSI=TSI+F2B
      IF(MDS.GE.0) GO TO 18
      DO 17 J=1,NDIM
      IAJ=IA(J)
17    D(IAJ,J)=D(IAJ,J)+F2B
18    K=NDIM
19    KG(K)=MOD(KG(K),NG)+1
      IF(KG(K).NE.1) GO TO 11
      K=K-1
      IF(K.GT.0) GO TO 19
      TI=TI/CALLS
      TSI=TSI*DV2G
      TI2=TI*TI
      WGT=TI2/TSI
      SI=SI+TI*WGT
      SI2=SI2+TI2
      SWGT=SWGT+WGT
      SCHI=SCHI+TI2*WGT
      SCALLS=SCALLS+CALLS
      AVGI=SI/SWGT
      SD=SWGT*IT/SI2
      CHI2A=0.D0
      IF(IT.GT.1)CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/(IT-1)
      SD=ONE/SD
      SD=SQRT(SD)
      IF(NPRN.EQ.0) GO TO 21
      TSI=SQRT(TSI)
      IF(NPRN.NE.10)PRINT 201,IPR,IT,TI,TSI,AVGI,SD,CHI2A
      RES=AVGI
      IF(NPRN.EQ.10)PRINT 203,IT,TI,TSI,AVGI,SD,CHI2A
      IF(NPRN.GE.0) GO TO 21
      DO 20 J=1,NDIM
      PRINT 202,J
20    PRINT 204,(XI(I,J),DI(I,J),D(I,J),I=1,ND)
21    IF(ABS(SD/AVGI).LE.ABS(ACC).OR.IT.GE.ITMX)NOW=2
      S1=AVGI
      S2=SD
      S3=TI
      S4=TSI
CS    IF(IGRAPH.GT.0)CALL PLOTIT(NOW,F1,W)
C      DO 23 J=1,NDIM
C      XO=D(1,J)
C      XN=D(2,J)
C      D(1,J)=(XO+XN)*0.5D0
C      DT(J)=D(1,J)
C      DO 22 I=2,NDM
C      D(I,J)=XO+XN
C      XO=XN
C      XN=D(I+1,J)
C      D(I,J)=(D(I,J)+XN)/3.D0
C22    DT(J)=DT(J)+D(I,J)
C      D(ND,J)=(XN+XO)*0.5D0
C23    DT(J)=DT(J)+D(ND,J)
C-----THIS PART OF THE VEGAS-ALGORITHM IS UNSTABLE
C-----IT SHOULD BE REPLACED BY
      DO 23 J=1,NDIM
      DT(J)=0.D0
      DO 23 I=1,ND
      IF(NXI(I,J).GT.0)D(I,J)=D(I,J)/NXI(I,J)
23    DT(J)=DT(J)+D(I,J)
      DO 28 J=1,NDIM
      RC=0.D0
      DO 24 I=1,ND
      R(I)=0.D0
      IF(D(I,J).LE.0.D0)GO TO 24
      XO=DT(J)/D(I,J)
      R(I)=((XO-ONE)/XO/LOG(XO))**ALPH
24    RC=RC+R(I)
      RC=RC/XND
      K=0
      XN=0.D0
      DR=XN
      I=K
25    K=K+1
      DR=DR+R(K)
      XO=XN
      XN=XI(K,J)
26    IF(RC.GT.DR) GO TO 25
      I=I+1
      DR=DR-RC
      XIN(I)=XN-(XN-XO)*DR/R(K)
      IF(I.LT.NDM) GO TO 26
      DO 27 I=1,NDM
27    XI(I,J)=XIN(I)
28    XI(ND,J)=ONE
      IF(IT.LT.ITMX.AND.ABS(ACC).LT.ABS(SD/AVGI))GO TO 9
200   FORMAT(35H0INPUT PARAMETERS FOR VEGAS   NDIM=,I3
     1,8H  NCALL=,F8.0/28X,5H  IT=,I5,8H  ITMX =,I5/28X
     2,6H  ACC=,G9.3/28X,6H  MDS=,I3,6H   ND=,I4//)
290   FORMAT(13H0VEGAS  NDIM=,I3,8H  NCALL=,F8.0,8H  ITMX =,I5
     1,6H  ACC=,G9.3,6H  MDS=,I3,6H   ND=,I4)
201   FORMAT(/I1,20HINTEGRATION BY VEGAS/13H0ITERATION NO,I3,
     114H.   INTEGRAL =,G14.8/20X,10HSTD DEV  =,G10.4/
     234H ACCUMULATED RESULTS.   INTEGRAL =,G14.8/
     324X,10HSTD DEV  =,G10.4 / 24X,18HCHI**2 PER ITN   =,G10.4)
202   FORMAT(14H0DATA FOR AXIS,I2 / 7X,1HX,7X,10H  DELT I  ,
     12X,11H CONVCE    ,11X,1HX,7X,10H  DELT I  ,2X,11H CONVCE
     2,11X,1HX,7X,10H  DELT I  ,2X,11H CONVCE     /)
204   FORMAT(1X,3G12.4,5X,3G12.4,5X,3G12.4)
203   FORMAT(1H ,I3,G20.8,G12.4,G20.8,G12.4,G12.4)
      S1=AVGI
      S2=SD
      S3=CHI2A
      RETURN
      END

C----------------------------------------------------------------------
C  A UNIVERSAL RANDOM NUMBER GENERATOR
 
        DOUBLE PRECISION FUNCTION RANDM(IDMY)
        IMPLICIT REAL*8(A-H,O-Z)
        REAL*4 UNIV
        RANDM=DBLE(UNIV(1))
        RETURN 
        END

C ---------------------------------------------------------------------

        FUNCTION UNIV(IDUM)
        REAL U(97)
        COMMON /SET1/ U,C,CD,CM,I,J
        UNIV=U(I)-U(J)
        IF(UNIV.LT.0.) UNIV=UNIV+1.
        U(I)=UNIV
        I=I-1
        IF(I.EQ.0) I=97
        J=J-1
        IF(J.EQ.0) J=97
        C=C-CD
        IF(C.LT.0.) C=C+CM
        UNIV=UNIV-C
        IF(UNIV.LT.0.) UNIV=UNIV+1
        RETURN
        END
 
C----------------------------------------------------------------------
C INITIALIZING THE RANDOM NUMBER GENERATOR
C TO INITIALIZE CALL RSTART(12,34,56,78)


        SUBROUTINE RSTART(I,J,K,L)
        REAL U(97)
        COMMON /SET1/ U,C,CD,CM,ISTART,JSTART
        IF ((I.LT.0).OR.(I.GT.178 )) STOP 'FIRST SEED .LT.0 OR .GT.178'
        IF ((J.LT.0).OR.(J.GT.178 )) STOP 'SECOND SEED .LT.0 OR .GT.178'
        IF ((K.LT.0).OR.(K.GT.178 )) STOP 'THIRD SEED .LT.0 OR .GT.178'
        IF ((L.LT.0).OR.(L.GT.168 )) STOP 'FOURTH SEED .LT.0 OR .GT.168'
        IF ( (I.EQ.1).AND.(J.EQ.1).AND.(K.EQ.1) ) STOP
     &     'FIRST, SECOND AND THIRD SEEDS ARE ALL EQUAL TO 1'
        ISTART=97
        JSTART=33
        IDUM=I
        JDUM=J
        KDUM=K
        LDUM=L
        DO 2 II=1,97
        S=0.
        T=.5
        DO 3 JJ=1,24
          M=MOD(MOD(IDUM*JDUM,179)*K,179)
          IDUM=JDUM
          JDUM=KDUM
          KDUM=M
          LDUM=MOD(53*LDUM+1,169)
          IF(MOD(LDUM*M,64).GE.32) S=S+T
3         T=.5*T
2         U(II)=S
        C=362436./16777216.
        CD=7654321./16777216.
        CM=16777213./16777216.
        RETURN
        END

C----------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION FINT1(X,XP,YP,NP)
C--ONE-DIMENSIONAL QUADRATIC INTERPOLATION
C--X  = WANTED POINT
C--XP = ARRAY OF DISCRETE X-VALUES
C--YP = ARRAY OF DISCRETE FUNCTION-VALUES
C--NP = NUMBER OF DISCTRETE POINTS
      PARAMETER(N0=999,N1=3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XP(N0),YP(N0),X0(N1),Y0(N1),XX(N1),YY(N1)
      DO 1 I=1,NP
       IF(X.GE.XP(I)) NX=I
1     CONTINUE
      IF(NX.EQ.NP)NX=NP-1
      IF(NX.EQ.1)NX=2
      DO 2 I=1,3
        X0(I)=XP(NX+I-2)
2     CONTINUE
      DO 3 I=1,3
        Y0(I)=YP(NX+I-2)
3     CONTINUE
      FINT1=FINT(X,X0,Y0)
      RETURN
      END

      DOUBLE PRECISION FUNCTION FINT2(X,XP,YP,NP)
C--TWO-DIMENSIONAL QUADRATIC INTERPOLATION
C--X  = WANTED POINT
C--XP = ARRAY OF DISCRETE X-VALUES
C--YP = ARRAY OF DISCRETE FUNCTION-VALUES
C--NP = NUMBER OF DISCTRETE POINTS
      PARAMETER(N0=500,NN=10,N1=2,N2=3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N1),XP(N0,N1),YP(N0,NN),NP(N1),NX(N1),X0(N2),
     .          Y0(N2),XX(N2,N1),YY(N2,N2)
      DO 1 I=1,2
       DO 1 J=1,NP(I)
        IF(X(I).GE.XP(J,I)) NX(I)=J
1     CONTINUE
      DO 2 I=1,2
       IF(NX(I).EQ.NP(I))NX(I)=NP(I)-1
       IF(NX(I).EQ.1)NX(I)=2
2     CONTINUE
      DO 3 I=1,3
       DO 3 J=1,2
        XX(I,J)=XP(NX(J)+I-2,J)
3     CONTINUE
      DO 4 I=1,3
       DO 4 J=1,3
        YY(I,J)=YP(NX(1)+I-2,NX(2)+J-2)
4     CONTINUE
      XF=X(1)
      X0(1)=XX(1,1)
      X0(2)=XX(2,1)
      X0(3)=XX(3,1)
      Y0(1)=YY(1,1)
      Y0(2)=YY(2,1)
      Y0(3)=YY(3,1)
      F0=FINT(XF,X0,Y0)
      Y0(1)=YY(1,2)
      Y0(2)=YY(2,2)
      Y0(3)=YY(3,2)
      F1=FINT(XF,X0,Y0)
      Y0(1)=YY(1,3)
      Y0(2)=YY(2,3)
      Y0(3)=YY(3,3)
      F2=FINT(XF,X0,Y0)
      XF=X(2)
      X0(1)=XX(1,2)
      X0(2)=XX(2,2)
      X0(3)=XX(3,2)
      Y0(1)=F0
      Y0(2)=F1
      Y0(3)=F2
      FINT2=FINT(XF,X0,Y0)
      RETURN
      END

      DOUBLE PRECISION FUNCTION FINT(X,XX,YY)
C--ONE-DIMENSIONAL QUADRATIC INTERPOLATION
C--X  = WANTED POINT
C--XX = ARRAY OF 3 DISCRETE X-VALUES AROUND X
C--YY = ARRAY OF 3 DISCRETE FUNCTION-VALUES AROUND X
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(3),YY(3)
      COMMON/INTPOL/KORD
      X0=XX(1)
      X1=XX(2)
      X2=XX(3)
      Y0=YY(1)
      Y1=YY(2)
      Y2=YY(3)
      D1F=(Y1-Y0)/(X1-X0)
      D2F=(Y2-Y1)/(X2-X1)
      IF(X.LT.X1)THEN
        DX=X-X0
        DF=D1F
      ELSE
        DX=X-X1
        DF=D2F
      ENDIF
      GINT=Y1+DF*DX
      IF(KORD.EQ.2) THEN
       A0=(X-X1)*(X-X2)/(X0-X1)/(X0-X2)
       A1=(X-X0)*(X-X2)/(X1-X0)/(X1-X2)
       A2=(X-X0)*(X-X1)/(X2-X0)/(X2-X1)
       GINT=A0*Y0+A1*Y1+A2*Y2
      ENDIF
      FINT=GINT
      RETURN
      END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      DOUBLE PRECISION FUNCTION TPTGG(X)
C--GG --> HG: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(3)
      COMMON/MASS/AMH,AMQ,S
      COMMON/PTY/PT,Y,NLOOP
      COMMON/PTCUT/PTMIN0,PTMAX0
      PTMAX = DSQRT((S+AMH**2)**2 - 4*AMH**2*S)/2.D0/DSQRT(S)
      PTMIN = PTMIN0
      PTMAX = DMIN1(PTMAX,PTMAX0)
      PT = PTMIN + (PTMAX-PTMIN)*X(1)
      TPTGG=FPTGG(X(2)) * (PTMAX-PTMIN)
      RETURN
      END

      DOUBLE PRECISION FUNCTION TPTGQ(X)
C--GQ --> HQ: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(2)
      COMMON/MASS/AMH,AMQ,S
      COMMON/PTY/PT,Y,NLOOP
      COMMON/PTCUT/PTMIN0,PTMAX0
      PTMAX = DSQRT((S+AMH**2)**2 - 4*AMH**2*S)/2.D0/DSQRT(S)
      PTMIN = PTMIN0
      PTMAX = DMIN1(PTMAX,PTMAX0)
      PT = PTMIN + (PTMAX-PTMIN)*X(1)
      TPTGQ=FPTGQ(X(2)) * (PTMAX-PTMIN)
      RETURN
      END

      DOUBLE PRECISION FUNCTION TPTQQ(X)
C--QQBAR --> HG: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(2)
      COMMON/MASS/AMH,AMQ,S
      COMMON/PTY/PT,Y,NLOOP
      COMMON/PTCUT/PTMIN0,PTMAX0
      PTMAX = DSQRT((S+AMH**2)**2 - 4*AMH**2*S)/2.D0/DSQRT(S)
      PTMIN = PTMIN0
      PTMAX = DMIN1(PTMAX,PTMAX0)
      PT = PTMIN + (PTMAX-PTMIN)*X(1)
      TPTQQ=FPTQQ(X(2)) * (PTMAX-PTMIN)
      RETURN
      END

      DOUBLE PRECISION FUNCTION FPTGG(X)
C--GG --> HG: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(2)
      COMMON/MASS/AMH,AMQ,S
      COMMON/PTY/PT,Y,NLOOP
      Y0 = (S+AMH**2)/2.D0/DSQRT(S*(AMH**2+PT**2))
      YMAX = DLOG(Y0 + DSQRT(Y0**2 -1.D0))
      Y = -YMAX + 2*YMAX*X(1)
      FPTGG=FPTYGG(X(2)) * 2*YMAX
      RETURN
      END

      DOUBLE PRECISION FUNCTION FPTGQ(X)
C--GQ --> HQ: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(2)
      COMMON/MASS/AMH,AMQ,S
      COMMON/PTY/PT,Y,NLOOP
      Y0 = (S+AMH**2)/2.D0/DSQRT(S*(AMH**2+PT**2))
      YMAX = DLOG(Y0 + DSQRT(Y0**2 -1.D0))
      Y = -YMAX + 2*YMAX*X(1)
      FPTGQ=FPTYGQ(X(2)) * 2*YMAX
      RETURN
      END

      DOUBLE PRECISION FUNCTION FPTQQ(X)
C--QQBAR --> HG: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(2)
      COMMON/MASS/AMH,AMQ,S
      COMMON/PTY/PT,Y,NLOOP
      Y0 = (S+AMH**2)/2.D0/DSQRT(S*(AMH**2+PT**2))
      YMAX = DLOG(Y0 + DSQRT(Y0**2 -1.D0))
      Y = -YMAX + 2*YMAX*X(1)
      FPTQQ=FPTYQQ(X(2)) * 2*YMAX
      RETURN
      END

      DOUBLE PRECISION FUNCTION FPTYGG(X)
C--GG --> HG: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(1)
      DIMENSION PDF(-6:6)
      COMMON/MASS/AMH,AMQ,S
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/ETACUT/ETA0
      COMMON/PTY/PT,Y,NLOOP
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/PARSC/XKAPM,XKAPQ
      PI = 4.D0*DATAN(1.D0)
      AMT2 = AMH**2+PT**2
      QM=XKAPM*DSQRT(AMT2)
      QQ=XKAPQ*DSQRT(AMT2)
C     QM=XKAPM*AMH
C     QQ=XKAPQ*AMH
      XCOF = DSQRT(S*AMT2)
      UU1 = -XCOF*DEXP(Y)
      TT1 = -XCOF*DEXP(-Y)
      X1M = -(AMH**2+UU1)/(S+TT1)
      ONE=1.D0
      X1 = X1M + (ONE-X1M)*X(1)
      X2 = -(X1*TT1+AMH**2)/(X1*S+UU1)
      if(x1.gt.1.d0.or.x2.gt.1.d0)
     .write(6,*)X1,X2,X1M,Y,AMH,DSQRT(S),TT1,UU1
      CALL STRUC(X1,QQ,PDF)
      GG=PDF(0)
      CALL STRUC(X2,QQ,PDF)
      DUM = GG*PDF(0)
      S0 = X1*X2*S
      T0 = X1*TT1 + AMH**2
      U0 = X2*UU1 + AMH**2
      ETAJ = DLOG(X1*T0/X2/U0)/2
c     ETAJ0 = DLOG(X1*(X1*TT1+AMH**2)/X2/(X2*UU1+AMH**2))/2
c     write(6,*)'eta: ',ETAJ,ETAJ0,ETAJ0/ETAJ,ETA0
      FPTYGG = 0
c     IF(DABS(ETAJ).LE.ETA0)THEN
       FPTYGG=DPTYGG(S0,T0)*DUM*S/(X1*S+UU1)*(ONE-X1M)
     .       * ALPHAS(QM,NLOOP)**3/PI
     .       * 2.D0*PT
c     ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION FPTYGQ(X)
C--GQ --> HQ: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(1)
      DIMENSION PDF(-6:6)
      COMMON/MASS/AMH,AMQ,S
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/ETACUT/ETA0
      COMMON/PTY/PT,Y,NLOOP
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/PARSC/XKAPM,XKAPQ
      PI = 4.D0*DATAN(1.D0)
      AMT2 = AMH**2+PT**2
      QM=XKAPM*DSQRT(AMT2)
      QQ=XKAPQ*DSQRT(AMT2)
C     QM=XKAPM*AMH
C     QQ=XKAPQ*AMH
      XCOF = DSQRT(S*AMT2)
      UU1 = -XCOF*DEXP(Y)
      TT1 = -XCOF*DEXP(-Y)
      X1M = -(AMH**2+UU1)/(S+TT1)
      ONE=1.D0
      X1 = X1M + (ONE-X1M)*X(1)
      X2 = -(X1*TT1+AMH**2)/(X1*S+UU1)
      CALL STRUC(X1,QQ,PDF)
      G1=PDF(0)
      Q1=0.D0
      DO 1 I=1,5
       Q1=Q1+PDF(I)+PDF(-I)
1     CONTINUE
      CALL STRUC(X2,QQ,PDF)
      G2=PDF(0)
      Q2=0.D0
      DO 2 I=1,5
       Q2=Q2+PDF(I)+PDF(-I)
2     CONTINUE
      S0 = X1*X2*S
      T0 = X1*TT1 + AMH**2
      U0 = X2*UU1 + AMH**2
      DUM=G1*Q2*DPTYGQ(S0,T0) + Q1*G2*DPTYGQ(S0,U0)
      DUM=G1*Q2*DPTYGQ(S0,T0) + Q1*G2*DPTYGQ(S0,U0)
      ETAJ = DLOG(X1*T0/X2/U0)/2
      FPTYGQ = 0
c     IF(DABS(ETAJ).LE.ETA0)THEN
       FPTYGQ=DUM*S/(X1*S+UU1)*(ONE-X1M)
     .       * ALPHAS(QM,NLOOP)**3/PI
     .       * 2.D0*PT
c     ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION FPTYQQ(X)
C--QQBAR --> HG: INTEGRAND FOR VEGAS-INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION X(1)
      DIMENSION PDF1(-6:6),PDF2(-6:6)
      COMMON/MASS/AMH,AMQ,S
      COMMON/CUT/EPST,EPSV,REPS
      COMMON/ETACUT/ETA0
      COMMON/PTY/PT,Y,NLOOP
      COMMON/CALLS/ICALL1,ICALL2,IFAIL,IF66
      COMMON/PARSC/XKAPM,XKAPQ
      COMMON/COLLI/ICOLL
      PI = 4.D0*DATAN(1.D0)
      AMT2 = AMH**2+PT**2
      QM=XKAPM*DSQRT(AMT2)
      QQ=XKAPQ*DSQRT(AMT2)
C     QM=XKAPM*AMH
C     QQ=XKAPQ*AMH
      XCOF = DSQRT(S*AMT2)
      UU1 = -XCOF*DEXP(Y)
      TT1 = -XCOF*DEXP(-Y)
      X1M = -(AMH**2+UU1)/(S+TT1)
      ONE=1.D0
      X1 = X1M + (ONE-X1M)*X(1)
      X2 = -(X1*TT1+AMH**2)/(X1*S+UU1)
      CALL STRUC(X1,QQ,PDF1)
      CALL STRUC(X2,QQ,PDF2)
      DUM1=0.D0
      DUM2=0.D0
      DO 1 I=1,5
       DUM1=DUM1+PDF1( I)*PDF2(-ICOLL*I)
       DUM2=DUM2+PDF1(-I)*PDF2( ICOLL*I)
1     CONTINUE
      S0 = X1*X2*S
      T0 = X1*TT1 + AMH**2
      U0 = X2*UU1 + AMH**2
      DUM=DUM1*DPTYQQ(S0,T0) + DUM2*DPTYQQ(S0,U0)
      ETAJ = DLOG(X1*T0/X2/U0)/2
      FPTYQQ = 0
c     IF(DABS(ETAJ).LE.ETA0)THEN
       FPTYQQ=DUM*S/(X1*S+UU1)*(ONE-X1M)
     .       * ALPHAS(QM,NLOOP)**3/PI
     .       * 2.D0*PT
c     ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DPTYQQ(SS,TT)
C--QQBAR --> HG: MASS-DEPENDENT INTEGRAND
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/HIGGS/IHIGGS
      COMMON/MASS/AMH,AMQ,S
      COMMON/SUSY/AMB,AMC,FACT,FACB,FACC,FACCG,FACCTG,ISILH
      COMMON/SUSY0/AMST1,AMST2,AMSB1,AMSB2,FACST1,FACST2,FACSB1,FACSB2
      COMMON/GRENZ/IGRENZ
      COMMON/GRENZSQ/IGRENZSQ
      COMMON/BORN/CCF0,CCF0T,CCF0B,CCF0C,CF0G
      COMMON/BORN00/CF0G00
C--MANDELSTAM-VARIABLES
      UU = AMH**2 - SS - TT
      TAUH = AMH**2/SS
      V = -TT/SS/(1.D0-TAUH)
      CSQQ=DCMPLX(0.D0,0.D0)
      CSQQP=DCMPLX(0.D0,0.D0)
C--FORM FACTOR
      IF(FACT.NE.0.D0)THEN
        IF(IGRENZ.EQ.1)THEN
         CSQQ=CSQQ+FACT*CLGQ(TT,SS,UU,AMQ)*CFBORN(AMH,AMQ)
         CSQQP=CSQQP+(FACT-1)*CLGQ(TT,SS,UU,AMQ)*CFBORN(AMH,AMQ)
        ELSE
         CSQQ=CSQQ+FACT*CIGQ(TT,SS,UU,AMQ)
         CSQQP=CSQQP+(FACT-1)*CIGQ(TT,SS,UU,AMQ)
        ENDIF
        IF(FACST1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          CSQQ=CSQQ+FACST1*CLGQ(TT,SS,UU,AMST1)*CSBORN(AMH,AMST1)
     .             +FACST2*CLGQ(TT,SS,UU,AMST2)*CSBORN(AMH,AMST2)
         ELSE
          CSQQ=CSQQ+FACST1*CIGQS(TT,SS,UU,AMST1)
     .             +FACST2*CIGQS(TT,SS,UU,AMST2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACB.NE.0.D0)THEN
        CSQQ=CSQQ+FACB*CIGQ(TT,SS,UU,AMB)
        CSQQP=CSQQP+(FACB-1)*CIGQ(TT,SS,UU,AMB)
        IF(FACSB1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          CSQQ=CSQQ+FACSB1*CLGQ(TT,SS,UU,AMSB1)*CSBORN(AMH,AMSB1)
     .             +FACSB2*CLGQ(TT,SS,UU,AMSB2)*CSBORN(AMH,AMSB2)
         ELSE
          CSQQ=CSQQ+FACSB1*CIGQS(TT,SS,UU,AMSB1)
     .             +FACSB2*CIGQS(TT,SS,UU,AMSB2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACC.NE.0.D0)THEN
        CSQQ=CSQQ+FACC*CIGQ(TT,SS,UU,AMC)
        CSQQP=CSQQP+(FACC-1)*CIGQ(TT,SS,UU,AMC)
      ENDIF
      IF(FACCG.NE.0.D0)THEN
        CSQQ=CSQQ-4*CF0G00/3
        CSQQP=CSQQP-4*CF0G00/3
      ENDIF
      IF(FACCTG.NE.0.D0)THEN
        CSQQ=CSQQ+FACCTG*CIGQT(TT,SS,UU,AMQ)
        CSQQP=CSQQP+FACCTG*CIGQT(TT,SS,UU,AMQ)
      ENDIF
C--CURRENT FACTOR
      HQQ=(UU**2+TT**2)/SS
C--LOWEST ORDER FORM FACTOR
      IF(IHIGGS.NE.1)THEN
C--SCALAR HIGGS
        CF0=FACT*CFBORN(AMH,AMQ)+FACB*CFBORN(AMH,AMB)
     .     +FACC*CFBORN(AMH,AMC)+CF0G
     .     +FACST1*CSBORN(AMH,AMST1)+FACST2*CSBORN(AMH,AMST2)
     .     +FACSB1*CSBORN(AMH,AMSB1)+FACSB2*CSBORN(AMH,AMSB2)
        CF0P=(FACT-1)*CFBORN(AMH,AMQ)+(FACB-1)*CFBORN(AMH,AMB)
     .      +(FACC-1)*CFBORN(AMH,AMC)+CF0G
      ELSE
C--PSEUDOSCALAR HIGGS
        CSQQP=0
        CF0P=0
        CF0=FACT*CFBORNA(AMH,AMQ)+FACB*CFBORNA(AMH,AMB)
     .     +FACC*CFBORNA(AMH,AMC)
      ENDIF
      DPTYQQ=CDABS(CSQQ)**2*HQQ/CDABS(CF0)**2
     .      *(1.D0-TAUH)/SS
     .      /(1.D0-TAUH)/SS
      IF(ISILH.NE.0) DPTYQQ=(CDABS(CSQQ)**2-CDABS(CSQQP)**2)*HQQ
     .                   / (CDABS(CF0)**2-CDABS(CF0P)**2)
     .      *(1.D0-TAUH)/SS
     .      /(1.D0-TAUH)/SS
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DPTYGQ(SS,TT)
C--GQ --> HQ: MASS-DEPENDENT INTEGRAND
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/HIGGS/IHIGGS
      COMMON/MASS/AMH,AMQ,S
      COMMON/SUSY/AMB,AMC,FACT,FACB,FACC,FACCG,FACCTG,ISILH
      COMMON/SUSY0/AMST1,AMST2,AMSB1,AMSB2,FACST1,FACST2,FACSB1,FACSB2
      COMMON/GRENZ/IGRENZ
      COMMON/GRENZSQ/IGRENZSQ
      COMMON/BORN/CCF0,CCF0T,CCF0B,CCF0C,CF0G
      COMMON/BORN00/CF0G00
C--MANDELSTAM-VARIABLES
      UU = AMH**2 - SS - TT
      TAUH = AMH**2/SS
      V = -TT/SS/(1.D0-TAUH)
      CSGQ=DCMPLX(0.D0,0.D0)
      CSGQP=DCMPLX(0.D0,0.D0)
C--FORM FACTOR
      IF(FACT.NE.0.D0)THEN
        IF(IGRENZ.EQ.1)THEN
         CSGQ=CSGQ+FACT*CLGQ(SS,TT,UU,AMQ)*CFBORN(AMH,AMQ)
         CSGQP=CSGQP+(FACT-1)*CLGQ(SS,TT,UU,AMQ)*CFBORN(AMH,AMQ)
        ELSE
         CSGQ=CSGQ+FACT*CIGQ(SS,TT,UU,AMQ)
         CSGQP=CSGQP+(FACT-1)*CIGQ(SS,TT,UU,AMQ)
        ENDIF
        IF(FACST1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          CSGQ=CSGQ+FACST1*CLGQ(SS,TT,UU,AMST1)*CSBORN(AMH,AMST1)
     .             +FACST2*CLGQ(SS,TT,UU,AMST2)*CSBORN(AMH,AMST2)
         ELSE
          CSGQ=CSGQ+FACST1*CIGQS(SS,TT,UU,AMST1)
     .             +FACST2*CIGQS(SS,TT,UU,AMST2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACB.NE.0.D0)THEN
        CSGQ=CSGQ+FACB*CIGQ(SS,TT,UU,AMB)
        CSGQP=CSGQP+(FACB-1)*CIGQ(SS,TT,UU,AMB)
        IF(FACSB1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          CSGQ=CSGQ+FACSB1*CLGQ(SS,TT,UU,AMSB1)*CSBORN(AMH,AMSB1)
     .             +FACSB2*CLGQ(SS,TT,UU,AMSB2)*CSBORN(AMH,AMSB2)
         ELSE
          CSGQ=CSGQ+FACSB1*CIGQS(SS,TT,UU,AMSB1)
     .             +FACSB2*CIGQS(SS,TT,UU,AMSB2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACC.NE.0.D0)THEN
        CSGQ=CSGQ+FACC*CIGQ(SS,TT,UU,AMC)
        CSGQP=CSGQP+(FACC-1)*CIGQ(SS,TT,UU,AMC)
      ENDIF
      IF(FACCG.NE.0.D0)THEN
        CSGQ=CSGQ-4*CF0G00/3
        CSGQP=CSGQP-4*CF0G00/3
      ENDIF
      IF(FACCTG.NE.0.D0)THEN
        CSGQ=CSGQ+FACCTG*CIGQT(SS,TT,UU,AMQ)
        CSGQP=CSGQP+FACCTG*CIGQT(SS,TT,UU,AMQ)
      ENDIF
C--CURRENT FACTOR
      HGQ=-(3.D0*(UU**2+SS**2))/(8.D0*TT)
C--LOWEST ORDER FORM FACTOR
      IF(IHIGGS.NE.1)THEN
C--SCALAR HIGGS
        CF0=FACT*CFBORN(AMH,AMQ)+FACB*CFBORN(AMH,AMB)
     .     +FACC*CFBORN(AMH,AMC)+CF0G
     .     +FACST1*CSBORN(AMH,AMST1)+FACST2*CSBORN(AMH,AMST2)
     .     +FACSB1*CSBORN(AMH,AMSB1)+FACSB2*CSBORN(AMH,AMSB2)
        CF0P=(FACT-1)*CFBORN(AMH,AMQ)+(FACB-1)*CFBORN(AMH,AMB)
     .      +(FACC-1)*CFBORN(AMH,AMC)+CF0G
      ELSE
C--PSEUDOSCALAR HIGGS
        CSGQP=0
        CF0P=0
        CF0=FACT*CFBORNA(AMH,AMQ)+FACB*CFBORNA(AMH,AMB)
     .     +FACC*CFBORNA(AMH,AMC)
      ENDIF
      DPTYGQ=(CDABS(CSGQ)**2*HGQ/CDABS(CF0)**2
     .       *(1.D0-TAUH)/SS)
     .       /(1.D0-TAUH)/SS
      IF(ISILH.NE.0) DPTYGQ=((CDABS(CSGQ)**2-CDABS(CSGQP)**2)*HGQ
     .                     / (CDABS(CF0)**2-CDABS(CF0P)**2)
     .                     *(1.D0-TAUH)/SS)
     .                     /(1.D0-TAUH)/SS
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DPTYGG(SS,TT)
C--GG --> HG: MASS-DEPENDENT INTEGRAND
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      COMMON/HIGGS/IHIGGS
      COMMON/MASS/AMH,AMQ,S
      COMMON/SUSY/AMB,AMC,FACT,FACB,FACC,FACCG,FACCTG,ISILH
      COMMON/SUSY0/AMST1,AMST2,AMSB1,AMSB2,FACST1,FACST2,FACSB1,FACSB2
      COMMON/GRENZ/IGRENZ
      COMMON/GRENZSQ/IGRENZSQ
      COMMON/BORN/CCF0,CCF0T,CCF0B,CCF0C,CF0G
      COMMON/BORN00/CF0G00
C--MANDELSTAM-VARIABLES
      UU = AMH**2 - SS - TT
      TAUH = AMH**2/SS
      V = -TT/SS/(1.D0-TAUH)
      C5=DCMPLX(0.D0,0.D0)
      C2=DCMPLX(0.D0,0.D0)
      C3=DCMPLX(0.D0,0.D0)
      C4=DCMPLX(0.D0,0.D0)
      C5P=DCMPLX(0.D0,0.D0)
      C2P=DCMPLX(0.D0,0.D0)
      C3P=DCMPLX(0.D0,0.D0)
      C4P=DCMPLX(0.D0,0.D0)
C--FORM FACTOR
      IF(FACT.NE.0.D0)THEN
        IF(IGRENZ.EQ.1)THEN
         C5=C5+FACT*CFGL5(SS,TT,UU,AMQ)*CFBORN(AMH,AMQ)/2
     .        +FACT*CFGL5(SS,TT,UU,AMQ)*CFBORN(AMH,AMQ)/2
         C2=C2+FACT*CFGL2(SS,TT,UU,AMQ)*CFBORN(AMH,AMQ)/2
     .        +FACT*CFGL2(SS,TT,UU,AMQ)*CFBORN(AMH,AMQ)/2
         C3=C3+FACT*CFGL2(TT,SS,UU,AMQ)*CFBORN(AMH,AMQ)/2
     .        +FACT*CFGL2(TT,SS,UU,AMQ)*CFBORN(AMH,AMQ)/2
         C4=C4+FACT*CFGL2(UU,TT,SS,AMQ)*CFBORN(AMH,AMQ)/2
     .        +FACT*CFGL2(UU,TT,SS,AMQ)*CFBORN(AMH,AMQ)/2
         C5P=C5P+(FACT-1)*CFGL5(SS,TT,UU,AMQ)*CFBORN(AMH,AMQ)/2
     .          +(FACT-1)*CFGL5(SS,TT,UU,AMQ)*CFBORN(AMH,AMQ)/2
         C2P=C2P+(FACT-1)*CFGL2(SS,TT,UU,AMQ)*CFBORN(AMH,AMQ)/2
     .          +(FACT-1)*CFGL2(SS,TT,UU,AMQ)*CFBORN(AMH,AMQ)/2
         C3P=C3P+(FACT-1)*CFGL2(TT,SS,UU,AMQ)*CFBORN(AMH,AMQ)/2
     .          +(FACT-1)*CFGL2(TT,SS,UU,AMQ)*CFBORN(AMH,AMQ)/2
         C4P=C4P+(FACT-1)*CFGL2(UU,TT,SS,AMQ)*CFBORN(AMH,AMQ)/2
     .          +(FACT-1)*CFGL2(UU,TT,SS,AMQ)*CFBORN(AMH,AMQ)/2
        ELSE
         C5=C5+FACT*CFGG5(SS,TT,UU,AMQ)
         C2=C2+FACT*CFGG2(SS,TT,UU,AMQ)
         C3=C3+FACT*CFGG2(TT,SS,UU,AMQ)
         C4=C4+FACT*CFGG2(UU,TT,SS,AMQ)
         C5P=C5P+(FACT-1)*CFGG5(SS,TT,UU,AMQ)
         C2P=C2P+(FACT-1)*CFGG2(SS,TT,UU,AMQ)
         C3P=C3P+(FACT-1)*CFGG2(TT,SS,UU,AMQ)
         C4P=C4P+(FACT-1)*CFGG2(UU,TT,SS,AMQ)
        ENDIF
        IF(FACST1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          C5=C5+FACST1*CFGL5(SS,TT,UU,AMST1)*CSBORN(AMH,AMST1)
     .         +FACST2*CFGL5(SS,TT,UU,AMST2)*CSBORN(AMH,AMST2)
          C2=C2+FACST1*CFGL2(SS,TT,UU,AMST1)*CSBORN(AMH,AMST1)
     .         +FACST2*CFGL2(SS,TT,UU,AMST2)*CSBORN(AMH,AMST2)
          C3=C3+FACST1*CFGL2(TT,SS,UU,AMST1)*CSBORN(AMH,AMST1)
     .         +FACST2*CFGL2(TT,SS,UU,AMST2)*CSBORN(AMH,AMST2)
          C4=C4+FACST1*CFGL2(UU,TT,SS,AMST1)*CSBORN(AMH,AMST1)
     .         +FACST2*CFGL2(UU,TT,SS,AMST2)*CSBORN(AMH,AMST2)
         ELSE
          C5=C5+FACST1*CFGG5S(SS,TT,UU,AMST1)
     .         +FACST2*CFGG5S(SS,TT,UU,AMST2)
          C2=C2+FACST1*CFGG2S(SS,TT,UU,AMST1)
     .         +FACST2*CFGG2S(SS,TT,UU,AMST2)
          C3=C3+FACST1*CFGG2S(TT,SS,UU,AMST1)
     .         +FACST2*CFGG2S(TT,SS,UU,AMST2)
          C4=C4+FACST1*CFGG2S(UU,TT,SS,AMST1)
     .         +FACST2*CFGG2S(UU,TT,SS,AMST2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACB.NE.0.D0)THEN
        C5=C5+FACB*CFGG5(SS,TT,UU,AMB)
        C2=C2+FACB*CFGG2(SS,TT,UU,AMB)
        C3=C3+FACB*CFGG2(TT,SS,UU,AMB)
        C4=C4+FACB*CFGG2(UU,TT,SS,AMB)
        C5P=C5P+(FACB-1)*CFGG5(SS,TT,UU,AMB)
        C2P=C2P+(FACB-1)*CFGG2(SS,TT,UU,AMB)
        C3P=C3P+(FACB-1)*CFGG2(TT,SS,UU,AMB)
        C4P=C4P+(FACB-1)*CFGG2(UU,TT,SS,AMB)
        IF(FACSB1.NE.0.D0)THEN
         IF(IGRENZSQ.EQ.1)THEN
          C5=C5+FACSB1*CFGL5(SS,TT,UU,AMSB1)*CSBORN(AMH,AMSB1)
     .         +FACSB2*CFGL5(SS,TT,UU,AMSB2)*CSBORN(AMH,AMSB2)
          C2=C2+FACSB1*CFGL2(SS,TT,UU,AMSB1)*CSBORN(AMH,AMSB1)
     .         +FACSB2*CFGL2(SS,TT,UU,AMSB2)*CSBORN(AMH,AMSB2)
          C3=C3+FACSB1*CFGL2(TT,SS,UU,AMSB1)*CSBORN(AMH,AMSB1)
     .         +FACSB2*CFGL2(TT,SS,UU,AMSB2)*CSBORN(AMH,AMSB2)
          C4=C4+FACSB1*CFGL2(UU,TT,SS,AMSB1)*CSBORN(AMH,AMSB1)
     .         +FACSB2*CFGL2(UU,TT,SS,AMSB2)*CSBORN(AMH,AMSB2)
         ELSE
          C5=C5+FACSB1*CFGG5S(SS,TT,UU,AMSB1)
     .         +FACSB2*CFGG5S(SS,TT,UU,AMSB2)
          C2=C2+FACSB1*CFGG2S(SS,TT,UU,AMSB1)
     .         +FACSB2*CFGG2S(SS,TT,UU,AMSB2)
          C3=C3+FACSB1*CFGG2S(TT,SS,UU,AMSB1)
     .         +FACSB2*CFGG2S(TT,SS,UU,AMSB2)
          C4=C4+FACSB1*CFGG2S(UU,TT,SS,AMSB1)
     .         +FACSB2*CFGG2S(UU,TT,SS,AMSB2)
         ENDIF
        ENDIF
      ENDIF
      IF(FACC.NE.0.D0)THEN
        C5=C5+FACC*CFGG5(SS,TT,UU,AMC)
        C2=C2+FACC*CFGG2(SS,TT,UU,AMC)
        C3=C3+FACC*CFGG2(TT,SS,UU,AMC)
        C4=C4+FACC*CFGG2(UU,TT,SS,AMC)
        C5P=C5P+(FACC-1)*CFGG5(SS,TT,UU,AMC)
        C2P=C2P+(FACC-1)*CFGG2(SS,TT,UU,AMC)
        C3P=C3P+(FACC-1)*CFGG2(TT,SS,UU,AMC)
        C4P=C4P+(FACC-1)*CFGG2(UU,TT,SS,AMC)
      ENDIF
      IF(FACCG.NE.0.D0)THEN
        C5=C5+CF0G00*(SS+TT+UU)**2
        C2=C2+CF0G00*SS**2
        C3=C3+CF0G00*TT**2
        C4=C4+CF0G00*UU**2
        C5P=C5P+CF0G00*(SS+TT+UU)**2
        C2P=C2P+CF0G00*SS**2
        C3P=C3P+CF0G00*TT**2
        C4P=C4P+CF0G00*UU**2
c       write(6,*)CDABS(CF0G00)/12/FACCG
      ENDIF
      IF(FACCTG.NE.0.D0)THEN
        C5=C5+FACCTG*CFGG5T(SS,TT,UU,AMQ)
        C2=C2+FACCTG*CFGG2T(SS,TT,UU,AMQ)
        C3=C3+FACCTG*CFGG2T(TT,SS,UU,AMQ)
        C4=C4+FACCTG*CFGG2T(UU,TT,SS,AMQ)
        C5P=C5P+FACCTG*CFGG5T(SS,TT,UU,AMQ)
        C2P=C2P+FACCTG*CFGG2T(SS,TT,UU,AMQ)
        C3P=C3P+FACCTG*CFGG2T(TT,SS,UU,AMQ)
        C4P=C4P+FACCTG*CFGG2T(UU,TT,SS,AMQ)
      ENDIF
      FGG=CDABS(C5)**2+CDABS(C2)**2+CDABS(C3)**2+CDABS(C4)**2
      FGGP=CDABS(C5P)**2+CDABS(C2P)**2+CDABS(C3P)**2+CDABS(C4P)**2
      IF(ISILH.NE.0)FGG = FGG - FGGP
      FGG=FGG/SS**4
C--LOWEST ORDER FORM FACTOR
      IF(IHIGGS.NE.1)THEN
C--SCALAR HIGGS
        CF0=FACT*CFBORN(AMH,AMQ)+FACB*CFBORN(AMH,AMB)
     .     +FACC*CFBORN(AMH,AMC)+CF0G
     .     +FACST1*CSBORN(AMH,AMST1)+FACST2*CSBORN(AMH,AMST2)
     .     +FACSB1*CSBORN(AMH,AMSB1)+FACSB2*CSBORN(AMH,AMSB2)
        CF0P=(FACT-1)*CFBORN(AMH,AMQ)+(FACB-1)*CFBORN(AMH,AMB)
     .      +(FACC-1)*CFBORN(AMH,AMC)+CF0G
      ELSE
C--PSEUDOSCALAR HIGGS
        FGGP=0
        CF0P=0
        CF0=FACT*CFBORNA(AMH,AMQ)+FACB*CFBORNA(AMH,AMB)
     .     +FACC*CFBORNA(AMH,AMC)
      ENDIF
      DPTYGG=3.D0*FGG/CDABS(CF0)**2
     .       /(1.D0-TAUH)/V/(1.D0-V)/2.D0
     .       /(1.D0-TAUH)/SS
      IF(ISILH.NE.0)
     .       DPTYGG=3.D0*(FGG-FGGP)/(CDABS(CF0)**2-CDABS(CF0P)**2)
     .             /(1.D0-TAUH)/V/(1.D0-V)/2.D0
     .             /(1.D0-TAUH)/SS
      RETURN
      END
 
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine hdecini(amx)
      implicit double precision (a-h,o-z)
c--higlu common blocks
      common/hmass/amsm,ama,amhl,amhh,amch,amar
      common/susy/amb1,amc1,fact,facb,facc,faccg,facctg,isilh
      common/susy0/amst1,amst2,amsb1,amsb2,facst1,facst2,facsb1,facsb2
      common/als/xlambda,amac,amab,amat,n0
      common/masses/ams,amc0,amb0,amt0
      common/strange/amsb
      common/susyp/gf,alpha,sw2,tgbet,amtq,amz,amsq
      common/break/amsq0,amur,amdr,au,ad,amu,am2
      common/breakglu/amglu
      common/sm4/amt4,amb4,fact4,facb4,ism4,iggelw
      common/param0/gf0,amw0,amz0
      common/coup/gat,gab,glt,glb,ght,ghb,gzah,gzal,
     .            ghhh,glll,ghll,glhh,ghaa,glaa,glvv,ghvv,
     .            glpm,ghpm,b,a
c--hdecay common blocks
      common/masses_hdec/ams_x,amc_x,amb_x,amt_x
      common/strange_hdec/amsb_x
      common/param_hdec/gf_x,alph_x,amtau_x,ammuon_x,amz_x,amw_x
      common/hmass_hdec/amsm_x,ama_x,aml_x,amh_x,amch_x,amar_x
      common/breakglu_hdec/amglu_x
      common/break_hdec/amel_x,amer_x,amsq_x,amur_x,amdr_x,
     .                  al_x,au_x,ad_x,amu_x,am2_x
      common/coup_hdec/gat_x,gab_x,glt_x,glb_x,ght_x,ghb_x,gzah_x,gzal_x
     .        ,ghhh_x,glll_x,ghll_x,glhh_x,ghaa_x,glaa_x,glvv_x,ghvv_x,
     .         glpm_x,ghpm_x,b_x,a_x
      common/sm4_hdec/amtp_x,ambp_x,amnup_x,amep_x,ism4_x,iggelw_x
      common/thdm_hdec/tgbet2hdm_x,alph2hdm_x,amhl2hdm_x,amhh2hdm_x,
     .     amha2hdm_x,amhc2hdm_x,am12sq_x,a1lam2hdm_x,a2lam2hdm_x,
     .     a3lam2hdm_x,a4lam2hdm_x,a5lam2hdm_x,itype2hdm_x,
     .     i2hdm_x,iparam2hdm_x
      common/model_hdec/imodel_x
      common/model/imodel,isusy,i2hdm
      common/flag_hdec/ihiggs_x,nnlo_x,ipole_x
      COMMON/SLHA_vals_HDEC/islhai,islhao
c----------------------------------------------------------------------
c--read hdecay input
      call read_hdec(tgbet,amabeg,amaend,nma)
      amar_x  = amx
      if(islhai.ne.0) amar_x = amabeg
      amsm_x  = amar_x
      ama_x   = amar_x
      ism4    = ism4_x
      iggelw  = iggelw_x
      i2hdm   = i2hdm_x
      imodel  = 0
      isusy   = 0
      if(ihiggs_x.ne.0.and.i2hdm.eq.0)then
       imodel = imodel_x
       isusy  = 1
      endif
c--call hdec to initialize everything
      call head_hdec(tgbet,ama_x)
      call hdec(tgbet)
      call write_hdec(tgbet)
      call close_hdec
      amsm    = amsm_x
      ama     = ama_x
      amhl    = aml_x
      amhh    = amh_x
      amch    = amch_x
      amb1    = amb_x
      amc1    = amc_x
      amst1   = amst1_x
      amst2   = amst2_x
      amsb1   = amsb1_x
      amsb2   = amsb2_x
      ams     = ams_x
      amc0    = amc_x
      amb0    = amb_x
      amt0    = amt_x
      amsb    = amsb_x
c     xlambda = xlambda_x
      amac    = amc_x
      amab    = amb_x
      amat    = amt_x
      n0      = n0_x
      gf      = gf_x
      alpha   = alpha_x
      sw2     = 1-amw_x**2/amz_x**2
      amtq    = amt_x
      amz     = amz_x
      amsq    = amsq_x
      amsq0   = amsq_x
      amur    = amur_x
      amdr    = amdr_x
      au      = au_x
      ad      = ad_x
      amu     = amu_x
      am2     = am2_x
      amglu   = amglu_x
      amt4    = amtp_x
      amb4    = ambp_x
      amw0    = amw_x
      amz0    = amz_x
      gat     = gat_x
      gab     = gab_x
      glt     = glt_x
      glb     = glb_x
      ght     = ght_x
      ghb     = ghb_x
      gzah    = gzah_x
      gzal    = gzal_x
      ghhh    = ghhh_x
      glll    = glll_x
      ghll    = ghll_x
      glhh    = glhh_x
      ghaa    = ghaa_x
      glaa    = glaa_x
      glvv    = glvv_x
      ghvv    = ghvv_x
      glpm    = glpm_x
      ghpm    = ghpm_x
      b       = b_x
      a       = a_x
      return
      end

      subroutine hdectrafo1
      implicit double precision (a-h,o-z)
c--higlu common blocks
      common/hmass/amsm,ama,amhl,amhh,amch,amar
c--hdecay common blocks
      common/hmass_hdec/amsm_x,ama_x,aml_x,amh_x,amch_x,amar_x
c----------------------------------------------------------------------
      amsm_x  = amsm
      ama_x   = ama
      amhl_x  = aml
      amhh_x  = amh
      amch_x  = amch
      return
      end

      subroutine hdectrafo2
      implicit double precision (a-h,o-z)
c--higlu common blocks
      common/hmass/amsm,ama,amhl,amhh,amch,amar
c--hdecay common blocks
      common/hmass_hdec/amsm_x,ama_x,aml_x,amh_x,amch_x,amar_x
c----------------------------------------------------------------------
      amsm    = amsm_x
      ama     = ama_x
      amhl    = aml_x
      amhh    = amh_x
      amch    = amch_x
      return
      end

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE PDGRV(ISET,X,Q2,XPDF)
 
C...ISET = 1 - LO,          Lambda_4=0.20 GeV, N_f=5
C...       2 - NLO, MS_bar, Lambda_4=0.20 GeV, N_f=5
C...X          - Bjorken x
C...Q2         - square of the momentum scale  (in GeV**2)
C...XPDF(-6:6) - matrix containing  x*p(x,Q2)
C...     IPDF = -6 ,  -5 ,  -4 ,  -3 ,  -2 ,  -1 ,0 ,1,2,3,4,5,6
C...          t_bar,b_bar,c_bar,s_bar,u_bar,d_bar,gl,d,u,s,c,b,t
C...range of validity:
C...     D-5  < X  < 1
C...      0.3 < Q2 < D8  GeV^2
C...REAL*8 version
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XPDF(-6:6)
 
C...LO PARAMETRIZATIONS :
      IF (ISET.EQ.1) THEN
        AMU2  = 0.25
        ALAM2 = 0.232 * 0.232
        S  = LOG (LOG(Q2/ALAM2)/LOG(AMU2/ALAM2))
        S2 = S * S
        S3 = S2 * S
 
C...X * (UV + DV) :
        DNUD = 0.663 + 0.191 * S - 0.041 * S2 + 0.031 * S3
        AKUD = 0.326
        AGUD = -1.97 +  6.74 * S -  1.96 * S2
        BUD  =  24.4 -  20.7 * S +  4.08 * S2
        DUD  =  2.86 +  0.70 * S -  0.02 * S2
        UDV  = PDFV (X, DNUD, AKUD, AGUD, BUD, DUD)
 
C...X * DV :
        DND  = 0.579 + 0.283 * S + 0.047 * S2
        AKD = 0.523 - 0.015 * S
        AGD =  2.22 -  0.59 * S -  0.27 * S2
        BD  =  5.95 -  6.19 * S +  1.55 * S2
        DD  =  3.57 +  0.94 * S -  0.16 * S2
        DV  = PDFV (X, DND, AKD, AGD, BD, DD)
 
C...X * G :
        ALG =  0.558
        BEG =  1.218
        AKG =   1.00 -  0.17 * S
        BKG =   0.0
        AGG =   0.0  + 4.879 * S - 1.383 * S2
        BGG =  25.92 - 28.97 * S + 5.596 * S2
        CG  = -25.69 + 23.68 * S - 1.975 * S2
        DG  =  2.537 + 1.718 * S + 0.353 * S2
        EG  =  0.595 + 2.138 * S
        ESG =  4.066
        GL =PDFW (X, S, ALG, BEG, AKG, BKG, AGG, BGG, CG, DG, EG, ESG)
 
C...X * UBAR = X * DBAR :
        ALU =  1.396
        BEU =  1.331
        AKU =  0.412 - 0.171 * S
        BKU =  0.566 - 0.496 * S
        AGU =  0.363
        BGU = -1.196
        CU  =  1.029 + 1.785 * S - 0.459 * S2
        DU  =  4.696 + 2.109 * S
        EU  =  3.838 + 1.944 * S
        ESU =  2.845
        UDB=PDFW (X, S, ALU, BEU, AKU, BKU, AGU, BGU, CU, DU, EU, ESU)
 
C...X * SBAR = X * S :
        SS  =   0.0
        ALS =  0.803
        BES =  0.563
        AKS =  2.082 - 0.577 * S
        AGS = -3.055 + 1.024 * S **  0.67
        BS  =   27.4 -  20.0 * S ** 0.154
        DS  =   6.22
        EST =   4.33 + 1.408 * S
        ESS =   8.27 - 0.437 * S
        SB =PDFWS (X, S, SS, ALS, BES, AKS, AGS, BS, DS, EST, ESS)
 
C...X * CBAR = X * C :
        SC  =  0.888
        ALC =   1.01
        BEC =   0.37
        AKC =   0.0
        AGC =   0.0
        BC  =   4.24 - 0.804 * S
        DC  =   3.46 + 1.076 * S
        EC  =   4.61 + 1.490 * S
        ESC =  2.555 + 1.961 * S
        CB =PDFWS (X, S, SC, ALC, BEC, AKC, AGC, BC, DC, EC, ESC)
 
C...X * BBAR = X * B :
        SBO =  1.351
        ALB =   1.00
        BEB =   0.51
        AKB =   0.0
        AGB =   0.0
        BBO =  1.848
        DB  =  2.929 + 1.396 * S
        EB  =   4.71 + 1.514 * S
        ESB =   4.02 + 1.239 * S
        BB =PDFWS (X, S, SBO, ALB, BEB, AKB, AGB, BBO, DB, EB, ESB)
 
C...HO parametrization:
      ELSEIF(ISET.EQ.2) THEN
        AMU2  = 0.3
        ALAM2 = 0.248 * 0.248
        S  = LOG (LOG(Q2/ALAM2)/LOG(AMU2/ALAM2))
        DS = SQRT (S)
        S2 = S * S
        S3 = S2 * S
 
C...X * (UV + DV) :
        DNUD  = 0.330 + 0.151 * S - 0.059 * S2 + 0.027 * S3
        AKUD = 0.285
        AGUD = -2.28 + 15.73 * S -  4.58 * S2
        BUD  =  56.7 -  53.6 * S + 11.21 * S2
        DUD  =  3.17 +  1.17 * S -  0.47 * S2 +  0.09 * S3
        UDV  = PDFV (X, DNUD, AKUD, AGUD, BUD, DUD)
 
C...X * DV :
        DND  = 0.459 + 0.315 * DS + 0.515 * S
        AKD = 0.624              - 0.031 * S
        AGD =  8.13 -  6.77 * DS +  0.46 * S
        BD  =  6.59 - 12.83 * DS +  5.65 * S
        DD  =  3.98              +  1.04 * S  -  0.34 * S2
        DV  = PDFV (X, DND, AKD, AGD, BD, DD)
 
C...X * G :
        ALG =  1.128
        BEG =  1.575
        AKG =  0.323 + 1.653 * S
        BKG =  0.811 + 2.044 * S
        AGG =   0.0  + 1.963 * S - 0.519 * S2
        BGG =  0.078 +  6.24 * S
        CG  =  30.77 - 24.19 * S
        DG  =  3.188 + 0.720 * S
        EG  = -0.881 + 2.687 * S
        ESG =  2.466
        GL =PDFW (X, S, ALG, BEG, AKG, BKG, AGG, BGG, CG, DG, EG, ESG)
 
C...X * UBAR = X * DBAR :
        ALU =  0.594
        BEU =  0.614
        AKU =  0.636 - 0.084 * S
        BKU =   0.0
        AGU =  1.121 - 0.193 * S
        BGU =  0.751 - 0.785 * S
        CU  =   8.57 - 1.763 * S
        DU  =  10.22 + 0.668 * S
        EU  =  3.784 + 1.280 * S
        ESU =  1.808 + 0.980 * S
        UDB=PDFW (X, S, ALU, BEU, AKU, BKU, AGU, BGU, CU, DU, EU, ESU)
 
C...X * SBAR = X * S :
        SS  =   0.0
        ALS =  0.756
        BES =  0.101
        AKS =  2.942 - 1.016 * S
        AGS =  -4.60 + 1.167 * S
        BS  =   9.31 - 1.324 * S
        DS  =  11.49 - 1.198 * S + 0.053 * S2
        EST =  2.630 + 1.729 * S
        ESS =   8.12
        SB =PDFWS (X, S, SS, ALS, BES, AKS, AGS, BS, DS, EST, ESS)
 
C...X * CBAR = X * C :
        SC  =  0.820
        ALC =   0.98
        BEC =   0.0
        AKC = -0.625 - 0.523 * S
        AGC =   0.0
        BC  =  1.896 + 1.616 * S
        DC  =   4.12 + 0.683 * S
        EC  =   4.36 + 1.328 * S
        ESC =  0.677 + 0.679 * S
        CB =PDFWS (X, S, SC, ALC, BEC, AKC, AGC, BC, DC, EC, ESC)
 
C...X * BBAR = X * B :
        SBO =  1.297
        ALB =   0.99
        BEB =   0.0
        AKB =   0.0  - 0.193 * S
        AGB =   0.0
        BBO =   0.0
        DB  =  3.447 + 0.927 * S
        EB  =   4.68 + 1.259 * S
        ESB =  1.892 + 2.199 * S
        BB =PDFWS (X, S, SBO, ALB, BEB, AKB, AGB, BBO, DB, EB, ESB)
      ELSE
       WRITE(*,*) ' error in PDGRV: wrong ISET value'
      ENDIF
 
C...final results
      XPDF(0)=GL
      XPDF(1)=DV+UDB
      XPDF(2)=UDV-DV+UDB
      XPDF(3)=SB
      XPDF(4)=CB
      XPDF(5)=BB
      XPDF(6)=0.
      XPDF(-1)=UDB
      XPDF(-2)=UDB
      XPDF(-3)=SB
      XPDF(-4)=CB
      XPDF(-5)=BB
      XPDF(-6)=0.
 
      RETURN
      END
C-----------------------------------------------------
      FUNCTION PDFV (X, DN, AK, AG, B, D)
 
C...functional forms for ho and lo parametrizations :
      IMPLICIT REAL*8 (A-H,O-Z)
       DX = SQRT (X)
       PDFV = DN * X**AK * (1.+ AG*DX + B*X) * (1.- X)**D
      RETURN
      END
C
      FUNCTION PDFW (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
      IMPLICIT REAL*8 (A-H,O-Z)
       ALX = LOG (1./X)
       PDFW = (X**AK * (AG + X * (BG + X*C)) *ALX**BK + S**AL
     1      * EXP (-E + SQRT (ES * S**BE *ALX))) * (1.- X)**D
      RETURN
      END
C-----------------------------------------------------
      FUNCTION PDFWS (X, S, ST, AL, BE, AK, AG, B, D, E, ES)
      IMPLICIT REAL*8 (A-H,O-Z)
       DX = SQRT (X)
       ALX = LOG (1./X)
       IF (S .LE. ST) THEN
         FWS = 0.0
       ELSE
         FWS = (S-ST)**AL / ALX**AK * (1.+ AG*DX + B*X) * (1.- X)**D
     1          * EXP (-E + SQRT (ES * S**BE *ALX))
       ENDIF
       PDFWS=FWS
      RETURN
      END
      
c     SUBROUTINE PDFSET(PARM,VALUE)
c     IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
c     DIMENSION VALUE(20)
c     CHARACTER*20 PARM(20)
c     CONTINUE
c     RETURN
c     END
 
c     SUBROUTINE PFTOPDG(X,Q,PDF)
c     IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
c     DIMENSION PDF(-6:6)
c     CONTINUE
c     RETURN
c     END

c     SUBROUTINE STRUC(X,Q,PDF)
c     IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
c     DIMENSION PDF(-6:6), VALUE(20)
c     DIMENSION PDF1(-6:6)
c     CHARACTER*20 PARM(20)
c     COMMON/STFU/ISET
c     COMMON/PDFLIB/IPDFLIB,NGROUP,NSET,ISETERR
C...X          - BJORKEN X
C...Q          - MOMENTUM SCALE  (IN GeV)
C...PDF(-6:6)  - MATRIX CONTAINING  X*P(X,Q)
C...    IPDF = -6 ,  -5 ,  -4 ,  -3 ,  -2 ,  -1 ,0 ,1,2,3,4,5,6
C...         T_BAR,B_BAR,C_BAR,S_BAR,U_BAR,D_BAR,GL,D,U,S,C,B,T
    
C--- CHOOSE PROTON STRUCTURE FUNCTIONS AND THEIR FACTORIZATION SCHEME
C--- ISCHEME:  0            1
C---           MSBAR        DIS
c     ICASE=IPDFLIB
c     IF(ICASE.EQ.1)THEN
c      Q2=Q**2
c      CALL PDGRV(ISET,X,Q2,PDF)
c     ELSE
c      PARM(1)='NPTYPE'
c      PARM(2)='NGROUP'
c      PARM(3)='NSET'
c      VALUE(1)=1.D0
c      VALUE(2)=NGROUP
c      VALUE(3)=NSET
c      CALL PDFSET(PARM,VALUE)
c      CALL PFTOPDG(X,Q,PDF)
c     ENDIF
c     RETURN
c     END

c     subroutine struc(x,q,pdf)
c     implicit double precision (a-h,o-z)
c     dimension pdf(-6:6), value(20)
c     character*20 parm(20)
c     character prefix*50
c     common/pdflib/ipdflib,ngroup,nset,iseterr
c     icase=ipdflib
c     ngroup = ipdflib
c     if(icase.eq.1)then
c      q2=q**2
c      call pdgrv(iset,x,q2,pdf)
c      return
c     endif
c     if(ngroup.gt.0)then
c      parm(1)='nptype'
c      parm(2)='ngroup'
c      parm(3)='nset'
c      value(1)=1.d0
c      value(2)=dfloat(ngroup)
c      value(3)=dfloat(nset)
c      call pdfset(parm,value)
c      call pftopdg(x,q,pdf)
c     elseif(ngroup.eq.-1)then
c      call SetCtq6(nset)
c      pdf(6)  = 0
c      pdf(-6) = 0
c      do i=-5,5
c       j = i
c       if(i.eq.1)j=2
c       if(i.eq.2)j=1
c       if(i.eq.-1)j=-2
c       if(i.eq.-2)j=-1
c       pdf(j) = x*Ctq6Pdf(i,x,q)
c      enddo
c     elseif(ngroup.eq.-2)then
c      mode = nset
c      call mrst2001(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
c      pdf(-6) = 0
c      pdf(-5) = bot
c      pdf(-4) = chm
c      pdf(-3) = str
c      pdf(-2) = usea
c      pdf(-1) = dsea
c      pdf(0) = glu
c      pdf(1) = dnv + dsea
c      pdf(2) = upv + usea
c      pdf(3) = str
c      pdf(4) = chm
c      pdf(5) = bot
c      pdf(6) = 0
c     else
c      if(iseterr.eq.0)then
c       iset = 0
c       if(nset.eq.0)then
c        prefix = "Grids/mstw2008lo"
c       elseif(nset.eq.1)then
c        prefix = "Grids/mstw2008nlo"
c       else
c        prefix = "Grids/mstw2008nnlo"
c       endif
c      elseif(iseterr.gt.0)then
c       iset = iseterr
c       if(nset.eq.0)then
c        prefix = "Grids/mstw2008lo.90cl"
c       else
c        prefix = "Grids/mstw2008nlo.90cl"
c       endif
c      else
c       iset =-iseterr
c       if(nset.eq.0)then
c        prefix = "Grids/mstw2008lo.68cl"
c       else
c        prefix = "Grids/mstw2008nlo.68cl"
c       endif
c      endif
C--   First the traditional MRST-like interface
C--   (but note the "sbar", "cbar", "bbar" and "phot").
c      CALL GetAllPDFs(prefix,iset,x,q,upv,dnv,usea,dsea,str,sbar,
c    &        chm,cbar,bot,bbar,glu,phot)
c      pdf(-6) = 0
c      pdf(-5) = bbar
c      pdf(-4) = cbar
c      pdf(-3) = sbar
c      pdf(-2) = usea
c      pdf(-1) = dsea
c      pdf(0) = glu
c      pdf(1) = dnv + dsea
c      pdf(2) = upv + usea
c      pdf(3) = str
c      pdf(4) = chm
c      pdf(5) = bot
c      pdf(6) = 0
c     endif
c     pdf( 6) = 0
c     pdf(-6) = 0
c     return
c     end

      subroutine struc(x,q,pdf)
      implicit double precision (a-h,o-z)
      dimension pdf(-6:6), value(20)
      common/pdflib/ipdflib,ngroup,nset,iseterr
      if(ipdflib.eq.1)then
       q2=q**2
       call pdgrv(nset,x,q2,pdf)
      elseif(ipdflib.eq.2)then
       call SetCtq6(nset)
       pdf(6)  = 0
       pdf(-6) = 0
       do i=-5,5
        j = i
        if(i.eq.1)j=2
        if(i.eq.2)j=1
        if(i.eq.-1)j=-2
        if(i.eq.-2)j=-1
        pdf(j) = x*Ctq6Pdf(i,x,q)
       enddo
      else
       xcut = 0.d-7
       if(x.ge.xcut)then
        call evolvePDF(x,q,pdf)
       else
        do i = -6,6
         pdf(i) = 0
        enddo
       endif
      endif
      pdf( 6) = 0
      pdf(-6) = 0
      return
      end

      subroutine pdfset(pathname,pdfname)
      implicit double precision (a-h,o-z)
      character*100 pdfname, pathname
      common/pdflib/ipdflib,ngroup,nset,iseterr
      if(ipdflib.eq.0)then
       call SetPDFpath(pathname)
       call InitPDFsetByName(pdfname)
       call InitPDF(nset)
      endif
      return
      end
