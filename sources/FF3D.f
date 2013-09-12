C ----------------------------------------------------------------------
C   IMSL ROUTINE NAME   - FFT3D                                         FFT30010
C                                                                       FFT30020
C-----------------------------------------------------------------------FFT30030
C                                                                       FFT30040
C   COMPUTER            - UNIVAC/SINGLE                                 FFT30050
C                                                                       FFT30060
C   LATEST REVISION     - JUNE 1, 1980                                  FFT30070
C                                                                       FFT30080
C   PURPOSE             - COMPUTE THE FAST FOURIER TRANSFORM OF         FFT30090
C                           A COMPLEX VALUED 1,2 OR 3 DIMENSIONAL       FFT30100
C                           ARRAY                                       FFT30110
C                                                                       FFT30120
C   USAGE               - CALL FFT3D (A,IA1,IA2,N1,N2,N3,IJOB,IWK,RWK,  FFT30130
C                           CWK)                                        FFT30140
C                                                                       FFT30150
C   ARGUMENTS    A      - COMPLEX ARRAY. A MAY BE A THREE               FFT30160
C                           DIMENSIONAL ARRAY OF DIMENSION N1 BY N2     FFT30170
C                           BY N3, A TWO DIMENSIONAL ARRAY OF           FFT30180
C                           DIMENSION N1 BY N2, OR A VECTOR OF          FFT30190
C                           LENGTH N1. ON INPUT A CONTAINS THE          FFT30200
C                           ARRAY TO BE TRANSFORMED. ON OUTPUT          FFT30210
C                           A IS REPLACED BY THE FOURIER OR             FFT30220
C                           INVERSE FOURIER TRANSFORM (DEPENDING ON     FFT30230
C                           THE VALUE OF INPUT PARAMETER IJOB).         FFT30240
C                IA1    - FIRST DIMENSION OF THE ARRAY A EXACTLY        FFT30250
C                           AS SPECIFIED IN THE DIMENSION STATEMENT     FFT30260
C                           IN THE CALLING PROGRAM. (INPUT)             FFT30270
C                IA2    - SECOND DIMENSION OF THE ARRAY A EXACTLY       FFT30280
C                           AS SPECIFIED IN THE DIMENSION STATEMENT     FFT30290
C                           IN THE CALLING PROGRAM. (INPUT)             FFT30300
C                N1     - LIMITS ON THE FIRST, SECOND, AND THIRD        FFT30310
C                N2         SUBSCRIPTS OF THE ARRAY A, RESPECTIVELY.    FFT30320
C                N3         (INPUT)                                     FFT30330
C                IJOB   - INPUT OPTION PARAMETER.                       FFT30340
C                           IF IJOB IS POSITIVE, THE FAST FOURIER       FFT30350
C                             TRANSFORM OF A IS TO BE CALCULATED.       FFT30360
C                           IF IJOB IS NEGATIVE, THE INVERSE            FFT30370
C                             FAST FOURIER TRANSFORM OF A IS TO BE      FFT30380
C                             CALCULATED.                               FFT30390
C                IWK    - INTEGER WORK VECTOR OF LENGTH                 FFT30400
C                           6*MAX(N1,N2,N3)+150.                        FFT30410
C                RWK    - REAL WORK VECTOR OF LENGTH                    FFT30420
C                           6*MAX(N1,N2,N3)+150.                        FFT30430
C                CWK    - COMPLEX WORK VECTOR OF LENGTH                 FFT30440
C                           MAX(N2,N3).                                 FFT30450
C                                                                       FFT30460
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         FFT30470
C                       - SINGLE/H32,H48,H60                            FFT30480
C                                                                       FFT30490
C   REQD. IMSL ROUTINES - FFTCC                                         FFT30500
C                                                                       FFT30510
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           FFT30520
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      FFT30530
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  FFT30540
C                                                                       FFT30550
C   REMARKS  1.  IF IJOB IS POSITIVE, FFT3D CALCULATES THE FOURIER      FFT30560
C                TRANSFORM, X, ACCORDING TO THE FOLLOWING FORMULA       FFT30570
C                                                                       FFT30580
C                  X(I+1,J+1,K+1)=TRIPLE SUM OF A(L+1,M+1,N+1)*         FFT30590
C                  EXP(2*PI*SQRT(-1)*(I*L/N1+J*M/N2+K*N/N3))            FFT30600
C                  WITH L=0...N1-1, M=0...N2-1, N=0...N3-1              FFT30610
C                  AND PI=3.1415...                                     FFT30620
C                                                                       FFT30630
C                IF IJOB IS NEGATIVE, FFT3D CALCULATES THE INVERSE      FFT30640
C                FOURIER TRANSFORM, X, ACCORDING TO THE FOLLOWING       FFT30650
C                FORMULA                                                FFT30660
C                                                                       FFT30670
C                  X(I+1,J+1,K+1)=1/(N1*N2*N3)*TRIPLE SUM OF            FFT30680
C                  A(L+1,M+1,N+1)*                                      FFT30690
C                  EXP(-2*PI*SQRT(-1)*(I*L/N1+J*M/N2+K*N/N3))           FFT30700
C                  WITH L=0...N1-1, M=0...N2-1, N=0...N3-1              FFT30710
C                  AND PI=3.1415...                                     FFT30720
C                                                                       FFT30730
C                NOTE THAT X OVERWRITES A ON OUTPUT.                    FFT30740
C            2.  IF A IS A TWO DIMENSIONAL ARRAY, SET N3 = 1.           FFT30750
C                IF A IS A ONE DIMENSIONAL ARRAY (VECTOR),              FFT30760
C                SET IA2 = N2 = N3 = 1.                                 FFT30770
C                                                                       FFT30780
C   COPYRIGHT           - 1980 BY IMSL, INC. ALL RIGHTS RESERVED.       FFT30790
C                                                                       FFT30800
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN FFT30810
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    FFT30820
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        FFT30830
C                                                                       FFT30840
C-----------------------------------------------------------------------FFT30850
C                                                                       FFT30860
      SUBROUTINE FFT3D  (A,IA1,IA2,N1,N2,N3,IJOB,IWK,RWK,CWK)           FFT30870
C                                  SPECIFICATIONS FOR ARGUMENTS         FFT30880
      INTEGER            IA1,IA2,N1,N2,N3,IJOB,IWK(1)                   FFT30890
      REAL               RWK(1)                                         FFT30900
      COMPLEX            A(IA1,IA2,N3),CWK(1)                           FFT30910
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   FFT30920
      INTEGER            I,J,K,L,M,N                                    FFT30930
      REAL               R123                                           FFT30940
      COMPLEX            C123                                           FFT30950
      
      Save 		 
C                                  FIRST EXECUTABLE STATEMENT           FFT30960
      IF (IJOB.GT.0) GO TO 10                                           FFT30970
C                                  INVERSE TRANSFORM                    FFT30980
      DO 5 I=1,N1                                                       FFT30990
      DO 5 J=1,N2                                                       FFT31000
      DO 5 K=1,N3                                                       FFT31010
         A(I,J,K) = CONJG(A(I,J,K))                                     FFT31020
    5 CONTINUE                                                          FFT31030
C                                  TRANSFORM THIRD SUBSCRIPT            FFT31040
   10 DO 25 L=1,N1                                                      FFT31050
      DO 25 M=1,N2                                                      FFT31060
         DO 15 N=1,N3                                                   FFT31070
            CWK(N) = A(L,M,N)                                           FFT31080
   15    CONTINUE                                                       FFT31090
         CALL FFTCC (CWK,N3,IWK,RWK)                                    FFT31100
         DO 20 K=1,N3                                                   FFT31110
            A(L,M,K) = CWK(K)                                           FFT31120
   20    CONTINUE                                                       FFT31130
   25 CONTINUE                                                          FFT31140
C                                  TRANSFORM SECOND SUBSCRIPT           FFT31150
      DO 40 L=1,N1                                                      FFT31160
      DO 40 K=1,N3                                                      FFT31170
         DO 30 M=1,N2                                                   FFT31180
            CWK(M) = A(L,M,K)                                           FFT31190
   30    CONTINUE                                                       FFT31200
         CALL FFTCC (CWK,N2,IWK,RWK)                                    FFT31210
         DO 35 J=1,N2                                                   FFT31220
            A(L,J,K) = CWK(J)                                           FFT31230
   35    CONTINUE                                                       FFT31240
   40 CONTINUE                                                          FFT31250
C                                  TRANSFORM FIRST SUBSCRIPT            FFT31260
      DO 45 J=1,N2                                                      FFT31270
      DO 45 K=1,N3                                                      FFT31280
         CALL FFTCC (A(1,J,K),N1,IWK,RWK)                               FFT31290
   45 CONTINUE                                                          FFT31300
      IF (IJOB.GT.0) GO TO 55                                           FFT31310
      R123 = N1*N2*N3                                                   FFT31320
      C123 = CMPLX(R123,0.0)                                            FFT31330
      DO 50 I=1,N1                                                      FFT31340
      DO 50 J=1,N2                                                      FFT31350
      DO 50 K=1,N3                                                      FFT31360
         A(I,J,K) = CONJG(A(I,J,K))/C123                                FFT31370
   50 CONTINUE                                                          FFT31380
   55 RETURN                                                            FFT31390
      END                                                               FFT31400
C   IMSL ROUTINE NAME   - FFTCC                                         FFTP0010
C                                                                       FFTP0020
C-----------------------------------------------------------------------FFTP0030
C                                                                       FFTP0040
C   COMPUTER            - UNIVAC/SINGLE                                 FFTP0050
C                                                                       FFTP0060
C   LATEST REVISION     - JANUARY 1, 1978                               FFTP0070
C                                                                       FFTP0080
C   PURPOSE             - COMPUTE THE FAST FOURIER TRANSFORM OF A       FFTP0090
C                           COMPLEX VALUED SEQUENCE                     FFTP0100
C                                                                       FFTP0110
C   USAGE               - CALL FFTCC (A,N,IWK,WK)                       FFTP0120
C                                                                       FFTP0130
C   ARGUMENTS    A      - COMPLEX VECTOR OF LENGTH N. ON INPUT A        FFTP0140
C                           CONTAINS THE COMPLEX VALUED SEQUENCE TO BE  FFTP0150
C                           TRANSFORMED. ON OUTPUT A IS REPLACED BY THE FFTP0160
C                           FOURIER TRANSFORM.                          FFTP0170
C                N      - INPUT NUMBER OF DATA POINTS TO BE             FFTP0180
C                           TRANSFORMED. N MAY BE ANY POSITIVE          FFTP0190
C                           INTEGER.                                    FFTP0200
C                IWK    - INTEGER WORK VECTOR OF LENGTH 6*N+150.        FFTP0210
C                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS) FFTP0220
C                WK     - REAL WORK VECTOR OF LENGTH 6*N+150.           FFTP0230
C                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS) FFTP0240
C                                                                       FFTP0250
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         FFTP0260
C                       - SINGLE/H36,H48,H60                            FFTP0270
C                                                                       FFTP0280
C   REQD. IMSL ROUTINES - NONE REQUIRED                                 FFTP0290
C                                                                       FFTP0300
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           FFTP0310
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      FFTP0320
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  FFTP0330
C                                                                       FFTP0340
C   REMARKS  1.  FFTCC COMPUTES THE FOURIER TRANSFORM, X, ACCORDING     FFTP0350
C                TO THE FOLLOWING FORMULA;                              FFTP0360
C                                                                       FFTP0370
C                  X(K+1) = SUM FROM J = 0 TO N-1 OF                    FFTP0380
C                           A(J+1)*CEXP((0.0,(2.0*PI*J*K)/N))           FFTP0390
C                  FOR K=0,1,...,N-1 AND PI=3.1415...                   FFTP0400
C                                                                       FFTP0410
C                NOTE THAT X OVERWRITES A ON OUTPUT.                    FFTP0420
C            2.  FFTCC CAN BE USED TO COMPUTE                           FFTP0430
C                                                                       FFTP0440
C                  X(K+1) = (1/N)*SUM FROM J = 0 TO N-1 OF              FFTP0450
C                           A(J+1)*CEXP((0.0,(-2.0*PI*J*K)/N))          FFTP0460
C                  FOR K=0,1,...,N-1 AND PI=3.1415...                   FFTP0470
C                                                                       FFTP0480
C                BY PERFORMING THE FOLLOWING STEPS;                     FFTP0490
C                                                                       FFTP0500
C                     DO 10 I=1,N                                       FFTP0510
C                        A(I) = CONJG(A(I))                             FFTP0520
C                  10 CONTINUE                                          FFTP0530
C                     CALL FFTCC (A,N,IWK,WK)                           FFTP0540
C                     DO 20 I=1,N                                       FFTP0550
C                        A(I) = CONJG(A(I))/N                           FFTP0560
C                  20 CONTINUE                                          FFTP0570
C                                                                       FFTP0580
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       FFTP0590
C                                                                       FFTP0600
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN FFTP0610
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    FFTP0620
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        FFTP0630
C                                                                       FFTP0640
C-----------------------------------------------------------------------FFTP0650
C                                                                       FFTP0660
      SUBROUTINE FFTCC (A,N,IWK,WK)                                     FFTP0670
C                                  SPECIFICATIONS FOR ARGUMENTS         FFTP0680
      INTEGER            N,IWK(1)                                       FFTP0690
      REAL               WK(1)                                          FFTP0700
      COMPLEX            A(N)                                           FFTP0710
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   FFTP0720
      INTEGER            I,IAM,IAP,IBM,IBP,IC,ICC,ICF,ICK,ID,IDM1,II,   FFTP0730
     1                   IJA,IKB,IKT,ILL,IM,IRD,ISF,ISK,ISP,ISS,ITA,ITB,FFTP0740
     2                   J,JA,JF,JJ,JK,K,K0,K1,K2,K3,KA,KB,KD2,KF,KH,KN,FFTP0750
     3                   KT,KTP,L,L1,M,MM,MM1,MP                        FFTP0760
      REAL               CM,SM,C1,C2,C3,S1,S2,S3,C30,RAD,A0,A1,A4,B4,   FFTP0770
     1                   A2,A3,B0,B1,B2,B3,ZERO,HALF,ONE,TWO,Z0(2),     FFTP0780
     2                   Z1(2),Z2(2),Z3(2),Z4(2)                        FFTP0790
      COMPLEX            ZA0,ZA1,ZA2,ZA3,ZA4,AK2                        FFTP0800
      EQUIVALENCE        (ZA0,Z0(1)),(ZA1,Z1(1)),(ZA2,Z2(1)),           FFTP0810
     1                   (ZA3,Z3(1)),(A0,Z0(1)),(B0,Z0(2)),(A1,Z1(1)),  FFTP0820
     2                   (B1,Z1(2)),(A2,Z2(1)),(B2,Z2(2)),(A3,Z3(1)),   FFTP0830
     3                   (B3,Z3(2)),(ZA4,Z4(1)),(Z4(1),A4),(Z4(2),B4)   FFTP0840
     
       save
     
      DATA               RAD/6.28318531/,                               FFTP0850
     1                   C30/.866025404/                                FFTP0860
      DATA               ZERO,HALF,ONE,TWO/0.0,0.5,1.0,2.0/             FFTP0870	
C                                  FIRST EXECUTABLE STATEMENT           FFTP0880
      IF (N .EQ. 1) GO TO 9005                                          FFTP0890
      K = N                                                             FFTP0900
      M = 0                                                             FFTP0910
      J = 2                                                             FFTP0920
      JJ = 4                                                            FFTP0930
      JF = 0                                                            FFTP0940
C                                  DETERMINE THE SQUARE FACTORS OF N    FFTP0950
      IWK(1) = 1                                                        FFTP0960
    5 I = K/JJ                                                          FFTP0970
      IF (I*JJ .NE. K) GO TO 10                                         FFTP0980
      M = M+1                                                           FFTP0990
      IWK(M+1) = J                                                      FFTP1000
      K = I                                                             FFTP1010
      GO TO 5                                                           FFTP1020
   10 J = J + 2                                                         FFTP1030
      IF (J .EQ. 4) J = 3                                               FFTP1040
      JJ = J * J                                                        FFTP1050
      IF (JJ .LE. K) GO TO 5                                            FFTP1060
      KT = M                                                            FFTP1070
C                                  DETERMINE THE REMAINING FACTORS OF N FFTP1080
      J = 2                                                             FFTP1090
   15 I = K / J                                                         FFTP1100
      IF (I*J .NE. K) GO TO 20                                          FFTP1110
      M = M + 1                                                         FFTP1120
      IWK(M+1) = J                                                      FFTP1130
      K = I                                                             FFTP1140
      GO TO 15                                                          FFTP1150
   20 J = J + 1                                                         FFTP1160
      IF (J .EQ. 3) GO TO 15                                            FFTP1170
      J = J + 1                                                         FFTP1180
      IF(J.LE.K) GO TO 15                                               FFTP1190
      K = IWK(M+1)                                                      FFTP1200
      IF (IWK(KT+1) .GT. IWK(M+1)) K = IWK(KT+1)                        FFTP1210
      IF(KT.LE.0) GO TO 30                                              FFTP1220
      KTP = KT + 2                                                      FFTP1230
      DO 25  I = 1,KT                                                   FFTP1240
         J = KTP - I                                                    FFTP1250
         M = M+1                                                        FFTP1260
         IWK(M+1) = IWK(J)                                              FFTP1270
   25 CONTINUE                                                          FFTP1280
   30 MP = M+1                                                          FFTP1290
      IC = MP+1                                                         FFTP1300
      ID = IC+MP                                                        FFTP1310
      ILL = ID+MP                                                       FFTP1320
      IRD = ILL+MP+1                                                    FFTP1330
      ICC = IRD+MP                                                      FFTP1340
      ISS = ICC+MP                                                      FFTP1350
      ICK = ISS+MP                                                      FFTP1360
      ISK = ICK+K                                                       FFTP1370
      ICF = ISK+K                                                       FFTP1380
      ISF = ICF+K                                                       FFTP1390
      IAP = ISF+K                                                       FFTP1400
      KD2 = (K-1) / 2 + 1                                               FFTP1410
      IBP = IAP + KD2                                                   FFTP1420
      IAM = IBP + KD2                                                   FFTP1430
      IBM = IAM + KD2                                                   FFTP1440
      MM1 = M-1                                                         FFTP1450
      I = 1                                                             FFTP1460
   35 L = MP - I                                                        FFTP1470
      J = IC - I                                                        FFTP1480
      IWK(ILL+L) = 0                                                    FFTP1490
      IF ((IWK(J-1) + IWK(J)) .EQ. 4) IWK(ILL+L) = 1                    FFTP1500
      IF (IWK(ILL+L) .EQ. 0) GO TO 40                                   FFTP1510
      I = I + 1                                                         FFTP1520
      L = L - 1                                                         FFTP1530
      IWK(ILL+L) = 0                                                    FFTP1540
   40 I = I + 1                                                         FFTP1550
      IF(I.LE.MM1) GO TO 35                                             FFTP1560
      IWK(ILL+1) = 0                                                    FFTP1570
      IWK(ILL+MP) = 0                                                   FFTP1580
      IWK(IC) = 1                                                       FFTP1590
      IWK(ID) = N                                                       FFTP1600
      DO 45  J = 1,M                                                    FFTP1610
         K = IWK(J+1)                                                   FFTP1620
         IWK(IC+J) = IWK(IC+J-1) * K                                    FFTP1630
         IWK(ID+J) = IWK(ID+J-1) / K                                    FFTP1640
         WK(IRD+J) = RAD/IWK(IC+J)                                      FFTP1650
         C1 = RAD/K                                                     FFTP1660
         IF (K .LE. 2) GO TO 45                                         FFTP1670
         WK(ICC+J) = COS(C1)                                            FFTP1680
         WK(ISS+J) = SIN(C1)                                            FFTP1690
   45 CONTINUE                                                          FFTP1700
      MM = M                                                            FFTP1710
      IF (IWK(ILL+M) .EQ. 1) MM = M - 1                                 FFTP1720
      IF (MM .LE. 1) GO TO 50                                           FFTP1730
      SM = IWK(IC+MM-2) * WK(IRD+M)                                     FFTP1740
      CM = COS(SM)                                                      FFTP1750
      SM = SIN(SM)                                                      FFTP1760
   50 KB = 0                                                            FFTP1770
      KN = N                                                            FFTP1780
      JJ = 0                                                            FFTP1790
      I = 1                                                             FFTP1800
      C1 = ONE                                                          FFTP1810
      S1 = ZERO                                                         FFTP1820
      L1 = 1                                                            FFTP1830
   55 IF (IWK(ILL+I+1) .EQ. 1) GO TO 60                                 FFTP1840
      KF = IWK(I+1)                                                     FFTP1850
      GO TO 65                                                          FFTP1860
   60 KF = 4                                                            FFTP1870
      I = I+1                                                           FFTP1880
   65 ISP = IWK(ID+I)                                                   FFTP1890
      IF (L1 .EQ. 1) GO TO 70                                           FFTP1900
      S1 = JJ * WK(IRD+I)                                               FFTP1910
      C1 = COS(S1)                                                      FFTP1920
      S1 = SIN(S1)                                                      FFTP1930
C                                  FACTORS OF 2, 3, AND 4 ARE           FFTP1940
C                                  HANDLED SEPARATELY.                  FFTP1950
   70 IF (KF .GT. 4) GO TO 140                                          FFTP1960
      GO TO (75,75,90,115), KF                                          FFTP1970
   75 K0 = KB + ISP                                                     FFTP1980
      K2 = K0 + ISP                                                     FFTP1990
      IF (L1 .EQ. 1) GO TO 85                                           FFTP2000
   80 K0 = K0 - 1                                                       FFTP2010
      IF (K0 .LT. KB) GO TO 190                                         FFTP2020
      K2 = K2 - 1                                                       FFTP2030
      ZA4 = A(K2+1)                                                     FFTP2040
      A0 = A4*C1-B4*S1                                                  FFTP2050
      B0 = A4*S1+B4*C1                                                  FFTP2060
      A(K2+1) = A(K0+1)-ZA0                                             FFTP2070
      A(K0+1) = A(K0+1)+ZA0                                             FFTP2080
      GO TO 80                                                          FFTP2090
   85 K0 = K0 - 1                                                       FFTP2100
      IF (K0 .LT. KB) GO TO 190                                         FFTP2110
      K2 = K2 - 1                                                       FFTP2120
      AK2 = A(K2+1)                                                     FFTP2130
      A(K2+1) = A(K0+1)-AK2                                             FFTP2140
      A(K0+1) = A(K0+1)+AK2                                             FFTP2150
      GO TO 85                                                          FFTP2160
   90 IF (L1 .EQ. 1) GO TO 95                                           FFTP2170
      C2 = C1 * C1 - S1 * S1                                            FFTP2180
      S2 = TWO * C1 * S1                                                FFTP2190
   95 JA = KB + ISP - 1                                                 FFTP2200
      KA = JA + KB                                                      FFTP2210
      IKB = KB+1                                                        FFTP2220
      IJA = JA+1                                                        FFTP2230
      DO 110 II = IKB,IJA                                               FFTP2240
         K0 = KA - II + 1                                               FFTP2250
         K1 = K0 + ISP                                                  FFTP2260
         K2 = K1 + ISP                                                  FFTP2270
         ZA0 = A(K0+1)                                                  FFTP2280
         IF (L1 .EQ. 1) GO TO 100                                       FFTP2290
         ZA4 = A(K1+1)                                                  FFTP2300
         A1 = A4*C1-B4*S1                                               FFTP2310
         B1 = A4*S1+B4*C1                                               FFTP2320
         ZA4 = A(K2+1)                                                  FFTP2330
         A2 = A4*C2-B4*S2                                               FFTP2340
         B2 = A4*S2+B4*C2                                               FFTP2350
         GO TO 105                                                      FFTP2360
  100    ZA1 = A(K1+1)                                                  FFTP2370
         ZA2 = A(K2+1)                                                  FFTP2380
  105    A(K0+1) = CMPLX(A0+A1+A2,B0+B1+B2)                             FFTP2390
         A0 = -HALF * (A1+A2) + A0                                      FFTP2400
         A1 = (A1-A2) * C30                                             FFTP2410
         B0 = -HALF * (B1+B2) + B0                                      FFTP2420
         B1 = (B1-B2) * C30                                             FFTP2430
         A(K1+1) = CMPLX(A0-B1,B0+A1)                                   FFTP2440
         A(K2+1) = CMPLX(A0+B1,B0-A1)                                   FFTP2450
  110 CONTINUE                                                          FFTP2460
      GO TO 190                                                         FFTP2470
  115 IF (L1 .EQ. 1) GO TO 120                                          FFTP2480
      C2 = C1 * C1 - S1 * S1                                            FFTP2490
      S2 = TWO * C1 * S1                                                FFTP2500
      C3 = C1 * C2 - S1 * S2                                            FFTP2510
      S3 = S1 * C2 + C1 * S2                                            FFTP2520
  120 JA = KB + ISP - 1                                                 FFTP2530
      KA = JA + KB                                                      FFTP2540
      IKB = KB+1                                                        FFTP2550
      IJA = JA+1                                                        FFTP2560
      DO 135 II = IKB,IJA                                               FFTP2570
         K0 = KA - II + 1                                               FFTP2580
         K1 = K0 + ISP                                                  FFTP2590
         K2 = K1 + ISP                                                  FFTP2600
         K3 = K2 + ISP                                                  FFTP2610
         ZA0 = A(K0+1)                                                  FFTP2620
         IF (L1 .EQ. 1) GO TO 125                                       FFTP2630
         ZA4 = A(K1+1)                                                  FFTP2640
         A1 = A4*C1-B4*S1                                               FFTP2650
         B1 = A4*S1+B4*C1                                               FFTP2660
         ZA4 = A(K2+1)                                                  FFTP2670
         A2 = A4*C2-B4*S2                                               FFTP2680
         B2 = A4*S2+B4*C2                                               FFTP2690
         ZA4 = A(K3+1)                                                  FFTP2700
         A3 = A4*C3-B4*S3                                               FFTP2710
         B3 = A4*S3+B4*C3                                               FFTP2720
         GO TO 130                                                      FFTP2730
  125    ZA1 = A(K1+1)                                                  FFTP2740
         ZA2 = A(K2+1)                                                  FFTP2750
         ZA3 = A(K3+1)                                                  FFTP2760
  130    A(K0+1) = CMPLX(A0+A2+A1+A3,B0+B2+B1+B3)                       FFTP2770
         A(K1+1) = CMPLX(A0+A2-A1-A3,B0+B2-B1-B3)                       FFTP2780
         A(K2+1) = CMPLX(A0-A2-B1+B3,B0-B2+A1-A3)                       FFTP2790
         A(K3+1) = CMPLX(A0-A2+B1-B3,B0-B2-A1+A3)                       FFTP2800
  135 CONTINUE                                                          FFTP2810
      GO TO 190                                                         FFTP2820
  140 JK = KF - 1                                                       FFTP2830
      KH = JK/2                                                         FFTP2840
      K3 = IWK(ID+I-1)                                                  FFTP2850
      K0 = KB + ISP                                                     FFTP2860
      IF (L1 .EQ. 1) GO TO 150                                          FFTP2870
      K = JK - 1                                                        FFTP2880
      WK(ICF+1) = C1                                                    FFTP2890
      WK(ISF+1) = S1                                                    FFTP2900
      DO 145 J = 1,K                                                    FFTP2910
         WK(ICF+J+1) = WK(ICF+J) * C1 - WK(ISF+J) * S1                  FFTP2920
         WK(ISF+J+1) = WK(ICF+J) * S1 + WK(ISF+J) * C1                  FFTP2930
  145 CONTINUE                                                          FFTP2940
  150 IF (KF .EQ. JF) GO TO 160                                         FFTP2950
      C2 = WK(ICC+I)                                                    FFTP2960
      WK(ICK+1) = C2                                                    FFTP2970
      WK(ICK+JK) = C2                                                   FFTP2980
      S2 = WK(ISS+I)                                                    FFTP2990
      WK(ISK+1) = S2                                                    FFTP3000
      WK(ISK+JK) = -S2                                                  FFTP3010
      DO 155 J = 1,KH                                                   FFTP3020
         K = JK - J                                                     FFTP3030
         WK(ICK+K) = WK(ICK+J) * C2 - WK(ISK+J) * S2                    FFTP3040
         WK(ICK+J+1) = WK(ICK+K)                                        FFTP3050
         WK(ISK+J+1) = WK(ICK+J) * S2 + WK(ISK+J) * C2                  FFTP3060
         WK(ISK+K) = -WK(ISK+J+1)                                       FFTP3070
  155 CONTINUE                                                          FFTP3080
  160 K0 = K0 - 1                                                       FFTP3090
      K1 = K0                                                           FFTP3100
      K2 = K0 + K3                                                      FFTP3110
      ZA0 = A(K0+1)                                                     FFTP3120
      A3 = A0                                                           FFTP3130
      B3 = B0                                                           FFTP3140
      DO 175 J = 1,KH                                                   FFTP3150
         K1 = K1 + ISP                                                  FFTP3160
         K2 = K2 - ISP                                                  FFTP3170
         IF (L1 .EQ. 1) GO TO 165                                       FFTP3180
         K = KF - J                                                     FFTP3190
         ZA4 = A(K1+1)                                                  FFTP3200
         A1 = A4*WK(ICF+J)-B4*WK(ISF+J)                                 FFTP3210
         B1 = A4*WK(ISF+J)+B4*WK(ICF+J)                                 FFTP3220
         ZA4 = A(K2+1)                                                  FFTP3230
         A2 = A4*WK(ICF+K)-B4*WK(ISF+K)                                 FFTP3240
         B2 = A4*WK(ISF+K)+B4*WK(ICF+K)                                 FFTP3250
         GO TO 170                                                      FFTP3260
  165    ZA1 = A(K1+1)                                                  FFTP3270
         ZA2 = A(K2+1)                                                  FFTP3280
  170    WK(IAP+J) = A1 + A2                                            FFTP3290
         WK(IAM+J) = A1 - A2                                            FFTP3300
         WK(IBP+J) = B1 + B2                                            FFTP3310
         WK(IBM+J) = B1 - B2                                            FFTP3320
         A3 = A1 + A2 + A3                                              FFTP3330
         B3 = B1 + B2 + B3                                              FFTP3340
  175 CONTINUE                                                          FFTP3350
      A(K0+1) = CMPLX(A3,B3)                                            FFTP3360
      K1 = K0                                                           FFTP3370
      K2 = K0 + K3                                                      FFTP3380
      DO 185 J = 1,KH                                                   FFTP3390
         K1 = K1 + ISP                                                  FFTP3400
         K2 = K2 - ISP                                                  FFTP3410
         JK = J                                                         FFTP3420
         A1 = A0                                                        FFTP3430
         B1 = B0                                                        FFTP3440
         A2 = ZERO                                                      FFTP3450
         B2 = ZERO                                                      FFTP3460
         DO 180  K = 1,KH                                               FFTP3470
            A1 = A1 + WK(IAP+K) * WK(ICK+JK)                            FFTP3480
            A2 = A2 + WK(IAM+K) * WK(ISK+JK)                            FFTP3490
            B1 = B1 + WK(IBP+K) * WK(ICK+JK)                            FFTP3500
            B2 = B2 + WK(IBM+K) * WK(ISK+JK)                            FFTP3510
            JK = JK + J                                                 FFTP3520
            IF (JK .GE. KF) JK = JK - KF                                FFTP3530
  180    CONTINUE                                                       FFTP3540
         A(K1+1) = CMPLX(A1-B2,B1+A2)                                   FFTP3550
         A(K2+1) = CMPLX(A1+B2,B1-A2)                                   FFTP3560
  185 CONTINUE                                                          FFTP3570
      IF (K0 .GT. KB) GO TO 160                                         FFTP3580
      JF = KF                                                           FFTP3590
  190 IF ( I .GE. MM ) GO TO 195                                        FFTP3600
      I = I + 1                                                         FFTP3610
      GO TO 55                                                          FFTP3620
  195 I = MM                                                            FFTP3630
      L1 = 0                                                            FFTP3640
      KB = IWK(ID+I-1) + KB                                             FFTP3650
      IF (KB .GE. KN) GO TO 215                                         FFTP3660
  200 JJ = IWK(IC+I-2) + JJ                                             FFTP3670
      IF (JJ .LT. IWK(IC+I-1)) GO TO 205                                FFTP3680
      I = I - 1                                                         FFTP3690
      JJ = JJ - IWK(IC+I)                                               FFTP3700
      GO TO 200                                                         FFTP3710
  205 IF (I .NE. MM) GO TO 210                                          FFTP3720
      C2 = C1                                                           FFTP3730
      C1 = CM * C1 - SM * S1                                            FFTP3740
      S1 = SM * C2 + CM * S1                                            FFTP3750
      GO TO 70                                                          FFTP3760
  210 IF (IWK(ILL+I) .EQ. 1) I = I + 1                                  FFTP3770
      GO TO 55                                                          FFTP3780
  215 I = 1                                                             FFTP3790
      JA = KT - 1                                                       FFTP3800
      KA = JA + 1                                                       FFTP3810
      IF(JA.LT.1) GO TO 225                                             FFTP3820
      DO 220  II = 1,JA                                                 FFTP3830
         J = KA - II                                                    FFTP3840
         IWK(J+1) = IWK(J+1) - 1                                        FFTP3850
         I = IWK(J+1) + I                                               FFTP3860
  220 CONTINUE                                                          FFTP3870
C                                  THE RESULT IS NOW PERMUTED TO        FFTP3880
C                                  NORMAL ORDER.                        FFTP3890
  225 IF (KT .LE. 0) GO TO 270                                          FFTP3900
      J = 1                                                             FFTP3910
      I = 0                                                             FFTP3920
      KB = 0                                                            FFTP3930
  230 K2 = IWK(ID+J) + KB                                               FFTP3940
      K3 = K2                                                           FFTP3950
      JJ = IWK(IC+J-1)                                                  FFTP3960
      JK = JJ                                                           FFTP3970
      K0 = KB + JJ                                                      FFTP3980
      ISP = IWK(IC+J) - JJ                                              FFTP3990
  235 K = K0 + JJ                                                       FFTP4000
  240 ZA4 = A(K0+1)                                                     FFTP4010
      A(K0+1) = A(K2+1)                                                 FFTP4020
      A(K2+1) = ZA4                                                     FFTP4030
      K0 = K0 + 1                                                       FFTP4040
      K2 = K2 + 1                                                       FFTP4050
      IF (K0 .LT. K) GO TO 240                                          FFTP4060
      K0 = K0 + ISP                                                     FFTP4070
      K2 = K2 + ISP                                                     FFTP4080
      IF (K0 .LT. K3) GO TO 235                                         FFTP4090
      IF (K0 .GE. K3 + ISP) GO TO 245                                   FFTP4100
      K0 = K0 - IWK(ID+J) + JJ                                          FFTP4110
      GO TO 235                                                         FFTP4120
  245 K3 = IWK(ID+J) + K3                                               FFTP4130
      IF (K3 - KB .GE. IWK(ID+J-1)) GO TO 250                           FFTP4140
      K2 = K3 + JK                                                      FFTP4150
      JK = JK + JJ                                                      FFTP4160
      K0 = K3 - IWK(ID+J) + JK                                          FFTP4170
      GO TO 235                                                         FFTP4180
  250 IF (J .GE. KT) GO TO 260                                          FFTP4190
      K = IWK(J+1) + I                                                  FFTP4200
      J = J + 1                                                         FFTP4210
  255 I = I + 1                                                         FFTP4220
      IWK(ILL+I) = J                                                    FFTP4230
      IF (I .LT. K) GO TO 255                                           FFTP4240
      GO TO 230                                                         FFTP4250
  260 KB = K3                                                           FFTP4260
      IF (I .LE. 0) GO TO 265                                           FFTP4270
      J = IWK(ILL+I)                                                    FFTP4280
      I = I - 1                                                         FFTP4290
      GO TO 230                                                         FFTP4300
  265 IF (KB .GE. N) GO TO 270                                          FFTP4310
      J = 1                                                             FFTP4320
      GO TO 230                                                         FFTP4330
  270 JK = IWK(IC+KT)                                                   FFTP4340
      ISP = IWK(ID+KT)                                                  FFTP4350
      M = M - KT                                                        FFTP4360
      KB = ISP/JK-2                                                     FFTP4370
      IF (KT .GE. M-1 ) GO TO 9005                                      FFTP4380
      ITA = ILL+KB+1                                                    FFTP4390
      ITB = ITA+JK                                                      FFTP4400
      IDM1 = ID-1                                                       FFTP4410
      IKT = KT+1                                                        FFTP4420
      IM = M+1                                                          FFTP4430
      DO 275 J = IKT,IM                                                 FFTP4440
         IWK(IDM1+J) = IWK(IDM1+J)/JK                                   FFTP4450
  275 CONTINUE                                                          FFTP4460
      JJ = 0                                                            FFTP4470
      DO 290 J = 1,KB                                                   FFTP4480
         K = KT                                                         FFTP4490
  280    JJ = IWK(ID+K+1) + JJ                                          FFTP4500
         IF (JJ .LT. IWK(ID+K)) GO TO 285                               FFTP4510
         JJ = JJ - IWK(ID+K)                                            FFTP4520
         K = K + 1                                                      FFTP4530
         GO TO 280                                                      FFTP4540
  285    IWK(ILL+J) = JJ                                                FFTP4550
         IF (JJ .EQ. J) IWK(ILL+J) = -J                                 FFTP4560
  290 CONTINUE                                                          FFTP4570
C                                  DETERMINE THE PERMUTATION CYCLES     FFTP4580
C                                  OF LENGTH GREATER THAN OR EQUAL      FFTP4590
C                                  TO TWO.                              FFTP4600
      DO 300  J = 1,KB                                                  FFTP4610
         IF (IWK(ILL+J) .LE. 0) GO TO 300                               FFTP4620
         K2 = J                                                         FFTP4630
  295    K2 = IABS(IWK(ILL+K2))                                         FFTP4640
         IF (K2 .EQ. J) GO TO 300                                       FFTP4650
         IWK(ILL+K2) = -IWK(ILL+K2)                                     FFTP4660
         GO TO 295                                                      FFTP4670
  300 CONTINUE                                                          FFTP4680
C                                  REORDER A FOLLOWING THE              FFTP4690
C                                  PERMUTATION CYCLES                   FFTP4700
      I = 0                                                             FFTP4710
      J = 0                                                             FFTP4720
      KB = 0                                                            FFTP4730
      KN = N                                                            FFTP4740
  305 J = J + 1                                                         FFTP4750
      IF (IWK(ILL+J) .LT. 0) GO TO 305                                  FFTP4760
      K = IWK(ILL+J)                                                    FFTP4770
      K0 = JK * K + KB                                                  FFTP4780
  310 ZA4 = A(K0+I+1)                                                   FFTP4790
      WK(ITA+I) = A4                                                    FFTP4800
      WK(ITB+I) = B4                                                    FFTP4810
      I = I + 1                                                         FFTP4820
      IF (I .LT. JK) GO TO 310                                          FFTP4830
      I = 0                                                             FFTP4840
  315 K = -IWK(ILL+K)                                                   FFTP4850
      JJ = K0                                                           FFTP4860
      K0 = JK * K + KB                                                  FFTP4870
  320 A(JJ+I+1) = A(K0+I+1)                                             FFTP4880
      I = I + 1                                                         FFTP4890
      IF (I .LT. JK) GO TO 320                                          FFTP4900
      I = 0                                                             FFTP4910
      IF (K .NE. J) GO TO 315                                           FFTP4920
  325 A(K0+I+1) = CMPLX(WK(ITA+I),WK(ITB+I))                            FFTP4930
      I = I + 1                                                         FFTP4940
      IF (I .LT. JK) GO TO 325                                          FFTP4950
      I = 0                                                             FFTP4960
      IF (J .LT. K2) GO TO 305                                          FFTP4970
      J = 0                                                             FFTP4980
      KB = KB + ISP                                                     FFTP4990
      IF (KB .LT. KN) GO TO 305                                         FFTP5000
 9005 RETURN                                                            FFTP5010
      END                                                               FFTP5020
      
