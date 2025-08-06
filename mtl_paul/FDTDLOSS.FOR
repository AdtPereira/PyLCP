C*************PROGRAM FDTDLOSS.FOR*********************************************
C
C     To Implement the Finite Difference-Time Domain Method for
C     Determining the Time-Domain Response of a (N+1)-Conductor
C     Transmission Line Having Resistive Loads.
C
C     Includes Conductor Loss via Skin Effect
C
C
C             Written by:
C                  Clayton R. Paul
C                  Department of Electrical Engineering
C                  University of Kentucky
C                  Lexington, KY 40506
C
C     Input Data are in FDTDLOSS.IN, VSVL.IN and PUL.DAT
C
C     Output Data are in FDTDLOSS.OUT
C
C     NDZMAX is the Maximum Number of Discretization Sections for the Line.
C     NMAX is the Maximum Number of Line Conductors (Excluding the Reference).
C
C******************************************************************************
C
      PARAMETER (NDZMAX=300)
      PARAMETER (NMAX=10)
      PARAMETER (MAXVSPTS=10)
      PARAMETER (NDZ1MAX=NDZMAX+1)
      REAL VOLD(NMAX,NDZ1MAX),VNEW(NMAX,NDZ1MAX),S(NMAX),
     1IOLD(NMAX,NDZMAX),INEW(NMAX,NDZMAX),D(NMAX,NMAX),E(NMAX,NMAX),
     2LEN,L(NMAX,NMAX),C(NMAX,NMAX),A(NMAX,NMAX),B(NMAX,NMAX),
     3RS(NMAX,NMAX),RL(NMAX,NMAX),
     4VSVLT(MAXVSPTS),VSM(MAXVSPTS),VLM(MAXVSPTS),VS(NMAX),VL(NMAX),
     5VS1(NMAX),VL1(NMAX),RDC(NMAX),RFREQ(NMAX),AA(NMAX,NMAX),
     6BB(NMAX,NMAX),AI(10),ALPHA(10),FI(NMAX,NMAX),FIB(NMAX,NMAX),
     7FIL(NMAX,NMAX),PSI(NMAX,NDZMAX),
     8M1(NMAX,NDZMAX),M2(NMAX,NDZMAX),
     9M3(NMAX,NDZMAX),M4(NMAX,NDZMAX),M5(NMAX,NDZMAX),M6(NMAX,NDZMAX),
     *M7(NMAX,NDZMAX),M8(NMAX,NDZMAX),M9(NMAX,NDZMAX),M10(NMAX,NDZMAX)
      INTEGER IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      REAL IDENT(NMAX,NMAX)
      OPEN (UNIT=1,FILE='FDTDLOSS.IN',STATUS='OLD')
      OPEN (UNIT=2,FILE='PUL.DAT',STATUS='OLD')
      OPEN (UNIT=3,FILE='VSVL.IN',STATUS='OLD')
      OPEN (UNIT=4,FILE='FDTDLOSS.OUT',STATUS='UNKNOWN')
      PI=4.*ATAN(1.)
      AI(1)=0.79098180E-1
      ALPHA(1)=-0.11484427E-2
      AI(2)=0.11543423
      ALPHA(2)=-0.13818329E-1
      AI(3)=0.13435380
      ALPHA(3)=-0.54037596E-1
      AI(4)=0.21870422
      ALPHA(4)=-0.14216494
      AI(5)=0.98229667E-1
      ALPHA(5)=-0.30128437
      AI(6)=0.51360484
      ALPHA(6)=-0.56142185
      AI(7)=-0.20962898
      ALPHA(7)=-0.97117126
      AI(8)=0.11974447E1
      ALPHA(8)=-0.16338433E1
      AI(9)=0.11225491E-1
      ALPHA(9)=-0.28951329E1
      AI(10)=-0.74425255
      ALPHA(10)=-0.50410969E1
C
C**********INPUT DATA*****************************************
C
      READ(1,*) NP1
      N=NP1-1
      READ(1,*) LEN
      READ(1,*) NDZ
      READ(1,*) NDT
      READ(1,*) FINAL
      READ(1,*) NPT
      NST=NPT
      READ(1,*) NVS
      READ(1,*) NVL
      READ(1,*) NOUT
      DZ=LEN/FLOAT(NDZ)
      DT=FINAL/FLOAT(NDT)
      READ(1,*) R0
      READ(1,*) R0FREQ
      DO 100 I=1,N
      READ(1,*) K,RDC(K)
      READ(1,*) K,RFREQ(K)
100   CONTINUE
      WRITE(4,1000) NP1,LEN,NDZ,NDT,DZ,DT,FINAL,NPT,NVS,NVL,NOUT
1000  FORMAT('Total Number of Conductors (N+1)=',I2/
     1'Line Length=',1PE12.5/
     2'NDeltaZ=',I5/
     3'NDeltaT=',I5/
     4'DeltaZ=',1PE12.5/
     5'DeltaT=',1PE12.5/
     6'Final Solution Time=',1PE12.5/
     7'Print Interval=',I4/
     8'Conductor Having Source Voltage, VS(t), (at near end)=',I2/
     9'Conductor Having Load Voltage, VL(t), (at far end)=',I2/
     9'Output Conductor=',I2/)
      WRITE(4,200) R0,R0FREQ
200   FORMAT('Reference Conductor dc Resistance=', 1PE12.5/
     1'Reference Conductor Break Frequency=',1PE12.5/)
      IF(R0FREQ.EQ.0.) R0FREQ=1.
      WRITE(4,300)
300   FORMAT('I',8X,'RDC',8X,'BREAK FREQUENCY')
      DO 400 I=1,N
      WRITE(4,500) I,RDC(I),RFREQ(I)
500   FORMAT(I2,2X,1PE12.5,5X,1PE12.5)
      IF(RFREQ(I).EQ.0.) RFREQ(I)=1.
400   CONTINUE
      WRITE(4,501)
501   FORMAT(/)
      F=DT/DZ
      NDZ1=NDZ+1
      DO 1 I=1,N
      DO 1 J=I,N
      READ(1,*) M,K,RS(M,K),RL(M,K)
      WRITE(4,2000) M,K,RS(M,K),M,K,RL(M,K)
2000  FORMAT('RS(',I2,',',I2,')=',1PE12.5,5X,
     1'RL(',I2,',',I2,')=',1PE12.5)
      RS(K,M)=RS(M,K)
1     RL(K,M)=RL(M,K)
      WRITE(4,3000)
3000  FORMAT(1X)
      DO 2 I=1,N
      DO 2 J=I,N
      READ(2,*) M,K,L(M,K)
      WRITE(4,4000) M,K,L(M,K)
4000  FORMAT('L(',I2,',',I2,')=',1PE12.5)
2     L(K,M)=L(M,K)
      DO 3 I=1,N
      DO 3 J=I,N
      READ(2,*) M,K,C(M,K)
      WRITE(4,5000) M,K,C(M,K)
5000  FORMAT('C(',I2,',',I2,')=',1PE12.5)
3     C(K,M)=C(M,K)
      DO 600 I=1,N
      DO 601 J=1,N
601   AA(I,J)=R0
600   AA(I,I)=AA(I,I)+RDC(I)
      DO 602 I=1,N
      DO 603 J=1,N
603   BB(I,J)=R0/SQRT(PI*R0FREQ)
602   BB(I,I)=BB(I,I)+RDC(I)/SQRT(PI*RFREQ(I))
      WRITE(4,6000)
6000  FORMAT(/10X,'VS(t),VL(t) Specification')
      WRITE(4,7000)
7000  FORMAT(1X,'NO.',8X,'TIME',14X,'VS',14X,'VL')
      READ(3,*) NVSPTS
      DO 4 I=1,NVSPTS
      READ(3,*) VSVLT(I),VSM(I),VLM(I)
      WRITE(4,8000) I,VSVLT(I),VSM(I),VLM(I)
8000  FORMAT(I3,5X,3(1PE12.5,5X))
4     CONTINUE
      IF(NDZ.GT.NDZMAX) GO TO 999
      WRITE(4,9000)
9000  FORMAT(/4X,'TIME',13X,'VNE',13X,'VFE',13X,'VS',13X,'VL')
C
C*****INITIAL COMPUTATIONS*****************************************************
C
      DO 5 I=1,N
      DO 5 J=1,NDZ1
5     VOLD(I,J)=0.
      DO 6 I=1,N
      DO 6 J=1,NDZ
      PSI(I,J)=0.
      M1(I,J)=0.
      M2(I,J)=0.
      M3(I,J)=0.
      M4(I,J)=0.
      M5(I,J)=0.
      M6(I,J)=0.
      M7(I,J)=0.
      M8(I,J)=0.
      M9(I,J)=0.
      M10(I,J)=0.
6     IOLD(I,J)=0.
      DO 7 I=1,N
      DO 8 J=1,N
      SUM1=0.
      SUM2=0.
      DO 9 K=1,N
      SUM1=SUM1+RS(I,K)*C(K,J)/F
9     SUM2=SUM2+RL(I,K)*C(K,J)/F
      A(I,J)=SUM1
      D(I,J)=SUM1
      E(I,J)=SUM2
8     B(I,J)=SUM2
      A(I,I)=A(I,I)+1.0
      D(I,I)=D(I,I)-1.0
      E(I,I)=E(I,I)-1.0
7     B(I,I)=B(I,I)+1.0
      SUM1=0.
      DO 700 I=1,10
700   SUM1=SUM1+AI(I)
      Z00=SUM1
      DO 705 I=1,N
      DO 705 J=1,N
      FI(I,J)=L(I,J)/F+AA(I,J)*DZ/2.+BB(I,J)*Z00*DZ/SQRT(PI*DT)
705   L(I,J)=FI(I,J)-AA(I,J)*DZ
      CALL IDENTITY(IDENT,N,NMAX)
      CALL GAUSSJ(A,N,NMAX,IDENT,N,NMAX,IPIV,INDXR,INDXC)
      CALL IDENTITY(IDENT,N,NMAX)
      CALL GAUSSJ(B,N,NMAX,IDENT,N,NMAX,IPIV,INDXR,INDXC)
      CALL IDENTITY(IDENT,N,NMAX)
      CALL GAUSSJ(FI,N,NMAX,IDENT,N,NMAX,IPIV,INDXR,INDXC)
      CALL IDENTITY(IDENT,N,NMAX)
      CALL GAUSSJ(C,N,NMAX,IDENT,N,NMAX,IPIV,INDXR,INDXC)
      DO 706 I=1,N
      DO 706 J=1,N
      SUM1=0.
      SUM2=0.
      DO 707 K=1,N
      SUM1=SUM1+FI(I,K)*L(K,J)
707   SUM2=SUM2+FI(I,K)*BB(K,J)
      FIL(I,J)=SUM1
706   FIB(I,J)=SUM2*DZ/SQRT(PI*DT)
C
C**********ITERATE POSITION THEN TIME**************************
C
      TIME=-DT
      NDT1=NDT+1
      DO 10 IDT=1,NDT1
      TIME=TIME+DT
      TIMEVS1=TIME
      TIMEVS=TIME-DT
      CALL VSVL(VS1,VL1,N,NMAX,NVS,NVL,TIMEVS,VSM,VLM,
     1VSVLT,NVSPTS,MAXVSPTS)
      CALL VSVL(VS,VL,N,NMAX,NVS,NVL,TIMEVS1,VSM,VLM,
     1VSVLT,NVSPTS,MAXVSPTS)
      DO 11 I=1,N
      SUM=0.
      DO 12 J=1,N
12    SUM=SUM-2.*RS(I,J)*IOLD(J,1)+D(I,J)*VOLD(J,1)
11    S(I)=SUM+VS1(I)+VS(I)
      DO 13 I=1,N
      SUM=0.
      DO 14 J=1,N
14    SUM=SUM+A(I,J)*S(J)
13    VNEW(I,1)=SUM
      DO 15 I=2,NDZ
      DO 16 J=1,N
      SUM=0.
      DO 17 K=1,N
17    SUM=SUM+F*C(J,K)*(IOLD(K,I)-IOLD(K,I-1))
16    VNEW(J,I)=VOLD(J,I)-SUM
15    CONTINUE
      DO 18 I=1,N
      SUM=0.
      DO 19 J=1,N
19    SUM=SUM+2.*RL(I,J)*IOLD(J,NDZ)+E(I,J)*VOLD(J,NDZ1)
18    S(I)=SUM+VL1(I)+VL(I)
      DO 20 I=1,N
      SUM=0.
      DO 21 J=1,N
21    SUM=SUM+B(I,J)*S(J)
20    VNEW(I,NDZ1)=SUM
      IF(IDT.GT.1) GO TO 33
      WRITE(4,10000) TIME,VNEW(NOUT,1),VNEW(NOUT,NDZ1),VS(NVS),VL(NVL)
      NST=NST+1
      GO TO 30
33    CONTINUE
      NTEST=IDT-NST
      IF(NTEST.LT.0) GO TO 30
      WRITE(4,10000) TIME,VNEW(NOUT,1),VNEW(NOUT,NDZ1),VS(NVS),VL(NVL)
10000 FORMAT(5(1PE12.5,4X))
      NST=NST+NPT
30    CONTINUE
      DO 800 I=1,N
      DO 800 J=1,NDZ
800   PSI(I,J)=M1(I,J)+M2(I,J)+M3(I,J)+M4(I,J)+M5(I,J)+M6(I,J)
     1+M7(I,J)+M8(I,J)+M9(I,J)+M10(I,J)
      DO 22 I=1,N
      DO 23 J=1,NDZ
      SUM=0.
      DO 24 K=1,N
24    SUM=SUM+FIL(I,K)*IOLD(K,J)-FI(I,K)*(VNEW(K,J+1)-VNEW(K,J))
     1-FIB(I,K)*PSI(K,J)
23    INEW(I,J)=SUM
22    CONTINUE
      IF(IDT.EQ.1) GO TO 900
      DO 801 I=1,N
      DO 801 J=1,NDZ
      M10(I,J)=AI(10)*EXP(ALPHA(10))*(INEW(I,J)-IOLD(I,J))
     1+EXP(ALPHA(10))*M10(I,J)
      M9(I,J)=AI(9)*EXP(ALPHA(9))*(INEW(I,J)-IOLD(I,J))
     1+EXP(ALPHA(9))*M9(I,J)
      M8(I,J)=AI(8)*EXP(ALPHA(8))*(INEW(I,J)-IOLD(I,J))
     1+EXP(ALPHA(8))*M8(I,J)
      M7(I,J)=AI(7)*EXP(ALPHA(7))*(INEW(I,J)-IOLD(I,J))
     1+EXP(ALPHA(7))*M7(I,J)
      M6(I,J)=AI(6)*EXP(ALPHA(6))*(INEW(I,J)-IOLD(I,J))
     1+EXP(ALPHA(6))*M6(I,J)
      M5(I,J)=AI(5)*EXP(ALPHA(5))*(INEW(I,J)-IOLD(I,J))
     1+EXP(ALPHA(5))*M5(I,J)
      M4(I,J)=AI(4)*EXP(ALPHA(4))*(INEW(I,J)-IOLD(I,J))
     1+EXP(ALPHA(4))*M4(I,J)
      M3(I,J)=AI(3)*EXP(ALPHA(3))*(INEW(I,J)-IOLD(I,J))
     1+EXP(ALPHA(3))*M3(I,J)
      M2(I,J)=AI(2)*EXP(ALPHA(2))*(INEW(I,J)-IOLD(I,J))
     1+EXP(ALPHA(2))*M2(I,J)
      M1(I,J)=AI(1)*EXP(ALPHA(1))*(INEW(I,J)-IOLD(I,J))
     1+EXP(ALPHA(1))*M1(I,J)
801   CONTINUE
900   CONTINUE
      DO 31 I=1,N
      DO 31 J=1,NDZ
      IOLD(I,J)=INEW(I,J)
31    VOLD(I,J)=VNEW(I,J)
      DO 32 I=1,N
32    VOLD(I,NDZ1)=VNEW(I,NDZ1)
10    CONTINUE
      GO TO 998
999   CONTINUE
      WRITE(4,20000)
20000 FORMAT(1X,///,10X,'NDZ EXCEEDS THE MAXIMUM',///)
998   CONTINUE
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
      STOP
      END
C
C*****SUBROUTINE VSVL TO FILL THE SOURCE AND LOAD VECTORS**********************
C
      SUBROUTINE VSVL(VS,VL,N,NMAX,NVS,NVL,TIME,VSM,VLM,
     1VSVLT,NVSPTS,MAXVSPTS)
      REAL VS(NMAX),VL(NMAX),VSVLT(MAXVSPTS),VSM(MAXVSPTS),VLM(MAXVSPTS)
      VSTEMP=0.
      VLTEMP=0.
      NVSPTSM1=NVSPTS-1
      DO 1 K=1,NVSPTSM1
      VS2=VSM(K+1)
      VS1=VSM(K)
      VL2=VLM(K+1)
      VL1=VLM(K)
      T2=VSVLT(K+1)
      T1=VSVLT(K)
      SLOPES=(VS2-VS1)/(T2-T1)
      SLOPEL=(VL2-VL1)/(T2-T1)
      IF(TIME.GE.T1.AND.TIME.LE.T2) VSTEMP=SLOPES*(TIME-T1)+VS1
      IF(TIME.GE.T1.AND.TIME.LE.T2) VLTEMP=SLOPEL*(TIME-T1)+VL1
1     CONTINUE
      DO 2 I=1,N
      VS(I)=0.
2     VL(I)=0.
      VS(NVS)=VSTEMP
      VL(NVL)=VLTEMP
      RETURN
      END
C
C*****SUBROUTINE IDENTITY TO PRODUCE THE NxN IDENTITY MATRIX*******************
C
      SUBROUTINE IDENTITY(IDENT,N,NMAX)
      REAL IDENT(NMAX,NMAX)
      DO 1 I=1,N
      DO 2 J=1,N
2     IDENT(I,J)=0.
1     IDENT(I,I)=1.
      RETURN
      END
C
C*****SUBROUTINE GAUSSJ TO INVERT A REAL NxN MATRIX****************************
C     Copyright (C) Numerical Recipes Software
C     Permission is granted for use only with the programs contained
C     within this book.
C
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP,IPIV,INDXR,INDXC)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NP),INDXR(NP),INDXC(NP)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) PAUSE 'Singular matrix.'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
