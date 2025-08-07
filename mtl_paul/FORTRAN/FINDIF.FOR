C*************PROGRAM FINDIF.FOR***********************************************
C
C     To Implement the Finite Difference-Time Domain Method for
C     Determining the Time-Domain Response of a (N+1)-Conductor
C     Transmission Line Having Resistive Loads.
C
C             Written by:
C                  Clayton R. Paul
C                  Department of Electrical Engineering
C                  University of Kentucky
C                  Lexington, KY 40506
C
C     Input Data are in FINDIF.IN, VSVL.IN and PUL.DAT
C
C     Output Data are in FINDIF.OUT
C
C     NDZMAX is the Maximum Number of Discretization Sections for the Line.
C     NMAX is the Maximum Number of Line Conductors (Excluding the Reference).
C
C******************************************************************************
C
      PARAMETER (NDZMAX=200)
      PARAMETER (NMAX=2)
      PARAMETER (MAXVSPTS=10)
      PARAMETER (NDZ1MAX=NDZMAX+1)
      REAL VOLD(NMAX,NDZ1MAX),VNEW(NMAX,NDZ1MAX),S(NMAX),
     1IOLD(NMAX,NDZMAX),INEW(NMAX,NDZMAX),D(NMAX,NMAX),E(NMAX,NMAX),
     2LEN,L(NMAX,NMAX),C(NMAX,NMAX),A(NMAX,NMAX),B(NMAX,NMAX),
     3RS(NMAX,NMAX),RL(NMAX,NMAX),
     4VSVLT(MAXVSPTS),VSM(MAXVSPTS),VLM(MAXVSPTS),VS(NMAX),VL(NMAX),
     5VS1(NMAX),VL1(NMAX)
      INTEGER IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      REAL IDENT(NMAX,NMAX)
      OPEN (UNIT=1,FILE='FINDIF.IN',STATUS='OLD')
      OPEN (UNIT=2,FILE='PUL.DAT',STATUS='OLD')
      OPEN (UNIT=3,FILE='VSVL.IN',STATUS='OLD')
      OPEN (UNIT=4,FILE='FINDIF.OUT',STATUS='UNKNOWN')
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
      CALL IDENTITY(IDENT,N,NMAX)
      CALL GAUSSJ(A,N,NMAX,IDENT,N,NMAX,IPIV,INDXR,INDXC)
      CALL IDENTITY(IDENT,N,NMAX)
      CALL GAUSSJ(B,N,NMAX,IDENT,N,NMAX,IPIV,INDXR,INDXC)
      CALL IDENTITY(IDENT,N,NMAX)
      CALL GAUSSJ(L,N,NMAX,IDENT,N,NMAX,IPIV,INDXR,INDXC)
      CALL IDENTITY(IDENT,N,NMAX)
      CALL GAUSSJ(C,N,NMAX,IDENT,N,NMAX,IPIV,INDXR,INDXC)
C
C**********ITERATE POSITION THEN TIME**************************
C
      NDT1=NDT+1
      DO 10 IDT=1,NDT1
      TIME=FLOAT(IDT-1)*DT
      TIMEVS1=TIME
      TIMEVS=FLOAT(IDT-2)*DT
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
      DO 22 I=1,NDZ
      DO 23 J=1,N
      SUM=0.
      DO 24 K=1,N
24    SUM=SUM+F*L(J,K)*(VNEW(K,I+1)-VNEW(K,I))
23    INEW(J,I)=IOLD(J,I)-SUM
22    CONTINUE
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
C     Copyright (C) 1987-1992 Numerical Recipes Software
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