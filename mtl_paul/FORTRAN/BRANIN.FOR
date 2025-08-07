C***************PROGRAM BRANIN.FOR*********************************************
C
C     To Implement Branin's Time-Domain Solution Method Extended to a
C     Lossless (N+1)-Conductor Transmission Line in a Homogeneous Medium
C     Having Resistive Loads.
C
C             Written by
C                  Clayton R. Paul
C                  Department of Electrical Engineering
C                  University of Kentucky
C                  Lexington, KY 40506
C
C     Input Data are in BRANIN.IN, VSVL.IN, PUL.DAT
C
C     Output is in BRANIN.OUT
C
C     NDTMAX is the Maximum Number of Time Points Per One-Way Delay
C     MAXVSPTS is the Maximum Number of Time Points Specifying VS(t) and VL(t)
C     NCONDMAX is the Maximum Number of Conductors (Excluding Reference)
C
C*****SET ARRAY DIMENSIONS AND OPEN FILES**************************************
C
      PARAMETER (NDTMAX=100)
      PARAMETER (MAXVSPTS=10)
      PARAMETER (NCONDMAX=2)
      REAL LEN,EIOLD(NCONDMAX,NDTMAX),EROLD(NCONDMAX,NDTMAX),
     1EINEW(NCONDMAX,NDTMAX),ERNEW(NCONDMAX,NDTMAX),
     2IDENT(NCONDMAX,NCONDMAX),TEMP1(NCONDMAX,NCONDMAX),
     3TEMP2(NCONDMAX,NCONDMAX),RS(NCONDMAX,NCONDMAX),
     4RL(NCONDMAX,NCONDMAX),L(NCONDMAX,NCONDMAX),
     5GS(NCONDMAX,NCONDMAX),GL(NCONDMAX,NCONDMAX),
     6VSLT(MAXVSPTS),VSM(MAXVSPTS),VLM(MAXVSPTS)
      INTEGER IPIV(NCONDMAX),INDXR(NCONDMAX),INDXC(NCONDMAX)
      OPEN (UNIT=1,FILE='BRANIN.IN',STATUS='OLD')
      OPEN (UNIT=2,FILE='PUL.DAT',STATUS='OLD')
      OPEN (UNIT=3,FILE='VSVL.IN',STATUS='OLD')
      OPEN (UNIT=4,FILE='BRANIN.OUT',STATUS='UNKNOWN')
C
C*****READ AND PRINT INPUT DATA************************************************
C
      READ(1,*) NP1
      N=NP1-1
      READ(1,*) LEN
      READ(1,*) ER
      READ(1,*) NTOT
      READ(1,*) NDT
      READ(1,*) NVS
      READ(1,*) NVL
      READ(1,*) NOUT
      WRITE(4,1000) NP1,LEN,ER,NTOT,NDT,NVS,NVL,NOUT
1000  FORMAT('Total Number of Conductors (N+1)=',I2/
     1'Line Length=',1PE12.5/'Dielectric Constant of Medium=',1PE12.5/
     2'Number of One-Way Time Delays=',I3/
     3'Number of Solution Points per One-Way Delay=',I3/
     4'Conductor Having Source Voltage, VS(t), (at near end)=',I2/
     5'Conductor Having Load Voltage, VL(t), (at far end)=',I2/
     6'Output Conductor =',I2/)
      V=2.997925E8/SQRT(ER)
      DO 1 I=1,N
      DO 1 J=I,N
      READ(1,*) M,K,RS(M,K),RL(M,K)
      WRITE(4,2000) M,K,RS(M,K),M,K,RL(M,K)
2000  FORMAT('RS(',I2,',',I2,')=',1PE12.5,5X,
     1'RL(',I2,',',I2,')=',1PE12.5)
      RS(K,M)=RS(M,K)
      RL(K,M)=RL(M,K)
1     CONTINUE
      WRITE(4,2001)
2001  FORMAT(1X)
      DO 2 I=1,N
      DO 2 J=I,N
      READ(2,*) M,K,L(M,K)
      WRITE(4,3000) M,K,L(M,K)
3000  FORMAT('L(',I2,',',I2,')=',1PE12.5)
      L(M,K)=V*L(M,K)
      L(K,M)=L(M,K)
2     CONTINUE
      WRITE(4,4000)
4000  FORMAT(/10X,'VS(t),VL(t) Specification')
      WRITE(4,5000)
5000  FORMAT(1X,'NO.',8X,'TIME',14X,'VS',14X,'VL')
      READ(3,*) NVSPTS
      DO 3 I=1,NVSPTS
      READ(3,*) VSLT(I),VSM(I),VLM(I)
      WRITE(4,6000) I,VSLT(I),VSM(I),VLM(I)
6000  FORMAT(I3,5X,3(1PE12.5,5X))
3     CONTINUE
      WRITE(4,7000)
7000  FORMAT(/4X,'TIME',13X,'VNE',13X,'VFE',13X,'VS',13X,'VL')
C
C***********INITIAL COMPUTATIONS******************************
C
      T=LEN/V
      XNDT=FLOAT(NDT)
      DT=T/XNDT
      DO 4 I=1,N
      DO 4 J=1,NDT
      EIOLD(I,J)=0.
4     EROLD(I,J)=0.
      DO 5 I=1,N
      DO 6 J=1,N
6     IDENT(I,J)=0.
5     IDENT(I,I)=1.
C*****COMPUTE GAMMAS, GAMMAL, MS, ML
      DO 7 I=1,N
      DO 7 J=1,N
      TEMP1(I,J)=RS(I,J)+L(I,J)
7     TEMP2(I,J)=RL(I,J)+L(I,J)
      CALL GAUSSJ(TEMP1,N,NCONDMAX,IDENT,N,NCONDMAX,IPIV,INDXR,INDXC)
      CALL GAUSSJ(TEMP2,N,NCONDMAX,IDENT,N,NCONDMAX,IPIV,INDXR,INDXC)
      DO 8 I=1,N
      DO 8 J=1,N
      SUM1=0.
      SUM2=0.
      DO 9 K=1,N
      SUM1=SUM1+(RS(I,K)-L(I,K))*TEMP1(K,J)
9     SUM2=SUM2+(RL(I,K)-L(I,K))*TEMP2(K,J)
      GS(I,J)=SUM1
8     GL(I,J)=SUM2
      DO 10 I=1,N
      DO 10 J=1,N
      SUM1=0.
      SUM2=0.
      DO 11 K=1,N
      SUM1=SUM1+L(I,K)*TEMP1(K,J)
11    SUM2=SUM2+L(I,K)*TEMP2(K,J)
      RS(I,J)=2.*SUM1
10    RL(I,J)=2.*SUM2
      TIME=0.
      INDEX=0
C
C*****START ITERATION**********************************************************
C
      DO 12 I=1,NTOT
      DO 13 J=1,NDT
C*************DEFINE VS(t) AND VL(t)***
      VS=0.
      VL=0.
      NVSPTSM1=NVSPTS-1
      DO 14 K=1,NVSPTSM1
      VS2=VSM(K+1)
      VS1=VSM(K)
      VL2=VLM(K+1)
      VL1=VLM(K)
      T2=VSLT(K+1)
      T1=VSLT(K)
      SLOPES=(VS2-VS1)/(T2-T1)
      SLOPEL=(VL2-VL1)/(T2-T1)
      IF(TIME.GE.T1.AND.TIME.LE.T2) VS=SLOPES*(TIME-T1)+VS1
      IF(TIME.GE.T1.AND.TIME.LE.T2) VL=SLOPEL*(TIME-T1)+VL1
14    CONTINUE
C***************************************
      DO 15 II=1,N
      SUM1=0.
      SUM2=0.
      DO 16 KK=1,N
      SUM1=SUM1+GS(II,KK)*EROLD(KK,J)
16    SUM2=SUM2+GL(II,KK)*EIOLD(KK,J)
      EINEW(II,J)=SUM1+RS(II,NVS)*VS
15    ERNEW(II,J)=SUM2+RL(II,NVL)*VL
      VNE=(EINEW(NOUT,J)+EROLD(NOUT,J))/2.
      VFE=(ERNEW(NOUT,J)+EIOLD(NOUT,J))/2.
      WRITE(4,8000) TIME,VNE,VFE,VS,VL
8000  FORMAT(5(1PE12.5,4X))
      INDEX=INDEX+1
      TIME=FLOAT(INDEX)*DT
13    CONTINUE
      DO 17 KK=1,N
      DO 18 JJ=1,NDT
      EIOLD(KK,JJ)=EINEW(KK,JJ)
      EROLD(KK,JJ)=ERNEW(KK,JJ)
18    CONTINUE
17    CONTINUE
12    CONTINUE
      CLOSE (1)
      CLOSE (2)
      CLOSE (3)
      CLOSE (4)
      STOP
      END
C********SUBROUTINE GAUSSJ FOR INVERTING REAL MATRICES*************************
C     Copyright (C) 1987-1992 Numerical Recipes Software
C     Permission is granted for use only with the programs contained
C     within this book
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