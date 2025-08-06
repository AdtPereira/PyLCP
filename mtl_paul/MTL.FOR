C*******************************************************************************
C             Program MTL.FOR
C             (FORTRAN)
C
C             Written by
C                  Clayton R. Paul
C                  Department of Electrical Engineering
C                  University of Kentucky
C                  Lexington, KY 40506
C
C     A digital computer program to compute the frequency response of a
C              (N+1)-conductor transmission line.
C
C        Input Data:
C             MTL.IN
C                  Contains (1) total number of conductors (N+1)
C                           (2) total line length
C                           (3) conductor number for which the
C                               terminal voltages are desired
C                           (4) dc resistances of conductors and break
C                               frequencies where dc resistance changes
C                               to skin effect (square root of frequency)
C                           (5) termination network characterizations as
C                               Generalized Thevenin Equivalents
C
C             FREQ.IN
C                  Contains solution frequencies sequentially
C
C             PUL.DAT
C                  Contains entries in per-unit-length inductance and
C                  capacitance matrices stored sequentially by rows
C                  for the upper diagonal
C
C     The Maximum Number of Conductors (Excluding the Reference) is MZ
C
C******************************************************************************
C
      PARAMETER (MZ= 2)
      PARAMETER (MZWK=2*MZ*(MZ+1))
      PARAMETER (NTMZ=2*MZ)
      REAL LEN,L(MZ,MZ),C(MZ,MZ),RDC(MZ),RFREQ(MZ),WK(MZWK),WA(NTMZ)
      COMPLEX XJ,ONEC,ZEROC,SUMC,ZREF,VS(MZ),VL(MZ),ZS(MZ,MZ),ZL(MZ,MZ),
     1Z(MZ,MZ),Y(MZ,MZ),SCR(MZ,MZ),T(MZ,MZ),GAM(MZ),
     2A(NTMZ,NTMZ),B(NTMZ),EP(MZ),EN(MZ),
     3ZST(MZ,MZ),ZLT(MZ,MZ),VNE,VFE
      OPEN (UNIT=1,FILE='MTL.IN',STATUS='OLD')
      OPEN (UNIT=2,FILE='PUL.DAT',STATUS='OLD')
      OPEN (UNIT=3,FILE='FREQ.IN',STATUS='OLD')
      OPEN (UNIT=4,FILE='MTL.OUT',STATUS='UNKNOWN')
      PI=4.0*ATAN(1.0)
      ONEC=CMPLX(1.E0,0.E0)
      ZEROC=CMPLX(0.E0,0.E0)
      XJ=CMPLX(0.E0,1.E0)
C
C*****READ AND PRINT INPUT DATA************************************************
C
      READ(1,*) NP1
      N=NP1-1
      READ(1,*) LEN
      READ(1,*) NOUT
      WRITE(4,100) NP1,LEN,NOUT
100   FORMAT('TOTAL NUMBER OF CONDUCTORS=',I2,/,'LINE LENGTH (METERS)='
     1,1PE10.3/,'OUTPUT CONDUCTOR=',I2/)
      NTWO=2*N
      READ(1,*) R0
      READ(1,*) R0FREQ
      WRITE(4,101) R0,R0FREQ
101   FORMAT('DC RESISTANCE OF REFERENCE CONDUCTOR (OHMS/M)=',1PE10.3/
     1'BREAK FREQUENCY (HZ)=',1PE10.3,/)
      WRITE(4,102)
102   FORMAT(' I',6X,'RDC',12X,'BREAK FREQUENCY')
      DO 1 I=1,N
      READ(1,*) K,RDC(K)
      READ(1,*) K,RFREQ(K)
      WRITE(4,103) K,RDC(K),RFREQ(K)
103   FORMAT(I2,2X,1PE10.3,10X,1PE10.3)
1     CONTINUE
      WRITE(4,104)
104   FORMAT(/,'TERMINAL NETWORK CHARACTERIZATION',/,
     1' I J',4X,'VSR',6X,'VSI',6X,'ZSR',6X,'ZSI',
     26X,'VLR',6X,'VLI',6X,'ZLR',6X,'ZLI')
      DO 2 I=1,N
      READ(1,*) K,VSR,VSI
2     VS(K)=CMPLX(VSR,VSI)
      DO 3 I=1,N
      DO 3 J=I,N
      READ(1,*) K,M,ZSR,ZSI
      ZS(K,M)=CMPLX(ZSR,ZSI)
3     ZS(M,K)=ZS(K,M)
      DO 4 I=1,N
      READ(1,*) K,VLR,VLI
4     VL(K)=CMPLX(VLR,VLI)
      DO 5 I=1,N
      DO 5 J=I,N
      READ(1,*) K,M,ZLR,ZLI
      ZL(K,M)=CMPLX(ZLR,ZLI)
5     ZL(M,K)=ZL(K,M)
      DO 6 I=1,N
      WRITE(4,105) I,I,VS(I),ZS(I,I),VL(I),ZL(I,I)
105   FORMAT(I2,I2,8(1X,1PE8.1))
6     CONTINUE
      DO 7 I=1,N
      DO 7 J=I+1,N
      WRITE(4,106) I,J,ZS(I,J),ZL(I,J)
106   FORMAT(I2,I2,18X,2(1X,1PE8.1),18X,2(1X,1PE8.1))
7     CONTINUE
C
C*****READ AND PRINT PER-UNIT-LENGTH PARAMETERS
C
      WRITE(4,107)
107   FORMAT(/,'PER-UNIT-LENGTH PARAMETERS')
      DO 8 I=1,N
      DO 8 J=I,N
      READ(2,*) M,K,L(M,K)
      L(J,I)=L(I,J)
      WRITE(4,108) I,J,L(I,J)
108   FORMAT('L(',I2,',',I2,')=',1PE12.5)
8     CONTINUE
      DO 9 I=1,N
      DO 9 J=I,N
      READ(2,*) M,K,C(M,K)
      C(J,I)=C(I,J)
      WRITE(4,109) I,J,C(I,J)
109   FORMAT('C(',I2,',',I2,')=',1PE12.5)
9     CONTINUE
      WRITE(4,110) NOUT
110   FORMAT(/,'OUTPUT CONDUCTOR=',I2)
      WRITE(4,111)
111   FORMAT(/,1X,'FREQUENCY',4X,'VNEMAG',6X,'VNEANG',6X,'VFEMAG',
     16X,'VFEANG')
C
C*****READ AND PRINT FREQUENCY*************************************************
C
1000  CONTINUE
      READ(3,*,END=2000) F
      OMEGA=2.*PI*F
C
C*****COMPUTE Z AND Y
C
      DO 10 I=1,N
      DO 10 J=1,N
      Z(I,J)=XJ*OMEGA*L(I,J)
10    Y(I,J)=XJ*OMEGA*C(I,J)
      DO 11 I=1,N
      IF (RFREQ(I).EQ.0.) RFREQ(I)=1.
      IF(F.GT.RFREQ(I)) GO TO 12
      Z(I,I)=Z(I,I)+RDC(I)*(1.+XJ*F/RFREQ(I))
      GO TO 13
12    Z(I,I)=Z(I,I)+RDC(I)*SQRT(F/RFREQ(I))*(1.+XJ)
13    CONTINUE
11    CONTINUE
      IF(R0FREQ.EQ.0.) R0FREQ=1.
      IF(F.GT.R0FREQ) GO TO 14
      ZREF=R0*(1.+XJ*F/R0FREQ)
      GO TO 15
14    ZREF=R0*SQRT(F/R0FREQ)*(1.+XJ)
15    CONTINUE
      DO 16 I=1,N
      DO 16 J=1,N
16    Z(I,J)=Z(I,J)+ZREF
C
C*****COMPUTE EIGENVALUES (GAM) AND EIGENVECTORS (T) OF Y*Z
C
      CALL MATMULT(SCR,Y,Z,N,MZ)
      CALL EIGCC(SCR,N,MZ,2,GAM,T,MZ,WK,IER)
C
C*****COMPUTE MATRIX EXPONENTIALS
C
      DO 17 I=1,N
      GAM(I)=CSQRT(GAM(I))
      EP(I)=CEXP(GAM(I)*LEN)
      EN(I)=CEXP(-GAM(I)*LEN)
17    GAM(I)=ONEC/GAM(I)
C
C*****COMPUTE Z*T/GAM AND ZS*T AND ZL*T
C
      CALL MATMULT(SCR,Z,T,N,MZ)
      CALL VECMULT(Z,SCR,GAM,N,MZ)
      CALL MATMULT(ZST,ZS,T,N,MZ)
      CALL MATMULT(ZLT,ZL,T,N,MZ)
C
C*****FILL EQUATIONS A*X=B FOR DETERMINING UNDETERMINED CONSTANTS Im+ AND Im-
C
      DO 18 I=1,N
      DO 19 J=1,N
      A(I,J)=Z(I,J)+ZST(I,J)
      A(I,J+N)=Z(I,J)-ZST(I,J)
      A(I+N,J)=(Z(I,J)-ZLT(I,J))*EN(J)
19    A(I+N,J+N)=(Z(I,J)+ZLT(I,J))*EP(J)
      B(I)=VS(I)
18    B(I+N)=VL(I)
      CALL LEQT1C(A,NTWO,NTMZ,B,1,NTMZ,0,WA,IER)
C
C*****FORM V(0) AND V(LEN) TO DETERMINE NEAR END AND FAR END SOLUTIONS FOR NCOND
C
      SUMC=ZEROC
      DO 20 J=1,N
20    SUMC=SUMC+Z(NOUT,J)*(B(J)+B(J+N))
      VNE=SUMC
      VNEM=CABS(VNE)
      VNER=REAL(VNE)
      VNEI=AIMAG(VNE)
      VNEA=ATAN2(VNEI,VNER)*180./PI
      SUMC=ZEROC
      DO 21 J=1,N
21    SUMC=SUMC+Z(NOUT,J)*(EN(J)*B(J)+EP(J)*B(J+N))
      VFE=SUMC
      VFEM=CABS(VFE)
      VFER=REAL(VFE)
      VFEI=AIMAG(VFE)
      VFEA=ATAN2(VFEI,VFER)*180./PI
      WRITE(4,112) F,VNEM,VNEA,VFEM,VFEA
112   FORMAT(5(1PE10.3,2X))
      GO TO 1000
2000  CONTINUE
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
      STOP
      END
C
C*****SUBROUTINE TO MULTIPLY TWO NxN MATRICES AS A=BxC
C
      SUBROUTINE MATMULT(A,B,C,N,MZ)
      COMPLEX A(MZ,MZ),B(MZ,MZ),C(MZ,MZ),SUMC,ZEROC
      ZEROC=CMPLX(0.,0.)
      DO 1 I=1,N
      DO 1 J=1,N
      SUMC=ZEROC
      DO 2 K=1,N
2     SUMC=SUMC+B(I,K)*C(K,J)
1     A(I,J)=SUMC
      RETURN
      END
C
C*****SUBROUTINE TO MULTIPLY NxN MATRIX B AND Nx1 VECTOR V AS A=B*V
C
      SUBROUTINE VECMULT(A,B,V,N,MZ)
      COMPLEX A(MZ,MZ),B(MZ,MZ),V(MZ)
      DO 1 I=1,N
      DO 1 J=1,N
1     A(I,J)=B(I,J)*V(J)
      RETURN
      END
C
C*****IMSL SUBROUTINES*********************************************************
C
C     SUBROUTINE EBALAC (AR,AI,N,IA,K,L,D)
C
C-EBALAC--------S/D-----LIBRARY 1---------------------------------------
C
C   FUNCTION            - BALANCES A COMPLEX GENERAL MATRIX AND ISOLATES
C                           EIGENVALUES WHENEVER POSSIBLE.
C   USAGE               - CALL EBALAC (AR,AI,N,IA,K,L,D)
C   PARAMETERS   AR     - INPUT/OUTPUT MATRICES OF DIMENSION N BY N.
C                AI         ON INPUT, AR AND AND AI CONTAIN THE REAL
C                           AND IMAGINARY PARTS, RESPECTIVELY, OF
C                           THE COMPLEX MATRIX OF ORDER N TO BE
C                           BALANCED. ON OUTPUT, AR AND AI CONTAIN THE
C                           REAL AND IMAGINARY PARTS OF THE
C                           TRANSFORMED MATRIX.
C                N      - INPUT VARIABLE CONTAINING THE ORDER
C                           OF THE MATRIX A = (AR,AI) TO BE BALANCED.
C                IA     - INPUT VARIABLE CONTAINING THE ROW DIMENSION OF
C                           AR AND AI IN THE CALLING PROGRAM.
C                K      - OUTPUT INTEGERS CONTAINING THE BOUNDARY
C                L          INDICES FOR THE BALANCED MATRIX A = (AR,AI)
C                           SUCH THAT
C                             AR(I,J) = 0. AND AI(I,J) = 0. IF
C                             (1) I IS GREATER THAN J AND
C                             (2) J = 1,...,K-1 OR
C                                 I = L+1,...,N
C                D      - OUTPUT VECTOR OF LENGTH N CONTAINING
C                           INFORMATION DETERMINING THE PERMUTATIONS
C                           USED AND THE SCALING FACTORS.
C   PRECISION           - SINGLE/DOUBLE
C   LANGUAGE            - FORTRAN
C-----------------------------------------------------------------------
C   LATEST REVISION     - MARCH 9, 1977
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C     4/19/79     MODIFICATIONS WERE MADE TO PORTIONS OF THE
C              IMSL ROUTINES REGARDING HEXADECIMAL DESIG-
C              NATIONS WITHIN DATA STATEMENTS; THIS ACTION WAS
C              TAKEN IN RESPONSE TO FATAL ERRORS ENCOUNTERED IN
C              PARTS OF THE LIBRARY PACKAGE DURING THE INITIAL
C              COMPILATION.
C
C     NOTE   THE ROUTINES AFFECTED ARE   EIGCC,ELRH1C,ELRH2C;
C            THE HEXADECIMAL Z31C00000 WAS CHANGED TO 1.0E-6
C            IN EACH ROUTINE.
C
C            ALL CHANGES WERE MADE VIA TEXT EDITOR AND ARE EN-
C            CLOSED WITHIN PER CENT (%) SIGNS FOR EASE IN IDEN-
C            IFICATION.
C
C                                                    ALSD
C
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE EBALAC (AR,AI,N,IA,K,L,D)
C
      DIMENSION          AR(IA,1),AI(IA,1),D(N)
      LOGICAL            NOCONV
C                                  RADIX IS A MACHINE DEPENDENT
C                                    PARAMETER SPECIFYING THE BASE OF
C                                    THE MACHINE FLOATING POINT REPRE-
C                                    SENTATION
C*    DOUBLE PRECISION   AR,AI,D,RADIX,ZERO,ONE,PT95,B2,F,C,G,R,S
C*    DOUBLE PRECISION   RRADIX,RB2
C*    DATA               RADIX/16.0D0/
C2    DATA               ZERO,ONE,PT95/0.0D0,1.0D0,0.95D0/
      DATA               RADIX/16.0/
      DATA               ZERO,ONE,PT95/0.0,1.0,0.95/
      B2 = RADIX*RADIX
      RRADIX = ONE/RADIX
      RB2 = RRADIX*RRADIX
      K = 1
      L = N
      GO TO 30
C                                  IN-LINE PROCEDURE FOR ROW AND COLUMN
C                                    EXCHANGE
    5 D(M) = J
      IF (J .EQ. M) GO TO 20
      DO 10 I = 1,L
         F = AR(I,J)
         AR(I,J) = AR(I,M)
         AR(I,M) = F
         F = AI(I,J)
         AI(I,J) = AI(I,M)
         AI(I,M) = F
   10 CONTINUE
      DO 15 I=K,N
         F = AR(J,I)
         AR(J,I) = AR(M,I)
         AR(M,I) = F
         F = AI(J,I)
         AI(J,I) = AI(M,I)
         AI(M,I) = F
   15 CONTINUE
   20 GO TO (25,45), IEXC
C                                  SEARCH FOR ROWS ISOLATING AN
C                                    EIGENVALUE AND PUSH THEM DOWN
   25 IF (L .EQ. 1) GO TO 115
      L = L-1
C                                  DO J=L,1,-1
   30 L1 = L+1
      DO 40 JJ = 1,L
         J = L1-JJ
         DO 35 I = 1,L
            IF (I .EQ. J) GO TO 35
            IF (AR(J,I) .NE. ZERO .OR. AI(J,I) .NE. ZERO) GO TO 40
   35    CONTINUE
         M = L
         IEXC = 1
         GO TO 5
   40 CONTINUE
      GO TO 50
C                                  SEARCH FOR COLUMNS ISOLATING AN
C                                    EIGENVALUE AND PUSH THEM LEFT
   45 K = K+1
   50 DO 60 J = K,L
         DO 55 I = K,L
            IF (I .EQ. J) GO TO 55
            IF (AR(I,J) .NE. ZERO .OR. AI(I,J) .NE. ZERO) GO TO 60
   55    CONTINUE
         M = K
         IEXC = 2
         GO TO 5
   60 CONTINUE
C                                  BALANCE THE SUBMATRIX IN ROWS
C                                    K TO L
      DO 65 I = K,L
         D(I) = ONE
   65 CONTINUE
C                                  ITERATIVE LOOP FOR NORM REDUCTION
   70 NOCONV = .FALSE.
      DO 110 I = K,L
         C = ZERO
         R = ZERO
         DO 75 J = K,L
            IF (J .EQ. I) GO TO 75
C1          C = C+DABS(AR(J,I))+DABS(AI(J,I))
            C = C+ABS(AR(J,I))+ABS(AI(J,I))
C1          R = R+DABS(AR(I,J))+DABS(AI(I,J))
            R = R+ABS(AR(I,J))+ABS(AI(I,J))
   75    CONTINUE
         G = R*RRADIX
         F = ONE
         S = C+R
   80    IF (C .GE. G) GO TO 85
         F = F*RADIX
         C = C*B2
         GO TO 80
   85    G = R*RADIX
   90    IF (C .LT. G) GO TO 95
         F = F*RRADIX
         C = C*RB2
         GO TO 90
C                                  BALANCE
   95    IF ((C+R)/F .GE. PT95*S) GO TO 110
         G = ONE/F
         D(I) = D(I)*F
         NOCONV = .TRUE.
         DO 100 J = K,N
            AR(I,J) = AR(I,J)*G
            AI(I,J) = AI(I,J)*G
  100    CONTINUE
         DO 105 J = 1,L
            AR(J,I) = AR(J,I)*F
            AI(J,I) = AI(J,I)*F
  105    CONTINUE
  110 CONTINUE
      IF (NOCONV) GO TO 70
  115 RETURN
      END
C     SUBROUTINE EBBCKC (ZR,ZI,N,IZ,K,L,M,D)
C
C-EBBCKC--------S/D-----LIBRARY 1---------------------------------------
C
C   FUNCTION            - BACKTRANSFORM THE EIGENVECTORS OF A BALANCED
C                           COMPLEX GENERAL MATRIX.
C   USAGE               - CALL EBBCKC (ZR,ZI,N,IZ,K,L,M,D)
C   PARAMETERS   ZR     - INPUT/OUTPUT MATRICES OF DIMENSION N BY M.
C                ZI         ON INPUT, THE FIRST M COLUMNS OF ZR AND
C                           ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C                           RESPECTIVELY, OF THE EIGENVECTORS TO BE
C                           BACK TRANSFORMED. ON OUTPUT, THESE M
C                           COLUMNS CONTAIN THE REAL AND IMAGINARY
C                           PARTS OF THE TRANSFORMED EIGENVECTORS.
C                N      - INPUT SCALAR CONTAINING THE NUMBER OF
C                           ROWS IN THE MATRIX Z = (ZR,ZI). N MUST
C                           NOT BE GREATER THAN IZ.
C                IZ     - INPUT SCALAR CONTAINING THE ROW DIMENSION
C                           OF MATRICES ZR AND ZI IN THE CALLING
C                           PROGRAM.
C                K      - INPUT SCALARS CONTAINING THE BOUNDARY
C                L          INDICES FOR THE BALANCED MATRIX. K AND L
C                           ARE TWO OUTPUT PARAMETERS FROM IMSL ROUTINE
C                           EBALAC.
C                M      - INPUT SCALAR CONTAINING THE NUMBER OF COLUMNS
C                           OF Z = (ZR,ZI) TO BE BACK TRANSFORMED.
C                D      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           DETAILS OF THE TRANSFORMATIONS PRODUCED
C                           BY IMSL ROUTINE EBALAC.
C   PRECISION           - SINGLE/DOUBLE
C   LANGUAGE            - FORTRAN
C-----------------------------------------------------------------------
C   LATEST REVISION     - MARCH 9, 1977
C
      SUBROUTINE EBBCKC (ZR,ZI,N,IZ,K,L,M,D)
C
      DIMENSION          ZR(IZ,1),ZI(IZ,1),D(1)
C*    DOUBLE PRECISION   ZR,ZI,D,S
      IF (L .EQ. K) GO TO 15
      DO 10 I = K,L
         S = D(I)
C                                  LEFT HAND EIGENVECTORS ARE BACK
C                                    TRANSFORMED IF THE ABOVE
C                                    STATEMENT IS REPLACED BY S=1.0/D(I)
         DO 5 J = 1,M
            ZR(I,J) = ZR(I,J)*S
            ZI(I,J) = ZI(I,J)*S
    5    CONTINUE
   10 CONTINUE
C                                  DO 25 I=K-1,1,-1 AND
C                                    DO I=L+1,N,1
   15 DO 25 II = 1,N
         I = II
         IF (I .GE. K .AND. I .LE. L) GO TO 25
         IF (I .LT. K) I = K-II
         KK = D(I)
         IF (KK .EQ. I) GO TO 25
         DO 20 J = 1,M
            S = ZR(I,J)
            ZR(I,J) = ZR(KK,J)
            ZR(KK,J) = S
            S = ZI(I,J)
            ZI(I,J) = ZI(KK,J)
            ZI(KK,J) = S
   20    CONTINUE
   25 CONTINUE
      RETURN
      END
C     SUBROUTINE EHESSC (AR,AI,K,L,N,IA,ID)
C
C-EHESSC--------S/D-----LIBRARY 1---------------------------------------
C
C   FUNCTION            - REDUCTION OF A COMPLEX MATRIX TO COMPLEX
C                           UPPER HESSENBERG FORM.
C   USAGE               - CALL EHESSC(AR,AI,K,L,N,IA,ID)
C   PARAMETERS   AR     - INPUT/OUTPUT MATRIX OF DIMENSION N BY N.
C                           ON INPUT CONTAINS THE REAL COMPONENTS OF
C                           THE MATRIX TO BE REDUCED.
C                           ON OUTPUT CONTAINS THE REAL COMPONENTS
C                           OF THE REDUCED HESSENBERG FORM IN THE
C                           UPPER TRIANGULAR PORTION (INCLUDING MAIN
C                           AND SUB-DIAGONAL) AND THE DETAILS OF
C                           THE REDUCTION IN THE LOWER TRIANGULAR
C                           PORTION.
C                AI     - INPUT/OUTPUT MATRIX OF DIMENSION N BY N
C                           CONTAINING THE IMAGINARY COUNTERPARTS
C                           TO AR, ABOVE.
C                K      - INPUT SCALAR CONTAINING THE ROW AND COLUMN
C                           INDEX OF THE STARTING ELEMENT TO BE
C                           REDUCED BY ROW SCALING.  FOR UNBALANCED
C                           MATRICES SET K = 1.
C                L      - INPUT SCALAR CONTAINING THE ROW AND
C                           COLUMN INDEX OF THE LAST ELEMENT TO BE
C                           REDUCED BY ROW SCALING. FOR UNBALANCED
C                           MATRICES SET L = N.
C                N      - INPUT SCALAR CONTAINING THE ORDER OF
C                           THE MATRIX TO BE REDUCED.
C                IA     - INPUT SCALAR CONTAINING ROW DIMENSION
C                           OF AR AND AI IN THE CALLING PROGRAM.
C                ID     - OUTPUT VECTOR OF LENGTH L CONTAINING
C                           DETAILS OF THE TRANSFORMATIONS.
C   PRECISION           - SINGLE/DOUBLE
C   LANGUAGE            - FORTRAN
C-----------------------------------------------------------------------
C   LATEST REVISION     - FEBRUARY 7, 1973
C
      SUBROUTINE EHESSC (AR,AI,K,L,N,IA,ID)
C
      DIMENSION          AR(IA,1),AI(IA,1),ID(1),T1(2),T2(2)
C*    DOUBLE PRECISION   AR,AI,XR,XI,YR,YI,T1,T2,ZERO
C1    COMPLEX*16         X,Y
      COMPLEX            X,Y
      EQUIVALENCE        (X,T1(1),XR),(T1(2),XI),(Y,T2(1),YR),
     1                   (T2(2),YI)
C1    DATA               ZERO/0.0D0/
      DATA               ZERO/0.0/
      LA=L-1
      KP1=K+1
      IF (LA .LT. KP1) GO TO 45
      DO 40 M=KP1,LA
         I=M
         XR=ZERO
         XI=ZERO
         DO 5 J=M,L
C*          IF (DABS(AR(J,M-1))+DABS(AI(J,M-1)) .LE. DABS(XR)+DABS(XI))
C2   1      GO TO 5
            IF ( ABS(AR(J,M-1))+ ABS(AI(J,M-1)) .LE.  ABS(XR)+ ABS(XI))
     1      GO TO 5
            XR=AR(J,M-1)
            XI=AI(J,M-1)
            I=J
    5    CONTINUE
         ID(M)=I
         IF (I .EQ. M) GO TO 20
C                                  INTERCHANGE ROWS AND COLUMNS OF
C                                    ARRAYS AR AND AI
         MM1=M-1
         DO 10 J=MM1,N
            YR=AR(I,J)
            AR(I,J)=AR(M,J)
            AR(M,J)=YR
            YI=AI(I,J)
            AI(I,J)=AI(M,J)
            AI(M,J)=YI
   10    CONTINUE
         DO 15 J=1,L
            YR=AR(J,I)
            AR(J,I)=AR(J,M)
            AR(J,M)=YR
            YI=AI(J,I)
            AI(J,I)=AI(J,M)
            AI(J,M)=YI
   15    CONTINUE
C                                  END INTERCHANGE
   20    IF (XR .EQ. ZERO .AND. XI .EQ. ZERO) GO TO 40
         MP1=M+1
         DO 35 I=MP1,L
            YR=AR(I,M-1)
            YI=AI(I,M-1)
            IF (YR .EQ. ZERO .AND. YI .EQ. ZERO) GO TO 35
            Y=Y/X
            AR(I,M-1)=YR
            AI(I,M-1)=YI
            DO 25 J=M,N
               AR(I,J)=AR(I,J)-YR*AR(M,J)+YI*AI(M,J)
               AI(I,J)=AI(I,J)-YR*AI(M,J)-YI*AR(M,J)
   25       CONTINUE
            DO 30 J=1,L
               AR(J,M)=AR(J,M)+YR*AR(J,I)-YI*AI(J,I)
               AI(J,M)=AI(J,M)+YR*AI(J,I)+YI*AR(J,I)
   30       CONTINUE
   35    CONTINUE
   40 CONTINUE
   45 RETURN
      END
C     SUBROUTINE EIGCC  (A,N,IA,IJOB,W,Z,IZ,WK,IER)
C
C-EIGCC---------S/D-----LIBRARY 1---------------------------------------
C
C   FUNCTION            - TO CALCULATE EIGENVALUES AND (OPTIONALLY)
C                           EIGENVECTORS OF A COMPLEX GENERAL MATRIX.
C   USAGE               - CALL EIGCC (A,N,IA,IJOB,W,Z,IZ,WK,IER)
C   PARAMETERS   A      - THE INPUT COMPLEX GENERAL MATRIX OF ORDER N
C                           WHOSE EIGENVALUES AND EIGENVECTORS ARE
C                           TO BE COMPUTED. INPUT A IS DESTROYED IF
C                           IJOB IS EQUAL TO 0 OR 1.
C                N      - THE ORDER OF THE MATRIX A AND MATRIX Z.
C                IA     - THE ROW DIMENSION OF THE MATRIX A IN THE
C                           CALLING PROGRAM. IA MUST BE GREATER THAN OR
C                           EQUAL TO N.
C                IJOB   - INPUT OPTION PARAMETER, WHEN
C                           IJOB = 0, COMPUTE EIGENVALUES ONLY
C                           IJOB = 1, COMPUTE EIGENVALUES AND EIGEN-
C                             VECTORS.
C                           IJOB = 2, COMPUTE EIGENVALUES, EIGENVECTORS
C                             AND PERFORMANCE INDEX.
C                           IJOB = 3, COMPUTE PERFORMANCE INDEX ONLY.
C                           IF THE PERFORMANCE INDEX IS COMPUTED, IT IS
C                           RETURNED IN WK(1). THE ROUTINES HAVE
C                           PERFORMED (WELL, SATISFACTORILY, POORLY) IF
C                           WK(1) IS (LESS THAN 1, BETWEEN 1 AND 100,
C                           GREATER THAN 100).
C                W      - THE OUTPUT COMPLEX VECTOR OF LENGTH N,
C                           CONTAINING THE EIGENVALUES OF A.
C                Z      - THE OUTPUT N BY N COMPLEX MATRIX CONTAINING
C                           THE EIGENVECTORS OF A.
C                           THE EIGENVECTOR IN COLUMN J OF Z CORRES-
C                           PONDS TO THE EIGENVALUE W(J).
C                           IF IJOB = 0, Z IS NOT USED.
C                IZ     - THE ROW DIMENSION OF THE MATRIX Z IN THE
C                           CALLING PROGRAM. IZ MUST BE GREATER THAN
C                           OR EQUAL TO N IF IJOB IS NOT EQUAL TO ZERO.
C                WK     - WORK AREA, THE LENGTH OF WK DEPENDS
C                           ON THE VALUE OF IJOB, WHEN
C                           IJOB = 0, THE LENGTH OF WK IS AT LEAST 2N.
C                           IJOB = 1, THE LENGTH OF WK IS AT LEAST 2N.
C                           IJOB = 2, THE LENGTH OF WK IS AT LEAST
C                             2(N*N+N).
C                           IJOB = 3, THE LENGTH OF WK IS AT LEAST 1.
C                IER    - ERROR PARAMETER
C                         TERMINAL ERROR
C                           IER = 128+J, INDICATES THAT ELRH1C OR ELRH2C
C                             FAILED TO CONVERGE ON EIGENVALUE J.
C                           EIGENVALUES J+1,J+2,...,N HAVE BEEN
C                           COMPUTED CORRECTLY.
C                           THE PERFORMANCE INDEX IS SET TO 1000.0.
C                         WARNING ERROR (WITH FIX)
C                           IER = 66, INDICATES IJOB IS LESS THAN 0 OR
C                             IJOB IS GREATER THAN 3. IJOB SET TO 1.
C                           IER = 67, INDICATES IJOB IS NOT EQUAL TO
C                             ZERO, AND IZ IS LESS THAN THE ORDER OF
C                             MATRIX A. IJOB IS SET TO ZERO.
C   PRECISION           - SINGLE/DOUBLE
C   REQD. IMSL ROUTINES - EBALAC,EBBCKC,EHESSC,ELRH1C,ELRH2C,UERTST
C   LANGUAGE            - FORTRAN
C-----------------------------------------------------------------------
C   LATEST REVISION     - APRIL 1, 1974
C
      SUBROUTINE EIGCC (A,N,IA,IJOB,W,Z,IZ,WK,IER)
C
      DIMENSION          A(1),W(1),WK(1),Z(1)
C*    DOUBLE PRECISION   A,W,WK,Z,ANORM,ASUM,PI,SUMZ,SUMR,SUMI,S
C*    DOUBLE PRECISION   TEN,RDELP,ZERO,ONE,THOUS,AN
C*    DATA               RDELP/Z3410000000000000/
C2    DATA               ZERO,ONE,TEN,THOUS/0.0D0,1.0D0,10.0D0,1000.0D0/
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     DATA               RDELP/Z3C100000/
      DATA               RDELP/1.0E-6/
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DATA               ZERO,ONE/0.0,1.0/,TEN/10.0/,THOUS/1000.0/
C                                  INITIALIZE ERROR PARAMETERS
      IER = 0
      JER = 0
      IAA = IA+IA
      IZZ = IZ+IZ
      IF (IJOB .GE. 0 .AND. IJOB .LE. 3) GO TO 5
C                                  WARNING ERROR - IJOB IS NOT IN THE
C                                    RANGE
      IER = 66
      IJOB = 1
      GO TO 10
    5 IF (IJOB .EQ. 0) GO TO 20
   10 IF (IZ .GE. N) GO TO 15
C                                  WARNING ERROR - IZ IS LESS THAN N
C                                    EIGENVECTORS CAN NOT BE COMPUTED,
C                                    IJOB SET TO ZERO
      IER = 67
      IJOB = 0
   15 IF (IJOB .EQ. 3) GO TO 100
C                                  PACK A INTO AN N BY N ARRAY
C                                    SEPARATE A INTO REAL AND IMAGINARY
C                                    PARTS AND PACK INTO AN N BY N
C                                    MATRIX
   20 NN = N+N
      N1 = 1
      IF (IJOB .EQ. 2) N1 = NN*N+1
      N2 = N1+N
      IF (IJOB .EQ. 0) N2 = 1
      L = 1
      M = 1
      NM = N1-1
      KA1 = 2
      DO 45 J = 1,N
         K = KA1
         DO 30 I = 1,N
            IF (IJOB .NE. 2) GO TO 25
            WK(M) = A(K-1)
            WK(M+1) = A(K)
            M = M+2
   25       WK(NM+I) = A(K)
            A(L) = A(K-1)
            K = K+2
            L = L+1
   30    CONTINUE
         DO 40 I = 1,N
   35       A(L) = WK(NM+I)
            L = L+1
   40    CONTINUE
         KA1 = KA1+IAA
   45 CONTINUE
C                                  BALANCE THE INPUT A
      CALL EBALAC (A(1),A(N+1),N,NN,K,L,WK(N1))
C                                  REDUCES A TO HESSENBERG FORM
      CALL EHESSC (A(1),A(N+1),K,L,N,NN,WK(N2))
      IF (IJOB .NE. 0) GO TO 50
C                                  COMPUTE EIGENVALUES ONLY
      CALL ELRH1C (A(1),A(N+1),K,L,N,NN,W(1),W(N+1),INFER,JER)
      GO TO 55
C                                  COMPUTE EIGENVALUES AND EIGENVECTORS
   50 CALL ELRH2C (A(1),A(N+1),K,L,N,NN,W(1),W(N+1),Z(1),Z(N+1),
     1   WK(N2),INFER,JER)
C                                  BACKTRANSFORM THE EIGENVECTORS
      CALL EBBCKC (Z(1),Z(N+1),N,NN,K,L,N,WK(N1))
   55 N2 = N2-1
      DO 60 I=1,N
         NPI = N+I
         N2PI = N2+I
         WK(N2PI) = W(NPI)
   60 CONTINUE
      JW = N+N
      J = N
      DO 65 I=1,N
         W(JW-1) = W(J)
         N2PJ = N2+J
         W(JW) = WK(N2PJ)
         JW = JW-2
         J = J-1
   65 CONTINUE
      IF (IJOB .EQ. 0) GO TO 9000
      IF (IJOB .NE. 2) GO TO 80
C                                  MOVE ORIGINAL MATRIX BACK TO A
      K = 1
      KA1 = 1
      DO 75 J = 1,N
         L = KA1
         DO 70 I = 1,N
            A(L) = WK(K)
            A(L+1) = WK(K+1)
            K = K+2
            L = L+2
   70    CONTINUE
         KA1 = KA1+IAA
   75 CONTINUE
C                                  CONVERT Z (EIGENVECTORS) TO COMPLEX
C                                    FORMAT Z(IZ,N)
   80 NP1 = N+1
      NM1 = N-1
      KA1 = NN*N-NM1
      KA2 = NN+IZZ*(N-1)
      KA3 = NN*N-N
      DO 95 JB = 1,N
         J = NP1-JB
         K = KA1
         DO 85 I = 1,N
            WK(I) = Z(K)
            K = K+1
   85    CONTINUE
         JZ = KA2
         JR = KA3
         JI = N
         DO 90 I = 1,N
            Z(JZ-1) = Z(JR)
            Z(JZ) = WK(JI)
            JZ = JZ-2
            JR = JR-1
            JI = JI-1
   90    CONTINUE
         KA1 = KA1-NN
         KA2 = KA2-IZZ
         KA3 = KA3-NN
   95 CONTINUE
C                                  Z IS NOW IN COMPLEX FORMAT Z(IZ,N).
      IF (IJOB .LE. 1) GO TO 9000
      WK(1) = THOUS
      IF (JER .NE. 0) GO TO 9000
C                                  COMPUTE 1-NORM OF A
  100 ANORM = ZERO
      KA1 = 1
      DO 110 J=1,N
         ASUM = ZERO
         KA = KA1
         DO 105 I=1,N
C1          ASUM = ASUM+CDABS(DCMPLX(A(KA),A(KA+1)))
            ASUM = ASUM+CABS(CMPLX(A(KA),A(KA+1)))
            KA = KA+2
  105    CONTINUE
C1       ANORM = DMAX1(ANORM,ASUM)
         ANORM = AMAX1(ANORM,ASUM)
         KA1 = KA1+IAA
  110 CONTINUE
      IF (ANORM .EQ. ZERO) ANORM = ONE
C                                  COMPUTE PERFORMANCE INDEX
      PI = ZERO
      LW = 1
      KA1 = 1
      DO 125 J=1,N
         S = ZERO
         SUMR = ZERO
         SUMI = ZERO
         SUMZ = ZERO
         LZ = KA1
         DO 120 L=1,N
C1          SUMZ = SUMZ+CDABS(DCMPLX(Z(LZ),Z(LZ+1)))
            SUMZ = SUMZ+CABS(CMPLX(Z(LZ),Z(LZ+1)))
            KA = L+L-1
            KZ = KA1
            DO 115 K=1,N
               SUMR =SUMR+A(KA)*Z(KZ)-A(KA+1)*Z(KZ+1)
               SUMI = SUMI+A(KA)*Z(KZ+1)+A(KA+1)*Z(KZ)
               KA = KA+IAA
               KZ = KZ+2
  115       CONTINUE
            SUMR = SUMR-W(LW)*Z(LZ)+W(LW+1)*Z(LZ+1)
            SUMI = SUMI-W(LW)*Z(LZ+1)-W(LW+1)*Z(LZ)
C1          S = S+CDABS(DCMPLX(SUMR,SUMI))
            S = S+CABS(CMPLX(SUMR,SUMI))
            LZ = LZ+2
  120    CONTINUE
C1       PI = DMAX1(PI,S/SUMZ)
         PI = AMAX1(PI,S/SUMZ)
         LW = LW+2
         KA1 = KA1+IZZ
  125 CONTINUE
      AN = N
      PI = PI/(ANORM*TEN*AN*RDELP)
      WK(1) = PI
 9000 CONTINUE
      IF (IER .NE. 0) CALL UERTST (IER,6HEIGCC )
      IF (JER .EQ. 0) GO TO 9005
      IER = JER+INFER
      CALL UERTST (IER,6HEIGCC )
 9005 RETURN
      END
C     SUBROUTINE ELRH1C (HR,HI,K,L,N,IH,WR,WI,INFER,IER)
C
C-ELRH1C--------S/D-----LIBRARY 1---------------------------------------
C
C   FUNCTION            - COMPUTATION OF ALL EIGENVALUES OF A COMPLEX
C                           UPPER HESSENBERG MATRIX.
C   USAGE               - CALL ELRH1C (HR,HI,K,L,N,IH,WR,WI,INFER,IER)
C   PARAMETERS   HR     - INPUT MATRIX OF DIMENSION N BY N CONTAINING
C                           THE REAL COMPONENTS OF THE COMPLEX
C                           HESSENBERG MATRIX. HR IS DESTROYED ON
C                           OUTPUT.
C                HI     - INPUT MATRIX OF DIMENSION N BY N CONTAINING
C                           THE IMAGINARY COUNTERPARTS TO HR, ABOVE.
C                           HI IS DESTROYED ON OUTPUT.
C                K      - INPUT SCALAR CONTAINING THE LOWER BOUNDARY
C                           INDEX FOR THE INPUT MATRIX.
C                           FOR UNBALANCED MATRICES SET K = 1.
C                L      - INPUT SCALAR CONTAINING THE UPPER BOUNDARY
C                           INDEX FOR THE INPUT MATRIX.
C                           FOR UNBALANCED MATRICES SET L = N.
C                N      - INPUT SCALAR CONTAINING THE ORDER OF THE
C                           MATRIX.
C                IH     - INPUT SCALAR CONTAINING THE ROW DIMENSION
C                           OF MATRICES HR AND HI, IN THE CALLING
C                           PROGRAM.
C                WR     - OUTPUT VECTOR OF LENGTH N CONTAINING THE REAL
C                           COMPONENTS OF THE EIGENVALUES.
C                WI     - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           IMAGINARY COMPONENTS OF THE EIGENVALUES.
C                INFER  - OUTPUT SCALAR CONTAINING THE INDEX OF THE
C                           EIGENVALUE WHICH GENERATED THE TERMINAL
C                           ERROR(SEE DESCRIPTION OF IER, BELOW).
C                IER    - ERROR PARAMETER
C                           TERMINAL ERROR = 128 + N
C                             N = 1 INDICATES THE EIGENVALUE RECORDED
C                                   IN THE OUTPUT PARAMETER, INFER,
C                                   COULD NOT BE DETERMINED AFTER 30
C                                   ITERATIONS. IF THE J-TH EIGENVALUE
C                                   COULD NOT BE SO DETERMINED,
C                                   THEN THE EIGENVALUES J+1,J+2,...,N
C                                   SHOULD BE CORRECT.
C   PRECISION           - SINGLE/DOUBLE
C   REQD. IMSL ROUTINES - UERTST
C   LANGUAGE            - FORTRAN
C-----------------------------------------------------------------------
C   LATEST REVISION     - MARCH 22, 1973
C
      SUBROUTINE ELRH1C (HR,HI,K,L,N,IH,WR,WI,INFER,IER)
C
      DIMENSION          HR(IH,1),HI(IH,1),WR(1),WI(1)
      DIMENSION          T1(2),T2(2),T3(2)
C1    COMPLEX*16         X,Y,Z
      COMPLEX            X,Y,Z
C*    DOUBLE PRECISION   HR,HI,WR,WI,EPS,ZR,ZI,T1,T2,T3,RDELP
C*    DOUBLE PRECISION   ZERO,ONE,TWO,SR,SI,XR,XI,YR,YI,TR,TI
      EQUIVALENCE        (X,T1(1),XR),(T1(2),XI),
     1                   (Y,T2(1),YR),(T2(2),YI),
     2                   (Z,T3(1),ZR),(T3(2),ZI)
C*    DATA               TWO/2.0D0/
C1    DATA               ZERO,ONE,RDELP/0.D0,1.D0,Z3410000000000000/
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     DATA               ZERO,ONE,TWO,RDELP/0.0,1.0,2.0,Z3C100000/
      DATA               ZERO,ONE,TWO,RDELP/0.0,1.0,2.0,1.0E-6/
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INFER=0
      IER=0
      DO 5 I=1,N
         IF (I .GE. K .AND. I .LE. L) GO TO 5
         WR(I)=HR(I,I)
         WI(I)=HI(I,I)
    5 CONTINUE
      NN=L
      TR=ZERO
      TI=ZERO
C                                  SEARCH FOR NEXT EIGENVALUE
   10 IF (NN .LT. K) GO TO 9005
      ITS=0
      NM1=NN-1
      IF (NN .EQ. K) GO TO 25
C                                  LOOK FOR SINGLE SMALL SUB-DIAGONAL
C                                    ELEMENTS
   15 NPL=NN+K
      DO 20 LL=K,NM1
         M=NPL-LL
         MM1=M-1
C*       IF (DABS(HR(M,MM1)) + DABS(HI(M,MM1)) .LE.
C*   1    RDELP*(DABS(HR(MM1,MM1)) + DABS(HI(MM1,MM1)) +
C3   2      DABS(HR(M,M)) + DABS(HI(M,M)))) GO TO 30
         IF ( ABS(HR(M,MM1)) +  ABS(HI(M,MM1)) .LE.
     1    RDELP*( ABS(HR(MM1,MM1)) +  ABS(HI(MM1,MM1)) +
     2       ABS(HR(M,M)) +  ABS(HI(M,M)))) GO TO 30
   20 CONTINUE
   25 M=K
   30 IF (M .EQ. NN) GO TO 110
      IF (ITS .EQ. 30) GO TO 115
C                                  FORM SHIFT
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 35
      SR=HR(NN,NN)
      SI=HI(NN,NN)
      XR=HR(NM1,NN)*HR(NN,NM1)-HI(NM1,NN)*HI(NN,NM1)
      XI=HR(NM1,NN)*HI(NN,NM1)+HI(NM1,NN)*HR(NN,NM1)
      IF (XR .EQ. ZERO .AND. XI .EQ. ZERO) GO TO 40
      YR=(HR(NM1,NM1)-SR)/TWO
      YI=(HI(NM1,NM1)-SI)/TWO
C1    Z=CDSQRT(DCMPLX(YR**2-YI**2+XR,TWO*YR*YI+XI))
      Z= CSQRT( CMPLX(YR**2-YI**2+XR,TWO*YR*YI+XI))
      IF (YR*ZR+YI*ZI .LT. ZERO) Z=-Z
      X=X/(Y+Z)
      SR=SR-XR
      SI=SI-XI
      GO TO 40
C* 35 SR=DABS(HR(NN,NM1))+DABS(HR(NM1,NN-2))
C2    SI=DABS(HI(NN,NM1))+DABS(HI(NM1,NN-2))
   35 SR= ABS(HR(NN,NM1))+ ABS(HR(NM1,NN-2))
      SI= ABS(HI(NN,NM1))+ ABS(HI(NM1,NN-2))
   40 DO 45 I=K,NN
         HR(I,I)=HR(I,I)-SR
         HI(I,I)=HI(I,I)-SI
   45 CONTINUE
      TR=TR+SR
      TI=TI+SI
      ITS=ITS+1
C                                  LOOK FOR TWO CONSECUTIVE SMALL
C                                    SUB-DIAGONAL ELEMENTS
C*    XR=DABS(HR(NM1,NM1))+DABS(HI(NM1,NM1))
C*    YR=DABS(HR(NN,NM1))+DABS(HI(NN,NM1))
C3    ZR=DABS(HR(NN,NN))+DABS(HI(NN,NN))
      XR= ABS(HR(NM1,NM1))+ ABS(HI(NM1,NM1))
      YR= ABS(HR(NN,NM1))+ ABS(HI(NN,NM1))
      ZR= ABS(HR(NN,NN))+ ABS(HI(NN,NN))
      NMJ=NM1-M
      IF (NMJ .EQ. 0) GO TO 55
C                                  FOR MM=NN-1 STEP -1 UNTIL M+1 DO
      DO 50 J=1,NMJ
         MM=NN-J
         M1=MM-1
         YI=YR
C1       YR=DABS(HR(MM,M1))+DABS(HI(MM,M1))
         YR= ABS(HR(MM,M1))+ ABS(HI(MM,M1))
         XI=ZR
         ZR=XR
C1       XR=DABS(HR(M1,M1))+DABS(HI(M1,M1))
         XR= ABS(HR(M1,M1))+ ABS(HI(M1,M1))
         IF (YR.LE.RDELP*ZR/YI*(ZR+XR+XI)) GO TO 60
   50 CONTINUE
   55 MM=M
C                                  TRIANGULAR DECOMPOSITION
   60 MP1=MM+1
      DO 85 I=MP1,NN
         IM1=I-1
         XR=HR(IM1,IM1)
         XI=HI(IM1,IM1)
         YR=HR(I,IM1)
         YI=HI(I,IM1)
C1       IF (DABS(XR)+DABS(XI) .GE. DABS(YR)+DABS(YI)) GO TO 70
         IF ( ABS(XR)+ ABS(XI) .GE.  ABS(YR)+ ABS(YI)) GO TO 70
C                                  INTERCHANGE ROWS OF HR AND HI
         DO 65 J=IM1,NN
            ZR=HR(IM1,J)
            HR(IM1,J)=HR(I,J)
            HR(I,J)=ZR
            ZI=HI(IM1,J)
            HI(IM1,J)=HI(I,J)
            HI(I,J)=ZI
   65    CONTINUE
         Z=X/Y
         WR(I)=ONE
         GO TO 75
   70    Z=Y/X
         WR(I)=-ONE
   75    HR(I,IM1)=ZR
         HI(I,IM1)=ZI
         DO 80 J=I,NN
            HR(I,J)=HR(I,J)-ZR*HR(IM1,J)+ZI*HI(IM1,J)
            HI(I,J)=HI(I,J)-ZR*HI(IM1,J)-ZI*HR(IM1,J)
   80    CONTINUE
   85 CONTINUE
C                                  COMPOSITION
      DO 105 J=MP1,NN
         JM1=J-1
         XR=HR(J,JM1)
         XI=HI(J,JM1)
         HR(J,JM1)=ZERO
         HI(J,JM1)=ZERO
C                                  INTERCHANGE COLUMNS OF HR AND HI IF
C                                    NECESSARY
         IF (WR(J) .LE. ZERO) GO TO 95
         DO 90 I=M,J
            ZR=HR(I,JM1)
            HR(I,JM1)=HR(I,J)
            HR(I,J)=ZR
            ZI=HI(I,JM1)
            HI(I,JM1)=HI(I,J)
            HI(I,J)=ZI
   90    CONTINUE
   95    DO 100 I=M,J
            HR(I,JM1)=HR(I,JM1)+XR*HR(I,J)-XI*HI(I,J)
            HI(I,JM1)=HI(I,JM1)+XR*HI(I,J)+XI*HR(I,J)
  100    CONTINUE
  105 CONTINUE
      GO TO 15
C                                  A ROOT FOUND
  110 WR(NN)=HR(NN,NN)+TR
      WI(NN)=HI(NN,NN)+TI
      NN=NM1
      GO TO 10
C                                  SET ERROR-NO CONVERGENCE TO AN
C                                    EIGENVALUE AFTER 30 ITERATIONS
  115 INFER=NN
      IER=129
 9000 CONTINUE
      CALL UERTST (IER,6HELRH1C)
 9005 RETURN
      END
C     SUBROUTINE ELRH2C (HR,HI,K,L,N,IH,WR,WI,ZR,ZI,ID,INFER,IER)
C
C-ELRH2C--------S/D-----LIBRARY 1---------------------------------------
C
C   FUNCTION            - COMPUTE THE EIGENVALUES AND EIGENVECTORS OF
C                           A COMPLEX UPPER HESSENBERG MATRIX AND
C                           BACK TRANSFORM THE EIGENVECTORS.
C   USAGE               - CALL ELRH2C (HR,HI,K,L,N,IH,WR,WI,ZR,ZI,ID,
C                             INFER,IER)
C   PARAMETERS   HR     - INPUT MATRIX OF DIMENSION N BY N CONTAINING
C                           THE REAL COMPONENTS OF THE COMPLEX
C                           HESSENBERG MATRIX. HR IS DESTROYED ON
C                           OUTPUT.
C                HI     - INPUT MATRIX OF DIMENSION N BY N CONTAINING
C                           THE IMAGINARY COUNTERPARTS TO HR, ABOVE.
C                           HI IS DESTROYED ON OUTPUT.
C                K      - INPUT SCALAR CONTAINING THE LOWER BOUNDARY
C                           INDEX FOR THE INPUT MATRIX.
C                           FOR UNBALANCED MATRICES SET K = 1.
C                L      - INPUT SCALAR CONTAINING THE UPPER BOUNDARY
C                           INDEX FOR THE INPUT MATRIX.
C                           FOR UNBALANCED MATRICES SET L = N.
C                N      - INPUT SCALAR CONTAINING THE ORDER OF THE
C                           HESSENBERG MATRIX AND THE EIGENVECTOR
C                           MATRIX.
C                IH     - INPUT SCALAR CONTAINING THE ROW DIMENSION
C                           OF MATRICES HR, HI, ZR AND ZI IN THE
C                           CALLING PROGRAM.
C                WR     - OUTPUT VECTOR OF LENGTH N CONTAINING THE REAL
C                           COMPONENTS OF THE EIGENVALUES.
C                WI     - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           IMAGINARY COMPONENTS OF THE EIGENVALUES.
C                ZR     - OUTPUT MATRIX OF DIMENSION N BY N CONTAINING
C                           THE REAL COMPONENTS OF THE EIGENVECTORS.
C                           THE EIGENVECTORS ARE NOT NORMALIZED.
C                ZI     - OUTPUT MATRIX OF DIMENSION N BY N CONTAINING
C                           THE IMAGINARY COUNTERPARTS TO ZR, ABOVE.
C                ID     - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           INFORMATION GENERATED BY IMSL ROUTINE
C                           EHESSC IDENTIFYING THE ROWS AND COLUMNS
C                           INTERCHANGED DURING THE REDUCTION TO
C                           HESSENBERG FORM. ONLY COMPONENTS K
C                           THROUGH L ARE USED.
C                INFER  - OUTPUT SCALAR CONTAINING THE INDEX OF THE
C                           EIGENVALUE WHICH GENERATED THE TERMINAL
C                           ERROR(SEE DESCRIPTION OF IER, BELOW).
C                IER    - ERROR PARAMETER
C                           TERMINAL ERROR = 128 + N
C                             N = 1 INDICATES THE EIGENVALUE RECORDED
C                                   IN THE OUTPUT PARAMETER, INFER,
C                                   COULD NOT BE DETERMINED AFTER 30
C                                   ITERATIONS. IF THE J-TH EIGENVALUE
C                                   COULD NOT BE SO DETERMINED,
C                                   THEN THE EIGENVALUES J+1,J+2,...,N
C                                   SHOULD BE CORRECT.
C   PRECISION           - SINGLE/DOUBLE
C   REQD. IMSL ROUTINES - UERTST
C   LANGUAGE            - FORTRAN
C-----------------------------------------------------------------------
C   LATEST REVISION     - APRIL 5, 1977
C
      SUBROUTINE ELRH2C (HR,HI,K,L,N,IH,WR,WI,ZR,ZI,ID,INFER,IER)
C
      DIMENSION          HR(IH,1),HI(IH,1),WR(1),WI(1),ZR(IH,1),ZI(IH,1)
      DIMENSION          T1(2),T2(2),T3(2),ID(1)
C1    COMPLEX*16         X,Y,Z
      COMPLEX            X,Y,Z
C*    DOUBLE PRECISION   HR,HI,WR,WI,ZR,ZI,EPS,FNORM
C*    DOUBLE PRECISION   T1,T2,T3,XR,XI,YR,YI,ZZR,ZZI,SR,SI,TR,TI
C*    DOUBLE PRECISION   ZERO,ONE,TWO
      EQUIVALENCE        (X,T1(1),XR),(T1(2),XI),
     1                   (Y,T2(1),YR),(T2(2),YI),
     2                   (Z,T3(1),ZZR),(T3(2),ZZI)
C*    DATA               ZERO,ONE,TWO/0.0D0,1.0D0,2.0D0/
C1    DATA               EPS/Z3410000000000000/
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     DATA               ZERO,ONE,TWO,EPS/0.0,1.0,2.0,Z3C100000/
      DATA               ZERO,ONE,TWO,EPS/0.0,1.0,2.0,1.0E-6/
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                  INITIALIZE IER
      IER=0
      INFER=0
      TR=ZERO
      TI=ZERO
      DO 5 I=1,N
         DO 3 J=1,N
            ZR(I,J)=ZERO
            ZI(I,J)=ZERO
    3    CONTINUE
         ZR(I,I) = ONE
    5 CONTINUE
C                                  FORM THE MATRIX OF ACCUMULATED
C                                    TRANSFORMATIONS FROM THE INFOR-
C                                    MATION LEFT BY ROUTINE "EHESSC"
      IEND=L-K-1
      IF (IEND .LE. 0) GO TO 25
C                                  DO I=L-1,K+1,-1
      DO 20 II=1,IEND
         I=L-II
         IP1=I+1
         IM1=I-1
         DO 10 M=IP1,L
            ZR(M,I)=HR(M,IM1)
            ZI(M,I)=HI(M,IM1)
   10    CONTINUE
         J=ID(I)
         IF (I .EQ. J) GO TO 20
         DO 15 M=I,L
            ZR(I,M)=ZR(J,M)
            ZI(I,M)=ZI(J,M)
            ZR(J,M)=ZERO
            ZI(J,M)=ZERO
   15    CONTINUE
         ZR(J,I)=ONE
   20 CONTINUE
   25 DO 30 I=1,N
         IF (I .GE. K .AND. I .LE. L) GO TO 30
         WR(I)=HR(I,I)
         WI(I)=HI(I,I)
   30 CONTINUE
      NN=L
C                                  SEARCH FOR NEXT EIGENVALUE
   35 IF (NN .LT. K) GO TO 150
      ITS=0
      NNM1=NN-1
      NNM2=NN-2
      IF (NN .EQ. K) GO TO 50
C                                  LOOK FOR SINGLE SMALL SUB-DIAGONAL
C                                    ELEMENT
   40 NPL=NN+K
C                                  DO M=NN,K+1,-1
      DO 45 KK=K,NNM1
         M=NPL-KK
         MM1=M-1
C*       IF (DABS(HR(M,MM1))+DABS(HI(M,MM1)) .LE. EPS*(DABS(HR(MM1,MM1))
C2   1      +DABS(HI(MM1,MM1))+DABS(HR(M,M))+DABS(HI(M,M)))) GO TO 55
         IF ( ABS(HR(M,MM1))+ ABS(HI(M,MM1)) .LE. EPS*( ABS(HR(MM1,MM1))
     1      + ABS(HI(MM1,MM1))+ ABS(HR(M,M))+ ABS(HI(M,M)))) GO TO 55
   45 CONTINUE
   50 M=K
   55 IF (M .EQ. NN) GO TO 145
      IF (ITS .EQ. 30) GO TO 205
C                                  FORM SHIFT
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 60
      SR=HR(NN,NN)
      SI=HI(NN,NN)
      XR=HR(NNM1,NN)*HR(NN,NNM1)-HI(NNM1,NN)*HI(NN,NNM1)
      XI=HR(NNM1,NN)*HI(NN,NNM1)+HI(NNM1,NN)*HR(NN,NNM1)
      IF (XR .EQ. ZERO .AND. XI .EQ. ZERO) GO TO 65
      YR=(HR(NNM1,NNM1)-SR)/TWO
      YI=(HI(NNM1,NNM1)-SI)/TWO
C1    Z=CDSQRT(DCMPLX(YR**2-YI**2+XR,TWO*YR*YI+XI))
      Z= CSQRT( CMPLX(YR**2-YI**2+XR,TWO*YR*YI+XI))
      IF (YR*ZZR+YI*ZZI .LT. ZERO) Z=-Z
      X=X/(Y+Z)
      SR=SR-XR
      SI=SI-XI
      GO TO 65
C* 60 SR=DABS(HR(NN,NNM1))+DABS(HR(NNM1,NNM2))
C2    SI=DABS(HI(NN,NNM1))+DABS(HI(NNM1,NNM2))
   60 SR= ABS(HR(NN,NNM1))+ ABS(HR(NNM1,NNM2))
      SI= ABS(HI(NN,NNM1))+ ABS(HI(NNM1,NNM2))
   65 DO 70 I=K,NN
         HR(I,I)=HR(I,I)-SR
         HI(I,I)=HI(I,I)-SI
   70 CONTINUE
      TR=TR+SR
      TI=TI+SI
      ITS=ITS+1
C                                  LOOK FOR TWO CONSECUTIVE SMALL
C                                    SUB-DIAGONAL ELEMENTS
C*    XR=DABS(HR(NNM1,NNM1))+DABS(HI(NNM1,NNM1))
C*    YR=DABS(HR(NN,NNM1))+DABS(HI(NN,NNM1))
C3    ZZR=DABS(HR(NN,NN))+DABS(HI(NN,NN))
      XR= ABS(HR(NNM1,NNM1))+ ABS(HI(NNM1,NNM1))
      YR= ABS(HR(NN,NNM1))+ ABS(HI(NN,NNM1))
      ZZR= ABS(HR(NN,NN))+ ABS(HI(NN,NN))
      NNMJ=NNM1-M
      IF (NNMJ .EQ. 0) GO TO 80
C                                  DO MM=NN-1,M+1,-1
      DO 75 NM=1,NNMJ
         MM=NN-NM
         MMM1=MM-1
         YI=YR
C1       YR=DABS(HR(MM,MMM1))+DABS(HI(MM,MMM1))
         YR= ABS(HR(MM,MMM1))+ ABS(HI(MM,MMM1))
         XI=ZZR
         ZZR=XR
C1       XR=DABS(HR(MMM1,MMM1))+DABS(HI(MMM1,MMM1))
         XR= ABS(HR(MMM1,MMM1))+ ABS(HI(MMM1,MMM1))
         IF (YR .LE. EPS*ZZR/YI*(ZZR+XR+XI)) GO TO 85
   75 CONTINUE
   80 MM=M
C                                  TRIANGULAR DECOMPOSITION H=L*R
   85 MP1=MM+1
      DO 110 I=MP1,NN
         IM1=I-1
         XR=HR(IM1,IM1)
         XI=HI(IM1,IM1)
         YR=HR(I,IM1)
         YI=HI(I,IM1)
C1       IF (DABS(XR)+DABS(XI) .GE. DABS(YR)+DABS(YI)) GO TO 95
         IF ( ABS(XR)+ ABS(XI) .GE.  ABS(YR)+ ABS(YI)) GO TO 95
C                                  INTERCHANGE ROWS OF HR AND HI
         DO 90 J=IM1,N
            ZZR=HR(IM1,J)
            HR(IM1,J)=HR(I,J)
            HR(I,J)=ZZR
            ZZI=HI(IM1,J)
            HI(IM1,J)=HI(I,J)
            HI(I,J)=ZZI
   90    CONTINUE
         Z=X/Y
         WR(I)=ONE
         GO TO 100
   95    Z=Y/X
         WR(I)=-ONE
  100    HR(I,IM1)=ZZR
         HI(I,IM1)=ZZI
         DO 105 J=I,N
            HR(I,J)=HR(I,J)-ZZR*HR(IM1,J)+ZZI*HI(IM1,J)
            HI(I,J)=HI(I,J)-ZZR*HI(IM1,J)-ZZI*HR(IM1,J)
  105    CONTINUE
  110 CONTINUE
C                                  COMPOSITION R*L=H
      DO 140 J=MP1,NN
         JM1=J-1
         XR=HR(J,JM1)
         XI=HI(J,JM1)
         HR(J,JM1)=ZERO
         HI(J,JM1)=ZERO
C                                  INTERCHANGE COLUMNS OF HR, HI,
C                                    ZR, AND ZI IF NECESSARY
         IF (WR(J) .LE. ZERO) GO TO 125
         DO 115 I=1,J
            ZZR=HR(I,JM1)
            HR(I,JM1)=HR(I,J)
            HR(I,J)=ZZR
            ZZI=HI(I,JM1)
            HI(I,JM1)=HI(I,J)
            HI(I,J)=ZZI
  115    CONTINUE
         DO 120 I=K,L
            ZZR=ZR(I,JM1)
            ZR(I,JM1)=ZR(I,J)
            ZR(I,J)=ZZR
            ZZI=ZI(I,JM1)
            ZI(I,JM1)=ZI(I,J)
            ZI(I,J)=ZZI
  120    CONTINUE
C                                  END INTERCHANGE COLUMNS
  125    DO 130 I=1,J
            HR(I,JM1)=HR(I,JM1)+XR*HR(I,J)-XI*HI(I,J)
            HI(I,JM1)=HI(I,JM1)+XR*HI(I,J)+XI*HR(I,J)
  130    CONTINUE
         DO 135 I=K,L
            ZR(I,JM1)=ZR(I,JM1)+XR*ZR(I,J)-XI*ZI(I,J)
            ZI(I,JM1)=ZI(I,JM1)+XR*ZI(I,J)+XI*ZR(I,J)
  135    CONTINUE
C                                  END ACCUMULATE TRANSFORMATIONS
  140 CONTINUE
      GO TO 40
C                                  A ROOT FOUND
  145 WR(NN)=HR(NN,NN)+TR
      WI(NN)=HI(NN,NN)+TI
      NN=NNM1
      GO TO 35
C                                  ALL ROOTS FOUND. BACKSUBSTITUTE TO
C                                    FIND VECTORS OF UPPER TRIANGULAR
C                                    FORM
  150 IF (N .EQ. 1) GO TO 9005
      FNORM=ZERO
      DO 160 I=1,N
C1       FNORM=FNORM+DABS(WR(I))+DABS(WI(I))
         FNORM=FNORM+ ABS(WR(I))+ ABS(WI(I))
         IF (I .EQ. N) GO TO 160
         IP1=I+1
         DO 155 J=IP1,N
C1          FNORM=FNORM+DABS(HR(I,J))+DABS(HI(I,J))
            FNORM=FNORM+ ABS(HR(I,J))+ ABS(HI(I,J))
  155    CONTINUE
  160 CONTINUE
      IF (FNORM .EQ. ZERO) GO TO 9005
C                                  DO NN=N,2,-1
      NP2 = N+2
      DO 180 NM=2,N
         NN=NP2-NM
         XR=WR(NN)
         XI=WI(NN)
         NNM1=NN-1
C                                  DO I=NN-1,1,-1
         DO 175 II=1,NNM1
            I=NN-II
            ZZR=HR(I,NN)
            ZZI=HI(I,NN)
            IF (I .EQ. NNM1) GO TO 170
            IP1=I+1
            DO 165 J=IP1,NNM1
               ZZR=ZZR+HR(I,J)*HR(J,NN)-HI(I,J)*HI(J,NN)
               ZZI=ZZI+HR(I,J)*HI(J,NN)+HI(I,J)*HR(J,NN)
  165       CONTINUE
  170       YR=XR-WR(I)
            YI=XI-WI(I)
            IF (YR .EQ. ZERO .AND. YI .EQ. ZERO) YR=EPS*FNORM
            Z=Z/Y
            HR(I,NN)=T3(1)
            HI(I,NN)=T3(2)
  175    CONTINUE
  180 CONTINUE
C                                  END BACKSUBSTITUTION
      NM1=N-1
C                                  VECTOR OF ISOLATED ROOTS
      DO 190 I=1,NM1
         IF (I .GE. K .AND. I .LE. L) GO TO 190
         IP1=I+1
         DO 185 J=IP1,N
            ZR(I,J)=HR(I,J)
            ZI(I,J)=HI(I,J)
  185    CONTINUE
  190 CONTINUE
      IF (L .EQ. 0) GO TO 9005
C                                  MULTIPLY BY TRANSFORMATION MATRIX
C                                    TO GIVE VECTORS OF ORIGINAL FULL
C                                    MATRIX
      NPL=N+K
C                                  DO J=N,K+1,-1
      DO 200 JJ=K,NM1
         J=NPL-JJ
         JM1=J-1
         DO 200 I=K,L
            ZZR=ZR(I,J)
            ZZI=ZI(I,J)
            MM=JM1
            IF (L .LT. J) MM=L
            DO 195 M=K,MM
               ZZR=ZZR+ZR(I,M)*HR(M,J)-ZI(I,M)*HI(M,J)
               ZZI=ZZI+ZR(I,M)*HI(M,J)+ZI(I,M)*HR(M,J)
  195       CONTINUE
            ZR(I,J)=ZZR
            ZI(I,J)=ZZI
  200 CONTINUE
      GO TO 9005
C                                  SET ERROR - NO CONVERGENCE TO AN
C                                    EIGENVALUE AFTER 30 ITERATIONS
  205 IER=129
      INFER=NN
 9000 CONTINUE
      CALL UERTST (IER,6HELRH2C)
 9005 RETURN
      END
C     SUBROUTINE LEQT1C (A,N,IA,B,M,IB,IJOB,WA,IER)
C
C-LEQT1C--------S/D-----LIBRARY 1---------------------------------------
C
C   FUNCTION            - MATRIX DECOMPOSITION, LINEAR EQUATION
C                           SOLUTION - SPACE ECONOMIZER SOLUTION -
C                           COMPLEX MATRICES
C   USAGE               - CALL LEQT1C (A,N,IA,B,M,IB,IJOB,WA,IER)
C   PARAMETERS   A      - INPUT COMPLEX MATRIX OF DIMENSION N BY N
C                           CONTAINING THE COMPLEX COEFFICIENTS OF THE
C                           EQUATION AX = B.
C                         ON OUTPUT, A CONTAINS THE L-U DECOMPOSITION OF
C                           A ROWWISE PERMUTATION OF THE INPUT MATRIX A.
C                N      - ORDER OF MATRIX A. (INPUT)
C                IA     - ROW DIMENSION OF A AS SPECIFIED IN THE CALLING
C                           PROGRAM. (INPUT)
C                B      - INPUT COMPLEX MATRIX OF DIMENSION N BY M
C                           CONTAINING THE M COMPLEX VALUED RIGHT HAND
C                           SIDES OF THE EQUATION AX = B.
C                         ON OUTPUT, THE SOLUTION MATRIX X REPLACES B.
C                           IF IJOB=1, B IS NOT USED.
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).
C                           (INPUT)
C                IB     - ROW DIMENSION OF B AS SPECIFIED IN THE
C                           CALLING PROGRAM. (INPUT)
C                IJOB   - INPUT OPTION PARAMETER.  IJOB=I IMPLIES WHEN
C                           I=0, FACTOR THE MATRIX AND SOLVE THE
C                             EQUATION AX=B.
C                           I=1, FACTOR THE MATRIX A.
C                           I=2, SOLVE THE EQUATION AX=B.  THIS
C                             OPTION IMPLIES THAT LEQT1C HAS ALREADY
C                             BEEN CALLED USING IJOB=0 OR 1 SO THAT
C                             THE MATRIX HAS ALREADY BEEN FACTORED.  IN
C                             THIS CASE OUTPUT MATRIX A MUST HAVE BEEN
C                             SAVED FOR REUSE IN THE CALL TO LEQT1C.
C                WA     - WORK AREA OF LENGTH N CONTAINING THE PIVOT
C                           INDICES.
C                IER    - ERROR PARAMETER.
C                         TERMINAL ERROR = 128+N
C                           N=1 INDICATES THAT MATRIX A IS
C                             ALGORITHMICALLY SINGULAR.  (SEE THE
C                             CHAPTER L PRELUDE).
C   PRECISION           - SINGLE/DOUBLE
C   REQD. IMSL ROUTINES - UERTST
C   LANGUAGE            - FORTRAN
C-----------------------------------------------------------------------
C   LATEST REVISION     - JANUARY 8,1973
C
      SUBROUTINE LEQT1C (A,N,IA,B,M,IB,IJOB,WA,IER)
C
      DIMENSION          A(IA,1),B(IB,1),WA(1),T(2)
C*    DOUBLE PRECISION   P,Q,ZERO,ONE,WA,T,RN,BIG
C1    COMPLEX*16         A,B,SUM,TEMP
      COMPLEX            A,B,SUM,TEMP
      EQUIVALENCE        (SUM,T(1))
C1    DATA               ZERO/0.D0/,ONE/1.D0/
      DATA               ZERO/0.0/,ONE/1.0/
C                                  INITIALIZATION
      IER = 0
      IF (IJOB .EQ. 2) GO TO 75
      RN = N
C                                  FIND EQUILIBRATION FACTORS
      DO 10 I=1,N
         BIG = ZERO
         DO 5 J=1,N
            TEMP = A(I,J)
C1          P = CDABS(TEMP)
            P = CABS(TEMP)
            IF (P .GT. BIG) BIG = P
    5    CONTINUE
         IF (BIG .EQ. ZERO) GO TO 105
         WA(I) = ONE/BIG
   10 CONTINUE
C                                  L-U DECOMPOSITION
      DO 70 J = 1,N
         JM1 = J-1
         IF (JM1 .LT. 1) GO TO 25
C                                  COMPUTE U(I,J), I=1,...,J-1
         DO 20 I=1,JM1
            SUM = A(I,J)
            IM1 = I-1
            IF (IM1 .LT. 1) GO TO 20
            DO 15 K=1,IM1
               SUM = SUM-A(I,K)*A(K,J)
   15       CONTINUE
            A(I,J) = SUM
   20    CONTINUE
   25    P = ZERO
C                                  COMPUTE U(J,J) AND L(I,J), I=J+1,...,
         DO 45 I=J,N
            SUM = A(I,J)
            IF (JM1 .LT. 1) GO TO 40
            DO 35 K=1,JM1
               SUM = SUM-A(I,K)*A(K,J)
   35       CONTINUE
            A(I,J) = SUM
C1 40       Q = WA(I)*CDABS(SUM)
   40       Q = WA(I)*CABS(SUM)
            IF (P .GE. Q) GO TO 45
            P = Q
            IMAX = I
   45    CONTINUE
C                                  TEST FOR ALGORITHMIC SINGULARITY
         IF (RN+P .EQ. RN) GO TO 105
         IF (J .EQ. IMAX) GO TO 60
C                                  INTERCHANGE ROWS J AND IMAX
         DO 50 K=1,N
            TEMP = A(IMAX,K)
            A(IMAX,K) = A(J,K)
            A(J,K) = TEMP
   50    CONTINUE
         WA(IMAX) = WA(J)
   60    WA(J) = IMAX
         JP1 = J+1
         IF (JP1 .GT. N) GO TO 70
C                                  DIVIDE BY PIVOT ELEMENT U(J,J)
         TEMP = A(J,J)
         DO 65 I = JP1,N
            A(I,J) = A(I,J)/TEMP
   65    CONTINUE
   70 CONTINUE
   75 IF (IJOB .EQ. 1) GO TO 9005
      DO 103 K = 1,M
C                                  SOLVE UX = Y FOR X
         IW = 0
         DO 90 I = 1,N
            IMAX = WA(I)
            SUM = B(IMAX,K)
            B(IMAX,K) = B(I,K)
            IF (IW .EQ. 0) GO TO 85
            IM1 = I-1
            DO 80 J = IW,IM1
               SUM = SUM-A(I,J)*B(J,K)
   80       CONTINUE
            GO TO 88
   85       IF (T(1) .NE. ZERO .OR. T(2) .NE. ZERO) IW = I
   88       B(I,K) = SUM
   90    CONTINUE
C                                  SOLVE LY = B FOR Y
         N1 = N+1
         DO 100 IW = 1,N
            I = N1-IW
            JP1 = I+1
            SUM = B(I,K)
            IF (JP1 .GT. N) GO TO 98
            DO 95 J = JP1,N
               SUM = SUM-A(I,J)*B(J,K)
   95       CONTINUE
   98       B(I,K) = SUM/A(I,I)
  100    CONTINUE
  103 CONTINUE
      GO TO 9005
C                                  ALGORITHMIC SINGULARITY
  105 IER = 129
 9000 CONTINUE
C                                  PRINT ERROR
      CALL UERTST(IER,6HLEQT1C)
 9005 RETURN
      END
C     SUBROUTINE UERTST (IER,NAME)
C
C-UERTST----------------LIBRARY 1---------------------------------------
C
C   FUNCTION            - ERROR MESSAGE GENERATION
C   USAGE               - CALL UERTST(IER,NAME)
C   PARAMETERS   IER    - ERROR PARAMETER. TYPE + N  WHERE
C                           TYPE= 128 IMPLIES TERMINAL ERROR
C                                  64 IMPLIES WARNING WITH FIX
C                                  32 IMPLIES WARNING
C                           N   = ERROR CODE RELEVANT TO CALLING ROUTINE
C                NAME   - INPUT VECTOR CONTAINING THE NAME OF THE
C                           CALLING ROUTINE AS A SIX CHARACTER LITERAL
C                           STRING.
C   LANGUAGE            - FORTRAN
C-----------------------------------------------------------------------
C   LATEST REVISION     - JANUARY 18, 1974
C
      SUBROUTINE UERTST(IER,NAME)
C
      DIMENSION          ITYP(5,4),IBIT(4)
C     ********************************************************
C     INTEGER*2          NAME(3)
      INTEGER            NAME
C     ***************************************************************
      INTEGER            WARN,WARF,TERM,PRINTR
      EQUIVALENCE        (IBIT(1),WARN),(IBIT(2),WARF),(IBIT(3),TERM)
      DATA     ITYP      /4HWARN,4HING ,4H    ,4H    ,4H    ,
     *                    4HWARN,4HING(,4HWITH,4H FIX,4H)   ,
     *                    4HTERM,4HINAL,4H    ,4H    ,4H    ,
     *                    4HNON-,4HDEFI,4HNED ,4H    ,4H    /,
     *         IBIT      / 32,64,128,0/
      DATA     PRINTR    / 6/
      IER2=IER
      IF (IER2 .GE. WARN) GO TO 5
C                                  NON-DEFINED
      IER1=4
      GO TO 20
   5  IF (IER2 .LT. TERM) GO TO 10
C                                  TERMINAL
      IER1=3
      GO TO 20
  10  IF (IER2 .LT. WARF) GO TO 15
C                                  WARNING(WITH FIX)
      IER1=2
      GO TO 20
C                                  WARNING
  15  IER1=1
C                                  EXTRACT "N"
  20  IER2=IER2-IBIT(IER1)
C                                  PRINT ERROR MESSAGE
      WRITE (PRINTR,25) (ITYP(I,IER1),I=1,5),NAME,IER2,IER
  25  FORMAT(26H *** I M S L(UERTST) ***  ,5A4,4X,3A2,4X,I2,
     *   8H (IER = ,I3,1H))
      RETURN
      END
