C*****PROGRAM TIMEFREQ.FOR****************************************************
C
C             Written by:
C                  Clayton R. Paul
C                  Department of Electrical Engineering
C                  University of Kentucky
C                  Lexington, KY 40506
C
C     To Implement the Time-Domain to Frequency-Domain Transformation
C     Method of Computing the Time-Domain Response of a Multi-Conductor
C     Transmission Line.
C
C     The Time-Domain Input Signal is Assumed to be a Periodic, Trapezoidal
C     Pulse Train with Maximum Voltage of VMAX, Period PERIOD, Duty Cycle
C     of DC, and Equal Rise/Fall Times of TR.
C
C     The Input Data are Contained in File TIMEFREQ.IN,
C     and the Frequency-Domain Transfer Function at the Various
C     Frequencies of Interest are Contained in File MTLFREQ.DAT.
C     The Time-Domain Response is Contained in File TIMEFREQ.OUT.
C
C     The Maximum Number of Frequencies is Parameter NMAX
C
C*****************************************************************************
      PARAMETER (NMAX=1000)
      REAL FREQ(NMAX),VM(NMAX),VP(NMAX),CN(NMAX),ANG(NMAX)
      OPEN (UNIT=1,FILE='TIMEFREQ.IN',STATUS='OLD')
      OPEN (UNIT=2,FILE='MTLFREQ.DAT',STATUS='OLD')
      OPEN (UNIT=3,FILE='TIMEFREQ.OUT',STATUS='UNKNOWN')
      PI=4.*ATAN(1.)
      RADEG=180./PI
      READ(1,*) NFREQ
      READ(1,*) VMAX
      READ(1,*) F0
      READ(1,*) DC
      READ(1,*) TR
      READ(1,*) FINAL
      READ(1,*) DT
      READ(1,*) DCOUT
      TIME=0.
      NTOT=IFIX(FINAL/DT)
      DO 1 I=1,NFREQ
      READ(2,*) FREQ(I),VM(I),VP(I)
1     CONTINUE
C*****COMPUTE FOURIER COEFFICIENTS FOR TRAPEZOIDAL PULSE TRAIN*********
      DO 2 I=1,NFREQ
      FI=FLOAT(I)
      CN(I)=2.*VMAX*DC*SIN(FI*DC*PI)*SIN(FI*PI*TR*F0)/((FI*PI*DC)
     1*(FI*PI*TR*F0))
      ANG(I)=-FI*PI*DC-FI*PI*TR*F0
2     CONTINUE
C**********************************************************************
      WRITE(3,6) NFREQ,VMAX,F0,DC,TR,FINAL,DT,DCOUT
6     FORMAT('NUMBER OF HARMONICS=',I6/
     1'PEAK PULSE LEVEL=',1PE10.3/
     2'REPETITION FREQUENCY OF PULSE TRAIN=',1PE10.3/
     3'DUTY CYCLE OF PULSE TRAIN=',1PE10.3/
     4'PULSE RISE/FALL TIME=',1PE10.3/
     5'FINAL SOLUTION TIME=',1PE10.3/
     6'SOLUTION TIME INCREMENT=',1PE10.3/
     7'DC COMMON-IMPEDANCE COUPLING LEVEL=',1PE10.3/)
      WRITE(3,7)
7     FORMAT(4X,'TIME',13X,'VOUT')
      DO 3 I=1,NTOT
      SUM=0.
      DO 4 N=1,NFREQ
      FN=FLOAT(N)
      SUM=SUM+CN(N)*VM(N)*COS(2.*PI*FN*F0*TIME+VP(N)/RADEG+ANG(N))
4     CONTINUE
      VOUT=SUM+DCOUT*VMAX*DC
      WRITE(3,5) TIME,VOUT
5     FORMAT(2(1PE12.5,5X))
      TIME=TIME+DT
3     CONTINUE
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      STOP
      END
