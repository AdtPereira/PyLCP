C******PROGRAM SPICELC.FOR*********************************************
C
C  Generates a SPICE subcircuit model of a Three-Conductor Transmission
C  Line using the Inductive-Capacitive Coupling Model.  Model is valid
C  for lines where (1) the line is weakly coupled, (2) the line is
C  electrically short at the frequency of interest, and (3) the
C  frequency is sufficiently small.  The input parameters are (1) the
C  line length, and (2) the per-unit-length parameters of mutual
C  inductance, lm, and mutual capacitance, cm.  These input parameters
C  are stored in SPICELC.IN.
C
C  The SPICE subcircuit is stored in SPICELC.OUT as output with the
C  following port node numbering.  Node S is the input end of the
C  Generator Circuit, node L is the output end of the Generator
C  Circuit, node NE is the input end of the Receptor Circuit, and
C  node FE is the output end of the Receptor Circuit.
C
C          Written by: Clayton R. Paul
C                      Department of Electrical Engineering
C                      University of Kentucky
C                      Lexington, KY 40506
C
C**********************************************************************
      REAL LEN,LM
      OPEN(UNIT=5,FILE='SPICELC.IN',STATUS='OLD')
      OPEN(UNIT=6,FILE='SPICELC.OUT',STATUS='UNKNOWN')
C
C*****READ INPUT DATA FROM SPICELC.IN
C
      READ(5,*) LEN
      READ(5,*) LM
      READ(5,*) CM
      WRITE(6,130) LEN,LM,CM
130   FORMAT('*SUBCIRCUIT MODEL OF A THREE-CONDUCTOR TRANSMISSION LINE*
     1'/'*INDUCTIVE-CAPACITIVE COUPLING MODEL*'/
     2'*      TOTAL LINE LENGTH (METERS)=',1PE12.5/
     3'*      PER-UNIT-LENGTH MUTUAL INDUCTANCE (H/m)=',1PE12.5/
     4'*      PER-UNIT-LENGTH MUTUAL CAPACITANCE (F/m)=',1PE12.5)
      LM=LM*LEN
      CM=CM*LEN
C
C*****PRINT SUBCKT CARDS TO SPICELC.OUT
C
      WRITE(6,103) LM,CM
103   FORMAT('.SUBCKT LC S NE L FE'/
     1'VG S L'/
     2'F1 0 1 VG 1'/
     3'E1 2 0 S 0 1'/
     4'VX 2 3'/
     5'LM 1 0',2X,1PE12.5/
     6'CM 3 0',2X,1PE12.5/
     7'E2 NE FE 1 0 1'/
     8'F2 0 FE VX 1'/
     9'.ENDS     LC')
      CLOSE(5)
      CLOSE(6)
      STOP
      END
