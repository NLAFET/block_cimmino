C *******************************************************************
C COPYRIGHT (c) 2002 ENSEEIHT-IRIT, Toulouse, France and
C  Council for the Central Laboratory of the Research Councils.
C All rights reserved.
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
C SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
C
C Please note that for an ACADEMIC Licence:
C
C 1. The Packages may only be used for academic research or teaching
C    purposes by the Licensee, and must not be copied by the Licensee for
C    use by any other persons. Use of the Packages in any commercial
C    application shall be subject to prior written agreement between
C    Hyprotech UK Limited and the Licensee on suitable terms and
C    conditions, which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C  Version 1.0.0 July 2004
C  Version 1.0.1 March 2008  Comments reflowed with length < 73
C  AUTHOR Daniel Ruiz (Daniel.Ruiz@enseeiht.fr)
C *** Copyright (c) 2004  Council for the Central Laboratory of the
C     Research Councils and Ecole Nationale Superieure
C     d'Electrotechnique, d'Electronique, d'Informatique,
C     d'Hydraulique et des Telecommunications de Toulouse.          ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC77 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***

C**********************************************************************
      SUBROUTINE MC77ID(ICNTL, CNTL)
      INTEGER LICNTL, LCNTL
      PARAMETER ( LICNTL=10, LCNTL=10 )
      INTEGER ICNTL(LICNTL)
      DOUBLE PRECISION CNTL(LCNTL)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = -1
      ICNTL(4) = 0
      ICNTL(5) = 0
      ICNTL(6) = 0
      ICNTL(7) = 10
      DO 10 I = 8,LICNTL
        ICNTL(I) = 0
   10 CONTINUE
      CNTL(1) = ZERO
      CNTL(2) = ONE
      DO 20 I = 3,LCNTL
        CNTL(I) = ZERO
   20 CONTINUE
      RETURN
      END
C**********************************************************************
C***           DRIVER FOR THE HARWELL-BOEING SPARSE FORMAT          ***
C**********************************************************************
      SUBROUTINE MC77AD(JOB,M,N,NNZ,JCST,IRN,A,IW,LIW,DW,LDW,
     &                  ICNTL,CNTL,INFO,RINFO)
      INTEGER LICNTL, LCNTL, LINFO, LRINFO
      PARAMETER ( LICNTL=10, LCNTL=10, LINFO=10, LRINFO=10 )
      INTEGER ICNTL(LICNTL),INFO(LINFO)
      DOUBLE PRECISION CNTL(LCNTL),RINFO(LRINFO)
      INTEGER JOB,M,N,NNZ,LIW,LDW
      INTEGER JCST(N+1),IRN(NNZ),IW(LIW)
      DOUBLE PRECISION A(NNZ),DW(LDW)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      DOUBLE PRECISION THRESH, DP
      INTEGER I, J, K, MAXIT, CHECK, SETUP
      EXTERNAL MC77ND,MC77OD,MC77PD,MC77QD
      INTRINSIC ABS,MAX
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      RINFO(1) = ZERO
      RINFO(2) = ZERO
      IF (JOB.LT.-1) THEN
        INFO(1) = -12
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
Cnew  IF (ICNTL(7).LT.0) THEN
      IF (ICNTL(7).LE.0) THEN
        INFO(1) = -10
        INFO(2) = 7
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9002) INFO(1),7,ICNTL(7)
        GO TO 99
      ENDIF
      IF ((JOB.EQ.-1) .AND. (CNTL(2).LT.ONE)) THEN
        INFO(1) = -11
        INFO(2) = 2
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9008) INFO(1),2,CNTL(2)
        GO TO 99
      ENDIF
      IF (M.LT.1) THEN
        INFO(1) = -1
        INFO(2) = M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
        GO TO 99
      ENDIF
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
      IF ( (ICNTL(6).NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -3
        INFO(2) = N-M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9003) INFO(1),(N-M)
        GO TO 99
      ENDIF
      IF ( (JOB.NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -13
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9009) INFO(1),JOB,M,N
        GO TO 99
      ENDIF
      IF (NNZ.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NNZ
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NNZ',NNZ
        GO TO 99
      ENDIF
      K = M
      IF (ICNTL(6).EQ.0)  K = M+N
      IF (LIW.LT.K) THEN
        INFO(1) = -5
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      IF ( (ICNTL(5).EQ.0) .OR. (JOB.GE.2) .OR. (JOB.EQ.-1) )
     &   K = K + NNZ
      IF (LDW.LT.K) THEN
        INFO(1) = -6
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
        GO TO 99
      ENDIF
      IF (ICNTL(4).EQ.0) THEN
        DO 3 I = 1,M
          IW(I) = 0
    3   CONTINUE
        DO 5 J = 1,N
          DO 4 K = JCST(J),JCST(J+1)-1
            I = IRN(K)
            IF (I.LT.1 .OR. I.GT.M) THEN
              INFO(1) = -7
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
            IF (ICNTL(6).NE.0 .AND. I.LT.J) THEN
              INFO(1) = -9
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
            IF (IW(I).EQ.J) THEN
              INFO(1) = -8
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
    4     CONTINUE
    5   CONTINUE
      ENDIF
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,M,N,NNZ,ICNTL(7),CNTL(1)
        IF (ICNTL(9).LT.1) THEN
          WRITE(ICNTL(3),9021) (JCST(J),J=1,N+1)
          WRITE(ICNTL(3),9023) (IRN(J),J=1,NNZ)
          WRITE(ICNTL(3),9024) (A(J),J=1,NNZ)
        ENDIF
      ENDIF
      DO 10 I=1,10
        INFO(I)  = 0
        RINFO(I) = ZERO
   10 CONTINUE
      THRESH = MAX( ZERO, CNTL(1) )
      MAXIT  = ICNTL(7)
Cnew  IF ( (MAXIT.LE.0) .AND. (CNTL(1).LE.ZERO) )  THEN
Cnew     INFO(1) = -14
Cnew     IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9010) INFO(1),MAXIT,CNTL(1)
Cnew     GO TO 99
Cnew  ENDIF
      IF (ICNTL(10).EQ.1) THEN
          DO I=1,M
            IF(DW(I).LE.ZERO) THEN
                INFO(1) = -14
                GOTO 99
            ENDIF
          ENDDO
          DO J=1,N
            IF(DW(J+M).LE.ZERO) THEN
                INFO(1) = -15
                GOTO 99
            ENDIF
          ENDDO
      ELSE
          DO I=1, M
            DW(I) = ONE
          ENDDO
          DO J=1, N
            DW(M+J) = ONE
          ENDDO
      ENDIF

      CHECK  = 1
      IF (CNTL(1).LE.ZERO)  CHECK = 0
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      SETUP = 0
      IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
        SETUP = 1
        IF (JOB.EQ.-1) THEN
          DP = CNTL(2)
        ELSE
          DP = REAL(JOB)
        ENDIF
        IF(ICNTL(10).EQ.1) THEN
          IF(ICNTL(6).EQ.0) THEN
            DO 18 J = 1,M
              DW(J) = DW(J)**DP
   18       CONTINUE
            DO 19 J = 1,N
              DW(M+J) = DW(M+J)**DP
   19       CONTINUE
          ELSE
            DO 17 J = 1,2*M
              DW(J) = DW(J)**DP
   17       CONTINUE
          ENDIF
        ENDIF
        DO 20 J = 1,NNZ
          DW(K+J) = ABS(A(J))**DP
   20   CONTINUE
        THRESH = ONE - (ABS(ONE-THRESH))**DP
      ELSE
        IF (ICNTL(5).EQ.0) THEN
          SETUP = 1
          DO 30 J = 1,NNZ
            DW(K+J) = ABS(A(J))
   30     CONTINUE
        ENDIF
      ENDIF
      
      IF (ICNTL(6).EQ.0)  THEN
        IF (JOB.EQ.0) THEN
          IF (SETUP.EQ.1) THEN
            CALL MC77ND(M,N,NNZ,JCST,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77ND(M,N,NNZ,JCST,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77OD(M,N,NNZ,JCST,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77OD(M,N,NNZ,JCST,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 40 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+N+I)**DP) )
   40       CONTINUE
            RINFO(2) = ZERO
            DO 50 J = 1,N
              DW(M+J) = DW(M+J)**DP
              IF (IW(M+J).NE.0)
     &           RINFO(2) = MAX( RINFO(2), ABS(ONE-DW(2*M+N+J)**DP) )
   50       CONTINUE
          ENDIF
        ENDIF
      ELSE
        IF (JOB.EQ.0) THEN
          IF (SETUP.EQ.1) THEN
            CALL MC77PD(M,NNZ,JCST,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77PD(M,NNZ,JCST,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77QD(M,NNZ,JCST,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77QD(M,NNZ,JCST,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 60 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+I)**DP) )
   60       CONTINUE
          ENDIF
        ENDIF
        RINFO(2) = RINFO(1)
      ENDIF
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,3),(RINFO(J),J=1,2)
        IF (ICNTL(9).LT.2) THEN
          WRITE(ICNTL(3),9031) (DW(I),I=1,M)
          IF (ICNTL(6).EQ.0)
     &      WRITE(ICNTL(3),9032) (DW(M+J),J=1,N)
        ENDIF
      ENDIF
   99 CONTINUE
      RETURN
 9001 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3,
     &        ' because ',(A),' = ',I10)
 9002 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Bad input control flag.',
     &        ' Value of ICNTL(',I1,') = ',I8)
 9003 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Input matrix is symmetric and N /= M.',
     &        ' Value of (N-M) = ',I8)
 9004 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Column ',I8,
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Column ',I8,
     &        ' contains two or more entries with row index ',I8)
 9008 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Bad input REAL control parameter.',
     &        ' Value of CNTL(',I1,') = ',1PD14.4)
 9009 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Only scaling in infinity norm is allowed',
     &        ' for rectangular matrices'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8)
C9010 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
Cnew &        '        Algorithm will loop indefinitely
Cnew &        '     Check for inconsistency'/
Cnew &        ' Maximum number of iterations = ',I8/
Cnew &        ' Threshold for convergence    = ',1PD14.4)
 9020 FORMAT (' ****** Input parameters for MC77A/AD:'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8/' NNZ = ',I8/
     &        ' Max n.b. Iters.  = ',I8/
     &        ' Cvgce. Threshold = ',1PD14.4)
 9021 FORMAT (' JCST(1:N+1) = ',8I8/(15X,8I8))
 9023 FORMAT (' IRN(1:NNZ)  = ',8I8/(15X,8I8))
 9024 FORMAT (' A(1:NNZ)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC77A/AD:'/
     &        ' INFO(1:3)   = ',3(2X,I8)/
     &        ' RINFO(1:2)  = ',2(1PD14.4))
 9031 FORMAT (' DW(1:M)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9032 FORMAT (' DW(M+1:M+N) = ',4(1PD14.4)/(15X,4(1PD14.4)))
c9031 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
c9032 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END
C**********************************************************************
      SUBROUTINE MC77ND(M,N,NNZ,JCST,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IW(M),JW(N)
      DOUBLE PRECISION THRESH,ERR(2)
      DOUBLE PRECISION A(NNZ),D(M),E(N),DW(M),EW(N)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
      INTRINSIC SQRT, MAX, ABS
      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO
      DO 5 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
!        D(I)  = ONE
    5 CONTINUE
      DO 10 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
!        E(J)  = ONE
   10 CONTINUE
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          I = IRN(K)
          IF (EW(J).LT.A(K)) THEN
            EW(J) = A(K)
            JW(J) = I
          ENDIF
          IF (DW(I).LT.A(K)) THEN
            DW(I) = A(K)
            IW(I) = J
          ENDIF
   30   CONTINUE
   20 CONTINUE
      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE
      DO 50 J=1,N
        I = JW(J)
        IF (I.GT.0) THEN
          IF (IW(I).EQ.J) THEN
            IW(I) = -IW(I)
            JW(J) = -JW(J)
          ENDIF
        ENDIF
   50 CONTINUE
      DO 60 I=1,M
        IF (IW(I).GT.0)  GOTO 99
   60 CONTINUE
      DO 65 J=1,N
        IF (JW(J).GT.0)  GOTO 99
   65 CONTINUE
      GOTO 200
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF
        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE
        DO 120 J=1,N
          IF (JW(J).GT.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (D(I)*E(J))
              IF (EW(J).LT.S) THEN
                EW(J) = S
                JW(J) = I
              ENDIF
              IF (IW(I).GT.0) THEN
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  130       CONTINUE
          ELSE
            DO 140 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              IF (IW(I).GT.0) THEN
                S = A(K) / (D(I)*E(J))
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  140       CONTINUE
          ENDIF
  120   CONTINUE
        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE
        DO 160 J=1,N
          I = JW(J)
          IF (I.GT.0) THEN
            IF (IW(I).EQ.J) THEN
              IW(I) = -IW(I)
              JW(J) = -JW(J)
            ENDIF
          ENDIF
  160   CONTINUE
        IF (CHECK.LE.0) GOTO 99
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) )  GOTO 200
      IF (CHECK.GT.0)  GOTO 99
  200 RETURN
      END
C**********************************************************************
      SUBROUTINE MC77OD(M,N,NNZ,JCST,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IW(M),JW(N)
      DOUBLE PRECISION THRESH,ERR(2)
      DOUBLE PRECISION A(NNZ),D(M),E(N),DW(M),EW(N)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
      INTRINSIC SQRT, MAX, ABS
      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO
      DO 10 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
!        D(I)  = ONE
   10 CONTINUE
      DO 15 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
!        E(J)  = ONE
   15 CONTINUE
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          IF (A(K).GT.ZERO) THEN
            I = IRN(K)
            EW(J) = EW(J) + A(K)
            IF (JW(J).EQ.0) THEN
               JW(J) = K
            ELSE
               JW(J) = -1
            ENDIF
            DW(I) = DW(I) + A(K)
            IF (IW(I).EQ.0) THEN
               IW(I) = K
            ELSE
               IW(I) = -1
            ENDIF
          ENDIF
   30   CONTINUE
   20 CONTINUE
      DO 40 K=1,M
        IF (IW(K).NE.0) D(K) = SQRT(DW(K))
   40 CONTINUE
      DO 45 K=1,N
        IF (JW(K).NE.0) E(K) = SQRT(EW(K))
   45 CONTINUE
      DO 50 J=1,N
        K = JW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IW(I).EQ.K) THEN
            IW(I) = 0
            JW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE
      DO 60 K=1,M
        IF ( (IW(K).NE.0) .OR. (JW(K).NE.0) )  GOTO 99
   60 CONTINUE
Cnew  DO 60 I=1,M
Cnew    IF (IW(I).NE.0)  GOTO 99
Cne60 CONTINUE
Cnew  DO 65 J=1,N
Cnew    IF (JW(J).NE.0)  GOTO 99
Cne65 CONTINUE
      GOTO 200
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF
        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE
        DO 120 J=1,N
          IF (JW(J).NE.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (D(I)*E(J))
              EW(J) = EW(J) + S
              DW(I) = DW(I) + S
  130       CONTINUE
          ENDIF
  120   CONTINUE
        DO 150 I=1,M
          IF (IW(I).NE.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).NE.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE
        IF (CHECK.LE.0) GOTO 99
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).NE.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).NE.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) ) GOTO 200
      IF (CHECK.GT.0)  GOTO 99
  200 RETURN
      END
C**********************************************************************
      SUBROUTINE MC77PD(N,NNZ,JCST,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IJW(N)
      DOUBLE PRECISION THRESH,ERR
      DOUBLE PRECISION A(NNZ),DE(N),DEW(N)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
      INTRINSIC SQRT, MAX, ABS
      INFO = 0
      NITER = 0
      ERR = ZERO
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          I = IRN(K)
          IF (DEW(J).LT.A(K)) THEN
            DEW(J) = A(K)
            IJW(J) = I
          ENDIF
          IF (DEW(I).LT.A(K)) THEN
            DEW(I) = A(K)
            IJW(I) = J
          ENDIF
   30   CONTINUE
   20 CONTINUE
      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE
      DO 50 J=1,N
        I = IJW(J)
        IF (I.GT.0) THEN
          IF (IJW(I).EQ.J) THEN
            IJW(I) = -IJW(I)
            IF (I.NE.J) IJW(J) = -IJW(J)
          ENDIF
        ENDIF
   50 CONTINUE
      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF
        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE
        DO 120 J=1,N
          IF (IJW(J).GT.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (DE(I)*DE(J))
              IF (DEW(J).LT.S) THEN
                DEW(J) = S
                IJW(J) = I
              ENDIF
              IF (IJW(I).GT.0) THEN
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
  130       CONTINUE
          ELSE
            DO 140 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              IF (IJW(I).GT.0) THEN
                S = A(K) / (DE(I)*DE(J))
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
  140       CONTINUE
          ENDIF
  120   CONTINUE
        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE
        DO 160 J=1,N
          I = IJW(J)
          IF (I.GT.0) THEN
            IF (IJW(I).EQ.J) THEN
              IJW(I) = -IJW(I)
              IF (I.NE.J) IJW(J) = -IJW(J)
            ENDIF
          ENDIF
  160   CONTINUE
        IF (CHECK.LE.0) GOTO 99
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200
      IF (CHECK.GT.0)  GOTO 99
  200 RETURN
      END
C**********************************************************************
      SUBROUTINE MC77QD(N,NNZ,JCST,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IJW(N)
      DOUBLE PRECISION THRESH,ERR
      DOUBLE PRECISION A(NNZ),DE(N),DEW(N)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
      INTRINSIC SQRT, MAX, ABS
      INFO = 0
      NITER = 0
      ERR = ZERO
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          IF (A(K).GT.ZERO) THEN
            DEW(J) = DEW(J) + A(K)
            IF (IJW(J).EQ.0) THEN
               IJW(J) = K
            ELSE
               IJW(J) = -1
            ENDIF
            I = IRN(K)
            IF (I.NE.J) THEN
              DEW(I) = DEW(I) + A(K)
              IJW(I) = IJW(I) + 1
              IF (IJW(I).EQ.0) THEN
                 IJW(I) = K
              ELSE
                 IJW(I) = -1
              ENDIF
            ENDIF
          ENDIF
   30   CONTINUE
   20 CONTINUE
      DO 40 K=1,N
        IF (IJW(K).NE.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE
      DO 50 J=1,N
        K = IJW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IJW(I).EQ.K) THEN
            IJW(I) = 0
            IJW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE
      DO 60 K=1,N
        IF (IJW(K).NE.0)  GOTO 99
   60 CONTINUE
      GOTO 200
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF
        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE
        DO 120 J=1,N
          IF (IJW(J).NE.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (DE(I)*DE(J))
              DEW(J) = DEW(J) + S
              IF (I.NE.J)  DEW(I) = DEW(I) + S
  130       CONTINUE
          ENDIF
  120   CONTINUE
        DO 150 K=1,N
          IF (IJW(K).NE.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE
        IF (CHECK.LE.0) GOTO 99
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).NE.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200
      IF (CHECK.GT.0)  GOTO 99
  200 RETURN
      END
C**********************************************************************
C***              DRIVER FOR THE GENERAL SPARSE FORMAT              ***
C**********************************************************************
      SUBROUTINE MC77BD(JOB,M,N,NNZ,IRN,JCN,A,IW,LIW,DW,LDW,
     &                  ICNTL,CNTL,INFO,RINFO)
      INTEGER LICNTL, LCNTL, LINFO, LRINFO
      PARAMETER ( LICNTL=10, LCNTL=10, LINFO=10, LRINFO=10 )
      INTEGER ICNTL(LICNTL),INFO(LINFO)
      DOUBLE PRECISION CNTL(LCNTL),RINFO(LRINFO)
      INTEGER JOB,M,N,NNZ,LIW,LDW
      INTEGER JCN(NNZ),IRN(NNZ),IW(LIW)
      DOUBLE PRECISION A(NNZ),DW(LDW)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      DOUBLE PRECISION THRESH, DP
      INTEGER I, J, K, MAXIT, CHECK, SETUP
      EXTERNAL MC77RD,MC77SD,MC77TD,MC77UD
      INTRINSIC ABS,MAX
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      RINFO(1) = ZERO
      RINFO(2) = ZERO
      IF (JOB.LT.-1) THEN
        INFO(1) = -12
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
Cnew  IF (ICNTL(7).LT.0) THEN
      IF (ICNTL(7).LE.0) THEN
        INFO(1) = -10
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9002) INFO(1),7,ICNTL(7)
        GO TO 99
      ENDIF
      IF ((JOB.EQ.-1) .AND. (CNTL(2).LT.ONE)) THEN
        INFO(1) = -11
        INFO(2) = 2
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9008) INFO(1),2,CNTL(2)
        GO TO 99
      ENDIF
      IF (M.LT.1) THEN
        INFO(1) = -1
        INFO(2) = M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
        GO TO 99
      ENDIF
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
      IF ( (ICNTL(6).NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -3
        INFO(2) = N-M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9003) INFO(1),(N-M)
        GO TO 99
      ENDIF
      IF ( (JOB.NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -13
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9009) INFO(1),JOB,M,N
        GO TO 99
      ENDIF
      IF (NNZ.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NNZ
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NNZ',NNZ
        GO TO 99
      ENDIF
      K = M
      IF (ICNTL(6).EQ.0)  K = M+N
      IF (LIW.LT.K) THEN
        INFO(1) = -5
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      IF ( (ICNTL(5).EQ.0) .OR. (JOB.GE.2) .OR. (JOB.EQ.-1) )
     &   K = K + NNZ
      IF (LDW.LT.K) THEN
        INFO(1) = -6
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
        GO TO 99
      ENDIF
      IF (ICNTL(4).EQ.0) THEN
        DO 6 K = 1,NNZ
          I = IRN(K)
          J = JCN(K)
          IF (I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N) THEN
            INFO(1) = -7
            INFO(2) = K
            IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),K,I,J
            GO TO 99
          ENDIF
          IF (ICNTL(6).NE.0 .AND. I.LT.J) THEN
            INFO(1) = -9
            INFO(2) = K
            IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),K,I,J
            GO TO 99
          ENDIF
    6   CONTINUE
        DO 7 I = 1,M
          IW(I) = 0
    7   CONTINUE
        J = 1
        DO 9 K = 1,NNZ
          IF (JCN(K).EQ.J) THEN
            I = IRN(K)
            IF (IW(I).EQ.J) THEN
              INFO(1) = -8
              INFO(2) = K
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),K,I,J
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
            J = J + 1
          ENDIF
    9   CONTINUE
C       DO 9 J = 1,N
C         DO 8 K = 1,NNZ
C         IF (JCN(K).EQ.J) THEN
C           I = IRN(K)
C           IF (IW(I).EQ.J) THEN
C             INFO(1) = -8
C             INFO(2) = K
C             IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),K,I,J
C             GO TO 99
C           ELSE
C             IW(I) = J
C           ENDIF
C         ENDIF
C   8     CONTINUE
C   9   CONTINUE
      ENDIF
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,M,N,NNZ,ICNTL(7),CNTL(1)
        IF (ICNTL(9).LT.1) THEN
          WRITE(ICNTL(3),9022) (JCN(J),J=1,NNZ)
          WRITE(ICNTL(3),9023) (IRN(J),J=1,NNZ)
          WRITE(ICNTL(3),9024) (A(J),J=1,NNZ)
        ENDIF
      ENDIF
      DO 10 I=1,10
        INFO(I)  = 0
        RINFO(I) = ZERO
   10 CONTINUE
      THRESH = MAX( ZERO, CNTL(1) )
      MAXIT  = ICNTL(7)
Cnew  IF ( (MAXIT.LE.0) .AND. (CNTL(1).LE.ZERO) )  THEN
Cnew     INFO(1) = -14
Cnew     IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9010) INFO(1),MAXIT,CNTL(1)
Cnew     GO TO 99
Cnew  ENDIF
      CHECK  = 1
      IF (CNTL(1).LE.ZERO)  CHECK = 0
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      SETUP = 0
      IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
        SETUP = 1
        IF (JOB.EQ.-1) THEN
          DP = CNTL(2)
        ELSE
          DP = REAL(JOB)
        ENDIF
        DO 20 J = 1,NNZ
          DW(K+J) = ABS(A(J))**DP
   20   CONTINUE
        THRESH = ONE - (ABS(ONE-THRESH))**DP
      ELSE
        IF (ICNTL(5).EQ.0) THEN
          SETUP = 1
          DO 30 J = 1,NNZ
            DW(K+J) = ABS(A(J))
   30     CONTINUE
        ENDIF
      ENDIF
      IF (ICNTL(6).EQ.0)  THEN
        IF (JOB.EQ.0) THEN
          IF (SETUP.EQ.1) THEN
            CALL MC77RD(M,N,NNZ,JCN,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77RD(M,N,NNZ,JCN,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77SD(M,N,NNZ,JCN,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77SD(M,N,NNZ,JCN,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 40 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+N+I)**DP) )
   40       CONTINUE
            RINFO(2) = ZERO
            DO 50 J = 1,N
              DW(M+J) = DW(M+J)**DP
              IF (IW(M+J).NE.0)
     &           RINFO(2) = MAX( RINFO(2), ABS(ONE-DW(2*M+N+J)**DP) )
   50       CONTINUE
          ENDIF
        ENDIF
      ELSE
        IF (JOB.EQ.0) THEN
          IF (SETUP.EQ.1) THEN
            CALL MC77TD(M,NNZ,JCN,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77TD(M,NNZ,JCN,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77UD(M,NNZ,JCN,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77UD(M,NNZ,JCN,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 60 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+I)**DP) )
   60       CONTINUE
          ENDIF
        ENDIF
        RINFO(2) = RINFO(1)
      ENDIF
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,3),(RINFO(J),J=1,2)
        IF (ICNTL(9).LT.2) THEN
          WRITE(ICNTL(3),9031) (DW(I),I=1,M)
          IF (ICNTL(6).EQ.0)
     &      WRITE(ICNTL(3),9032) (DW(M+J),J=1,N)
        ENDIF
      ENDIF
   99 CONTINUE
      RETURN
 9001 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3,
     &        ' because ',(A),' = ',I10)
 9002 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Bad input control flag.',
     &        ' Value of ICNTL(',I1,') = ',I8)
 9003 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Input matrix is symmetric and N /= M.',
     &        ' Value of (N-M) = ',I8)
 9004 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Entry ',I8,
     &        ' has invalid row index ',I8, ' or column index ',I8)
 9007 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Duplicate entry ',I8, '   with row index ',I8/
     &        '                                 and column index ',I8)
 9008 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Bad input REAL control parameter.',
     &        ' Value of CNTL(',I1,') = ',1PD14.4)
 9009 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Only scaling in infinity norm is allowed',
     &        ' for rectangular matrices'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8)
C9010 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
Cnew &        '        Algorithm will loop indefinitely
Cnew &        '     Check for inconsistency'/
Cnew &        ' Maximum number of iterations = ',I8/
Cnew &        ' Threshold for convergence    = ',1PD14.4)
 9020 FORMAT (' ****** Input parameters for MC77B/BD:'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8/' NNZ = ',I8/
     &        ' Max n.b. Iters.  = ',I8/
     &        ' Cvgce. Threshold = ',1PD14.4)
 9022 FORMAT (' JCN(1:NNZ)  = ',8I8/(15X,8I8))
 9023 FORMAT (' IRN(1:NNZ)  = ',8I8/(15X,8I8))
 9024 FORMAT (' A(1:NNZ)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC77B/BD:'/
     &        ' INFO(1:3)   = ',3(2X,I8)/
     &        ' RINFO(1:2)  = ',2(1PD14.4))
 9031 FORMAT (' DW(1:M)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9032 FORMAT (' DW(M+1:M+N) = ',4(1PD14.4)/(15X,4(1PD14.4)))
c9031 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
c9032 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END
C**********************************************************************
      SUBROUTINE MC77RD(M,N,NNZ,JCN,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IW(M),JW(N)
      DOUBLE PRECISION THRESH,ERR(2)
      DOUBLE PRECISION A(NNZ),D(M),E(N),DW(M),EW(N)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
      INTRINSIC SQRT, MAX, ABS
      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO
      DO 5 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
    5 CONTINUE
      DO 10 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   10 CONTINUE
      DO 20 K=1,NNZ
        I = IRN(K)
        J = JCN(K)
        IF (EW(J).LT.A(K)) THEN
          EW(J) = A(K)
          JW(J) = I
        ENDIF
        IF (DW(I).LT.A(K)) THEN
          DW(I) = A(K)
          IW(I) = J
        ENDIF
   20 CONTINUE
      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE
      DO 50 J=1,N
        I = JW(J)
        IF (I.GT.0) THEN
          IF (IW(I).EQ.J) THEN
            IW(I) = -IW(I)
            JW(J) = -JW(J)
          ENDIF
        ENDIF
   50 CONTINUE
      DO 60 I=1,M
        IF (IW(I).GT.0)  GOTO 99
   60 CONTINUE
      DO 65 J=1,N
        IF (JW(J).GT.0)  GOTO 99
   65 CONTINUE
      GOTO 200
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF
        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE
        DO 120 K=1,NNZ
          I = IRN(K)
          J = JCN(K)
          IF (JW(J).GT.0) THEN
            S = A(K) / (D(I)*E(J))
            IF (EW(J).LT.S) THEN
              EW(J) = S
              JW(J) = I
            ENDIF
            IF (IW(I).GT.0) THEN
              IF (DW(I).LT.S) THEN
                DW(I) = S
                IW(I) = J
              ENDIF
            ENDIF
          ELSE
            IF (IW(I).GT.0) THEN
              S = A(K) / (D(I)*E(J))
              IF (DW(I).LT.S) THEN
                DW(I) = S
                IW(I) = J
              ENDIF
            ENDIF
          ENDIF
  120   CONTINUE
        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE
        DO 160 J=1,N
          I = JW(J)
          IF (I.GT.0) THEN
            IF (IW(I).EQ.J) THEN
              IW(I) = -IW(I)
              JW(J) = -JW(J)
            ENDIF
          ENDIF
  160   CONTINUE
        IF (CHECK.LE.0) GOTO 99
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) )  GOTO 200
      IF (CHECK.GT.0)  GOTO 99
  200 RETURN
      END
C**********************************************************************
      SUBROUTINE MC77SD(M,N,NNZ,JCN,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IW(M),JW(N)
      DOUBLE PRECISION THRESH,ERR(2)
      DOUBLE PRECISION A(NNZ),D(M),E(N),DW(M),EW(N)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
      INTRINSIC SQRT, MAX, ABS
      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO
      DO 10 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
   10 CONTINUE
      DO 15 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   15 CONTINUE
      DO 20 K=1,NNZ
        IF (A(K).GT.ZERO) THEN
          I = IRN(K)
          J = JCN(K)
          EW(J) = EW(J) + A(K)
          IF (JW(J).EQ.0) THEN
             JW(J) = K
          ELSE
             JW(J) = -1
          ENDIF
          DW(I) = DW(I) + A(K)
          IF (IW(I).EQ.0) THEN
             IW(I) = K
          ELSE
             IW(I) = -1
          ENDIF
        ENDIF
   20 CONTINUE
      DO 40 I=1,M
        IF (IW(I).NE.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).NE.0) E(J) = SQRT(EW(J))
   45 CONTINUE
      DO 50 J=1,N
        K = JW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IW(I).EQ.K) THEN
            IW(I) = 0
            JW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE
      DO 60 K=1,M
        IF ( (IW(K).NE.0) .OR. (JW(K).NE.0) )  GOTO 99
   60 CONTINUE
Cnew  DO 60 I=1,M
Cnew    IF (IW(I).NE.0)  GOTO 99
Cne60 CONTINUE
Cnew  DO 65 J=1,N
Cnew    IF (JW(J).NE.0)  GOTO 99
Cne65 CONTINUE
      GOTO 200
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF
        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE
        DO 120 K=1,NNZ
          J = JCN(K)
          I = IRN(K)
          S = A(K) / (D(I)*E(J))
          EW(J) = EW(J) + S
          DW(I) = DW(I) + S
  120   CONTINUE
        DO 150 I=1,M
          IF (IW(I).NE.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).NE.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE
        IF (CHECK.LE.0) GOTO 99
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).NE.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).NE.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) ) GOTO 200
      IF (CHECK.GT.0)  GOTO 99
  200 RETURN
      END
C**********************************************************************
      SUBROUTINE MC77TD(N,NNZ,JCN,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IJW(N)
      DOUBLE PRECISION THRESH,ERR
      DOUBLE PRECISION A(NNZ),DE(N),DEW(N)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
      INTRINSIC SQRT, MAX, ABS
      INFO = 0
      NITER = 0
      ERR = ZERO
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE
      DO 20 K=1,NNZ
        I = IRN(K)
        J = JCN(K)
        IF (DEW(J).LT.A(K)) THEN
          DEW(J) = A(K)
          IJW(J) = I
        ENDIF
        IF (DEW(I).LT.A(K)) THEN
          DEW(I) = A(K)
          IJW(I) = J
        ENDIF
   20 CONTINUE
      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE
      DO 50 J=1,N
        I = IJW(J)
        IF (I.GT.0) THEN
          IF (IJW(I).EQ.J) THEN
            IJW(I) = -IJW(I)
            IF (I.NE.J) IJW(J) = -IJW(J)
          ENDIF
        ENDIF
   50 CONTINUE
      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF
        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE
        DO 120 K=1,NNZ
          I = IRN(K)
          J = JCN(K)
          IF (IJW(J).GT.0) THEN
            S = A(K) / (DE(I)*DE(J))
            IF (DEW(J).LT.S) THEN
              DEW(J) = S
              IJW(J) = I
            ENDIF
            IF (IJW(I).GT.0) THEN
              IF (DEW(I).LT.S) THEN
                DEW(I) = S
                IJW(I) = J
              ENDIF
            ENDIF
          ELSE
            IF (IJW(I).GT.0) THEN
              S = A(K) / (DE(I)*DE(J))
              IF (DEW(I).LT.S) THEN
                DEW(I) = S
                IJW(I) = J
              ENDIF
            ENDIF
          ENDIF
  120   CONTINUE
        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE
        DO 160 J=1,N
          I = IJW(J)
          IF (I.GT.0) THEN
            IF (IJW(I).EQ.J) THEN
              IJW(I) = -IJW(I)
              IF (I.NE.J) IJW(J) = -IJW(J)
            ENDIF
          ENDIF
  160   CONTINUE
        IF (CHECK.LE.0) GOTO 99
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200
      IF (CHECK.GT.0)  GOTO 99
  200 RETURN
      END
C**********************************************************************
      SUBROUTINE MC77UD(N,NNZ,JCN,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IJW(N)
      DOUBLE PRECISION THRESH,ERR
      DOUBLE PRECISION A(NNZ),DE(N),DEW(N)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
      INTRINSIC SQRT, MAX, ABS
      INFO = 0
      NITER = 0
      ERR = ZERO
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE
      DO 20 K=1,NNZ
        IF (A(K).GT.ZERO) THEN
          J = JCN(K)
          DEW(J) = DEW(J) + A(K)
          IF (IJW(J).EQ.0) THEN
             IJW(J) = K
          ELSE
             IJW(J) = -1
          ENDIF
          I = IRN(K)
          IF (I.NE.J) THEN
            DEW(I) = DEW(I) + A(K)
            IF (IJW(I).EQ.0) THEN
               IJW(I) = K
            ELSE
               IJW(I) = -1
            ENDIF
          ENDIF
        ENDIF
   20 CONTINUE
      DO 40 K=1,N
        IF (IJW(K).NE.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE
      DO 50 J=1,N
        K = IJW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IJW(I).EQ.K) THEN
            IJW(I) = 0
            IJW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE
      DO 60 K=1,N
        IF (IJW(K).NE.0)  GOTO 99
   60 CONTINUE
      GOTO 200
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF
        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE
        DO 120 K=1,NNZ
          J = JCN(K)
          I = IRN(K)
          S = A(K) / (DE(I)*DE(J))
          DEW(J) = DEW(J) + S
          IF (I.NE.J)  DEW(I) = DEW(I) + S
  120   CONTINUE
        DO 150 K=1,N
          IF (IJW(K).NE.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE
        IF (CHECK.LE.0) GOTO 99
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).NE.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200
      IF (CHECK.GT.0)  GOTO 99
  200 RETURN
      END
C**********************************************************************
C***             DRIVER FOR THE CASE OF DENSE MATRICES              ***
C**********************************************************************
      SUBROUTINE MC77CD(JOB,M,N,A,LDA,IW,LIW,DW,LDW,
     &                  ICNTL,CNTL,INFO,RINFO)
      INTEGER LICNTL, LCNTL, LINFO, LRINFO
      PARAMETER ( LICNTL=10, LCNTL=10, LINFO=10, LRINFO=10 )
      INTEGER ICNTL(LICNTL),INFO(LINFO)
      DOUBLE PRECISION CNTL(LCNTL),RINFO(LRINFO)
      INTEGER JOB,M,N,LDA,LIW,LDW
      INTEGER IW(LIW)
      DOUBLE PRECISION A(LDA,*),DW(LDW)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      DOUBLE PRECISION THRESH, DP
      INTEGER I, J, K, NNZ, MAXIT, CHECK, SETUP
      EXTERNAL MC77JD,MC77KD,MC77LD,MC77MD
      INTRINSIC ABS,MAX
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      RINFO(1) = ZERO
      RINFO(2) = ZERO
      IF (JOB.LT.-1) THEN
        INFO(1) = -12
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
Cnew  IF (ICNTL(7).LT.0) THEN
      IF (ICNTL(7).LE.0) THEN
        INFO(1) = -10
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9002) INFO(1),7,ICNTL(7)
        GO TO 99
      ENDIF
      IF ((JOB.EQ.-1) .AND. (CNTL(2).LT.ONE)) THEN
        INFO(1) = -11
        INFO(2) = 2
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9008) INFO(1),2,CNTL(2)
        GO TO 99
      ENDIF
      IF (M.LT.1) THEN
        INFO(1) = -1
        INFO(2) = M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
        GO TO 99
      ENDIF
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
      IF ( (ICNTL(6).NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -3
        INFO(2) = N-M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9003) INFO(1),(N-M)
        GO TO 99
      ENDIF
      IF ( (JOB.NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -13
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9009) INFO(1),JOB,M,N
        GO TO 99
      ENDIF
      IF (ICNTL(6).EQ.0) THEN
        IF (M.GT.LDA) THEN
          INFO(1) = -4
          INFO(2) = LDA-M
          IF (ICNTL(1).GE.0)
     &      WRITE(ICNTL(1),9001) INFO(1),'LDA < M; (LDA-M)',INFO(2)
          GO TO 99
        ENDIF
      ELSE
        IF (((M*(M+1))/2).GT.LDA) THEN
          INFO(1) = -4
          INFO(2) = LDA-((M*(M+1))/2)
          IF (ICNTL(1).GE.0)
     &      WRITE(ICNTL(1),9001) INFO(1),
     &        'LDA < M*(M+1)/2; (LDA-(M*(M+1)/2))',INFO(2)
          GO TO 99
        ENDIF
      ENDIF
      K = M
      IF (ICNTL(6).EQ.0)  K = M+N
      IF (LIW.LT.K) THEN
        INFO(1) = -5
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
      K = 2*M
      NNZ = (M*(M+1))/2
      IF (ICNTL(6).EQ.0)  THEN
        K = 2*(M+N)
        NNZ = LDA*N
      ENDIF
      IF ( (ICNTL(5).EQ.0) .OR. (JOB.GE.2) .OR. (JOB.EQ.-1) )
     &   K = K + NNZ
      IF (LDW.LT.K) THEN
        INFO(1) = -6
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
        GO TO 99
      ENDIF
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,M,N,ICNTL(7),CNTL(1)
        IF (ICNTL(9).LT.1) THEN
          K = 1
          DO 5 I=1,M
            IF (ICNTL(6).EQ.0)
     &        WRITE(ICNTL(3),9021) I,(A(I,J),J=1,N)
            IF (ICNTL(6).NE.0) THEN
              WRITE(ICNTL(3),9022) I,(A(J,1),J=K,K+M-I)
              K = K+M-I+1
            ENDIF
    5     CONTINUE
        ENDIF
      ENDIF
      DO 10 I=1,10
        INFO(I)  = 0
        RINFO(I) = ZERO
   10 CONTINUE
      THRESH = MAX( ZERO, CNTL(1) )
      MAXIT  = ICNTL(7)
Cnew  IF ( (MAXIT.LE.0) .AND. (CNTL(1).LE.ZERO) )  THEN
Cnew     INFO(1) = -14
Cnew     IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9010) INFO(1),MAXIT,CNTL(1)
Cnew     GO TO 99
Cnew  ENDIF
      CHECK  = 1
      IF (CNTL(1).LE.ZERO)  CHECK = 0
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      SETUP = 0
      IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
        SETUP = 1
        IF (JOB.EQ.-1) THEN
          DP = CNTL(2)
        ELSE
          DP = REAL(JOB)
        ENDIF
        IF (ICNTL(6).EQ.0) THEN
          DO 25 J = 1,N
            DO 20 I = 1,M
              DW(K+(J-1)*LDA+I) = ABS(A(I,J))**DP
   20       CONTINUE
   25     CONTINUE
        ELSE
          NNZ = (M*(M+1))/2
          DO 30 J = 1,NNZ
            DW(K+J) = ABS(A(J,1))**DP
   30     CONTINUE
        ENDIF
        THRESH = ONE - (ABS(ONE-THRESH))**DP
      ELSE
        IF (ICNTL(5).EQ.0) THEN
          SETUP = 1
          IF (ICNTL(6).EQ.0) THEN
            DO 40 J = 1,N
              DO 35 I = 1,M
                DW(K+(J-1)*LDA+I) = ABS(A(I,J))
   35         CONTINUE
   40       CONTINUE
          ELSE
            NNZ = (M*(M+1))/2
            DO 45 J = 1,NNZ
              DW(K+J) = ABS(A(J,1))
   45       CONTINUE
          ENDIF
        ENDIF
      ENDIF
      IF (ICNTL(6).EQ.0)  THEN
        IF (JOB.EQ.0) THEN
          IF (SETUP.EQ.1) THEN
            CALL MC77JD(M,N,DW(2*(M+N)+1),LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77JD(M,N,A,LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77KD(M,N,DW(2*(M+N)+1),LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77KD(M,N,A,LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 50 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+N+I)**DP) )
   50       CONTINUE
            RINFO(2) = ZERO
            DO 55 J = 1,N
              DW(M+J) = DW(M+J)**DP
              IF (IW(M+J).NE.0)
     &           RINFO(2) = MAX( RINFO(2), ABS(ONE-DW(2*M+N+J)**DP) )
   55       CONTINUE
          ENDIF
        ENDIF
      ELSE
        IF (JOB.EQ.0) THEN
          IF (SETUP.EQ.1) THEN
            CALL MC77LD(M,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77LD(M,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77MD(M,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77MD(M,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 60 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+I)**DP) )
   60       CONTINUE
          ENDIF
        ENDIF
        RINFO(2) = RINFO(1)
      ENDIF
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,3),(RINFO(J),J=1,2)
        IF (ICNTL(9).LT.2) THEN
          WRITE(ICNTL(3),9031) (DW(I),I=1,M)
          IF (ICNTL(6).EQ.0)
     &      WRITE(ICNTL(3),9032) (DW(M+J),J=1,N)
        ENDIF
      ENDIF
   99 CONTINUE
      RETURN
 9001 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3,
     &        ' because ',(A),' = ',I10)
 9002 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Bad input control flag.',
     &        ' Value of ICNTL(',I1,') = ',I8)
 9003 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Input matrix is symmetric and N /= M.',
     &        ' Value of (N-M) = ',I8)
 9004 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        LDW too small, must be at least ',I8)
 9008 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Bad input REAL control parameter.',
     &        ' Value of CNTL(',I1,') = ',1PD14.4)
 9009 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Only scaling in infinity norm is allowed',
     &        ' for rectangular matrices'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8)
C9010 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
Cnew &        '        Algorithm will loop indefinitely
Cnew &        '     Check for inconsistency'/
Cnew &        ' Maximum number of iterations = ',I8/
Cnew &        ' Threshold for convergence    = ',1PD14.4)
 9020 FORMAT (' ****** Input parameters for MC77C/CD:'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8/
     &        ' Max n.b. Iters.  = ',I8/
     &        ' Cvgce. Threshold = ',1PD14.4)
 9021 FORMAT (' ROW     (I) = ',I8/
     &        ' A(I,1:N)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9022 FORMAT (' ROW     (I) = ',I8/
     &        ' A(I,I:N)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC77C/CD:'/
     &        ' INFO(1:3)   = ',3(2X,I8)/
     &        ' RINFO(1:2)  = ',2(1PD14.4))
 9031 FORMAT (' DW(1:M)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9032 FORMAT (' DW(M+1:M+N) = ',4(1PD14.4)/(15X,4(1PD14.4)))
c9031 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
c9032 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END
C**********************************************************************
      SUBROUTINE MC77JD(M,N,A,LDA,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
      INTEGER M,N,LDA,MAXIT,NITER,CHECK,INFO
      INTEGER IW(M),JW(N)
      DOUBLE PRECISION THRESH,ERR(2)
      DOUBLE PRECISION A(LDA,N),D(M),E(N),DW(M),EW(N)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J
      DOUBLE PRECISION S
      INTRINSIC SQRT, MAX, ABS
      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO
      DO 5 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
    5 CONTINUE
      DO 10 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   10 CONTINUE
      DO 20 J=1,N
        DO 30 I=1,M
          IF (EW(J).LT.A(I,J)) THEN
            EW(J) = A(I,J)
            JW(J) = I
          ENDIF
          IF (DW(I).LT.A(I,J)) THEN
            DW(I) = A(I,J)
            IW(I) = J
          ENDIF
   30   CONTINUE
   20 CONTINUE
      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE
      DO 50 J=1,N
        I = JW(J)
        IF (I.GT.0) THEN
          IF (IW(I).EQ.J) THEN
            IW(I) = -IW(I)
            JW(J) = -JW(J)
          ENDIF
        ENDIF
   50 CONTINUE
      DO 60 I=1,M
        IF (IW(I).GT.0)  GOTO 99
   60 CONTINUE
      DO 65 J=1,N
        IF (JW(J).GT.0)  GOTO 99
   65 CONTINUE
      GOTO 200
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF
        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE
        DO 120 J=1,N
          IF (JW(J).GT.0) THEN
            DO 130 I=1,M
              S = A(I,J) / (D(I)*E(J))
              IF (EW(J).LT.S) THEN
                EW(J) = S
                JW(J) = I
              ENDIF
              IF (IW(I).GT.0) THEN
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  130       CONTINUE
          ELSE
            DO 140 I=1,M
              IF (IW(I).GT.0) THEN
                S = A(I,J) / (D(I)*E(J))
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  140       CONTINUE
          ENDIF
  120   CONTINUE
        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE
        DO 160 J=1,N
          I = JW(J)
          IF (I.GT.0) THEN
            IF (IW(I).EQ.J) THEN
              IW(I) = -IW(I)
              JW(J) = -JW(J)
            ENDIF
          ENDIF
  160   CONTINUE
        IF (CHECK.LE.0) GOTO 99
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) )  GOTO 200
      IF (CHECK.GT.0)  GOTO 99
  200 RETURN
      END
C**********************************************************************
      SUBROUTINE MC77KD(M,N,A,LDA,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
      INTEGER M,N,LDA,MAXIT,NITER,CHECK,INFO
      INTEGER IW(M),JW(N)
      DOUBLE PRECISION THRESH,ERR(2)
      DOUBLE PRECISION A(LDA,N),D(M),E(N),DW(M),EW(N)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
      INTRINSIC SQRT, MAX, ABS
      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO
      DO 10 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
   10 CONTINUE
      DO 15 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   15 CONTINUE
      DO 20 J=1,N
        DO 30 I=1,M
          IF (A(I,J).GT.ZERO) THEN
            DW(I) = DW(I) + A(I,J)
            IW(I) = IW(I) + 1
            EW(J) = EW(J) + A(I,J)
            JW(J) = JW(J) + 1
          ENDIF
   30   CONTINUE
   20 CONTINUE
      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE
      DO 60 K=1,M
        IF ( (IW(K).GT.0) .OR. (JW(K).GT.0) )  GOTO 99
   60 CONTINUE
Cnew  DO 60 I=1,M
Cnew    IF (IW(I).GT.0)  GOTO 99
Cne60 CONTINUE
Cnew  DO 65 J=1,N
Cnew    IF (JW(J).GT.0)  GOTO 99
Cne65 CONTINUE
      GOTO 200
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF
        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE
        DO 120 J=1,N
          IF (JW(J).GT.0) THEN
            DO 130 I=1,M
              S = A(I,J) / (D(I)*E(J))
              EW(J) = EW(J) + S
              DW(I) = DW(I) + S
  130       CONTINUE
          ENDIF
  120   CONTINUE
        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE
        IF (CHECK.LE.0) GOTO 99
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) ) GOTO 200
      IF (CHECK.GT.0)  GOTO 99
  200 RETURN
      END
C**********************************************************************
      SUBROUTINE MC77LD(N,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
      INTEGER N,MAXIT,NITER,CHECK,INFO
      INTEGER IJW(N)
      DOUBLE PRECISION THRESH,ERR
      DOUBLE PRECISION A(*),DE(N),DEW(N)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
      INTRINSIC SQRT, MAX, ABS
      INFO = 0
      NITER = 0
      ERR = ZERO
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE
      K = 1
      DO 20 J=1,N
        DO 30 I=J,N
          IF (DEW(J).LT.A(K)) THEN
            DEW(J) = A(K)
            IJW(J) = I
          ENDIF
          IF (DEW(I).LT.A(K)) THEN
            DEW(I) = A(K)
            IJW(I) = J
          ENDIF
          K = K+1
   30   CONTINUE
   20 CONTINUE
      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE
      DO 50 J=1,N
        I = IJW(J)
        IF (I.GT.0) THEN
          IF (IJW(I).EQ.J) THEN
            IJW(I) = -IJW(I)
            IF (I.NE.J)  IJW(J) = -IJW(J)
          ENDIF
        ENDIF
   50 CONTINUE
      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF
        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE
        K = 1
        DO 120 J=1,N
          IF (IJW(J).GT.0) THEN
            DO 130 I=J,N
              S = A(K) / (DE(I)*DE(J))
              IF (DEW(J).LT.S) THEN
                DEW(J) = S
                IJW(J) = I
              ENDIF
              IF (IJW(I).GT.0) THEN
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
              K = K+1
  130       CONTINUE
          ELSE
            DO 140 I=J,N
              IF (IJW(I).GT.0) THEN
                S = A(K) / (DE(I)*DE(J))
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
              K = K+1
  140       CONTINUE
          ENDIF
  120   CONTINUE
        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE
        DO 160 J=1,N
          I = IJW(J)
          IF (I.GT.0) THEN
            IF (IJW(I).EQ.J) THEN
              IJW(I) = -IJW(I)
              IF (I.NE.J)  IJW(J) = -IJW(J)
            ENDIF
          ENDIF
  160   CONTINUE
        IF (CHECK.LE.0) GOTO 99
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200
      IF (CHECK.GT.0)  GOTO 99
  200 RETURN
      END
C**********************************************************************
      SUBROUTINE MC77MD(N,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
      INTEGER N,MAXIT,NITER,CHECK,INFO
      INTEGER IJW(N)
      DOUBLE PRECISION THRESH,ERR
      DOUBLE PRECISION A(*),DE(N),DEW(N)
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE=1.0D+00, ZERO=0.0D+00 )
      INTEGER I, J, K
      DOUBLE PRECISION S
      INTRINSIC SQRT, MAX, ABS
      INFO = 0
      NITER = 0
      ERR = ZERO
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE
      K = 1
      DO 20 J=1,N
        DO 30 I=J,N
          IF (A(K).GT.ZERO) THEN
            DEW(J) = DEW(J) + A(K)
            IJW(J) = IJW(J) + 1
            IF (I.NE.J) THEN
              DEW(I) = DEW(I) + A(K)
              IJW(I) = IJW(I) + 1
            ENDIF
          ENDIF
          K = K+1
   30   CONTINUE
   20 CONTINUE
      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE
      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF
        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE
        K = 1
        DO 120 J=1,N
          IF (IJW(J).GT.0) THEN
            DO 130 I=J,N
              S = A(K) / (DE(I)*DE(J))
              DEW(J) = DEW(J) + S
              IF (I.NE.J)  DEW(I) = DEW(I) + S
              K = K+1
  130       CONTINUE
          ENDIF
  120   CONTINUE
        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE
        IF (CHECK.LE.0) GOTO 99
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200
      IF (CHECK.GT.0)  GOTO 99
  200 RETURN
      END
